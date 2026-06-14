from __future__ import annotations

import os
import platform
import random
import ssl
import time
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional
import hashlib
import json
import socket
import Bio
from Bio import Entrez
import xml.etree.ElementTree as ET

from src.pago_pipeline.storage import (
    read_json_file,
    sha256_of_file,
    sha256_of_lines,
    write_bytes_atomic,
    write_json_atomic,
)


NCBI_SSL_CERT_FILE_ENV_VARS = (
    "NCBI_SSL_CERT_FILE",
    "SSL_CERT_FILE",
)
NCBI_SSL_CERT_DIR_ENV_VARS = (
    "NCBI_SSL_CERT_DIR",
    "SSL_CERT_DIR",
)
_active_ncbi_ssl_context: ContextVar[Optional[ssl.SSLContext]] = ContextVar(
    "_active_ncbi_ssl_context",
    default=None,
)
_ncbi_entrez_delegate_urlopen = getattr(Entrez, "urlopen", None)
DEFAULT_XML_SOCKET_IDLE_TIMEOUT_SECONDS = 90.0
DEFAULT_XML_BATCH_DEADLINE_SECONDS = 300.0
DEFAULT_XML_RETRY_DELAY_SECONDS = (0.0, 5.0, 15.0, 45.0, 120.0)
DEFAULT_XML_CIRCUIT_BREAKER_FAILURE_THRESHOLD = 3
DEFAULT_XML_CIRCUIT_BREAKER_COOLDOWN_SECONDS = 180.0
DEFAULT_XML_CHUNK_SIZE_BYTES = 64 * 1024


class XMLBatchReadError(RuntimeError):
    def __init__(self, message: str, *, bytes_received: int) -> None:
        super().__init__(message)
        self.bytes_received = bytes_received


class XMLBatchReadTimeoutError(TimeoutError):
    def __init__(self, message: str, *, bytes_received: int) -> None:
        super().__init__(message)
        self.bytes_received = bytes_received


def _first_configured_environment_variable(
    variable_names: tuple[str, ...],
) -> Optional[str]:
    """
    Return the first non-empty configured environment variable value.
    """
    for variable_name in variable_names:
        variable_value = os.getenv(variable_name)
        if variable_value and variable_value.strip():
            return variable_value.strip()

    return None


def _resolve_ncbi_ssl_ca_configuration(
    *,
    ssl_ca_file: Optional[str],
    ssl_ca_directory: Optional[str],
) -> tuple[Optional[str], Optional[str]]:
    """
    Resolve explicit and environment-provided CA configuration for NCBI HTTPS calls.
    """
    resolved_ssl_ca_file = ssl_ca_file
    if resolved_ssl_ca_file is None:
        resolved_ssl_ca_file = _first_configured_environment_variable(
            NCBI_SSL_CERT_FILE_ENV_VARS,
        )

    resolved_ssl_ca_directory = ssl_ca_directory
    if resolved_ssl_ca_directory is None:
        resolved_ssl_ca_directory = _first_configured_environment_variable(
            NCBI_SSL_CERT_DIR_ENV_VARS,
        )

    return resolved_ssl_ca_file, resolved_ssl_ca_directory


def _validate_ncbi_ssl_ca_configuration(
    *,
    ssl_ca_file: Optional[str],
    ssl_ca_directory: Optional[str],
) -> None:
    """
    Validate CA bundle paths before attempting an HTTPS request.
    """
    if ssl_ca_file:
        resolved_ssl_ca_file_path = Path(ssl_ca_file).expanduser()
        if not resolved_ssl_ca_file_path.exists():
            raise ValueError(
                "Invalid NCBI CA configuration: ssl_ca_file does not exist: "
                f"{resolved_ssl_ca_file_path}. Provide a readable CA bundle "
                "file with ssl_ca_file or set NCBI_SSL_CERT_FILE."
            )

        if not resolved_ssl_ca_file_path.is_file():
            raise ValueError(
                "Invalid NCBI CA configuration: ssl_ca_file is not a file: "
                f"{resolved_ssl_ca_file_path}. Provide a readable CA bundle "
                "file with ssl_ca_file or set NCBI_SSL_CERT_FILE."
            )

    if ssl_ca_directory:
        resolved_ssl_ca_directory_path = Path(ssl_ca_directory).expanduser()
        if not resolved_ssl_ca_directory_path.exists():
            raise ValueError(
                "Invalid NCBI CA configuration: ssl_ca_directory does not exist: "
                f"{resolved_ssl_ca_directory_path}. Provide a valid CA "
                "directory with ssl_ca_directory or set NCBI_SSL_CERT_DIR."
            )

        if not resolved_ssl_ca_directory_path.is_dir():
            raise ValueError(
                "Invalid NCBI CA configuration: ssl_ca_directory is not a directory: "
                f"{resolved_ssl_ca_directory_path}. Provide a valid CA "
                "directory with ssl_ca_directory or set NCBI_SSL_CERT_DIR."
            )


def _build_ncbi_tls_configuration_error_message(
    *,
    error: Exception,
    ssl_ca_file: Optional[str],
    ssl_ca_directory: Optional[str],
) -> str:
    """
    Build a clear guidance message for TLS trust failures against NCBI.
    """
    configured_inputs: List[str] = []

    if ssl_ca_file:
        configured_inputs.append(f"ssl_ca_file={ssl_ca_file}")
    if ssl_ca_directory:
        configured_inputs.append(f"ssl_ca_directory={ssl_ca_directory}")

    configured_inputs_summary = (
        ", ".join(configured_inputs)
        if configured_inputs
        else "no explicit CA path was configured"
    )

    return (
        "NCBI HTTPS request failed during TLS certificate verification. "
        "Configure a trusted corporate/root CA bundle with ssl_ca_file and/or "
        "ssl_ca_directory, or set NCBI_SSL_CERT_FILE / NCBI_SSL_CERT_DIR. "
        "TLS verification remains enabled and must not be disabled. "
        f"Resolved CA configuration: {configured_inputs_summary}. "
        f"Original error: {error}"
    )


def _ncbi_entrez_urlopen_with_context(*args, **kwargs):
    """
    Delegate Entrez HTTPS calls while injecting the active SSL context for this flow.
    """
    active_ssl_context = _active_ncbi_ssl_context.get()
    if active_ssl_context is not None and "context" not in kwargs:
        kwargs["context"] = active_ssl_context

    if _ncbi_entrez_delegate_urlopen is None:
        raise RuntimeError(
            "Bio.Entrez.urlopen is unavailable; unable to configure TLS for NCBI requests."
        )

    return _ncbi_entrez_delegate_urlopen(*args, **kwargs)


def _ensure_ncbi_entrez_urlopen_hook_installed() -> None:
    """
    Install a stable Entrez urlopen hook once and route per-call SSL via ContextVar.
    """
    global _ncbi_entrez_delegate_urlopen

    current_urlopen = getattr(Entrez, "urlopen", None)
    if current_urlopen is _ncbi_entrez_urlopen_with_context:
        return

    if current_urlopen is None:
        raise RuntimeError(
            "Bio.Entrez.urlopen is unavailable; unable to configure TLS for NCBI requests."
        )

    _ncbi_entrez_delegate_urlopen = current_urlopen
    Entrez.urlopen = _ncbi_entrez_urlopen_with_context


def _is_ncbi_tls_configuration_error(error: Exception) -> bool:
    """
    Return True when an exception indicates TLS trust-store misconfiguration.
    """
    pending_errors: List[BaseException] = [error]
    visited_error_ids: set[int] = set()

    while pending_errors:
        current_error = pending_errors.pop()
        current_error_id = id(current_error)
        if current_error_id in visited_error_ids:
            continue
        visited_error_ids.add(current_error_id)

        if isinstance(current_error, ssl.SSLError):
            return True

        current_error_message = str(current_error)
        if (
            "CERTIFICATE_VERIFY_FAILED" in current_error_message
            or "certificate verify failed" in current_error_message.lower()
            or "self-signed certificate" in current_error_message.lower()
            or "unable to get local issuer certificate" in current_error_message.lower()
        ):
            return True

        nested_reason = getattr(current_error, "reason", None)
        if isinstance(nested_reason, BaseException):
            pending_errors.append(nested_reason)

        nested_cause = getattr(current_error, "__cause__", None)
        if isinstance(nested_cause, BaseException):
            pending_errors.append(nested_cause)

        nested_context = getattr(current_error, "__context__", None)
        if isinstance(nested_context, BaseException):
            pending_errors.append(nested_context)

    return False


@contextmanager
def _configured_ncbi_entrez_urlopen(
    *,
    ssl_ca_file: Optional[str],
    ssl_ca_directory: Optional[str],
):
    """
    Temporarily configure Bio.Entrez to use a custom SSL trust store.
    """
    resolved_ssl_ca_file, resolved_ssl_ca_directory = (
        _resolve_ncbi_ssl_ca_configuration(
            ssl_ca_file=ssl_ca_file,
            ssl_ca_directory=ssl_ca_directory,
        )
    )

    _validate_ncbi_ssl_ca_configuration(
        ssl_ca_file=resolved_ssl_ca_file,
        ssl_ca_directory=resolved_ssl_ca_directory,
    )

    if not resolved_ssl_ca_file and not resolved_ssl_ca_directory:
        yield resolved_ssl_ca_file, resolved_ssl_ca_directory
        return

    try:
        ssl_context = ssl.create_default_context(
            cafile=resolved_ssl_ca_file,
            capath=resolved_ssl_ca_directory,
        )
    except Exception as error:
        raise ValueError(
            "Invalid NCBI CA configuration. Provide a readable CA bundle with "
            "ssl_ca_file, a valid CA directory with ssl_ca_directory, or set "
            "NCBI_SSL_CERT_FILE / NCBI_SSL_CERT_DIR. "
            f"Original error: {error}"
        ) from error

    _ensure_ncbi_entrez_urlopen_hook_installed()
    ssl_context_token = _active_ncbi_ssl_context.set(ssl_context)

    try:
        yield resolved_ssl_ca_file, resolved_ssl_ca_directory
    finally:
        _active_ncbi_ssl_context.reset(ssl_context_token)


@dataclass(frozen=True)
class NCBIProteinUidFetchResult:
    """
    Immutable payload representing one NCBI protein-UID retrieval event.

    This object is designed to become the basis of a reproducible local snapshot.
    """
    database_name: str
    search_query: str
    translated_query: Optional[str]
    identifier_type: str
    retrieved_at_utc: str
    ncbi_reported_result_count: int
    protein_uids: List[str]
    raw_protein_uid_count: int
    normalized_protein_uid_count: int
    deduplicate_uids: bool
    sort_uids: bool
    protein_uids_sha256: str
    page_size: int
    request_delay_seconds: float
    max_retry_attempts: int
    history_web_env: Optional[str]
    history_query_key: Optional[str]
    python_version: str
    biopython_version: str
    source_type: Optional[str] = None
    source_file_path: Optional[str] = None
    source_file_sha256: Optional[str] = None
    source_sheet_name: Optional[str] = None
    source_identifier_column: Optional[str] = None


def _normalize_protein_uid_list(
    *,
    protein_uids: List[str],
    deduplicate_uids: bool,
    sort_uids: bool,
) -> List[str]:
    """
    Normalize protein UIDs by:
    - casting to str
    - stripping surrounding whitespace
    - removing empty values
    - optionally deduplicating while preserving first occurrence order
    - optionally sorting canonically
    """
    normalized_protein_uids = [
        str(protein_uid).strip()
        for protein_uid in protein_uids
        if str(protein_uid).strip()
    ]

    if deduplicate_uids:
        unique_protein_uids: List[str] = []
        seen_protein_uids = set()

        for protein_uid in normalized_protein_uids:
            if protein_uid not in seen_protein_uids:
                seen_protein_uids.add(protein_uid)
                unique_protein_uids.append(protein_uid)

        normalized_protein_uids = unique_protein_uids

    if sort_uids:
        normalized_protein_uids = sorted(normalized_protein_uids)

    return normalized_protein_uids


def fetch_ncbi_protein_uid_snapshot(
    *,
    ncbi_email: str,
    ncbi_api_key: Optional[str],
    query: str,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    deduplicate_uids: bool = True,
    sort_uids: bool = True,
    page_size: int = 1000,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
) -> NCBIProteinUidFetchResult:
    """
    Retrieve protein UIDs from NCBI Entrez and return a snapshot-ready payload.

    Important:
    - This function does NOT provide long-term reproducibility by itself.
    - Reproducibility comes from persisting the returned payload locally
      as an immutable snapshot artifact.

    Args:
        ncbi_email:
            Required by NCBI usage policies / best practice.
        ncbi_api_key:
            Optional NCBI API key. If provided, a higher request rate is allowed.
        query:
            Entrez search term.
        ssl_ca_file:
            Optional CA bundle file used to validate HTTPS connections to NCBI.
            If omitted, NCBI_SSL_CERT_FILE and then SSL_CERT_FILE are checked.
        ssl_ca_directory:
            Optional CA directory used to validate HTTPS connections to NCBI.
            If omitted, NCBI_SSL_CERT_DIR and then SSL_CERT_DIR are checked.
        deduplicate_uids:
            Whether to deduplicate UIDs in the returned payload.
        sort_uids:
            Whether to sort UIDs canonically in the returned payload.
        page_size:
            Number of protein UIDs requested per paginated Entrez call.
        max_retry_attempts:
            Maximum retry attempts for each failed paginated request.
        request_delay_seconds:
            Optional custom delay between successful requests.
            If None, defaults to:
            - 0.10 seconds when an API key is provided
            - 0.34 seconds otherwise

    Returns:
        NCBIProteinUidFetchResult:
            Snapshot-ready result including UIDs and retrieval metadata.
    """
    if not ncbi_email:
        raise ValueError(
            "No NCBI email credential found. Please add .env file with "
            "NCBI_EMAIL. NCBI_API_KEY is optional."
        )

    if not query or not query.strip():
        raise ValueError("Query must be a non-empty string.")

    if page_size <= 0:
        raise ValueError("page_size must be a positive integer.")

    if max_retry_attempts <= 0:
        raise ValueError("max_retry_attempts must be a positive integer.")

    Entrez.email = ncbi_email
    Entrez.api_key = ncbi_api_key

    database_name = "protein"
    identifier_type = "uid"

    resolved_request_delay_seconds: float
    if request_delay_seconds is None:
        resolved_request_delay_seconds = 0.10 if ncbi_api_key else 0.34
    else:
        resolved_request_delay_seconds = request_delay_seconds

    retrieved_at_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    with _configured_ncbi_entrez_urlopen(
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
    ) as (resolved_ssl_ca_file, resolved_ssl_ca_directory):
        try:
            initial_search_handle = Entrez.esearch(
                db=database_name,
                term=query,
                retmax=0,
                usehistory="y",
            )
            try:
                initial_search_response = Entrez.read(initial_search_handle)
            finally:
                initial_search_handle.close()
        except Exception as error:
            if _is_ncbi_tls_configuration_error(error):
                raise RuntimeError(
                    _build_ncbi_tls_configuration_error_message(
                        error=error,
                        ssl_ca_file=resolved_ssl_ca_file,
                        ssl_ca_directory=resolved_ssl_ca_directory,
                    )
                ) from error

            raise

        ncbi_reported_result_count = int(initial_search_response["Count"])
        translated_query = initial_search_response.get("QueryTranslation")
        history_web_env = initial_search_response.get("WebEnv")
        history_query_key = initial_search_response.get("QueryKey")

        print(f"Found {ncbi_reported_result_count} protein UIDs.")
        print(
            f"History session: WebEnv={history_web_env}, "
            f"QueryKey={history_query_key}."
        )

        if ncbi_reported_result_count == 0:
            print("NCBI returned zero protein UIDs for the provided query.")

        if ncbi_reported_result_count > 0 and (
            not history_web_env or not history_query_key
        ):
            raise RuntimeError(
                "NCBI search history metadata is required to fetch "
                "protein UIDs, but WebEnv/query_key was missing "
                "from the ESearch response."
            )

        raw_protein_uids: List[str] = []
        page_start_index = 0

        while page_start_index < ncbi_reported_result_count:
            for retry_attempt_index in range(max_retry_attempts):
                try:
                    paginated_search_handle = Entrez.esearch(
                        db=database_name,
                        term=query,
                        usehistory="y",
                        retmax=page_size,
                        retstart=page_start_index,
                    )
                    try:
                        paginated_search_response = Entrez.read(
                            paginated_search_handle
                        )
                    finally:
                        paginated_search_handle.close()

                    if "IdList" not in paginated_search_response:
                        raise RuntimeError("NCBI response missing 'IdList'.")

                    current_batch_protein_uids = paginated_search_response["IdList"]

                    if not current_batch_protein_uids:
                        raise RuntimeError(
                            "NCBI returned an empty IdList before the reported "
                            "result set was fully paginated."
                        )

                    raw_protein_uids.extend(current_batch_protein_uids)
                    page_start_index += len(current_batch_protein_uids)

                    retrieval_progress = min(
                        page_start_index / ncbi_reported_result_count,
                        1.0,
                    )
                    print(
                        f"Extracted {len(raw_protein_uids)} protein UIDs "
                        f"({retrieval_progress:.2%})."
                    )

                    jitter_seconds = random.uniform(0.01, 0.05)
                    time.sleep(resolved_request_delay_seconds + jitter_seconds)
                    break

                except Exception as error:
                    if _is_ncbi_tls_configuration_error(error):
                        raise RuntimeError(
                            _build_ncbi_tls_configuration_error_message(
                                error=error,
                                ssl_ca_file=resolved_ssl_ca_file,
                                ssl_ca_directory=resolved_ssl_ca_directory,
                            )
                        ) from error

                    if retry_attempt_index == max_retry_attempts - 1:
                        raise RuntimeError(
                            f"Failed to extract protein UIDs after "
                            f"{max_retry_attempts} attempts at "
                            f"page_start_index={page_start_index}: {error}"
                        ) from error

                    retry_backoff_seconds = (
                        2 ** retry_attempt_index
                    ) * resolved_request_delay_seconds

                    print(
                        f"Retry {retry_attempt_index + 1}/{max_retry_attempts} "
                        f"after error: {error}"
                    )
                    time.sleep(retry_backoff_seconds)

        normalized_protein_uids = _normalize_protein_uid_list(
            protein_uids=raw_protein_uids,
            deduplicate_uids=deduplicate_uids,
            sort_uids=sort_uids,
        )

        protein_uids_sha256 = sha256_of_lines(
            text_lines=normalized_protein_uids,
            deduplicate_lines_preserving_order=False,
            sort_lines=False,
        )

        fetch_result = NCBIProteinUidFetchResult(
            database_name=database_name,
            search_query=query,
            translated_query=translated_query,
            identifier_type=identifier_type,
            retrieved_at_utc=retrieved_at_utc,
            ncbi_reported_result_count=ncbi_reported_result_count,
            protein_uids=normalized_protein_uids,
            raw_protein_uid_count=len(raw_protein_uids),
            normalized_protein_uid_count=len(normalized_protein_uids),
            deduplicate_uids=deduplicate_uids,
            sort_uids=sort_uids,
            protein_uids_sha256=protein_uids_sha256,
            page_size=page_size,
            request_delay_seconds=resolved_request_delay_seconds,
            max_retry_attempts=max_retry_attempts,
            history_web_env=history_web_env,
            history_query_key=history_query_key,
            python_version=platform.python_version(),
            biopython_version=getattr(Bio, "__version__", "unknown"),
        )

        print(f"Final raw UID count: {fetch_result.raw_protein_uid_count}")
        print(
            "Final normalized UID count: "
            f"{fetch_result.normalized_protein_uid_count}"
        )
        print(f"Dataset SHA-256: {fetch_result.protein_uids_sha256}")

        return fetch_result


@dataclass(frozen=True)
class NCBIProteinXmlBatchFetchResult:
    """
    Immutable payload representing one XML batch fetched from NCBI Entrez.
    """
    batch_index: int
    batch_start_index: int
    batch_end_index: int
    protein_uids: List[str]
    protein_uid_count: int
    xml_payload_bytes: bytes
    xml_payload_sha256: str
    xml_payload_file_path: Optional[str] = None
    open_latency_seconds: Optional[float] = None
    read_latency_seconds: Optional[float] = None
    total_latency_seconds: Optional[float] = None
    response_byte_count: Optional[int] = None
    average_throughput_bytes_per_second: Optional[float] = None
    attempt_count: Optional[int] = None


@dataclass(frozen=True)
class NCBIProteinXmlFetchResult:
    """
    Immutable payload representing one NCBI protein-XML retrieval event.

    Important:
    This object describes the API retrieval step only.
    It does not persist any local snapshot by itself.
    """
    database_name: str
    identifier_type: str
    retrieved_at_utc: str
    requested_protein_uid_count: int
    normalized_protein_uid_count: int
    protein_uids_sha256: str
    batch_size: int
    batch_count: int
    retmode: str
    request_delay_seconds: float
    max_retry_attempts: int
    python_version: str
    biopython_version: str
    xml_batches: List[NCBIProteinXmlBatchFetchResult]
    batch_workspace_directory: Optional[str] = None
    socket_idle_timeout_seconds: Optional[float] = None
    batch_deadline_seconds: Optional[float] = None
    circuit_breaker_failure_threshold: Optional[int] = None
    circuit_breaker_cooldown_seconds: Optional[float] = None


def _coerce_payload_to_bytes(
    *,
    payload: str | bytes,
    encoding: str = "utf-8",
) -> bytes:
    """
    Normalize a response payload into raw bytes.
    """
    if isinstance(payload, bytes):
        return payload

    return payload.encode(encoding)


def _utc_now_string() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _as_path(path_like: str | Path) -> Path:
    return Path(path_like)


def _build_xml_batch_workspace_id(
    *,
    protein_uids_sha256: str,
    batch_size: int,
    retmode: str,
) -> str:
    return f"uids_{protein_uids_sha256[:16]}__batch_{batch_size}__{retmode}"


def _build_xml_batch_plan(
    *,
    normalized_protein_uids: List[str],
    protein_uids_sha256: str,
    batch_size: int,
    retmode: str,
) -> dict:
    batches: list[dict[str, object]] = []
    total_batch_count = (
        len(normalized_protein_uids) + batch_size - 1
    ) // batch_size

    for batch_index, batch_start_index in enumerate(
        range(0, len(normalized_protein_uids), batch_size),
        start=1,
    ):
        current_batch_protein_uids = normalized_protein_uids[
            batch_start_index: batch_start_index + batch_size
        ]
        current_batch_end_index = (
            batch_start_index + len(current_batch_protein_uids) - 1
        )
        batches.append(
            {
                "batch_index": batch_index,
                "batch_start_index": batch_start_index,
                "batch_end_index": current_batch_end_index,
                "protein_uid_count": len(current_batch_protein_uids),
                "protein_uids_sha256": sha256_of_lines(
                    text_lines=current_batch_protein_uids,
                    deduplicate_lines_preserving_order=False,
                    sort_lines=False,
                ),
            }
        )

    return {
        "plan_format_version": "1.0",
        "database_name": "protein",
        "identifier_type": "uid",
        "retmode": retmode,
        "batch_size": batch_size,
        "batch_count": total_batch_count,
        "normalized_protein_uid_count": len(normalized_protein_uids),
        "protein_uids_sha256": protein_uids_sha256,
        "batches": batches,
    }


def _write_or_validate_xml_batch_plan(
    *,
    workspace_directory: Path,
    batch_plan: dict,
) -> None:
    plan_file_path = workspace_directory / "batch_plan.json"
    workspace_directory.mkdir(parents=True, exist_ok=True)

    if not plan_file_path.exists():
        write_json_atomic(payload=batch_plan, output_file_path=plan_file_path)
        return

    saved_plan = read_json_file(input_file_path=plan_file_path)
    comparable_keys = (
        "retmode",
        "batch_size",
        "batch_count",
        "normalized_protein_uid_count",
        "protein_uids_sha256",
    )
    for comparable_key in comparable_keys:
        if saved_plan.get(comparable_key) != batch_plan.get(comparable_key):
            raise RuntimeError(
                "Existing XML batch workspace plan does not match the current "
                f"request. Mismatch key: {comparable_key}."
            )


def _xml_batch_file_paths(
    *,
    workspace_directory: Path,
    batch_index: int,
) -> tuple[Path, Path]:
    batches_directory = workspace_directory / "batches"
    batch_file_stem = f"batch_{batch_index:06d}"
    return (
        batches_directory / f"{batch_file_stem}.xml",
        batches_directory / f"{batch_file_stem}.json",
    )


def _count_xml_batch_records(
    *,
    xml_payload_bytes: bytes,
    expected_root_tag: str = "GBSet",
    expected_record_tag: str = "GBSeq",
) -> int:
    try:
        root_element = ET.fromstring(xml_payload_bytes)
    except ET.ParseError as error:
        raise RuntimeError(f"XML batch is not well-formed: {error}") from error

    if root_element.tag != expected_root_tag:
        raise RuntimeError(
            "XML batch root tag mismatch. "
            f"Expected {expected_root_tag}, got {root_element.tag}."
        )

    child_elements = list(root_element)
    invalid_child_tags = sorted(
        {child.tag for child in child_elements if child.tag != expected_record_tag}
    )
    if invalid_child_tags:
        raise RuntimeError(
            "XML batch contains unexpected child tags under "
            f"{expected_root_tag}: {invalid_child_tags}."
        )

    return len(child_elements)


def _validate_persisted_xml_batch(
    *,
    batch_xml_file_path: Path,
    batch_metadata_file_path: Path,
    expected_batch_index: int,
    expected_batch_start_index: int,
    expected_batch_end_index: int,
    expected_protein_uids: List[str],
) -> Optional[dict]:
    if not batch_xml_file_path.exists() or not batch_metadata_file_path.exists():
        return None

    try:
        batch_metadata = read_json_file(input_file_path=batch_metadata_file_path)
    except Exception:
        return None

    expected_uid_sha256 = sha256_of_lines(
        text_lines=expected_protein_uids,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )
    if (
        batch_metadata.get("batch_index") != expected_batch_index
        or batch_metadata.get("batch_start_index") != expected_batch_start_index
        or batch_metadata.get("batch_end_index") != expected_batch_end_index
        or batch_metadata.get("protein_uid_count") != len(expected_protein_uids)
        or batch_metadata.get("protein_uids_sha256") != expected_uid_sha256
    ):
        return None

    actual_xml_sha256 = sha256_of_file(input_file_path=batch_xml_file_path)
    if batch_metadata.get("xml_payload_sha256") != actual_xml_sha256:
        return None

    try:
        xml_payload_bytes = batch_xml_file_path.read_bytes()
        record_count = _count_xml_batch_records(xml_payload_bytes=xml_payload_bytes)
    except Exception:
        return None

    if record_count != len(expected_protein_uids):
        return None

    return batch_metadata


def _record_xml_latency_event(
    *,
    workspace_directory: Optional[Path],
    event_payload: dict,
) -> None:
    if workspace_directory is None:
        return

    event_payload = {
        **event_payload,
        "recorded_at_utc": _utc_now_string(),
    }
    latency_events_file_path = workspace_directory / "latency_events.jsonl"
    latency_events_file_path.parent.mkdir(parents=True, exist_ok=True)
    with latency_events_file_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(event_payload, sort_keys=True) + "\n")


def _read_response_body_with_deadline(
    *,
    fetch_handle,
    batch_deadline_seconds: float,
    deadline_started_at: float,
    chunk_size_bytes: int = DEFAULT_XML_CHUNK_SIZE_BYTES,
) -> tuple[bytes, float, int]:
    response_chunks: list[bytes] = []
    read_started_at = time.perf_counter()
    bytes_received = 0

    while True:
        elapsed_seconds = time.perf_counter() - deadline_started_at
        if elapsed_seconds > batch_deadline_seconds:
            raise XMLBatchReadTimeoutError(
                "Batch exceeded deadline of "
                f"{batch_deadline_seconds:.1f} seconds after receiving "
                f"{bytes_received} bytes.",
                bytes_received=bytes_received,
            )

        try:
            payload_chunk = fetch_handle.read(chunk_size_bytes)
        except Exception as error:
            raise XMLBatchReadError(
                "Failed while reading XML response body after receiving "
                f"{bytes_received} bytes: {error}",
                bytes_received=bytes_received,
            ) from error
        if not payload_chunk:
            break

        payload_chunk_bytes = _coerce_payload_to_bytes(payload=payload_chunk)
        bytes_received += len(payload_chunk_bytes)
        response_chunks.append(payload_chunk_bytes)

    read_latency_seconds = time.perf_counter() - read_started_at
    return b"".join(response_chunks), read_latency_seconds, bytes_received


def _classify_xml_fetch_error(error: BaseException) -> str:
    error_message = str(error)
    if isinstance(error, TimeoutError) or "timed out" in error_message.lower():
        return "timeout"
    if isinstance(error, ConnectionResetError) or "10054" in error_message:
        return "remote_connection_reset"
    return error.__class__.__name__


def fetch_ncbi_protein_xml_batches(
    *,
    ncbi_email: str,
    ncbi_api_key: Optional[str],
    protein_uids: List[str],
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    batch_size: int = 50,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
    batch_workspace_directory: Optional[str | Path] = None,
    socket_idle_timeout_seconds: float = DEFAULT_XML_SOCKET_IDLE_TIMEOUT_SECONDS,
    batch_deadline_seconds: float = DEFAULT_XML_BATCH_DEADLINE_SECONDS,
    retry_delay_seconds: tuple[float, ...] = DEFAULT_XML_RETRY_DELAY_SECONDS,
    circuit_breaker_failure_threshold: int = (
        DEFAULT_XML_CIRCUIT_BREAKER_FAILURE_THRESHOLD
    ),
    circuit_breaker_cooldown_seconds: float = (
        DEFAULT_XML_CIRCUIT_BREAKER_COOLDOWN_SECONDS
    ),
) -> NCBIProteinXmlFetchResult:
    """
    Fetch raw protein XML batches from NCBI Entrez using a frozen list of protein UIDs.

    This function is intended to run downstream of a previously frozen UID snapshot.
    When batch_workspace_directory is provided, each batch is persisted
    immediately and reusable on subsequent runs.

    Args:
        ncbi_email:
            Required by NCBI usage policies / best practice.
        ncbi_api_key:
            Optional NCBI API key. If provided, a higher request rate is allowed.
        protein_uids:
            Frozen list of protein UIDs that should be fetched as XML.
        ssl_ca_file:
            Optional CA bundle file used to validate HTTPS connections to NCBI.
            If omitted, NCBI_SSL_CERT_FILE and then SSL_CERT_FILE are checked.
        ssl_ca_directory:
            Optional CA directory used to validate HTTPS connections to NCBI.
            If omitted, NCBI_SSL_CERT_DIR and then SSL_CERT_DIR are checked.
        batch_size:
            Number of protein UIDs requested per EFetch call.
        max_retry_attempts:
            Maximum retry attempts for each failed batch request.
        request_delay_seconds:
            Optional custom delay between successful requests.
            If None, defaults to:
            - 0.10 seconds when an API key is provided
            - 0.34 seconds otherwise

    Returns:
        NCBIProteinXmlFetchResult:
            Retrieval metadata plus raw XML payload bytes for each batch.
    """
    if not ncbi_email:
        raise ValueError(
            "No NCBI email credential found. Please add .env file with "
            "NCBI_EMAIL. NCBI_API_KEY is optional."
        )

    if batch_size <= 0:
        raise ValueError("batch_size must be a positive integer.")

    if max_retry_attempts <= 0:
        raise ValueError("max_retry_attempts must be a positive integer.")

    if socket_idle_timeout_seconds <= 0:
        raise ValueError("socket_idle_timeout_seconds must be positive.")

    if batch_deadline_seconds <= 0:
        raise ValueError("batch_deadline_seconds must be positive.")

    if circuit_breaker_failure_threshold <= 0:
        raise ValueError("circuit_breaker_failure_threshold must be positive.")

    if circuit_breaker_cooldown_seconds < 0:
        raise ValueError("circuit_breaker_cooldown_seconds must be non-negative.")

    normalized_protein_uids = _normalize_protein_uid_list(
        protein_uids=protein_uids,
        deduplicate_uids=False,
        sort_uids=False,
    )

    if not normalized_protein_uids:
        raise ValueError("protein_uids must contain at least one non-empty UID.")

    Entrez.email = ncbi_email
    Entrez.api_key = ncbi_api_key
    previous_entrez_max_tries = getattr(Entrez, "max_tries", None)
    Entrez.max_tries = 1

    database_name = "protein"
    identifier_type = "uid"
    retmode = "xml"

    resolved_request_delay_seconds: float
    if request_delay_seconds is None:
        resolved_request_delay_seconds = 0.10 if ncbi_api_key else 0.34
    else:
        resolved_request_delay_seconds = request_delay_seconds

    retrieved_at_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    protein_uids_sha256 = sha256_of_lines(
        text_lines=normalized_protein_uids,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )

    total_batch_count = (
        len(normalized_protein_uids) + batch_size - 1
    ) // batch_size
    batch_plan = _build_xml_batch_plan(
        normalized_protein_uids=normalized_protein_uids,
        protein_uids_sha256=protein_uids_sha256,
        batch_size=batch_size,
        retmode=retmode,
    )

    workspace_directory: Optional[Path] = None
    if batch_workspace_directory is not None:
        workspace_directory = _as_path(batch_workspace_directory)
    elif batch_workspace_directory is None:
        workspace_directory = None

    if workspace_directory is not None:
        _write_or_validate_xml_batch_plan(
            workspace_directory=workspace_directory,
            batch_plan=batch_plan,
        )

    previous_socket_default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(socket_idle_timeout_seconds)

    try:
        with _configured_ncbi_entrez_urlopen(
            ssl_ca_file=ssl_ca_file,
            ssl_ca_directory=ssl_ca_directory,
        ) as (resolved_ssl_ca_file, resolved_ssl_ca_directory):
            xml_batches: List[NCBIProteinXmlBatchFetchResult] = []
            consecutive_degraded_batch_count = 0
            slow_batch_threshold_seconds = min(180.0, batch_deadline_seconds * 0.75)

            for batch_index, batch_start_index in enumerate(
                range(0, len(normalized_protein_uids), batch_size),
                start=1,
            ):
                current_batch_protein_uids = normalized_protein_uids[
                    batch_start_index: batch_start_index + batch_size
                ]
                current_batch_end_index = (
                    batch_start_index + len(current_batch_protein_uids) - 1
                )
                current_batch_uid_sha256 = sha256_of_lines(
                    text_lines=current_batch_protein_uids,
                    deduplicate_lines_preserving_order=False,
                    sort_lines=False,
                )
                batch_xml_file_path: Optional[Path] = None
                batch_metadata_file_path: Optional[Path] = None

                if workspace_directory is not None:
                    batch_xml_file_path, batch_metadata_file_path = (
                        _xml_batch_file_paths(
                            workspace_directory=workspace_directory,
                            batch_index=batch_index,
                        )
                    )
                    batch_metadata = _validate_persisted_xml_batch(
                        batch_xml_file_path=batch_xml_file_path,
                        batch_metadata_file_path=batch_metadata_file_path,
                        expected_batch_index=batch_index,
                        expected_batch_start_index=batch_start_index,
                        expected_batch_end_index=current_batch_end_index,
                        expected_protein_uids=current_batch_protein_uids,
                    )

                    if batch_metadata is not None:
                        print(
                            f"Skipping validated XML batch "
                            f"{batch_index}/{total_batch_count}."
                        )
                        xml_batches.append(
                            NCBIProteinXmlBatchFetchResult(
                                batch_index=batch_index,
                                batch_start_index=batch_start_index,
                                batch_end_index=current_batch_end_index,
                                protein_uids=current_batch_protein_uids,
                                protein_uid_count=len(current_batch_protein_uids),
                                xml_payload_bytes=b"",
                                xml_payload_sha256=batch_metadata[
                                    "xml_payload_sha256"
                                ],
                                xml_payload_file_path=str(batch_xml_file_path),
                                response_byte_count=batch_metadata.get(
                                    "response_byte_count"
                                ),
                                average_throughput_bytes_per_second=(
                                    batch_metadata.get(
                                        "average_throughput_bytes_per_second"
                                    )
                                ),
                                attempt_count=batch_metadata.get("attempt_count"),
                            )
                        )
                        consecutive_degraded_batch_count = 0
                        continue

                if (
                    consecutive_degraded_batch_count
                    >= circuit_breaker_failure_threshold
                ):
                    print(
                        "Circuit breaker cooldown after "
                        f"{consecutive_degraded_batch_count} consecutive "
                        f"degraded XML batches: "
                        f"{circuit_breaker_cooldown_seconds:.1f} seconds."
                    )
                    time.sleep(circuit_breaker_cooldown_seconds)
                    consecutive_degraded_batch_count = 0

                batch_succeeded = False
                for retry_attempt_index in range(max_retry_attempts):
                    attempt_index = retry_attempt_index + 1
                    if retry_attempt_index > 0:
                        base_delay_seconds = retry_delay_seconds[
                            min(retry_attempt_index, len(retry_delay_seconds) - 1)
                        ]
                        sleep_seconds = base_delay_seconds + random.uniform(0.0, 2.0)
                        print(
                            f"Retry {attempt_index}/{max_retry_attempts} for XML "
                            f"batch {batch_index}/{total_batch_count} after "
                            f"{sleep_seconds:.3f} seconds."
                        )
                        time.sleep(sleep_seconds)

                    fetch_handle = None
                    open_latency_seconds: Optional[float] = None
                    read_latency_seconds: Optional[float] = None
                    request_started_at: Optional[float] = None
                    bytes_received_before_failure = 0

                    try:
                        print(
                            f"Starting XML request for batch "
                            f"{batch_index}/{total_batch_count} "
                            f"(attempt {attempt_index}/{max_retry_attempts}) "
                            f"with {len(current_batch_protein_uids)} protein UIDs."
                        )

                        request_started_at = time.perf_counter()
                        fetch_handle = Entrez.efetch(
                            db=database_name,
                            id=",".join(current_batch_protein_uids),
                            retmode=retmode,
                        )
                        response_opened_at = time.perf_counter()
                        open_latency_seconds = response_opened_at - request_started_at
                        response_url = getattr(fetch_handle, "url", None)
                        response_headers = getattr(fetch_handle, "headers", None)

                        xml_payload_bytes, read_latency_seconds, byte_count = (
                            _read_response_body_with_deadline(
                                fetch_handle=fetch_handle,
                                batch_deadline_seconds=batch_deadline_seconds,
                                deadline_started_at=request_started_at,
                            )
                        )
                        total_latency_seconds = time.perf_counter() - request_started_at

                        if not xml_payload_bytes.strip():
                            raise RuntimeError(
                                "NCBI returned an empty XML payload for the current batch."
                            )

                        record_count = _count_xml_batch_records(
                            xml_payload_bytes=xml_payload_bytes,
                        )
                        if record_count != len(current_batch_protein_uids):
                            raise RuntimeError(
                                "NCBI XML batch record count does not match the "
                                "requested identifier count. "
                                f"Expected {len(current_batch_protein_uids)}, "
                                f"got {record_count}."
                            )

                        xml_payload_sha256 = hashlib.sha256(
                            xml_payload_bytes
                        ).hexdigest()
                        average_throughput = (
                            byte_count / read_latency_seconds
                            if read_latency_seconds > 0
                            else float(byte_count)
                        )

                        if batch_xml_file_path is not None:
                            write_bytes_atomic(
                                binary_payload=xml_payload_bytes,
                                output_file_path=batch_xml_file_path,
                            )
                            batch_metadata_payload = {
                                "batch_index": batch_index,
                                "batch_start_index": batch_start_index,
                                "batch_end_index": current_batch_end_index,
                                "protein_uid_count": len(
                                    current_batch_protein_uids
                                ),
                                "protein_uids_sha256": current_batch_uid_sha256,
                                "xml_payload_sha256": xml_payload_sha256,
                                "response_byte_count": byte_count,
                                "record_count": record_count,
                                "open_latency_seconds": open_latency_seconds,
                                "read_latency_seconds": read_latency_seconds,
                                "total_latency_seconds": total_latency_seconds,
                                "average_throughput_bytes_per_second": (
                                    average_throughput
                                ),
                                "attempt_count": attempt_index,
                                "recorded_at_utc": _utc_now_string(),
                            }
                            write_json_atomic(
                                payload=batch_metadata_payload,
                                output_file_path=batch_metadata_file_path,
                            )
                            returned_xml_payload_bytes = b""
                            returned_xml_payload_file_path = str(batch_xml_file_path)
                        else:
                            returned_xml_payload_bytes = xml_payload_bytes
                            returned_xml_payload_file_path = None

                        print(
                            f"Received XML response for batch "
                            f"{batch_index}/{total_batch_count} "
                            f"in {total_latency_seconds:.3f} seconds "
                            f"({byte_count} bytes, "
                            f"{average_throughput:.1f} bytes/s)."
                        )

                        if response_url is not None:
                            print(f"NCBI XML request URL: {response_url}")

                        if response_headers is not None:
                            print(
                                f"NCBI XML response headers: {dict(response_headers)}"
                            )

                        if (
                            b"<ERROR>" in xml_payload_bytes
                            or b"<ErrorList>" in xml_payload_bytes
                            or b"<ERRORS>" in xml_payload_bytes
                            or b"<error>" in xml_payload_bytes
                            or b"<Error>" in xml_payload_bytes
                        ):
                            print(
                                f"Warning: XML batch {batch_index}/"
                                f"{total_batch_count} contains server-side "
                                "error markers despite a successful response."
                            )

                        _record_xml_latency_event(
                            workspace_directory=workspace_directory,
                            event_payload={
                                "batch_index": batch_index,
                                "attempt_index": attempt_index,
                                "status": "success",
                                "open_latency_seconds": open_latency_seconds,
                                "read_latency_seconds": read_latency_seconds,
                                "total_latency_seconds": total_latency_seconds,
                                "bytes_received_before_failure": None,
                                "response_byte_count": byte_count,
                                "average_throughput_bytes_per_second": (
                                    average_throughput
                                ),
                            },
                        )

                        xml_batches.append(
                            NCBIProteinXmlBatchFetchResult(
                                batch_index=batch_index,
                                batch_start_index=batch_start_index,
                                batch_end_index=current_batch_end_index,
                                protein_uids=current_batch_protein_uids,
                                protein_uid_count=len(current_batch_protein_uids),
                                xml_payload_bytes=returned_xml_payload_bytes,
                                xml_payload_sha256=xml_payload_sha256,
                                xml_payload_file_path=returned_xml_payload_file_path,
                                open_latency_seconds=open_latency_seconds,
                                read_latency_seconds=read_latency_seconds,
                                total_latency_seconds=total_latency_seconds,
                                response_byte_count=byte_count,
                                average_throughput_bytes_per_second=average_throughput,
                                attempt_count=attempt_index,
                            )
                        )

                        if total_latency_seconds >= slow_batch_threshold_seconds:
                            consecutive_degraded_batch_count += 1
                        else:
                            consecutive_degraded_batch_count = 0

                        jitter_seconds = random.uniform(0.01, 0.05)
                        time.sleep(resolved_request_delay_seconds + jitter_seconds)
                        batch_succeeded = True
                        break

                    except Exception as error:
                        if _is_ncbi_tls_configuration_error(error):
                            raise RuntimeError(
                                _build_ncbi_tls_configuration_error_message(
                                    error=error,
                                    ssl_ca_file=resolved_ssl_ca_file,
                                    ssl_ca_directory=resolved_ssl_ca_directory,
                                )
                            ) from error

                        bytes_received_before_failure = getattr(
                            error,
                            "bytes_received",
                            bytes_received_before_failure,
                        )
                        total_latency_seconds = (
                            time.perf_counter() - request_started_at
                            if request_started_at is not None
                            else None
                        )
                        error_status = _classify_xml_fetch_error(error)
                        _record_xml_latency_event(
                            workspace_directory=workspace_directory,
                            event_payload={
                                "batch_index": batch_index,
                                "attempt_index": attempt_index,
                                "status": error_status,
                                "open_latency_seconds": open_latency_seconds,
                                "read_latency_seconds": read_latency_seconds,
                                "total_latency_seconds": total_latency_seconds,
                                "bytes_received_before_failure": (
                                    bytes_received_before_failure
                                ),
                                "error_type": error.__class__.__name__,
                                "error_message": str(error),
                            },
                        )

                        if retry_attempt_index == max_retry_attempts - 1:
                            raise RuntimeError(
                                "Failed to fetch XML after "
                                f"{max_retry_attempts} attempts for "
                                f"batch_index={batch_index}, "
                                f"batch_start_index={batch_start_index}: {error}"
                            ) from error

                        consecutive_degraded_batch_count += 1
                    finally:
                        if fetch_handle is not None:
                            fetch_handle.close()

                if not batch_succeeded:
                    raise RuntimeError(
                        f"Failed to fetch XML batch {batch_index}/{total_batch_count}."
                    )

            fetch_result = NCBIProteinXmlFetchResult(
                database_name=database_name,
                identifier_type=identifier_type,
                retrieved_at_utc=retrieved_at_utc,
                requested_protein_uid_count=len(protein_uids),
                normalized_protein_uid_count=len(normalized_protein_uids),
                protein_uids_sha256=protein_uids_sha256,
                batch_size=batch_size,
                batch_count=len(xml_batches),
                retmode=retmode,
                request_delay_seconds=resolved_request_delay_seconds,
                max_retry_attempts=max_retry_attempts,
                python_version=platform.python_version(),
                biopython_version=getattr(Bio, "__version__", "unknown"),
                xml_batches=xml_batches,
                batch_workspace_directory=(
                    str(workspace_directory) if workspace_directory is not None else None
                ),
                socket_idle_timeout_seconds=socket_idle_timeout_seconds,
                batch_deadline_seconds=batch_deadline_seconds,
                circuit_breaker_failure_threshold=(
                    circuit_breaker_failure_threshold
                ),
                circuit_breaker_cooldown_seconds=circuit_breaker_cooldown_seconds,
            )

            print(
                f"Fetched {fetch_result.batch_count} XML batches for "
                f"{fetch_result.normalized_protein_uid_count} protein UIDs."
            )

            return fetch_result
    finally:
        socket.setdefaulttimeout(previous_socket_default_timeout)
        if previous_entrez_max_tries is None:
            try:
                delattr(Entrez, "max_tries")
            except AttributeError:
                pass
        else:
            Entrez.max_tries = previous_entrez_max_tries
