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
import Bio
from Bio import Entrez

from src.pago_pipeline.storage import sha256_of_lines


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


def fetch_ncbi_protein_xml_batches(
    *,
    ncbi_email: str,
    ncbi_api_key: Optional[str],
    protein_uids: List[str],
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    batch_size: int = 100,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
) -> NCBIProteinXmlFetchResult:
    """
    Fetch raw protein XML batches from NCBI Entrez using a frozen list of protein UIDs.

    This function is intended to run downstream of a previously frozen UID snapshot.
    It performs remote retrieval only and returns in-memory batch payloads as raw bytes.

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

    normalized_protein_uids = _normalize_protein_uid_list(
        protein_uids=protein_uids,
        deduplicate_uids=False,
        sort_uids=False,
    )

    if not normalized_protein_uids:
        raise ValueError("protein_uids must contain at least one non-empty UID.")

    Entrez.email = ncbi_email
    Entrez.api_key = ncbi_api_key

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

    with _configured_ncbi_entrez_urlopen(
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
    ) as (resolved_ssl_ca_file, resolved_ssl_ca_directory):
        xml_batches: List[NCBIProteinXmlBatchFetchResult] = []
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

            for retry_attempt_index in range(max_retry_attempts):
                try:
                    print(
                        f"Starting XML request for batch {batch_index}/{total_batch_count} "
                        f"(attempt {retry_attempt_index + 1}/{max_retry_attempts}) "
                        f"with {len(current_batch_protein_uids)} protein UIDs."
                    )

                    request_started_at = time.perf_counter()
                    fetch_handle = Entrez.efetch(
                        db=database_name,
                        id=",".join(current_batch_protein_uids),
                        retmode=retmode,
                    )
                    response_url = getattr(fetch_handle, "url", None)
                    response_headers = getattr(fetch_handle, "headers", None)
                    try:
                        raw_xml_payload = fetch_handle.read()
                    finally:
                        fetch_handle.close()

                    round_trip_latency_seconds = (
                        time.perf_counter() - request_started_at
                    )

                    xml_payload_bytes = _coerce_payload_to_bytes(
                        payload=raw_xml_payload,
                    )

                    if not xml_payload_bytes.strip():
                        raise RuntimeError(
                            "NCBI returned an empty XML payload for the current batch."
                        )

                    xml_payload_sha256 = hashlib.sha256(
                        xml_payload_bytes
                    ).hexdigest()

                    print(
                        f"Received XML response for batch {batch_index}/{total_batch_count} "
                        f"in {round_trip_latency_seconds:.3f} seconds "
                        f"({len(xml_payload_bytes)} bytes)."
                    )

                    if response_url is not None:
                        print(f"NCBI XML request URL: {response_url}")

                    if response_headers is not None:
                        print(f"NCBI XML response headers: {dict(response_headers)}")

                    if (
                        b"<ERROR>" in xml_payload_bytes
                        or b"<ErrorList>" in xml_payload_bytes
                        or b"<ERRORS>" in xml_payload_bytes
                        or b"<error>" in xml_payload_bytes
                        or b"<Error>" in xml_payload_bytes
                    ):
                        print(
                            f"Warning: XML batch {batch_index}/{total_batch_count} "
                            f"contains server-side error markers despite a successful response."
                        )

                    xml_batches.append(
                        NCBIProteinXmlBatchFetchResult(
                            batch_index=batch_index,
                            batch_start_index=batch_start_index,
                            batch_end_index=current_batch_end_index,
                            protein_uids=current_batch_protein_uids,
                            protein_uid_count=len(current_batch_protein_uids),
                            xml_payload_bytes=xml_payload_bytes,
                            xml_payload_sha256=xml_payload_sha256,
                        )
                    )

                    print(
                        f"Fetched XML batch {batch_index}/{total_batch_count} "
                        f"with {len(current_batch_protein_uids)} protein UIDs."
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
                            f"Failed to fetch XML after {max_retry_attempts} attempts "
                            f"for batch_index={batch_index}, "
                            f"batch_start_index={batch_start_index}: {error}"
                        ) from error

                    retry_backoff_seconds = (
                        2 ** retry_attempt_index
                    ) * resolved_request_delay_seconds

                    print(
                        f"Retry {retry_attempt_index + 1}/{max_retry_attempts} "
                        f"after error: {error}"
                    )
                    time.sleep(retry_backoff_seconds)

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
        )

        print(
            f"Fetched {fetch_result.batch_count} XML batches for "
            f"{fetch_result.normalized_protein_uid_count} protein UIDs."
        )

        return fetch_result
