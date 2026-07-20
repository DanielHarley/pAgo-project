from __future__ import annotations

import os
import platform
import random
import socket
import ssl
import time
import json
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlopen as _urllib_urlopen
import hashlib
import Bio
from Bio import Entrez

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
_active_ncbi_network_timeout_seconds: ContextVar[Optional[float]] = ContextVar(
    "_active_ncbi_network_timeout_seconds",
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
    Delegate Entrez HTTPS calls while injecting active request controls.
    """
    active_ssl_context = _active_ncbi_ssl_context.get()
    if active_ssl_context is not None and "context" not in kwargs:
        kwargs["context"] = active_ssl_context

    active_network_timeout_seconds = _active_ncbi_network_timeout_seconds.get()
    if active_network_timeout_seconds is not None and "timeout" not in kwargs:
        kwargs["timeout"] = active_network_timeout_seconds

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

    _ncbi_entrez_delegate_urlopen = (
        current_urlopen if current_urlopen is not None else _urllib_urlopen
    )
    Entrez.urlopen = _ncbi_entrez_urlopen_with_context


@contextmanager
def _ncbi_entrez_request_timeout(
    *,
    network_timeout_seconds: float,
):
    """
    Temporarily apply a socket/network idle timeout to Bio.Entrez urlopen calls.
    """
    _ensure_ncbi_entrez_urlopen_hook_installed()
    network_timeout_token = _active_ncbi_network_timeout_seconds.set(
        network_timeout_seconds,
    )

    try:
        yield
    finally:
        _active_ncbi_network_timeout_seconds.reset(network_timeout_token)


@contextmanager
def _temporary_ncbi_entrez_retry_policy(
    *,
    max_tries: int,
    sleep_between_tries: float,
):
    """
    Temporarily make Bio.Entrez retry behavior explicit for this fetch layer.
    """
    original_max_tries = getattr(Entrez, "max_tries", None)
    original_sleep_between_tries = getattr(Entrez, "sleep_between_tries", None)

    if original_max_tries is not None:
        Entrez.max_tries = max_tries
    if original_sleep_between_tries is not None:
        Entrez.sleep_between_tries = sleep_between_tries

    try:
        yield
    finally:
        if original_max_tries is not None:
            Entrez.max_tries = original_max_tries
        if original_sleep_between_tries is not None:
            Entrez.sleep_between_tries = original_sleep_between_tries


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
    xml_payload_file_path: Optional[str] = None
    response_byte_count: Optional[int] = None
    round_trip_latency_seconds: Optional[float] = None
    attempt_count: Optional[int] = None


class NCBITransientFetchError(RuntimeError):
    """
    Error raised when a fetch attempt failed in a way that may recover on retry.
    """


class NCBIXmlBatchDeadlineExceeded(RuntimeError):
    """
    Error raised when a complete XML batch operation exceeds its total deadline.
    """


class NCBICircuitBreakerOpen(RuntimeError):
    """
    Error raised when recent NCBI failures have opened the circuit breaker.
    """


class NCBIXmlBatchFetchError(RuntimeError):
    """
    Error raised when an XML batch cannot be fetched within the configured policy.
    """


@dataclass
class NCBIXmlCircuitBreaker:
    """
    Track consecutive failed XML batches and temporarily block further requests.
    """
    failure_threshold: int
    cooldown_seconds: float
    consecutive_failure_count: int = 0
    opened_at_seconds: Optional[float] = None

    def __post_init__(self) -> None:
        _validate_positive_integer(
            parameter_name="circuit_breaker_failure_threshold",
            parameter_value=self.failure_threshold,
        )
        _validate_non_negative_float(
            parameter_name="circuit_breaker_cooldown_seconds",
            parameter_value=self.cooldown_seconds,
        )

    def assert_request_allowed(
        self,
        *,
        current_time_seconds: float,
        batch_context: str,
    ) -> None:
        if self.opened_at_seconds is None:
            return

        elapsed_cooldown_seconds = current_time_seconds - self.opened_at_seconds
        if elapsed_cooldown_seconds < self.cooldown_seconds:
            remaining_cooldown_seconds = (
                self.cooldown_seconds - elapsed_cooldown_seconds
            )
            raise NCBICircuitBreakerOpen(
                "NCBI XML circuit breaker is open after "
                f"{self.consecutive_failure_count} consecutive failed batches; "
                f"cooldown has {remaining_cooldown_seconds:.3f} seconds remaining "
                f"before {batch_context} can be requested."
            )

    def record_success(self) -> None:
        self.consecutive_failure_count = 0
        self.opened_at_seconds = None

    def record_failure(self, *, current_time_seconds: float) -> None:
        self.consecutive_failure_count += 1
        if self.consecutive_failure_count >= self.failure_threshold:
            self.opened_at_seconds = current_time_seconds


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
    fetch_timeout_seconds: Optional[float] = None
    batch_deadline_seconds: Optional[float] = None
    retry_backoff_initial_seconds: Optional[float] = None
    retry_backoff_multiplier: Optional[float] = None
    retry_backoff_max_seconds: Optional[float] = None
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


def _xml_batch_file_paths(
    *,
    workspace_directory: Path,
    batch_index: int,
) -> tuple[Path, Path]:
    batch_file_stem = f"batch_{batch_index:06d}"
    batch_directory = workspace_directory / "batches"
    return (
        batch_directory / f"{batch_file_stem}.xml",
        batch_directory / f"{batch_file_stem}.json",
    )


def _count_xml_batch_records(
    *,
    xml_payload_bytes: bytes,
) -> int:
    batch_root = ET.fromstring(xml_payload_bytes)
    return len(list(batch_root))


def _record_xml_latency_event(
    *,
    workspace_directory: Optional[Path],
    event_payload: dict,
) -> None:
    if workspace_directory is None:
        return

    latency_events_file_path = workspace_directory / "latency_events.jsonl"
    latency_events_file_path.parent.mkdir(parents=True, exist_ok=True)
    event_payload = {
        **event_payload,
        "recorded_at_utc": _utc_now_string(),
    }

    with latency_events_file_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(event_payload, sort_keys=True) + "\n")


def _validate_persisted_xml_batch(
    *,
    batch_xml_file_path: Path,
    batch_metadata_file_path: Path,
    expected_batch_index: int,
    expected_batch_start_index: int,
    expected_batch_end_index: int,
    expected_protein_uids: List[str],
) -> Optional[tuple[dict, bytes]]:
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

    try:
        xml_payload_sha256 = sha256_of_file(input_file_path=batch_xml_file_path)
        if batch_metadata.get("xml_payload_sha256") != xml_payload_sha256:
            return None

        xml_payload_bytes = batch_xml_file_path.read_bytes()
        if _xml_payload_contains_server_error_marker(
            xml_payload_bytes=xml_payload_bytes,
        ):
            return None

        record_count = _count_xml_batch_records(
            xml_payload_bytes=xml_payload_bytes,
        )
    except Exception:
        return None

    if record_count != len(expected_protein_uids):
        return None

    return batch_metadata, xml_payload_bytes


def _validate_positive_integer(
    *,
    parameter_name: str,
    parameter_value: int,
) -> None:
    if parameter_value <= 0:
        raise ValueError(f"{parameter_name} must be a positive integer.")


def _validate_positive_float(
    *,
    parameter_name: str,
    parameter_value: float,
) -> None:
    if parameter_value <= 0:
        raise ValueError(f"{parameter_name} must be a positive number.")


def _validate_non_negative_float(
    *,
    parameter_name: str,
    parameter_value: float,
) -> None:
    if parameter_value < 0:
        raise ValueError(f"{parameter_name} must be a non-negative number.")


def _validate_xml_fetch_controls(
    *,
    fetch_timeout_seconds: float,
    batch_deadline_seconds: float,
    request_delay_seconds: float,
    retry_backoff_initial_seconds: float,
    retry_backoff_multiplier: float,
    retry_backoff_max_seconds: float,
    circuit_breaker_failure_threshold: int,
    circuit_breaker_cooldown_seconds: float,
) -> None:
    _validate_positive_float(
        parameter_name="fetch_timeout_seconds",
        parameter_value=fetch_timeout_seconds,
    )
    _validate_positive_float(
        parameter_name="batch_deadline_seconds",
        parameter_value=batch_deadline_seconds,
    )
    _validate_non_negative_float(
        parameter_name="request_delay_seconds",
        parameter_value=request_delay_seconds,
    )
    _validate_positive_float(
        parameter_name="retry_backoff_initial_seconds",
        parameter_value=retry_backoff_initial_seconds,
    )
    if retry_backoff_multiplier < 1:
        raise ValueError(
            "retry_backoff_multiplier must be greater than or equal to 1."
        )
    _validate_positive_float(
        parameter_name="retry_backoff_max_seconds",
        parameter_value=retry_backoff_max_seconds,
    )
    if retry_backoff_max_seconds < retry_backoff_initial_seconds:
        raise ValueError(
            "retry_backoff_max_seconds must be greater than or equal to "
            "retry_backoff_initial_seconds."
        )
    _validate_positive_integer(
        parameter_name="circuit_breaker_failure_threshold",
        parameter_value=circuit_breaker_failure_threshold,
    )
    _validate_non_negative_float(
        parameter_name="circuit_breaker_cooldown_seconds",
        parameter_value=circuit_breaker_cooldown_seconds,
    )


def _build_xml_batch_context(
    *,
    batch_index: int,
    total_batch_count: int,
    batch_start_index: int,
    current_batch_protein_uids: List[str],
) -> str:
    batch_end_index = batch_start_index + len(current_batch_protein_uids) - 1
    return (
        f"batch_index={batch_index}/{total_batch_count}, "
        f"batch_start_index={batch_start_index}, "
        f"batch_end_index={batch_end_index}, "
        f"protein_uids={current_batch_protein_uids}"
    )


def _is_transient_ncbi_fetch_error(error: Exception) -> bool:
    if isinstance(error, NCBITransientFetchError):
        return True

    if isinstance(error, HTTPError):
        return error.code in {408, 429, 500, 502, 503, 504}

    if isinstance(error, (TimeoutError, socket.timeout, URLError)):
        return True

    error_message = str(error).lower()
    return (
        "timed out" in error_message
        or "timeout" in error_message
        or "temporarily unavailable" in error_message
        or "connection reset" in error_message
        or "connection aborted" in error_message
    )


def _xml_payload_contains_server_error_marker(
    *,
    xml_payload_bytes: bytes,
) -> bool:
    return (
        b"<ERROR>" in xml_payload_bytes
        or b"<ErrorList>" in xml_payload_bytes
        or b"<ERRORS>" in xml_payload_bytes
        or b"<error>" in xml_payload_bytes
        or b"<Error>" in xml_payload_bytes
    )


def _calculate_retry_backoff_seconds(
    *,
    retry_attempt_index: int,
    retry_backoff_initial_seconds: float,
    retry_backoff_multiplier: float,
    retry_backoff_max_seconds: float,
) -> float:
    retry_backoff_seconds = retry_backoff_initial_seconds * (
        retry_backoff_multiplier ** retry_attempt_index
    )
    return min(retry_backoff_seconds, retry_backoff_max_seconds)


def _raise_xml_batch_deadline_exceeded(
    *,
    batch_context: str,
    batch_deadline_seconds: float,
    elapsed_seconds: float,
    last_error: Optional[Exception] = None,
) -> None:
    message = (
        "NCBI XML batch deadline exhausted after "
        f"{elapsed_seconds:.3f} seconds with a configured deadline of "
        f"{batch_deadline_seconds:.3f} seconds for {batch_context}."
    )
    if last_error is not None:
        message += f" Last error: {last_error}"

    raise NCBIXmlBatchDeadlineExceeded(message)


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
    batch_workspace_directory: Optional[str | Path] = None,
    fetch_timeout_seconds: float = 30.0,
    batch_deadline_seconds: float = 300.0,
    retry_backoff_initial_seconds: Optional[float] = None,
    retry_backoff_multiplier: float = 2.0,
    retry_backoff_max_seconds: float = 30.0,
    circuit_breaker_failure_threshold: int = 3,
    circuit_breaker_cooldown_seconds: float = 60.0,
    circuit_breaker: Optional[NCBIXmlCircuitBreaker] = None,
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
        batch_workspace_directory:
            Optional directory used to persist and validate individual XML
            batches before reusing them in later runs.
        fetch_timeout_seconds:
            Socket/network idle timeout applied to each NCBI EFetch request.
        batch_deadline_seconds:
            Maximum total time allowed for one complete XML batch operation,
            including retries and backoff delays.
        retry_backoff_initial_seconds:
            Delay before the first retry. If None, uses the resolved request
            delay, with a 0.10 second minimum.
        retry_backoff_multiplier:
            Multiplier applied to each subsequent retry delay.
        retry_backoff_max_seconds:
            Upper bound for any single retry backoff delay.
        circuit_breaker_failure_threshold:
            Number of consecutive failed batches that opens the circuit breaker.
        circuit_breaker_cooldown_seconds:
            Cooldown interval before another request is allowed after opening.
        circuit_breaker:
            Optional shared circuit breaker state. Supplying the same instance
            across calls lets sustained failures block later fetch operations.

    Returns:
        NCBIProteinXmlFetchResult:
            Retrieval metadata plus raw XML payload bytes for each batch.
    """
    if not ncbi_email:
        raise ValueError(
            "No NCBI email credential found. Please add .env file with "
            "NCBI_EMAIL. NCBI_API_KEY is optional."
        )

    _validate_positive_integer(
        parameter_name="batch_size",
        parameter_value=batch_size,
    )

    _validate_positive_integer(
        parameter_name="max_retry_attempts",
        parameter_value=max_retry_attempts,
    )

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

    resolved_retry_backoff_initial_seconds = (
        max(resolved_request_delay_seconds, 0.10)
        if retry_backoff_initial_seconds is None
        else retry_backoff_initial_seconds
    )

    _validate_xml_fetch_controls(
        fetch_timeout_seconds=fetch_timeout_seconds,
        batch_deadline_seconds=batch_deadline_seconds,
        request_delay_seconds=resolved_request_delay_seconds,
        retry_backoff_initial_seconds=resolved_retry_backoff_initial_seconds,
        retry_backoff_multiplier=retry_backoff_multiplier,
        retry_backoff_max_seconds=retry_backoff_max_seconds,
        circuit_breaker_failure_threshold=circuit_breaker_failure_threshold,
        circuit_breaker_cooldown_seconds=circuit_breaker_cooldown_seconds,
    )

    active_circuit_breaker = circuit_breaker
    if active_circuit_breaker is None:
        active_circuit_breaker = NCBIXmlCircuitBreaker(
            failure_threshold=circuit_breaker_failure_threshold,
            cooldown_seconds=circuit_breaker_cooldown_seconds,
        )

    retrieved_at_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    protein_uids_sha256 = sha256_of_lines(
        text_lines=normalized_protein_uids,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )
    workspace_directory = (
        Path(batch_workspace_directory).expanduser()
        if batch_workspace_directory is not None
        else None
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
            current_batch_uid_sha256 = sha256_of_lines(
                text_lines=current_batch_protein_uids,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            )
            batch_xml_file_path: Optional[Path] = None
            batch_metadata_file_path: Optional[Path] = None
            if workspace_directory is not None:
                batch_xml_file_path, batch_metadata_file_path = _xml_batch_file_paths(
                    workspace_directory=workspace_directory,
                    batch_index=batch_index,
                )
                persisted_batch = _validate_persisted_xml_batch(
                    batch_xml_file_path=batch_xml_file_path,
                    batch_metadata_file_path=batch_metadata_file_path,
                    expected_batch_index=batch_index,
                    expected_batch_start_index=batch_start_index,
                    expected_batch_end_index=current_batch_end_index,
                    expected_protein_uids=current_batch_protein_uids,
                )
                if persisted_batch is not None:
                    batch_metadata, xml_payload_bytes = persisted_batch
                    print(
                        f"Skipping validated XML batch "
                        f"{batch_index}/{total_batch_count}."
                    )
                    _record_xml_latency_event(
                        workspace_directory=workspace_directory,
                        event_payload={
                            "batch_index": batch_index,
                            "attempt_index": 0,
                            "status": "reused_validated_batch",
                            "response_byte_count": batch_metadata.get(
                                "response_byte_count"
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
                            xml_payload_bytes=xml_payload_bytes,
                            xml_payload_sha256=batch_metadata[
                                "xml_payload_sha256"
                            ],
                            xml_payload_file_path=str(batch_xml_file_path),
                            response_byte_count=batch_metadata.get(
                                "response_byte_count"
                            ),
                            round_trip_latency_seconds=batch_metadata.get(
                                "round_trip_latency_seconds"
                            ),
                            attempt_count=batch_metadata.get("attempt_count"),
                        )
                    )
                    active_circuit_breaker.record_success()
                    continue

            batch_context = _build_xml_batch_context(
                batch_index=batch_index,
                total_batch_count=total_batch_count,
                batch_start_index=batch_start_index,
                current_batch_protein_uids=current_batch_protein_uids,
            )
            batch_started_at_seconds = time.monotonic()
            batch_deadline_at_seconds = (
                batch_started_at_seconds + batch_deadline_seconds
            )
            active_circuit_breaker.assert_request_allowed(
                current_time_seconds=batch_started_at_seconds,
                batch_context=batch_context,
            )
            last_transient_error: Optional[Exception] = None

            for retry_attempt_index in range(max_retry_attempts):
                try:
                    attempt_started_at_seconds = time.monotonic()
                    if attempt_started_at_seconds >= batch_deadline_at_seconds:
                        _raise_xml_batch_deadline_exceeded(
                            batch_context=batch_context,
                            batch_deadline_seconds=batch_deadline_seconds,
                            elapsed_seconds=(
                                attempt_started_at_seconds - batch_started_at_seconds
                            ),
                            last_error=last_transient_error,
                        )

                    remaining_batch_seconds = (
                        batch_deadline_at_seconds - attempt_started_at_seconds
                    )
                    effective_fetch_timeout_seconds = min(
                        fetch_timeout_seconds,
                        remaining_batch_seconds,
                    )

                    print(
                        f"Starting XML request for batch {batch_index}/{total_batch_count} "
                        f"(attempt {retry_attempt_index + 1}/{max_retry_attempts}) "
                        f"with {len(current_batch_protein_uids)} protein UIDs."
                    )

                    request_started_at = time.perf_counter()
                    with _ncbi_entrez_request_timeout(
                        network_timeout_seconds=effective_fetch_timeout_seconds,
                    ):
                        with _temporary_ncbi_entrez_retry_policy(
                            max_tries=1,
                            sleep_between_tries=0.0,
                        ):
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
                    attempt_finished_at_seconds = time.monotonic()
                    if attempt_finished_at_seconds > batch_deadline_at_seconds:
                        _raise_xml_batch_deadline_exceeded(
                            batch_context=batch_context,
                            batch_deadline_seconds=batch_deadline_seconds,
                            elapsed_seconds=(
                                attempt_finished_at_seconds
                                - batch_started_at_seconds
                            ),
                            last_error=last_transient_error,
                        )

                    xml_payload_bytes = _coerce_payload_to_bytes(
                        payload=raw_xml_payload,
                    )

                    if not xml_payload_bytes.strip():
                        raise NCBITransientFetchError(
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

                    if _xml_payload_contains_server_error_marker(
                        xml_payload_bytes=xml_payload_bytes,
                    ):
                        raise NCBITransientFetchError(
                            "NCBI XML payload contains server-side error markers "
                            "despite a successful HTTP response."
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

                    response_byte_count = len(xml_payload_bytes)
                    attempt_count = retry_attempt_index + 1
                    returned_xml_payload_file_path: Optional[str] = None

                    if (
                        batch_xml_file_path is not None
                        and batch_metadata_file_path is not None
                    ):
                        write_bytes_atomic(
                            binary_payload=xml_payload_bytes,
                            output_file_path=batch_xml_file_path,
                        )
                        batch_metadata_payload = {
                            "batch_index": batch_index,
                            "batch_start_index": batch_start_index,
                            "batch_end_index": current_batch_end_index,
                            "protein_uid_count": len(current_batch_protein_uids),
                            "protein_uids_sha256": current_batch_uid_sha256,
                            "xml_payload_sha256": xml_payload_sha256,
                            "record_count": record_count,
                            "response_byte_count": response_byte_count,
                            "round_trip_latency_seconds": round_trip_latency_seconds,
                            "attempt_count": attempt_count,
                            "fetch_timeout_seconds": effective_fetch_timeout_seconds,
                            "batch_deadline_seconds": batch_deadline_seconds,
                            "recorded_at_utc": _utc_now_string(),
                        }
                        write_json_atomic(
                            payload=batch_metadata_payload,
                            output_file_path=batch_metadata_file_path,
                        )
                        returned_xml_payload_file_path = str(batch_xml_file_path)

                    _record_xml_latency_event(
                        workspace_directory=workspace_directory,
                        event_payload={
                            "batch_index": batch_index,
                            "attempt_index": attempt_count,
                            "status": "success",
                            "round_trip_latency_seconds": round_trip_latency_seconds,
                            "response_byte_count": response_byte_count,
                        },
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
                            xml_payload_file_path=returned_xml_payload_file_path,
                            response_byte_count=response_byte_count,
                            round_trip_latency_seconds=round_trip_latency_seconds,
                            attempt_count=attempt_count,
                        )
                    )
                    active_circuit_breaker.record_success()

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

                    if isinstance(error, NCBIXmlBatchDeadlineExceeded):
                        failure_status = "deadline_exhausted"
                    elif _is_transient_ncbi_fetch_error(error):
                        failure_status = "transient_failure"
                    else:
                        failure_status = "permanent_failure"

                    _record_xml_latency_event(
                        workspace_directory=workspace_directory,
                        event_payload={
                            "batch_index": batch_index,
                            "attempt_index": retry_attempt_index + 1,
                            "status": failure_status,
                            "error_type": type(error).__name__,
                            "error_message": str(error),
                        },
                    )

                    if isinstance(error, NCBIXmlBatchDeadlineExceeded):
                        active_circuit_breaker.record_failure(
                            current_time_seconds=time.monotonic(),
                        )
                        raise error

                    if not _is_transient_ncbi_fetch_error(error):
                        raise RuntimeError(
                            "Permanent NCBI XML batch failure; not retrying "
                            f"{batch_context}: {error}"
                        ) from error

                    last_transient_error = error

                    if retry_attempt_index == max_retry_attempts - 1:
                        active_circuit_breaker.record_failure(
                            current_time_seconds=time.monotonic(),
                        )
                        raise NCBIXmlBatchFetchError(
                            f"Failed to fetch NCBI XML batch after "
                            f"{max_retry_attempts} attempts for {batch_context}: "
                            f"{error}"
                        ) from error

                    retry_backoff_seconds = _calculate_retry_backoff_seconds(
                        retry_attempt_index=retry_attempt_index,
                        retry_backoff_initial_seconds=(
                            resolved_retry_backoff_initial_seconds
                        ),
                        retry_backoff_multiplier=retry_backoff_multiplier,
                        retry_backoff_max_seconds=retry_backoff_max_seconds,
                    )
                    retry_scheduled_at_seconds = time.monotonic()
                    remaining_batch_seconds = (
                        batch_deadline_at_seconds - retry_scheduled_at_seconds
                    )
                    if retry_backoff_seconds >= remaining_batch_seconds:
                        active_circuit_breaker.record_failure(
                            current_time_seconds=retry_scheduled_at_seconds,
                        )
                        _raise_xml_batch_deadline_exceeded(
                            batch_context=batch_context,
                            batch_deadline_seconds=batch_deadline_seconds,
                            elapsed_seconds=(
                                retry_scheduled_at_seconds
                                - batch_started_at_seconds
                            ),
                            last_error=error,
                        )

                    print(
                        f"Retry {retry_attempt_index + 1}/{max_retry_attempts} "
                        f"after transient error: {error}"
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
            batch_workspace_directory=(
                str(workspace_directory) if workspace_directory is not None else None
            ),
            fetch_timeout_seconds=fetch_timeout_seconds,
            batch_deadline_seconds=batch_deadline_seconds,
            retry_backoff_initial_seconds=resolved_retry_backoff_initial_seconds,
            retry_backoff_multiplier=retry_backoff_multiplier,
            retry_backoff_max_seconds=retry_backoff_max_seconds,
            circuit_breaker_failure_threshold=circuit_breaker_failure_threshold,
            circuit_breaker_cooldown_seconds=circuit_breaker_cooldown_seconds,
        )

        print(
            f"Fetched {fetch_result.batch_count} XML batches for "
            f"{fetch_result.normalized_protein_uid_count} protein UIDs."
        )

        return fetch_result
