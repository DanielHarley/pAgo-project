from __future__ import annotations

import hashlib
import os
import platform
import random
import socket
import ssl
import threading
import time
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from contextvars import ContextVar, copy_context
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, List, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen as _urllib_urlopen

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
_active_ncbi_network_timeout_seconds: ContextVar[Optional[float]] = ContextVar(
    "_active_ncbi_network_timeout_seconds",
    default=None,
)
_ncbi_entrez_delegate_urlopen = getattr(Entrez, "urlopen", None)
_default_ncbi_xml_circuit_breakers_lock = threading.RLock()
_default_ncbi_xml_circuit_breakers: dict[
    tuple[int, float],
    "NCBIXmlCircuitBreaker",
] = {}


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


def _open_ncbi_entrez_efetch_once(
    *,
    database_name: str,
    protein_uids: List[str],
    retmode: str,
):
    """
    Open one EFetch request without Bio.Entrez's module-global retry policy.

    Bio.Entrez.efetch delegates to a private retry loop configured through the
    module-global ``max_tries`` and ``sleep_between_tries`` values. Building the
    equivalent EFetch request here and opening it once keeps retry ownership in
    this layer without mutating process-global Bio.Entrez state.
    """
    request_parameters = {
        "db": database_name,
        "id": ",".join(protein_uids),
        "retmode": retmode,
        "tool": "biopython",
        "email": getattr(Entrez, "email", None),
        "api_key": getattr(Entrez, "api_key", None),
    }
    encoded_parameters = urlencode(
        {
            parameter_name: parameter_value
            for parameter_name, parameter_value in request_parameters.items()
            if parameter_value is not None
        }
    )
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    use_post_request = (
        len(encoded_parameters) > 1000 or len(protein_uids) >= 200
    )
    if use_post_request:
        request = Request(
            efetch_url,
            data=encoded_parameters.encode("utf-8"),
            method="POST",
        )
    else:
        request = Request(
            f"{efetch_url}?{encoded_parameters}",
            method="GET",
        )

    return Entrez.urlopen(request)


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


class NCBIXmlResponseValidationError(RuntimeError):
    """
    Error raised when an XML response is complete enough to read but invalid.
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
    _state_lock: threading.RLock = field(
        default_factory=threading.RLock,
        init=False,
        repr=False,
        compare=False,
    )

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
        with self._state_lock:
            if self.opened_at_seconds is None:
                return

            elapsed_cooldown_seconds = (
                current_time_seconds - self.opened_at_seconds
            )
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
        with self._state_lock:
            self.consecutive_failure_count = 0
            self.opened_at_seconds = None

    def record_failure(self, *, current_time_seconds: float) -> None:
        with self._state_lock:
            self.consecutive_failure_count += 1
            if self.consecutive_failure_count >= self.failure_threshold:
                self.opened_at_seconds = current_time_seconds


def get_default_ncbi_xml_circuit_breaker(
    *,
    failure_threshold: int,
    cooldown_seconds: float,
) -> NCBIXmlCircuitBreaker:
    """
    Return process-persistent circuit-breaker state for one fetch policy.

    The default state is shared across fetch calls, so repeated workflow
    invocations accumulate failures without requiring application callers to
    construct and retain a circuit-breaker instance themselves.
    """
    breaker_key = (failure_threshold, float(cooldown_seconds))
    with _default_ncbi_xml_circuit_breakers_lock:
        circuit_breaker = _default_ncbi_xml_circuit_breakers.get(breaker_key)
        if circuit_breaker is None:
            circuit_breaker = NCBIXmlCircuitBreaker(
                failure_threshold=failure_threshold,
                cooldown_seconds=cooldown_seconds,
            )
            _default_ncbi_xml_circuit_breakers[breaker_key] = circuit_breaker

    return circuit_breaker


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


def _count_xml_batch_records(
    *,
    xml_payload_bytes: bytes,
) -> int:
    batch_root = ET.fromstring(xml_payload_bytes)
    return len(list(batch_root))


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


def _fetch_ncbi_xml_attempt_with_deadline(
    *,
    database_name: str,
    protein_uids: List[str],
    retmode: str,
    effective_fetch_timeout_seconds: float,
    batch_deadline_at_seconds: float,
    batch_started_at_seconds: float,
    batch_deadline_seconds: float,
    batch_context: str,
) -> tuple[str | bytes, Any, Any]:
    """
    Open, read, and close one EFetch response within an absolute deadline.

    ``urllib`` timeouts bound individual blocking network operations, but they
    do not provide a scheduler-level deadline for the complete response. The
    request therefore runs in a daemon worker while the calling thread waits
    only until the remaining batch deadline. On deadline exhaustion, the
    handle is closed asynchronously when available so a blocked read is asked
    to stop without delaying the caller beyond the configured total deadline.
    """
    attempt_completed = threading.Event()
    cancellation_requested = threading.Event()
    handle_closed = threading.Event()
    attempt_state_lock = threading.Lock()
    handle_close_lock = threading.Lock()
    attempt_state: dict[str, Any] = {}

    def close_fetch_handle_once() -> None:
        with handle_close_lock:
            if handle_closed.is_set():
                return

            with attempt_state_lock:
                fetch_handle = attempt_state.get("fetch_handle")

            if fetch_handle is None:
                return

            try:
                fetch_handle.close()
            finally:
                handle_closed.set()

    def execute_attempt() -> None:
        try:
            with _ncbi_entrez_request_timeout(
                network_timeout_seconds=effective_fetch_timeout_seconds,
            ):
                fetch_handle = _open_ncbi_entrez_efetch_once(
                    database_name=database_name,
                    protein_uids=protein_uids,
                    retmode=retmode,
                )
                with attempt_state_lock:
                    attempt_state["fetch_handle"] = fetch_handle
                    attempt_state["response_url"] = getattr(
                        fetch_handle,
                        "url",
                        None,
                    )
                    attempt_state["response_headers"] = getattr(
                        fetch_handle,
                        "headers",
                        None,
                    )

                if cancellation_requested.is_set():
                    return

                raw_xml_payload = fetch_handle.read()
                with attempt_state_lock:
                    attempt_state["raw_xml_payload"] = raw_xml_payload
        except BaseException as error:
            with attempt_state_lock:
                attempt_state["error"] = error
        finally:
            try:
                close_fetch_handle_once()
            except BaseException as close_error:
                with attempt_state_lock:
                    if "error" not in attempt_state:
                        attempt_state["error"] = close_error
            finally:
                attempt_completed.set()

    def cancel_attempt() -> None:
        try:
            close_fetch_handle_once()
        except BaseException:
            # Deadline exhaustion is the primary failure. The worker records
            # close errors when it owns normal cleanup; asynchronous
            # cancellation must not emit an unhandled background exception.
            pass

    attempt_context = copy_context()
    attempt_thread = threading.Thread(
        target=attempt_context.run,
        args=(execute_attempt,),
        name="ncbi-xml-fetch-attempt",
        daemon=True,
    )
    attempt_thread.start()

    remaining_batch_seconds = (
        batch_deadline_at_seconds - time.monotonic()
    )
    attempt_finished_before_deadline = (
        remaining_batch_seconds > 0
        and attempt_completed.wait(timeout=remaining_batch_seconds)
    )
    if not attempt_finished_before_deadline:
        cancellation_requested.set()
        threading.Thread(
            target=cancel_attempt,
            name="ncbi-xml-fetch-canceller",
            daemon=True,
        ).start()
        deadline_detected_at_seconds = time.monotonic()
        _raise_xml_batch_deadline_exceeded(
            batch_context=batch_context,
            batch_deadline_seconds=batch_deadline_seconds,
            elapsed_seconds=(
                deadline_detected_at_seconds - batch_started_at_seconds
            ),
        )

    with attempt_state_lock:
        attempt_error = attempt_state.get("error")
        raw_xml_payload = attempt_state.get("raw_xml_payload")
        response_url = attempt_state.get("response_url")
        response_headers = attempt_state.get("response_headers")

    if attempt_error is not None:
        raise attempt_error

    return raw_xml_payload, response_url, response_headers


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
        fetch_timeout_seconds:
            Socket/network idle timeout applied to opening and reading each
            NCBI EFetch response. A blocked response read is interrupted by
            the underlying socket timeout.
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
            Optional explicitly scoped circuit breaker. When omitted, a
            process-persistent breaker for the configured threshold and
            cooldown is reused automatically across calls.

    Partial responses:
        NCBI may omit an invalid, removed, or otherwise unavailable UID while
        returning syntactically valid XML. A batch with fewer or more direct
        XML records than requested UIDs is treated as a permanent response
        validation failure and is not retried because retrying cannot establish
        which UID was omitted.

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
        active_circuit_breaker = get_default_ncbi_xml_circuit_breaker(
            failure_threshold=circuit_breaker_failure_threshold,
            cooldown_seconds=circuit_breaker_cooldown_seconds,
        )

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
                    (
                        raw_xml_payload,
                        response_url,
                        response_headers,
                    ) = _fetch_ncbi_xml_attempt_with_deadline(
                        database_name=database_name,
                        protein_uids=current_batch_protein_uids,
                        retmode=retmode,
                        effective_fetch_timeout_seconds=(
                            effective_fetch_timeout_seconds
                        ),
                        batch_deadline_at_seconds=batch_deadline_at_seconds,
                        batch_started_at_seconds=batch_started_at_seconds,
                        batch_deadline_seconds=batch_deadline_seconds,
                        batch_context=batch_context,
                    )

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
                        raise NCBIXmlResponseValidationError(
                            "NCBI returned a partial or mismatched XML batch; "
                            "the direct XML record count does not match the "
                            "requested protein UID count. "
                            f"Expected {len(current_batch_protein_uids)}, "
                            f"got {record_count}; {batch_context}."
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
