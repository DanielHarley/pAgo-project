from __future__ import annotations

import hashlib
import http.client
import io
import json
import os
import platform
import random
import socket
import ssl
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from contextvars import ContextVar, copy_context
from dataclasses import dataclass, field, replace
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Union
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen as _urllib_urlopen

import Bio
from Bio import Entrez

from src.pago_pipeline.ncbi_batch_workspace import (
    NCBIXmlBatchWorkspace,
    PlannedXmlBatch,
    build_xml_batch_plan,
)
from src.pago_pipeline.ncbi_telemetry import (
    DEADLINE_FAILURE_KIND,
    HTTP_429_FAILURE_KIND,
    HTTP_5XX_FAILURE_KIND,
    TRUNCATED_RESPONSE_FAILURE_KIND,
    NCBIAdaptiveConcurrencyGovernor,
    NCBIRequestStartRateLimiter,
    NCBIRetrievalTelemetry,
    OTHER_FAILURE_KIND,
    RATE_LIMIT_FAILURE_KIND,
    TIMEOUT_FAILURE_KIND,
    UID_PAGE_STAGE_NAME,
    UID_SEARCH_STAGE_NAME,
    VALIDATION_FAILURE_KIND,
    XML_BATCH_STAGE_NAME,
)
from src.pago_pipeline.ncbi_xml_stream import validate_xml_batch_payload
from src.pago_pipeline.storage import sha256_of_lines

PathLike = Union[str, Path]

NCBI_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EFETCH_URL = f"{NCBI_EUTILS_BASE_URL}/efetch.fcgi"

# NBK25499 documents ``rettype=gp`` with ``retmode=xml`` as the GBSeq XML
# combination for db=protein. Omitting rettype currently resolves to the same
# payload, but that resolution is undocumented and may change.
NCBI_PROTEIN_GBSEQ_XML_RETTYPE = "gp"
NCBI_PROTEIN_UID_LIST_RETTYPE = "uilist"

# EFetch caps retstart/retmax paging at 10,000 identifiers per request.
NCBI_MAX_EFETCH_RETMAX = 10_000

# NCBI's published allowance is 10 requests/second with an API key and 3
# without. The defaults below stay inside those ceilings with margin because
# the limit aggregates every request issued under the same key.
NCBI_MAX_REQUEST_STARTS_PER_SECOND_WITH_API_KEY = 8.0
NCBI_MAX_REQUEST_STARTS_PER_SECOND_WITHOUT_API_KEY = 2.5

_NCBI_RATE_LIMIT_BODY_MARKERS = (
    b"API rate limit exceeded",
    b"api rate limit exceeded",
    b"Too Many Requests",
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


class _NCBIPersistentHttpsResponse:
    """
    Expose one keep-alive HTTP response through the urlopen handle contract.

    The retrieval layer only requires ``read``, ``close``, ``url``, and
    ``headers``. Tracking whether the body was fully consumed matters because a
    keep-alive connection whose response was abandoned mid-stream, as happens on
    deadline cancellation, can no longer be reused for the next request.
    """

    def __init__(
        self,
        *,
        http_response,
        transport: "_NCBIPersistentHttpsTransport",
        request_url: str,
        connect_seconds: Optional[float],
    ) -> None:
        self._http_response = http_response
        self._transport = transport
        self._response_fully_read = False
        self.url = request_url
        self.headers = http_response.headers
        self.status = http_response.status
        self.ncbi_connect_seconds = connect_seconds

    def read(self, *args, **kwargs):
        payload = self._http_response.read(*args, **kwargs)
        if not args and not kwargs:
            self._response_fully_read = True

        return payload

    def close(self) -> None:
        try:
            if not self._response_fully_read:
                self._transport.discard_current_connection()
        finally:
            self._http_response.close()


class _NCBIPersistentHttpsTransport:
    """
    Reuse one HTTPS connection per thread across consecutive NCBI requests.

    ``urllib`` opens and discards a connection for every call, so each request
    pays a full TCP and TLS handshake. Holding one ``http.client`` connection
    per thread removes that cost while keeping the same SSL context, and
    therefore the same custom certificate-authority behavior, that the urllib
    path uses.
    """

    def __init__(self) -> None:
        self._thread_state = threading.local()
        self._open_connections_lock = threading.Lock()
        self._open_connections: list[Any] = []

    def _current_connection(self):
        return getattr(self._thread_state, "connection", None)

    @staticmethod
    def _close_connection_quietly(connection) -> None:
        try:
            connection.close()
        except Exception:
            pass

    def _forget_connection(self, connection) -> None:
        with self._open_connections_lock:
            if connection in self._open_connections:
                self._open_connections.remove(connection)

    def discard_current_connection(self) -> None:
        connection = self._current_connection()
        self._thread_state.connection = None
        if connection is None:
            return

        self._forget_connection(connection)
        self._close_connection_quietly(connection)

    def discard_all_connections(self) -> None:
        """
        Close every connection this transport opened, on any thread.

        Worker threads keep their connection in thread-local state, so the
        thread that owns a run cannot reach them. Without a central registry
        those sockets would stay open until the connection objects were
        collected, which is not a guarantee worth relying on.
        """
        with self._open_connections_lock:
            connections_to_close = list(self._open_connections)
            self._open_connections.clear()

        self._thread_state.connection = None
        for connection in connections_to_close:
            self._close_connection_quietly(connection)

    def _acquire_connection(
        self,
        *,
        host_name: str,
        port_number: Optional[int],
        timeout_seconds: Optional[float],
        ssl_context: Optional[ssl.SSLContext],
    ) -> tuple[Any, Optional[float], bool]:
        connection = self._current_connection()
        connection_signature = (host_name, port_number, timeout_seconds)

        if (
            connection is not None
            and getattr(self._thread_state, "signature", None)
            == connection_signature
            and getattr(self._thread_state, "ssl_context", None) is ssl_context
        ):
            return connection, None, True

        self.discard_current_connection()

        connection = http.client.HTTPSConnection(
            host_name,
            port=port_number,
            timeout=timeout_seconds,
            context=ssl_context,
        )
        connect_started_at_seconds = time.perf_counter()
        connection.connect()
        connect_seconds = time.perf_counter() - connect_started_at_seconds

        self._thread_state.connection = connection
        self._thread_state.signature = connection_signature
        self._thread_state.ssl_context = ssl_context
        with self._open_connections_lock:
            self._open_connections.append(connection)

        return connection, connect_seconds, False

    def open_request(
        self,
        *,
        request_url: str,
        http_method: str,
        request_body: Optional[bytes],
        timeout_seconds: Optional[float],
        ssl_context: Optional[ssl.SSLContext],
    ) -> _NCBIPersistentHttpsResponse:
        split_url = urlsplit(request_url)
        request_target = split_url.path or "/"
        if split_url.query:
            request_target = f"{request_target}?{split_url.query}"

        request_headers = {
            "Connection": "keep-alive",
            "Accept-Encoding": "identity",
            "User-Agent": "pago-pipeline/ncbi-retrieval",
        }
        if request_body is not None:
            request_headers["Content-Type"] = (
                "application/x-www-form-urlencoded"
            )

        last_error: Optional[Exception] = None
        for connection_attempt_index in range(2):
            connection, connect_seconds, connection_was_reused = (
                self._acquire_connection(
                    host_name=split_url.hostname,
                    port_number=split_url.port,
                    timeout_seconds=timeout_seconds,
                    ssl_context=ssl_context,
                )
            )

            try:
                connection.request(
                    http_method,
                    request_target,
                    body=request_body,
                    headers=request_headers,
                )
                http_response = connection.getresponse()
            except (http.client.HTTPException, OSError) as error:
                self.discard_current_connection()
                last_error = error
                if connection_was_reused and connection_attempt_index == 0:
                    # A reused keep-alive connection can be closed by the peer
                    # between requests. One reconnect distinguishes that from a
                    # genuine transport failure.
                    continue
                raise

            if http_response.status >= 400:
                response_headers = http_response.headers
                try:
                    http_response.read()
                finally:
                    http_response.close()

                raise HTTPError(
                    request_url,
                    http_response.status,
                    http_response.reason,
                    response_headers,
                    None,
                )

            return _NCBIPersistentHttpsResponse(
                http_response=http_response,
                transport=self,
                request_url=request_url,
                connect_seconds=connect_seconds,
            )

        raise last_error if last_error is not None else RuntimeError(
            "Failed to open a persistent NCBI HTTPS request."
        )


_ncbi_persistent_https_transport = _NCBIPersistentHttpsTransport()


def _build_ncbi_efetch_request_parameters(
    *,
    database_name: str,
    retmode: str,
    rettype: Optional[str] = None,
    protein_uids: Optional[Sequence[str]] = None,
    history_web_env: Optional[str] = None,
    history_query_key: Optional[str] = None,
    retstart: Optional[int] = None,
    retmax: Optional[int] = None,
) -> Dict[str, str]:
    """
    Build the EFetch query parameters for either an explicit UID list or History.
    """
    request_parameters: Dict[str, Any] = {
        "db": database_name,
        "retmode": retmode,
        "rettype": rettype,
        "tool": "biopython",
        "email": getattr(Entrez, "email", None),
        "api_key": getattr(Entrez, "api_key", None),
    }

    if protein_uids is not None:
        request_parameters["id"] = ",".join(protein_uids)

    if history_web_env is not None:
        request_parameters["WebEnv"] = history_web_env

    if history_query_key is not None:
        request_parameters["query_key"] = history_query_key

    if retstart is not None:
        request_parameters["retstart"] = str(retstart)

    if retmax is not None:
        request_parameters["retmax"] = str(retmax)

    return {
        parameter_name: parameter_value
        for parameter_name, parameter_value in request_parameters.items()
        if parameter_value is not None
    }


def _open_ncbi_entrez_efetch_once(
    *,
    database_name: str,
    retmode: str,
    rettype: Optional[str] = None,
    protein_uids: Optional[Sequence[str]] = None,
    history_web_env: Optional[str] = None,
    history_query_key: Optional[str] = None,
    retstart: Optional[int] = None,
    retmax: Optional[int] = None,
    reuse_http_connection: bool = False,
):
    """
    Open one EFetch request without Bio.Entrez's module-global retry policy.

    Bio.Entrez.efetch delegates to a private retry loop configured through the
    module-global ``max_tries`` and ``sleep_between_tries`` values. Building the
    equivalent EFetch request here and opening it once keeps retry ownership in
    this layer without mutating process-global Bio.Entrez state.
    """
    request_parameters = _build_ncbi_efetch_request_parameters(
        database_name=database_name,
        retmode=retmode,
        rettype=rettype,
        protein_uids=protein_uids,
        history_web_env=history_web_env,
        history_query_key=history_query_key,
        retstart=retstart,
        retmax=retmax,
    )
    encoded_parameters = urlencode(request_parameters)
    requested_uid_count = len(protein_uids) if protein_uids is not None else 0
    use_post_request = (
        len(encoded_parameters) > 1000 or requested_uid_count >= 200
    )

    if reuse_http_connection:
        return _ncbi_persistent_https_transport.open_request(
            request_url=(
                NCBI_EFETCH_URL
                if use_post_request
                else f"{NCBI_EFETCH_URL}?{encoded_parameters}"
            ),
            http_method="POST" if use_post_request else "GET",
            request_body=(
                encoded_parameters.encode("utf-8") if use_post_request else None
            ),
            timeout_seconds=_active_ncbi_network_timeout_seconds.get(),
            ssl_context=_active_ncbi_ssl_context.get(),
        )

    if use_post_request:
        request = Request(
            NCBI_EFETCH_URL,
            data=encoded_parameters.encode("utf-8"),
            method="POST",
        )
    else:
        request = Request(
            f"{NCBI_EFETCH_URL}?{encoded_parameters}",
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

        # SSLEOFError and SSLZeroReturnError mean the TLS session ended while
        # data was still expected. They are transport failures, and reporting
        # them as trust-store problems sends the reader to configure a
        # certificate authority that was never involved.
        if isinstance(current_error, ssl.SSLError) and not isinstance(
            current_error,
            (ssl.SSLEOFError, ssl.SSLZeroReturnError),
        ):
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
    uid_retrieval_strategy: str = "esearch_history_efetch_uilist"
    esearch_request_count: int = 0
    efetch_request_count: int = 0
    fetch_timeout_seconds: Optional[float] = None
    request_deadline_seconds: Optional[float] = None
    retrieval_telemetry: Optional[Dict[str, Any]] = None


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


def _parse_ncbi_uid_list_payload(
    *,
    payload_bytes: bytes,
) -> List[str]:
    """
    Parse one newline-delimited EFetch ``uilist`` page into protein UIDs.
    """
    decoded_payload = payload_bytes.decode("utf-8", "replace")
    page_protein_uids: List[str] = []

    for payload_line in decoded_payload.splitlines():
        candidate_protein_uid = payload_line.strip()
        if not candidate_protein_uid:
            continue

        if not candidate_protein_uid.isdigit():
            raise NCBIXmlResponseValidationError(
                "NCBI returned a non-numeric value in a protein UID list page: "
                f"{candidate_protein_uid!r}."
            )

        page_protein_uids.append(candidate_protein_uid)

    return page_protein_uids


def fetch_ncbi_protein_uid_snapshot(
    *,
    ncbi_email: str,
    ncbi_api_key: Optional[str],
    query: str,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    deduplicate_uids: bool = True,
    sort_uids: bool = True,
    page_size: int = NCBI_MAX_EFETCH_RETMAX,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
    fetch_timeout_seconds: float = 30.0,
    request_deadline_seconds: float = 300.0,
    retry_backoff_initial_seconds: Optional[float] = None,
    retry_backoff_multiplier: float = 2.0,
    retry_backoff_max_seconds: float = 30.0,
    rate_limit_backoff_seconds: float = 5.0,
    reuse_http_connection: bool = False,
    verbose_logging: bool = False,
    telemetry: Optional[NCBIRetrievalTelemetry] = None,
) -> NCBIProteinUidFetchResult:
    """
    Retrieve protein UIDs from NCBI Entrez and return a snapshot-ready payload.

    The query is executed exactly once. The resulting set is retained on the
    NCBI History server and every page is then read back from that frozen set
    with EFetch. Re-running the query per page, as ``retstart`` pagination over
    ESearch does, offers no snapshot isolation: an index change between pages
    can duplicate or silently drop UIDs, and an omission is undetectable after
    the fact.

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
            Number of protein UIDs requested per History-backed EFetch call.
            EFetch paging is capped at 10,000 identifiers per request.
        max_retry_attempts:
            Maximum retry attempts for each failed paginated request.
        request_delay_seconds:
            Optional custom delay between successful requests.
            If None, defaults to:
            - 0.10 seconds when an API key is provided
            - 0.34 seconds otherwise
        fetch_timeout_seconds:
            Socket/network idle timeout applied to each UID request.
        request_deadline_seconds:
            Maximum total time allowed for one complete UID request, including
            retries and backoff delays.
        retry_backoff_initial_seconds:
            Delay before the first retry. If None, uses the resolved request
            delay, with a 0.10 second minimum.
        retry_backoff_multiplier:
            Multiplier applied to each subsequent retry delay.
        retry_backoff_max_seconds:
            Upper bound for any single retry backoff delay.
        rate_limit_backoff_seconds:
            Minimum delay applied before retrying a request that NCBI rejected
            for exceeding the request-rate allowance.
        reuse_http_connection:
            Reuse one keep-alive HTTPS connection instead of reconnecting per
            request.
        verbose_logging:
            Include request URLs and response headers in logs.
        telemetry:
            Optional shared telemetry collector. One is created when omitted.

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

    if page_size > NCBI_MAX_EFETCH_RETMAX:
        raise ValueError(
            "page_size must not exceed the EFetch paging ceiling of "
            f"{NCBI_MAX_EFETCH_RETMAX} identifiers per request."
        )

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

    resolved_retry_backoff_initial_seconds = (
        max(resolved_request_delay_seconds, 0.10)
        if retry_backoff_initial_seconds is None
        else retry_backoff_initial_seconds
    )

    _validate_xml_fetch_controls(
        fetch_timeout_seconds=fetch_timeout_seconds,
        batch_deadline_seconds=request_deadline_seconds,
        request_delay_seconds=resolved_request_delay_seconds,
        retry_backoff_initial_seconds=resolved_retry_backoff_initial_seconds,
        retry_backoff_multiplier=retry_backoff_multiplier,
        retry_backoff_max_seconds=retry_backoff_max_seconds,
        circuit_breaker_failure_threshold=1,
        circuit_breaker_cooldown_seconds=0.0,
        rate_limit_backoff_seconds=rate_limit_backoff_seconds,
    )

    active_telemetry = telemetry if telemetry is not None else NCBIRetrievalTelemetry()
    retrieved_at_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    esearch_request_count = 0
    efetch_request_count = 0

    with _configured_ncbi_entrez_urlopen(
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
    ) as (resolved_ssl_ca_file, resolved_ssl_ca_directory):
        request_controls = _NCBIRequestControls(
            max_retry_attempts=max_retry_attempts,
            fetch_timeout_seconds=fetch_timeout_seconds,
            request_deadline_seconds=request_deadline_seconds,
            request_delay_seconds=resolved_request_delay_seconds,
            retry_backoff_initial_seconds=resolved_retry_backoff_initial_seconds,
            retry_backoff_multiplier=retry_backoff_multiplier,
            retry_backoff_max_seconds=retry_backoff_max_seconds,
            rate_limit_backoff_seconds=rate_limit_backoff_seconds,
            resolved_ssl_ca_file=resolved_ssl_ca_file,
            resolved_ssl_ca_directory=resolved_ssl_ca_directory,
            telemetry=active_telemetry,
            verbose_logging=verbose_logging,
        )

        def open_initial_search_request():
            return Entrez.esearch(
                db=database_name,
                term=query,
                retmax=0,
                usehistory="y",
            )

        with active_telemetry.measure_stage(stage_name=UID_SEARCH_STAGE_NAME):
            initial_search_result = _execute_ncbi_request_with_controls(
                stage_name=UID_SEARCH_STAGE_NAME,
                request_context=f"initial ESearch for query {query!r}",
                operation_label="protein UID search",
                open_request=open_initial_search_request,
                controls=request_controls,
                permanent_failure_message_prefix=(
                    "Permanent NCBI protein UID search failure; not retrying"
                ),
                build_retry_exhausted_error=lambda error: RuntimeError(
                    "Failed to run the initial NCBI protein UID search after "
                    f"{max_retry_attempts} attempts: {error}"
                ),
            )
        esearch_request_count += 1

        initial_search_response = Entrez.read(
            io.BytesIO(initial_search_result.payload_bytes)
        )

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
        page_number = 0
        uid_progress_reporter = NCBIRetrievalProgressReporter(
            total_unit_count=ncbi_reported_result_count,
            unit_label="UIDs",
        )
        request_controls = replace(
            request_controls,
            progress_reporter=uid_progress_reporter,
        )

        with active_telemetry.measure_stage(stage_name=UID_PAGE_STAGE_NAME):
            while page_start_index < ncbi_reported_result_count:
                current_page_start_index = page_start_index
                page_number += 1
                page_context = (
                    f"protein UID page retstart={current_page_start_index}, "
                    f"retmax={page_size}"
                )
                parsed_page_protein_uids: List[List[str]] = []

                def open_uid_page_request():
                    return _open_ncbi_entrez_efetch_once(
                        database_name=database_name,
                        retmode="text",
                        rettype=NCBI_PROTEIN_UID_LIST_RETTYPE,
                        history_web_env=history_web_env,
                        history_query_key=history_query_key,
                        retstart=current_page_start_index,
                        retmax=page_size,
                        reuse_http_connection=reuse_http_connection,
                    )

                def validate_uid_page_payload(payload_bytes: bytes) -> None:
                    page_protein_uids = _parse_ncbi_uid_list_payload(
                        payload_bytes=payload_bytes,
                    )

                    if not page_protein_uids:
                        raise NCBITransientFetchError(
                            "NCBI returned an empty protein UID page before the "
                            "reported result set was fully paginated at "
                            f"retstart={current_page_start_index}."
                        )

                    if len(page_protein_uids) > page_size:
                        raise NCBIXmlResponseValidationError(
                            "NCBI returned more protein UIDs than requested for "
                            f"{page_context}. Expected at most {page_size}, got "
                            f"{len(page_protein_uids)}."
                        )

                    parsed_page_protein_uids.append(page_protein_uids)

                uid_page_result = _execute_ncbi_request_with_controls(
                    stage_name=UID_PAGE_STAGE_NAME,
                    request_context=page_context,
                    operation_label="protein UID page",
                    open_request=open_uid_page_request,
                    controls=request_controls,
                    permanent_failure_message_prefix=(
                        "Permanent NCBI protein UID page failure; not retrying"
                    ),
                    build_retry_exhausted_error=lambda error: RuntimeError(
                        "Failed to extract protein UIDs after "
                        f"{max_retry_attempts} attempts at "
                        f"page_start_index={current_page_start_index}: {error}"
                    ),
                    validate_payload=validate_uid_page_payload,
                )
                efetch_request_count += 1

                current_page_protein_uids = parsed_page_protein_uids[-1]
                raw_protein_uids.extend(current_page_protein_uids)
                page_start_index += len(current_page_protein_uids)

                uid_progress_reporter.report_completed_units(
                    completed_unit_delta=len(current_page_protein_uids),
                    item_label=f"page {page_number}",
                    response_byte_count=len(uid_page_result.payload_bytes),
                    round_trip_latency_seconds=(
                        uid_page_result.round_trip_latency_seconds
                    ),
                )

        uid_progress_reporter.finish()
        uid_retrieval_seconds = uid_progress_reporter.elapsed_seconds
        print(
            f"Retrieved {len(raw_protein_uids):,} protein UIDs in "
            f"{_format_duration_seconds(uid_retrieval_seconds)}."
        )
        if uid_retrieval_seconds > 0 and uid_progress_reporter.network_byte_count:
            print(
                "Transferred "
                f"{uid_progress_reporter.network_byte_count / (1024 * 1024):.1f} MB "
                "at "
                f"{uid_progress_reporter.network_byte_count / 1024 / uid_retrieval_seconds:,.0f}"
                " KB/s."
            )

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
            uid_retrieval_strategy="esearch_history_efetch_uilist",
            esearch_request_count=esearch_request_count,
            efetch_request_count=efetch_request_count,
            fetch_timeout_seconds=fetch_timeout_seconds,
            request_deadline_seconds=request_deadline_seconds,
            retrieval_telemetry=active_telemetry.build_summary(),
        )

        print(
            f"Issued {esearch_request_count} ESearch and "
            f"{efetch_request_count} EFetch requests."
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

    A batch carries either its raw payload in memory or the path of the file it
    was streamed to. Streaming persistence is the default because retaining
    every payload until consolidation is what makes peak memory scale with the
    size of the whole result set rather than with one batch.
    """
    batch_index: int
    batch_start_index: int
    batch_end_index: int
    protein_uids: List[str]
    protein_uid_count: int
    xml_payload_sha256: str
    xml_payload_bytes: Optional[bytes] = None
    xml_payload_file_path: Optional[str] = None
    record_count: Optional[int] = None
    response_byte_count: Optional[int] = None
    round_trip_latency_seconds: Optional[float] = None
    attempt_count: Optional[int] = None
    reused_from_workspace: bool = False


class NCBITransientFetchError(RuntimeError):
    """
    Error raised when a fetch attempt failed in a way that may recover on retry.
    """


class NCBIRateLimitError(NCBITransientFetchError):
    """
    Error raised when NCBI reports that the request rate allowance was exceeded.

    NCBI answers a rate-limit violation with a small JSON body that may arrive
    with HTTP 200. Without explicit detection that body fails XML parsing and is
    classified as a permanent failure, which would abort an otherwise healthy
    run at exactly the moment concurrency is introduced.
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
    rettype: Optional[str] = None
    rate_limit_backoff_seconds: Optional[float] = None
    max_concurrent_requests: int = 1
    max_request_starts_per_second: Optional[float] = None
    reuse_http_connection: bool = False
    batch_workspace_directory: Optional[str] = None
    batch_workspace_is_ephemeral: bool = False
    reused_batch_count: int = 0
    fetched_batch_count: int = 0
    retrieval_telemetry: Optional[Dict[str, Any]] = None


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
    rate_limit_backoff_seconds: float = 0.0,
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
    _validate_non_negative_float(
        parameter_name="rate_limit_backoff_seconds",
        parameter_value=rate_limit_backoff_seconds,
    )


def _build_xml_batch_context(
    *,
    batch_index: int,
    total_batch_count: int,
    batch_start_index: int,
    current_batch_protein_uids: Sequence[str],
    include_full_protein_uid_list: bool = False,
) -> str:
    """
    Describe one XML batch compactly enough to appear in every log line.

    Embedding the whole UID list made a single failure message hundreds of
    identifiers long and grew linearly with batch size. The bounded summary
    keeps the identifying range, and the full list stays available behind the
    verbose flag when a specific batch has to be reconstructed by hand.
    """
    batch_end_index = batch_start_index + len(current_batch_protein_uids) - 1

    if include_full_protein_uid_list:
        protein_uid_summary = f"{list(current_batch_protein_uids)}"
    elif current_batch_protein_uids:
        protein_uid_summary = (
            f"{current_batch_protein_uids[0]}..{current_batch_protein_uids[-1]} "
            f"(n={len(current_batch_protein_uids)})"
        )
    else:
        protein_uid_summary = "(n=0)"

    return (
        f"batch_index={batch_index}/{total_batch_count}, "
        f"batch_start_index={batch_start_index}, "
        f"batch_end_index={batch_end_index}, "
        f"protein_uids={protein_uid_summary}"
    )


def _is_transient_ncbi_fetch_error(error: Exception) -> bool:
    if isinstance(error, NCBITransientFetchError):
        return True

    if isinstance(error, HTTPError):
        return error.code in {408, 429, 500, 502, 503, 504}

    # A response that ends early is the textbook retryable failure: the request
    # was accepted and the server began answering it. http.client reports these
    # as HTTPException, which is not an OSError and carries no wording the text
    # heuristics below would recognize, so without this branch a truncated
    # response is classified permanent and aborts an otherwise healthy run.
    if isinstance(error, http.client.HTTPException):
        return True

    # ConnectionResetError, ConnectionAbortedError, and BrokenPipeError are
    # OSError subclasses rather than URLError, and their messages do not
    # reliably contain the words matched below.
    if isinstance(error, ConnectionError):
        return True

    # A TLS session that ends mid-stream is a transport failure. Trust-store
    # misconfiguration is classified earlier and separately.
    if isinstance(error, (ssl.SSLEOFError, ssl.SSLZeroReturnError)):
        return True

    if isinstance(error, (TimeoutError, socket.timeout, URLError)):
        return True

    error_message = str(error).lower()
    return (
        "timed out" in error_message
        or "timeout" in error_message
        or "temporarily unavailable" in error_message
        or "connection reset" in error_message
        or "connection aborted" in error_message
        or "incompleteread" in error_message
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


def _payload_indicates_ncbi_rate_limit(
    *,
    payload_bytes: bytes,
) -> bool:
    """
    Return True when a response body is an NCBI rate-limit notice.

    The documented shape is ``{"error":"API rate limit exceeded","count":"11"}``
    delivered instead of the requested payload. The body is matched both as raw
    text and as JSON so that a reworded message with the same structure is still
    recognized.
    """
    inspected_payload = payload_bytes[:4096]

    for rate_limit_marker in _NCBI_RATE_LIMIT_BODY_MARKERS:
        if rate_limit_marker in inspected_payload:
            return True

    stripped_payload = inspected_payload.strip()
    if not stripped_payload.startswith(b"{"):
        return False

    try:
        decoded_payload = json.loads(payload_bytes.decode("utf-8", "replace"))
    except Exception:
        return False

    if not isinstance(decoded_payload, dict):
        return False

    reported_error = decoded_payload.get("error")
    if not isinstance(reported_error, str):
        return False

    normalized_error = reported_error.lower()
    return "rate limit" in normalized_error or "too many requests" in normalized_error


def _classify_ncbi_failure_kind(error: BaseException) -> str:
    """
    Map one failed attempt onto a telemetry failure category.
    """
    if isinstance(error, NCBIRateLimitError):
        return RATE_LIMIT_FAILURE_KIND

    if isinstance(error, NCBIXmlBatchDeadlineExceeded):
        return DEADLINE_FAILURE_KIND

    if isinstance(error, NCBIXmlResponseValidationError):
        return VALIDATION_FAILURE_KIND

    if isinstance(error, HTTPError):
        if error.code == 429:
            return HTTP_429_FAILURE_KIND
        if 500 <= error.code < 600:
            return HTTP_5XX_FAILURE_KIND
        return OTHER_FAILURE_KIND

    if isinstance(error, (TimeoutError, socket.timeout)):
        return TIMEOUT_FAILURE_KIND

    if isinstance(
        error,
        (
            http.client.HTTPException,
            ConnectionError,
            ssl.SSLEOFError,
            ssl.SSLZeroReturnError,
        ),
    ):
        return TRUNCATED_RESPONSE_FAILURE_KIND

    error_message = str(error).lower()
    if "timed out" in error_message or "timeout" in error_message:
        return TIMEOUT_FAILURE_KIND

    if "incompleteread" in error_message:
        return TRUNCATED_RESPONSE_FAILURE_KIND

    return OTHER_FAILURE_KIND


def _failure_kind_indicates_overload(failure_kind: str) -> bool:
    """
    Return True for failures that should reduce request pressure on NCBI.
    """
    return failure_kind in {
        RATE_LIMIT_FAILURE_KIND,
        HTTP_429_FAILURE_KIND,
        HTTP_5XX_FAILURE_KIND,
        TIMEOUT_FAILURE_KIND,
        TRUNCATED_RESPONSE_FAILURE_KIND,
        DEADLINE_FAILURE_KIND,
    }


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
    operation_label: str = "XML batch",
) -> None:
    message = (
        f"NCBI {operation_label} deadline exhausted after "
        f"{elapsed_seconds:.3f} seconds with a configured deadline of "
        f"{batch_deadline_seconds:.3f} seconds for {batch_context}."
    )
    if last_error is not None:
        message += f" Last error: {last_error}"

    raise NCBIXmlBatchDeadlineExceeded(message)


@dataclass(frozen=True)
class NCBIFetchAttemptOutcome:
    """
    One completed request attempt with its transport timing breakdown.

    ``time_to_first_byte_seconds`` covers everything up to the response headers
    being available, which includes connection setup. ``connect_seconds`` is
    reported separately only when the transport in use can measure it, so the
    two can be compared to tell "NCBI is still rendering" apart from "the link
    is slow".
    """
    payload: Any
    response_url: Any
    response_headers: Any
    connect_seconds: Optional[float] = None
    time_to_first_byte_seconds: Optional[float] = None
    response_read_seconds: Optional[float] = None


def _fetch_ncbi_payload_attempt_with_deadline(
    *,
    open_request: Callable[[], Any],
    effective_fetch_timeout_seconds: float,
    request_deadline_at_seconds: float,
    request_started_at_seconds: float,
    request_deadline_seconds: float,
    request_context: str,
) -> NCBIFetchAttemptOutcome:
    """
    Open, read, and close one EFetch response within an absolute deadline.

    ``urllib`` timeouts bound individual blocking network operations, but they
    do not provide a scheduler-level deadline for the complete response. The
    request therefore runs in a daemon worker while the calling thread waits
    only until the remaining request deadline. On deadline exhaustion, the
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
                request_opened_at_seconds = time.perf_counter()
                fetch_handle = open_request()
                time_to_first_byte_seconds = (
                    time.perf_counter() - request_opened_at_seconds
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
                    attempt_state["connect_seconds"] = getattr(
                        fetch_handle,
                        "ncbi_connect_seconds",
                        None,
                    )
                    attempt_state["time_to_first_byte_seconds"] = (
                        time_to_first_byte_seconds
                    )

                if cancellation_requested.is_set():
                    return

                response_read_started_at_seconds = time.perf_counter()
                raw_xml_payload = fetch_handle.read()
                with attempt_state_lock:
                    attempt_state["raw_xml_payload"] = raw_xml_payload
                    attempt_state["response_read_seconds"] = (
                        time.perf_counter() - response_read_started_at_seconds
                    )
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

    remaining_request_seconds = (
        request_deadline_at_seconds - time.monotonic()
    )
    attempt_finished_before_deadline = (
        remaining_request_seconds > 0
        and attempt_completed.wait(timeout=remaining_request_seconds)
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
            batch_context=request_context,
            batch_deadline_seconds=request_deadline_seconds,
            elapsed_seconds=(
                deadline_detected_at_seconds - request_started_at_seconds
            ),
        )

    with attempt_state_lock:
        attempt_error = attempt_state.get("error")
        raw_xml_payload = attempt_state.get("raw_xml_payload")
        response_url = attempt_state.get("response_url")
        response_headers = attempt_state.get("response_headers")
        connect_seconds = attempt_state.get("connect_seconds")
        time_to_first_byte_seconds = attempt_state.get(
            "time_to_first_byte_seconds"
        )
        response_read_seconds = attempt_state.get("response_read_seconds")

    if attempt_error is not None:
        raise attempt_error

    return NCBIFetchAttemptOutcome(
        payload=raw_xml_payload,
        response_url=response_url,
        response_headers=response_headers,
        connect_seconds=connect_seconds,
        time_to_first_byte_seconds=time_to_first_byte_seconds,
        response_read_seconds=response_read_seconds,
    )


@dataclass(frozen=True)
class _NCBIRequestControls:
    """
    The failure-control policy shared by every NCBI request in one run.
    """
    max_retry_attempts: int
    fetch_timeout_seconds: float
    request_deadline_seconds: float
    request_delay_seconds: float
    retry_backoff_initial_seconds: float
    retry_backoff_multiplier: float
    retry_backoff_max_seconds: float
    rate_limit_backoff_seconds: float
    resolved_ssl_ca_file: Optional[str]
    resolved_ssl_ca_directory: Optional[str]
    telemetry: Optional[NCBIRetrievalTelemetry] = None
    verbose_logging: bool = False
    progress_reporter: Optional["NCBIRetrievalProgressReporter"] = None

    def write_log_line(self, message: str) -> None:
        """
        Emit a message without corrupting an in-place progress line.
        """
        if self.progress_reporter is not None:
            self.progress_reporter.write_line_above_progress(message)
            return

        print(message)


@dataclass(frozen=True)
class _NCBIControlledRequestResult:
    """
    One successfully retrieved and validated NCBI response.
    """
    payload_bytes: bytes
    payload_sha256: str
    attempt_count: int
    round_trip_latency_seconds: float
    response_url: Any
    response_headers: Any


def _execute_ncbi_request_with_controls(
    *,
    stage_name: str,
    request_context: str,
    operation_label: str,
    open_request: Callable[[], Any],
    controls: _NCBIRequestControls,
    permanent_failure_message_prefix: str,
    build_retry_exhausted_error: Callable[[Exception], Exception],
    validate_payload: Optional[Callable[[bytes], None]] = None,
    circuit_breaker: Optional[NCBIXmlCircuitBreaker] = None,
    rate_limiter: Optional[NCBIRequestStartRateLimiter] = None,
    concurrency_governor: Optional[NCBIAdaptiveConcurrencyGovernor] = None,
    on_attempt_start: Optional[Callable[[int], None]] = None,
    on_response_received: Optional[
        Callable[[NCBIFetchAttemptOutcome, bytes, float], None]
    ] = None,
    on_response_validated: Optional[Callable[[], None]] = None,
) -> _NCBIControlledRequestResult:
    """
    Run one NCBI request under the full timeout, deadline, and retry policy.

    Both retrieval stages share this loop so that the UID path cannot silently
    drift away from the failure controls the XML path already enforces. Stage
    specific behavior enters only through the injected request opener, payload
    validator, and logging callbacks.
    """
    telemetry = controls.telemetry
    request_started_at_seconds = time.monotonic()
    request_deadline_at_seconds = (
        request_started_at_seconds + controls.request_deadline_seconds
    )

    if circuit_breaker is not None:
        circuit_breaker.assert_request_allowed(
            current_time_seconds=request_started_at_seconds,
            batch_context=request_context,
        )

    last_transient_error: Optional[Exception] = None

    for retry_attempt_index in range(controls.max_retry_attempts):
        try:
            attempt_started_at_seconds = time.monotonic()
            if attempt_started_at_seconds >= request_deadline_at_seconds:
                _raise_xml_batch_deadline_exceeded(
                    batch_context=request_context,
                    batch_deadline_seconds=controls.request_deadline_seconds,
                    elapsed_seconds=(
                        attempt_started_at_seconds - request_started_at_seconds
                    ),
                    last_error=last_transient_error,
                    operation_label=operation_label,
                )

            remaining_request_seconds = (
                request_deadline_at_seconds - attempt_started_at_seconds
            )
            effective_fetch_timeout_seconds = min(
                controls.fetch_timeout_seconds,
                remaining_request_seconds,
            )

            if on_attempt_start is not None:
                on_attempt_start(retry_attempt_index)

            request_started_at = time.perf_counter()
            with (
                concurrency_governor.request_slot()
                if concurrency_governor is not None
                else _null_request_slot()
            ):
                if rate_limiter is not None:
                    rate_limiter_wait_seconds = (
                        rate_limiter.acquire_request_start_slot()
                    )
                    if rate_limiter_wait_seconds and telemetry is not None:
                        telemetry.record_sleep(
                            stage_name=stage_name,
                            sleep_seconds=rate_limiter_wait_seconds,
                        )
                    request_started_at = time.perf_counter()

                attempt_outcome = _fetch_ncbi_payload_attempt_with_deadline(
                    open_request=open_request,
                    effective_fetch_timeout_seconds=(
                        effective_fetch_timeout_seconds
                    ),
                    request_deadline_at_seconds=request_deadline_at_seconds,
                    request_started_at_seconds=request_started_at_seconds,
                    request_deadline_seconds=controls.request_deadline_seconds,
                    request_context=request_context,
                )

            round_trip_latency_seconds = time.perf_counter() - request_started_at
            attempt_finished_at_seconds = time.monotonic()
            if attempt_finished_at_seconds > request_deadline_at_seconds:
                _raise_xml_batch_deadline_exceeded(
                    batch_context=request_context,
                    batch_deadline_seconds=controls.request_deadline_seconds,
                    elapsed_seconds=(
                        attempt_finished_at_seconds - request_started_at_seconds
                    ),
                    last_error=last_transient_error,
                    operation_label=operation_label,
                )

            payload_bytes = _coerce_payload_to_bytes(
                payload=attempt_outcome.payload,
            )

            if telemetry is not None:
                telemetry.record_request(
                    stage_name=stage_name,
                    total_seconds=round_trip_latency_seconds,
                    response_byte_count=len(payload_bytes),
                    connect_seconds=attempt_outcome.connect_seconds,
                    time_to_first_byte_seconds=(
                        attempt_outcome.time_to_first_byte_seconds
                    ),
                    response_read_seconds=attempt_outcome.response_read_seconds,
                )

            if not payload_bytes.strip():
                raise NCBITransientFetchError(
                    f"NCBI returned an empty payload for {operation_label}."
                )

            if _payload_indicates_ncbi_rate_limit(payload_bytes=payload_bytes):
                raise NCBIRateLimitError(
                    "NCBI reported an API rate-limit violation in the response "
                    f"body for {request_context}."
                )

            payload_sha256 = hashlib.sha256(payload_bytes).hexdigest()

            if on_response_received is not None:
                on_response_received(
                    attempt_outcome,
                    payload_bytes,
                    round_trip_latency_seconds,
                )

            if _xml_payload_contains_server_error_marker(
                xml_payload_bytes=payload_bytes,
            ):
                raise NCBITransientFetchError(
                    "NCBI XML payload contains server-side error markers "
                    "despite a successful HTTP response."
                )

            if validate_payload is not None:
                validate_payload(payload_bytes)

            if circuit_breaker is not None:
                circuit_breaker.record_success()

            if concurrency_governor is not None:
                concurrency_governor.record_successful_request()

            if on_response_validated is not None:
                on_response_validated()

            if rate_limiter is None:
                jitter_seconds = random.uniform(0.01, 0.05)
                regulatory_sleep_seconds = (
                    controls.request_delay_seconds + jitter_seconds
                )
                time.sleep(regulatory_sleep_seconds)
                if telemetry is not None:
                    telemetry.record_sleep(
                        stage_name=stage_name,
                        sleep_seconds=regulatory_sleep_seconds,
                    )

            return _NCBIControlledRequestResult(
                payload_bytes=payload_bytes,
                payload_sha256=payload_sha256,
                attempt_count=retry_attempt_index + 1,
                round_trip_latency_seconds=round_trip_latency_seconds,
                response_url=attempt_outcome.response_url,
                response_headers=attempt_outcome.response_headers,
            )

        except Exception as error:
            if _is_ncbi_tls_configuration_error(error):
                raise RuntimeError(
                    _build_ncbi_tls_configuration_error_message(
                        error=error,
                        ssl_ca_file=controls.resolved_ssl_ca_file,
                        ssl_ca_directory=controls.resolved_ssl_ca_directory,
                    )
                ) from error

            failure_kind = _classify_ncbi_failure_kind(error)
            if telemetry is not None:
                telemetry.record_failure(
                    stage_name=stage_name,
                    failure_kind=failure_kind,
                )

            if concurrency_governor is not None and _failure_kind_indicates_overload(
                failure_kind,
            ):
                concurrency_governor.record_overload_signal()

            if isinstance(error, NCBIXmlBatchDeadlineExceeded):
                if circuit_breaker is not None:
                    circuit_breaker.record_failure(
                        current_time_seconds=time.monotonic(),
                    )
                raise error

            if not _is_transient_ncbi_fetch_error(error):
                raise RuntimeError(
                    f"{permanent_failure_message_prefix} "
                    f"{request_context}: {error}"
                ) from error

            last_transient_error = error

            if retry_attempt_index == controls.max_retry_attempts - 1:
                if circuit_breaker is not None:
                    circuit_breaker.record_failure(
                        current_time_seconds=time.monotonic(),
                    )
                raise build_retry_exhausted_error(error) from error

            retry_backoff_seconds = _calculate_retry_backoff_seconds(
                retry_attempt_index=retry_attempt_index,
                retry_backoff_initial_seconds=(
                    controls.retry_backoff_initial_seconds
                ),
                retry_backoff_multiplier=controls.retry_backoff_multiplier,
                retry_backoff_max_seconds=controls.retry_backoff_max_seconds,
            )
            if isinstance(error, NCBIRateLimitError):
                # A rate-limit response means the allowance is already spent.
                # Backing off by the ordinary transient delay would re-enter the
                # same violation, so the wait is extended explicitly.
                retry_backoff_seconds = max(
                    retry_backoff_seconds,
                    controls.rate_limit_backoff_seconds,
                )

            retry_scheduled_at_seconds = time.monotonic()
            remaining_request_seconds = (
                request_deadline_at_seconds - retry_scheduled_at_seconds
            )
            if retry_backoff_seconds >= remaining_request_seconds:
                if circuit_breaker is not None:
                    circuit_breaker.record_failure(
                        current_time_seconds=retry_scheduled_at_seconds,
                    )
                _raise_xml_batch_deadline_exceeded(
                    batch_context=request_context,
                    batch_deadline_seconds=controls.request_deadline_seconds,
                    elapsed_seconds=(
                        retry_scheduled_at_seconds - request_started_at_seconds
                    ),
                    last_error=error,
                    operation_label=operation_label,
                )

            controls.write_log_line(
                f"Retry {retry_attempt_index + 1}/{controls.max_retry_attempts} "
                f"for {request_context} after transient error "
                f"({failure_kind}): {error}"
            )
            if controls.progress_reporter is not None:
                controls.progress_reporter.record_retry()
            if telemetry is not None:
                telemetry.record_retry(stage_name=stage_name)
                telemetry.record_sleep(
                    stage_name=stage_name,
                    sleep_seconds=retry_backoff_seconds,
                )
            time.sleep(retry_backoff_seconds)

    raise build_retry_exhausted_error(
        last_transient_error
        if last_transient_error is not None
        else RuntimeError("no attempt was executed")
    )


@contextmanager
def _null_request_slot():
    yield


def _format_duration_seconds(total_seconds: float) -> str:
    """
    Format a duration for a human watching a long retrieval scroll past.
    """
    if total_seconds < 60:
        return f"{total_seconds:.1f}s"

    total_whole_seconds = int(round(total_seconds))
    elapsed_hours, remaining_seconds = divmod(total_whole_seconds, 3600)
    elapsed_minutes, elapsed_seconds = divmod(remaining_seconds, 60)

    if elapsed_hours:
        return f"{elapsed_hours}h{elapsed_minutes:02d}m{elapsed_seconds:02d}s"

    return f"{elapsed_minutes}m{elapsed_seconds:02d}s"


PROGRESS_BAR_WIDTH = 24


def _resolve_progress_bar_characters() -> tuple[str, str]:
    """
    Pick bar glyphs the active output stream can actually encode.
    """
    output_encoding = getattr(sys.stdout, "encoding", None) or "ascii"

    try:
        "█░".encode(output_encoding)
    except (LookupError, UnicodeEncodeError):
        return "#", "-"

    return "█", "░"


def _supports_inline_progress_updates() -> bool:
    """
    Return True when the output stream can be rewritten in place.

    A notebook front end accepts a carriage return even though it is not a
    terminal, while a redirected stream would only accumulate control
    characters, so both are detected rather than assumed.
    """
    try:
        if sys.stdout.isatty():
            return True
    except Exception:
        return False

    return "ipykernel" in sys.modules


class NCBIRetrievalProgressReporter:
    """
    Render one self-rewriting progress line for a whole retrieval stage.

    A long retrieval is a decision point while it is still running: the operator
    needs to know whether to let it finish. Under concurrency, per-request log
    lines from several workers interleave and arrive out of order, which buries
    that answer in noise, so the per-request lines move behind the verbose flag
    and one line carries the stage: how far along it is, what the unit that just
    landed cost, the aggregate transfer rate, and the projected finish time.

    Progress is counted in whatever unit the stage advances by. The XML stage
    advances one batch at a time; the UID stage advances by however many
    identifiers a History page returned. Both share the same line so a reader
    does not have to learn two formats.

    The projection is derived from completion throughput, not from per-request
    latency, because the two diverge as soon as requests overlap. Units that did
    not cost a request, such as batches reused from the workspace, are excluded
    from the rate so the projection cannot collapse to zero on a resumed run.
    """

    # Below this, elapsed time is too short to divide by: the derived rate and
    # the projection would be dominated by measurement noise and would print
    # numbers that are confidently wrong.
    MINIMUM_ELAPSED_SECONDS_FOR_ESTIMATES = 1.0

    def __init__(
        self,
        *,
        total_unit_count: int,
        unit_label: str = "batches",
        inline_updates_enabled: Optional[bool] = None,
        monotonic_function=time.monotonic,
    ) -> None:
        self._total_unit_count = max(1, total_unit_count)
        self._unit_label = unit_label
        self._state_lock = threading.RLock()
        self._monotonic_function = monotonic_function
        self._started_at_seconds = monotonic_function()
        self._completed_unit_count = 0
        self._network_unit_count = 0
        self._network_byte_count = 0
        self._retry_count = 0
        self._last_rendered_line = ""
        self._last_rendered_length = 0
        self._inline_updates_enabled = (
            _supports_inline_progress_updates()
            if inline_updates_enabled is None
            else inline_updates_enabled
        )
        self._filled_bar_character, self._empty_bar_character = (
            _resolve_progress_bar_characters()
        )

    @property
    def elapsed_seconds(self) -> float:
        return self._monotonic_function() - self._started_at_seconds

    @property
    def network_byte_count(self) -> int:
        with self._state_lock:
            return self._network_byte_count

    def _render_progress_bar(self, *, completion_fraction: float) -> str:
        filled_segment_count = int(completion_fraction * PROGRESS_BAR_WIDTH)
        filled_segment_count = max(0, min(PROGRESS_BAR_WIDTH, filled_segment_count))
        return (
            self._filled_bar_character * filled_segment_count
            + self._empty_bar_character
            * (PROGRESS_BAR_WIDTH - filled_segment_count)
        )

    def _build_progress_line(
        self,
        *,
        completed_unit_count: int,
        network_unit_count: int,
        network_byte_count: int,
        retry_count: int,
        elapsed_seconds: float,
        completed_item_summary: str,
    ) -> str:
        completion_fraction = min(
            1.0,
            completed_unit_count / self._total_unit_count,
        )
        remaining_unit_count = self._total_unit_count - completed_unit_count
        unit_count_width = len(f"{self._total_unit_count:,}")
        estimates_are_meaningful = (
            elapsed_seconds >= self.MINIMUM_ELAPSED_SECONDS_FOR_ESTIMATES
        )

        progress_line_segments = [
            f"[{self._render_progress_bar(completion_fraction=completion_fraction)}] "
            f"{completed_unit_count:>{unit_count_width},}/"
            f"{self._total_unit_count:,} {completion_fraction:>6.1%}",
            completed_item_summary,
        ]

        if network_byte_count:
            transferred_megabytes = network_byte_count / (1024 * 1024)
            if estimates_are_meaningful:
                progress_line_segments.append(
                    f"{transferred_megabytes:.1f} MB at "
                    f"{network_byte_count / 1024 / elapsed_seconds:,.0f} KB/s"
                )
            else:
                progress_line_segments.append(f"{transferred_megabytes:.1f} MB")

        if retry_count:
            progress_line_segments.append(f"{retry_count} retries")

        progress_line_segments.append(
            f"{_format_duration_seconds(elapsed_seconds)} elapsed"
        )

        if (
            network_unit_count
            and estimates_are_meaningful
            and remaining_unit_count > 0
        ):
            projected_remaining_seconds = remaining_unit_count * (
                elapsed_seconds / network_unit_count
            )
            projected_finish_at = datetime.now().astimezone() + timedelta(
                seconds=projected_remaining_seconds,
            )
            progress_line_segments.append(
                f"~{_format_duration_seconds(projected_remaining_seconds)} left"
            )
            progress_line_segments.append(
                f"ends {projected_finish_at.strftime('%H:%M:%S')}"
            )

        return " | ".join(progress_line_segments)

    def _emit(self, progress_line: str) -> None:
        if not self._inline_updates_enabled:
            print(progress_line)
            return

        erase_padding = " " * max(
            0,
            self._last_rendered_length - len(progress_line),
        )
        print(f"\r{progress_line}{erase_padding}", end="", flush=True)
        self._last_rendered_line = progress_line
        self._last_rendered_length = len(progress_line)

    def record_retry(self) -> None:
        with self._state_lock:
            self._retry_count += 1

    def report_completed_units(
        self,
        *,
        completed_unit_delta: int,
        item_label: str,
        response_byte_count: Optional[int] = None,
        round_trip_latency_seconds: Optional[float] = None,
        required_a_request: bool = True,
    ) -> None:
        """
        Record one completed unit of work and redraw the progress line.
        """
        with self._state_lock:
            self._completed_unit_count += completed_unit_delta
            if required_a_request:
                self._network_unit_count += completed_unit_delta
                self._network_byte_count += int(response_byte_count or 0)

            completed_item_summary = item_label
            if not required_a_request:
                completed_item_summary += " reused"
            else:
                if round_trip_latency_seconds is not None:
                    completed_item_summary += f" {round_trip_latency_seconds:.1f}s"
                if response_byte_count:
                    completed_item_summary += (
                        f" {response_byte_count / (1024 * 1024):.2f} MB"
                    )

            self._emit(
                self._build_progress_line(
                    completed_unit_count=self._completed_unit_count,
                    network_unit_count=self._network_unit_count,
                    network_byte_count=self._network_byte_count,
                    retry_count=self._retry_count,
                    elapsed_seconds=self.elapsed_seconds,
                    completed_item_summary=completed_item_summary,
                )
            )

    def write_line_above_progress(self, message: str) -> None:
        """
        Emit a standalone message without leaving the progress line corrupted.
        """
        with self._state_lock:
            if not self._inline_updates_enabled or not self._last_rendered_length:
                print(message)
                return

            print(
                "\r" + " " * self._last_rendered_length + "\r",
                end="",
            )
            print(message)
            print(self._last_rendered_line, end="", flush=True)

    def finish(self) -> None:
        """
        Terminate the progress line so later output starts on its own row.
        """
        with self._state_lock:
            if self._inline_updates_enabled and self._last_rendered_length:
                print(flush=True)
                self._last_rendered_length = 0
                self._last_rendered_line = ""


def read_ncbi_xml_batch_payload_bytes(
    *,
    xml_batch: NCBIProteinXmlBatchFetchResult,
) -> bytes:
    """
    Return the raw payload of one fetched batch from memory or from disk.
    """
    if xml_batch.xml_payload_bytes is not None:
        return xml_batch.xml_payload_bytes

    if xml_batch.xml_payload_file_path is None:
        raise ValueError(
            "XML batch "
            f"{xml_batch.batch_index} carries neither an in-memory payload nor "
            "a persisted payload file path."
        )

    return Path(xml_batch.xml_payload_file_path).read_bytes()


def _resolve_max_request_starts_per_second(
    *,
    max_request_starts_per_second: Optional[float],
    ncbi_api_key: Optional[str],
) -> float:
    if max_request_starts_per_second is not None:
        return max_request_starts_per_second

    return (
        NCBI_MAX_REQUEST_STARTS_PER_SECOND_WITH_API_KEY
        if ncbi_api_key
        else NCBI_MAX_REQUEST_STARTS_PER_SECOND_WITHOUT_API_KEY
    )


def _fetch_single_ncbi_xml_batch(
    *,
    planned_batch: PlannedXmlBatch,
    total_batch_count: int,
    database_name: str,
    retmode: str,
    rettype: Optional[str],
    controls: _NCBIRequestControls,
    circuit_breaker: NCBIXmlCircuitBreaker,
    workspace: Optional[NCBIXmlBatchWorkspace],
    rate_limiter: Optional[NCBIRequestStartRateLimiter],
    concurrency_governor: Optional[NCBIAdaptiveConcurrencyGovernor],
    validate_batch_protein_uids: bool,
    reuse_http_connection: bool,
    progress_reporter: Optional[NCBIRetrievalProgressReporter] = None,
) -> NCBIProteinXmlBatchFetchResult:
    """
    Retrieve, validate, and persist one XML batch, reusing prior work if valid.
    """
    telemetry = controls.telemetry
    batch_context = _build_xml_batch_context(
        batch_index=planned_batch.batch_index,
        total_batch_count=total_batch_count,
        batch_start_index=planned_batch.batch_start_index,
        current_batch_protein_uids=planned_batch.protein_uids,
        include_full_protein_uid_list=controls.verbose_logging,
    )

    if workspace is not None and workspace.resume_enabled:
        reusable_batch = workspace.load_reusable_batch(
            planned_batch=planned_batch,
            validate_protein_uids=validate_batch_protein_uids,
        )
        if reusable_batch is not None:
            if controls.verbose_logging:
                controls.write_log_line(
                    f"Reusing validated XML batch {planned_batch.batch_index}/"
                    f"{total_batch_count} from the batch workspace."
                )
            if telemetry is not None:
                telemetry.record_reused_batch(stage_name=XML_BATCH_STAGE_NAME)

            if progress_reporter is not None:
                progress_reporter.report_completed_units(
                    completed_unit_delta=1,
                    item_label=f"#{planned_batch.batch_index}",
                    response_byte_count=reusable_batch.response_byte_count,
                    required_a_request=False,
                )

            return NCBIProteinXmlBatchFetchResult(
                batch_index=planned_batch.batch_index,
                batch_start_index=planned_batch.batch_start_index,
                batch_end_index=planned_batch.batch_end_index,
                protein_uids=planned_batch.protein_uids,
                protein_uid_count=planned_batch.protein_uid_count,
                xml_payload_sha256=reusable_batch.xml_payload_sha256,
                xml_payload_file_path=str(reusable_batch.payload_file_path),
                record_count=reusable_batch.record_count,
                response_byte_count=reusable_batch.response_byte_count,
                round_trip_latency_seconds=(
                    reusable_batch.round_trip_latency_seconds
                ),
                attempt_count=reusable_batch.attempt_count,
                reused_from_workspace=True,
            )

    validated_record_count: Dict[str, int] = {}

    def open_batch_request():
        return _open_ncbi_entrez_efetch_once(
            database_name=database_name,
            protein_uids=planned_batch.protein_uids,
            retmode=retmode,
            rettype=rettype,
            reuse_http_connection=reuse_http_connection,
        )

    def announce_attempt(retry_attempt_index: int) -> None:
        if not controls.verbose_logging:
            return

        controls.write_log_line(
            f"Starting XML request for batch {planned_batch.batch_index}/"
            f"{total_batch_count} "
            f"(attempt {retry_attempt_index + 1}/{controls.max_retry_attempts}) "
            f"with {planned_batch.protein_uid_count} protein UIDs."
        )

    def announce_response(
        attempt_outcome: NCBIFetchAttemptOutcome,
        payload_bytes: bytes,
        round_trip_latency_seconds: float,
    ) -> None:
        if not controls.verbose_logging:
            return

        controls.write_log_line(
            f"Received XML response for batch {planned_batch.batch_index}/"
            f"{total_batch_count} in {round_trip_latency_seconds:.3f} seconds "
            f"({len(payload_bytes)} bytes)."
        )

        if attempt_outcome.response_url is not None:
            controls.write_log_line(
                f"NCBI XML request URL: {attempt_outcome.response_url}"
            )

        if attempt_outcome.response_headers is not None:
            controls.write_log_line(
                "NCBI XML response headers: "
                f"{dict(attempt_outcome.response_headers)}"
            )

    def validate_batch_payload(payload_bytes: bytes) -> None:
        try:
            batch_validation_result = validate_xml_batch_payload(
                xml_payload_bytes=payload_bytes,
                expected_protein_uids=planned_batch.protein_uids,
                validate_protein_uids=validate_batch_protein_uids,
            )
        except ET.ParseError as error:
            raise NCBITransientFetchError(
                "NCBI returned a payload that is not well-formed XML for "
                f"{batch_context}: {error}"
            ) from error
        except RuntimeError as error:
            raise NCBIXmlResponseValidationError(
                "NCBI returned a partial or mismatched XML batch; "
                f"{error} {batch_context}."
            ) from error

        if batch_validation_result.record_count != planned_batch.protein_uid_count:
            raise NCBIXmlResponseValidationError(
                "NCBI returned a partial or mismatched XML batch; "
                "the direct XML record count does not match the "
                "requested protein UID count. "
                f"Expected {planned_batch.protein_uid_count}, "
                f"got {batch_validation_result.record_count}; {batch_context}."
            )

        validated_record_count["record_count"] = (
            batch_validation_result.record_count
        )

    def announce_validated_batch() -> None:
        if not controls.verbose_logging:
            return

        controls.write_log_line(
            f"Fetched XML batch {planned_batch.batch_index}/{total_batch_count} "
            f"with {planned_batch.protein_uid_count} protein UIDs."
        )

    controlled_request_result = _execute_ncbi_request_with_controls(
        stage_name=XML_BATCH_STAGE_NAME,
        request_context=batch_context,
        operation_label="XML batch",
        open_request=open_batch_request,
        controls=controls,
        permanent_failure_message_prefix=(
            "Permanent NCBI XML batch failure; not retrying"
        ),
        build_retry_exhausted_error=lambda error: NCBIXmlBatchFetchError(
            f"Failed to fetch NCBI XML batch after "
            f"{controls.max_retry_attempts} attempts for {batch_context}: "
            f"{error}"
        ),
        validate_payload=validate_batch_payload,
        circuit_breaker=circuit_breaker,
        rate_limiter=rate_limiter,
        concurrency_governor=concurrency_governor,
        on_attempt_start=announce_attempt,
        on_response_received=announce_response,
        on_response_validated=announce_validated_batch,
    )

    persisted_payload_file_path: Optional[Path] = None
    if workspace is not None:
        persisted_payload_file_path = workspace.store_batch(
            planned_batch=planned_batch,
            xml_payload_bytes=controlled_request_result.payload_bytes,
            xml_payload_sha256=controlled_request_result.payload_sha256,
            record_count=validated_record_count.get(
                "record_count",
                planned_batch.protein_uid_count,
            ),
            round_trip_latency_seconds=(
                controlled_request_result.round_trip_latency_seconds
            ),
            attempt_count=controlled_request_result.attempt_count,
            fetch_timeout_seconds=controls.fetch_timeout_seconds,
            batch_deadline_seconds=controls.request_deadline_seconds,
        )

    if progress_reporter is not None:
        progress_reporter.report_completed_units(
            completed_unit_delta=1,
            item_label=f"#{planned_batch.batch_index}",
            response_byte_count=len(controlled_request_result.payload_bytes),
            round_trip_latency_seconds=(
                controlled_request_result.round_trip_latency_seconds
            ),
        )

    return NCBIProteinXmlBatchFetchResult(
        batch_index=planned_batch.batch_index,
        batch_start_index=planned_batch.batch_start_index,
        batch_end_index=planned_batch.batch_end_index,
        protein_uids=planned_batch.protein_uids,
        protein_uid_count=planned_batch.protein_uid_count,
        xml_payload_sha256=controlled_request_result.payload_sha256,
        xml_payload_bytes=(
            None
            if persisted_payload_file_path is not None
            else controlled_request_result.payload_bytes
        ),
        xml_payload_file_path=(
            str(persisted_payload_file_path)
            if persisted_payload_file_path is not None
            else None
        ),
        record_count=validated_record_count.get("record_count"),
        response_byte_count=len(controlled_request_result.payload_bytes),
        round_trip_latency_seconds=(
            controlled_request_result.round_trip_latency_seconds
        ),
        attempt_count=controlled_request_result.attempt_count,
        reused_from_workspace=False,
    )


def _fetch_planned_xml_batches_sequentially(
    *,
    planned_batches: Sequence[PlannedXmlBatch],
    fetch_one_batch: Callable[[PlannedXmlBatch], NCBIProteinXmlBatchFetchResult],
) -> List[NCBIProteinXmlBatchFetchResult]:
    return [
        fetch_one_batch(planned_batch)
        for planned_batch in planned_batches
    ]


def _fetch_planned_xml_batches_concurrently(
    *,
    planned_batches: Sequence[PlannedXmlBatch],
    fetch_one_batch: Callable[[PlannedXmlBatch], NCBIProteinXmlBatchFetchResult],
    max_concurrent_requests: int,
) -> List[NCBIProteinXmlBatchFetchResult]:
    """
    Fetch planned batches concurrently while preserving plan order.

    Results are written into a slot addressed by the batch index rather than
    appended, because completion order under concurrency is not plan order and
    the consolidated XML hash depends on plan order alone. The worker context is
    copied at submit time so that a custom certificate authority configured in
    the caller's context still applies inside the executor threads.
    """
    ordered_batch_results: List[Optional[NCBIProteinXmlBatchFetchResult]] = [
        None
    ] * len(planned_batches)
    abort_requested = threading.Event()
    first_fetch_error: Optional[BaseException] = None

    def fetch_planned_batch(
        planned_batch: PlannedXmlBatch,
    ) -> NCBIProteinXmlBatchFetchResult:
        if abort_requested.is_set():
            raise NCBIXmlBatchFetchError(
                "NCBI XML batch "
                f"{planned_batch.batch_index} was skipped because an earlier "
                "batch failed."
            )

        return fetch_one_batch(planned_batch)

    with ThreadPoolExecutor(
        max_workers=max_concurrent_requests,
        thread_name_prefix="ncbi-xml-batch",
    ) as batch_executor:
        submitted_futures = {}
        for planned_batch in planned_batches:
            worker_context = copy_context()
            submitted_futures[
                batch_executor.submit(
                    worker_context.run,
                    fetch_planned_batch,
                    planned_batch,
                )
            ] = planned_batch

        for completed_future in as_completed(submitted_futures):
            completed_planned_batch = submitted_futures[completed_future]
            try:
                ordered_batch_results[
                    completed_planned_batch.batch_index - 1
                ] = completed_future.result()
            except BaseException as error:
                abort_requested.set()
                if first_fetch_error is None:
                    first_fetch_error = error

    if first_fetch_error is not None:
        raise first_fetch_error

    return [
        batch_result
        for batch_result in ordered_batch_results
        if batch_result is not None
    ]


def fetch_ncbi_protein_xml_batches(
    *,
    ncbi_email: str,
    ncbi_api_key: Optional[str],
    protein_uids: List[str],
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    batch_size: int = 500,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
    fetch_timeout_seconds: float = 30.0,
    batch_deadline_seconds: float = 300.0,
    retry_backoff_initial_seconds: Optional[float] = None,
    retry_backoff_multiplier: float = 2.0,
    retry_backoff_max_seconds: float = 30.0,
    rate_limit_backoff_seconds: float = 5.0,
    circuit_breaker_failure_threshold: int = 3,
    circuit_breaker_cooldown_seconds: float = 60.0,
    circuit_breaker: Optional[NCBIXmlCircuitBreaker] = None,
    rettype: Optional[str] = NCBI_PROTEIN_GBSEQ_XML_RETTYPE,
    max_concurrent_requests: int = 1,
    max_request_starts_per_second: Optional[float] = None,
    reuse_http_connection: bool = False,
    batch_workspace_directory: Optional[PathLike] = None,
    enable_batch_resume: bool = True,
    spill_batch_payloads_to_disk: bool = True,
    batch_spill_parent_directory: Optional[PathLike] = None,
    validate_batch_protein_uids: bool = True,
    verbose_batch_logging: bool = False,
    telemetry: Optional[NCBIRetrievalTelemetry] = None,
) -> NCBIProteinXmlFetchResult:
    """
    Fetch raw protein XML batches from NCBI Entrez using a frozen list of protein UIDs.

    This function is intended to run downstream of a previously frozen UID snapshot.
    It performs remote retrieval only and does not publish any local snapshot.

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
            Number of protein UIDs requested per EFetch call. EFetch paging is
            capped at 10,000 identifiers per request.
        max_retry_attempts:
            Maximum retry attempts for each failed batch request.
        request_delay_seconds:
            Optional custom delay between successful requests.
            If None, defaults to:
            - 0.10 seconds when an API key is provided
            - 0.34 seconds otherwise
            The delay is not applied when bounded concurrency is enabled,
            because a shared request-start rate limiter governs pacing instead.
        fetch_timeout_seconds:
            Socket/network idle timeout applied to opening and reading each
            NCBI EFetch response.
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
        rate_limit_backoff_seconds:
            Minimum delay applied before retrying a request that NCBI rejected
            for exceeding the request-rate allowance.
        circuit_breaker_failure_threshold:
            Number of consecutive failed batches that opens the circuit breaker.
        circuit_breaker_cooldown_seconds:
            Cooldown interval before another request is allowed after opening.
        circuit_breaker:
            Optional explicitly scoped circuit breaker. When omitted, a
            process-persistent breaker for the configured threshold and
            cooldown is reused automatically across calls.
        rettype:
            EFetch record type. Defaults to the documented GBSeq XML value
            ``gp``; pass None to reproduce the previous undocumented default.
        max_concurrent_requests:
            Number of XML batches fetched in parallel. Defaults to 1 so that
            batch-size effects can be measured without concurrency present.
        max_request_starts_per_second:
            Ceiling on request *start* frequency shared by all workers. If None,
            a conservative default is chosen from API key presence.
        reuse_http_connection:
            Reuse one keep-alive HTTPS connection per worker thread instead of
            opening a new connection, and therefore a new TLS handshake, per
            request.
        batch_workspace_directory:
            Directory holding fetched batch payloads. When provided, a validated
            batch from an earlier run is reused instead of being re-fetched.
        enable_batch_resume:
            Whether a provided workspace may be reused across runs.
        spill_batch_payloads_to_disk:
            Persist each batch payload as soon as it is validated instead of
            retaining every payload in memory until consolidation.
        batch_spill_parent_directory:
            Parent directory for the ephemeral spill directory used when no
            workspace directory is provided.
        validate_batch_protein_uids:
            Validate the exact protein UID set of every batch, not only its
            record count.
        verbose_batch_logging:
            Include full UID lists, request URLs, and response headers in logs.
        telemetry:
            Optional shared telemetry collector. One is created when omitted.

    Partial responses:
        NCBI may omit an invalid, removed, or otherwise unavailable UID while
        returning syntactically valid XML. A batch whose record count or exact
        UID set does not match the request is treated as a permanent response
        validation failure and is not retried, because retrying cannot establish
        which UID was omitted. The failing UIDs are named in the error so the
        offending batch can be isolated and excluded deliberately.

    Returns:
        NCBIProteinXmlFetchResult:
            Retrieval metadata plus, for each batch, either the raw payload
            bytes or the path of the persisted payload file.
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

    if batch_size > NCBI_MAX_EFETCH_RETMAX:
        raise ValueError(
            "batch_size must not exceed the EFetch paging ceiling of "
            f"{NCBI_MAX_EFETCH_RETMAX} identifiers per request."
        )

    _validate_positive_integer(
        parameter_name="max_retry_attempts",
        parameter_value=max_retry_attempts,
    )

    _validate_positive_integer(
        parameter_name="max_concurrent_requests",
        parameter_value=max_concurrent_requests,
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
        rate_limit_backoff_seconds=rate_limit_backoff_seconds,
    )

    active_circuit_breaker = circuit_breaker
    if active_circuit_breaker is None:
        active_circuit_breaker = get_default_ncbi_xml_circuit_breaker(
            failure_threshold=circuit_breaker_failure_threshold,
            cooldown_seconds=circuit_breaker_cooldown_seconds,
        )

    active_telemetry = telemetry if telemetry is not None else NCBIRetrievalTelemetry()
    retrieved_at_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    protein_uids_sha256 = sha256_of_lines(
        text_lines=normalized_protein_uids,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )

    planned_batches = build_xml_batch_plan(
        protein_uids=normalized_protein_uids,
        batch_size=batch_size,
    )
    total_batch_count = len(planned_batches)

    active_workspace: Optional[NCBIXmlBatchWorkspace] = None
    if batch_workspace_directory is not None or spill_batch_payloads_to_disk:
        active_workspace = NCBIXmlBatchWorkspace(
            workspace_directory=batch_workspace_directory,
            plan_identity={
                "protein_uids_sha256": protein_uids_sha256,
                "protein_uid_count": len(normalized_protein_uids),
                "batch_size": batch_size,
                "database_name": database_name,
                "retmode": retmode,
                "rettype": rettype,
            },
            request_policy={
                "max_retry_attempts": max_retry_attempts,
                "fetch_timeout_seconds": fetch_timeout_seconds,
                "batch_deadline_seconds": batch_deadline_seconds,
                "request_delay_seconds": resolved_request_delay_seconds,
                "retry_backoff_initial_seconds": (
                    resolved_retry_backoff_initial_seconds
                ),
                "retry_backoff_multiplier": retry_backoff_multiplier,
                "retry_backoff_max_seconds": retry_backoff_max_seconds,
                "max_concurrent_requests": max_concurrent_requests,
            },
            enable_resume=enable_batch_resume,
            ephemeral_parent_directory=batch_spill_parent_directory,
        )
        active_workspace.open(planned_batches=planned_batches)

    active_rate_limiter: Optional[NCBIRequestStartRateLimiter] = None
    active_concurrency_governor: Optional[NCBIAdaptiveConcurrencyGovernor] = None
    if max_concurrent_requests > 1:
        active_rate_limiter = NCBIRequestStartRateLimiter(
            max_request_starts_per_second=_resolve_max_request_starts_per_second(
                max_request_starts_per_second=max_request_starts_per_second,
                ncbi_api_key=ncbi_api_key,
            ),
        )
        active_concurrency_governor = NCBIAdaptiveConcurrencyGovernor(
            max_concurrent_requests=max_concurrent_requests,
            rate_limiter=active_rate_limiter,
        )

    batch_progress_reporter = NCBIRetrievalProgressReporter(
        total_unit_count=total_batch_count,
        unit_label="batches",
    )

    try:
        with _configured_ncbi_entrez_urlopen(
            ssl_ca_file=ssl_ca_file,
            ssl_ca_directory=ssl_ca_directory,
        ) as (resolved_ssl_ca_file, resolved_ssl_ca_directory):
            request_controls = _NCBIRequestControls(
                max_retry_attempts=max_retry_attempts,
                fetch_timeout_seconds=fetch_timeout_seconds,
                request_deadline_seconds=batch_deadline_seconds,
                request_delay_seconds=resolved_request_delay_seconds,
                retry_backoff_initial_seconds=(
                    resolved_retry_backoff_initial_seconds
                ),
                retry_backoff_multiplier=retry_backoff_multiplier,
                retry_backoff_max_seconds=retry_backoff_max_seconds,
                rate_limit_backoff_seconds=rate_limit_backoff_seconds,
                resolved_ssl_ca_file=resolved_ssl_ca_file,
                resolved_ssl_ca_directory=resolved_ssl_ca_directory,
                telemetry=active_telemetry,
                verbose_logging=verbose_batch_logging,
                progress_reporter=batch_progress_reporter,
            )

            def fetch_one_planned_batch(
                planned_batch: PlannedXmlBatch,
            ) -> NCBIProteinXmlBatchFetchResult:
                return _fetch_single_ncbi_xml_batch(
                    planned_batch=planned_batch,
                    total_batch_count=total_batch_count,
                    database_name=database_name,
                    retmode=retmode,
                    rettype=rettype,
                    controls=request_controls,
                    circuit_breaker=active_circuit_breaker,
                    workspace=active_workspace,
                    rate_limiter=active_rate_limiter,
                    concurrency_governor=active_concurrency_governor,
                    validate_batch_protein_uids=validate_batch_protein_uids,
                    reuse_http_connection=reuse_http_connection,
                    progress_reporter=batch_progress_reporter,
                )

            with active_telemetry.measure_stage(stage_name=XML_BATCH_STAGE_NAME):
                if max_concurrent_requests == 1:
                    xml_batches = _fetch_planned_xml_batches_sequentially(
                        planned_batches=planned_batches,
                        fetch_one_batch=fetch_one_planned_batch,
                    )
                else:
                    xml_batches = _fetch_planned_xml_batches_concurrently(
                        planned_batches=planned_batches,
                        fetch_one_batch=fetch_one_planned_batch,
                        max_concurrent_requests=max_concurrent_requests,
                    )
    except BaseException:
        if active_workspace is not None and active_workspace.is_ephemeral:
            active_workspace.purge()
        raise
    finally:
        batch_progress_reporter.finish()
        if reuse_http_connection:
            _ncbi_persistent_https_transport.discard_all_connections()

    reused_batch_count = sum(
        1 for xml_batch in xml_batches if xml_batch.reused_from_workspace
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
        fetch_timeout_seconds=fetch_timeout_seconds,
        batch_deadline_seconds=batch_deadline_seconds,
        retry_backoff_initial_seconds=resolved_retry_backoff_initial_seconds,
        retry_backoff_multiplier=retry_backoff_multiplier,
        retry_backoff_max_seconds=retry_backoff_max_seconds,
        circuit_breaker_failure_threshold=circuit_breaker_failure_threshold,
        circuit_breaker_cooldown_seconds=circuit_breaker_cooldown_seconds,
        rettype=rettype,
        rate_limit_backoff_seconds=rate_limit_backoff_seconds,
        max_concurrent_requests=max_concurrent_requests,
        max_request_starts_per_second=(
            active_rate_limiter.max_request_starts_per_second
            if active_rate_limiter is not None
            else None
        ),
        reuse_http_connection=reuse_http_connection,
        batch_workspace_directory=(
            str(active_workspace.workspace_directory)
            if active_workspace is not None
            else None
        ),
        batch_workspace_is_ephemeral=(
            active_workspace.is_ephemeral if active_workspace is not None else False
        ),
        reused_batch_count=reused_batch_count,
        fetched_batch_count=len(xml_batches) - reused_batch_count,
        retrieval_telemetry=active_telemetry.build_summary(),
    )

    total_retrieval_seconds = batch_progress_reporter.elapsed_seconds
    print(
        f"Fetched {fetch_result.batch_count} XML batches for "
        f"{fetch_result.normalized_protein_uid_count} protein UIDs "
        f"in {_format_duration_seconds(total_retrieval_seconds)}."
    )
    if total_retrieval_seconds > 0 and batch_progress_reporter.network_byte_count:
        print(
            "Transferred "
            f"{batch_progress_reporter.network_byte_count / (1024 * 1024):.1f} MB "
            "at "
            f"{batch_progress_reporter.network_byte_count / 1024 / total_retrieval_seconds:,.0f}"
            " KB/s."
        )
    if reused_batch_count:
        print(
            f"Reused {reused_batch_count} previously validated XML batches "
            "from the batch workspace."
        )

    return fetch_result
