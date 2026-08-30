from __future__ import annotations

import io
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Optional

from Bio import Entrez

from src.pago_pipeline.ncbi_api import _configured_ncbi_entrez_urlopen
from src.pago_pipeline.ncbi_xml_stream import extract_protein_uid_from_gbseq_element

DEFAULT_MAX_UID_COUNT = 150_000
DEFAULT_SAMPLE_SIZE = 200
DEFAULT_MAX_RETRY_ATTEMPTS = 5
DEFAULT_RETRY_BACKOFF_SECONDS = 5.0
NCBI_PROTEIN_DATABASE = "protein"


@dataclass(frozen=True)
class EsearchPreflightResult:
    """
    Cheap-query summary of what a full NCBI protein retrieval would download.
    """

    search_query: str
    translated_query: Optional[str]
    result_count: int
    history_web_env: Optional[str]
    history_query_key: Optional[str]
    retrieved_at_utc: str
    max_uid_count: int
    exceeds_max_uid_count: bool
    sample_requested_count: int
    sample_uid_list: list[str] = field(default_factory=list)
    sample_record_count: int = 0
    sample_records_with_sequence: int = 0
    sample_records_missing_sequence: int = 0
    sample_records_with_extractable_uid: int = 0
    sample_fetch_error: Optional[str] = None
    python_version: str = ""
    biopython_version: str = ""

    def as_report_payload(self) -> dict[str, object]:
        return asdict(self)


def parse_esearch_history_response(*, response_payload_bytes: bytes) -> dict[str, object]:
    """
    Extract Count, QueryTranslation and History handles from an ESearch response.
    """
    parsed = Entrez.read(io.BytesIO(response_payload_bytes))
    return {
        "result_count": int(parsed["Count"]),
        "translated_query": parsed.get("QueryTranslation"),
        "history_web_env": parsed.get("WebEnv"),
        "history_query_key": parsed.get("QueryKey"),
    }


def parse_uilist_text(*, payload_text: str) -> list[str]:
    return [line.strip() for line in payload_text.splitlines() if line.strip()]


def summarize_sample_xml(*, xml_payload_bytes: bytes) -> dict[str, int]:
    """
    Count records, sequence presence and extractable protein_uid in a sample
    GBSet XML payload.
    """
    root_element = ET.fromstring(xml_payload_bytes)
    gbseq_elements = list(root_element.iter("GBSeq"))
    with_sequence = 0
    with_extractable_uid = 0
    for gbseq_element in gbseq_elements:
        sequence_element = gbseq_element.find("GBSeq_sequence")
        if sequence_element is not None and (sequence_element.text or "").strip():
            with_sequence += 1
        try:
            extract_protein_uid_from_gbseq_element(gbseq_element=gbseq_element)
            with_extractable_uid += 1
        except Exception:  # noqa: BLE001 - a sample QC counter, never fatal here
            pass
    return {
        "sample_record_count": len(gbseq_elements),
        "sample_records_with_sequence": with_sequence,
        "sample_records_missing_sequence": len(gbseq_elements) - with_sequence,
        "sample_records_with_extractable_uid": with_extractable_uid,
    }


def build_esearch_preflight_result(
    *,
    search_query: str,
    retrieved_at_utc: str,
    result_count: int,
    translated_query: Optional[str],
    history_web_env: Optional[str],
    history_query_key: Optional[str],
    max_uid_count: int,
    sample_requested_count: int,
    sample_uid_list: Optional[list[str]] = None,
    sample_summary: Optional[dict[str, int]] = None,
    sample_fetch_error: Optional[str] = None,
) -> EsearchPreflightResult:
    resolved_sample_summary = sample_summary or {}
    return EsearchPreflightResult(
        search_query=search_query,
        translated_query=translated_query,
        result_count=int(result_count),
        history_web_env=history_web_env,
        history_query_key=history_query_key,
        retrieved_at_utc=retrieved_at_utc,
        max_uid_count=int(max_uid_count),
        exceeds_max_uid_count=int(result_count) > int(max_uid_count),
        sample_requested_count=int(sample_requested_count),
        sample_uid_list=list(sample_uid_list or []),
        sample_record_count=int(resolved_sample_summary.get("sample_record_count", 0)),
        sample_records_with_sequence=int(
            resolved_sample_summary.get("sample_records_with_sequence", 0)
        ),
        sample_records_missing_sequence=int(
            resolved_sample_summary.get("sample_records_missing_sequence", 0)
        ),
        sample_records_with_extractable_uid=int(
            resolved_sample_summary.get("sample_records_with_extractable_uid", 0)
        ),
        sample_fetch_error=sample_fetch_error,
        python_version=sys.version.split()[0],
        biopython_version=getattr(
            __import__("Bio"), "__version__", "unknown"
        ),
    )


def _read_entrez_handle_bytes(handle) -> bytes:
    try:
        payload = handle.read()
    finally:
        handle.close()
    if isinstance(payload, str):
        return payload.encode("utf-8")
    return payload


def _run_with_retries(
    *,
    operation_label: str,
    open_request,
    max_retry_attempts: int,
    retry_backoff_seconds: float,
) -> bytes:
    last_error: Optional[Exception] = None
    for attempt_index in range(max_retry_attempts):
        try:
            return _read_entrez_handle_bytes(open_request())
        except Exception as error:  # noqa: BLE001 - retried transport failures
            last_error = error
            if attempt_index + 1 < max_retry_attempts:
                time.sleep(retry_backoff_seconds * (2 ** attempt_index))
    raise RuntimeError(
        f"NCBI {operation_label} failed after {max_retry_attempts} attempts: "
        f"{last_error}"
    )


def run_ncbi_esearch_preflight(
    *,
    search_query: str,
    ncbi_email: str,
    ncbi_api_key: Optional[str] = None,
    max_uid_count: int = DEFAULT_MAX_UID_COUNT,
    sample_size: int = DEFAULT_SAMPLE_SIZE,
    max_retry_attempts: int = DEFAULT_MAX_RETRY_ATTEMPTS,
    retry_backoff_seconds: float = DEFAULT_RETRY_BACKOFF_SECONDS,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
) -> EsearchPreflightResult:
    """
    Run a single ESearch (with History) plus one small sample fetch to describe
    what a full retrieval for search_query would cost, before committing to it.
    """
    if not search_query or not search_query.strip():
        raise ValueError("search_query must be a non-empty string.")
    if not ncbi_email:
        raise ValueError("ncbi_email is required for an NCBI ESearch preflight.")
    if max_uid_count <= 0:
        raise ValueError("max_uid_count must be a positive integer.")
    if sample_size < 0:
        raise ValueError("sample_size must be a non-negative integer.")

    Entrez.email = ncbi_email
    Entrez.api_key = ncbi_api_key

    retrieved_at_utc = (
        datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
            "+00:00", "Z"
        )
    )

    with _configured_ncbi_entrez_urlopen(
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
    ):
        esearch_bytes = _run_with_retries(
            operation_label="ESearch preflight",
            open_request=lambda: Entrez.esearch(
                db=NCBI_PROTEIN_DATABASE,
                term=search_query,
                retmax=0,
                usehistory="y",
            ),
            max_retry_attempts=max_retry_attempts,
            retry_backoff_seconds=retry_backoff_seconds,
        )
        history = parse_esearch_history_response(
            response_payload_bytes=esearch_bytes
        )
        result_count = int(history["result_count"])
        history_web_env = history["history_web_env"]
        history_query_key = history["history_query_key"]

        sample_requested_count = min(sample_size, result_count) if result_count else 0
        sample_uid_list: list[str] = []
        sample_summary: Optional[dict[str, int]] = None
        sample_fetch_error: Optional[str] = None

        if (
            sample_requested_count > 0
            and history_web_env
            and history_query_key
        ):
            try:
                uilist_bytes = _run_with_retries(
                    operation_label="ESearch preflight sample UID fetch",
                    open_request=lambda: Entrez.efetch(
                        db=NCBI_PROTEIN_DATABASE,
                        rettype="uilist",
                        retmode="text",
                        retstart=0,
                        retmax=sample_requested_count,
                        webenv=history_web_env,
                        query_key=history_query_key,
                    ),
                    max_retry_attempts=max_retry_attempts,
                    retry_backoff_seconds=retry_backoff_seconds,
                )
                sample_uid_list = parse_uilist_text(
                    payload_text=uilist_bytes.decode("utf-8")
                )
                if sample_uid_list:
                    sample_xml_bytes = _run_with_retries(
                        operation_label="ESearch preflight sample XML fetch",
                        open_request=lambda: Entrez.efetch(
                            db=NCBI_PROTEIN_DATABASE,
                            id=",".join(sample_uid_list),
                            rettype="gb",
                            retmode="xml",
                        ),
                        max_retry_attempts=max_retry_attempts,
                        retry_backoff_seconds=retry_backoff_seconds,
                    )
                    sample_summary = summarize_sample_xml(
                        xml_payload_bytes=sample_xml_bytes
                    )
            except Exception as error:  # noqa: BLE001 - sample QC is best-effort
                sample_fetch_error = str(error)

    return build_esearch_preflight_result(
        search_query=search_query,
        retrieved_at_utc=retrieved_at_utc,
        result_count=result_count,
        translated_query=history["translated_query"],
        history_web_env=history_web_env,
        history_query_key=history_query_key,
        max_uid_count=max_uid_count,
        sample_requested_count=sample_requested_count,
        sample_uid_list=sample_uid_list,
        sample_summary=sample_summary,
        sample_fetch_error=sample_fetch_error,
    )
