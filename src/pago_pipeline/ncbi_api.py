from __future__ import annotations

import platform
import random
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import List, Optional

import Bio
from Bio import Entrez

from src.pago_pipeline.storage import sha256_of_lines


@dataclass(frozen=True)
class NCBIProteinIdFetchResult:
    """
    Immutable payload representing one NCBI protein-ID retrieval event.

    This object is designed to become the basis of a reproducible local snapshot.
    """
    database_name: str
    search_query: str
    translated_query: Optional[str]
    identifier_type: str
    retrieved_at_utc: str
    ncbi_reported_result_count: int
    protein_ids: List[str]
    raw_protein_id_count: int
    normalized_protein_id_count: int
    deduplicate_ids: bool
    sort_ids: bool
    protein_ids_sha256: str
    page_size: int
    request_delay_seconds: float
    max_retry_attempts: int
    history_web_env: Optional[str]
    history_query_key: Optional[str]
    python_version: str
    biopython_version: str


def _normalize_protein_id_list(
    *,
    protein_ids: List[str],
    deduplicate_ids: bool,
    sort_ids: bool,
) -> List[str]:
    """
    Normalize protein IDs by:
    - casting to str
    - stripping surrounding whitespace
    - removing empty values
    - optionally deduplicating while preserving first occurrence order
    - optionally sorting canonically
    """
    normalized_protein_ids = [
        str(protein_id).strip()
        for protein_id in protein_ids
        if str(protein_id).strip()
    ]

    if deduplicate_ids:
        unique_protein_ids: List[str] = []
        seen_protein_ids = set()

        for protein_id in normalized_protein_ids:
            if protein_id not in seen_protein_ids:
                seen_protein_ids.add(protein_id)
                unique_protein_ids.append(protein_id)

        normalized_protein_ids = unique_protein_ids

    if sort_ids:
        normalized_protein_ids = sorted(normalized_protein_ids)

    return normalized_protein_ids


def fetch_ncbi_protein_id_snapshot(
    *,
    ncbi_email: str,
    ncbi_api_key: Optional[str],
    query: str,
    deduplicate_ids: bool = True,
    sort_ids: bool = True,
    page_size: int = 1000,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
) -> NCBIProteinIdFetchResult:
    """
    Retrieve protein IDs from NCBI Entrez and return a snapshot-ready payload.

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
        deduplicate_ids:
            Whether to deduplicate IDs in the returned payload.
        sort_ids:
            Whether to sort IDs canonically in the returned payload.
        page_size:
            Number of protein IDs requested per paginated Entrez call.
        max_retry_attempts:
            Maximum retry attempts for each failed paginated request.
        request_delay_seconds:
            Optional custom delay between successful requests.
            If None, defaults to:
            - 0.10 seconds when an API key is provided
            - 0.34 seconds otherwise

    Returns:
        NCBIProteinIdFetchResult:
            Snapshot-ready result including IDs and retrieval metadata.
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
    identifier_type = "acc"

    if request_delay_seconds is None:
        request_delay_seconds = 0.10 if ncbi_api_key else 0.34

    retrieved_at_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    initial_search_handle = Entrez.esearch(
        db=database_name,
        term=query,
        retmax=0,
        idtype=identifier_type,
        usehistory="y",
    )
    try:
        initial_search_response = Entrez.read(initial_search_handle)
    finally:
        initial_search_handle.close()

    ncbi_reported_result_count = int(initial_search_response["Count"])
    translated_query = initial_search_response.get("QueryTranslation")
    history_web_env = initial_search_response.get("WebEnv")
    history_query_key = initial_search_response.get("QueryKey")

    print(f"Found {ncbi_reported_result_count} protein IDs.")
    print(
        f"History session: WebEnv={history_web_env}, "
        f"QueryKey={history_query_key}."
    )

    if ncbi_reported_result_count == 0:
        print("NCBI returned zero protein IDs for the provided query.")

    raw_protein_ids: List[str] = []
    page_start_index = 0

    while page_start_index < ncbi_reported_result_count:
        for retry_attempt_index in range(max_retry_attempts):
            try:
                paginated_search_handle = Entrez.esearch(
                    db=database_name,
                    term=query,
                    idtype=identifier_type,
                    usehistory="y",
                    retmax=page_size,
                    retstart=page_start_index,
                    WebEnv=history_web_env,
                    query_key=history_query_key,
                )
                try:
                    search_response = Entrez.read(paginated_search_handle)
                finally:
                    paginated_search_handle.close()

                if "IdList" not in search_response:
                    raise RuntimeError("NCBI response missing 'IdList'.")

                current_batch_protein_ids = search_response["IdList"]

                if not current_batch_protein_ids:
                    raise RuntimeError(
                        "NCBI returned an empty IdList before the reported "
                        "result set was fully paginated."
                    )

                raw_protein_ids.extend(current_batch_protein_ids)
                page_start_index += len(current_batch_protein_ids)

                retrieval_progress = min(
                    page_start_index / ncbi_reported_result_count,
                    1.0,
                )
                print(
                    f"Extracted {len(raw_protein_ids)} protein IDs "
                    f"({retrieval_progress:.2%})."
                )

                jitter_seconds = random.uniform(0.01, 0.05)
                time.sleep(request_delay_seconds + jitter_seconds)
                break

            except Exception as error:
                if retry_attempt_index == max_retry_attempts - 1:
                    raise RuntimeError(
                        f"Failed to extract protein IDs after "
                        f"{max_retry_attempts} attempts at "
                        f"page_start_index={page_start_index}: {error}"
                    ) from error

                retry_backoff_seconds = (
                    2 ** retry_attempt_index
                ) * request_delay_seconds

                print(
                    f"Retry {retry_attempt_index + 1}/{max_retry_attempts} "
                    f"after error: {error}"
                )
                time.sleep(retry_backoff_seconds)

    normalized_protein_ids = _normalize_protein_id_list(
        protein_ids=raw_protein_ids,
        deduplicate_ids=deduplicate_ids,
        sort_ids=sort_ids,
    )

    # The protein ID list has already been normalized according to the
    # selected deduplication/sorting policy, so hashing must operate on
    # the list exactly as received here.
    protein_ids_sha256 = sha256_of_lines(
        text_lines=normalized_protein_ids,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )

    fetch_result = NCBIProteinIdFetchResult(
        database_name=database_name,
        search_query=query,
        translated_query=translated_query,
        identifier_type=identifier_type,
        retrieved_at_utc=retrieved_at_utc,
        ncbi_reported_result_count=ncbi_reported_result_count,
        protein_ids=normalized_protein_ids,
        raw_protein_id_count=len(raw_protein_ids),
        normalized_protein_id_count=len(normalized_protein_ids),
        deduplicate_ids=deduplicate_ids,
        sort_ids=sort_ids,
        protein_ids_sha256=protein_ids_sha256,
        page_size=page_size,
        request_delay_seconds=request_delay_seconds,
        max_retry_attempts=max_retry_attempts,
        history_web_env=history_web_env,
        history_query_key=history_query_key,
        python_version=platform.python_version(),
        biopython_version=getattr(Bio, "__version__", "unknown"),
    )

    print(f"Final raw ID count: {fetch_result.raw_protein_id_count}")
    print(
        "Final normalized ID count: "
        f"{fetch_result.normalized_protein_id_count}"
    )
    print(f"Dataset SHA-256: {fetch_result.protein_ids_sha256}")

    return fetch_result