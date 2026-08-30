from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass

import pandas as pd

DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN = "protein_uid"
DEFAULT_RETRIEVED_ACCESSION_COLUMN = "gbseq__accession_version"

REQUIRED_REFERENCE_COLUMNS: tuple[str, ...] = (
    "accession",
    "ago_family",
    "clade",
    "reference_label_source",
    "reference_label_evidence",
)

# Recall is reported overall and per MID-PIWI clade. PIWI-RE is its own family,
# not a clade, but a stratum is still useful for auditing text-query coverage.
RECALL_STRATA: tuple[tuple[str, str], ...] = (
    ("overall", "overall_reference_recall"),
    ("LONG_A", "long_a_reference_recall"),
    ("LONG_B", "long_b_reference_recall"),
    ("SHORT", "short_reference_recall"),
    ("PIWI_RE", "piwi_re_reference_recall"),
)

_VERSION_SUFFIX_PATTERN = re.compile(r"\.\d+$")


@dataclass(frozen=True)
class QueryReferenceRecallResult:
    summary_dataframe: pd.DataFrame
    detail_dataframe: pd.DataFrame
    stratum_recall: dict[str, float]
    reference_count: int
    recovered_count: int


def _normalize_accession(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def _bare_accession(accession: str) -> str:
    return _VERSION_SUFFIX_PATTERN.sub("", accession)


def _validate_reference_columns(*, reference_dataframe: pd.DataFrame) -> None:
    missing_columns = [
        column_name
        for column_name in REQUIRED_REFERENCE_COLUMNS
        if column_name not in reference_dataframe.columns
    ]
    if missing_columns:
        raise RuntimeError(
            "reference_dataframe is missing required columns: "
            f"{missing_columns}."
        )


def compute_query_reference_recall(
    *,
    reference_dataframe: pd.DataFrame,
    retrieved_metadata_dataframe: pd.DataFrame,
    retrieved_accession_column: str = DEFAULT_RETRIEVED_ACCESSION_COLUMN,
    retrieved_protein_uid_column: str = DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN,
    progress_callback: Callable[[str], None] | None = None,
) -> QueryReferenceRecallResult:
    """
    Measure how many curated reference pAgos the retrieved candidate set
    recovered, overall and stratified by MID-PIWI clade.

    Matching is done on accession.version and its version-stripped form, so a
    reference listed as ``WP_000000001.1`` is recovered by a retrieved
    ``WP_000000001.2`` and vice versa.
    """
    _validate_reference_columns(reference_dataframe=reference_dataframe)

    if retrieved_accession_column not in retrieved_metadata_dataframe.columns:
        raise RuntimeError(
            "retrieved_metadata_dataframe is missing the accession column "
            f"{retrieved_accession_column!r}."
        )

    retrieved_accessions = (
        retrieved_metadata_dataframe[retrieved_accession_column]
        .map(_normalize_accession)
    )
    retrieved_versioned_set = {value for value in retrieved_accessions if value}
    retrieved_bare_set = {
        _bare_accession(value) for value in retrieved_versioned_set
    }

    retrieved_protein_uid_by_bare: dict[str, str] = {}
    if retrieved_protein_uid_column in retrieved_metadata_dataframe.columns:
        for accession_value, uid_value in zip(
            retrieved_accessions,
            retrieved_metadata_dataframe[retrieved_protein_uid_column].astype(str),
        ):
            if accession_value:
                retrieved_protein_uid_by_bare.setdefault(
                    _bare_accession(accession_value), uid_value
                )

    detail_rows: list[dict[str, object]] = []
    for _, reference_row in reference_dataframe.iterrows():
        reference_accession = _normalize_accession(reference_row["accession"])
        reference_bare = _bare_accession(reference_accession)
        matched_versioned: str | None = None
        recovered = False
        if reference_accession and reference_accession in retrieved_versioned_set:
            recovered = True
            matched_versioned = reference_accession
        elif reference_bare and reference_bare in retrieved_bare_set:
            recovered = True
            matched_versioned = reference_bare

        detail_rows.append(
            {
                "accession": reference_accession,
                "clade": str(reference_row["clade"]).strip().upper(),
                "ago_family": str(reference_row["ago_family"]).strip().upper(),
                "recovered": bool(recovered),
                "matched_accession": matched_versioned or "",
                "matched_protein_uid": (
                    retrieved_protein_uid_by_bare.get(reference_bare, "")
                    if recovered
                    else ""
                ),
                "reference_label_source": reference_row["reference_label_source"],
                "reference_label_evidence": reference_row[
                    "reference_label_evidence"
                ],
            }
        )

    detail_dataframe = pd.DataFrame(detail_rows)

    summary_rows: list[dict[str, object]] = []
    stratum_recall: dict[str, float] = {}
    for stratum_name, metric_name in RECALL_STRATA:
        if stratum_name == "overall":
            stratum_frame = detail_dataframe
        else:
            stratum_frame = detail_dataframe[
                detail_dataframe["clade"] == stratum_name
            ]
        reference_count = int(len(stratum_frame))
        recovered_count = int(stratum_frame["recovered"].sum())
        recall_value = (
            float(recovered_count / reference_count) if reference_count else 0.0
        )
        stratum_recall[metric_name] = recall_value
        summary_rows.append(
            {
                "stratum": stratum_name,
                "metric_name": metric_name,
                "reference_count": reference_count,
                "recovered_count": recovered_count,
                "recall": recall_value,
            }
        )

    summary_dataframe = pd.DataFrame(summary_rows)

    if progress_callback is not None:
        progress_callback(
            "query reference recall: "
            f"{int(detail_dataframe['recovered'].sum())}/"
            f"{len(detail_dataframe)} references recovered"
        )

    return QueryReferenceRecallResult(
        summary_dataframe=summary_dataframe,
        detail_dataframe=detail_dataframe,
        stratum_recall=stratum_recall,
        reference_count=int(len(detail_dataframe)),
        recovered_count=int(detail_dataframe["recovered"].sum()),
    )
