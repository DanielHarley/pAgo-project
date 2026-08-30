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
# not a pAgo clade, so its stratum is selected on ago_family rather than clade
# (a PIWI-RE reference row leaves clade empty / NA). Each stratum is
# (display_name, metric_name, detail_column, detail_value); a None column means
# "all rows".
RECALL_STRATA: tuple[tuple[str, str, str | None, str | None], ...] = (
    ("overall", "overall_reference_recall", None, None),
    ("LONG_A", "long_a_reference_recall", "clade", "LONG_A"),
    ("LONG_B", "long_b_reference_recall", "clade", "LONG_B"),
    ("SHORT", "short_reference_recall", "clade", "SHORT"),
    ("PIWI_RE", "piwi_re_reference_recall", "ago_family", "PIWI_RE"),
)

_VERSION_SUFFIX_PATTERN = re.compile(r"\.\d+$")


RECALL_STATUS_EVALUABLE = "EVALUABLE"
RECALL_STATUS_NOT_EVALUABLE = "NOT_EVALUABLE"


@dataclass(frozen=True)
class QueryReferenceRecallResult:
    summary_dataframe: pd.DataFrame
    detail_dataframe: pd.DataFrame
    # None (not 0.0) for a stratum with zero reference sequences.
    stratum_recall: dict[str, float | None]
    stratum_recall_status: dict[str, str]
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
    stratum_recall: dict[str, float | None] = {}
    stratum_recall_status: dict[str, str] = {}
    for stratum_name, metric_name, detail_column, detail_value in RECALL_STRATA:
        if detail_column is None:
            stratum_frame = detail_dataframe
        else:
            stratum_frame = detail_dataframe[
                detail_dataframe[detail_column] == detail_value
            ]
        reference_count = int(len(stratum_frame))
        recovered_count = int(stratum_frame["recovered"].sum())
        if reference_count == 0:
            recall_value: float | None = None
            recall_status = RECALL_STATUS_NOT_EVALUABLE
        else:
            recall_value = float(recovered_count / reference_count)
            recall_status = RECALL_STATUS_EVALUABLE
        stratum_recall[metric_name] = recall_value
        stratum_recall_status[metric_name] = recall_status
        summary_rows.append(
            {
                "stratum": stratum_name,
                "metric_name": metric_name,
                "reference_count": reference_count,
                "recovered_count": recovered_count,
                "recall": (
                    float("nan") if recall_value is None else recall_value
                ),
                "recall_status": recall_status,
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
        stratum_recall_status=stratum_recall_status,
        reference_count=int(len(detail_dataframe)),
        recovered_count=int(detail_dataframe["recovered"].sum()),
    )
