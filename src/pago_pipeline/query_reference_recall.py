from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

import pandas as pd

DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN = "protein_uid"
DEFAULT_RETRIEVED_ACCESSION_COLUMN = "gbseq__accession_version"
DEFAULT_RETRIEVED_SEQUENCE_COLUMN = "gbseq__sequence"

REQUIRED_REFERENCE_COLUMNS: tuple[str, ...] = (
    "accession",
    "ago_family",
    "clade",
    "reference_label_source",
    "reference_label_evidence",
)
# Optional but strongly recommended: enables the SEQUENCE_SHA256 match method.
REFERENCE_SEQUENCE_SHA256_COLUMN = "sequence_sha256"

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

# Documented, irreversible-free protein-sequence normalization applied before
# hashing on BOTH the reference side (during curation) and the retrieved side.
SEQUENCE_NORMALIZATION = "strip_all_whitespace_then_uppercase"


class ReferenceMatchMethod(str, Enum):
    """
    How a curated reference was matched to the retrieved candidate set,
    ordered by decreasing identifier strictness.
    """

    EXACT_ACCESSION_VERSION = "EXACT_ACCESSION_VERSION"
    SAME_BASE_ACCESSION = "SAME_BASE_ACCESSION"
    SEQUENCE_SHA256 = "SEQUENCE_SHA256"
    NONE = "NONE"


# Identifier-strict methods: the retrieval "found the same accession".
_STRICT_ACCESSION_METHODS: frozenset[str] = frozenset(
    {
        ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value,
        ReferenceMatchMethod.SAME_BASE_ACCESSION.value,
    }
)


@dataclass(frozen=True)
class QueryReferenceRecallResult:
    summary_dataframe: pd.DataFrame
    detail_dataframe: pd.DataFrame
    # Strict-identifier reading (EXACT + SAME_BASE only). None (not 0.0) for a
    # stratum with zero reference sequences.
    stratum_exact_recall: dict[str, float | None]
    # Retrieval-equivalent reading: also counts a reference whose exact protein
    # sequence was retrieved under a different accession (SEQUENCE_SHA256).
    stratum_equivalent_recall: dict[str, float | None]
    stratum_recall_status: dict[str, str]
    reference_count: int
    exact_recovered_count: int
    equivalent_recovered_count: int


def normalize_protein_sequence(value: object) -> str:
    """
    Strip all whitespace and uppercase. No other transformation is applied.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return "".join(str(value).split()).upper()


def protein_sequence_sha256(value: object) -> str:
    normalized = normalize_protein_sequence(value)
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def build_matching_strategy_payload() -> dict[str, object]:
    return {
        "strategy_kind": "query_reference_recall_matching",
        "strategy_version": "2.0",
        "match_hierarchy": [
            ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value,
            ReferenceMatchMethod.SAME_BASE_ACCESSION.value,
            ReferenceMatchMethod.SEQUENCE_SHA256.value,
        ],
        "sequence_normalization": SEQUENCE_NORMALIZATION,
        "sequence_hash": "sha256",
        "multi_sequence_hit_representative": "min_accession_lexicographic",
        "recall_readings": [
            "exact_accession_recall",
            "retrieval_equivalent_recall",
        ],
    }


def build_matching_strategy_sha256() -> str:
    serialized = json.dumps(
        build_matching_strategy_payload(),
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


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


def _recall(reference_count: int, recovered_count: int) -> tuple[float | None, str]:
    if reference_count == 0:
        return None, RECALL_STATUS_NOT_EVALUABLE
    return float(recovered_count / reference_count), RECALL_STATUS_EVALUABLE


def compute_query_reference_recall(
    *,
    reference_dataframe: pd.DataFrame,
    retrieved_metadata_dataframe: pd.DataFrame,
    retrieved_accession_column: str = DEFAULT_RETRIEVED_ACCESSION_COLUMN,
    retrieved_protein_uid_column: str = DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN,
    retrieved_sequence_column: str = DEFAULT_RETRIEVED_SEQUENCE_COLUMN,
    progress_callback: Callable[[str], None] | None = None,
) -> QueryReferenceRecallResult:
    """
    Measure how many curated reference pAgos the retrieved candidate set
    recovered, overall and stratified by clade / ago_family.

    A reference is matched by, in order:
      1. EXACT_ACCESSION_VERSION - identical accession.version;
      2. SAME_BASE_ACCESSION    - same base accession, different version;
      3. SEQUENCE_SHA256        - a different accession whose normalized protein
                                  sequence hashes identically (requires the
                                  reference row to carry sequence_sha256).
    Two recall readings are reported per stratum: exact_accession_recall
    (methods 1-2 only) and retrieval_equivalent_recall (methods 1-3).
    """
    _validate_reference_columns(reference_dataframe=reference_dataframe)

    if retrieved_accession_column not in retrieved_metadata_dataframe.columns:
        raise RuntimeError(
            "retrieved_metadata_dataframe is missing the accession column "
            f"{retrieved_accession_column!r}."
        )

    retrieved_accessions = retrieved_metadata_dataframe[
        retrieved_accession_column
    ].map(_normalize_accession)
    if retrieved_protein_uid_column in retrieved_metadata_dataframe.columns:
        retrieved_protein_uids = retrieved_metadata_dataframe[
            retrieved_protein_uid_column
        ].astype("string").fillna("")
    else:
        retrieved_protein_uids = pd.Series(
            [""] * len(retrieved_metadata_dataframe),
            index=retrieved_metadata_dataframe.index,
        )

    retrieved_versioned_set: set[str] = set()
    # bare accession -> deterministic (min) versioned accession + its protein_uid
    retrieved_by_bare: dict[str, tuple[str, str]] = {}
    versioned_to_uid: dict[str, str] = {}
    for accession_value, uid_value in zip(retrieved_accessions, retrieved_protein_uids):
        if not accession_value:
            continue
        retrieved_versioned_set.add(accession_value)
        versioned_to_uid.setdefault(accession_value, str(uid_value))
        bare = _bare_accession(accession_value)
        existing = retrieved_by_bare.get(bare)
        if existing is None or accession_value < existing[0]:
            retrieved_by_bare[bare] = (accession_value, str(uid_value))

    # sequence sha256 -> sorted list of (accession, protein_uid)
    sequence_sha256_to_hits: dict[str, list[tuple[str, str]]] = {}
    has_sequence_column = (
        retrieved_sequence_column in retrieved_metadata_dataframe.columns
    )
    if has_sequence_column:
        for accession_value, uid_value, sequence_value in zip(
            retrieved_accessions,
            retrieved_protein_uids,
            retrieved_metadata_dataframe[retrieved_sequence_column],
        ):
            normalized_sequence = normalize_protein_sequence(sequence_value)
            if not normalized_sequence:
                continue
            digest = hashlib.sha256(
                normalized_sequence.encode("utf-8")
            ).hexdigest()
            sequence_sha256_to_hits.setdefault(digest, []).append(
                (str(accession_value), str(uid_value))
            )
        for digest, hits in sequence_sha256_to_hits.items():
            hits.sort()

    reference_has_sequence_sha256 = (
        REFERENCE_SEQUENCE_SHA256_COLUMN in reference_dataframe.columns
    )

    detail_rows: list[dict[str, object]] = []
    for _, reference_row in reference_dataframe.iterrows():
        reference_accession = _normalize_accession(reference_row["accession"])
        reference_bare = _bare_accession(reference_accession)
        reference_sequence_sha256 = (
            str(reference_row[REFERENCE_SEQUENCE_SHA256_COLUMN]).strip().lower()
            if reference_has_sequence_sha256
            else ""
        )

        match_method = ReferenceMatchMethod.NONE
        matched_accession = ""
        matched_protein_uid = ""
        sequence_match_count = 0

        if reference_accession and reference_accession in retrieved_versioned_set:
            match_method = ReferenceMatchMethod.EXACT_ACCESSION_VERSION
            matched_accession = reference_accession
            matched_protein_uid = versioned_to_uid.get(reference_accession, "")
        elif reference_bare and reference_bare in retrieved_by_bare:
            match_method = ReferenceMatchMethod.SAME_BASE_ACCESSION
            matched_accession, matched_protein_uid = retrieved_by_bare[
                reference_bare
            ]
        elif (
            reference_sequence_sha256
            and reference_sequence_sha256 in sequence_sha256_to_hits
        ):
            match_method = ReferenceMatchMethod.SEQUENCE_SHA256
            hits = sequence_sha256_to_hits[reference_sequence_sha256]
            matched_accession, matched_protein_uid = hits[0]
            sequence_match_count = len(hits)

        recovered = match_method != ReferenceMatchMethod.NONE
        recovered_exact_accession = (
            match_method.value in _STRICT_ACCESSION_METHODS
        )

        detail_rows.append(
            {
                "accession": reference_accession,
                "clade": str(reference_row["clade"]).strip().upper(),
                "ago_family": str(reference_row["ago_family"]).strip().upper(),
                "sequence_sha256": reference_sequence_sha256,
                "recovered": bool(recovered),
                "recovered_exact_accession": bool(recovered_exact_accession),
                "match_method": match_method.value,
                "matched_accession": matched_accession,
                "matched_protein_uid": matched_protein_uid,
                "sequence_match_count": int(sequence_match_count),
                "reference_label_source": reference_row["reference_label_source"],
                "reference_label_evidence": reference_row[
                    "reference_label_evidence"
                ],
            }
        )

    detail_dataframe = pd.DataFrame(detail_rows)

    summary_rows: list[dict[str, object]] = []
    stratum_exact_recall: dict[str, float | None] = {}
    stratum_equivalent_recall: dict[str, float | None] = {}
    stratum_recall_status: dict[str, str] = {}
    for stratum_name, metric_name, detail_column, detail_value in RECALL_STRATA:
        if detail_column is None:
            stratum_frame = detail_dataframe
        else:
            stratum_frame = detail_dataframe[
                detail_dataframe[detail_column] == detail_value
            ]
        reference_count = int(len(stratum_frame))
        exact_recovered_count = int(
            stratum_frame["recovered_exact_accession"].sum()
        )
        equivalent_recovered_count = int(stratum_frame["recovered"].sum())

        exact_recall_value, recall_status = _recall(
            reference_count, exact_recovered_count
        )
        equivalent_recall_value, _ = _recall(
            reference_count, equivalent_recovered_count
        )

        stratum_exact_recall[metric_name] = exact_recall_value
        stratum_equivalent_recall[metric_name] = equivalent_recall_value
        stratum_recall_status[metric_name] = recall_status
        summary_rows.append(
            {
                "stratum": stratum_name,
                "metric_name": metric_name,
                "reference_count": reference_count,
                "exact_recovered_count": exact_recovered_count,
                "equivalent_recovered_count": equivalent_recovered_count,
                "exact_accession_recall": (
                    float("nan")
                    if exact_recall_value is None
                    else exact_recall_value
                ),
                "retrieval_equivalent_recall": (
                    float("nan")
                    if equivalent_recall_value is None
                    else equivalent_recall_value
                ),
                "recall_status": recall_status,
            }
        )

    summary_dataframe = pd.DataFrame(summary_rows)

    exact_recovered_total = int(detail_dataframe["recovered_exact_accession"].sum())
    equivalent_recovered_total = int(detail_dataframe["recovered"].sum())

    if progress_callback is not None:
        progress_callback(
            "query reference recall: "
            f"exact {exact_recovered_total}/{len(detail_dataframe)}, "
            f"retrieval-equivalent {equivalent_recovered_total}/"
            f"{len(detail_dataframe)}"
        )

    return QueryReferenceRecallResult(
        summary_dataframe=summary_dataframe,
        detail_dataframe=detail_dataframe,
        stratum_exact_recall=stratum_exact_recall,
        stratum_equivalent_recall=stratum_equivalent_recall,
        stratum_recall_status=stratum_recall_status,
        reference_count=int(len(detail_dataframe)),
        exact_recovered_count=exact_recovered_total,
        equivalent_recovered_count=equivalent_recovered_total,
    )
