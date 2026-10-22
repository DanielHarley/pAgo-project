"""Compute NCBI query recall against the curated query-recall panel (ADR-0005).

This module answers one narrow question: how many rows of the versioned
query-recall reference panel does a given NCBI retrieval recover? It reports
two readings per stratum:

- ``same_accession_recall`` — recovery demonstrated by the admissible
  identifiers available in the retrieved metadata (tiers
  ``EXACT_ACCESSION_VERSION`` and ``SAME_BASE_ACCESSION``). A negative is not
  proof that the record was absent from the retrieval: a record that was in
  fact retrieved can go unrecognized when its identifier is missing or
  inadmissible. ``retrieved_inadmissible_accession_count`` is the qualifier.
- ``retrieval_equivalent_recall`` — the union of all three admitted tiers. A
  value of 1.0 does **not** establish that the retrieved sequences are
  identical to the reference sequences: a row may be counted through accession
  identity alone, under a different record version, carrying a different
  sequence. ``match_method`` is the only field that states how a row was
  recovered.

The denominator is the number of **panel rows** in the stratum — not distinct
accessions and not distinct protein sequences. Two panel rows may legitimately
share one protein sequence, in which case one retrieved sequence recovers both;
such groups are reported rather than rejected.

The panel is a query-coverage audit instrument only (``intended_use =
NCBI_QUERY_COVERAGE_AUDIT``, ``detector_validation_eligible = false``). This
module refuses to run against a panel whose lock declares otherwise. That is a
contract the resource declares and this module checks for self-consistency; it
cannot prevent a future consumer from misusing either the panel or these
outputs.

The panel's 21 sequence hashes are historically derived and frozen, not
re-derived here. A ``SEQUENCE_SHA256`` match is only as good as the stored
reference hash.

No network access, no persistence: see ``query_reference_recall_snapshot``.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from enum import StrEnum

import pandas as pd

from src.pago_pipeline.query_recall_reference_panel import (
    QueryRecallReferencePanelValidation,
)

DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN = "protein_uid"
DEFAULT_RETRIEVED_ACCESSION_COLUMN = "gbseq__accession_version"
DEFAULT_RETRIEVED_SEQUENCE_COLUMN = "gbseq__sequence"

# An accession is admissible to the accession tiers only in explicit
# accession.version form. Without the version suffix the identifier denotes
# "the current version" rather than one immutable record, so EXACT_ACCESSION_
# VERSION could never fire and SAME_BASE_ACCESSION would silently answer in its
# place (see ADR on evidence admissibility).
ADMISSIBLE_ACCESSION = re.compile(r"^[A-Za-z0-9_]+\.\d+$")
# An NCBI protein UID is a decimal integer. `str.isdigit()` is not usable here:
# it accepts "²" (which int() rejects) and "٣" (which int() accepts as 3).
ADMISSIBLE_UID = re.compile(r"^[0-9]+$")
_VERSION_SUFFIX_PATTERN = re.compile(r"\.\d+$")

SEQUENCE_NORMALIZATION = "strip_all_whitespace_then_uppercase"
STRATEGY_KIND = "query_reference_recall_matching"
STRATEGY_VERSION = "3.0"
MATCHING_POLICY_VERSION = "1.0"
REFERENCE_SCHEMA = "query_recall_reference_panel_lock/1.0"

RECALL_STATUS_EVALUABLE = "EVALUABLE"
RECALL_STATUS_NOT_EVALUABLE = "NOT_EVALUABLE"

REQUIRED_INTENDED_USE = "NCBI_QUERY_COVERAGE_AUDIT"


class ReferenceMatchMethod(StrEnum):
    """How a panel row was recovered, ordered by decreasing identifier strictness."""

    EXACT_ACCESSION_VERSION = "EXACT_ACCESSION_VERSION"
    SAME_BASE_ACCESSION = "SAME_BASE_ACCESSION"
    SEQUENCE_SHA256 = "SEQUENCE_SHA256"
    NONE = "NONE"


# The two identifier-strict tiers: "the retrieval returned the same accession",
# which deliberately allows a different version of that accession.
_SAME_ACCESSION_METHODS = frozenset(
    {
        ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value,
        ReferenceMatchMethod.SAME_BASE_ACCESSION.value,
    }
)

# (display_name, metric_name, required_ago_family, required_pago_clade); None
# means "do not filter on this field". The clade strata also require
# ago_family == PAGO because pago_clade is meaningful only within PAGO
# (ADR-0012), and PIWI_RE rows carry pago_clade = UNRESOLVED, so a clade filter
# alone would mix them with the quarantined PAGO rows (ADR-0011).
RECALL_STRATA: tuple[tuple[str, str, str | None, str | None], ...] = (
    ("overall", "overall_reference_recall", None, None),
    ("LONG_A", "long_a_reference_recall", "PAGO", "LONG_A"),
    ("LONG_B", "long_b_reference_recall", "PAGO", "LONG_B"),
    ("SHORT", "short_reference_recall", "PAGO", "SHORT"),
    (
        "PAGO_UNRESOLVED_CLADE",
        "pago_unresolved_clade_reference_recall",
        "PAGO",
        "UNRESOLVED",
    ),
    ("PIWI_RE", "piwi_re_reference_recall", "PIWI_RE", None),
)

SUMMARY_FIELDNAMES: tuple[str, ...] = (
    "stratum",
    "metric_name",
    "reference_count",
    "exact_version_match_count",
    "same_base_match_count",
    "same_accession_recovered_count",
    "equivalent_recovered_count",
    "same_accession_recall",
    "retrieval_equivalent_recall",
    "recall_status",
    "retrieval_equivalent_recall_is_lower_bound",
    "equivalent_recall_counts_shared_sequences",
)

DETAIL_FIELDNAMES: tuple[str, ...] = (
    "accession",
    "protein_short_name",
    "ago_family",
    "pago_clade",
    "verification_status",
    "accession_curation_status",
    "sequence_sha256",
    "recovered",
    "recovered_same_accession",
    "match_method",
    "matched_accession",
    "matched_protein_uid",
    "accession_version_distinct_uid_count",
    "sequence_match_count",
    "reference_label_source",
    "reference_label_evidence",
)


@dataclass(frozen=True)
class QueryReferenceRecallResult:
    summary_rows: tuple[dict[str, object], ...]
    detail_rows: tuple[dict[str, object], ...]
    stratum_same_accession_recall: dict[str, float | None]
    stratum_equivalent_recall: dict[str, float | None]
    stratum_recall_status: dict[str, str]
    stratum_equivalent_is_lower_bound: dict[str, bool]
    reference_count: int
    same_accession_recovered_count: int
    equivalent_recovered_count: int
    retrieved_row_count: int
    retrieved_rows_with_sequence_count: int
    sequence_evidence_complete: bool
    retrieved_inadmissible_accession_count: int
    retrieved_unattributable_row_count: int
    duplicate_reference_sequence_groups: tuple[tuple[str, ...], ...]


# --------------------------------------------------------------------------
# Cell contract: reject non-scalars, then classify missingness, then require str
# --------------------------------------------------------------------------


def _require_scalar_cell(value: object, *, column_name: str) -> None:
    if not pd.api.types.is_scalar(value):
        raise RuntimeError(
            f"Retrieved metadata column {column_name!r} contains a non-scalar "
            f"cell of type {type(value).__name__}; one scalar per cell is "
            "required."
        )


def _is_missing(value: object) -> bool:
    """Missingness for a value already established to be scalar."""
    return value is None or bool(pd.isna(value))


def _cell_text(value: object, *, column_name: str) -> str:
    """Return the cell's text, or "" when it is missing.

    A container is rejected before missingness is considered, so neither
    ``[None]`` nor an all-NaN array can be mistaken for an absent cell. A
    present value must already be a string: the documented input is
    ``read_csv(..., dtype=str)``, so a third type means the caller is outside
    the contract, and coercing it would turn ``1.0`` into an accession.
    """
    _require_scalar_cell(value, column_name=column_name)
    if _is_missing(value):
        return ""
    if not isinstance(value, str):
        raise RuntimeError(
            f"Retrieved metadata column {column_name!r} contains a "
            f"{type(value).__name__} cell ({value!r}); string cells are "
            "required. Read the metadata CSV with dtype=str."
        )
    return value


def normalize_protein_sequence(
    value: object,
    *,
    column_name: str = DEFAULT_RETRIEVED_SEQUENCE_COLUMN,
) -> str:
    """Strip every whitespace character, then uppercase. Nothing else."""
    return "".join(_cell_text(value, column_name=column_name).split()).upper()


def protein_sequence_sha256(
    value: object,
    *,
    column_name: str = DEFAULT_RETRIEVED_SEQUENCE_COLUMN,
) -> str:
    normalized = normalize_protein_sequence(value, column_name=column_name)
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def _normalize_accession(value: object, *, column_name: str) -> str:
    return _cell_text(value, column_name=column_name).strip()


def _accession_is_admissible(accession: str) -> bool:
    return bool(ADMISSIBLE_ACCESSION.fullmatch(accession))


def _uid_is_admissible(uid: str) -> bool:
    return bool(ADMISSIBLE_UID.fullmatch(uid))


def _bare_accession(accession: str) -> str:
    return _VERSION_SUFFIX_PATTERN.sub("", accession)


def _uid_sort_key(uid: str) -> tuple[int, int, str]:
    """Total order over UID tokens.

    Admissible UIDs order by numeric value, with the raw text breaking a tie
    only when the numeric values are equal — so "01" and "1" are distinguished
    without either being canonicalized away. Present-but-inadmissible tokens
    and missing tokens remain orderable, in that order, so the key is total
    over every input; orderability is not admissibility (see
    ``_is_attributable``).
    """
    if _uid_is_admissible(uid):
        return (0, int(uid), uid)
    if uid:
        return (1, 0, uid)
    return (2, 0, "")


def _is_attributable(accession: str, uid: str) -> bool:
    """Whether a retrieved observation can be attributed to an NCBI record."""
    return _accession_is_admissible(accession) or _uid_is_admissible(uid)


def _attribution_sort_key(accession: str, uid: str) -> tuple[int, str, tuple[int, int, str]]:
    """Representative order for sequence hits: admissible accessions lead.

    Leading on admissibility rather than on emptiness matters: byte order alone
    would rank "!bad" above "A.1".
    """
    accession_admissible = _accession_is_admissible(accession)
    return (
        0 if accession_admissible else 1,
        accession if accession_admissible else "",
        _uid_sort_key(uid),
    )


def _emitted_accession(accession: str) -> str:
    return accession if _accession_is_admissible(accession) else ""


def _emitted_uid(uid: str) -> str:
    return uid if _uid_is_admissible(uid) else ""


# --------------------------------------------------------------------------
# Methodological identity
# --------------------------------------------------------------------------


def build_matching_strategy_payload() -> dict[str, object]:
    """Every policy that can change a result or an evaluation state.

    A field is listed here precisely so that ``matching_strategy_sha256``
    covers the methodology it claims to cover: a stored snapshot computed under
    a different policy must not look compatible.
    """
    return {
        "strategy_kind": STRATEGY_KIND,
        "strategy_version": STRATEGY_VERSION,
        "matching_policy_version": MATCHING_POLICY_VERSION,
        "reference_schema": REFERENCE_SCHEMA,
        "match_hierarchy": [
            ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value,
            ReferenceMatchMethod.SAME_BASE_ACCESSION.value,
            ReferenceMatchMethod.SEQUENCE_SHA256.value,
        ],
        "admissible_retrieved_accession_pattern": ADMISSIBLE_ACCESSION.pattern,
        "admissible_retrieved_uid_pattern": ADMISSIBLE_UID.pattern,
        "inadmissible_accession_policy": (
            "excluded_from_accession_tiers_sequence_still_eligible"
        ),
        "attribution_rule": "admissible_accession_or_admissible_uid",
        "unattributable_sequence_observation_policy": (
            "excluded_from_index_and_from_usable_evidence"
        ),
        "emitted_identifier_policy": "inadmissible_token_emitted_as_empty",
        "duplicate_reference_base_accession_policy": "fatal_instrument_policy",
        "duplicate_reference_sequence_hash_policy": "permitted_and_reported",
        "conflicting_sequence_under_one_accession_version_policy": (
            "fatal_integrity_error"
        ),
        "missing_scalar_policy": "scalar_check_then_isna_gate_then_str_only_else_reject",
        "sequence_normalization": SEQUENCE_NORMALIZATION,
        "sequence_hash": "sha256",
        "incomplete_sequence_evidence_policy": (
            "tracked_counts_and_per_stratum_lower_bound_flag"
        ),
        "lower_bound_saturation_rule": (
            "not_complete_and_equivalent_recovered_lt_reference_count"
        ),
        "exact_accession_representative": "min_uid_numeric_then_raw_text",
        "base_accession_representative": (
            "min_versioned_accession_lexicographic_then_min_uid"
        ),
        "multi_sequence_hit_representative": (
            "admissible_accession_first_then_min_accession_then_min_uid"
        ),
        "recall_readings": ["same_accession_recall", "retrieval_equivalent_recall"],
        "strata": [
            {"stratum": name, "ago_family": family, "pago_clade": clade}
            for name, _metric, family, clade in RECALL_STRATA
        ],
    }


def build_matching_strategy_sha256() -> str:
    serialized = json.dumps(
        build_matching_strategy_payload(), sort_keys=True, separators=(",", ":")
    )
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


# --------------------------------------------------------------------------
# Panel preconditions
# --------------------------------------------------------------------------


def _validate_panel_use_contract(lock_payload: Mapping[str, object]) -> None:
    if lock_payload.get("intended_use") != REQUIRED_INTENDED_USE:
        raise RuntimeError(
            "Query-recall panel lock intended_use must be "
            f"{REQUIRED_INTENDED_USE!r} to compute a query-coverage audit, got "
            f"{lock_payload.get('intended_use')!r}."
        )
    if lock_payload.get("detector_validation_eligible") is not False:
        raise RuntimeError(
            "Query-recall panel lock detector_validation_eligible must be "
            f"false, got {lock_payload.get('detector_validation_eligible')!r}. "
            "This panel is a query-coverage audit instrument, not a "
            "detector-validation set."
        )


def _validate_reference_identity(
    panel_validation: QueryRecallReferencePanelValidation,
) -> tuple[tuple[str, ...], ...]:
    """Check the matcher's preconditions on the panel and group shared hashes.

    Versioned accessions and unique base accessions are preconditions of the
    tier-1/tier-2 distinction. Duplicate sequence hashes are *not* an error:
    distinct accessions can legitimately encode identical proteins, so they are
    grouped and reported instead.
    """
    rows = panel_validation.rows
    for row in rows:
        if not _accession_is_admissible(row.accession):
            raise RuntimeError(
                f"Reference accession {row.accession!r} ({row.protein_short_name}) "
                "is not an explicit accession.version identifier matching "
                f"{ADMISSIBLE_ACCESSION.pattern!r}. Without the version the "
                "exact-version tier can never apply and the same-base tier "
                "would answer in its place."
            )

    base_counts = Counter(_bare_accession(row.accession) for row in rows)
    shared_bases = sorted(base for base, count in base_counts.items() if count > 1)
    if shared_bases:
        raise RuntimeError(
            "Reference rows share a base accession: "
            f"{shared_bases}. This panel's policy requires one row per base "
            "accession, because a single retrieved record would otherwise "
            "satisfy the same-base tier for two rows and be counted twice in "
            "the strict reading. Permitting such pairs requires a decision "
            "that redefines the same-base tier to bind at most one reference "
            "per base accession."
        )

    hash_counts = Counter(row.sequence_sha256 for row in rows)
    return tuple(
        tuple(
            sorted(row.accession for row in rows if row.sequence_sha256 == digest)
        )
        for digest in sorted(
            digest for digest, count in hash_counts.items() if count > 1
        )
    )


# --------------------------------------------------------------------------
# Retrieved-side indices
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class _RetrievedIndices:
    exact_uids_by_accession: dict[str, tuple[str, ...]]
    representative_by_base: dict[str, tuple[str, str]]
    hits_by_sequence_sha256: dict[str, tuple[tuple[str, str], ...]]
    row_count: int
    rows_with_sequence_count: int
    inadmissible_accession_count: int
    unattributable_row_count: int


def _build_retrieved_indices(
    *,
    retrieved_metadata_dataframe: pd.DataFrame,
    accession_column: str,
    protein_uid_column: str,
    sequence_column: str,
) -> _RetrievedIndices:
    for required_column in (accession_column, sequence_column):
        if required_column not in retrieved_metadata_dataframe.columns:
            raise RuntimeError(
                f"Retrieved metadata is missing the required column "
                f"{required_column!r}. Both {accession_column!r} and "
                f"{sequence_column!r} are required: the sequence tier is part "
                "of the decided matching hierarchy (ADR-0005), so a run "
                "without it does not produce the documented metric."
            )

    accession_values = retrieved_metadata_dataframe[accession_column]
    sequence_values = retrieved_metadata_dataframe[sequence_column]
    if protein_uid_column in retrieved_metadata_dataframe.columns:
        uid_values = retrieved_metadata_dataframe[protein_uid_column]
    else:
        uid_values = pd.Series([None] * len(retrieved_metadata_dataframe))

    uids_by_accession: dict[str, list[str]] = {}
    sequence_hashes_by_accession: dict[str, set[str]] = {}
    representative_by_base: dict[str, tuple[str, str]] = {}
    hits_by_sequence: dict[str, list[tuple[str, str]]] = {}
    row_count = 0
    rows_with_sequence_count = 0
    inadmissible_accession_count = 0
    unattributable_row_count = 0

    for raw_accession, raw_uid, raw_sequence in zip(
        accession_values, uid_values, sequence_values, strict=True
    ):
        row_count += 1
        accession = _normalize_accession(raw_accession, column_name=accession_column)
        uid = _cell_text(raw_uid, column_name=protein_uid_column).strip()
        normalized_sequence = normalize_protein_sequence(
            raw_sequence, column_name=sequence_column
        )

        accession_admissible = _accession_is_admissible(accession)
        if not accession_admissible:
            inadmissible_accession_count += 1
        attributable = _is_attributable(accession, uid)
        if not attributable:
            unattributable_row_count += 1

        if accession_admissible:
            uids_by_accession.setdefault(accession, []).append(uid)
            if normalized_sequence:
                sequence_hashes_by_accession.setdefault(accession, set()).add(
                    hashlib.sha256(normalized_sequence.encode("utf-8")).hexdigest()
                )
            base = _bare_accession(accession)
            candidate = (accession, uid)
            existing = representative_by_base.get(base)
            if existing is None or (
                accession,
                _uid_sort_key(uid),
            ) < (existing[0], _uid_sort_key(existing[1])):
                representative_by_base[base] = candidate

        if normalized_sequence and attributable:
            rows_with_sequence_count += 1
            digest = hashlib.sha256(normalized_sequence.encode("utf-8")).hexdigest()
            hits_by_sequence.setdefault(digest, []).append((accession, uid))

    conflicting = sorted(
        accession
        for accession, digests in sequence_hashes_by_accession.items()
        if len(digests) > 1
    )
    if conflicting:
        raise RuntimeError(
            "Retrieved metadata carries conflicting protein sequences under a "
            f"single accession.version: {conflicting}. One accession.version "
            "denotes one immutable record, so this is a metadata integrity "
            "error rather than a tie to be broken."
        )

    return _RetrievedIndices(
        exact_uids_by_accession={
            accession: tuple(sorted(set(uids), key=_uid_sort_key))
            for accession, uids in uids_by_accession.items()
        },
        representative_by_base=representative_by_base,
        hits_by_sequence_sha256={
            digest: tuple(
                sorted(hits, key=lambda hit: _attribution_sort_key(hit[0], hit[1]))
            )
            for digest, hits in hits_by_sequence.items()
        },
        row_count=row_count,
        rows_with_sequence_count=rows_with_sequence_count,
        inadmissible_accession_count=inadmissible_accession_count,
        unattributable_row_count=unattributable_row_count,
    )


def _recall(reference_count: int, recovered_count: int) -> tuple[float | None, str]:
    if reference_count == 0:
        return None, RECALL_STATUS_NOT_EVALUABLE
    return float(recovered_count / reference_count), RECALL_STATUS_EVALUABLE


# --------------------------------------------------------------------------
# Public entry point
# --------------------------------------------------------------------------


def compute_query_reference_recall(
    *,
    panel_validation: QueryRecallReferencePanelValidation,
    retrieved_metadata_dataframe: pd.DataFrame,
    retrieved_accession_column: str = DEFAULT_RETRIEVED_ACCESSION_COLUMN,
    retrieved_protein_uid_column: str = DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN,
    retrieved_sequence_column: str = DEFAULT_RETRIEVED_SEQUENCE_COLUMN,
    progress_callback: Callable[[str], None] | None = None,
) -> QueryReferenceRecallResult:
    """Measure panel-row coverage of a retrieval, overall and per stratum.

    Raises ``RuntimeError`` on a panel precondition violation, a missing
    required retrieved column, a non-scalar or non-string retrieved cell, or
    conflicting sequences under one retrieved accession.version. Performs no
    network access and writes nothing.
    """
    _validate_panel_use_contract(panel_validation.lock_payload)
    duplicate_groups = _validate_reference_identity(panel_validation)

    indices = _build_retrieved_indices(
        retrieved_metadata_dataframe=retrieved_metadata_dataframe,
        accession_column=retrieved_accession_column,
        protein_uid_column=retrieved_protein_uid_column,
        sequence_column=retrieved_sequence_column,
    )
    sequence_evidence_complete = (
        indices.rows_with_sequence_count == indices.row_count
    )

    detail_rows: list[dict[str, object]] = []
    for reference_row in panel_validation.rows:
        accession = reference_row.accession
        match_method = ReferenceMatchMethod.NONE
        matched_accession = ""
        matched_uid = ""
        sequence_match_count = 0

        if accession in indices.exact_uids_by_accession:
            match_method = ReferenceMatchMethod.EXACT_ACCESSION_VERSION
            matched_accession = accession
            uids = indices.exact_uids_by_accession[accession]
            matched_uid = uids[0] if uids else ""
        else:
            base = _bare_accession(accession)
            if base in indices.representative_by_base:
                match_method = ReferenceMatchMethod.SAME_BASE_ACCESSION
                matched_accession, matched_uid = indices.representative_by_base[base]
            elif reference_row.sequence_sha256 in indices.hits_by_sequence_sha256:
                match_method = ReferenceMatchMethod.SEQUENCE_SHA256
                hits = indices.hits_by_sequence_sha256[
                    reference_row.sequence_sha256
                ]
                matched_accession, matched_uid = hits[0]
                sequence_match_count = len(hits)

        distinct_uid_count = len(
            {
                uid
                for uid in indices.exact_uids_by_accession.get(matched_accession, ())
                if uid
            }
        )
        recovered = match_method is not ReferenceMatchMethod.NONE
        detail_rows.append(
            {
                "accession": accession,
                "protein_short_name": reference_row.protein_short_name,
                "ago_family": reference_row.ago_family,
                "pago_clade": reference_row.pago_clade,
                "verification_status": reference_row.verification_status,
                "accession_curation_status": reference_row.accession_curation_status,
                "sequence_sha256": reference_row.sequence_sha256,
                "recovered": recovered,
                "recovered_same_accession": match_method.value
                in _SAME_ACCESSION_METHODS,
                "match_method": match_method.value,
                "matched_accession": _emitted_accession(matched_accession),
                "matched_protein_uid": _emitted_uid(matched_uid),
                "accession_version_distinct_uid_count": distinct_uid_count,
                "sequence_match_count": sequence_match_count,
                "reference_label_source": reference_row.reference_label_source,
                "reference_label_evidence": reference_row.reference_label_evidence,
            }
        )

    duplicated_accessions = {
        accession for group in duplicate_groups for accession in group
    }

    summary_rows: list[dict[str, object]] = []
    stratum_same_accession_recall: dict[str, float | None] = {}
    stratum_equivalent_recall: dict[str, float | None] = {}
    stratum_recall_status: dict[str, str] = {}
    stratum_equivalent_is_lower_bound: dict[str, bool] = {}
    partitioned_reference_count = 0

    for stratum_name, metric_name, required_family, required_clade in RECALL_STRATA:
        stratum_detail = [
            row
            for row in detail_rows
            if (required_family is None or row["ago_family"] == required_family)
            and (required_clade is None or row["pago_clade"] == required_clade)
        ]
        reference_count = len(stratum_detail)
        if stratum_name != "overall":
            partitioned_reference_count += reference_count

        exact_version_match_count = sum(
            1
            for row in stratum_detail
            if row["match_method"]
            == ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value
        )
        same_base_match_count = sum(
            1
            for row in stratum_detail
            if row["match_method"] == ReferenceMatchMethod.SAME_BASE_ACCESSION.value
        )
        same_accession_recovered_count = (
            exact_version_match_count + same_base_match_count
        )
        equivalent_recovered_count = sum(
            1 for row in stratum_detail if row["recovered"]
        )

        same_accession_recall, recall_status = _recall(
            reference_count, same_accession_recovered_count
        )
        equivalent_recall, _ = _recall(reference_count, equivalent_recovered_count)

        # Saturation: a stratum whose every reference is already recovered has a
        # determined equivalent reading, however incomplete the sequence
        # evidence is, because the metric cannot rise further.
        is_lower_bound = bool(
            reference_count
            and not sequence_evidence_complete
            and equivalent_recovered_count < reference_count
        )
        recovered_shared = Counter(
            row["sequence_sha256"]
            for row in stratum_detail
            if row["recovered"] and row["accession"] in duplicated_accessions
        )
        counts_shared_sequences = any(
            count >= 2 for count in recovered_shared.values()
        )

        stratum_same_accession_recall[metric_name] = same_accession_recall
        stratum_equivalent_recall[metric_name] = equivalent_recall
        stratum_recall_status[metric_name] = recall_status
        stratum_equivalent_is_lower_bound[metric_name] = is_lower_bound
        summary_rows.append(
            {
                "stratum": stratum_name,
                "metric_name": metric_name,
                "reference_count": reference_count,
                "exact_version_match_count": exact_version_match_count,
                "same_base_match_count": same_base_match_count,
                "same_accession_recovered_count": same_accession_recovered_count,
                "equivalent_recovered_count": equivalent_recovered_count,
                "same_accession_recall": same_accession_recall,
                "retrieval_equivalent_recall": equivalent_recall,
                "recall_status": recall_status,
                "retrieval_equivalent_recall_is_lower_bound": is_lower_bound,
                "equivalent_recall_counts_shared_sequences": counts_shared_sequences,
            }
        )

    if partitioned_reference_count != len(detail_rows):
        raise RuntimeError(
            "The stratum definitions do not partition the panel: the "
            f"non-overall strata cover {partitioned_reference_count} rows of "
            f"{len(detail_rows)}. Every reference must fall in exactly one "
            "stratum, or stratified recall silently omits references."
        )

    same_accession_total = sum(
        1 for row in detail_rows if row["recovered_same_accession"]
    )
    equivalent_total = sum(1 for row in detail_rows if row["recovered"])

    if progress_callback is not None:
        # The message reports the overall stratum, so its lower-bound note must
        # come from that stratum's computed state. Incomplete sequence evidence
        # alone is not enough: a saturated reading is determined.
        overall_is_lower_bound = stratum_equivalent_is_lower_bound[
            "overall_reference_recall"
        ]
        progress_callback(
            "query reference recall: same-accession "
            f"{same_accession_total}/{len(detail_rows)}, retrieval-equivalent "
            f"{equivalent_total}/{len(detail_rows)}"
            + (" (lower bound)" if overall_is_lower_bound else "")
        )

    return QueryReferenceRecallResult(
        summary_rows=tuple(summary_rows),
        detail_rows=tuple(detail_rows),
        stratum_same_accession_recall=stratum_same_accession_recall,
        stratum_equivalent_recall=stratum_equivalent_recall,
        stratum_recall_status=stratum_recall_status,
        stratum_equivalent_is_lower_bound=stratum_equivalent_is_lower_bound,
        reference_count=len(detail_rows),
        same_accession_recovered_count=same_accession_total,
        equivalent_recovered_count=equivalent_total,
        retrieved_row_count=indices.row_count,
        retrieved_rows_with_sequence_count=indices.rows_with_sequence_count,
        sequence_evidence_complete=sequence_evidence_complete,
        retrieved_inadmissible_accession_count=indices.inadmissible_accession_count,
        retrieved_unattributable_row_count=indices.unattributable_row_count,
        duplicate_reference_sequence_groups=duplicate_groups,
    )
