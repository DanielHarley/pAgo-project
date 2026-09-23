"""Load and validate the curated query-recall reference panel (ADR-0005).

The panel is a **query-coverage audit instrument**: it measures whether an
NCBI text query recovers known pAgo/PIWI-RE proteins. It is not a classifier
gold standard, detector benchmark, detector CALIBRATION set, detector
FINAL_HOLDOUT, or generic pAgo ground truth (ADR-0005). The PIWI-RE rows
carry an additional, pre-existing circularity restriction (ADR-0014): most
were selected using PIWI-RE-specific profile evidence, so they must never be
used to validate a profile-based PIWI-RE detector. Both restrictions are
enforced here as machine-checkable invariants (`intended_use`,
`forbidden_uses`, `detector_validation_eligible` in the lock), not only in
prose.

Each row separates two independent judgments: `verification_status` (evidence
strength for the protein's biology) and `accession_curation_status` (whether
this specific accession.version is the settled, canonical identifier for the
reference). TtAgo is `verification_status = verified` — its biology is not in
question — but `accession_curation_status = UNRESOLVED`, since the panel's
`WP_011174533.1` is not confirmed as the crystallographic HB8 reference
(SCIENTIFIC_STATE.md, Open scientific question 2).

This module only loads and validates the panel. It does not compute
query-recall (matching against retrieved records, recall metrics) — that is
separate, future work building on this resource.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

from src.pago_pipeline.storage import read_json_file, sha256_of_file

PathLike: TypeAlias = str | Path

PANEL_LOCK_ARTIFACT_TYPE = "query_recall_reference_panel_lock"
PANEL_LOCK_FORMAT_VERSION = "1.0"
DEFAULT_PANEL_DIRECTORY = (
    Path(__file__).resolve().parent / "resources" / "query_recall"
)
LOCK_FILE_NAME = "panel_lock.json"

EXPECTED_REFERENCE_COUNT = 21
EXPECTED_COUNTS_BY_AGO_FAMILY: Mapping[str, int] = {"PAGO": 14, "PIWI_RE": 7}

REQUIRED_SCHEMA_FIELDS = (
    "accession",
    "protein_short_name",
    "organism",
    "ago_family",
    "pago_clade",
    "historical_query_recall_clade",
    "sequence_sha256",
    "sequence_length",
    "uniprot_accession",
    "reference_label_source",
    "reference_label_evidence",
    "verification_status",
    "accession_curation_status",
    "notes",
)
ALLOWED_AGO_FAMILY_VALUES = frozenset({"PAGO", "PIWI_RE"})
ALLOWED_PAGO_CLADE_VALUES = frozenset({"LONG_A", "LONG_B", "SHORT", "UNRESOLVED"})
ALLOWED_EVIDENCE_TIER_VALUES = frozenset(
    {
        "EXPERIMENTAL",
        "LITERATURE_PHYLOGENETIC",
        "CURATED_COMPUTATIONAL",
        "DATABASE_ANNOTATION",
    }
)
ALLOWED_VERIFICATION_STATUS_VALUES = frozenset({"verified", "provisional"})
# accession_curation_status is deliberately independent of verification_status:
# verification_status describes evidence strength for the protein's biology,
# while accession_curation_status describes whether *this accession* is the
# settled, canonical identifier for the reference. A row can be biologically
# `verified` while its accession identity is still `UNRESOLVED` (TtAgo).
ALLOWED_ACCESSION_CURATION_STATUS_VALUES = frozenset({"RESOLVED", "UNRESOLVED"})
REQUIRED_FORBIDDEN_USES = frozenset(
    {
        "DETECTOR_VALIDATION",
        "DETECTOR_CALIBRATION",
        "DETECTOR_FINAL_HOLDOUT",
        "CLASSIFIER_GOLD_STANDARD",
    }
)
REQUIRED_INTENDED_USE = "NCBI_QUERY_COVERAGE_AUDIT"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")

# Known reconciliation invariants (ADR-0011 / SCIENTIFIC_STATE.md open
# question 2) that must never silently revert to their historical values.
KNOWN_PAGO_CLADE_RECONCILIATIONS: Mapping[str, str] = {
    "WP_005580376.1": "LONG_A",  # NgAgo: historical LONG_B was a probable curation error
    "WP_010878815.1": "UNRESOLVED",  # AfAgo: quarantined (architecture vs. Ryazansky conflict)
    "WP_012735993.1": "UNRESOLVED",  # SiAgo: quarantined (architecture vs. Ryazansky conflict)
}
TTAGO_ACCESSION = "WP_011174533.1"
# TtAgo's verification_status (biological evidence strength) is unchanged
# from the historical panel ("verified") — no authoritative decision changes
# that field's value for this row. The open question is accession identity
# only (SCIENTIFIC_STATE.md, Open scientific question 2), represented by the
# independent accession_curation_status field below.
TTAGO_REQUIRED_ACCESSION_CURATION_STATUS = "UNRESOLVED"

# Lock fields that only *describe* the committed CSV are cross-checked against
# the parsed rows, so a declared count, schema, or allowed-value set cannot
# silently drift from the panel it describes. `piwi_re_circularity_note` stays
# descriptive: the machine-readable contract is `intended_use`,
# `detector_validation_eligible`, and `forbidden_uses`, already validated above.
_DECLARED_COUNT_FIELDS: Mapping[str, tuple[str, frozenset[str]]] = {
    "counts_by_ago_family": ("ago_family", ALLOWED_AGO_FAMILY_VALUES),
    "counts_by_evidence_tier": (
        "reference_label_evidence",
        ALLOWED_EVIDENCE_TIER_VALUES,
    ),
    "counts_by_verification_status": (
        "verification_status",
        ALLOWED_VERIFICATION_STATUS_VALUES,
    ),
    "counts_by_accession_curation_status": (
        "accession_curation_status",
        ALLOWED_ACCESSION_CURATION_STATUS_VALUES,
    ),
}
_DECLARED_ALLOWED_VALUE_FIELDS: Mapping[str, frozenset[str]] = {
    "allowed_ago_family_values": ALLOWED_AGO_FAMILY_VALUES,
    "allowed_pago_clade_values": ALLOWED_PAGO_CLADE_VALUES,
    "allowed_evidence_tier_values": ALLOWED_EVIDENCE_TIER_VALUES,
    "allowed_verification_status_values": ALLOWED_VERIFICATION_STATUS_VALUES,
    "allowed_accession_curation_status_values": (
        ALLOWED_ACCESSION_CURATION_STATUS_VALUES
    ),
}


@dataclass(frozen=True)
class QueryRecallReferencePanelRow:
    accession: str
    protein_short_name: str
    organism: str
    ago_family: str
    pago_clade: str
    historical_query_recall_clade: str
    sequence_sha256: str
    sequence_length: int
    uniprot_accession: str
    reference_label_source: str
    reference_label_evidence: str
    verification_status: str
    accession_curation_status: str
    notes: str


@dataclass(frozen=True)
class QueryRecallReferencePanelValidation:
    panel_directory: Path
    lock_payload: Mapping[str, object]
    rows: tuple[QueryRecallReferencePanelRow, ...]
    panel_sha256: str
    curation_notes_sha256: str


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _require_sha256(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not _SHA256_PATTERN.match(value):
        raise RuntimeError(
            f"{field_name} must be a 64-character lowercase hex SHA-256, got {value!r}."
        )
    return value


def _parse_panel_rows(panel_file_path: Path) -> tuple[QueryRecallReferencePanelRow, ...]:
    with panel_file_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or tuple(reader.fieldnames) != REQUIRED_SCHEMA_FIELDS:
            raise RuntimeError(
                "Query-recall reference panel CSV header does not match the "
                f"required schema. Expected {REQUIRED_SCHEMA_FIELDS!r}, got "
                f"{tuple(reader.fieldnames) if reader.fieldnames else None!r}."
            )
        raw_rows = list(reader)

    rows: list[QueryRecallReferencePanelRow] = []
    for row_index, raw_row in enumerate(raw_rows, start=2):  # header is line 1
        accession = raw_row["accession"].strip()
        if not accession:
            raise RuntimeError(f"Panel row {row_index}: accession must not be empty.")

        ago_family = raw_row["ago_family"].strip()
        if ago_family not in ALLOWED_AGO_FAMILY_VALUES:
            raise RuntimeError(
                f"Panel row {row_index} ({accession}): ago_family {ago_family!r} "
                f"is not one of {sorted(ALLOWED_AGO_FAMILY_VALUES)}."
            )

        pago_clade = raw_row["pago_clade"].strip()
        if pago_clade not in ALLOWED_PAGO_CLADE_VALUES:
            raise RuntimeError(
                f"Panel row {row_index} ({accession}): pago_clade {pago_clade!r} "
                f"is not one of {sorted(ALLOWED_PAGO_CLADE_VALUES)}."
            )
        if ago_family == "PIWI_RE" and pago_clade != "UNRESOLVED":
            raise RuntimeError(
                f"Panel row {row_index} ({accession}): PIWI_RE rows must carry "
                f"pago_clade = UNRESOLVED (ADR-0012), got {pago_clade!r}."
            )

        evidence_tier = raw_row["reference_label_evidence"].strip()
        if evidence_tier not in ALLOWED_EVIDENCE_TIER_VALUES:
            raise RuntimeError(
                f"Panel row {row_index} ({accession}): reference_label_evidence "
                f"{evidence_tier!r} is not one of {sorted(ALLOWED_EVIDENCE_TIER_VALUES)}."
            )

        verification_status = raw_row["verification_status"].strip()
        if verification_status not in ALLOWED_VERIFICATION_STATUS_VALUES:
            raise RuntimeError(
                f"Panel row {row_index} ({accession}): verification_status "
                f"{verification_status!r} is not one of "
                f"{sorted(ALLOWED_VERIFICATION_STATUS_VALUES)}."
            )

        accession_curation_status = raw_row["accession_curation_status"].strip()
        if accession_curation_status not in ALLOWED_ACCESSION_CURATION_STATUS_VALUES:
            raise RuntimeError(
                f"Panel row {row_index} ({accession}): accession_curation_status "
                f"{accession_curation_status!r} is not one of "
                f"{sorted(ALLOWED_ACCESSION_CURATION_STATUS_VALUES)}."
            )

        sequence_sha256 = _require_sha256(
            raw_row["sequence_sha256"].strip(),
            field_name=f"panel row {row_index} ({accession}) sequence_sha256",
        )

        raw_length = raw_row["sequence_length"].strip()
        if not raw_length.isdigit() or int(raw_length) <= 0:
            raise RuntimeError(
                f"Panel row {row_index} ({accession}): sequence_length must be a "
                f"positive integer, got {raw_length!r}."
            )

        rows.append(
            QueryRecallReferencePanelRow(
                accession=accession,
                protein_short_name=raw_row["protein_short_name"].strip(),
                organism=raw_row["organism"].strip(),
                ago_family=ago_family,
                pago_clade=pago_clade,
                historical_query_recall_clade=raw_row[
                    "historical_query_recall_clade"
                ].strip(),
                sequence_sha256=sequence_sha256,
                sequence_length=int(raw_length),
                uniprot_accession=raw_row["uniprot_accession"].strip(),
                reference_label_source=raw_row["reference_label_source"].strip(),
                reference_label_evidence=evidence_tier,
                verification_status=verification_status,
                accession_curation_status=accession_curation_status,
                notes=raw_row["notes"].strip(),
            )
        )
    return tuple(rows)


def _validate_semantic_invariants(
    rows: tuple[QueryRecallReferencePanelRow, ...],
) -> None:
    if len(rows) != EXPECTED_REFERENCE_COUNT:
        raise RuntimeError(
            f"Query-recall reference panel must contain exactly "
            f"{EXPECTED_REFERENCE_COUNT} rows, got {len(rows)}."
        )

    counts_by_accession = Counter(row.accession for row in rows)
    duplicates = sorted(
        accession for accession, count in counts_by_accession.items() if count > 1
    )
    if duplicates:
        raise RuntimeError(
            f"Query-recall reference panel contains duplicate accessions: {duplicates}."
        )

    counts_by_family = dict(Counter(row.ago_family for row in rows))
    if counts_by_family != dict(EXPECTED_COUNTS_BY_AGO_FAMILY):
        raise RuntimeError(
            "Query-recall reference panel ago_family composition changed. "
            f"Expected {dict(EXPECTED_COUNTS_BY_AGO_FAMILY)}, got {counts_by_family}. "
            "Changing panel composition requires a new decision (see "
            "documentation/decisions/0005-query-recall-reference-panel-and-"
            "sequence-identity-matching.md); it must not happen silently."
        )

    rows_by_accession = {row.accession: row for row in rows}
    for accession, expected_pago_clade in KNOWN_PAGO_CLADE_RECONCILIATIONS.items():
        row = rows_by_accession.get(accession)
        if row is None:
            raise RuntimeError(
                f"Known reconciliation accession {accession!r} is missing from the panel."
            )
        if row.pago_clade != expected_pago_clade:
            raise RuntimeError(
                f"Panel row {accession} ({row.protein_short_name}): pago_clade "
                f"must remain {expected_pago_clade!r} per ADR-0011, got "
                f"{row.pago_clade!r}. This is a known reconciliation invariant "
                "and must not silently revert to its historical value."
            )

    ttago_row = rows_by_accession.get(TTAGO_ACCESSION)
    if ttago_row is None:
        raise RuntimeError(f"TtAgo accession {TTAGO_ACCESSION!r} is missing from the panel.")
    if ttago_row.accession_curation_status != TTAGO_REQUIRED_ACCESSION_CURATION_STATUS:
        raise RuntimeError(
            f"TtAgo ({TTAGO_ACCESSION}) accession_curation_status must remain "
            f"{TTAGO_REQUIRED_ACCESSION_CURATION_STATUS!r} until the open "
            "accession-identity question in SCIENTIFIC_STATE.md is resolved; "
            f"got {ttago_row.accession_curation_status!r}. Do not mark the "
            "accession identity resolved for convenience. This is independent "
            "of verification_status, which describes evidence strength for "
            "the protein's biology, not accession identity."
        )


def _validate_declared_lock_descriptions(
    lock_payload: Mapping[str, object],
    rows: tuple[QueryRecallReferencePanelRow, ...],
) -> None:
    """Cross-check the lock's descriptive fields against the parsed CSV.

    These fields describe the committed panel and are recomputable offline, so
    they must not be able to drift from it silently. Counts are compared over
    the full allowed-value set, so a declared zero (for example
    ``DATABASE_ANNOTATION``) is checked too. This adds no scientific invariant:
    it only proves the lock describes the CSV it ships with.
    """
    for lock_field, (row_field, allowed_values) in _DECLARED_COUNT_FIELDS.items():
        counts = Counter(getattr(row, row_field) for row in rows)
        actual = {value: counts[value] for value in sorted(allowed_values)}
        if lock_payload.get(lock_field) != actual:
            raise RuntimeError(
                f"Query-recall panel lock {lock_field} does not match the panel "
                f"CSV. Lock declares {lock_payload.get(lock_field)!r}, CSV has "
                f"{actual!r}."
            )

    declared_schema_fields = lock_payload.get("schema_fields")
    if tuple(declared_schema_fields or ()) != REQUIRED_SCHEMA_FIELDS:
        raise RuntimeError(
            "Query-recall panel lock schema_fields must match the validated CSV "
            f"schema {list(REQUIRED_SCHEMA_FIELDS)!r}, got "
            f"{declared_schema_fields!r}."
        )

    for lock_field, allowed_values in _DECLARED_ALLOWED_VALUE_FIELDS.items():
        declared = lock_payload.get(lock_field)
        if set(declared or ()) != set(allowed_values):
            raise RuntimeError(
                f"Query-recall panel lock {lock_field} must match the values the "
                f"loader enforces {sorted(allowed_values)!r}, got {declared!r}."
            )


def validate_query_recall_reference_panel(
    *,
    panel_directory: PathLike = DEFAULT_PANEL_DIRECTORY,
) -> QueryRecallReferencePanelValidation:
    """Validate the panel's lock, schema, and scientific invariants offline.

    Raises ``RuntimeError`` (or a subclass) on any violation. Performs no
    network access and no query-recall computation.
    """
    resolved_panel_directory = _as_path(panel_directory)
    lock_file_path = resolved_panel_directory / LOCK_FILE_NAME
    lock_payload = read_json_file(input_file_path=lock_file_path)
    if not isinstance(lock_payload, Mapping):
        raise RuntimeError("Query-recall panel lock must deserialize into a JSON object.")

    if lock_payload.get("artifact_type") != PANEL_LOCK_ARTIFACT_TYPE:
        raise RuntimeError(
            "Query-recall panel lock artifact_type mismatch. Expected "
            f"{PANEL_LOCK_ARTIFACT_TYPE!r}, got {lock_payload.get('artifact_type')!r}."
        )
    if lock_payload.get("format_version") != PANEL_LOCK_FORMAT_VERSION:
        raise RuntimeError(
            "Query-recall panel lock format_version mismatch. Expected "
            f"{PANEL_LOCK_FORMAT_VERSION!r}, got {lock_payload.get('format_version')!r}."
        )

    if lock_payload.get("intended_use") != REQUIRED_INTENDED_USE:
        raise RuntimeError(
            "Query-recall panel lock intended_use must be "
            f"{REQUIRED_INTENDED_USE!r}, got {lock_payload.get('intended_use')!r}."
        )
    declared_forbidden_uses = set(lock_payload.get("forbidden_uses") or [])
    if not REQUIRED_FORBIDDEN_USES.issubset(declared_forbidden_uses):
        raise RuntimeError(
            "Query-recall panel lock forbidden_uses must include at least "
            f"{sorted(REQUIRED_FORBIDDEN_USES)}, got {sorted(declared_forbidden_uses)}."
        )
    if lock_payload.get("detector_validation_eligible") is not False:
        raise RuntimeError(
            "Query-recall panel lock detector_validation_eligible must be "
            f"false, got {lock_payload.get('detector_validation_eligible')!r}."
        )

    panel_file_name = lock_payload.get("panel_file")
    curation_notes_file_name = lock_payload.get("curation_notes_file")
    if not isinstance(panel_file_name, str) or not isinstance(
        curation_notes_file_name, str
    ):
        raise RuntimeError(
            "Query-recall panel lock must declare panel_file and curation_notes_file."
        )

    panel_file_path = resolved_panel_directory / panel_file_name
    curation_notes_file_path = resolved_panel_directory / curation_notes_file_name

    expected_panel_sha256 = _require_sha256(
        lock_payload.get("panel_file_sha256"), field_name="panel_file_sha256"
    )
    actual_panel_sha256 = sha256_of_file(input_file_path=panel_file_path)
    if actual_panel_sha256 != expected_panel_sha256:
        raise RuntimeError(
            f"Query-recall panel CSV hash mismatch. Expected "
            f"{expected_panel_sha256}, got {actual_panel_sha256}."
        )

    expected_notes_sha256 = _require_sha256(
        lock_payload.get("curation_notes_sha256"), field_name="curation_notes_sha256"
    )
    actual_notes_sha256 = sha256_of_file(input_file_path=curation_notes_file_path)
    if actual_notes_sha256 != expected_notes_sha256:
        raise RuntimeError(
            f"Query-recall curation notes hash mismatch. Expected "
            f"{expected_notes_sha256}, got {actual_notes_sha256}."
        )

    declared_reference_count = lock_payload.get("reference_count")
    if declared_reference_count != EXPECTED_REFERENCE_COUNT:
        raise RuntimeError(
            f"Query-recall panel lock reference_count must be "
            f"{EXPECTED_REFERENCE_COUNT}, got {declared_reference_count!r}."
        )

    rows = _parse_panel_rows(panel_file_path)
    _validate_semantic_invariants(rows)
    _validate_declared_lock_descriptions(lock_payload, rows)

    return QueryRecallReferencePanelValidation(
        panel_directory=resolved_panel_directory,
        lock_payload=lock_payload,
        rows=rows,
        panel_sha256=actual_panel_sha256,
        curation_notes_sha256=actual_notes_sha256,
    )
