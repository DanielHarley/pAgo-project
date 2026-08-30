from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

import pandas as pd

DEFAULT_PROTEIN_UID_COLUMN = "protein_uid"
DEFAULT_SEQUENCE_COLUMN = "gbseq__sequence"
DEFAULT_LENGTH_COLUMN = "gbseq__length"

REQUIRED_TECHNICAL_PREFILTER_INPUT_COLUMNS: tuple[str, ...] = (
    DEFAULT_PROTEIN_UID_COLUMN,
    DEFAULT_SEQUENCE_COLUMN,
)

# Canonical amino acids.
DEFAULT_ALLOWED_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
# Ambiguity and non-canonical residue codes that are still biologically
# interpretable downstream (B=Asx, Z=Glx, J=Xle, X=any, U=selenocysteine,
# O=pyrrolysine). Kept, never a reason to exclude.
DEFAULT_TOLERATED_AMBIGUOUS_RESIDUES = "BZJXUO"
# Symbols that appear in NCBI protein sequences without making the record
# technically unusable (`*` terminator, `-` gap).
DEFAULT_TOLERATED_SEQUENCE_SYMBOLS = "*-"

# Length band used ONLY to raise a non-blocking warning flag. It never removes a
# record: the domain scan, not the prefilter, decides whether a length is
# compatible with a pAgo.
DEFAULT_LENGTH_WARNING_MIN = 200
DEFAULT_LENGTH_WARNING_MAX = 2000

TECHNICAL_PREFILTER_DECISION_COLUMN = "technical_prefilter_decision"
TECHNICAL_PREFILTER_REASON_COLUMN = "technical_prefilter_reason"
LENGTH_WARNING_COLUMN = "length_warning"

TECHNICAL_PREFILTER_COUNT_COLUMNS: tuple[str, ...] = (
    "decision",
    "count",
    "fraction",
    "description",
)


class PagoTechnicalPrefilterDecision(str, Enum):
    RETAIN = "retain"
    DROP_UNPROCESSABLE_RECORD = "drop_unprocessable_record"
    DROP_TECHNICAL_DUPLICATE = "drop_technical_duplicate"
    DROP_MISSING_SEQUENCE = "drop_missing_sequence"
    DROP_INVALID_SEQUENCE_CHARACTERS = "drop_invalid_sequence_characters"


DECISION_DESCRIPTIONS: dict[str, str] = {
    PagoTechnicalPrefilterDecision.RETAIN.value: (
        "Technically usable record kept for downstream analysis. No biological "
        "or annotation-based judgement is applied here."
    ),
    PagoTechnicalPrefilterDecision.DROP_UNPROCESSABLE_RECORD.value: (
        "protein_uid is missing or blank, so the record cannot be tracked "
        "through the pipeline."
    ),
    PagoTechnicalPrefilterDecision.DROP_TECHNICAL_DUPLICATE.value: (
        "A record with the same protein_uid already appeared; later copies are "
        "dropped so protein_uid stays unique."
    ),
    PagoTechnicalPrefilterDecision.DROP_MISSING_SEQUENCE.value: (
        "The amino-acid sequence column is missing or empty."
    ),
    PagoTechnicalPrefilterDecision.DROP_INVALID_SEQUENCE_CHARACTERS.value: (
        "The amino-acid sequence contains characters outside the configured "
        "residue policy after normalization."
    ),
}


@dataclass(frozen=True)
class PagoTechnicalPrefilterResult:
    """
    Partition returned after applying the technical-only prefilter to one
    flattened metadata DataFrame.
    """

    retained_records: pd.DataFrame
    excluded_records: pd.DataFrame
    counts_by_decision: dict[str, int]
    retained_protein_uids: tuple[str, ...]
    input_record_count: int


def _as_json_serializable_char_set(value: str) -> list[str]:
    return sorted(set(value))


def build_technical_prefilter_policy_payload(
    *,
    allowed_residues: str = DEFAULT_ALLOWED_RESIDUES,
    tolerated_ambiguous_residues: str = DEFAULT_TOLERATED_AMBIGUOUS_RESIDUES,
    tolerated_sequence_symbols: str = DEFAULT_TOLERATED_SEQUENCE_SYMBOLS,
    length_warning_min: int = DEFAULT_LENGTH_WARNING_MIN,
    length_warning_max: int = DEFAULT_LENGTH_WARNING_MAX,
    protein_uid_column: str = DEFAULT_PROTEIN_UID_COLUMN,
    sequence_column: str = DEFAULT_SEQUENCE_COLUMN,
    length_column: str = DEFAULT_LENGTH_COLUMN,
) -> dict[str, object]:
    return {
        "policy_kind": "pago_technical_prefilter",
        "policy_version": "1.0",
        "allowed_residues": _as_json_serializable_char_set(allowed_residues),
        "tolerated_ambiguous_residues": _as_json_serializable_char_set(
            tolerated_ambiguous_residues
        ),
        "tolerated_sequence_symbols": _as_json_serializable_char_set(
            tolerated_sequence_symbols
        ),
        "length_warning_min": int(length_warning_min),
        "length_warning_max": int(length_warning_max),
        "protein_uid_column": protein_uid_column,
        "sequence_column": sequence_column,
        "length_column": length_column,
        "exclusion_reasons": [
            PagoTechnicalPrefilterDecision.DROP_UNPROCESSABLE_RECORD.value,
            PagoTechnicalPrefilterDecision.DROP_TECHNICAL_DUPLICATE.value,
            PagoTechnicalPrefilterDecision.DROP_MISSING_SEQUENCE.value,
            PagoTechnicalPrefilterDecision.DROP_INVALID_SEQUENCE_CHARACTERS.value,
        ],
        "notes": (
            "Technical prefilter only. It never excludes a record by sequence "
            "length or by any text in the NCBI annotation. Length outside the "
            "warning band sets length_warning=True but keeps the record."
        ),
    }


def build_technical_prefilter_policy_sha256(
    *,
    allowed_residues: str = DEFAULT_ALLOWED_RESIDUES,
    tolerated_ambiguous_residues: str = DEFAULT_TOLERATED_AMBIGUOUS_RESIDUES,
    tolerated_sequence_symbols: str = DEFAULT_TOLERATED_SEQUENCE_SYMBOLS,
    length_warning_min: int = DEFAULT_LENGTH_WARNING_MIN,
    length_warning_max: int = DEFAULT_LENGTH_WARNING_MAX,
    protein_uid_column: str = DEFAULT_PROTEIN_UID_COLUMN,
    sequence_column: str = DEFAULT_SEQUENCE_COLUMN,
    length_column: str = DEFAULT_LENGTH_COLUMN,
) -> str:
    serialized_payload = json.dumps(
        build_technical_prefilter_policy_payload(
            allowed_residues=allowed_residues,
            tolerated_ambiguous_residues=tolerated_ambiguous_residues,
            tolerated_sequence_symbols=tolerated_sequence_symbols,
            length_warning_min=length_warning_min,
            length_warning_max=length_warning_max,
            protein_uid_column=protein_uid_column,
            sequence_column=sequence_column,
            length_column=length_column,
        ),
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(serialized_payload.encode("utf-8")).hexdigest()


def _validate_required_columns(*, metadata_dataframe: pd.DataFrame) -> None:
    missing_columns = [
        column_name
        for column_name in REQUIRED_TECHNICAL_PREFILTER_INPUT_COLUMNS
        if column_name not in metadata_dataframe.columns
    ]
    if missing_columns:
        raise RuntimeError(
            "metadata_dataframe is missing required technical prefilter input "
            f"columns: {missing_columns}."
        )


def _normalize_sequence_text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return "".join(str(value).split()).upper()


def build_pago_technical_prefilter_partition(
    *,
    metadata_dataframe: pd.DataFrame,
    allowed_residues: str = DEFAULT_ALLOWED_RESIDUES,
    tolerated_ambiguous_residues: str = DEFAULT_TOLERATED_AMBIGUOUS_RESIDUES,
    tolerated_sequence_symbols: str = DEFAULT_TOLERATED_SEQUENCE_SYMBOLS,
    length_warning_min: int = DEFAULT_LENGTH_WARNING_MIN,
    length_warning_max: int = DEFAULT_LENGTH_WARNING_MAX,
    protein_uid_column: str = DEFAULT_PROTEIN_UID_COLUMN,
    sequence_column: str = DEFAULT_SEQUENCE_COLUMN,
    length_column: str = DEFAULT_LENGTH_COLUMN,
    progress_callback: Callable[[str], None] | None = None,
) -> PagoTechnicalPrefilterResult:
    """
    Split one flattened metadata DataFrame into technically-usable and
    technically-unusable records. First matching rule wins.
    """
    _validate_required_columns(metadata_dataframe=metadata_dataframe)

    if length_warning_min < 0 or length_warning_max < length_warning_min:
        raise ValueError(
            "length_warning_min must be >= 0 and length_warning_max must be "
            ">= length_warning_min."
        )

    working_dataframe = metadata_dataframe.copy().reset_index(drop=True)
    input_record_count = int(len(working_dataframe))

    normalized_protein_uid = (
        working_dataframe[protein_uid_column].astype("string").str.strip()
    )
    working_dataframe[protein_uid_column] = normalized_protein_uid

    allowed_character_set = set(allowed_residues.upper())
    allowed_character_set |= set(tolerated_ambiguous_residues.upper())
    allowed_character_set |= set(tolerated_sequence_symbols.upper())

    normalized_sequence = working_dataframe[sequence_column].map(_normalize_sequence_text)

    missing_protein_uid_mask = normalized_protein_uid.isna() | (
        normalized_protein_uid.fillna("") == ""
    )
    duplicate_protein_uid_mask = normalized_protein_uid.duplicated(keep="first") & (
        ~missing_protein_uid_mask
    )
    missing_sequence_mask = normalized_sequence == ""

    def _sequence_has_invalid_characters(sequence_text: str) -> bool:
        if not sequence_text:
            return False
        return not set(sequence_text).issubset(allowed_character_set)

    invalid_sequence_mask = normalized_sequence.map(_sequence_has_invalid_characters)

    decision_series = pd.Series(
        PagoTechnicalPrefilterDecision.RETAIN.value,
        index=working_dataframe.index,
        dtype="object",
    )
    # Applied in reverse priority so earlier rules overwrite later ones.
    decision_series[invalid_sequence_mask] = (
        PagoTechnicalPrefilterDecision.DROP_INVALID_SEQUENCE_CHARACTERS.value
    )
    decision_series[missing_sequence_mask] = (
        PagoTechnicalPrefilterDecision.DROP_MISSING_SEQUENCE.value
    )
    decision_series[duplicate_protein_uid_mask] = (
        PagoTechnicalPrefilterDecision.DROP_TECHNICAL_DUPLICATE.value
    )
    decision_series[missing_protein_uid_mask] = (
        PagoTechnicalPrefilterDecision.DROP_UNPROCESSABLE_RECORD.value
    )

    working_dataframe[TECHNICAL_PREFILTER_DECISION_COLUMN] = decision_series
    working_dataframe[TECHNICAL_PREFILTER_REASON_COLUMN] = decision_series.map(
        DECISION_DESCRIPTIONS
    )

    numeric_length = pd.to_numeric(
        working_dataframe.get(length_column, pd.Series(index=working_dataframe.index)),
        errors="coerce",
    )
    length_warning_series = numeric_length.isna() | (
        numeric_length < length_warning_min
    ) | (numeric_length > length_warning_max)
    working_dataframe[LENGTH_WARNING_COLUMN] = length_warning_series.astype(bool)

    retained_mask = (
        decision_series == PagoTechnicalPrefilterDecision.RETAIN.value
    )
    retained_records = (
        working_dataframe.loc[retained_mask].copy().reset_index(drop=True)
    )
    excluded_records = (
        working_dataframe.loc[~retained_mask].copy().reset_index(drop=True)
    )

    if len(retained_records) + len(excluded_records) != input_record_count:
        raise RuntimeError(
            "Technical prefilter must assign each record to exactly one bucket. "
            f"Expected {input_record_count}, got "
            f"{len(retained_records) + len(excluded_records)}."
        )

    counts_by_decision = {
        decision.value: int((decision_series == decision.value).sum())
        for decision in PagoTechnicalPrefilterDecision
    }

    if progress_callback is not None:
        progress_callback(
            f"technical prefilter: {counts_by_decision[PagoTechnicalPrefilterDecision.RETAIN.value]} "
            f"retained of {input_record_count}"
        )

    return PagoTechnicalPrefilterResult(
        retained_records=retained_records,
        excluded_records=excluded_records,
        counts_by_decision=counts_by_decision,
        retained_protein_uids=tuple(
            retained_records[protein_uid_column].astype(str).tolist()
        ),
        input_record_count=input_record_count,
    )


def build_technical_prefilter_counts_dataframe(
    *,
    result: PagoTechnicalPrefilterResult,
) -> pd.DataFrame:
    total_count = result.input_record_count
    rows: list[dict[str, object]] = []
    for decision in PagoTechnicalPrefilterDecision:
        decision_count = result.counts_by_decision.get(decision.value, 0)
        rows.append(
            {
                "decision": decision.value,
                "count": int(decision_count),
                "fraction": (
                    float(decision_count / total_count) if total_count else 0.0
                ),
                "description": DECISION_DESCRIPTIONS[decision.value],
            }
        )
    return pd.DataFrame(rows, columns=list(TECHNICAL_PREFILTER_COUNT_COLUMNS))
