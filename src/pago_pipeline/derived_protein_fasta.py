from __future__ import annotations

import hashlib
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import Enum

import pandas as pd

DEFAULT_PROTEIN_UID_COLUMN = "protein_uid"
DEFAULT_SEQUENCE_COLUMN = "gbseq__sequence"

REQUIRED_DERIVED_FASTA_INPUT_COLUMNS: tuple[str, ...] = (
    DEFAULT_PROTEIN_UID_COLUMN,
    DEFAULT_SEQUENCE_COLUMN,
)


class DerivedFastaRecordOrder(str, Enum):
    AS_SELECTED = "as_selected"
    SORTED_BY_UID = "sorted_by_uid"


@dataclass(frozen=True)
class DerivedFastaSelectionResult:
    """
    Result of resolving a protein_uid selection against a flattened metadata
    DataFrame, ready to be written as a provenance-tracked derived FASTA.
    """

    selected_metadata: pd.DataFrame
    resolved_protein_uids: tuple[str, ...]
    requested_uid_count: int
    resolved_uid_count: int
    missing_uids: tuple[str, ...]
    duplicate_requested_uids: tuple[str, ...]
    empty_sequence_uids: tuple[str, ...]
    record_order: str
    source_record_ids_sha256: str


def _normalize_uid(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def _compute_record_ids_sha256(resolved_protein_uids: Sequence[str]) -> str:
    serialized = "\n".join(resolved_protein_uids)
    if serialized:
        serialized += "\n"
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def build_derived_fasta_selection(
    *,
    metadata_dataframe: pd.DataFrame,
    selected_protein_uids: Sequence[str],
    protein_uid_column: str = DEFAULT_PROTEIN_UID_COLUMN,
    sequence_column: str = DEFAULT_SEQUENCE_COLUMN,
    record_order: str = DerivedFastaRecordOrder.AS_SELECTED.value,
    drop_missing_uids: bool = False,
    progress_callback: Callable[[str], None] | None = None,
) -> DerivedFastaSelectionResult:
    """
    Resolve an ordered protein_uid selection against one flattened metadata
    DataFrame. Raises on a selection that cannot be materialized faithfully
    unless drop_missing_uids is set.
    """
    missing_columns = [
        column_name
        for column_name in (protein_uid_column, sequence_column)
        if column_name not in metadata_dataframe.columns
    ]
    if missing_columns:
        raise RuntimeError(
            "metadata_dataframe is missing required derived FASTA input columns: "
            f"{missing_columns}."
        )

    try:
        resolved_record_order = DerivedFastaRecordOrder(record_order)
    except ValueError as error:
        raise ValueError(
            "record_order must be one of "
            f"{[member.value for member in DerivedFastaRecordOrder]}, got "
            f"{record_order!r}."
        ) from error

    metadata_dataframe_copy = metadata_dataframe.copy()
    metadata_dataframe_copy[protein_uid_column] = metadata_dataframe_copy[
        protein_uid_column
    ].map(_normalize_uid)

    metadata_duplicate_count = int(
        metadata_dataframe_copy[protein_uid_column]
        .loc[metadata_dataframe_copy[protein_uid_column] != ""]
        .duplicated()
        .sum()
    )
    if metadata_duplicate_count:
        raise RuntimeError(
            "metadata_dataframe protein_uid values must be unique for a derived "
            f"FASTA selection. Duplicate count: {metadata_duplicate_count}."
        )

    normalized_requested = [_normalize_uid(value) for value in selected_protein_uids]
    normalized_requested = [value for value in normalized_requested if value]
    requested_uid_count = len(normalized_requested)

    seen: set[str] = set()
    duplicate_requested: list[str] = []
    ordered_unique_requested: list[str] = []
    for uid in normalized_requested:
        if uid in seen:
            duplicate_requested.append(uid)
            continue
        seen.add(uid)
        ordered_unique_requested.append(uid)

    if duplicate_requested and not drop_missing_uids:
        raise RuntimeError(
            "selected_protein_uids contains duplicates: "
            f"{sorted(set(duplicate_requested))}."
        )

    metadata_by_uid = metadata_dataframe_copy.set_index(protein_uid_column, drop=False)
    available_uids = set(metadata_by_uid.index)

    missing_uids = [
        uid for uid in ordered_unique_requested if uid not in available_uids
    ]
    if missing_uids and not drop_missing_uids:
        raise RuntimeError(
            f"selected_protein_uids has {len(missing_uids)} uid(s) absent from the "
            f"metadata: {missing_uids[:10]}."
        )

    resolved_uids = [
        uid for uid in ordered_unique_requested if uid in available_uids
    ]
    if resolved_record_order == DerivedFastaRecordOrder.SORTED_BY_UID:
        resolved_uids = sorted(resolved_uids)

    selected_metadata = metadata_by_uid.loc[resolved_uids].reset_index(drop=True)

    empty_sequence_mask = (
        selected_metadata[sequence_column]
        .map(lambda value: "".join(str(value or "").split()))
        .eq("")
    )
    empty_sequence_uids = tuple(
        selected_metadata.loc[empty_sequence_mask, protein_uid_column].astype(str)
    )
    if empty_sequence_uids:
        raise RuntimeError(
            "Derived FASTA selection contains records with an empty amino-acid "
            f"sequence: {list(empty_sequence_uids)[:10]}."
        )

    if progress_callback is not None:
        progress_callback(
            f"derived FASTA selection: {len(resolved_uids)} of "
            f"{requested_uid_count} requested uid(s) resolved"
        )

    return DerivedFastaSelectionResult(
        selected_metadata=selected_metadata,
        resolved_protein_uids=tuple(resolved_uids),
        requested_uid_count=requested_uid_count,
        resolved_uid_count=len(resolved_uids),
        missing_uids=tuple(missing_uids),
        duplicate_requested_uids=tuple(dict.fromkeys(duplicate_requested)),
        empty_sequence_uids=empty_sequence_uids,
        record_order=resolved_record_order.value,
        source_record_ids_sha256=_compute_record_ids_sha256(resolved_uids),
    )


def parse_protein_uids_from_fasta_deflines(
    *,
    fasta_file_path: str,
) -> list[str]:
    """
    Extract the protein_uid= token from each defline, in file order.
    """
    protein_uids: list[str] = []
    with open(fasta_file_path, "r", encoding="utf-8") as fasta_file_handle:
        for line in fasta_file_handle:
            if not line.startswith(">"):
                continue
            header_core = line[1:].strip().split(" ", 1)[0]
            for field in header_core.split("|"):
                if field.startswith("protein_uid="):
                    protein_uids.append(field.split("=", 1)[1])
                    break
    return protein_uids


def compute_record_ids_sha256_for_order(
    *,
    protein_uids: Sequence[str],
    record_order: str,
) -> str:
    resolved_record_order = DerivedFastaRecordOrder(record_order)
    ordered = list(protein_uids)
    if resolved_record_order == DerivedFastaRecordOrder.SORTED_BY_UID:
        ordered = sorted(ordered)
    return _compute_record_ids_sha256(ordered)
