from __future__ import annotations

import csv
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, TypeAlias

from src.pago_pipeline.storage import sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

DEFAULT_PROTEIN_UID_COLUMN = "protein_uid"
DEFAULT_ACCESSION_COLUMN = "gbseq__accession_version"
DEFAULT_LENGTH_COLUMN = "gbseq__length"
DEFAULT_ORGANISM_COLUMN = "gbseq__organism"
DEFAULT_DESCRIPTION_COLUMN = "gbseq__definition"
DEFAULT_SEQUENCE_COLUMN = "gbseq__sequence"
DEFAULT_REQUIRED_COLUMNS = (
    DEFAULT_PROTEIN_UID_COLUMN,
    DEFAULT_ACCESSION_COLUMN,
    DEFAULT_LENGTH_COLUMN,
    DEFAULT_ORGANISM_COLUMN,
    DEFAULT_DESCRIPTION_COLUMN,
    DEFAULT_SEQUENCE_COLUMN,
)


@dataclass(frozen=True)
class MetadataFastaExportResult:
    """
    Summary returned after converting one flattened metadata CSV into FASTA.
    """

    metadata_csv_file_path: Path
    fasta_file_path: Path
    manifest_file_path: Optional[Path]
    row_count: int
    fasta_record_count: int
    skipped_missing_sequence_count: int
    sequence_line_width: int
    required_columns: list[str]
    metadata_csv_file_sha256: str
    fasta_file_sha256: str


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _normalize_header_value(
    value: str | None,
    *,
    max_len: int | None = None,
) -> str:
    """
    Normalize FASTA header key-value fields into one parser-friendly token.
    """

    normalized_value = " ".join(str(value or "").split())
    normalized_value = normalized_value.replace("|", "/").replace(";", ",")
    normalized_value = normalized_value.replace(" ", "_")

    if max_len is not None and len(normalized_value) > max_len:
        if max_len <= 3:
            return normalized_value[:max_len]
        return normalized_value[: max_len - 3] + "..."

    return normalized_value


def _normalize_description_text(value: str | None) -> str:
    """
    Normalize the human-readable description placed after the FASTA identifier.
    """

    return " ".join(str(value or "").split()).replace("|", "/")


def _normalize_sequence_text(sequence_text: str | None) -> str:
    return "".join(str(sequence_text or "").split())


def _wrap_sequence_lines(
    sequence_text: str,
    *,
    width: int,
) -> list[str]:
    if width <= 0:
        raise ValueError("sequence_line_width must be a positive integer.")

    return [
        sequence_text[start_index : start_index + width]
        for start_index in range(0, len(sequence_text), width)
    ]


def _validate_required_columns(
    *,
    fieldnames: list[str] | None,
    required_columns: tuple[str, ...],
) -> None:
    if fieldnames is None:
        raise RuntimeError("Metadata CSV is missing a header row.")

    missing_columns = sorted(
        required_column
        for required_column in required_columns
        if required_column not in fieldnames
    )
    if missing_columns:
        raise RuntimeError(
            "Metadata CSV is missing required columns for FASTA export: "
            f"{missing_columns}."
        )


def _build_fasta_defline(
    *,
    row: Mapping[str, str],
    protein_uid_column: str,
    accession_column: str,
    length_column: str,
    organism_column: str,
    description_column: str,
) -> str:
    protein_uid = _normalize_header_value(row.get(protein_uid_column))
    accession = _normalize_header_value(row.get(accession_column))
    sequence_length = _normalize_header_value(row.get(length_column))
    organism = _normalize_header_value(row.get(organism_column), max_len=120)
    description = _normalize_description_text(row.get(description_column))

    header_core_parts = [
        f"protein_uid={protein_uid}" if protein_uid else None,
        f"accession={accession}" if accession else None,
        f"length={sequence_length}" if sequence_length else None,
        f"organism={organism}" if organism else None,
    ]
    header_core = "|".join(
        header_core_part
        for header_core_part in header_core_parts
        if header_core_part is not None
    )
    if not header_core:
        header_core = "record"

    if description:
        return f">{header_core} {description}"

    return f">{header_core}"


def _build_manifest_payload(
    *,
    metadata_csv_file_path: Path,
    fasta_file_path: Path,
    metadata_csv_file_sha256: str,
    fasta_file_sha256: str,
    row_count: int,
    fasta_record_count: int,
    skipped_missing_sequence_count: int,
    sequence_line_width: int,
    required_columns: tuple[str, ...],
    protein_uid_column: str,
    accession_column: str,
    length_column: str,
    organism_column: str,
    description_column: str,
    sequence_column: str,
    source_metadata_manifest_payload: Optional[Mapping[str, object]],
    source_metadata_manifest_file_path: Optional[Path],
) -> dict[str, object]:
    manifest_payload: dict[str, object] = {
        "artifact_type": "ncbi_protein_fasta_export",
        "snapshot_format_version": "1.0",
        "metadata_csv_file_name": metadata_csv_file_path.name,
        "metadata_csv_file_path": str(metadata_csv_file_path),
        "metadata_csv_file_sha256": metadata_csv_file_sha256,
        "fasta_file_name": fasta_file_path.name,
        "fasta_file_path": str(fasta_file_path),
        "fasta_file_sha256": fasta_file_sha256,
        "row_count": row_count,
        "fasta_record_count": fasta_record_count,
        "skipped_missing_sequence_count": skipped_missing_sequence_count,
        "sequence_line_width": sequence_line_width,
        "required_columns": list(required_columns),
        "protein_uid_column": protein_uid_column,
        "accession_column": accession_column,
        "length_column": length_column,
        "organism_column": organism_column,
        "description_column": description_column,
        "sequence_column": sequence_column,
        "defline_primary_fields": [
            "protein_uid",
            "accession",
            "length",
            "organism",
        ],
    }

    if source_metadata_manifest_payload is not None:
        manifest_payload["source_metadata_snapshot_relative_path"] = (
            source_metadata_manifest_payload.get("immutable_snapshot_relative_path")
        )
        manifest_payload["source_metadata_snapshot_directory_name"] = (
            source_metadata_manifest_payload.get("immutable_snapshot_directory_name")
        )
        manifest_payload["source_metadata_row_count"] = source_metadata_manifest_payload.get(
            "row_count"
        )
        manifest_payload["source_metadata_csv_sha256"] = source_metadata_manifest_payload.get(
            "csv_file_sha256"
        )
        manifest_payload["source_xml_snapshot_relative_path"] = (
            source_metadata_manifest_payload.get("source_xml_snapshot_relative_path")
        )
        manifest_payload["source_xml_snapshot_directory_name"] = (
            source_metadata_manifest_payload.get("source_xml_snapshot_directory_name")
        )
        manifest_payload["source_xml_file_sha256"] = source_metadata_manifest_payload.get(
            "source_xml_file_sha256"
        )
        manifest_payload["search_query"] = source_metadata_manifest_payload.get(
            "search_query"
        )
        manifest_payload["translated_query"] = source_metadata_manifest_payload.get(
            "translated_query"
        )

    if source_metadata_manifest_file_path is not None:
        manifest_payload["source_metadata_manifest_file_path"] = str(
            source_metadata_manifest_file_path
        )
        manifest_payload["source_metadata_manifest_file_sha256"] = sha256_of_file(
            input_file_path=source_metadata_manifest_file_path
        )

    return manifest_payload


def export_metadata_csv_to_fasta(
    *,
    metadata_csv_file_path: PathLike,
    output_fasta_file_path: PathLike,
    output_manifest_file_path: Optional[PathLike] = None,
    sequence_line_width: int = 60,
    protein_uid_column: str = DEFAULT_PROTEIN_UID_COLUMN,
    accession_column: str = DEFAULT_ACCESSION_COLUMN,
    length_column: str = DEFAULT_LENGTH_COLUMN,
    organism_column: str = DEFAULT_ORGANISM_COLUMN,
    description_column: str = DEFAULT_DESCRIPTION_COLUMN,
    sequence_column: str = DEFAULT_SEQUENCE_COLUMN,
    source_metadata_manifest_payload: Optional[Mapping[str, object]] = None,
    source_metadata_manifest_file_path: Optional[PathLike] = None,
) -> MetadataFastaExportResult:
    """
    Convert one flattened protein metadata CSV into a multi-FASTA file.
    """

    resolved_metadata_csv_file_path = _as_path(metadata_csv_file_path)
    resolved_output_fasta_file_path = _as_path(output_fasta_file_path)
    resolved_output_fasta_file_path.parent.mkdir(parents=True, exist_ok=True)

    required_columns = (
        protein_uid_column,
        accession_column,
        length_column,
        organism_column,
        description_column,
        sequence_column,
    )

    with resolved_metadata_csv_file_path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as metadata_file_handle:
        metadata_csv_reader = csv.DictReader(metadata_file_handle)
        _validate_required_columns(
            fieldnames=metadata_csv_reader.fieldnames,
            required_columns=required_columns,
        )

        row_count = 0
        fasta_record_count = 0
        skipped_missing_sequence_count = 0

        with tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            dir=resolved_output_fasta_file_path.parent,
            encoding="utf-8",
            newline="\n",
        ) as temporary_fasta_file:
            resolved_temporary_fasta_file_path = Path(temporary_fasta_file.name)

            for row in metadata_csv_reader:
                row_count += 1
                normalized_sequence = _normalize_sequence_text(row.get(sequence_column))
                if not normalized_sequence:
                    skipped_missing_sequence_count += 1
                    continue

                fasta_defline = _build_fasta_defline(
                    row=row,
                    protein_uid_column=protein_uid_column,
                    accession_column=accession_column,
                    length_column=length_column,
                    organism_column=organism_column,
                    description_column=description_column,
                )
                wrapped_sequence_lines = _wrap_sequence_lines(
                    normalized_sequence,
                    width=sequence_line_width,
                )

                temporary_fasta_file.write(f"{fasta_defline}\n")
                temporary_fasta_file.write("\n".join(wrapped_sequence_lines))
                temporary_fasta_file.write("\n")
                fasta_record_count += 1

    resolved_temporary_fasta_file_path.replace(resolved_output_fasta_file_path)

    resolved_output_manifest_file_path: Optional[Path] = None
    resolved_source_metadata_manifest_file_path: Optional[Path] = None
    if source_metadata_manifest_file_path is not None:
        resolved_source_metadata_manifest_file_path = _as_path(
            source_metadata_manifest_file_path
        )

    metadata_csv_file_sha256 = sha256_of_file(
        input_file_path=resolved_metadata_csv_file_path
    )
    fasta_file_sha256 = sha256_of_file(input_file_path=resolved_output_fasta_file_path)

    if output_manifest_file_path is not None:
        resolved_output_manifest_file_path = _as_path(output_manifest_file_path)
        manifest_payload = _build_manifest_payload(
            metadata_csv_file_path=resolved_metadata_csv_file_path,
            fasta_file_path=resolved_output_fasta_file_path,
            metadata_csv_file_sha256=metadata_csv_file_sha256,
            fasta_file_sha256=fasta_file_sha256,
            row_count=row_count,
            fasta_record_count=fasta_record_count,
            skipped_missing_sequence_count=skipped_missing_sequence_count,
            sequence_line_width=sequence_line_width,
            required_columns=required_columns,
            protein_uid_column=protein_uid_column,
            accession_column=accession_column,
            length_column=length_column,
            organism_column=organism_column,
            description_column=description_column,
            sequence_column=sequence_column,
            source_metadata_manifest_payload=source_metadata_manifest_payload,
            source_metadata_manifest_file_path=resolved_source_metadata_manifest_file_path,
        )
        write_json_atomic(
            payload=manifest_payload,
            output_file_path=resolved_output_manifest_file_path,
        )

    return MetadataFastaExportResult(
        metadata_csv_file_path=resolved_metadata_csv_file_path,
        fasta_file_path=resolved_output_fasta_file_path,
        manifest_file_path=resolved_output_manifest_file_path,
        row_count=row_count,
        fasta_record_count=fasta_record_count,
        skipped_missing_sequence_count=skipped_missing_sequence_count,
        sequence_line_width=sequence_line_width,
        required_columns=list(required_columns),
        metadata_csv_file_sha256=metadata_csv_file_sha256,
        fasta_file_sha256=fasta_file_sha256,
    )
