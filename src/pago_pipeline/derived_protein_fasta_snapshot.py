from __future__ import annotations

import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

import pandas as pd

from src.pago_pipeline.derived_protein_fasta import (
    DEFAULT_PROTEIN_UID_COLUMN,
    DerivedFastaRecordOrder,
    build_derived_fasta_selection,
    compute_record_ids_sha256_for_order,
    parse_protein_uids_from_fasta_deflines,
)
from src.pago_pipeline.metadata_to_fasta import export_metadata_csv_to_fasta
from src.pago_pipeline.ncbi_metadata_snapshot import load_latest_metadata_snapshot
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.storage import (
    read_json_file,
    read_text_lines_from_file,
    sha256_of_file,
    write_json_atomic,
)

PathLike: TypeAlias = str | Path

ARTIFACT_TYPE = "derived_protein_fasta_snapshot"
SNAPSHOT_FORMAT_VERSION = "1.0"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
DEFAULT_FASTA_FILE_NAME = "protein_sequences.fasta"
DEFAULT_SELECTION_REPORT_FILE_NAME = "selection_report.json"


@dataclass(frozen=True)
class DerivedProteinFastaSnapshotResult:
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    fasta_file_path: Path
    selection_report_file_path: Path
    dataset_kind: str
    requested_uid_count: int
    resolved_uid_count: int
    fasta_record_count: int


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _load_selection_snapshot(
    *,
    snapshot_root_directory: Path,
    expected_selection_artifact_type: str,
    selection_uid_list_file_name: str,
) -> dict[str, object]:
    latest_directory = snapshot_root_directory / "latest"
    manifest_file_path = latest_directory / DEFAULT_MANIFEST_FILE_NAME
    uid_list_file_path = latest_directory / selection_uid_list_file_name
    if not manifest_file_path.exists():
        raise FileNotFoundError(
            f"Selection snapshot manifest not found: {manifest_file_path}."
        )
    if not uid_list_file_path.exists():
        raise FileNotFoundError(
            f"Selection snapshot uid list file not found: {uid_list_file_path}."
        )
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Selection snapshot manifest must deserialize into a JSON object."
        )
    selection_artifact_type = manifest_payload.get("artifact_type")
    if selection_artifact_type != expected_selection_artifact_type:
        raise RuntimeError(
            "Selection snapshot artifact_type mismatch. Expected "
            f"{expected_selection_artifact_type!r}, got {selection_artifact_type!r}."
        )
    uid_list = read_text_lines_from_file(input_file_path=uid_list_file_path)
    return {
        "manifest": dict(manifest_payload),
        "manifest_file_path": manifest_file_path,
        "uid_list_file_path": uid_list_file_path,
        "uid_list": uid_list,
    }


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    dataset_kind: str,
    record_selection_rule: str,
    record_selection_config_sha256: str,
    record_order: str,
    source_record_ids_sha256: str,
    requested_uid_count: int,
    resolved_uid_count: int,
    missing_uid_count: int,
    fasta_file_path: Path,
    fasta_record_count: int,
    skipped_missing_sequence_count: int,
    sequence_line_width: int,
    selection_report_file_path: Path,
    source_metadata_manifest: Mapping[str, object],
    source_metadata_manifest_file_path: Path,
    source_metadata_csv_file_path: Path,
    selection_manifest: Mapping[str, object],
    selection_manifest_file_path: Path,
    selection_artifact_type: str,
) -> dict[str, object]:
    return {
        "artifact_type": ARTIFACT_TYPE,
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "dataset_kind": dataset_kind,
        # Provenance of the transformation that produced this FASTA.
        "derived_from_artifact_type": selection_artifact_type,
        "derived_from_manifest_sha256": sha256_of_file(
            input_file_path=selection_manifest_file_path
        ),
        "record_selection_rule": record_selection_rule,
        "record_selection_config_sha256": record_selection_config_sha256,
        "record_order": record_order,
        "source_record_ids_sha256": source_record_ids_sha256,
        "requested_uid_count": int(requested_uid_count),
        "resolved_uid_count": int(resolved_uid_count),
        "missing_uid_count": int(missing_uid_count),
        # FASTA artifact facts.
        "fasta_file_name": fasta_file_path.name,
        "fasta_file_sha256": sha256_of_file(input_file_path=fasta_file_path),
        "fasta_record_count": int(fasta_record_count),
        "row_count": int(resolved_uid_count),
        "skipped_missing_sequence_count": int(skipped_missing_sequence_count),
        "sequence_line_width": int(sequence_line_width),
        "selection_report_file_name": selection_report_file_path.name,
        "selection_report_file_sha256": sha256_of_file(
            input_file_path=selection_report_file_path
        ),
        # Source metadata snapshot provenance (also used by downstream SWeeP).
        "source_metadata_snapshot_relative_path": source_metadata_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_metadata_snapshot_directory_name": source_metadata_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_metadata_manifest_sha256": sha256_of_file(
            input_file_path=source_metadata_manifest_file_path
        ),
        "source_metadata_csv_sha256": sha256_of_file(
            input_file_path=source_metadata_csv_file_path
        ),
        "source_xml_snapshot_relative_path": source_metadata_manifest.get(
            "source_xml_snapshot_relative_path"
        ),
        "search_query": source_metadata_manifest.get("search_query"),
        "translated_query": source_metadata_manifest.get("translated_query"),
        # Selection snapshot provenance.
        "source_selection_artifact_type": selection_artifact_type,
        "source_selection_snapshot_relative_path": selection_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_selection_snapshot_directory_name": selection_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_selection_manifest_sha256": sha256_of_file(
            input_file_path=selection_manifest_file_path
        ),
    }


def save_derived_protein_fasta_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    source_selection_snapshot_root_directory: PathLike,
    selection_artifact_type: str,
    selection_uid_list_file_name: str,
    record_selection_rule: str,
    record_selection_config_sha256: str,
    dataset_kind: str,
    sequence_line_width: int = 60,
    record_order: str = DerivedFastaRecordOrder.AS_SELECTED.value,
    drop_missing_uids: bool = False,
    update_latest_directory: bool = True,
) -> DerivedProteinFastaSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )
    resolved_source_selection_snapshot_root_directory = _as_path(
        source_selection_snapshot_root_directory
    )

    try:
        source_metadata_snapshot_payload = load_latest_metadata_snapshot(
            snapshot_root_directory=resolved_source_metadata_snapshot_root_directory
        )
    except FileNotFoundError as error:
        raise FileNotFoundError(
            "No reusable metadata snapshot was found for the derived FASTA "
            "stage. Run the upstream metadata extraction stage first."
        ) from error

    source_metadata_manifest = source_metadata_snapshot_payload.get("manifest")
    source_metadata_manifest_file_path = source_metadata_snapshot_payload.get(
        "manifest_file_path"
    )
    source_metadata_csv_file_path = source_metadata_snapshot_payload.get("csv_file_path")
    if (
        not isinstance(source_metadata_manifest, Mapping)
        or not isinstance(source_metadata_manifest_file_path, (str, Path))
        or not isinstance(source_metadata_csv_file_path, (str, Path))
    ):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing manifest / csv path."
        )
    source_metadata_manifest_file_path = _as_path(source_metadata_manifest_file_path)
    source_metadata_csv_file_path = _as_path(source_metadata_csv_file_path)

    selection_snapshot = _load_selection_snapshot(
        snapshot_root_directory=resolved_source_selection_snapshot_root_directory,
        expected_selection_artifact_type=selection_artifact_type,
        selection_uid_list_file_name=selection_uid_list_file_name,
    )

    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query=f"derived_protein_fasta::{dataset_kind}",
    )
    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    immutable_snapshot_relative_path = str(
        Path("snapshots") / snapshot_directory_name
    )
    immutable_snapshot_complete = False

    try:
        metadata_dataframe = pd.read_csv(
            source_metadata_csv_file_path,
            low_memory=False,
            dtype={DEFAULT_PROTEIN_UID_COLUMN: "string"},
        )
        selection_result = build_derived_fasta_selection(
            metadata_dataframe=metadata_dataframe,
            selected_protein_uids=list(selection_snapshot["uid_list"]),
            record_order=record_order,
            drop_missing_uids=drop_missing_uids,
        )

        fasta_file_path = immutable_snapshot_directory / DEFAULT_FASTA_FILE_NAME
        with tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            dir=immutable_snapshot_directory,
            encoding="utf-8",
            newline="",
            suffix=".csv",
        ) as selected_csv_file:
            selected_csv_path = Path(selected_csv_file.name)
        selection_result.selected_metadata.to_csv(selected_csv_path, index=False)
        try:
            fasta_export_result = export_metadata_csv_to_fasta(
                metadata_csv_file_path=selected_csv_path,
                output_fasta_file_path=fasta_file_path,
                sequence_line_width=sequence_line_width,
            )
        finally:
            selected_csv_path.unlink(missing_ok=True)

        selection_report_file_path = (
            immutable_snapshot_directory / DEFAULT_SELECTION_REPORT_FILE_NAME
        )
        write_json_atomic(
            payload={
                "dataset_kind": dataset_kind,
                "record_selection_rule": record_selection_rule,
                "record_order": selection_result.record_order,
                "requested_uid_count": selection_result.requested_uid_count,
                "resolved_uid_count": selection_result.resolved_uid_count,
                "missing_uids": list(selection_result.missing_uids),
                "duplicate_requested_uids": list(
                    selection_result.duplicate_requested_uids
                ),
                "empty_sequence_uids": list(selection_result.empty_sequence_uids),
                "source_record_ids_sha256": selection_result.source_record_ids_sha256,
                "fasta_record_count": fasta_export_result.fasta_record_count,
                "skipped_missing_sequence_count": (
                    fasta_export_result.skipped_missing_sequence_count
                ),
            },
            output_file_path=selection_report_file_path,
        )

        manifest_file_path = (
            immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
        )
        manifest_payload = _build_manifest_payload(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            dataset_kind=dataset_kind,
            record_selection_rule=record_selection_rule,
            record_selection_config_sha256=record_selection_config_sha256,
            record_order=selection_result.record_order,
            source_record_ids_sha256=selection_result.source_record_ids_sha256,
            requested_uid_count=selection_result.requested_uid_count,
            resolved_uid_count=selection_result.resolved_uid_count,
            missing_uid_count=len(selection_result.missing_uids),
            fasta_file_path=fasta_file_path,
            fasta_record_count=fasta_export_result.fasta_record_count,
            skipped_missing_sequence_count=(
                fasta_export_result.skipped_missing_sequence_count
            ),
            sequence_line_width=sequence_line_width,
            selection_report_file_path=selection_report_file_path,
            source_metadata_manifest=source_metadata_manifest,
            source_metadata_manifest_file_path=source_metadata_manifest_file_path,
            source_metadata_csv_file_path=source_metadata_csv_file_path,
            selection_manifest=selection_snapshot["manifest"],
            selection_manifest_file_path=selection_snapshot["manifest_file_path"],
            selection_artifact_type=selection_artifact_type,
        )
        write_json_atomic(
            payload=manifest_payload,
            output_file_path=manifest_file_path,
        )
        immutable_snapshot_complete = True

        if update_latest_directory:
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=[
                    (fasta_file_path, DEFAULT_FASTA_FILE_NAME),
                    (
                        selection_report_file_path,
                        DEFAULT_SELECTION_REPORT_FILE_NAME,
                    ),
                    (manifest_file_path, DEFAULT_MANIFEST_FILE_NAME),
                ],
            )
    except Exception:
        if (
            not immutable_snapshot_complete
            and immutable_snapshot_directory.exists()
        ):
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return DerivedProteinFastaSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        fasta_file_path=fasta_file_path,
        selection_report_file_path=selection_report_file_path,
        dataset_kind=dataset_kind,
        requested_uid_count=selection_result.requested_uid_count,
        resolved_uid_count=selection_result.resolved_uid_count,
        fasta_record_count=fasta_export_result.fasta_record_count,
    )


def _validate_loaded_derived_protein_fasta_payload(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
) -> dict[str, Path]:
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != ARTIFACT_TYPE:
        raise RuntimeError(
            "Saved derived protein FASTA snapshot manifest artifact_type "
            f"mismatch. Expected {ARTIFACT_TYPE!r}, got {artifact_type!r}."
        )

    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved derived protein FASTA snapshot snapshot_format_version "
            f"mismatch. Expected {SNAPSHOT_FORMAT_VERSION!r}, got "
            f"{snapshot_format_version!r}."
        )

    fasta_file_name = manifest_payload.get("fasta_file_name")
    if not isinstance(fasta_file_name, str) or not fasta_file_name.strip():
        raise RuntimeError(
            "Saved derived protein FASTA snapshot must define fasta_file_name."
        )
    fasta_file_path = snapshot_directory / fasta_file_name
    if not fasta_file_path.exists():
        raise FileNotFoundError(
            f"Saved derived protein FASTA file not found: {fasta_file_path}."
        )

    expected_fasta_sha256 = manifest_payload.get("fasta_file_sha256")
    if isinstance(expected_fasta_sha256, str):
        actual_fasta_sha256 = sha256_of_file(input_file_path=fasta_file_path)
        if actual_fasta_sha256 != expected_fasta_sha256:
            raise RuntimeError(
                "Saved derived protein FASTA file hash mismatch. "
                f"Expected {expected_fasta_sha256}, got {actual_fasta_sha256}."
            )

    fasta_protein_uids = parse_protein_uids_from_fasta_deflines(
        fasta_file_path=str(fasta_file_path)
    )

    expected_record_count = manifest_payload.get("fasta_record_count")
    if isinstance(expected_record_count, int) and len(
        fasta_protein_uids
    ) != expected_record_count:
        raise RuntimeError(
            "Saved derived protein FASTA record count mismatch. Expected "
            f"{expected_record_count}, found {len(fasta_protein_uids)} deflines."
        )

    expected_record_ids_sha256 = manifest_payload.get("source_record_ids_sha256")
    record_order = manifest_payload.get("record_order")
    if isinstance(expected_record_ids_sha256, str) and isinstance(record_order, str):
        actual_record_ids_sha256 = compute_record_ids_sha256_for_order(
            protein_uids=fasta_protein_uids,
            record_order=record_order,
        )
        if actual_record_ids_sha256 != expected_record_ids_sha256:
            raise RuntimeError(
                "Saved derived protein FASTA record-id set does not match "
                "source_record_ids_sha256; the FASTA content or order changed."
            )

    selection_report_file_name = manifest_payload.get("selection_report_file_name")
    selection_report_file_path: Optional[Path] = None
    if isinstance(selection_report_file_name, str) and selection_report_file_name:
        selection_report_file_path = snapshot_directory / selection_report_file_name
        if not selection_report_file_path.exists():
            raise FileNotFoundError(
                "Saved derived protein FASTA selection report not found: "
                f"{selection_report_file_path}."
            )
        expected_report_sha256 = manifest_payload.get("selection_report_file_sha256")
        if isinstance(expected_report_sha256, str):
            actual_report_sha256 = sha256_of_file(
                input_file_path=selection_report_file_path
            )
            if actual_report_sha256 != expected_report_sha256:
                raise RuntimeError(
                    "Saved derived protein FASTA selection report hash mismatch."
                )

    return {
        "fasta_file": fasta_file_path,
        **(
            {"selection_report_file": selection_report_file_path}
            if selection_report_file_path is not None
            else {}
        ),
    }


def load_derived_protein_fasta_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved derived protein FASTA snapshot manifest must deserialize into "
            "a JSON object."
        )
    resolved_file_path_by_key = _validate_loaded_derived_protein_fasta_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )
    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "fasta_file_path": resolved_file_path_by_key["fasta_file"],
    }


def load_latest_derived_protein_fasta_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_derived_protein_fasta_snapshot_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest"
    )


def _source_identity_matches(
    *,
    manifest_payload: Mapping[str, object],
    source_metadata_snapshot_root_directory: Optional[Path],
    source_selection_snapshot_root_directory: Optional[Path],
    selection_artifact_type: str,
    selection_uid_list_file_name: str,
) -> bool:
    if source_metadata_snapshot_root_directory is not None:
        source_metadata_snapshot_payload = load_latest_metadata_snapshot(
            snapshot_root_directory=source_metadata_snapshot_root_directory
        )
        source_metadata_manifest_file_path = source_metadata_snapshot_payload.get(
            "manifest_file_path"
        )
        if not isinstance(source_metadata_manifest_file_path, (str, Path)):
            return False
        if manifest_payload.get("source_metadata_manifest_sha256") != sha256_of_file(
            input_file_path=_as_path(source_metadata_manifest_file_path)
        ):
            return False

    if source_selection_snapshot_root_directory is not None:
        selection_snapshot = _load_selection_snapshot(
            snapshot_root_directory=source_selection_snapshot_root_directory,
            expected_selection_artifact_type=selection_artifact_type,
            selection_uid_list_file_name=selection_uid_list_file_name,
        )
        if manifest_payload.get("source_selection_manifest_sha256") != sha256_of_file(
            input_file_path=_as_path(selection_snapshot["manifest_file_path"])
        ):
            return False

    return True


def latest_derived_protein_fasta_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: Optional[PathLike] = None,
    source_selection_snapshot_root_directory: Optional[PathLike] = None,
    selection_artifact_type: Optional[str] = None,
    selection_uid_list_file_name: Optional[str] = None,
) -> bool:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / DEFAULT_MANIFEST_FILE_NAME
    if not latest_directory.exists() or not latest_manifest_file_path.exists():
        return False
    try:
        manifest_payload = read_json_file(input_file_path=latest_manifest_file_path)
        if not isinstance(manifest_payload, Mapping):
            return False
        _validate_loaded_derived_protein_fasta_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        if not _source_identity_matches(
            manifest_payload=manifest_payload,
            source_metadata_snapshot_root_directory=(
                _as_path(source_metadata_snapshot_root_directory)
                if source_metadata_snapshot_root_directory is not None
                else None
            ),
            source_selection_snapshot_root_directory=(
                _as_path(source_selection_snapshot_root_directory)
                if source_selection_snapshot_root_directory is not None
                else None
            ),
            selection_artifact_type=(
                selection_artifact_type
                or str(manifest_payload.get("source_selection_artifact_type"))
            ),
            selection_uid_list_file_name=(
                selection_uid_list_file_name or ""
            ),
        ):
            return False
        return True
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False


def resolve_derived_protein_fasta_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    source_selection_snapshot_root_directory: PathLike,
    selection_artifact_type: str,
    selection_uid_list_file_name: str,
    record_selection_rule: str,
    record_selection_config_sha256: str,
    dataset_kind: str,
    sequence_line_width: int = 60,
    record_order: str = DerivedFastaRecordOrder.AS_SELECTED.value,
    drop_missing_uids: bool = False,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )
    resolved_source_selection_snapshot_root_directory = _as_path(
        source_selection_snapshot_root_directory
    )

    latest_is_available = latest_derived_protein_fasta_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            resolved_source_metadata_snapshot_root_directory
        ),
        source_selection_snapshot_root_directory=(
            resolved_source_selection_snapshot_root_directory
        ),
        selection_artifact_type=selection_artifact_type,
        selection_uid_list_file_name=selection_uid_list_file_name,
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest derived protein FASTA snapshot was found for the "
                "current metadata and selection snapshots. Run the workflow "
                "once with snapshot_mode='create_new' before using 'reuse_latest'."
            )
        return load_latest_derived_protein_fasta_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print(
            "Latest derived protein FASTA snapshot is available. Reusing frozen "
            "snapshot."
        )
        return load_latest_derived_protein_fasta_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: "
            "'create_new', 'reuse_latest', 'reuse_latest_or_create'."
        )

    saved_result = save_derived_protein_fasta_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            resolved_source_metadata_snapshot_root_directory
        ),
        source_selection_snapshot_root_directory=(
            resolved_source_selection_snapshot_root_directory
        ),
        selection_artifact_type=selection_artifact_type,
        selection_uid_list_file_name=selection_uid_list_file_name,
        record_selection_rule=record_selection_rule,
        record_selection_config_sha256=record_selection_config_sha256,
        dataset_kind=dataset_kind,
        sequence_line_width=sequence_line_width,
        record_order=record_order,
        drop_missing_uids=drop_missing_uids,
        update_latest_directory=update_latest_directory,
    )
    return load_derived_protein_fasta_snapshot_by_directory(
        snapshot_directory=saved_result.snapshot_directory
    )


def list_saved_derived_protein_fasta_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory
    )


def get_most_recent_derived_protein_fasta_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
