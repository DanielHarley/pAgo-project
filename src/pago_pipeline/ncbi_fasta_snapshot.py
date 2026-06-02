from __future__ import annotations

import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias, cast

from src.pago_pipeline.metadata_to_fasta import (
    MetadataFastaExportResult,
    export_metadata_csv_to_fasta,
)
from src.pago_pipeline.ncbi_metadata_snapshot import (
    load_latest_metadata_snapshot,
    load_metadata_snapshot_by_directory,
)
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path


@dataclass(frozen=True)
class NCBIProteinFastaSnapshotResult:
    """
    Summary returned after creating one immutable FASTA snapshot.
    """

    snapshot_directory: Path
    manifest_file_path: Path
    fasta_file_path: Path
    source_metadata_snapshot_directory: Path
    source_metadata_csv_file_path: Path
    fasta_export_result: MetadataFastaExportResult


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _build_fasta_snapshot_manifest(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    export_result: MetadataFastaExportResult,
    source_metadata_snapshot_payload: Mapping[str, object],
) -> dict[str, object]:
    source_metadata_manifest = source_metadata_snapshot_payload.get("manifest")
    if not isinstance(source_metadata_manifest, Mapping):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing a valid source manifest."
        )

    source_metadata_manifest_file_path = source_metadata_snapshot_payload.get(
        "manifest_file_path"
    )

    manifest_payload: dict[str, object] = {
        "snapshot_format_version": "1.0",
        "artifact_type": "ncbi_protein_fasta_snapshot",
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": "manifest.json",
        "fasta_file_name": export_result.fasta_file_path.name,
        "fasta_file_sha256": export_result.fasta_file_sha256,
        "row_count": export_result.row_count,
        "fasta_record_count": export_result.fasta_record_count,
        "skipped_missing_sequence_count": export_result.skipped_missing_sequence_count,
        "sequence_line_width": export_result.sequence_line_width,
        "required_columns": export_result.required_columns,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_metadata_snapshot_relative_path": source_metadata_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_metadata_snapshot_directory_name": source_metadata_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_metadata_csv_file_name": source_metadata_manifest.get("csv_file_name"),
        "source_metadata_csv_file_sha256": source_metadata_manifest.get(
            "csv_file_sha256"
        ),
        "source_metadata_row_count": source_metadata_manifest.get("row_count"),
        "source_xml_snapshot_relative_path": source_metadata_manifest.get(
            "source_xml_snapshot_relative_path"
        ),
        "source_xml_snapshot_directory_name": source_metadata_manifest.get(
            "source_xml_snapshot_directory_name"
        ),
        "source_xml_file_name": source_metadata_manifest.get("source_xml_file_name"),
        "source_xml_file_sha256": source_metadata_manifest.get("source_xml_file_sha256"),
        "source_xml_retrieved_at_utc": source_metadata_manifest.get(
            "source_xml_retrieved_at_utc"
        ),
        "search_query": source_metadata_manifest.get("search_query"),
        "translated_query": source_metadata_manifest.get("translated_query"),
    }

    if isinstance(source_metadata_manifest_file_path, (str, Path)):
        manifest_payload["source_metadata_manifest_sha256"] = sha256_of_file(
            input_file_path=source_metadata_manifest_file_path
        )

    return manifest_payload


def save_ncbi_protein_fasta_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_directory: PathLike,
    sequence_line_width: int = 60,
    update_latest_directory: bool = True,
) -> NCBIProteinFastaSnapshotResult:
    """
    Persist one immutable protein FASTA snapshot derived from a metadata snapshot.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    source_metadata_snapshot_payload = load_metadata_snapshot_by_directory(
        snapshot_directory=source_metadata_snapshot_directory,
    )
    source_metadata_csv_file_path = source_metadata_snapshot_payload.get("csv_file_path")
    if not isinstance(source_metadata_csv_file_path, Path):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing a valid csv_file_path."
        )

    source_metadata_manifest = source_metadata_snapshot_payload.get("manifest")
    if not isinstance(source_metadata_manifest, Mapping):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing a valid source manifest."
        )

    snapshot_created_at_utc = _current_utc_timestamp()
    search_query = str(
        source_metadata_manifest.get("search_query")
        or source_metadata_manifest.get("translated_query")
        or source_metadata_csv_file_path.name
    )
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query=search_query,
    )

    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    immutable_snapshot_complete = False

    try:
        immutable_snapshot_relative_path = str(
            Path("snapshots") / snapshot_directory_name
        )
        fasta_file_path = immutable_snapshot_directory / "protein_sequences.fasta"
        manifest_file_path = immutable_snapshot_directory / "manifest.json"

        fasta_export_result = export_metadata_csv_to_fasta(
            metadata_csv_file_path=source_metadata_csv_file_path,
            output_fasta_file_path=fasta_file_path,
            output_manifest_file_path=manifest_file_path,
            sequence_line_width=sequence_line_width,
            source_metadata_manifest_payload=cast(
                Mapping[str, object], source_metadata_snapshot_payload["manifest"]
            ),
            source_metadata_manifest_file_path=cast(
                Path, source_metadata_snapshot_payload["manifest_file_path"]
            ),
        )

        manifest_payload = _build_fasta_snapshot_manifest(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            export_result=fasta_export_result,
            source_metadata_snapshot_payload=cast(
                Mapping[str, object], source_metadata_snapshot_payload
            ),
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
                    (fasta_file_path, "protein_sequences.fasta"),
                    (manifest_file_path, "manifest.json"),
                ],
            )

    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return NCBIProteinFastaSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        manifest_file_path=manifest_file_path,
        fasta_file_path=fasta_file_path,
        source_metadata_snapshot_directory=_as_path(source_metadata_snapshot_directory),
        source_metadata_csv_file_path=source_metadata_csv_file_path,
        fasta_export_result=fasta_export_result,
    )


def _validate_loaded_fasta_snapshot_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Mapping[str, object],
) -> Path:
    resolved_snapshot_directory = _as_path(snapshot_directory)

    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "ncbi_protein_fasta_snapshot":
        raise RuntimeError(
            "Saved FASTA snapshot manifest artifact_type mismatch. "
            f"Expected 'ncbi_protein_fasta_snapshot', got {artifact_type!r}."
        )

    fasta_file_name = manifest_payload.get("fasta_file_name")
    if not isinstance(fasta_file_name, str) or not fasta_file_name.strip():
        raise RuntimeError(
            "Saved FASTA snapshot manifest must define a non-empty fasta_file_name."
        )

    resolved_fasta_file_path = resolved_snapshot_directory / fasta_file_name
    if not resolved_fasta_file_path.exists():
        raise FileNotFoundError(
            f"Saved FASTA snapshot file not found: {resolved_fasta_file_path}."
        )

    expected_fasta_file_sha256 = manifest_payload.get("fasta_file_sha256")
    if expected_fasta_file_sha256 is not None:
        actual_fasta_file_sha256 = sha256_of_file(
            input_file_path=resolved_fasta_file_path
        )
        if actual_fasta_file_sha256 != expected_fasta_file_sha256:
            raise RuntimeError(
                "Saved FASTA snapshot file hash mismatch. "
                f"Expected {expected_fasta_file_sha256}, got {actual_fasta_file_sha256}."
            )

    return resolved_fasta_file_path


def load_fasta_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    """
    Load and validate one saved FASTA snapshot directory.
    """

    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved FASTA snapshot manifest must deserialize into a JSON object."
        )

    fasta_file_path = _validate_loaded_fasta_snapshot_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "fasta_file_path": fasta_file_path,
    }


def load_latest_fasta_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    """
    Load and validate the convenience latest FASTA snapshot copy.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    return load_fasta_snapshot_by_directory(snapshot_directory=latest_directory)


def get_snapshot_fasta_file_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the FASTA file path inside a specific FASTA snapshot directory.
    """

    resolved_snapshot_directory = _as_path(snapshot_directory)
    return resolved_snapshot_directory / "protein_sequences.fasta"


def get_latest_fasta_file_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the FASTA file path for the convenience latest FASTA snapshot copy.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return resolved_snapshot_root_directory / "latest" / "protein_sequences.fasta"


def list_saved_fasta_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    """
    List saved immutable FASTA snapshot directories in lexical order.
    """

    return list_saved_snapshot_directories(snapshot_root_directory=snapshot_root_directory)


def get_most_recent_fasta_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    """
    Return the most recent immutable FASTA snapshot directory, or None.
    """

    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )


def latest_fasta_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: Optional[PathLike] = None,
) -> bool:
    """
    Return True only if the convenience latest FASTA snapshot copy is complete.

    When source_metadata_snapshot_root_directory is provided, the latest FASTA
    snapshot must also match the current latest metadata input snapshot.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / "manifest.json"
    if not latest_directory.exists() or not latest_manifest_file_path.exists():
        return False

    try:
        manifest_payload = read_json_file(input_file_path=latest_manifest_file_path)
        if not isinstance(manifest_payload, Mapping):
            return False
        _validate_loaded_fasta_snapshot_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        if source_metadata_snapshot_root_directory is not None:
            source_metadata_snapshot_payload = load_latest_metadata_snapshot(
                snapshot_root_directory=source_metadata_snapshot_root_directory,
            )
            if not _source_metadata_snapshot_identity_matches(
                manifest_payload=manifest_payload,
                source_metadata_snapshot_payload=source_metadata_snapshot_payload,
            ):
                return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False

    return True


def _source_metadata_snapshot_identity_matches(
    *,
    manifest_payload: Mapping[str, object],
    source_metadata_snapshot_payload: Mapping[str, object],
) -> bool:
    expected_manifest_sha256 = manifest_payload.get(
        "source_metadata_manifest_sha256"
    )
    expected_csv_sha256 = manifest_payload.get("source_metadata_csv_file_sha256")
    source_manifest_file_path = source_metadata_snapshot_payload.get(
        "manifest_file_path"
    )
    source_csv_file_path = source_metadata_snapshot_payload.get("csv_file_path")

    if not isinstance(expected_manifest_sha256, str) or not expected_manifest_sha256:
        return False
    if not isinstance(expected_csv_sha256, str) or not expected_csv_sha256:
        return False
    if not isinstance(source_manifest_file_path, Path):
        return False
    if not isinstance(source_csv_file_path, Path):
        return False

    return (
        sha256_of_file(input_file_path=source_manifest_file_path)
        == expected_manifest_sha256
        and sha256_of_file(input_file_path=source_csv_file_path)
        == expected_csv_sha256
    )


def resolve_ncbi_protein_fasta_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    sequence_line_width: int = 60,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    """
    Resolve the active FASTA snapshot payload according to the requested mode.
    """

    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )

    latest_is_available = latest_fasta_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            resolved_source_metadata_snapshot_root_directory
        ),
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest FASTA snapshot directory was found for the current "
                "source metadata snapshot. Run the workflow once with "
                "snapshot_mode='create_new' before using 'reuse_latest'."
            )

        return load_latest_fasta_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print("Latest FASTA snapshot is available. Reusing frozen snapshot.")
        return load_latest_fasta_snapshot(
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

    try:
        source_metadata_snapshot_payload = load_latest_metadata_snapshot(
            snapshot_root_directory=resolved_source_metadata_snapshot_root_directory
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No reusable source metadata snapshot was found for the FASTA workflow. "
            "Create the upstream metadata snapshot before running the FASTA step."
        ) from exc

    saved_snapshot = save_ncbi_protein_fasta_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_directory=cast(
            Path, source_metadata_snapshot_payload["snapshot_directory"]
        ),
        sequence_line_width=sequence_line_width,
        update_latest_directory=update_latest_directory,
    )

    return load_fasta_snapshot_by_directory(
        snapshot_directory=saved_snapshot.snapshot_directory
    )
