from __future__ import annotations

import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias, cast

from src.pago_pipeline.ncbi_metadata_csv import (
    NCBIProteinMetadataCsvExportResult,
    export_ncbi_protein_metadata_csv,
)
from src.pago_pipeline.ncbi_metadata_qc import (
    NCBIProteinMetadataCsvQcResult,
    run_ncbi_protein_metadata_csv_qc,
)
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    load_latest_xml_snapshot,
    load_xml_snapshot_by_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path


@dataclass(frozen=True)
class NCBIProteinMetadataSnapshotResult:
    """
    Summary returned after creating one immutable metadata snapshot.
    """

    snapshot_directory: Path
    manifest_file_path: Path
    csv_file_path: Path
    qc_report_file_path: Path
    source_xml_snapshot_directory: Path
    source_xml_file_path: Path
    metadata_export_result: NCBIProteinMetadataCsvExportResult
    metadata_qc_result: NCBIProteinMetadataCsvQcResult


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _build_metadata_snapshot_manifest(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    export_result: NCBIProteinMetadataCsvExportResult,
    qc_result: NCBIProteinMetadataCsvQcResult,
    qc_report_file_path: PathLike,
    source_xml_snapshot_payload: Mapping[str, object],
    reference_metadata_manifest_file_path: Optional[PathLike],
) -> dict[str, object]:
    source_xml_manifest = source_xml_snapshot_payload.get("manifest")
    if not isinstance(source_xml_manifest, Mapping):
        raise RuntimeError(
            "Resolved XML snapshot payload is missing a valid source manifest."
        )

    source_xml_manifest_file_path = source_xml_snapshot_payload.get("manifest_file_path")
    resolved_qc_report_file_path = _as_path(qc_report_file_path)

    manifest_payload: dict[str, object] = {
        "snapshot_format_version": "1.0",
        "artifact_type": "ncbi_protein_metadata_snapshot",
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": "manifest.json",
        "csv_file_name": export_result.csv_file_path.name,
        "csv_file_sha256": export_result.csv_file_sha256,
        "qc_report_file_name": resolved_qc_report_file_path.name,
        "qc_report_file_sha256": sha256_of_file(
            input_file_path=resolved_qc_report_file_path
        ),
        "row_count": export_result.row_count,
        "column_count": export_result.column_count,
        "columns": export_result.columns,
        "max_taxonomy_depth": export_result.max_taxonomy_depth,
        "observed_feature_keys": export_result.observed_feature_keys,
        "observed_feature_qualifiers": export_result.observed_feature_qualifiers,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_xml_snapshot_relative_path": source_xml_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_xml_snapshot_directory_name": source_xml_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_xml_file_name": export_result.xml_file_path.name,
        "source_xml_file_sha256": export_result.xml_file_sha256,
        "source_xml_retrieved_at_utc": source_xml_manifest.get("retrieved_at_utc"),
        "search_query": source_xml_manifest.get("search_query"),
        "translated_query": source_xml_manifest.get("translated_query"),
        "qc_summary": {
            "row_count": qc_result.row_count,
            "column_count": qc_result.column_count,
            "empty_protein_uid_count": qc_result.empty_protein_uid_count,
            "duplicate_protein_uid_count": qc_result.duplicate_protein_uid_count,
            "fully_empty_column_count": qc_result.fully_empty_column_count,
            "schema_match_with_manifest": qc_result.schema_match_with_manifest,
            "row_count_match_with_manifest": qc_result.row_count_match_with_manifest,
            "row_count_match_with_source_xml": qc_result.row_count_match_with_source_xml,
            "normalization_collision_count": qc_result.normalization_collision_count,
        },
    }

    if isinstance(source_xml_manifest_file_path, (str, Path)):
        manifest_payload["source_xml_snapshot_manifest_sha256"] = sha256_of_file(
            input_file_path=source_xml_manifest_file_path
        )

    if reference_metadata_manifest_file_path is not None:
        resolved_reference_manifest_file_path = _as_path(
            reference_metadata_manifest_file_path
        )
        manifest_payload["reference_metadata_manifest_file_path"] = str(
            resolved_reference_manifest_file_path
        )
        manifest_payload["reference_metadata_manifest_sha256"] = sha256_of_file(
            input_file_path=resolved_reference_manifest_file_path
        )

    return manifest_payload


def save_ncbi_protein_metadata_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_xml_snapshot_directory: PathLike,
    reference_metadata_manifest_file_path: Optional[PathLike] = None,
    update_latest_directory: bool = True,
) -> NCBIProteinMetadataSnapshotResult:
    """
    Persist one immutable metadata CSV snapshot plus its QC report.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    source_xml_snapshot_payload = load_xml_snapshot_by_directory(
        snapshot_directory=source_xml_snapshot_directory,
    )
    source_xml_file_path = source_xml_snapshot_payload.get("xml_file_path")
    if not isinstance(source_xml_file_path, Path):
        raise RuntimeError(
            "Resolved XML snapshot payload is missing a valid xml_file_path."
        )

    source_xml_manifest = source_xml_snapshot_payload.get("manifest")
    if not isinstance(source_xml_manifest, Mapping):
        raise RuntimeError(
            "Resolved XML snapshot payload is missing a valid source manifest."
        )

    snapshot_created_at_utc = _current_utc_timestamp()
    search_query = str(
        source_xml_manifest.get("search_query")
        or source_xml_manifest.get("translated_query")
        or source_xml_file_path.name
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
        csv_file_path = immutable_snapshot_directory / "protein_metadata.csv"
        manifest_file_path = immutable_snapshot_directory / "manifest.json"
        qc_report_file_path = immutable_snapshot_directory / "qc_report.json"

        metadata_export_result = export_ncbi_protein_metadata_csv(
            xml_file_path=source_xml_file_path,
            output_csv_file_path=csv_file_path,
            output_manifest_file_path=manifest_file_path,
            source_snapshot_payload=cast(
                Mapping[str, object], source_xml_snapshot_payload
            ),
        )
        metadata_qc_result = run_ncbi_protein_metadata_csv_qc(
            csv_file_path=csv_file_path,
            metadata_manifest_file_path=manifest_file_path,
            source_xml_file_path=source_xml_file_path,
            source_xml_manifest_payload=cast(
                Mapping[str, object], source_xml_manifest
            ),
            reference_metadata_manifest_file_path=reference_metadata_manifest_file_path,
            output_report_file_path=qc_report_file_path,
        )

        manifest_payload = _build_metadata_snapshot_manifest(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            export_result=metadata_export_result,
            qc_result=metadata_qc_result,
            qc_report_file_path=qc_report_file_path,
            source_xml_snapshot_payload=cast(
                Mapping[str, object], source_xml_snapshot_payload
            ),
            reference_metadata_manifest_file_path=reference_metadata_manifest_file_path,
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
                    (csv_file_path, "protein_metadata.csv"),
                    (manifest_file_path, "manifest.json"),
                    (qc_report_file_path, "qc_report.json"),
                ],
            )

    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return NCBIProteinMetadataSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        manifest_file_path=manifest_file_path,
        csv_file_path=csv_file_path,
        qc_report_file_path=qc_report_file_path,
        source_xml_snapshot_directory=_as_path(source_xml_snapshot_directory),
        source_xml_file_path=source_xml_file_path,
        metadata_export_result=metadata_export_result,
        metadata_qc_result=metadata_qc_result,
    )


def _validate_loaded_metadata_snapshot_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Mapping[str, object],
) -> tuple[Path, Path]:
    resolved_snapshot_directory = _as_path(snapshot_directory)

    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "ncbi_protein_metadata_snapshot":
        raise RuntimeError(
            "Saved metadata snapshot manifest artifact_type mismatch. "
            f"Expected 'ncbi_protein_metadata_snapshot', got {artifact_type!r}."
        )

    csv_file_name = manifest_payload.get("csv_file_name")
    if not isinstance(csv_file_name, str) or not csv_file_name.strip():
        raise RuntimeError(
            "Saved metadata snapshot manifest must define a non-empty csv_file_name."
        )

    qc_report_file_name = manifest_payload.get("qc_report_file_name")
    if not isinstance(qc_report_file_name, str) or not qc_report_file_name.strip():
        raise RuntimeError(
            "Saved metadata snapshot manifest must define a non-empty "
            "qc_report_file_name."
        )

    resolved_csv_file_path = resolved_snapshot_directory / csv_file_name
    resolved_qc_report_file_path = resolved_snapshot_directory / qc_report_file_name

    if not resolved_csv_file_path.exists():
        raise FileNotFoundError(
            f"Saved metadata snapshot CSV file not found: {resolved_csv_file_path}."
        )

    if not resolved_qc_report_file_path.exists():
        raise FileNotFoundError(
            "Saved metadata snapshot QC report file not found: "
            f"{resolved_qc_report_file_path}."
        )

    expected_csv_file_sha256 = manifest_payload.get("csv_file_sha256")
    if expected_csv_file_sha256 is not None:
        actual_csv_file_sha256 = sha256_of_file(
            input_file_path=resolved_csv_file_path,
        )
        if actual_csv_file_sha256 != expected_csv_file_sha256:
            raise RuntimeError(
                "Saved metadata snapshot CSV file hash mismatch. "
                f"Expected {expected_csv_file_sha256}, got {actual_csv_file_sha256}."
            )

    expected_qc_report_file_sha256 = manifest_payload.get("qc_report_file_sha256")
    if expected_qc_report_file_sha256 is not None:
        actual_qc_report_file_sha256 = sha256_of_file(
            input_file_path=resolved_qc_report_file_path,
        )
        if actual_qc_report_file_sha256 != expected_qc_report_file_sha256:
            raise RuntimeError(
                "Saved metadata snapshot QC report file hash mismatch. "
                f"Expected {expected_qc_report_file_sha256}, got "
                f"{actual_qc_report_file_sha256}."
            )

    return resolved_csv_file_path, resolved_qc_report_file_path


def load_metadata_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    """
    Load and validate one saved metadata snapshot directory.
    """

    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved metadata snapshot manifest must deserialize into a JSON object."
        )

    csv_file_path, qc_report_file_path = _validate_loaded_metadata_snapshot_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "csv_file_path": csv_file_path,
        "qc_report_file_path": qc_report_file_path,
    }


def load_latest_metadata_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    """
    Load and validate the convenience latest metadata snapshot copy.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    return load_metadata_snapshot_by_directory(snapshot_directory=latest_directory)


def get_snapshot_metadata_csv_file_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the metadata CSV path inside a specific metadata snapshot directory.
    """

    resolved_snapshot_directory = _as_path(snapshot_directory)
    return resolved_snapshot_directory / "protein_metadata.csv"


def get_snapshot_metadata_qc_report_file_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the QC report path inside a specific metadata snapshot directory.
    """

    resolved_snapshot_directory = _as_path(snapshot_directory)
    return resolved_snapshot_directory / "qc_report.json"


def get_latest_metadata_csv_file_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the metadata CSV path for the convenience latest snapshot copy.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return resolved_snapshot_root_directory / "latest" / "protein_metadata.csv"


def get_latest_metadata_qc_report_file_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the QC report path for the convenience latest snapshot copy.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return resolved_snapshot_root_directory / "latest" / "qc_report.json"


def list_saved_metadata_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    """
    List saved immutable metadata snapshot directories in lexical order.
    """

    return list_saved_snapshot_directories(snapshot_root_directory=snapshot_root_directory)


def get_most_recent_metadata_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    """
    Return the most recent immutable metadata snapshot directory, or None.
    """

    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )


def latest_metadata_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    """
    Return True only if the convenience latest metadata snapshot copy is complete.
    """

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / "manifest.json"
    latest_csv_file_path = latest_directory / "protein_metadata.csv"
    latest_qc_report_file_path = latest_directory / "qc_report.json"

    return (
        latest_directory.exists()
        and latest_manifest_file_path.exists()
        and latest_csv_file_path.exists()
        and latest_qc_report_file_path.exists()
    )


def _metadata_snapshot_matches_source_xml(
    *,
    metadata_snapshot_payload: dict[str, object],
    source_xml_snapshot_payload: dict[str, object],
) -> bool:
    metadata_manifest = metadata_snapshot_payload.get("manifest")
    source_xml_manifest = source_xml_snapshot_payload.get("manifest")
    if not isinstance(metadata_manifest, Mapping) or not isinstance(
        source_xml_manifest, Mapping
    ):
        return False

    return (
        metadata_manifest.get("source_xml_file_sha256")
        == source_xml_manifest.get("xml_file_sha256")
    )


def resolve_ncbi_protein_metadata_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_xml_snapshot_root_directory: PathLike,
    reference_metadata_manifest_file_path: Optional[PathLike] = None,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    """
    Resolve the active metadata snapshot payload according to the requested mode.
    """

    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_xml_snapshot_root_directory = _as_path(
        source_xml_snapshot_root_directory
    )

    latest_is_available = latest_metadata_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            latest_directory = resolved_snapshot_root_directory / "latest"
            latest_manifest_file_path = latest_directory / "manifest.json"
            latest_csv_file_path = latest_directory / "protein_metadata.csv"
            latest_qc_report_file_path = latest_directory / "qc_report.json"

            if not latest_directory.exists():
                raise FileNotFoundError(
                    "No latest metadata snapshot directory was found. Run the "
                    "workflow once with snapshot_mode='create_new' before using "
                    "'reuse_latest'."
                )

            if not latest_manifest_file_path.exists():
                raise FileNotFoundError(
                    f"Latest metadata snapshot manifest not found: "
                    f"{latest_manifest_file_path}. Run the workflow once with "
                    f"snapshot_mode='create_new' to create it."
                )

            if not latest_csv_file_path.exists():
                raise FileNotFoundError(
                    f"Latest metadata snapshot CSV file not found: "
                    f"{latest_csv_file_path}. Run the workflow once with "
                    f"snapshot_mode='create_new' to create it."
                )

            raise FileNotFoundError(
                "Latest metadata snapshot QC report file not found: "
                f"{latest_qc_report_file_path}. Run the workflow once with "
                f"snapshot_mode='create_new' to create it."
            )

        return load_latest_metadata_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
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
        source_xml_snapshot_payload = load_latest_xml_snapshot(
            snapshot_root_directory=resolved_source_xml_snapshot_root_directory,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No reusable source XML snapshot was found for the metadata workflow. "
            "Create the upstream XML snapshot before running the metadata step."
        ) from exc

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        latest_metadata_snapshot_payload = load_latest_metadata_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )
        if _metadata_snapshot_matches_source_xml(
            metadata_snapshot_payload=latest_metadata_snapshot_payload,
            source_xml_snapshot_payload=source_xml_snapshot_payload,
        ):
            print("Latest metadata snapshot is available. Reusing frozen snapshot.")
            return latest_metadata_snapshot_payload

        print(
            "Latest metadata snapshot is available but does not match the "
            "latest source XML snapshot. Creating a new frozen snapshot."
        )

    saved_snapshot = save_ncbi_protein_metadata_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_xml_snapshot_directory=cast(
            Path, source_xml_snapshot_payload["snapshot_directory"]
        ),
        reference_metadata_manifest_file_path=reference_metadata_manifest_file_path,
        update_latest_directory=update_latest_directory,
    )

    return load_metadata_snapshot_by_directory(
        snapshot_directory=saved_snapshot.snapshot_directory,
    )
