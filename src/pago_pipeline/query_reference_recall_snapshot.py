from __future__ import annotations

import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

import pandas as pd

from src.pago_pipeline.ncbi_metadata_snapshot import load_latest_metadata_snapshot
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.query_reference_recall import (
    compute_query_reference_recall,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

ARTIFACT_TYPE = "query_reference_recall"
SNAPSHOT_FORMAT_VERSION = "1.0"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
DEFAULT_SUMMARY_FILE_NAME = "reference_recall_summary.csv"
DEFAULT_DETAIL_FILE_NAME = "reference_recall_detail.csv"

_SNAPSHOT_DIRECTORY_QUERY_LITERAL = "query_reference_recall"


@dataclass(frozen=True)
class QueryReferenceRecallSnapshotResult:
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    summary_file_path: Path
    detail_file_path: Path
    reference_count: int
    recovered_count: int
    stratum_recall: dict[str, float | None]


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _write_dataframe_csv_atomic(
    *,
    dataframe: pd.DataFrame,
    output_file_path: PathLike,
) -> Path:
    resolved_output_file_path = _as_path(output_file_path)
    resolved_output_file_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=resolved_output_file_path.parent,
        encoding="utf-8",
        newline="",
        suffix=".csv",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)
    dataframe.to_csv(resolved_temporary_file_path, index=False)
    resolved_temporary_file_path.replace(resolved_output_file_path)
    return resolved_output_file_path


def _resolve_metadata_snapshot_paths(
    *,
    metadata_snapshot_payload: Mapping[str, object],
) -> tuple[Mapping[str, object], Path, Path]:
    manifest = metadata_snapshot_payload.get("manifest")
    manifest_file_path = metadata_snapshot_payload.get("manifest_file_path")
    csv_file_path = metadata_snapshot_payload.get("csv_file_path")
    if (
        not isinstance(manifest, Mapping)
        or not isinstance(manifest_file_path, (str, Path))
        or not isinstance(csv_file_path, (str, Path))
    ):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing manifest / csv path."
        )
    return manifest, _as_path(manifest_file_path), _as_path(csv_file_path)


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    source_metadata_manifest: Mapping[str, object],
    source_metadata_manifest_file_path: Path,
    source_metadata_csv_file_path: Path,
    query_recall_reference_set_csv_file_path: Path,
    output_file_path_by_key: Mapping[str, Path],
    reference_count: int,
    recovered_count: int,
    stratum_recall: Mapping[str, float | None],
    stratum_recall_status: Mapping[str, str],
) -> dict[str, object]:
    output_files = {
        key: {
            "file_name": file_path.name,
            "path": str(file_path),
            "sha256": sha256_of_file(input_file_path=file_path),
        }
        for key, file_path in output_file_path_by_key.items()
    }
    return {
        "artifact_type": ARTIFACT_TYPE,
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "reference_count": int(reference_count),
        "recovered_count": int(recovered_count),
        # None (JSON null), not 0.0, for a stratum with zero reference sequences.
        "stratum_recall": {
            str(key): (None if value is None else float(value))
            for key, value in stratum_recall.items()
        },
        "stratum_recall_status": {
            str(key): str(value) for key, value in stratum_recall_status.items()
        },
        "query_recall_reference_set_csv_path": str(
            query_recall_reference_set_csv_file_path
        ),
        "query_recall_reference_set_csv_sha256": sha256_of_file(
            input_file_path=query_recall_reference_set_csv_file_path
        ),
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
        "search_query": source_metadata_manifest.get("search_query"),
        "translated_query": source_metadata_manifest.get("translated_query"),
        "output_files": output_files,
    }


def save_query_reference_recall_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    query_recall_reference_set_csv_path: PathLike,
    update_latest_directory: bool = True,
) -> QueryReferenceRecallSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )
    resolved_reference_set_csv_path = _as_path(query_recall_reference_set_csv_path)
    if not resolved_reference_set_csv_path.exists():
        raise FileNotFoundError(
            "query_recall_reference_set CSV not found: "
            f"{resolved_reference_set_csv_path}."
        )

    try:
        metadata_snapshot_payload = load_latest_metadata_snapshot(
            snapshot_root_directory=resolved_source_metadata_snapshot_root_directory
        )
    except FileNotFoundError as error:
        raise FileNotFoundError(
            "No reusable metadata snapshot was found for query reference recall. "
            "Run the upstream metadata extraction stage first."
        ) from error

    (
        source_metadata_manifest,
        source_metadata_manifest_file_path,
        source_metadata_csv_file_path,
    ) = _resolve_metadata_snapshot_paths(
        metadata_snapshot_payload=metadata_snapshot_payload
    )

    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query=_SNAPSHOT_DIRECTORY_QUERY_LITERAL,
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
        reference_dataframe = pd.read_csv(
            resolved_reference_set_csv_path, dtype=str
        ).fillna("")
        retrieved_metadata_dataframe = pd.read_csv(
            source_metadata_csv_file_path,
            low_memory=False,
            usecols=lambda column_name: column_name
            in {"protein_uid", "gbseq__accession_version"},
            dtype=str,
        )
        recall_result = compute_query_reference_recall(
            reference_dataframe=reference_dataframe,
            retrieved_metadata_dataframe=retrieved_metadata_dataframe,
        )

        output_file_path_by_key = {
            "summary_file": immutable_snapshot_directory / DEFAULT_SUMMARY_FILE_NAME,
            "detail_file": immutable_snapshot_directory / DEFAULT_DETAIL_FILE_NAME,
        }
        _write_dataframe_csv_atomic(
            dataframe=recall_result.summary_dataframe,
            output_file_path=output_file_path_by_key["summary_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=recall_result.detail_dataframe,
            output_file_path=output_file_path_by_key["detail_file"],
        )

        manifest_file_path = (
            immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
        )
        manifest_payload = _build_manifest_payload(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            source_metadata_manifest=source_metadata_manifest,
            source_metadata_manifest_file_path=source_metadata_manifest_file_path,
            source_metadata_csv_file_path=source_metadata_csv_file_path,
            query_recall_reference_set_csv_file_path=resolved_reference_set_csv_path,
            output_file_path_by_key=output_file_path_by_key,
            reference_count=recall_result.reference_count,
            recovered_count=recall_result.recovered_count,
            stratum_recall=recall_result.stratum_recall,
            stratum_recall_status=recall_result.stratum_recall_status,
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
                    (
                        output_file_path_by_key["summary_file"],
                        DEFAULT_SUMMARY_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["detail_file"],
                        DEFAULT_DETAIL_FILE_NAME,
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

    return QueryReferenceRecallSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        summary_file_path=output_file_path_by_key["summary_file"],
        detail_file_path=output_file_path_by_key["detail_file"],
        reference_count=recall_result.reference_count,
        recovered_count=recall_result.recovered_count,
        stratum_recall=dict(recall_result.stratum_recall),
    )


def _validate_loaded_query_reference_recall_payload(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
) -> dict[str, Path]:
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != ARTIFACT_TYPE:
        raise RuntimeError(
            "Saved query reference recall snapshot artifact_type mismatch. "
            f"Expected {ARTIFACT_TYPE!r}, got {artifact_type!r}."
        )
    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved query reference recall snapshot snapshot_format_version "
            f"mismatch. Expected {SNAPSHOT_FORMAT_VERSION!r}, got "
            f"{snapshot_format_version!r}."
        )
    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError(
            "Saved query reference recall snapshot manifest must define "
            "output_files."
        )
    resolved: dict[str, Path] = {}
    for key, entry in output_files.items():
        if not isinstance(entry, Mapping):
            raise RuntimeError(f"output_files entry {key!r} must be an object.")
        file_name = entry.get("file_name")
        expected_sha256 = entry.get("sha256")
        if not isinstance(file_name, str) or not file_name.strip():
            raise RuntimeError(
                f"output_files entry {key!r} must define file_name."
            )
        if not isinstance(expected_sha256, str) or not expected_sha256.strip():
            raise RuntimeError(f"output_files entry {key!r} must define sha256.")
        file_path = snapshot_directory / file_name
        if not file_path.exists():
            raise FileNotFoundError(
                f"Saved query reference recall snapshot file not found: {file_path}."
            )
        actual_sha256 = sha256_of_file(input_file_path=file_path)
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Saved query reference recall snapshot hash mismatch for {key!r}."
            )
        resolved[key] = file_path
    return resolved


def load_query_reference_recall_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved query reference recall snapshot manifest must deserialize into "
            "a JSON object."
        )
    resolved = _validate_loaded_query_reference_recall_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )
    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "summary_file_path": resolved["summary_file"],
        "detail_file_path": resolved["detail_file"],
        "summary": pd.read_csv(resolved["summary_file"]),
        "detail": pd.read_csv(resolved["detail_file"]),
    }


def load_latest_query_reference_recall_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_query_reference_recall_snapshot_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest"
    )


def latest_query_reference_recall_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: Optional[PathLike] = None,
    query_recall_reference_set_csv_path: Optional[PathLike] = None,
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
        _validate_loaded_query_reference_recall_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        if query_recall_reference_set_csv_path is not None:
            resolved_reference_set_csv_path = _as_path(
                query_recall_reference_set_csv_path
            )
            if not resolved_reference_set_csv_path.exists():
                return False
            if manifest_payload.get(
                "query_recall_reference_set_csv_sha256"
            ) != sha256_of_file(input_file_path=resolved_reference_set_csv_path):
                return False
        if source_metadata_snapshot_root_directory is not None:
            metadata_snapshot_payload = load_latest_metadata_snapshot(
                snapshot_root_directory=_as_path(
                    source_metadata_snapshot_root_directory
                )
            )
            metadata_manifest_file_path = metadata_snapshot_payload.get(
                "manifest_file_path"
            )
            if not isinstance(metadata_manifest_file_path, (str, Path)):
                return False
            if manifest_payload.get(
                "source_metadata_manifest_sha256"
            ) != sha256_of_file(
                input_file_path=_as_path(metadata_manifest_file_path)
            ):
                return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False
    return True


def resolve_query_reference_recall_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    query_recall_reference_set_csv_path: PathLike,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )
    resolved_reference_set_csv_path = _as_path(query_recall_reference_set_csv_path)

    latest_is_available = latest_query_reference_recall_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            resolved_source_metadata_snapshot_root_directory
        ),
        query_recall_reference_set_csv_path=resolved_reference_set_csv_path,
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest query reference recall snapshot was found for the "
                "current metadata snapshot and reference set. Run the workflow "
                "once with snapshot_mode='create_new' before using 'reuse_latest'."
            )
        return load_latest_query_reference_recall_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print(
            "Latest query reference recall snapshot is available. Reusing frozen "
            "snapshot."
        )
        return load_latest_query_reference_recall_snapshot(
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

    saved_result = save_query_reference_recall_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            resolved_source_metadata_snapshot_root_directory
        ),
        query_recall_reference_set_csv_path=resolved_reference_set_csv_path,
        update_latest_directory=update_latest_directory,
    )
    return load_query_reference_recall_snapshot_by_directory(
        snapshot_directory=saved_result.snapshot_directory
    )


def list_saved_query_reference_recall_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory
    )


def get_most_recent_query_reference_recall_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
