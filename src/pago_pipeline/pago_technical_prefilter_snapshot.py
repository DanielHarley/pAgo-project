from __future__ import annotations

import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

import pandas as pd

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
from src.pago_pipeline.pago_technical_prefilter import (
    DEFAULT_LENGTH_COLUMN,
    DEFAULT_PROTEIN_UID_COLUMN,
    DEFAULT_SEQUENCE_COLUMN,
    build_pago_technical_prefilter_partition,
    build_technical_prefilter_counts_dataframe,
    build_technical_prefilter_policy_payload,
    build_technical_prefilter_policy_sha256,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

ARTIFACT_TYPE = "pago_technical_prefilter"
SNAPSHOT_FORMAT_VERSION = "1.0"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
DEFAULT_RETAINED_PROTEIN_UIDS_FILE_NAME = "retained_protein_uids.txt"
DEFAULT_RETAINED_RECORDS_FILE_NAME = "retained_records.csv"
DEFAULT_EXCLUDED_TECHNICAL_FILE_NAME = "excluded_technical.csv"
DEFAULT_PREFILTER_COUNTS_FILE_NAME = "prefilter_counts.csv"

_SNAPSHOT_DIRECTORY_QUERY_LITERAL = "pago_technical_prefilter"


@dataclass(frozen=True)
class PagoTechnicalPrefilterSnapshotResult:
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    retained_protein_uids_file_path: Path
    retained_records_file_path: Path
    excluded_technical_file_path: Path
    prefilter_counts_file_path: Path
    input_record_count: int
    retained_record_count: int
    excluded_record_count: int


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


def _write_text_lines_atomic(
    *,
    text_lines: list[str],
    output_file_path: PathLike,
) -> Path:
    resolved_output_file_path = _as_path(output_file_path)
    resolved_output_file_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=resolved_output_file_path.parent,
        encoding="utf-8",
        newline="\n",
        suffix=".txt",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)
        for text_line in text_lines:
            temporary_file.write(f"{text_line}\n")
    resolved_temporary_file_path.replace(resolved_output_file_path)
    return resolved_output_file_path


def _resolve_source_metadata_manifest(
    *,
    source_metadata_snapshot_payload: Mapping[str, object],
) -> tuple[Mapping[str, object], Path, Path]:
    source_metadata_manifest = source_metadata_snapshot_payload.get("manifest")
    if not isinstance(source_metadata_manifest, Mapping):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing a valid source manifest."
        )
    source_metadata_manifest_file_path = source_metadata_snapshot_payload.get(
        "manifest_file_path"
    )
    source_metadata_csv_file_path = source_metadata_snapshot_payload.get("csv_file_path")
    if not isinstance(source_metadata_manifest_file_path, (str, Path)):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing manifest_file_path."
        )
    if not isinstance(source_metadata_csv_file_path, (str, Path)):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing csv_file_path."
        )
    return (
        source_metadata_manifest,
        _as_path(source_metadata_manifest_file_path),
        _as_path(source_metadata_csv_file_path),
    )


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    source_metadata_manifest: Mapping[str, object],
    source_metadata_manifest_file_path: Path,
    source_metadata_csv_file_path: Path,
    output_file_path_by_key: Mapping[str, Path],
    input_record_count: int,
    retained_record_count: int,
    excluded_record_count: int,
    counts_by_decision: Mapping[str, int],
    policy_payload: Mapping[str, object],
    policy_sha256: str,
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
        "input_record_count": int(input_record_count),
        "retained_record_count": int(retained_record_count),
        "excluded_record_count": int(excluded_record_count),
        "counts_by_decision": {
            str(key): int(value) for key, value in counts_by_decision.items()
        },
        "technical_prefilter_policy": dict(policy_payload),
        "technical_prefilter_policy_sha256": policy_sha256,
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
        "source_metadata_row_count": source_metadata_manifest.get("row_count"),
        "search_query": source_metadata_manifest.get("search_query"),
        "translated_query": source_metadata_manifest.get("translated_query"),
        "source_xml_snapshot_relative_path": source_metadata_manifest.get(
            "source_xml_snapshot_relative_path"
        ),
        "output_files": output_files,
    }


def save_pago_technical_prefilter_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> PagoTechnicalPrefilterSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )

    try:
        source_metadata_snapshot_payload = load_latest_metadata_snapshot(
            snapshot_root_directory=resolved_source_metadata_snapshot_root_directory
        )
    except FileNotFoundError as error:
        raise FileNotFoundError(
            "No reusable metadata snapshot was found for the technical prefilter. "
            "Run the upstream metadata extraction stage before running the "
            "technical prefilter."
        ) from error

    (
        source_metadata_manifest,
        source_metadata_manifest_file_path,
        source_metadata_csv_file_path,
    ) = _resolve_source_metadata_manifest(
        source_metadata_snapshot_payload=source_metadata_snapshot_payload
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
        metadata_dataframe = pd.read_csv(
            source_metadata_csv_file_path,
            low_memory=False,
            dtype={DEFAULT_PROTEIN_UID_COLUMN: "string"},
        )
        partition_result = build_pago_technical_prefilter_partition(
            metadata_dataframe=metadata_dataframe,
        )
        counts_dataframe = build_technical_prefilter_counts_dataframe(
            result=partition_result,
        )

        output_file_path_by_key = {
            "retained_protein_uids_file": (
                immutable_snapshot_directory
                / DEFAULT_RETAINED_PROTEIN_UIDS_FILE_NAME
            ),
            "retained_records_file": (
                immutable_snapshot_directory / DEFAULT_RETAINED_RECORDS_FILE_NAME
            ),
            "excluded_technical_file": (
                immutable_snapshot_directory / DEFAULT_EXCLUDED_TECHNICAL_FILE_NAME
            ),
            "prefilter_counts_file": (
                immutable_snapshot_directory / DEFAULT_PREFILTER_COUNTS_FILE_NAME
            ),
        }

        _write_text_lines_atomic(
            text_lines=list(partition_result.retained_protein_uids),
            output_file_path=output_file_path_by_key["retained_protein_uids_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=partition_result.retained_records,
            output_file_path=output_file_path_by_key["retained_records_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=partition_result.excluded_records,
            output_file_path=output_file_path_by_key["excluded_technical_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=counts_dataframe,
            output_file_path=output_file_path_by_key["prefilter_counts_file"],
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
            output_file_path_by_key=output_file_path_by_key,
            input_record_count=partition_result.input_record_count,
            retained_record_count=int(len(partition_result.retained_records)),
            excluded_record_count=int(len(partition_result.excluded_records)),
            counts_by_decision=partition_result.counts_by_decision,
            policy_payload=build_technical_prefilter_policy_payload(),
            policy_sha256=build_technical_prefilter_policy_sha256(),
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
                        output_file_path_by_key["retained_protein_uids_file"],
                        DEFAULT_RETAINED_PROTEIN_UIDS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["retained_records_file"],
                        DEFAULT_RETAINED_RECORDS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["excluded_technical_file"],
                        DEFAULT_EXCLUDED_TECHNICAL_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["prefilter_counts_file"],
                        DEFAULT_PREFILTER_COUNTS_FILE_NAME,
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

    return PagoTechnicalPrefilterSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        retained_protein_uids_file_path=output_file_path_by_key[
            "retained_protein_uids_file"
        ],
        retained_records_file_path=output_file_path_by_key["retained_records_file"],
        excluded_technical_file_path=output_file_path_by_key["excluded_technical_file"],
        prefilter_counts_file_path=output_file_path_by_key["prefilter_counts_file"],
        input_record_count=partition_result.input_record_count,
        retained_record_count=int(len(partition_result.retained_records)),
        excluded_record_count=int(len(partition_result.excluded_records)),
    )


def _validate_loaded_pago_technical_prefilter_payload(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
) -> dict[str, Path]:
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != ARTIFACT_TYPE:
        raise RuntimeError(
            "Saved technical prefilter snapshot manifest artifact_type mismatch. "
            f"Expected {ARTIFACT_TYPE!r}, got {artifact_type!r}."
        )

    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved technical prefilter snapshot manifest snapshot_format_version "
            f"mismatch. Expected {SNAPSHOT_FORMAT_VERSION!r}, got "
            f"{snapshot_format_version!r}."
        )

    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError(
            "Saved technical prefilter snapshot manifest must define output_files."
        )

    resolved_output_file_path_by_key: dict[str, Path] = {}
    for key, entry in output_files.items():
        if not isinstance(entry, Mapping):
            raise RuntimeError(
                f"Saved technical prefilter output_files entry {key!r} must be an "
                "object."
            )
        file_name = entry.get("file_name")
        expected_sha256 = entry.get("sha256")
        if not isinstance(file_name, str) or not file_name.strip():
            raise RuntimeError(
                f"Saved technical prefilter output_files entry {key!r} must define "
                "a non-empty file_name."
            )
        if not isinstance(expected_sha256, str) or not expected_sha256.strip():
            raise RuntimeError(
                f"Saved technical prefilter output_files entry {key!r} must define "
                "a non-empty sha256."
            )
        resolved_file_path = snapshot_directory / file_name
        if not resolved_file_path.exists():
            raise FileNotFoundError(
                f"Saved technical prefilter snapshot file not found: "
                f"{resolved_file_path}."
            )
        actual_sha256 = sha256_of_file(input_file_path=resolved_file_path)
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Saved technical prefilter snapshot hash mismatch for {key!r}. "
                f"Expected {expected_sha256}, got {actual_sha256}."
            )
        resolved_output_file_path_by_key[key] = resolved_file_path

    return resolved_output_file_path_by_key


def load_pago_technical_prefilter_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved technical prefilter snapshot manifest must deserialize into a "
            "JSON object."
        )

    resolved_output_file_path_by_key = (
        _validate_loaded_pago_technical_prefilter_payload(
            snapshot_directory=resolved_snapshot_directory,
            manifest_payload=manifest_payload,
        )
    )

    retained_records_dataframe = pd.read_csv(
        resolved_output_file_path_by_key["retained_records_file"],
        low_memory=False,
        dtype={DEFAULT_PROTEIN_UID_COLUMN: "string"},
    )
    prefilter_counts_dataframe = pd.read_csv(
        resolved_output_file_path_by_key["prefilter_counts_file"]
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "retained_protein_uids_file_path": resolved_output_file_path_by_key[
            "retained_protein_uids_file"
        ],
        "retained_records_file_path": resolved_output_file_path_by_key[
            "retained_records_file"
        ],
        "excluded_technical_file_path": resolved_output_file_path_by_key[
            "excluded_technical_file"
        ],
        "prefilter_counts_file_path": resolved_output_file_path_by_key[
            "prefilter_counts_file"
        ],
        "retained_records": retained_records_dataframe,
        "prefilter_counts": prefilter_counts_dataframe,
    }


def load_latest_pago_technical_prefilter_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_pago_technical_prefilter_snapshot_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest"
    )


def _source_metadata_snapshot_identity_matches(
    *,
    manifest_payload: Mapping[str, object],
    source_metadata_snapshot_root_directory: Path,
) -> bool:
    source_metadata_snapshot_payload = load_latest_metadata_snapshot(
        snapshot_root_directory=source_metadata_snapshot_root_directory
    )
    source_metadata_manifest_file_path = source_metadata_snapshot_payload.get(
        "manifest_file_path"
    )
    if not isinstance(source_metadata_manifest_file_path, (str, Path)):
        return False
    current_manifest_sha256 = sha256_of_file(
        input_file_path=_as_path(source_metadata_manifest_file_path)
    )
    return (
        manifest_payload.get("source_metadata_manifest_sha256")
        == current_manifest_sha256
    )


def latest_pago_technical_prefilter_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: Optional[PathLike] = None,
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
        _validate_loaded_pago_technical_prefilter_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        if (
            manifest_payload.get("technical_prefilter_policy_sha256")
            != build_technical_prefilter_policy_sha256()
        ):
            return False
        if source_metadata_snapshot_root_directory is not None:
            if not _source_metadata_snapshot_identity_matches(
                manifest_payload=manifest_payload,
                source_metadata_snapshot_root_directory=_as_path(
                    source_metadata_snapshot_root_directory
                ),
            ):
                return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False

    return True


def resolve_pago_technical_prefilter_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )

    latest_is_available = latest_pago_technical_prefilter_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            resolved_source_metadata_snapshot_root_directory
        ),
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest technical prefilter snapshot was found for the current "
                "metadata snapshot. Run the workflow once with "
                "snapshot_mode='create_new' before using 'reuse_latest'."
            )
        return load_latest_pago_technical_prefilter_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print(
            "Latest technical prefilter snapshot is available. Reusing frozen "
            "snapshot."
        )
        return load_latest_pago_technical_prefilter_snapshot(
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

    saved_result = save_pago_technical_prefilter_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            resolved_source_metadata_snapshot_root_directory
        ),
        update_latest_directory=update_latest_directory,
    )
    return load_pago_technical_prefilter_snapshot_by_directory(
        snapshot_directory=saved_result.snapshot_directory
    )


def list_saved_pago_technical_prefilter_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory
    )


def get_most_recent_pago_technical_prefilter_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
