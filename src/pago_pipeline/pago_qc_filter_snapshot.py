from __future__ import annotations

import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

import pandas as pd

from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.pago_qc_filter import (
    CLASSIC_PAGO_HIGH_PRECISION_DATASET_NAME,
    CLASSIC_PAGO_REVIEW_DATASET_NAME,
    EXCLUDED_DATASET_NAME,
    FILTERED_DATASET_DEFINITIONS,
    PIWI_RE_DATASET_NAME,
    build_filtered_dataset_counts_dataframe,
    build_filter_policy_payload,
    build_filter_policy_sha256,
    build_pago_qc_filtered_datasets,
)
from src.pago_pipeline.pago_qc_snapshot import (
    DEFAULT_MANIFEST_FILE_NAME,
    load_latest_pago_qc_evidence_inventory,
    load_pago_qc_evidence_inventory_by_directory,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

DEFAULT_FILTERED_DATASET_COUNTS_FILE_NAME = "filtered_dataset_counts.csv"
SNAPSHOT_FORMAT_VERSION = "1.0"
EXCLUDED_RECORDS_SEMANTICS = (
    "Records excluded from the conservative classic pAgo positive dataset; "
    "not a biologically validated negative class."
)
DATASET_DEFINITION_BY_NAME = {
    dataset_definition.dataset_name: dataset_definition
    for dataset_definition in FILTERED_DATASET_DEFINITIONS
}


@dataclass(frozen=True)
class PagoQcFilteredDatasetsResult:
    source_qc_snapshot_directory: Path
    source_qc_manifest_file_path: Path
    source_labelled_records_file_path: Path
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    classic_pago_high_precision_records_file_path: Path
    classic_pago_review_records_file_path: Path
    piwi_re_records_file_path: Path
    excluded_records_file_path: Path
    filtered_dataset_counts_file_path: Path
    source_record_count: int


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


def _build_output_file_payload_by_key(
    *,
    output_file_path_by_key: Mapping[str, Path],
) -> dict[str, dict[str, str]]:
    return {
        key: {
            "file_name": file_path.name,
            "path": str(file_path),
            "sha256": sha256_of_file(input_file_path=file_path),
        }
        for key, file_path in output_file_path_by_key.items()
    }


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    source_qc_snapshot_directory: Path,
    source_qc_manifest_file_path: Path,
    source_qc_manifest_payload: Mapping[str, object],
    source_labelled_records_file_path: Path,
    output_file_path_by_key: Mapping[str, Path],
    source_record_count: int,
) -> dict[str, object]:
    source_qc_snapshot_format_version = source_qc_manifest_payload.get(
        "snapshot_format_version"
    )
    if (
        not isinstance(source_qc_snapshot_format_version, str)
        or not source_qc_snapshot_format_version.strip()
    ):
        raise RuntimeError(
            "Source pAgo QC manifest must define a non-empty "
            "snapshot_format_version."
        )

    return {
        "artifact_type": "pago_qc_filtered_datasets",
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_qc_snapshot_directory": str(source_qc_snapshot_directory),
        "source_qc_artifact_type": source_qc_manifest_payload.get("artifact_type"),
        "source_qc_manifest_path": str(source_qc_manifest_file_path),
        "source_qc_manifest_sha256": sha256_of_file(
            input_file_path=source_qc_manifest_file_path,
        ),
        "source_qc_snapshot_format_version": source_qc_snapshot_format_version.strip(),
        "source_labelled_records_path": str(source_labelled_records_file_path),
        "source_labelled_records_sha256": sha256_of_file(
            input_file_path=source_labelled_records_file_path,
        ),
        "source_record_count": int(source_record_count),
        "filter_policy_sha256": build_filter_policy_sha256(),
        "filtered_dataset_definitions": build_filter_policy_payload(),
        "excluded_records_semantics": EXCLUDED_RECORDS_SEMANTICS,
        "output_files": _build_output_file_payload_by_key(
            output_file_path_by_key=output_file_path_by_key,
        ),
    }


def _resolve_immutable_source_qc_snapshot_directory(
    *,
    source_qc_snapshot_directory: Path,
    source_qc_manifest_payload: Mapping[str, object],
) -> Path:
    immutable_snapshot_relative_path = source_qc_manifest_payload.get(
        "immutable_snapshot_relative_path"
    )
    if not isinstance(immutable_snapshot_relative_path, str):
        return source_qc_snapshot_directory

    source_qc_snapshot_root_directory = source_qc_snapshot_directory.parent
    if source_qc_snapshot_directory.name == "latest":
        source_qc_snapshot_root_directory = source_qc_snapshot_directory.parent
    elif source_qc_snapshot_directory.parent.name == "snapshots":
        source_qc_snapshot_root_directory = source_qc_snapshot_directory.parent.parent

    immutable_snapshot_directory = (
        source_qc_snapshot_root_directory / immutable_snapshot_relative_path
    )
    if immutable_snapshot_directory.exists():
        return immutable_snapshot_directory

    return source_qc_snapshot_directory


def save_pago_qc_filtered_datasets(
    *,
    source_qc_snapshot_directory: PathLike,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> PagoQcFilteredDatasetsResult:
    requested_source_qc_snapshot_directory = _as_path(source_qc_snapshot_directory)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    source_qc_payload = load_pago_qc_evidence_inventory_by_directory(
        snapshot_directory=requested_source_qc_snapshot_directory,
    )
    source_qc_manifest_payload = source_qc_payload["manifest"]
    if not isinstance(source_qc_manifest_payload, Mapping):
        raise RuntimeError("Loaded source QC payload must include a manifest object.")

    resolved_source_qc_snapshot_directory = (
        _resolve_immutable_source_qc_snapshot_directory(
            source_qc_snapshot_directory=requested_source_qc_snapshot_directory,
            source_qc_manifest_payload=source_qc_manifest_payload,
        )
    )
    if resolved_source_qc_snapshot_directory != requested_source_qc_snapshot_directory:
        source_qc_payload = load_pago_qc_evidence_inventory_by_directory(
            snapshot_directory=resolved_source_qc_snapshot_directory,
        )
        source_qc_manifest_payload = source_qc_payload["manifest"]
        if not isinstance(source_qc_manifest_payload, Mapping):
            raise RuntimeError(
                "Loaded immutable source QC payload must include a manifest object."
            )

    source_labelled_records_dataframe = source_qc_payload["labelled_records"]
    if not isinstance(source_labelled_records_dataframe, pd.DataFrame):
        raise RuntimeError("Loaded source QC payload must include labelled_records.")

    filtered_dataset_by_name = build_pago_qc_filtered_datasets(
        labelled_records_dataframe=source_labelled_records_dataframe,
    )
    filtered_dataset_counts_dataframe = build_filtered_dataset_counts_dataframe(
        filtered_dataset_by_name=filtered_dataset_by_name,
        total_record_count=int(len(source_labelled_records_dataframe)),
    )

    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query="pago_qc_filtered_datasets",
    )
    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    immutable_snapshot_relative_path = str(Path("snapshots") / snapshot_directory_name)

    output_file_path_by_key = {
        "classic_pago_high_precision_records_file": (
            immutable_snapshot_directory
            / DATASET_DEFINITION_BY_NAME[
                CLASSIC_PAGO_HIGH_PRECISION_DATASET_NAME
            ].file_name
        ),
        "classic_pago_review_records_file": (
            immutable_snapshot_directory
            / DATASET_DEFINITION_BY_NAME[CLASSIC_PAGO_REVIEW_DATASET_NAME].file_name
        ),
        "piwi_re_records_file": (
            immutable_snapshot_directory
            / DATASET_DEFINITION_BY_NAME[PIWI_RE_DATASET_NAME].file_name
        ),
        "excluded_records_file": (
            immutable_snapshot_directory
            / DATASET_DEFINITION_BY_NAME[EXCLUDED_DATASET_NAME].file_name
        ),
        "filtered_dataset_counts_file": (
            immutable_snapshot_directory / DEFAULT_FILTERED_DATASET_COUNTS_FILE_NAME
        ),
    }

    _write_dataframe_csv_atomic(
        dataframe=filtered_dataset_by_name["classic_pago_high_precision"],
        output_file_path=output_file_path_by_key[
            "classic_pago_high_precision_records_file"
        ],
    )
    _write_dataframe_csv_atomic(
        dataframe=filtered_dataset_by_name["classic_pago_review"],
        output_file_path=output_file_path_by_key["classic_pago_review_records_file"],
    )
    _write_dataframe_csv_atomic(
        dataframe=filtered_dataset_by_name["piwi_re"],
        output_file_path=output_file_path_by_key["piwi_re_records_file"],
    )
    _write_dataframe_csv_atomic(
        dataframe=filtered_dataset_by_name[EXCLUDED_DATASET_NAME],
        output_file_path=output_file_path_by_key["excluded_records_file"],
    )
    _write_dataframe_csv_atomic(
        dataframe=filtered_dataset_counts_dataframe,
        output_file_path=output_file_path_by_key["filtered_dataset_counts_file"],
    )

    manifest_file_path = immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = _build_manifest_payload(
        snapshot_created_at_utc=snapshot_created_at_utc,
        immutable_snapshot_directory_name=snapshot_directory_name,
        immutable_snapshot_relative_path=immutable_snapshot_relative_path,
        source_qc_snapshot_directory=resolved_source_qc_snapshot_directory,
        source_qc_manifest_file_path=source_qc_payload["manifest_file_path"],
        source_qc_manifest_payload=source_qc_manifest_payload,
        source_labelled_records_file_path=source_qc_payload[
            "labelled_records_file_path"
        ],
        output_file_path_by_key=output_file_path_by_key,
        source_record_count=int(len(source_labelled_records_dataframe)),
    )
    write_json_atomic(payload=manifest_payload, output_file_path=manifest_file_path)

    if update_latest_directory:
        _replace_latest_directory(
            latest_directory=resolved_snapshot_root_directory / "latest",
            files_to_copy=[
                (
                    output_file_path_by_key[
                        "classic_pago_high_precision_records_file"
                    ],
                    DATASET_DEFINITION_BY_NAME[
                        CLASSIC_PAGO_HIGH_PRECISION_DATASET_NAME
                    ].file_name,
                ),
                (
                    output_file_path_by_key["classic_pago_review_records_file"],
                    DATASET_DEFINITION_BY_NAME[
                        CLASSIC_PAGO_REVIEW_DATASET_NAME
                    ].file_name,
                ),
                (
                    output_file_path_by_key["piwi_re_records_file"],
                    DATASET_DEFINITION_BY_NAME[PIWI_RE_DATASET_NAME].file_name,
                ),
                (
                    output_file_path_by_key["excluded_records_file"],
                    DATASET_DEFINITION_BY_NAME[EXCLUDED_DATASET_NAME].file_name,
                ),
                (
                    output_file_path_by_key["filtered_dataset_counts_file"],
                    DEFAULT_FILTERED_DATASET_COUNTS_FILE_NAME,
                ),
                (manifest_file_path, DEFAULT_MANIFEST_FILE_NAME),
            ],
        )

    return PagoQcFilteredDatasetsResult(
        source_qc_snapshot_directory=resolved_source_qc_snapshot_directory,
        source_qc_manifest_file_path=source_qc_payload["manifest_file_path"],
        source_labelled_records_file_path=source_qc_payload[
            "labelled_records_file_path"
        ],
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        classic_pago_high_precision_records_file_path=output_file_path_by_key[
            "classic_pago_high_precision_records_file"
        ],
        classic_pago_review_records_file_path=output_file_path_by_key[
            "classic_pago_review_records_file"
        ],
        piwi_re_records_file_path=output_file_path_by_key["piwi_re_records_file"],
        excluded_records_file_path=output_file_path_by_key["excluded_records_file"],
        filtered_dataset_counts_file_path=output_file_path_by_key[
            "filtered_dataset_counts_file"
        ],
        source_record_count=int(len(source_labelled_records_dataframe)),
    )


def _get_manifest_output_file_path(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
    output_file_key: str,
) -> Path:
    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError(
            "Saved pAgo QC filtered datasets manifest must define output_files."
        )

    file_payload = output_files.get(output_file_key)
    if not isinstance(file_payload, Mapping):
        raise RuntimeError(
            "Saved pAgo QC filtered datasets manifest is missing "
            f"{output_file_key!r}."
        )

    file_name = file_payload.get("file_name")
    if not isinstance(file_name, str) or not file_name.strip():
        raise RuntimeError(
            "Saved pAgo QC filtered datasets manifest output entry must define "
            f"a non-empty file_name for {output_file_key!r}."
        )

    resolved_file_path = snapshot_directory / file_name
    if not resolved_file_path.exists():
        raise FileNotFoundError(
            "Saved pAgo QC filtered datasets file not found: "
            f"{resolved_file_path}."
        )
    return resolved_file_path


def _validate_loaded_pago_qc_filtered_datasets_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Mapping[str, object],
) -> dict[str, Path]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "pago_qc_filtered_datasets":
        raise RuntimeError(
            "Saved pAgo QC filtered datasets manifest artifact_type mismatch. "
            f"Expected 'pago_qc_filtered_datasets', got {artifact_type!r}."
        )
    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved pAgo QC filtered datasets snapshot_format_version mismatch. "
            f"Expected {SNAPSHOT_FORMAT_VERSION!r}, got "
            f"{snapshot_format_version!r}."
        )
    filter_policy_sha256 = manifest_payload.get("filter_policy_sha256")
    expected_filter_policy_sha256 = build_filter_policy_sha256()
    if filter_policy_sha256 != expected_filter_policy_sha256:
        raise RuntimeError(
            "Saved pAgo QC filtered datasets filter_policy_sha256 mismatch. "
            f"Expected {expected_filter_policy_sha256!r}, got "
            f"{filter_policy_sha256!r}."
        )

    required_output_file_keys = (
        "classic_pago_high_precision_records_file",
        "classic_pago_review_records_file",
        "piwi_re_records_file",
        "excluded_records_file",
        "filtered_dataset_counts_file",
    )
    resolved_file_path_by_key: dict[str, Path] = {}
    for output_file_key in required_output_file_keys:
        resolved_file_path_by_key[output_file_key] = _get_manifest_output_file_path(
            snapshot_directory=resolved_snapshot_directory,
            manifest_payload=manifest_payload,
            output_file_key=output_file_key,
        )

    output_files = manifest_payload.get("output_files")
    if isinstance(output_files, Mapping):
        for output_file_key, resolved_file_path in resolved_file_path_by_key.items():
            file_payload = output_files.get(output_file_key)
            if not isinstance(file_payload, Mapping):
                continue
            expected_sha256 = file_payload.get("sha256")
            if not isinstance(expected_sha256, str) or not expected_sha256.strip():
                raise RuntimeError(
                    "Saved pAgo QC filtered datasets manifest output entry "
                    f"must define a non-empty sha256 for {output_file_key!r}."
                )
            actual_sha256 = sha256_of_file(input_file_path=resolved_file_path)
            if actual_sha256 != expected_sha256.strip():
                raise RuntimeError(
                    "Saved pAgo QC filtered datasets hash mismatch for "
                    f"{output_file_key}. Expected {expected_sha256}, got "
                    f"{actual_sha256}."
                )

    return resolved_file_path_by_key


def load_pago_qc_filtered_datasets_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved pAgo QC filtered datasets manifest must deserialize into a "
            "JSON object."
        )

    resolved_file_path_by_key = _validate_loaded_pago_qc_filtered_datasets_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "classic_pago_high_precision_records_file_path": resolved_file_path_by_key[
            "classic_pago_high_precision_records_file"
        ],
        "classic_pago_review_records_file_path": resolved_file_path_by_key[
            "classic_pago_review_records_file"
        ],
        "piwi_re_records_file_path": resolved_file_path_by_key[
            "piwi_re_records_file"
        ],
        "excluded_records_file_path": resolved_file_path_by_key[
            "excluded_records_file"
        ],
        "filtered_dataset_counts_file_path": resolved_file_path_by_key[
            "filtered_dataset_counts_file"
        ],
        "classic_pago_high_precision_records": pd.read_csv(
            resolved_file_path_by_key["classic_pago_high_precision_records_file"],
            low_memory=False,
            dtype={"protein_uid": "string"},
        ),
        "classic_pago_review_records": pd.read_csv(
            resolved_file_path_by_key["classic_pago_review_records_file"],
            low_memory=False,
            dtype={"protein_uid": "string"},
        ),
        "piwi_re_records": pd.read_csv(
            resolved_file_path_by_key["piwi_re_records_file"],
            low_memory=False,
            dtype={"protein_uid": "string"},
        ),
        "excluded_records": pd.read_csv(
            resolved_file_path_by_key["excluded_records_file"],
            low_memory=False,
            dtype={"protein_uid": "string"},
        ),
        "filtered_dataset_counts": pd.read_csv(
            resolved_file_path_by_key["filtered_dataset_counts_file"],
        ),
    }


def load_latest_pago_qc_filtered_datasets(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_pago_qc_filtered_datasets_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest",
    )


def _source_qc_snapshot_identity_matches(
    *,
    manifest_payload: Mapping[str, object],
    source_qc_payload: Mapping[str, object],
) -> bool:
    expected_manifest_sha256 = manifest_payload.get("source_qc_manifest_sha256")
    expected_labelled_records_sha256 = manifest_payload.get(
        "source_labelled_records_sha256"
    )
    expected_format_version = manifest_payload.get("source_qc_snapshot_format_version")

    if not isinstance(expected_manifest_sha256, str) or not expected_manifest_sha256:
        return False
    if (
        not isinstance(expected_labelled_records_sha256, str)
        or not expected_labelled_records_sha256
    ):
        return False

    source_manifest_file_path = source_qc_payload.get("manifest_file_path")
    source_labelled_records_file_path = source_qc_payload.get(
        "labelled_records_file_path"
    )
    source_manifest_payload = source_qc_payload.get("manifest")
    if not isinstance(source_manifest_file_path, Path):
        return False
    if not isinstance(source_labelled_records_file_path, Path):
        return False
    if not isinstance(source_manifest_payload, Mapping):
        return False

    actual_manifest_sha256 = sha256_of_file(
        input_file_path=source_manifest_file_path,
    )
    actual_labelled_records_sha256 = sha256_of_file(
        input_file_path=source_labelled_records_file_path,
    )

    return (
        expected_manifest_sha256 == actual_manifest_sha256
        and expected_labelled_records_sha256 == actual_labelled_records_sha256
        and expected_format_version
        == source_manifest_payload.get("snapshot_format_version")
    )


def latest_pago_qc_filtered_datasets_is_available(
    *,
    snapshot_root_directory: PathLike,
    source_qc_snapshot_root_directory: PathLike,
) -> bool:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    manifest_file_path = latest_directory / DEFAULT_MANIFEST_FILE_NAME
    if not latest_directory.exists() or not manifest_file_path.exists():
        return False

    try:
        manifest_payload = read_json_file(input_file_path=manifest_file_path)
        if not isinstance(manifest_payload, Mapping):
            return False
        _validate_loaded_pago_qc_filtered_datasets_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        source_qc_payload = load_latest_pago_qc_evidence_inventory(
            snapshot_root_directory=source_qc_snapshot_root_directory,
        )
        if not _source_qc_snapshot_identity_matches(
            manifest_payload=manifest_payload,
            source_qc_payload=source_qc_payload,
        ):
            return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False

    return True


def resolve_pago_qc_filtered_datasets(
    *,
    snapshot_mode: SnapshotMode | str,
    source_qc_snapshot_root_directory: PathLike,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_qc_snapshot_root_directory = _as_path(
        source_qc_snapshot_root_directory
    )

    latest_is_available = latest_pago_qc_filtered_datasets_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_qc_snapshot_root_directory=resolved_source_qc_snapshot_root_directory,
    )
    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest pAgo QC filtered datasets snapshot directory was found "
                "for the current source QC snapshot. Run the workflow once with "
                "snapshot_mode='create_new' before using 'reuse_latest'."
            )
        return load_latest_pago_qc_filtered_datasets(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print(
            "Latest pAgo QC filtered datasets snapshot is available. "
            "Reusing frozen snapshot."
        )
        return load_latest_pago_qc_filtered_datasets(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: 'create_new', 'reuse_latest', "
            "'reuse_latest_or_create'."
        )

    source_qc_payload = load_latest_pago_qc_evidence_inventory(
        snapshot_root_directory=resolved_source_qc_snapshot_root_directory,
    )
    saved_snapshot = save_pago_qc_filtered_datasets(
        source_qc_snapshot_directory=source_qc_payload["snapshot_directory"],
        snapshot_root_directory=resolved_snapshot_root_directory,
        update_latest_directory=update_latest_directory,
    )
    return load_pago_qc_filtered_datasets_by_directory(
        snapshot_directory=saved_snapshot.snapshot_directory,
    )


def run_pago_qc_filtered_datasets(
    *,
    source_qc_snapshot_directory: PathLike,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> PagoQcFilteredDatasetsResult:
    return save_pago_qc_filtered_datasets(
        source_qc_snapshot_directory=source_qc_snapshot_directory,
        snapshot_root_directory=snapshot_root_directory,
        update_latest_directory=update_latest_directory,
    )


def list_saved_pago_qc_filtered_datasets_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory,
    )


def get_most_recent_pago_qc_filtered_datasets_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory,
    )
