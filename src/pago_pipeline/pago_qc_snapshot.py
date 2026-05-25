from __future__ import annotations

import shutil
import tempfile
from collections.abc import Mapping, Sequence
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
from src.pago_pipeline.pago_qc import (
    DEPRECATED_ALIAS_BY_COLUMN,
    FLAG_COLUMNS,
    build_filter_decision_counts_dataframe,
    build_pago_qc_evidence_flags,
    build_pago_qc_labelled_records,
    build_evidence_counts_dataframe,
    build_label_counts_dataframe,
    build_suspicious_terms_dataframe,
    build_top_value_count_dataframe,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

DEFAULT_EVIDENCE_FLAGS_FILE_NAME = "evidence_flags.csv"
DEFAULT_EVIDENCE_COUNTS_FILE_NAME = "evidence_counts.csv"
DEFAULT_LABELLED_RECORDS_FILE_NAME = "labelled_records.csv"
DEFAULT_LABEL_COUNTS_FILE_NAME = "label_counts.csv"
DEFAULT_FILTER_DECISION_COUNTS_FILE_NAME = "filter_decision_counts.csv"
DEFAULT_TOP_REGION_NAMES_FILE_NAME = "top_region_names.csv"
DEFAULT_TOP_PRODUCTS_FILE_NAME = "top_products.csv"
DEFAULT_SUSPICIOUS_TERMS_FILE_NAME = "suspicious_terms.csv"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
SNAPSHOT_FORMAT_VERSION = "2.0"
LABEL_COLUMNS: tuple[str, ...] = (
    "primary_label",
    "qc_decision",
    "confidence_score",
    "rationale",
)
LENGTH_BIN_COLUMN = "length_bin"


@dataclass(frozen=True)
class PagoQcEvidenceInventoryResult:
    metadata_csv_file_path: Path
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    evidence_flags_file_path: Path
    evidence_counts_file_path: Path
    labelled_records_file_path: Path
    label_counts_file_path: Path
    filter_decision_counts_file_path: Path
    top_region_names_file_path: Path
    top_products_file_path: Path
    suspicious_terms_file_path: Path
    metadata_row_count: int
    metadata_column_count: int


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


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    metadata_csv_file_path: Path,
    output_file_path_by_key: Mapping[str, Path],
    metadata_dataframe: pd.DataFrame,
    flag_columns: Sequence[str],
) -> dict[str, object]:
    output_files = {
        key: {
            "file_name": file_path.name,
            "path": str(file_path),
            "sha256": sha256_of_file(input_file_path=file_path),
        }
        for key, file_path in output_file_path_by_key.items()
    }

    manifest_payload: dict[str, object] = {
        "artifact_type": "pago_qc_evidence_inventory",
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "metadata_csv_file_path": str(metadata_csv_file_path),
        "metadata_csv_sha256": sha256_of_file(input_file_path=metadata_csv_file_path),
        "metadata_row_count": int(len(metadata_dataframe)),
        "metadata_column_count": int(len(metadata_dataframe.columns)),
        "flag_columns": list(flag_columns),
        "label_columns": list(LABEL_COLUMNS),
        "length_bin_column": LENGTH_BIN_COLUMN,
        "deprecated_aliases": dict(DEPRECATED_ALIAS_BY_COLUMN),
        "output_files": output_files,
    }

    return manifest_payload


def save_pago_qc_evidence_inventory(
    *,
    metadata_csv_file_path: PathLike,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> PagoQcEvidenceInventoryResult:
    resolved_metadata_csv_file_path = _as_path(metadata_csv_file_path)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query="pago_qc_evidence_inventory",
    )
    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    immutable_snapshot_relative_path = str(Path("snapshots") / snapshot_directory_name)
    immutable_snapshot_complete = False

    try:
        metadata_dataframe = pd.read_csv(
            resolved_metadata_csv_file_path,
            low_memory=False,
            dtype={"protein_uid": "string"},
        )
        evidence_dataframe = build_pago_qc_evidence_flags(
            metadata_dataframe=metadata_dataframe,
        )
        labelled_records_dataframe = build_pago_qc_labelled_records(
            evidence_dataframe=evidence_dataframe,
        )
        evidence_counts_dataframe = build_evidence_counts_dataframe(
            evidence_dataframe=evidence_dataframe,
        )
        label_counts_dataframe = build_label_counts_dataframe(
            labelled_dataframe=labelled_records_dataframe,
        )
        filter_decision_counts_dataframe = build_filter_decision_counts_dataframe(
            labelled_dataframe=labelled_records_dataframe,
        )
        top_region_names_dataframe = build_top_value_count_dataframe(
            dataframe=metadata_dataframe,
            column_name="feature__region__qual__region_name",
        )
        top_products_dataframe = build_top_value_count_dataframe(
            dataframe=metadata_dataframe,
            column_name="feature__protein__qual__product",
        )
        suspicious_terms_dataframe = build_suspicious_terms_dataframe(
            evidence_dataframe=evidence_dataframe,
        )

        output_file_path_by_key = {
            "evidence_flags_file": (
                immutable_snapshot_directory / DEFAULT_EVIDENCE_FLAGS_FILE_NAME
            ),
            "evidence_counts_file": (
                immutable_snapshot_directory / DEFAULT_EVIDENCE_COUNTS_FILE_NAME
            ),
            "labelled_records_file": (
                immutable_snapshot_directory / DEFAULT_LABELLED_RECORDS_FILE_NAME
            ),
            "label_counts_file": (
                immutable_snapshot_directory / DEFAULT_LABEL_COUNTS_FILE_NAME
            ),
            "filter_decision_counts_file": (
                immutable_snapshot_directory / DEFAULT_FILTER_DECISION_COUNTS_FILE_NAME
            ),
            "top_region_names_file": (
                immutable_snapshot_directory / DEFAULT_TOP_REGION_NAMES_FILE_NAME
            ),
            "top_products_file": (
                immutable_snapshot_directory / DEFAULT_TOP_PRODUCTS_FILE_NAME
            ),
            "suspicious_terms_file": (
                immutable_snapshot_directory / DEFAULT_SUSPICIOUS_TERMS_FILE_NAME
            ),
        }

        _write_dataframe_csv_atomic(
            dataframe=evidence_dataframe,
            output_file_path=output_file_path_by_key["evidence_flags_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=evidence_counts_dataframe,
            output_file_path=output_file_path_by_key["evidence_counts_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=labelled_records_dataframe,
            output_file_path=output_file_path_by_key["labelled_records_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=label_counts_dataframe,
            output_file_path=output_file_path_by_key["label_counts_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=filter_decision_counts_dataframe,
            output_file_path=output_file_path_by_key["filter_decision_counts_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=top_region_names_dataframe,
            output_file_path=output_file_path_by_key["top_region_names_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=top_products_dataframe,
            output_file_path=output_file_path_by_key["top_products_file"],
        )
        _write_dataframe_csv_atomic(
            dataframe=suspicious_terms_dataframe,
            output_file_path=output_file_path_by_key["suspicious_terms_file"],
        )

        manifest_file_path = immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
        manifest_payload = _build_manifest_payload(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            metadata_csv_file_path=resolved_metadata_csv_file_path,
            output_file_path_by_key=output_file_path_by_key,
            metadata_dataframe=metadata_dataframe,
            flag_columns=FLAG_COLUMNS,
        )
        write_json_atomic(payload=manifest_payload, output_file_path=manifest_file_path)
        immutable_snapshot_complete = True

        if update_latest_directory:
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=[
                    (
                        output_file_path_by_key["evidence_flags_file"],
                        DEFAULT_EVIDENCE_FLAGS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["evidence_counts_file"],
                        DEFAULT_EVIDENCE_COUNTS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["labelled_records_file"],
                        DEFAULT_LABELLED_RECORDS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["label_counts_file"],
                        DEFAULT_LABEL_COUNTS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["filter_decision_counts_file"],
                        DEFAULT_FILTER_DECISION_COUNTS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["top_region_names_file"],
                        DEFAULT_TOP_REGION_NAMES_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["top_products_file"],
                        DEFAULT_TOP_PRODUCTS_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["suspicious_terms_file"],
                        DEFAULT_SUSPICIOUS_TERMS_FILE_NAME,
                    ),
                    (manifest_file_path, DEFAULT_MANIFEST_FILE_NAME),
                ],
            )
    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return PagoQcEvidenceInventoryResult(
        metadata_csv_file_path=resolved_metadata_csv_file_path,
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        evidence_flags_file_path=output_file_path_by_key["evidence_flags_file"],
        evidence_counts_file_path=output_file_path_by_key["evidence_counts_file"],
        labelled_records_file_path=output_file_path_by_key["labelled_records_file"],
        label_counts_file_path=output_file_path_by_key["label_counts_file"],
        filter_decision_counts_file_path=output_file_path_by_key[
            "filter_decision_counts_file"
        ],
        top_region_names_file_path=output_file_path_by_key["top_region_names_file"],
        top_products_file_path=output_file_path_by_key["top_products_file"],
        suspicious_terms_file_path=output_file_path_by_key["suspicious_terms_file"],
        metadata_row_count=int(len(metadata_dataframe)),
        metadata_column_count=int(len(metadata_dataframe.columns)),
    )


def _get_manifest_output_file_path(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
    output_file_key: str,
    required: bool = True,
) -> Optional[Path]:
    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError(
            "Saved pAgo QC evidence inventory manifest must define output_files."
        )

    file_payload = output_files.get(output_file_key)
    if file_payload is None:
        if required:
            raise RuntimeError(
                "Saved pAgo QC evidence inventory manifest is missing "
                f"{output_file_key!r}."
            )
        return None
    if not isinstance(file_payload, Mapping):
        raise RuntimeError(
            "Saved pAgo QC evidence inventory manifest output entry must be an "
            f"object for {output_file_key!r}."
        )

    file_name = file_payload.get("file_name")
    if not isinstance(file_name, str) or not file_name.strip():
        raise RuntimeError(
            "Saved pAgo QC evidence inventory manifest output entry must define "
            f"a non-empty file_name for {output_file_key!r}."
        )

    resolved_file_path = snapshot_directory / file_name
    if not resolved_file_path.exists():
        raise FileNotFoundError(
            "Saved pAgo QC evidence inventory file not found: "
            f"{resolved_file_path}."
        )
    return resolved_file_path


def _validate_loaded_pago_qc_evidence_inventory_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Mapping[str, object],
) -> dict[str, Optional[Path]]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "pago_qc_evidence_inventory":
        raise RuntimeError(
            "Saved pAgo QC evidence inventory manifest artifact_type mismatch. "
            f"Expected 'pago_qc_evidence_inventory', got {artifact_type!r}."
        )
    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved pAgo QC evidence inventory snapshot_format_version mismatch. "
            f"Expected {SNAPSHOT_FORMAT_VERSION!r}, got "
            f"{snapshot_format_version!r}."
        )

    required_output_file_keys = (
        "evidence_flags_file",
        "evidence_counts_file",
        "labelled_records_file",
        "label_counts_file",
        "filter_decision_counts_file",
        "top_region_names_file",
        "top_products_file",
        "suspicious_terms_file",
    )
    resolved_file_path_by_key: dict[str, Optional[Path]] = {}
    for output_file_key in required_output_file_keys:
        resolved_file_path_by_key[output_file_key] = _get_manifest_output_file_path(
            snapshot_directory=resolved_snapshot_directory,
            manifest_payload=manifest_payload,
            output_file_key=output_file_key,
        )

    output_files = manifest_payload.get("output_files")
    if isinstance(output_files, Mapping):
        for output_file_key, resolved_file_path in resolved_file_path_by_key.items():
            if resolved_file_path is None:
                continue
            file_payload = output_files.get(output_file_key)
            if not isinstance(file_payload, Mapping):
                continue
            expected_sha256 = file_payload.get("sha256")
            if not isinstance(expected_sha256, str) or not expected_sha256.strip():
                raise RuntimeError(
                    "Saved pAgo QC evidence inventory manifest output entry "
                    f"must define a non-empty sha256 for {output_file_key!r}."
                )
            expected_sha256 = expected_sha256.strip()
            actual_sha256 = sha256_of_file(input_file_path=resolved_file_path)
            if actual_sha256 != expected_sha256:
                raise RuntimeError(
                    "Saved pAgo QC evidence inventory hash mismatch for "
                    f"{output_file_key}. Expected {expected_sha256}, got "
                    f"{actual_sha256}."
                )

    return resolved_file_path_by_key


def load_pago_qc_evidence_inventory_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved pAgo QC evidence inventory manifest must deserialize into a "
            "JSON object."
        )

    resolved_file_path_by_key = _validate_loaded_pago_qc_evidence_inventory_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "evidence_flags_file_path": resolved_file_path_by_key["evidence_flags_file"],
        "evidence_counts_file_path": resolved_file_path_by_key["evidence_counts_file"],
        "labelled_records_file_path": resolved_file_path_by_key[
            "labelled_records_file"
        ],
        "label_counts_file_path": resolved_file_path_by_key["label_counts_file"],
        "filter_decision_counts_file_path": resolved_file_path_by_key[
            "filter_decision_counts_file"
        ],
        "top_region_names_file_path": resolved_file_path_by_key[
            "top_region_names_file"
        ],
        "top_products_file_path": resolved_file_path_by_key["top_products_file"],
        "suspicious_terms_file_path": resolved_file_path_by_key[
            "suspicious_terms_file"
        ],
        "evidence_flags": pd.read_csv(
            resolved_file_path_by_key["evidence_flags_file"],
            low_memory=False,
            dtype={"protein_uid": "string"},
        ),
        "evidence_counts": pd.read_csv(
            resolved_file_path_by_key["evidence_counts_file"]
        ),
        "labelled_records": pd.read_csv(
            resolved_file_path_by_key["labelled_records_file"],
            low_memory=False,
            dtype={"protein_uid": "string"},
        ),
        "label_counts": pd.read_csv(resolved_file_path_by_key["label_counts_file"]),
        "filter_decision_counts": pd.read_csv(
            resolved_file_path_by_key["filter_decision_counts_file"]
        ),
        "top_region_names": pd.read_csv(
            resolved_file_path_by_key["top_region_names_file"]
        ),
        "top_products": pd.read_csv(resolved_file_path_by_key["top_products_file"]),
        "suspicious_terms": pd.read_csv(
            resolved_file_path_by_key["suspicious_terms_file"]
        ),
    }


def load_latest_pago_qc_evidence_inventory(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_pago_qc_evidence_inventory_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest",
    )


def latest_pago_qc_evidence_inventory_is_available(
    *,
    snapshot_root_directory: PathLike,
    metadata_csv_file_path: Optional[PathLike] = None,
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
        _validate_loaded_pago_qc_evidence_inventory_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        if metadata_csv_file_path is not None:
            expected_metadata_csv_sha256 = manifest_payload.get("metadata_csv_sha256")
            if (
                not isinstance(expected_metadata_csv_sha256, str)
                or not expected_metadata_csv_sha256.strip()
            ):
                return False
            actual_metadata_csv_sha256 = sha256_of_file(
                input_file_path=metadata_csv_file_path
            )
            if expected_metadata_csv_sha256.strip() != actual_metadata_csv_sha256:
                return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False

    return True


def resolve_pago_qc_evidence_inventory(
    *,
    snapshot_mode: SnapshotMode | str,
    metadata_csv_file_path: PathLike,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    latest_is_available = latest_pago_qc_evidence_inventory_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        metadata_csv_file_path=metadata_csv_file_path,
    )
    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest pAgo QC evidence inventory snapshot directory was found. "
                "Run the workflow once with snapshot_mode='create_new' before using "
                "'reuse_latest'."
            )
        return load_latest_pago_qc_evidence_inventory(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print(
            "Latest pAgo QC evidence inventory snapshot is available. "
            "Reusing frozen snapshot."
        )
        return load_latest_pago_qc_evidence_inventory(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: 'create_new', 'reuse_latest', "
            "'reuse_latest_or_create'."
        )

    saved_snapshot = save_pago_qc_evidence_inventory(
        metadata_csv_file_path=metadata_csv_file_path,
        snapshot_root_directory=resolved_snapshot_root_directory,
        update_latest_directory=update_latest_directory,
    )
    return load_pago_qc_evidence_inventory_by_directory(
        snapshot_directory=saved_snapshot.snapshot_directory
    )


def run_pago_qc_evidence_inventory(
    *,
    metadata_csv_file_path: PathLike,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> PagoQcEvidenceInventoryResult:
    return save_pago_qc_evidence_inventory(
        metadata_csv_file_path=metadata_csv_file_path,
        snapshot_root_directory=snapshot_root_directory,
        update_latest_directory=update_latest_directory,
    )


def list_saved_pago_qc_evidence_inventory_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(snapshot_root_directory=snapshot_root_directory)


def get_most_recent_pago_qc_evidence_inventory_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
