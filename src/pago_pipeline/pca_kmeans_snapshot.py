from __future__ import annotations

import re
import shutil
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, TypeAlias

import numpy as np
import pandas as pd
import psutil

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
from src.pago_pipeline.pca_kmeans import (
    DEFAULT_KMEANS_INITIALIZATION_REPEAT_COUNT,
    DEFAULT_KMEANS_N_INIT,
    DEFAULT_MINIMUM_ACCEPTABLE_INIT_ARI_MIN,
    DEFAULT_MINIMUM_ACCEPTABLE_SUBSAMPLE_ARI_MIN,
    DEFAULT_PCA_RANDOM_STATE,
    DEFAULT_PCA_SVD_SOLVER,
    DEFAULT_SILHOUETTE_RANDOM_STATE,
    DEFAULT_SILHOUETTE_SAMPLE_SIZE,
    DEFAULT_SUBSAMPLE_FRACTION,
    DEFAULT_SUBSAMPLE_RANDOM_STATE,
    DEFAULT_SUBSAMPLE_REPEAT_COUNT,
    build_cluster_assignment_dataframe,
    choose_best_m_and_k_from_stability_grid,
    compute_pca_coordinates_and_model,
    fit_kmeans_and_compute_sampled_silhouette,
    get_effective_pca_component_count_grid,
    run_pca_kmeans_stability_grid_sweep,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic
from src.pago_pipeline.sweep_genes_snapshot import (
    load_latest_sweep_genes_snapshot,
    load_sweep_genes_snapshot_by_directory,
)

PathLike: TypeAlias = str | Path

DEFAULT_PCA_COMPONENT_COUNT_GRID: tuple[int, ...] = (10, 20, 50, 100, 200)
DEFAULT_KMEANS_CLUSTER_COUNT_GRID: tuple[int, ...] = tuple(range(2, 11))
DEFAULT_EXPORT_PROJECTION_COMPONENT_COUNT = 3
DEFAULT_CLUSTER_ASSIGNMENTS_FILE_NAME = "cluster_assignments.csv"
DEFAULT_STABILITY_GRID_FILE_NAME = "stability_grid.csv"
DEFAULT_PROFILING_LOG_FILE_NAME = "profiling_log.csv"
DEFAULT_ALIGNMENT_REPORT_FILE_NAME = "alignment_report.json"

PROTEIN_UID_PATTERN = re.compile(r"(?:^|\|)protein_uid=([^|]+)")


@dataclass(frozen=True)
class PcaKMeansSnapshotResult:
    snapshot_directory: Path
    manifest_file_path: Path
    pca_coordinates_file_path: Path
    explained_variance_ratio_file_path: Path
    cluster_assignments_file_path: Path
    stability_grid_file_path: Path
    profiling_log_file_path: Path
    alignment_report_file_path: Path
    source_sweep_snapshot_directory: Path
    source_metadata_snapshot_directory: Path


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _get_process_cpu_seconds(process: psutil.Process) -> float:
    cpu_times = process.cpu_times()
    return float(cpu_times.user + cpu_times.system)


def _get_process_rss_memory_mb(process: psutil.Process) -> float:
    return float(process.memory_info().rss / (1024**2))


def _finalize_profiling_row(*, row_payload: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: (round(value, 6) if isinstance(value, float) else value)
        for key, value in row_payload.items()
    }


def _build_numpy_file_name(prefix: str, component_count: int) -> str:
    if component_count <= 0:
        raise ValueError("component_count must be a positive integer.")
    return f"{prefix}_{component_count}D.npy"


def _extract_protein_uid_from_record_id(record_id: object) -> Optional[str]:
    if record_id is None:
        return None
    match = PROTEIN_UID_PATTERN.search(str(record_id))
    if match is None:
        return None
    protein_uid = str(match.group(1)).strip()
    return protein_uid or None


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


def _write_numpy_array_atomic(
    *,
    array: np.ndarray,
    output_file_path: PathLike,
) -> Path:
    resolved_output_file_path = _as_path(output_file_path)
    resolved_output_file_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb",
        delete=False,
        dir=resolved_output_file_path.parent,
        suffix=".npy",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)
    np.save(resolved_temporary_file_path, array)
    resolved_temporary_file_path.replace(resolved_output_file_path)
    return resolved_output_file_path


def _align_sequence_metadata_with_ncbi_metadata(
    *,
    sequence_metadata_dataframe: pd.DataFrame,
    metadata_dataframe: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    sequence_dataframe = sequence_metadata_dataframe.copy()
    metadata_dataframe_copy = metadata_dataframe.copy()

    if "record_id" not in sequence_dataframe.columns:
        raise RuntimeError(
            "The SWeeP sequence metadata must contain a 'record_id' column."
        )
    if "sequence_index" not in sequence_dataframe.columns:
        raise RuntimeError(
            "The SWeeP sequence metadata must contain a 'sequence_index' column."
        )
    if "protein_uid" not in metadata_dataframe_copy.columns:
        raise RuntimeError(
            "The NCBI metadata snapshot must contain a 'protein_uid' column."
        )

    sequence_dataframe["protein_uid"] = sequence_dataframe["record_id"].map(
        _extract_protein_uid_from_record_id
    )
    metadata_dataframe_copy["protein_uid"] = (
        metadata_dataframe_copy["protein_uid"].astype(str).str.strip()
    )

    missing_protein_uid_in_sequence_count = int(
        sequence_dataframe["protein_uid"].isna().sum()
    )
    duplicate_protein_uid_in_sequence_count = int(
        sequence_dataframe["protein_uid"].duplicated().sum()
    )
    duplicate_protein_uid_in_metadata_count = int(
        metadata_dataframe_copy["protein_uid"].duplicated().sum()
    )

    if missing_protein_uid_in_sequence_count > 0:
        raise RuntimeError(
            "Could not extract protein_uid from all SWeeP sequence metadata rows. "
            f"Missing count: {missing_protein_uid_in_sequence_count}."
        )
    if duplicate_protein_uid_in_sequence_count > 0:
        raise RuntimeError(
            "SWeeP sequence metadata protein_uid values must be unique for one-to-one "
            "alignment. Duplicate count: "
            f"{duplicate_protein_uid_in_sequence_count}."
        )
    if duplicate_protein_uid_in_metadata_count > 0:
        raise RuntimeError(
            "NCBI metadata protein_uid values must be unique for one-to-one alignment. "
            f"Duplicate count: {duplicate_protein_uid_in_metadata_count}."
        )

    aligned_dataframe = sequence_dataframe.merge(
        metadata_dataframe_copy,
        how="left",
        on="protein_uid",
        validate="one_to_one",
        suffixes=("", "__metadata"),
    )
    aligned_dataframe = aligned_dataframe.sort_values("sequence_index").reset_index(
        drop=True
    )

    metadata_columns = [
        column_name
        for column_name in metadata_dataframe_copy.columns
        if column_name != "protein_uid"
    ]
    missing_metadata_row_mask = aligned_dataframe[metadata_columns].isna().all(axis=1)
    matched_metadata_row_count = int((~missing_metadata_row_mask).sum())
    missing_metadata_row_count = int(missing_metadata_row_mask.sum())

    alignment_report = {
        "sequence_row_count": int(len(sequence_dataframe)),
        "metadata_row_count": int(len(metadata_dataframe_copy)),
        "matched_metadata_row_count": matched_metadata_row_count,
        "missing_metadata_row_count": missing_metadata_row_count,
        "missing_protein_uid_in_sequence_count": missing_protein_uid_in_sequence_count,
        "duplicate_protein_uid_in_sequence_count": duplicate_protein_uid_in_sequence_count,
        "duplicate_protein_uid_in_metadata_count": duplicate_protein_uid_in_metadata_count,
    }
    return aligned_dataframe, alignment_report


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    pca_coordinates_file_path: Path,
    explained_variance_ratio_file_path: Path,
    cluster_assignments_file_path: Path,
    stability_grid_file_path: Path,
    profiling_log_file_path: Path,
    alignment_report_file_path: Path,
    source_sweep_snapshot_payload: Mapping[str, object],
    source_metadata_snapshot_payload: Mapping[str, object],
    alignment_report: Mapping[str, Any],
    profiling_rows: Sequence[Mapping[str, Any]],
    selected_configuration_row: Mapping[str, Any],
    selection_reason: str,
    pca_component_count_grid: Sequence[int],
    kmeans_cluster_count_grid: Sequence[int],
    kmeans_initialization_repeat_count: int,
    subsample_repeat_count: int,
    subsample_fraction: float,
    pca_svd_solver: str,
    pca_random_state: int,
    kmeans_n_init: int | str,
    silhouette_sample_size: int,
    silhouette_random_state: int,
    subsample_random_state: int,
    minimum_acceptable_init_ari_min: float,
    minimum_acceptable_subsample_ari_min: float,
    export_projection_component_count: int,
) -> dict[str, object]:
    source_sweep_manifest = source_sweep_snapshot_payload.get("manifest")
    source_metadata_manifest = source_metadata_snapshot_payload.get("manifest")
    if not isinstance(source_sweep_manifest, Mapping):
        raise RuntimeError(
            "Resolved SWeeP snapshot payload is missing a valid source manifest."
        )
    if not isinstance(source_metadata_manifest, Mapping):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing a valid source manifest."
        )

    snapshot_total_row = next(
        (
            row
            for row in profiling_rows
            if row.get("work_unit_kind") == "snapshot_total"
        ),
        None,
    )

    return {
        "snapshot_format_version": "1.0",
        "artifact_type": "pca_kmeans_snapshot",
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": "manifest.json",
        "pca_coordinates_file_name": pca_coordinates_file_path.name,
        "explained_variance_ratio_file_name": explained_variance_ratio_file_path.name,
        "cluster_assignments_file_name": cluster_assignments_file_path.name,
        "stability_grid_file_name": stability_grid_file_path.name,
        "profiling_log_file_name": profiling_log_file_path.name,
        "alignment_report_file_name": alignment_report_file_path.name,
        "pca_coordinates_file_sha256": sha256_of_file(
            input_file_path=pca_coordinates_file_path
        ),
        "explained_variance_ratio_file_sha256": sha256_of_file(
            input_file_path=explained_variance_ratio_file_path
        ),
        "cluster_assignments_file_sha256": sha256_of_file(
            input_file_path=cluster_assignments_file_path
        ),
        "stability_grid_file_sha256": sha256_of_file(
            input_file_path=stability_grid_file_path
        ),
        "profiling_log_file_sha256": sha256_of_file(
            input_file_path=profiling_log_file_path
        ),
        "alignment_report_file_sha256": sha256_of_file(
            input_file_path=alignment_report_file_path
        ),
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_sweep_snapshot_relative_path": source_sweep_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_sweep_snapshot_directory_name": source_sweep_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_sweep_manifest_sha256": sha256_of_file(
            input_file_path=Path(source_sweep_snapshot_payload["manifest_file_path"])
        ),
        "source_metadata_snapshot_relative_path": source_metadata_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_metadata_snapshot_directory_name": source_metadata_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_metadata_manifest_sha256": sha256_of_file(
            input_file_path=Path(source_metadata_snapshot_payload["manifest_file_path"])
        ),
        "sequence_count": int(alignment_report["sequence_row_count"]),
        "selected_pca_component_count": int(
            selected_configuration_row["pca_component_count"]
        ),
        "selected_cluster_count_k": int(selected_configuration_row["k"]),
        "selected_variance_explained_fraction": float(
            selected_configuration_row["variance_explained_fraction"]
        ),
        "selected_silhouette_best_sampled": float(
            selected_configuration_row["silhouette_best_sampled"]
        ),
        "selected_init_ari_min": float(selected_configuration_row["init_ari_min"]),
        "selected_subsample_ari_min": float(
            selected_configuration_row["subsample_ari_min"]
        ),
        "selection_reason": str(selection_reason),
        "alignment_summary": dict(alignment_report),
        "parameter_search_configuration": {
            "pca_component_count_grid": [int(value) for value in pca_component_count_grid],
            "kmeans_cluster_count_grid": [
                int(value) for value in kmeans_cluster_count_grid
            ],
            "kmeans_initialization_repeat_count": int(
                kmeans_initialization_repeat_count
            ),
            "subsample_repeat_count": int(subsample_repeat_count),
            "subsample_fraction": float(subsample_fraction),
            "pca_svd_solver": str(pca_svd_solver),
            "pca_random_state": int(pca_random_state),
            "kmeans_n_init": kmeans_n_init,
            "silhouette_sample_size": int(silhouette_sample_size),
            "silhouette_random_state": int(silhouette_random_state),
            "subsample_random_state": int(subsample_random_state),
            "minimum_acceptable_init_ari_min": float(
                minimum_acceptable_init_ari_min
            ),
            "minimum_acceptable_subsample_ari_min": float(
                minimum_acceptable_subsample_ari_min
            ),
            "export_projection_component_count": int(
                export_projection_component_count
            ),
        },
        "profiling_summary": {
            "row_count": int(len(profiling_rows)),
            "snapshot_total_elapsed_seconds": (
                snapshot_total_row.get("elapsed_seconds")
                if snapshot_total_row is not None
                else None
            ),
            "snapshot_total_cpu_seconds": (
                snapshot_total_row.get("cpu_seconds")
                if snapshot_total_row is not None
                else None
            ),
            "snapshot_total_rss_memory_delta_mb": (
                snapshot_total_row.get("rss_memory_delta_mb")
                if snapshot_total_row is not None
                else None
            ),
        },
    }


def save_pca_kmeans_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_sweep_snapshot_directory: PathLike,
    source_metadata_snapshot_directory: PathLike,
    pca_component_count_grid: Sequence[int] = DEFAULT_PCA_COMPONENT_COUNT_GRID,
    kmeans_cluster_count_grid: Sequence[int] = DEFAULT_KMEANS_CLUSTER_COUNT_GRID,
    pca_svd_solver: str = DEFAULT_PCA_SVD_SOLVER,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
    kmeans_n_init: int | str = DEFAULT_KMEANS_N_INIT,
    silhouette_sample_size: int = DEFAULT_SILHOUETTE_SAMPLE_SIZE,
    silhouette_random_state: int = DEFAULT_SILHOUETTE_RANDOM_STATE,
    kmeans_initialization_repeat_count: int = DEFAULT_KMEANS_INITIALIZATION_REPEAT_COUNT,
    subsample_repeat_count: int = DEFAULT_SUBSAMPLE_REPEAT_COUNT,
    subsample_fraction: float = DEFAULT_SUBSAMPLE_FRACTION,
    subsample_random_state: int = DEFAULT_SUBSAMPLE_RANDOM_STATE,
    minimum_acceptable_init_ari_min: float = DEFAULT_MINIMUM_ACCEPTABLE_INIT_ARI_MIN,
    minimum_acceptable_subsample_ari_min: float = DEFAULT_MINIMUM_ACCEPTABLE_SUBSAMPLE_ARI_MIN,
    export_projection_component_count: int = DEFAULT_EXPORT_PROJECTION_COMPONENT_COUNT,
    update_latest_directory: bool = True,
) -> PcaKMeansSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_sweep_snapshot_directory = _as_path(
        source_sweep_snapshot_directory
    )
    resolved_source_metadata_snapshot_directory = _as_path(
        source_metadata_snapshot_directory
    )

    source_sweep_snapshot_payload = load_sweep_genes_snapshot_by_directory(
        snapshot_directory=resolved_source_sweep_snapshot_directory,
    )
    source_metadata_snapshot_payload = load_metadata_snapshot_by_directory(
        snapshot_directory=resolved_source_metadata_snapshot_directory,
    )

    embeddings = source_sweep_snapshot_payload.get("embeddings")
    sequence_metadata_dataframe = source_sweep_snapshot_payload.get("sequence_metadata")
    metadata_csv_file_path = source_metadata_snapshot_payload.get("csv_file_path")
    if not isinstance(sequence_metadata_dataframe, pd.DataFrame):
        raise RuntimeError(
            "Resolved SWeeP snapshot payload is missing a valid sequence metadata DataFrame."
        )
    if not isinstance(metadata_csv_file_path, Path):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing a valid csv_file_path."
        )
    metadata_dataframe = pd.read_csv(metadata_csv_file_path, low_memory=False)

    if embeddings is None:
        raise RuntimeError(
            "Resolved SWeeP snapshot payload is missing the embeddings matrix."
        )

    process = psutil.Process()
    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_started_at_utc = _current_utc_timestamp()
    snapshot_elapsed_started_at = time.perf_counter()
    snapshot_cpu_seconds_before = _get_process_cpu_seconds(process)
    snapshot_rss_memory_before_mb = _get_process_rss_memory_mb(process)

    source_sweep_manifest = source_sweep_snapshot_payload.get("manifest")
    if not isinstance(source_sweep_manifest, Mapping):
        raise RuntimeError(
            "Resolved SWeeP snapshot payload is missing a valid source manifest."
        )
    search_query = str(
        source_sweep_manifest.get("search_query")
        or source_sweep_manifest.get("translated_query")
        or "pca_kmeans"
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
        profiling_rows: list[dict[str, Any]] = []

        def _profile_step(
            *,
            work_unit_kind: str,
            started_at_utc: str,
            elapsed_started_at: float,
            cpu_seconds_before: float,
            rss_memory_before_mb: float,
            extra_fields: Optional[Mapping[str, Any]] = None,
        ) -> None:
            completed_at_utc = _current_utc_timestamp()
            cpu_seconds_after = _get_process_cpu_seconds(process)
            rss_memory_after_mb = _get_process_rss_memory_mb(process)
            row_payload: dict[str, Any] = {
                "work_unit_kind": work_unit_kind,
                "started_at_utc": started_at_utc,
                "completed_at_utc": completed_at_utc,
                "elapsed_seconds": float(time.perf_counter() - elapsed_started_at),
                "cpu_seconds": float(cpu_seconds_after - cpu_seconds_before),
                "rss_memory_before_mb": rss_memory_before_mb,
                "rss_memory_after_mb": rss_memory_after_mb,
                "rss_memory_delta_mb": float(
                    rss_memory_after_mb - rss_memory_before_mb
                ),
            }
            if extra_fields is not None:
                row_payload.update(dict(extra_fields))
            profiling_rows.append(_finalize_profiling_row(row_payload=row_payload))

        align_started_at_utc = _current_utc_timestamp()
        align_elapsed_started_at = time.perf_counter()
        align_cpu_seconds_before = _get_process_cpu_seconds(process)
        align_rss_memory_before_mb = _get_process_rss_memory_mb(process)
        aligned_metadata_dataframe, alignment_report = (
            _align_sequence_metadata_with_ncbi_metadata(
                sequence_metadata_dataframe=sequence_metadata_dataframe,
                metadata_dataframe=metadata_dataframe,
            )
        )
        if int(alignment_report["missing_metadata_row_count"]) > 0:
            raise RuntimeError(
                "The aligned PCA input still has rows without metadata. Missing row "
                f"count: {alignment_report['missing_metadata_row_count']}."
            )
        _profile_step(
            work_unit_kind="align_inputs",
            started_at_utc=align_started_at_utc,
            elapsed_started_at=align_elapsed_started_at,
            cpu_seconds_before=align_cpu_seconds_before,
            rss_memory_before_mb=align_rss_memory_before_mb,
            extra_fields=alignment_report,
        )

        embeddings_array = np.asarray(embeddings)
        if embeddings_array.shape[0] != len(aligned_metadata_dataframe):
            raise RuntimeError(
                "Embeddings row count does not match aligned metadata row count. "
                f"Embeddings rows: {embeddings_array.shape[0]}; aligned rows: "
                f"{len(aligned_metadata_dataframe)}."
            )

        effective_pca_component_count_grid = get_effective_pca_component_count_grid(
            data_matrix=embeddings_array,
            requested_pca_component_count_grid=pca_component_count_grid,
        )

        grid_started_at_utc = _current_utc_timestamp()
        grid_elapsed_started_at = time.perf_counter()
        grid_cpu_seconds_before = _get_process_cpu_seconds(process)
        grid_rss_memory_before_mb = _get_process_rss_memory_mb(process)
        stability_grid_dataframe, pca_coordinate_cache_by_component_count = (
            run_pca_kmeans_stability_grid_sweep(
                data_matrix=embeddings_array,
                pca_component_count_grid=effective_pca_component_count_grid,
                kmeans_cluster_count_grid=kmeans_cluster_count_grid,
                kmeans_initialization_repeat_count=kmeans_initialization_repeat_count,
                subsample_repeat_count=subsample_repeat_count,
                subsample_fraction=subsample_fraction,
                silhouette_sample_size=silhouette_sample_size,
                pca_svd_solver=pca_svd_solver,
                pca_random_state=pca_random_state,
                kmeans_n_init=kmeans_n_init,
                silhouette_random_state=silhouette_random_state,
                subsample_random_state=subsample_random_state,
                progress_callback=print,
            )
        )
        _profile_step(
            work_unit_kind="grid_sweep",
            started_at_utc=grid_started_at_utc,
            elapsed_started_at=grid_elapsed_started_at,
            cpu_seconds_before=grid_cpu_seconds_before,
            rss_memory_before_mb=grid_rss_memory_before_mb,
            extra_fields={
                "grid_row_count": int(len(stability_grid_dataframe)),
                "effective_pca_component_count_grid": ",".join(
                    str(value) for value in effective_pca_component_count_grid
                ),
            },
        )

        (
            selected_configuration_row,
            selection_reason,
            scored_grid_dataframe,
        ) = choose_best_m_and_k_from_stability_grid(
            stability_grid_dataframe=stability_grid_dataframe,
            minimum_acceptable_init_ari_min=minimum_acceptable_init_ari_min,
            minimum_acceptable_subsample_ari_min=minimum_acceptable_subsample_ari_min,
        )
        selected_pca_component_count = int(
            selected_configuration_row["pca_component_count"]
        )
        selected_cluster_count_k = int(selected_configuration_row["k"])

        final_pca_started_at_utc = _current_utc_timestamp()
        final_pca_elapsed_started_at = time.perf_counter()
        final_pca_cpu_seconds_before = _get_process_cpu_seconds(process)
        final_pca_rss_memory_before_mb = _get_process_rss_memory_mb(process)
        final_pca_coordinate_matrix = pca_coordinate_cache_by_component_count.get(
            selected_pca_component_count
        )
        final_explained_variance_ratio_vector: np.ndarray
        if final_pca_coordinate_matrix is None:
            (
                final_pca_coordinate_matrix,
                final_explained_variance_ratio_vector,
                _,
            ) = compute_pca_coordinates_and_model(
                data_matrix=embeddings_array,
                requested_pca_component_count=selected_pca_component_count,
                pca_svd_solver=pca_svd_solver,
                pca_random_state=pca_random_state,
            )
        else:
            (
                _recomputed_coordinate_matrix,
                final_explained_variance_ratio_vector,
                _,
            ) = compute_pca_coordinates_and_model(
                data_matrix=embeddings_array,
                requested_pca_component_count=selected_pca_component_count,
                pca_svd_solver=pca_svd_solver,
                pca_random_state=pca_random_state,
            )
        _profile_step(
            work_unit_kind="final_pca_fit",
            started_at_utc=final_pca_started_at_utc,
            elapsed_started_at=final_pca_elapsed_started_at,
            cpu_seconds_before=final_pca_cpu_seconds_before,
            rss_memory_before_mb=final_pca_rss_memory_before_mb,
            extra_fields={
                "selected_pca_component_count": selected_pca_component_count,
                "selected_variance_explained_fraction": float(
                    np.sum(final_explained_variance_ratio_vector)
                ),
            },
        )

        final_kmeans_started_at_utc = _current_utc_timestamp()
        final_kmeans_elapsed_started_at = time.perf_counter()
        final_kmeans_cpu_seconds_before = _get_process_cpu_seconds(process)
        final_kmeans_rss_memory_before_mb = _get_process_rss_memory_mb(process)
        (
            cluster_label_vector,
            final_sampled_silhouette_value,
        ) = fit_kmeans_and_compute_sampled_silhouette(
            coordinate_matrix_for_clustering=final_pca_coordinate_matrix,
            requested_cluster_count_k=selected_cluster_count_k,
            kmeans_random_state=42,
            silhouette_sample_size=silhouette_sample_size,
            silhouette_random_state=silhouette_random_state,
            kmeans_n_init=kmeans_n_init,
        )
        _profile_step(
            work_unit_kind="final_kmeans_fit",
            started_at_utc=final_kmeans_started_at_utc,
            elapsed_started_at=final_kmeans_elapsed_started_at,
            cpu_seconds_before=final_kmeans_cpu_seconds_before,
            rss_memory_before_mb=final_kmeans_rss_memory_before_mb,
            extra_fields={
                "selected_cluster_count_k": selected_cluster_count_k,
                "final_sampled_silhouette_value": float(
                    final_sampled_silhouette_value
                ),
            },
        )

        cluster_assignments_dataframe = build_cluster_assignment_dataframe(
            aligned_metadata_dataframe=aligned_metadata_dataframe,
            cluster_label_vector=cluster_label_vector,
            pca_coordinate_matrix=final_pca_coordinate_matrix,
            projection_dimension_for_table=export_projection_component_count,
        )

        immutable_snapshot_relative_path = str(
            Path("snapshots") / snapshot_directory_name
        )
        pca_coordinates_file_path = immutable_snapshot_directory / _build_numpy_file_name(
            "pca_coordinates",
            selected_pca_component_count,
        )
        explained_variance_ratio_file_path = (
            immutable_snapshot_directory
            / _build_numpy_file_name(
                "explained_variance_ratio",
                selected_pca_component_count,
            )
        )
        cluster_assignments_file_path = (
            immutable_snapshot_directory / DEFAULT_CLUSTER_ASSIGNMENTS_FILE_NAME
        )
        stability_grid_file_path = (
            immutable_snapshot_directory / DEFAULT_STABILITY_GRID_FILE_NAME
        )
        profiling_log_file_path = (
            immutable_snapshot_directory / DEFAULT_PROFILING_LOG_FILE_NAME
        )
        alignment_report_file_path = (
            immutable_snapshot_directory / DEFAULT_ALIGNMENT_REPORT_FILE_NAME
        )
        manifest_file_path = immutable_snapshot_directory / "manifest.json"

        _write_numpy_array_atomic(
            array=np.asarray(final_pca_coordinate_matrix, dtype=np.float32),
            output_file_path=pca_coordinates_file_path,
        )
        _write_numpy_array_atomic(
            array=np.asarray(final_explained_variance_ratio_vector, dtype=np.float64),
            output_file_path=explained_variance_ratio_file_path,
        )
        _write_dataframe_csv_atomic(
            dataframe=cluster_assignments_dataframe,
            output_file_path=cluster_assignments_file_path,
        )
        _write_dataframe_csv_atomic(
            dataframe=scored_grid_dataframe,
            output_file_path=stability_grid_file_path,
        )

        snapshot_completed_at_utc = _current_utc_timestamp()
        snapshot_cpu_seconds_after = _get_process_cpu_seconds(process)
        snapshot_rss_memory_after_mb = _get_process_rss_memory_mb(process)
        profiling_rows.append(
            _finalize_profiling_row(
                row_payload={
                    "work_unit_kind": "snapshot_total",
                    "started_at_utc": snapshot_started_at_utc,
                    "completed_at_utc": snapshot_completed_at_utc,
                    "elapsed_seconds": float(
                        time.perf_counter() - snapshot_elapsed_started_at
                    ),
                    "cpu_seconds": float(
                        snapshot_cpu_seconds_after - snapshot_cpu_seconds_before
                    ),
                    "rss_memory_before_mb": snapshot_rss_memory_before_mb,
                    "rss_memory_after_mb": snapshot_rss_memory_after_mb,
                    "rss_memory_delta_mb": float(
                        snapshot_rss_memory_after_mb - snapshot_rss_memory_before_mb
                    ),
                }
            )
        )
        _write_dataframe_csv_atomic(
            dataframe=pd.DataFrame(profiling_rows),
            output_file_path=profiling_log_file_path,
        )
        write_json_atomic(
            payload=dict(alignment_report),
            output_file_path=alignment_report_file_path,
        )

        manifest_payload = _build_manifest_payload(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            pca_coordinates_file_path=pca_coordinates_file_path,
            explained_variance_ratio_file_path=explained_variance_ratio_file_path,
            cluster_assignments_file_path=cluster_assignments_file_path,
            stability_grid_file_path=stability_grid_file_path,
            profiling_log_file_path=profiling_log_file_path,
            alignment_report_file_path=alignment_report_file_path,
            source_sweep_snapshot_payload=source_sweep_snapshot_payload,
            source_metadata_snapshot_payload=source_metadata_snapshot_payload,
            alignment_report=alignment_report,
            profiling_rows=profiling_rows,
            selected_configuration_row=selected_configuration_row.to_dict(),
            selection_reason=selection_reason,
            pca_component_count_grid=effective_pca_component_count_grid,
            kmeans_cluster_count_grid=kmeans_cluster_count_grid,
            kmeans_initialization_repeat_count=kmeans_initialization_repeat_count,
            subsample_repeat_count=subsample_repeat_count,
            subsample_fraction=subsample_fraction,
            pca_svd_solver=pca_svd_solver,
            pca_random_state=pca_random_state,
            kmeans_n_init=kmeans_n_init,
            silhouette_sample_size=silhouette_sample_size,
            silhouette_random_state=silhouette_random_state,
            subsample_random_state=subsample_random_state,
            minimum_acceptable_init_ari_min=minimum_acceptable_init_ari_min,
            minimum_acceptable_subsample_ari_min=minimum_acceptable_subsample_ari_min,
            export_projection_component_count=export_projection_component_count,
        )
        manifest_payload["final_sampled_silhouette_value"] = float(
            final_sampled_silhouette_value
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
                    (pca_coordinates_file_path, pca_coordinates_file_path.name),
                    (
                        explained_variance_ratio_file_path,
                        explained_variance_ratio_file_path.name,
                    ),
                    (
                        cluster_assignments_file_path,
                        cluster_assignments_file_path.name,
                    ),
                    (stability_grid_file_path, stability_grid_file_path.name),
                    (profiling_log_file_path, profiling_log_file_path.name),
                    (alignment_report_file_path, alignment_report_file_path.name),
                    (manifest_file_path, manifest_file_path.name),
                ],
            )

    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return PcaKMeansSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        manifest_file_path=manifest_file_path,
        pca_coordinates_file_path=pca_coordinates_file_path,
        explained_variance_ratio_file_path=explained_variance_ratio_file_path,
        cluster_assignments_file_path=cluster_assignments_file_path,
        stability_grid_file_path=stability_grid_file_path,
        profiling_log_file_path=profiling_log_file_path,
        alignment_report_file_path=alignment_report_file_path,
        source_sweep_snapshot_directory=resolved_source_sweep_snapshot_directory,
        source_metadata_snapshot_directory=resolved_source_metadata_snapshot_directory,
    )


def _validate_loaded_pca_kmeans_snapshot_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Mapping[str, object],
) -> tuple[Path, Path, Path, Path, Path, Path]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "pca_kmeans_snapshot":
        raise RuntimeError(
            "Saved PCA/KMeans snapshot manifest artifact_type mismatch. "
            f"Expected 'pca_kmeans_snapshot', got {artifact_type!r}."
        )

    required_file_name_keys = (
        "pca_coordinates_file_name",
        "explained_variance_ratio_file_name",
        "cluster_assignments_file_name",
        "stability_grid_file_name",
        "profiling_log_file_name",
        "alignment_report_file_name",
    )
    resolved_file_paths: dict[str, Path] = {}
    for key in required_file_name_keys:
        file_name = manifest_payload.get(key)
        if not isinstance(file_name, str) or not file_name.strip():
            raise RuntimeError(
                f"Saved PCA/KMeans snapshot manifest must define a non-empty {key}."
            )
        resolved_file_path = resolved_snapshot_directory / file_name
        if not resolved_file_path.exists():
            raise FileNotFoundError(
                "Saved PCA/KMeans snapshot file not found: "
                f"{resolved_file_path}."
            )
        resolved_file_paths[key] = resolved_file_path

    hash_key_map = {
        "pca_coordinates_file_name": "pca_coordinates_file_sha256",
        "explained_variance_ratio_file_name": "explained_variance_ratio_file_sha256",
        "cluster_assignments_file_name": "cluster_assignments_file_sha256",
        "stability_grid_file_name": "stability_grid_file_sha256",
        "profiling_log_file_name": "profiling_log_file_sha256",
        "alignment_report_file_name": "alignment_report_file_sha256",
    }
    for file_name_key, hash_key in hash_key_map.items():
        expected_sha256 = manifest_payload.get(hash_key)
        if expected_sha256 is None:
            continue
        actual_sha256 = sha256_of_file(
            input_file_path=resolved_file_paths[file_name_key]
        )
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Saved PCA/KMeans snapshot hash mismatch for {file_name_key}. "
                f"Expected {expected_sha256}, got {actual_sha256}."
            )

    return (
        resolved_file_paths["pca_coordinates_file_name"],
        resolved_file_paths["explained_variance_ratio_file_name"],
        resolved_file_paths["cluster_assignments_file_name"],
        resolved_file_paths["stability_grid_file_name"],
        resolved_file_paths["profiling_log_file_name"],
        resolved_file_paths["alignment_report_file_name"],
    )


def load_pca_kmeans_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved PCA/KMeans snapshot manifest must deserialize into a JSON object."
        )

    (
        pca_coordinates_file_path,
        explained_variance_ratio_file_path,
        cluster_assignments_file_path,
        stability_grid_file_path,
        profiling_log_file_path,
        alignment_report_file_path,
    ) = _validate_loaded_pca_kmeans_snapshot_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "pca_coordinates_file_path": pca_coordinates_file_path,
        "explained_variance_ratio_file_path": explained_variance_ratio_file_path,
        "cluster_assignments_file_path": cluster_assignments_file_path,
        "stability_grid_file_path": stability_grid_file_path,
        "profiling_log_file_path": profiling_log_file_path,
        "alignment_report_file_path": alignment_report_file_path,
        "pca_coordinates": np.load(pca_coordinates_file_path, mmap_mode="r"),
        "explained_variance_ratio": np.load(
            explained_variance_ratio_file_path,
            mmap_mode="r",
        ),
        "cluster_assignments": pd.read_csv(
            cluster_assignments_file_path,
            low_memory=False,
        ),
        "stability_grid": pd.read_csv(stability_grid_file_path),
        "profiling_log": pd.read_csv(profiling_log_file_path),
        "alignment_report": read_json_file(input_file_path=alignment_report_file_path),
    }


def load_latest_pca_kmeans_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_pca_kmeans_snapshot_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest",
    )


def latest_pca_kmeans_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    manifest_file_path = latest_directory / "manifest.json"
    if not latest_directory.exists() or not manifest_file_path.exists():
        return False

    try:
        manifest_payload = read_json_file(input_file_path=manifest_file_path)
        if not isinstance(manifest_payload, Mapping):
            return False
        _validate_loaded_pca_kmeans_snapshot_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False

    return True


def resolve_pca_kmeans_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_sweep_snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    pca_component_count_grid: Sequence[int] = DEFAULT_PCA_COMPONENT_COUNT_GRID,
    kmeans_cluster_count_grid: Sequence[int] = DEFAULT_KMEANS_CLUSTER_COUNT_GRID,
    pca_svd_solver: str = DEFAULT_PCA_SVD_SOLVER,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
    kmeans_n_init: int | str = DEFAULT_KMEANS_N_INIT,
    silhouette_sample_size: int = DEFAULT_SILHOUETTE_SAMPLE_SIZE,
    silhouette_random_state: int = DEFAULT_SILHOUETTE_RANDOM_STATE,
    kmeans_initialization_repeat_count: int = DEFAULT_KMEANS_INITIALIZATION_REPEAT_COUNT,
    subsample_repeat_count: int = DEFAULT_SUBSAMPLE_REPEAT_COUNT,
    subsample_fraction: float = DEFAULT_SUBSAMPLE_FRACTION,
    subsample_random_state: int = DEFAULT_SUBSAMPLE_RANDOM_STATE,
    minimum_acceptable_init_ari_min: float = DEFAULT_MINIMUM_ACCEPTABLE_INIT_ARI_MIN,
    minimum_acceptable_subsample_ari_min: float = DEFAULT_MINIMUM_ACCEPTABLE_SUBSAMPLE_ARI_MIN,
    export_projection_component_count: int = DEFAULT_EXPORT_PROJECTION_COMPONENT_COUNT,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_sweep_snapshot_root_directory = _as_path(
        source_sweep_snapshot_root_directory
    )
    resolved_source_metadata_snapshot_root_directory = _as_path(
        source_metadata_snapshot_root_directory
    )

    latest_is_available = latest_pca_kmeans_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory
    )
    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest PCA/KMeans snapshot directory was found. Run the workflow "
                "once with snapshot_mode='create_new' before using 'reuse_latest'."
            )
        return load_latest_pca_kmeans_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print("Latest PCA/KMeans snapshot is available. Reusing frozen snapshot.")
        return load_latest_pca_kmeans_snapshot(
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

    try:
        source_sweep_snapshot_payload = load_latest_sweep_genes_snapshot(
            snapshot_root_directory=resolved_source_sweep_snapshot_root_directory
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No reusable SWeeP Genes snapshot was found for the PCA/KMeans workflow. "
            "Create the upstream SWeeP snapshot before running this step."
        ) from exc

    try:
        source_metadata_snapshot_payload = load_latest_metadata_snapshot(
            snapshot_root_directory=resolved_source_metadata_snapshot_root_directory
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No reusable metadata snapshot was found for the PCA/KMeans workflow. "
            "Create the upstream metadata snapshot before running this step."
        ) from exc

    saved_snapshot = save_pca_kmeans_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_sweep_snapshot_directory=Path(
            source_sweep_snapshot_payload["snapshot_directory"]
        ),
        source_metadata_snapshot_directory=Path(
            source_metadata_snapshot_payload["snapshot_directory"]
        ),
        pca_component_count_grid=pca_component_count_grid,
        kmeans_cluster_count_grid=kmeans_cluster_count_grid,
        pca_svd_solver=pca_svd_solver,
        pca_random_state=pca_random_state,
        kmeans_n_init=kmeans_n_init,
        silhouette_sample_size=silhouette_sample_size,
        silhouette_random_state=silhouette_random_state,
        kmeans_initialization_repeat_count=kmeans_initialization_repeat_count,
        subsample_repeat_count=subsample_repeat_count,
        subsample_fraction=subsample_fraction,
        subsample_random_state=subsample_random_state,
        minimum_acceptable_init_ari_min=minimum_acceptable_init_ari_min,
        minimum_acceptable_subsample_ari_min=minimum_acceptable_subsample_ari_min,
        export_projection_component_count=export_projection_component_count,
        update_latest_directory=update_latest_directory,
    )
    return load_pca_kmeans_snapshot_by_directory(
        snapshot_directory=saved_snapshot.snapshot_directory
    )


def list_saved_pca_kmeans_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(snapshot_root_directory=snapshot_root_directory)


def get_most_recent_pca_kmeans_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
