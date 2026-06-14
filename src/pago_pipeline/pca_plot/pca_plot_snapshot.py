from __future__ import annotations

import shutil
import tempfile
import time
import webbrowser
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, TypeAlias

import pandas as pd
import plotly.graph_objects as go
import psutil

from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.pca_plot.pca_plot_dataset import (
    DEFAULT_PLOT_DIMENSION_MODE,
    DEFAULT_PLOT_MIRROR_X_AXIS,
    DEFAULT_PLOT_ROTATION_DEGREES,
    build_pca_plot_dataframe,
)
from src.pago_pipeline.pca_plot.pca_plot_render import (
    DEFAULT_SCHEME_NAME_TO_COLUMN_NAME,
    build_pca_plot_figure,
    build_pca_plot_title,
)
from src.pago_pipeline.pca_kmeans_snapshot import (
    load_latest_pca_kmeans_snapshot,
    load_pca_kmeans_snapshot_by_directory,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

DEFAULT_PLOT_DATA_FILE_NAME = "plot_data.csv"
DEFAULT_PLOT_HTML_FILE_NAME = "interactive_plot.html"
DEFAULT_PROFILING_LOG_FILE_NAME = "profiling_log.csv"
VALID_PLOT_SNAPSHOT_ARTIFACT_TYPES = frozenset(
    {"pca_plot_snapshot", "pca_3d_plot_snapshot"}
)


@dataclass(frozen=True)
class PcaPlotSnapshotResult:
    snapshot_directory: Path
    manifest_file_path: Path
    plot_data_file_path: Path
    plot_html_file_path: Path
    profiling_log_file_path: Path
    source_pca_kmeans_snapshot_directory: Path


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


def _write_plotly_html_atomic(
    *,
    figure_object: go.Figure,
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
        suffix=".html",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)
    figure_object.write_html(
        str(resolved_temporary_file_path),
        include_plotlyjs="cdn",
        full_html=True,
    )
    resolved_temporary_file_path.replace(resolved_output_file_path)
    return resolved_output_file_path


def _resolve_rendered_plot_dimension_count(plot_dataframe: pd.DataFrame) -> int:
    if "plot_dimension_count" in plot_dataframe.columns and not plot_dataframe.empty:
        unique_dimension_count_list = (
            pd.Series(plot_dataframe["plot_dimension_count"])
            .dropna()
            .astype(int)
            .unique()
            .tolist()
        )
        if len(unique_dimension_count_list) == 1:
            resolved_plot_dimension_count = int(unique_dimension_count_list[0])
            if resolved_plot_dimension_count in {2, 3}:
                return resolved_plot_dimension_count
        if len(unique_dimension_count_list) > 1:
            raise RuntimeError(
                "Plot dataframe contains inconsistent plot_dimension_count values: "
                f"{unique_dimension_count_list}."
            )

    if "plot_pc3" in plot_dataframe.columns:
        return 3
    return 2


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    plot_data_file_path: Path,
    plot_html_file_path: Path,
    profiling_log_file_path: Path,
    source_pca_kmeans_snapshot_payload: Mapping[str, object],
    profiling_rows: list[dict[str, Any]],
    plot_dimension_mode: str,
    rendered_plot_dimension_count: int,
    plot_rotation_degrees: float,
    plot_mirror_x_axis: bool,
    open_html_in_browser: bool,
) -> dict[str, object]:
    source_pca_kmeans_manifest = source_pca_kmeans_snapshot_payload.get("manifest")
    if not isinstance(source_pca_kmeans_manifest, Mapping):
        raise RuntimeError(
            "Resolved PCA/KMeans snapshot payload is missing a valid source manifest."
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
        "artifact_type": "pca_plot_snapshot",
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": "manifest.json",
        "plot_data_file_name": plot_data_file_path.name,
        "plot_html_file_name": plot_html_file_path.name,
        "profiling_log_file_name": profiling_log_file_path.name,
        "plot_data_file_sha256": sha256_of_file(input_file_path=plot_data_file_path),
        "plot_html_file_sha256": sha256_of_file(input_file_path=plot_html_file_path),
        "profiling_log_file_sha256": sha256_of_file(
            input_file_path=profiling_log_file_path
        ),
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_pca_kmeans_snapshot_relative_path": source_pca_kmeans_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_pca_kmeans_snapshot_directory_name": source_pca_kmeans_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_pca_kmeans_manifest_sha256": sha256_of_file(
            input_file_path=Path(source_pca_kmeans_snapshot_payload["manifest_file_path"])
        ),
        "selected_pca_component_count": source_pca_kmeans_manifest.get(
            "selected_pca_component_count"
        ),
        "selected_cluster_count_k": source_pca_kmeans_manifest.get(
            "selected_cluster_count_k"
        ),
        "selected_variance_explained_fraction": source_pca_kmeans_manifest.get(
            "selected_variance_explained_fraction"
        ),
        "selected_silhouette_best_sampled": source_pca_kmeans_manifest.get(
            "selected_silhouette_best_sampled"
        ),
        "final_sampled_silhouette_value": source_pca_kmeans_manifest.get(
            "final_sampled_silhouette_value"
        ),
        "selection_reason": source_pca_kmeans_manifest.get("selection_reason"),
        "rendered_plot_dimension_count": int(rendered_plot_dimension_count),
        "plot_configuration": {
            "plot_dimension_mode": str(plot_dimension_mode),
            "plot_rotation_degrees": float(plot_rotation_degrees),
            "plot_mirror_x_axis": bool(plot_mirror_x_axis),
            "open_html_in_browser": bool(open_html_in_browser),
            "color_schemes": list(DEFAULT_SCHEME_NAME_TO_COLUMN_NAME.keys()),
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


def save_pca_plot_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_pca_kmeans_snapshot_directory: PathLike,
    plot_dimension_mode: str = DEFAULT_PLOT_DIMENSION_MODE,
    plot_rotation_degrees: float = DEFAULT_PLOT_ROTATION_DEGREES,
    plot_mirror_x_axis: bool = DEFAULT_PLOT_MIRROR_X_AXIS,
    open_html_in_browser: bool = False,
    update_latest_directory: bool = True,
) -> PcaPlotSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_pca_kmeans_snapshot_directory = _as_path(
        source_pca_kmeans_snapshot_directory
    )

    source_pca_kmeans_snapshot_payload = load_pca_kmeans_snapshot_by_directory(
        snapshot_directory=resolved_source_pca_kmeans_snapshot_directory
    )
    cluster_assignments_dataframe = source_pca_kmeans_snapshot_payload.get(
        "cluster_assignments"
    )
    explained_variance_ratio = source_pca_kmeans_snapshot_payload.get(
        "explained_variance_ratio"
    )
    source_manifest = source_pca_kmeans_snapshot_payload.get("manifest")
    if not isinstance(cluster_assignments_dataframe, pd.DataFrame):
        raise RuntimeError(
            "Resolved PCA/KMeans snapshot payload is missing a valid cluster assignments DataFrame."
        )
    if explained_variance_ratio is None:
        raise RuntimeError(
            "Resolved PCA/KMeans snapshot payload is missing explained variance ratios."
        )
    if not isinstance(source_manifest, Mapping):
        raise RuntimeError(
            "Resolved PCA/KMeans snapshot payload is missing a valid source manifest."
        )

    process = psutil.Process()
    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_started_at_utc = _current_utc_timestamp()
    snapshot_elapsed_started_at = time.perf_counter()
    snapshot_cpu_seconds_before = _get_process_cpu_seconds(process)
    snapshot_rss_memory_before_mb = _get_process_rss_memory_mb(process)

    search_query = str(
        source_manifest.get("search_query")
        or source_manifest.get("translated_query")
        or "pca_plot"
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

        plot_data_started_at_utc = _current_utc_timestamp()
        plot_data_elapsed_started_at = time.perf_counter()
        plot_data_cpu_seconds_before = _get_process_cpu_seconds(process)
        plot_data_rss_memory_before_mb = _get_process_rss_memory_mb(process)
        plot_dataframe = build_pca_plot_dataframe(
            cluster_assignments_dataframe=cluster_assignments_dataframe,
            plot_dimension_mode=plot_dimension_mode,
            rotation_degrees=plot_rotation_degrees,
            mirror_x_axis=plot_mirror_x_axis,
        )
        rendered_plot_dimension_count = _resolve_rendered_plot_dimension_count(
            plot_dataframe
        )
        _profile_step(
            work_unit_kind="build_plot_dataframe",
            started_at_utc=plot_data_started_at_utc,
            elapsed_started_at=plot_data_elapsed_started_at,
            cpu_seconds_before=plot_data_cpu_seconds_before,
            rss_memory_before_mb=plot_data_rss_memory_before_mb,
            extra_fields={"row_count": int(len(plot_dataframe))},
        )

        variance_explained_ratio_array = pd.Series(explained_variance_ratio, dtype=float)
        figure_title = build_pca_plot_title(
            selected_pca_component_count=int(
                source_manifest["selected_pca_component_count"]
            ),
            selected_cluster_count_k=int(source_manifest["selected_cluster_count_k"]),
            rendered_plot_dimension_count=rendered_plot_dimension_count,
            variance_explained_rendered_fraction=float(
                variance_explained_ratio_array.iloc[
                    :rendered_plot_dimension_count
                ].sum()
            ),
            variance_explained_total_fraction=float(
                variance_explained_ratio_array.sum()
            ),
            final_sampled_silhouette_value=float(
                source_manifest["final_sampled_silhouette_value"]
            ),
            selection_reason=str(source_manifest["selection_reason"]),
        )

        figure_started_at_utc = _current_utc_timestamp()
        figure_elapsed_started_at = time.perf_counter()
        figure_cpu_seconds_before = _get_process_cpu_seconds(process)
        figure_rss_memory_before_mb = _get_process_rss_memory_mb(process)
        figure_object = build_pca_plot_figure(
            plot_dataframe=plot_dataframe,
            figure_title=figure_title,
        )
        _profile_step(
            work_unit_kind="build_plot_figure",
            started_at_utc=figure_started_at_utc,
            elapsed_started_at=figure_elapsed_started_at,
            cpu_seconds_before=figure_cpu_seconds_before,
            rss_memory_before_mb=figure_rss_memory_before_mb,
        )

        immutable_snapshot_relative_path = str(
            Path("snapshots") / snapshot_directory_name
        )
        plot_data_file_path = immutable_snapshot_directory / DEFAULT_PLOT_DATA_FILE_NAME
        plot_html_file_path = immutable_snapshot_directory / DEFAULT_PLOT_HTML_FILE_NAME
        profiling_log_file_path = (
            immutable_snapshot_directory / DEFAULT_PROFILING_LOG_FILE_NAME
        )
        manifest_file_path = immutable_snapshot_directory / "manifest.json"

        _write_dataframe_csv_atomic(
            dataframe=plot_dataframe,
            output_file_path=plot_data_file_path,
        )
        _write_plotly_html_atomic(
            figure_object=figure_object,
            output_file_path=plot_html_file_path,
        )
        if open_html_in_browser:
            webbrowser.open_new_tab(plot_html_file_path.resolve().as_uri())

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

        manifest_payload = _build_manifest_payload(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            plot_data_file_path=plot_data_file_path,
            plot_html_file_path=plot_html_file_path,
            profiling_log_file_path=profiling_log_file_path,
            source_pca_kmeans_snapshot_payload=source_pca_kmeans_snapshot_payload,
            profiling_rows=profiling_rows,
            plot_dimension_mode=plot_dimension_mode,
            rendered_plot_dimension_count=rendered_plot_dimension_count,
            plot_rotation_degrees=plot_rotation_degrees,
            plot_mirror_x_axis=plot_mirror_x_axis,
            open_html_in_browser=open_html_in_browser,
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
                    (plot_data_file_path, plot_data_file_path.name),
                    (plot_html_file_path, plot_html_file_path.name),
                    (profiling_log_file_path, profiling_log_file_path.name),
                    (manifest_file_path, manifest_file_path.name),
                ],
            )

    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return PcaPlotSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        manifest_file_path=manifest_file_path,
        plot_data_file_path=plot_data_file_path,
        plot_html_file_path=plot_html_file_path,
        profiling_log_file_path=profiling_log_file_path,
        source_pca_kmeans_snapshot_directory=resolved_source_pca_kmeans_snapshot_directory,
    )


def _validate_loaded_pca_plot_snapshot_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Mapping[str, object],
) -> tuple[Path, Path, Path]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type not in VALID_PLOT_SNAPSHOT_ARTIFACT_TYPES:
        raise RuntimeError(
            "Saved PCA plot snapshot manifest artifact_type mismatch. "
            f"Expected one of {sorted(VALID_PLOT_SNAPSHOT_ARTIFACT_TYPES)}, "
            f"got {artifact_type!r}."
        )

    required_file_name_keys = (
        "plot_data_file_name",
        "plot_html_file_name",
        "profiling_log_file_name",
    )
    resolved_file_paths: dict[str, Path] = {}
    for key in required_file_name_keys:
        file_name = manifest_payload.get(key)
        if not isinstance(file_name, str) or not file_name.strip():
            raise RuntimeError(
                f"Saved PCA plot snapshot manifest must define a non-empty {key}."
            )
        resolved_file_path = resolved_snapshot_directory / file_name
        if not resolved_file_path.exists():
            raise FileNotFoundError(
                f"Saved PCA plot snapshot file not found: {resolved_file_path}."
            )
        resolved_file_paths[key] = resolved_file_path

    hash_key_map = {
        "plot_data_file_name": "plot_data_file_sha256",
        "plot_html_file_name": "plot_html_file_sha256",
        "profiling_log_file_name": "profiling_log_file_sha256",
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
                f"Saved PCA plot snapshot hash mismatch for {file_name_key}. "
                f"Expected {expected_sha256}, got {actual_sha256}."
            )

    return (
        resolved_file_paths["plot_data_file_name"],
        resolved_file_paths["plot_html_file_name"],
        resolved_file_paths["profiling_log_file_name"],
    )


def load_pca_plot_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved PCA plot snapshot manifest must deserialize into a JSON object."
        )

    (
        plot_data_file_path,
        plot_html_file_path,
        profiling_log_file_path,
    ) = _validate_loaded_pca_plot_snapshot_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "plot_data_file_path": plot_data_file_path,
        "plot_html_file_path": plot_html_file_path,
        "profiling_log_file_path": profiling_log_file_path,
        "plot_data": pd.read_csv(plot_data_file_path, low_memory=False),
        "profiling_log": pd.read_csv(profiling_log_file_path),
    }


def load_latest_pca_plot_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_pca_plot_snapshot_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest"
    )


def latest_pca_plot_snapshot_is_available(
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
        _validate_loaded_pca_plot_snapshot_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False

    return True


def _pca_plot_snapshot_matches_source_pca_kmeans(
    *,
    pca_plot_snapshot_payload: dict[str, object],
    source_pca_kmeans_snapshot_payload: dict[str, object],
) -> bool:
    pca_plot_manifest = pca_plot_snapshot_payload.get("manifest")
    source_manifest_file_path = source_pca_kmeans_snapshot_payload.get(
        "manifest_file_path"
    )
    if not isinstance(pca_plot_manifest, Mapping) or not isinstance(
        source_manifest_file_path,
        (str, Path),
    ):
        return False

    return (
        pca_plot_manifest.get("source_pca_kmeans_manifest_sha256")
        == sha256_of_file(input_file_path=Path(source_manifest_file_path))
    )


def resolve_pca_plot_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_pca_kmeans_snapshot_root_directory: PathLike,
    plot_dimension_mode: str = DEFAULT_PLOT_DIMENSION_MODE,
    plot_rotation_degrees: float = DEFAULT_PLOT_ROTATION_DEGREES,
    plot_mirror_x_axis: bool = DEFAULT_PLOT_MIRROR_X_AXIS,
    open_html_in_browser: bool = False,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_pca_kmeans_snapshot_root_directory = _as_path(
        source_pca_kmeans_snapshot_root_directory
    )

    latest_is_available = latest_pca_plot_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory
    )
    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest PCA plot snapshot directory was found. Run the workflow "
                "once with snapshot_mode='create_new' before using 'reuse_latest'."
            )
        return load_latest_pca_plot_snapshot(
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
        source_pca_kmeans_snapshot_payload = load_latest_pca_kmeans_snapshot(
            snapshot_root_directory=resolved_source_pca_kmeans_snapshot_root_directory
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No reusable PCA/KMeans snapshot was found for the PCA plot workflow. "
            "Create the upstream PCA/KMeans snapshot before running this step."
        ) from exc

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        latest_snapshot_payload = load_latest_pca_plot_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )
        if _pca_plot_snapshot_matches_source_pca_kmeans(
            pca_plot_snapshot_payload=latest_snapshot_payload,
            source_pca_kmeans_snapshot_payload=source_pca_kmeans_snapshot_payload,
        ):
            print("Latest PCA plot snapshot is available. Reusing frozen snapshot.")
            return latest_snapshot_payload

        print(
            "Latest PCA plot snapshot is available but does not match the latest "
            "source PCA/KMeans snapshot. Creating a new frozen snapshot."
        )

    saved_snapshot = save_pca_plot_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_pca_kmeans_snapshot_directory=Path(
            source_pca_kmeans_snapshot_payload["snapshot_directory"]
        ),
        plot_dimension_mode=plot_dimension_mode,
        plot_rotation_degrees=plot_rotation_degrees,
        plot_mirror_x_axis=plot_mirror_x_axis,
        open_html_in_browser=open_html_in_browser,
        update_latest_directory=update_latest_directory,
    )
    return load_pca_plot_snapshot_by_directory(
        snapshot_directory=saved_snapshot.snapshot_directory
    )


def list_saved_pca_plot_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(snapshot_root_directory=snapshot_root_directory)


def get_most_recent_pca_plot_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
