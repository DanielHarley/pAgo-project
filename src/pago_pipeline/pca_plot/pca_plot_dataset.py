from __future__ import annotations

import re
from collections.abc import Sequence

import numpy as np
import pandas as pd

DEFAULT_PLOT_DIMENSION_MODE = "auto"
VALID_PLOT_DIMENSION_MODES = frozenset({"auto", "2d", "3d"})
DEFAULT_PLOT_ROTATION_DEGREES = 225.0
DEFAULT_PLOT_MIRROR_X_AXIS = True
UNKNOWN_LABEL = "Unknown"
PAGO_REGEX = re.compile(r"argonaute|piwi|pago", flags=re.IGNORECASE)

DEFAULT_DESCRIPTION_SOURCE_COLUMNS: tuple[str, ...] = (
    "description",
    "gbseq__definition",
    "feature__protein__qual__product",
    "feature__protein__qual__name",
)
DEFAULT_ORGANISM_SOURCE_COLUMNS: tuple[str, ...] = (
    "gbseq__organism",
    "feature__source__qual__organism",
)
DEFAULT_PAGO_SOURCE_COLUMNS: tuple[str, ...] = (
    "description",
    "gbseq__definition",
    "feature__protein__qual__product",
    "feature__protein__qual__name",
    "feature__protein__qual__note",
    "feature__region__qual__note",
    "feature__region__qual__region_name",
)


def _coalesce_first_nonempty_text_value(
    *,
    row: pd.Series,
    candidate_columns: Sequence[str],
    default_value: str = UNKNOWN_LABEL,
) -> str:
    for column_name in candidate_columns:
        if column_name not in row.index:
            continue
        value = row.get(column_name)
        if pd.isna(value):
            continue
        normalized_value = str(value).strip()
        if normalized_value:
            return normalized_value
    return default_value


def _build_text_mask_from_columns(
    *,
    dataframe: pd.DataFrame,
    candidate_columns: Sequence[str],
    regex_pattern: re.Pattern[str],
) -> pd.Series:
    result_mask = pd.Series(False, index=dataframe.index)
    for column_name in candidate_columns:
        if column_name not in dataframe.columns:
            continue
        current_column_mask = (
            dataframe[column_name]
            .fillna("")
            .astype(str)
            .str.contains(regex_pattern, regex=True, na=False)
        )
        result_mask = result_mask | current_column_mask
    return result_mask


def _coerce_plot_dimension_mode(plot_dimension_mode: str) -> str:
    resolved_plot_dimension_mode = str(plot_dimension_mode).strip().lower()
    if resolved_plot_dimension_mode not in VALID_PLOT_DIMENSION_MODES:
        raise ValueError(
            "plot_dimension_mode must be one of "
            f"{sorted(VALID_PLOT_DIMENSION_MODES)}, got {plot_dimension_mode!r}."
        )
    return resolved_plot_dimension_mode


def _resolve_plot_dimension_count(
    *,
    cluster_assignments_dataframe: pd.DataFrame,
    plot_dimension_mode: str,
) -> int:
    resolved_plot_dimension_mode = _coerce_plot_dimension_mode(plot_dimension_mode)
    available_coordinate_columns = {
        column_name
        for column_name in cluster_assignments_dataframe.columns
        if column_name in {"pc1", "pc2", "pc3"}
    }

    if resolved_plot_dimension_mode == "3d":
        required_columns = {"pc1", "pc2", "pc3"}
        missing_required_columns = sorted(
            required_columns.difference(available_coordinate_columns)
        )
        if missing_required_columns:
            raise RuntimeError(
                "3D plotting requires PCA coordinates pc1, pc2, and pc3. "
                f"Missing columns: {missing_required_columns}."
            )
        return 3

    required_2d_columns = {"pc1", "pc2"}
    missing_required_2d_columns = sorted(
        required_2d_columns.difference(available_coordinate_columns)
    )
    if missing_required_2d_columns:
        raise RuntimeError(
            "PCA plotting requires at least pc1 and pc2 coordinates. "
            f"Missing columns: {missing_required_2d_columns}."
        )

    if resolved_plot_dimension_mode == "2d":
        return 2
    if "pc3" in available_coordinate_columns:
        return 3
    return 2


def transform_pca_coordinates_for_plot(
    *,
    coordinate_matrix: np.ndarray,
    rotation_degrees: float = DEFAULT_PLOT_ROTATION_DEGREES,
    mirror_x_axis: bool = DEFAULT_PLOT_MIRROR_X_AXIS,
) -> np.ndarray:
    resolved_coordinate_matrix = np.asarray(coordinate_matrix, dtype=np.float64)
    if resolved_coordinate_matrix.ndim != 2 or resolved_coordinate_matrix.shape[1] not in {
        2,
        3,
    }:
        raise ValueError(
            "coordinate_matrix must be a two-dimensional array with 2 or 3 columns."
        )

    coordinate_dimension_count = resolved_coordinate_matrix.shape[1]

    rotation_angle_radians = np.deg2rad(float(rotation_degrees))
    rotation_matrix = np.eye(coordinate_dimension_count, dtype=np.float64)
    rotation_matrix[:2, :2] = np.array(
        [
            [np.cos(rotation_angle_radians), -np.sin(rotation_angle_radians)],
            [np.sin(rotation_angle_radians), np.cos(rotation_angle_radians)],
        ],
        dtype=np.float64,
    )
    transformed_coordinate_matrix = resolved_coordinate_matrix @ rotation_matrix

    if mirror_x_axis:
        mirror_matrix = np.eye(coordinate_dimension_count, dtype=np.float64)
        mirror_matrix[0, 0] = -1.0
        transformed_coordinate_matrix = transformed_coordinate_matrix @ mirror_matrix

    return transformed_coordinate_matrix


def build_pca_plot_dataframe(
    *,
    cluster_assignments_dataframe: pd.DataFrame,
    plot_dimension_mode: str = DEFAULT_PLOT_DIMENSION_MODE,
    rotation_degrees: float = DEFAULT_PLOT_ROTATION_DEGREES,
    mirror_x_axis: bool = DEFAULT_PLOT_MIRROR_X_AXIS,
) -> pd.DataFrame:
    required_columns = ("cluster_label", "pc1", "pc2")
    missing_required_columns = [
        column_name
        for column_name in required_columns
        if column_name not in cluster_assignments_dataframe.columns
    ]
    if missing_required_columns:
        raise RuntimeError(
            "Cluster assignments dataframe is missing required columns for PCA plotting: "
            f"{missing_required_columns}."
        )

    plot_dimension_count = _resolve_plot_dimension_count(
        cluster_assignments_dataframe=cluster_assignments_dataframe,
        plot_dimension_mode=plot_dimension_mode,
    )
    coordinate_column_names = ["pc1", "pc2"]
    if plot_dimension_count == 3:
        coordinate_column_names.append("pc3")

    dataframe_copy = cluster_assignments_dataframe.copy()
    transformed_coordinate_matrix = transform_pca_coordinates_for_plot(
        coordinate_matrix=dataframe_copy[coordinate_column_names].to_numpy(),
        rotation_degrees=rotation_degrees,
        mirror_x_axis=mirror_x_axis,
    )

    dataframe_copy["plot_pc1"] = transformed_coordinate_matrix[:, 0]
    dataframe_copy["plot_pc2"] = transformed_coordinate_matrix[:, 1]
    if plot_dimension_count == 3:
        dataframe_copy["plot_pc3"] = transformed_coordinate_matrix[:, 2]
    dataframe_copy["plot_dimension_count"] = plot_dimension_count
    dataframe_copy["plot_dimension_label"] = f"{plot_dimension_count}d"
    dataframe_copy["description_label"] = dataframe_copy.apply(
        lambda row: _coalesce_first_nonempty_text_value(
            row=row,
            candidate_columns=DEFAULT_DESCRIPTION_SOURCE_COLUMNS,
            default_value="Protein record",
        ),
        axis=1,
    )
    dataframe_copy["organism_label"] = dataframe_copy.apply(
        lambda row: _coalesce_first_nonempty_text_value(
            row=row,
            candidate_columns=DEFAULT_ORGANISM_SOURCE_COLUMNS,
        ),
        axis=1,
    )

    taxonomy_column_to_label_map = {
        "taxonomy__03": "phylum_label",
        "taxonomy__04": "class_label",
        "taxonomy__06": "family_label",
        "taxonomy__07": "genus_label",
    }
    for source_column_name, target_column_name in taxonomy_column_to_label_map.items():
        if source_column_name in dataframe_copy.columns:
            dataframe_copy[target_column_name] = (
                dataframe_copy[source_column_name]
                .fillna(UNKNOWN_LABEL)
                .astype(str)
                .replace("", UNKNOWN_LABEL)
            )
        else:
            dataframe_copy[target_column_name] = UNKNOWN_LABEL

    dataframe_copy["is_pago"] = _build_text_mask_from_columns(
        dataframe=dataframe_copy,
        candidate_columns=DEFAULT_PAGO_SOURCE_COLUMNS,
        regex_pattern=PAGO_REGEX,
    )
    dataframe_copy["pago_label"] = np.where(
        dataframe_copy["is_pago"],
        "pAgo",
        "Other",
    )

    return dataframe_copy
