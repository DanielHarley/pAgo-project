from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
import plotly.colors as plotly_colors
import plotly.graph_objects as go

UNKNOWN_LABEL = "Unknown"
DEFAULT_COLOR_SCALE_NAME = "turbo"
DEFAULT_SCHEME_NAME_TO_COLUMN_NAME: dict[str, str] = {
    "Cluster": "cluster_label",
    "Phylum": "phylum_label",
    "Class": "class_label",
    "Family": "family_label",
    "Genus": "genus_label",
    "pAgo?": "pago_label",
    "No colors": "__all_points__",
}


def _resolve_plot_dimension_count(plot_dataframe: pd.DataFrame) -> int:
    if "plot_dimension_count" in plot_dataframe.columns:
        unique_dimension_count_list = [
            int(value)
            for value in pd.Series(plot_dataframe["plot_dimension_count"])
            .dropna()
            .astype(int)
            .unique()
            .tolist()
        ]
        if len(unique_dimension_count_list) == 1:
            resolved_plot_dimension_count = unique_dimension_count_list[0]
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


def build_distinct_categorical_color_vector(
    *,
    category_series: pd.Series | Sequence[object],
    colorscale_name: str = DEFAULT_COLOR_SCALE_NAME,
) -> list[str]:
    category_as_categorical = (
        pd.Series(category_series)
        .fillna(UNKNOWN_LABEL)
        .astype(str)
        .replace("", UNKNOWN_LABEL)
        .astype("category")
    )
    unique_category_list = category_as_categorical.cat.categories.tolist()
    category_count = len(unique_category_list)
    if category_count <= 0:
        return ["lightgray"] * len(category_as_categorical)

    sampling_points_in_open_interval = [
        (index + 0.5) / category_count for index in range(category_count)
    ]
    base_colorscale = plotly_colors.get_colorscale(colorscale_name)
    sampled_color_palette = plotly_colors.sample_colorscale(
        base_colorscale,
        sampling_points_in_open_interval,
    )
    category_to_color_map = dict(zip(unique_category_list, sampled_color_palette))
    return category_as_categorical.map(category_to_color_map).tolist()


def build_pca_plot_title(
    *,
    selected_pca_component_count: int,
    selected_cluster_count_k: int,
    rendered_plot_dimension_count: int,
    variance_explained_rendered_fraction: float,
    variance_explained_total_fraction: float,
    final_sampled_silhouette_value: float,
    selection_reason: str,
) -> str:
    return (
        f"PCA {rendered_plot_dimension_count}D - PCA(m={selected_pca_component_count}) "
        f"+ KMeans(k={selected_cluster_count_k}) "
        f"- Var(PC1-{rendered_plot_dimension_count})="
        f"{variance_explained_rendered_fraction * 100:.2f}% "
        f"- Var(total m)={variance_explained_total_fraction * 100:.2f}% "
        f"- Silhouette(sampled)={final_sampled_silhouette_value:.3f} "
        f"- Selection={selection_reason}"
    )


def build_pca_plot_figure(
    *,
    plot_dataframe: pd.DataFrame,
    figure_title: str,
    scheme_name_to_column_name: Mapping[str, str] = DEFAULT_SCHEME_NAME_TO_COLUMN_NAME,
    colorscale_name: str = DEFAULT_COLOR_SCALE_NAME,
) -> go.Figure:
    required_columns = (
        "plot_pc1",
        "plot_pc2",
        "cluster_label",
        "description_label",
        "organism_label",
    )
    missing_required_columns = [
        column_name for column_name in required_columns if column_name not in plot_dataframe.columns
    ]
    if missing_required_columns:
        raise RuntimeError(
            "Plot dataframe is missing required columns for rendering: "
            f"{missing_required_columns}."
        )

    plot_dimension_count = _resolve_plot_dimension_count(plot_dataframe)
    if plot_dimension_count == 3 and "plot_pc3" not in plot_dataframe.columns:
        raise RuntimeError(
            "Plot dataframe resolved to 3D rendering but is missing plot_pc3."
        )

    plot_traces_list: list[go.Scatter | go.Scatter3d] = []
    scheme_name_list = list(scheme_name_to_column_name.keys())

    for scheme_display_name in scheme_name_list:
        source_column_name = scheme_name_to_column_name[scheme_display_name]
        if source_column_name == "__all_points__":
            marker_color_vector: str | list[str] = "lightgray"
            hover_value_series = pd.Series(
                ["All points"] * len(plot_dataframe),
                index=plot_dataframe.index,
            )
        else:
            source_series = plot_dataframe.get(
                source_column_name,
                pd.Series([UNKNOWN_LABEL] * len(plot_dataframe), index=plot_dataframe.index),
            )
            marker_color_vector = build_distinct_categorical_color_vector(
                category_series=source_series,
                colorscale_name=colorscale_name,
            )
            hover_value_series = (
                pd.Series(source_series, index=plot_dataframe.index)
                .fillna(UNKNOWN_LABEL)
                .astype(str)
                .replace("", UNKNOWN_LABEL)
            )

        trace_kwargs = {
            "x": plot_dataframe["plot_pc1"],
            "y": plot_dataframe["plot_pc2"],
            "mode": "markers",
            "marker": dict(
                size=2.5,
                color=marker_color_vector,
                opacity=0.8,
                line=dict(width=0),
            ),
            "name": scheme_display_name,
            "text": plot_dataframe["description_label"],
            "customdata": np.stack(
                [
                    plot_dataframe["organism_label"].astype(str).to_numpy(),
                    plot_dataframe["cluster_label"].astype(int).to_numpy(),
                    hover_value_series.to_numpy(),
                ],
                axis=-1,
            ),
            "hovertemplate": (
                "<b>%{text}</b><br>"
                "Organism: %{customdata[0]}<br>"
                "Cluster: %{customdata[1]}<br>"
                f"{scheme_display_name}: "
                "%{customdata[2]}<extra></extra>"
            ),
            "visible": (scheme_display_name == "Cluster"),
        }
        if plot_dimension_count == 3:
            scatter_trace = go.Scatter3d(
                z=plot_dataframe["plot_pc3"],
                **trace_kwargs,
            )
        else:
            scatter_trace = go.Scatter(**trace_kwargs)
        plot_traces_list.append(scatter_trace)

    dropdown_button_list: list[dict[str, Any]] = []
    for scheme_index, scheme_display_name in enumerate(scheme_name_list):
        visibility_mask_list = [False] * len(scheme_name_list)
        visibility_mask_list[scheme_index] = True
        dropdown_button_list.append(
            {
                "label": scheme_display_name,
                "method": "update",
                "args": [
                    {"visible": visibility_mask_list},
                    {
                        "title": (
                            f"{figure_title} - Coloring by {scheme_display_name}"
                        )
                    },
                ],
            }
        )

    figure_object = go.Figure(data=plot_traces_list)
    layout_kwargs: dict[str, Any] = {
        "title": figure_title,
        "updatemenus": [
            dict(
                buttons=dropdown_button_list,
                direction="down",
                x=0.02,
                y=1.12,
                showactive=True,
            )
        ],
        "margin": dict(l=0, r=0, b=0, t=50),
    }
    if plot_dimension_count == 3:
        layout_kwargs["scene"] = dict(
            xaxis_title="PC1",
            yaxis_title="PC2",
            zaxis_title="PC3",
        )
    else:
        layout_kwargs["xaxis"] = dict(title="PC1")
        layout_kwargs["yaxis"] = dict(title="PC2")
    figure_object.update_layout(**layout_kwargs)
    return figure_object
