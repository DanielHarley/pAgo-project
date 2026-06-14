from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypeAlias

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    silhouette_score,
)

from src.pago_pipeline.ncbi_snapshot import (
    _replace_latest_directory,
    build_snapshot_directory_name,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

DEFAULT_PC_BIOLOGY_LATEST_DIRECTORY = Path(
    "data/04-analysis/pc_biology_interpretation/latest"
)
DEFAULT_OUTPUT_ROOT = Path("data/04-analysis/cluster_comparison")

INTEGRATED_TABLE_FILE_NAME = "pc_biology_integrated_table.csv"
CLUSTER_ENRICHMENT_FILE_NAME = "cluster_biology_enrichment.csv"
MANIFEST_FILE_NAME = "manifest.json"

GEOMETRIC_METRICS_FILE_NAME = "cluster_geometric_metrics.csv"
ENRICHMENT_COMPARISON_FILE_NAME = "ag_vs_kmeans_cramers_v.csv"
PURITY_COMPARISON_FILE_NAME = "ag_vs_kmeans_purity_entropy.csv"
KMEANS_GA_MATRIX_FILE_NAME = "kmeans_by_ga_cluster_matrix.csv"
KMEANS_SUBDIVISION_FILE_NAME = "kmeans_to_ga_subdivision_summary.csv"
SUMMARY_FILE_NAME = "ag_vs_kmeans_cluster_comparison_summary.md"
PLOT_HTML_FILE_NAME = "ag_vs_kmeans_cluster_comparison.html"

GA_CLUSTER_COLUMN = "ga_cluster_label"
KMEANS_CLUSTER_COLUMN = "kmeans_cluster_label"
DEFAULT_CLUSTER_COLUMNS = (GA_CLUSTER_COLUMN, KMEANS_CLUSTER_COLUMN)
DEFAULT_PC_COLUMNS = tuple(f"pc{component_index}" for component_index in range(1, 11))
DEFAULT_BIOLOGICAL_VARIABLES = (
    "taxonomy__04",
    "taxonomy__03",
    "paper_length_bin",
    "paper_ago_type_family",
    "paper_ago_type_raw",
    "paper_domain_architecture",
    "paper_paz_type",
    "paper_mid_5p_type",
    "paper_mid_5oh_type",
    "paper_has_piwi_catalytic_tetrad",
    "qc__qc_decision",
    "qc__primary_label",
)


@dataclass(frozen=True)
class ClusterComparisonResult:
    snapshot_directory: Path
    latest_directory: Path
    manifest_file_path: Path
    geometric_metrics_file_path: Path
    enrichment_comparison_file_path: Path
    purity_comparison_file_path: Path
    kmeans_ga_matrix_file_path: Path
    kmeans_subdivision_file_path: Path
    summary_file_path: Path
    plot_html_file_path: Path


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
        temporary_file_path = Path(temporary_file.name)
    dataframe.to_csv(temporary_file_path, index=False)
    temporary_file_path.replace(resolved_output_file_path)
    return resolved_output_file_path


def _write_text_atomic(
    *,
    text: str,
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
        suffix=".md",
    ) as temporary_file:
        temporary_file_path = Path(temporary_file.name)
        temporary_file.write(text)
        if not text.endswith("\n"):
            temporary_file.write("\n")
    temporary_file_path.replace(resolved_output_file_path)
    return resolved_output_file_path


def _write_html_atomic(
    *,
    html_text: str,
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
        temporary_file_path = Path(temporary_file.name)
        temporary_file.write(html_text)
    temporary_file_path.replace(resolved_output_file_path)
    return resolved_output_file_path


def _clean_category(value: object, *, missing_label: str = "missing") -> str:
    if pd.isna(value):
        return missing_label
    normalized_value = str(value).strip()
    return normalized_value if normalized_value else missing_label


def _collapse_rare_biological_categories(
    series: pd.Series,
    *,
    minimum_category_count: int,
) -> pd.Series:
    cleaned_series = series.map(_clean_category)
    value_counts = cleaned_series.value_counts(dropna=False)
    rare_values = set(value_counts[value_counts < minimum_category_count].index)
    if not rare_values:
        return cleaned_series
    return cleaned_series.map(
        lambda value: "other_rare" if value in rare_values else value
    )


def _entropy_bits_from_counts(counts: np.ndarray) -> float:
    total = float(counts.sum())
    if total <= 0:
        return np.nan
    probabilities = counts.astype(float) / total
    probabilities = probabilities[probabilities > 0]
    return float(-(probabilities * np.log2(probabilities)).sum())


def _safe_metric_value(function, coordinates: np.ndarray, labels: np.ndarray) -> float:
    try:
        unique_label_count = len(np.unique(labels))
        if unique_label_count < 2 or unique_label_count >= len(labels):
            return np.nan
        return float(function(coordinates, labels))
    except Exception:
        return np.nan


def compute_geometric_metrics(
    *,
    integrated_dataframe: pd.DataFrame,
    cluster_columns: tuple[str, ...] = DEFAULT_CLUSTER_COLUMNS,
    pc_columns: tuple[str, ...] = DEFAULT_PC_COLUMNS,
) -> pd.DataFrame:
    available_pc_columns = [
        column_name
        for column_name in pc_columns
        if column_name in integrated_dataframe.columns
    ]
    if len(available_pc_columns) < 2:
        raise ValueError("At least two PC columns are required for geometric metrics.")

    coordinates = integrated_dataframe[available_pc_columns].to_numpy(dtype=float)
    rows: list[dict[str, Any]] = []
    for cluster_column in cluster_columns:
        if cluster_column not in integrated_dataframe.columns:
            continue
        labels = integrated_dataframe[cluster_column].map(_clean_category).to_numpy()
        cluster_sizes = pd.Series(labels).value_counts().sort_index()
        rows.append(
            {
                "cluster_column": cluster_column,
                "sample_count": int(len(labels)),
                "cluster_count": int(cluster_sizes.size),
                "minimum_cluster_size": int(cluster_sizes.min()),
                "maximum_cluster_size": int(cluster_sizes.max()),
                "largest_cluster_fraction": float(cluster_sizes.max() / len(labels)),
                "silhouette": _safe_metric_value(
                    silhouette_score,
                    coordinates,
                    labels,
                ),
                "davies_bouldin": _safe_metric_value(
                    davies_bouldin_score,
                    coordinates,
                    labels,
                ),
                "calinski_harabasz": _safe_metric_value(
                    calinski_harabasz_score,
                    coordinates,
                    labels,
                ),
            }
        )

    if all(
        column_name in integrated_dataframe.columns
        for column_name in (GA_CLUSTER_COLUMN, KMEANS_CLUSTER_COLUMN)
    ):
        ga_labels = integrated_dataframe[GA_CLUSTER_COLUMN].map(_clean_category)
        kmeans_labels = integrated_dataframe[KMEANS_CLUSTER_COLUMN].map(_clean_category)
        rows.append(
            {
                "cluster_column": "ag_vs_kmeans",
                "sample_count": int(len(integrated_dataframe)),
                "cluster_count": np.nan,
                "minimum_cluster_size": np.nan,
                "maximum_cluster_size": np.nan,
                "largest_cluster_fraction": np.nan,
                "silhouette": np.nan,
                "davies_bouldin": np.nan,
                "calinski_harabasz": np.nan,
                "adjusted_rand_index": float(
                    adjusted_rand_score(kmeans_labels, ga_labels)
                ),
                "normalized_mutual_information": float(
                    normalized_mutual_info_score(kmeans_labels, ga_labels)
                ),
            }
        )

    return pd.DataFrame(rows)


def build_enrichment_comparison(
    *,
    cluster_enrichment_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    required_columns = {
        "cluster_column",
        "biological_variable",
        "effect_size",
        "q_value_bh",
    }
    missing_columns = required_columns - set(cluster_enrichment_dataframe.columns)
    if missing_columns:
        raise ValueError(f"Missing enrichment columns: {sorted(missing_columns)}")

    subset = cluster_enrichment_dataframe[
        cluster_enrichment_dataframe["cluster_column"].isin(DEFAULT_CLUSTER_COLUMNS)
    ].copy()
    effect_pivot = subset.pivot_table(
        index="biological_variable",
        columns="cluster_column",
        values="effect_size",
        aggfunc="first",
    )
    q_pivot = subset.pivot_table(
        index="biological_variable",
        columns="cluster_column",
        values="q_value_bh",
        aggfunc="first",
    )

    comparison = pd.DataFrame(index=effect_pivot.index)
    comparison["cramers_v_ag"] = effect_pivot.get(GA_CLUSTER_COLUMN)
    comparison["cramers_v_kmeans"] = effect_pivot.get(KMEANS_CLUSTER_COLUMN)
    comparison["q_value_bh_ag"] = q_pivot.get(GA_CLUSTER_COLUMN)
    comparison["q_value_bh_kmeans"] = q_pivot.get(KMEANS_CLUSTER_COLUMN)
    comparison["delta_cramers_v_ag_minus_kmeans"] = (
        comparison["cramers_v_ag"] - comparison["cramers_v_kmeans"]
    )
    comparison["better_by_cramers_v"] = np.select(
        [
            comparison["delta_cramers_v_ag_minus_kmeans"] > 0,
            comparison["delta_cramers_v_ag_minus_kmeans"] < 0,
        ],
        [
            "AG",
            "KMeans",
        ],
        default="tie",
    )
    return (
        comparison.reset_index()
        .sort_values(
            ["delta_cramers_v_ag_minus_kmeans", "cramers_v_ag"],
            ascending=[False, False],
        )
        .reset_index(drop=True)
    )


def _cluster_variable_metrics(
    *,
    dataframe: pd.DataFrame,
    cluster_column: str,
    biological_variable: str,
    minimum_biological_category_count: int,
) -> dict[str, Any]:
    working = dataframe[[cluster_column, biological_variable]].dropna().copy()
    working[cluster_column] = working[cluster_column].map(_clean_category)
    working[biological_variable] = _collapse_rare_biological_categories(
        working[biological_variable],
        minimum_category_count=minimum_biological_category_count,
    )
    contingency = pd.crosstab(working[cluster_column], working[biological_variable])
    counts = contingency.to_numpy(dtype=float)
    sample_count = int(counts.sum())
    if sample_count == 0 or contingency.shape[0] == 0 or contingency.shape[1] == 0:
        return {
            "cluster_column": cluster_column,
            "biological_variable": biological_variable,
            "sample_count": sample_count,
        }

    cluster_totals = counts.sum(axis=1)
    majority_counts = counts.max(axis=1)
    cluster_purity = majority_counts / cluster_totals
    cluster_entropy = np.array(
        [_entropy_bits_from_counts(cluster_counts) for cluster_counts in counts],
        dtype=float,
    )
    maximum_entropy = math.log2(contingency.shape[1]) if contingency.shape[1] > 1 else 0
    weighted_entropy = float(np.sum(cluster_entropy * cluster_totals) / sample_count)
    normalized_entropy = (
        float(weighted_entropy / maximum_entropy)
        if maximum_entropy > 0
        else 0.0
    )

    return {
        "cluster_column": cluster_column,
        "biological_variable": biological_variable,
        "sample_count": sample_count,
        "cluster_count": int(contingency.shape[0]),
        "biological_category_count": int(contingency.shape[1]),
        "weighted_majority_purity": float(majority_counts.sum() / sample_count),
        "macro_majority_purity": float(np.mean(cluster_purity)),
        "cluster_size_weighted_entropy_bits": weighted_entropy,
        "cluster_size_weighted_normalized_entropy": normalized_entropy,
        "macro_entropy_bits": float(np.mean(cluster_entropy)),
    }


def compute_purity_entropy_metrics(
    *,
    integrated_dataframe: pd.DataFrame,
    biological_variables: tuple[str, ...] = DEFAULT_BIOLOGICAL_VARIABLES,
    cluster_columns: tuple[str, ...] = DEFAULT_CLUSTER_COLUMNS,
    minimum_biological_category_count: int = 5,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cluster_column in cluster_columns:
        if cluster_column not in integrated_dataframe.columns:
            continue
        for biological_variable in biological_variables:
            if biological_variable not in integrated_dataframe.columns:
                continue
            rows.append(
                _cluster_variable_metrics(
                    dataframe=integrated_dataframe,
                    cluster_column=cluster_column,
                    biological_variable=biological_variable,
                    minimum_biological_category_count=minimum_biological_category_count,
                )
            )
    metrics = pd.DataFrame(rows)
    if metrics.empty:
        return metrics

    index_columns = ["biological_variable"]
    value_columns = [
        "weighted_majority_purity",
        "macro_majority_purity",
        "cluster_size_weighted_entropy_bits",
        "cluster_size_weighted_normalized_entropy",
    ]
    comparison = metrics.pivot_table(
        index=index_columns,
        columns="cluster_column",
        values=value_columns,
        aggfunc="first",
    )
    flattened_columns = []
    for metric_name, cluster_column in comparison.columns:
        method_name = "ag" if cluster_column == GA_CLUSTER_COLUMN else "kmeans"
        flattened_columns.append(f"{metric_name}_{method_name}")
    comparison.columns = flattened_columns
    comparison = comparison.reset_index()

    comparison["delta_weighted_majority_purity_ag_minus_kmeans"] = (
        comparison.get("weighted_majority_purity_ag")
        - comparison.get("weighted_majority_purity_kmeans")
    )
    comparison["delta_normalized_entropy_kmeans_minus_ag"] = (
        comparison.get("cluster_size_weighted_normalized_entropy_kmeans")
        - comparison.get("cluster_size_weighted_normalized_entropy_ag")
    )
    comparison["better_by_weighted_purity"] = np.select(
        [
            comparison["delta_weighted_majority_purity_ag_minus_kmeans"] > 0,
            comparison["delta_weighted_majority_purity_ag_minus_kmeans"] < 0,
        ],
        ["AG", "KMeans"],
        default="tie",
    )
    comparison["better_by_weighted_entropy"] = np.select(
        [
            comparison["delta_normalized_entropy_kmeans_minus_ag"] > 0,
            comparison["delta_normalized_entropy_kmeans_minus_ag"] < 0,
        ],
        ["AG", "KMeans"],
        default="tie",
    )
    return comparison.sort_values(
        [
            "delta_weighted_majority_purity_ag_minus_kmeans",
            "delta_normalized_entropy_kmeans_minus_ag",
        ],
        ascending=[False, False],
    ).reset_index(drop=True)


def build_kmeans_ga_matrix(
    *,
    integrated_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    matrix = pd.crosstab(
        integrated_dataframe[KMEANS_CLUSTER_COLUMN].map(_clean_category),
        integrated_dataframe[GA_CLUSTER_COLUMN].map(_clean_category),
    )
    matrix.index.name = KMEANS_CLUSTER_COLUMN
    matrix.columns.name = GA_CLUSTER_COLUMN
    return matrix.reset_index()


def build_kmeans_subdivision_summary(
    *,
    integrated_dataframe: pd.DataFrame,
    minimum_ga_subcluster_size: int = 5,
) -> pd.DataFrame:
    matrix = pd.crosstab(
        integrated_dataframe[KMEANS_CLUSTER_COLUMN].map(_clean_category),
        integrated_dataframe[GA_CLUSTER_COLUMN].map(_clean_category),
    )
    rows: list[dict[str, Any]] = []
    for kmeans_cluster, counts in matrix.iterrows():
        cluster_total = int(counts.sum())
        nonzero_counts = counts[counts > 0].sort_values(ascending=False)
        dominant_ga_cluster = str(nonzero_counts.index[0])
        dominant_count = int(nonzero_counts.iloc[0])
        entropy_bits = _entropy_bits_from_counts(nonzero_counts.to_numpy(dtype=float))
        maximum_entropy = math.log2(len(nonzero_counts)) if len(nonzero_counts) > 1 else 0
        normalized_entropy = (
            float(entropy_bits / maximum_entropy) if maximum_entropy > 0 else 0.0
        )
        rows.append(
            {
                "kmeans_cluster_label": kmeans_cluster,
                "kmeans_cluster_size": cluster_total,
                "ga_subcluster_count": int(len(nonzero_counts)),
                "ga_subcluster_count_ge_minimum": int(
                    (nonzero_counts >= minimum_ga_subcluster_size).sum()
                ),
                "dominant_ga_cluster_label": dominant_ga_cluster,
                "dominant_ga_cluster_count": dominant_count,
                "dominant_ga_cluster_fraction": float(dominant_count / cluster_total),
                "ga_distribution_entropy_bits": float(entropy_bits),
                "ga_distribution_normalized_entropy": normalized_entropy,
                "top_ga_subclusters": "; ".join(
                    f"{cluster_label}:{int(count)} ({count / cluster_total:.1%})"
                    for cluster_label, count in nonzero_counts.head(5).items()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        [
            "ga_subcluster_count_ge_minimum",
            "ga_distribution_normalized_entropy",
            "kmeans_cluster_size",
        ],
        ascending=[False, False, False],
    )


def build_summary_markdown(
    *,
    geometric_metrics: pd.DataFrame,
    enrichment_comparison: pd.DataFrame,
    purity_comparison: pd.DataFrame,
    subdivision_summary: pd.DataFrame,
) -> str:
    lines: list[str] = [
        "# AG vs KMeans cluster comparison",
        "",
        "This report compares the genetic algorithm (AG) clustering against the KMeans baseline using geometry, biological enrichment, homogeneity metrics, and the KMeans-to-AG subdivision matrix.",
        "",
        "Statistical caution: protein sequences are evolutionarily related, so observations are not fully independent. Significance values should be interpreted cautiously; conclusions prioritize effect sizes, homogeneity metrics, and coherent biological subdivisions.",
        "",
        "## Geometric and agreement metrics",
        "",
    ]

    for _, row in geometric_metrics.iterrows():
        cluster_column = row["cluster_column"]
        if cluster_column == "ag_vs_kmeans":
            lines.append(
                f"- AG vs KMeans: ARI={row.get('adjusted_rand_index', np.nan):.3f}, "
                f"NMI={row.get('normalized_mutual_information', np.nan):.3f}"
            )
            continue
        lines.append(
            f"- {cluster_column}: clusters={int(row['cluster_count'])}, "
            f"silhouette={row['silhouette']:.3f}, "
            f"Davies-Bouldin={row['davies_bouldin']:.3f}, "
            f"Calinski-Harabasz={row['calinski_harabasz']:.1f}, "
            f"largest cluster={row['largest_cluster_fraction']:.1%}"
        )

    lines.extend(["", "## Biological enrichment: Cramer's V", ""])
    for _, row in enrichment_comparison.head(8).iterrows():
        lines.append(
            f"- {row['biological_variable']}: AG={row['cramers_v_ag']:.3f}, "
            f"KMeans={row['cramers_v_kmeans']:.3f}, "
            f"delta={row['delta_cramers_v_ag_minus_kmeans']:+.3f}, "
            f"better={row['better_by_cramers_v']}"
        )

    lines.extend(["", "## Homogeneity: purity and entropy", ""])
    for _, row in purity_comparison.head(8).iterrows():
        lines.append(
            f"- {row['biological_variable']}: weighted purity AG="
            f"{row['weighted_majority_purity_ag']:.3f}, KMeans="
            f"{row['weighted_majority_purity_kmeans']:.3f}, delta="
            f"{row['delta_weighted_majority_purity_ag_minus_kmeans']:+.3f}; "
            f"normalized entropy delta(KMeans-AG)="
            f"{row['delta_normalized_entropy_kmeans_minus_ag']:+.3f}"
        )

    lines.extend(["", "## KMeans clusters subdivided by AG", ""])
    for _, row in subdivision_summary.head(8).iterrows():
        lines.append(
            f"- KMeans {row['kmeans_cluster_label']} (n={int(row['kmeans_cluster_size'])}): "
            f"{int(row['ga_subcluster_count_ge_minimum'])} AG subclusters >= minimum size; "
            f"dominant AG={row['dominant_ga_cluster_label']} "
            f"({row['dominant_ga_cluster_fraction']:.1%}); "
            f"top={row['top_ga_subclusters']}"
        )

    ag_v_wins = int((enrichment_comparison["better_by_cramers_v"] == "AG").sum())
    kmeans_v_wins = int((enrichment_comparison["better_by_cramers_v"] == "KMeans").sum())
    ag_purity_wins = int((purity_comparison["better_by_weighted_purity"] == "AG").sum())
    kmeans_purity_wins = int(
        (purity_comparison["better_by_weighted_purity"] == "KMeans").sum()
    )
    lines.extend(
        [
            "",
            "## Conclusion",
            "",
            f"- AG has higher Cramer's V for {ag_v_wins} variables; KMeans is higher for {kmeans_v_wins}.",
            f"- AG has higher weighted majority purity for {ag_purity_wins} variables; KMeans is higher for {kmeans_purity_wins}.",
            "- The main question for the final report is not whether AG is globally different from KMeans, because ARI is high; it is whether the extra AG clusters split KMeans groups into biologically coherent subgroups.",
        ]
    )
    return "\n".join(lines)


def build_html_report(
    *,
    enrichment_comparison: pd.DataFrame,
    purity_comparison: pd.DataFrame,
    kmeans_ga_matrix: pd.DataFrame,
) -> str:
    enrichment_figure = go.Figure()
    enrichment_figure.add_trace(
        go.Bar(
            x=enrichment_comparison["biological_variable"],
            y=enrichment_comparison["delta_cramers_v_ag_minus_kmeans"],
            marker_color=np.where(
                enrichment_comparison["delta_cramers_v_ag_minus_kmeans"] >= 0,
                "#2ca02c",
                "#d62728",
            ),
            name="AG - KMeans",
        )
    )
    enrichment_figure.update_layout(
        title="Cramer's V difference by biological variable",
        xaxis_title="Biological variable",
        yaxis_title="Delta Cramer's V (AG - KMeans)",
    )

    purity_figure = go.Figure()
    purity_figure.add_trace(
        go.Bar(
            x=purity_comparison["biological_variable"],
            y=purity_comparison["delta_weighted_majority_purity_ag_minus_kmeans"],
            marker_color=np.where(
                purity_comparison[
                    "delta_weighted_majority_purity_ag_minus_kmeans"
                ]
                >= 0,
                "#2ca02c",
                "#d62728",
            ),
            name="AG - KMeans",
        )
    )
    purity_figure.update_layout(
        title="Weighted majority-purity difference by biological variable",
        xaxis_title="Biological variable",
        yaxis_title="Delta weighted purity (AG - KMeans)",
    )

    heatmap_matrix = kmeans_ga_matrix.set_index(KMEANS_CLUSTER_COLUMN)
    heatmap_figure = go.Figure(
        data=go.Heatmap(
            z=heatmap_matrix.to_numpy(dtype=float),
            x=[str(column_name) for column_name in heatmap_matrix.columns],
            y=[str(index_name) for index_name in heatmap_matrix.index],
            colorscale="Blues",
            colorbar={"title": "Count"},
        )
    )
    heatmap_figure.update_layout(
        title="KMeans x AG cluster overlap",
        xaxis_title="AG cluster",
        yaxis_title="KMeans cluster",
    )

    fragments = [
        "<!doctype html><html><head><meta charset='utf-8'><title>AG vs KMeans cluster comparison</title></head><body>",
        "<h1>AG vs KMeans cluster comparison</h1>",
    ]
    for figure_index, figure in enumerate(
        [enrichment_figure, purity_figure, heatmap_figure]
    ):
        fragments.append(
            pio.to_html(
                figure,
                include_plotlyjs="cdn" if figure_index == 0 else False,
                full_html=False,
            )
        )
    fragments.append("</body></html>")
    return "\n".join(fragments)


def _build_manifest(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    pc_biology_latest_directory: Path,
    output_paths: dict[str, Path],
    minimum_biological_category_count: int,
    minimum_ga_subcluster_size: int,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {
        "snapshot_format_version": "1.0",
        "artifact_type": "cluster_comparison_snapshot",
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "pc_biology_latest_directory": str(pc_biology_latest_directory.resolve()),
        "minimum_biological_category_count": int(minimum_biological_category_count),
        "minimum_ga_subcluster_size": int(minimum_ga_subcluster_size),
    }
    for key, output_path in output_paths.items():
        manifest[f"{key}_file_name"] = output_path.name
        manifest[f"{key}_file_sha256"] = sha256_of_file(input_file_path=output_path)
    source_manifest_path = pc_biology_latest_directory / MANIFEST_FILE_NAME
    if source_manifest_path.exists():
        manifest["source_pc_biology_manifest_sha256"] = sha256_of_file(
            input_file_path=source_manifest_path
        )
        manifest["source_pc_biology_manifest"] = read_json_file(
            input_file_path=source_manifest_path
        )
    return manifest


def create_cluster_comparison_snapshot(
    *,
    pc_biology_latest_directory: PathLike = DEFAULT_PC_BIOLOGY_LATEST_DIRECTORY,
    output_root_directory: PathLike = DEFAULT_OUTPUT_ROOT,
    minimum_biological_category_count: int = 5,
    minimum_ga_subcluster_size: int = 5,
    update_latest_directory: bool = True,
    verbose: bool = True,
) -> ClusterComparisonResult:
    resolved_pc_biology_latest_directory = _as_path(pc_biology_latest_directory)
    resolved_output_root_directory = _as_path(output_root_directory)

    integrated_table_path = (
        resolved_pc_biology_latest_directory / INTEGRATED_TABLE_FILE_NAME
    )
    cluster_enrichment_path = (
        resolved_pc_biology_latest_directory / CLUSTER_ENRICHMENT_FILE_NAME
    )

    if verbose:
        print("Loading PC biology interpretation artifacts...")
    integrated_dataframe = pd.read_csv(integrated_table_path, low_memory=False)
    cluster_enrichment_dataframe = pd.read_csv(cluster_enrichment_path)

    if verbose:
        print("Computing AG vs KMeans geometric metrics...")
    geometric_metrics = compute_geometric_metrics(
        integrated_dataframe=integrated_dataframe
    )
    enrichment_comparison = build_enrichment_comparison(
        cluster_enrichment_dataframe=cluster_enrichment_dataframe
    )

    if verbose:
        print("Computing purity, entropy, and KMeans-to-AG subdivisions...")
    purity_comparison = compute_purity_entropy_metrics(
        integrated_dataframe=integrated_dataframe,
        minimum_biological_category_count=minimum_biological_category_count,
    )
    kmeans_ga_matrix = build_kmeans_ga_matrix(
        integrated_dataframe=integrated_dataframe
    )
    subdivision_summary = build_kmeans_subdivision_summary(
        integrated_dataframe=integrated_dataframe,
        minimum_ga_subcluster_size=minimum_ga_subcluster_size,
    )

    if verbose:
        print("Rendering comparison Markdown and HTML report...")
    summary_markdown = build_summary_markdown(
        geometric_metrics=geometric_metrics,
        enrichment_comparison=enrichment_comparison,
        purity_comparison=purity_comparison,
        subdivision_summary=subdivision_summary,
    )
    html_report = build_html_report(
        enrichment_comparison=enrichment_comparison,
        purity_comparison=purity_comparison,
        kmeans_ga_matrix=kmeans_ga_matrix,
    )

    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query="ag_vs_kmeans_cluster_comparison",
    )
    snapshot_directory = (
        resolved_output_root_directory / "snapshots" / snapshot_directory_name
    )
    snapshot_directory.mkdir(parents=True, exist_ok=False)

    geometric_metrics_file_path = snapshot_directory / GEOMETRIC_METRICS_FILE_NAME
    enrichment_comparison_file_path = (
        snapshot_directory / ENRICHMENT_COMPARISON_FILE_NAME
    )
    purity_comparison_file_path = snapshot_directory / PURITY_COMPARISON_FILE_NAME
    kmeans_ga_matrix_file_path = snapshot_directory / KMEANS_GA_MATRIX_FILE_NAME
    kmeans_subdivision_file_path = snapshot_directory / KMEANS_SUBDIVISION_FILE_NAME
    summary_file_path = snapshot_directory / SUMMARY_FILE_NAME
    plot_html_file_path = snapshot_directory / PLOT_HTML_FILE_NAME
    manifest_file_path = snapshot_directory / MANIFEST_FILE_NAME

    _write_dataframe_csv_atomic(
        dataframe=geometric_metrics,
        output_file_path=geometric_metrics_file_path,
    )
    _write_dataframe_csv_atomic(
        dataframe=enrichment_comparison,
        output_file_path=enrichment_comparison_file_path,
    )
    _write_dataframe_csv_atomic(
        dataframe=purity_comparison,
        output_file_path=purity_comparison_file_path,
    )
    _write_dataframe_csv_atomic(
        dataframe=kmeans_ga_matrix,
        output_file_path=kmeans_ga_matrix_file_path,
    )
    _write_dataframe_csv_atomic(
        dataframe=subdivision_summary,
        output_file_path=kmeans_subdivision_file_path,
    )
    _write_text_atomic(text=summary_markdown, output_file_path=summary_file_path)
    _write_html_atomic(html_text=html_report, output_file_path=plot_html_file_path)

    output_paths = {
        "geometric_metrics": geometric_metrics_file_path,
        "enrichment_comparison": enrichment_comparison_file_path,
        "purity_comparison": purity_comparison_file_path,
        "kmeans_ga_matrix": kmeans_ga_matrix_file_path,
        "kmeans_subdivision": kmeans_subdivision_file_path,
        "summary": summary_file_path,
        "plot_html": plot_html_file_path,
    }
    manifest_payload = _build_manifest(
        snapshot_created_at_utc=snapshot_created_at_utc,
        immutable_snapshot_directory_name=snapshot_directory_name,
        immutable_snapshot_relative_path=str(
            snapshot_directory.relative_to(resolved_output_root_directory)
        ),
        pc_biology_latest_directory=resolved_pc_biology_latest_directory,
        output_paths=output_paths,
        minimum_biological_category_count=minimum_biological_category_count,
        minimum_ga_subcluster_size=minimum_ga_subcluster_size,
    )
    write_json_atomic(payload=manifest_payload, output_file_path=manifest_file_path)

    latest_directory = resolved_output_root_directory / "latest"
    if update_latest_directory:
        _replace_latest_directory(
            latest_directory=latest_directory,
            files_to_copy=[
                (geometric_metrics_file_path, GEOMETRIC_METRICS_FILE_NAME),
                (enrichment_comparison_file_path, ENRICHMENT_COMPARISON_FILE_NAME),
                (purity_comparison_file_path, PURITY_COMPARISON_FILE_NAME),
                (kmeans_ga_matrix_file_path, KMEANS_GA_MATRIX_FILE_NAME),
                (kmeans_subdivision_file_path, KMEANS_SUBDIVISION_FILE_NAME),
                (summary_file_path, SUMMARY_FILE_NAME),
                (plot_html_file_path, PLOT_HTML_FILE_NAME),
                (manifest_file_path, MANIFEST_FILE_NAME),
            ],
        )

    return ClusterComparisonResult(
        snapshot_directory=snapshot_directory,
        latest_directory=latest_directory,
        manifest_file_path=manifest_file_path,
        geometric_metrics_file_path=geometric_metrics_file_path,
        enrichment_comparison_file_path=enrichment_comparison_file_path,
        purity_comparison_file_path=purity_comparison_file_path,
        kmeans_ga_matrix_file_path=kmeans_ga_matrix_file_path,
        kmeans_subdivision_file_path=kmeans_subdivision_file_path,
        summary_file_path=summary_file_path,
        plot_html_file_path=plot_html_file_path,
    )
