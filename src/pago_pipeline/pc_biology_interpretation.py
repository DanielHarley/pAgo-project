from __future__ import annotations

import math
import re
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypeAlias

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from scipy import stats
from sklearn.linear_model import LinearRegression

from src.pago_pipeline.ncbi_snapshot import (
    _replace_latest_directory,
    build_snapshot_directory_name,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

DEFAULT_EXCEL_SHEET_NAME = "novel_pago_fixed_july2018"
DEFAULT_OUTPUT_ROOT = Path("data/04-analysis/pc_biology_interpretation")
DEFAULT_PERMUTATION_COUNT = 999
DEFAULT_RANDOM_STATE = 42

INTEGRATED_TABLE_FILE_NAME = "pc_biology_integrated_table.csv"
ASSOCIATIONS_FILE_NAME = "pc_variable_associations.csv"
CLUSTER_ENRICHMENT_FILE_NAME = "cluster_biology_enrichment.csv"
SUMMARY_FILE_NAME = "pc_axis_interpretation_summary.md"
PLOT_HTML_FILE_NAME = "pc_biology_3d.html"
MANIFEST_FILE_NAME = "manifest.json"

DOMAIN_INTERVAL_PATTERN = re.compile(r"(?P<start>\d+)\s*[-:]\s*(?P<end>\d+)")


@dataclass(frozen=True)
class PcBiologyInterpretationResult:
    snapshot_directory: Path
    latest_directory: Path
    manifest_file_path: Path
    integrated_table_file_path: Path
    associations_file_path: Path
    cluster_enrichment_file_path: Path
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


def _normalize_accession(value: object) -> str | None:
    if pd.isna(value):
        return None
    normalized_value = str(value).strip()
    return normalized_value or None


def _clean_category(value: object, *, missing_label: str = "missing") -> str:
    if pd.isna(value):
        return missing_label
    normalized_value = str(value).strip()
    return normalized_value if normalized_value else missing_label


def _simplify_ago_type(value: object) -> str:
    normalized_value = _clean_category(value).lower()
    if normalized_value.startswith("longa"):
        return "longA"
    if normalized_value.startswith("longb"):
        return "longB"
    if normalized_value.startswith("short"):
        return "short"
    return "unknown"


def _simplify_putative_motif(value: object) -> str:
    normalized_value = _clean_category(value)
    if normalized_value == "missing":
        return normalized_value
    return normalized_value.replace("-putative", "").replace("_putative", "").strip()


def _parse_domain_interval(value: object) -> tuple[float, float, float]:
    if pd.isna(value):
        return (np.nan, np.nan, np.nan)
    match = DOMAIN_INTERVAL_PATTERN.search(str(value))
    if match is None:
        return (np.nan, np.nan, np.nan)
    start = float(match.group("start"))
    end = float(match.group("end"))
    if end < start:
        start, end = end, start
    return (start, end, end - start + 1.0)


def _add_domain_interval_features(
    *,
    dataframe: pd.DataFrame,
    source_column_name: str,
    output_prefix: str,
    length_column_name: str,
) -> pd.DataFrame:
    interval_dataframe = dataframe[source_column_name].map(_parse_domain_interval).apply(
        pd.Series
    )
    interval_dataframe.columns = [
        f"{output_prefix}_start",
        f"{output_prefix}_end",
        f"{output_prefix}_length",
    ]
    dataframe = pd.concat([dataframe, interval_dataframe], axis=1)
    dataframe[f"has_{output_prefix}_domain_interval"] = dataframe[
        f"{output_prefix}_length"
    ].notna()

    protein_length = pd.to_numeric(dataframe[length_column_name], errors="coerce")
    for coordinate_name in ("start", "end", "length"):
        absolute_column = f"{output_prefix}_{coordinate_name}"
        relative_column = f"{absolute_column}_relative"
        dataframe[relative_column] = dataframe[absolute_column] / protein_length

    return dataframe


def _build_domain_architecture(row: pd.Series) -> str:
    domain_names = []
    if bool(row.get("has_paper_paz_domain_interval")):
        domain_names.append("PAZ")
    if bool(row.get("has_paper_mid_domain_interval")):
        domain_names.append("MID")
    if bool(row.get("has_paper_piwi_domain_interval")):
        domain_names.append("PIWI")
    return "+".join(domain_names) if domain_names else "unknown"


def load_paper_annotation_table(
    *,
    excel_file_path: PathLike,
    sheet_name: str = DEFAULT_EXCEL_SHEET_NAME,
) -> pd.DataFrame:
    paper_dataframe = pd.read_excel(excel_file_path, sheet_name=sheet_name)
    required_columns = {
        "Accession",
        "Description",
        "Species",
        "NCBI taxon_id",
        "length",
        "Ago_type",
        "PAZ",
        "MID",
        "PIWI",
        "PAZ_type",
        "MID_type-5P",
        "MID_type-5OH",
        "PIWI_catalytic tetrad",
    }
    missing_columns = sorted(required_columns.difference(paper_dataframe.columns))
    if missing_columns:
        raise RuntimeError(
            "The paper Excel table is missing required columns: "
            f"{missing_columns}."
        )

    renamed_dataframe = paper_dataframe.rename(
        columns={
            "Accession": "paper_accession",
            "Description": "paper_description",
            "Species": "paper_species",
            "NCBI taxon_id": "paper_taxon_id",
            "length": "paper_length",
            "Ago_type": "paper_ago_type_raw",
            "PAZ": "paper_paz_interval_raw",
            "MID": "paper_mid_interval_raw",
            "PIWI": "paper_piwi_interval_raw",
            "PAZ_type": "paper_paz_type_raw",
            "MID_type-5P": "paper_mid_5p_type_raw",
            "MID_type-5OH": "paper_mid_5oh_type_raw",
            "PIWI_catalytic tetrad": "paper_piwi_catalytic_tetrad_raw",
        }
    ).copy()

    renamed_dataframe["paper_accession"] = renamed_dataframe["paper_accession"].map(
        _normalize_accession
    )
    renamed_dataframe["paper_length"] = pd.to_numeric(
        renamed_dataframe["paper_length"], errors="coerce"
    )
    renamed_dataframe["paper_ago_type_family"] = renamed_dataframe[
        "paper_ago_type_raw"
    ].map(_simplify_ago_type)
    renamed_dataframe["paper_is_truncated"] = renamed_dataframe[
        "paper_ago_type_raw"
    ].map(lambda value: "trun" in _clean_category(value).lower())
    renamed_dataframe["paper_paz_type"] = renamed_dataframe["paper_paz_type_raw"].map(
        _clean_category
    )
    renamed_dataframe["paper_mid_5p_type"] = renamed_dataframe[
        "paper_mid_5p_type_raw"
    ].map(_simplify_putative_motif)
    renamed_dataframe["paper_mid_5oh_type"] = renamed_dataframe[
        "paper_mid_5oh_type_raw"
    ].map(_simplify_putative_motif)
    renamed_dataframe["paper_has_piwi_catalytic_tetrad"] = renamed_dataframe[
        "paper_piwi_catalytic_tetrad_raw"
    ].notna()
    renamed_dataframe["paper_length_bin"] = pd.cut(
        renamed_dataframe["paper_length"],
        bins=[0, 299, 599, 900, 1300, np.inf],
        labels=["lt_300", "300_599", "600_900", "901_1300", "gt_1300"],
        include_lowest=True,
    ).astype("string")

    for source_column_name, output_prefix in (
        ("paper_paz_interval_raw", "paper_paz"),
        ("paper_mid_interval_raw", "paper_mid"),
        ("paper_piwi_interval_raw", "paper_piwi"),
    ):
        renamed_dataframe = _add_domain_interval_features(
            dataframe=renamed_dataframe,
            source_column_name=source_column_name,
            output_prefix=output_prefix,
            length_column_name="paper_length",
        )

    renamed_dataframe["paper_domain_architecture"] = renamed_dataframe.apply(
        _build_domain_architecture,
        axis=1,
    )
    renamed_dataframe["paper_has_all_core_domains"] = (
        renamed_dataframe["has_paper_paz_domain_interval"]
        & renamed_dataframe["has_paper_mid_domain_interval"]
        & renamed_dataframe["has_paper_piwi_domain_interval"]
    )

    if renamed_dataframe["paper_accession"].duplicated().any():
        duplicated_accessions = (
            renamed_dataframe.loc[
                renamed_dataframe["paper_accession"].duplicated(), "paper_accession"
            ]
            .dropna()
            .head(10)
            .tolist()
        )
        raise RuntimeError(
            "The paper Excel accession column must be unique. Example duplicates: "
            f"{duplicated_accessions}."
        )

    return renamed_dataframe


def load_pca_kmeans_table(*, pca_kmeans_latest_directory: PathLike) -> pd.DataFrame:
    latest_directory = _as_path(pca_kmeans_latest_directory)
    cluster_assignments_file_path = latest_directory / "cluster_assignments.csv"
    pca_coordinates_file_path = latest_directory / "pca_coordinates_10D.npy"
    if not cluster_assignments_file_path.exists():
        raise FileNotFoundError(cluster_assignments_file_path)
    if not pca_coordinates_file_path.exists():
        raise FileNotFoundError(pca_coordinates_file_path)

    pca_dataframe = pd.read_csv(cluster_assignments_file_path, low_memory=False)
    pca_coordinates = np.load(pca_coordinates_file_path)
    if pca_coordinates.shape[0] != len(pca_dataframe):
        raise RuntimeError(
            "PCA coordinates row count does not match cluster assignments row count: "
            f"{pca_coordinates.shape[0]} != {len(pca_dataframe)}."
        )
    if pca_coordinates.shape[1] < 10:
        raise RuntimeError(
            "Expected at least 10 PCA components in pca_coordinates_10D.npy."
        )

    for component_index in range(10):
        pca_dataframe[f"pc{component_index + 1}"] = pca_coordinates[:, component_index]

    if "cluster_label" in pca_dataframe.columns:
        pca_dataframe["kmeans_cluster_label"] = pca_dataframe["cluster_label"]

    pca_dataframe["gbseq__accession_version"] = pca_dataframe[
        "gbseq__accession_version"
    ].map(_normalize_accession)

    return pca_dataframe


def _load_ga_cluster_labels(
    *,
    ga_latest_directory: PathLike | None,
) -> pd.DataFrame | None:
    if ga_latest_directory is None:
        return None
    assignments_file_path = _as_path(ga_latest_directory) / "ga_cluster_assignments.csv"
    if not assignments_file_path.exists():
        return None

    ga_dataframe = pd.read_csv(assignments_file_path, low_memory=False)
    required_columns = {"gbseq__accession_version", "ga_cluster_label"}
    if not required_columns.issubset(ga_dataframe.columns):
        return None

    selected_columns = [
        column_name
        for column_name in (
            "gbseq__accession_version",
            "ga_cluster_label",
            "kmeans_cluster_label",
        )
        if column_name in ga_dataframe.columns
    ]
    ga_dataframe = ga_dataframe[selected_columns].copy()
    ga_dataframe["gbseq__accession_version"] = ga_dataframe[
        "gbseq__accession_version"
    ].map(_normalize_accession)
    return ga_dataframe.drop_duplicates(subset=["gbseq__accession_version"])


def _load_qc_annotations(*, qc_latest_directory: PathLike | None) -> pd.DataFrame | None:
    if qc_latest_directory is None:
        return None
    labelled_records_file_path = _as_path(qc_latest_directory) / "labelled_records.csv"
    if not labelled_records_file_path.exists():
        return None

    qc_dataframe = pd.read_csv(labelled_records_file_path, low_memory=False)
    selected_columns = [
        column_name
        for column_name in (
            "protein_uid",
            "gbseq__accession_version",
            "taxonomy__03",
            "taxonomy__04",
            "length_bin",
            "has_paz_region",
            "has_mid_region",
            "has_piwi_region",
            "has_active_site_annotation",
            "has_cdd_region",
            "primary_label",
            "qc_decision",
            "confidence_score",
        )
        if column_name in qc_dataframe.columns
    ]
    qc_dataframe = qc_dataframe[selected_columns].copy()
    qc_dataframe["gbseq__accession_version"] = qc_dataframe[
        "gbseq__accession_version"
    ].map(_normalize_accession)
    qc_dataframe = qc_dataframe.add_prefix("qc__")
    return qc_dataframe.drop_duplicates(subset=["qc__gbseq__accession_version"])


def build_integrated_pc_biology_table(
    *,
    excel_file_path: PathLike,
    pca_kmeans_latest_directory: PathLike,
    ga_latest_directory: PathLike | None = None,
    qc_latest_directory: PathLike | None = None,
    excel_sheet_name: str = DEFAULT_EXCEL_SHEET_NAME,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    paper_dataframe = load_paper_annotation_table(
        excel_file_path=excel_file_path,
        sheet_name=excel_sheet_name,
    )
    pca_dataframe = load_pca_kmeans_table(
        pca_kmeans_latest_directory=pca_kmeans_latest_directory
    )

    integrated_dataframe = pca_dataframe.merge(
        paper_dataframe,
        how="left",
        left_on="gbseq__accession_version",
        right_on="paper_accession",
        validate="one_to_one",
    )

    ga_dataframe = _load_ga_cluster_labels(ga_latest_directory=ga_latest_directory)
    if ga_dataframe is not None:
        integrated_dataframe = integrated_dataframe.drop(
            columns=[
                column_name
                for column_name in ("ga_cluster_label",)
                if column_name in integrated_dataframe.columns
            ],
            errors="ignore",
        ).merge(
            ga_dataframe,
            how="left",
            on="gbseq__accession_version",
            validate="one_to_one",
            suffixes=("", "_ga"),
        )
        if "kmeans_cluster_label_ga" in integrated_dataframe.columns:
            integrated_dataframe["kmeans_cluster_label"] = integrated_dataframe[
                "kmeans_cluster_label_ga"
            ].combine_first(integrated_dataframe.get("kmeans_cluster_label"))
            integrated_dataframe = integrated_dataframe.drop(
                columns=["kmeans_cluster_label_ga"]
            )

    qc_dataframe = _load_qc_annotations(qc_latest_directory=qc_latest_directory)
    if qc_dataframe is not None:
        integrated_dataframe = integrated_dataframe.merge(
            qc_dataframe,
            how="left",
            left_on="gbseq__accession_version",
            right_on="qc__gbseq__accession_version",
            validate="one_to_one",
        )

    paper_match_count = int(integrated_dataframe["paper_accession"].notna().sum())
    if paper_match_count != len(integrated_dataframe):
        missing_examples = (
            integrated_dataframe.loc[
                integrated_dataframe["paper_accession"].isna(),
                "gbseq__accession_version",
            ]
            .dropna()
            .head(10)
            .tolist()
        )
        raise RuntimeError(
            "PCA/KMeans records do not fully match the paper Excel accessions. "
            f"Matched {paper_match_count}/{len(integrated_dataframe)}. "
            f"Examples missing from paper table: {missing_examples}."
        )

    integrated_dataframe["paper_length_minus_ncbi_sequence_length"] = (
        pd.to_numeric(integrated_dataframe["paper_length"], errors="coerce")
        - pd.to_numeric(integrated_dataframe["sequence_length"], errors="coerce")
    )

    integration_metadata = {
        "pca_record_count": int(len(pca_dataframe)),
        "paper_record_count": int(len(paper_dataframe)),
        "paper_match_count": paper_match_count,
        "ga_labels_loaded": ga_dataframe is not None,
        "qc_annotations_loaded": qc_dataframe is not None,
    }

    return integrated_dataframe, integration_metadata


def _benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    numeric_p_values = pd.to_numeric(p_values, errors="coerce")
    valid_mask = numeric_p_values.notna()
    q_values = pd.Series(np.nan, index=p_values.index, dtype="float64")
    if not valid_mask.any():
        return q_values

    valid_p_values = numeric_p_values.loc[valid_mask].to_numpy(dtype=float)
    order = np.argsort(valid_p_values)
    sorted_p_values = valid_p_values[order]
    rank = np.arange(1, len(sorted_p_values) + 1, dtype=float)
    sorted_q_values = sorted_p_values * len(sorted_p_values) / rank
    sorted_q_values = np.minimum.accumulate(sorted_q_values[::-1])[::-1]
    sorted_q_values = np.clip(sorted_q_values, 0.0, 1.0)

    original_order_q_values = np.empty_like(sorted_q_values)
    original_order_q_values[order] = sorted_q_values
    q_values.loc[valid_mask] = original_order_q_values
    return q_values


def _spearman_permutation_p_value(
    *,
    x_values: np.ndarray,
    y_values: np.ndarray,
    observed_rho: float,
    permutation_count: int,
    rng: np.random.Generator,
) -> float:
    if permutation_count <= 0 or not np.isfinite(observed_rho):
        return np.nan
    x_ranks = stats.rankdata(x_values)
    y_ranks = stats.rankdata(y_values)
    if np.nanstd(x_ranks) == 0 or np.nanstd(y_ranks) == 0:
        return np.nan

    exceedance_count = 0
    observed_abs = abs(observed_rho)
    for _ in range(permutation_count):
        permuted_y_ranks = rng.permutation(y_ranks)
        permuted_rho = float(np.corrcoef(x_ranks, permuted_y_ranks)[0, 1])
        if abs(permuted_rho) >= observed_abs:
            exceedance_count += 1
    return float((exceedance_count + 1) / (permutation_count + 1))


def _collapse_rare_categories(
    series: pd.Series,
    *,
    minimum_category_count: int,
    include_missing: bool = True,
) -> pd.Series:
    if include_missing:
        cleaned_series = series.map(_clean_category)
    else:
        cleaned_series = series.map(
            lambda value: np.nan if pd.isna(value) else _clean_category(value)
        )
    value_counts = cleaned_series.value_counts(dropna=False)
    rare_values = set(value_counts[value_counts < minimum_category_count].index)
    if not rare_values:
        return cleaned_series
    return cleaned_series.map(
        lambda value: "other_rare" if value in rare_values else value
    )


def _compute_continuous_associations(
    *,
    dataframe: pd.DataFrame,
    pc_columns: list[str],
    continuous_columns: list[str],
    permutation_count: int,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pc_column in pc_columns:
        pc_values = pd.to_numeric(dataframe[pc_column], errors="coerce")
        for variable_name in continuous_columns:
            variable_values = pd.to_numeric(dataframe[variable_name], errors="coerce")
            valid_mask = pc_values.notna() & variable_values.notna()
            sample_count = int(valid_mask.sum())
            if sample_count < 4:
                continue
            x_values = variable_values.loc[valid_mask].to_numpy(dtype=float)
            y_values = pc_values.loc[valid_mask].to_numpy(dtype=float)
            if np.nanstd(x_values) == 0 or np.nanstd(y_values) == 0:
                continue

            statistic = stats.spearmanr(x_values, y_values)
            rho = float(statistic.statistic)
            p_value = _spearman_permutation_p_value(
                x_values=x_values,
                y_values=y_values,
                observed_rho=rho,
                permutation_count=permutation_count,
                rng=rng,
            )
            rows.append(
                {
                    "analysis_scope": "raw_pc",
                    "pc_name": pc_column,
                    "variable_name": variable_name,
                    "variable_kind": "continuous",
                    "test_name": "spearman_permutation",
                    "sample_count": sample_count,
                    "group_count": np.nan,
                    "statistic_name": "spearman_rho",
                    "statistic": rho,
                    "effect_size_name": "abs_spearman_rho",
                    "effect_size": abs(rho),
                    "p_value": p_value,
                    "direction": "positive" if rho >= 0 else "negative",
                    "top_group_high": np.nan,
                    "top_group_low": np.nan,
                }
            )
    return rows


def _categorical_test_rows_for_columns(
    *,
    dataframe: pd.DataFrame,
    pc_columns: list[str],
    categorical_columns: list[str],
    analysis_scope: str,
    minimum_category_count: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pc_column in pc_columns:
        pc_values = pd.to_numeric(dataframe[pc_column], errors="coerce")
        for variable_name in categorical_columns:
            category_values = _collapse_rare_categories(
                dataframe[variable_name],
                minimum_category_count=minimum_category_count,
            )
            valid_mask = pc_values.notna() & category_values.notna()
            valid_pc_values = pc_values.loc[valid_mask]
            valid_category_values = category_values.loc[valid_mask]
            grouped_values = [
                valid_pc_values.loc[valid_category_values == category].to_numpy(
                    dtype=float
                )
                for category in sorted(valid_category_values.unique())
            ]
            grouped_values = [values for values in grouped_values if len(values) > 1]
            group_count = len(grouped_values)
            sample_count = int(sum(len(values) for values in grouped_values))
            if group_count < 2 or sample_count <= group_count:
                continue

            statistic = stats.kruskal(*grouped_values)
            h_value = float(statistic.statistic)
            epsilon_squared = (h_value - group_count + 1.0) / (
                sample_count - group_count
            )
            epsilon_squared = float(np.clip(epsilon_squared, 0.0, 1.0))
            group_medians = (
                pd.DataFrame(
                    {
                        "pc_value": valid_pc_values,
                        "category": valid_category_values,
                    }
                )
                .groupby("category", dropna=False)["pc_value"]
                .median()
                .sort_values()
            )
            rows.append(
                {
                    "analysis_scope": analysis_scope,
                    "pc_name": pc_column.replace("_length_residual", ""),
                    "variable_name": variable_name,
                    "variable_kind": "categorical",
                    "test_name": "kruskal_wallis",
                    "sample_count": sample_count,
                    "group_count": group_count,
                    "statistic_name": "kruskal_h",
                    "statistic": h_value,
                    "effect_size_name": "epsilon_squared",
                    "effect_size": epsilon_squared,
                    "p_value": float(statistic.pvalue),
                    "direction": "category_shift",
                    "top_group_high": str(group_medians.index[-1]),
                    "top_group_low": str(group_medians.index[0]),
                }
            )
    return rows


def _add_length_residual_pc_columns(
    *,
    dataframe: pd.DataFrame,
    pc_columns: list[str],
    length_column_name: str,
) -> tuple[pd.DataFrame, list[str]]:
    output_dataframe = dataframe.copy()
    length_values = pd.to_numeric(
        output_dataframe[length_column_name], errors="coerce"
    ).to_numpy(dtype=float)
    valid_length_mask = np.isfinite(length_values)
    if int(valid_length_mask.sum()) < 3:
        return output_dataframe, []

    residual_columns: list[str] = []
    x_values = length_values[valid_length_mask].reshape(-1, 1)
    for pc_column in pc_columns:
        pc_values = pd.to_numeric(output_dataframe[pc_column], errors="coerce").to_numpy(
            dtype=float
        )
        valid_mask = valid_length_mask & np.isfinite(pc_values)
        if int(valid_mask.sum()) < 3:
            continue
        model = LinearRegression()
        model.fit(length_values[valid_mask].reshape(-1, 1), pc_values[valid_mask])
        predictions = np.full(len(output_dataframe), np.nan, dtype=float)
        predictions[valid_length_mask] = model.predict(x_values)
        residual_column = f"{pc_column}_length_residual"
        output_dataframe[residual_column] = pc_values - predictions
        residual_columns.append(residual_column)

    return output_dataframe, residual_columns


def compute_pc_biology_associations(
    *,
    integrated_dataframe: pd.DataFrame,
    permutation_count: int = DEFAULT_PERMUTATION_COUNT,
    random_state: int = DEFAULT_RANDOM_STATE,
    minimum_category_count: int = 5,
) -> pd.DataFrame:
    pc_columns = [f"pc{component_index}" for component_index in range(1, 11)]
    missing_pc_columns = [
        column_name
        for column_name in pc_columns
        if column_name not in integrated_dataframe.columns
    ]
    if missing_pc_columns:
        raise RuntimeError(f"Missing PCA columns: {missing_pc_columns}.")

    continuous_columns = [
        column_name
        for column_name in (
            "paper_length",
            "sequence_length",
            "paper_paz_length",
            "paper_mid_length",
            "paper_piwi_length",
            "paper_paz_start_relative",
            "paper_mid_start_relative",
            "paper_piwi_start_relative",
            "paper_paz_length_relative",
            "paper_mid_length_relative",
            "paper_piwi_length_relative",
            "qc__confidence_score",
        )
        if column_name in integrated_dataframe.columns
    ]
    categorical_columns = [
        column_name
        for column_name in (
            "paper_ago_type_raw",
            "paper_ago_type_family",
            "paper_is_truncated",
            "paper_domain_architecture",
            "paper_paz_type",
            "paper_mid_5p_type",
            "paper_mid_5oh_type",
            "paper_has_piwi_catalytic_tetrad",
            "paper_length_bin",
            "taxonomy__03",
            "taxonomy__04",
            "qc__length_bin",
            "qc__has_paz_region",
            "qc__has_mid_region",
            "qc__has_piwi_region",
            "qc__has_active_site_annotation",
            "qc__primary_label",
            "qc__qc_decision",
        )
        if column_name in integrated_dataframe.columns
    ]

    rng = np.random.default_rng(random_state)
    rows = _compute_continuous_associations(
        dataframe=integrated_dataframe,
        pc_columns=pc_columns,
        continuous_columns=continuous_columns,
        permutation_count=permutation_count,
        rng=rng,
    )
    rows.extend(
        _categorical_test_rows_for_columns(
            dataframe=integrated_dataframe,
            pc_columns=pc_columns,
            categorical_columns=categorical_columns,
            analysis_scope="raw_pc",
            minimum_category_count=minimum_category_count,
        )
    )

    residual_dataframe, residual_pc_columns = _add_length_residual_pc_columns(
        dataframe=integrated_dataframe,
        pc_columns=pc_columns,
        length_column_name="paper_length",
    )
    if residual_pc_columns:
        rows.extend(
            _categorical_test_rows_for_columns(
                dataframe=residual_dataframe,
                pc_columns=residual_pc_columns,
                categorical_columns=categorical_columns,
                analysis_scope="length_residual_pc",
                minimum_category_count=minimum_category_count,
            )
        )

    associations_dataframe = pd.DataFrame(rows)
    if associations_dataframe.empty:
        return associations_dataframe

    associations_dataframe["q_value_bh"] = _benjamini_hochberg(
        associations_dataframe["p_value"]
    )
    associations_dataframe = associations_dataframe.sort_values(
        ["analysis_scope", "pc_name", "q_value_bh", "effect_size"],
        ascending=[True, True, True, False],
    ).reset_index(drop=True)
    return associations_dataframe


def compute_cluster_biology_enrichment(
    *,
    integrated_dataframe: pd.DataFrame,
    minimum_category_count: int = 5,
) -> pd.DataFrame:
    cluster_columns = [
        column_name
        for column_name in ("ga_cluster_label", "kmeans_cluster_label")
        if column_name in integrated_dataframe.columns
    ]
    biological_columns = [
        column_name
        for column_name in (
            "paper_ago_type_raw",
            "paper_ago_type_family",
            "paper_domain_architecture",
            "paper_paz_type",
            "paper_mid_5p_type",
            "paper_mid_5oh_type",
            "paper_has_piwi_catalytic_tetrad",
            "paper_length_bin",
            "taxonomy__03",
            "taxonomy__04",
            "qc__primary_label",
            "qc__qc_decision",
        )
        if column_name in integrated_dataframe.columns
    ]

    rows: list[dict[str, Any]] = []
    for cluster_column in cluster_columns:
        cluster_values = _collapse_rare_categories(
            integrated_dataframe[cluster_column],
            minimum_category_count=minimum_category_count,
        )
        for biological_column in biological_columns:
            biological_values = _collapse_rare_categories(
                integrated_dataframe[biological_column],
                minimum_category_count=minimum_category_count,
            )
            contingency_table = pd.crosstab(cluster_values, biological_values)
            if contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
                continue
            chi2_value, p_value, dof, _ = stats.chi2_contingency(contingency_table)
            sample_count = int(contingency_table.to_numpy().sum())
            denominator = sample_count * min(
                contingency_table.shape[0] - 1,
                contingency_table.shape[1] - 1,
            )
            cramer_v = math.sqrt(chi2_value / denominator) if denominator > 0 else np.nan
            rows.append(
                {
                    "cluster_column": cluster_column,
                    "biological_variable": biological_column,
                    "test_name": "chi_square_independence",
                    "sample_count": sample_count,
                    "cluster_group_count": int(contingency_table.shape[0]),
                    "biological_group_count": int(contingency_table.shape[1]),
                    "statistic_name": "chi2",
                    "statistic": float(chi2_value),
                    "degrees_of_freedom": int(dof),
                    "effect_size_name": "cramers_v",
                    "effect_size": float(cramer_v),
                    "p_value": float(p_value),
                }
            )

    enrichment_dataframe = pd.DataFrame(rows)
    if enrichment_dataframe.empty:
        return enrichment_dataframe
    enrichment_dataframe["q_value_bh"] = _benjamini_hochberg(
        enrichment_dataframe["p_value"]
    )
    return enrichment_dataframe.sort_values(
        ["cluster_column", "q_value_bh", "effect_size"],
        ascending=[True, True, False],
    ).reset_index(drop=True)


def build_pc_axis_summary_markdown(
    *,
    associations_dataframe: pd.DataFrame,
    cluster_enrichment_dataframe: pd.DataFrame,
    integrated_dataframe: pd.DataFrame,
) -> str:
    lines = [
        "# PC biological interpretation summary",
        "",
        "Interpretation rule: PCs were computed from SWeeP embeddings, not from the "
        "biological annotations below. Therefore, a PC should be described as "
        "associated with a property, not as directly measuring that property.",
        "",
        f"Record count: {len(integrated_dataframe)}",
        "",
    ]

    raw_associations = associations_dataframe[
        associations_dataframe["analysis_scope"] == "raw_pc"
    ].copy()
    residual_associations = associations_dataframe[
        associations_dataframe["analysis_scope"] == "length_residual_pc"
    ].copy()

    for pc_name in [f"pc{component_index}" for component_index in range(1, 11)]:
        lines.append(f"## {pc_name.upper()}")
        pc_rows = raw_associations[raw_associations["pc_name"] == pc_name]
        if pc_rows.empty:
            lines.append("")
            lines.append("No valid association tests were available.")
            lines.append("")
            continue

        continuous_top = pc_rows[pc_rows["variable_kind"] == "continuous"].head(3)
        categorical_top = pc_rows[pc_rows["variable_kind"] == "categorical"].head(3)
        if not continuous_top.empty:
            lines.append("")
            lines.append("Top continuous associations:")
            for _, row in continuous_top.iterrows():
                direction = row["direction"]
                lines.append(
                    "- "
                    f"{row['variable_name']}: rho={row['statistic']:.3f}, "
                    f"q={row['q_value_bh']:.4g}, direction={direction}"
                )
        if not categorical_top.empty:
            lines.append("")
            lines.append("Top categorical associations:")
            for _, row in categorical_top.iterrows():
                lines.append(
                    "- "
                    f"{row['variable_name']}: epsilon_squared="
                    f"{row['effect_size']:.3f}, q={row['q_value_bh']:.4g}, "
                    f"high={row['top_group_high']}, low={row['top_group_low']}"
                )

        pc_residual_rows = residual_associations[
            residual_associations["pc_name"] == pc_name
        ].head(3)
        if not pc_residual_rows.empty:
            lines.append("")
            lines.append("Top categorical associations after length residualization:")
            for _, row in pc_residual_rows.iterrows():
                lines.append(
                    "- "
                    f"{row['variable_name']}: epsilon_squared="
                    f"{row['effect_size']:.3f}, q={row['q_value_bh']:.4g}, "
                    f"high={row['top_group_high']}, low={row['top_group_low']}"
                )
        lines.append("")

    if not cluster_enrichment_dataframe.empty:
        lines.append("## Cluster enrichment")
        for _, row in cluster_enrichment_dataframe.head(10).iterrows():
            lines.append(
                "- "
                f"{row['cluster_column']} vs {row['biological_variable']}: "
                f"Cramer's V={row['effect_size']:.3f}, "
                f"q={row['q_value_bh']:.4g}"
            )
        lines.append("")

    return "\n".join(lines)


def _build_effect_size_heatmap(
    *,
    associations_dataframe: pd.DataFrame,
    variable_kind: str,
    title: str,
) -> go.Figure | None:
    filtered_dataframe = associations_dataframe[
        (associations_dataframe["analysis_scope"] == "raw_pc")
        & (associations_dataframe["variable_kind"] == variable_kind)
    ].copy()
    if filtered_dataframe.empty:
        return None

    pivot_dataframe = filtered_dataframe.pivot_table(
        index="variable_name",
        columns="pc_name",
        values="effect_size",
        aggfunc="max",
    )
    pc_columns = [f"pc{component_index}" for component_index in range(1, 11)]
    pivot_dataframe = pivot_dataframe.reindex(columns=pc_columns)
    figure = go.Figure(
        data=go.Heatmap(
            z=pivot_dataframe.to_numpy(),
            x=[column.upper() for column in pivot_dataframe.columns],
            y=pivot_dataframe.index.tolist(),
            colorscale="Viridis",
            colorbar={"title": "Effect"},
        )
    )
    figure.update_layout(
        title=title,
        margin={"l": 260, "r": 30, "t": 70, "b": 60},
        height=max(450, 28 * len(pivot_dataframe.index) + 180),
    )
    return figure


def build_pc_biology_html_report(
    *,
    integrated_dataframe: pd.DataFrame,
    associations_dataframe: pd.DataFrame,
    cluster_enrichment_dataframe: pd.DataFrame,
) -> str:
    hover_columns = [
        column_name
        for column_name in (
            "gbseq__accession_version",
            "paper_length",
            "paper_ago_type_raw",
            "paper_domain_architecture",
            "paper_paz_type",
            "paper_mid_5p_type",
            "paper_has_piwi_catalytic_tetrad",
            "ga_cluster_label",
            "kmeans_cluster_label",
        )
        if column_name in integrated_dataframe.columns
    ]
    plot_dataframe = integrated_dataframe.copy()
    for column_name in ("ga_cluster_label", "kmeans_cluster_label"):
        if column_name in plot_dataframe.columns:
            plot_dataframe[column_name] = plot_dataframe[column_name].map(_clean_category)

    figure_specs = [
        ("paper_length", "PC1/PC2/PC3 colored by paper protein length"),
        ("paper_ago_type_family", "PC1/PC2/PC3 colored by pAgo family"),
        ("paper_paz_type", "PC1/PC2/PC3 colored by PAZ type"),
        (
            "paper_has_piwi_catalytic_tetrad",
            "PC1/PC2/PC3 colored by PIWI catalytic tetrad annotation",
        ),
        ("paper_mid_5p_type", "PC1/PC2/PC3 colored by MID 5P motif"),
        ("ga_cluster_label", "PC1/PC2/PC3 colored by genetic algorithm cluster"),
        ("kmeans_cluster_label", "PC1/PC2/PC3 colored by KMeans cluster"),
    ]

    figures: list[go.Figure] = []
    for color_column, title in figure_specs:
        if color_column not in plot_dataframe.columns:
            continue
        figure = px.scatter_3d(
            plot_dataframe,
            x="pc1",
            y="pc2",
            z="pc3",
            color=color_column,
            hover_data=hover_columns,
            title=title,
            opacity=0.82,
            height=760,
        )
        figure.update_traces(marker={"size": 4})
        figure.update_layout(margin={"l": 0, "r": 0, "t": 70, "b": 0})
        figures.append(figure)

    for variable_kind, title in (
        ("continuous", "Continuous variable association effect sizes"),
        ("categorical", "Categorical variable association effect sizes"),
    ):
        heatmap_figure = _build_effect_size_heatmap(
            associations_dataframe=associations_dataframe,
            variable_kind=variable_kind,
            title=title,
        )
        if heatmap_figure is not None:
            figures.append(heatmap_figure)

    if not cluster_enrichment_dataframe.empty:
        top_enrichment = cluster_enrichment_dataframe.head(20).copy()
        top_enrichment["comparison"] = (
            top_enrichment["cluster_column"]
            + " vs "
            + top_enrichment["biological_variable"]
        )
        enrichment_figure = px.bar(
            top_enrichment.sort_values("effect_size"),
            x="effect_size",
            y="comparison",
            color="cluster_column",
            orientation="h",
            title="Top cluster/biology enrichment effects",
            height=720,
        )
        enrichment_figure.update_layout(margin={"l": 260, "r": 30, "t": 70, "b": 60})
        figures.append(enrichment_figure)

    html_fragments = []
    for figure_index, figure in enumerate(figures):
        html_fragments.append(
            pio.to_html(
                figure,
                include_plotlyjs="cdn" if figure_index == 0 else False,
                full_html=False,
            )
        )

    return "\n".join(
        [
            "<!doctype html>",
            "<html>",
            "<head>",
            '<meta charset="utf-8">',
            "<title>PC biological interpretation</title>",
            "<style>",
            "body { font-family: Arial, sans-serif; margin: 24px; }",
            "h1, h2 { color: #1f2933; }",
            ".note { max-width: 980px; line-height: 1.45; color: #364152; }",
            ".figure-block { margin-bottom: 42px; }",
            "</style>",
            "</head>",
            "<body>",
            "<h1>PC biological interpretation</h1>",
            '<p class="note">PCs were computed from SWeeP embeddings. The plots '
            "and tests below should be interpreted as biological associations "
            "with PCs, not as direct definitions of what each PC measures.</p>",
            *[
                f'<div class="figure-block">{fragment}</div>'
                for fragment in html_fragments
            ],
            "</body>",
            "</html>",
        ]
    )


def _build_manifest(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    excel_file_path: Path,
    excel_sheet_name: str,
    pca_kmeans_latest_directory: Path,
    ga_latest_directory: Path | None,
    qc_latest_directory: Path | None,
    integrated_table_file_path: Path,
    associations_file_path: Path,
    cluster_enrichment_file_path: Path,
    summary_file_path: Path,
    plot_html_file_path: Path,
    integration_metadata: dict[str, Any],
    permutation_count: int,
    random_state: int,
    minimum_category_count: int,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {
        "snapshot_format_version": "1.0",
        "artifact_type": "pc_biology_interpretation_snapshot",
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "manifest_file_name": MANIFEST_FILE_NAME,
        "integrated_table_file_name": integrated_table_file_path.name,
        "associations_file_name": associations_file_path.name,
        "cluster_enrichment_file_name": cluster_enrichment_file_path.name,
        "summary_file_name": summary_file_path.name,
        "plot_html_file_name": plot_html_file_path.name,
        "integrated_table_file_sha256": sha256_of_file(
            input_file_path=integrated_table_file_path
        ),
        "associations_file_sha256": sha256_of_file(
            input_file_path=associations_file_path
        ),
        "cluster_enrichment_file_sha256": sha256_of_file(
            input_file_path=cluster_enrichment_file_path
        ),
        "summary_file_sha256": sha256_of_file(input_file_path=summary_file_path),
        "plot_html_file_sha256": sha256_of_file(input_file_path=plot_html_file_path),
        "excel_file_path": str(excel_file_path),
        "excel_file_sha256": sha256_of_file(input_file_path=excel_file_path),
        "excel_sheet_name": excel_sheet_name,
        "pca_kmeans_latest_directory": str(pca_kmeans_latest_directory),
        "permutation_count": permutation_count,
        "random_state": random_state,
        "minimum_category_count": minimum_category_count,
        "integration_metadata": integration_metadata,
    }

    source_manifest_specs = {
        "source_pca_kmeans_manifest_sha256": pca_kmeans_latest_directory
        / MANIFEST_FILE_NAME,
        "source_ga_clustering_manifest_sha256": (
            ga_latest_directory / MANIFEST_FILE_NAME if ga_latest_directory else None
        ),
        "source_qc_evidence_manifest_sha256": (
            qc_latest_directory / MANIFEST_FILE_NAME if qc_latest_directory else None
        ),
    }
    for manifest_key, source_manifest_path in source_manifest_specs.items():
        if source_manifest_path is not None and source_manifest_path.exists():
            manifest[manifest_key] = sha256_of_file(
                input_file_path=source_manifest_path
            )
            try:
                manifest[manifest_key.replace("_sha256", "")] = read_json_file(
                    input_file_path=source_manifest_path
                )
            except Exception:
                pass

    return manifest


def create_pc_biology_interpretation_snapshot(
    *,
    project_root: PathLike,
    excel_file_path: PathLike,
    pca_kmeans_latest_directory: PathLike,
    output_root_directory: PathLike = DEFAULT_OUTPUT_ROOT,
    ga_latest_directory: PathLike | None = None,
    qc_latest_directory: PathLike | None = None,
    excel_sheet_name: str = DEFAULT_EXCEL_SHEET_NAME,
    permutation_count: int = DEFAULT_PERMUTATION_COUNT,
    random_state: int = DEFAULT_RANDOM_STATE,
    minimum_category_count: int = 5,
    update_latest_directory: bool = True,
    verbose: bool = True,
) -> PcBiologyInterpretationResult:
    resolved_project_root = _as_path(project_root).resolve()
    resolved_excel_file_path = (resolved_project_root / excel_file_path).resolve()
    resolved_pca_kmeans_latest_directory = (
        resolved_project_root / pca_kmeans_latest_directory
    ).resolve()
    resolved_output_root_directory = (
        resolved_project_root / output_root_directory
    ).resolve()
    resolved_ga_latest_directory = (
        (resolved_project_root / ga_latest_directory).resolve()
        if ga_latest_directory is not None
        else None
    )
    resolved_qc_latest_directory = (
        (resolved_project_root / qc_latest_directory).resolve()
        if qc_latest_directory is not None
        else None
    )
    if verbose:
        print("Building integrated PC/biology table...")
    integrated_dataframe, integration_metadata = build_integrated_pc_biology_table(
        excel_file_path=resolved_excel_file_path,
        excel_sheet_name=excel_sheet_name,
        pca_kmeans_latest_directory=resolved_pca_kmeans_latest_directory,
        ga_latest_directory=resolved_ga_latest_directory,
        qc_latest_directory=resolved_qc_latest_directory,
    )

    if verbose:
        print("Computing PC/variable association tests...")
    associations_dataframe = compute_pc_biology_associations(
        integrated_dataframe=integrated_dataframe,
        permutation_count=permutation_count,
        random_state=random_state,
        minimum_category_count=minimum_category_count,
    )

    if verbose:
        print("Computing AG/KMeans cluster enrichment tests...")
    cluster_enrichment_dataframe = compute_cluster_biology_enrichment(
        integrated_dataframe=integrated_dataframe,
        minimum_category_count=minimum_category_count,
    )

    if verbose:
        print("Rendering markdown and HTML report...")
    summary_markdown = build_pc_axis_summary_markdown(
        associations_dataframe=associations_dataframe,
        cluster_enrichment_dataframe=cluster_enrichment_dataframe,
        integrated_dataframe=integrated_dataframe,
    )
    html_report = build_pc_biology_html_report(
        integrated_dataframe=integrated_dataframe,
        associations_dataframe=associations_dataframe,
        cluster_enrichment_dataframe=cluster_enrichment_dataframe,
    )

    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query="pc_biology_interpretation",
    )
    snapshot_directory = (
        resolved_output_root_directory / "snapshots" / snapshot_directory_name
    )
    if snapshot_directory.exists():
        shutil.rmtree(snapshot_directory)
    snapshot_directory.mkdir(parents=True, exist_ok=True)

    integrated_table_file_path = snapshot_directory / INTEGRATED_TABLE_FILE_NAME
    associations_file_path = snapshot_directory / ASSOCIATIONS_FILE_NAME
    cluster_enrichment_file_path = snapshot_directory / CLUSTER_ENRICHMENT_FILE_NAME
    summary_file_path = snapshot_directory / SUMMARY_FILE_NAME
    plot_html_file_path = snapshot_directory / PLOT_HTML_FILE_NAME
    manifest_file_path = snapshot_directory / MANIFEST_FILE_NAME

    _write_dataframe_csv_atomic(
        dataframe=integrated_dataframe,
        output_file_path=integrated_table_file_path,
    )
    _write_dataframe_csv_atomic(
        dataframe=associations_dataframe,
        output_file_path=associations_file_path,
    )
    _write_dataframe_csv_atomic(
        dataframe=cluster_enrichment_dataframe,
        output_file_path=cluster_enrichment_file_path,
    )
    _write_text_atomic(text=summary_markdown, output_file_path=summary_file_path)
    _write_html_atomic(html_text=html_report, output_file_path=plot_html_file_path)

    manifest_payload = _build_manifest(
        snapshot_created_at_utc=snapshot_created_at_utc,
        immutable_snapshot_directory_name=snapshot_directory_name,
        immutable_snapshot_relative_path=str(
            Path("snapshots") / snapshot_directory_name
        ),
        excel_file_path=resolved_excel_file_path,
        excel_sheet_name=excel_sheet_name,
        pca_kmeans_latest_directory=resolved_pca_kmeans_latest_directory,
        ga_latest_directory=resolved_ga_latest_directory,
        qc_latest_directory=resolved_qc_latest_directory,
        integrated_table_file_path=integrated_table_file_path,
        associations_file_path=associations_file_path,
        cluster_enrichment_file_path=cluster_enrichment_file_path,
        summary_file_path=summary_file_path,
        plot_html_file_path=plot_html_file_path,
        integration_metadata=integration_metadata,
        permutation_count=permutation_count,
        random_state=random_state,
        minimum_category_count=minimum_category_count,
    )
    write_json_atomic(payload=manifest_payload, output_file_path=manifest_file_path)

    latest_directory = resolved_output_root_directory / "latest"
    if update_latest_directory:
        _replace_latest_directory(
            latest_directory=latest_directory,
            files_to_copy=[
                (integrated_table_file_path, INTEGRATED_TABLE_FILE_NAME),
                (associations_file_path, ASSOCIATIONS_FILE_NAME),
                (cluster_enrichment_file_path, CLUSTER_ENRICHMENT_FILE_NAME),
                (summary_file_path, SUMMARY_FILE_NAME),
                (plot_html_file_path, PLOT_HTML_FILE_NAME),
                (manifest_file_path, MANIFEST_FILE_NAME),
            ],
        )

    if verbose:
        print(f"Snapshot directory: {snapshot_directory}")
        if update_latest_directory:
            print(f"Latest directory: {latest_directory}")

    return PcBiologyInterpretationResult(
        snapshot_directory=snapshot_directory,
        latest_directory=latest_directory,
        manifest_file_path=manifest_file_path,
        integrated_table_file_path=integrated_table_file_path,
        associations_file_path=associations_file_path,
        cluster_enrichment_file_path=cluster_enrichment_file_path,
        summary_file_path=summary_file_path,
        plot_html_file_path=plot_html_file_path,
    )
