from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score

DEFAULT_PCA_SVD_SOLVER = "randomized"
DEFAULT_PCA_RANDOM_STATE = 42
DEFAULT_KMEANS_N_INIT: int | str = "auto"
DEFAULT_SILHOUETTE_SAMPLE_SIZE = 8000
DEFAULT_SILHOUETTE_RANDOM_STATE = 42
DEFAULT_KMEANS_INITIALIZATION_REPEAT_COUNT = 6
DEFAULT_SUBSAMPLE_REPEAT_COUNT = 6
DEFAULT_SUBSAMPLE_FRACTION = 0.80
DEFAULT_SUBSAMPLE_RANDOM_STATE = 123
DEFAULT_MINIMUM_ACCEPTABLE_INIT_ARI_MIN = 0.70
DEFAULT_MINIMUM_ACCEPTABLE_SUBSAMPLE_ARI_MIN = 0.60


def get_effective_pca_component_count_grid(
    *,
    data_matrix: np.ndarray,
    requested_pca_component_count_grid: Sequence[int],
) -> tuple[int, ...]:
    maximum_valid_component_count = int(min(data_matrix.shape[0], data_matrix.shape[1]))
    effective_grid = tuple(
        int(requested_component_count)
        for requested_component_count in requested_pca_component_count_grid
        if 1 <= int(requested_component_count) <= maximum_valid_component_count
    )
    if not effective_grid:
        raise ValueError(
            "No valid PCA component counts remain after filtering by matrix shape "
            f"{data_matrix.shape}. Original grid was "
            f"{tuple(requested_pca_component_count_grid)}."
        )
    return effective_grid


def compute_pca_coordinates_and_model(
    *,
    data_matrix: np.ndarray,
    requested_pca_component_count: int,
    pca_svd_solver: str = DEFAULT_PCA_SVD_SOLVER,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
) -> tuple[np.ndarray, np.ndarray, PCA]:
    pca_model_object = PCA(
        n_components=int(requested_pca_component_count),
        svd_solver=pca_svd_solver,
        random_state=int(pca_random_state),
    )
    pca_coordinate_matrix = pca_model_object.fit_transform(data_matrix)
    pca_explained_variance_ratio_vector = pca_model_object.explained_variance_ratio_
    return pca_coordinate_matrix, pca_explained_variance_ratio_vector, pca_model_object


def compute_pairwise_ari_statistics_from_label_vectors(
    label_vectors: Sequence[np.ndarray],
) -> tuple[float, float, float, int]:
    if len(label_vectors) <= 1:
        return 1.0, 1.0, 1.0, 0

    pairwise_ari_values: list[float] = []
    for first_index in range(len(label_vectors)):
        for second_index in range(first_index + 1, len(label_vectors)):
            pairwise_ari_values.append(
                float(
                    adjusted_rand_score(
                        label_vectors[first_index],
                        label_vectors[second_index],
                    )
                )
            )

    return (
        float(np.mean(pairwise_ari_values)),
        float(np.min(pairwise_ari_values)),
        float(np.max(pairwise_ari_values)),
        int(len(pairwise_ari_values)),
    )


def compute_pairwise_ari_statistics_from_subsample_labelings(
    subsample_labelings: Sequence[tuple[np.ndarray, np.ndarray]],
) -> tuple[float, float, float, int, float]:
    if len(subsample_labelings) <= 1:
        return 1.0, 1.0, 1.0, 0, 0.0

    pairwise_ari_values: list[float] = []
    intersection_size_values: list[int] = []

    for first_index in range(len(subsample_labelings)):
        first_subsample_index_vector, first_subsample_label_vector = subsample_labelings[
            first_index
        ]

        for second_index in range(first_index + 1, len(subsample_labelings)):
            (
                second_subsample_index_vector,
                second_subsample_label_vector,
            ) = subsample_labelings[second_index]

            common_index_vector, first_positions, second_positions = np.intersect1d(
                first_subsample_index_vector,
                second_subsample_index_vector,
                return_indices=True,
            )
            if common_index_vector.size < 2:
                continue

            first_common_labels = first_subsample_label_vector[first_positions]
            second_common_labels = second_subsample_label_vector[second_positions]

            pairwise_ari_values.append(
                float(adjusted_rand_score(first_common_labels, second_common_labels))
            )
            intersection_size_values.append(int(common_index_vector.size))

    if not pairwise_ari_values:
        return 0.0, 0.0, 0.0, 0, 0.0

    return (
        float(np.mean(pairwise_ari_values)),
        float(np.min(pairwise_ari_values)),
        float(np.max(pairwise_ari_values)),
        int(len(pairwise_ari_values)),
        float(np.mean(intersection_size_values)) if intersection_size_values else 0.0,
    )


def fit_kmeans_and_compute_sampled_silhouette(
    *,
    coordinate_matrix_for_clustering: np.ndarray,
    requested_cluster_count_k: int,
    kmeans_random_state: int,
    silhouette_sample_size: int,
    silhouette_random_state: int,
    kmeans_n_init: int | str = DEFAULT_KMEANS_N_INIT,
) -> tuple[np.ndarray, float]:
    kmeans_model_object = KMeans(
        n_clusters=int(requested_cluster_count_k),
        random_state=int(kmeans_random_state),
        n_init=kmeans_n_init,
    )
    label_vector = kmeans_model_object.fit_predict(coordinate_matrix_for_clustering)
    total_point_count = int(coordinate_matrix_for_clustering.shape[0])
    effective_sample_size = min(int(silhouette_sample_size), total_point_count)
    unique_cluster_count = int(np.unique(label_vector).size)
    if unique_cluster_count < 2 or effective_sample_size < 2:
        return label_vector, float("nan")

    sampled_silhouette_value = float(
        silhouette_score(
            coordinate_matrix_for_clustering,
            label_vector,
            sample_size=effective_sample_size,
            random_state=int(silhouette_random_state),
        )
    )
    return label_vector, sampled_silhouette_value


def run_pca_kmeans_stability_grid_sweep(
    *,
    data_matrix: np.ndarray,
    pca_component_count_grid: Sequence[int],
    kmeans_cluster_count_grid: Sequence[int],
    kmeans_initialization_repeat_count: int = DEFAULT_KMEANS_INITIALIZATION_REPEAT_COUNT,
    subsample_repeat_count: int = DEFAULT_SUBSAMPLE_REPEAT_COUNT,
    subsample_fraction: float = DEFAULT_SUBSAMPLE_FRACTION,
    silhouette_sample_size: int = DEFAULT_SILHOUETTE_SAMPLE_SIZE,
    pca_svd_solver: str = DEFAULT_PCA_SVD_SOLVER,
    pca_random_state: int = DEFAULT_PCA_RANDOM_STATE,
    kmeans_n_init: int | str = DEFAULT_KMEANS_N_INIT,
    silhouette_random_state: int = DEFAULT_SILHOUETTE_RANDOM_STATE,
    subsample_random_state: int = DEFAULT_SUBSAMPLE_RANDOM_STATE,
    progress_callback: Callable[[str], None] | None = None,
) -> tuple[pd.DataFrame, dict[int, np.ndarray]]:
    stability_grid_row_list: list[dict[str, Any]] = []
    pca_coordinate_cache_by_component_count: dict[int, np.ndarray] = {}

    for requested_pca_component_count in pca_component_count_grid:
        (
            pca_coordinate_matrix_m,
            pca_explained_variance_ratio_vector,
            _,
        ) = compute_pca_coordinates_and_model(
            data_matrix=data_matrix,
            requested_pca_component_count=int(requested_pca_component_count),
            pca_svd_solver=pca_svd_solver,
            pca_random_state=int(pca_random_state),
        )
        explained_variance_fraction_for_m = float(
            np.sum(pca_explained_variance_ratio_vector)
        )
        pca_coordinate_cache_by_component_count[int(requested_pca_component_count)] = (
            pca_coordinate_matrix_m
        )

        total_point_count = int(pca_coordinate_matrix_m.shape[0])
        subsample_point_count = int(round(float(subsample_fraction) * total_point_count))

        for requested_cluster_count_k in kmeans_cluster_count_grid:
            if int(requested_cluster_count_k) < 2:
                continue
            if int(requested_cluster_count_k) >= total_point_count:
                continue

            (
                baseline_label_vector,
                baseline_sampled_silhouette_value,
            ) = fit_kmeans_and_compute_sampled_silhouette(
                coordinate_matrix_for_clustering=pca_coordinate_matrix_m,
                requested_cluster_count_k=int(requested_cluster_count_k),
                kmeans_random_state=42,
                silhouette_sample_size=int(silhouette_sample_size),
                silhouette_random_state=int(silhouette_random_state),
                kmeans_n_init=kmeans_n_init,
            )

            label_vectors_across_initializations: list[np.ndarray] = [
                baseline_label_vector
            ]
            for seed_value in range(1, int(kmeans_initialization_repeat_count)):
                kmeans_model_object = KMeans(
                    n_clusters=int(requested_cluster_count_k),
                    random_state=int(seed_value),
                    n_init=kmeans_n_init,
                )
                label_vectors_across_initializations.append(
                    kmeans_model_object.fit_predict(pca_coordinate_matrix_m)
                )

            (
                init_ari_mean,
                init_ari_min,
                init_ari_max,
                init_ari_pair_count,
            ) = compute_pairwise_ari_statistics_from_label_vectors(
                label_vectors_across_initializations
            )

            subsample_labeling_list: list[tuple[np.ndarray, np.ndarray]] = []
            for subsample_iteration_index in range(int(subsample_repeat_count)):
                iteration_rng = np.random.default_rng(
                    int(subsample_random_state)
                    + int(subsample_iteration_index) * 10007
                )
                subsample_index_vector = iteration_rng.choice(
                    total_point_count,
                    size=subsample_point_count,
                    replace=False,
                )
                subsample_index_vector_sorted = np.sort(subsample_index_vector)
                subsample_coordinate_matrix = pca_coordinate_matrix_m[
                    subsample_index_vector_sorted, :
                ]
                kmeans_model_object_subsample = KMeans(
                    n_clusters=int(requested_cluster_count_k),
                    random_state=42,
                    n_init=kmeans_n_init,
                )
                subsample_label_vector = kmeans_model_object_subsample.fit_predict(
                    subsample_coordinate_matrix
                )
                subsample_labeling_list.append(
                    (subsample_index_vector_sorted, subsample_label_vector)
                )

            (
                subsample_ari_mean,
                subsample_ari_min,
                subsample_ari_max,
                subsample_ari_pair_count,
                mean_intersection_size,
            ) = compute_pairwise_ari_statistics_from_subsample_labelings(
                subsample_labeling_list
            )

            stability_grid_row_list.append(
                {
                    "pca_component_count": int(requested_pca_component_count),
                    "variance_explained_fraction": float(
                        explained_variance_fraction_for_m
                    ),
                    "k": int(requested_cluster_count_k),
                    "silhouette_best_sampled": float(
                        baseline_sampled_silhouette_value
                    ),
                    "init_ari_mean": float(init_ari_mean),
                    "init_ari_min": float(init_ari_min),
                    "init_ari_max": float(init_ari_max),
                    "init_ari_pair_count": int(init_ari_pair_count),
                    "subsample_ari_mean": float(subsample_ari_mean),
                    "subsample_ari_min": float(subsample_ari_min),
                    "subsample_ari_max": float(subsample_ari_max),
                    "subsample_ari_pair_count": int(subsample_ari_pair_count),
                    "subsample_fraction": float(subsample_fraction),
                    "subsample_point_count": int(subsample_point_count),
                    "mean_pairwise_intersection_size": float(
                        mean_intersection_size
                    ),
                }
            )

            if progress_callback is not None:
                progress_callback(
                    "m="
                    f"{int(requested_pca_component_count):3d}, "
                    f"k={int(requested_cluster_count_k):2d} | "
                    f"sil={baseline_sampled_silhouette_value:.4f} | "
                    f"initARI(min/mean)={init_ari_min:.3f}/{init_ari_mean:.3f} | "
                    f"subARI(min/mean)={subsample_ari_min:.3f}/{subsample_ari_mean:.3f}"
                )

    stability_grid_dataframe = pd.DataFrame(stability_grid_row_list)
    if stability_grid_dataframe.empty:
        raise RuntimeError("The PCA/KMeans stability sweep did not produce any rows.")

    stability_grid_dataframe = stability_grid_dataframe.sort_values(
        ["pca_component_count", "k"]
    ).reset_index(drop=True)
    return stability_grid_dataframe, pca_coordinate_cache_by_component_count


def choose_best_m_and_k_from_stability_grid(
    *,
    stability_grid_dataframe: pd.DataFrame,
    minimum_acceptable_init_ari_min: float = DEFAULT_MINIMUM_ACCEPTABLE_INIT_ARI_MIN,
    minimum_acceptable_subsample_ari_min: float = DEFAULT_MINIMUM_ACCEPTABLE_SUBSAMPLE_ARI_MIN,
) -> tuple[pd.Series, str, pd.DataFrame]:
    if stability_grid_dataframe.empty:
        raise ValueError("stability_grid_dataframe must not be empty.")

    dataframe_copy = stability_grid_dataframe.copy()
    dataframe_copy["silhouette_best_sampled_filled"] = dataframe_copy[
        "silhouette_best_sampled"
    ].fillna(-1.0)
    dataframe_copy["init_ari_min_clipped"] = dataframe_copy["init_ari_min"].clip(
        lower=0.0
    )
    dataframe_copy["subsample_ari_min_clipped"] = dataframe_copy[
        "subsample_ari_min"
    ].clip(lower=0.0)
    dataframe_copy["composite_score_silhouette_times_min_ari"] = (
        dataframe_copy["silhouette_best_sampled_filled"]
        * dataframe_copy["init_ari_min_clipped"]
        * dataframe_copy["subsample_ari_min_clipped"]
    )

    threshold_filtered_dataframe = dataframe_copy[
        (
            dataframe_copy["init_ari_min"]
            >= float(minimum_acceptable_init_ari_min)
        )
        & (
            dataframe_copy["subsample_ari_min"]
            >= float(minimum_acceptable_subsample_ari_min)
        )
    ]
    if not threshold_filtered_dataframe.empty:
        best_row = threshold_filtered_dataframe.sort_values(
            [
                "silhouette_best_sampled_filled",
                "variance_explained_fraction",
                "pca_component_count",
            ],
            ascending=[False, False, True],
        ).iloc[0]
        return best_row, "threshold_pass_then_max_silhouette", dataframe_copy

    best_row = dataframe_copy.sort_values(
        [
            "composite_score_silhouette_times_min_ari",
            "silhouette_best_sampled_filled",
            "variance_explained_fraction",
        ],
        ascending=[False, False, False],
    ).iloc[0]
    return best_row, "fallback_max_composite_score", dataframe_copy


def build_cluster_assignment_dataframe(
    *,
    aligned_metadata_dataframe: pd.DataFrame,
    cluster_label_vector: np.ndarray,
    pca_coordinate_matrix: np.ndarray,
    projection_dimension_for_table: int = 3,
) -> pd.DataFrame:
    dataframe_copy = aligned_metadata_dataframe.copy()
    dataframe_copy["cluster_label"] = cluster_label_vector.astype(int)
    table_projection_dimension = min(
        int(projection_dimension_for_table),
        int(pca_coordinate_matrix.shape[1]),
    )
    for component_index in range(table_projection_dimension):
        dataframe_copy[f"pc{component_index + 1}"] = pca_coordinate_matrix[
            :, component_index
        ]
    return dataframe_copy
