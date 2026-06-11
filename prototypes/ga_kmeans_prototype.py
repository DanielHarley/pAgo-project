from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.metrics import silhouette_score


@dataclass(frozen=True)
class Evaluation:
    k: int
    silhouette: float
    fitness: float
    labels: np.ndarray
    inertia: float


@dataclass(frozen=True)
class GeneticKMeansResult:
    best_k: int
    best_silhouette: float
    best_fitness: float
    labels: np.ndarray
    history: pd.DataFrame
    evaluated_candidates: pd.DataFrame


def _validate(matrix: np.ndarray, minimum_k: int, maximum_k: int) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] < 3 or matrix.shape[1] < 1:
        raise ValueError("data_matrix must be a finite 2D matrix with at least 3 rows.")
    if not np.isfinite(matrix).all():
        raise ValueError("data_matrix must contain only finite values.")
    if minimum_k < 2 or maximum_k < minimum_k or maximum_k >= matrix.shape[0]:
        raise ValueError("Require 2 <= minimum_k <= maximum_k < number of samples.")
    return matrix


def run_genetic_kmeans_search(
    *,
    data_matrix: np.ndarray,
    minimum_k: int = 2,
    maximum_k: int = 15,
    population_size: int = 20,
    generation_count: int = 15,
    elite_count: int = 2,
    tournament_size: int = 3,
    mutation_probability: float = 0.25,
    complexity_penalty_weight: float = 0.0,
    silhouette_sample_size: int | None = 8000,
    random_state: int = 42,
    kmeans_random_state: int = 42,
    silhouette_random_state: int = 42,
    kmeans_n_init: int | str = "auto",
) -> GeneticKMeansResult:
    """Search for the KMeans cluster count with a one-gene integer-coded GA."""
    matrix = _validate(data_matrix, int(minimum_k), int(maximum_k))
    if population_size < 2 or not 1 <= elite_count < population_size:
        raise ValueError("Require population_size >= 2 and 1 <= elite_count < population_size.")
    if generation_count < 1 or tournament_size < 2:
        raise ValueError("Require generation_count >= 1 and tournament_size >= 2.")
    if not 0.0 <= mutation_probability <= 1.0 or complexity_penalty_weight < 0.0:
        raise ValueError("Invalid mutation probability or complexity penalty.")

    rng = np.random.default_rng(int(random_state))
    population = rng.integers(minimum_k, maximum_k + 1, size=population_size)
    cache: dict[int, Evaluation] = {}
    history: list[dict[str, Any]] = []

    def evaluate(k: int) -> Evaluation:
        k = int(k)
        if k not in cache:
            model = KMeans(n_clusters=k, random_state=kmeans_random_state, n_init=kmeans_n_init)
            labels = model.fit_predict(matrix)
            sample_size = None if silhouette_sample_size is None else min(silhouette_sample_size, len(matrix))
            silhouette = float(
                silhouette_score(matrix, labels, sample_size=sample_size, random_state=silhouette_random_state)
            )
            span = max(1, maximum_k - minimum_k)
            fitness = silhouette - complexity_penalty_weight * ((k - minimum_k) / span)
            cache[k] = Evaluation(k, silhouette, float(fitness), labels.astype(int), float(model.inertia_))
        return cache[k]

    def select(fitness_by_k: dict[int, float]) -> int:
        contenders = population[rng.integers(0, len(population), size=tournament_size)]
        return int(max(contenders, key=lambda k: fitness_by_k[int(k)]))

    for generation in range(generation_count):
        evaluations = [evaluate(k) for k in population]
        fitness_by_k = {item.k: item.fitness for item in evaluations}
        ranked = sorted(population.tolist(), key=lambda k: fitness_by_k[int(k)], reverse=True)
        best = evaluate(ranked[0])
        history.append(
            {
                "generation": generation,
                "best_k": best.k,
                "best_silhouette": best.silhouette,
                "best_fitness": best.fitness,
                "mean_fitness": float(np.mean([item.fitness for item in evaluations])),
                "unique_k_count": int(np.unique(population).size),
            }
        )
        next_population = [int(k) for k in ranked[:elite_count]]
        while len(next_population) < population_size:
            first_parent, second_parent = select(fitness_by_k), select(fitness_by_k)
            child = first_parent if rng.random() < 0.5 else second_parent
            if rng.random() < mutation_probability:
                child = int(np.clip(child + int(rng.choice([-1, 1])), minimum_k, maximum_k))
            next_population.append(child)
        population = np.asarray(next_population, dtype=int)

    best = max(cache.values(), key=lambda item: item.fitness)
    evaluated_candidates = pd.DataFrame(
        [{"k": item.k, "silhouette": item.silhouette, "fitness": item.fitness, "inertia": item.inertia}
         for item in sorted(cache.values(), key=lambda item: item.k)]
    )
    return GeneticKMeansResult(
        best_k=best.k,
        best_silhouette=best.silhouette,
        best_fitness=best.fitness,
        labels=best.labels.copy(),
        history=pd.DataFrame(history),
        evaluated_candidates=evaluated_candidates,
    )


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Minimal GA search for a KMeans cluster count.")
    parser.add_argument("--input", type=Path, default=None, help="Optional .npy matrix. Synthetic data is used when omitted.")
    parser.add_argument("--minimum-k", type=int, default=2)
    parser.add_argument("--maximum-k", type=int, default=12)
    parser.add_argument("--output-directory", type=Path, default=Path("data/04-analysis/ga_kmeans_prototype"))
    args = parser.parse_args()

    matrix = np.load(args.input) if args.input else make_blobs(n_samples=600, centers=4, cluster_std=0.70, random_state=42)[0]
    result = run_genetic_kmeans_search(data_matrix=matrix, minimum_k=args.minimum_k, maximum_k=args.maximum_k)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    result.history.to_csv(args.output_directory / "ga_history.csv", index=False)
    result.evaluated_candidates.to_csv(args.output_directory / "evaluated_candidates.csv", index=False)
    pd.DataFrame({"cluster_label": result.labels}).to_csv(args.output_directory / "cluster_labels.csv", index=False)
    print(f"best_k={result.best_k} | silhouette={result.best_silhouette:.6f} | fitness={result.best_fitness:.6f}")


if __name__ == "__main__":
    main()
