# GA + KMeans clustering prototype

This directory contains an intentionally minimal proof of concept for the DS804 genetic-algorithm clustering assignment.

## Objective

The prototype uses a one-gene integer-coded genetic algorithm. Each individual stores a candidate number of clusters `k`. For each candidate, the program fits `sklearn.cluster.KMeans`, computes the silhouette score, and uses that score as the default fitness value.

The resulting search is pedagogical. For a small interval of candidate `k` values, an exhaustive sweep remains simpler and computationally stronger. This prototype exists to establish the genetic-algorithm mechanics before extending the chromosome and the fitness function.

## Genetic operators

* Representation: one integer gene, `k`
* Selection: tournament selection
* Crossover: inheritance of `k` from one of two selected parents
* Mutation: increment or decrement of `k` by one unit
* Elitism: direct preservation of the best individuals
* Evaluation cache: repeated values of `k` reuse previous KMeans results

## Fitness

By default:

```text
fitness = silhouette_score
```

An optional complexity penalty can discourage excessive fragmentation:

```text
fitness = silhouette_score - complexity_penalty_weight * normalized_k
```

## Run with synthetic data

From the repository root:

```powershell
python prototypes/ga_kmeans_prototype.py
```

The default synthetic dataset contains four Gaussian clusters. Expected output is similar to:

```text
best_k=4 | silhouette=0.852434 | fitness=0.852434
```

## Run with an existing numeric matrix

Provide a NumPy `.npy` matrix with shape `(n_samples, n_features)`:

```powershell
python prototypes/ga_kmeans_prototype.py `
  --input path/to/matrix.npy `
  --minimum-k 2 `
  --maximum-k 15 `
  --output-directory data/04-analysis/ga_kmeans_prototype
```

The script writes:

* `ga_history.csv`
* `evaluated_candidates.csv`
* `cluster_labels.csv`

## Immediate next extension

The next prototype should evolve a richer chromosome, such as `(pca_component_count, k)`, and incorporate cluster stability across KMeans initializations and subsamples. The existing module `src/pago_pipeline/pca_kmeans.py` already implements silhouette and ARI-based stability measurements that can be reused.
