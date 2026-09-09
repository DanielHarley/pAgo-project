# ADR-0006 — SWeeP → PCA → KMeans is exploratory, not validated biological classification

## Lifecycle status

`CURRENT`

## Integration status

`INTEGRATED` — the pipeline (notebooks 05–07, PRs #9–#11) is in `master`. The
*framing* — exploratory only — is the decision recorded here.

## Epistemic status

`PARTIALLY_SUPPORTED` — the methods do what they do; they do not, by
construction, assign biological meaning to a cluster.

## Type

`SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (pipeline PRs dated
2026-04)

## Context

The project produces a numeric representation of the retained proteome and
explores its structure.

- **SWeeP** is an alignment-free vector representation of biological sequences.
- **PCA** (Principal Component Analysis) is a linear dimensionality-reduction
  technique that describes axes of variance.
- **KMeans** partitions points into *k* groups by minimising within-cluster
  squared distance to centroids.

## Problem

Prevent an exploratory clustering from being read as a validated biological
classification.

## Evidence

- PRs #9 (`2593376`, SWeeP), #10 (`99946c2`, PCA/KMeans), #11 (`05661a0`,
  PCA plot) — merged to `master`.
- `documentation/pago_qc.md` "Known Boundaries"; the repository `README.md`
  ("exploratory PCA and KMeans analyses").
- `tmp/siepe_report_analysis/summary.json` (local historical artifact,
  gitignored) — a run over the **earlier** filtered dataset of **41,345**
  records: `selected_cluster_count_k = 3`, `final_sampled_silhouette_value ≈
  0.656`, `selected_init_ari_min ≈ 0.998`, `selected_subsample_ari_min ≈
  0.993`, Cramér's V (cluster vs primary label) ≈ 0.386, Spearman ρ(PC2, length)
  ≈ 0.878. Its `dataset_counts` (2,629 / 6,703 / 15,056 / 16,957; total 41,345)
  match `documentation/pago_qc.md`.

## Previous state

No explicit statement that the clustering is not a classification.

## Decision

- SWeeP → PCA → KMeans is an **exploratory analysis of the sequence space**. A
  cluster is not a biological class. Silhouette measures partition cohesion /
  separation; ARI measures agreement between partitions with chance correction;
  neither validates biology.
- The `summary.json` numbers belong to the **earlier 41,345-record dataset**
  (a different query and snapshot) and must **not** be mixed with the Phase A
  52,473-record candidate set (ADR-0004).

## Rationale

Stated in `pago_qc.md` "Known Boundaries": the labels and any structure found
are not equivalent to external domain validation (HMMER / Pfam / CDD /
InterProScan / motif scans). The methods' own definitions bound what they can
claim.

## Alternatives considered

`HISTORICAL_EVIDENCE_INSUFFICIENT` — no design document weighs alternative
representations or clustering methods for this stage. SWeeP is developed with
UFPR participation (De Pierri et al. 2020) and chosen for scalability; the
functional interpretation of pAgo clusters remains a hypothesis.

## Consequences

- Downstream analysis (Phase E in the roadmap) is scoped as exploratory over
  (a) the pre-filtered proteome and (b) classified pAgos — not as classification
  itself.
- The 41,345 vs 52,473 distinction must be stated wherever these numbers appear.

## Limitations

Length is a visible confounder (Spearman ρ(PC2, length) ≈ 0.878 in the 41,345
run). No biological validation of any cluster has been done.

## Supersedes

Nothing formal.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

PR #9, #10, #11.

## Related commits

`2593376`, `99946c2`, `05661a0`.

## Related data / artifacts

`data/03-features/sweep_genes/**`, `data/03-features/pca_kmeans/**` (local,
gitignored); `tmp/siepe_report_analysis/summary.json` (local, gitignored,
historical).

## Validation

The clustering has internal-stability evidence (ARI) for the 41,345 run. It has
**no** biological validation, and none is claimed.

## Scientific impact

Bounds the claims of the exploratory analysis; separates the historical 41,345
dataset from the Phase A 52,473 set.

## Data impact

None on labels or partitions.

## Reproducibility impact

The `summary.json` run is a local artifact; a clean clone does not reproduce it.

## Historical reconstruction note

- Directly demonstrated: the merged PRs, `pago_qc.md`, `summary.json`.
- Inferred: that the exploratory framing is a deliberate decision rather than an
  omission (supported by the explicit "Known Boundaries" text).
- Not recoverable: the deliberation over representation / clustering choices.
