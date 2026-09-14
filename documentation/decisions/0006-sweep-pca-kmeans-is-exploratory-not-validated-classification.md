# ADR-0006 — Treat SWeeP → PCA → KMeans as exploratory analysis, not validated biological classification

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (pipeline PRs dated
2026-04)

## Context

The project represents protein sequences numerically with SWeeP, reduces those
vectors with PCA, and partitions the reduced space with KMeans.

These methods reveal geometric structure in sequence space, but they do not by
themselves assign biological meaning to a cluster.

## Decision

Treat SWeeP → PCA → KMeans as exploratory analysis only.

A cluster is not a biological class. Silhouette measures geometric
cohesion/separation of a partition; ARI measures agreement between partitions
with chance correction. Neither validates a biological interpretation.

Keep historical results tied to the dataset that produced them. In particular,
the local clustering summary based on 41,345 records must not be mixed with the
later 52,473-record annotation-enriched candidate set.

## Evidence

- PRs #9, #10, and #11 introduced the SWeeP, PCA/KMeans, and PCA-plot stages
- `documentation/pago_qc.md` explicitly describes the metadata-derived labels
  and downstream analyses as bounded, non-ground-truth evidence
- the historical local summary for the 41,345-record run recorded approximately
  `k = 3`, sampled silhouette ≈ 0.656, selected-init ARI minimum ≈ 0.998,
  selected-subsample ARI minimum ≈ 0.993, Cramér's V ≈ 0.386, and Spearman
  ρ(PC2, length) ≈ 0.878

## Rationale

The methods' definitions bound what they can establish. Geometric separation and
clustering stability are useful exploratory properties, but biological
classification requires independent biological evidence.

## Consequences

- downstream cluster interpretation must remain explicitly exploratory unless
  independently validated
- the strong association between PC2 and sequence length must be treated as a
  potential confounding signal
- historical results must retain their dataset provenance

## Limitations

No biological validation of the historical clusters has been demonstrated.
SWeeP is also an external dependency whose exact historical execution is not
fully reproduced by a clean clone.

## Validation

The historical run has internal clustering-stability evidence. That is not a
biological validation result.

## Related records

- PRs #9, #10, #11
- commits `2593376`, `99946c2`, `05661a0`
- `documentation/pago_qc.md`
- historical local `tmp/siepe_report_analysis/summary.json`

## Historical reconstruction note

- Directly demonstrated: merged PRs, `pago_qc.md`, and the historical summary
  artifact.
- Not recoverable: the original deliberation over alternative representations
  or clustering methods.
