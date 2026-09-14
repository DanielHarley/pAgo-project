# ADR-0007 — Require provenance, defensible labels, redundancy control, and a protected holdout before supervised ML

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09`

## Context

A supervised pAgo classifier depends on labels that are biologically defensible
and on an evaluation design that prevents information from leaking from
validation data into model development.

The earlier pipeline reached exploratory clustering from text-query candidates
and metadata-derived QC labels. Those inputs are useful for exploration but are
not sufficient ground truth for final supervised learning.

## Decision

Before training or evaluating a final supervised model, require:

1. traceable provenance for reference data and labels
2. explicit evidence supporting biological labels
3. redundancy/high-similarity control
4. development and evaluation partitions that keep high-similarity groups
   together
5. a protected holdout that does not influence model, threshold, or strategy
   selection

Do not treat metadata-derived QC labels or query-recall references as automatic
supervised ground truth.

## Evidence

- the Phase A audit explicitly warns against treating text recovery as biological
  identity
- APAZ and clade-reference curation notes record evidence tiers and unresolved
  cases instead of forcing convenient labels
- the historical ontology specifies BUILD/CALIBRATION/FINAL_HOLDOUT and
  tree/placement equivalents for future validation
- the reference-layer work introduced high-similarity grouping specifically to
  prevent near-identical sequences from crossing development/evaluation
  boundaries

## Rationale

The main scientific risk is a plausible but weakly supported label becoming
training truth and then being evaluated against data that are not genuinely
independent of development.

The full literature argument that originally motivated every part of this design
was not preserved as one versioned source. Where that rationale cannot be
reconstructed directly, it should not be invented.

## Consequences

- reference construction and label curation precede final supervised modelling
- final holdouts must not be consulted during development or calibration
- evaluation must report the actual evidence and partitioning contract rather
  than only a performance score
- a reproducible pipeline does not by itself make the labels scientifically
  valid

## Limitations

The historical record supports the methodological boundary, but not every part
of the broader ML-literature rationale has been reconstructed from primary
sources inside the repository.

## Related records

- `documentation/history/2026-08-30-phase-a-audit.md`
- historical ontology staging record
- APAZ and clade-reference curation notes
- historical reference-layer commits `3de5f1b` … `49651f1`

## Historical reconstruction note

- Directly demonstrated: the audit's cautions, curation practices, historical
  validation protocols, and reference-layer commit sequence.
- Inferred: these pieces form one methodological requirement that precedes final
  supervised ML.
- Not fully recoverable: the complete contemporaneous literature justification.
