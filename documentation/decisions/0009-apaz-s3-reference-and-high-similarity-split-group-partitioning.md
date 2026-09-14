# ADR-0009 — Partition APAZ references by high-similarity split groups

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commits
`4d8aa99`, `620ff85`, `45cc637`, 2026-08-31)

## Context

Ryazansky Data Set S3 contains 481 representative APAZ domains derived from a
larger set by a reported 90% UCLUST reduction. Building an APAZ classifier
requires development and evaluation partitions that do not place near-identical
references on opposite sides of the split.

## Decision

Use a **split group** as the indivisible partition unit.

A split group is a connected component of an all-vs-all sequence-similarity
graph where an edge requires identity ≥ 0.90 and coverage ≥ 0.80, with identity
and coverage recomputed from the aligned strings and exact-sequence duplicates
joined first.

This 90/80 criterion is an operational high-similarity/redundancy-control rule,
not a definition of homology.

For the 481 APAZ representatives, the historical derivation produced 460 split
groups. Assign whole groups to FINAL_HOLDOUT first, then CALIBRATION, then BUILD,
using the deterministic partition policy. Keep HisG (`PF01634`) and EIIB
(`PF00367`) hard-negative families out of BUILD.

Once selected, the final holdout must not be redistributed in response to HMM
performance.

## Evidence

- historical APAZ curation notes describe the 481 S3 representatives, the v1
  leakage problem, the 90/80 split-group method, and the holdout-first partition
  policy
- historical manifests record 481 sequences, 460 split groups, 20 multi-member
  groups, and maximum group size 3
- a documented v1 leak placed two 100%-identical/0.816-coverage references in
  different development partitions
- commit `45cc637` records that the known v1 cross-partition high-similarity
  leaks were resolved and that v1 had not been used for predictive evaluation

## Rationale

Keeping a whole high-similarity component within one partition removes the
specific leakage mode in which near-identical sequences appear in both
development and evaluation.

Freezing the final holdout before model evaluation prevents later
performance-driven redistribution.

## Consequences

- APAZ HMMs are constructed from BUILD references only
- high-similarity grouping reduces an important leakage risk but does not prove
  complete statistical independence between partitions
- evaluation against HisG/EIIB tests generalisation to new sequences from known
  negative families, not necessarily to unseen negative families

## Limitations

The exact historical procedure used to enumerate every reported v1 leakage case
was not fully preserved. One worked example and the final count survive.

The origin of the exact 90/80 threshold is only partially recoverable.

## Supersedes

The earlier accession-hash partition unit that allowed high-similarity groups to
cross partition boundaries.

## Validation

Historical reference verification checks partition and seed invariants.
Historical regeneration also reproduced the split groups deterministically and
recorded byte-identical results across independent MMseqs2 runs.

## Related records

- historical commits `4d8aa99`, `620ff85`, `45cc637`
- historical `src/pago_pipeline/resources/apaz_seed/**`
- historical roadmap reference: B2

## Historical reconstruction note

- Directly demonstrated: curation notes, manifests, commits, and the documented
  leakage example.
- Not recoverable: the full historical enumeration procedure for every v1
  leakage case.
