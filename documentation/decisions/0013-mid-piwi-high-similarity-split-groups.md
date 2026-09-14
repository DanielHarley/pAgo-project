# ADR-0013 — Use MID-PIWI high-similarity groups as future partition units

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commits
`18c35e4`, `49651f1`, 2026-08-31)

## Context

The 1,002 MID-PIWI regions extracted from the Ryazansky catalog will eventually
feed tree-building and placement-evaluation partitions. Near-identical regions
should not be allowed to cross those partitions.

## Decision

Apply the same high-similarity grouping method used for APAZ to the 1,002
MID-PIWI regions:

- use the MID-PIWI region as the comparison unit
- join exact region duplicates first
- build connected components using the declared 90% identity / 80% coverage
  criterion
- treat those components as indivisible future partition units

Run lower-identity diagnostic passes only as exploratory similarity checks.
Do not promote a diagnostic cross-clade similarity component into phylogenetic
evidence, relabeling, or automatic quarantine.

## Evidence

Historical manifests and commits record:

- 1,002 MID-PIWI regions and full all-vs-all pair coverage
- 625 qualifying edges at the operational 90/80 criterion
- 701 split groups
- 697 groups eligible for later partitioning
- 556 singleton groups
- 145 multi-member groups
- maximum group size 23
- 0 resolved cross-clade groups at 90/80
- 4 unresolved-only groups
- 2 groups requiring manual review because of AfAgo/SiAgo ambiguity
- 14 groups containing aliases or exact duplicates
- 984 unique MID-PIWI region sequences

A cross-clade component appears only in a lower 50%-identity diagnostic pass.
The historical record explicitly keeps that observation diagnostic rather than
using it to rewrite the source clade label.

## Rationale

Using one high-similarity rule across APAZ and MID-PIWI keeps the leakage-control
method consistent. Keeping a whole high-similarity group within one future
partition reduces a specific near-duplicate leakage risk.

The absence of cross-clade groups at 90/80 does not establish statistical
independence, absence of evolutionary relatedness, or correctness of the clade
labels themselves.

## Consequences

- the 697 eligible groups provide candidate indivisible units for later
  tree-build/calibration/holdout design
- AfAgo and SiAgo remain unresolved/manual-review references rather than being
  forced into a partition label
- lower-threshold similarity structure remains diagnostic only

## Limitations

The 90/80 rule is a redundancy-control convention, not a phylogenetic method or
a definition of homology. Future partition quality still requires independent
validation of the tree/placement design.

## Validation

Historical verification records byte-identical group outputs across independent
MMseqs2 runs and checks the committed split-group manifests.

## Related records

- historical commits `18c35e4`, `49651f1`
- historical MID-PIWI split-group resources and manifest
- ADR-0009 for the shared high-similarity grouping method
- historical roadmap reference: B4.3

## Historical reconstruction note

- Directly demonstrated: historical manifests, commit messages, and grouping
  workflow.
- All quoted group counts come from the preserved historical manifest.
