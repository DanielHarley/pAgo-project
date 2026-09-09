# ADR-0013 — MID-PIWI high-similarity split groups over the Ryazansky catalog

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — Phase B (milestone B4.3) on `(feat)-siepe-ready-project`.

## Epistemic status

`SUPPORTED` — deterministic derivation, byte-identical across two independent
runs.

## Type

`TECHNICAL + SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commits
`18c35e4`, `49651f1`, 2026-08-31)

## Context

The 1002 MID-PIWI regions from the Ryazansky catalog (ADR-0011) will eventually
be split between tree-build and placement-evaluation partitions. Near-identical
regions must not span those partitions.

## Problem

Identify the high-similarity structure among the MID-PIWI regions with the same
method used for APAZ (ADR-0009), and check whether high similarity ever crosses
a resolved clade boundary.

## Evidence

- Staging commit `49651f1` and
  `src/pago_pipeline/resources/clade_seed/split_groups/midpiwi/manifest.json`:
  - `n_regions = 1002`; `all_vs_all_observed_unordered_pairs = 501501`
    (= C(1002,2), 100 % combinatorial coverage; `pairs_raw_tsv_row_count =
    1004004`, both orientations). An **unordered pair** treats A–B and B–A as
    the same pair.
  - `pairs.filtered.tsv` has **626** rows = 1 header + **625 edges** at the
    90/80 criterion.
  - `n_split_groups = 701`; `partition_eligible_groups = 697`
    (LONG_A 165 / LONG_B 139 / SHORT 393); `n_singletons = 556`;
    `n_multi_member_groups = 145`; `max_group_size = 23`.
  - `label_conflict_groups = 0` at 90/80 (0 cross-clade groups);
    `unresolved_only_groups = 4`; `groups_requiring_manual_review = 2`
    (AfAgo, SiAgo — not partition-eligible);
    `groups_with_alias_or_exact_duplicate = 14`; `n_unique_regions = 984`.
- Staging commit `49651f1` message: a cross-clade component appears **only** at
  a 50 % identity **diagnostic** threshold (one 38-member group,
  37 LONG_B + `WP_076701676.1` LONG_A). "Recorded as diagnostic evidence only —
  no quarantine, no exclusion; `WP_076701676.1` stays LONG_A per Table S1."
- Same MMseqs2 workflow as ADR-0009 (`src/pago_pipeline/apaz_split_groups.py`).

## Previous state

No MID-PIWI similarity analysis.

## Decision

- Run the ADR-0009 MMseqs2 90/80 connected-component workflow on the 1002
  MID-PIWI regions (clustering unit = the region, not the full protein), with a
  prior union of exact `midpiwi_region_sha256` duplicates.
- Record the result as high-similarity split groups: 701 groups, 697
  partition-eligible, **0 label-conflict groups at 90/80**.
- Run diagnostic passes at 1.00 / 0.90 / 0.70 / 0.50 identity (coverage 0.80)
  and record the members of any cross-clade group **as diagnostic evidence
  only** — never as a partition rule, never as a leakage claim, and **never as
  phylogenetic evidence**.

## Rationale

Stated in the commits: using the same, already-validated workflow keeps the
methodology consistent. Zero resolved cross-clade groups at the operational
90/80 threshold supports using these high-similarity split groups as indivisible
units to prevent 90/80-defined high-similarity cross-partition leakage. It does
**not** establish statistical independence, absence of evolutionary
relatedness, or absence of other forms of leakage. The 50 % cross-clade
component is a similarity observation only, with no phylogenetic
interpretation.

## Alternatives considered

Treating the 50 % cross-clade component as a reason to quarantine or exclude
`WP_076701676.1` (rejected — it stays LONG_A per Table S1; the observation is
diagnostic). `HISTORICAL_EVIDENCE_INSUFFICIENT` for other thresholds or units.

## Consequences

- The 697 partition-eligible groups are the indivisible units for a future
  TREE_BUILD / PLACEMENT_CALIBRATION / PLACEMENT_HOLDOUT split (not yet built).
- AfAgo and SiAgo remain quarantined (not partition-eligible).

## Limitations

Not integrated. "0 cross-clade at 90/80" does not imply statistical independence
between future partitions. The diagnostic 50 % component carries no
phylogenetic weight.

## Supersedes

Nothing.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`18c35e4`, `49651f1`.

## Related data / artifacts

`src/pago_pipeline/resources/clade_seed/split_groups/midpiwi/**` (staging):
`split_groups.csv`, `groups.csv`, `pairs.filtered.tsv`,
`diagnostic_thresholds.json`, `manifest.json`, `mmseqs_search.log`.

## Validation

`scripts/verify_reference_data.py --scope clade` (staging);
`tests/test_clade_midpiwi_split_groups.py`; byte-identical across two
independent MMseqs2 runs (`split_groups.csv` SHA-256 `329b98c8…`).

## Scientific impact

Provides the indivisible units for a future clade-placement split and shows no
high-similarity clade conflict at 90/80.

## Data impact

Adds committed split-group resources. No candidate-protein label assigned.

## Reproducibility impact

Positive: deterministic, byte-identical across independent runs, offline
verifier.

## Historical reconstruction note

- Directly demonstrated: the manifest counts, the commit messages, the split
  workflow.
- Inferred: none material.
- Not recoverable: none — all quoted numbers are in the committed manifest.
