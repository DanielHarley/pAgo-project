# ADR-0009 — APAZ reference from Ryazansky Data Set S3; high-similarity split-group partitioning (v2), superseding v1

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — Phase B (milestone B2) on `(feat)-siepe-ready-project`.

## Epistemic status

`PARTIALLY_SUPPORTED` — the partition is deterministic and offline-verifiable;
statistical independence between partitions is not claimed, and predictive
performance is not evaluated.

## Type

`TECHNICAL + SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commits
`4d8aa99`, `620ff85`, `45cc637`, 2026-08-31)

## Context

**APAZ** is a class of domains associated with some prokaryotic Argonaute
systems. Building an APAZ profile classifier needs a reference set partitioned so
that near-identical sequences do not span development and evaluation.

## Problem

The v1 curation keyed the indivisible partition unit on a hash of the accession
(one singleton per accession). Near-identical sequences under different
accessions could then be split across partitions — **data leakage**: information
from the evaluation set appearing in development, which can inflate apparent
performance.

## Evidence

- `src/pago_pipeline/resources/apaz_seed/curation_notes.md` (staging):
  - §2 — Ryazansky et al. 2018 **Data Set S3** contains **481 representative
    APAZ domains** (not 481 complete pAgo proteins), reduced from 632 with
    UCLUST 4.2 at 90 % identity, longest member per cluster. Groups (S3
    representatives): **Ia 83, Ib 109, IIa 98, IIb 179, III 12**; the paper's
    original group counts are 117 / 150 / 122 / 225 / 18.
  - §3 — the 481 are `source_redundancy_representative`s; the operational split
    unit is a fresh MMseqs2 90/80 connected component. "The v1 curation … is
    **superseded** and was never used for predictive evaluation."
  - §4 — hard negatives HisG (Pfam **PF01634**) and EIIB (**PF00367**), from
    Ryazansky's own search artefacts; never enter BUILD.
  - §5 — split-group method; "**not a claim that homology is absent below 90 %
    identity**".
  - §6 — partition policy `apaz_partition_v2_mmseqs90_80`: FINAL_HOLDOUT first,
    then CALIBRATION, then BUILD; deterministic subset-sum; **aborts** on an
    unreachable target.
  - §8 — one worked leak: `WP_077684787.1` (BUILD in v1) / `WP_081684175.1`
    (CALIBRATION in v1), identity 1.000, coverage 0.816 — inseparable, so one
    split group in v2.
- `src/pago_pipeline/resources/apaz_seed/split_groups/apaz/manifest.json`:
  `n_sequences = 481`, `n_split_groups = 460`, `multi_member_groups = 20`,
  `max_group_size = 3`; MMseqs2 `18.8cc5c`, `--prefilter-mode 2`, WSL2 /
  Ubuntu 24.04.4.
- Commit `45cc637` — *"All eight v1 cross-partition homology leaks are
  resolved … v1 was never used for predictive evaluation … The FINAL_HOLDOUT
  partition frozen here is definitive: no later redistribution may be motivated
  by HMM performance."*

## Previous state

APAZ partition v1: split unit = hash of the accession; homology clusters could
cross BUILD / CALIBRATION.

## Decision

- The indivisible partition unit is a **split group**: a connected component of
  an MMseqs2 all-vs-all graph with an edge iff identity ≥ 0.90 **and** coverage
  ≥ 0.80, with identity and coverage **recomputed from the aligned strings**
  (never from MMseqs `pident` / `alnlen`), plus a prior union of
  `sequence_sha256` duplicates. **Coverage** here is the fraction of the shorter
  compared region that participates in the alignment.
- **This 90/80 rule is a high-similarity / redundancy control, not a definition
  of homology.** (Homology is common evolutionary ancestry; there is no
  universal sequence-identity threshold that defines it.)
- 481 sequences → **460** split groups (20 multi-member, max size 3).
- Partition `apaz_partition_v2_mmseqs90_80`: whole split groups only;
  FINAL_HOLDOUT filled first, then CALIBRATION, then BUILD; the FINAL_HOLDOUT
  partition is **frozen** — no later redistribution may be motivated by HMM
  performance.
- HisG (PF01634) and EIIB (PF00367) hard negatives are partitioned by family
  into CALIBRATION / FINAL_HOLDOUT only; never BUILD.

## Rationale

Stated in the curation notes and commit `45cc637`: keeping a whole
high-similarity group in one partition removes the leakage that v1 allowed;
the FINAL_HOLDOUT freeze prevents the holdout being tuned against.

## Alternatives considered

The v1 accession-hash unit (rejected — leaked). The curation notes also record
the **epistemological limitation** of splitting a single negative family
(HisG / EIIB) between CALIBRATION and FINAL_HOLDOUT: it measures generalisation
to new sequences of the *same* negative family, not to unseen negative families.

## Consequences

- APAZ profile HMMs (ADR-0010) are built only from the BUILD partition.
- Keeping split groups whole **does not** guarantee complete statistical
  independence between partitions.

## Limitations

Not integrated. The exact procedure by which all **eight** v1 leaks were
enumerated is `HISTORICAL_EVIDENCE_INSUFFICIENT` — the count is stated and one
example is documented (§8), but there is no versioned v1 artefact.

## Supersedes

APAZ reference partitioning v1 (accession-hash split unit).

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`4d8aa99`, `620ff85`, `45cc637`.

## Related data / artifacts

`src/pago_pipeline/resources/apaz_seed/**` (staging): `apaz_partitions.csv`,
`split_groups/{apaz,hisg,eiib}/**`, `seeds_lock.json`,
`apaz_validation_sequences.fasta`.

## Validation

`scripts/verify_reference_data.py --scope apaz` (staging): 14 partition
invariants + 7 seed-consistency checks; byte-for-byte regeneration with
`scripts/regenerate_apaz_reference.py`; split groups byte-identical across two
independent MMseqs2 runs.

## Scientific impact

Defines the indivisible train/test unit for the APAZ work and rules out
performance-motivated holdout redistribution.

## Data impact

Adds committed partition tables and split-group resources. No candidate-protein
label assigned.

## Reproducibility impact

Positive: deterministic split-group and partition derivation, pinned tool
version and environment, offline verifier.

## Historical reconstruction note

- Directly demonstrated: the curation notes, the manifest counts, the commits.
- Inferred: none material.
- Not recoverable: how all eight v1 leaks were enumerated (one example only).
