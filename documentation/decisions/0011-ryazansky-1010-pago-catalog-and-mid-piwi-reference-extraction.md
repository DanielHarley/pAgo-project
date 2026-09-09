# ADR-0011 — Ryazansky 1010-pAgo catalog reconciliation and MID-PIWI reference extraction

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — Phase B (milestone B4.2) on `(feat)-siepe-ready-project`.

## Epistemic status

`SUPPORTED` — the catalog is a deterministic parse of a published supplementary
table, with the coordinate convention proven from the data and every count
recorded.

## Type

`TECHNICAL + SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commits
`b7490e4`, `481d0c2`, 2026-08-31)

## Context

The clade reference layer needs a curated catalog of known pAgos. **MID** and
**PIWI** are structural domains of Argonaute proteins: MID mainly anchors the
guide 5′ end; PIWI carries the catalytic core in catalytically active
Argonautes. Family, phylogenetic clade, domain architecture and catalytic
activity are **different concepts**.

## Problem

Ryazansky et al. 2018 Table S1 lists the pAgos, but the reusable inputs for
their phylogenetic tree are not published, and some accessions are now
suppressed or replaced.

## Evidence

- `src/pago_pipeline/resources/clade_seed/ryazansky_s1_catalog_notes.md`
  (staging):
  - Table S1 = **1010 full pAgos** found by PSI-BLAST over RefSeq; it
    **precedes** the redundancy reduction.
  - Ryazansky **reports** a **721**-protein nonredundant representative set
    (UCLUST 4.2, 90 % identity) used for the Fig. 1B tree. The reusable
    721-member tree-input mapping, the MID-PIWI tree MSA, and the Newick tree
    are **not** available as published supplementary artifacts (only the figure
    images).
  - Catalog (rebuild SHA-256 `dde09c03…`): **1010 rows**; **1010 / 1010**
    sequences recovered; NCBI status **731 LIVE / 272 SUPPRESSED / 7
    DEAD_REPLACED** (all 7 replacements byte-identical); effective distinct
    proteins ≈ **1004** (`n_unique_full_sequence_sha256 = 1004`; 6
    obsolete-accession double-listings).
  - **MID-PIWI extracted for 1002 / 1010** (8 rows have no MID coordinate);
    **984 unique** MID-PIWI regions; **972 canonical / 38 truncated**.
  - Coordinate convention (PAZ/MID/PIWI columns 1-based, inclusive, MID·PIWI
    contiguous) **proven from all 1010 rows** before any slice; the builder
    raises if it cannot be re-proven.
  - AfAgo (`WP_010878815.1`) and SiAgo (`WP_012735993.1`) → **QUARANTINE**,
    `curated_pago_clade = UNRESOLVED` (architecture SHORT vs Ryazansky
    truncated-long); NgAgo kept as Ryazansky `longA` (the recall panel's LONG_B
    label contradicts the source it cites — treated as a probable recall-panel
    error; the recall-panel artifact is not modified); the 2 `unkn` rows →
    `UNRESOLVED`.
  - Nothing is excluded.
- Staging commits `b7490e4` (`clade_reference_catalog.py`), `481d0c2` (frozen
  sources + `ryazansky_s1_catalog.csv`).
- Primary source: Ryazansky, Kulbachinskiy & Aravin 2018, *mBio*
  9(6):e01935-18 (PMC6299218, CC BY).

## Previous state

An earlier note claimed TtAgo and RsAgo were "absent from Table S1 because
removed by the 90 % filter". That claim is **withdrawn** in the catalog notes
(Correction 1): Table S1 is the pre-reduction 1010 set; RsAgo is present as
`WP_011910606.1` and the crystallographic HB8 TtAgo is present as `YP_145307.1`.

## Decision

- Build a deterministic 1010-row catalog from Table S1; recover every sequence
  by `accession.version`; snapshot NCBI record status; keep
  `source_ago_type` / `source_phylogenetic_clade` / `architecture_status` /
  `curated_pago_clade` / `curation_status` as **separate** fields — architecture
  never sets clade.
- Extract the MID-PIWI region as `sequence[MID_start-1 : PIWI_end]` after
  proving the coordinate convention.
- Exclude nothing; quarantine only where the architecture and the source clade
  genuinely conflict.

## Rationale

Stated in the catalog notes: Table S1 is pre-reduction, so nothing should be
dropped by a redundancy assumption; the coordinate convention must be proven,
not assumed; family / clade / architecture / catalysis are distinct, so an
architecture observation cannot override a phylogenetic clade label.

## Alternatives considered

Forcing AfAgo / SiAgo / NgAgo to a single clade (rejected — the conflict is
genuine and is recorded, not resolved). `HISTORICAL_EVIDENCE_INSUFFICIENT` for
other parse or reconciliation strategies.

## Consequences

- The 1002 MID-PIWI regions feed the split-group derivation (ADR-0013).
- The project's own reference tree, alignment and substitution model still have
  to be built (not yet decided or executed).
- Which redundant homolog stood in for RsAgo / TtAgo in Ryazansky's 721-set
  tree is **not reconstructable** from the published supplements.

## Limitations

Not integrated. The 721-set / MSA / tree are unavailable, so the original tree
cannot be reused directly. The TtAgo accession for the recall panel is a pending
curation decision (see `SCIENTIFIC_STATE.md` §4).

## Supersedes

The withdrawn "TtAgo / RsAgo removed by the 90 % filter" claim.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`b7490e4`, `481d0c2`.

## Related data / artifacts

`src/pago_pipeline/resources/clade_seed/**` (staging): `ryazansky_s1_catalog.csv`
(SHA-256 `dde09c03…`), `source_data/`, `ryazansky_s1_catalog_notes.md`.

## Validation

`scripts/verify_reference_data.py --scope clade` (staging);
`tests/test_clade_reference_catalog.py` (off-by-one guards on the MID-PIWI
slice); QC checks with 0 hits for duplicate accession, coordinate out of range,
length drift.

## Scientific impact

Provides the pAgo reference catalog and the MID-PIWI regions; enforces the
family / clade / architecture / catalysis distinction.

## Data impact

Adds a committed 1010-row catalog. No candidate-protein label assigned.

## Reproducibility impact

Positive: `--reuse-frozen` rebuilds the catalog bit-for-bit offline.

## Historical reconstruction note

- Directly demonstrated: the catalog notes, the counts, the commits, the
  primary paper.
- Inferred: none material.
- Not recoverable: the specific redundant homologs used for RsAgo / TtAgo in
  the original 721-set tree.
