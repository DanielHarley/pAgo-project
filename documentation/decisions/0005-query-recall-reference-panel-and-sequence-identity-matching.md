# ADR-0005 — Query-recall reference panel and sequence-identity equivalence matching

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — Phase A work on `(feat)-siepe-ready-project`.

## Epistemic status

`SUPPORTED` for auditing whether the text query recovers known references. It is
**not** a gold standard and **not** a validation holdout.

## Type

`SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commits dated
2026-08-30)

## Context

The candidate set (ADR-0004) is a coverage of a text query. That coverage needs
to be measured against known pAgo / PIWI-RE proteins.

## Problem

Measure query recall honestly, stratified, with the strength of each reference
label made explicit, and without a benchmark artefact making a real recovery
look like a miss.

## Evidence

- Staging commits `599a111` (initial recall stage), `d21f25a` (`NOT_EVALUABLE`
  for empty strata), `0fd283a` (PsPIWI-RE; stratify PIWI_RE by `ago_family`),
  `21793bc` (+6 PIWI-RE `CURATED_COMPUTATIONAL`), `98e66c6` (test drops
  `PIWI_RE` from allowed `clade` values), `34c0ca6` (document the RsAgo miss),
  `0013c6d` (`SEQUENCE_SHA256` matcher; two readings).
- `tests/fixtures/query_recall_reference_set_curation_notes.md` (staging) — the
  21-row panel, the offline hash derivation, the RsAgo case, the historical GI
  resolution, and the explicit "Not a validation holdout" section.
- `documentation/history/2026-08-30-phase-a-audit.md` D6 (panel), D7
  (`SEQUENCE_SHA256`), D8 (two readings), P1–P4.

## Previous state

Matching only on `accession.version` (and its version-stripped form).

## Decision

- A curated panel of **21** references: 14 pAgo across `LONG_A` / `LONG_B` /
  `SHORT`, 7 PIWI-RE. Each row carries a `reference_label_evidence` tier
  (`EXPERIMENTAL` / `LITERATURE_PHYLOGENETIC` / `CURATED_COMPUTATIONAL` /
  `DATABASE_ANNOTATION`).
- Matching hierarchy: exact `accession.version` → same base accession →
  normalized-sequence SHA-256 (`SEQUENCE_SHA256`); two recall readings are kept
  (`exact_accession_recall`, `retrieval_equivalent_recall`); the matching
  strategy is pinned by `matching_strategy_sha256`.
- An empty stratum is `NOT_EVALUABLE`, never `0.0`.
- The panel **must not** be used as a validation gold standard or as a
  `FINAL_HOLDOUT`: 6 of the 7 PIWI-RE references are `CURATED_COMPUTATIONAL`,
  selected *using* PIWI-RE profile models, so scoring a profile-based detector
  against them would be circular.

## Rationale

Stated in the curation notes and audit: the tiers make label strength explicit;
`SEQUENCE_SHA256` records a **sequence-equivalent** recovery — the reference's
exact normalized protein sequence appears in the retrieved set under some
accession (RsAgo `ABP72561.1` under the byte-identical IPG alias `A4WYU7.1`) —
offline, deterministically, and without falsifying the fixture; keeping both
readings avoids hiding either the exact-accession miss or the
sequence-equivalent recovery.

Equal normalized-sequence SHA-256 demonstrates **sequence equivalence for the
recall report** (equal hashes ⇒ byte-identical normalized sequence with
extremely high probability). It does **not**, on its own, demonstrate locus
identity, record provenance, or biological context. The reported recall is that
the 21 selected references were recovered by sequence equivalence — not a
statement about all known pAgos.

## Alternatives considered

Per D7 of the audit: swap the fixture accession to `A4WYU7.1` (rejected by the
project owner — `ABP72561.1` is a correct RsAgo accession); resolve IPG online
at runtime (rejected — introduces a network dependency). Per P1–P3: several
historical Burroughs 2013 GI identifiers were resolved and 3 excluded (dead /
suppressed record, or the associated restriction endonuclease).

## Consequences

- All five recall strata become `EVALUABLE`.
- `LONG_B` has only 2 references; 6 of 7 PIWI-RE references carry the
  circularity caveat.
- The panel is a query-coverage instrument, not a classifier benchmark.

## Limitations

Small (21) and annotation-enriched — the best case for a text search; it does
not estimate recall over divergent or mis-annotated pAgos.

## Supersedes

- Exact-accession-only recall matching.
- `PIWI_RE` treated as a value of `pago_clade` in the early fixtures/tests
  (removed in commit `98e66c6`; see ADR-0012).

## Superseded by

None. The **independent** PIWI-RE reference arm (ADR-0014) is `PROPOSED` and
distinct from this panel.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`599a111`, `d21f25a`, `0fd283a`, `21793bc`, `98e66c6`, `34c0ca6`, `0013c6d`.

## Related data / artifacts

`tests/fixtures/query_recall_reference_set.csv` and `_curation_notes.md`
(staging); `data/02-intermediate/query_reference_recall__…/` (local).

## Validation

`test_query_reference_recall*` (staging suite): match-priority tests,
`test_empty_stratum_is_not_evaluable_not_zero`,
`test_rsago_is_recovered_by_sequence_identity_under_alias_accession`,
`matching_strategy_sha256` pinned. Recall *figures* are `[EXEC-LOCAL]`.

## Scientific impact

Defines what "recall" means here (panel coverage of a text query) and rules the
panel out as a classifier benchmark.

## Data impact

Adds the 21-row reference fixture. No candidate-protein label is assigned.

## Reproducibility impact

The matcher is offline and deterministic; the panel hashes were derived from a
specific local metadata CSV.

## Historical reconstruction note

- Directly demonstrated: commits, the curation notes, the phase-A audit.
- Inferred: none material.
- Not recoverable: the full deliberation over panel size and composition.
