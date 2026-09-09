# ADR-0014 — Independent PIWI-RE reference arm (proposed, not executed)

## Lifecycle status

`PROPOSED`

## Integration status

`NOT_INTEGRATED`

## Epistemic status

`UNRESOLVED` — no independent PIWI-RE reference set exists yet.

## Type

`SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09`

## Context

A future PIWI-RE detector (a profile HMM, or a placement arm) would need an
evaluation set of PIWI-RE references that is independent of the models used to
build it. **Circularity** means using the same kind of evidence or method to
choose the examples and then using those examples to "validate" the method.

## Problem

The current query-recall panel's PIWI-RE references cannot serve as an
independent validation set.

## Evidence

- `tests/fixtures/query_recall_reference_set_curation_notes.md` (staging),
  "PIWI-RE set (7 entries)" and "Not a validation holdout": 1 `EXPERIMENTAL`
  (PsPIWI-RE `WP_014597637.1`, Huang et al. 2022) and **6
  `CURATED_COMPUTATIONAL`** — "NCBI RefSeq/CDD representatives that appear in
  two or three of the PIWI-RE-specific CDD/Pfam models"; "the
  CURATED_COMPUTATIONAL members were selected *using* PIWI-RE profile models, so
  scoring a profile-based detector against them would be circular"; "This set …
  must **not** become the `FINAL_HOLDOUT` for a future PIWI-RE HMM".
- ADR-0005 (the recall panel) and ADR-0012 (PIWI-RE is a family).
- Primary literature: Burroughs, Iyer & Aravind 2013 (*Biology Direct* 8:13,
  PMC3702460) — the founding PIWI-RE description; Huang et al. 2022 — PsPIWI-RE
  experimental characterisation.

## Previous state

No independent PIWI-RE reference arm; the recall panel's 7 PIWI-RE references
exist only for query-recall auditing.

## Decision

- An independent PIWI-RE reference arm, separate from Ryazansky, **must** be
  curated from primary PIWI-RE literature — starting with Burroughs et al. 2013
  and subsequent experimental sources — with explicit circularity auditing of
  every included reference.
- The current recall-panel PIWI-RE references are **not** automatically eligible
  as an independent validation set, because 6 of the 7 are
  `CURATED_COMPUTATIONAL`, selected using PIWI-RE profile models.
- This work is **not started** (roadmap reference B4.4).

## Rationale

Stated in the curation notes: using profile-model-selected examples to score a
profile-based detector is circular; an independent arm must trace to primary
experimental / literature evidence.

## Alternatives considered

Reusing the recall panel's PIWI-RE rows as a holdout — explicitly rejected in
the curation notes. `HISTORICAL_EVIDENCE_INSUFFICIENT` for the specific
composition of a future independent arm.

## Consequences

- Any future PIWI-RE detector validation depends on this arm being built first.
- Until then, no PIWI-RE detection performance claim is admissible.

## Limitations

Not started; no evidence yet supports any concrete composition.

## Supersedes

Nothing.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`0fd283a`, `21793bc` (the recall-panel PIWI-RE additions — **not** this arm).

## Related data / artifacts

`tests/fixtures/query_recall_reference_set_curation_notes.md` (staging).

## Validation

None — proposed only.

## Scientific impact

Prevents a circular PIWI-RE validation.

## Data impact

None yet.

## Reproducibility impact

None yet.

## Historical reconstruction note

- Directly demonstrated: the curation notes' explicit circularity statement and
  "not a validation holdout" section.
- Inferred: that this constitutes a standing `PROPOSED` decision rather than a
  passing remark (supported by the emphatic wording).
- Not recoverable: the future arm's composition — deliberately left open.
