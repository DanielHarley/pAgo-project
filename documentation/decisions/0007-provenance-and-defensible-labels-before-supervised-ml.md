# ADR-0007 — Provenance, defensible labels, redundancy control and a protected holdout come before supervised ML

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — the reference / validation layer that implements this
(Phase B) is on `(feat)-siepe-ready-project`.

## Epistemic status

`PARTIALLY_SUPPORTED` — the in-repository framing is verifiable; parts of the
broader machine-learning-methodology justification are literature-dependent and
not all primary sources have been re-verified here.

## Type

`SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09`

## Context

The original aim is a computational approach to identify prokaryotic Argonautes,
with interest in mesophilic organisms, eventually with supervised machine
learning (**ML** — machine learning). A **holdout** is a reserved evaluation set
that must not be used to develop, calibrate or choose a model.

## Problem

Training a supervised classifier directly on NCBI records risks producing
plausible but conceptually wrong labels that are then treated downstream as
experimental truth and as ML ground truth.

## Evidence

- `documentation/history/2026-08-30-phase-a-audit.md` §17.2 ("we cannot
  claim…") and §17.3 (future phases), and D2 ("epistemological integrity — do
  not confuse *recoverable by text* with *existing*").
- `src/pago_pipeline/resources/apaz_seed/curation_notes.md` and
  `.../clade_seed/ryazansky_s1_catalog_notes.md` (staging) — every curation
  decision is recorded with its evidence tier and its limitations.
- `documentation/history/2026-08-31-pago-annotation-ontology-staging.md` — the
  BUILD / CALIBRATION / FINAL_HOLDOUT protocol for custom HMMs and the
  TREE_BUILD / PLACEMENT_CALIBRATION / PLACEMENT_HOLDOUT protocol for placement;
  the reference-label evidence tiers.
- Primary literature already used in the project (cited in the ontology and
  curation notes): Ryazansky et al. 2018 (mBio 9:e01935-18) — a broad
  PIWI-related search recovers a heterogeneous universe of ~1,010 pAgos across
  diverse architectures.

## Previous state

The pipeline reached exploratory clustering (ADR-0006) directly from a
text-annotation-derived candidate set with metadata-derived QC labels
(ADR-0004, `pago_qc.md`).

## Decision

Before any final supervised model, prioritise: (1) traceable provenance;
(2) biologically defensible labels with explicit evidence tiers; (3) redundancy
control; (4) high-similarity control between development and evaluation;
(5) a cluster-separated development / evaluation split; (6) a protected holdout
consulted once. No final model has been trained or validated.

## Rationale

Stated in the phase-A audit and the curation notes: the overriding risk is
conceptually-wrong-but-plausible labels becoming downstream ground truth.
`RATIONALE_NOT_RECOVERABLE` for the full literature argument as a versioned
statement — the broader data-leakage / negative-set / annotation-bias literature
appears only in local research notes (`tmp/…`), which are a map of leads, not an
authority (see PROJECT_LEDGER, "Local bibliographic material").

## Alternatives considered

Proceeding to a supervised classifier on public annotations as ground truth —
rejected, per the audit's "we cannot claim" list. `HISTORICAL_EVIDENCE_INSUFFICIENT`
for other sequencings of the work.

## Consequences

- Phase B builds evidence instruments (Pfam bundle, APAZ, clade catalog,
  MID-PIWI split groups) and freezes calibration policies and one-shot
  validation reports *before* any candidate scan.
- `FINAL_HOLDOUT` / `PLACEMENT_HOLDOUT` must not be consulted early.

## Limitations

The decision is reconstructed from framing text and curation practice, not from
a single design document in the repository. B5, B6 and final validation are not
integrated.

## Supersedes

The implicit path "candidate set → exploratory clustering → supervised model".

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found` — Phase B staging commits `3de5f1b` … `49651f1`.

## Related commits

`3de5f1b`, `69ca40e`, `4d8aa99`, `620ff85`, `45cc637`, `9e76e78`, `c11732a`,
`b7490e4`, `481d0c2`, `18c35e4`, `49651f1`.

## Related data / artifacts

`src/pago_pipeline/resources/{pfam_hmm,apaz_seed,clade_seed}/**` (staging).

## Validation

`scripts/verify_reference_data.py --scope {pfam,apaz,clade}` on the staging
branch. No model validation exists.

## Scientific impact

Reframes the project as a reproducible prospecting + evidence-tiered curation +
exploratory-analysis pipeline, not a classifier.

## Data impact

Motivates the reference layer; no candidate-protein label is assigned in it.

## Reproducibility impact

Central: every reference resource is SHA-256-pinned and offline-verifiable.

## Historical reconstruction note

- Directly demonstrated: the audit's conclusions, the curation notes, the
  ontology protocol sections, the Phase B commit sequence.
- Inferred: that these together constitute a single deliberate methodological
  priority.
- Not recoverable / not confirmed: the full ML-methodology literature
  justification as a versioned artefact — flagged `PARTIALLY_SUPPORTED`.
