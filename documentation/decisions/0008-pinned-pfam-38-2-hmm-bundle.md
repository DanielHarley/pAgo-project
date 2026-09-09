# ADR-0008 — Pinned Pfam 38.2 HMM bundle (ten profiles, versioned, SHA-256-locked)

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — Phase B (milestone B1) on `(feat)-siepe-ready-project`.

## Epistemic status

`PARTIALLY_SUPPORTED` — the bundle is pinned and offline-verified; the
predictive performance of scanning with it has not been evaluated.

## Type

`TECHNICAL + SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commit
`69ca40e`, 2026-08-31)

## Context

**Pfam** is a database of protein families and domains represented by profile
models. An **HMM** (Hidden Markov Model) is a probabilistic model used here to
recognise conserved sequence patterns of a domain. The pAgo work needs a fixed
set of core-Argonaute and system-associated domain profiles.

## Problem

A domain scan must use a fixed, reproducible set of profiles; a silent Pfam
release or profile-version change would alter results without a trace.

## Evidence

- Staging commit `69ca40e` — *"(feat) B1 pinned Pfam HMM reference bundle"*:
  ten version-pinned HMMs (PIWI, PAZ, ArgoN, ArgoL1, ArgoL2, ArgoMid, SIR2,
  TIR_2, TIR, Mrr_cat), a SHA-256 lock, structural validation against PyHMMER
  (single model, amino alphabet, gathering cutoffs, exact accession/name),
  canonical concatenation, and `hmmpress`.
- The commit adds a rule that an **unversioned** Pfam accession is rejected,
  with a positive and a negative test.
- `documentation/history/2026-08-31-pago-annotation-ontology-staging.md` —
  "Reference integrity": committed Phase B resources pin source releases, file
  inventories, SHA-256 values.
- Staging: `src/pago_pipeline/resources/pfam_hmm/` +
  `pfam_hmm_bundle_lock.json`; `scripts/materialize_pfam_hmm_bundle.py`.

## Previous state

No committed domain-profile set.

## Decision

Commit ten Pfam 38.2 profiles as `<Name>__<ACC.version>.hmm` with a SHA-256
lock; reject an unversioned accession; concatenate in canonical order and press
via `pyhmmer.hmmer.hmmpress`; verify the whole chain offline with
`scripts/verify_reference_data.py --scope pfam`.

## Rationale

Stated in the commit and the ontology's reference-integrity section: pinning
release + version + hash makes the instrument reproducible and prevents a silent
upstream change.

## Alternatives considered

`RATIONALE_NOT_RECOVERABLE` for **which** ten profiles and **why Pfam 38.2** —
the commit lists the names ("core Argonaute and system-associated domains") but
no design document weighs the selection or the release.

## Consequences

- Any domain scan (Phase C) uses exactly this bundle.
- Changing a profile version or the Pfam release requires a new lock and an
  explicit decision.

## Limitations

Not integrated. The profile selection rationale is not documented. Scanning
performance (sensitivity / specificity) is not evaluated.

## Supersedes

Nothing.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`3de5f1b` (pin `pyhmmer==0.12.3`), `69ca40e` (B1 bundle).

## Related data / artifacts

`src/pago_pipeline/resources/pfam_hmm/**`, `pfam_hmm_bundle_lock.json`
(staging).

## Validation

`scripts/verify_reference_data.py --scope pfam` on the staging branch;
`tests/test_pfam_hmm_bundle*`.

## Scientific impact

Fixes the domain vocabulary used to describe pAgo architecture.

## Data impact

Adds committed reference resources; no candidate-protein label assigned.

## Reproducibility impact

Positive: release + version + hash pinned; offline verifier.

## Historical reconstruction note

- Directly demonstrated: the commit, the lock file, the ontology text.
- Inferred: none material.
- Not recoverable: the rationale for the specific ten profiles and the release
  choice.
