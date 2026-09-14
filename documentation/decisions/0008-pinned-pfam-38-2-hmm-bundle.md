# ADR-0008 — Pin the Pfam 38.2 HMM reference bundle

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commit
`69ca40e`, 2026-08-31)

## Context

The pAgo reference layer needs a fixed set of Pfam profile HMMs for core
Argonaute and associated-system domains. A silent change in Pfam release,
profile accession version, or file contents would change downstream scans
without a trace.

## Decision

Pin ten Pfam 38.2 profiles — PIWI, PAZ, ArgoN, ArgoL1, ArgoL2, ArgoMid, SIR2,
TIR_2, TIR, and Mrr_cat — using versioned accessions and SHA-256 identities.

Reject unversioned Pfam accessions, concatenate the profiles in a canonical
order, press the bundle with PyHMMER, and provide an offline verifier for the
reference chain.

## Evidence

Historical commit `69ca40e` added the ten HMM resources, a lock file with
versioned accessions and SHA-256 values, canonical concatenation, structural
validation, `hmmpress`, and tests that reject unversioned accessions.

The historical ontology's reference-integrity section also requires pinned
source releases, file inventories, and hashes for committed reference resources.

## Rationale

Pinning release, profile version, and file hash makes the reference instrument
identifiable and guards against silent upstream changes.

`RATIONALE_NOT_RECOVERABLE` for why exactly these ten profiles were selected and
why Pfam 38.2 was chosen over another release.

## Consequences

- downstream domain scans can refer to an exact, reproducible reference bundle
- changing the Pfam release or a profile version requires an explicit reference
  update rather than silently changing results
- structural parsing and byte-level integrity do **not** establish predictive
  sensitivity, specificity, or biological validity

## Limitations

The original selection rationale for the ten profiles and the release choice is
not documented. Predictive performance has not been established by the bundle's
integrity checks.

## Validation

The historical reference verifier checks the pinned files offline, and the
historical tests cover expected profile identity and rejection of unversioned
accessions.

## Related records

- historical commits `3de5f1b` and `69ca40e`
- historical `src/pago_pipeline/resources/pfam_hmm/**`
- historical `pfam_hmm_bundle_lock.json`
- historical roadmap reference: B1

## Historical reconstruction note

- Directly demonstrated: historical commit, lock file, tests, and ontology text.
- Not recoverable: the rationale for the exact ten-profile selection and Pfam
  release choice.
