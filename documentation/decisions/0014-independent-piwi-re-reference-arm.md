# ADR-0014 — Require an independent PIWI-RE reference arm for detector validation

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09`

## Context

A future PIWI-RE detector requires evaluation references that were not selected
using the same profile-model evidence that the detector itself would later be
asked to validate.

The existing query-recall panel contains seven PIWI-RE references, but six were
selected computationally using PIWI-RE-specific profile-model evidence.

## Decision

Before making a detector-performance claim for PIWI-RE, curate an independent
PIWI-RE reference arm from primary PIWI-RE literature and experimental evidence,
with explicit circularity auditing for every included reference.

Do not automatically reuse the existing query-recall PIWI-RE rows as an
independent validation set for a profile-based detector.

## Evidence

The historical query-recall curation notes state that:

- one PIWI-RE panel reference has experimental support
- six are `CURATED_COMPUTATIONAL` references selected using PIWI-RE-specific
  CDD/Pfam profile evidence
- scoring a profile-based detector against those computationally selected
  references would be circular
- the recall panel must not become the final holdout for a future PIWI-RE HMM

Burroughs, Iyer & Aravind 2013 provides the founding PIWI-RE family context, and
later experimental literature such as the PsPIWI-RE work provides a starting
point for independent curation.

## Rationale

Using profile-selected examples to validate a profile-based detector risks
measuring agreement with the selection method rather than independent detector
performance.

## Consequences

- no independent PIWI-RE detector-performance claim should be made until this
  reference arm exists
- the query-recall panel remains useful for query-coverage auditing but does not
  become a detector-validation gold standard
- future reference curation must preserve the source and evidence supporting
  each PIWI-RE label

## Limitations

No concrete composition for the independent arm has yet been established. The
historical record supports the need for independence, not a final reference
list.

## Related records

- historical PIWI-RE query-recall curation notes
- ADR-0005 for the query-recall panel
- ADR-0012 for PIWI-RE family semantics
- Burroughs, Iyer & Aravind 2013 and subsequent experimental PIWI-RE literature
- historical roadmap reference: B4.4

## Historical reconstruction note

- Directly demonstrated: the historical curation notes explicitly reject the
  recall-panel PIWI-RE set as a profile-detector holdout.
- Not recoverable because it was never defined: the final composition of the
  independent PIWI-RE reference arm.
