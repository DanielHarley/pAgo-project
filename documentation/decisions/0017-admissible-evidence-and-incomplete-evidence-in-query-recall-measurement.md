# ADR-0017 — Define admissible evidence and incomplete-evidence reporting for query-recall measurement

## Date

2026-10-22

## Context

ADR-0005 decides the query-recall matching hierarchy (exact `accession.version`,
then the same base accession, then normalized-sequence SHA-256 equivalence) and
that both the strict-accession and retrieval-equivalent readings are kept. It
presupposes well-formed identifiers and complete sequence data.

Real flattened NCBI metadata does not guarantee either. Cells can be missing,
version suffixes can be lost, identifier columns can carry malformed tokens,
and a retrieval can arrive without usable protein sequences. Measuring recall
requires deciding what counts as evidence for each tier, and what to report
when the evidence needed for a tier is not there.

Two concrete failure modes made the decision necessary rather than theoretical.
A retrieved accession that had lost its version suffix satisfied the
same-base tier for a reference whose sequence was entirely different, and was
counted as strict-accession recovery. A retrieval without a usable sequence
column produced a retrieval-equivalent recall of `0.0` marked `EVALUABLE`, so
an unevaluated tier read as a completely evaluated zero.

## Decision

**Admissibility.** A retrieved accession participates in the accession tiers
only in explicit `accession.version` form, matching
`^[A-Za-z0-9_]+\.\d+$`. A protein UID is admissible only as `^[0-9]+$`. An
inadmissible identifier is excluded from the tier that depends on it, counted,
and left eligible for sequence matching through its own sequence.

A reference accession must be admissible by the same pattern; an inadmissible
one is a fatal error, because the tier-1 and tier-2 distinction has no meaning
without an explicit version.

**Attribution.** A retrieved sequence observation is usable only when it can be
attributed to an NCBI record by an admissible accession or an admissible UID.
An observation with neither is excluded from the sequence index and from the
count of usable sequence evidence. An inadmissible token is never emitted as
identification in the output; it survives only as a count.

**Cell contract.** A retrieved cell must be scalar; a container is rejected
before missingness is considered. Admitted missing scalars are treated as
absent. A present value must be a string; any other type is rejected rather
than coerced.

**Incomplete evidence.** A missing required column is fatal. A present column
with missing values is not: the number of retrieved rows and the number with
usable sequence evidence are recorded, and per stratum the retrieval-equivalent
reading is marked a lower bound when the evidence is incomplete *and* the
stratum is not already saturated. A stratum whose references are all recovered
has a determined value regardless of evidence missing elsewhere.

**Panel policy.** Reference base accessions must be unique. Duplicate reference
sequence hashes are permitted and reported.

**Integrity.** Conflicting non-missing sequences under one retrieved
`accession.version` are a fatal integrity error.

## Rationale

Rejecting an entire retrieval because one upstream cell is malformed would make
the audit hostage to metadata noise; accepting a truncated identifier
manufactures recovery that was never demonstrated. Excluding the token from the
tier that depends on it, while keeping the independent sequence evidence that
row still carries, preserves both robustness and correctness.

Requiring an explicit version on both sides is what gives the two accession
tiers distinct meanings. Without it the exact-version tier can never apply, and
the same-base tier answers in its place while the output still claims
same-accession recovery.

Reporting incompleteness, rather than folding it into a complete-looking value,
keeps the error in the direction that can be read correctly: the figure can
become a declared lower bound, never a silently inflated value. The saturation
rule keeps that honest in the other direction too — it does not label a value a
bound when it is in fact determined.

Base-accession uniqueness is a property of this instrument rather than an
arithmetic necessity. Recall counts panel rows, so two rows sharing a base
accession each legitimately ask whether that record was retrieved; a single
retrieved record answering both would be counted twice in the strict reading,
which defeats a coverage instrument. Duplicate sequence hashes are different:
distinct accessions can legitimately encode identical proteins, so forbidding
them would impose a biological claim the evidence does not support.

One `accession.version` denotes one immutable record, so two different
sequences under it mean the metadata is corrupt rather than merely ambiguous.

## Consequences

- A negative at the strict reading is not proof that the record was absent from
  the retrieval: an identifier may have been excluded. The count of
  inadmissible retrieved accessions is what qualifies that reading, and it is
  recorded with every measurement.
- A retrieval-equivalent value of 1.0 does not establish sequence identity for
  every reference: it is the union of the three admitted tiers, and a reference
  may be counted through accession identity alone. The per-reference match
  method is the only field that states how a reference was recovered.
- Duplicate reference sequence hashes remain loadable. When they occur, one
  retrieved sequence recovers every row in the group, and the group is reported
  so the figure is not read as several distinct proteins.
- Base-accession uniqueness stays fatal until a decision redefines the same-base
  tier to bind at most one reference per base accession.
- Changing any of these policies changes what the metric counts, so it changes
  the matching-strategy identity and invalidates stored measurements even when
  the data inputs are unchanged.

## Related records

- ADR-0005 for the panel, the matching hierarchy, and the two readings
- ADR-0001 for the protein UID as the canonical retrieval identifier
- ADR-0002 for the immutable snapshot and provenance architecture
- ADR-0012 for the family/clade separation used by the recall strata
- ADR-0014 for the PIWI-RE circularity boundary the panel carries
