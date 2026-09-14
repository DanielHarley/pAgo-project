# ADR-0005 — Use an evidence-tiered query-recall panel with sequence-identity equivalence matching

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commits
dated 2026-08-30)

## Context

The annotation-enriched candidate set is produced by a text query. That query
needs an audit against known pAgo and PIWI-RE references without turning the
audit panel into a classifier benchmark.

## Decision

Use a curated panel of 21 references: 14 pAgo references across LONG_A,
LONG_B, and SHORT, plus 7 PIWI-RE references. Each row records the evidence tier
supporting its reference label.

Match retrieved records in this order:

1. exact `accession.version`
2. the same base accession
3. normalized-sequence SHA-256 equivalence

Keep both the exact-accession and retrieval-equivalent recall readings. Treat an
empty stratum as not evaluable rather than as zero recall.

The panel is a query-coverage instrument, not a classifier gold standard and not
a final holdout. In particular, the six PIWI-RE references selected using
PIWI-RE profile models cannot be reused uncritically to validate a profile-based
PIWI-RE detector.

## Evidence

- historical commits `599a111`, `d21f25a`, `0fd283a`, `21793bc`, `98e66c6`,
  `34c0ca6`, `0013c6d`
- historical curation notes for the 21-reference panel and the RsAgo accession
  case
- Phase A audit decisions covering the panel, sequence-hash matcher, and the two
  recall readings

## Rationale

Evidence tiers make label strength explicit. Sequence-hash equivalence allows a
real sequence recovery under another accession to be recorded without rewriting
the reference fixture or hiding an exact-accession miss.

Equal normalized-sequence hashes establish sequence equivalence for this recall
report. They do not establish identical locus identity, record provenance, or
biological context.

## Consequences

- all recall strata can be reported consistently
- the LONG_B stratum remains small
- the PIWI-RE portion carries an explicit circularity boundary
- query recall must not be promoted into detector-validation performance

## Limitations

The panel is small and annotation-enriched. It does not estimate recall over all
divergent, poorly annotated, or unknown pAgos.

## Supersedes

Exact-accession-only recall matching and the early fixture convention that
represented PIWI-RE as a `pago_clade` value.

## Validation

Historical tests cover match priority, empty-stratum handling, sequence-equivalent
RsAgo recovery, and the pinned matching-strategy identity.

## Related records

- historical query-recall fixture and curation notes
- `documentation/history/2026-08-30-phase-a-audit.md`
- ADR-0012 for family-versus-clade semantics
- ADR-0014 for the independent PIWI-RE validation reference requirement

## Historical reconstruction note

- Directly demonstrated: historical commits, curation notes, and the Phase A
  audit.
- Not recoverable: the full deliberation over the final panel size and
  composition.
