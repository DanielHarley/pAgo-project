# ADR-0001 — Use the protein UID as the canonical NCBI retrieval identifier

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (implementation dated
2026-03-31, PR #5)

## Context

Early retrieval used the NCBI `accession.version` string as the record
identifier. NCBI protein records also carry a numeric UID (a GI number in the
`db=protein` context).

The pipeline needs one canonical identifier to key ESearch → EFetch → XML
validation and snapshot manifests.

## Decision

Use the protein UID as the canonical retrieval and XML-validation identifier.

## Evidence

- PR #5 merge `71ff9b2` — *"(feat) switch NCBI protein snapshot retrieval from
  accession version IDs to UIDs"* — and implementation commit `a10d89b`
- current `ncbi_snapshot.py` extracts the UID from `GBSeqid` fields and validates
  XML record UIDs against the requested UIDs
- the historical branch named for reverting to `accession.version` predates PR
  #5 and contains no implemented reversal

## Rationale

`RATIONALE_NOT_RECOVERABLE`. The surviving commit and PR record the switch but
do not preserve a design comparison between the identifiers.

NCBI documentation indicates that sequence revision changes can assign a new GI
and increment the `accession.version` suffix, so no claim is made that either
identifier is invariant across sequence revisions.

## Consequences

- raw snapshots and manifests key on `protein_uid`
- the Phase A query-recall matcher may still use `accession.version`, base
  accession, and normalized-sequence SHA-256 as benchmark-matching strategies;
  that does not change the canonical retrieval identifier
- the recorded UID list, frozen XML, and manifest hashes together identify the
  records actually analysed

## Supersedes

`accession.version` as the canonical retrieval identifier.

## Validation

The repository test suite covers UID extraction and XML-vs-UID validation in the
NCBI snapshot workflow.

## Related records

- PR #5, merge `71ff9b2`
- commit `a10d89b`
- `data/01-raw/**/manifest.json`
- ADR-0005 for query-recall matching semantics

## Historical reconstruction note

- Directly demonstrated: PR #5 merged, current code uses UIDs, and the later
  reversal-named branch contains no reversal.
- Not recoverable: the original rationale for choosing UID over
  `accession.version`.
