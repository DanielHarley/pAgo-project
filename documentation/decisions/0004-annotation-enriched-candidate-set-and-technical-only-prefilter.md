# ADR-0004 — Treat the annotation-enriched query result as a candidate set and keep the prefilter technical only

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commits
dated 2026-08-30)

## Context

A broad NCBI Protein text query can recover a heterogeneous set of records whose
annotations mention PIWI or Argonaute. That query result is useful for candidate
discovery but does not establish that every retrieved protein is a pAgo.

## Decision

Use the query

`(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])`

and name the result `annotation_enriched_candidate_set`.

Describe it as candidate proteins enriched by textual annotation associated with
PIWI/Argonaute, not as a universe of confirmed pAgos.

Apply only a technical prefilter at this stage. Exclude records for technical
reasons such as missing/invalid sequence, missing `protein_uid`, or technical
duplication. Do not exclude candidates solely because of annotation text or
sequence length; length outside the declared range produces a warning rather
than a biological negative decision.

## Evidence

- historical commits `f7e91ec`, `8629d20`, `b54c698`, `175a12a`, `8a200ff`,
  `c66e467`
- the Phase A audit records the dataset naming and the explicit rules against
  annotation-text and length filtering
- local execution manifests recorded the query string, 52,473 retrieved
  records, and matching UID/XML hashes for the audited 2026-08-30/31 run

## Rationale

Text annotation is a discovery aid, not sufficient evidence of pAgo identity.
Biological discrimination should happen later using appropriate sequence/domain
or other defensible evidence rather than a text or length heuristic at
acquisition time.

## Consequences

- the audited execution retained all 52,473 retrieved records through the
  technical prefilter
- 8,299 records carried `length_warning` rather than being removed
- downstream work must not interpret the candidate-set size as the number of
  pAgos
- a sequence-driven discovery route remains necessary if the project later aims
  to recover poorly annotated or divergent pAgos outside text-query coverage

## Supersedes

The earlier narrow query design as the intended candidate-set strategy, once the
broader design is adopted in the active pipeline.

## Validation

Historical Phase A tests cover the technical-only contract, including retention
regardless of annotation text or sequence length. The 52,473-record counts are
execution evidence from local manifests rather than versioned result artifacts.

## Related records

- `documentation/history/2026-08-30-phase-a-audit.md`
- historical commits `f7e91ec`, `8629d20`, `b54c698`, `175a12a`, `8a200ff`,
  `c66e467`
- `notebooks/10_dataset_audit.ipynb` in the historical work

## Historical reconstruction note

- Directly demonstrated: historical commits, query/count manifests, and the
  Phase A audit.
- Not recoverable: the full deliberation over alternative query formulations.
