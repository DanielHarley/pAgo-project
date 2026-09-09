# ADR-0004 — Annotation-enriched candidate set; the prefilter is technical, not biological

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — Phase A work on `(feat)-siepe-ready-project`.

## Epistemic status

`SUPPORTED` — a sound design with an audited execution. It is explicitly **not**
a set of confirmed pAgos.

## Type

`TECHNICAL + SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commits dated
2026-08-30)

## Context

A supervised pAgo classifier would need a candidate set. A broad text query over
NCBI Protein returns a biologically heterogeneous universe of PIWI/Argonaute-
associated proteins.

## Problem

Acquire a candidate set without implying that a text hit is a confirmed pAgo,
and without discarding real but atypical proteins at the acquisition stage.

## Evidence

- Staging commits: `b54c698` (notebook 10 + README section), `f7e91ec`
  (ESearch preflight), `175a12a` (technical prefilter), `8a200ff` / `c66e467`
  (derived FASTA + SWeeP input).
- Untracked local manifests (now gitignored, preserved on disk):
  `search_query = "(PIWI[All Fields] OR Argonaute[All Fields]) AND
  (Bacteria[Organism] OR Archaea[Organism])"`, `result_count = 52473`,
  `consolidated_record_count = 52473`, `xml_file_sha256 = 77739043…6169`,
  `protein_uids_sha256 = 949bcaab…b253`.
- `documentation/history/2026-08-30-phase-a-audit.md` — decisions D1
  (dataset name), D2 (not a universe), D3 (no annotation-text filter), D4 (no
  length filter); §13.2 (execution numbers); §17.1–§17.2 (what can / cannot be
  concluded).

## Previous state

The earlier narrow query `PIWI[All Fields] AND Bacteria[Organism]` (still the
integrated query in `master`).

## Decision

- Query `(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism]
  OR Archaea[Organism])`; name the result `annotation_enriched_candidate_set`.
- The set is **candidate proteins enriched by textual annotation associated
  with PIWI/Argonaute**, never "confirmed pAgos".
- An ESearch preflight records `Count`, the translated query and a sample QC,
  and blocks a full retrieval above `max_uid_count` until reviewed.
- The **technical prefilter** excludes only technically unusable records
  (missing / invalid sequence, missing `protein_uid`, technical duplicates). It
  **never** filters by annotation text and **never** by sequence length; length
  outside a warning band sets `length_warning = True` and the record is kept.

## Rationale

Stated in the phase-A audit (D2–D4): the name carries the annotation bias;
separating pAgo from non-pAgo is a domain / HMM task for a later phase, decided
on structural evidence, not a text or length heuristic at acquisition.

## Alternatives considered

Per D3/D4 of the audit: a term blacklist or an annotation classifier at the
prefilter stage (rejected — noisy and circular); a length cut outside a band
(rejected — fragments and fusions can be real). `HISTORICAL_EVIDENCE_INSUFFICIENT`
for other query formulations.

## Consequences

- One execution (2026-08-30/31) returned **52,473** records; the technical
  prefilter retained all 52,473 (0 exclusions) with **8,299** carrying
  `length_warning`.
- The retained proteome contains textual false positives (transposases,
  methyltransferases, restriction endonucleases) to be removed by a later
  domain scan.
- The 52,473-record execution artifacts are regenerable but not versioned.

## Limitations

The candidate set is a coverage of a text query, not a pAgo universe. A
sequence-based discovery route (HMM / PSI-BLAST over RefSeq) is out of scope and
documented as future work.

## Supersedes

(pending reintegration) the narrow `PIWI[All Fields] AND Bacteria[Organism]`
query and its 41,345-record filtered set (see ADR-0007 and `pago_qc.md`).

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`f7e91ec`, `8629d20`, `b54c698`, `175a12a`, `8a200ff`, `c66e467`.

## Related data / artifacts

`data/01-raw/*__annotation_enriched_candidate_set/**` (local, gitignored);
`notebooks/10_dataset_audit.ipynb`; `tmp/nb10_out_of_scope_preserved/` (the
executed working copy, preserved for provenance).

## Validation

51 Phase A unit tests (part of the 203 in the staging suite) cover the
prefilter's "technical only" contract, including
`test_retains_records_regardless_of_annotation_text_or_length`. The 52,473
figure is `[EXEC-LOCAL]` — from local manifests, not a versioned check.

## Scientific impact

Establishes the epistemic boundary: text recovery ≠ existence ≠ identity. Feeds
every later phase.

## Data impact

Defines the candidate set and the retained proteome. No biological label is
assigned.

## Reproducibility impact

Conditional: a clean clone does not reproduce the 52,473-record chain without
re-running notebook 10 against NCBI.

## Historical reconstruction note

- Directly demonstrated: commits, the query string and counts in local
  manifests, the phase-A audit.
- Inferred: none material.
- Not recoverable: the exact deliberation over alternative query strings.
