# ADR-0003 — NCBI retrieval performance and correctness rework (F1–F12 / P0–P5)

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — exists only on `(feat)-siepe-ready-project`.

## Epistemic status

`SUPPORTED` — the two staging documents record a measured baseline, an
`rettype=gp` invariance test, and byte-identical consolidation at full scale.

## Type

`TECHNICAL`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commits dated
2026-08-15)

## Context

Retrieving the full candidate set from NCBI E-utilities was slow and fragile:
UID pagination discarded the History handle, consolidation held the whole XML in
memory, a rate-limit response could kill a run, and partial responses had no
resume path.

## Problem

Retrieve tens of thousands of records reliably and reproducibly, and prove the
payload contract is stable.

## Evidence

- `documentation/ncbi_retrieval_performance_plan.md` and
  `documentation/ncbi_retrieval_performance_implementation.md` — **staging
  copies on `(feat)-siepe-ready-project`, not in `master`**. They contain the
  measured baseline, cost model, findings F1–F12, phased plan P0–P5, an
  acceptance-criteria section, and a "two defects found while running the real
  workflow" section.
- Staging commits `f5671b8` … `0b4ed5a` (2026-08-15): streaming XML parser,
  request telemetry + start-rate limiting + concurrency governor, History
  paging for UIDs, resumable batch workspace, `rettype=gp` invariance check.
- A parallel, rewritten history of the same work exists on the branches
  `backup/original-messages` and `(perf)-improve-ncbi-fetch-performance`; the
  motivation for the rewrite is not documented.

## Previous state

UID pagination without History; in-memory consolidation; brittle partial-response
handling; no bounded concurrency.

## Decision

Keep the History handle for UID retrieval; stream persistence; larger XML
batches; selective resume; bounded concurrency (default off); reclassify a
truncated response as a transient failure with retry; order-safe `latest/`
publish; add an `rettype=gp` payload-invariance test.

## Rationale

Stated in the plan document (findings F1–F12 each name a concrete defect and its
fix). Verifiable against the staging code.

## Alternatives considered

The plan's §8 "open decisions" lists options that were left open. Beyond that,
`HISTORICAL_EVIDENCE_INSUFFICIENT`.

## Consequences

- The Phase A run (ADR-0004) used this retrieval path.
- Its relationship to the PR #27 failure controls already in `master`
  (ADR-0002) must be reconciled at reintegration — overlap vs supersession.
- The plan's "open decisions" were never explicitly closed.

## Limitations

Not integrated; not reviewed via a PR; the branch-rewrite provenance is unclear.

## Supersedes

The pre-rework UID-pagination and consolidation behaviour, **on the staging
branch only**. `master` still runs the PR #27 controls.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found` — committed directly to the staging branch.

## Related commits

`f5671b8`, `50016f0`, `24cba76`, `c2837df`, `22a967e`, `4598a25`, `ba62ceb`,
`df4eed4`, `0b4ed5a`.

## Related data / artifacts

`scripts/ncbi_rettype_invariance_check.py`; the two staging performance
documents.

## Validation

Per the implementation document: the `rettype=gp` invariance test (2 NCBI
requests) and byte-identical consolidation at full scale. Not independently
re-verified here.

## Scientific impact

None. Retrieval mechanics only.

## Data impact

Enables the annotation-enriched acquisition (ADR-0004); no label or partition
change.

## Reproducibility impact

Positive on the staging branch: History paging + streaming + resume make a large
retrieval reproducible and restartable.

## Historical reconstruction note

- Directly demonstrated: the staging commits and the two staging documents.
- Inferred: that this is one coherent effort distinct from the integrated
  PR #27 controls.
- Not recoverable: why the parallel branch history was rewritten; the closure of
  the "open decisions".
