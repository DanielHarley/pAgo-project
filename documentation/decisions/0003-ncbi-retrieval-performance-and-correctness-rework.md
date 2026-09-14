# ADR-0003 — Rework NCBI retrieval for measured performance and correctness constraints

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commits
dated 2026-08-15)

## Context

Large NCBI E-utilities retrievals were slow and fragile: UID paging discarded
the Entrez History handle, XML consolidation held large payloads in memory,
partial responses had weak recovery paths, and request behaviour needed clearer
bounds.

## Decision

Rework large-record retrieval to:

- retain the Entrez History handle for UID paging
- persist XML incrementally instead of consolidating everything in memory
- use larger batches where measured behaviour supports them
- support selective resume
- provide bounded concurrency with conservative defaults
- treat truncated responses as transient failures eligible for retry
- publish `latest/` only after a complete result is ready
- verify the `rettype=gp` payload contract with an invariance check

## Evidence

Historical commits `f5671b8` … `0b4ed5a` implement the measured retrieval
rework, including streaming XML handling, request telemetry, rate limiting,
History paging, resumable batch workspaces, and the `rettype=gp` check.

Two historical performance documents preserved on the staging branch recorded
the measured baseline, findings F1–F12, phased plan P0–P5, acceptance criteria,
and defects found while exercising the real workflow.

## Rationale

The historical plan ties the individual changes to measured bottlenecks and
failure modes. The code provides direct evidence for the implemented mechanics.

## Consequences

- later large acquisitions can reuse the History-backed, resumable retrieval
  path
- the historical implementation overlaps conceptually with failure controls
  introduced later in PR #27, so code should not be reused blindly without
  reconciling duplicated or superseding behaviour

## Limitations

The original plan left some design choices explicitly open, and the reason for a
parallel rewritten branch history was not documented.

## Validation

The historical implementation records a `rettype=gp` invariance check and
byte-identical consolidation at full scale. These results are historical
evidence and were not independently rerun during the retrospective ADR
reconstruction.

## Related records

- commits `f5671b8`, `50016f0`, `24cba76`, `c2837df`, `22a967e`, `4598a25`,
  `ba62ceb`, `df4eed4`, `0b4ed5a`
- historical roadmap reference: pre-Phase-A
- `scripts/ncbi_rettype_invariance_check.py`
- historical staging documents `ncbi_retrieval_performance_plan.md` and
  `ncbi_retrieval_performance_implementation.md`

## Historical reconstruction note

- Directly demonstrated: historical commits and staging performance documents.
- Inferred: these changes form one coherent retrieval-rework effort distinct
  from later PR #27 controls.
- Not recoverable: why the parallel branch history was rewritten and how every
  open decision in the plan was ultimately resolved.
