# ADR-NNNN — <semantic title: what was decided>

This template is a guide, not a mandatory schema.

Delete any section that does not materially help explain the decision. Do not
add empty, trivial, redundant, speculative, or `N/A` sections merely to make an
ADR look complete.

For many decisions, `Context`, `Decision`, `Rationale`, and `Consequences` are
sufficient.

## Date

Include the decision date when it is useful and recoverable.

For a retrospective record, say that the date is reconstructed if the original
date is uncertain.

## Context

What situation, constraint, or scientific/technical need made the decision
necessary?

## Decision

State the decision clearly enough that a future contributor can tell whether a
proposed change follows or contradicts it.

## Rationale

Why was this option chosen?

Do not invent rationale. If the original rationale cannot be recovered, write
`RATIONALE_NOT_RECOVERABLE` and preserve only what the surviving evidence
supports.

## Consequences

What important constraints, requirements, or trade-offs follow from the
decision?

---

The sections below are optional. Add only the ones that materially improve the
record.

## Problem

Use when the specific problem is not already clear from `Context`.

## Evidence

Use when evidence is central to understanding or defending the decision.
Separate observed or executed evidence from interpretation.

## Previous state

Use when the earlier behaviour or interpretation is necessary to understand the
change.

## Alternatives considered

Use only for alternatives that are actually supported by the historical or
contemporaneous record.

If historical evidence cannot reconstruct them, write
`HISTORICAL_EVIDENCE_INSUFFICIENT` rather than manufacturing alternatives.

## Limitations

Use when the decision has important scientific, technical, or methodological
boundaries that a future reader could otherwise overinterpret.

## Validation

Use when executed validation or an explicit validation boundary is important to
the decision.

Never describe planned validation as executed validation.

## Supersedes

Use only when this decision actually replaces a specific earlier decision.
Link the earlier ADR or describe the earlier rule precisely.

## Related records

Use only when Issues, PRs, commits, artifacts, reports, or other records
materially improve traceability.

Do not duplicate repository state merely to maintain a second copy of what Git
already records.

## Historical reconstruction note

Use for retrospective ADRs when it helps preserve the difference between:

- what is directly demonstrated by surviving evidence
- what is inferred
- what could not be recovered

Do not rewrite historical Git history or silently turn an inference into a
historical fact.
