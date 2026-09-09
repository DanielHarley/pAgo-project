# ADR-NNNN — <semantic title: what was decided>

The three status fields are independent. No field mixes two axes.

## Lifecycle status

Decision lifecycle / governance status:

`CURRENT` | `PROPOSED` | `HISTORICAL` | `SUPERSEDED` | `EXPERIMENTAL` | `DEPRECATED`

## Integration status

A fact about Git:

`INTEGRATED (via PR #N)` | `NOT_INTEGRATED` | `WAS_INTEGRATED_REMOVED` | `N/A`

## Epistemic status

Whether the decision itself is sound. For `SCIENTIFIC` and `TECHNICAL +
SCIENTIFIC` decisions only; otherwise `N/A`.

`SUPPORTED` | `PARTIALLY_SUPPORTED` | `UNRESOLVED` | `N/A`

A staging-branch decision can be `CURRENT` (lifecycle), `NOT_INTEGRATED` (Git),
and `SUPPORTED` (evidence) at the same time.

## Type

`TECHNICAL` | `SCIENTIFIC` | `TECHNICAL + SCIENTIFIC`

## Date

Original date of the decision, if recoverable (ISO `YYYY-MM-DD`).
Otherwise: `UNKNOWN — reconstructed retrospectively on YYYY-MM-DD`.

## Context

The situation and constraints in which the decision was taken.

## Problem

What problem existed. When it was identified.

## Evidence

What evidence demonstrated the problem (tests, artifacts, manifests, hashes,
literature, logs, prior audits). Cite commits / files / line ranges.

## Previous state

How the project behaved before this decision.

## Decision

What was changed, precisely.

## Rationale

Why this option was chosen over the others.

## Alternatives considered

Each alternative and why it was not chosen. If none are recoverable from the
historical record, write `HISTORICAL_EVIDENCE_INSUFFICIENT` — do not invent
alternatives.

## Consequences

What became different. Positive and negative.

## Limitations

Known limits of the decision as made.

## Supersedes

Earlier decision(s) this one replaces. Preserve, do not delete, the old
decision: OLD STATE -> EVIDENCE -> NEW DECISION -> SUPERSEDED.

## Superseded by

Later decision(s) that replace this one, once they exist.

## Related Issue

Issue link, or `No historical Issue found`.

## Related PR

PR link, or `No historical PR found`.

## Related commits

SHAs (preserved; history is not rewritten).

## Related data / artifacts

Reference resources, snapshots, locks, fixtures affected.

## Historical source

Only for decisions reconstructed from a historical / staging branch.

```
Branch:
Commit(s):
Historical roadmap reference:   (e.g. B4.3 — secondary metadata only)
```

## Validation

How we know it worked (tests, reproducibility checks, holdout reports,
verifier runs).

## Scientific impact

Did it change a biological interpretation, an ontology field, a label, a
partition, or a leakage risk? If not, say so.

## Data impact

Did it change data, labels, partitions, or artifacts?

## Reproducibility impact

Did it change what can be reconstructed, and how?

## Historical reconstruction note

**Mandatory for retrospective ADRs.** State explicitly:

- which parts of this record are directly demonstrated by the Git / artifact
  history;
- which parts are inferred by the author;
- which parts could not be recovered (`RATIONALE_NOT_RECOVERABLE`,
  `HISTORICAL_EVIDENCE_INSUFFICIENT`).
