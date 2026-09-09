# Decision records (ADRs)

An ADR (Architecture Decision Record) captures **one** significant decision:
the problem, the evidence, the options, the choice, the rationale, and the
provenance. In this project ADRs also cover **scientific and methodological**
decisions (label definitions, ontology, phylogenetic reference choices,
partition / holdout policy, evidence criteria), not only software architecture.

## When to write an ADR

Write an ADR when the project changed its behaviour, interpretation,
architecture, methodology, or a scientific rule in a way that a future
maintainer would need explained.

Do **not** write an ADR for:

- typos, formatting, mechanical renames, comments, pure refactors — Git is
  enough;
- small technical changes with an obvious rationale already in the commit and
  the PR — a `PROJECT_LEDGER.md` line pointing at the PR is enough.

One ADR per decision. Never one ADR per commit. Group commits that together
implement a single conceptual change.

## Numbering and files

`NNNN-semantic-kebab-title.md`, four digits, zero-padded, allocated in order.
`0000-adr-template.md` is the template and is not a decision.

The title states the decision, not the roadmap step. Prefer
`0017-split-group-as-indivisible-train-test-unit.md` over `0017-b2-partitions.md`.

## Status axes

ADRs carry up to **three independent** fields. No field mixes two axes. See
[`README.md`](../README.md#status-axes) for the full definitions.

`Lifecycle status` is the decision lifecycle / governance status — where the
decision stands, independently of Git:

| Value | Meaning |
| --- | --- |
| `CURRENT` | The operative decision for its topic. |
| `PROPOSED` | Drafted, not yet acted on. |
| `HISTORICAL` | Made and recorded, no longer driving current work, no direct successor. |
| `SUPERSEDED` | Replaced by a specific later decision (link it under *Superseded by*). |
| `EXPERIMENTAL` | A trial that was not adopted as the operative decision. |
| `DEPRECATED` | Withdrawn, not replaced. |

`Integration status` is a fact about Git:

| Value | Meaning |
| --- | --- |
| `INTEGRATED` | Merged into `master` (add `(via PR #N)`). |
| `NOT_INTEGRATED` | The implementation exists only on a staging / historical branch. |
| `WAS_INTEGRATED_REMOVED` | Merged into `master` and later removed. |
| `N/A` | Documentation-only or process decision. |

`Epistemic status` (for `SCIENTIFIC` and `TECHNICAL + SCIENTIFIC` ADRs only) is
whether the decision itself is sound, regardless of Git:
`SUPPORTED` / `PARTIALLY_SUPPORTED` / `UNRESOLVED` / `N/A`.

A retrospectively reconstructed decision that currently lives only on the
staging branch and is well evidenced is:

```
## Lifecycle status
CURRENT

## Integration status
NOT_INTEGRATED

## Epistemic status
SUPPORTED
```

After the work is reintegrated through a reviewed PR, only the integration axis
moves:

```
## Integration status
INTEGRATED (via PR #N)
```

The `## Historical source` section is preserved across that transition so the
decision that *existed on the staging branch* stays distinguishable from the
decision *accepted on the integration line*.

## Required sections

See [`0000-adr-template.md`](0000-adr-template.md). Every ADR must include the
**Historical reconstruction note** when it is written retrospectively: it must
say which parts are directly demonstrated by history, which are inferred, and
which could not be recovered. Use `RATIONALE_NOT_RECOVERABLE` or
`HISTORICAL_EVIDENCE_INSUFFICIENT` rather than inventing a rationale, an
alternative, or an author.

## Do not rewrite history

Retrospective ADRs describe decisions that were already made. Do not attribute
a decision to a person without evidence (`git blame` shows who committed a
line, not who made the decision). Do not turn a consequence observed later into
an original justification. Preserve the difference between "it was decided this
way at the time" and "we now know this is the correct interpretation".

## Relationship to Issues and PRs

ADRs link to Issues and PRs; they do not copy them. For decisions with no
historical Issue or PR (much of the pre-framework history), write
`No historical Issue found` / `No historical PR found` rather than
manufacturing one. An audit Issue created now is allowed, but must be labelled
as retrospective.
