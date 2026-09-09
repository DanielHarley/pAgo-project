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

## Status and integration status

ADRs carry **two** independent fields.

`Status` is the epistemic lifecycle of the decision:

| Status | Meaning |
| --- | --- |
| `Proposed` | Under discussion. |
| `Historical decision` | The decision was genuinely made and acted on, on a historical / staging branch, and has not yet been re-accepted onto `master` through review. |
| `Accepted` | Accepted on the integration line. |
| `Superseded` | Replaced by a later decision (link it under *Superseded by*). |
| `Deprecated` | No longer relevant, not replaced. |

`Integration status` is a statement about Git:

| Integration status | Meaning |
| --- | --- |
| `NOT INTEGRATED INTO MASTER` | The decision's implementation exists only on a historical / staging branch. |
| `Integrated into master via PR #N` | Merged. |
| `N/A` | Documentation-only or process decision. |

A retrospectively reconstructed decision that currently lives only on
`(feat)-siepe-ready-project` therefore starts as:

```
## Status
Historical decision

## Integration status
NOT INTEGRATED INTO MASTER
```

and, after the work is reintegrated through a reviewed PR, becomes:

```
## Status
Accepted

## Integration status
Integrated into master via PR #N
```

The `## Historical source` section is preserved across that transition so the
decision that *existed* stays distinguishable from the decision that has been
*accepted on the integration line*.

`SCIENTIFIC` and `TECHNICAL + SCIENTIFIC` ADRs also carry a third, independent
field, `Epistemic status` (`SUPPORTED` / `PARTIALLY_SUPPORTED` / `UNRESOLVED` /
`N/A`): is the decision itself scientifically sound, regardless of whether it is
merged? A staging-branch decision can be `Historical decision`,
`NOT INTEGRATED INTO MASTER`, and `SUPPORTED` at once.

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
