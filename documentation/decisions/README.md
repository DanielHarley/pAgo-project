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
| `N/A` | No meaningful Git integration state applies to this record. |

`Epistemic status` is whether the decision itself is sound, regardless of Git:
`SUPPORTED` / `PARTIALLY_SUPPORTED` / `UNRESOLVED` / `N/A`. Use it for
`SCIENTIFIC` and `TECHNICAL + SCIENTIFIC` ADRs, and also for `TECHNICAL` ADRs
when there is an evidence-bearing claim about technical soundness or validation.
`N/A` when no epistemic / evidential judgment is meaningful.

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

When the work is reintegrated through a reviewed PR: if the same decision is
accepted unchanged, only the integration axis moves —

```
## Integration status
INTEGRATED (via PR #N)
```

If review changes the decision, its rationale, or its evidential support,
update the relevant axes and record a revised or superseding decision
explicitly (a new ADR under *Superseded by* / *Supersedes*).

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

## Index

Retrospectively reconstructed decisions. `Lc` = lifecycle, `Int` = integration,
`Ep` = epistemic status.

| ADR | Decision | Type | Lc | Int | Ep |
| --- | --- | --- | --- | --- | --- |
| [0001](0001-canonical-ncbi-retrieval-identifier-is-the-protein-uid.md) | Canonical NCBI retrieval identifier is the protein UID | T | CURRENT | INTEGRATED | SUPPORTED |
| [0002](0002-immutable-snapshot-and-provenance-architecture.md) | Immutable snapshot + manifest + SHA-256 provenance architecture | T | CURRENT | INTEGRATED | SUPPORTED |
| [0003](0003-ncbi-retrieval-performance-and-correctness-rework.md) | NCBI retrieval performance & correctness rework (F1–F12 / P0–P5) | T | CURRENT | NOT_INTEGRATED | SUPPORTED |
| [0004](0004-annotation-enriched-candidate-set-and-technical-only-prefilter.md) | Annotation-enriched candidate set; technical-only prefilter | T+S | CURRENT | NOT_INTEGRATED | SUPPORTED |
| [0005](0005-query-recall-reference-panel-and-sequence-identity-matching.md) | Query-recall reference panel; sequence-identity matching | S | CURRENT | NOT_INTEGRATED | SUPPORTED |
| [0006](0006-sweep-pca-kmeans-is-exploratory-not-validated-classification.md) | SWeeP → PCA → KMeans is exploratory, not validated classification | S | CURRENT | INTEGRATED | PARTIALLY_SUPPORTED |
| [0007](0007-provenance-and-defensible-labels-before-supervised-ml.md) | Provenance / defensible labels / protected holdout before supervised ML | S | CURRENT | NOT_INTEGRATED | PARTIALLY_SUPPORTED |
| [0008](0008-pinned-pfam-38-2-hmm-bundle.md) | Pinned Pfam 38.2 HMM bundle (10 profiles, versioned, SHA-256) | T+S | CURRENT | NOT_INTEGRATED | PARTIALLY_SUPPORTED |
| [0009](0009-apaz-s3-reference-and-high-similarity-split-group-partitioning.md) | APAZ S3 reference; high-similarity split-group partitioning (v2 supersedes v1) | T+S | CURRENT | NOT_INTEGRATED | PARTIALLY_SUPPORTED |
| [0010](0010-apaz-profile-hmms-built-from-build-only.md) | Six APAZ profile HMMs built from BUILD only | T+S | CURRENT | NOT_INTEGRATED | PARTIALLY_SUPPORTED |
| [0011](0011-ryazansky-1010-pago-catalog-and-mid-piwi-reference-extraction.md) | Ryazansky 1010-pAgo catalog reconciliation; MID-PIWI extraction | T+S | CURRENT | NOT_INTEGRATED | SUPPORTED |
| [0012](0012-annotation-ontology-family-versus-clade.md) | Ontology: `ago_family` vs `pago_clade`; PIWI-RE is a family | S | CURRENT | NOT_INTEGRATED | SUPPORTED |
| [0013](0013-mid-piwi-high-similarity-split-groups.md) | MID-PIWI high-similarity split groups over the Ryazansky catalog | T+S | CURRENT | NOT_INTEGRATED | SUPPORTED |
| [0014](0014-independent-piwi-re-reference-arm.md) | Independent PIWI-RE reference arm (proposed, not executed) | S | PROPOSED | NOT_INTEGRATED | UNRESOLVED |
