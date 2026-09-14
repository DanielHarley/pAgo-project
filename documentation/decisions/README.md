# Decision records (ADRs)

ADRs inherit the documentation principles in [`../README.md`](../README.md).

An ADR (Architecture Decision Record) captures one significant decision whose
reasoning is worth preserving. In this project that includes technical,
methodological, and scientific decisions — not only software architecture.

The goal is not to classify every decision. The goal is to preserve enough
context that a future reader can understand what was decided and why.

Write ADR content in English.

## When to write an ADR

Write an ADR when a future maintainer or researcher would reasonably need to
know why the project chose one important approach over another.

Typical examples include:

- a canonical identifier or persistence contract
- a dataset-construction rule
- a leakage-control or holdout policy
- a scientific annotation boundary
- a reference-resource choice
- a methodology that constrains later work
- a change in how an existing result may be interpreted

Do not write an ADR merely because there is a commit or PR.

Do not write one for typos, formatting, mechanical renames, routine dependency
updates, ordinary tests, small refactors, repository hygiene, or other changes
whose meaning is already clear from Git history.

One ADR records one conceptual decision, not one commit.

## Minimal sufficient ADRs — MUST

ADR sections are not mandatory merely because they exist in a template.

Use only the sections needed to explain the specific decision.

For many decisions, this is enough:

- Context
- Decision
- Rationale
- Consequences

Add `Evidence`, `Problem`, `Previous state`, `Alternatives considered`,
`Validation`, `Limitations`, `Supersedes`, `Related records`, or a historical
reconstruction note only when they contain information that materially improves
the record.

Never add an empty, trivial, redundant, speculative, or `N/A` section just to
satisfy a documentation shape.

## Do not maintain ADR status taxonomies

ADRs do not carry lifecycle, integration, or epistemic status fields.

Git records repository state. Scientific confidence and limitations should be
stated directly in the relevant prose when they matter.

Do not add substitute fields that recreate the same state machine under new
names.

For example, instead of a generic label such as `PARTIALLY_SUPPORTED`, state the
actual boundary: a resource may be byte-reproducible while predictive
sensitivity and specificity remain unevaluated.

## Titles and filenames

Use:

`NNNN-semantic-kebab-title.md`

Numbers are four-digit, zero-padded, and allocated in order.
`0000-adr-template.md` is a template and is not a decision.

Do not renumber existing ADRs. Do not reuse a number that has already
appeared in repository history, even for a deleted, historical, superseded,
or abandoned ADR.

The title must state the decision itself, not a roadmap step or vague topic.

Prefer:

`0017-keep-high-similarity-groups-within-one-evaluation-partition.md`

instead of:

`0017-b2-partitions.md`

Historical roadmap identifiers may appear inside the record as provenance when
useful, but never as the primary semantic identity.

## Evidence and scientific claims

Evidence should be as specific as necessary to support the decision.

Do not invent evidence, alternatives, historical intent, validation results,
or decision provenance. If evidence is genuinely missing, say so explicitly
(see `RATIONALE_NOT_RECOVERABLE` and `HISTORICAL_EVIDENCE_INSUFFICIENT`
below) rather than reconstructing it from guesswork.

Distinguish what was observed, computed, inferred, decided, validated, and left
unresolved. Do not silently collapse those categories.

Examples:

- a profile HMM parsing successfully does not establish predictive performance
- a byte-identical sequence under two accessions does not establish identical
  biological provenance
- a high-similarity split group is a redundancy-control unit, not a definition
  of homology
- reproducibility does not by itself establish scientific validity

If the original rationale cannot be recovered, use
`RATIONALE_NOT_RECOVERABLE` rather than inventing one.

If the historical evidence is insufficient to reconstruct an alternative,
motivation, or event, use `HISTORICAL_EVIDENCE_INSUFFICIENT`.

## Retrospective ADRs

A retrospective ADR reconstructs a decision that predates the ADR system.

Preserve the distinction between:

- what surviving evidence directly demonstrates
- what is a reasonable inference
- what cannot be recovered

Do not rewrite historical Git history. Do not attribute a decision to a person
or tool without evidence. Do not turn a later consequence into an original
rationale.

A historical branch, commit SHA, roadmap identifier, or old document may be
recorded as a source when it is genuinely useful for reconstructing the
decision. It does not need to be repeated as ongoing repository-state metadata.

## Supersession

When a later decision replaces an earlier one, link the two decisions and
explain the replacement where that relationship matters.

Do not add `Supersedes: None` or `Superseded by: None` sections.

Do not delete the earlier ADR merely because a later decision replaces it.

## Relationship to Git

ADRs may link to Issues, PRs, commits, or artifacts when those links materially
help trace the decision.

They do not need to continuously state whether an implementation is currently
merged, staged, or removed. Git is the source for that information.

Do not create a follow-up documentation PR merely to update an ADR after a
normal merge.

## Index

| ADR | Decision |
| --- | --- |
| [0001](0001-canonical-ncbi-retrieval-identifier-is-the-protein-uid.md) | Use the protein UID as the canonical NCBI retrieval identifier |
| [0002](0002-immutable-snapshot-and-provenance-architecture.md) | Use immutable snapshots, manifests, and SHA-256 provenance |
| [0003](0003-ncbi-retrieval-performance-and-correctness-rework.md) | Rework NCBI retrieval for measured performance and correctness constraints |
| [0004](0004-annotation-enriched-candidate-set-and-technical-only-prefilter.md) | Treat the annotation-enriched query result as a candidate set and keep the prefilter technical only |
| [0005](0005-query-recall-reference-panel-and-sequence-identity-matching.md) | Use an evidence-tiered query-recall panel with sequence-identity equivalence matching |
| [0006](0006-sweep-pca-kmeans-is-exploratory-not-validated-classification.md) | Treat SWeeP → PCA → KMeans as exploratory analysis, not validated biological classification |
| [0007](0007-provenance-and-defensible-labels-before-supervised-ml.md) | Require provenance, defensible labels, redundancy control, and a protected holdout before supervised ML |
| [0008](0008-pinned-pfam-38-2-hmm-bundle.md) | Pin the Pfam 38.2 HMM reference bundle |
| [0009](0009-apaz-s3-reference-and-high-similarity-split-group-partitioning.md) | Partition APAZ references by high-similarity split groups |
| [0010](0010-apaz-profile-hmms-built-from-build-only.md) | Build APAZ profile HMMs from BUILD references only |
| [0011](0011-ryazansky-1010-pago-catalog-and-mid-piwi-reference-extraction.md) | Reconcile the Ryazansky 1010-pAgo catalog and extract MID-PIWI references |
| [0012](0012-annotation-ontology-family-versus-clade.md) | Separate Argonaute family identity from pAgo phylogenetic clade |
| [0013](0013-mid-piwi-high-similarity-split-groups.md) | Use MID-PIWI high-similarity split groups as future partition units |
| [0014](0014-independent-piwi-re-reference-arm.md) | Require an independent PIWI-RE reference arm for detector validation |
