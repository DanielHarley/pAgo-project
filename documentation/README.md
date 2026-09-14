# pAgo project documentation

This directory contains the project's tracked documentation.

Its purpose is to preserve the information a future researcher or developer
would actually need to understand the project: important decisions, current
scientific interpretation, significant historical evidence, and exceptional
operational events.

The documentation is not a second state machine for the repository. Git records
repository state. Documentation explains meaning.

## Documentation principles

### Minimal sufficient documentation — MUST

Documentation MUST contain only information that is materially useful to a
future reader.

Do not add fields, sections, statuses, classifications, metadata, or records
merely because the documentation system allows them.

If removing a field or section would not meaningfully reduce understanding,
traceability, reproducibility, or scientific interpretation, omit it.

### Document information, not documentation state — MUST

Do not maintain parallel lifecycle, integration, or epistemic state machines
inside documentation.

Repository state belongs to Git. Scientific interpretation should be stated in
plain language in the scientific documentation itself.

Do not duplicate facts that can be obtained reliably from Git unless the fact
is necessary to understand the document.

### Structure is optional unless meaning requires it — MUST

Templates and section lists are guides, not mandatory schemas.

Use only the sections needed to explain the specific subject.

Never add empty, trivial, redundant, speculative, or `N/A` sections merely to
satisfy a template.

### Preserve semantic distinctions — MUST

Being merged into Git does not establish scientific validity.

Historical evidence does not automatically define current project behaviour.

Exploratory results must not be described as validated classification.

Incomplete evidence must be described as incomplete.

A reproducible artifact is not automatically a scientifically valid model or
biological conclusion.

### Prefer durable information — MUST

Documentation should preferentially record information that remains useful as
the repository evolves.

Avoid transient operational facts unless they are necessary to reconstruct an
important event.

An ADR should document a decision, not continuously describe the repository
state around that decision.

### Do not create documentation recursively — MUST

Do not create or update documentation merely to record that another
documentation-only change was merged or updated.

Git history is sufficient for ordinary documentation bookkeeping.

### Use semantic names — MUST

Names must describe what a change, decision, record, or artifact means.

Do not use roadmap identifiers such as `Phase A`, `Phase B`, `B1`, `B4.3`, or
`B5` as the primary semantic identity of branches, commits, PRs, ADRs, operation
records, or other tracked artifacts.

Historical roadmap identifiers may appear as secondary provenance metadata when
that genuinely helps reconstruct earlier work.

## Map

| Path | Purpose |
| --- | --- |
| [`SCIENTIFIC_STATE.md`](SCIENTIFIC_STATE.md) | Current scientific interpretation, important boundaries, unresolved questions, and known limitations. |
| [`PROJECT_LEDGER.md`](PROJECT_LEDGER.md) | Compact chronological record of important decision events: how the project got here. |
| [`decisions/`](decisions/) | ADRs for significant technical, methodological, or scientific decisions whose rationale is worth preserving. |
| [`history/`](history/) | Point-in-time historical evidence preserved for provenance. It is not rewritten to match later project state. |
| [`operations/`](operations/) | Detailed records of unusual operational events whose forensic detail is worth preserving beyond normal Git history. |
| [`pago_qc.md`](pago_qc.md) | Documentation of the integrated legacy metadata-derived QC workflow and its scientific boundaries. |
| [`git_conventions.md`](git_conventions.md) | Authority for branch, commit, Pull Request, and Issue naming and writing conventions. |

There is currently no authoritative `pago_annotation_ontology.md` in
`documentation/`. The staging draft is preserved under `history/`, while the
relevant semantic decisions are recorded in the ADRs.

## How to decide what to update

Use the smallest documentation surface that preserves the information that
matters.

- If a significant decision needs its reasoning preserved, use an ADR.
- If an important event belongs in the project's chronology, add a compact
  `PROJECT_LEDGER.md` entry.
- If the project's current scientific interpretation changes, update
  `SCIENTIFIC_STATE.md`.
- If an unusual operational event needs detailed forensic provenance, use an
  Operation Record.
- If the information is already adequately preserved by the commit and PR, do
  not create another document.
- When creating a branch, commit, PR, or Issue, consult
  `git_conventions.md` for naming and writing conventions.

Most ordinary code changes do not require changes to every documentation layer.

## ADRs

ADRs inherit the principles in this README.

The ADR template is intentionally minimal. `Context`, `Decision`, `Rationale`,
and `Consequences` are often sufficient. Add `Evidence`, `Alternatives`,
`Validation`, `Limitations`, `Supersedes`, historical reconstruction notes, or
other sections only when they materially improve the record.

See [`decisions/README.md`](decisions/README.md).

## History

Files under `history/` are historical evidence only.

Do not treat a statement from the middle of a historical document as current
project authority without checking current code and current documentation.

Do not rewrite a historical body merely because later terminology or decisions
changed. Later corrections belong in current documents.

See [`history/README.md`](history/README.md).

## Operations

Operation Records are exceptional, not routine.

Create one only when an operational event contains provenance, recovery,
backup, migration, incident, or verification detail that is important enough to
preserve beyond the commit and PR.

See [`operations/README.md`](operations/README.md).

## Git and scientific meaning

Git answers questions such as which commit exists on `master`, which PR merged,
and what a branch contains.

Scientific documentation answers different questions: what a dataset means,
what a result does and does not establish, why a methodological decision was
made, and what remains unresolved.

Do not reproduce Git state as documentation metadata merely to keep a second
copy of the same fact synchronized.

## Naming convention

Branch, commit, PR, and Issue naming and writing conventions are defined in
[`git_conventions.md`](git_conventions.md).
