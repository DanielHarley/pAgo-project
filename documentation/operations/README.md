# Operation records

Operation Records inherit the documentation principles in
[`../README.md`](../README.md).

An Operation Record preserves the factual detail of an unusual repository or
provenance operation when that detail is important enough to outlive the commit
and PR.

Typical examples include:

- a backup or recovery operation
- a repository migration
- an incident with meaningful forensic detail
- protection of important local files from accidental staging
- an administrative change whose before/after state and verification matter

Most day-to-day work does **not** need an Operation Record.

## When to write one

Create an Operation Record only when the event contains details such as paths,
hashes, before/after state, recovery evidence, or verification steps that a
future maintainer may genuinely need.

Do not create one for a routine `git status`, an ordinary commit, a normal PR,
or any event already explained adequately by Git history and a compact Project
Ledger entry.

An Operation Record is not an ADR. It records what happened during an operation;
it does not manufacture a methodological or scientific decision.

An Operation Record is not a substitute for scientific interpretation. If an
event changes what the project can scientifically claim, update the appropriate
scientific documentation separately.

## Content

Use only the sections required by the event.

Useful content may include:

- purpose
- relevant initial state
- what was found
- what was done
- verification performed
- important paths or hashes
- Git provenance when it materially helps reconstruction
- explicit scope limits when a future reader could otherwise overinterpret the
  operation

Do not add empty sections merely to follow a template.

## Naming

Use:

`YYYY-MM-DD-semantic-kebab-title.md`

The name must describe the operation itself. Historical roadmap identifiers may
appear inside the record as provenance metadata when useful, but not as the
record's primary identity.

## Sensitive or unnecessary payloads

Reference paths and hashes when they are sufficient. Do not copy holdout data,
large backed-up payloads, or unrelated scientific results into an Operation
Record merely for completeness.

## Relationship to the rest of the documentation

- [`PROJECT_LEDGER.md`](../PROJECT_LEDGER.md) may contain a compact chronological
  pointer to an important operation.
- [`decisions/`](../decisions/) records significant decisions and their
  rationale.
- [`SCIENTIFIC_STATE.md`](../SCIENTIFIC_STATE.md) records current scientific
  interpretation and boundaries.

Use the smallest set of documents that preserves the information that matters.
