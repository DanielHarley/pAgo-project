# Operation records

An Operation Record documents a factual, non-scientific, non-architectural
operation on the repository or its provenance — a backup, a migration, a
recovery, an incident, a staging-safety change, or an administrative action
whose operational detail (before/after state, hashes, paths, verification
steps) is worth preserving but does not belong in a more compact document.

## When to write one

Write an Operation Record when the operation:

- involves concrete forensic detail (hashes, paths, before/after state,
  verification steps) that would make `PROJECT_LEDGER.md` unwieldy if inlined
  there in full;
- is not itself a decision about project methodology, architecture, or a
  scientific conclusion — that is an ADR (`decisions/`);
- does not, by itself, change what the project currently considers
  scientifically valid — that is `SCIENTIFIC_STATE.md`;
- is not a scientific or computational execution whose result matters for a
  conclusion — that is a Run Record plus its manifest.

Do not write one for a routine `git status` check, a routine commit with an
obvious diff, or anything a single `PROJECT_LEDGER.md` line already covers on
its own. Most day-to-day work needs no Operation Record.

## How it differs from the other documents

| Document | Answers |
| --- | --- |
| `decisions/*.md` (ADR) | Why was a significant technical / scientific / architectural decision made, over what alternatives? |
| `PROJECT_LEDGER.md` | When did this happen, and what decision or event does it belong to — one compact, chronological entry? |
| `operations/*.md` (this directory) | What exactly happened during an important operational event — the full forensic detail: before/after state, discoveries, verification, hashes? |
| `SCIENTIFIC_STATE.md` | What does the project currently consider scientifically valid, right now? |
| A Run Record + manifest | What execution produced a specific artifact or result, from what inputs, verified how? |

An Operation Record never replaces an ADR, a `PROJECT_LEDGER.md` entry, a
commit, or a Pull Request — it sits alongside them, and the corresponding
`PROJECT_LEDGER.md` entry links to it.

## Format

- File name: `YYYY-MM-DD-semantic-kebab-title.md` — semantic, not a roadmap
  step (same naming rule as branches/commits/PRs; see
  `documentation/README.md#naming-convention`).
- Plain factual record: purpose, initial state, what was found, what was
  done, verification performed, Git provenance (branch, commit), and explicit
  scope limits — what the operation does **not** claim (e.g. it does not
  integrate, does not validate, does not itself constitute a backup or a
  scientific result).
- Never include the content of files being backed up or protected, holdout
  data, or scientific results in the record itself — reference paths and
  hashes, not payloads.
- A historical roadmap reference (e.g. `B5/B6`) may appear as metadata only,
  never as the record's primary identity.

## Related

See `documentation/README.md` for how this fits the rest of the decision
chain, and `documentation/decisions/README.md` for when something is an ADR
instead.
