# pAgo project documentation

This directory holds the project's tracked, versioned documentation: the
decision history, the current scientific state, the annotation ontology, and
the architecture / methodology decision records.

The goal is that a clone of this repository contains the full chain needed to
answer, for any significant change:

1. what problem existed
2. when it was identified
3. what the previous state was
4. what evidence demonstrated the problem
5. what alternatives were considered (where recoverable)
6. what was decided, and why
7. how it was implemented and validated
8. what earlier result it invalidated or made historical
9. which decision is authoritative today

## Map

| Document | Question it answers |
| --- | --- |
| [`SCIENTIFIC_STATE.md`](SCIENTIFIC_STATE.md) | What is considered scientifically valid **right now**, and what is only staged or historical? |
| [`PROJECT_LEDGER.md`](PROJECT_LEDGER.md) | How did we get here? A dated list of decision events, each pointing to more detailed sources. |
| [`decisions/`](decisions/) | Architecture / methodology Decision Records (ADRs) — one per decision, with full rationale and provenance. Start at [`decisions/README.md`](decisions/README.md). |
| `pago_annotation_ontology.md` | The formal semantics of every annotation field (family, clade, detection, architecture, evidence tiers). **Not yet present as an authoritative document** — see [The staging ontology](#the-staging-ontology). |
| [`pago_qc.md`](pago_qc.md) | The legacy metadata-derived QC labelling layer that is integrated in `master`. Its labels are audit features, not ground truth. |
| `ncbi_retrieval_performance_plan.md` / `ncbi_retrieval_performance_implementation.md` | The measured baseline, findings, and phased rework of NCBI retrieval. Staging sources not yet in `master`; preserved / integrated by the retrospective backfill, which turns these into real links. |
| `history/` | Point-in-time audits and staging documents preserved verbatim as provenance. Never authoritative for current behaviour. Populated by the retrospective backfill PR. |

## How Git and this documentation relate

The decision comes before the implementation. An ADR may be written on the same
branch and approved in the same PR as the code, but it records a choice that was
reasoned through — not a rationalisation written after the diff already exists.

```
problem / Issue          what problem exists, and the evidence for it
  |
alternatives             the options weighed
  |
decision                 what was chosen, and why
  |
ADR (when required)      the decision recorded
  |
branch + implementation  the work
  |
commit(s)                each logical change
  |
validation               tests, reproducibility checks, holdout reports
  |
Pull Request / review    the complete decision proposed for integration
  |
merge into master        integration
  |
PROJECT_LEDGER           where it sits in the project's evolution
  |
SCIENTIFIC_STATE         what we consider valid now (if affected)
```

Git history plus this documentation is the decision chain. Reading a diff
should not be required in order to understand what a change means, or why it
was made.

## Status axes

Every documented decision is described on up to three **independent** axes. No
single field mixes two of them. The age of a branch never implies a status —
evidence does.

### Integration status — a fact about Git

| Value | Meaning |
| --- | --- |
| `INTEGRATED` | Merged into `master`. |
| `NOT_INTEGRATED` | Exists only on a staging / historical branch. |
| `WAS_INTEGRATED_REMOVED` | Merged into `master` and later removed. |
| `N/A` | No meaningful Git integration state applies to this record. |

### Lifecycle status — where the decision stands

| Value | Meaning |
| --- | --- |
| `CURRENT` | The operative decision for its topic. |
| `PROPOSED` | Drafted, not yet acted on. |
| `HISTORICAL` | Made and recorded, no longer driving current work, no direct successor. |
| `SUPERSEDED` | Replaced by a specific later decision (linked). |
| `EXPERIMENTAL` | A trial that was not adopted as the operative decision. |
| `DEPRECATED` | Withdrawn, not replaced. |

### Epistemic status — is it scientifically / technically sound

Use it for `SCIENTIFIC` and `TECHNICAL + SCIENTIFIC` decisions, and also for
`TECHNICAL` decisions when there is an evidence-bearing claim about technical
soundness or validation (a passing test suite, a reproducibility check, a
measured baseline). Use `N/A` when no epistemic / evidential judgment is
meaningful.

| Value | Meaning |
| --- | --- |
| `SUPPORTED` | The available evidence supports the decision. |
| `PARTIALLY_SUPPORTED` | Some parts are evidenced; others are proposed or not yet evaluated. |
| `UNRESOLVED` | The evidence does not yet settle it. |
| `N/A` | Not a scientific / methodological decision. |

`master` is the authoritative integration line. `INTEGRATED` is a fact about
Git, not a claim of scientific truth: a decision can be `INTEGRATED` and
`SUPERSEDED`, or `SUPPORTED` and `NOT_INTEGRATED`. `SCIENTIFIC_STATE.md` keeps
the axes apart.

## Branches, and where authority lives

- `master` is the authoritative integration line. Work is accepted onto it only
  through review plus a Pull Request plus a merge.
- `(feat)-siepe-ready-project` is a historical / staging branch. It contains
  real work and useful evidence, but it has **no epistemic or semantic
  authority**: its commit organisation and its roadmap labels (`Phase B`,
  `B1`, `B4.3`, ...) are not the project's final structure. Its relevant
  content is being reintegrated onto `master` as separate, semantically
  described changes — not as one monolithic merge.
- Other branches (`backup/original-messages`, `(perf)-improve-ncbi-fetch-performance`,
  `(fix)-change-canon-id-from-UID-to-accession-version`, ...) are historical or
  experimental and do not orient new work.

## Naming convention

Branches, commits, and PRs use the project's parenthetical-scope pattern and
must describe **what the change does**, not which roadmap step it belongs to.

| Kind | Format | Example |
| --- | --- | --- |
| branch | `(<scope>)-<semantic-kebab-description>` | `(reference)-reproducible-apaz-profile-hmms` |
| commit | `(<scope>) <semantic imperative description>` | `(ontology) separate PIWI-RE family labels from pAgo clade labels` |
| PR | `(<scope>) <Semantic description>` | `(data) Prevent highly similar MID-PIWI references from crossing dataset splits` |

Roadmap identifiers (`Phase B`, `B4.3`, ...) may appear only as secondary
metadata — in a PR body or a ledger entry as
`Historical roadmap reference: B4.3` — never as the primary meaning of a title.

## The staging ontology

`pago_annotation_ontology.md` was drafted during staging work that is not yet
in `master`. Until the corresponding work is reintegrated through a reviewed
PR, no authoritative ontology document exists here; the draft is preserved
under `history/` as provenance. See `PROJECT_LEDGER.md` and the ontology ADR
once the retrospective backfill lands.

## Contributing a decision

1. Open an Issue describing the problem and the evidence for it
   (template: *Relevant change*).
2. Decide the approach: weigh the alternatives, choose one, and state why. If
   this is a decision event, draft the ADR now — before the implementation.
3. Implement on a semantically named branch.
4. Open a PR using the pull request template; state the decision type, whether
   it touches labels / ontology / holdouts, and the reproducibility impact.
5. Include the ADR and a `PROJECT_LEDGER.md` entry in the same PR; update
   `SCIENTIFIC_STATE.md` in the same PR if the change alters what the project
   considers valid.
