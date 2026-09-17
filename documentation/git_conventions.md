# Git and GitHub conventions

This document defines how project changes are named and described in Git and
GitHub. It applies to human contributors and to automated coding agents
working on this repository.

## General principles

Infer the real change, problem, limitation, or required repository behavior
from the diff, staged changes, branch diff, filenames, symbols, surrounding
code context, inspected artifacts, logs, and explicit instructions. Do not
invent changes, motivations, validations, side effects, risk statements, root
causes, reproduction results, or test-execution claims that are not directly
supported by the available evidence.

Write every branch name, commit message, PR description, and Issue in
English, regardless of the language used to request it.

Never include an AI/tool attribution footer (for example
`Co-authored-by: <tool>`, `Generated with <tool>`, or `Made with <tool>`)
unless explicitly requested.

Preserve semantic distinctions material to truthful writing, for example:
protein UID vs accession identity, structural validation vs
semantic/scientific validation, snapshot creation vs snapshot reuse, and a
planned validation step vs an executed test.

Be precise and specific. Use imperative mood in commit titles. Name the exact
mechanism instead of vague claims such as "improves reliability" or "adds
robustness". Avoid vague wording such as "several updates" or "cleanup"
unless the change genuinely is broad maintenance work.

## Truthfulness constraints

- Never claim tests ran unless the diff or the request proves that.
- Never claim a fix, refactor, validation, hardening, or performance
  improvement unless the diff clearly supports it.
- Never describe a validation plan as an executed test result.
- Never infer successful end-to-end execution from code changes alone.
- Never invent filenames, artifacts, helpers, or functions not present in the
  diff or the request.
- Never claim an issue's root cause is identified unless code, artifacts,
  logs, or an explicit instruction directly establish it.
- Never claim an issue reproduces consistently unless repeated reproduction
  is evidenced.

## Naming convention

Branches, commits, and PRs use the project's parenthetical-scope pattern and
must describe what the change does.

| Kind | Format | Example |
| --- | --- | --- |
| branch | `(<scope>)-<semantic-kebab-description>` | `(reference)-pin-pfam-38-2-hmm-bundle` |
| commit | `(<scope>) <semantic imperative description>` | `(docs) simplify documentation governance` |
| PR | `(<scope>) <Semantic description>` | `(docs) Simplify documentation governance` |

A title such as `Phase B`, `B5 files`, or `misc changes` is not sufficiently
semantic. Roadmap identifiers such as `Phase A`, `B1`, or `B4.3` may appear
only as secondary provenance metadata inside a record, never as the primary
identity of a branch, commit, or PR.

## Scopes

Choose the scope by the primary semantic effect of the change, not by file
count, code size, or how much work it took.

Decision order:

1. If the change only affects documentation, use `docs`.
2. If the change only affects tests or test fixtures, use `test`.
3. If the change only affects CI/CD workflow or automation config, use `ci`.
4. If the change only affects packaging, dependency management, build
   tooling, or environment/build configuration, use `build`.
5. If the main purpose is measurable runtime or memory improvement, use
   `perf`.
6. If the change corrects incorrect behavior, missing validation that should
   already have existed, unsafe acceptance of invalid state, wrong output,
   broken logic, or a real defect, use `fix`.
7. If the change introduces a new capability, new workflow, new persisted
   artifact, new validation surface, or materially expands what the system
   can do, use `feat`.
8. If the change restructures existing code without intended behavioral
   change, use `refactor`.
9. If the change is internal maintenance, hygiene, consistency work,
   non-user-facing operational tightening, ignore rules, artifact cleanup, or
   validation tightening that is not best described as a defect correction or
   new feature, use `chore`.

Tie-break rules:

- Prefer `fix` over `chore` when the previous behavior was actually incorrect
  or unsafe.
- Prefer `feat` over `fix` when the branch adds a genuinely new capability
  rather than correcting an existing one.
- Prefer `refactor` only if behavior is intended to stay the same.
- Prefer `chore` over `refactor` when the change is mostly
  maintenance/hygiene rather than structural redesign.
- If a change mixes multiple categories, choose the scope that best matches
  the main reviewer-visible outcome, not a small secondary change such as a
  `.gitignore` edit or formatting cleanup.

For this repository specifically: use `feat` for a new snapshot type,
resolver, materialization flow, persistence path, or end-to-end pipeline
capability; use `fix` for correcting snapshot validation semantics, artifact
integrity checks, reuse correctness, UID matching logic, XML parsing
correctness, or publication safety for already existing behavior; use `chore`
for repository hygiene, ignore rules, or non-user-facing maintenance; use
`refactor` only when restructuring snapshot or API code without intended
semantic change.

## Commit messages

A commit message has exactly these parts, in order:

1. A one-line title: `(<scope>) <concise imperative summary>`
2. A blank line
3. One short paragraph stating what the change does at a high level
4. A blank line
5. A section titled exactly `Changes:` followed by a bullet list of concrete
   changes
6. A blank line
7. A section titled exactly `Rationale:` followed by one short paragraph
   explaining why the change was necessary

Do not add extra sections. Do not add markdown headings. Do not wrap the
message in code fences unless explicitly asked to.

Example:

```
(chore) verify XML record UIDs against saved protein UID snapshots

Ensure consolidated XML snapshots are validated not only by structure and
counts but also by the exact protein UID set encoded in the XML records,
both when saving new snapshots and when reusing existing ones.

Changes:
- add a helper to recover a single protein UID from each `GBSeq` record via
  its `GBSeqid` fields
- add a validator that parses the consolidated XML file and compares
  extracted record UIDs against the expected protein UID list
- reject XML snapshots when a record is missing a UID, exposes multiple
  conflicting UID candidates, or contains unexpected XML child tags
- run the same validation during snapshot creation and snapshot reuse

Rationale:
Matching record counts alone is not enough to prove that a saved XML
snapshot corresponds to the intended protein UID snapshot. Validating the
exact UID set inside the XML records prevents mismatched, duplicated, or
semantically corrupted XML payloads from being silently accepted by the
pipeline.
```

## Pull requests

Do not add a PR title unless explicitly asked for it; when asked, use the PR
naming convention above.

A PR description uses exactly this section order and heading style, with no
other sections:

- `## Summary` — 1–2 short paragraphs: what the branch adds or changes at a
  high level, then its architectural, workflow, validation, persistence, or
  behavioral impact.
- `## What changed` — a bullet list of concrete changes derived from the
  diff; nested bullets only when needed for grouped artifacts or checks.
- `## Why` — one short paragraph explaining the underlying problem,
  limitation, or gap before this branch, causally (not merely descriptively),
  optionally followed by a short bullet list of what the branch now
  guarantees or prevents.
- `## Notable implementation details` — a bullet list of implementation
  details, invariants, stricter validation rules, reuse semantics,
  publication guarantees, or failure-handling behavior reviewers should
  understand.
- `## Files changed` — a bullet list of the most relevant changed files
  (repository-relative paths), omitting files not materially relevant to
  review.
- `## Test plan` — a bullet list that distinguishes executed verification
  supported by evidence from an intended, not-yet-run validation step;
  describe the latter as a plan, never as completed verification.
- `## Risks / review focus` — a bullet list of the most important real
  review risks, assumptions, sharp edges, or compatibility implications
  suggested by the diff; do not invent generic risks.

Do not wrap the PR description in code fences unless explicitly asked to.

This document is the structural authority for PR descriptions.
`.github/pull_request_template.md` must reflect exactly the seven sections
above and must not introduce additional `##` headings or a competing
structure. The documentation impact assessment required by
`documentation/README.md` remains mandatory for a substantive PR; its result
belongs inside one of the seven sections above (typically `Notable
implementation details`), not in a section of its own.

## Issues

An Issue body uses exactly this section order and heading style, with no
other sections, and no labels/assignees/milestones/priority unless
explicitly asked for:

- `## Problem` — one short paragraph describing the observed problem,
  missing capability, limitation, or repository requirement, distinguishing
  an observed incorrect behavior from a requested new capability. Do not
  assert a root cause unless directly evidenced.
- `## Evidence / current behavior` — a bullet list of concrete, observed or
  explicitly supplied evidence (affected module/file/artifact, error
  messages, reproduction context); if evidence is incomplete, state exactly
  what is known without filling gaps.
- `## Expected behavior` — a bullet list stating the correct behavior as
  explicit, verifiable invariants or workflow guarantees, not vague
  expectations such as "should work correctly".
- `## Scope` — a bullet list of the smallest review-relevant implementation
  surface justified by the evidence; do not turn the issue into a
  speculative implementation design.
- `## Acceptance criteria` — a `- [ ]` checklist of concrete, reviewable,
  testable criteria; do not claim any criterion is already satisfied unless
  proven.
- `## Validation plan` — a bullet list of concrete steps that could verify
  the resolution, written as planned verification unless execution is
  explicitly evidenced; include negative/rejection cases when the issue
  concerns invalid state.
- `## Risks / review focus` — a bullet list of real assumptions, semantic
  ambiguities, or compatibility implications suggested by the issue; do not
  add generic risks detached from it.

Do not describe a proposed fix as already implemented. Do not wrap the issue
body in code fences unless explicitly asked to.
