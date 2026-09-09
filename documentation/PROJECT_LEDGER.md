# Project ledger

A dated list of decision events — the moments the project changed its
behaviour, interpretation, architecture, methodology, or a scientific rule.

Each entry is short and points to more detail (an ADR, a PR, commits, an audit).
This file answers **"how did we get here?"**. For **"what is valid now?"** see
[`SCIENTIFIC_STATE.md`](SCIENTIFIC_STATE.md).

## Conventions

- Entries are chronological.
- Each entry states, on **independent** axes: **lifecycle status**
  (`CURRENT` / `PROPOSED` / `HISTORICAL` / `SUPERSEDED` / `EXPERIMENTAL` /
  `DEPRECATED`), **integration status** (`INTEGRATED` / `NOT_INTEGRATED` /
  `WAS_INTEGRATED_REMOVED` / `N/A`), and, for scientific decisions,
  **epistemic status** (`SUPPORTED` / `PARTIALLY_SUPPORTED` / `UNRESOLVED` /
  `N/A`). No field mixes two axes. See [`README.md`](README.md#status-axes).
- Historical Git commits keep their original SHAs and messages. Where a commit
  message names only a roadmap step, the entry adds a **semantic
  interpretation** and keeps the roadmap label as
  `Historical roadmap reference:`.
- Where the historical record cannot support a field, the entry says so
  (`HISTORICAL_EVIDENCE_INSUFFICIENT`, `RATIONALE_NOT_RECOVERABLE`) instead of
  guessing.
- An entry whose integration status says `INTEGRATED (on merge of this PR)`
  describes the state that becomes true when the Pull Request introducing this
  framework is merged into `master`.

## Entry template

```
## YYYY-MM-DD — <semantic title>

### Type
TECHNICAL | SCIENTIFIC | TECHNICAL + SCIENTIFIC

### Lifecycle status
CURRENT | PROPOSED | HISTORICAL | SUPERSEDED | EXPERIMENTAL | DEPRECATED

### Integration status
INTEGRATED | NOT_INTEGRATED | WAS_INTEGRATED_REMOVED | N/A

### Epistemic status
SUPPORTED | PARTIALLY_SUPPORTED | UNRESOLVED | N/A
(for SCIENTIFIC / TECHNICAL + SCIENTIFIC decisions, and TECHNICAL decisions with
an evidence-bearing soundness or validation claim; N/A when no such judgment
applies)

### Summary
One or two sentences: problem -> decision.

### Semantic interpretation
(only when the historical commit/PR is labelled by roadmap step)

### Sources
ADR:
Issue:
PR:
Commits:
Historical roadmap reference:
Audit / notes:

### Supersedes
```

---

# Governance baseline (2026-09-09)

## 2026-09-09 — `master` is the authoritative integration line

### Type
TECHNICAL (process)

### Lifecycle status
CURRENT

### Integration status
INTEGRATED (on merge of this PR)

### Epistemic status
N/A

### Summary
`master` is where accepted changes are integrated, through review + PR + merge.
Being in `master` is a fact about integration, not a claim of scientific truth.
The documentation describes decisions on three independent axes — lifecycle
status, integration status, and epistemic status — and never merges two into
one field. A branch's age never implies any of them; evidence does.

### Sources
ADR: (governance ADR, retrospective backfill PR)
PR: this PR — `(docs) Establish a tracked decision-traceability framework`

## 2026-09-09 — `(feat)-siepe-ready-project` is a historical / staging branch

### Type
TECHNICAL (process)

### Lifecycle status
HISTORICAL (as a branch)

### Integration status
NOT_INTEGRATED (its content)

### Epistemic status
N/A (per-item epistemic status is in `SCIENTIFIC_STATE.md` §2)

### Summary
The branch was created to prepare project state for the SIEPE report. It holds
real work and useful evidence but has **no epistemic or semantic authority**:
its commit organisation and its roadmap labels (`Phase A/B`, `B1`, `B4.3`, ...)
are not the project's final structure. Its relevant content will be reintegrated
onto `master` as separate, semantically described PRs — **not** as a single
monolithic merge. Per change: audit -> identify the problem -> semantic Issue ->
semantic branch -> semantic commits -> semantic PR -> validation -> merge, with
cherry-pick / reimplementation / refactor decided case by case.

### Sources
Owner clarification, 2026-09-09.
Retrospective audit (added under `documentation/history/` by the backfill PR).

## 2026-09-09 — Agent scratch directories stay out of the versioned project

### Type
TECHNICAL (process)

### Lifecycle status
CURRENT

### Integration status
INTEGRATED (on merge of this PR)

### Epistemic status
N/A

### Summary
Local agent scratch and configuration directories, and local agent plan files,
are not part of the scientific or technical source that gets published. They
may be used as historical evidence during reconstruction, but a confirmed
decision is migrated into tracked documentation (with a confidence level)
rather than the directory being versioned. If a rationale exists only there,
extract the decision; do not commit the directory.

As of `master` `9d1e55a` these directories were not all listed in `.gitignore`.
This PR adds ignore rules for local agent scratch and configuration directories
in an isolated commit so the policy is enforced, not merely stated. Nothing was
already tracked under those paths (`git ls-files` returned empty).

### Sources
Owner clarification, 2026-09-09.
Commit: `(chore) ignore local agent scratch directories`

## 2026-09-09 — Naming convention: semantics first, roadmap IDs as metadata

### Type
TECHNICAL (process)

### Lifecycle status
CURRENT

### Integration status
INTEGRATED (on merge of this PR)

### Epistemic status
N/A

### Summary
Branches `(<scope>)-<semantic-kebab-description>`, commits
`(<scope>) <semantic imperative description>`, PRs
`(<scope>) <Semantic description>`. Git must explain **what** a change does.
Roadmap identifiers (`Phase B`, `B4.3`) appear only as secondary metadata in a
PR body or a ledger entry (`Historical roadmap reference: B4.3`), never as the
primary meaning of a title. Published history is **not** rewritten; the
translation from roadmap label to meaning is done here, in documentation.

### Sources
Owner clarification, 2026-09-09; existing project branch naming.

---

# Decision events

> The detailed reconstruction of the events below is added by the follow-up PR
> `(docs) Reconstruct the rationale and provenance of major historical decisions`.
> The retrospective reconstruction adds ADRs as needed, one per significant
> conceptual decision — a broad historical ledger entry may map to multiple
> ADRs. The entries here establish the timeline and the status axes.

## 2026-02-03 — Repository scaffold

### Type
TECHNICAL

### Lifecycle status
HISTORICAL

### Integration status
INTEGRATED

### Epistemic status
N/A

### Summary
Initial pipeline scaffold.

### Sources
Commits: `83ba528`, `2c0db28`

## 2026-02 – 2026-07 — NCBI acquisition, snapshot system, exploratory analysis, textual QC, snapshot robustness

### Type
TECHNICAL + SCIENTIFIC

### Lifecycle status
CURRENT (the integrated pipeline; the narrow search query within it is
`SUPERSEDED` — see `SCIENTIFIC_STATE.md` §3)

### Integration status
INTEGRATED

### Epistemic status
PARTIALLY_SUPPORTED — acquisition and snapshot behaviour are verified; the
exploratory SWeeP / PCA / KMeans analysis is exploratory only and the textual
QC labels are audit features, not ground truth

### Summary
The integrated pipeline was built over this period: NCBI protein retrieval
(canonical identifier = protein UID; ESearch via Entrez History), the immutable
snapshot + manifest + `reuse_latest_or_create` system, Git-LFS-tracked raw XML,
SWeeP / PCA / KMeans exploratory analysis, the metadata-derived textual QC
labelling layer (`excluded` is not a biological negative; PIWI-RE handled as a
separate class; labels are audit features), custom CA / TLS support, hash-based
`latest/` invalidation, Windows `WinError 32` avoidance in snapshot reuse, and
NCBI XML failure controls.

### Sources
PRs: #1–#27 (merged to `master`).
ADRs: to be added by the backfill PR.
Audit / notes: `documentation/pago_qc.md`, `data/01-raw/README.md`,
`.github/workflows/ci.yml`.

## 2026-08-15 — NCBI retrieval performance and correctness rework

### Type
TECHNICAL

### Lifecycle status
CURRENT

### Integration status
NOT_INTEGRATED

### Epistemic status
SUPPORTED — technically validated by the available evidence

### Summary
Measured baseline plus a phased rework (findings F1–F12, phases P0–P5): keep the
History handle for UID retrieval, streaming persistence, larger XML batches,
selective resume, bounded concurrency (default off), an `rettype=gp` invariance
test, truncated responses reclassified as transient, order-safe `latest/`
publish.

### Semantic interpretation
Retrieve NCBI records via History paging with bounded concurrency and resumable
batches, and prove the payload contract is stable.

### Sources
Commits: `f5671b8` … `0b4ed5a`
PR: No historical PR found (committed directly to the staging branch)
Historical roadmap reference: pre-Phase-A
Audit / notes: `documentation/ncbi_retrieval_performance_plan.md`,
`documentation/ncbi_retrieval_performance_implementation.md`
(staging source; preserved / integrated by the retrospective backfill)

### Open at reintegration
Relationship to the PR #27 failure controls already in `master`
(overlap vs supersession).

## 2026-08-30 — Annotation-enriched pAgo candidate set, recall panel, technical prefilter, derived FASTA

### Type
TECHNICAL + SCIENTIFIC

### Lifecycle status
CURRENT

### Integration status
NOT_INTEGRATED

### Epistemic status
SUPPORTED — sound design, audited execution; the candidate set is explicitly
not a pAgo universe and the recall panel is not a validation gold standard

### Summary
Second NCBI acquisition with a broad query
`(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])`,
named `annotation_enriched_candidate_set` and explicitly **not** the pAgo
universe. A strictly technical prefilter (never filters by annotation text or
length). A 21-reference query-recall panel with per-row evidence tiers, matched
by accession then by normalized-sequence identity, reported as two readings; the
panel is not a global gold standard and not a validation holdout. A
provenance-honest derived FASTA. `ago_family` (PIWI_RE) separated from
`pago_clade`.

### Semantic interpretation
Build an annotation-enriched pAgo candidate set with a stratified,
evidence-tiered query-recall report; recognise references retrieved under a
different accession by identical sequence.

### Sources
Commits: `175a12a` … `0013c6d`
PR: No historical PR found
Historical roadmap reference: Phase A
Audit / notes: `documentation/history/2026-08-30-phase-a-audit.md`
(added by the backfill PR)

### Supersedes
Exact-accession-only recall matching; `PIWI_RE` as a `pago_clade` value;
(pending reintegration) the narrow `PIWI[All Fields] AND Bacteria[Organism]`
query.

## 2026-08-31 — pAgo reference layer and annotation ontology

### Type
TECHNICAL + SCIENTIFIC

### Lifecycle status
CURRENT

### Integration status
NOT_INTEGRATED

### Epistemic status
PARTIALLY_SUPPORTED — resources, construction, and integrity checks are in
place; predictive performance is not yet evaluated, and the ontology's future
phylogenetic-placement protocol is `PROPOSED` and not executed

### Summary
Pinned Pfam 38.2 HMM bundle (versioned accessions, SHA-256 lock). A
high-similarity `split_group` unit (MMseqs2 90 % identity / 80 % coverage
connected components) as the indivisible train/test unit — explicitly **not** a
definition of homology. APAZ partitioning `apaz_partition_v2_mmseqs90_80` with a
frozen FINAL_HOLDOUT, superseding a v1 in which high-similarity split groups
could cross partition boundaries; HisG / EIIB hard negatives kept out of BUILD.
Six APAZ profile HMMs
built only from BUILD, bit-for-bit reproducible, with generated `.hmm` files
gitignored. A deterministic catalog of Ryazansky Table S1 (1010 proteins,
pre-reduction) with proven MID-PIWI coordinates, NCBI record status, and
quarantine of conflicting references (AfAgo, SiAgo). MID-PIWI high-similarity
split groups with 0 cross-clade groups at 90/80 and a 50 %-identity cross-clade
structure recorded as diagnostic only. A formal annotation ontology.

### Semantic interpretation
- `(reference) Pin a versioned Pfam 38.2 HMM bundle with SHA-256 locks`
- `(data) Keep high-similarity APAZ groups within a single evaluation partition`
- `(reference) Build reproducible APAZ profile HMMs from training references only`
- `(reference) Reconcile Ryazansky pAgo accessions and quarantine conflicting clade labels`
- `(data) Prevent highly similar MID-PIWI references from crossing future dataset splits`
- `(ontology) Separate PIWI-RE family identity from pAgo phylogenetic clades`

### Sources
Commits: `3de5f1b` (B0), `69ca40e` (B1), `4d8aa99` `620ff85` `45cc637` (B2),
`9e76e78` `c11732a` (B3), `b7490e4` `481d0c2` (B4.2), `18c35e4` `49651f1` (B4.3)
PR: No historical PR found
Historical roadmap reference: Phase B, milestones B0–B4.3
Audit / notes: `src/pago_pipeline/resources/apaz_seed/curation_notes.md`,
`src/pago_pipeline/resources/clade_seed/ryazansky_s1_catalog_notes.md`,
`documentation/history/2026-08-31-pago-annotation-ontology-staging.md`
(added by the backfill PR — preserved verbatim, with a header terminology note
that the staging text uses "homology cluster" for the operational 90/80
grouping now called a "high-similarity split group"; the historical content is
not rewritten)

### Supersedes
APAZ reference partitioning v1 (accession-hash split unit; high-similarity split
groups could cross BUILD / CALIBRATION boundaries).

## 2026-09-09 — Decision-traceability framework established

### Type
TECHNICAL (process)

### Lifecycle status
CURRENT

### Integration status
INTEGRATED (on merge of this PR)

### Epistemic status
N/A

### Summary
Added `documentation/README.md`, this ledger, `SCIENTIFIC_STATE.md`,
`decisions/` with an ADR convention and template, and PR / Issue templates, so
that the rationale and provenance of decisions live in the repository rather
than in `tmp/`, agent memory, or a private plan. A follow-up commit separates
the lifecycle, integration, and epistemic axes into independent fields.

### Sources
PR: this PR — `(docs) Establish a tracked decision-traceability framework`
