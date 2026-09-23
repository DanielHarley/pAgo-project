# Project ledger

A dated list of important decision events. This file answers **"how did we get
here?"** by keeping a compact chronology and pointing to the records that carry
the detail.

For current scientific interpretation, see
[`SCIENTIFIC_STATE.md`](SCIENTIFIC_STATE.md).

This ledger inherits the documentation principles in
[`README.md`](README.md). It does not maintain a parallel classification system
for lifecycle, Git integration, or scientific confidence.

## Conventions

- Entries are chronological.
- Keep each entry as small as the event allows.
- Record the semantic meaning of a change, not only a historical roadmap label.
- Preserve historical commit messages and SHAs; do not rewrite history to make
  old naming cleaner.
- Use `HISTORICAL_EVIDENCE_INSUFFICIENT` or `RATIONALE_NOT_RECOVERABLE` when the
  surviving evidence cannot support a stronger reconstruction.
- Link to ADRs, PRs, commits, audits, or Operation Records only when those links
  materially help a future reader reconstruct the event.
- Do not update an entry merely to mirror ordinary Git state that Git already
  records reliably.

## Entry template

```text
## YYYY-MM-DD — <semantic title>

### Summary
One or two sentences describing what changed and why it mattered.

### Sources
Only the ADRs, PRs, commits, audits, records, or historical metadata that
materially help reconstruct the event.

### Notes
Optional. Include only if a boundary, unresolved question, or historical naming
translation matters.
```

---

# Governance baseline

## 2026-09-09 — Treat `(feat)-siepe-ready-project` as historical/staging evidence

### Summary
The SIEPE branch contains real work and useful evidence, but its roadmap labels
and commit organisation do not define the project's semantic structure. Relevant
work is reviewed and reintegrated as separate, semantically named changes rather
than by a monolithic merge.

### Sources
Owner clarification, 2026-09-09.
`documentation/history/2026-08-30-phase-a-audit.md`.

## 2026-09-09 — Keep agent scratch directories out of the versioned project

### Summary
Local agent scratch/configuration directories and local plan files are not
published as project source. When they contain evidence worth preserving, the
relevant decision or fact is migrated into tracked documentation instead of
versioning the scratch directory itself.

The same change added ignore rules after confirming that nothing under those
paths was already tracked.

### Sources
Owner clarification, 2026-09-09.
Commit: `(chore) ignore local agent scratch directories`.
PR #28.

---

# Decision events

## 2026-02-03 — Repository scaffold

### Summary
Initial project and pipeline scaffold.

### Sources
Commits: `83ba528`, `2c0db28`.

## 2026-02 – 2026-07 — Build NCBI acquisition, snapshots, exploratory analysis, and legacy textual QC

### Summary
The integrated pipeline grew to include NCBI protein retrieval using protein
UIDs, Entrez History pagination, immutable snapshot/manifests, Git-LFS-tracked
raw XML, SWeeP/PCA/KMeans exploratory analysis, metadata-derived QC labels,
custom CA/TLS support, hash-based snapshot invalidation, Windows-safe snapshot
reuse, and NCBI XML failure controls.

The important scientific boundaries are preserved separately: the clustering is
exploratory, and the legacy QC labels are audit features rather than biological
ground truth.

### Sources
PRs #1–#27.
ADRs: [0001](decisions/0001-canonical-ncbi-retrieval-identifier-is-the-protein-uid.md),
[0002](decisions/0002-immutable-snapshot-and-provenance-architecture.md),
[0006](decisions/0006-sweep-pca-kmeans-is-exploratory-not-validated-classification.md).
`documentation/pago_qc.md`.
`data/01-raw/README.md`.
`.github/workflows/ci.yml`.

## 2026-08-15 — Rework NCBI retrieval performance and correctness

### Summary
A measured retrieval rework kept the Entrez History handle for UID paging,
introduced streaming persistence, larger XML batches, selective resume, bounded
concurrency, a `rettype=gp` invariance check, transient handling for truncated
responses, and order-safe `latest/` publication.

### Sources
ADR: [0003](decisions/0003-ncbi-retrieval-performance-and-correctness-rework.md).
Commits: `f5671b8` … `0b4ed5a`.
Historical roadmap reference: pre-Phase-A.

### Notes
The historical implementation overlaps conceptually with later failure controls
from PR #27; that overlap must be understood before reusing the old code.

## 2026-08-30 — Build the annotation-enriched candidate set and query-recall audit

### Summary
A broader NCBI text query produced the
`annotation_enriched_candidate_set`, explicitly a text-annotation-enriched
candidate set rather than a universe of confirmed pAgos. The same work added an
ESearch preflight, a technical-only prefilter, a provenance-preserving derived
FASTA, and an evidence-tiered 21-reference query-recall panel.

Recall matching evolved from accession-only matching to a three-tier strategy:
exact `accession.version`, base accession, then normalized-sequence SHA-256.
PIWI-RE family identity was separated from pAgo clade labels.

### Sources
ADRs: [0004](decisions/0004-annotation-enriched-candidate-set-and-technical-only-prefilter.md),
[0005](decisions/0005-query-recall-reference-panel-and-sequence-identity-matching.md),
[0012](decisions/0012-annotation-ontology-family-versus-clade.md).
Commits: `175a12a` … `0013c6d`.
Historical roadmap reference: Phase A.
Audit: `documentation/history/2026-08-30-phase-a-audit.md`.

### Notes
One audited execution returned 52,473 records. Those execution artifacts were
local rather than versioned. The query-recall panel is a query-coverage
instrument, not a classifier gold standard or final holdout.

## 2026-08-31 — Build the pAgo reference layer and annotation semantics

### Summary
The historical reference-layer work pinned a Pfam 38.2 HMM bundle, defined
high-similarity split groups for leakage control, rebuilt the APAZ reference
partition, built deterministic APAZ HMMs from BUILD references only, reconciled
the Ryazansky 1,010-pAgo catalog, extracted MID-PIWI regions, derived MID-PIWI
high-similarity groups, and formalized the family/clade semantic distinction.

The work also established important boundaries that remain useful independent of
implementation state: the 90/80 split-group rule is a redundancy-control rule,
not a definition of homology; reproducible HMM construction does not establish
predictive performance; and PIWI-RE is a family, not a pAgo clade.

### Sources
ADRs:
[0007](decisions/0007-provenance-and-defensible-labels-before-supervised-ml.md),
[0008](decisions/0008-pinned-pfam-38-2-hmm-bundle.md),
[0009](decisions/0009-apaz-s3-reference-and-high-similarity-split-group-partitioning.md),
[0010](decisions/0010-apaz-profile-hmms-built-from-build-only.md),
[0011](decisions/0011-ryazansky-1010-pago-catalog-and-mid-piwi-reference-extraction.md),
[0012](decisions/0012-annotation-ontology-family-versus-clade.md),
[0013](decisions/0013-mid-piwi-high-similarity-split-groups.md).
Commits: `3de5f1b`, `69ca40e`, `4d8aa99`, `620ff85`, `45cc637`, `9e76e78`,
`c11732a`, `b7490e4`, `481d0c2`, `18c35e4`, `49651f1`.
Historical roadmap reference: Phase B, B0–B4.3.
`documentation/history/2026-08-31-pago-annotation-ontology-staging.md`.

## 2026-09-09 — Establish tracked decision traceability

### Summary
Added the documentation map, Project Ledger, Scientific State, ADR convention,
and GitHub templates so important reasoning and provenance no longer depended on
agent memory, local `tmp/` material, or an untracked private plan.

The first version of this framework also introduced lifecycle, integration, and
epistemic status classifications. Those classifications were later removed when
they proved more costly to maintain than the information they added; the useful
distinctions remain expressed directly in prose.

### Sources
PR #28 — `(docs) Establish a tracked decision-traceability framework`.

## 2026-09-09 — Harden local-artifact ignore rules

### Summary
Changed `.gitignore` so `data/01-raw/` is default-deny except for its tracked
README, and ignored `tmp/`, preventing broad staging from accidentally including
large local acquisition artifacts. Already tracked files were not removed.

### Sources
PR #29.
Commit: `7fc2efe`.

## 2026-09-09 — Reconstruct historical project decisions

### Summary
Preserved the Phase A audit and staging ontology under `history/`, added
retrospective ADRs 0001–0014, expanded the ledger, and aligned current
scientific documentation with the reconstructed record.

The reconstruction used explicit markers such as
`RATIONALE_NOT_RECOVERABLE` and `HISTORICAL_EVIDENCE_INSUFFICIENT` rather than
inventing missing rationale or alternatives. No scientific code, data,
notebooks, HMMs, results, partitions, labels, or holdouts from the staging branch
were introduced by the documentation reconstruction itself.

### Sources
Direct merge commit: `3ab5079`.
PR #30 — post-merge provenance reconciliation.
Branch: `(docs)-reconstruct-historical-project-decisions`.
Reconstruction commits: `442a442`, `7622986`, `7ef7cbf`, `e0284d6`, `6e58e13`.
Provenance-reconciliation commits: `4ed0462`, `e605c37`, `6c8409f8`.

### Notes
The reconstruction reached `master` through a direct `--no-ff` merge rather
than the normal reviewed-PR path. PR #30 documented that exception afterward.

## 2026-09-14 — Protect local reference and validation files from accidental staging

### Summary
Fourteen local, untracked reference-tree and profile/placement-validation files
could have been swept into an unrelated `git add -A`. A verified external
byte-for-byte backup was created, then fourteen exact-path `.gitignore` rules
were added so ordinary broad staging would not select them.

This was a repository-safety operation only. Being ignored does not establish
anything about the scientific validity of those files, and the `.gitignore`
rules are not themselves the backup.

### Sources
PR #31 — `(chore) Protect local reference and validation staging files from
accidental staging`.
Merge commit: `06c52140dfc0b2ce24240df76f7f9acd7cf2ed8c`.
Branch: `(chore)-protect-local-reference-validation-staging-files`.
Operation Record:
[`operations/2026-09-14-local-reference-validation-staging-protection.md`](operations/2026-09-14-local-reference-validation-staging-protection.md).
Historical roadmap reference: B5/B6.

## 2026-09-17 — Reintegrate the Ryazansky pAgo catalog and MID-PIWI extraction

### Summary
The Ryazansky 1,010-pAgo catalog and its MID-PIWI region extraction — the
first scientific unit of the historical `(feat)-siepe-ready-project` branch
reintegrated to `master` through the semantic-reintegration plan — were
rebuilt on current `master`. The catalog's frozen sources, the resolved NCBI
record-status snapshot, the frozen sequence FASTA, the provenance manifest,
and the reconstruction script are now versioned, and the build was
demonstrated to be reproducible offline from those frozen sources. The
historical dependency on the (not-yet-reintegrated) query-recall reference
panel was deliberately removed; the catalog's schema no longer carries the
panel-dependent `experimental_anchor`/`experimental_anchor_name` columns — see
ADR-0015. Remaining reference-layer work, including MID-PIWI high-similarity
grouping, stays separate.

### Sources
ADRs: [0011](decisions/0011-ryazansky-1010-pago-catalog-and-mid-piwi-reference-extraction.md),
[0015](decisions/0015-decouple-the-ryazansky-reference-catalog-from-the-query-recall-panel.md).
PR #35 — `(feat) Reconcile the Ryazansky pAgo catalog and extract MID-PIWI
regions`.
Merge commit: `c6e3562fd4877ba023de24728e174e8b9114b121`.
Historical commits `b7490e4`, `481d0c2`.
Historical roadmap reference: B4.2.

### Notes
`experimental_anchor`, `experimental_anchor_name`, and `proposed_partition`
from the historical catalog schema are not reproduced; see ADR-0015.

## 2026-09-17 — Reintegrate the MID-PIWI high-similarity split groups

### Summary
The MID-PIWI high-similarity split groups — the second scientific unit of
the historical `(feat)-siepe-ready-project` branch reintegrated to `master`
— were rebuilt from the reintegrated Ryazansky catalog's 1,002 MID-PIWI
regions using MMseqs2 (pinned at 18.8cc5c, bioconda) in the declared WSL/
conda environment. The per-region and per-group assignments, the diagnostic-
threshold report, and a provenance manifest are now versioned. The
reconstruction reproduced every historically documented count exactly (701
split groups, 697 eligible, 0 resolved cross-clade groups at 90/80, the
AfAgo/SiAgo manual-review groups, and the 50%-identity diagnostic cross-clade
observation) and was byte-identical across two independent MMseqs2 runs. The
connected-component reducer is a small, deliberate, temporary duplicate of
the historical APAZ reducer, which is not yet reintegrated; no partition
(BUILD/CALIBRATION/FINAL_HOLDOUT), tree, placement, or HMM work was
introduced.

### Sources
ADR-0013.
PR #37 — `(feat) Reintegrate the MID-PIWI high-similarity split groups`.
Historical commits `18c35e4`, `49651f1`.
Historical roadmap reference: B4.3.

### Notes
The connected-component reducer duplicates (rather than imports) the
historical `src/pago_pipeline/apaz_split_groups.py` logic, since that module
is not yet reintegrated; see the module docstring in
`src/pago_pipeline/clade_midpiwi_split_groups.py`. No new ADR was written for
this duplication — it is an integration-ordering/engineering detail, not a
scientific or dataset-construction decision, and ADR-0009/ADR-0013 already
establish that APAZ and MID-PIWI share one high-similarity grouping method.

## 2026-09-23 — Reintegrate the Pfam 38.2 reference bundle

### Summary
The Pfam 38.2 HMM reference bundle — a third, independent scientific unit of
the historical `(feat)-siepe-ready-project` branch, with no code or data
dependency on the APAZ reference lineage — was reinstated to `master`. The ten
version-pinned profiles (PIWI, PAZ, ArgoN, ArgoL1, ArgoL2, ArgoMid, SIR2,
TIR_2, TIR, Mrr_cat) and their SHA-256 lock were ported byte-for-byte from the
historical bundle and verified SHA-256-identical to it. `pyhmmer==0.12.3` was
restored as a pinned runtime dependency, and `scripts/verify_reference_data.py
--scope pfam` validates accession versioning, hash identity, and structural
HMM shape fully offline. No claim of predictive sensitivity, specificity, or
biological validity is made; the APAZ reference set and profile HMMs remain
not yet reintegrated.

### Sources
ADR-0008.
Historical commits `3de5f1b`, `69ca40e`.
Historical roadmap reference: B0 (partial: the `pyhmmer` pin only), B1.

### Notes
Historical commit `3de5f1b` (B0) also touched unrelated `.gitignore` scratch-
directory rules and other housekeeping not needed by any module currently on
`master`; only its `pyhmmer==0.12.3` dependency pin was reinstated, since
that is the only part of B0 with a real consumer in this reintegration.
`scripts/verify_reference_data.py` was rewritten rather than ported
byte-for-byte: the historical file had grown speculative verifiers for
`apaz`/`clade`/`tree` scopes whose backing modules did not exist even at the
commit that added them; only the `pfam` scope is implemented here.

## 2026-09-23 — Reintegrate the APAZ reference set and partitions

### Summary
The APAZ reference-curation layer — historical commits `4d8aa99`, `620ff85`,
and `45cc637`, reintegrated as one coherent unit because they share a single
lock/schema (`apaz_partitions.csv`, `seeds_lock.json`) — was reinstated to
`master`. The frozen Ryazansky Data Set S3 source (481 representatives),
the HisG/EIIB Pfam hard-negative seed alignments, the canonical 90%
identity/80% coverage split-group reducer, the 460/509/517 frozen split
groups, and the deterministic FINAL_HOLDOUT-then-CALIBRATION-then-BUILD
partition (positives 337/72/72; HisG 0/255/254; EIIB 0/259/258) were ported
byte-for-byte and verified SHA-256-identical to the historical freeze. Two
independent offline regenerations from the frozen sources (including the
PyFAMSA 0.7.0 global BUILD alignment) reproduced every derived artifact
byte-for-byte, both against each other and against the committed resources.
Whole 90/80 groups never cross a partition, and HisG/EIIB never enter BUILD.
No APAZ profile HMM was built and no predictive claim is made; that
construction step (B3) remains separate and not yet reintegrated.

### Sources
ADR-0009.
Historical commits `4d8aa99`, `620ff85`, `45cc637`.
Historical roadmap reference: B2.

### Notes
Independent regeneration of the MMseqs2 90/80 split-group edges themselves
(`scripts/derive_apaz_split_groups.py`, which needs MMseqs2 18.8cc5c in the
declared WSL/conda environment) was not executed in this reintegration: no
MMseqs2 installation was available in the session's environment. This does
not affect the acceptance evidence, since the committed `split_groups/`
resources are verified byte-identical to the historical commit and the
frozen edge list is independently proven, by the ported test suite, to
reproduce the committed `split_groups.csv` through the pure-Python reducer.
The full historical procedure that enumerated every one of the eight
reported v1 cross-partition leaks is not recoverable beyond the one worked
example already in `curation_notes.md`; this is treated as an informational
historical limitation, not a reintegration blocker, since v1 was never used
for predictive evaluation and v2 is independently reproducible.

## 2026-09-23 — Add deterministic APAZ profile HMM construction

### Summary
The repository gained the ability to deterministically construct, structurally
validate, and snapshot the six APAZ profile HMMs (global, Ia, Ib, IIa, IIb,
III) exclusively from the frozen BUILD seed alignments, under the pinned
`pyhmmer==0.12.3` / `pyhmmer.plan7.Builder.build_msa` contract. Every written
model is reopened from disk and structurally proven (single model, amino
alphabet, expected name) before its SHA-256 is recorded, and CALIBRATION/
FINAL_HOLDOUT accessions are rejected if found in a seed alignment. Two
independent real builds from the committed BUILD seeds produced
byte-identical output for all six models. The generated `.hmm` files are
derived build artifacts, not versioned source-of-truth resources. No
calibration, threshold selection, or final-holdout evaluation was performed.

### Sources
ADR-0010.
PR — `(feat) Add deterministic APAZ profile HMM construction`.
Historical commits `9e76e78`, `c11732a` (secondary provenance only — this
capability is documented and named by what it currently adds to `master`,
not by its historical roadmap step).

### Notes
`src/pago_pipeline/apaz_hmm_build.py` and `src/pago_pipeline/apaz_hmm_build_snapshot.py`
ported from that historical evidence with no code adaptation required: both
already matched the current APAZ reference schema (`apaz_partitions.csv`,
`seeds_lock.json` format 2.0) and the current shared snapshot-directory
helpers in `src/pago_pipeline/ncbi_snapshot.py`.

---

## Local bibliographic material

Local research notes under `tmp/` are bibliographic leads produced during
report preparation. They are not scientific authority. Claims that depend on
those leads must be checked against the primary paper or an appropriate official
source before being promoted into scientific documentation.

## ADR index

See [`decisions/README.md#index`](decisions/README.md#index).
