# Scientific state

What the project currently considers valid, and where each item stands relative
to the `master` integration line.

This document separates two axes that must not be conflated:

- **epistemic status** — is the result / decision scientifically sound?
- **Git integration status** — is it merged into `master`?

"It is in `master`" does not mean "it is scientifically true". "It is on the
staging branch" does not mean "it is authoritative". See
[`README.md`](README.md#integration-states).

Reference commits: `master` at `9d1e55a`; staging branch
`(feat)-siepe-ready-project` at `49651f1` (as of 2026-09-09).

Detailed rationale and provenance for the items below live in
[`decisions/`](decisions/) and [`PROJECT_LEDGER.md`](PROJECT_LEDGER.md); the
retrospective ADRs are added by a follow-up backfill PR.

---

## 1. Integrated authoritative state

Integrated in `master`. `INTEGRATED`.

### 1.1 Data acquisition

- **Canonical NCBI record identifier: protein UID (GI).** Retrieval and XML
  validation operate on protein UIDs (PR #5). A later branch named for reverting
  this to `accession.version` (`(fix)-change-canon-id-from-UID-to-accession-version`)
  never implemented or integrated a reversal; UID remains authoritative.
- **ESearch pagination via Entrez History (WebEnv / QueryKey)** for
  reproducibility (PR #3).
- **Search query (integrated):** `PIWI[All Fields] AND Bacteria[Organism]`.
  This is the narrow, earlier query. It has been superseded on the staging
  branch by a broader annotation-enriched query (see §2.2); that replacement is
  not yet in `master`.
- **Consolidated XML snapshot is the reproducibility source of record.**
  Re-fetching the same UIDs from NCBI is not guaranteed to reproduce the inputs
  because NCBI records change. Raw XML is tracked with Git LFS; snapshot runs
  under `snapshots/` are published one at a time together with the `latest/`
  pointer that names them.
- **Custom CA / TLS configuration** is supported without disabling certificate
  verification (PR #13).
- **NCBI failure controls** — retry / backoff, circuit breaker, deadlines,
  resumable XML batch workspace, batch audit cache (PR #27).

### 1.2 Snapshot and provenance system

- Immutable `snapshots/<timestamp>__<hash>/` directories plus a replaceable
  `latest/`, each with a `manifest.json` recording SHA-256 hashes of outputs
  and of parent snapshots.
- `reuse_latest_or_create`: a `latest/` snapshot is reused only when every
  input hash it declares still matches; upstream input changes invalidate it
  (PR #20).
- Manifest-only validation happens before any `.npy` array is loaded, and a
  compatible immutable snapshot is republished as `latest/` before the old
  `latest/` is touched — this avoids Windows `PermissionError` / `WinError 32`
  from memory-mapped `.npy` files (PR #21).
- Snapshot reuse requires declared artifact hashes and, for SWeeP, a matching
  package-root identity (PR #21).
- CI runs the snapshot-reuse regressions on `windows-latest` as a required job.

### 1.3 Exploratory analysis (integrated, exploratory only)

- SWeeP Genes embeddings (SWeeP 2.1.3.0 is **not** committed — intellectual
  property restriction; reproducibility of this step is conditional), PCA,
  KMeans, and a 3-D PCA plot (PRs #9–#11).

### 1.4 Legacy textual QC layer

Implemented by `pago_qc*` and documented in [`pago_qc.md`](pago_qc.md).

- Three separated contracts: `evidence_flags` -> `labelled_records` ->
  `filtered_datasets`.
- **`excluded_records.csv` is not a biologically validated negative class** —
  it means "excluded from the conservative classic-pAgo positive set".
- **PIWI-RE is handled as a separate class** (`separate_dataset`), not as a
  doubtful classic pAgo.
- **QC labels are deterministic audit features, not ground truth.** They are
  not equivalent to HMMER / Pfam / CDD / InterProScan validation.
- The filtered-dataset counts in `pago_qc.md`
  (`classic_pago_high_precision` 2,629 / `classic_pago_review` 6,703 /
  `piwi_re` 15,056 / `excluded` 16,957 / total **41,345**) describe the run of
  the **earlier narrow query** only. They are not an estimate of the number of
  pAgos and are not the staging-branch 52,473-record candidate set (§2.2).

### 1.5 Notebooks

`notebooks/00`–`09` (acquisition -> metadata -> QC -> FASTA -> SWeeP -> PCA ->
KMeans -> 3-D plot -> QC evidence -> QC filtered datasets).

---

## 2. Validated work awaiting integration

Present only on `(feat)-siepe-ready-project`. Each item has strong evidence but
has **not** been reintegrated onto `master` through a reviewed PR. The two axes
are independent: the work is scientifically / technically supported by the
evidence on the staging branch, and it is separately not yet part of the
integrated project state.

| Field | Value |
| --- | --- |
| Epistemic status | `VALIDATED` — supported by the evidence available on the staging branch |
| Integration status | `NOT INTEGRATED INTO MASTER` |
| Source branch | `(feat)-siepe-ready-project` |
| Reintegration | separate, semantically described PRs — **not** one monolithic merge |

### 2.1 NCBI retrieval performance and correctness rework

- **Evidence:** [`ncbi_retrieval_performance_plan.md`](ncbi_retrieval_performance_plan.md)
  and [`ncbi_retrieval_performance_implementation.md`](ncbi_retrieval_performance_implementation.md)
  (staging copies), measured baselines, `rettype=gp` invariance test, byte-identical
  consolidation at full scale.
- **Content:** UID retrieval keeps the History handle, streaming persistence,
  larger XML batches, selective resume, bounded concurrency (default off),
  truncated responses reclassified as transient, order-safe `latest/` publish.
- **Historical commits:** `f5671b8` … `0b4ed5a`.
- **Historical roadmap reference:** pre-Phase-A.
- **Open at reintegration:** relationship to the PR #27 failure controls already
  in `master` must be reconciled (overlap vs supersession).

### 2.2 Annotation-enriched pAgo candidate set

- **Evidence:** `documentation/history/2026-08-30-phase-a-audit.md` (added by the
  backfill PR), notebook 10, execution manifests on the author's disk (not
  versioned).
- **Content:** query
  `(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])`;
  dataset named `annotation_enriched_candidate_set`, explicitly **not** the pAgo
  universe (a text query only recovers proteins already annotated with the
  terminology). One execution (2026-08-30/31) returned **52,473** protein
  records. The 52,473-record execution artifacts are regenerable but **not
  versioned**.
- **Historical commits:** `175a12a` … `0013c6d`.
- **Historical roadmap reference:** Phase A.

### 2.3 Technical prefilter

- **Content:** excludes only technically unusable records (missing / invalid
  sequence, missing `protein_uid`, technical duplicates). **Never** excludes by
  annotation text or by sequence length; length outside a band sets
  `length_warning` and the record is kept. On the 52,473-record run: 0
  exclusions, 8,299 `length_warning`.
- **Rationale:** the pAgo / non-pAgo separation is a domain / HMM task for a
  later phase, decided on structural evidence.

### 2.4 Query-recall reference panel

- **Content:** 21 curated references (14 pAgo across LONG_A / LONG_B / SHORT,
  7 PIWI-RE), each with a `reference_label_evidence` tier
  (`EXPERIMENTAL` / `LITERATURE_PHYLOGENETIC` / `CURATED_COMPUTATIONAL` /
  `DATABASE_ANNOTATION`). Matching hierarchy: exact `accession.version` ->
  base accession -> normalized-sequence SHA-256. Two recall readings are kept
  (`exact_accession_recall`, `retrieval_equivalent_recall`).
- **This panel is not a global gold standard and not a validation holdout.**
  The 6 `CURATED_COMPUTATIONAL` PIWI-RE members were selected *using* PIWI-RE
  profile models — scoring a profile detector against them would be circular.
- **Evidence:** `tests/fixtures/query_recall_reference_set_curation_notes.md`
  (staging).

### 2.5 Sequence-identity equivalence in recall matching

- **Content:** a third matching tier, `SEQUENCE_SHA256` (normalized
  whitespace-stripped uppercased sequence hash), offline and deterministic,
  pinned by `matching_strategy_sha256`. Motivated by RsAgo (`ABP72561.1`),
  which is retrieved under the byte-identical IPG alias `A4WYU7.1`. The
  reference fixture accession was **not** changed.
- **Supersedes:** exact-accession-only matching (§3).

### 2.6 pAgo reference layer (HMM / phylogenetics instruments)

- **Pfam 38.2 bundle** — 10 version-pinned HMMs (PIWI, PAZ, ArgoN, ArgoL1,
  ArgoL2, ArgoMid, SIR2, TIR_2, TIR, Mrr_cat), SHA-256 lock, unversioned Pfam
  accessions rejected. Historical roadmap reference: B1.
- **`split_group` as the indivisible train/test unit** — a connected component
  of a MMseqs2 graph with edges at >= 90 % identity **and** >= 80 % coverage,
  identity and coverage recomputed from the aligned strings (never MMseqs
  `pident` / `alnlen`), with a prior union of exact-sequence duplicates. **This
  is a high-similarity / redundancy control unit. It is not a definition of
  homology and does not assert absence of homology below 90 %.** The 481
  Ryazansky S3 representatives (UCLUST 90 %) are an upstream reduction, not this
  unit. Historical roadmap reference: B2.
- **APAZ partition `apaz_partition_v2_mmseqs90_80`** — whole split groups
  assigned to `FINAL_HOLDOUT` first, then `CALIBRATION`, then `BUILD`, by a
  deterministic subset-sum search that aborts rather than silently adjusting an
  unreachable target. Hard negatives HisG (`PF01634`) and EIIB (`PF00367`)
  never enter `BUILD`. **The FINAL_HOLDOUT partition is frozen — no
  redistribution may be motivated by HMM performance.** Supersedes v1 (§3).
  Historical roadmap reference: B2.
- **6 APAZ profile HMMs** (global + Ia / Ib / IIa / IIb / III) built **only from
  BUILD** sequences, bit-for-bit reproducible (volatile `COM` / `DATE` lines
  stripped), re-validated structurally from the on-disk artifact. Generated
  `.hmm` files are gitignored — a generated artifact is not the versioned
  source. Historical roadmap reference: B3.
- **Ryazansky Table S1 clade catalog** — the 1010-protein pre-redundancy-reduction
  set (the 721 nonredundant representatives, their MID-PIWI alignment, and the
  tree are **not published**). MID-PIWI coordinate convention proven from the
  data before slicing. NCBI record status recorded (LIVE / SUPPRESSED /
  DEAD_REPLACED); nothing excluded. Architecture is kept separate from clade.
  Historical roadmap reference: B4.2.
- **Quarantine of conflicting references** — AfAgo (`WP_010878815.1`) and SiAgo
  (`WP_012735993.1`) held at `curated_pago_clade = UNRESOLVED` (architecture
  SHORT vs Ryazansky truncated-long). NgAgo kept as Ryazansky `longA` (the
  recall panel's LONG_B label contradicts the source it cites and is treated as
  a probable recall-panel error; the Phase A artifact is not modified). The two
  Ryazansky `unkn` proteins -> `UNRESOLVED`. Historical roadmap reference: B4.2.
- **MID-PIWI high-similarity split groups** — the B2 workflow applied to the
  1002 extracted MID-PIWI regions: 701 groups, 697 partition-eligible, **0
  cross-clade groups at 90 / 80**. Cross-clade structure appears only at 50 %
  identity and is recorded as **diagnostic evidence only** — no quarantine, no
  exclusion. Byte-identical across two independent runs. Historical roadmap
  reference: B4.3.
- **Formal annotation ontology** — see §2.7.

### 2.7 Annotation ontology

- **Content:** `ago_family` in {PAGO, PIWI_RE, UNRESOLVED} — **PIWI-RE is a
  family, not a pAgo clade**; `pago_clade` in {LONG_A, LONG_B, SHORT,
  UNRESOLVED}, defined only within PAGO and only by phylogenetic placement plus
  three calibrated thresholds; `ago_detection_status` /
  `detection_evidence_status` as an epistemic axis where a triage-only MID-PIWI
  signal never becomes `DETECTED`; disjoint evidence namespaces (Pfam / APAZ /
  clade-HMM); `catalytic_site_status = UNKNOWN` in version 1; the
  BUILD / CALIBRATION / FINAL_HOLDOUT protocol for custom HMMs and the
  TREE_BUILD / PLACEMENT_CALIBRATION / PLACEMENT_HOLDOUT protocol for placement;
  reference-label evidence tiers as provenance, not interchangeable ground
  truth.
- **Status:** drafted on the staging branch as
  `documentation/pago_annotation_ontology.md`; preserved under `history/` by the
  backfill PR; becomes authoritative only when reintegrated.

---

## 3. Superseded, historical, and currently integrated legacy states

Two axes, kept separate: **epistemic / methodological status** and **integration
status**. An item can be scientifically superseded and still be the behaviour
that runs in `master`.

| Item | Epistemic / methodological status | Integration status | Replacement / note |
| --- | --- | --- | --- |
| Query `PIWI[All Fields] AND Bacteria[Organism]` and its 41,345-record filtered set | `SUPERSEDED` by the validated annotation-enriched design | `INTEGRATED` — **current `master` behaviour** | The annotation-enriched query and candidate set (§2.2), on the staging branch, replace it once reintegrated. |
| APAZ partition v1 (split unit keyed on an accession hash; one singleton per accession; high-similarity groups could cross BUILD / CALIBRATION) | `SUPERSEDED` | `NEVER INTEGRATED` (staging only) | `apaz_partition_v2_mmseqs90_80` (§2.6). v1 was never used for predictive evaluation. |
| Exact-accession-only recall matching | `SUPERSEDED` | `NEVER INTEGRATED` (staging only) | Three-tier matching with `SEQUENCE_SHA256` (§2.5). |
| `PIWI_RE` as a value of `pago_clade` (early fixtures / tests) | `SUPERSEDED` | `NEVER INTEGRATED` (staging only) | `PIWI_RE` as a value of `ago_family`; PIWI-RE rows carry `pago_clade = UNRESOLVED` (§2.7). |
| Branch `(fix)-change-canon-id-from-UID-to-accession-version` | `HISTORICAL` — a reversal that was never designed or argued | `NEVER INTEGRATED` | UID remains canonical (§1.1). Kept as a historical pointer, not a competing decision. |
| Minimal GA-KMeans clustering prototype | `HISTORICAL` / abandoned | `WAS INTEGRATED, later removed` from `master` | — |
| Branches `backup/original-messages`, `(perf)-improve-ncbi-fetch-performance` | `HISTORICAL` — a rewritten parallel history of the retrieval rework (§2.1) | `NEVER INTEGRATED` | Kept for provenance; the motivation for the rewrite is not documented. |
| `HisG PF00815` (as written in the private architecture plan) | `CORRECTED` to `PF01634` | `N/A` (plan text, not code) | Fixed in the staging code and curation notes before use. |

---

## 4. Open scientific questions

Not resolved. Do not assume an answer.

1. **Origin and justification of the 90 / 80 threshold.** The 90 % identity
   echoes Ryazansky's UCLUST 90 % reduction; the 80 % coverage requirement
   appears operational. Whether 90 / 80 was chosen deliberately to match the
   literature or independently is only partially recoverable.
2. **Canonical TtAgo accession for the recall panel.** The panel's
   `WP_011174533.1` is a non-matching *T. thermophilus* isolate; the
   crystallographic HB8 TtAgo is Table S1 `YP_145307.1` (dead ->
   `WP_011229221.1`). The panel entry is a pending curation decision.
3. **Construction of the MID-PIWI reference tree.** Ryazansky's 721-representative
   set, its MID-PIWI MSA, and the tree are unpublished. The reference tree,
   alignment, and substitution model have not been built; the method is not
   yet decided or executed.
4. **Independent PIWI-RE reference arm.** A PIWI-RE curation arm sourced from
   Burroughs 2013 plus the recall set, independent of Ryazansky, is planned but
   not started (historical roadmap reference: B4.4).
5. **Enumeration of the APAZ v1 cross-partition leaks.** Only one worked example
   is documented; the full set and the exact discovery procedure could not be
   reconstructed (no versioned v1 artifact).
6. **Sequence-based discovery route** (HMM / PSI-BLAST over RefSeq) to reach
   pAgos not annotated with "PIWI" / "Argonaute" — out of current scope,
   documented as future work.

---

## 5. Standing limitations

- **Execution artifacts are not versioned.** The 52,473-record acquisition
  chain and the Phase B HMM builds are regenerable but absent from the
  repository; a clean clone does not reproduce them without re-running against
  NCBI / external tools.
- **Conditional external tools.** SWeeP 2.1.3.0 (not committed), and — when
  their phases arrive — MMseqs2 and EPA-ng, have no pip-wheel contract; those
  steps are conditionally reproducible.
- **The legacy textual QC layer** (§1.4) predates the ontology (§2.7) and its
  labels must not be promoted to ontology fields.
