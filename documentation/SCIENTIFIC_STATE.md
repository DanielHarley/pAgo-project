# Scientific state

This document summarizes the project's current scientific interpretation: what
the datasets and analyses mean, what they do not establish, and which questions
remain unresolved.

It does not mirror Git state. Use Git, PRs, and commits to determine where an
implementation lives or whether it has been merged.

Detailed rationale and provenance live in [`decisions/`](decisions/) and
[`PROJECT_LEDGER.md`](PROJECT_LEDGER.md).

## Dataset acquisition and meaning

### Text-query datasets are candidate sets, not confirmed pAgos

NCBI text queries recover proteins whose records contain the searched
terminology. They do not define the biological universe of pAgos.

Two different acquisition contexts must not be conflated:

- the earlier narrow query `PIWI[All Fields] AND Bacteria[Organism]`, whose
  downstream legacy QC run contained **41,345** records
- the later annotation-enriched query
  `(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR Archaea[Organism])`,
  whose audited 2026-08-30/31 execution returned **52,473** records

The 52,473-record set is appropriately described as **candidate proteins
enriched by textual annotation associated with PIWI/Argonaute**. It is not an
estimate of the total number of pAgos.

A strictly technical prefilter may remove technically unusable records, but
annotation text and sequence length alone are not sufficient biological reasons
to declare a candidate a non-pAgo. In the audited 52,473-record execution, the
technical prefilter removed 0 records and marked 8,299 with a length warning.

See ADR-0004 and the Phase A historical audit.

## Query-recall panel

The query-recall panel contains 21 curated references: 14 pAgo references across
LONG_A, LONG_B, and SHORT, plus 7 PIWI-RE references.

It is a **query-coverage instrument**, not a classifier gold standard and not a
final holdout.

Reference matching keeps three distinct possibilities:

1. exact `accession.version`
2. the same base accession
3. normalized-sequence SHA-256 equivalence

Sequence-hash equivalence means that the same normalized protein sequence was
retrieved under another accession. It does not by itself establish identical
locus identity, record provenance, or biological context.

Six of the seven PIWI-RE panel references were selected using PIWI-RE profile
models. They therefore cannot be reused uncritically as an independent
validation set for a profile-based PIWI-RE detector; that would be circular.

See ADR-0005 and ADR-0014.

## Legacy metadata-derived QC

The workflow documented in [`pago_qc.md`](pago_qc.md) derives deterministic
flags and labels from NCBI metadata.

Those labels are **audit features, not biological ground truth**.

Important boundaries:

- `excluded_records.csv` means excluded from the conservative classic-pAgo
  positive set; it is not a biologically validated negative class
- PIWI-RE is handled separately rather than treated as a doubtful classic pAgo
- metadata-derived labels are not equivalent to domain validation with HMMER,
  Pfam, CDD, InterProScan, or experimentally supported biological labels
- the historical 2,629 / 6,703 / 15,056 / 16,957 filtered counts belong only to
  the earlier 41,345-record query

These QC labels must not be promoted directly into the formal annotation
ontology.

## SWeeP, PCA, and KMeans are exploratory

SWeeP converts biological sequences into fixed-length numerical vectors. PCA
(Principal Component Analysis) linearly reduces dimensionality, and KMeans
partitions points around centroids.

The SWeeP → PCA → KMeans workflow is an **exploratory analysis of sequence
space**, not a validated biological classifier.

Silhouette scores describe geometric cohesion/separation of a partition. ARI
(Adjusted Rand Index) measures agreement between clusterings with chance
correction. Neither establishes that a cluster is a biological class.

For the historical 41,345-record run, the local summary recorded approximately:

- selected `k = 3`
- sampled silhouette ≈ 0.656
- selected-init ARI minimum ≈ 0.998
- selected-subsample ARI minimum ≈ 0.993
- Cramér's V between cluster and legacy primary label ≈ 0.386
- Spearman ρ between PC2 and sequence length ≈ 0.878

The strong PC2/length association is an important confounding signal. No
biological validation of those clusters has been demonstrated.

See ADR-0006.

## Reference-layer boundaries

### Pfam 38.2 bundle

The selected Pfam 38.2 reference bundle contains ten version-pinned profiles:
PIWI, PAZ, ArgoN, ArgoL1, ArgoL2, ArgoMid, SIR2, TIR_2, TIR, and Mrr_cat.

Pinned accession versions and SHA-256 hashes establish reference identity and
support reproducible integrity checks. Structural parsing or byte-level
reproducibility does **not** establish sensitivity, specificity, or predictive
performance for pAgo detection.

The original rationale for exactly these ten profiles and for choosing Pfam
38.2 has not been recovered.

See ADR-0008.

### High-similarity split groups

A `split_group` is an operational redundancy-control unit built as a connected
component of a sequence-similarity graph using the declared 90% identity and
80% coverage criterion, with exact-sequence duplicates joined first.

This 90/80 rule is **not a definition of homology**. Homology means shared
evolutionary ancestry and has no universal sequence-identity threshold.

Keeping a whole high-similarity group within one development/evaluation
partition reduces one important leakage mechanism. It does not prove complete
statistical independence between partitions.

The origin of the exact 90/80 choice is only partially recoverable: the 90%
identity value resembles Ryazansky's UCLUST reduction, while the 80% coverage
requirement appears operational.

See ADR-0009 and ADR-0013.

### APAZ references and profile HMMs

Ryazansky Data Set S3 contains 481 representative APAZ domains derived from 632
sequences by a reported UCLUST 90% reduction.

The project's later high-similarity grouping of those 481 representatives
produced 460 split groups. Whole groups are kept together when partitioning the
reference set. HisG (`PF01634`) and EIIB (`PF00367`) are used as hard-negative
families and do not enter the HMM BUILD references.

Six APAZ profile HMMs — global plus Ia, Ib, IIa, IIb, and III — are constructed
from BUILD references only. Their deterministic construction and byte-level
reproducibility do not constitute a performance evaluation.

No sensitivity, specificity, F1, or frozen predictive threshold result should
be claimed from HMM construction alone.

See ADR-0009 and ADR-0010.

### Ryazansky pAgo catalog and MID-PIWI regions

Ryazansky Table S1 contains a 1,010-pAgo catalog before the paper's reported
90%-identity redundancy reduction. The catalog is reintegrated and versioned:
`src/pago_pipeline/resources/clade_seed/` carries the frozen Table S1
workbook, the Data Set S1/S2 diagnostic alignments, the PMC metadata, the
resolved NCBI record-status snapshot, the frozen sequence FASTA, a source
manifest with per-file md5/sha256, the built catalog CSV, and its summary;
`scripts/build_clade_reference_catalog.py` rebuilds the catalog offline and
deterministically from those frozen sources (`--reuse-frozen`).

The reconstruction recovers sequences for all 1,010 accessions (1,010/1,010)
and extracts MID-PIWI regions for 1,002 of them (984 unique regions; 972
canonical and 38 truncated architectures) after proving the coordinate
convention from the table data. This reproducibility is computational: it
demonstrates that the catalog can be rebuilt bit-for-bit from its declared
sources, not that any clade or curation label it carries has been
phylogenetically, classifier-, or HMM-validated.

The reusable mapping of the reported 721 nonredundant tree representatives, the
MID-PIWI tree alignment, and the original Newick tree are not available as
published supplementary artifacts. The project's own validated reference-tree
construction and placement calibration therefore remain separate, unresolved
work (see Open scientific question 3); this catalog does not build or imply a
tree.

AfAgo and SiAgo are retained as unresolved (`curated_pago_clade=UNRESOLVED`,
`curation_status=QUARANTINE`) where source clade and architecture conflict.
NgAgo follows the Ryazansky `longA` source label (`LONG_A`/`OK`); the older
query-recall-panel LONG_B label is treated as a probable curation error rather
than as a reason to rewrite the source catalog. These three outcomes are
implemented as a small, explicit, ADR-0011-cited exception table rather than a
live cross-check against the query-recall panel, which belongs to a separate,
not-yet-reintegrated unit — see ADR-0015.

See ADR-0011 and ADR-0015.

### MID-PIWI high-similarity groups

For the 1,002 extracted MID-PIWI regions, the historical 90/80 workflow produced
701 groups, 697 of them eligible for later partitioning, with no resolved
cross-clade group at the operational 90/80 threshold.

A cross-clade component appears at a 50% identity diagnostic threshold. That
observation is a similarity diagnostic only; it is not a phylogenetic inference
and is not a reason by itself to relabel or quarantine the affected reference.

See ADR-0013.

## Annotation semantics

The project keeps several biological concepts separate.

- `ago_family ∈ {PAGO, PIWI_RE, UNRESOLVED}` describes family identity.
- `pago_clade ∈ {LONG_A, LONG_B, SHORT, UNRESOLVED}` describes pAgo
  phylogenetic clade and is meaningful only within the PAGO family.
- **PIWI-RE is a family, not a pAgo clade.**
- domain architecture does not by itself determine phylogenetic clade
- a triage signal is not automatically a confirmed Argonaute detection
- absence of one profile hit is not automatically a biological negative
- catalytic-site annotation is separate from family, clade, architecture, and
  experimentally demonstrated catalytic activity

The historical staging ontology contains additional proposed placement and
validation machinery. Its historical body is not current authority merely
because it is detailed; current semantic decisions should be checked against the
ADRs and current implementation.

See ADR-0012 and `documentation/history/2026-08-31-pago-annotation-ontology-staging.md`.

## Superseded or corrected interpretations

The following older interpretations should not be reused as current scientific
claims:

- `PIWI_RE` as a value of `pago_clade` → PIWI-RE is represented at the family
  level, while `pago_clade` remains unresolved for PIWI-RE
- exact-accession-only query-recall matching → sequence-equivalent recovery is
  also recorded explicitly while preserving the exact-accession reading
- APAZ accession-hash partition units → high-similarity split groups are kept
  intact across evaluation partitions
- the claim that TtAgo and RsAgo were absent from Ryazansky Table S1 because of
  the 90% reduction → Table S1 is the pre-reduction 1,010-protein catalog and
  contains the relevant entries
- `HisG PF00815` in an old private architecture plan → the curated hard-negative
  family uses `PF01634`

## Open scientific questions

1. **Origin and justification of the 90/80 threshold.** The precise rationale
   for combining 90% identity with 80% coverage has not been fully recovered.
2. **Canonical TtAgo accession for the recall panel.** The panel's
   `WP_011174533.1` is not the crystallographic HB8 reference represented in
   Table S1 as `YP_145307.1` (later `WP_011229221.1`). The panel entry still
   requires curation.
3. **Construction of the MID-PIWI reference tree.** The original reusable tree
   inputs are not available from the paper's supplements, and the project's own
   tree-building/placement method has not yet been scientifically finalized and
   validated.
4. **Independent PIWI-RE reference arm.** A detector-validation reference set
   must be curated independently from primary PIWI-RE literature with explicit
   circularity auditing. The existing recall-panel PIWI-RE rows are not
   automatically eligible.
5. **Complete reconstruction of APAZ v1 leakage.** One worked example is
   documented, but the full historical procedure used to enumerate the reported
   v1 cross-partition leaks has not been recovered.
6. **Sequence-based discovery beyond text annotation.** A sequence-driven route
   such as profile-HMM or PSI-BLAST discovery over an appropriate protein
   database remains future work for finding poorly annotated or divergent
   pAgos.
7. **Mesophily.** Mesophily is a property of the source organism's ecology,
   commonly described using its optimal growth temperature (OGT). It is distinct
   from the catalytic temperature optimum or working range of a specific pAgo
   protein, which requires biochemical evidence for that protein. The two must
   not be conflated.

## Standing limitations

- The 52,473-record acquisition execution remains local rather than versioned;
  that execution is not reproduced by a clean clone without rerunning external
  data or tools. The Ryazansky pAgo catalog and MID-PIWI extraction are
  versioned and reproducible offline from frozen sources (see above); other
  historical reference-layer components (the Pfam HMM bundle, APAZ references
  and profile HMMs, MID-PIWI high-similarity grouping) remain either not yet
  reintegrated to `master` or dependent on external tools/data not distributed
  with the repository.
- SWeeP is not distributed with the repository, and later tools such as MMseqs2
  or EPA-ng require external environments. Reproducibility of those stages is
  conditional on the declared tool/environment contract.
- Legacy textual QC labels must not be promoted directly into formal biological
  ontology fields.
- No final predictive-performance result has been established for the pAgo
  reference/profile/placement instruments described above.
- Protected final holdouts must not be used to make development or calibration
  decisions before the corresponding method and thresholds are frozen.
