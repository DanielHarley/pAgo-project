# pAgo QC Curation

This document describes the pAgo QC workflow added after the initial metadata,
FASTA, embedding, PCA, and KMeans stages. The goal is to turn a broad NCBI query
result into auditable evidence flags, conservative biological labels, and
materialized filtered datasets for downstream analysis.

The workflow is intentionally split into three contracts:

```text
1. Evidence layer
   metadata.csv -> evidence_flags.csv

2. Label and decision layer
   evidence_flags.csv -> labelled_records.csv

3. Filtered dataset layer
   labelled_records.csv -> curated subset CSVs
```

This separation keeps raw observations, biological interpretation, and
downstream materialization from being mixed in one notebook.

## Evidence and Label Snapshot

Implementation:

- `src/pago_pipeline/pago_qc.py`
- `src/pago_pipeline/pago_qc_snapshot.py`
- `notebooks/08_pago_qc_evidence_inventory.ipynb`

Output root:

```text
data/03-features/pago_qc/evidence_inventory/
  latest/
  snapshots/
```

Snapshot manifest contract:

```text
artifact_type = "pago_qc_evidence_inventory"
snapshot_format_version = "2.0"
```

Main output files:

- `evidence_flags.csv`
- `evidence_counts.csv`
- `labelled_records.csv`
- `label_counts.csv`
- `filter_decision_counts.csv`
- `suspicious_terms.csv`
- `top_region_names.csv`
- `top_products.csv`
- `manifest.json`

The manifest records the source metadata CSV path and SHA-256, output file
SHA-256 hashes, label columns, the length bin column, deprecated aliases, and
the snapshot format version. Loading or reusing a snapshot validates the
manifest, required output files, output hashes, source metadata hash when
provided, and `snapshot_format_version`.

## Evidence Flags

Evidence flags are deterministic boolean or categorical columns derived from
NCBI metadata fields. They are observations, not final labels.

Important PIWI-related flags include:

- `has_classic_piwi_region`
- `has_ppiwi_re_region`
- `has_any_piwi_related_region`
- `has_classic_piwi_text_anywhere`
- `has_ppiwi_re_text_anywhere`
- `has_any_piwi_related_text_anywhere`
- `has_any_piwi_related_evidence`

Important methyltransferase and conflict flags include:

- `has_ubig_term`
- `has_sam_methyltransferase_term`
- `is_probable_methyltransferase_noise`
- `has_conflicting_piwi_and_methyltransferase_evidence`
- `is_methyltransferase_noise_or_conflict`

Length and annotation quality flags include:

- `length_bin`
- `is_short_lt_300`
- `is_possible_partial_300_599`
- `is_compatible_long_pago_600_900`
- `is_large_or_fused_901_1300`
- `is_long_outlier_gt_1300`
- `has_partial_or_fragment_keyword`

Deprecated aliases are kept for compatibility while new code uses clearer
names:

- `has_piwi_region` -> `has_any_piwi_related_region`
- `has_piwi_text_anywhere` -> `has_classic_piwi_text_anywhere`
- `has_active_site_annotation` -> `has_active_or_functional_site_annotation`

## Labels and Decisions

`labelled_records.csv` adds interpretive and operational columns:

- `primary_label`
- `qc_decision`
- `confidence_score`
- `rationale`

Primary labels:

- `classic_piwi_candidate`
- `piwi_re_candidate`
- `methyltransferase_noise_or_unresolved_conflict`
- `low_evidence_or_unrelated`

QC decisions:

- `include`
- `review`
- `exclude`
- `separate_dataset`

Decision semantics:

- `include` is reserved for conservative classic pAgo candidates.
- `review` is reserved for classic candidates with weaker or less typical
  evidence, such as possible partials or atypical length bins.
- `separate_dataset` is used for PIWI-RE candidates. PIWI-RE is treated as a
  separate class rather than a doubtful classic pAgo.
- `exclude` means excluded from the conservative classic pAgo positive dataset.
  It does not mean biologically validated negative.

Rule ordering is conservative:

1. Methyltransferase noise or unresolved conflict blocks positive assignment.
2. PIWI-RE evidence is assigned before classic PIWI evidence.
3. Classic PIWI evidence becomes `classic_piwi_candidate`.
4. Remaining records become `low_evidence_or_unrelated`.

`confidence_score` is an ordinal audit score, not a calibrated probability.

## Filtered Dataset Snapshot

Implementation:

- `src/pago_pipeline/pago_qc_filter.py`
- `src/pago_pipeline/pago_qc_filter_snapshot.py`
- `notebooks/09_pago_qc_filtered_datasets.ipynb`

Output root:

```text
data/03-features/pago_qc/filtered_datasets/
  latest/
  snapshots/
```

Snapshot manifest contract:

```text
artifact_type = "pago_qc_filtered_datasets"
snapshot_format_version = "1.0"
```

Main output files:

- `classic_pago_high_precision_records.csv`
- `classic_pago_review_records.csv`
- `piwi_re_records.csv`
- `excluded_records.csv`
- `filtered_dataset_counts.csv`
- `manifest.json`

Dataset definitions:

```text
classic_pago_high_precision:
  primary_label == "classic_piwi_candidate"
  and qc_decision == "include"

classic_pago_review:
  primary_label == "classic_piwi_candidate"
  and qc_decision == "review"

piwi_re:
  primary_label == "piwi_re_candidate"
  and qc_decision == "separate_dataset"

excluded:
  qc_decision == "exclude"
```

Every labelled record must be assigned to exactly one filtered dataset. The
filter module validates this invariant and fails if any record is unassigned or
multiply assigned.

The filter snapshot manifest records:

- source QC snapshot directory
- source QC artifact type
- source QC manifest path and SHA-256
- source QC snapshot format version
- source `labelled_records.csv` path and SHA-256
- filter policy payload and `filter_policy_sha256`
- output file names and SHA-256 hashes
- the explicit semantics of `excluded_records.csv`

Changing the filter policy changes `filter_policy_sha256`, which prevents
accidental reuse of an outdated `latest/` filtered snapshot.

## Current Local Counts

The current local filtered snapshot produced:

```text
classic_pago_high_precision   2,629
classic_pago_review           6,703
piwi_re                      15,056
excluded                     16,957
total                        41,345
```

These counts are useful for audit, but downstream code should depend on the
snapshot manifests and CSVs rather than hard-coded counts.

## Known Boundaries

The QC labels are based on existing NCBI metadata annotations. They are
deterministic and auditable, but they are not equivalent to external domain
validation from HMMER, Pfam, CDD re-search, InterProScan, AGODB, or sequence
motif scans.

`classic_pago_high_precision` prioritizes a practical conservative positive
set. It may include a small number of text-supported candidates without
explicit classic PIWI region annotation when their length is compatible with
long pAgo candidates.

`excluded_records.csv` must not be used directly as a supervised negative class
without additional negative-class curation.

## Next Stage

The next logical pipeline stage is a FASTA export snapshot:

```text
filtered_datasets + protein_metadata.csv -> filtered FASTA files
```

That stage should keep its own manifest and provenance, including the source
filtered dataset snapshot and the source metadata CSV hash.
