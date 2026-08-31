# pAgo Project

This repository contains a reproducible pipeline for collecting NCBI protein
records related to pAgo candidates, materializing raw XML snapshots, extracting
protein metadata, curating pAgo evidence and labels, materializing filtered
datasets, generating sequence embeddings, and running exploratory PCA and
KMeans analyses.

## Recreating the environment

Recreating the environment means rebuilding the system, Python, and project
conditions needed to run the notebooks and validation scripts from a clean
machine.

The current public contract targets Python 3.14.0, matching the environment
used to generate `requirements.lock.txt`.

### 1. Clone the repository

```powershell
git clone <repo-url>
cd pAgo-project
```

### 2. Enable Git LFS

Raw NCBI XML snapshots are tracked with Git LFS.

```powershell
git lfs install
git lfs pull
```

### 3. Create and activate a virtual environment

```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\Activate.ps1
```

### 4. Upgrade packaging tools

```powershell
python -m pip install --upgrade pip setuptools wheel
```

### 5. Install locked dependencies

```powershell
python -m pip install -r requirements.lock.txt
```

`requirements.in` lists the direct project dependencies. `requirements.lock.txt`
captures the exact installed package versions used for environment recreation.

### 6. Configure environment variables

```powershell
Copy-Item .env.example .env
```

Edit `.env` and set at least:

```text
NCBI_EMAIL=your.email@example.com
```

`NCBI_API_KEY` is optional. TLS certificate variables are only needed on
networks that require custom certificate authority configuration.

### 7. Validate the environment

```powershell
python scripts/print_environment.py
python scripts/verify_raw_data.py
python -m unittest discover -s tests -q
```

The raw-data verification checks local `data/01-raw` files against the SHA-256
values recorded in snapshot manifests.

## Project configuration

Public pipeline defaults are stored in `config/config.yaml`. Secrets and local
machine overrides should stay out of Git and should be placed in `.env` or a
local ignored file.

## Pipeline order

Run notebooks in numeric order:

1. `notebooks/00_ncbi_protein_ids_query.ipynb`
2. `notebooks/01_ncbi_protein_xml_snapshot.ipynb`
3. `notebooks/02_protein_metadata_extraction_to_csv.ipynb`
4. `notebooks/03_protein_metadata_csv_qc.ipynb`
5. `notebooks/04_metadata_to_fasta.ipynb`
6. `notebooks/05_sweep_genes_embeddings.ipynb`
7. `notebooks/06_pca_kmeans_analysis.ipynb`
8. `notebooks/07_pca_3d_plot.ipynb`
9. `notebooks/08_pago_qc_evidence_inventory.ipynb`
10. `notebooks/09_pago_qc_filtered_datasets.ipynb`
11. `notebooks/10_dataset_audit.ipynb`

Notebooks are orchestration layers only. Reusable logic lives under
`src/pago_pipeline/`, where each snapshot-producing stage owns its validation,
manifest, and reuse behavior.

## Annotation-enriched candidate set (from notebook 10)

`notebooks/10_dataset_audit.ipynb` starts a second NCBI retrieval with a broader
text query
(`(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR
Archaea[Organism])`) and materializes the audit artifacts around it:

- `ncbi_esearch_preflight` — one ESearch (with History) plus a small sample
  fetch, recording `Count`, `QueryTranslation`, and a sample-record QC. The full
  retrieval is blocked when `Count` exceeds `max_uid_count` until reviewed.
- UID -> XML -> flattened metadata snapshots, reusing the existing NCBI modules
  with `annotation_enriched_candidate_set`-suffixed snapshot roots.
- `query_reference_recall` — how many curated reference pAgos
  (`tests/fixtures/query_recall_reference_set.csv`, see the sibling
  `_curation_notes.md`) the text query recovered, overall, per MID-PIWI clade
  (LONG_A / LONG_B / SHORT), and for the PIWI-RE family. A reference is matched
  by accession.version, then by base accession, then by identical protein
  sequence (`sequence_sha256`, whitespace-stripped + uppercased) — so a protein
  retrieved under a different accession still counts. Two readings are reported:
  `exact_accession_recall` (strict) and `retrieval_equivalent_recall` (also
  sequence identity). A stratum with zero references reports `NOT_EVALUABLE`,
  never `0.0`. The matching methodology is pinned by `matching_strategy_sha256`.
  Each reference row carries a `reference_label_evidence` level (EXPERIMENTAL /
  LITERATURE_PHYLOGENETIC / CURATED_COMPUTATIONAL / DATABASE_ANNOTATION).
- `pago_technical_prefilter` — excludes only technically unusable records
  (missing/invalid sequence, missing `protein_uid`, technical duplicates). It
  never removes a record by sequence length or by NCBI annotation text; length
  outside a warning band sets `length_warning=True` and the record is kept.
- `derived_protein_fasta` — the retained proteome as a FASTA snapshot with its
  own `artifact_type` (`derived_protein_fasta_snapshot`), full transformation
  provenance (`derived_from_*`, `record_selection_rule`, `source_record_ids_sha256`,
  `record_order`) and load-time re-validation of the record multiset and order.

The candidate set is *annotation-enriched*, not a full pAgo universe: a text
query recovers proteins already annotated with the terminology. A
sequence-based discovery route (HMM / PSI-BLAST over RefSeq) is future work.

## pAgo QC curation

The pAgo QC workflow separates observation, interpretation, and operational
filtering:

```text
metadata.csv
  -> evidence_flags.csv
  -> labelled_records.csv
  -> filtered dataset CSVs
```

The evidence and label layer is implemented by:

- `src/pago_pipeline/pago_qc.py`
- `src/pago_pipeline/pago_qc_snapshot.py`
- `notebooks/08_pago_qc_evidence_inventory.ipynb`

It writes snapshot artifacts under:

```text
data/03-features/pago_qc/evidence_inventory/
  latest/
  snapshots/
```

The filtered dataset layer is implemented by:

- `src/pago_pipeline/pago_qc_filter.py`
- `src/pago_pipeline/pago_qc_filter_snapshot.py`
- `notebooks/09_pago_qc_filtered_datasets.ipynb`

It writes snapshot artifacts under:

```text
data/03-features/pago_qc/filtered_datasets/
  latest/
  snapshots/
```

The filtered datasets currently materialized are:

- `classic_pago_high_precision_records.csv`
- `classic_pago_review_records.csv`
- `piwi_re_records.csv`
- `excluded_records.csv`
- `filtered_dataset_counts.csv`

`excluded_records.csv` means excluded from the conservative classic pAgo
positive dataset. It is not a biologically validated negative class.

See `docs/pago_qc.md` for the QC labels, decisions, output contracts, and
snapshot integrity rules.

## Generated feature outputs

Feature outputs are generated under `data/03-features/` and are intentionally
ignored by Git. Each generated feature stage uses a `latest/` directory for the
most convenient current artifact and a `snapshots/` directory for immutable
historical runs.

The SWeeP embedding notebook writes to:

```text
data/03-features/sweep_genes/
  latest/
  snapshots/
```

The pAgo QC notebooks write to:

```text
data/03-features/pago_qc/
  evidence_inventory/
  filtered_datasets/
```

## SWeeP dependency

The embedding step depends on SWeeP 2.1.3.0 at
`scripts/sweep-2.1.3.0/`. That dependency is intentionally not committed because
of intellectual property restrictions.

Public end-to-end reproducibility is therefore conditional for the SWeeP-backed
embedding step: users must obtain SWeeP 2.1.3.0 through an authorized channel and
place it at the configured path before running the embedding notebook.
