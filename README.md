# pAgo Project

**Reproducible research software for prokaryotic Argonaute (pAgo) candidate discovery, evidence curation, and exploratory sequence analysis.**

This project investigates how to construct scientifically defensible protein candidate and reference datasets before attempting functional prediction. It combines NCBI protein retrieval, sequence and metadata processing, curated reference resources, similarity-aware data partitioning, exploratory protein embeddings, and reproducible computational workflows. A final validated functional or catalytic-activity predictor is **not** provided.

## Project at a glance

| Area | Implemented in this repository | Scientific boundary |
| --- | --- | --- |
| Acquisition | NCBI protein UID and XML retrieval, metadata extraction, and FASTA materialization with immutable snapshots | Text-query hits are candidates, not confirmed pAgos |
| Curation | Metadata-based evidence inventory and conservative filtered datasets; pinned Pfam and APAZ reference resources | Legacy QC labels are audit features, not ground-truth classes |
| Reference analysis | Reconstructed 1,010-entry literature pAgo catalog, MID-PIWI extraction, similarity-group resources, and a 21-entry query-recall panel | Query recall measures retrieval coverage, not detector accuracy |
| Exploratory ML | SWeeP sequence embeddings, PCA, and KMeans analysis | Clusters and silhouette scores do not establish functional classes |
| Software engineering | Versioned manifests, SHA-256 integrity checks, reusable Python modules, regression tests, and Linux/Windows CI | Computational reproducibility does not establish biological validity |

The audited annotation-enriched NCBI search described in the scientific documentation retrieved **52,473 records**; that execution dataset remains local and is not committed as a complete public artifact. The frozen Ryazansky reference reconstruction includes **1,010 entries**, with MID-PIWI regions extracted for **1,002**. These counts describe distinct datasets and should not be combined as a classifier benchmark.

## What a reviewer can inspect

1. **Problem framing and scientific limitations:** [`documentation/SCIENTIFIC_STATE.md`](documentation/SCIENTIFIC_STATE.md) describes dataset meaning, evidence boundaries, unresolved questions, and where predictive validation remains pending.
2. **Reusable implementation:** [`src/pago_pipeline/`](src/pago_pipeline/) contains the Python modules for retrieval, snapshots, curation, references, and exploratory analysis.
3. **Tests and continuous integration:** [`tests/`](tests/) and [`.github/workflows/ci.yml`](.github/workflows/ci.yml) show automated checks and the declared Linux/Windows test environments.
4. **Traceability:** [`documentation/decisions/`](documentation/decisions/) records methodological and engineering choices, while [`documentation/PROJECT_LEDGER.md`](documentation/PROJECT_LEDGER.md) records major decision events.
5. **Notebook execution order:** [`notebooks/`](notebooks/) contains the numbered orchestration notebooks; reusable logic is implemented in Python modules rather than duplicated across notebooks.

## Pipeline architecture

~~~text
NCBI search / protein UIDs
  -> XML snapshots and metadata
  -> sequence FASTA and QC evidence
  -> conservative candidate subsets
  -> exploratory embeddings -> PCA / KMeans / visualization

Independent reference/evaluation instruments:
  Pfam / APAZ resources -> profile-HMM construction and split-group controls
  literature pAgo catalog -> MID-PIWI regions and similarity groups
  curated query-recall panel -> retrieval-coverage measurement
~~~

Each arrow represents a data-processing relationship, not a claim that the complete sequence is runnable without external dependencies or that exploratory outputs provide validated biological labels.

## Research status and attribution

**Implemented and inspectable:** the acquisition, data-quality, reference-resource, snapshot/provenance, query-recall, and exploratory-analysis components described above. **Pending as validated research outputs:** a final supervised pAgo functional classifier, calibrated catalytic-activity prediction, experimentally validated guide/target specificity, and a validated wet-lab activity result.

This repository is maintained by **Daniel Harley** in the context of research at the **Federal University of Paraná (UFPR)**, supervised by **Roberto Tadeu Raittz**, with scientific collaboration involving **Diego Bonatto (UFRGS)**. Git commit and pull-request history provide evidence of repository changes, but commit authorship alone does not establish individual authorship of every scientific idea or experimental contribution. Additional contributor-specific responsibilities should be attributed from verified project records.

The source code is publicly readable. No license file is currently included, so redistribution or reuse permissions should not be assumed. The SWeeP dependency is intentionally absent from Git because of intellectual-property restrictions; the limitations of public reproduction are detailed below.

## Recreating the environment

Recreating the environment means rebuilding the system, Python, and project
conditions needed to run the notebooks and validation scripts from a clean
machine.

The current public contract targets Python 3.14.0, matching the environment
used to generate `requirements.lock.txt`.

### 1. Clone the repository

```powershell
git clone https://github.com/DanielHarley/pAgo-project.git
cd pAgo-project
```

### 2. Materialize raw XML snapshots only when needed

Raw NCBI XML snapshots are tracked with Git LFS. Routine code development,
tests, and CI can operate with the small LFS pointer files. Materialize the
full raw XML bytes only when a workflow needs to inspect or reprocess them:

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

The raw-data verification checks local `data/01-raw` identities against the
SHA-256 values recorded in snapshot manifests. For a materialized XML file it
hashes the file bytes. For an unmaterialized Git LFS pointer it compares the
pointer's content-addressed SHA-256 OID with the manifest hash, so routine CI
does not need to download the 363 MB raw XML object.

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

Notebooks are orchestration layers only. Reusable logic lives under
`src/pago_pipeline/`, where each snapshot-producing stage owns its validation,
manifest, and reuse behavior.

The import cell of notebooks `01` to `09` removes every cached
`src.pago_pipeline` module before importing, so rerunning that cell in a
running kernel picks up edits to any file under `src/pago_pipeline/`, including
modules that are only imported transitively. Objects created before the rerun
keep the old code, so restart the kernel after changing a class, dataclass
fields, or `Enum` members, after editing
`src/pago_pipeline/notebook_reload.py` itself, or whenever a notebook raises an
unexpected `ImportError` or `AttributeError` after a code change.

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

See [`documentation/pago_qc.md`](documentation/pago_qc.md) for the QC labels,
decisions, output contracts, and snapshot integrity rules, and
[`documentation/README.md`](documentation/README.md) for the project's decision
history, scientific state, and architecture decision records.

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
