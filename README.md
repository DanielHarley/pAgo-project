# pAgo Project

This repository contains a reproducible pipeline for collecting NCBI protein
records related to pAgo candidates, materializing raw XML snapshots, extracting
protein metadata, generating sequence embeddings, and running exploratory PCA
and KMeans analyses.

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

## SWeeP dependency

The embedding step depends on SWeeP 2.1.3.0 at
`scripts/sweep-2.1.3.0/`. That dependency is intentionally not committed because
of intellectual property restrictions.

Public end-to-end reproducibility is therefore conditional for the SWeeP-backed
embedding step: users must obtain SWeeP 2.1.3.0 through an authorized channel and
place it at the configured path before running the embedding notebook.
