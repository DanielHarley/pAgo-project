# Raw Data Reproducibility

This directory contains the raw NCBI snapshots needed to reproduce the analysis
from the same source records.

## Tracking policy

- Commit UID snapshots and manifests with normal Git.
- Track raw XML snapshots with Git LFS.
- Keep SHA-256 hashes in the manifests as the authority for dataset identity.

The XML files are required for strong reproducibility. Re-fetching records from
NCBI using the same UIDs is useful, but it may not reproduce the exact input
data because NCBI records can change after the original retrieval date.

## Expected layout

```text
data/01-raw/
  protein_uid_snapshots/
    latest/
    snapshots/
  protein_xml_snapshots/
    latest/
    snapshots/
```

## Restoring data from a fresh clone

Install Git LFS before cloning or pulling the repository data:

```powershell
git lfs install
git lfs pull
```

Then verify that the local files match the recorded manifests:

```powershell
python scripts/verify_raw_data.py
```

The verification checks the recorded SHA-256 values for `protein_uids.txt` and
`protein_records.xml` wherever those files are referenced by a raw snapshot
manifest.
