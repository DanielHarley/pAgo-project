# Raw Data Reproducibility

This directory contains the raw NCBI snapshots needed to reproduce the analysis
from the same source records.

## Tracking policy

- `latest/` is tracked with normal Git and always reflects the published dataset.
- Runs under `snapshots/` are ignored by default: each one is a regenerable
  multi-hundred-MB artifact, and the repository publishes only curated runs.
  Publishing one is deliberate -- add a `!` line for its directory in
  `.gitignore`, then commit it together with the `latest/` update that points
  at it. The two must land in the same commit, because `latest/manifest.json`
  names its source run in `immutable_snapshot_relative_path`.
- Track raw XML snapshots with Git LFS. Pushed LFS storage is never reclaimed,
  which is the other reason snapshot runs are published one at a time.
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
