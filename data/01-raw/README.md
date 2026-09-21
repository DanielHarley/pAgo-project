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

A normal checkout may leave `protein_records.xml` as a Git LFS pointer. That
is sufficient for repository-level identity verification because the pointer
contains the SHA-256 OID of the exact raw object.

Verify the checked-out identities against the recorded manifests:

```powershell
python scripts/verify_raw_data.py
```

For an unmaterialized XML pointer, the verifier compares the LFS OID with the
manifest's `xml_file_sha256`. To inspect, parse, or perform a full byte-level
verification of the raw XML itself, materialize the object explicitly and run
the same verifier again:

```powershell
git lfs install
git lfs pull
python scripts/verify_raw_data.py
```

The verification covers the recorded SHA-256 identities for `protein_uids.txt`
and `protein_records.xml` wherever those files are referenced by a raw snapshot
manifest.
