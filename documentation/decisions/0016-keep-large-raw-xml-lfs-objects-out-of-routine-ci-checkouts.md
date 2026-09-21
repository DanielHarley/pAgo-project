# ADR-0016 — Keep large raw XML LFS objects out of routine CI checkouts

## Date

2026-09-21

## Context

The repository tracks the frozen NCBI XML snapshot under
`data/01-raw/protein_xml_snapshots/**/*.xml` with Git LFS. The current raw
object is 363,034,184 bytes and its LFS OID is the same SHA-256 recorded in the
snapshot manifests.

The CI workflow previously configured `actions/checkout@v4` with
`lfs: true` in both the Linux `environment` job and the Windows
`windows-snapshot-reuse` job. Each fresh job was therefore configured to fetch
the large raw object before running tests, including jobs whose tests never read
that snapshot. A sampled successful run for PR #38, workflow run
`35238280184`, confirms that checkout executed
`git lfs fetch origin refs/remotes/pull/38/merge`. Because CI runs on both
`push` and `pull_request`, ordinary branch and PR activity could repeatedly
consume Git LFS bandwidth for the same immutable object.

## Decision

Routine CI checkouts do not hydrate Git LFS objects. The workflow sets
`GIT_LFS_SKIP_SMUDGE=1` and checks out with `lfs: false`.

`scripts/verify_raw_data.py` accepts either representation of an XML snapshot:

1. If `protein_records.xml` is materialized, verify the SHA-256 of its bytes
   against `xml_file_sha256` in the manifest.
2. If the checkout contains a Git LFS pointer, verify that the pointer's
   `oid sha256:<digest>` equals the same manifest SHA-256.

Full raw-byte materialization remains an explicit operation through
`git lfs pull` when a researcher needs to inspect, parse, or independently
rehash the XML payload.

## Rationale

A Git LFS OID is content addressed by SHA-256. The repository manifest already
uses SHA-256 as the authority for raw dataset identity. Comparing the pointer
OID with the manifest therefore verifies that Git references the intended
immutable raw object without downloading hundreds of megabytes on every
ephemeral CI runner.

This separates routine repository validation from deliberate raw-data
materialization while preserving the exact artifact identity required for
reproducibility.

## Consequences

Default CI validates the committed identity of the raw XML object and can run
without Git LFS bandwidth. It does not prove during every run that the remote
LFS object is currently retrievable or rehash its full bytes.

A full byte-level check remains available by running `git lfs pull` followed
by `python scripts/verify_raw_data.py`. Existing Git history and the LFS object
are left intact; no history rewrite or scientific dataset migration is part of
this decision.
