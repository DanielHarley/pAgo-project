# ADR-0018 - Support explicit UID and accession.version domains in NCBI acquisition

## Context

The ESearch workflow uses protein UIDs for retrieval, XML validation, and
snapshot provenance. Acquisition from a supplied list also needs to preserve
the exact revision of each protein record identified by `accession.version`.
See [ADR-0001](0001-canonical-ncbi-retrieval-identifier-is-the-protein-uid.md)
and [issue #26](https://github.com/DanielHarley/pAgo-project/issues/26).

## Decision

Each acquisition snapshot declares one `identifier_type`, either `uid` or
`accession.version`. ESearch produces UID snapshots. A supplied list can
declare either domain. EFetch and XML validation use the declared domain and
check record identity and multiplicity against the saved list. All supported
`accession.version` candidates in a record must agree, including versioned
accessions in `GBSeqid`.

| Artifact in format 1.1 | Required integrity metadata |
| --- | --- |
| Identifier list | `protein_identifiers_sha256`, `normalized_protein_identifier_count` |
| Consolidated XML | The identifier-list fields, `xml_file_sha256`, `consolidated_record_count` |

XML creation selects the source through `latest`, checks its manifest against
the referenced immutable `snapshots/<run>/manifest.json`, and loads the
identifiers from that run. The resolver captures the source manifest and its
SHA-256 before EFetch. Publication uses this captured provenance throughout
the run, including when another run replaces `latest` during the fetch. The
identifier list, fetch result, and persisted XML must agree before promotion.

Identifier and XML paths declared in a manifest are relative to the owning
snapshot. The loader and auditor reject paths that escape that directory,
including through a symbolic link. Directory names use the UTC rendering of
an ISO 8601 retrieval timestamp with an explicit offset.

`reuse_latest_or_create` with a source root checks the selected source
manifest hash, domain, and identifier list before opening and validating
the cached XML. A changed source triggers XML creation from the newly
selected run. Without a source root, an existing XML snapshot can still be
reused offline. `reuse_latest` always selects the existing XML snapshot
directly.

## Compatibility

Unversioned and `1.0` manifests default to UID when `identifier_type` is
absent. Reuse retains legacy field fallbacks and accepts missing optional
counts. XML creation requires a matching source list hash and count. For a
legacy source without an immutable directory reference, XML creation uses
the manifest bytes and SHA-256 captured from `latest` at selection. A
declared reference must identify the matching immutable run. When matching
cached XML provenance, a legacy source manifest hash may reflect LF or
CRLF line endings introduced by a Git checkout.

Format `1.1` requires an explicit domain and its canonical integrity fields.
An explicit version must be exactly `1.0` or `1.1`.
New writers also provide the legacy UID-named aliases and file copy for
existing readers. These aliases carry identifiers in the declared domain.
The loaded list's hash preserves its order and duplicate occurrences.

The auditor and loader share validation of materialized `1.1` snapshots.
For an XML LFS pointer, the auditor checks the OID against the XML hash and
checks the manifest and identifier list. Record validation requires the
materialized XML. See [ADR-0016](0016-keep-large-raw-xml-lfs-objects-out-of-routine-ci-checkouts.md).

## Rationale

An explicit domain preserves versioned accessions and allows exact record
matching while keeping the UID workflow compatible. Capturing the manifest
identity before the request binds the XML to the source that produced it.
Contained paths and canonical UTC directory names make snapshots portable
and keep their inputs within the recorded run. See
[ADR-0002](0002-immutable-snapshot-and-provenance-architecture.md).
