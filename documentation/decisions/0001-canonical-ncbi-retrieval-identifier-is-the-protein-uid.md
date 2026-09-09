# ADR-0001 — The canonical NCBI retrieval identifier is the protein UID

## Lifecycle status

`CURRENT`

## Integration status

`INTEGRATED` (via PR #5)

## Epistemic status

`SUPPORTED` — the choice is a technical convention; the code on `master` operates
on protein UIDs consistently.

## Type

`TECHNICAL`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (implementation dated
2026-03-31, PR #5)

## Context

Early retrieval used the NCBI `accession.version` string as the record
identifier. NCBI protein records also carry a numeric UID (a GI number in the
`db=protein` context).

## Problem

The pipeline needs one stable identifier to key ESearch → EFetch → XML
validation and the snapshot manifests.

## Evidence

- PR #5 merge `71ff9b2` — *"(feat) switch NCBI protein snapshot retrieval from
  accession version IDs to UIDs"* — is an ancestor of `master`.
- Implementation commit `a10d89b`.
- `master` code: `src/pago_pipeline/ncbi_snapshot.py` extracts the UID from
  `GBSeqid` fields with `re.fullmatch(r"gi\|(\d+)\|?", ...)` and validates the
  XML record UIDs against the requested UIDs; `fetch_ncbi_protein_uid_snapshot`
  is the retrieval entry point.
- The branch `(fix)-change-canon-id-from-UID-to-accession-version` exists on the
  remote. Its tip `0177e37` is the *"Merge pull request #4"* commit, which
  **predates** PR #5, is ~53 commits behind `master`, and contains **no**
  reversal. It never implemented or integrated a change back to
  `accession.version`.

## Previous state

`accession.version` was the retrieval / query identifier.

## Decision

Use the protein UID as the canonical retrieval and XML-validation identifier.

## Rationale

`RATIONALE_NOT_RECOVERABLE`. The commit and PR record only state the switch; no
design document weighs the two identifiers. NCBI's own documentation states that
when a sequence changes, a **new** GI is assigned and the `accession.version`
suffix is incremented — so neither identifier is invariant across a sequence
revision. No stability advantage is claimed here.

## Alternatives considered

`accession.version` — it was the previous choice and the target of a
never-completed reversal branch. No design document weighs the two.
`HISTORICAL_EVIDENCE_INSUFFICIENT` for any other alternative.

## Consequences

- All raw snapshots and manifests key on `protein_uid`.
- The Phase A query-recall matcher (see ADR-0005) matches references by
  `accession.version` → base accession → normalized-sequence SHA-256. That is a
  **benchmark matching strategy for a specific report**, not the canonical
  retrieval identifier — the two coexist without conflict.

## Limitations

The rationale for the switch is not documented.

## Supersedes

`accession.version` as the canonical retrieval identifier.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

PR #5 (`71ff9b2`). The unrelated, incomplete branch
`(fix)-change-canon-id-from-UID-to-accession-version` had no merged PR that
reverted this.

## Related commits

`a10d89b`, `71ff9b2`.

## Related data / artifacts

`data/01-raw/**/manifest.json` (`protein_uids_sha256`, per-record `protein_uid`).

## Validation

`master` test suite covers `ncbi_snapshot` UID extraction and XML-vs-UID
validation.

## Scientific impact

None directly. It does not change any biological interpretation.

## Data impact

Identifier convention only; no label or partition change.

## Reproducibility impact

The recorded UID list, together with the frozen consolidated XML and the
manifest hashes, identifies the set of records that was actually analysed.
(No claim is made that a UID is invariant across a sequence revision — it is
not; see *Rationale*.)

## Historical reconstruction note

- Directly demonstrated: PR #5 merged into `master`; current code uses UIDs; the
  reversal branch predates PR #5 and contains no reversal.
- Inferred: nothing material.
- Not recoverable: the original rationale for the switch.
