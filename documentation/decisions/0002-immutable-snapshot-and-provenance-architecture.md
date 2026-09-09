# ADR-0002 — Immutable snapshot, manifest and SHA-256 provenance architecture

## Lifecycle status

`CURRENT`

## Integration status

`INTEGRATED` (PRs #4, #6, #13, #20, #21, #27)

## Epistemic status

`SUPPORTED` — the reuse / invalidation / Windows-safety behaviours are covered
by tests, including a dedicated `windows-snapshot-reuse` CI job on
`windows-latest`. `master` has no branch protection or required status checks,
so nothing enforces that the job stays green before a merge.

## Type

`TECHNICAL`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (built incrementally
2026-03 to 2026-07)

## Context

Every pipeline stage that fetches or derives data needs a reproducible,
auditable output that can be safely reused without re-fetching.

## Problem

Re-running a notebook must not silently reuse a stale artifact, must not
re-download when inputs are unchanged, and must survive interrupted or partial
writes and Windows file-locking.

## Evidence

- PR #4 (`3835181`) / PR #6 (`e58c54e`, `1eb8a7f`, `1520b41`) — reproducible
  snapshot workflow: `snapshots/<timestamp>__<hash>/` immutable directory plus a
  replaceable `latest/`, each with a `manifest.json` recording SHA-256 hashes of
  outputs and of parent snapshots.
- `data/01-raw/README.md` — the raw-data tracking policy (LFS for XML;
  `latest/` tracked; snapshot runs published one at a time with an explicit
  `!` line).
- PR #20 (`c871ecc`) — a `latest/` snapshot is reused only when every declared
  upstream input hash still matches.
- PR #21 (`05ebc81`, `2bf341c`, `b73aa27`) — manifest-only validation before any
  `.npy` array is loaded; a compatible immutable snapshot is republished as
  `latest/` before the old one is touched (avoids Windows
  `PermissionError` / `WinError 32` from memory-mapped `.npy`); artifact hashes
  and SWeeP package-root identity are required for reuse; a `windows-latest` CI
  job runs the reuse regressions.
- PR #13 (`33e7113`) — custom CA / TLS trust without disabling verification.
- PR #27 (`b335fcd`, `9a73c71`, `d5190e6`) — NCBI XML failure controls
  (retry / backoff, circuit breaker, deadlines, resumable batch workspace,
  batch audit cache) across snapshot runs.

## Previous state

Ad hoc outputs without a manifest, hash chain, or reuse contract.

## Decision

Adopt one snapshot contract for every stage: immutable run directory +
replaceable `latest/` + `manifest.json` with output and parent-snapshot
SHA-256s; `reuse_latest_or_create` semantics gated on upstream hash equality;
rollback on partial write; manifest-only checks before loading large arrays.

## Rationale

Verifiable in the code and `data/01-raw/README.md`: hashes make dataset identity
explicit and let reuse be conditional on *everything* that defines a result;
immutable directories give an audit trail; the Windows-specific ordering is
forced by `.npy` memory-mapping behaviour (stated in the PR #21 commit).

## Alternatives considered

`HISTORICAL_EVIDENCE_INSUFFICIENT` — the early PRs (#4, #6) carry no design
document or Issue. The Windows-ordering fix (PR #21) explicitly rejects
"`np.load` the stale `latest/` then replace it".

## Consequences

- Every artifact carries a full provenance chain to its inputs.
- Re-execution with valid snapshots performs no network I/O.
- A changed input, policy hash, or methodology hash invalidates `latest/`.

## Limitations

Storage grows with each published immutable run; which runs to version is a
manual decision (`.gitignore` `!` lines). LFS storage is never reclaimed.

## Supersedes

Unversioned stage outputs.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

PR #4, #6, #13, #20, #21, #27.

## Related commits

`3835181`, `e58c54e`, `1eb8a7f`, `1520b41`, `33e7113`, `c871ecc`, `05ebc81`,
`2bf341c`, `b73aa27`, `b335fcd`, `9a73c71`, `d5190e6`.

## Related data / artifacts

`data/01-raw/**`, `data/03-features/**` layout; every `manifest.json`.

## Validation

`master` test suite + the `windows-snapshot-reuse` CI job on `windows-latest`.
(Note: `master` currently has no branch protection enforcing that job — see
`SCIENTIFIC_STATE.md` §1.2.)

## Scientific impact

None. It is provenance infrastructure.

## Data impact

Defines how data artifacts are stored and reused; no label or partition change.

## Reproducibility impact

Central to it: hash-chained manifests plus immutable runs let a specific
analysed dataset be identified and, given the frozen inputs, reconstructed.

## Historical reconstruction note

- Directly demonstrated: the merged PRs, the committed code, and
  `data/01-raw/README.md`.
- Inferred: that these PRs form one coherent architecture rather than
  independent decisions (grouped here for that reason).
- Not recoverable: rationale for the earliest snapshot PRs beyond what the code
  and README state.
