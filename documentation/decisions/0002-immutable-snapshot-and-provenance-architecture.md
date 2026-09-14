# ADR-0002 — Use immutable snapshots, manifests, and SHA-256 provenance

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (built incrementally
2026-03 to 2026-07)

## Context

Pipeline stages fetch or derive data that may later be reused. Re-running a
notebook must not silently reuse stale outputs, partially publish a failed run,
or lose the provenance needed to identify what inputs produced a result.

## Decision

Use one snapshot contract across pipeline stages:

- immutable `snapshots/<timestamp>__<hash>/` run directories
- a replaceable `latest/`
- `manifest.json` with output and parent-input SHA-256 identities
- `reuse_latest_or_create` semantics gated on declared input identity
- cleanup/rollback when publication fails
- manifest-level compatibility checks before loading large array artifacts

For Windows-safe reuse, compatible immutable snapshots are republished as
`latest/` before stale memory-mapped `.npy` artifacts are opened or replaced.

## Evidence

- PRs #4 and #6 introduced reproducible snapshot storage, manifests, and parent
  hashes
- PR #20 made reuse conditional on upstream input identity
- PR #21 added manifest-first validation, stricter artifact identity checks, and
  Windows-safe publication ordering
- PR #13 added configurable CA/TLS trust without disabling certificate
  verification
- PR #27 added retry/backoff, deadlines, circuit breaking, resumable XML batch
  workspaces, and batch-audit caching
- `data/01-raw/README.md` documents the raw-data tracking policy

## Rationale

Hashes make artifact identity explicit. Immutable run directories preserve an
audit trail. Conditional reuse avoids unnecessary work while preventing stale
results from being silently accepted. The Windows-specific ordering avoids
replacing files that are still held open through memory maps.

## Consequences

- persisted artifacts carry a provenance chain to their inputs
- unchanged inputs can reuse a valid snapshot without repeating network work
- changes to declared inputs, policies, or methodology identities invalidate
  reuse
- immutable published runs accumulate storage and require an explicit retention
  policy

## Validation

The repository test suite covers snapshot creation/reuse, integrity checks, and
Windows reuse behaviour, including the dedicated Windows CI path introduced for
the regression.

## Supersedes

Ad hoc mutable stage outputs without an explicit manifest/provenance contract.

## Related records

- PRs #4, #6, #13, #20, #21, #27
- commits `3835181`, `e58c54e`, `1eb8a7f`, `1520b41`, `33e7113`, `c871ecc`,
  `05ebc81`, `2bf341c`, `b73aa27`, `b335fcd`, `9a73c71`, `d5190e6`
- `data/01-raw/README.md`

## Historical reconstruction note

- Directly demonstrated: merged PRs, committed code, tests, and raw-data README.
- Inferred: these incremental changes form one coherent snapshot/provenance
  architecture.
- Not recoverable: the full rationale for the earliest snapshot changes beyond
  what the code and README preserve.
