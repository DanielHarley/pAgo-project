# ADR-0010 — Build APAZ profile HMMs from BUILD references only

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commits
`9e76e78`, `c11732a`, 2026-08-31)

## Context

Once APAZ references are partitioned, profile construction must not use
calibration or final-evaluation sequences.

## Decision

Build six APAZ profile HMMs — global plus Ia, Ib, IIa, IIb, and III —
exclusively from BUILD references.

Make the build deterministic by removing volatile `COM` and `DATE` provenance
lines, rejecting artifacts that retain them, and reopening the written HMM for
structural validation before recording its SHA-256 identity.

Treat generated `.hmm` files as reproducible build artifacts rather than as the
versioned source of truth. Defer model strategy, thresholds, and performance
claims to later calibration/evaluation work.

## Evidence

- historical commit `9e76e78` implements the six BUILD-only HMMs with
  `pyhmmer==0.12.3`, volatile-line removal, and post-write structural checks
- historical commit `c11732a` adds validation tests
- local historical build iterations reproduced `apaz_global.hmm` byte-for-byte
  across four builds and `apaz_Ia.hmm` across multiple builds
- the historical ontology explicitly states that BUILD contains only sequences
  used to construct the HMM

## Rationale

BUILD-only construction keeps calibration and final-holdout information out of
the classifier. Deterministic generation makes an HMM a reproducible function of
its declared BUILD references rather than an opaque local artifact.

## Consequences

- later calibration can compare model strategies and thresholds without having
  leaked evaluation sequences into HMM construction
- deterministic construction and structural validity do not establish
  sensitivity, specificity, F1, or a valid decision threshold
- subgroup III has comparatively few BUILD references, so subgroup-specific
  performance must not be assumed from construction alone

## Validation

Historical tests check the six real HMMs, global/subgroup reference consistency,
absence of volatile provenance lines, deterministic output, and structural
validity of written artifacts.

## Related records

- historical commits `9e76e78`, `c11732a`
- historical APAZ BUILD seed alignments
- local historical generated HMM iterations
- historical roadmap reference: B3

## Historical reconstruction note

- Directly demonstrated: historical commits, tests, ontology protocol, and
  byte-identical local artifacts.
- Not recoverable: the full deliberation over alternative alignment/build
  strategies.
