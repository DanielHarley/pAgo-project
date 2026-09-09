# ADR-0010 — Six APAZ profile HMMs built from the BUILD partition only

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — Phase B (milestone B3) on `(feat)-siepe-ready-project`.

## Epistemic status

`PARTIALLY_SUPPORTED` — the build is deterministic and reproduced byte-for-byte
locally; no sensitivity, specificity, F1 or threshold result exists.

## Type

`TECHNICAL + SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (staging commits
`9e76e78`, `c11732a`, 2026-08-31)

## Context

With the APAZ partition frozen (ADR-0009), the profile models must be built
without seeing calibration or holdout data.

## Problem

An HMM built from any sequence outside BUILD would leak evaluation information
into the classifier.

## Evidence

- Staging commit `9e76e78` — *"(feat) B3 add reproducible APAZ HMM builder"*:
  six HMMs (global + Ia / Ib / IIa / IIb / III) built with
  `pyhmmer.plan7.Builder`, pinned `pyhmmer==0.12.3`, **from the frozen BUILD
  seeds only**; `creation_time` and `command_line` cleared and the payload
  rejected if a `COM` / `DATE` line survives (bit-for-bit reproducibility);
  every written HMM reopened from its closed on-disk artifact and re-validated
  (one model, amino alphabet, expected name, Plan7-valid, non-empty) before its
  SHA-256 is taken. Build output → `data/00-reference/phase_b_builds/
  apaz_hmm_build/` (**gitignored**).
- Staging commit `c11732a` — validation tests.
- `documentation/history/2026-08-31-pago-annotation-ontology-staging.md` — "HMM
  classifier validation protocol": `BUILD` contains only sequences used to
  construct the HMM.
- Local evidence (`tmp/phase_b/apaz_hmms{,2,3,4}/`, gitignored): `apaz_global.hmm`
  is byte-identical (SHA-256 `9b3dd8de…`) across four separate local build
  iterations; `apaz_Ia.hmm` is byte-identical between two of them.

## Previous state

No APAZ profile models.

## Decision

- Build six APAZ HMMs (global + Ia / Ib / IIa / IIb / III) **exclusively from
  the BUILD partition**, deterministically (volatile `COM` / `DATE` lines
  stripped and rejected), with mandatory post-build structural re-validation
  from the on-disk artifact.
- The generated `.hmm` files are **local, gitignored artifacts** — they are
  reproduced from the committed seeds, not frozen in Git.
- The global-vs-subgroup choice and any thresholds are deferred to a later
  calibration step (`pago_domain_profile_validation`) and were **not** run.

## Rationale

Stated in the commit and the ontology protocol: BUILD-only construction keeps
calibration and holdout information out of the classifier; determinism makes the
build a function of the committed seeds; a generated artifact is not a source of
truth and so is not versioned.

## Alternatives considered

`HISTORICAL_EVIDENCE_INSUFFICIENT` for alignment / builder alternatives. A
non-determinism source was found and neutralised: `pyhmmer` writes `COM` / `DATE`
lines from `sys.argv` and a timestamp.

## Consequences

- The APAZ classifier can be calibrated later without leakage.
- No performance number for the APAZ HMMs exists or may be quoted.

## Limitations

Not integrated. Subgroup III has only 12 representatives (BUILD 8), so
III-specific performance cannot be estimated with useful precision.

## Supersedes

Nothing.

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`9e76e78`, `c11732a`.

## Related data / artifacts

`src/pago_pipeline/resources/apaz_seed/apaz_*.sto` (committed seeds, staging);
`data/00-reference/phase_b_builds/apaz_hmm_build/**` (local, gitignored);
`tmp/phase_b/apaz_hmms*/**` (local dev iterations, gitignored).

## Validation

`tests/test_apaz_hmm_build*` (staging): global == union of subgroups, six real
HMMs built and reproducible, no volatile provenance lines, structural validation
of a written artifact. Byte-identical global HMM across four local iterations.

## Scientific impact

Establishes that the APAZ profile instruments are leakage-controlled by
construction. Makes no detection claim.

## Data impact

No candidate-protein label assigned. Generated `.hmm` files are not versioned.

## Reproducibility impact

Positive: bit-for-bit reproducible from the committed BUILD seeds.

## Historical reconstruction note

- Directly demonstrated: the commits, the tests, the ontology protocol, the
  byte-identical local artifacts.
- Inferred: none material.
- Not recoverable: alternative alignment/builder options weighed.
