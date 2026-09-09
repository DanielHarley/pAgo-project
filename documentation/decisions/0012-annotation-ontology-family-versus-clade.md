# ADR-0012 — Annotation ontology: `ago_family` vs `pago_clade`; PIWI-RE is a family, never a clade

## Lifecycle status

`CURRENT`

## Integration status

`NOT_INTEGRATED` — the ontology document is staging (see
`documentation/history/2026-08-31-pago-annotation-ontology-staging.md`).

## Epistemic status

`SUPPORTED` for the semantic separation. The future phylogenetic-placement
protocol referenced by the ontology is `PROPOSED` and not executed.

## Type

`SCIENTIFIC`

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09`

## Context

pAgo annotations mix several orthogonal dimensions: whether a protein is an
Argonaute at all, which biological family it belongs to, which phylogenetic
clade (within the PAGO family), its domain architecture, and its catalytic
status.

## Problem

A computationally plausible label (a text hit, an architecture, a competitive-
HMM score) can be silently promoted to biological identity, phylogenetic clade,
or catalytic activity.

## Evidence

- `documentation/history/2026-08-31-pago-annotation-ontology-staging.md`:
  `ago_family` has the vocabulary `PAGO`, `PIWI_RE`, `UNRESOLVED`; **"PIWI-RE is
  a family-level assignment. It is not a pAgo clade."** `pago_clade` has
  `LONG_A`, `LONG_B`, `SHORT`, `UNRESOLVED`, defined only within PAGO and only
  by phylogenetic placement against the frozen reference tree plus three
  calibrated likelihood-weight-ratio thresholds. `catalytic_site_status` is
  fixed to `UNKNOWN` in pipeline version 1. "An architecture such as `MID_PIWI`
  does not imply `pago_clade = SHORT`."
- Staging commit `98e66c6` — *"(test) drop PIWI_RE from allowed clade values;
  assert PIWI-RE rows are clade=UNRESOLVED"*: the test no longer accepts
  `"PIWI_RE"` as a `clade` value and asserts every `ago_family == "PIWI_RE"`
  row has `clade == "UNRESOLVED"`.
- `documentation/history/2026-08-30-phase-a-audit.md` D5 — separate
  `ago_family` (`PIWI_RE`) from `pago_clade`.
- Primary literature: Burroughs, Iyer & Aravind 2013 (*Biology Direct* 8:13,
  PMC3702460) — the PIWI-RE family; Ryazansky et al. 2018 — the LONG_A /
  LONG_B / SHORT MID-PIWI clades. PIWI-RE is absent from Ryazansky.

## Previous state

Early fixtures / tests treated `PIWI_RE` as a value of `pago_clade`.

## Decision

- `ago_family ∈ {PAGO, PIWI_RE, UNRESOLVED}` (biological identity).
- `pago_clade ∈ {LONG_A, LONG_B, SHORT, UNRESOLVED}` (phylogenetic property,
  only within PAGO, only by placement + three thresholds).
- **PIWI-RE is a family, never a pAgo clade.** Every PIWI-RE row carries
  `pago_clade = UNRESOLVED`.
- Detection is a separate epistemic axis; a triage-only MID-PIWI / competitive-
  HMM signal never becomes `DETECTED`; absence of a Pfam PIWI hit never sets
  `NOT_DETECTED`.
- Domain architecture does not imply a clade. `catalytic_site_status = UNKNOWN`
  in v1.

## Rationale

Stated in the ontology and the audit: PIWI-RE is a divergent *family*
(Burroughs 2013), not a clade inside the long/short pAgo MID-PIWI tree; keeping
the axes disjoint prevents a convenient collapse of ambiguous evidence.

## Alternatives considered

`clade = "PIWI_RE"` (the earlier fixture convention, rejected by `98e66c6`);
ignoring PIWI-RE entirely (rejected — it is a real family in the recall panel).
`HISTORICAL_EVIDENCE_INSUFFICIENT` for other vocabularies.

## Consequences

- The recall panel selects the PIWI-RE stratum on `ago_family == "PIWI_RE"`,
  not on `clade`.
- The placement protocol must distinguish a family error (SHORT → PIWI-RE) from
  an intra-pAgo error (SHORT → LONG_B).

## Limitations

Not integrated. The placement protocol (three thresholds, reference tree) is
`PROPOSED`; no placement has been run. The ontology document is preserved under
`documentation/history/` and is **not** the authoritative current ontology —
that document will be written when the corresponding work is reintegrated.

## Supersedes

`PIWI_RE` treated as a value of `pago_clade` (early fixtures / tests).

## Superseded by

None.

## Related Issue

`No historical Issue found`.

## Related PR

`No historical PR found`.

## Related commits

`98e66c6`, `0fd283a`, `21793bc`.

## Related data / artifacts

`documentation/history/2026-08-31-pago-annotation-ontology-staging.md`;
`tests/fixtures/query_recall_reference_set.csv` (staging).

## Validation

`test_committed_reference_set_is_well_formed` (staging) — rejects `"PIWI_RE"`
as a `clade` value; asserts PIWI-RE rows are `clade == "UNRESOLVED"`.

## Scientific impact

Core of the annotation semantics: identity, clade, architecture and catalysis
are kept distinct; PIWI-RE is fixed as a family.

## Data impact

Defines the label vocabularies. No candidate-protein label assigned yet.

## Reproducibility impact

None directly — it is a semantic contract.

## Historical reconstruction note

- Directly demonstrated: the ontology text, the `98e66c6` test change, the
  audit, the primary papers.
- Inferred: none material.
- Not recoverable: the full deliberation over the vocabularies.
