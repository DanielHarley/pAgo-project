# ADR-0012 — Separate Argonaute family identity from pAgo phylogenetic clade

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09`

## Context

pAgo annotation mixes several distinct concepts: whether a protein is an
Argonaute at all, which family it belongs to, which phylogenetic clade applies
within PAGO, its domain architecture, and its catalytic status.

Collapsing those concepts into one convenient label can turn weak or indirect
evidence into an unsupported biological claim.

## Decision

Keep the annotation dimensions separate:

- `ago_family ∈ {PAGO, PIWI_RE, UNRESOLVED}`
- `pago_clade ∈ {LONG_A, LONG_B, SHORT, UNRESOLVED}` and only applies within
  the PAGO family
- PIWI-RE is a family, never a pAgo clade
- PIWI-RE references therefore carry `pago_clade = UNRESOLVED`
- domain architecture does not by itself determine clade
- a triage-only competitive-HMM signal does not automatically establish
  Argonaute detection
- absence of one Pfam PIWI hit does not automatically establish a biological
  negative
- catalytic-site annotation remains separate from family, clade, architecture,
  and experimentally demonstrated catalytic activity

The historical placement protocol described in the staging ontology should not
be treated as experimentally or computationally validated merely because it was
written down.

## Evidence

- the historical ontology explicitly defines PIWI-RE as a family-level
  assignment rather than a pAgo clade
- historical commit `98e66c6` removes `PIWI_RE` from allowed clade values and
  requires PIWI-RE rows to use `clade = UNRESOLVED`
- the Phase A audit records the same family/clade separation
- Burroughs et al. 2013 describes the PIWI-RE family, while Ryazansky et al.
  2018 provides the LONG_A/LONG_B/SHORT pAgo clade context

## Rationale

Family identity, phylogenetic clade, architecture, detection evidence, and
catalytic activity answer different biological questions. Keeping them separate
prevents one evidence source from silently being promoted into another kind of
claim.

## Consequences

- query-recall stratification for PIWI-RE uses `ago_family`, not `pago_clade`
- later placement/evaluation must distinguish a family-level error from an
  intra-pAgo clade error
- architecture-only evidence cannot be used to force a clade label

## Supersedes

The early fixture convention that represented `PIWI_RE` as a value of
`pago_clade`.

## Limitations

The historical staging ontology also describes a future phylogenetic-placement
procedure and likelihood-weight thresholds. Writing that procedure did not by
itself validate the reference tree, thresholds, or placement performance.

## Validation

Historical tests reject `PIWI_RE` as a clade value and assert that PIWI-RE rows
carry an unresolved pAgo clade.

## Related records

- historical `documentation/history/2026-08-31-pago-annotation-ontology-staging.md`
- historical commit `98e66c6` and related recall-panel commits
- `documentation/history/2026-08-30-phase-a-audit.md`
- Burroughs et al. 2013 and Ryazansky et al. 2018

## Historical reconstruction note

- Directly demonstrated: historical ontology, test change, audit, and cited
  primary literature.
- Not recoverable: the complete historical deliberation over every ontology
  vocabulary choice.
