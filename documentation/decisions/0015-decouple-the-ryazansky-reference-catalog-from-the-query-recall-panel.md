# ADR-0015 — Decouple the Ryazansky reference catalog from the query-recall panel

## Date

2026-09-17

## Context

The historical implementation of the Ryazansky reference-catalog builder
(staging branch `(feat)-siepe-ready-project`, commits `b7490e4`/`481d0c2`)
cross-referenced every catalog row against the Phase A query-recall
reference panel (`tests/fixtures/query_recall_reference_set.csv`) to
populate `experimental_anchor`/`experimental_anchor_name` and to decide the
curation outcome for AfAgo, SiAgo, and NgAgo (ADR-0011). That panel belongs
to a separate reintegration unit not yet on `master` when the catalog was
reintegrated (PR #35).

## Decision

The Ryazansky reference catalog (`src/pago_pipeline/clade_reference_catalog.py`,
`scripts/build_clade_reference_catalog.py`) does not depend on the
query-recall panel. The three curation exceptions ADR-0011 already documents
by name — AfAgo and SiAgo quarantined to `UNRESOLVED`; NgAgo's Ryazansky
`LONG_A` label retained with a note — are implemented as a small, explicit,
ADR-0011-cited exception table (`apply_known_architecture_clade_conflicts`)
keyed by accession, not as a live cross-check against the panel. The
historical `experimental_anchor`, `experimental_anchor_name`, and
`proposed_partition` catalog columns, and the panel-dependent
`experimental_anchors_in_table_s1`/`review` summary fields, are not part of
the catalog's schema.

## Rationale

A reference-source catalog (evidence about what Ryazansky Table S1 says) and
a query-coverage audit instrument (evidence about how well a text query
recovers known pAgos) answer different questions and are independently
reintegrable units. Coupling the catalog's schema to the panel would make the
catalog's output depend on which other reintegration units happen to already
be on `master`, and would let a future, unrelated change to the panel
silently change the catalog's schema or curation outcome.

## Consequences

When the query-recall panel is reintegrated, it must not automatically
restore `experimental_anchor`/`experimental_anchor_name` to the Ryazansky
catalog. Any future decision to reconnect the two — for example, to
generalize the AfAgo/SiAgo/NgAgo exception table into a live cross-check
against the full panel — requires its own explicit decision, not an
incidental side effect of reintegrating the panel.

## Related records

- ADR-0011 (documents the AfAgo/SiAgo/NgAgo outcome this decision implements
  without the panel)
- PR #35 — `(feat) Reconcile the Ryazansky pAgo catalog and extract MID-PIWI
  regions`
- commit `a00eb99` (`apply_known_architecture_clade_conflicts`,
  `src/pago_pipeline/clade_reference_catalog.py`)
- historical commits `b7490e4`, `481d0c2` (the coupled implementation this
  decision departs from)
