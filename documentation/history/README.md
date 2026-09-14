# Historical records

Historical records inherit the documentation principles in
[`../README.md`](../README.md).

This directory preserves point-in-time evidence: snapshots of what the project
looked like, what a branch contained, or what a staging document said at a
specific moment.

Nothing in this directory should be treated as a description of current project
state merely because it is detailed or searchable.

For current scientific interpretation, read
[`../SCIENTIFIC_STATE.md`](../SCIENTIFIC_STATE.md) and the relevant current
ADRs.

Historical bodies are preserved without semantic rewriting. If terminology,
interpretation, or decisions change later, record the correction in current
documentation rather than editing the old body to make history look cleaner.

This matters especially for agent/tool search: a search result may return a
fragment from the middle of a historical file without its header. Always check
the file path and current documentation before treating a historical fragment
as current authority.

| File | What it is |
| --- | --- |
| `2026-08-30-phase-a-audit.md` | Read-only technical audit of the Phase A staging work as it existed on 2026-08-30. |
| `2026-08-31-pago-annotation-ontology-staging.md` | Annotation-ontology draft from staging work, preserved as historical provenance rather than rewritten as a current ontology. |
