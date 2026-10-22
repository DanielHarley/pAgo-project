# `query_recall_reference_panel.csv` — curation notes

## 1. Purpose and scope

This panel exists to answer one narrow question: **does a given NCBI text
query recover these known pAgo / PIWI-RE proteins?** It is a
**query-coverage audit instrument**, per ADR-0005.

It is **not**:

- a classifier gold standard
- a detector benchmark
- a detector CALIBRATION set
- a detector FINAL_HOLDOUT
- generic pAgo ground truth
- a PIWI-RE detector-validation set

This restriction is machine-verifiable, not only documented: see
`panel_lock.json`'s `intended_use`/`forbidden_uses`/
`detector_validation_eligible` fields and
`scripts/verify_reference_data.py --scope query_recall`.

This PR adds the versioned panel, its provenance, and its offline integrity
verification only. It does **not** implement query-recall computation
(matching against retrieved records, recall metrics, or a recall snapshot);
that is separate, future work.

## 2. Columns

`accession` (accession.version), `protein_short_name`, `organism`,
`ago_family` (`PAGO` / `PIWI_RE`, per ADR-0012), `pago_clade` (`LONG_A` /
`LONG_B` / `SHORT` / `UNRESOLVED`, meaningful only within `PAGO`, per
ADR-0012 — every `PIWI_RE` row carries `pago_clade = UNRESOLVED`),
`historical_query_recall_clade` (populated only where the historical Phase A
panel's `clade` column diverges from the current `pago_clade`; empty
otherwise — see §4), `sequence_sha256`, `sequence_length`,
`uniprot_accession`, `reference_label_source`, `reference_label_evidence`
(`EXPERIMENTAL` / `LITERATURE_PHYLOGENETIC` / `CURATED_COMPUTATIONAL` /
`DATABASE_ANNOTATION`), `verification_status` (`verified` / `provisional`),
`accession_curation_status` (`RESOLVED` / `UNRESOLVED`), `notes`.

`verification_status` and `accession_curation_status` are deliberately
independent columns, not one merged field: `verification_status` describes
evidence strength for the protein's *biology*; `accession_curation_status`
describes whether *this accession.version* is the settled, canonical
identifier for the reference. A row can be biologically `verified` while its
accession identity is still `UNRESOLVED` — see TtAgo in §4.

Evidence tiers are not interchangeable: `EXPERIMENTAL` means direct
biochemical/structural evidence for the labeled protein; `LITERATURE_PHYLOGENETIC`
means a literature-supported phylogenetic placement; `CURATED_COMPUTATIONAL`
means selection via computational family/profile evidence (see §5 for why
this matters for the PIWI-RE rows specifically); `DATABASE_ANNOTATION` (not
used by any current row) would mean unreviewed database annotation alone.

## 3. `sequence_sha256` — protein-sequence identity column

`sequence_sha256` is `sha256(normalize(protein_sequence))` where `normalize`
= **strip every whitespace character, then uppercase**. Equal hashes
establish sequence equivalence between two records; they do **not**
establish identical genomic locus, identical accession provenance, or
identical biological context (ADR-0005). Matching a retrieved record against
this panel by sequence hash is the responsibility of the future
query-recall-computation PR, not of this resource.

**These 21 hashes are historically derived and frozen, not independently
regenerated in this PR.** The original derivation (documented in the
historical Phase A curation notes) computed each hash from the local
metadata CSV of the first real Phase A retrieval run, which is not present
on `master`. Re-deriving them would require either that local run's
artifacts or a live NCBI/UniProt lookup; neither was performed here (this PR
is fully offline, per its scope). This is stated explicitly rather than
implying a reproducibility this resource does not (yet) have.

## 4. Reconciliation with current project ontology (ADR-0012, ADR-0011)

The historical Phase A panel used a single `clade` column that conflated the
historical Ryazansky-based clade label with today's curated `pago_clade`
semantics. Reconciling the panel against current governance (ADR-0011,
ADR-0012, ADR-0015) required exactly four row-level decisions; every other
row is unchanged from the historical panel.

| Accession | Protein | Historical `clade` | Current `pago_clade` | Reconciliation | Authority |
|---|---|---|---|---|---|
| `WP_005580376.1` | NgAgo | `LONG_B` | `LONG_A` | **Current-label correction.** The historical panel recorded `LONG_B` while citing Ryazansky, but Ryazansky Table S1 records `LONG_A` for this accession; treated as a probable historical recall-panel curation error. The Ryazansky classification is retained. | ADR-0011 (`clade_reference_catalog.py::KNOWN_ARCHITECTURE_CLADE_CONFLICTS["WP_005580376.1"]`, `curation_status = OK`) |
| `WP_010878815.1` | AfAgo | `SHORT` | `UNRESOLVED` (quarantine) | **Current-label correction.** Architecture-based `SHORT` conflicts with the Ryazansky phylogenetic clade `longB_trun`; pending a specific audit against primary AfAgo sources. | ADR-0011 (`KNOWN_ARCHITECTURE_CLADE_CONFLICTS["WP_010878815.1"]`, `curation_status = QUARANTINE`) |
| `WP_012735993.1` | SiAgo | `SHORT` | `UNRESOLVED` (quarantine) | **Current-label correction.** Architecture-based `SHORT` conflicts with the Ryazansky phylogenetic clade `longA_trun`; pending a specific audit against primary SiAgo sources. | ADR-0011 (`KNOWN_ARCHITECTURE_CLADE_CONFLICTS["WP_012735993.1"]`, `curation_status = QUARANTINE`) |
| `WP_011174533.1` | TtAgo | `LONG_A` (unchanged) | `LONG_A` (unchanged) | **Currently unresolved — `accession_curation_status` only, not `verification_status`.** `verification_status` stays historical `verified`: no authoritative decision changes evidence strength for TtAgo's biology. The open question is accession *identity*: this accession is not the crystallographic HB8 reference (Table S1 `YP_145307.1`, later `WP_011229221.1`). Represented by a new `accession_curation_status = UNRESOLVED`, independent of `verification_status`. Not resolved in this PR; `WP_011174533.1` remains the panel accession and is not silently replaced. | `SCIENTIFIC_STATE.md`, Open scientific question 2 |

All other 17 rows are **unchanged**: `pago_clade` equals the historical
`clade` value, `historical_query_recall_clade` is left empty (nothing to
reconcile), `accession_curation_status = RESOLVED`, and no other field
differs from the historical panel.

This is a reconciliation of *current biological labels against current
project governance*, not a re-derivation of biology from data in this PR:
every changed value traces to an ADR or to `SCIENTIFIC_STATE.md`, never to a
new heuristic invented here.

## 5. RsAgo — accession preserved, not swapped

`ABP72561.1` (RsAgo, GenBank, no RefSeq `WP_` record) is kept as the
reference accession. Its byte-identical Swiss-Prot alias `A4WYU7.1` (same
NCBI Identical Protein Group, independently confirmed byte-identical, 777
aa) is **not** substituted in as the reference accession — doing so would
erase the distinction between "the same accession was retrieved" and "the
same sequence was retrieved under a different accession," which is exactly
the distinction the future query-recall computation exists to measure via
the `sequence_sha256` matching tier. `ABP72561.1` stays as-is; it is a
correct accession for RsAgo.

## 6. PIWI-RE set (7 entries) and the circularity boundary

- 1 `EXPERIMENTAL`: `WP_014597637.1` PsPIWI-RE (*Pseudomonas stutzeri* DSM
  4166), Huang et al. 2022.
- 6 `CURATED_COMPUTATIONAL`: NCBI RefSeq/CDD representatives selected from
  the PIWI-RE-specific CDD/Pfam models (`pPIWI_RE_X` PF13111,
  `MID_pPIWI_RE` PF18157, `RNaseH_pPIWI_RE` PF13032), from the family
  established by Burroughs, Iyer & Aravind 2013 (Biol Direct 8:13).

**This restriction is not new to this PR.** The historical Phase A curation
notes already stated, verbatim: *"This is a coverage test of the query, not
a biochemical-activity claim,"* and closed with a section titled *"Not a
validation holdout": "This set may be used for query_reference_recall
(Phase A). It must not become the FINAL_HOLDOUT for a future PIWI-RE HMM:
the CURATED_COMPUTATIONAL members were selected using PIWI-RE profile
models, so scoring a profile-based detector against them would be
circular."* ADR-0014 formalizes this pre-existing constraint as project-wide
governance; it does not correct a historical design error. This panel's
`panel_lock.json` and offline verifier encode the same boundary as a
machine-checkable invariant, for the whole panel, not only its PIWI-RE rows
(§1 above) — the 7 PIWI-RE rows carry an *additional* circularity reason on
top of the panel-wide restriction, not a separate rule.

## 7. Provenance

Recovered from historical branch `(feat)-siepe-ready-project`, commit
`0013c6d24d90070ee297ba6af7d12d125ebd1732` (`(feat) recognize reference
pAgos recovered under a different accession by identical sequence`), the
last commit to touch either historical source file. Historical file
identities at that commit (verified byte-for-byte against the branch tip):

- `tests/fixtures/query_recall_reference_set.csv` —
  `sha256 f16c2a136e362ee2f6483249af34ebf613be22b5f7fafe9297c91ec9be3885d6`
- `tests/fixtures/query_recall_reference_set_curation_notes.md` —
  `sha256 34eead28aad2628c3dc7583cd9d8c52af137a030312af85ab6f84060c0704ffc`

This is secondary provenance for traceability. It is not this resource's
public identity — see `documentation/PROJECT_LEDGER.md` for the current
integration event.

The offline verifier does **not** re-verify these two historical hashes. The
historical branch is not required to be present in a checkout of `master`, so
making verification depend on reading it would break on a shallow clone or a
source export. They are recorded provenance, checked once during curation —
not a current dependency of `scripts/verify_reference_data.py`.

See `panel_lock.json` for the machine-readable identity, counts, and
integrity hashes of this resource. Its declared counts
(`counts_by_ago_family`, `counts_by_evidence_tier`,
`counts_by_verification_status`, `counts_by_accession_curation_status`), its
`schema_fields`, and its `allowed_*_values` sets are cross-checked against
this CSV by the loader and the offline verifier, so they cannot drift from
the panel they describe. `piwi_re_circularity_note` is descriptive prose; the
machine-readable use contract is `intended_use`,
`detector_validation_eligible`, and `forbidden_uses`.
