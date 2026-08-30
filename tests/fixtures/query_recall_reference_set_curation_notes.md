# `query_recall_reference_set.csv` — curation notes

Purpose: measure whether the notebook 10 text query
`(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR
Archaea[Organism])` recovers *known* pAgo / PIWI-RE proteins, stratified by
MID-PIWI clade (`LONG_A` / `LONG_B` / `SHORT`) and by `ago_family` for
`PIWI_RE`. This is a **coverage test of the query**, not a biochemical-activity
claim — hence `reference_label_evidence` per row.

Columns: `accession` (accession.version), `protein_short_name`, `organism`,
`ago_family` (`PAGO` / `PIWI_RE`), `clade` (`LONG_A` / `LONG_B` / `SHORT` /
`UNRESOLVED`), `uniprot_accession`, `reference_label_source`,
`reference_label_evidence` (`EXPERIMENTAL` / `LITERATURE_PHYLOGENETIC` /
`CURATED_COMPUTATIONAL` / `DATABASE_ANNOTATION`), `verification_status`
(`verified` / `provisional`), `notes`.

`PIWI_RE` rows leave `clade = UNRESOLVED`: PIWI-RE is a family, not a pAgo clade.
`compute_query_reference_recall` selects the `PIWI_RE` stratum on
`ago_family == "PIWI_RE"`.

## PIWI-RE set (7 entries)

- 1 EXPERIMENTAL: `WP_014597637.1` PsPIWI-RE (*Pseudomonas stutzeri* DSM 4166),
  Huang et al. 2022.
- 6 CURATED_COMPUTATIONAL: NCBI RefSeq/CDD representatives that appear in two or
  three of the PIWI-RE-specific CDD/Pfam models
  (`pPIWI_RE_X` PF13111, `MID_pPIWI_RE` PF18157, `RNaseH_pPIWI_RE` PF13032),
  from the family established by Burroughs, Iyer & Aravind 2013 (Biol Direct
  8:13).

## Historical Burroughs 2013 identifiers — pending / excluded

The founding PIWI-RE paper cites GenInfo Identifier (GI) numbers. GI was retired
as a primary identifier, so each must be resolved to a current accession.version.

| GI | organism | resolution | outcome |
|----|----------|-----------|---------|
| 269125748 | *Thermomonospora curvata* | GI -> YP_003299118.1 -> **WP_012851864.1** (live) | **included** (CURATED_COMPUTATIONAL) |
| 228927677 | *Bacillus thuringiensis* serovar pondicheriensis | ZP_04090728.1 — dead record, **removed, no replacement** | **excluded** — no resolvable current accession |
| 119855142 | *Mycobacterium* sp. KMS | GI -> YP_935747.1 -> WP_011767947.1 (`argonaute/piwi family protein`, 477 aa) but WP_011767947.1 is **suppressed / no longer annotated on any genome** | **excluded** — not an unambiguous live record; also short (477 aa) |

Do **not** add `GI 158336201` (*Acaryochloris marina*): per Burroughs et al. 2013
it is the associated restriction endonuclease (REase) neighbouring the PIWI-RE
gene, not the PIWI-RE protein itself.

## Not a validation holdout

This set may be used for `query_reference_recall` (Phase A). It must **not**
become the `FINAL_HOLDOUT` for a future PIWI-RE HMM: the CURATED_COMPUTATIONAL
members were selected *using* PIWI-RE profile models, so scoring a profile-based
detector against them would be circular. That circularity does not apply to the
Phase A question (does an NCBI *text* query recover them?).
