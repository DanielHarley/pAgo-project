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

## Known limitation: exact-accession matching can undercount recall

`compute_query_reference_recall` matches a reference to the retrieved set on
`accession.version` (and its version-stripped form). If the *same protein
sequence* is retrieved under a **different accession** than the one in this CSV,
it is scored as a miss even though the retrieval biologically succeeded.

### Investigated case: `ABP72561.1` / RsAgo (first real Phase A run, 52,473 records)

`ABP72561.1` (GenBank, `hypothetical protein Rsph17025_3694 (plasmid)`,
*Cereibacter sphaeroides* ATCC 17025, 777 aa; UniProt reviewed A4WYU7 "Protein
argonaute", gene `Rsph17025_3694`) was reported as the single miss
(`long_b_reference_recall = 0.5`, `overall = 0.952`).

Mechanical investigation (not an accession swap):

- NCBI Identical Protein Group for `ABP72561.1` contains `A4WYU7.1` (Swiss-Prot)
  and `XLG71013.1` (patent) — no RefSeq `WP_` record exists for this sequence.
- `A4WYU7.1` **is present** in the 52,473 retrieved records (protein_uid
  `2500461169`) and its `gbseq__sequence` is **byte-identical** to `ABP72561.1`
  (both 777 aa, `sha256 cbdb6bb64718c9e8ca78a34ac8445eff1556cb87b5ad687026373ed401c5fb36`).
- RsAgo is also present as PDB chains `5AWH_A/B`, `6D8A`, `6D8F`, `6D8P`, `6D92`
  (*C. sphaeroides* ATCC 17025).

**Conclusion: A** — RsAgo was recovered by the query (via `A4WYU7.1`, same IPG,
identical sequence). The reported miss is a benchmark artifact of exact-accession
matching, not a query gap.

**Consequence for the recall figures** from that run: the biologically correct
recall on this panel is `overall = 21/21 = 1.000` and `LONG_B = 2/2 = 1.000`;
the reported `0.952` / `0.500` reflect the matcher, not the retrieval.

Recommended fix (pending approval, not applied): match a reference by protein
sequence identity (e.g. a `sequence_sha256` column, or NCBI IPG resolution),
falling back to accession. `ABP72561.1` should stay as-is in the CSV — it is a
correct accession for RsAgo.

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
