# `query_recall_reference_set.csv` — curation notes

Purpose: measure whether the notebook 10 text query
`(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism] OR
Archaea[Organism])` recovers *known* pAgo / PIWI-RE proteins, stratified by
MID-PIWI clade (`LONG_A` / `LONG_B` / `SHORT`) and by `ago_family` for
`PIWI_RE`. This is a **coverage test of the query**, not a biochemical-activity
claim — hence `reference_label_evidence` per row.

Columns: `accession` (accession.version), `protein_short_name`, `organism`,
`ago_family` (`PAGO` / `PIWI_RE`), `clade` (`LONG_A` / `LONG_B` / `SHORT` /
`UNRESOLVED`), `sequence_sha256`, `sequence_length`, `uniprot_accession`,
`reference_label_source`, `reference_label_evidence` (`EXPERIMENTAL` /
`LITERATURE_PHYLOGENETIC` / `CURATED_COMPUTATIONAL` / `DATABASE_ANNOTATION`),
`verification_status` (`verified` / `provisional`), `notes`.

## `sequence_sha256` — protein-sequence identity column

`sequence_sha256` is `sha256(normalize(protein_sequence))` where `normalize` =
**strip every whitespace character, then uppercase** (`SEQUENCE_NORMALIZATION =
"strip_all_whitespace_then_uppercase"` in `query_reference_recall.py`). The same
normalization is applied to the retrieved `gbseq__sequence` before hashing, so a
reference is recognised when its exact protein sequence appears in the retrieved
set under *any* accession. `sequence_length` is the residue count after
normalization (informational; not used for matching).

### How the 21 hashes were derived (offline, deterministic, reproducible)

Every hash was computed from the **local metadata CSV of the first real Phase A
run** —
`data/02-intermediate/protein_metadata_csv__annotation_enriched_candidate_set/latest/protein_metadata.csv`
— by looking up each reference accession's `gbseq__sequence` and hashing it with
the normalization above. No online NCBI / IPG call is made, during curation or
during normal pipeline execution. For `ABP72561.1` (RsAgo, GenBank, no `WP_`
RefSeq record and not itself in the retrieved set) the byte-identical Swiss-Prot
alias `A4WYU7.1` — same NCBI Identical Protein Group, independently confirmed
byte-identical, 777 aa — was used as the offline sequence source. The resulting
RsAgo hash is
`cbdb6bb64718c9e8ca78a34ac8445eff1556cb87b5ad687026373ed401c5fb36`.

`compute_query_reference_recall` matches each reference in this order:
`EXACT_ACCESSION_VERSION` -> `SAME_BASE_ACCESSION` -> `SEQUENCE_SHA256` ->
`NONE`, and reports **two** recall readings per stratum: `exact_accession_recall`
(strict-accession methods only) and `retrieval_equivalent_recall` (also counts a
`SEQUENCE_SHA256` hit). When several retrieved records share one sequence hash,
the lexicographically smallest accession is recorded as the representative and
`sequence_match_count` holds the number of colliding records. The matching
methodology is pinned by `matching_strategy_sha256`; changing it invalidates the
`query_reference_recall` snapshot (and nothing else — no downstream artifact
consumes recall).

`PIWI_RE` rows leave `clade = UNRESOLVED`: PIWI-RE is a family, not a pAgo clade.
`compute_query_reference_recall` selects the `PIWI_RE` stratum on
`ago_family == "PIWI_RE"`.

## Resolved: sequence-identity matching (was: exact-accession undercount)

Earlier versions of `compute_query_reference_recall` matched only on
`accession.version` (and its version-stripped form), so the *same protein
sequence* retrieved under a **different accession** scored as a miss even though
the retrieval biologically succeeded. This is now fixed: `SEQUENCE_SHA256` is the
third matching tier (see the `sequence_sha256` section above), and the
`retrieval_equivalent_recall` reading counts such recoveries. The strict
`exact_accession_recall` reading is retained alongside it. The RsAgo case below
is what motivated the change; it is now recovered via `A4WYU7.1`.

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
identical sequence). The reported miss was a benchmark artifact of exact-accession
matching, not a query gap.

**Recall figures** for this panel with the sequence-identity matcher:
`retrieval_equivalent_recall` `overall = 21/21 = 1.000`, `LONG_B = 2/2 = 1.000`;
the strict `exact_accession_recall` stays `overall = 20/21 = 0.952`,
`LONG_B = 1/2 = 0.500` — RsAgo's `ABP72561.1` is genuinely absent from the
retrieved accessions, and that reading is reported honestly rather than hidden.
`ABP72561.1` stays as-is in the CSV — it is a correct accession for RsAgo.

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
