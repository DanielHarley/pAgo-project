# `ryazansky_s1_catalog.csv` — audit notes (Phase B, sub-step B4.2)

Deterministic catalog of the **1010 prokaryotic Argonaute (pAgo) proteins** listed
in Ryazansky, Kulbachinskiy & Aravin 2018, *mBio* 9(6):e01935-18
(doi 10.1128/mBio.01935-18, PMC6299218.1, CC BY), built by
`scripts/build_clade_reference_catalog.py`.

This is an **intermediate working table**, not a frozen reference. It contains
no partition assignments (`proposed_partition` is `UNASSIGNED` on every row), no
alignment, no split groups, no HMM. `clade_partitions.csv` is deliberately not
created.

## Provenance

Sources materialized under `source_data/` (all md5-verified against the official
PMC Open Access metadata, which is byte-identical to the copy frozen in B2 at
`apaz_seed/source_data/PMC6299218.1.metadata.json`, sha256
`89aa43f27ebb729d388abfe6e6098f8e02d7dc9095cb63b4f6e5e3bc3f8385cb`):

| file | role | md5 | bytes |
|---|---|---|---:|
| `mbo006184236st1.xls` | Table S1 — the 1010 pAgos, clade label (`Ago_type`), PAZ/MID/PIWI coordinates | `36c9e5b00dfd327b856bdc69de9d70f0` | 309 760 |
| `mbo006184236sd1.txt` | Data Set S1 — MID\* 5′-motif fragment alignment, 75 non-canonical long pAgos (diagnostic; **not** the tree alignment) | `4a151e924f85810a34f2ad5cbcc3278f` | 17 848 |
| `mbo006184236sd2.txt` | Data Set S2 — PAZ / PAZ\* domain alignment of long pAgos (diagnostic; **not** the tree alignment) | `3a160c67945eff04ad161662ae781a77` | 143 035 |
| `ryazansky_s1_sequences.fasta` | frozen full protein sequences for the 1010 accessions + 7 replacement accessions (fetched by `accession.version` from NCBI protein) | — | — |
| `ryazansky_s1_ncbi_resolution.json` | NCBI `esummary` record status snapshot for the 1010 accessions | — | — |

**Not materialized** (recorded here for provenance only — genomic-context data,
not a clade or sequence resource):

| file | md5 | sha256 |
|---|---|---|
| `mbo006184236st2.xlsx` (Table S2, operon-neighbour orthogroups) | `ffc634d19468f4ab872fc71cb35ed649` | `f8731f8012ad316a106b17aaab4b607876750d03223a0fe83fc2d4d8f48aa59d` |
| `mbo006184236st3.xls` (Table S3, operon size / orthogroup frequency) | `fd61523055f660ebcd75bc9cedefb9ed` | `c0976419a6917068c7fe6fa10538d47df91de4fc1fd1bdddcb8febba8d0776cb` |

`PMC6299218.1.txt` (article full text, sha256
`789d073026bf3c31511bfe7ea8effb49b7617d1a3ef88b03c53034469cb5d15b`) was read as
the interpretive source for the notes below; it is not materialized.

## What Table S1 is (from the article Methods / Results)

- **1010 pAgos** found by PSI-BLAST (query = PIWI-MID domains of known pAgos)
  over RefSeq, in 1385 genomes of 1248 strains. Table S1 is this full set and
  **precedes** the redundancy reduction.
- The **721 nonredundant representatives** (UCLUST 4.2, 90% identity, longest
  sequence per cluster) are what the phylogenetic tree (Fig. 1B) was built on:
  MID-PIWI regions → MAFFT (`-ep 0 --genafpair --maxiterate 1000`) → trimAl
  (drop columns with ≥0.5 gaps) → FastTree (WAG, discrete gamma 20 categories),
  bootstrap n=100. **The 721-set accession list, the MID-PIWI MSA, and the tree
  itself are not published** — only Fig. 1B (image) and Fig. S2 (collapsed
  triangles).
- The tree resolved **three clades** — short, long-A, long-B — plus **two
  proteins that could not be classified** (high divergence): from
  *Thermoproteus uzoniensis* and *Vulcanisaeta moutnovskia*
  (`Ago_type = unkn`).
- `Ago_type` is the phylogenetic clade label. `_trun` = truncated / partial
  variant of that clade (usually missing the N-terminal / PAZ region).
- Data Set S1 = MID\* motif fragment for 75 long pAgos with non-canonical
  5′-binding residues (Fig. 3B logo). Data Set S2 = PAZ/PAZ\* R3-region
  alignment of long pAgos (Fig. 4). Neither is the MID-PIWI tree alignment.

## TtAgo / RsAgo (Correction 1)

The earlier claim that TtAgo and RsAgo are "absent from Table S1 because removed
by the 90% filter" is **withdrawn**. Table S1 is the pre-reduction 1010 set.
Investigated by exact accession and by sequence sha256:

- **RsAgo is present** in Table S1 as **`WP_011910606.1`** (*Rhodobacter
  sphaeroides*, `Ago_type = longB`, 777 aa). Its sequence is **byte-identical**
  (sha256 `cbdb6bb64718c9e8ca78a34ac8445eff1556cb87b5ad687026373ed401c5fb36`) to
  the `query_recall_reference_set.csv` accession `ABP72561.1` and to Swiss-Prot
  `A4WYU7.1`. The recall-set note "no RefSeq WP_ record exists for this sequence"
  is incorrect; `WP_011910606.1` is exactly that record. Cause: the recall set
  curated the strain-specific GenBank accession (*C. sphaeroides* ATCC 17025)
  and did not locate the multispecies WP_ collapse.
- **The canonical TtAgo** (*T. thermophilus* HB8, gene TTHB068, the
  crystallography workhorse, UniProt Q746M7) **is present** in Table S1 as
  **`YP_145307.1`** (`Ago_type = longA`, 685 aa). `YP_145307.1` is now a dead
  accession → NCBI replacement `WP_011229221.1` (byte-identical, sha256
  `535945486b3530a4925ded522620ed1ddc91c92a4bce85e3c583b93410d173fd`);
  `WP_011229221.1` is **not** separately a Table S1 row.
- The recall-set TtAgo accession **`WP_011174533.1`** (also 685 aa, "longA") is
  **not** in Table S1 and is **not byte-identical** to the HB8 sequence
  (sha256 `1a654ce2…`). It is a divergent *T. thermophilus* Ago from a different
  isolate. **Curation decision pending**: which accession represents the
  crystallized TtAgo for the reference set.
- Which specific redundant homolog stood in for RsAgo / TtAgo in the **721-set
  tree** is **not reconstructable** from the published supplements.

## Coordinate convention (proven, then asserted before slicing)

The PAZ/MID/PIWI columns of Table S1 are **1-based, inclusive, with MID and PIWI
contiguous** (MID then PIWI). Evidence over all 1010 rows:

- `PIWI_end > protein_length`: **0** (never exceeds);
- `MID_end + 1 == PIWI_start`: **1002 / 1002** (every row with both);
- `MID_end == PIWI_start` (a 0-based half-open signature): **0**;
- `PAZ_end < MID_start`: **393 / 393**.

MID-PIWI region = `sequence[MID_start-1 : PIWI_end]` (Python slice). The builder
raises if this convention cannot be re-proven, and
`tests/test_clade_reference_catalog.py` guards the slice against an off-by-one
regression.

## Catalog contents & QC (rebuild sha256 `dde09c035919ad6ebd23199200e13d78c73afa5ef842193d9c39114cff95d631`)

- **1010 rows**, one per Table S1 entry. **1010 / 1010 sequences recovered.**
- NCBI record status: **731 LIVE · 272 SUPPRESSED · 7 DEAD_REPLACED**.
  Suppression is concentrated in `_trun` rows (32–67%) vs ~23–27% for canonical.
  A suppressed record's sequence is still retrievable; suppressed rows are kept
  with `ncbi_record_status = SUPPRESSED`, never auto-excluded.
- **7 dead accessions** (`NP_*` / `YP_*`), each with an NCBI replacement; all 7
  replacements are **byte-identical** to the source (`replacement_equivalent =
  True`). For **6 of the 7**, the replacement `WP_*` is **also a separate Table
  S1 row** → those 6 proteins are **double-listed** (obsolete + current
  accession). Effective distinct-protein count ≈ **1004**
  (`n_unique_full_sequence_sha256 = 1004`).
- **MID-PIWI extracted for 1002 / 1010**; 8 `MISSING_MID` (rows with no MID
  coordinate in Table S1 — all short-truncated fragments).
- **984 unique MID-PIWI regions** among the 1002 extracted → 18 duplicate
  groups: 6 are the obsolete-accession double-listings, **12 are genuine**
  (distinct proteins / distinct full sequences that share a byte-identical
  MID-PIWI region).
- `length_matches_source`: **1010 / 1010** — `accession.version` is immutable;
  no sequence drift.
- 30 `_trun` rows nonetheless have a **complete** MID-PIWI region (usable
  despite the truncated N-terminus); 8 `_trun` rows do not.
- QC checks with **0 hits**: duplicate accession string, coordinate out of
  range, `MID_start > MID_end`, `PIWI_end > current_length`, non-equivalent
  replacement, MID|PIWI gap/overlap, length drift.

### Clade counts

| `source_phylogenetic_clade` | count | `curated_pago_clade` | count |
|---|---:|---|---:|
| LONG_A | 226 | LONG_A | 225 |
| LONG_B | 199 | LONG_B | 198 |
| SHORT | 583 | SHORT | 583 |
| UNRESOLVED | 2 | UNRESOLVED | 4 |

`curated_pago_clade` differs from `source_phylogenetic_clade` only for the 2
`unkn` rows and the 2 QUARANTINE rows (all → `UNRESOLVED`). Canonical 972 /
truncated 38.

## Conflicts and quarantine (Correction 2 — architecture and clade kept separate)

`source_ago_type` and `source_phylogenetic_clade` are Ryazansky's phylogenetic
evidence and are always preserved. MID-PIWI-only architecture does **not** by
itself define SHORT.

| protein | recall-set `clade` | Ryazansky `Ago_type` | `curation_status` | handling |
|---|---|---|---|---|
| `WP_005580376.1` NgAgo | LONG_B (cites Ryazansky) | `longA` | **OK** | Recall set contradicts the source it cites → **probable recall-set error**. Ryazansky classification retained (`LONG_A`). The Phase A recall artifact is **not** modified. |
| `WP_010878815.1` AfAgo | SHORT (Parker 2004 / Ma 2005) | `longB_trun` | **QUARANTINE** | Architecture is MID-PIWI-only ("half Argonaute"); Ryazansky places it as truncated long-B. Genuine ambiguity. `curated_pago_clade = UNRESOLVED` until a specific audit against the primary AfAgo literature. Not forced either way. |
| `WP_012735993.1` SiAgo (strain M.16.4) | SHORT | `longA_trun` | **QUARANTINE** | Recall set also has the wrong strain (M.16.4 ≠ characterized REY15A). `curated_pago_clade = UNRESOLVED` pending audit. Not forced. |
| `WP_010942012.1` GsAgo | SHORT | `short` | OK | agree. |

The **2 `unkn` proteins** (`WP_013604390.1` *V. moutnovskia*, live;
`WP_052886194.1` *T. uzoniensis*, suppressed) stay `ago_family = PAGO`,
`curated_pago_clade = UNRESOLVED`, and are never LONG_A/LONG_B/SHORT seeds.

## Experimental anchors present in Table S1 (matched by sequence sha256)

`AaAgo`, `AfAgo`, `GsAgo`, `MhAgo`, `MjAgo`, `MpAgo`, `NgAgo`, `PfAgo`, `RsAgo`,
`SeAgo`, `SiAgo` (11 of the 14 recall-set pAgo anchors). Absent from Table S1:
**TtAgo** (recall-set accession is a non-matching isolate — see Correction 1
above; the HB8 TtAgo is present under `YP_145307.1`), **TpsAgo**
(`WP_060384876.1`), **MapAgo** (`WP_109649955.1`, Koopal 2022 — post-dates the
2018 study). `experimental_anchor` is a flag only and does **not** set the
partition.

## Potential-exclusion inventory (no exclusions applied)

| reason | count |
|---|---:|
| SUPPRESSED | 272 |
| SUPPRESSED and truncated | 16 |
| DEAD_REPLACED | 7 |
| obsolete accession double-listed (same protein twice) | 6 |
| missing MID coordinate | 8 |
| MID-PIWI extraction failed | 8 |
| `curated_pago_clade == UNRESOLVED` | 4 |
| QUARANTINE | 2 |

Nothing is excluded in B4.2. These are candidates for a later, explicit
inclusion policy.
