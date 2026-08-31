# APAZ reference curation notes

## 1. Primary source

Ryazansky S, Hilbig JK, Kolchenko M, Aravind L, Rechavi O.
"The Expanded Universe of Prokaryotic Argonaute Proteins."
mBio 9(6):e01935-18 (2018). DOI 10.1128/mBio.01935-18. PMC6299218. PMID 30563906.
License **CC BY**.

The APAZ alignment is Data Set S3 (`source_data/mbo006184236sd3.txt`),
md5 `ac30768ca799feb85165cea9d11aa3e2`, sha256
`a591e80ba5ea79eec0d7f8736122b8cfb67592496e0033196d8381766c151247`, obtained from
the PMC OA Open Data bucket and pinned in `seeds_lock.json`. No network access is
made during pipeline execution or validation.

## 2. Groups Ia–III

Ryazansky et al. define the **five APAZ groups by a phylogenetic analysis of the
APAZ domain itself**, then describe the architecture of the proteins that host
each group. The architectures are **predominant tendencies, with exceptions** —
they are not the definition of the group and must not be treated as
detection/family/clade state.

| group | described host architecture (predominant) | original count | S3 representatives |
|---|---|---|---|
| Ia | SIR2–APAZ–MID–PIWI | 117 | 83 |
| Ib | SIR2–APAZ | 150 | 109 |
| IIa | TIR–APAZ | 122 | 98 |
| IIb | mixed: 112 DUF4365–APAZ, 13 SIR2–APAZ, 100 with no additional domain | 225 | 179 |
| III | solo APAZ in Archaea | 18 | 12 |
| **total** | | **632** | **481** |

(Group / architecture counts are from the paper; cross-check against Table S1 when
this document is finalised.)

Do **not** write "every IIb protein is DUF4365–APAZ": IIb is a mixed set (see the
table).

## 3. The 481 representatives — `source_redundancy_representative` vs `split_group`

Ryazansky et al. reduced the 632-sequence APAZ set with **UCLUST 4.2 at 90 %
identity** and kept **the longest member of each cluster** to build their
phylogenetic representative set of **481**. Each of the 481 is therefore a
`source_redundancy_representative` from that upstream reduction.

That representative is **not** our operational split unit. Our indivisible unit
for the Phase B partitions is the **`split_group`** — a connected component of an
MMseqs2 90 % identity / 80 % coverage graph computed fresh over the 481 sequences
(see §5). The v1 curation, which keyed the split unit on a hash of the accession
(one singleton per accession), is **superseded** and was never used for
predictive evaluation.

`apaz_partitions.csv` carries both `split_group_id` and, for backward
compatibility only, `homology_cluster_id`; the two columns are required to be
equal.

## 4. Hard negatives — HisG (PF01634) and EIIB (PF00367)

During its APAZ searches Ryazansky et al. recovered hits to the **HisG** family
(ATP phosphoribosyltransferase, C-terminal domain, Pfam **PF01634**) and to the
**EIIB** family (PTS system EIIB component, Pfam **PF00367**). After **reverse
DELTA-BLAST** they concluded these were **search artefacts, not true APAZ**. We
use their Pfam 38.2 seed alignments as *hard negatives*: homologous-looking but
established non-APAZ families.

- PF01634.24 seed: 509 sequences, sha256 pinned in `seeds_lock.json`.
- PF00367.26 seed (`#=GF ID PTS_EIIB`): 517 sequences, sha256 pinned.

Hard negatives never enter BUILD. They are partitioned **by family**, by whole
`split_group`, into CALIBRATION and FINAL_HOLDOUT only.

**Epistemological limitation.** Splitting HisG between CALIBRATION and
FINAL_HOLDOUT measures generalisation to **new sequences of the same negative
family**, not generalisation to **entirely new negative families**. The same
applies to EIIB. A future calibration should not read a low FINAL_HOLDOUT false
positive rate on HisG/EIIB as evidence that the model rejects unseen non-APAZ
families.

## 5. Split-group method (MMseqs2 90 / 80)

- Deterministic ungapped FASTA rebuilt from the frozen source (S3 for APAZ, Pfam
  seed for HisG/EIIB).
- MMseqs2 **18.8cc5c** (bioconda, build `hd6d6fdc_0`), all-vs-all, **nofilter**:
  `mmseqs search DB DB alnDB tmp --prefilter-mode 2 --max-seqs 100000 -e 1e100
  --min-seq-id 0 -c 0 -a --alignment-mode 3`. Environment: WSL 2.7.12.0 /
  Ubuntu 24.04.4 / conda 26.3.2 / env `pago-linux`.
- MMseqs2 `pident` / `alnlen` are **not** used. For every candidate pair, identity
  and coverage are recomputed from the aligned strings `qaln` / `taln`:
  - `comparable` = columns where both sequences carry a residue (no gap);
  - `identity`   = identical residue columns / comparable columns;
  - `coverage`   = comparable columns / min(query ungapped length, target ungapped length).
- An unordered pair is an **edge** iff `identity >= 0.90 AND coverage >= 0.80`
  over at least one reported orientation.
- **Connected components** by single-linkage union-find (transitive: A~B, B~C
  puts A, B, C in one group even if A≁C). The exact-duplicate rule (`sequence_sha256`
  union) is applied **before** any alignment edge.
- `split_group_id = <PREFIX>_<sha256 of the newline-joined sorted member ids>[:16]`
  — derived from the sorted member list, never from input order.

A 90 / 80 split group is a **high-similarity / redundancy control unit. It is not
a claim that homology is absent below 90 % identity** — HisG and EIIB are each a
single homologous family whose members mostly sit well below 90 % pairwise
identity, yet every member is genuinely homologous.

Result: APAZ 460 split groups (440 singletons, 19 pairs, 1 triple);
HisG 509 singletons; EIIB 517 singletons.

## 6. Partition policy (`apaz_partition_v2_mmseqs90_80`)

`PARTITION_SALT = pago-phase-b-apaz-partitions-v2-mmseqs90_80`.

Per stratum (APAZ subgroup Ia/Ib/IIa/IIb/III for positives; Pfam family for
negatives), fill **FINAL_HOLDOUT first, then CALIBRATION, then BUILD gets the
remainder**. Whole `split_group`s only. Groups are ordered for each partition by
`sha256(PARTITION_SALT | stratum | split_group_id | partition)`; a deterministic
subset-sum search (DFS, include-before-exclude, prune on overshoot and on
insufficient suffix sum) selects the first combination summing exactly to the
target. If no exact combination exists the run **aborts** — targets are never
adjusted silently. Two sequences with the same `sequence_sha256` are always in the
same `split_group`, hence the same partition.

Target counts (70 / 15 / 15 for positives):

| stratum | BUILD | CALIBRATION | FINAL_HOLDOUT |
|---|---|---|---|
| Ia | 59 | 12 | 12 |
| Ib | 77 | 16 | 16 |
| IIa | 68 | 15 | 15 |
| IIb | 125 | 27 | 27 |
| III | 8 | 2 | 2 |
| **positives** | **337** | **72** | **72** |
| HisG | 0 | 255 | 254 |
| EIIB | 0 | 259 | 258 |

## 7. Subgroup III limitation

Group III has only 12 representatives → BUILD 8, CALIBRATION 2, FINAL_HOLDOUT 2.
Two calibration positives cannot estimate III-specific performance with useful
precision. III is kept in the candidate construction (its BUILD HMM is built from
the 8 BUILD sequences), but:

- the GLOBAL vs SUBGROUP_MAX strategy choice is made in Phase B calibration by the
  documented protocol, not tuned for III on two examples;
- no III-specific threshold is set from those two examples;
- when the FINAL_HOLDOUT is evaluated (Phase B11), III-specific results are
  reported as statistically very uncertain and are not used to change the method.

## 8. Fragment `WP_077684787.1`

`WP_077684787.1` (group Ib) is a **38-residue** APAZ fragment. MMseqs2 found a
local alignment of **31 identical residues** to `WP_081684175.1` (137 aa) —
identity 1.000, coverage 0.816 of the shorter sequence — so the two are one
`split_group` and are inseparable across the partitions. It is **not excluded, not
moved by hand, and no length-based partition rule is applied**; length never
determines the partition. This pair is one of the eight cross-partition leaks in
v1 that v2 resolves (v1 had `WP_077684787.1` in BUILD and `WP_081684175.1` in
CALIBRATION). The frozen-column preview over the S3 alignment missed this pair
because the global S3 MSA does not place the fragment's residues on the homologous
columns; the de-novo local alignment does.

## 9. Reproducibility

`seeds_lock.json` (`lock_format_version` 2.0) pins: the S3 source (md5 + sha256),
the PMC metadata, both Pfam negative seeds (versioned accession + sha256), the six
`.sto` seeds, `apaz_partitions.csv`, `apaz_validation_sequences.fasta`, the three
`split_groups/<dataset>/{split_groups.csv,pairs.filtered.tsv,manifest.json,mmseqs_search.log}`,
the canonical library `apaz_split_groups.py` sha256, the MMseqs2 version/build,
the environment lock hashes, the identity/coverage/connected-components
definitions, `PARTITION_SALT`, the target counts, and the subset-sum / tie-break
policy. `python scripts/verify_reference_data.py --scope apaz` proves the whole
chain offline.

The global BUILD MSA (`apaz_global.sto`) is independently reconstructed from the
BUILD ungapped sequences with **PyFAMSA 0.7.0** (FAMSA 2, UPGMA guide tree). The
five subgroup seeds keep the original Ryazansky S3 alignment columns.
