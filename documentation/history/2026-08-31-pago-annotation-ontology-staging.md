<!--
================================================================================
POINT-IN-TIME HISTORICAL RECORD
Date: 2026-08-31 (approximate — drafted during the Phase B staging work)
Historical branch: (feat)-siepe-ready-project
Status when written: STAGING DOCUMENT, never integrated into master
Authoritative current state: NO
Preserved for provenance and as a source for the future authoritative ontology.
Current decisions must be read from documentation/SCIENTIFIC_STATE.md and the
current documentation/decisions/ ADRs.
The body below is the original staging document, preserved without semantic
rewriting.

Historical terminology note: This staging document uses "homology cluster" in
several places for the operational 90% identity / 80% coverage grouping. Current
project terminology calls these "high-similarity split groups". The 90/80 rule
is not a definition of homology. This historical copy is preserved without
semantic rewriting.
================================================================================
-->

# pAgo annotation ontology and evidence boundaries

This document defines the semantic contract for protein and locus annotations in
the pAgo pipeline. It prevents a computationally plausible label from being
silently promoted to biological identity, phylogenetic clade, catalytic
activity, or complete defense-system membership.

The Phase B reference workflow builds and validates evidence instruments. It
does not scan the annotation-enriched proteome and does not assign biological
labels to its proteins. Those operations begin in Phase C.

## Orthogonal annotation dimensions

### Detection state

`ago_detection_status` records the current epistemic conclusion about whether a
protein is an Argonaute-family protein.

| Value | Meaning |
| --- | --- |
| `DETECTED` | Validated core Argonaute evidence is present. |
| `NOT_DETECTED` | Explicit, sufficiently strong negative evidence is present. |
| `AMBIGUOUS` | Material evidence sources are incompatible. |
| `UNRESOLVED` | Available evidence cannot support a positive or negative decision. |

`detection_evidence_status` records the evidence class used for that state.

| Value | Meaning |
| --- | --- |
| `CORE_HMM` | A validated core PIWI or MID Argonaute profile supports detection. |
| `TRIAGE_ONLY` | A MID–PIWI clade profile produced a candidate signal without validated core detection evidence. |
| `NO_CORE_HMM_EVIDENCE` | No validated core signal is available. Absence alone is not a biological negative. |
| `NEGATIVE_EVIDENCE` | An explicit negative-evidence rule supports `NOT_DETECTED`. |

The containment rules are:

1. A competitive clade HMM result cannot set `ago_detection_status=DETECTED`.
2. Absence of a Pfam PIWI hit cannot set `ago_detection_status=NOT_DETECTED`.
3. A triage-only signal remains `UNRESOLVED` until independent evidence supports
   a family assignment.
4. Conflicting evidence is retained as `AMBIGUOUS` rather than collapsed to a
   convenient positive or negative label.

### Family identity

`ago_family` is a biological identity field with the vocabulary `PAGO`,
`PIWI_RE`, and `UNRESOLVED`. PIWI-RE is a family-level assignment. It is not a
pAgo clade.

`ago_family_evidence` identifies the decision surface, such as `PLACEMENT`,
`PIWI_RE_RESIDUES`, `HMM`, or `NONE`. A family assignment must retain its
confidence and rationale rather than relying on a bare categorical value.

### pAgo clade

`pago_clade` has the vocabulary `LONG_A`, `LONG_B`, `SHORT`, and `UNRESOLVED`.
It is a phylogenetic property defined only within the PAGO family.

The competitive MID–PIWI HMM result is stored separately as
`pago_clade_hmm_candidate`. It is a triage feature and never becomes
`pago_clade` by itself. A resolved clade requires phylogenetic placement against
the frozen MID–PIWI reference tree and all three calibrated decision criteria:

1. `best_clade_lwr_sum >= t1`
2. `best_clade_lwr_sum - second_clade_lwr_sum >= t2`
3. `unresolved_lwr <= t3`

Failure of any criterion yields `pago_clade=UNRESOLVED`. Placement mass assigned
to PIWI-RE can support `ago_family=PIWI_RE`, while `pago_clade` remains
`UNRESOLVED`.

### Architecture

`architecture_pfam` is derived only from resolved Pfam-domain hits.
`architecture_extended` can additionally include an APAZ event that passed the
frozen APAZ validation policy. Raw Pfam, APAZ, and competitive clade-HMM results
remain in disjoint namespaces:

* `domain_hits_raw` and `domain_hits_resolved`
* `apaz_hits_raw` and `apaz_hits`
* `clade_hmm_scores`

An overlap resolver may transform Pfam hits into `domain_hits_resolved`; it must
not destroy `domain_hits_raw` or combine unrelated evidence namespaces.

An architecture such as `MID_PIWI` does not imply `pago_clade=SHORT`.

### Catalytic-site status

`catalytic_site_status` is fixed to `UNKNOWN` in pipeline version 1. Catalytic
motifs or residue annotations do not decide clade, family, detection status, or
system class. Literature-derived catalytic annotations can be retained in
separate audit fields without being promoted to a computed activity claim.

### Redundancy and supervised-learning splits

`redundancy_cluster_id_90` describes sequence redundancy under a declared
identity and coverage policy. It is not an ML split assignment.

`ml_split_group` remains null in version 1 until a supervised generalization
design is specified. Any future supervised evaluation must keep complete
homology clusters in a single split. Random row-level splitting is outside this
contract because it can leak near-identical proteins across training and test
sets.

## Protein and locus units

Protein-sequence properties belong in `pago_protein_annotation.csv`, with one
row per `protein_uid`. Genomic-context properties belong in
`pago_locus_annotation.csv`, with zero or more rows per protein because a protein
may map to multiple genomic occurrences through IPG. Neighbor observations
belong in `pago_locus_neighbors.csv`, with one row per locus and neighbor.

`pago_protein_system_summary.csv` aggregates resolved locus evidence while
preserving two distinct states:

* absence of a resolved locus conclusion
* disagreement among resolved locus conclusions

The aggregation vocabulary therefore includes `MIXED`, `NO_CONSENSUS`, and
`UNRESOLVED` where appropriate.

`SPARTA`, `SPARSA`, and `SPARDA` require a complete fused architecture or
supported genomic context. Partial evidence produces an associated-candidate or
unresolved state. `OTHER` requires a complete incompatible system; missing
evidence cannot produce `OTHER`.

## Reference-label evidence

Every curated reference label records both its source and evidence tier.

| Evidence tier | Intended interpretation |
| --- | --- |
| `EXPERIMENTAL` | Direct experimental evidence relevant to the recorded label. |
| `LITERATURE_PHYLOGENETIC` | A phylogenetic label reported in the literature. |
| `CURATED_COMPUTATIONAL` | A reproducible computational curation decision. |
| `DATABASE_ANNOTATION` | A database annotation without stronger evidence in this dataset. |

These tiers are provenance metadata. They must not be treated as mutually
interchangeable ground truth.

## HMM classifier validation protocol

Every custom HMM classifier uses immutable, cluster-separated partitions.

1. `BUILD` contains only sequences used to construct the HMM.
2. `CALIBRATION` selects the strategy, thresholds, coverage constraints, and
   confidence rules.
3. `FINAL_HOLDOUT` is consulted exactly once after the strategy and thresholds
   are frozen.

Positive and hard-negative homology clusters are assigned wholly to one
partition. A sequence or homology cluster crossing partition boundaries is a
validation error. The APAZ hard negatives include HisG (`PF01634`) and EIIB
(`PF00367`) reference families because related analyses identified them as
important false-homology controls.

Calibration writes `apaz_strategy.json` and `apaz_thresholds.json`, each with a
content identity. Final evaluation writes an immutable report keyed by the
frozen policy and holdout identities. Reusing the same scientific evaluation
identity is rejected. A method changed after final evaluation requires a new,
previously unconsulted holdout and a new evaluation identity.

Final scientific holdout reports are versioned evidence artifacts. They are not
routine CI assertions. CI validates code behavior and committed reference
integrity without repeatedly optimizing against or reinterpreting final
holdouts.

## Phylogenetic placement validation protocol

Clade evaluation uses cluster-separated `TREE_BUILD`,
`PLACEMENT_CALIBRATION`, and `PLACEMENT_HOLDOUT` partitions. Calibration chooses
the three likelihood-weight-ratio thresholds. The final report preserves:

1. family-level PAGO versus PIWI-RE performance
2. clade-level LONG_A versus LONG_B versus SHORT performance within PAGO
3. the raw four-group confusion matrix

This hierarchy distinguishes a family error such as SHORT to PIWI-RE from an
intra-pAgo error such as SHORT to LONG_B.

Reference sequences and query sequences must be aligned with the same canonical
`reference.hmm`. The pre-placement invariant is exact equality between
`reference_profile_sha256` and `query_alignment_profile_sha256`, together with
an identical alignment-column contract. Shared ancestral or backbone edges are
mapped to `UNRESOLVED` in `edge_clade_map.csv`.

## Reference integrity

Committed Phase B resources are under `src/pago_pipeline/resources/`. Their
locks pin source releases, file inventories, SHA-256 values, partition tables,
alignment dimensions, sequence identities, and tree metadata where applicable.

Run the read-only integrity verifier from the repository root:

```powershell
python scripts/verify_reference_data.py
```

The verifier validates the Pfam bundle inputs, APAZ seeds and partitions, clade
seeds and partitions, and the canonical MID–PIWI profile, alignment, tree,
model, metadata, and edge map. It performs no network request and writes no
artifact.

## Phase boundary

Phase B ends after the reference resources, HMM builds, frozen calibration
policies, and one-shot validation reports have been validated and recorded.
Phase B does not perform the following operations:

* scan the 52,473-protein annotation-enriched dataset
* assign detection, family, architecture, or clade labels to candidate proteins
* reconstruct genomic neighborhoods
* assign SPARTA, SPARSA, or SPARDA system classes
* generate SWeeP embeddings, PCA, clustering, or three-dimensional plots

Candidate scanning and protein annotation begin in Phase C. Locus inference and
integrated exploratory analyses belong to later reviewed phases.

## Literature and database anchors

The principal clade and APAZ reference source is Ryazansky et al. (2018),
*Prokaryotic Argonaute Proteins: Phylogenetic Distribution, Pathways of
Evolution, and the Role of Catalytic Activities in Biological Functions*,
[mBio 9:e01935-18](https://pmc.ncbi.nlm.nih.gov/articles/PMC6299218/).
PIWI-RE family context is anchored in Burroughs et al. (2013),
[Biology Direct 8:13](https://pmc.ncbi.nlm.nih.gov/articles/PMC3702460/).
Pfam profile data follow the release pinned in
`pfam_hmm_bundle_lock.json` and the
[InterPro data license](https://www.ebi.ac.uk/interpro/about/license/).
