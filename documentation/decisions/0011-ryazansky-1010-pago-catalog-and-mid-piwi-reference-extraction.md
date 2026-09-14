# ADR-0011 — Reconcile the Ryazansky 1010-pAgo catalog and extract MID-PIWI references

## Date

`UNKNOWN — reconstructed retrospectively on 2026-09-09` (historical commits
`b7490e4`, `481d0c2`, 2026-08-31)

## Context

The clade-reference layer needs a reproducible catalog of known pAgos and a
well-defined way to extract the MID-PIWI region used in later similarity and
phylogenetic work.

Ryazansky et al. 2018 Table S1 lists 1,010 pAgos before the paper's reported
90%-identity redundancy reduction, but the reusable 721-member tree-input
mapping, MID-PIWI tree alignment, and Newick tree are not available as
published supplementary artifacts.

## Decision

Build a deterministic 1,010-row catalog from Table S1.

For each row:

- recover the protein sequence by `accession.version`
- record NCBI record status rather than silently dropping suppressed or replaced
  entries
- keep source Ago type, source phylogenetic clade, architecture status, curated
  clade, and curation status as separate concepts
- prove the MID/PIWI coordinate convention from the table before slicing
- extract MID-PIWI as `sequence[MID_start-1 : PIWI_end]`
- exclude nothing solely because an accession is suppressed/replaced
- quarantine references whose architecture and source clade genuinely conflict
  instead of forcing a convenient label

## Evidence

Historical catalog notes and commits record:

- 1,010/1,010 sequences recovered
- 731 LIVE, 272 SUPPRESSED, and 7 DEAD_REPLACED records, with the seven
  replacements byte-identical to their historical sequences
- approximately 1,004 distinct full sequences after obsolete-accession
  duplication
- MID-PIWI extracted for 1,002/1,010 proteins
- 984 unique MID-PIWI regions
- 972 canonical and 38 truncated architectures
- AfAgo and SiAgo retained as unresolved because architecture and source clade
  conflict
- NgAgo retained as Ryazansky `longA`; the conflicting older recall-panel label
  is treated as a probable curation error rather than rewriting the source

## Rationale

Table S1 is the pre-reduction catalog, so entries should not be removed based on
an assumed redundancy filter. Coordinate semantics must be demonstrated before
sequence slicing. Architecture, phylogenetic clade, family, and catalytic
activity are distinct biological properties and should not overwrite one
another.

## Consequences

- 1,002 MID-PIWI regions become the input to the high-similarity grouping in
  ADR-0013
- the project's own reference-tree construction must be defined independently
  because the original reusable 721-member tree inputs are unavailable
- the exact redundant representatives used for some published tree positions
  cannot be reconstructed from the available supplements

## Supersedes

The earlier claim that TtAgo and RsAgo were absent from Table S1 because they
had been removed by the 90% redundancy reduction.

## Limitations

The original 721-member mapping, MID-PIWI tree alignment, and Newick tree are not
available from the published supplementary material. The canonical TtAgo
accession used by the separate recall panel still requires curation.

## Validation

Historical catalog verification checks record/coordinate invariants and guards
the MID-PIWI slice against off-by-one errors. The catalog can be rebuilt from
the frozen historical sources.

## Related records

- historical commits `b7490e4`, `481d0c2`
- historical Ryazansky Table S1 catalog and curation notes
- Ryazansky, Kulbachinskiy & Aravin 2018, *mBio* 9:e01935-18
- historical roadmap reference: B4.2

## Historical reconstruction note

- Directly demonstrated: historical catalog notes, counts, commits, and primary
  paper.
- Not recoverable: the exact 721-member tree-input mapping and the redundant
  homolog choices used in the published tree.
