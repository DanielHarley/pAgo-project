"""Pure parsing/coordinate logic for the Ryazansky 2018 pAgo catalog (B4.2).

No network, no filesystem side effects beyond reading the Table S1 workbook.
``scripts/build_clade_reference_catalog.py`` is the orchestration wrapper that
adds source fetching, NCBI sequence retrieval, QC and CSV output.

Nothing here clusters, aligns, partitions, or builds an HMM.
"""
from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TypeAlias

PathLike: TypeAlias = str | Path

TABLE_S1_HEADER = [
    "Accession", "Description", "Species", "NCBI taxon_id", "length",
    "Ago_type", "PAZ", "MID", "PIWI", "PAZ_type", "MID_type-5P",
    "MID_type-5OH", "PIWI_catalytic tetrad",
]

# Ryazansky's Ago_type -> pipeline phylogenetic-clade vocabulary. ``_trun`` is a
# truncated variant of the same clade; ``unkn`` is the two proteins the paper
# itself could not classify.
AGO_TYPE_TO_CLADE: Mapping[str, str] = {
    "longA": "LONG_A", "longA_trun": "LONG_A",
    "longB": "LONG_B", "longB_trun": "LONG_B",
    "short": "SHORT", "short_trun": "SHORT",
    "unkn": "UNRESOLVED",
}

MIDPIWI_COORDINATE_CONVENTION = "one_based_inclusive_contiguous_mid_then_piwi"

_RANGE = re.compile(r"^(\d+)-(\d+)$")


def parse_coordinate(value: object) -> tuple[int | None, int | None]:
    """Parse a Table S1 ``start-end`` cell; ``NULL``/blank -> ``(None, None)``."""
    match = _RANGE.match(str(value).strip())
    if match is None:
        return (None, None)
    return (int(match.group(1)), int(match.group(2)))


def parse_table_s1_rows(rows: Sequence[Sequence[object]]) -> list[dict]:
    """Parse Table S1 body rows (header already stripped) into record dicts."""
    records: list[dict] = []
    for index, row in enumerate(rows, start=2):
        accession = str(row[0]).strip()
        ago_type = str(row[5]).strip()
        if ago_type not in AGO_TYPE_TO_CLADE:
            raise RuntimeError(f"Unknown Ago_type {ago_type!r} at Table S1 row {index}.")
        paz_start, paz_end = parse_coordinate(row[6])
        mid_start, mid_end = parse_coordinate(row[7])
        piwi_start, piwi_end = parse_coordinate(row[8])
        length_value = row[4]
        records.append({
            "source_accession": accession,
            "species": str(row[2]).strip(),
            "taxon_id": (
                int(row[3]) if isinstance(row[3], float) else str(row[3])
            ),
            "source_length": (
                int(length_value) if isinstance(length_value, float) else None
            ),
            "description": str(row[1]).strip(),
            "source_ago_type": ago_type,
            "ago_family": "PAGO",
            "source_phylogenetic_clade": AGO_TYPE_TO_CLADE[ago_type],
            "truncated_flag": ago_type.endswith("_trun"),
            "source_paz_start": paz_start, "source_paz_end": paz_end,
            "source_mid_start": mid_start, "source_mid_end": mid_end,
            "source_piwi_start": piwi_start, "source_piwi_end": piwi_end,
            "source_paz_type": str(row[9]).strip(),
            "source_mid_5p_motif": str(row[10]).strip(),
            "source_piwi_tetrad": str(row[12]).strip(),
        })
    return records


def read_table_s1(path: PathLike) -> list[dict]:
    """Read the Table S1 ``.xls`` workbook into record dicts."""
    import xlrd

    book = xlrd.open_workbook(str(path))
    sheet = book.sheet_by_index(0)
    header = [str(sheet.cell_value(0, c)).strip() for c in range(sheet.ncols)]
    if header != TABLE_S1_HEADER:
        raise RuntimeError(f"Unexpected Table S1 header: {header}.")
    body = [
        [sheet.cell_value(r, c) for c in range(sheet.ncols)]
        for r in range(1, sheet.nrows)
    ]
    records = parse_table_s1_rows(body)
    if len(records) != 1010:
        raise RuntimeError(f"Expected 1010 Table S1 rows, parsed {len(records)}.")
    if len({r["source_accession"] for r in records}) != len(records):
        raise RuntimeError("Table S1 contains a duplicate accession string.")
    return records


def coordinate_convention_evidence(records: Sequence[Mapping]) -> dict:
    """Tally the signatures that fix the PAZ/MID/PIWI coordinate convention."""
    piwi_over_length = piwi_equals_length = 0
    mid_end_plus_one_eq_piwi_start = mid_end_eq_piwi_start = 0
    mid_piwi_gap_or_overlap = 0
    paz_before_mid = paz_pairs = 0
    for record in records:
        length = record["source_length"]
        mid_end = record["source_mid_end"]
        piwi_start = record["source_piwi_start"]
        piwi_end = record["source_piwi_end"]
        if piwi_end is not None and length is not None:
            piwi_over_length += piwi_end > length
            piwi_equals_length += piwi_end == length
        if mid_end is not None and piwi_start is not None:
            if mid_end + 1 == piwi_start:
                mid_end_plus_one_eq_piwi_start += 1
            elif mid_end == piwi_start:
                mid_end_eq_piwi_start += 1
            else:
                mid_piwi_gap_or_overlap += 1
        if record["source_paz_end"] is not None and record["source_mid_start"] is not None:
            paz_pairs += 1
            paz_before_mid += record["source_paz_end"] < record["source_mid_start"]
    return {
        "piwi_end_gt_length": piwi_over_length,
        "piwi_end_eq_length": piwi_equals_length,
        "mid_end_plus_one_eq_piwi_start": mid_end_plus_one_eq_piwi_start,
        "mid_end_eq_piwi_start": mid_end_eq_piwi_start,
        "mid_piwi_gap_or_overlap": mid_piwi_gap_or_overlap,
        "paz_end_lt_mid_start": paz_before_mid,
        "paz_pairs": paz_pairs,
    }


def prove_coordinate_convention(records: Sequence[Mapping]) -> dict:
    """Prove the coordinates are 1-based, inclusive, contiguous (MID then PIWI).

    Raises ``RuntimeError`` if any signature disagrees. Returns the evidence
    dict with a ``convention`` key on success.
    """
    evidence = coordinate_convention_evidence(records)
    ok = (
        evidence["piwi_end_gt_length"] == 0
        and evidence["mid_end_plus_one_eq_piwi_start"] > 0
        and evidence["mid_end_eq_piwi_start"] == 0
        and evidence["mid_piwi_gap_or_overlap"] == 0
        and evidence["paz_end_lt_mid_start"] == evidence["paz_pairs"]
    )
    if not ok:
        raise RuntimeError(
            "Ryazansky coordinate convention is not "
            f"1-based/inclusive/contiguous: {evidence}."
        )
    return {**evidence, "convention": MIDPIWI_COORDINATE_CONVENTION}


def slice_midpiwi(sequence: str, *, mid_start: int, piwi_end: int) -> str:
    """Return the MID-PIWI region under the 1-based inclusive convention.

    ``mid_start`` and ``piwi_end`` are 1-based inclusive residue positions, so
    the region is ``sequence[mid_start - 1 : piwi_end]``.
    """
    if mid_start < 1:
        raise ValueError(f"mid_start must be >= 1, got {mid_start}.")
    if piwi_end > len(sequence):
        raise ValueError(
            f"piwi_end {piwi_end} exceeds sequence length {len(sequence)}."
        )
    if mid_start > piwi_end:
        raise ValueError(f"mid_start {mid_start} > piwi_end {piwi_end}.")
    return sequence[mid_start - 1:piwi_end]
