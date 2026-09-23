"""Fetch and validate the external APAZ reference sources (maintenance, off-CI).

Responsibility boundary:

  * THIS script (network) fetches the Ryazansky Data Set S3 and the Pfam HisG /
    EIIB seed alignments, validates their identity (md5 / versioned Pfam
    accession / InterPro release / PMC licence / Stockholm completeness / the
    exact APAZ group counts), and writes ONLY the frozen inputs:

        resources/apaz_seed/source_data/
            mbo006184236sd3.txt
            PMC6299218.1.metadata.json
            PF01634.seed.sto
            PF00367.seed.sto

  * ``scripts/derive_apaz_split_groups.py`` (needs MMseqs2) turns those sources
    into the frozen, validated ``resources/apaz_seed/split_groups/`` resources.

  * ``scripts/regenerate_apaz_reference.py`` (offline, needs pyfamsa) turns the
    sources + split groups into ``apaz_partitions.csv``, the five subgroup seeds,
    ``apaz_global.sto``, ``apaz_validation_sequences.fasta`` and
    ``seeds_lock.json``.

All three share the single canonical library
``src.pago_pipeline.apaz_split_groups``. This script performs no clustering and
no partitioning.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import requests

SOURCE_ARTICLE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6299218/"
SOURCE_METADATA_URL = "https://pmc-oa-opendata.s3.amazonaws.com/metadata/PMC6299218.1.json"
SOURCE_ALIGNMENT_URL = (
    "https://pmc-oa-opendata.s3.amazonaws.com/PMC6299218.1/mbo006184236sd3.txt"
)
SOURCE_ALIGNMENT_MD5 = "ac30768ca799feb85165cea9d11aa3e2"
INTERPRO_ENTRY_API = "https://www.ebi.ac.uk/interpro/api/entry/pfam"
INTERPRO_RELEASE = "109.0"
PFAM_RELEASE = "38.2"

APAZ_SOURCE_COUNTS = {"Ia": 83, "Ib": 109, "IIa": 98, "IIb": 179, "III": 12}
APAZ_SUBGROUP_ORDER = ("Ia", "Ib", "IIa", "IIb", "III")
NEGATIVE_FAMILIES = {"PF01634": "HisG", "PF00367": "EIIB"}

DEFAULT_SOURCE_DATA_DIRECTORY = (
    Path(__file__).resolve().parents[1]
    / "src" / "pago_pipeline" / "resources" / "apaz_seed" / "source_data"
)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _fetch(session: requests.Session, url: str, timeout_seconds: float) -> bytes:
    response = session.get(url, timeout=timeout_seconds)
    response.raise_for_status()
    return response.content


def _parse_apaz_alignment(payload: bytes) -> None:
    """Validate the Ryazansky S3 alignment: group counts, unique accessions, width."""
    text = payload.decode("utf-8-sig")
    current_subgroup: str | None = None
    counts: dict[str, int] = {subgroup: 0 for subgroup in APAZ_SUBGROUP_ORDER}
    accessions: list[str] = []
    widths: set[int] = set()
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        heading = [g for g in APAZ_SUBGROUP_ORDER if ("Group %s:" % g) in line]
        if heading:
            current_subgroup = heading[0]
            continue
        if line.startswith("#"):
            continue
        if current_subgroup is None:
            raise RuntimeError(f"APAZ sequence row precedes a subgroup heading at line {line_number}.")
        fields = line.split()
        if len(fields) != 2:
            raise RuntimeError(f"Malformed APAZ source row at line {line_number}.")
        identifier, aligned = fields
        identifier_fields = identifier.split("|")
        if len(identifier_fields) < 3 or not aligned:
            raise RuntimeError(f"Malformed APAZ source identifier at line {line_number}.")
        counts[current_subgroup] += 1
        accessions.append(identifier_fields[0])
        widths.add(len(aligned))
        if not set(aligned) <= set("ABCDEFGHIJKLMNOPQRSTUVWXYZ-."):
            raise RuntimeError(f"Invalid APAZ aligned sequence at line {line_number}.")
    if counts != APAZ_SOURCE_COUNTS:
        raise RuntimeError(f"APAZ source counts drifted: expected {APAZ_SOURCE_COUNTS}, observed {counts}.")
    if len(accessions) != len(set(accessions)):
        raise RuntimeError("APAZ source contains duplicate accessions.")
    if widths != {173}:
        raise RuntimeError(f"Unexpected APAZ source alignment widths: {sorted(widths)}.")


def _stockholm_family_identity(payload: bytes) -> tuple[str, str]:
    text = payload.decode("utf-8-sig")
    accession_match = re.search(r"^#=GF AC\s+(\S+)", text, re.MULTILINE)
    name_match = re.search(r"^#=GF ID\s+(\S+)", text, re.MULTILINE)
    if accession_match is None or name_match is None:
        raise RuntimeError("Pfam seed alignment lacks #=GF AC or #=GF ID metadata.")
    return accession_match.group(1), name_match.group(1)


def _validate_stockholm(payload: bytes) -> int:
    text = payload.decode("utf-8-sig")
    fragments: dict[str, list[str]] = defaultdict(list)
    first_content: str | None = None
    terminated = False
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if first_content is None:
            first_content = line
        if line == "//":
            terminated = True
            continue
        if line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 2:
            raise RuntimeError(f"Malformed Stockholm row at line {line_number}.")
        fragments[fields[0]].append(fields[1])
    if first_content != "# STOCKHOLM 1.0" or not terminated or not fragments:
        raise RuntimeError("Pfam seed response is not a complete Stockholm alignment.")
    aligned = {identifier: "".join(parts) for identifier, parts in fragments.items()}
    if len({len(sequence) for sequence in aligned.values()}) != 1:
        raise RuntimeError("Pfam seed alignment contains unequal sequence widths.")
    return len(aligned)


def fetch_sources(*, source_data_directory: Path, timeout_seconds: float) -> dict[str, str]:
    source_data_directory.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update(
        {"User-Agent": "pAgo-project APAZ reference source fetcher (offline reference preparation)"}
    )

    interpro_response = session.get(
        f"{INTERPRO_ENTRY_API}/{next(iter(NEGATIVE_FAMILIES))}/", timeout=timeout_seconds
    )
    interpro_response.raise_for_status()
    observed_release = interpro_response.headers.get("InterPro-Version")
    if observed_release != INTERPRO_RELEASE:
        raise RuntimeError(
            f"InterPro release drift. Expected {INTERPRO_RELEASE}, got {observed_release!r}."
        )

    metadata_payload = _fetch(session, SOURCE_METADATA_URL, timeout_seconds)
    metadata = json.loads(metadata_payload)
    if metadata.get("license_code") != "CC BY":
        raise RuntimeError("Unexpected PMC licence for the APAZ source article.")

    alignment_payload = _fetch(session, SOURCE_ALIGNMENT_URL, timeout_seconds)
    if hashlib.md5(alignment_payload).hexdigest() != SOURCE_ALIGNMENT_MD5:
        raise RuntimeError("Ryazansky Data Set S3 md5 mismatch.")
    _parse_apaz_alignment(alignment_payload)

    written: dict[str, str] = {}
    metadata_path = source_data_directory / "PMC6299218.1.metadata.json"
    alignment_path = source_data_directory / "mbo006184236sd3.txt"
    _write_bytes(metadata_path, metadata_payload)
    _write_bytes(alignment_path, alignment_payload)
    written[metadata_path.name] = _sha256(metadata_payload)
    written[alignment_path.name] = _sha256(alignment_payload)

    for family_accession, family_name in NEGATIVE_FAMILIES.items():
        seed_url = f"{INTERPRO_ENTRY_API}/{family_accession}/?annotation=alignment:seed"
        compressed = _fetch(session, seed_url, timeout_seconds)
        seed_payload = (
            gzip.decompress(compressed) if compressed.startswith(b"\x1f\x8b") else compressed
        )
        versioned_accession, seed_name = _stockholm_family_identity(seed_payload)
        if versioned_accession.split(".", 1)[0] != family_accession:
            raise RuntimeError(
                f"Pfam seed accession mismatch for {family_accession}: {versioned_accession}."
            )
        sequence_count = _validate_stockholm(seed_payload)
        seed_path = source_data_directory / f"{family_accession}.seed.sto"
        _write_bytes(seed_path, seed_payload)
        written[seed_path.name] = _sha256(seed_payload)
        print(
            f"  {family_accession} ({family_name}) -> {seed_path.name}  "
            f"{versioned_accession} / {seed_name}  {sequence_count} sequences"
        )

    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-data-directory", type=Path, default=DEFAULT_SOURCE_DATA_DIRECTORY)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    arguments = parser.parse_args(argv)
    written = fetch_sources(
        source_data_directory=arguments.source_data_directory,
        timeout_seconds=arguments.timeout_seconds,
    )
    print("\nwrote frozen sources:")
    for name, digest in sorted(written.items()):
        print("  %-30s %s" % (name, digest))
    print(
        "\nnext: scripts/derive_apaz_split_groups.py  (needs MMseqs2)\n"
        "then: scripts/regenerate_apaz_reference.py  (needs pyfamsa==0.7.0)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
