"""Build the deterministic Ryazansky 2018 pAgo catalog (Phase B, sub-step B4.2).

Responsibility boundary:

  * Fetches and freezes the three Ryazansky 2018 supplements needed for the
    MID-PIWI clade reference (Table S1, Data Set S1, Data Set S2) plus the
    official PMC metadata, verifying every md5 against that metadata, and writes
    them under ``resources/clade_seed/source_data/`` with a ``source_manifest``.

  * Fetches every protein sequence named in Table S1 by ``accession.version``
    (plus the current replacement for each dead accession), freezing them in
    ``ryazansky_s1_sequences.fasta`` so the catalog rebuilds offline.

  * Parses Table S1 into ``ryazansky_s1_catalog.csv`` (1010 rows, one per
    Table S1 entry). Preserves ``source_ago_type`` / ``source_phylogenetic_clade``
    and only sets ``curated_pago_clade`` = UNRESOLVED for the two proteins
    Ryazansky itself declared unclassifiable and for QUARANTINE cases.

This script performs NO clustering, NO split-group construction, NO
partitioning, NO HMM construction, NO alignment. It never writes
``clade_partitions.csv`` and never touches the APAZ (B2) or APAZ-HMM (B3)
resources.

The MID-PIWI region is extracted purely from the PAZ/MID/PIWI coordinate columns
of Table S1, whose convention (1-based, inclusive, contiguous MID|PIWI) is
proven from the data and asserted before any slicing.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import requests

PROJECT_ROOT = Path(__file__).resolve().parents[1]
import sys

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pago_pipeline.clade_reference_catalog import (  # noqa: E402
    AGO_TYPE_TO_CLADE,
    prove_coordinate_convention,
    read_table_s1,
    slice_midpiwi,
)
DEFAULT_CLADE_SEED_DIRECTORY = (
    PROJECT_ROOT / "src" / "pago_pipeline" / "resources" / "clade_seed"
)

PMC_ID = "PMC6299218.1"
SOURCE_ARTICLE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6299218/"
SOURCE_METADATA_URL = (
    f"https://pmc-oa-opendata.s3.amazonaws.com/metadata/{PMC_ID}.json"
)
PMC_OA_BASE_URL = f"https://pmc-oa-opendata.s3.amazonaws.com/{PMC_ID}/"
DOI = "10.1128/mBio.01935-18"
PMID = 30563906

# Files materialized as B4 sources, with the role each plays.
MATERIALIZED_SUPPLEMENTS = {
    "mbo006184236st1.xls": (
        "table_s1",
        "The list of 1010 identified pAgo proteins: accession, species, taxon, "
        "length, Ago_type clade label, and PAZ/MID/PIWI domain coordinates.",
    ),
    "mbo006184236sd1.txt": (
        "data_set_s1",
        "Alignment of the MID* 5'-end-binding-motif fragment for 75 long pAgos "
        "with a non-canonical MID motif (diagnostic; not the tree alignment).",
    ),
    "mbo006184236sd2.txt": (
        "data_set_s2",
        "Alignment of PAZ and PAZ* domains of long pAgos (diagnostic; not the "
        "tree alignment).",
    ),
}
# Recorded in the audit notes only, not materialized under clade_seed.
NON_MATERIALIZED_SUPPLEMENTS = (
    "mbo006184236st2.xlsx",  # operon-neighbour orthogroups (genomic context)
    "mbo006184236st3.xls",   # operon size / orthogroup frequency
)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = {
    "User-Agent": "pAgo-project clade reference catalog (offline reference preparation)"
}

# Recall-set anchor overrides that must not be silently "corrected".
QUARANTINE_ANCHORS = {"AfAgo", "SiAgo"}


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def _md5_bytes(payload: bytes) -> str:
    return hashlib.md5(payload).hexdigest()


def _normalize_sequence(sequence: str) -> str:
    return re.sub(r"\s", "", sequence).upper()


def _session() -> requests.Session:
    session = requests.Session()
    session.headers.update(USER_AGENT)
    return session


def _get(session: requests.Session, url: str, *, params=None, timeout=120) -> bytes:
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            response = session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.content
        except Exception as error:  # noqa: BLE001
            last_error = error
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET failed after retries: {url}: {last_error}")


def fetch_sources(*, source_data_directory: Path, session: requests.Session) -> dict:
    source_data_directory.mkdir(parents=True, exist_ok=True)

    metadata_payload = _get(session, SOURCE_METADATA_URL)
    metadata = json.loads(metadata_payload)
    if metadata.get("license_code") != "CC BY":
        raise RuntimeError("Unexpected PMC licence for the Ryazansky 2018 article.")
    frozen_b2_metadata = (
        PROJECT_ROOT / "src" / "pago_pipeline" / "resources" / "apaz_seed"
        / "source_data" / f"{PMC_ID}.metadata.json"
    )
    if frozen_b2_metadata.is_file():
        if _sha256_bytes(metadata_payload) != _sha256_bytes(
            frozen_b2_metadata.read_bytes()
        ):
            raise RuntimeError(
                "PMC metadata for PMC6299218.1 no longer matches the copy frozen "
                "in B2 (apaz_seed/source_data); provenance chain broken."
            )
    declared_md5 = {}
    for media_url in metadata["media_urls"]:
        name = media_url.split("/")[-1].split("?")[0]
        declared_md5[name] = media_url.split("md5=")[-1]

    manifest_files = []
    metadata_path = source_data_directory / f"{PMC_ID}.metadata.json"
    metadata_path.write_bytes(metadata_payload)
    manifest_files.append({
        "file_name": metadata_path.name,
        "role": "pmc_oa_metadata",
        "description": "Official PMC Open Access metadata record for the article.",
        "url": SOURCE_METADATA_URL,
        "bytes": len(metadata_payload),
        "md5": _md5_bytes(metadata_payload),
        "sha256": _sha256_bytes(metadata_payload),
    })

    for name, (role, description) in MATERIALIZED_SUPPLEMENTS.items():
        payload = _get(session, PMC_OA_BASE_URL + name)
        observed_md5 = _md5_bytes(payload)
        if observed_md5 != declared_md5.get(name):
            raise RuntimeError(
                f"md5 mismatch for {name}: metadata declares "
                f"{declared_md5.get(name)}, downloaded {observed_md5}."
            )
        (source_data_directory / name).write_bytes(payload)
        manifest_files.append({
            "file_name": name,
            "role": role,
            "description": description,
            "url": PMC_OA_BASE_URL + name,
            "bytes": len(payload),
            "md5": observed_md5,
            "sha256": _sha256_bytes(payload),
        })

    # Non-materialized supplements: record identity only.
    non_materialized = []
    for name in NON_MATERIALIZED_SUPPLEMENTS:
        payload = _get(session, PMC_OA_BASE_URL + name)
        observed_md5 = _md5_bytes(payload)
        if observed_md5 != declared_md5.get(name):
            raise RuntimeError(f"md5 mismatch for {name}.")
        non_materialized.append({
            "file_name": name,
            "url": PMC_OA_BASE_URL + name,
            "bytes": len(payload),
            "md5": observed_md5,
            "sha256": _sha256_bytes(payload),
            "note": "genomic-context data (operon neighbour orthogroups); not a "
                    "clade or sequence resource; recorded for provenance only.",
        })

    manifest = {
        "artifact_type": "clade_reference_source_manifest",
        "source_article_url": SOURCE_ARTICLE_URL,
        "pmc_id": PMC_ID,
        "doi": DOI,
        "pmid": PMID,
        "license_code": "CC BY",
        "pmc_oa_metadata_url": SOURCE_METADATA_URL,
        "pmc_oa_base_url": PMC_OA_BASE_URL,
        "materialized_files": manifest_files,
        "non_materialized_supplements": non_materialized,
    }
    (source_data_directory / "source_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def resolve_records(
    *, accessions: list[str], session: requests.Session
) -> dict[str, dict]:
    result: dict[str, dict] = {}
    for start in range(0, len(accessions), 50):
        chunk = accessions[start:start + 50]
        payload = _get(
            session, f"{EUTILS}/esummary.fcgi",
            params={"db": "protein", "id": ",".join(chunk), "retmode": "json"},
        )
        body = json.loads(payload).get("result", {})
        for uid in body.get("uids", []):
            record = body[uid]
            accession = record.get("accessionversion") or record.get("caption")
            result[accession] = {
                "status": record.get("status") or "live",
                "slen": record.get("slen"),
                "replacedby": record.get("replacedby") or "",
            }
        time.sleep(0.34)
    return result


def _read_frozen_fasta(fasta_path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = None
    for line in fasta_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            current = line[1:].strip()
            records[current] = ""
        elif current is not None:
            records[current] += line.strip()
    return records


def fetch_sequences(
    *,
    accessions: list[str],
    fasta_path: Path,
    session: requests.Session,
    reuse_frozen: bool = False,
) -> dict[str, str]:
    unique = sorted(set(accessions))
    if reuse_frozen and fasta_path.is_file():
        frozen = _read_frozen_fasta(fasta_path)
        if set(frozen) >= set(unique):
            return {a: _normalize_sequence(frozen[a]) for a in unique}
    records: dict[str, str] = {}
    for start in range(0, len(unique), 50):
        chunk = unique[start:start + 50]
        payload = _get(
            session, f"{EUTILS}/efetch.fcgi",
            params={"db": "protein", "id": ",".join(chunk),
                    "rettype": "fasta", "retmode": "text"},
        )
        current = None
        for line in payload.decode("utf-8").splitlines():
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = ""
            elif current is not None:
                records[current] += line.strip()
        time.sleep(0.34)

    def key_for(accession: str) -> str | None:
        if accession in records:
            return accession
        base = accession.split(".")[0]
        for candidate in records:
            tokens = candidate.replace("|", " ").split()
            if base in tokens or candidate.split(".")[0] == base:
                return candidate
        return None

    ordered = {}
    with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in unique:
            key = key_for(accession)
            if key is None:
                raise RuntimeError(f"NCBI returned no sequence for {accession}.")
            sequence = _normalize_sequence(records[key])
            ordered[accession] = sequence
            handle.write(f">{accession}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset:offset + 60] + "\n")
    return ordered


def build_catalog(
    *,
    records: list[dict],
    status_by_accession: dict[str, dict],
    sequence_by_accession: dict[str, str],
    recall_reference_path: Path,
) -> tuple[list[dict], dict]:
    convention = prove_coordinate_convention(records)

    recall_by_sha: dict[str, dict] = {}
    for recall_row in csv.DictReader(
        recall_reference_path.read_text(encoding="utf-8").splitlines()
    ):
        if recall_row["ago_family"] == "PAGO":
            recall_by_sha[recall_row["sequence_sha256"]] = recall_row

    catalog: list[dict] = []
    for record in records:
        accession = record["source_accession"]
        resolution = status_by_accession.get(accession, {"status": "unknown"})
        raw_status = resolution["status"]
        replacement = resolution.get("replacedby", "")
        ncbi_status = {
            "live": "LIVE", "suppressed": "SUPPRESSED", "dead": "DEAD_REPLACED",
        }.get(raw_status, raw_status.upper())

        sequence = sequence_by_accession.get(accession)
        if sequence is None:
            raise RuntimeError(f"Missing frozen sequence for {accession}.")
        full_sha = _sha256_text(sequence)
        length_current = len(sequence)
        length_matches_source = length_current == record["source_length"]

        replacement_sha = ""
        replacement_equivalent = ""
        if replacement:
            replacement_sequence = sequence_by_accession.get(replacement)
            if replacement_sequence is not None:
                replacement_sha = _sha256_text(replacement_sequence)
                replacement_equivalent = replacement_sha == full_sha

        mid_start = record["source_mid_start"]
        piwi_end = record["source_piwi_end"]
        midpiwi_sequence = ""
        midpiwi_length = ""
        midpiwi_region_sha256 = ""
        if mid_start is None:
            midpiwi_status = "MISSING_MID"
        elif mid_start < 1 or piwi_end > length_current or mid_start > piwi_end:
            midpiwi_status = "COORD_OUT_OF_RANGE"
        else:
            midpiwi_sequence = slice_midpiwi(
                sequence, mid_start=mid_start, piwi_end=piwi_end
            )
            midpiwi_length = len(midpiwi_sequence)
            midpiwi_region_sha256 = _sha256_text(midpiwi_sequence)
            midpiwi_status = "OK"

        anchor_row = recall_by_sha.get(full_sha)
        experimental_anchor = anchor_row is not None
        anchor_name = anchor_row["protein_short_name"] if anchor_row else ""

        curated_clade = record["source_phylogenetic_clade"]
        curation_status = "OK"
        quarantine_reason = ""
        notes: list[str] = []

        if record["source_ago_type"] == "unkn":
            curated_clade = "UNRESOLVED"
            notes.append(
                "Ryazansky declared this protein unclassifiable (high divergence); "
                "curated_pago_clade=UNRESOLVED; never a LONG_A/LONG_B/SHORT seed"
            )

        if anchor_row is not None:
            recall_clade = anchor_row["clade"]
            notes.append(
                f"experimental anchor {anchor_name} "
                f"(recall-set clade={recall_clade or 'UNRESOLVED'}, "
                f"evidence={anchor_row['reference_label_evidence']})"
            )
            if (
                recall_clade
                and recall_clade != "UNRESOLVED"
                and recall_clade != record["source_phylogenetic_clade"]
            ):
                if anchor_name in QUARANTINE_ANCHORS:
                    curation_status = "QUARANTINE"
                    curated_clade = "UNRESOLVED"
                    quarantine_reason = (
                        f"literature/architecture clade {recall_clade} vs "
                        f"Ryazansky phylogenetic {record['source_ago_type']}; "
                        "pending specific audit against the primary sources"
                    )
                elif anchor_name == "NgAgo":
                    notes.append(
                        "recall-set assigned LONG_B while citing Ryazansky, but "
                        "Table S1 records longA -> probable recall-set curation "
                        "error; Ryazansky classification retained; the Phase A "
                        "recall artifact is NOT modified here"
                    )
                else:
                    curation_status = "REVIEW"
                    notes.append(
                        f"recall-set clade {recall_clade} disagrees with "
                        f"Ryazansky {record['source_phylogenetic_clade']}; "
                        "flagged for review"
                    )

        catalog.append({
            "source_accession": accession,
            "resolved_accession": (
                replacement if ncbi_status == "DEAD_REPLACED" and replacement
                else accession
            ),
            "species": record["species"],
            "taxon_id": record["taxon_id"],
            "source_length": record["source_length"],
            "source_ago_type": record["source_ago_type"],
            "ago_family": "PAGO",
            "source_phylogenetic_clade": record["source_phylogenetic_clade"],
            "curated_pago_clade": curated_clade,
            "architecture_status": (
                "TRUNCATED" if record["truncated_flag"]
                else ("LONG_ARCHITECTURE" if record["source_paz_start"] is not None
                      else "SHORT_ARCHITECTURE")
            ),
            "truncated_flag": record["truncated_flag"],
            "source_paz_start": record["source_paz_start"],
            "source_paz_end": record["source_paz_end"],
            "source_mid_start": record["source_mid_start"],
            "source_mid_end": record["source_mid_end"],
            "source_piwi_start": record["source_piwi_start"],
            "source_piwi_end": record["source_piwi_end"],
            "source_paz_type": record["source_paz_type"],
            "source_mid_5p_motif": record["source_mid_5p_motif"],
            "source_piwi_tetrad": record["source_piwi_tetrad"],
            "ncbi_record_status": ncbi_status,
            "replacement_accession": replacement,
            "replacement_sequence_sha256": replacement_sha,
            "replacement_equivalent": replacement_equivalent,
            "full_sequence_sha256": full_sha,
            "sequence_length_current": length_current,
            "length_matches_source": length_matches_source,
            "midpiwi_sequence": midpiwi_sequence,
            "midpiwi_length": midpiwi_length,
            "midpiwi_region_sha256": midpiwi_region_sha256,
            "midpiwi_extraction_status": midpiwi_status,
            "experimental_anchor": experimental_anchor,
            "experimental_anchor_name": anchor_name,
            "curation_status": curation_status,
            "quarantine_reason": quarantine_reason,
            "curation_notes": " | ".join(notes),
            "proposed_partition": "UNASSIGNED",
        })

    quality = run_quality_checks(catalog)
    quality["coordinate_convention"] = convention
    return catalog, quality


def run_quality_checks(catalog: list[dict]) -> dict:
    accessions = [c["source_accession"] for c in catalog]
    full_by_sha: dict[str, list[str]] = defaultdict(list)
    region_by_sha: dict[str, list[str]] = defaultdict(list)
    for entry in catalog:
        if entry["full_sequence_sha256"]:
            full_by_sha[entry["full_sequence_sha256"]].append(entry["source_accession"])
        if entry["midpiwi_region_sha256"]:
            region_by_sha[entry["midpiwi_region_sha256"]].append(entry["source_accession"])

    checks = {
        "duplicate_accession": sorted(
            a for a, n in Counter(accessions).items() if n > 1
        ),
        "duplicate_full_sequence_groups": {
            k: v for k, v in full_by_sha.items() if len(v) > 1
        },
        "duplicate_midpiwi_region_groups": {
            k: v for k, v in region_by_sha.items() if len(v) > 1
        },
        "coord_out_of_range": [
            c["source_accession"] for c in catalog
            if c["midpiwi_extraction_status"] == "COORD_OUT_OF_RANGE"
        ],
        "mid_start_gt_mid_end": [
            c["source_accession"] for c in catalog
            if c["source_mid_start"] and c["source_mid_end"]
            and c["source_mid_start"] > c["source_mid_end"]
        ],
        "piwi_end_gt_current_length": [
            (c["source_accession"], c["source_piwi_end"], c["sequence_length_current"])
            for c in catalog
            if c["source_piwi_end"] and c["source_piwi_end"] > c["sequence_length_current"]
        ],
        "missing_mid": [
            c["source_accession"] for c in catalog
            if c["source_mid_start"] is None
        ],
        "length_current_ne_source": [
            (c["source_accession"], c["source_length"], c["sequence_length_current"])
            for c in catalog if c["length_matches_source"] is False
        ],
        "replacement_not_equivalent": [
            (c["source_accession"], c["replacement_accession"])
            for c in catalog
            if c["replacement_accession"] and c["replacement_equivalent"] is False
        ],
        "mid_then_piwi_gap_or_overlap": [
            (c["source_accession"], c["source_mid_end"], c["source_piwi_start"])
            for c in catalog
            if c["source_mid_end"] and c["source_piwi_start"]
            and c["source_mid_end"] + 1 != c["source_piwi_start"]
        ],
        "trun_with_complete_midpiwi": [
            c["source_accession"] for c in catalog
            if c["truncated_flag"] and c["midpiwi_extraction_status"] == "OK"
        ],
        "trun_with_incomplete_midpiwi": [
            c["source_accession"] for c in catalog
            if c["truncated_flag"] and c["midpiwi_extraction_status"] != "OK"
        ],
    }
    return checks


def summarize(catalog: list[dict], quality: dict) -> dict:
    def by(key: str) -> dict:
        return dict(Counter(c[key] for c in catalog))

    recovered = [c for c in catalog if c["full_sequence_sha256"]]
    midpiwi_ok = [c for c in catalog if c["midpiwi_extraction_status"] == "OK"]
    return {
        "n_total": len(catalog),
        "n_sequences_recovered": len(recovered),
        "n_sequences_not_recovered": len(catalog) - len(recovered),
        "ncbi_record_status": by("ncbi_record_status"),
        "n_unique_full_sequence_sha256": len(
            {c["full_sequence_sha256"] for c in recovered}
        ),
        "n_midpiwi_extracted": len(midpiwi_ok),
        "n_unique_midpiwi_region_sha256": len(
            {c["midpiwi_region_sha256"] for c in midpiwi_ok}
        ),
        "by_source_phylogenetic_clade": by("source_phylogenetic_clade"),
        "by_curated_pago_clade": by("curated_pago_clade"),
        "by_source_ago_type": by("source_ago_type"),
        "canonical_vs_truncated": {
            "canonical": sum(not c["truncated_flag"] for c in catalog),
            "truncated": sum(c["truncated_flag"] for c in catalog),
        },
        "by_architecture_status": by("architecture_status"),
        "by_curation_status": by("curation_status"),
        "by_midpiwi_extraction_status": by("midpiwi_extraction_status"),
        "experimental_anchors_in_table_s1": sorted(
            {c["experimental_anchor_name"] for c in catalog if c["experimental_anchor"]}
        ),
        "quarantine": [
            {"source_accession": c["source_accession"],
             "anchor": c["experimental_anchor_name"],
             "source_ago_type": c["source_ago_type"],
             "quarantine_reason": c["quarantine_reason"]}
            for c in catalog if c["curation_status"] == "QUARANTINE"
        ],
        "review": [
            {"source_accession": c["source_accession"],
             "anchor": c["experimental_anchor_name"],
             "curation_notes": c["curation_notes"]}
            for c in catalog if c["curation_status"] == "REVIEW"
        ],
        "replacements": [
            {"source": c["source_accession"],
             "replacement": c["replacement_accession"],
             "equivalent": c["replacement_equivalent"]}
            for c in catalog if c["replacement_accession"]
        ],
        "duplicate_full_sequence_pairs": [
            v for v in quality["duplicate_full_sequence_groups"].values()
        ],
        "duplicate_midpiwi_region_groups": [
            v for v in quality["duplicate_midpiwi_region_groups"].values()
        ],
        "coordinate_convention": quality["coordinate_convention"],
        "quality_checks": {
            k: v for k, v in quality.items()
            if k not in {
                "duplicate_full_sequence_groups",
                "duplicate_midpiwi_region_groups",
                "coordinate_convention",
            }
        },
        "potential_exclusion_reasons": {
            "suppressed": sum(c["ncbi_record_status"] == "SUPPRESSED" for c in catalog),
            "suppressed_and_truncated": sum(
                c["ncbi_record_status"] == "SUPPRESSED" and c["truncated_flag"]
                for c in catalog
            ),
            "dead_replaced": sum(
                c["ncbi_record_status"] == "DEAD_REPLACED" for c in catalog
            ),
            "obsolete_accession_double_listed": len([
                v for v in quality["duplicate_full_sequence_groups"].values()
            ]),
            "missing_mid_coordinate": len(quality["missing_mid"]),
            "midpiwi_extraction_failed": sum(
                c["midpiwi_extraction_status"] != "OK" for c in catalog
            ),
            "curated_clade_unresolved": sum(
                c["curated_pago_clade"] == "UNRESOLVED" for c in catalog
            ),
            "quarantine": sum(c["curation_status"] == "QUARANTINE" for c in catalog),
        },
    }


CATALOG_FIELDNAMES = [
    "source_accession", "resolved_accession", "species", "taxon_id",
    "source_length", "source_ago_type", "ago_family",
    "source_phylogenetic_clade", "curated_pago_clade", "architecture_status",
    "truncated_flag", "source_paz_start", "source_paz_end", "source_mid_start",
    "source_mid_end", "source_piwi_start", "source_piwi_end", "source_paz_type",
    "source_mid_5p_motif", "source_piwi_tetrad", "ncbi_record_status",
    "replacement_accession", "replacement_sequence_sha256", "replacement_equivalent",
    "full_sequence_sha256", "sequence_length_current", "length_matches_source",
    "midpiwi_sequence", "midpiwi_length", "midpiwi_region_sha256",
    "midpiwi_extraction_status", "experimental_anchor", "experimental_anchor_name",
    "curation_status", "quarantine_reason", "curation_notes", "proposed_partition",
]


def write_catalog_csv(*, catalog: list[dict], output_path: Path) -> str:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CATALOG_FIELDNAMES, extrasaction="ignore"
        )
        writer.writeheader()
        for entry in sorted(catalog, key=lambda z: z["source_accession"]):
            writer.writerow(entry)
    return _sha256_bytes(output_path.read_bytes())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--clade-seed-directory", type=Path, default=DEFAULT_CLADE_SEED_DIRECTORY
    )
    parser.add_argument(
        "--recall-reference",
        type=Path,
        default=PROJECT_ROOT / "tests" / "fixtures" / "query_recall_reference_set.csv",
    )
    parser.add_argument("--skip-fetch-sources", action="store_true")
    parser.add_argument(
        "--reuse-frozen", action="store_true",
        help="reuse the frozen NCBI resolution + sequence FASTA if present",
    )
    arguments = parser.parse_args(argv)

    clade_seed_directory = arguments.clade_seed_directory
    source_data_directory = clade_seed_directory / "source_data"
    session = _session()

    if not arguments.skip_fetch_sources:
        fetch_sources(source_data_directory=source_data_directory, session=session)

    records = read_table_s1(source_data_directory / "mbo006184236st1.xls")

    resolution_path = source_data_directory / "ryazansky_s1_ncbi_resolution.json"
    if arguments.reuse_frozen and resolution_path.is_file():
        status_by_accession = json.loads(resolution_path.read_text())
    else:
        status_by_accession = resolve_records(
            accessions=[r["source_accession"] for r in records], session=session
        )
        resolution_path.write_text(
            json.dumps(status_by_accession, indent=1, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    replacement_accessions = sorted(
        {v["replacedby"] for v in status_by_accession.values() if v.get("replacedby")}
    )
    sequence_by_accession = fetch_sequences(
        accessions=[r["source_accession"] for r in records] + replacement_accessions,
        fasta_path=source_data_directory / "ryazansky_s1_sequences.fasta",
        session=session,
        reuse_frozen=arguments.reuse_frozen,
    )

    catalog, quality = build_catalog(
        records=records,
        status_by_accession=status_by_accession,
        sequence_by_accession=sequence_by_accession,
        recall_reference_path=arguments.recall_reference,
    )
    catalog_path = clade_seed_directory / "ryazansky_s1_catalog.csv"
    catalog_sha256 = write_catalog_csv(catalog=catalog, output_path=catalog_path)
    summary = summarize(catalog, quality)
    summary["ryazansky_s1_catalog_sha256"] = catalog_sha256
    (clade_seed_directory / "ryazansky_s1_catalog_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
