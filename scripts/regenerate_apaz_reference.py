"""Deterministically regenerate the derived APAZ Phase B reference resources.

Offline, no network, no HMMs, no predictive metrics. Consumes exclusively:

  * the frozen Ryazansky Data Set S3 and the frozen Pfam HisG / EIIB seed
    alignments under ``<source>/source_data/``;
  * the frozen, validated split groups under ``<source>/split_groups/`` (produced
    once by ``scripts/derive_apaz_split_groups.py`` from an MMseqs2 all-vs-all run);
  * the canonical partition library ``src.pago_pipeline.apaz_split_groups``.

Produces, deterministically, into ``<output>``:

  * apaz_partitions.csv
  * apaz_Ia.sto apaz_Ib.sto apaz_IIa.sto apaz_IIb.sto apaz_III.sto
  * apaz_global.sto            (PyFAMSA 0.7.0 reconstruction of the BUILD set)
  * apaz_validation_sequences.fasta
  * seeds_lock.json

``pyfamsa==0.7.0`` is required (curation/development extra, not a runtime
dependency). Run ``python scripts/regenerate_apaz_reference.py`` to rewrite the
committed resources in place, or pass ``--output-directory`` to regenerate into
an empty directory for a clean-room reproduction check.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pago_pipeline.apaz_split_groups import (  # noqa: E402
    CONNECTED_COMPONENTS_DEFINITION,
    COVERAGE_DEFINITION,
    DEFAULT_COVERAGE_THRESHOLD,
    DEFAULT_IDENTITY_THRESHOLD,
    DUPLICATE_RULE,
    IDENTITY_DEFINITION,
    SPLIT_GROUP_METHOD,
    build_split_groups,
    partition_split_groups,
    validate_apaz_partition_invariants,
)

DEFAULT_RESOURCE_DIRECTORY = (
    PROJECT_ROOT / "src" / "pago_pipeline" / "resources" / "apaz_seed"
)

PARTITION_PROTOCOL = "apaz_partition_v2_mmseqs90_80"
PARTITION_SALT = "pago-phase-b-apaz-partitions-v2-mmseqs90_80"
CURATION_PROTOCOL_VERSION = "apaz_reference_curation_v2"
PFAM_RELEASE = "38.2"
SUBGROUPS = ("Ia", "Ib", "IIa", "IIb", "III")

# per subgroup: (FINAL_HOLDOUT, CALIBRATION, BUILD)  -> 70 / 15 / 15
APAZ_TARGETS = {
    "Ia": (12, 12, 59), "Ib": (16, 16, 77), "IIa": (15, 15, 68),
    "IIb": (27, 27, 125), "III": (2, 2, 8),
}
# per Pfam family: (FINAL_HOLDOUT, CALIBRATION)  BUILD is always 0
#   HisG: CALIBRATION 255 / FINAL_HOLDOUT 254 ; EIIB: CALIBRATION 259 / FINAL_HOLDOUT 258
NEGATIVE_TARGETS = {"PF01634": (254, 255), "PF00367": (258, 259)}
NEGATIVE_FAMILY_NAME = {"PF01634": "HisG", "PF00367": "EIIB"}
NEGATIVE_SEED_FILE = {"PF01634": "PF01634.seed.sto", "PF00367": "PF00367.seed.sto"}
NEGATIVE_DATASET = {"PF01634": "hisg", "PF00367": "eiib"}
SPLIT_GROUP_PREFIX = {"apaz": "APAZ", "hisg": "HISG", "eiib": "EIIB"}

ENVIRONMENT = {
    "ubuntu": "24.04.4", "wsl": "2.7.12.0", "conda": "26.3.2", "conda_env": "pago-linux",
    "pago_linux_explicit_lock_sha256":
        "aa362f3eabfadf86445061a2e01c6f903015138b4434b10aa18b18b2fe98ffcc",
    "pago_linux_env_locked_yml_sha256":
        "8a52c5f27d8c7eff016f90ccf5d08502b67b3ed234f8224b5e8ac809add04375",
}
MMSEQS2 = {"version": "18.8cc5c", "build": "hd6d6fdc_0", "channel": "bioconda"}
MMSEQS_COMMAND = (
    "mmseqs createdb in.fasta DB ; "
    "mmseqs search DB DB alnDB tmp --prefilter-mode 2 --max-seqs 100000 -e 1e100 "
    "--min-seq-id 0 -c 0 -a --alignment-mode 3 ; "
    "mmseqs convertalis DB DB alnDB pairs.raw.tsv --format-output "
    "query,target,qaln,taln,qlen,tlen,alnlen,nident,pident,fident,qstart,qend,tstart,tend,evalue,bits"
)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(text.encode("utf-8"))


def _ungap(sequence: str) -> str:
    return sequence.replace("-", "").replace(".", "")


# --------------------------------------------------------------------------- #
# frozen inputs                                                              #
# --------------------------------------------------------------------------- #
def _read_positive_records(source_directory: Path) -> dict[str, dict[str, str]]:
    text = (source_directory / "source_data" / "mbo006184236sd3.txt").read_bytes().decode(
        "utf-8-sig"
    )
    current_subgroup: str | None = None
    records: dict[str, dict[str, str]] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        heading = [g for g in SUBGROUPS if ("Group %s:" % g) in line]
        if heading:
            current_subgroup = heading[0]
            continue
        if line.startswith("#"):
            continue
        fields = line.split()
        identifier_parts = fields[0].split("|")
        records[identifier_parts[0]] = {
            "taxid": identifier_parts[1],
            "annotation": "|".join(identifier_parts[2:]),
            "subgroup": current_subgroup or "",
            "aligned": fields[1],
        }
    if len(records) != 481:
        raise RuntimeError(
            f"Ryazansky Data Set S3 must contain 481 APAZ records, found {len(records)}."
        )
    return records


def _read_pfam_seed(path: Path) -> dict[str, str]:
    fragments: dict[str, list[str]] = defaultdict(list)
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line or line == "//" or line.startswith("#"):
            continue
        fields = line.split()
        fragments[fields[0]].append(fields[1])
    return {identifier: "".join(parts) for identifier, parts in fragments.items()}


def _read_split_groups_resource(
    split_groups_directory: Path, dataset: str
) -> tuple[dict[str, str], list[tuple[str, str]]]:
    sequence_sha256_by_id: dict[str, str] = {}
    with (split_groups_directory / dataset / "split_groups.csv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle):
            sequence_sha256_by_id[row["sequence_id"]] = row["sequence_sha256"]
    edges: list[tuple[str, str]] = []
    with (split_groups_directory / dataset / "pairs.filtered.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            edges.append(tuple(sorted((row["seq_a"], row["seq_b"]))))
    return sequence_sha256_by_id, sorted(set(edges))


def _committed_split_group_membership(
    split_groups_directory: Path, dataset: str
) -> dict[str, tuple[str, ...]]:
    members: dict[str, list[str]] = defaultdict(list)
    with (split_groups_directory / dataset / "split_groups.csv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle):
            members[row["split_group_id"]].append(row["sequence_id"])
    return {group: tuple(sorted(ids)) for group, ids in members.items()}


# --------------------------------------------------------------------------- #
# writers                                                                    #
# --------------------------------------------------------------------------- #
def _stockholm_text(model_id: str, aligned_by_accession: dict[str, str], description: str) -> str:
    widths = {len(sequence) for sequence in aligned_by_accession.values()}
    if len(widths) != 1:
        raise RuntimeError(f"Seed {model_id} has unequal alignment widths: {sorted(widths)}.")
    identifier_width = max(len(accession) for accession in aligned_by_accession)
    rows = [
        "# STOCKHOLM 1.0",
        "#=GF ID %s" % model_id,
        "#=GF DE %s" % description,
        "#=GF AU Ryazansky et al. 2018 source curation; pAgo-project Phase B",
    ]
    for accession in sorted(aligned_by_accession):
        rows.append("%s %s" % (accession.ljust(identifier_width), aligned_by_accession[accession]))
    rows.append("//")
    return "\n".join(rows) + "\n"


def _fasta_text(records: list[tuple[str, str, str]]) -> str:
    rows: list[str] = []
    for sequence_id, description, sequence in sorted(records, key=lambda record: record[0]):
        rows.append((">%s %s" % (sequence_id, description)).rstrip())
        for start in range(0, len(sequence), 80):
            rows.append(sequence[start:start + 80])
    return "\n".join(rows) + "\n"


def _famsa_global_alignment(accession_to_ungapped: dict[str, str]) -> tuple[dict[str, str], str]:
    try:
        import pyfamsa
    except ImportError as error:
        raise RuntimeError(
            "pyfamsa==0.7.0 is required to reconstruct apaz_global.sto. It is a "
            "curation/development extra, not a runtime dependency: "
            "pip install 'pago-project[curation]' (or pip install pyfamsa==0.7.0)."
        ) from error
    sequences = [
        pyfamsa.Sequence(accession.encode("ascii"), accession_to_ungapped[accession].encode("ascii"))
        for accession in sorted(accession_to_ungapped)
    ]
    alignment = pyfamsa.Aligner(guide_tree="upgma").align(sequences)
    aligned = {
        sequence.id.decode("ascii"): sequence.sequence.decode("ascii")
        for sequence in alignment
    }
    if set(aligned) != set(accession_to_ungapped):
        raise RuntimeError("PyFAMSA changed the accession inventory of the global BUILD alignment.")
    return aligned, str(pyfamsa.__version__)


# --------------------------------------------------------------------------- #
# main                                                                       #
# --------------------------------------------------------------------------- #
def regenerate(*, source_directory: Path, output_directory: Path) -> dict[str, str]:
    source_directory = source_directory.resolve()
    output_directory = output_directory.resolve()
    source_data = source_directory / "source_data"
    split_groups_directory = source_directory / "split_groups"
    output_directory.mkdir(parents=True, exist_ok=True)

    positives = _read_positive_records(source_directory)
    subgroup_by_accession = {a: positives[a]["subgroup"] for a in positives}
    ungapped_by_accession = {a: _ungap(positives[a]["aligned"]) for a in positives}
    positive_sha256 = {
        a: _sha256_bytes(ungapped_by_accession[a].encode("ascii")) for a in positives
    }

    negative_seeds = {
        family: _read_pfam_seed(source_data / NEGATIVE_SEED_FILE[family])
        for family in NEGATIVE_SEED_FILE
    }

    # ---- split groups: rebuild from the frozen edge list, prove self-consistency ----
    split_group_result: dict[str, object] = {}
    for dataset in ("apaz", "hisg", "eiib"):
        sequence_sha256_by_id, edges = _read_split_groups_resource(
            split_groups_directory, dataset
        )
        result = build_split_groups(
            sequence_ids=sorted(sequence_sha256_by_id),
            sequence_sha256_by_id=sequence_sha256_by_id,
            edge_pairs=edges,
            split_group_id_prefix=SPLIT_GROUP_PREFIX[dataset],
        )
        committed = _committed_split_group_membership(split_groups_directory, dataset)
        if result.members_by_split_group_id != committed:
            raise RuntimeError(
                f"{dataset}: the committed split_groups.csv does not match what the "
                "canonical library rebuilds from pairs.filtered.tsv."
            )
        split_group_result[dataset] = (result, sequence_sha256_by_id)

    apaz_result, _ = split_group_result["apaz"]
    if any(
        split_group_result["apaz"][1][a] != positive_sha256[a] for a in positives
    ):
        raise RuntimeError("APAZ split-group sequence_sha256 disagrees with the S3 source.")

    # ---- deterministic stratified v2 partition ----
    apaz_partition = partition_split_groups(
        members_by_split_group_id=apaz_result.members_by_split_group_id,
        stratum_by_sequence_id=subgroup_by_accession,
        fill_order_by_stratum={
            s: [
                ("FINAL_HOLDOUT", APAZ_TARGETS[s][0]),
                ("CALIBRATION", APAZ_TARGETS[s][1]),
                ("BUILD", APAZ_TARGETS[s][2]),
            ]
            for s in SUBGROUPS
        },
        partition_salt=PARTITION_SALT,
    )
    negative_partition: dict[str, object] = {}
    for family, (holdout_target, calibration_target) in NEGATIVE_TARGETS.items():
        dataset = NEGATIVE_DATASET[family]
        result, sequence_sha256_by_id = split_group_result[dataset]
        negative_partition[family] = partition_split_groups(
            members_by_split_group_id=result.members_by_split_group_id,
            stratum_by_sequence_id={sid: family for sid in sequence_sha256_by_id},
            fill_order_by_stratum={
                family: [
                    ("FINAL_HOLDOUT", holdout_target),
                    ("CALIBRATION", calibration_target),
                ]
            },
            partition_salt=PARTITION_SALT,
        )

    # ---- apaz_partitions.csv ----
    rows: list[dict[str, str]] = []
    for accession in sorted(positives):
        group_id = apaz_result.split_group_id_by_sequence_id[accession]
        rows.append({
            "accession": accession,
            "homology_cluster_id": group_id,
            "split_group_id": group_id,
            "partition": apaz_partition.partition_by_sequence_id[accession],
            "class_label": "APAZ",
            "subgroup": subgroup_by_accession[accession],
            "evidence_level": "LITERATURE_PHYLOGENETIC_SET",
            "source": "Ryazansky_et_al_2018_Data_Set_S3",
            "source_sequence_id": "%s|%s|%s" % (
                accession, positives[accession]["taxid"], positives[accession]["annotation"]
            ),
            "sequence_sha256": positive_sha256[accession],
        })
    for family in NEGATIVE_TARGETS:
        dataset = NEGATIVE_DATASET[family]
        result, sequence_sha256_by_id = split_group_result[dataset]
        partition = negative_partition[family]
        for seed_id in sorted(sequence_sha256_by_id):
            group_id = result.split_group_id_by_sequence_id[seed_id]
            negative_accession = "NEG_%s_%s" % (
                family, hashlib.sha256(seed_id.encode("utf-8")).hexdigest()[:20]
            )
            rows.append({
                "accession": negative_accession,
                "homology_cluster_id": group_id,
                "split_group_id": group_id,
                "partition": partition.partition_by_sequence_id[seed_id],
                "class_label": "HARD_NEGATIVE",
                "subgroup": "",
                "evidence_level": "LITERATURE_EXCLUDED_FALSE_HOMOLOG_FAMILY",
                "source": "Pfam_%s_%s_seed" % (PFAM_RELEASE, family),
                "source_sequence_id": seed_id,
                "sequence_sha256": sequence_sha256_by_id[seed_id],
            })
    rows.sort(key=lambda row: row["accession"])
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=[
            "accession", "homology_cluster_id", "split_group_id", "partition",
            "class_label", "subgroup", "evidence_level", "source",
            "source_sequence_id", "sequence_sha256",
        ],
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    _write(output_directory / "apaz_partitions.csv", buffer.getvalue())

    # ---- subgroup seeds (S3 columns preserved) ----
    build_alignment_by_subgroup: dict[str, dict[str, str]] = defaultdict(dict)
    for accession in positives:
        if apaz_partition.partition_by_sequence_id[accession] == "BUILD":
            build_alignment_by_subgroup[subgroup_by_accession[accession]][accession] = (
                positives[accession]["aligned"]
            )
    for subgroup in SUBGROUPS:
        _write(
            output_directory / ("apaz_%s.sto" % subgroup),
            _stockholm_text(
                "APAZ_%s" % subgroup.upper(),
                build_alignment_by_subgroup[subgroup],
                "APAZ subgroup %s BUILD alignment from Ryazansky et al. 2018" % subgroup,
            ),
        )

    # ---- global BUILD alignment (independent PyFAMSA reconstruction) ----
    build_ungapped = {
        accession: ungapped_by_accession[accession]
        for accession in positives
        if apaz_partition.partition_by_sequence_id[accession] == "BUILD"
    }
    global_alignment, pyfamsa_version = _famsa_global_alignment(build_ungapped)
    _write(
        output_directory / "apaz_global.sto",
        _stockholm_text(
            "APAZ_GLOBAL",
            global_alignment,
            "Global APAZ BUILD alignment reconstructed with PyFAMSA from the five "
            "Ryazansky 2018 subgroup alignments",
        ),
    )

    # ---- validation FASTA (non-BUILD positives + all hard negatives) ----
    validation_records: list[tuple[str, str, str]] = []
    for accession in positives:
        partition = apaz_partition.partition_by_sequence_id[accession]
        if partition != "BUILD":
            validation_records.append((
                accession,
                "class=APAZ subgroup=%s partition=%s" % (subgroup_by_accession[accession], partition),
                ungapped_by_accession[accession],
            ))
    for family in NEGATIVE_TARGETS:
        dataset = NEGATIVE_DATASET[family]
        _, sequence_sha256_by_id = split_group_result[dataset]
        partition = negative_partition[family]
        for seed_id in sequence_sha256_by_id:
            negative_accession = "NEG_%s_%s" % (
                family, hashlib.sha256(seed_id.encode("utf-8")).hexdigest()[:20]
            )
            validation_records.append((
                negative_accession,
                "class=HARD_NEGATIVE family=%s partition=%s"
                % (family, partition.partition_by_sequence_id[seed_id]),
                _ungap(negative_seeds[family][seed_id]),
            ))
    validation_fasta = _fasta_text(validation_records)
    _write(output_directory / "apaz_validation_sequences.fasta", validation_fasta)

    # ---- seeds_lock.json ----
    lock = _build_lock(
        source_directory=source_directory,
        split_groups_directory=split_groups_directory,
        output_directory=output_directory,
        rows=rows,
        negative_seeds=negative_seeds,
        pyfamsa_version=pyfamsa_version,
        validation_record_count=len(validation_records),
    )
    _write(
        output_directory / "seeds_lock.json",
        json.dumps(lock, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )

    return {
        name: _sha256_file(output_directory / name)
        for name in (
            "apaz_partitions.csv", "apaz_global.sto", "apaz_Ia.sto", "apaz_Ib.sto",
            "apaz_IIa.sto", "apaz_IIb.sto", "apaz_III.sto",
            "apaz_validation_sequences.fasta", "seeds_lock.json",
        )
    }


def _build_lock(
    *,
    source_directory: Path,
    split_groups_directory: Path,
    output_directory: Path,
    rows: list[dict[str, str]],
    negative_seeds: dict[str, dict[str, str]],
    pyfamsa_version: str,
    validation_record_count: int,
) -> dict[str, object]:
    source_data = source_directory / "source_data"

    seed_entries = []
    lock_subgroup = {
        "apaz_global.sto": "GLOBAL", "apaz_Ia.sto": "IA", "apaz_Ib.sto": "IB",
        "apaz_IIa.sto": "IIA", "apaz_IIb.sto": "IIB", "apaz_III.sto": "III",
    }
    for file_name, subgroup in sorted(lock_subgroup.items()):
        path = output_directory / file_name
        sequence_rows = [
            line for line in path.read_text().splitlines()
            if line and not line.startswith("#") and line != "//"
        ]
        seed_entries.append({
            "path": file_name, "subgroup": subgroup, "alphabet": "amino",
            "sequence_count": len(sequence_rows),
            "alignment_length": len(sequence_rows[0].split()[1]),
            "sha256": _sha256_file(path),
        })

    negative_source_entries = []
    for family, seed_file in NEGATIVE_SEED_FILE.items():
        path = source_data / seed_file
        text = path.read_text(encoding="utf-8-sig")
        versioned_accession = re.search(r"^#=GF AC\s+(\S+)", text, re.MULTILINE).group(1)
        seed_name = re.search(r"^#=GF ID\s+(\S+)", text, re.MULTILINE).group(1)
        negative_source_entries.append({
            "family_accession": family,
            "family_name": NEGATIVE_FAMILY_NAME[family],
            "pfam_seed_name": seed_name,
            "versioned_family_accession": versioned_accession,
            "path": "source_data/%s" % seed_file,
            "sequence_count": len(negative_seeds[family]),
            "sha256": _sha256_file(path),
        })

    split_group_resources = {}
    for dataset in ("apaz", "hisg", "eiib"):
        directory = split_groups_directory / dataset
        split_group_resources[dataset] = {
            key: {
                "path": "split_groups/%s/%s" % (dataset, file_name),
                "sha256": _sha256_file(directory / file_name),
            }
            for key, file_name in (
                ("split_groups_csv", "split_groups.csv"),
                ("pairs_filtered_tsv", "pairs.filtered.tsv"),
                ("manifest_json", "manifest.json"),
                ("mmseqs_search_log", "mmseqs_search.log"),
            )
        }

    partition_counts_by_class: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in rows:
        partition_counts_by_class[row["class_label"]][row["partition"]] += 1

    return {
        "artifact_type": "apaz_seed_lock",
        "lock_format_version": "2.0",
        "curation_protocol_version": CURATION_PROTOCOL_VERSION,
        "partition_protocol": PARTITION_PROTOCOL,
        "partition_salt": PARTITION_SALT,
        "supersedes": {
            "curation_protocol_version": "apaz_reference_curation_v1",
            "partition_protocol": "apaz_partition_v1",
            "note": "v1 keyed the split unit on a hash of the accession and was never "
                    "used for predictive evaluation; superseded by the MMseqs2 90/80 "
                    "split-group unit.",
        },
        "split_group_method": SPLIT_GROUP_METHOD,
        "identity_threshold": DEFAULT_IDENTITY_THRESHOLD,
        "coverage_threshold": DEFAULT_COVERAGE_THRESHOLD,
        "identity_definition": IDENTITY_DEFINITION,
        "coverage_definition": COVERAGE_DEFINITION,
        "connected_components_definition": CONNECTED_COMPONENTS_DEFINITION,
        "duplicate_sha256_rule": DUPLICATE_RULE,
        "high_similarity_control_note": (
            "A 90/80 split group is a high-similarity / redundancy control unit, "
            "NOT a claim that homology is absent below 90% identity."
        ),
        "mmseqs2": MMSEQS2,
        "mmseqs_command": MMSEQS_COMMAND,
        "environment": ENVIRONMENT,
        "canonical_library": {
            "module": "src/pago_pipeline/apaz_split_groups.py",
            "sha256": _sha256_file(PROJECT_ROOT / "src/pago_pipeline/apaz_split_groups.py"),
        },
        "regeneration_script": {
            "script": "scripts/regenerate_apaz_reference.py",
            "sha256": _sha256_file(Path(__file__).resolve()),
        },
        "partition_policy": {
            "positive_partition_unit": "mmseqs2_90_80_split_group",
            "positive_partition_fractions": {"BUILD": 0.70, "CALIBRATION": 0.15, "FINAL_HOLDOUT": 0.15},
            "positive_target_counts_by_subgroup": {
                s: {"FINAL_HOLDOUT": APAZ_TARGETS[s][0], "CALIBRATION": APAZ_TARGETS[s][1],
                    "BUILD": APAZ_TARGETS[s][2]}
                for s in SUBGROUPS
            },
            "negative_target_counts_by_family": {
                family: {"FINAL_HOLDOUT": NEGATIVE_TARGETS[family][0],
                         "CALIBRATION": NEGATIVE_TARGETS[family][1], "BUILD": 0}
                for family in NEGATIVE_TARGETS
            },
            "stratification": (
                "positives by Ryazansky APAZ subgroup Ia/Ib/IIa/IIb/III; negatives by Pfam family"
            ),
            "fill_order": ["FINAL_HOLDOUT", "CALIBRATION", "BUILD"],
            "selection_order_key": "sha256(PARTITION_SALT | stratum | split_group_id | partition)",
            "subset_sum_policy": (
                "deterministic DFS over whole split groups in selection-order, "
                "INCLUDE-before-EXCLUDE, prune on overshoot and on insufficient suffix "
                "sum; first exact hit wins; abort (no silent adjustment) if no exact "
                "combination exists"
            ),
            "tie_break": "lexicographic on the selection_order_key sha256",
            "duplicate_grouping": "sequence_sha256 union applied before any alignment edge",
            "row_order_independence": (
                "build_split_groups and partition_split_groups are invariant to input order"
            ),
            "positive_source_alignment_columns": (
                "subgroup seeds keep the Ryazansky S3 columns; the global BUILD MSA is "
                "independently reconstructed with PyFAMSA"
            ),
        },
        "partition_counts": {
            label: dict(sorted(counts.items()))
            for label, counts in sorted(partition_counts_by_class.items())
        },
        "global_alignment_builder": {
            "algorithm": "FAMSA 2 through PyFAMSA", "guide_tree": "upgma",
            "pyfamsa_version": pyfamsa_version,
        },
        "files": sorted(seed_entries, key=lambda entry: entry["path"]),
        "negative_source_files": negative_source_entries,
        "source_alignment": {
            "article_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6299218/",
            "license": "CC-BY-4.0",
            "md5": "ac30768ca799feb85165cea9d11aa3e2",
            "path": "source_data/mbo006184236sd3.txt",
            "sha256": _sha256_file(source_data / "mbo006184236sd3.txt"),
        },
        "source_metadata": {
            "path": "source_data/PMC6299218.1.metadata.json",
            "sha256": _sha256_file(source_data / "PMC6299218.1.metadata.json"),
        },
        "pfam_release": PFAM_RELEASE,
        "curation_notes": {
            "path": "curation_notes.md",
            "sha256": _sha256_file(source_directory / "curation_notes.md"),
        },
        "partitions_file": {
            "path": "apaz_partitions.csv",
            "sha256": _sha256_file(output_directory / "apaz_partitions.csv"),
        },
        "validation_sequences_file": {
            "path": "apaz_validation_sequences.fasta",
            "sequence_count": validation_record_count,
            "sha256": _sha256_file(output_directory / "apaz_validation_sequences.fasta"),
        },
        "split_group_resources": split_group_resources,
    }


def _summary(output_directory: Path) -> None:
    with (output_directory / "apaz_partitions.csv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    positives = [r for r in rows if r["class_label"] == "APAZ"]
    negatives = [r for r in rows if r["class_label"] == "HARD_NEGATIVE"]
    print("apaz_partitions.csv rows:", len(rows))
    print("  positives:", dict(sorted(Counter(r["partition"] for r in positives).items())))
    print("  negatives:", dict(sorted(Counter(r["partition"] for r in negatives).items())))
    print("  split groups:", len({r["split_group_id"] for r in rows}))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", type=Path, default=DEFAULT_RESOURCE_DIRECTORY,
                        help="frozen inputs: source_data/, split_groups/, curation_notes.md")
    parser.add_argument("--output-directory", type=Path, default=None,
                        help="where to write the derived resources (default: --source-directory)")
    arguments = parser.parse_args(argv)
    output_directory = arguments.output_directory or arguments.source_directory
    hashes = regenerate(
        source_directory=arguments.source_directory,
        output_directory=output_directory,
    )
    validate_apaz_partition_invariants(
        partitions_csv_path=output_directory / "apaz_partitions.csv",
        split_groups_directory=arguments.source_directory / "split_groups",
    )
    _summary(output_directory)
    print("\nregenerated file sha256:")
    for name, digest in hashes.items():
        print("  %-34s %s" % (name, digest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
