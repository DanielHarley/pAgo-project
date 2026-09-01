"""Derive the MID-PIWI split groups for the clade reference layer (B4.3).

Same authoritative workflow as B2 (``scripts/derive_apaz_split_groups.py``):

  * MMseqs2 all-vs-all with ``--prefilter-mode 2`` (nofilter), full alignment,
    exporting ``qaln`` / ``taln`` (never ``pident`` / ``alnlen``);
  * identity and coverage recomputed by the canonical reducer
    (``src.pago_pipeline.apaz_split_groups``);
  * edge rule: ``identity >= 0.90 AND coverage >= 0.80``;
  * connected components, additionally unioning any sequences that share a
    ``midpiwi_region_sha256``.

The clustering unit is the **MID-PIWI region**, not the full protein.

Diagnostic-only passes at 1.00 / 0.90 / 0.70 / 0.50 identity (coverage fixed at
0.80) are reported but do NOT define split groups or leakage. Cross-clade at
0.70 / 0.50 is not treated as a conflict.

No partitioning. No TREE_BUILD / PLACEMENT_CALIBRATION / PLACEMENT_HOLDOUT. No
HMM. No tree. Does not touch B2 or B3.

Inputs  : resources/clade_seed/ryazansky_s1_catalog.csv  (rows with
          midpiwi_extraction_status == OK; 1002 regions)
Outputs : resources/clade_seed/split_groups/midpiwi/{split_groups.csv,
          pairs.filtered.tsv, diagnostic_thresholds.json, manifest.json,
          mmseqs_search.log}
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
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
    AlignedPair,
    build_split_groups,
    score_pairs,
)
from src.pago_pipeline.clade_midpiwi_split_groups import (  # noqa: E402
    RESOLVED_CLADES,
    assert_complete_all_vs_all,
    expected_unordered_pairs,
    label_split_group,
)

DEFAULT_CLADE_SEED_DIRECTORY = (
    PROJECT_ROOT / "src" / "pago_pipeline" / "resources" / "clade_seed"
)
SPLIT_GROUP_PREFIX = "MIDPIWI"
FORMAT_OUTPUT = (
    "query,target,qaln,taln,qlen,tlen,alnlen,nident,pident,fident,"
    "qstart,qend,tstart,tend,evalue,bits"
)
FORMAT_COLUMNS = FORMAT_OUTPUT.split(",")
DIAGNOSTIC_IDENTITY_THRESHOLDS = (1.00, 0.90, 0.70, 0.50)
ENVIRONMENT = {
    "ubuntu": "24.04.4", "wsl": "2.7.12.0", "conda": "26.3.2", "conda_env": "pago-linux",
    "pago_linux_explicit_lock_sha256":
        "aa362f3eabfadf86445061a2e01c6f903015138b4434b10aa18b18b2fe98ffcc",
    "pago_linux_env_locked_yml_sha256":
        "8a52c5f27d8c7eff016f90ccf5d08502b67b3ed234f8224b5e8ac809add04375",
}
_NON_DETERMINISTIC_LOG_MARKERS = (
    "Time for ", " ETA ", "k-mers per position", "k-mers per split",
    "Estimated memory", "Index table", "Process ", "rss ", "sys.", "user.",
)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _deterministic_mmseqs_log(log_text: str) -> bytes:
    kept: list[str] = []
    for raw_line in log_text.splitlines():
        line = raw_line.rstrip()
        if not line or line.lstrip().startswith(("[", "=")):
            continue
        if any(marker in line for marker in _NON_DETERMINISTIC_LOG_MARKERS):
            continue
        line = re.sub(r"tmp/\d+/", "tmp/<run>/", line)
        line = re.sub(r"--threads \d+", "--threads <n>", line)
        kept.append(line)
    return ("\n".join(kept) + "\n").encode("utf-8")


def read_midpiwi_regions(catalog_path: Path) -> tuple[dict[str, str], list[dict]]:
    rows = list(csv.DictReader(catalog_path.read_text(encoding="utf-8").splitlines()))
    ok_rows = [r for r in rows if r["midpiwi_extraction_status"] == "OK"]
    regions: dict[str, str] = {}
    for row in ok_rows:
        region = row["midpiwi_sequence"]
        if _sha256_bytes(region.encode("ascii")) != row["midpiwi_region_sha256"]:
            raise RuntimeError(
                f"midpiwi_region_sha256 mismatch for {row['source_accession']}."
            )
        regions[row["source_accession"]] = region
    return regions, rows


def _write_fasta(path: Path, identifier_to_sequence: dict[str, str]) -> None:
    with path.open("w", newline="\n") as handle:
        for identifier in sorted(identifier_to_sequence):
            handle.write(f">{identifier}\n")
            sequence = identifier_to_sequence[identifier]
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def _run_mmseqs_via_wsl(
    *, wsl_distro: str, fasta_path: Path, work_subdir: str
) -> tuple[str, str]:
    """Run createdb/search/convertalis in the pinned WSL bioconda env.

    Returns (pairs_raw_tsv_text, mmseqs_search_log_text).
    """
    posix_fasta = "/mnt/" + str(fasta_path).replace("\\", "/").replace(":", "").lower()[0:1] \
        + str(fasta_path).replace("\\", "/")[2:]
    script = f"""
set -e
source ~/miniforge3/etc/profile.d/conda.sh
conda activate pago-linux
W=~/pago-work/{work_subdir}
rm -rf "$W"; mkdir -p "$W"; cd "$W"
cp "{posix_fasta}" in.fasta
mmseqs createdb in.fasta DB > createdb.log 2>&1
mmseqs search DB DB alnDB tmp --prefilter-mode 2 --max-seqs 100000 -e 1e100 \
    --min-seq-id 0 -c 0 -a --alignment-mode 3 > mmseqs_search.log 2>&1
mmseqs convertalis DB DB alnDB pairs.raw.tsv --format-output "{FORMAT_OUTPUT}" > convertalis.log 2>&1
echo "===PAIRS==="
cat pairs.raw.tsv
echo "===LOG==="
cat mmseqs_search.log
"""
    completed = subprocess.run(
        ["wsl", "-d", wsl_distro, "-e", "bash", "-lc", script],
        capture_output=True, text=True,
    )
    output = completed.stdout.replace("\x00", "")
    if completed.returncode != 0 or "===PAIRS===" not in output:
        sys.stderr.write(output + completed.stderr)
        raise SystemExit("MMseqs2 run under WSL failed")
    pairs_text, _, log_text = output.partition("===LOG===")
    pairs_text = pairs_text.split("===PAIRS===", 1)[1]
    return pairs_text.strip("\n") + "\n", log_text.strip("\n") + "\n"


def _diagnostic_pass(
    scored_by_pair: dict, regions: dict[str, str], clade_by_id: dict[str, str],
    *, identity_threshold: float, coverage_threshold: float,
) -> dict:
    edges = [
        key for key, score in scored_by_pair.items()
        if score.pair_identity >= identity_threshold
        and score.pair_coverage >= coverage_threshold
    ]
    result = build_split_groups(
        sequence_ids=sorted(regions),
        sequence_sha256_by_id={
            i: _sha256_bytes(regions[i].encode("ascii")) for i in regions
        },
        edge_pairs=edges,
        split_group_id_prefix="DIAG",
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )
    cross_clade_members: list[dict] = []
    for members in result.members_by_split_group_id.values():
        resolved = {clade_by_id[m] for m in members} & set(RESOLVED_CLADES)
        if len(resolved) > 1:
            cross_clade_members.append({
                "n_members": len(members),
                "clade_breakdown": dict(sorted(
                    Counter(clade_by_id[m] for m in members).items()
                )),
                "members": {m: clade_by_id[m] for m in sorted(members)},
            })
    sizes = Counter(len(m) for m in result.members_by_split_group_id.values())
    return {
        "identity_threshold": identity_threshold,
        "coverage_threshold": coverage_threshold,
        "n_edges": len(edges),
        "n_split_groups": len(result.members_by_split_group_id),
        "n_singletons": sum(1 for m in result.members_by_split_group_id.values() if len(m) == 1),
        "n_multi_member": sum(1 for m in result.members_by_split_group_id.values() if len(m) > 1),
        "max_group_size": max(sizes),
        "size_distribution": dict(sorted(sizes.items())),
        "cross_clade_groups": len(cross_clade_members),
        # diagnostic evidence only - no quarantine, no exclusion
        "cross_clade_group_members": cross_clade_members,
    }


def derive(
    *, clade_seed_directory: Path, wsl_distro: str, reuse_pairs: Path | None = None,
) -> dict:
    catalog_path = clade_seed_directory / "ryazansky_s1_catalog.csv"
    regions, catalog_rows = read_midpiwi_regions(catalog_path)
    catalog_by_id = {r["source_accession"]: r for r in catalog_rows}
    clade_by_id = {a: catalog_by_id[a]["curated_pago_clade"] for a in regions}

    output_directory = clade_seed_directory / "split_groups" / "midpiwi"
    output_directory.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="clade_midpiwi_") as tmp:
        fasta_path = Path(tmp) / "midpiwi_regions.fasta"
        _write_fasta(fasta_path, regions)
        input_fasta_sha256 = _sha256_file(fasta_path)

        if reuse_pairs is not None:
            pairs_text = reuse_pairs.read_text(encoding="utf-8")
            log_text = (reuse_pairs.parent / "mmseqs_search.log").read_text(
                encoding="utf-8", errors="replace"
            )
        else:
            pairs_text, log_text = _run_mmseqs_via_wsl(
                wsl_distro=wsl_distro, fasta_path=fasta_path,
                work_subdir="clade_midpiwi",
            )

    pairs: list[AlignedPair] = []
    seen_unordered: set[tuple[str, str]] = set()
    raw_lines: list[str] = []
    for line in pairs_text.splitlines():
        stripped = line.rstrip("\n")
        if not stripped:
            continue
        raw_lines.append(stripped)
        row = dict(zip(FORMAT_COLUMNS, stripped.split("\t")))
        if row["query"] == row["target"]:
            continue
        seen_unordered.add(tuple(sorted((row["query"], row["target"]))))
        pairs.append(AlignedPair(
            row["query"], row["target"], row["qaln"], row["taln"],
            int(row["qlen"]), int(row["tlen"]),
        ))

    expected_unordered = expected_unordered_pairs(len(regions))
    combinatorial_coverage = len(seen_unordered) / expected_unordered
    # Refuse to build split groups on an incomplete all-vs-all.
    assert_complete_all_vs_all(
        observed_unordered_pairs=len(seen_unordered), n_sequences=len(regions)
    )

    pairs_raw_canonical_sha256 = _sha256_bytes(
        ("\n".join(sorted(set(raw_lines))) + "\n").encode("utf-8")
    )
    (output_directory / "mmseqs_search.log").write_bytes(
        _deterministic_mmseqs_log(log_text)
    )

    # best score per unordered pair (canonical reducer, never MMseqs pident/alnlen)
    best: dict[tuple[str, str], object] = {}
    for score in score_pairs(pairs):
        key = tuple(sorted((score.query_id, score.target_id)))
        current = best.get(key)
        if current is None or (score.pair_identity, score.pair_coverage) > (
            current.pair_identity, current.pair_coverage
        ):
            best[key] = score

    edges_90_80 = sorted(k for k, s in best.items() if s.is_edge)
    sequence_sha256_by_id = {
        i: _sha256_bytes(regions[i].encode("ascii")) for i in regions
    }
    result = build_split_groups(
        sequence_ids=sorted(regions),
        sequence_sha256_by_id=sequence_sha256_by_id,
        edge_pairs=edges_90_80,
        split_group_id_prefix=SPLIT_GROUP_PREFIX,
    )

    # per-group clade analysis
    region_sha_members: dict[str, list[str]] = defaultdict(list)
    for identifier, sha in sequence_sha256_by_id.items():
        region_sha_members[sha].append(identifier)
    groups_with_alias = sum(
        1 for members in result.members_by_split_group_id.values()
        if any(len(region_sha_members[sequence_sha256_by_id[m]]) > 1 for m in members)
    )

    group_rows = []
    counters = Counter()
    for group_id, members in result.members_by_split_group_id.items():
        member_clades = {clade_by_id[m] for m in members}
        label = label_split_group(
            member_ids=members, curated_pago_clade_by_id=clade_by_id
        )
        counters[label.group_label] += 1
        counters[f"status_{label.group_status}"] += 1
        if label.contains_unresolved:
            counters["contains_unresolved"] += 1
        if label.requires_manual_review:
            counters["requires_manual_review"] += 1
        group_rows.append({
            "split_group_id": group_id,
            "n_members": len(members),
            "members": " ".join(members),
            "member_clades": " ".join(sorted(member_clades)),
            "resolved_clades": " ".join(label.resolved_clades),
            "contains_unresolved": label.contains_unresolved,
            "quarantine_anchors": " ".join(label.quarantine_anchors),
            "requires_manual_review": label.requires_manual_review,
            "group_label": label.group_label,
            "group_status": label.group_status,
            "partition_eligible": label.partition_eligible,
        })

    # write split_groups.csv (per-sequence) + groups.csv (per-group) + pairs.filtered
    with (output_directory / "split_groups.csv").open("w", newline="\n") as handle:
        handle.write("sequence_id,midpiwi_region_sha256,split_group_id,curated_pago_clade\n")
        for identifier in sorted(regions):
            handle.write("%s,%s,%s,%s\n" % (
                identifier, sequence_sha256_by_id[identifier],
                result.split_group_id_by_sequence_id[identifier],
                clade_by_id[identifier],
            ))
    with (output_directory / "groups.csv").open("w", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(group_rows[0]))
        writer.writeheader()
        for row in sorted(group_rows, key=lambda z: z["split_group_id"]):
            writer.writerow(row)
    with (output_directory / "pairs.filtered.tsv").open("w", newline="\n") as handle:
        handle.write("seq_a\tseq_b\tidentity\tcoverage\tcomparable\tidentical\n")
        for left, right in edges_90_80:
            score = best[(left, right)]
            handle.write("%s\t%s\t%.6f\t%.6f\t%d\t%d\n" % (
                left, right, score.pair_identity, score.pair_coverage,
                score.comparable_columns, score.identical_columns,
            ))

    diagnostics = [
        _diagnostic_pass(
            best, regions, clade_by_id,
            identity_threshold=t, coverage_threshold=DEFAULT_COVERAGE_THRESHOLD,
        )
        for t in DIAGNOSTIC_IDENTITY_THRESHOLDS
    ]
    (output_directory / "diagnostic_thresholds.json").write_bytes(
        (json.dumps({
            "note": "diagnostic only; ONLY 0.90/0.80 defines split groups and "
                    "leakage; cross-clade at 0.70/0.50 is not a conflict",
            "coverage_threshold": DEFAULT_COVERAGE_THRESHOLD,
            "passes": diagnostics,
        }, indent=2, sort_keys=True) + "\n").encode("utf-8")
    )

    sizes = Counter(len(m) for m in result.members_by_split_group_id.values())
    clade_group_counts = Counter()
    for row in group_rows:
        if row["group_label"].startswith("SINGLE_") and row["group_label"] != "SINGLE_RESOLVED_PLUS_UNRESOLVED":
            clade_group_counts[row["group_label"].removeprefix("SINGLE_")] += 1

    manifest = {
        "artifact_type": "clade_midpiwi_split_groups",
        "clustering_unit": "mid_piwi_region",
        "n_regions": len(regions),
        "n_unique_regions": len(set(regions.values())),
        "n_split_groups": len(result.members_by_split_group_id),
        "n_singletons": sum(1 for m in result.members_by_split_group_id.values() if len(m) == 1),
        "n_multi_member_groups": sum(1 for m in result.members_by_split_group_id.values() if len(m) > 1),
        "max_group_size": max(sizes),
        "size_distribution": dict(sorted(sizes.items())),
        "groups_with_alias_or_exact_duplicate": groups_with_alias,
        "clade_group_counts": dict(clade_group_counts),
        "label_conflict_groups": counters["LABEL_CONFLICT"],
        "single_resolved_plus_unresolved_groups": counters["SINGLE_RESOLVED_PLUS_UNRESOLVED"],
        "unresolved_only_groups": counters["UNRESOLVED_ONLY"],
        "groups_containing_unresolved": counters["contains_unresolved"],
        "groups_requiring_manual_review": counters["requires_manual_review"],
        "partition_eligible_groups": sum(1 for r in group_rows if r["partition_eligible"]),
        "group_status_counts": {
            k.removeprefix("status_"): v for k, v in counters.items()
            if k.startswith("status_")
        },
        "edges_90_80": len(edges_90_80),
        "identity_threshold": DEFAULT_IDENTITY_THRESHOLD,
        "coverage_threshold": DEFAULT_COVERAGE_THRESHOLD,
        "split_group_method": SPLIT_GROUP_METHOD,
        "identity_definition": IDENTITY_DEFINITION,
        "coverage_definition": COVERAGE_DEFINITION,
        "connected_components_definition": CONNECTED_COMPONENTS_DEFINITION,
        "duplicate_rule": DUPLICATE_RULE,
        "all_vs_all_expected_unordered_pairs": expected_unordered,
        "all_vs_all_observed_unordered_pairs": len(seen_unordered),
        "all_vs_all_combinatorial_coverage": combinatorial_coverage,
        "mmseqs2": {"version": "18.8cc5c", "build": "hd6d6fdc_0", "channel": "bioconda"},
        "mmseqs_command": (
            "mmseqs createdb in.fasta DB ; "
            "mmseqs search DB DB alnDB tmp --prefilter-mode 2 --max-seqs 100000 -e 1e100 "
            "--min-seq-id 0 -c 0 -a --alignment-mode 3 ; "
            "mmseqs convertalis DB DB alnDB pairs.raw.tsv --format-output \"%s\"" % FORMAT_OUTPUT
        ),
        "environment": ENVIRONMENT,
        "input_fasta_sha256": input_fasta_sha256,
        "pairs_raw_tsv_canonical_sha256": pairs_raw_canonical_sha256,
        "pairs_raw_tsv_row_count": len(raw_lines),
        "split_groups_csv_sha256": _sha256_file(output_directory / "split_groups.csv"),
        "groups_csv_sha256": _sha256_file(output_directory / "groups.csv"),
        "pairs_filtered_tsv_sha256": _sha256_file(output_directory / "pairs.filtered.tsv"),
        "diagnostic_thresholds_json_sha256": _sha256_file(
            output_directory / "diagnostic_thresholds.json"
        ),
        "mmseqs_search_log_sha256": _sha256_file(output_directory / "mmseqs_search.log"),
        "source_catalog": "ryazansky_s1_catalog.csv",
        "source_catalog_sha256": _sha256_file(catalog_path),
        "canonical_library": "src/pago_pipeline/apaz_split_groups.py",
        "canonical_library_sha256": _sha256_file(
            PROJECT_ROOT / "src/pago_pipeline/apaz_split_groups.py"
        ),
        "derivation_script": "scripts/derive_clade_midpiwi_split_groups.py",
        "derivation_script_sha256": _sha256_file(Path(__file__).resolve()),
    }
    (output_directory / "manifest.json").write_bytes(
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clade-seed-directory", type=Path, default=DEFAULT_CLADE_SEED_DIRECTORY)
    parser.add_argument("--wsl-distro", type=str, default="Ubuntu-24.04")
    parser.add_argument("--reuse-pairs", type=Path, default=None,
                        help="path to a frozen pairs.raw.tsv (skips the WSL MMseqs2 run)")
    arguments = parser.parse_args(argv)
    manifest = derive(
        clade_seed_directory=arguments.clade_seed_directory,
        wsl_distro=arguments.wsl_distro,
        reuse_pairs=arguments.reuse_pairs,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
