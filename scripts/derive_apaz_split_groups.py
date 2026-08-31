"""Derive the frozen APAZ / HisG / EIIB split-group resources from the frozen
sources with an MMseqs2 all-vs-all alignment and the canonical library.

Maintenance / curation step, run outside CI in the pinned external-tool
environment (WSL 2 + bioconda MMseqs2, see
``src/pago_pipeline/resources/apaz_seed/seeds_lock.json``). It performs no
partitioning and builds no HMMs.

Inputs  : resources/apaz_seed/source_data/{mbo006184236sd3.txt, PF01634.seed.sto,
          PF00367.seed.sto}
Outputs : resources/apaz_seed/split_groups/<dataset>/{split_groups.csv,
          pairs.filtered.tsv, manifest.json, mmseqs_search.log}

The edge rule (identity over comparable columns, coverage over the shorter
sequence, 0.90 / 0.80) and the connected-component construction come exclusively
from ``src.pago_pipeline.apaz_split_groups``. MMseqs2 ``pident`` / ``alnlen`` are
never used to decide an edge.
"""
from __future__ import annotations

import argparse
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

DEFAULT_RESOURCE_DIRECTORY = (
    PROJECT_ROOT / "src" / "pago_pipeline" / "resources" / "apaz_seed"
)
FORMAT_OUTPUT = (
    "query,target,qaln,taln,qlen,tlen,alnlen,nident,pident,fident,"
    "qstart,qend,tstart,tend,evalue,bits"
)
FORMAT_COLUMNS = FORMAT_OUTPUT.split(",")
SUBGROUPS = ("Ia", "Ib", "IIa", "IIb", "III")
DATASETS = (
    ("apaz", "APAZ", 481),
    ("hisg", "HISG", 509),
    ("eiib", "EIIB", 517),
)
NEGATIVE_SEED_FILE = {"hisg": "PF01634.seed.sto", "eiib": "PF00367.seed.sto"}
ENVIRONMENT = {
    "ubuntu": "24.04.4", "wsl": "2.7.12.0", "conda": "26.3.2", "conda_env": "pago-linux",
    "pago_linux_explicit_lock_sha256":
        "aa362f3eabfadf86445061a2e01c6f903015138b4434b10aa18b18b2fe98ffcc",
    "pago_linux_env_locked_yml_sha256":
        "8a52c5f27d8c7eff016f90ccf5d08502b67b3ed234f8224b5e8ac809add04375",
}


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _ungap(sequence: str) -> str:
    return sequence.replace("-", "").replace(".", "")


_NON_DETERMINISTIC_LOG_MARKERS = (
    "Time for ", " ETA ", "k-mers per position", "k-mers per split",
    "Estimated memory", "Index table", "Process ", "rss ", "sys.", "user.",
)


def _deterministic_mmseqs_log(log_text: str) -> bytes:
    """Keep the reproducible MMseqs2 parameter echo; drop timings / progress bars.

    Preserves the ``key<tab>value`` configuration block (MMseqs Version,
    Prefilter mode, Alignment mode, ...) that documents how the alignment was
    run, and drops wall-clock timing, ETAs and progress bars so the committed
    log has a stable SHA-256.
    """
    kept: list[str] = []
    for raw_line in log_text.splitlines():
        line = raw_line.rstrip()
        if not line:
            continue
        if line.lstrip().startswith(("[", "=")):
            continue
        if any(marker in line for marker in _NON_DETERMINISTIC_LOG_MARKERS):
            continue
        # normalise machine/run-specific tokens that do not affect the result
        line = re.sub(r"tmp/\d+/", "tmp/<run>/", line)
        line = re.sub(r"--threads \d+", "--threads <n>", line)
        kept.append(line)
    return ("\n".join(kept) + "\n").encode("utf-8")


def _read_positive_ungapped(resource_directory: Path) -> dict[str, str]:
    text = (resource_directory / "source_data" / "mbo006184236sd3.txt").read_bytes().decode(
        "utf-8-sig"
    )
    ungapped: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "Group " in line:
            continue
        fields = line.split()
        ungapped[fields[0].split("|")[0]] = _ungap(fields[1])
    return ungapped


def _read_seed_ungapped(path: Path) -> dict[str, str]:
    fragments: dict[str, list[str]] = defaultdict(list)
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line or line == "//" or line.startswith("#"):
            continue
        fields = line.split()
        fragments[fields[0]].append(fields[1])
    return {identifier: _ungap("".join(parts)) for identifier, parts in fragments.items()}


def _write_fasta(path: Path, identifier_to_sequence: dict[str, str]) -> None:
    with path.open("w", newline="\n") as handle:
        for identifier in sorted(identifier_to_sequence):
            handle.write(">%s\n" % identifier)
            sequence = identifier_to_sequence[identifier]
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def _run_mmseqs(mmseqs: str, work_directory: Path, fasta_path: Path) -> tuple[Path, Path]:
    def run(arguments: list[str], log_path: Path | None = None) -> None:
        result = subprocess.run(
            [mmseqs, *arguments], cwd=work_directory, capture_output=True, text=True
        )
        if log_path is not None:
            log_path.write_text(result.stdout + result.stderr)
        if result.returncode != 0:
            sys.stderr.write(result.stdout + result.stderr)
            raise SystemExit("mmseqs %s failed" % arguments[0])

    (work_directory / "in.fasta").write_bytes(fasta_path.read_bytes())
    search_log = work_directory / "mmseqs_search.log"
    run(["createdb", "in.fasta", "DB"])
    run(
        [
            "search", "DB", "DB", "alnDB", "tmp",
            "--prefilter-mode", "2", "--max-seqs", "100000", "-e", "1e100",
            "--min-seq-id", "0", "-c", "0", "-a", "--alignment-mode", "3",
        ],
        log_path=search_log,
    )
    run(["convertalis", "DB", "DB", "alnDB", "pairs.raw.tsv", "--format-output", FORMAT_OUTPUT])
    return work_directory / "pairs.raw.tsv", search_log


def derive_dataset(
    *,
    dataset: str,
    prefix: str,
    expected_count: int,
    ungapped_by_id: dict[str, str],
    resource_directory: Path,
    mmseqs: str,
) -> dict[str, object]:
    if len(ungapped_by_id) != expected_count:
        raise RuntimeError(
            f"{dataset}: expected {expected_count} sequences, found {len(ungapped_by_id)}."
        )
    sequence_ids = sorted(ungapped_by_id)
    sequence_sha256_by_id = {
        identifier: _sha256_bytes(ungapped_by_id[identifier].encode("ascii"))
        for identifier in sequence_ids
    }

    output_directory = resource_directory / "split_groups" / dataset
    output_directory.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"mmseqs_{dataset}_") as temporary_name:
        work_directory = Path(temporary_name)
        fasta_path = work_directory / "input.fasta"
        _write_fasta(fasta_path, ungapped_by_id)
        input_fasta_sha256 = _sha256_file(fasta_path)
        pairs_raw_path, search_log_path = _run_mmseqs(mmseqs, work_directory, fasta_path)

        pairs: list[AlignedPair] = []
        seen_unordered: set[tuple[str, str]] = set()
        raw_lines: list[str] = []
        with pairs_raw_path.open("r", encoding="utf-8") as handle:
            for line in handle:
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
        expected_unordered = expected_count * (expected_count - 1) // 2
        missing = expected_unordered - len(seen_unordered)
        if missing:
            raise RuntimeError(
                f"{dataset}: MMseqs2 all-vs-all is incomplete "
                f"({len(seen_unordered)}/{expected_unordered} unordered pairs)."
            )
        # MMseqs2 is multithreaded: the pairs.raw.tsv row *order* is not stable,
        # so hash the canonicalised (sorted, deduplicated) row set for provenance.
        pairs_raw_canonical_sha256 = _sha256_bytes(
            ("\n".join(sorted(set(raw_lines))) + "\n").encode("utf-8")
        )
        pairs_raw_row_count = len(raw_lines)
        (output_directory / "mmseqs_search.log").write_bytes(
            _deterministic_mmseqs_log(search_log_path.read_text(encoding="utf-8", errors="replace"))
        )

    best: dict[tuple[str, str], object] = {}
    for score in score_pairs(pairs):
        key = tuple(sorted((score.query_id, score.target_id)))
        current = best.get(key)
        if current is None or (score.pair_identity, score.pair_coverage) > (
            current.pair_identity, current.pair_coverage
        ):
            best[key] = score
    edges = sorted(key for key, score in best.items() if score.is_edge)

    result = build_split_groups(
        sequence_ids=sequence_ids,
        sequence_sha256_by_id=sequence_sha256_by_id,
        edge_pairs=edges,
        split_group_id_prefix=prefix,
    )

    with (output_directory / "split_groups.csv").open("w", newline="\n") as handle:
        handle.write("sequence_id,sequence_sha256,split_group_id\n")
        for identifier in sequence_ids:
            handle.write("%s,%s,%s\n" % (
                identifier, sequence_sha256_by_id[identifier],
                result.split_group_id_by_sequence_id[identifier],
            ))
    with (output_directory / "pairs.filtered.tsv").open("w", newline="\n") as handle:
        handle.write("seq_a\tseq_b\tidentity\tcoverage\tcomparable\tidentical\n")
        for left, right in edges:
            score = best[(left, right)]
            handle.write("%s\t%s\t%.6f\t%.6f\t%d\t%d\n" % (
                left, right, score.pair_identity, score.pair_coverage,
                score.comparable_columns, score.identical_columns,
            ))

    sizes = Counter(len(members) for members in result.members_by_split_group_id.values())
    manifest = {
        "dataset": dataset,
        "n_sequences": len(sequence_ids),
        "n_split_groups": len(result.members_by_split_group_id),
        "size_distribution": dict(sorted(sizes.items())),
        "max_group_size": max(sizes),
        "multi_member_groups": sum(
            1 for members in result.members_by_split_group_id.values() if len(members) > 1
        ),
        "edges_90_80": len(edges),
        "identity_threshold": DEFAULT_IDENTITY_THRESHOLD,
        "coverage_threshold": DEFAULT_COVERAGE_THRESHOLD,
        "split_group_method": SPLIT_GROUP_METHOD,
        "identity_definition": IDENTITY_DEFINITION,
        "coverage_definition": COVERAGE_DEFINITION,
        "connected_components_definition": CONNECTED_COMPONENTS_DEFINITION,
        "duplicate_rule": DUPLICATE_RULE,
        "mmseqs2": {"version": "18.8cc5c", "build": "hd6d6fdc_0", "channel": "bioconda"},
        "mmseqs_command": (
            "mmseqs createdb in.fasta DB ; "
            "mmseqs search DB DB alnDB tmp --prefilter-mode 2 --max-seqs 100000 -e 1e100 "
            "--min-seq-id 0 -c 0 -a --alignment-mode 3 ; "
            "mmseqs convertalis DB DB alnDB pairs.raw.tsv --format-output \"%s\"" % FORMAT_OUTPUT
        ),
        "environment": ENVIRONMENT,
        "input_fasta_sha256": input_fasta_sha256,
        "input_fasta_construction": (
            "deterministic ungapped FASTA from the frozen source: id order sorted, "
            "60-col wrap, LF, trailing newline"
        ),
        "pairs_raw_tsv_canonical_sha256": pairs_raw_canonical_sha256,
        "pairs_raw_tsv_row_count": pairs_raw_row_count,
        "split_groups_csv_sha256": _sha256_file(output_directory / "split_groups.csv"),
        "pairs_filtered_tsv_sha256": _sha256_file(output_directory / "pairs.filtered.tsv"),
        "mmseqs_search_log_sha256": _sha256_file(output_directory / "mmseqs_search.log"),
        "canonical_library": "src/pago_pipeline/apaz_split_groups.py",
        "canonical_library_sha256": _sha256_file(
            PROJECT_ROOT / "src/pago_pipeline/apaz_split_groups.py"
        ),
        "derivation_script": "scripts/derive_apaz_split_groups.py",
        "derivation_script_sha256": _sha256_file(Path(__file__).resolve()),
    }
    (output_directory / "manifest.json").write_bytes(
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resource-directory", type=Path, default=DEFAULT_RESOURCE_DIRECTORY)
    parser.add_argument("--mmseqs", type=str, default="mmseqs")
    arguments = parser.parse_args(argv)

    positives = _read_positive_ungapped(arguments.resource_directory)
    ungapped_by_dataset = {
        "apaz": positives,
        "hisg": _read_seed_ungapped(
            arguments.resource_directory / "source_data" / NEGATIVE_SEED_FILE["hisg"]
        ),
        "eiib": _read_seed_ungapped(
            arguments.resource_directory / "source_data" / NEGATIVE_SEED_FILE["eiib"]
        ),
    }
    for dataset, prefix, expected_count in DATASETS:
        manifest = derive_dataset(
            dataset=dataset,
            prefix=prefix,
            expected_count=expected_count,
            ungapped_by_id=ungapped_by_dataset[dataset],
            resource_directory=arguments.resource_directory,
            mmseqs=arguments.mmseqs,
        )
        print(
            "%-5s %d sequences | %d split groups | %d edges | max %d | %d multi-member"
            % (
                dataset, manifest["n_sequences"], manifest["n_split_groups"],
                manifest["edges_90_80"], manifest["max_group_size"],
                manifest["multi_member_groups"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
