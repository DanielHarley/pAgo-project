from __future__ import annotations

import argparse
import sys
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pago_pipeline.storage import read_json_file, sha256_of_file  # noqa: E402

# Per-layer validators are imported lazily inside each verifier so that a single
# approved scope (for example ``--scope pfam``) runs without importing the
# reference layers that have not been curated or scientifically reviewed yet.


@dataclass(frozen=True)
class VerificationSummary:
    reference_layer: str
    details: tuple[str, ...]


Verifier = Callable[[], VerificationSummary]


def _declared_file_hashes(lock_payload: object) -> dict[str, str]:
    declared_hashes: dict[str, str] = {}

    def visit(value: object) -> None:
        if isinstance(value, Mapping):
            path_value = value.get("path")
            sha256_value = value.get("sha256")
            if isinstance(path_value, str) and isinstance(sha256_value, str):
                previous_sha256 = declared_hashes.setdefault(path_value, sha256_value)
                if previous_sha256 != sha256_value:
                    raise RuntimeError(
                        f"A reference lock declares conflicting hashes for {path_value!r}."
                    )
            for nested_value in value.values():
                visit(nested_value)
        elif isinstance(value, list):
            for nested_value in value:
                visit(nested_value)

    visit(lock_payload)
    return declared_hashes


def verify_every_declared_lock_file(*, lock_file_path: Path) -> int:
    lock_payload = read_json_file(input_file_path=lock_file_path)
    declared_hashes = _declared_file_hashes(lock_payload)
    if not declared_hashes:
        raise RuntimeError(
            f"Reference lock declares no path and SHA-256 entries: {lock_file_path}."
        )
    physical_relative_paths = {
        file_path.relative_to(lock_file_path.parent).as_posix()
        for file_path in lock_file_path.parent.rglob("*")
        if file_path.is_file() and file_path != lock_file_path
    }
    declared_relative_paths = set(declared_hashes)
    unexpected_relative_paths = sorted(
        physical_relative_paths.difference(declared_relative_paths)
    )
    if unexpected_relative_paths:
        raise RuntimeError(
            "Reference directory contains files omitted from its lock: "
            + ", ".join(unexpected_relative_paths)
        )
    for relative_path, expected_sha256 in sorted(declared_hashes.items()):
        pure_path = PurePosixPath(relative_path)
        if (
            pure_path.is_absolute()
            or ".." in pure_path.parts
            or pure_path.as_posix() != relative_path
        ):
            raise RuntimeError(
                f"Reference lock path is not canonical and relative: {relative_path!r}."
            )
        file_path = lock_file_path.parent.joinpath(*pure_path.parts)
        if not file_path.is_file():
            raise FileNotFoundError(f"Locked reference file was not found: {file_path}.")
        actual_sha256 = sha256_of_file(input_file_path=file_path)
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Locked reference file hash mismatch for {relative_path!r}. "
                f"Expected {expected_sha256}, got {actual_sha256}."
            )
    return len(declared_hashes)


def verify_pfam_reference_data() -> VerificationSummary:
    from src.pago_pipeline.pfam_hmm_bundle import validate_pfam_hmm_reference_data

    result = validate_pfam_hmm_reference_data()
    return VerificationSummary(
        reference_layer="pfam_hmm",
        details=(
            f"source={result.registry.source_database} {result.registry.source_release}",
            f"models={len(result.validated_models)}",
            f"registry_sha256={result.registry.registry_sha256}",
            f"pyhmmer={result.pyhmmer_version}",
        ),
    )


def verify_apaz_reference_data() -> VerificationSummary:
    # B2 verification depends only on the canonical split-group library and the
    # committed apaz_seed resources; it must not import any B3 module.
    import csv as _csv
    from collections import Counter as _Counter

    from src.pago_pipeline.apaz_split_groups import (
        validate_apaz_partition_invariants,
        validate_apaz_seed_consistency,
    )

    apaz_seed_directory = (
        PROJECT_ROOT / "src" / "pago_pipeline" / "resources" / "apaz_seed"
    )
    partitions_csv_path = apaz_seed_directory / "apaz_partitions.csv"
    locked_file_count = verify_every_declared_lock_file(
        lock_file_path=apaz_seed_directory / "seeds_lock.json"
    )
    partition_invariants = validate_apaz_partition_invariants(
        partitions_csv_path=partitions_csv_path,
        split_groups_directory=apaz_seed_directory / "split_groups",
    )
    seed_checks = validate_apaz_seed_consistency(resource_directory=apaz_seed_directory)

    with partitions_csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(_csv.DictReader(handle))
    partition_counts = _Counter(row["partition"] for row in rows)
    formatted_partitions = ",".join(
        f"{partition}={partition_counts[partition]}"
        for partition in sorted(partition_counts)
    )
    lock_sha256 = sha256_of_file(input_file_path=apaz_seed_directory / "seeds_lock.json")
    partitions_sha256 = sha256_of_file(input_file_path=partitions_csv_path)
    return VerificationSummary(
        reference_layer="apaz_seed",
        details=(
            f"locked_files={locked_file_count}",
            f"partition_rows={len(rows)}",
            f"split_groups={len({row['split_group_id'] for row in rows})}",
            f"partitions={formatted_partitions}",
            f"partition_invariants_passed={len(partition_invariants)}",
            f"seed_consistency_checks_passed={len(seed_checks)}",
            f"seed_lock_sha256={lock_sha256}",
            f"partitions_sha256={partitions_sha256}",
        ),
    )


def verify_clade_reference_data() -> VerificationSummary:
    from src.pago_pipeline.clade_hmm_build import validate_clade_hmm_build_inputs

    result = validate_clade_hmm_build_inputs()
    locked_file_count = verify_every_declared_lock_file(
        lock_file_path=result.seed_lock_file_path
    )
    partition_counts = Counter(result.partition_by_accession.values())
    formatted_partitions = ",".join(
        f"{partition}={partition_counts[partition]}"
        for partition in sorted(partition_counts)
    )
    return VerificationSummary(
        reference_layer="clade_seed",
        details=(
            f"seed_alignments={len(result.seed_artifacts)}",
            f"locked_files={locked_file_count}",
            f"partition_rows={len(result.partition_by_accession)}",
            f"homology_clusters={len(set(result.homology_cluster_by_accession.values()))}",
            f"partitions={formatted_partitions}",
            f"seed_lock_sha256={result.seed_lock_sha256}",
            f"partitions_sha256={result.partitions_sha256}",
        ),
    )


def verify_midpiwi_reference_data() -> VerificationSummary:
    from src.pago_pipeline.midpiwi_reference_tree import (
        validate_midpiwi_reference_tree,
    )

    result = validate_midpiwi_reference_tree()
    return VerificationSummary(
        reference_layer="midpiwi_reference",
        details=(
            f"sequences={result.reference_sequence_count}",
            f"alignment_columns={result.reference_alignment_length}",
            f"edges={result.reference_edge_count}",
            f"groups={','.join(result.reference_groups)}",
            f"reference_profile_sha256={result.reference_profile_sha256}",
            f"tree_lock_sha256={result.tree_lock_sha256}",
        ),
    )


VERIFIERS: dict[str, Verifier] = {
    "pfam": verify_pfam_reference_data,
    "apaz": verify_apaz_reference_data,
    "clade": verify_clade_reference_data,
    "tree": verify_midpiwi_reference_data,
}


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate committed Phase B reference locks, SHA-256 hashes, "
            "inventories, partitions, and phylogenetic-reference invariants "
            "without modifying any artifact."
        )
    )
    parser.add_argument(
        "--scope",
        action="append",
        choices=tuple(VERIFIERS),
        help=(
            "Reference layer to validate. Repeat for multiple layers. "
            "Omission validates every layer."
        ),
    )
    return parser


def verify_selected_reference_layers(
    *,
    scopes: Sequence[str] | None = None,
) -> tuple[VerificationSummary, ...]:
    selected_scopes = tuple(scopes) if scopes else tuple(VERIFIERS)
    if len(set(selected_scopes)) != len(selected_scopes):
        raise ValueError("Each verification scope may be selected only once.")
    return tuple(VERIFIERS[scope]() for scope in selected_scopes)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    summaries = verify_selected_reference_layers(scopes=arguments.scope)
    for summary in summaries:
        print(f"{summary.reference_layer}: OK")
        for detail in summary.details:
            print(f"  {detail}")
    print(f"Verified reference layers: {len(summaries)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
