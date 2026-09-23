"""Canonical APAZ split-group library — the single source of truth for the
indivisible unit of the Phase B curation partitions.

Two responsibilities, both pure and deterministic:

1. ``build_split_groups`` — turn MMseqs2 all-vs-all pairwise alignments into
   connected components under an explicit edge rule, plus a hard exact-duplicate
   invariant. MMseqs2 is never called here; it only produces the candidate
   alignments consumed as ``AlignedPair`` rows.
2. ``partition_split_groups`` — assign whole split groups to
   FINAL_HOLDOUT / CALIBRATION / BUILD, stratified by APAZ subgroup, hitting
   exact per-stratum sequence-count targets via a deterministic subset-sum
   search. It aborts rather than silently adjusting a target it cannot meet.

Edge rule (identical to the B2 sensitivity analysis semantics)::

    comparable_columns = columns where BOTH aligned sequences carry a residue
    pair_identity      = identical_residue_columns / comparable_columns
    pair_coverage      = comparable_columns / min(query_len, target_len)   # smaller member
    edge  iff  pair_identity >= identity_threshold AND pair_coverage >= coverage_threshold

``query_len`` / ``target_len`` are the full ungapped residue counts.

Hard invariant, independent of any alignment::

    one sequence_sha256  ->  exactly one split_group_id

MMseqs2 ``pident`` / ``alnlen`` are never used to decide an edge; identity and
coverage are recomputed from the aligned strings.

A 90 %/80 % split group is a HIGH-SIMILARITY / REDUNDANCY control unit. It is
**not** a claim that homology is absent below 90 % identity.
"""
from __future__ import annotations

import csv
import hashlib
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

GAP_CHARACTERS = frozenset("-.")

DEFAULT_IDENTITY_THRESHOLD = 0.90
DEFAULT_COVERAGE_THRESHOLD = 0.80

SPLIT_GROUP_METHOD = "mmseqs2_all_vs_all_then_explicit_90_80_connected_components"
IDENTITY_DEFINITION = "identical_residue_columns / comparable_residue_columns"
COVERAGE_DEFINITION = "comparable_residue_columns / min(query_ungapped_len, target_ungapped_len)"
CONNECTED_COMPONENTS_DEFINITION = (
    "single-linkage union-find over 90/80 edges plus exact sequence_sha256 edges; "
    "component id = <prefix>_<sha256 of the newline-joined sorted member ids>[:16]"
)
DUPLICATE_RULE = "sequences sharing sequence_sha256 are unioned regardless of alignment"


# --------------------------------------------------------------------------- #
# 1. split groups                                                             #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class AlignedPair:
    query_id: str
    target_id: str
    query_aligned: str
    target_aligned: str
    query_length: int
    target_length: int


@dataclass(frozen=True)
class PairScore:
    query_id: str
    target_id: str
    comparable_columns: int
    identical_columns: int
    pair_identity: float
    pair_coverage: float
    is_edge: bool


@dataclass(frozen=True)
class SplitGroupResult:
    split_group_id_by_sequence_id: Mapping[str, str]
    members_by_split_group_id: Mapping[str, tuple[str, ...]]
    edge_pairs: tuple[tuple[str, str], ...]
    identity_threshold: float
    coverage_threshold: float


def _require_equal_alignment_widths(query_aligned: str, target_aligned: str) -> None:
    if len(query_aligned) != len(target_aligned):
        raise ValueError(
            "Aligned query and target strings must have equal length; got "
            f"{len(query_aligned)} and {len(target_aligned)}."
        )


def score_aligned_pair(
    *,
    query_aligned: str,
    target_aligned: str,
    query_length: int,
    target_length: int,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
) -> tuple[int, int, float, float, bool]:
    """Return (comparable, identical, identity, coverage, is_edge) for one pair."""
    _require_equal_alignment_widths(query_aligned, target_aligned)
    if query_length <= 0 or target_length <= 0:
        raise ValueError("query_length and target_length must be positive.")
    query_upper = query_aligned.upper()
    target_upper = target_aligned.upper()
    comparable = 0
    identical = 0
    for query_residue, target_residue in zip(query_upper, target_upper):
        if query_residue in GAP_CHARACTERS or target_residue in GAP_CHARACTERS:
            continue
        comparable += 1
        if query_residue == target_residue:
            identical += 1
    if comparable == 0:
        return 0, 0, 0.0, 0.0, False
    identity = identical / comparable
    coverage = comparable / min(query_length, target_length)
    is_edge = identity >= identity_threshold and coverage >= coverage_threshold
    return comparable, identical, identity, coverage, is_edge


def score_pairs(
    pairs: Iterable[AlignedPair],
    *,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
) -> list[PairScore]:
    scored: list[PairScore] = []
    for pair in pairs:
        if pair.query_id == pair.target_id:
            continue
        comparable, identical, identity, coverage, is_edge = score_aligned_pair(
            query_aligned=pair.query_aligned,
            target_aligned=pair.target_aligned,
            query_length=pair.query_length,
            target_length=pair.target_length,
            identity_threshold=identity_threshold,
            coverage_threshold=coverage_threshold,
        )
        scored.append(
            PairScore(
                query_id=pair.query_id,
                target_id=pair.target_id,
                comparable_columns=comparable,
                identical_columns=identical,
                pair_identity=identity,
                pair_coverage=coverage,
                is_edge=is_edge,
            )
        )
    return scored


class _DisjointSet:
    def __init__(self, items: Sequence[str]) -> None:
        self._parent = {item: item for item in items}

    def contains(self, item: str) -> bool:
        return item in self._parent

    def find(self, item: str) -> str:
        root = item
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[item] != root:
            self._parent[item], item = root, self._parent[item]
        return root

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        low, high = sorted((left_root, right_root))
        self._parent[high] = low


def build_split_groups(
    *,
    sequence_ids: Sequence[str],
    sequence_sha256_by_id: Mapping[str, str],
    edge_pairs: Iterable[tuple[str, str]],
    split_group_id_prefix: str,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD,
) -> SplitGroupResult:
    """Connected components over (alignment edges) + (exact sequence_sha256 edges).

    The result is invariant to the order of ``sequence_ids`` and ``edge_pairs``.
    """
    material = list(sequence_ids)
    unique_ids = sorted(set(material))
    if len(unique_ids) != len(material):
        raise ValueError("sequence_ids must be unique.")
    missing_hashes = [i for i in unique_ids if i not in sequence_sha256_by_id]
    if missing_hashes:
        raise ValueError(f"sequence_sha256 missing for: {missing_hashes}.")

    disjoint = _DisjointSet(unique_ids)

    ids_by_hash: dict[str, list[str]] = {}
    for sequence_id in unique_ids:
        ids_by_hash.setdefault(sequence_sha256_by_id[sequence_id], []).append(sequence_id)
    for members in ids_by_hash.values():
        for other in members[1:]:
            disjoint.union(members[0], other)

    normalized_edges: list[tuple[str, str]] = []
    for left, right in edge_pairs:
        if not disjoint.contains(left) or not disjoint.contains(right):
            raise ValueError(f"edge references an unknown sequence id: {(left, right)}.")
        if left == right:
            continue
        disjoint.union(left, right)
        normalized_edges.append(tuple(sorted((left, right))))

    components: dict[str, list[str]] = {}
    for sequence_id in unique_ids:
        components.setdefault(disjoint.find(sequence_id), []).append(sequence_id)

    split_group_id_by_sequence_id: dict[str, str] = {}
    members_by_split_group_id: dict[str, tuple[str, ...]] = {}
    for members in components.values():
        ordered_members = tuple(sorted(members))
        digest = hashlib.sha256("\n".join(ordered_members).encode("utf-8")).hexdigest()[:16]
        split_group_id = f"{split_group_id_prefix}_{digest}"
        if split_group_id in members_by_split_group_id:
            raise RuntimeError(f"split_group_id collision for {split_group_id!r}.")
        members_by_split_group_id[split_group_id] = ordered_members
        for member in ordered_members:
            split_group_id_by_sequence_id[member] = split_group_id

    _assert_hash_maps_to_single_group(
        sequence_sha256_by_id=sequence_sha256_by_id,
        split_group_id_by_sequence_id=split_group_id_by_sequence_id,
    )

    return SplitGroupResult(
        split_group_id_by_sequence_id=dict(sorted(split_group_id_by_sequence_id.items())),
        members_by_split_group_id=dict(sorted(members_by_split_group_id.items())),
        edge_pairs=tuple(sorted(set(normalized_edges))),
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )


def _assert_hash_maps_to_single_group(
    *,
    sequence_sha256_by_id: Mapping[str, str],
    split_group_id_by_sequence_id: Mapping[str, str],
) -> None:
    groups_by_hash: dict[str, set[str]] = {}
    for sequence_id, split_group_id in split_group_id_by_sequence_id.items():
        groups_by_hash.setdefault(sequence_sha256_by_id[sequence_id], set()).add(split_group_id)
    violations = {h: sorted(g) for h, g in groups_by_hash.items() if len(g) != 1}
    if violations:
        raise RuntimeError(
            "Invariant violated: a sequence_sha256 maps to more than one "
            f"split_group_id: {violations}."
        )


# --------------------------------------------------------------------------- #
# 2. deterministic stratified partition                                       #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class StratumPartitionPlan:
    stratum: str
    fill_order: tuple[tuple[str, int], ...]
    selection_order_by_partition: Mapping[str, tuple[str, ...]]
    groups_by_partition: Mapping[str, tuple[str, ...]]


@dataclass(frozen=True)
class StratifiedPartitionResult:
    partition_by_sequence_id: Mapping[str, str]
    partition_by_split_group_id: Mapping[str, str]
    stratum_by_split_group_id: Mapping[str, str]
    partition_salt: str
    plans_by_stratum: Mapping[str, StratumPartitionPlan]


def selection_key_sha256(
    *, partition_salt: str, stratum: str, split_group_id: str, partition: str
) -> str:
    payload = "|".join((partition_salt, stratum, split_group_id, partition))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _subset_sum_first(
    items: Sequence[tuple[str, int]], target: int
) -> list[str] | None:
    """First subset (in the given deterministic item order, preferring to
    INCLUDE earlier items) whose sizes sum exactly to ``target``. None if no
    subset reaches the target. ``items`` is [(split_group_id, size), ...]."""
    if target < 0:
        return None
    n = len(items)
    suffix = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suffix[i] = suffix[i + 1] + items[i][1]
    result: list[str] | None = None
    chosen: list[str] = []

    def dfs(index: int, current: int) -> bool:
        nonlocal result
        if current == target:
            result = list(chosen)
            return True
        if index == n or current > target:
            return False
        if current + suffix[index] < target:
            return False
        group_id, size = items[index]
        chosen.append(group_id)
        if dfs(index + 1, current + size):
            return True
        chosen.pop()
        return dfs(index + 1, current)

    dfs(0, 0)
    return result


def partition_split_groups(
    *,
    members_by_split_group_id: Mapping[str, Sequence[str]],
    stratum_by_sequence_id: Mapping[str, str],
    fill_order_by_stratum: Mapping[str, Sequence[tuple[str, int]]],
    partition_salt: str,
) -> StratifiedPartitionResult:
    """Assign whole split groups to partitions, stratified, hitting exact counts.

    ``fill_order_by_stratum[stratum]`` is an ordered list of
    ``(partition_name, target_sequence_count)``. Every partition except the last
    is filled by a deterministic subset-sum over the not-yet-assigned groups of
    that stratum, ordered by ``selection_key_sha256``. The last partition
    receives the remainder and its size is asserted against its declared target.
    A subset-sum that cannot be met raises RuntimeError (no silent adjustment).
    """
    if not partition_salt or partition_salt != partition_salt.strip():
        raise ValueError("partition_salt must be a non-empty trimmed string.")

    size_by_group = {g: len(tuple(m)) for g, m in members_by_split_group_id.items()}
    stratum_by_group: dict[str, str] = {}
    for group_id, member_ids in members_by_split_group_id.items():
        strata = {stratum_by_sequence_id[m] for m in member_ids}
        if len(strata) != 1:
            raise RuntimeError(
                f"split_group {group_id!r} crosses strata {sorted(strata)}; "
                "a split group must be contained in a single stratum."
            )
        stratum_by_group[group_id] = next(iter(strata))

    groups_by_stratum: dict[str, list[str]] = {}
    for group_id, stratum in stratum_by_group.items():
        groups_by_stratum.setdefault(stratum, []).append(group_id)

    declared_strata = set(fill_order_by_stratum)
    observed_strata = set(groups_by_stratum)
    if declared_strata != observed_strata:
        raise RuntimeError(
            "fill_order strata do not match the split-group strata. "
            f"Declared={sorted(declared_strata)}; observed={sorted(observed_strata)}."
        )

    partition_by_group: dict[str, str] = {}
    plans: dict[str, StratumPartitionPlan] = {}

    for stratum in sorted(fill_order_by_stratum):
        fill_order = tuple((str(p), int(t)) for p, t in fill_order_by_stratum[stratum])
        if len(fill_order) < 2:
            raise ValueError("Each stratum needs at least two partitions in fill_order.")
        stratum_total = sum(size_by_group[g] for g in groups_by_stratum[stratum])
        if sum(t for _, t in fill_order) != stratum_total:
            raise RuntimeError(
                f"stratum {stratum!r}: targets {[t for _, t in fill_order]} sum to "
                f"{sum(t for _, t in fill_order)} but the stratum holds {stratum_total} "
                "sequences."
            )

        remaining = set(groups_by_stratum[stratum])
        selection_order_by_partition: dict[str, tuple[str, ...]] = {}
        groups_by_partition: dict[str, tuple[str, ...]] = {}

        for partition, target in fill_order[:-1]:
            ordered = sorted(
                remaining,
                key=lambda g: selection_key_sha256(
                    partition_salt=partition_salt,
                    stratum=stratum,
                    split_group_id=g,
                    partition=partition,
                ),
            )
            selection_order_by_partition[partition] = tuple(ordered)
            chosen = _subset_sum_first(
                [(g, size_by_group[g]) for g in ordered], target
            )
            if chosen is None:
                raise RuntimeError(
                    f"stratum {stratum!r}: no combination of whole split groups sums "
                    f"to the {partition} target of {target} sequences. Aborting; "
                    "targets are not adjusted silently."
                )
            for group_id in chosen:
                partition_by_group[group_id] = partition
            groups_by_partition[partition] = tuple(sorted(chosen))
            remaining.difference_update(chosen)

        last_partition, last_target = fill_order[-1]
        remainder_size = sum(size_by_group[g] for g in remaining)
        if remainder_size != last_target:
            raise RuntimeError(
                f"stratum {stratum!r}: remainder for {last_partition} is "
                f"{remainder_size} sequences, expected {last_target}."
            )
        for group_id in remaining:
            partition_by_group[group_id] = last_partition
        groups_by_partition[last_partition] = tuple(sorted(remaining))
        selection_order_by_partition[last_partition] = tuple(sorted(remaining))

        plans[stratum] = StratumPartitionPlan(
            stratum=stratum,
            fill_order=fill_order,
            selection_order_by_partition=dict(selection_order_by_partition),
            groups_by_partition=dict(groups_by_partition),
        )

    partition_by_sequence_id: dict[str, str] = {}
    for group_id, member_ids in members_by_split_group_id.items():
        for member in member_ids:
            partition_by_sequence_id[member] = partition_by_group[group_id]

    return StratifiedPartitionResult(
        partition_by_sequence_id=dict(sorted(partition_by_sequence_id.items())),
        partition_by_split_group_id=dict(sorted(partition_by_group.items())),
        stratum_by_split_group_id=dict(sorted(stratum_by_group.items())),
        partition_salt=partition_salt,
        plans_by_stratum=dict(sorted(plans.items())),
    )


# --------------------------------------------------------------------------- #
# 3. offline invariant validation                                            #
# --------------------------------------------------------------------------- #
APAZ_SUBGROUPS = ("Ia", "Ib", "IIa", "IIb", "III")
NEGATIVE_FAMILY_SOURCE_SUFFIX = {"PF01634": "hisg", "PF00367": "eiib"}
APAZ_POSITIVE_COUNTS = {"BUILD": 337, "CALIBRATION": 72, "FINAL_HOLDOUT": 72}
NEGATIVE_FAMILY_COUNTS = {
    "PF01634": {"CALIBRATION": 255, "FINAL_HOLDOUT": 254, "BUILD": 0},
    "PF00367": {"CALIBRATION": 259, "FINAL_HOLDOUT": 258, "BUILD": 0},
}


def _read_split_group_csv(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    sha_by_id: dict[str, str] = {}
    group_by_id: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            sha_by_id[row["sequence_id"]] = row["sequence_sha256"]
            group_by_id[row["sequence_id"]] = row["split_group_id"]
    return sha_by_id, group_by_id


def validate_apaz_partition_invariants(
    *,
    partitions_csv_path: Path,
    split_groups_directory: Path,
) -> dict[str, bool]:
    """Prove every B2 partition invariant offline. Raises on any violation."""
    rows: list[dict[str, str]] = []
    with Path(partitions_csv_path).open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    apaz_rows = [r for r in rows if r["class_label"] == "APAZ"]
    neg_rows = [r for r in rows if r["class_label"] == "HARD_NEGATIVE"]

    checks: dict[str, bool] = {}

    accessions = [r["accession"] for r in rows]
    checks["no_duplicate_accession"] = len(accessions) == len(set(accessions))

    checks["homology_cluster_id_equals_split_group_id"] = all(
        r["homology_cluster_id"] == r["split_group_id"] for r in rows
    )

    part_by_sha: dict[str, set[str]] = {}
    part_by_group: dict[str, set[str]] = {}
    for r in rows:
        part_by_sha.setdefault(r["sequence_sha256"], set()).add(r["partition"])
        part_by_group.setdefault(r["split_group_id"], set()).add(r["partition"])
    checks["sequence_sha256_one_partition"] = all(len(v) == 1 for v in part_by_sha.values())
    checks["split_group_id_one_partition"] = all(len(v) == 1 for v in part_by_group.values())

    subgroups_by_group: dict[str, set[str]] = {}
    for r in apaz_rows:
        subgroups_by_group.setdefault(r["split_group_id"], set()).add(r["subgroup"])
    checks["no_split_group_crosses_subgroup"] = all(
        len(v) == 1 for v in subgroups_by_group.values()
    )

    build_pos = {r["accession"] for r in apaz_rows if r["partition"] == "BUILD"}
    nonbuild_pos = {r["accession"] for r in apaz_rows if r["partition"] != "BUILD"}
    checks["no_nonbuild_positive_in_build"] = build_pos.isdisjoint(nonbuild_pos)
    checks["no_hard_negative_in_build"] = all(r["partition"] != "BUILD" for r in neg_rows)

    pos_sha = {r["sequence_sha256"] for r in apaz_rows}
    neg_sha = {r["sequence_sha256"] for r in neg_rows}
    checks["positives_negatives_disjoint_sequence_sha256"] = pos_sha.isdisjoint(neg_sha)

    pos_counts: dict[str, int] = {}
    for r in apaz_rows:
        pos_counts[r["partition"]] = pos_counts.get(r["partition"], 0) + 1
    checks["positive_counts_match_targets"] = pos_counts == APAZ_POSITIVE_COUNTS

    for family, expected in NEGATIVE_FAMILY_COUNTS.items():
        fam_rows = [r for r in neg_rows if r["source"].endswith("_%s_seed" % family)]
        got: dict[str, int] = {"BUILD": 0, "CALIBRATION": 0, "FINAL_HOLDOUT": 0}
        for r in fam_rows:
            got[r["partition"]] += 1
        checks["negative_counts_%s" % family] = got == expected

    # cross-check split_group_id and sequence_sha256 against the committed
    # split_groups/ resources.
    sg_dir = Path(split_groups_directory)
    apaz_sha, apaz_group = _read_split_group_csv(sg_dir / "apaz" / "split_groups.csv")
    resource_group_by_row_accession_ok = True
    resource_sha_ok = True
    for r in apaz_rows:
        acc = r["accession"]
        if apaz_group.get(acc) != r["split_group_id"]:
            resource_group_by_row_accession_ok = False
        if apaz_sha.get(acc) != r["sequence_sha256"]:
            resource_sha_ok = False
    checks["apaz_rows_match_split_groups_resource"] = (
        resource_group_by_row_accession_ok and resource_sha_ok
    )

    for family, dataset in NEGATIVE_FAMILY_SOURCE_SUFFIX.items():
        n_sha, n_group = _read_split_group_csv(sg_dir / dataset / "split_groups.csv")
        fam_rows = [r for r in neg_rows if r["source"].endswith("_%s_seed" % family)]
        ok = True
        for r in fam_rows:
            src = r["source_sequence_id"]
            if n_group.get(src) != r["split_group_id"] or n_sha.get(src) != r["sequence_sha256"]:
                ok = False
        checks["%s_rows_match_split_groups_resource" % dataset] = ok

    failed = sorted(name for name, ok in checks.items() if not ok)
    if failed:
        raise RuntimeError("APAZ partition invariant(s) failed: " + ", ".join(failed))
    return checks


APAZ_SUBGROUP_SEED_FILE = {
    "Ia": "apaz_Ia.sto", "Ib": "apaz_Ib.sto", "IIa": "apaz_IIa.sto",
    "IIb": "apaz_IIb.sto", "III": "apaz_III.sto",
}
_STOCKHOLM_COORDINATE_SUFFIX = re.compile(r"/\d+-\d+$")


def _read_stockholm_alignment(path: Path) -> dict[str, str]:
    fragments: dict[str, list[str]] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line == "//" or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 2:
            raise RuntimeError(f"Malformed Stockholm row in {path}.")
        fragments.setdefault(fields[0], []).append(fields[1])
    return {identifier: "".join(parts) for identifier, parts in fragments.items()}


def _ungapped(sequence: str) -> str:
    return sequence.replace("-", "").replace(".", "").upper()


def validate_apaz_seed_consistency(
    *,
    resource_directory: Path,
) -> dict[str, bool]:
    """Validate the committed APAZ seed alignments against apaz_partitions.csv.

    Pure data checks on the B2 resources only (no HMM construction, no B3
    module). Raises on any violation.
    """
    resource_directory = Path(resource_directory)
    with (resource_directory / "apaz_partitions.csv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rows = [r for r in csv.DictReader(handle) if r["class_label"] == "APAZ"]
    partition_by_accession = {r["accession"]: r["partition"] for r in rows}
    subgroup_by_accession = {r["accession"]: r["subgroup"] for r in rows}
    sha256_by_accession = {r["accession"]: r["sequence_sha256"] for r in rows}

    checks: dict[str, bool] = {}

    subgroup_alignments: dict[str, dict[str, str]] = {}
    for subgroup, file_name in APAZ_SUBGROUP_SEED_FILE.items():
        aligned_by_raw_id = _read_stockholm_alignment(resource_directory / file_name)
        aligned_by_accession = {
            _STOCKHOLM_COORDINATE_SUFFIX.sub("", raw_id): sequence
            for raw_id, sequence in aligned_by_raw_id.items()
        }
        subgroup_alignments[subgroup] = aligned_by_accession

    global_alignment = {
        _STOCKHOLM_COORDINATE_SUFFIX.sub("", raw_id): sequence
        for raw_id, sequence in _read_stockholm_alignment(
            resource_directory / "apaz_global.sto"
        ).items()
    }

    subgroup_accession_sets = {sg: set(a) for sg, a in subgroup_alignments.items()}
    all_subgroup_accessions: set[str] = set()
    overlap = False
    for accessions in subgroup_accession_sets.values():
        if all_subgroup_accessions & accessions:
            overlap = True
        all_subgroup_accessions |= accessions
    checks["subgroup_seeds_are_accession_disjoint"] = not overlap

    seed_partition_ok = True
    seed_subgroup_ok = True
    for subgroup, accessions in subgroup_accession_sets.items():
        for accession in accessions:
            if partition_by_accession.get(accession) != "BUILD":
                seed_partition_ok = False
            if subgroup_by_accession.get(accession) != subgroup:
                seed_subgroup_ok = False
    checks["subgroup_seed_members_are_all_build"] = seed_partition_ok
    checks["subgroup_seed_members_match_declared_subgroup"] = seed_subgroup_ok

    build_accessions = {a for a, p in partition_by_accession.items() if p == "BUILD"}
    checks["subgroup_seed_union_equals_build_set"] = (
        all_subgroup_accessions == build_accessions
    )
    checks["global_seed_accessions_equal_subgroup_union"] = (
        set(global_alignment) == all_subgroup_accessions
    )

    ungapped_consistent = True
    sha256_consistent = True
    for subgroup, aligned_by_accession in subgroup_alignments.items():
        for accession, aligned in aligned_by_accession.items():
            subgroup_ungapped = _ungapped(aligned)
            global_ungapped = _ungapped(global_alignment.get(accession, ""))
            if subgroup_ungapped != global_ungapped:
                ungapped_consistent = False
            if (
                hashlib.sha256(subgroup_ungapped.encode("ascii")).hexdigest()
                != sha256_by_accession.get(accession)
            ):
                sha256_consistent = False
    checks["seed_ungapped_sequences_match_between_subgroup_and_global"] = ungapped_consistent
    checks["seed_ungapped_sequence_sha256_matches_partitions_csv"] = sha256_consistent

    failed = sorted(name for name, ok in checks.items() if not ok)
    if failed:
        raise RuntimeError("APAZ seed consistency check(s) failed: " + ", ".join(failed))
    return checks
