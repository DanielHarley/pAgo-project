"""MID-PIWI high-similarity split-group labelling and validation (ADR-0013).

Groups the 1,002 MID-PIWI regions extracted by the Ryazansky catalog
(``clade_reference_catalog.py``, ADR-0011) into high-similarity split groups
under the same 90% identity / 80% coverage redundancy-control rule ADR-0009
established for APAZ. A split group is a redundancy-control unit, not a
phylogenetic clade and not a claim of homology (ADR-0013).

Duplication note: the connected-component construction and the identity/
coverage edge rule below (``score_aligned_pair``, ``score_pairs``,
``build_split_groups``, and the constants describing them) are a minimal,
deliberate duplicate of the same logic in the historical
``src/pago_pipeline/apaz_split_groups.py`` (APAZ reference layer, not yet
reintegrated to `master`). This module does not depend on that file — it
depends on nothing outside the standard library. When the APAZ lineage is
reintegrated, whether to extract one shared reducer module is a separate,
later decision; duplicating this ~140-line pure function here now avoids
coupling this PR to unintegrated code.

This module adds, beyond the reducer:

* ``label_split_group`` — classify a finished 90/80 component by the
  ``curated_pago_clade`` of its members (LABEL_CONFLICT, manual-review
  propagation, partition eligibility);
* ``assert_complete_all_vs_all`` — refuse to build split groups unless the
  MMseqs2 all-vs-all covered every unordered pair;
* ``validate_clade_reference_data`` — revalidate the committed Ryazansky
  catalog and MID-PIWI split-group resources without importing any B5/B6
  module.

Nothing here builds a tree, a placement, or an HMM, and nothing here assigns
BUILD/CALIBRATION/FINAL_HOLDOUT partitions.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

PathLike: TypeAlias = str | Path

# --------------------------------------------------------------------------- #
# 90/80 connected-component reducer (duplicated from apaz_split_groups.py —
# see the module docstring).                                                 #
# --------------------------------------------------------------------------- #
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
# MID-PIWI-specific labelling and validation.                                 #
# --------------------------------------------------------------------------- #
RESOLVED_CLADES: tuple[str, ...] = ("LONG_A", "LONG_B", "SHORT")

# Accessions parked in QUARANTINE by the Ryazansky catalog (architecture vs
# Ryazansky phylogeny, ADR-0011). Any split group containing one of these
# needs a human decision before it can be considered for a partition.
QUARANTINE_ANCHOR_ACCESSIONS: Mapping[str, str] = {
    "WP_010878815.1": "AfAgo",
    "WP_012735993.1": "SiAgo",
}

GROUP_STATUS_OK = "OK"
GROUP_STATUS_REVIEW = "REVIEW"
GROUP_STATUS_QUARANTINE = "QUARANTINE"


@dataclass(frozen=True)
class GroupLabel:
    group_label: str
    group_status: str
    requires_manual_review: bool
    partition_eligible: bool
    resolved_clades: tuple[str, ...]
    contains_unresolved: bool
    quarantine_anchors: tuple[str, ...]


def label_split_group(
    *,
    member_ids: Sequence[str],
    curated_pago_clade_by_id: Mapping[str, str],
) -> GroupLabel:
    """Classify one finished 90/80 split group.

    * more than one resolved clade among {LONG_A, LONG_B, SHORT}
      -> ``LABEL_CONFLICT`` / QUARANTINE (never forced);
    * exactly one resolved clade + at least one UNRESOLVED member
      -> ``SINGLE_RESOLVED_PLUS_UNRESOLVED`` / REVIEW (not an automatic conflict);
    * exactly one resolved clade -> ``SINGLE_<clade>`` / OK;
    * only UNRESOLVED members -> ``UNRESOLVED_ONLY`` / REVIEW.

    A member in ``QUARANTINE_ANCHOR_ACCESSIONS`` forces
    ``requires_manual_review`` and status QUARANTINE, and removes partition
    eligibility.
    """
    member_clades = {curated_pago_clade_by_id[m] for m in member_ids}
    resolved = tuple(c for c in RESOLVED_CLADES if c in member_clades)
    contains_unresolved = "UNRESOLVED" in member_clades
    anchors = tuple(sorted(
        QUARANTINE_ANCHOR_ACCESSIONS[m]
        for m in member_ids
        if m in QUARANTINE_ANCHOR_ACCESSIONS
    ))
    requires_manual_review = bool(anchors)

    if len(resolved) > 1:
        group_label, group_status = "LABEL_CONFLICT", GROUP_STATUS_QUARANTINE
    elif len(resolved) == 1 and contains_unresolved:
        group_label = "SINGLE_RESOLVED_PLUS_UNRESOLVED"
        group_status = GROUP_STATUS_REVIEW
    elif len(resolved) == 1:
        group_label = f"SINGLE_{resolved[0]}"
        group_status = GROUP_STATUS_OK
    else:
        group_label, group_status = "UNRESOLVED_ONLY", GROUP_STATUS_REVIEW

    if requires_manual_review:
        group_status = GROUP_STATUS_QUARANTINE

    partition_eligible = (
        group_status == GROUP_STATUS_OK and not requires_manual_review
    )
    return GroupLabel(
        group_label=group_label,
        group_status=group_status,
        requires_manual_review=requires_manual_review,
        partition_eligible=partition_eligible,
        resolved_clades=resolved,
        contains_unresolved=contains_unresolved,
        quarantine_anchors=anchors,
    )


def expected_unordered_pairs(n_sequences: int) -> int:
    return n_sequences * (n_sequences - 1) // 2


def assert_complete_all_vs_all(
    *, observed_unordered_pairs: int, n_sequences: int
) -> None:
    """Raise unless every unordered pair was aligned by the MMseqs2 all-vs-all."""
    expected = expected_unordered_pairs(n_sequences)
    if observed_unordered_pairs != expected:
        raise RuntimeError(
            "MID-PIWI all-vs-all is incomplete: "
            f"{observed_unordered_pairs}/{expected} unordered pairs aligned."
        )


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(Path(path).read_bytes())


def _read_csv(path: Path) -> list[dict]:
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8").splitlines()))


def validate_clade_reference_data(*, clade_seed_directory: PathLike) -> dict:
    """Revalidate the committed Ryazansky catalog and MID-PIWI split groups.

    Raises ``RuntimeError`` on any inconsistency. B5/B6-free.
    """
    directory = Path(clade_seed_directory)
    source_data = directory / "source_data"
    checks: dict[str, bool] = {}

    # 1. source manifest hashes
    manifest = json.loads((source_data / "source_manifest.json").read_text("utf-8"))
    for entry in manifest["materialized_files"]:
        path = source_data / entry["file_name"]
        actual_md5 = hashlib.md5(path.read_bytes()).hexdigest()
        actual_sha = _sha256_file(path)
        if actual_md5 != entry["md5"] or actual_sha != entry["sha256"]:
            raise RuntimeError(
                f"source_manifest hash mismatch for {entry['file_name']}."
            )
    checks["source_manifest_hashes_match"] = True

    # 2. catalog identity + internal consistency
    catalog_path = directory / "ryazansky_s1_catalog.csv"
    summary = json.loads(
        (directory / "ryazansky_s1_catalog_summary.json").read_text("utf-8")
    )
    if _sha256_file(catalog_path) != summary["ryazansky_s1_catalog_sha256"]:
        raise RuntimeError("ryazansky_s1_catalog.csv sha256 != summary.")
    checks["catalog_sha256_matches_summary"] = True

    catalog_rows = _read_csv(catalog_path)
    if len(catalog_rows) != summary["n_total"]:
        raise RuntimeError("catalog row count != summary n_total.")
    clade_by_id = {r["source_accession"]: r["curated_pago_clade"] for r in catalog_rows}
    for row in catalog_rows:
        if row["midpiwi_extraction_status"] == "OK":
            expected = _sha256_bytes(row["midpiwi_sequence"].encode("ascii"))
            if expected != row["midpiwi_region_sha256"]:
                raise RuntimeError(
                    f"catalog midpiwi_region_sha256 mismatch for "
                    f"{row['source_accession']}."
                )
    checks["catalog_midpiwi_region_hashes_match"] = True

    ok_ids = {
        r["source_accession"] for r in catalog_rows
        if r["midpiwi_extraction_status"] == "OK"
    }

    # 3. split-group manifest hashes
    sg_dir = directory / "split_groups" / "midpiwi"
    sg_manifest = json.loads((sg_dir / "manifest.json").read_text("utf-8"))
    for file_name, key in (
        ("split_groups.csv", "split_groups_csv_sha256"),
        ("groups.csv", "groups_csv_sha256"),
        ("pairs.filtered.tsv", "pairs_filtered_tsv_sha256"),
        ("diagnostic_thresholds.json", "diagnostic_thresholds_json_sha256"),
        ("mmseqs_search.log", "mmseqs_search_log_sha256"),
    ):
        if _sha256_file(sg_dir / file_name) != sg_manifest[key]:
            raise RuntimeError(f"split-group manifest hash mismatch for {file_name}.")
    if _sha256_file(catalog_path) != sg_manifest["source_catalog_sha256"]:
        raise RuntimeError("split-group manifest source_catalog_sha256 stale.")
    checks["split_group_manifest_hashes_match"] = True

    # 4. split_groups.csv covers exactly the OK regions
    sg_rows = _read_csv(sg_dir / "split_groups.csv")
    sg_ids = {r["sequence_id"] for r in sg_rows}
    if sg_ids != ok_ids:
        raise RuntimeError("split_groups.csv sequence set != catalog OK regions.")
    if any(
        r["curated_pago_clade"] != clade_by_id[r["sequence_id"]] for r in sg_rows
    ):
        raise RuntimeError("split_groups.csv clade disagrees with the catalog.")
    checks["split_groups_cover_ok_regions"] = True

    # 5. re-label every group from split_groups.csv + catalog -> must match groups.csv
    members_by_group: dict[str, list[str]] = {}
    for row in sg_rows:
        members_by_group.setdefault(row["split_group_id"], []).append(row["sequence_id"])
    group_rows = {r["split_group_id"]: r for r in _read_csv(sg_dir / "groups.csv")}
    if set(group_rows) != set(members_by_group):
        raise RuntimeError("groups.csv group set != split_groups.csv.")
    label_conflicts = 0
    partition_eligible = 0
    unresolved_groups = 0
    quarantine_groups = 0
    clade_group_counts: Counter = Counter()
    for group_id, members in members_by_group.items():
        label = label_split_group(
            member_ids=members, curated_pago_clade_by_id=clade_by_id
        )
        recorded = group_rows[group_id]
        if (
            recorded["group_label"] != label.group_label
            or recorded["group_status"] != label.group_status
            or (recorded["requires_manual_review"] == "True")
            != label.requires_manual_review
            or (recorded["partition_eligible"] == "True") != label.partition_eligible
        ):
            raise RuntimeError(
                f"groups.csv label disagrees with recomputation for {group_id}."
            )
        if int(recorded["n_members"]) != len(members):
            raise RuntimeError(f"groups.csv n_members wrong for {group_id}.")
        label_conflicts += label.group_label == "LABEL_CONFLICT"
        partition_eligible += label.partition_eligible
        unresolved_groups += label.contains_unresolved
        quarantine_groups += label.group_status == GROUP_STATUS_QUARANTINE
        if label.group_label.startswith("SINGLE_") and len(label.resolved_clades) == 1:
            clade_group_counts[label.resolved_clades[0]] += 1
    checks["groups_csv_matches_recomputation"] = True

    if label_conflicts != 0:
        raise RuntimeError(f"{label_conflicts} LABEL_CONFLICT group(s) present.")
    checks["no_label_conflict_groups"] = True

    # 6. counts vs manifest
    if (
        sg_manifest["n_split_groups"] != len(members_by_group)
        or sg_manifest["partition_eligible_groups"] != partition_eligible
        or sg_manifest["label_conflict_groups"] != label_conflicts
        or sg_manifest["groups_containing_unresolved"] != unresolved_groups
        or sg_manifest["clade_group_counts"] != dict(clade_group_counts)
        or sg_manifest["all_vs_all_combinatorial_coverage"] != 1.0
    ):
        raise RuntimeError("split-group manifest counts disagree with recomputation.")
    checks["manifest_counts_match_recomputation"] = True

    failed = sorted(name for name, ok in checks.items() if not ok)
    if failed:
        raise RuntimeError("clade reference validation failed: " + ", ".join(failed))

    return {
        "catalog_rows": len(catalog_rows),
        "midpiwi_regions": len(ok_ids),
        "split_groups": len(members_by_group),
        "partition_eligible_groups": partition_eligible,
        "clade_group_counts": dict(clade_group_counts),
        "label_conflict_groups": label_conflicts,
        "groups_containing_unresolved": unresolved_groups,
        "quarantine_groups": quarantine_groups,
        "all_vs_all_combinatorial_coverage": sg_manifest[
            "all_vs_all_combinatorial_coverage"
        ],
        "catalog_sha256": _sha256_file(catalog_path),
        "split_groups_csv_sha256": sg_manifest["split_groups_csv_sha256"],
        "checks": checks,
    }
