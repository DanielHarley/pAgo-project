from __future__ import annotations

import hashlib
import unittest
from pathlib import Path

from src.pago_pipeline.apaz_split_groups import (
    AlignedPair,
    build_split_groups,
    partition_split_groups,
    score_aligned_pair,
    score_pairs,
    selection_key_sha256,
    validate_apaz_partition_invariants,
    validate_apaz_seed_consistency,
    _subset_sum_first,
)

APAZ_SEED_DIRECTORY = (
    Path(__file__).resolve().parents[1]
    / "src" / "pago_pipeline" / "resources" / "apaz_seed"
)


def _sha(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


class ScoreRuleTests(unittest.TestCase):
    def test_identity_is_over_comparable_columns_and_coverage_over_smaller(self):
        # 6/10 columns comparable, 6 identical -> identity 1.0
        # coverage 6 / min(10, 10) = 0.6  -> below 0.80 -> not an edge
        comparable, identical, identity, coverage, is_edge = score_aligned_pair(
            query_aligned="ACDEFGHIKL",
            target_aligned="ACDEFG----",
            query_length=10,
            target_length=10,
        )
        self.assertEqual((comparable, identical), (6, 6))
        self.assertEqual(identity, 1.0)
        self.assertAlmostEqual(coverage, 0.6)
        self.assertFalse(is_edge)

    def test_coverage_is_against_the_shorter_sequence(self):
        # target is a 6-residue fragment fully aligned -> coverage 6/6 = 1.0
        *_, is_edge = score_aligned_pair(
            query_aligned="ACDEFGHIKL",
            target_aligned="ACDEFG----",
            query_length=10,
            target_length=6,
        )
        self.assertTrue(is_edge)

    def test_below_identity_threshold_is_not_an_edge(self):
        *_, identity, _, is_edge = score_aligned_pair(
            query_aligned="ACDEFGHIKL",
            target_aligned="ACDEFGHIKW",  # 9/10
            query_length=10,
            target_length=10,
        )
        self.assertEqual(identity, 0.9)
        self.assertTrue(is_edge)  # 0.9 >= 0.90 and coverage 1.0
        *_, identity2, _, is_edge2 = score_aligned_pair(
            query_aligned="ACDEFGHIKL",
            target_aligned="ACDEFGHWWW",  # 7/10
            query_length=10,
            target_length=10,
        )
        self.assertEqual(identity2, 0.7)
        self.assertFalse(is_edge2)

    def test_mismatched_alignment_widths_raise(self):
        with self.assertRaises(ValueError):
            score_aligned_pair(query_aligned="ACD", target_aligned="ACDE",
                               query_length=3, target_length=4)


class BuildSplitGroupsTests(unittest.TestCase):
    def _sha_map(self, ids):
        return {i: _sha(i) for i in ids}

    def test_connected_components_are_transitive_not_cliques(self):
        ids = ["A", "B", "C", "D"]
        # A-B and B-C are edges; A-C is not.  D is isolated.
        res = build_split_groups(
            sequence_ids=ids, sequence_sha256_by_id=self._sha_map(ids),
            edge_pairs=[("A", "B"), ("B", "C")], split_group_id_prefix="T")
        groups = {frozenset(m) for m in res.members_by_split_group_id.values()}
        self.assertEqual(groups, {frozenset({"A", "B", "C"}), frozenset({"D"})})

    def test_exact_duplicates_are_unioned_without_any_edge(self):
        ids = ["A", "A_dup", "B"]
        sha = {"A": _sha("same"), "A_dup": _sha("same"), "B": _sha("other")}
        res = build_split_groups(
            sequence_ids=ids, sequence_sha256_by_id=sha,
            edge_pairs=[], split_group_id_prefix="T")
        self.assertEqual(res.split_group_id_by_sequence_id["A"],
                         res.split_group_id_by_sequence_id["A_dup"])
        self.assertNotEqual(res.split_group_id_by_sequence_id["A"],
                            res.split_group_id_by_sequence_id["B"])

    def test_result_is_invariant_to_input_order(self):
        ids = ["S1", "S2", "S3", "S4", "S5"]
        sha = self._sha_map(ids)
        edges = [("S1", "S3"), ("S3", "S5")]
        a = build_split_groups(sequence_ids=ids, sequence_sha256_by_id=sha,
                               edge_pairs=edges, split_group_id_prefix="T")
        b = build_split_groups(sequence_ids=list(reversed(ids)),
                               sequence_sha256_by_id=sha,
                               edge_pairs=list(reversed(edges)),
                               split_group_id_prefix="T")
        self.assertEqual(a.split_group_id_by_sequence_id, b.split_group_id_by_sequence_id)
        self.assertEqual(a.members_by_split_group_id, b.members_by_split_group_id)

    def test_split_group_id_derives_from_sorted_members_not_input_order(self):
        ids = ["B", "A"]
        sha = {"A": _sha("x"), "B": _sha("x")}  # duplicates -> one group
        res = build_split_groups(sequence_ids=ids, sequence_sha256_by_id=sha,
                                 edge_pairs=[], split_group_id_prefix="T")
        expected = "T_" + hashlib.sha256("A\nB".encode()).hexdigest()[:16]
        self.assertEqual(list(res.members_by_split_group_id), [expected])

    def test_edge_to_unknown_id_raises(self):
        ids = ["A", "B"]
        with self.assertRaises(ValueError):
            build_split_groups(sequence_ids=ids, sequence_sha256_by_id=self._sha_map(ids),
                               edge_pairs=[("A", "Z")], split_group_id_prefix="T")


class SubsetSumTests(unittest.TestCase):
    def test_first_solution_prefers_earlier_items(self):
        items = [("g1", 2), ("g2", 1), ("g3", 1), ("g4", 2)]
        self.assertEqual(_subset_sum_first(items, 3), ["g1", "g2"])

    def test_none_when_no_subset_reaches_target(self):
        self.assertIsNone(_subset_sum_first([("a", 2), ("b", 2)], 3))

    def test_zero_target(self):
        self.assertEqual(_subset_sum_first([("a", 1)], 0), [])


class PartitionTests(unittest.TestCase):
    def _members(self):
        # stratum X: 4 singletons + 1 pair (size 6); stratum Y: 6 singletons
        return {
            "X_s1": ("x1",), "X_s2": ("x2",), "X_s3": ("x3",), "X_s4": ("x4",),
            "X_p1": ("x5", "x6"),
            "Y_s1": ("y1",), "Y_s2": ("y2",), "Y_s3": ("y3",),
            "Y_s4": ("y4",), "Y_s5": ("y5",), "Y_s6": ("y6",),
        }

    def _stratum(self):
        return {**{f"x{i}": "X" for i in range(1, 7)},
                **{f"y{i}": "Y" for i in range(1, 7)}}

    def test_exact_targets_stratified_and_whole_groups(self):
        res = partition_split_groups(
            members_by_split_group_id=self._members(),
            stratum_by_sequence_id=self._stratum(),
            fill_order_by_stratum={
                "X": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)],
                "Y": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)],
            },
            partition_salt="salt-v2",
        )
        from collections import Counter
        cx = Counter(res.partition_by_sequence_id[s] for s in self._stratum() if s.startswith("x"))
        cy = Counter(res.partition_by_sequence_id[s] for s in self._stratum() if s.startswith("y"))
        self.assertEqual(dict(cx), {"FINAL_HOLDOUT": 2, "CALIBRATION": 2, "BUILD": 2})
        self.assertEqual(dict(cy), {"FINAL_HOLDOUT": 2, "CALIBRATION": 2, "BUILD": 2})
        # the pair x5/x6 is never split across partitions
        self.assertEqual(res.partition_by_sequence_id["x5"], res.partition_by_sequence_id["x6"])

    def test_partition_is_invariant_to_group_and_sequence_order(self):
        members = self._members()
        a = partition_split_groups(
            members_by_split_group_id=members, stratum_by_sequence_id=self._stratum(),
            fill_order_by_stratum={"X": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)],
                                   "Y": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)]},
            partition_salt="salt-v2")
        shuffled = {k: tuple(reversed(v)) for k, v in reversed(list(members.items()))}
        b = partition_split_groups(
            members_by_split_group_id=shuffled, stratum_by_sequence_id=self._stratum(),
            fill_order_by_stratum={"Y": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)],
                                   "X": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)]},
            partition_salt="salt-v2")
        self.assertEqual(a.partition_by_sequence_id, b.partition_by_sequence_id)
        self.assertEqual(a.partition_by_split_group_id, b.partition_by_split_group_id)

    def test_salt_changes_the_assignment(self):
        kw = dict(members_by_split_group_id=self._members(),
                  stratum_by_sequence_id=self._stratum(),
                  fill_order_by_stratum={"X": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)],
                                         "Y": [("FINAL_HOLDOUT", 2), ("CALIBRATION", 2), ("BUILD", 2)]})
        a = partition_split_groups(partition_salt="salt-a", **kw)
        b = partition_split_groups(partition_salt="salt-b", **kw)
        self.assertNotEqual(a.partition_by_split_group_id, b.partition_by_split_group_id)

    def test_aborts_when_a_target_is_not_subset_sum_reachable(self):
        members = {"X_p1": ("x1", "x2"), "X_p2": ("x3", "x4")}  # only sizes {2, 2}
        with self.assertRaisesRegex(RuntimeError, "no combination of whole split groups"):
            partition_split_groups(
                members_by_split_group_id=members,
                stratum_by_sequence_id={f"x{i}": "X" for i in range(1, 5)},
                fill_order_by_stratum={"X": [("FINAL_HOLDOUT", 1), ("BUILD", 3)]},
                partition_salt="salt-v2")

    def test_split_group_crossing_strata_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "crosses strata"):
            partition_split_groups(
                members_by_split_group_id={"g": ("x1", "y1")},
                stratum_by_sequence_id={"x1": "X", "y1": "Y"},
                fill_order_by_stratum={"X": [("A", 1), ("B", 0)]},
                partition_salt="salt-v2")


@unittest.skipUnless(
    (APAZ_SEED_DIRECTORY / "apaz_partitions.csv").is_file(),
    "committed APAZ B2 resources are not present",
)
class CommittedResourceTests(unittest.TestCase):
    def test_partition_invariants_hold_for_the_committed_resources(self):
        checks = validate_apaz_partition_invariants(
            partitions_csv_path=APAZ_SEED_DIRECTORY / "apaz_partitions.csv",
            split_groups_directory=APAZ_SEED_DIRECTORY / "split_groups",
        )
        self.assertEqual(len(checks), 14)
        self.assertTrue(all(checks.values()))

    def test_seed_consistency_holds_for_the_committed_resources(self):
        checks = validate_apaz_seed_consistency(resource_directory=APAZ_SEED_DIRECTORY)
        self.assertEqual(len(checks), 7)
        self.assertTrue(all(checks.values()))

    def test_committed_split_groups_rebuild_from_the_frozen_edge_list(self):
        import csv

        for dataset, prefix in (("apaz", "APAZ"), ("hisg", "HISG"), ("eiib", "EIIB")):
            directory = APAZ_SEED_DIRECTORY / "split_groups" / dataset
            sha_by_id = {}
            committed = {}
            with (directory / "split_groups.csv").open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    sha_by_id[row["sequence_id"]] = row["sequence_sha256"]
                    committed.setdefault(row["split_group_id"], []).append(row["sequence_id"])
            edges = []
            with (directory / "pairs.filtered.tsv").open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    edges.append(tuple(sorted((row["seq_a"], row["seq_b"]))))
            result = build_split_groups(
                sequence_ids=sorted(sha_by_id), sequence_sha256_by_id=sha_by_id,
                edge_pairs=sorted(set(edges)), split_group_id_prefix=prefix)
            self.assertEqual(
                result.members_by_split_group_id,
                {group: tuple(sorted(ids)) for group, ids in committed.items()},
            )


class SelectionKeyTests(unittest.TestCase):
    def test_selection_key_is_stable_and_field_sensitive(self):
        base = dict(partition_salt="s", stratum="Ia", split_group_id="g1", partition="BUILD")
        self.assertEqual(selection_key_sha256(**base), selection_key_sha256(**base))
        self.assertNotEqual(selection_key_sha256(**base),
                            selection_key_sha256(**{**base, "partition": "CALIBRATION"}))


if __name__ == "__main__":
    unittest.main()
