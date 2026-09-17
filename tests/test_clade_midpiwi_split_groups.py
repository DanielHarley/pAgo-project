from __future__ import annotations

import hashlib
import unittest

from src.pago_pipeline.clade_midpiwi_split_groups import (
    assert_complete_all_vs_all,
    build_split_groups,
    expected_unordered_pairs,
    label_split_group,
    score_aligned_pair,
)


def _sha(token: str) -> str:
    return hashlib.sha256(token.encode("ascii")).hexdigest()


class MidpiwiEdgeRuleTests(unittest.TestCase):
    """The 90/80 edge rule ADR-0013 applies to MID-PIWI regions."""

    def test_valid_90_80_edge(self) -> None:
        query = "ACDEFGHIKL" * 10
        target = ("ACDEFGHIKL" * 9) + "ACDEFGHIKM"  # 99/100 identical
        comparable, identical, identity, coverage, is_edge = score_aligned_pair(
            query_aligned=query, target_aligned=target,
            query_length=100, target_length=100,
        )
        self.assertEqual((comparable, identical), (100, 99))
        self.assertAlmostEqual(identity, 0.99)
        self.assertAlmostEqual(coverage, 1.0)
        self.assertTrue(is_edge)

    def test_rejected_by_identity(self) -> None:
        query = "A" * 100
        target = ("A" * 85) + ("C" * 15)  # 85% identical over full overlap
        *_, identity, coverage, is_edge = score_aligned_pair(
            query_aligned=query, target_aligned=target,
            query_length=100, target_length=100,
        )
        self.assertAlmostEqual(identity, 0.85)
        self.assertAlmostEqual(coverage, 1.0)
        self.assertFalse(is_edge)

    def test_rejected_by_coverage(self) -> None:
        # 100 fully identical comparable columns, but the shorter member is 150 aa
        query = "A" * 100
        target = "A" * 100
        *_, identity, coverage, is_edge = score_aligned_pair(
            query_aligned=query, target_aligned=target,
            query_length=400, target_length=150,
        )
        self.assertAlmostEqual(identity, 1.0)
        self.assertAlmostEqual(coverage, 100 / 150)
        self.assertFalse(is_edge)


class MidpiwiConnectedComponentTests(unittest.TestCase):
    def test_exact_duplicate_regions_are_unioned_without_an_edge(self) -> None:
        result = build_split_groups(
            sequence_ids=["WP_a.1", "WP_b.1", "WP_c.1"],
            sequence_sha256_by_id={
                "WP_a.1": _sha("REGION_ONE"),
                "WP_b.1": _sha("REGION_ONE"),  # identical MID-PIWI region
                "WP_c.1": _sha("REGION_TWO"),
            },
            edge_pairs=[],
            split_group_id_prefix="MIDPIWI",
        )
        self.assertEqual(
            result.split_group_id_by_sequence_id["WP_a.1"],
            result.split_group_id_by_sequence_id["WP_b.1"],
        )
        self.assertNotEqual(
            result.split_group_id_by_sequence_id["WP_a.1"],
            result.split_group_id_by_sequence_id["WP_c.1"],
        )
        self.assertEqual(len(result.members_by_split_group_id), 2)

    def test_transitive_edges_form_one_component(self) -> None:
        result = build_split_groups(
            sequence_ids=["WP_a.1", "WP_b.1", "WP_c.1", "WP_d.1"],
            sequence_sha256_by_id={i: _sha(i) for i in
                                   ("WP_a.1", "WP_b.1", "WP_c.1", "WP_d.1")},
            edge_pairs=[("WP_a.1", "WP_b.1"), ("WP_b.1", "WP_c.1")],  # A-B-C, no A-C
            split_group_id_prefix="MIDPIWI",
        )
        groups = {frozenset(m) for m in result.members_by_split_group_id.values()}
        self.assertIn(frozenset({"WP_a.1", "WP_b.1", "WP_c.1"}), groups)
        self.assertIn(frozenset({"WP_d.1"}), groups)

    def test_no_sequence_lost_and_groups_are_mutually_exclusive(self) -> None:
        ids = [f"WP_{i}.1" for i in range(20)]
        sha = {i: _sha(i) for i in ids}
        result = build_split_groups(
            sequence_ids=ids, sequence_sha256_by_id=sha,
            edge_pairs=[("WP_0.1", "WP_1.1"), ("WP_2.1", "WP_3.1")],
            split_group_id_prefix="MIDPIWI",
        )
        covered = [m for members in result.members_by_split_group_id.values() for m in members]
        self.assertEqual(sorted(covered), sorted(ids))
        self.assertEqual(len(covered), len(set(covered)))


class LabelSplitGroupTests(unittest.TestCase):
    def _clades(self, mapping):
        return {"member_ids": list(mapping), "curated_pago_clade_by_id": mapping}

    def test_two_resolved_clades_is_label_conflict_and_not_forced(self) -> None:
        label = label_split_group(**self._clades({
            "WP_1.1": "LONG_A", "WP_2.1": "LONG_B", "WP_3.1": "LONG_A",
        }))
        self.assertEqual(label.group_label, "LABEL_CONFLICT")
        self.assertEqual(label.group_status, "QUARANTINE")
        self.assertFalse(label.partition_eligible)
        self.assertEqual(set(label.resolved_clades), {"LONG_A", "LONG_B"})

    def test_single_resolved_plus_unresolved_is_review_not_conflict(self) -> None:
        label = label_split_group(**self._clades({
            "WP_1.1": "LONG_B", "WP_2.1": "UNRESOLVED",
        }))
        self.assertEqual(label.group_label, "SINGLE_RESOLVED_PLUS_UNRESOLVED")
        self.assertEqual(label.group_status, "REVIEW")
        self.assertFalse(label.partition_eligible)
        self.assertNotEqual(label.group_label, "LABEL_CONFLICT")

    def test_single_resolved_clade_is_ok_and_partition_eligible(self) -> None:
        label = label_split_group(**self._clades({
            "WP_1.1": "SHORT", "WP_2.1": "SHORT", "WP_3.1": "SHORT",
        }))
        self.assertEqual(label.group_label, "SINGLE_SHORT")
        self.assertEqual(label.group_status, "OK")
        self.assertTrue(label.partition_eligible)

    def test_unresolved_only_group_is_review(self) -> None:
        label = label_split_group(**self._clades({"WP_1.1": "UNRESOLVED"}))
        self.assertEqual(label.group_label, "UNRESOLVED_ONLY")
        self.assertEqual(label.group_status, "REVIEW")
        self.assertFalse(label.partition_eligible)

    def test_afago_quarantine_anchor_propagates_requires_manual_review(self) -> None:
        # AfAgo (WP_010878815.1) sitting inside an otherwise single-clade group
        label = label_split_group(**self._clades({
            "WP_010878815.1": "UNRESOLVED", "WP_x.1": "LONG_B",
        }))
        self.assertTrue(label.requires_manual_review)
        self.assertEqual(label.quarantine_anchors, ("AfAgo",))
        self.assertEqual(label.group_status, "QUARANTINE")
        self.assertFalse(label.partition_eligible)

    def test_siago_quarantine_anchor_alone_is_quarantined(self) -> None:
        label = label_split_group(**self._clades({"WP_012735993.1": "UNRESOLVED"}))
        self.assertEqual(label.quarantine_anchors, ("SiAgo",))
        self.assertEqual(label.group_status, "QUARANTINE")
        self.assertFalse(label.partition_eligible)


class AllVsAllCoverageTests(unittest.TestCase):
    def test_complete_all_vs_all_passes(self) -> None:
        assert_complete_all_vs_all(observed_unordered_pairs=6, n_sequences=4)
        assert_complete_all_vs_all(
            observed_unordered_pairs=expected_unordered_pairs(1002),
            n_sequences=1002,
        )

    def test_incomplete_all_vs_all_aborts(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "incomplete"):
            assert_complete_all_vs_all(observed_unordered_pairs=5, n_sequences=4)
        with self.assertRaisesRegex(RuntimeError, "incomplete"):
            assert_complete_all_vs_all(
                observed_unordered_pairs=501500, n_sequences=1002
            )


class DiagnosticThresholdIndependenceTests(unittest.TestCase):
    """A diagnostic pass at another threshold never changes the 90/80 grouping,
    and never changes any curated_pago_clade/curation_status label (ADR-0013:
    diagnostic passes are evidence only, never automatic relabeling)."""

    def _edges_at(self, scored, identity_threshold, coverage_threshold=0.80):
        return [
            key for key, (identity, coverage) in scored.items()
            if identity >= identity_threshold and coverage >= coverage_threshold
        ]

    def test_authoritative_grouping_is_independent_of_diagnostic_passes(self) -> None:
        ids = ["WP_a.1", "WP_b.1", "WP_c.1", "WP_d.1"]
        sha = {i: _sha(i) for i in ids}
        # a-b are 95% identical (edge at 90); b-c are 75% (edge only at <=70);
        # d is unrelated.
        scored = {
            ("WP_a.1", "WP_b.1"): (0.95, 1.0),
            ("WP_b.1", "WP_c.1"): (0.75, 1.0),
            ("WP_a.1", "WP_c.1"): (0.72, 1.0),
            ("WP_a.1", "WP_d.1"): (0.20, 1.0),
            ("WP_b.1", "WP_d.1"): (0.20, 1.0),
            ("WP_c.1", "WP_d.1"): (0.20, 1.0),
        }
        authoritative = build_split_groups(
            sequence_ids=ids, sequence_sha256_by_id=sha,
            edge_pairs=self._edges_at(scored, 0.90),
            split_group_id_prefix="MIDPIWI",
        )
        authoritative_partition = {
            frozenset(m) for m in authoritative.members_by_split_group_id.values()
        }
        clades_before = {"WP_a.1": "LONG_A", "WP_b.1": "LONG_A", "WP_c.1": "LONG_B", "WP_d.1": "SHORT"}
        labels_before = {
            m: label_split_group(member_ids=[m], curated_pago_clade_by_id=clades_before)
            for m in ids
        }
        # run the diagnostic passes
        for threshold in (1.00, 0.90, 0.70, 0.50):
            build_split_groups(
                sequence_ids=ids, sequence_sha256_by_id=sha,
                edge_pairs=self._edges_at(scored, threshold),
                split_group_id_prefix="DIAG",
            )
        # recompute the authoritative grouping: unchanged
        again = build_split_groups(
            sequence_ids=ids, sequence_sha256_by_id=sha,
            edge_pairs=self._edges_at(scored, 0.90),
            split_group_id_prefix="MIDPIWI",
        )
        self.assertEqual(
            authoritative_partition,
            {frozenset(m) for m in again.members_by_split_group_id.values()},
        )
        # and it is {a,b} | {c} | {d}, not the coarser {a,b,c} the 70% pass gives
        self.assertIn(frozenset({"WP_a.1", "WP_b.1"}), authoritative_partition)
        self.assertIn(frozenset({"WP_c.1"}), authoritative_partition)
        coarser = build_split_groups(
            sequence_ids=ids, sequence_sha256_by_id=sha,
            edge_pairs=self._edges_at(scored, 0.70),
            split_group_id_prefix="DIAG",
        )
        self.assertIn(
            frozenset({"WP_a.1", "WP_b.1", "WP_c.1"}),
            {frozenset(m) for m in coarser.members_by_split_group_id.values()},
        )
        # curated_pago_clade/curation_status labels are untouched by any of the
        # diagnostic passes above (the diagnostic result was never fed back in)
        labels_after = {
            m: label_split_group(member_ids=[m], curated_pago_clade_by_id=clades_before)
            for m in ids
        }
        self.assertEqual(labels_before, labels_after)


if __name__ == "__main__":
    unittest.main()
