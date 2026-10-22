from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from src.pago_pipeline.query_recall_reference_panel import (
    DEFAULT_PANEL_DIRECTORY,
    QueryRecallReferencePanelRow,
    QueryRecallReferencePanelValidation,
    validate_query_recall_reference_panel,
)
from src.pago_pipeline.query_reference_recall import (
    RECALL_STATUS_EVALUABLE,
    RECALL_STATUS_NOT_EVALUABLE,
    ReferenceMatchMethod,
    _accession_is_admissible,
    _bare_accession,
    _recall,
    _uid_is_admissible,
    _uid_sort_key,
    build_matching_strategy_sha256,
    compute_query_reference_recall,
    normalize_protein_sequence,
    protein_sequence_sha256,
)

# Pinned so any change to the matching methodology — hierarchy, admissibility,
# attribution, normalization, strata, representatives, completeness policy —
# forces this test and the snapshot-reuse invariant to be revisited.
_EXPECTED_MATCHING_STRATEGY_SHA256 = (
    "e44ed7e341a25014e57729b1c054d2696b4aa7df77d845321580387ae77b76c7"
)

_RSAGO_ACCESSION = "ABP72561.1"
_TTAGO_ACCESSION = "WP_011174533.1"
_AFAGO_ACCESSION = "WP_010878815.1"
_SIAGO_ACCESSION = "WP_012735993.1"


def _reference_row(
    accession: str,
    *,
    sequence: str = "MREFERENCESEQ",
    ago_family: str = "PAGO",
    pago_clade: str = "LONG_B",
    protein_short_name: str = "Xx",
    verification_status: str = "verified",
    accession_curation_status: str = "RESOLVED",
) -> QueryRecallReferencePanelRow:
    return QueryRecallReferencePanelRow(
        accession=accession,
        protein_short_name=protein_short_name,
        organism="Synthetic organism",
        ago_family=ago_family,
        pago_clade=pago_clade,
        historical_query_recall_clade="",
        sequence_sha256=protein_sequence_sha256(sequence),
        sequence_length=len(sequence),
        uniprot_accession="",
        reference_label_source="synthetic",
        reference_label_evidence="EXPERIMENTAL",
        verification_status=verification_status,
        accession_curation_status=accession_curation_status,
        notes="",
    )


def _panel(
    rows: list[QueryRecallReferencePanelRow],
    **lock_overrides: object,
) -> QueryRecallReferencePanelValidation:
    """A validation object for the matcher's own boundary checks.

    Built directly rather than through the loader: the loader already rejects
    everything it owns (row count, composition, enums, reconciliations), so the
    only way to exercise P11's preconditions is to hand it a validation object.
    The committed P10 resource is never modified.
    """
    lock_payload: dict[str, object] = {
        "intended_use": "NCBI_QUERY_COVERAGE_AUDIT",
        "detector_validation_eligible": False,
    }
    lock_payload.update(lock_overrides)
    return QueryRecallReferencePanelValidation(
        panel_directory=DEFAULT_PANEL_DIRECTORY,
        lock_payload=lock_payload,
        rows=tuple(rows),
        panel_sha256="0" * 64,
        curation_notes_sha256="1" * 64,
    )


def _retrieved(
    rows: list[tuple[object, object, object]],
    *,
    columns: tuple[str, ...] = (
        "protein_uid",
        "gbseq__accession_version",
        "gbseq__sequence",
    ),
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            columns[0]: [row[0] for row in rows],
            columns[1]: [row[1] for row in rows],
            columns[2]: [row[2] for row in rows],
        }
    )
    return frame


def _detail(result, accession: str) -> dict[str, object]:
    return next(row for row in result.detail_rows if row["accession"] == accession)


def _summary(result, stratum: str) -> dict[str, object]:
    return next(row for row in result.summary_rows if row["stratum"] == stratum)


class CellContractTests(unittest.TestCase):
    """T-U1, T-U2, T-U2b — normalization, missingness, and type rejection."""

    def test_normalization_strips_all_whitespace_and_uppercases(self) -> None:
        self.assertEqual(normalize_protein_sequence("  m a\ta\na  "), "MAAA")
        self.assertEqual(
            protein_sequence_sha256(" maaa "), protein_sequence_sha256("MAAA")
        )

    def test_every_admitted_missing_scalar_normalizes_to_absent(self) -> None:
        string_series = pd.Series(["MAAA", None], dtype="string")
        object_series = pd.Series(["MAAA", None], dtype=object)
        missing_values = [
            None,
            float("nan"),
            np.float32("nan"),
            np.float64("nan"),
            pd.NA,
            pd.NaT,
            "",
            "   ",
            string_series.iloc[1],
            object_series.iloc[1],
        ]
        for value in missing_values:
            with self.subTest(value=repr(value)):
                normalized = normalize_protein_sequence(value)
                self.assertEqual(normalized, "")
                for leaked in ("<NA>", "NAN", "NAT", "NONE"):
                    self.assertNotIn(leaked, normalized)

    def test_missing_sequence_never_counts_as_evidence(self) -> None:
        # pd.NA used to stringify to "<NA>", get hashed, and count as present
        # sequence evidence. It must not.
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("5", "WP_9.1", pd.NA), ("6", "WP_8.1", pd.NaT)]
            ),
        )
        self.assertEqual(result.retrieved_row_count, 2)
        self.assertEqual(result.retrieved_rows_with_sequence_count, 0)
        self.assertFalse(result.sequence_evidence_complete)
        self.assertEqual(
            _detail(result, "WP_1.1")["match_method"], ReferenceMatchMethod.NONE.value
        )

    def test_non_scalar_cell_is_rejected(self) -> None:
        panel = _panel([_reference_row("WP_1.1")])
        for container in ([None], np.array([float("nan")])):
            with self.subTest(container=type(container).__name__):
                frame = pd.DataFrame(
                    {
                        "protein_uid": ["1"],
                        "gbseq__accession_version": ["WP_1.1"],
                        "gbseq__sequence": pd.Series([container], dtype=object),
                    }
                )
                with self.assertRaisesRegex(RuntimeError, "non-scalar"):
                    compute_query_reference_recall(
                        panel_validation=panel, retrieved_metadata_dataframe=frame
                    )

    def test_non_string_scalar_cell_is_rejected(self) -> None:
        panel = _panel([_reference_row("WP_1.1")])
        for column, value in (
            ("gbseq__sequence", 1.0),
            ("gbseq__accession_version", 3),
            ("protein_uid", True),
        ):
            with self.subTest(column=column):
                rows = {
                    "protein_uid": ["1"],
                    "gbseq__accession_version": ["WP_1.1"],
                    "gbseq__sequence": ["MAAA"],
                }
                rows[column] = [value]
                with self.assertRaisesRegex(RuntimeError, column):
                    compute_query_reference_recall(
                        panel_validation=panel,
                        retrieved_metadata_dataframe=pd.DataFrame(rows),
                    )


class IdentifierAdmissibilityTests(unittest.TestCase):
    """T-U3, T-U4 — accession and UID admissibility, and the total UID order."""

    def test_accession_admissibility(self) -> None:
        for accession in ("WP_1.2", "WP_1.10", "ABP72561.1", "A4WYU7.1"):
            self.assertTrue(_accession_is_admissible(accession), accession)
        for accession in ("ABP72561", "", "WP 1.1", "WP_1.", "!bad", "WP_1.1.1"):
            self.assertFalse(_accession_is_admissible(accession), accession)

    def test_bare_accession_strips_only_a_trailing_version(self) -> None:
        self.assertEqual(_bare_accession("WP_1.10"), "WP_1")
        self.assertEqual(_bare_accession("ABP72561.1"), "ABP72561")
        self.assertEqual(_bare_accession("ABP72561"), "ABP72561")

    def test_uid_admissibility_rejects_non_ascii_digits(self) -> None:
        for uid in ("1", "01", "9", "10", "2500461169"):
            self.assertTrue(_uid_is_admissible(uid), uid)
        # str.isdigit() accepts both of these; int() rejects the first and
        # accepts the second as 3. Neither is an NCBI UID.
        for uid in ("²", "٣", "abc", "", "1a"):
            self.assertFalse(_uid_is_admissible(uid), uid)

    def test_uid_sort_key_is_total_and_numeric(self) -> None:
        self.assertLess(_uid_sort_key("9"), _uid_sort_key("10"))
        # Numeric value first; raw text breaks a tie only when values are equal.
        self.assertLess(_uid_sort_key("01"), _uid_sort_key("1"))
        self.assertLess(_uid_sort_key("7"), _uid_sort_key("0012"))
        self.assertLess(_uid_sort_key("10"), _uid_sort_key("abc"))
        self.assertLess(_uid_sort_key("abc"), _uid_sort_key(""))
        for uid in ("²", "٣"):
            # Never routed through int(), so no exception and no numeric band.
            self.assertEqual(_uid_sort_key(uid)[0], 1)

    def test_recall_of_an_empty_denominator(self) -> None:
        self.assertEqual(_recall(0, 0), (None, RECALL_STATUS_NOT_EVALUABLE))
        self.assertEqual(_recall(2, 1), (0.5, RECALL_STATUS_EVALUABLE))


class MatchingStrategyIdentityTests(unittest.TestCase):
    def test_matching_strategy_sha256_is_pinned(self) -> None:
        self.assertEqual(
            build_matching_strategy_sha256(), _EXPECTED_MATCHING_STRATEGY_SHA256
        )


class CommittedPanelStratumTests(unittest.TestCase):
    """T-S1, T-S2, T-S4, T-S6 — invariants against the real committed panel."""

    def setUp(self) -> None:
        self.panel = validate_query_recall_reference_panel()
        self.result = compute_query_reference_recall(
            panel_validation=self.panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", _TTAGO_ACCESSION, "MAAA")]
            ),
        )

    def test_stratum_sizes_and_exhaustive_partition(self) -> None:
        expected = {
            "overall": 21,
            "LONG_A": 9,
            "LONG_B": 1,
            "SHORT": 2,
            "PAGO_UNRESOLVED_CLADE": 2,
            "PIWI_RE": 7,
        }
        for stratum, reference_count in expected.items():
            self.assertEqual(
                _summary(self.result, stratum)["reference_count"],
                reference_count,
                stratum,
            )
        partitioned = sum(
            row["reference_count"]
            for row in self.result.summary_rows
            if row["stratum"] != "overall"
        )
        self.assertEqual(partitioned, 21)

    def test_no_reference_falls_in_two_non_overall_strata(self) -> None:
        seen: dict[str, str] = {}
        for stratum, _metric, family, clade in (
            ("LONG_A", "", "PAGO", "LONG_A"),
            ("LONG_B", "", "PAGO", "LONG_B"),
            ("SHORT", "", "PAGO", "SHORT"),
            ("PAGO_UNRESOLVED_CLADE", "", "PAGO", "UNRESOLVED"),
            ("PIWI_RE", "", "PIWI_RE", None),
        ):
            for row in self.result.detail_rows:
                if row["ago_family"] == family and (
                    clade is None or row["pago_clade"] == clade
                ):
                    self.assertNotIn(row["accession"], seen)
                    seen[str(row["accession"])] = stratum
        self.assertEqual(len(seen), 21)

    def test_clade_strata_never_contain_a_piwi_re_row(self) -> None:
        for row in self.result.detail_rows:
            if row["ago_family"] == "PIWI_RE":
                self.assertEqual(row["pago_clade"], "UNRESOLVED")

    def test_pago_unresolved_clade_is_exactly_afago_and_siago(self) -> None:
        quarantined = sorted(
            str(row["accession"])
            for row in self.result.detail_rows
            if row["ago_family"] == "PAGO" and row["pago_clade"] == "UNRESOLVED"
        )
        self.assertEqual(quarantined, sorted([_AFAGO_ACCESSION, _SIAGO_ACCESSION]))

    def test_ttago_curation_uncertainty_reaches_the_detail(self) -> None:
        row = _detail(self.result, _TTAGO_ACCESSION)
        self.assertEqual(row["accession_curation_status"], "UNRESOLVED")
        self.assertEqual(row["verification_status"], "verified")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value
        )

    def test_tier_counts_sum_to_the_same_accession_numerator(self) -> None:
        for row in self.result.summary_rows:
            self.assertEqual(
                row["exact_version_match_count"] + row["same_base_match_count"],
                row["same_accession_recovered_count"],
                row["stratum"],
            )

    def test_recovered_same_accession_matches_the_tier_exactly(self) -> None:
        for row in self.result.detail_rows:
            expected = row["match_method"] in {
                ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value,
                ReferenceMatchMethod.SAME_BASE_ACCESSION.value,
            }
            self.assertEqual(row["recovered_same_accession"], expected)
            if row["match_method"] == ReferenceMatchMethod.SEQUENCE_SHA256.value:
                self.assertFalse(row["recovered_same_accession"])

    def test_no_duplicate_reference_sequence_groups_in_the_committed_panel(
        self,
    ) -> None:
        self.assertEqual(self.result.duplicate_reference_sequence_groups, ())


class MetricInterpretationTests(unittest.TestCase):
    """T-S3, T-S9 — what each reading does and does not demonstrate."""

    def test_same_base_match_does_not_establish_sequence_identity(self) -> None:
        # Both readings reach 1.0 while no sequence identity was established
        # anywhere: the equivalent reading is the union of the three admitted
        # tiers, not evidence that retrieved sequences equal reference
        # sequences. match_method is the only field that says how.
        panel = _panel([_reference_row("WP_X.1", sequence="TRUE_REFERENCE_SEQ")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", "WP_X.2", "COMPLETELY_DIFFERENT_SEQUENCE")]
            ),
        )
        row = _detail(result, "WP_X.1")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SAME_BASE_ACCESSION.value
        )
        self.assertEqual(
            result.stratum_same_accession_recall["long_b_reference_recall"], 1.0
        )
        self.assertEqual(
            result.stratum_equivalent_recall["long_b_reference_recall"], 1.0
        )

    def test_reference_recovered_by_sequence_under_an_alias_accession(self) -> None:
        # The RsAgo shape: the reference accession is absent from the retrieval
        # but its sequence is present under a different accession. Recovered
        # for the equivalent reading only, and the strict reading records the
        # accession miss.
        panel = _panel([_reference_row(_RSAGO_ACCESSION, sequence="MRSAGOSEQ")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("2500461169", "A4WYU7.1", "MRSAGOSEQ")]
            ),
        )
        row = _detail(result, _RSAGO_ACCESSION)
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SEQUENCE_SHA256.value
        )
        self.assertEqual(row["matched_accession"], "A4WYU7.1")
        self.assertEqual(row["matched_protein_uid"], "2500461169")
        self.assertEqual(row["sequence_match_count"], 1)
        self.assertTrue(row["recovered"])
        self.assertFalse(row["recovered_same_accession"])
        self.assertEqual(
            result.stratum_same_accession_recall["long_b_reference_recall"], 0.0
        )
        self.assertEqual(
            result.stratum_equivalent_recall["long_b_reference_recall"], 1.0
        )


class DuplicateReferenceSequenceTests(unittest.TestCase):
    """T-S8 — duplicate hashes are permitted and reported, never fatal."""

    def test_two_rows_sharing_one_sequence_are_both_recovered_and_reported(
        self,
    ) -> None:
        panel = _panel(
            [
                _reference_row("AAA111.1", sequence="SHAREDSEQ"),
                _reference_row("BBB222.1", sequence="SHAREDSEQ"),
            ]
        )
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("9", "CCC333.1", "SHAREDSEQ")]),
        )
        for accession in ("AAA111.1", "BBB222.1"):
            row = _detail(result, accession)
            self.assertTrue(row["recovered"])
            self.assertEqual(
                row["match_method"], ReferenceMatchMethod.SEQUENCE_SHA256.value
            )
            self.assertEqual(row["matched_accession"], "CCC333.1")
        self.assertEqual(
            result.duplicate_reference_sequence_groups,
            (("AAA111.1", "BBB222.1"),),
        )
        self.assertEqual(
            result.stratum_equivalent_recall["long_b_reference_recall"], 1.0
        )
        self.assertTrue(
            _summary(result, "LONG_B")["equivalent_recall_counts_shared_sequences"]
        )


class EmptyStratumTests(unittest.TestCase):
    """T-S7 — one representation for a stratum with no references."""

    def test_stratum_without_references_is_not_evaluable(self) -> None:
        panel = _panel([_reference_row("WP_1.1", pago_clade="LONG_A")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("1", "WP_1.1", "MAAA")]),
        )
        short = _summary(result, "SHORT")
        self.assertEqual(short["reference_count"], 0)
        self.assertEqual(short["exact_version_match_count"], 0)
        self.assertEqual(short["same_base_match_count"], 0)
        self.assertEqual(short["same_accession_recovered_count"], 0)
        self.assertEqual(short["equivalent_recovered_count"], 0)
        self.assertIsNone(short["same_accession_recall"])
        self.assertIsNone(short["retrieval_equivalent_recall"])
        self.assertEqual(short["recall_status"], RECALL_STATUS_NOT_EVALUABLE)
        self.assertIs(short["retrieval_equivalent_recall_is_lower_bound"], False)
        self.assertIs(short["equivalent_recall_counts_shared_sequences"], False)
        self.assertIsNone(result.stratum_equivalent_recall["short_reference_recall"])
        self.assertIs(
            result.stratum_equivalent_is_lower_bound["short_reference_recall"], False
        )
        # A present-but-unrecovered stratum stays EVALUABLE at 0.0.
        self.assertEqual(
            result.stratum_recall_status["long_a_reference_recall"],
            RECALL_STATUS_EVALUABLE,
        )


class RejectionTests(unittest.TestCase):
    """T-N1 … T-N7."""

    def test_rejects_unversioned_reference_accession(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "ABP72561"):
            compute_query_reference_recall(
                panel_validation=_panel([_reference_row("ABP72561")]),
                retrieved_metadata_dataframe=_retrieved([("1", "WP_1.1", "MAAA")]),
            )

    def test_rejects_duplicate_reference_base_accessions(self) -> None:
        panel = _panel(
            [_reference_row("WP_1.1"), _reference_row("WP_1.2", sequence="OTHERSEQ")]
        )
        with self.assertRaisesRegex(RuntimeError, "share a base accession"):
            compute_query_reference_recall(
                panel_validation=panel,
                retrieved_metadata_dataframe=_retrieved([("1", "WP_1.5", "MAAA")]),
            )

    def test_rejects_missing_sequence_column(self) -> None:
        frame = pd.DataFrame(
            {"protein_uid": ["1"], "gbseq__accession_version": ["WP_1.1"]}
        )
        with self.assertRaisesRegex(RuntimeError, "gbseq__sequence"):
            compute_query_reference_recall(
                panel_validation=_panel([_reference_row("WP_1.1")]),
                retrieved_metadata_dataframe=frame,
            )

    def test_rejects_missing_accession_column(self) -> None:
        frame = pd.DataFrame({"protein_uid": ["1"], "gbseq__sequence": ["MAAA"]})
        with self.assertRaisesRegex(RuntimeError, "gbseq__accession_version"):
            compute_query_reference_recall(
                panel_validation=_panel([_reference_row("WP_1.1")]),
                retrieved_metadata_dataframe=frame,
            )

    def test_rejects_conflicting_sequences_under_one_accession_version(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "WP_9.1"):
            compute_query_reference_recall(
                panel_validation=_panel([_reference_row("WP_1.1")]),
                retrieved_metadata_dataframe=_retrieved(
                    [("1", "WP_9.1", "MAAA"), ("2", "WP_9.1", "MZZZ")]
                ),
            )

    def test_rejects_a_panel_declared_eligible_for_detector_validation(self) -> None:
        panel = _panel([_reference_row("WP_1.1")], detector_validation_eligible=True)
        with self.assertRaisesRegex(RuntimeError, "detector_validation_eligible"):
            compute_query_reference_recall(
                panel_validation=panel,
                retrieved_metadata_dataframe=_retrieved([("1", "WP_1.1", "MAAA")]),
            )

    def test_rejects_a_panel_with_another_intended_use(self) -> None:
        panel = _panel([_reference_row("WP_1.1")], intended_use="SOMETHING_ELSE")
        with self.assertRaisesRegex(RuntimeError, "intended_use"):
            compute_query_reference_recall(
                panel_validation=panel,
                retrieved_metadata_dataframe=_retrieved([("1", "WP_1.1", "MAAA")]),
            )

    def test_truncated_retrieved_accession_does_not_recover_a_reference(self) -> None:
        # The historical defect: "ABP72561" satisfied the same-base tier for
        # "ABP72561.1" while carrying an unrelated sequence.
        panel = _panel([_reference_row(_RSAGO_ACCESSION, sequence="TRUE_RSAGO_SEQ")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("5", "ABP72561", "SOMETHING_ELSE")]
            ),
        )
        row = _detail(result, _RSAGO_ACCESSION)
        self.assertEqual(row["match_method"], ReferenceMatchMethod.NONE.value)
        self.assertFalse(row["recovered"])
        self.assertEqual(
            result.stratum_same_accession_recall["long_b_reference_recall"], 0.0
        )
        self.assertEqual(result.retrieved_inadmissible_accession_count, 1)


class DeterminismTests(unittest.TestCase):
    """T-D1, T-D3 … T-D10 — representatives and details independent of row order."""

    def _both_orders(self, panel, rows):
        forward = compute_query_reference_recall(
            panel_validation=panel, retrieved_metadata_dataframe=_retrieved(rows)
        )
        backward = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(list(reversed(rows))),
        )
        return forward, backward

    def test_shuffling_rows_changes_neither_summary_nor_detail(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        rows = [("900", "WP_1.1", "MAAA"), ("111", "WP_1.1", "MAAA")]
        forward, backward = self._both_orders(panel, rows)
        self.assertEqual(forward.summary_rows, backward.summary_rows)
        self.assertEqual(forward.detail_rows, backward.detail_rows)
        row = _detail(forward, "WP_1.1")
        self.assertEqual(row["matched_protein_uid"], "111")
        self.assertEqual(row["accession_version_distinct_uid_count"], 2)

    def test_exact_accession_takes_priority_over_the_sequence_tier(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="TRUESEQ")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", "WP_1.1", "MUTATEDSEQ"), ("2", "DECOY_1.1", "TRUESEQ")]
            ),
        )
        row = _detail(result, "WP_1.1")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value
        )
        self.assertTrue(row["recovered_same_accession"])

    def test_same_base_accession_takes_priority_over_the_sequence_tier(self) -> None:
        panel = _panel([_reference_row("WP_2.2", sequence="TRUESEQ")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", "WP_2.5", "MUTATEDSEQ"), ("2", "DECOY_2.1", "TRUESEQ")]
            ),
        )
        row = _detail(result, "WP_2.2")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SAME_BASE_ACCESSION.value
        )
        self.assertEqual(row["matched_accession"], "WP_2.5")

    def test_three_accessions_sharing_one_sequence_resolve_identically(self) -> None:
        panel = _panel([_reference_row("WP_5.1", sequence="SHARED")])
        rows = [
            ("30", "ZZ_9.1", "SHARED"),
            ("31", "AA_1.1", "SHARED"),
            ("32", "MM_5.1", "SHARED"),
        ]
        forward, backward = self._both_orders(panel, rows)
        for result in (forward, backward):
            row = _detail(result, "WP_5.1")
            self.assertEqual(row["matched_accession"], "AA_1.1")
            self.assertEqual(row["matched_protein_uid"], "31")
            self.assertEqual(row["sequence_match_count"], 3)

    def test_duplicated_accession_and_uid_pair_collapses(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("7", "WP_1.1", "MAAA"), ("7", "WP_1.1", "MAAA")]
            ),
        )
        row = _detail(result, "WP_1.1")
        self.assertEqual(row["accession_version_distinct_uid_count"], 1)
        self.assertEqual(row["matched_protein_uid"], "7")

    def test_duplicate_with_one_missing_sequence_is_not_a_conflict(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", "WP_9.1", "MAAA"), ("2", "WP_9.1", None)]
            ),
        )
        row = _detail(result, "WP_1.1")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SEQUENCE_SHA256.value
        )
        self.assertEqual(row["matched_accession"], "WP_9.1")

    def test_leading_zero_uids_are_distinguished_without_canonicalization(
        self,
    ) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        rows = [("1", "WP_1.1", "MAAA"), ("01", "WP_1.1", "MAAA")]
        forward, backward = self._both_orders(panel, rows)
        for result in (forward, backward):
            row = _detail(result, "WP_1.1")
            self.assertEqual(row["matched_protein_uid"], "01")
            self.assertEqual(row["accession_version_distinct_uid_count"], 2)

    def test_sequence_representative_prefers_an_admissible_accession(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="SHARED")])
        rows = [("5", "!bad", "SHARED"), ("7", "A.1", "SHARED"), ("9", "", "SHARED")]
        forward, backward = self._both_orders(panel, rows)
        for result in (forward, backward):
            row = _detail(result, "WP_1.1")
            self.assertEqual(row["matched_accession"], "A.1")
            self.assertEqual(row["matched_protein_uid"], "7")
            self.assertEqual(row["sequence_match_count"], 3)

    def test_uid_only_representative_orders_numerically(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="SHARED")])
        rows = [("7", "", "SHARED"), ("0012", "", "SHARED")]
        forward, backward = self._both_orders(panel, rows)
        for result in (forward, backward):
            row = _detail(result, "WP_1.1")
            self.assertEqual(row["matched_accession"], "")
            self.assertEqual(row["matched_protein_uid"], "7")


class SequenceEvidenceCompletenessTests(unittest.TestCase):
    """T-E1 … T-E8 — empty recovery, absent evidence, saturation."""

    def test_zero_retrieved_rows_is_a_valid_empty_retrieval(self) -> None:
        panel = _panel([_reference_row("WP_1.1")])
        result = compute_query_reference_recall(
            panel_validation=panel, retrieved_metadata_dataframe=_retrieved([])
        )
        self.assertEqual(result.retrieved_row_count, 0)
        self.assertTrue(result.sequence_evidence_complete)
        summary = _summary(result, "LONG_B")
        self.assertEqual(summary["retrieval_equivalent_recall"], 0.0)
        self.assertEqual(summary["recall_status"], RECALL_STATUS_EVALUABLE)
        self.assertIs(summary["retrieval_equivalent_recall_is_lower_bound"], False)

    def test_no_usable_sequences_flags_an_unsaturated_stratum(self) -> None:
        panel = _panel(
            [
                _reference_row("WP_1.1"),
                _reference_row("WP_2.1", sequence="OTHERSEQ"),
            ]
        )
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("1", "WP_1.1", None)]),
        )
        self.assertFalse(result.sequence_evidence_complete)
        summary = _summary(result, "LONG_B")
        self.assertEqual(summary["equivalent_recovered_count"], 1)
        self.assertEqual(summary["reference_count"], 2)
        self.assertIs(summary["retrieval_equivalent_recall_is_lower_bound"], True)

    def test_a_saturated_stratum_is_determined_despite_missing_sequences(self) -> None:
        panel = _panel([_reference_row("WP_1.1")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", "WP_1.1", None), ("2", "WP_7.1", None)]
            ),
        )
        self.assertFalse(result.sequence_evidence_complete)
        summary = _summary(result, "LONG_B")
        self.assertEqual(summary["retrieval_equivalent_recall"], 1.0)
        self.assertIs(summary["retrieval_equivalent_recall_is_lower_bound"], False)

    def test_partial_sequence_availability_flags_only_unsaturated_strata(
        self,
    ) -> None:
        panel = _panel(
            [
                _reference_row("WP_1.1", pago_clade="LONG_A"),
                _reference_row("WP_2.1", sequence="OTHERSEQ", pago_clade="SHORT"),
            ]
        )
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", "WP_1.1", "MAAA"), ("2", "WP_8.1", None)]
            ),
        )
        self.assertFalse(result.sequence_evidence_complete)
        self.assertIs(
            _summary(result, "LONG_A")["retrieval_equivalent_recall_is_lower_bound"],
            False,
        )
        self.assertIs(
            _summary(result, "SHORT")["retrieval_equivalent_recall_is_lower_bound"],
            True,
        )

    def test_unattributable_sequence_observation_is_missing_evidence(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("", "", "MAAA")]),
        )
        self.assertEqual(result.retrieved_unattributable_row_count, 1)
        self.assertEqual(result.retrieved_rows_with_sequence_count, 0)
        self.assertFalse(result.sequence_evidence_complete)
        self.assertEqual(
            _detail(result, "WP_1.1")["match_method"], ReferenceMatchMethod.NONE.value
        )

    def test_a_filled_token_is_not_automatically_a_uid(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("abc", "!bad", "MAAA")]),
        )
        self.assertEqual(result.retrieved_unattributable_row_count, 1)
        self.assertEqual(result.retrieved_rows_with_sequence_count, 0)
        self.assertEqual(
            _detail(result, "WP_1.1")["match_method"], ReferenceMatchMethod.NONE.value
        )

    def test_an_admissible_uid_alone_attributes_a_sequence_match(self) -> None:
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("123", "!bad", "MAAA")]),
        )
        row = _detail(result, "WP_1.1")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SEQUENCE_SHA256.value
        )
        self.assertEqual(row["matched_accession"], "")
        self.assertEqual(row["matched_protein_uid"], "123")
        self.assertEqual(result.retrieved_rows_with_sequence_count, 1)
        self.assertEqual(result.retrieved_unattributable_row_count, 0)
        self.assertEqual(result.retrieved_inadmissible_accession_count, 1)

    def test_progress_message_does_not_claim_a_bound_for_a_saturated_result(
        self,
    ) -> None:
        # Every reference recovered by accession, no sequence evidence at all:
        # the equivalent reading is determined, so the message must not call it
        # a lower bound even though the evidence is incomplete.
        panel = _panel([_reference_row("WP_1.1")])
        messages: list[str] = []
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved(
                [("1", "WP_1.1", None), ("2", "WP_7.1", None)]
            ),
            progress_callback=messages.append,
        )
        self.assertFalse(result.sequence_evidence_complete)
        self.assertIs(
            result.stratum_equivalent_is_lower_bound["overall_reference_recall"],
            False,
        )
        self.assertEqual(len(messages), 1)
        self.assertNotIn("lower bound", messages[0])
        self.assertIn("retrieval-equivalent 1/1", messages[0])

    def test_progress_message_reports_a_bound_for_an_unsaturated_result(
        self,
    ) -> None:
        panel = _panel(
            [
                _reference_row("WP_1.1"),
                _reference_row("WP_2.1", sequence="OTHERSEQ"),
            ]
        )
        messages: list[str] = []
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("1", "WP_1.1", None)]),
            progress_callback=messages.append,
        )
        self.assertIs(
            result.stratum_equivalent_is_lower_bound["overall_reference_recall"],
            True,
        )
        self.assertIn("lower bound", messages[0])

    def test_progress_message_omits_the_bound_when_evidence_is_complete(
        self,
    ) -> None:
        panel = _panel(
            [
                _reference_row("WP_1.1"),
                _reference_row("WP_2.1", sequence="OTHERSEQ"),
            ]
        )
        messages: list[str] = []
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("1", "WP_1.1", "MAAA")]),
            progress_callback=messages.append,
        )
        self.assertTrue(result.sequence_evidence_complete)
        self.assertNotIn("lower bound", messages[0])
        self.assertIn("same-accession 1/2", messages[0])

    def test_an_unrecovered_reference_is_qualified_by_the_inadmissible_count(
        self,
    ) -> None:
        # A negative at the strict reading is not proof the record was absent:
        # the retrieval carried an identifier that could not be used.
        panel = _panel([_reference_row("WP_1.1", sequence="MAAA")])
        result = compute_query_reference_recall(
            panel_validation=panel,
            retrieved_metadata_dataframe=_retrieved([("5", "WP_1", "OTHER")]),
        )
        self.assertEqual(
            result.stratum_same_accession_recall["long_b_reference_recall"], 0.0
        )
        self.assertEqual(result.retrieved_inadmissible_accession_count, 1)


if __name__ == "__main__":
    unittest.main()
