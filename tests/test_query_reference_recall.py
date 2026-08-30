from __future__ import annotations

import unittest
from pathlib import Path

import pandas as pd

from src.pago_pipeline.query_reference_recall import compute_query_reference_recall

_FIXTURE_CSV = (
    Path(__file__).resolve().parent / "fixtures" / "query_recall_reference_set.csv"
)


def _reference_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accession": "WP_001.1",
                "ago_family": "PAGO",
                "clade": "LONG_A",
                "reference_label_source": "src",
                "reference_label_evidence": "EXPERIMENTAL",
            },
            {
                "accession": "WP_002.2",
                "ago_family": "PAGO",
                "clade": "LONG_A",
                "reference_label_source": "src",
                "reference_label_evidence": "EXPERIMENTAL",
            },
            {
                "accession": "WP_003.1",
                "ago_family": "PAGO",
                "clade": "LONG_B",
                "reference_label_source": "src",
                "reference_label_evidence": "EXPERIMENTAL",
            },
            {
                "accession": "WP_004.1",
                "ago_family": "PAGO",
                "clade": "SHORT",
                "reference_label_source": "src",
                "reference_label_evidence": "EXPERIMENTAL",
            },
            {
                "accession": "WP_005.1",
                "ago_family": "PIWI_RE",
                "clade": "PIWI_RE",
                "reference_label_source": "src",
                "reference_label_evidence": "LITERATURE_PHYLOGENETIC",
            },
        ]
    )


class QueryReferenceRecallLogicTests(unittest.TestCase):
    def test_version_tolerant_matching_and_stratified_recall(self) -> None:
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["11", "12", "13"],
                # WP_001.1 exact; WP_002.9 matches WP_002.2 via bare accession;
                # WP_003 absent; WP_004 absent; WP_005 absent.
                "gbseq__accession_version": ["WP_001.1", "WP_002.9", "WP_999.1"],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        self.assertEqual(result.reference_count, 5)
        self.assertEqual(result.recovered_count, 2)
        self.assertAlmostEqual(
            result.stratum_recall["overall_reference_recall"], 2 / 5
        )
        self.assertAlmostEqual(
            result.stratum_recall["long_a_reference_recall"], 2 / 2
        )
        # Every clade has references in this fixture, so 0 recovered means an
        # EVALUABLE recall of 0.0 (never confused with NOT_EVALUABLE).
        self.assertAlmostEqual(
            result.stratum_recall["long_b_reference_recall"], 0.0
        )
        self.assertAlmostEqual(result.stratum_recall["short_reference_recall"], 0.0)
        self.assertAlmostEqual(
            result.stratum_recall["piwi_re_reference_recall"], 0.0
        )
        for metric_name in (
            "overall_reference_recall",
            "long_a_reference_recall",
            "long_b_reference_recall",
            "short_reference_recall",
            "piwi_re_reference_recall",
        ):
            self.assertEqual(
                result.stratum_recall_status[metric_name], "EVALUABLE"
            )

        detail_by_accession = dict(
            zip(
                result.detail_dataframe["accession"],
                result.detail_dataframe["recovered"],
            )
        )
        self.assertTrue(detail_by_accession["WP_001.1"])
        self.assertTrue(detail_by_accession["WP_002.2"])
        self.assertFalse(detail_by_accession["WP_003.1"])

    def test_empty_clade_is_not_evaluable_not_zero(self) -> None:
        reference_dataframe = _reference_dataframe()
        reference_dataframe = reference_dataframe[
            reference_dataframe["clade"] != "SHORT"
        ]
        result = compute_query_reference_recall(
            reference_dataframe=reference_dataframe,
            retrieved_metadata_dataframe=pd.DataFrame(
                {
                    "protein_uid": ["1"],
                    "gbseq__accession_version": ["WP_001.1"],
                }
            ),
        )
        self.assertIsNone(result.stratum_recall["short_reference_recall"])
        self.assertEqual(
            result.stratum_recall_status["short_reference_recall"],
            "NOT_EVALUABLE",
        )
        short_row = result.summary_dataframe[
            result.summary_dataframe["stratum"] == "SHORT"
        ].iloc[0]
        self.assertEqual(short_row["reference_count"], 0)
        self.assertTrue(pd.isna(short_row["recall"]))
        self.assertEqual(short_row["recall_status"], "NOT_EVALUABLE")
        # A present-but-unrecovered clade is still EVALUABLE at 0.0.
        self.assertEqual(
            result.stratum_recall_status["long_b_reference_recall"], "EVALUABLE"
        )
        self.assertEqual(result.stratum_recall["long_b_reference_recall"], 0.0)

    def test_missing_reference_columns_raise(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "missing required columns"):
            compute_query_reference_recall(
                reference_dataframe=pd.DataFrame([{"accession": "x"}]),
                retrieved_metadata_dataframe=pd.DataFrame(
                    {"gbseq__accession_version": ["x"]}
                ),
            )

    def test_committed_reference_set_is_well_formed(self) -> None:
        reference_dataframe = pd.read_csv(_FIXTURE_CSV, dtype=str).fillna("")
        for column_name in (
            "accession",
            "ago_family",
            "clade",
            "reference_label_source",
            "reference_label_evidence",
            "verification_status",
        ):
            self.assertIn(column_name, reference_dataframe.columns)
        self.assertGreaterEqual(len(reference_dataframe), 20)
        self.assertTrue((reference_dataframe["accession"].str.len() > 0).all())
        self.assertTrue((reference_dataframe["accession"].str.contains(r"\.\d+$")).all())
        clades = set(reference_dataframe["clade"].str.upper())
        self.assertTrue(
            clades <= {"LONG_A", "LONG_B", "SHORT", "PIWI_RE", "UNRESOLVED", "NA", ""}
        )
        # PIWI-RE stratum keys on ago_family and must be non-empty (EVALUABLE).
        self.assertGreaterEqual(
            int((reference_dataframe["ago_family"].str.upper() == "PIWI_RE").sum()),
            1,
        )
        # LONG_A / LONG_B / SHORT must each have >= 1 reference (recall evaluable).
        for required_clade in ("LONG_A", "LONG_B", "SHORT"):
            self.assertIn(required_clade, clades)
        self.assertTrue(
            set(reference_dataframe["reference_label_evidence"].str.upper())
            <= {
                "EXPERIMENTAL",
                "LITERATURE_PHYLOGENETIC",
                "CURATED_COMPUTATIONAL",
                "DATABASE_ANNOTATION",
            }
        )
        self.assertTrue(
            set(reference_dataframe["verification_status"].str.lower())
            <= {"verified", "provisional"}
        )


if __name__ == "__main__":
    unittest.main()
