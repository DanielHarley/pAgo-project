from __future__ import annotations

import unittest

import pandas as pd

from src.pago_pipeline.pago_technical_prefilter import (
    LENGTH_WARNING_COLUMN,
    PagoTechnicalPrefilterDecision,
    build_pago_technical_prefilter_partition,
    build_technical_prefilter_counts_dataframe,
    build_technical_prefilter_policy_sha256,
)


def _metadata_dataframe(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows)


class PagoTechnicalPrefilterLogicTests(unittest.TestCase):
    def test_retains_records_regardless_of_annotation_text_or_length(self) -> None:
        metadata_dataframe = _metadata_dataframe(
            [
                {
                    "protein_uid": "1",
                    "gbseq__sequence": "MKTAYIAKQR",
                    "gbseq__length": 10,
                    "gbseq__definition": "SAM-dependent methyltransferase",
                },
                {
                    "protein_uid": "2",
                    "gbseq__sequence": "M" * 5000,
                    "gbseq__length": 5000,
                    "gbseq__definition": "hypothetical protein",
                },
                {
                    "protein_uid": "3",
                    "gbseq__sequence": "MPIWIARGONAUTE",
                    "gbseq__length": 14,
                    "gbseq__definition": "transposase",
                },
            ]
        )

        result = build_pago_technical_prefilter_partition(
            metadata_dataframe=metadata_dataframe
        )

        self.assertEqual(
            result.retained_records["protein_uid"].tolist(), ["1", "2", "3"]
        )
        self.assertEqual(len(result.excluded_records), 0)
        # Length outside [200, 2000] only warns, never excludes.
        self.assertEqual(
            result.retained_records[LENGTH_WARNING_COLUMN].tolist(),
            [True, True, True],
        )

    def test_excludes_only_technical_problems(self) -> None:
        metadata_dataframe = _metadata_dataframe(
            [
                {"protein_uid": "1", "gbseq__sequence": "MKTAYIAK", "gbseq__length": 8},
                {"protein_uid": "", "gbseq__sequence": "MKTAYIAK", "gbseq__length": 8},
                {"protein_uid": "1", "gbseq__sequence": "MKTAYIAK", "gbseq__length": 8},
                {"protein_uid": "3", "gbseq__sequence": "   ", "gbseq__length": 0},
                {"protein_uid": "4", "gbseq__sequence": "MKT@YIAK1", "gbseq__length": 9},
            ]
        )

        result = build_pago_technical_prefilter_partition(
            metadata_dataframe=metadata_dataframe
        )

        decision_by_uid = dict(
            zip(
                result.excluded_records["protein_uid"].tolist(),
                result.excluded_records["technical_prefilter_decision"].tolist(),
            )
        )
        self.assertEqual(result.retained_records["protein_uid"].tolist(), ["1"])
        self.assertEqual(
            decision_by_uid.get(""),
            PagoTechnicalPrefilterDecision.DROP_UNPROCESSABLE_RECORD.value,
        )
        self.assertEqual(
            result.excluded_records[
                result.excluded_records["technical_prefilter_decision"]
                == PagoTechnicalPrefilterDecision.DROP_TECHNICAL_DUPLICATE.value
            ]["protein_uid"].tolist(),
            ["1"],
        )
        self.assertEqual(
            decision_by_uid.get("3"),
            PagoTechnicalPrefilterDecision.DROP_MISSING_SEQUENCE.value,
        )
        self.assertEqual(
            decision_by_uid.get("4"),
            PagoTechnicalPrefilterDecision.DROP_INVALID_SEQUENCE_CHARACTERS.value,
        )

    def test_tolerated_ambiguous_and_symbol_residues_are_kept(self) -> None:
        metadata_dataframe = _metadata_dataframe(
            [
                {"protein_uid": "1", "gbseq__sequence": "MKTXBZUO*-", "gbseq__length": 9},
            ]
        )
        result = build_pago_technical_prefilter_partition(
            metadata_dataframe=metadata_dataframe
        )
        self.assertEqual(result.retained_records["protein_uid"].tolist(), ["1"])

    def test_partition_is_total_and_disjoint(self) -> None:
        metadata_dataframe = _metadata_dataframe(
            [
                {"protein_uid": str(i), "gbseq__sequence": "MKT", "gbseq__length": 3}
                for i in range(7)
            ]
            + [{"protein_uid": "", "gbseq__sequence": "MKT", "gbseq__length": 3}]
        )
        result = build_pago_technical_prefilter_partition(
            metadata_dataframe=metadata_dataframe
        )
        self.assertEqual(
            len(result.retained_records) + len(result.excluded_records),
            result.input_record_count,
        )
        self.assertEqual(
            sum(result.counts_by_decision.values()), result.input_record_count
        )

    def test_missing_required_columns_raises(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "required technical prefilter"):
            build_pago_technical_prefilter_partition(
                metadata_dataframe=pd.DataFrame([{"protein_uid": "1"}])
            )

    def test_policy_sha256_is_stable_and_parameter_sensitive(self) -> None:
        self.assertEqual(
            build_technical_prefilter_policy_sha256(),
            build_technical_prefilter_policy_sha256(),
        )
        self.assertNotEqual(
            build_technical_prefilter_policy_sha256(),
            build_technical_prefilter_policy_sha256(length_warning_max=3000),
        )

    def test_counts_dataframe_sums_to_input(self) -> None:
        metadata_dataframe = _metadata_dataframe(
            [
                {"protein_uid": "1", "gbseq__sequence": "MKT", "gbseq__length": 3},
                {"protein_uid": "2", "gbseq__sequence": "", "gbseq__length": 0},
            ]
        )
        result = build_pago_technical_prefilter_partition(
            metadata_dataframe=metadata_dataframe
        )
        counts_dataframe = build_technical_prefilter_counts_dataframe(result=result)
        self.assertEqual(int(counts_dataframe["count"].sum()), 2)


if __name__ == "__main__":
    unittest.main()
