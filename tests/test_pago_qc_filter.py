from __future__ import annotations

import unittest

import pandas as pd

from src.pago_pipeline.pago_qc_filter import (
    build_filtered_dataset_counts_dataframe,
    build_pago_qc_filtered_datasets,
)


class PagoQcFilteredDatasetsTests(unittest.TestCase):
    def test_build_pago_qc_filtered_datasets_partitions_labelled_records(
        self,
    ) -> None:
        labelled_records_dataframe = pd.DataFrame(
            [
                {
                    "protein_uid": "1",
                    "primary_label": "classic_piwi_candidate",
                    "qc_decision": "include",
                },
                {
                    "protein_uid": "2",
                    "primary_label": "classic_piwi_candidate",
                    "qc_decision": "review",
                },
                {
                    "protein_uid": "3",
                    "primary_label": "piwi_re_candidate",
                    "qc_decision": "separate_dataset",
                },
                {
                    "protein_uid": "4",
                    "primary_label": (
                        "methyltransferase_noise_or_unresolved_conflict"
                    ),
                    "qc_decision": "exclude",
                },
                {
                    "protein_uid": "5",
                    "primary_label": "low_evidence_or_unrelated",
                    "qc_decision": "exclude",
                },
            ]
        )

        filtered_dataset_by_name = build_pago_qc_filtered_datasets(
            labelled_records_dataframe=labelled_records_dataframe,
        )
        filtered_dataset_counts_dataframe = build_filtered_dataset_counts_dataframe(
            filtered_dataset_by_name=filtered_dataset_by_name,
            total_record_count=len(labelled_records_dataframe),
        )

        self.assertEqual(
            filtered_dataset_by_name["classic_pago_high_precision"][
                "protein_uid"
            ].tolist(),
            ["1"],
        )
        self.assertEqual(
            filtered_dataset_by_name["classic_pago_review"][
                "protein_uid"
            ].tolist(),
            ["2"],
        )
        self.assertEqual(
            filtered_dataset_by_name["piwi_re"]["protein_uid"].tolist(),
            ["3"],
        )
        self.assertEqual(
            filtered_dataset_by_name["excluded"]["protein_uid"].tolist(),
            ["4", "5"],
        )

        count_by_dataset_name = dict(
            zip(
                filtered_dataset_counts_dataframe["dataset_name"],
                filtered_dataset_counts_dataframe["count"],
                strict=True,
            )
        )
        self.assertEqual(
            count_by_dataset_name,
            {
                "classic_pago_high_precision": 1,
                "classic_pago_review": 1,
                "piwi_re": 1,
                "excluded": 2,
            },
        )
        excluded_description = filtered_dataset_counts_dataframe.loc[
            filtered_dataset_counts_dataframe["dataset_name"] == "excluded",
            "description",
        ].iloc[0]
        self.assertIn(
            "not a biologically validated negative class",
            excluded_description,
        )

    def test_build_pago_qc_filtered_datasets_reports_missing_input_columns(
        self,
    ) -> None:
        with self.assertRaisesRegex(
            RuntimeError,
            "missing required filter input columns",
        ):
            build_pago_qc_filtered_datasets(
                labelled_records_dataframe=pd.DataFrame(),
            )

    def test_build_pago_qc_filtered_datasets_rejects_unassigned_records(
        self,
    ) -> None:
        labelled_records_dataframe = pd.DataFrame(
            [
                {
                    "protein_uid": "1",
                    "primary_label": "piwi_re_candidate",
                    "qc_decision": "review",
                },
            ]
        )

        with self.assertRaisesRegex(
            RuntimeError,
            "Unassigned records: 1",
        ):
            build_pago_qc_filtered_datasets(
                labelled_records_dataframe=labelled_records_dataframe,
            )

    def test_build_filtered_dataset_counts_dataframe_requires_complete_counts(
        self,
    ) -> None:
        filtered_dataset_by_name = {
            "classic_pago_high_precision": pd.DataFrame({"protein_uid": ["1"]}),
            "classic_pago_review": pd.DataFrame({"protein_uid": []}),
            "piwi_re": pd.DataFrame({"protein_uid": []}),
            "excluded": pd.DataFrame({"protein_uid": []}),
        }

        with self.assertRaisesRegex(
            RuntimeError,
            "must sum to total_record_count",
        ):
            build_filtered_dataset_counts_dataframe(
                filtered_dataset_by_name=filtered_dataset_by_name,
                total_record_count=2,
            )
