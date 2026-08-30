from __future__ import annotations

import unittest

import pandas as pd

from src.pago_pipeline.derived_protein_fasta import (
    DerivedFastaRecordOrder,
    build_derived_fasta_selection,
    compute_record_ids_sha256_for_order,
)


def _metadata_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"protein_uid": "10", "gbseq__sequence": "MKAAA"},
            {"protein_uid": "11", "gbseq__sequence": "MKBBB"},
            {"protein_uid": "12", "gbseq__sequence": "MKCCC"},
        ]
    )


class DerivedProteinFastaSelectionTests(unittest.TestCase):
    def test_selection_preserves_requested_order(self) -> None:
        result = build_derived_fasta_selection(
            metadata_dataframe=_metadata_dataframe(),
            selected_protein_uids=["12", "10"],
        )
        self.assertEqual(list(result.resolved_protein_uids), ["12", "10"])
        self.assertEqual(
            result.selected_metadata["protein_uid"].astype(str).tolist(),
            ["12", "10"],
        )
        self.assertEqual(result.requested_uid_count, 2)
        self.assertEqual(result.resolved_uid_count, 2)

    def test_sorted_by_uid_order(self) -> None:
        result = build_derived_fasta_selection(
            metadata_dataframe=_metadata_dataframe(),
            selected_protein_uids=["12", "10", "11"],
            record_order=DerivedFastaRecordOrder.SORTED_BY_UID.value,
        )
        self.assertEqual(list(result.resolved_protein_uids), ["10", "11", "12"])

    def test_record_ids_sha256_matches_helper(self) -> None:
        result = build_derived_fasta_selection(
            metadata_dataframe=_metadata_dataframe(),
            selected_protein_uids=["12", "10"],
        )
        self.assertEqual(
            result.source_record_ids_sha256,
            compute_record_ids_sha256_for_order(
                protein_uids=["12", "10"],
                record_order=DerivedFastaRecordOrder.AS_SELECTED.value,
            ),
        )

    def test_missing_uid_raises_unless_dropped(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "absent from the metadata"):
            build_derived_fasta_selection(
                metadata_dataframe=_metadata_dataframe(),
                selected_protein_uids=["10", "999"],
            )
        result = build_derived_fasta_selection(
            metadata_dataframe=_metadata_dataframe(),
            selected_protein_uids=["10", "999"],
            drop_missing_uids=True,
        )
        self.assertEqual(list(result.resolved_protein_uids), ["10"])
        self.assertEqual(list(result.missing_uids), ["999"])

    def test_duplicate_request_raises(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "duplicates"):
            build_derived_fasta_selection(
                metadata_dataframe=_metadata_dataframe(),
                selected_protein_uids=["10", "10"],
            )

    def test_duplicate_metadata_protein_uid_raises(self) -> None:
        metadata_dataframe = pd.concat(
            [_metadata_dataframe(), _metadata_dataframe().head(1)], ignore_index=True
        )
        with self.assertRaisesRegex(RuntimeError, "must be unique"):
            build_derived_fasta_selection(
                metadata_dataframe=metadata_dataframe,
                selected_protein_uids=["10"],
            )

    def test_empty_sequence_raises(self) -> None:
        metadata_dataframe = pd.DataFrame(
            [{"protein_uid": "10", "gbseq__sequence": "   "}]
        )
        with self.assertRaisesRegex(RuntimeError, "empty amino-acid sequence"):
            build_derived_fasta_selection(
                metadata_dataframe=metadata_dataframe,
                selected_protein_uids=["10"],
            )


if __name__ == "__main__":
    unittest.main()
