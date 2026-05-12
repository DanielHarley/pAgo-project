from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.pago_qc_snapshot import (
    DEFAULT_EVIDENCE_COUNTS_FILE_NAME,
    DEFAULT_EVIDENCE_FLAGS_FILE_NAME,
    DEFAULT_FILTER_DECISION_COUNTS_FILE_NAME,
    DEFAULT_LABEL_COUNTS_FILE_NAME,
    DEFAULT_LABELLED_RECORDS_FILE_NAME,
    DEFAULT_MANIFEST_FILE_NAME,
    DEFAULT_SUSPICIOUS_TERMS_FILE_NAME,
    DEFAULT_TOP_PRODUCTS_FILE_NAME,
    DEFAULT_TOP_REGION_NAMES_FILE_NAME,
    SNAPSHOT_FORMAT_VERSION,
    _validate_loaded_pago_qc_evidence_inventory_payload,
    latest_pago_qc_evidence_inventory_is_available,
    resolve_pago_qc_evidence_inventory,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic


class PagoQcSnapshotValidationTests(unittest.TestCase):
    def _write_minimal_latest_snapshot(
        self,
        *,
        latest_directory: Path,
        metadata_csv_file_path: Path,
    ) -> None:
        latest_directory.mkdir(parents=True, exist_ok=True)
        output_file_names_by_key = {
            "evidence_flags_file": DEFAULT_EVIDENCE_FLAGS_FILE_NAME,
            "evidence_counts_file": DEFAULT_EVIDENCE_COUNTS_FILE_NAME,
            "labelled_records_file": DEFAULT_LABELLED_RECORDS_FILE_NAME,
            "label_counts_file": DEFAULT_LABEL_COUNTS_FILE_NAME,
            "filter_decision_counts_file": DEFAULT_FILTER_DECISION_COUNTS_FILE_NAME,
            "top_region_names_file": DEFAULT_TOP_REGION_NAMES_FILE_NAME,
            "top_products_file": DEFAULT_TOP_PRODUCTS_FILE_NAME,
            "suspicious_terms_file": DEFAULT_SUSPICIOUS_TERMS_FILE_NAME,
        }
        output_files = {}
        for output_file_key, output_file_name in output_file_names_by_key.items():
            output_file_path = latest_directory / output_file_name
            output_file_path.write_text("header\n", encoding="utf-8")
            output_files[output_file_key] = {
                "file_name": output_file_name,
                "sha256": sha256_of_file(input_file_path=output_file_path),
            }

        write_json_atomic(
            payload={
                "artifact_type": "pago_qc_evidence_inventory",
                "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
                "metadata_csv_sha256": sha256_of_file(
                    input_file_path=metadata_csv_file_path
                ),
                "output_files": output_files,
            },
            output_file_path=latest_directory / DEFAULT_MANIFEST_FILE_NAME,
        )

    def test_pago_qc_inventory_validation_requires_output_file_sha256(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            output_file_names_by_key = {
                "evidence_flags_file": DEFAULT_EVIDENCE_FLAGS_FILE_NAME,
                "evidence_counts_file": DEFAULT_EVIDENCE_COUNTS_FILE_NAME,
                "labelled_records_file": DEFAULT_LABELLED_RECORDS_FILE_NAME,
                "label_counts_file": DEFAULT_LABEL_COUNTS_FILE_NAME,
                "filter_decision_counts_file": (
                    DEFAULT_FILTER_DECISION_COUNTS_FILE_NAME
                ),
                "top_region_names_file": DEFAULT_TOP_REGION_NAMES_FILE_NAME,
                "top_products_file": DEFAULT_TOP_PRODUCTS_FILE_NAME,
                "suspicious_terms_file": DEFAULT_SUSPICIOUS_TERMS_FILE_NAME,
            }
            for output_file_name in output_file_names_by_key.values():
                (temporary_directory / output_file_name).write_text(
                    "header\n",
                    encoding="utf-8",
                )

            manifest_payload = {
                "artifact_type": "pago_qc_evidence_inventory",
                "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
                "output_files": {
                    output_file_key: {
                        "file_name": output_file_name,
                        "sha256": sha256_of_file(
                            input_file_path=temporary_directory / output_file_name
                        ),
                    }
                    for output_file_key, output_file_name in (
                        output_file_names_by_key.items()
                    )
                },
            }
            del manifest_payload["output_files"]["evidence_counts_file"]["sha256"]

            with self.assertRaisesRegex(RuntimeError, "non-empty sha256"):
                _validate_loaded_pago_qc_evidence_inventory_payload(
                    snapshot_directory=temporary_directory,
                    manifest_payload=manifest_payload,
                )

    def test_pago_qc_inventory_validation_rejects_snapshot_version_mismatch(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            output_file_names_by_key = {
                "evidence_flags_file": DEFAULT_EVIDENCE_FLAGS_FILE_NAME,
                "evidence_counts_file": DEFAULT_EVIDENCE_COUNTS_FILE_NAME,
                "labelled_records_file": DEFAULT_LABELLED_RECORDS_FILE_NAME,
                "label_counts_file": DEFAULT_LABEL_COUNTS_FILE_NAME,
                "filter_decision_counts_file": (
                    DEFAULT_FILTER_DECISION_COUNTS_FILE_NAME
                ),
                "top_region_names_file": DEFAULT_TOP_REGION_NAMES_FILE_NAME,
                "top_products_file": DEFAULT_TOP_PRODUCTS_FILE_NAME,
                "suspicious_terms_file": DEFAULT_SUSPICIOUS_TERMS_FILE_NAME,
            }
            for output_file_name in output_file_names_by_key.values():
                (temporary_directory / output_file_name).write_text(
                    "header\n",
                    encoding="utf-8",
                )

            manifest_payload = {
                "artifact_type": "pago_qc_evidence_inventory",
                "snapshot_format_version": "1.0",
                "output_files": {
                    output_file_key: {
                        "file_name": output_file_name,
                        "sha256": sha256_of_file(
                            input_file_path=temporary_directory / output_file_name
                        ),
                    }
                    for output_file_key, output_file_name in (
                        output_file_names_by_key.items()
                    )
                },
            }

            with self.assertRaisesRegex(RuntimeError, "snapshot_format_version"):
                _validate_loaded_pago_qc_evidence_inventory_payload(
                    snapshot_directory=temporary_directory,
                    manifest_payload=manifest_payload,
                )

    def test_latest_pago_qc_inventory_availability_checks_metadata_csv_hash(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "pago_qc"
            metadata_csv_file_path = temporary_directory / "metadata.csv"
            changed_metadata_csv_file_path = temporary_directory / "changed.csv"
            metadata_csv_file_path.write_text(
                "protein_uid,gbseq__length\n1,750\n",
                encoding="utf-8",
            )
            changed_metadata_csv_file_path.write_text(
                "protein_uid,gbseq__length\n2,750\n",
                encoding="utf-8",
            )

            self._write_minimal_latest_snapshot(
                latest_directory=snapshot_root_directory / "latest",
                metadata_csv_file_path=metadata_csv_file_path,
            )

            self.assertTrue(
                latest_pago_qc_evidence_inventory_is_available(
                    snapshot_root_directory=snapshot_root_directory,
                    metadata_csv_file_path=metadata_csv_file_path,
                )
            )
            self.assertFalse(
                latest_pago_qc_evidence_inventory_is_available(
                    snapshot_root_directory=snapshot_root_directory,
                    metadata_csv_file_path=changed_metadata_csv_file_path,
                )
            )

    def test_resolve_pago_qc_inventory_recreates_latest_for_changed_metadata(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "pago_qc"
            old_metadata_csv_file_path = temporary_directory / "old_metadata.csv"
            new_metadata_csv_file_path = temporary_directory / "new_metadata.csv"
            old_metadata_csv_file_path.write_text(
                "protein_uid,gbseq__length\n1,750\n",
                encoding="utf-8",
            )
            new_metadata_csv_file_path.write_text(
                "protein_uid,gbseq__length\n2,982\n",
                encoding="utf-8",
            )
            self._write_minimal_latest_snapshot(
                latest_directory=snapshot_root_directory / "latest",
                metadata_csv_file_path=old_metadata_csv_file_path,
            )

            inventory_payload = resolve_pago_qc_evidence_inventory(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                metadata_csv_file_path=new_metadata_csv_file_path,
                snapshot_root_directory=snapshot_root_directory,
            )

            self.assertEqual(
                inventory_payload["manifest"]["metadata_csv_sha256"],
                sha256_of_file(input_file_path=new_metadata_csv_file_path),
            )
            pd.testing.assert_series_equal(
                inventory_payload["evidence_flags"]["protein_uid"],
                pd.Series(["2"], name="protein_uid", dtype="string"),
            )
            self.assertIn("labelled_records", inventory_payload)
            self.assertIn("label_counts", inventory_payload)
            self.assertIn("filter_decision_counts", inventory_payload)
            self.assertEqual(
                inventory_payload["manifest"]["snapshot_format_version"],
                "2.0",
            )
            self.assertEqual(
                inventory_payload["manifest"]["label_columns"],
                ["primary_label", "qc_decision", "confidence_score", "rationale"],
            )
            self.assertEqual(
                inventory_payload["manifest"]["length_bin_column"],
                "length_bin",
            )
