from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.pago_qc_filter_snapshot import (
    DEFAULT_FILTERED_DATASET_COUNTS_FILE_NAME,
    EXCLUDED_RECORDS_SEMANTICS,
    SNAPSHOT_FORMAT_VERSION,
    _validate_loaded_pago_qc_filtered_datasets_payload,
    latest_pago_qc_filtered_datasets_is_available,
    resolve_pago_qc_filtered_datasets,
)
from src.pago_pipeline.pago_qc_filter import build_filter_policy_sha256
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
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic


class PagoQcFilteredDatasetsSnapshotTests(unittest.TestCase):
    def _write_source_qc_latest_snapshot(
        self,
        *,
        snapshot_root_directory: Path,
        labelled_records_dataframe: pd.DataFrame,
        snapshot_directory_name: str = "source_snapshot",
    ) -> Path:
        immutable_snapshot_directory = (
            snapshot_root_directory / "snapshots" / snapshot_directory_name
        )
        immutable_snapshot_directory.mkdir(parents=True, exist_ok=True)

        labelled_records_file_path = (
            immutable_snapshot_directory / DEFAULT_LABELLED_RECORDS_FILE_NAME
        )
        labelled_records_dataframe.to_csv(labelled_records_file_path, index=False)

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
        for output_file_name in output_file_names_by_key.values():
            output_file_path = immutable_snapshot_directory / output_file_name
            if output_file_path.exists():
                continue
            output_file_path.write_text("header\n", encoding="utf-8")

        output_files = {
            output_file_key: {
                "file_name": output_file_name,
                "sha256": sha256_of_file(
                    input_file_path=immutable_snapshot_directory / output_file_name
                ),
            }
            for output_file_key, output_file_name in output_file_names_by_key.items()
        }
        write_json_atomic(
            payload={
                "artifact_type": "pago_qc_evidence_inventory",
                "snapshot_format_version": "2.0",
                "immutable_snapshot_relative_path": str(
                    Path("snapshots") / snapshot_directory_name
                ),
                "output_files": output_files,
            },
            output_file_path=immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME,
        )

        latest_directory = snapshot_root_directory / "latest"
        if latest_directory.exists():
            shutil.rmtree(latest_directory)
        shutil.copytree(immutable_snapshot_directory, latest_directory)
        return immutable_snapshot_directory

    def test_resolve_pago_qc_filtered_datasets_persists_outputs_and_manifest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_qc_root_directory = temporary_directory / "pago_qc"
            filtered_root_directory = temporary_directory / "filtered"
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
                ]
            )
            source_snapshot_directory = self._write_source_qc_latest_snapshot(
                snapshot_root_directory=source_qc_root_directory,
                labelled_records_dataframe=labelled_records_dataframe,
            )

            filtered_payload = resolve_pago_qc_filtered_datasets(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                source_qc_snapshot_root_directory=source_qc_root_directory,
                snapshot_root_directory=filtered_root_directory,
            )

            self.assertEqual(
                filtered_payload["manifest"]["artifact_type"],
                "pago_qc_filtered_datasets",
            )
            self.assertEqual(
                filtered_payload["manifest"]["snapshot_format_version"],
                SNAPSHOT_FORMAT_VERSION,
            )
            self.assertEqual(
                filtered_payload["manifest"]["filter_policy_sha256"],
                build_filter_policy_sha256(),
            )
            self.assertEqual(
                filtered_payload["manifest"]["source_qc_artifact_type"],
                "pago_qc_evidence_inventory",
            )
            self.assertEqual(
                filtered_payload["manifest"]["source_qc_snapshot_directory"],
                str(source_snapshot_directory),
            )
            self.assertEqual(
                filtered_payload["manifest"]["source_qc_snapshot_format_version"],
                "2.0",
            )
            self.assertEqual(
                filtered_payload["manifest"]["excluded_records_semantics"],
                EXCLUDED_RECORDS_SEMANTICS,
            )
            self.assertTrue(
                filtered_payload["filtered_dataset_counts_file_path"]
                .name
                .endswith(DEFAULT_FILTERED_DATASET_COUNTS_FILE_NAME)
            )

            self.assertEqual(
                filtered_payload["classic_pago_high_precision_records"][
                    "protein_uid"
                ].tolist(),
                ["1"],
            )
            self.assertEqual(
                filtered_payload["classic_pago_review_records"][
                    "protein_uid"
                ].tolist(),
                ["2"],
            )
            self.assertEqual(
                filtered_payload["piwi_re_records"]["protein_uid"].tolist(),
                ["3"],
            )
            self.assertEqual(
                filtered_payload["excluded_records"]["protein_uid"].tolist(),
                ["4"],
            )

    def test_filtered_dataset_validation_rejects_snapshot_version_mismatch(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            records_file_path = temporary_directory / "records.csv"
            counts_file_path = temporary_directory / DEFAULT_FILTERED_DATASET_COUNTS_FILE_NAME
            records_file_path.write_text("protein_uid\n1\n", encoding="utf-8")
            counts_file_path.write_text("dataset_name,count\nx,1\n", encoding="utf-8")
            manifest_payload = {
                "artifact_type": "pago_qc_filtered_datasets",
                "snapshot_format_version": "0.9",
                "filter_policy_sha256": build_filter_policy_sha256(),
                "output_files": {
                    "classic_pago_high_precision_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "classic_pago_review_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "piwi_re_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "excluded_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "filtered_dataset_counts_file": {
                        "file_name": counts_file_path.name,
                        "sha256": sha256_of_file(input_file_path=counts_file_path),
                    },
                },
            }

            with self.assertRaisesRegex(RuntimeError, "snapshot_format_version"):
                _validate_loaded_pago_qc_filtered_datasets_payload(
                    snapshot_directory=temporary_directory,
                    manifest_payload=manifest_payload,
                )

    def test_filtered_dataset_validation_rejects_filter_policy_mismatch(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            records_file_path = temporary_directory / "records.csv"
            counts_file_path = temporary_directory / DEFAULT_FILTERED_DATASET_COUNTS_FILE_NAME
            records_file_path.write_text("protein_uid\n1\n", encoding="utf-8")
            counts_file_path.write_text("dataset_name,count\nx,1\n", encoding="utf-8")
            manifest_payload = {
                "artifact_type": "pago_qc_filtered_datasets",
                "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
                "filter_policy_sha256": "not-current-policy",
                "output_files": {
                    "classic_pago_high_precision_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "classic_pago_review_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "piwi_re_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "excluded_records_file": {
                        "file_name": records_file_path.name,
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    },
                    "filtered_dataset_counts_file": {
                        "file_name": counts_file_path.name,
                        "sha256": sha256_of_file(input_file_path=counts_file_path),
                    },
                },
            }

            with self.assertRaisesRegex(RuntimeError, "filter_policy_sha256"):
                _validate_loaded_pago_qc_filtered_datasets_payload(
                    snapshot_directory=temporary_directory,
                    manifest_payload=manifest_payload,
                )

    def test_latest_pago_qc_filtered_datasets_rejects_changed_source_qc(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_qc_root_directory = temporary_directory / "pago_qc"
            filtered_root_directory = temporary_directory / "filtered"
            first_labelled_records_dataframe = pd.DataFrame(
                [
                    {
                        "protein_uid": "1",
                        "primary_label": "classic_piwi_candidate",
                        "qc_decision": "include",
                    },
                ]
            )
            second_labelled_records_dataframe = pd.DataFrame(
                [
                    {
                        "protein_uid": "2",
                        "primary_label": "piwi_re_candidate",
                        "qc_decision": "separate_dataset",
                    },
                ]
            )
            self._write_source_qc_latest_snapshot(
                snapshot_root_directory=source_qc_root_directory,
                labelled_records_dataframe=first_labelled_records_dataframe,
                snapshot_directory_name="source_snapshot_1",
            )
            resolve_pago_qc_filtered_datasets(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                source_qc_snapshot_root_directory=source_qc_root_directory,
                snapshot_root_directory=filtered_root_directory,
            )

            self.assertTrue(
                latest_pago_qc_filtered_datasets_is_available(
                    snapshot_root_directory=filtered_root_directory,
                    source_qc_snapshot_root_directory=source_qc_root_directory,
                )
            )

            self._write_source_qc_latest_snapshot(
                snapshot_root_directory=source_qc_root_directory,
                labelled_records_dataframe=second_labelled_records_dataframe,
                snapshot_directory_name="source_snapshot_2",
            )

            self.assertFalse(
                latest_pago_qc_filtered_datasets_is_available(
                    snapshot_root_directory=filtered_root_directory,
                    source_qc_snapshot_root_directory=source_qc_root_directory,
                )
            )
