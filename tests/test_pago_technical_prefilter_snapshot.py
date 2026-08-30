from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.pago_technical_prefilter import (
    build_technical_prefilter_policy_sha256,
)
from src.pago_pipeline.pago_technical_prefilter_snapshot import (
    ARTIFACT_TYPE,
    DEFAULT_MANIFEST_FILE_NAME,
    SNAPSHOT_FORMAT_VERSION,
    _validate_loaded_pago_technical_prefilter_payload,
    latest_pago_technical_prefilter_snapshot_is_available,
    load_pago_technical_prefilter_snapshot_by_directory,
    resolve_pago_technical_prefilter_snapshot,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic


def _write_metadata_latest_snapshot(
    *,
    snapshot_root_directory: Path,
    metadata_dataframe: pd.DataFrame,
    snapshot_directory_name: str = "metadata_snapshot_1",
    search_query: str = "(PIWI[All Fields] OR Argonaute[All Fields])",
) -> Path:
    immutable_snapshot_directory = (
        snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=True)

    csv_file_path = immutable_snapshot_directory / "protein_metadata.csv"
    metadata_dataframe.to_csv(csv_file_path, index=False)
    qc_report_file_path = immutable_snapshot_directory / "qc_report.json"
    write_json_atomic(payload={"checks": {}}, output_file_path=qc_report_file_path)

    manifest_file_path = immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    write_json_atomic(
        payload={
            "artifact_type": "ncbi_protein_metadata_snapshot",
            "snapshot_format_version": "1.0",
            "csv_file_name": "protein_metadata.csv",
            "csv_file_sha256": sha256_of_file(input_file_path=csv_file_path),
            "qc_report_file_name": "qc_report.json",
            "immutable_snapshot_directory_name": snapshot_directory_name,
            "immutable_snapshot_relative_path": str(
                Path("snapshots") / snapshot_directory_name
            ),
            "row_count": int(len(metadata_dataframe)),
            "column_count": int(len(metadata_dataframe.columns)),
            "search_query": search_query,
            "translated_query": search_query,
            "source_xml_snapshot_relative_path": "snapshots/xml_1",
        },
        output_file_path=manifest_file_path,
    )

    latest_directory = snapshot_root_directory / "latest"
    if latest_directory.exists():
        shutil.rmtree(latest_directory)
    shutil.copytree(immutable_snapshot_directory, latest_directory)
    return immutable_snapshot_directory


def _sample_metadata_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"protein_uid": "10", "gbseq__sequence": "MKTAYIAK", "gbseq__length": 8},
            {"protein_uid": "11", "gbseq__sequence": "", "gbseq__length": 0},
            {"protein_uid": "12", "gbseq__sequence": "MKT@1", "gbseq__length": 5},
            {"protein_uid": "13", "gbseq__sequence": "M" * 900, "gbseq__length": 900},
        ]
    )


class PagoTechnicalPrefilterSnapshotTests(unittest.TestCase):
    def test_resolve_persists_outputs_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            metadata_root = temporary_directory / "metadata"
            prefilter_root = temporary_directory / "prefilter"
            source_snapshot_directory = _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                metadata_dataframe=_sample_metadata_dataframe(),
            )

            payload = resolve_pago_technical_prefilter_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                snapshot_root_directory=prefilter_root,
                source_metadata_snapshot_root_directory=metadata_root,
            )

            manifest = payload["manifest"]
            self.assertEqual(manifest["artifact_type"], ARTIFACT_TYPE)
            self.assertEqual(
                manifest["snapshot_format_version"], SNAPSHOT_FORMAT_VERSION
            )
            self.assertEqual(
                manifest["technical_prefilter_policy_sha256"],
                build_technical_prefilter_policy_sha256(),
            )
            self.assertEqual(manifest["input_record_count"], 4)
            self.assertEqual(manifest["retained_record_count"], 2)
            self.assertEqual(manifest["excluded_record_count"], 2)
            self.assertEqual(
                manifest["source_metadata_snapshot_directory_name"],
                source_snapshot_directory.name,
            )
            self.assertEqual(
                manifest["search_query"],
                "(PIWI[All Fields] OR Argonaute[All Fields])",
            )
            self.assertEqual(
                payload["retained_records"]["protein_uid"].astype(str).tolist(),
                ["10", "13"],
            )
            retained_uids_file = Path(payload["retained_protein_uids_file_path"])
            self.assertEqual(
                retained_uids_file.read_text(encoding="utf-8").split(),
                ["10", "13"],
            )

    def test_reuse_latest_or_create_reuses_frozen_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            metadata_root = temporary_directory / "metadata"
            prefilter_root = temporary_directory / "prefilter"
            _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                metadata_dataframe=_sample_metadata_dataframe(),
            )
            first_payload = resolve_pago_technical_prefilter_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                snapshot_root_directory=prefilter_root,
                source_metadata_snapshot_root_directory=metadata_root,
            )
            with mock.patch(
                "src.pago_pipeline.pago_technical_prefilter_snapshot"
                ".save_pago_technical_prefilter_snapshot"
            ) as save_mock:
                second_payload = resolve_pago_technical_prefilter_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=prefilter_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                )
                save_mock.assert_not_called()
            self.assertEqual(
                first_payload["manifest"]["immutable_snapshot_directory_name"],
                second_payload["manifest"]["immutable_snapshot_directory_name"],
            )

    def test_latest_is_unavailable_when_metadata_snapshot_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            metadata_root = temporary_directory / "metadata"
            prefilter_root = temporary_directory / "prefilter"
            _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                metadata_dataframe=_sample_metadata_dataframe(),
                snapshot_directory_name="metadata_snapshot_1",
            )
            resolve_pago_technical_prefilter_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                snapshot_root_directory=prefilter_root,
                source_metadata_snapshot_root_directory=metadata_root,
            )
            self.assertTrue(
                latest_pago_technical_prefilter_snapshot_is_available(
                    snapshot_root_directory=prefilter_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                )
            )

            changed_metadata_dataframe = pd.concat(
                [
                    _sample_metadata_dataframe(),
                    pd.DataFrame(
                        [
                            {
                                "protein_uid": "99",
                                "gbseq__sequence": "MKT",
                                "gbseq__length": 3,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                metadata_dataframe=changed_metadata_dataframe,
                snapshot_directory_name="metadata_snapshot_2",
            )
            self.assertFalse(
                latest_pago_technical_prefilter_snapshot_is_available(
                    snapshot_root_directory=prefilter_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                )
            )

    def test_validation_rejects_artifact_type_and_hash_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            records_file_path = temporary_directory / "retained_records.csv"
            records_file_path.write_text("protein_uid\n1\n", encoding="utf-8")
            base_manifest = {
                "artifact_type": ARTIFACT_TYPE,
                "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
                "output_files": {
                    "retained_records_file": {
                        "file_name": "retained_records.csv",
                        "sha256": sha256_of_file(input_file_path=records_file_path),
                    }
                },
            }

            with self.assertRaisesRegex(RuntimeError, "artifact_type"):
                _validate_loaded_pago_technical_prefilter_payload(
                    snapshot_directory=temporary_directory,
                    manifest_payload={**base_manifest, "artifact_type": "wrong"},
                )

            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                _validate_loaded_pago_technical_prefilter_payload(
                    snapshot_directory=temporary_directory,
                    manifest_payload={
                        **base_manifest,
                        "output_files": {
                            "retained_records_file": {
                                "file_name": "retained_records.csv",
                                "sha256": "0" * 64,
                            }
                        },
                    },
                )


if __name__ == "__main__":
    unittest.main()
