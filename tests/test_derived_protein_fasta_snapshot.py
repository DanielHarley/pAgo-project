from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from src.pago_pipeline.derived_protein_fasta import (
    parse_protein_uids_from_fasta_deflines,
)
from src.pago_pipeline.derived_protein_fasta_snapshot import (
    ARTIFACT_TYPE,
    DEFAULT_MANIFEST_FILE_NAME,
    _validate_loaded_derived_protein_fasta_payload,
    latest_derived_protein_fasta_snapshot_is_available,
    resolve_derived_protein_fasta_snapshot,
)
from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic

_SELECTION_ARTIFACT_TYPE = "pago_technical_prefilter"
_SELECTION_UID_LIST_FILE_NAME = "retained_protein_uids.txt"


def _metadata_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "protein_uid": str(10 + i),
                "gbseq__accession_version": f"WP_{i:09d}.1",
                "gbseq__length": 5,
                "gbseq__organism": "Escherichia coli",
                "gbseq__definition": "hypothetical protein",
                "gbseq__sequence": "MK" + "ACD"[i % 3] + "EF",
            }
            for i in range(5)
        ]
    )


def _write_metadata_latest_snapshot(*, snapshot_root_directory: Path) -> Path:
    directory = snapshot_root_directory / "snapshots" / "metadata_1"
    directory.mkdir(parents=True, exist_ok=True)
    csv_file_path = directory / "protein_metadata.csv"
    _metadata_dataframe().to_csv(csv_file_path, index=False)
    write_json_atomic(payload={"checks": {}}, output_file_path=directory / "qc_report.json")
    write_json_atomic(
        payload={
            "artifact_type": "ncbi_protein_metadata_snapshot",
            "snapshot_format_version": "1.0",
            "csv_file_name": "protein_metadata.csv",
            "csv_file_sha256": sha256_of_file(input_file_path=csv_file_path),
            "qc_report_file_name": "qc_report.json",
            "immutable_snapshot_directory_name": "metadata_1",
            "immutable_snapshot_relative_path": "snapshots/metadata_1",
            "row_count": 5,
            "search_query": "(PIWI[All Fields] OR Argonaute[All Fields])",
            "translated_query": "translated",
            "source_xml_snapshot_relative_path": "snapshots/xml_1",
        },
        output_file_path=directory / DEFAULT_MANIFEST_FILE_NAME,
    )
    latest = snapshot_root_directory / "latest"
    if latest.exists():
        shutil.rmtree(latest)
    shutil.copytree(directory, latest)
    return directory


def _write_selection_latest_snapshot(
    *,
    snapshot_root_directory: Path,
    retained_uids: list[str],
    snapshot_directory_name: str = "selection_1",
) -> Path:
    directory = snapshot_root_directory / "snapshots" / snapshot_directory_name
    directory.mkdir(parents=True, exist_ok=True)
    uid_file_path = directory / _SELECTION_UID_LIST_FILE_NAME
    uid_file_path.write_text("\n".join(retained_uids) + "\n", encoding="utf-8")
    write_json_atomic(
        payload={
            "artifact_type": _SELECTION_ARTIFACT_TYPE,
            "snapshot_format_version": "1.0",
            "immutable_snapshot_directory_name": snapshot_directory_name,
            "immutable_snapshot_relative_path": (
                f"snapshots/{snapshot_directory_name}"
            ),
            "technical_prefilter_policy_sha256": "policy-sha",
        },
        output_file_path=directory / DEFAULT_MANIFEST_FILE_NAME,
    )
    latest = snapshot_root_directory / "latest"
    if latest.exists():
        shutil.rmtree(latest)
    shutil.copytree(directory, latest)
    return directory


class DerivedProteinFastaSnapshotTests(unittest.TestCase):
    def _resolve(self, temporary_directory: Path, retained_uids: list[str]):
        metadata_root = temporary_directory / "metadata"
        selection_root = temporary_directory / "selection"
        derived_root = temporary_directory / "derived"
        _write_metadata_latest_snapshot(snapshot_root_directory=metadata_root)
        _write_selection_latest_snapshot(
            snapshot_root_directory=selection_root, retained_uids=retained_uids
        )
        payload = resolve_derived_protein_fasta_snapshot(
            snapshot_mode=SnapshotMode.reuse_latest_or_create,
            snapshot_root_directory=derived_root,
            source_metadata_snapshot_root_directory=metadata_root,
            source_selection_snapshot_root_directory=selection_root,
            selection_artifact_type=_SELECTION_ARTIFACT_TYPE,
            selection_uid_list_file_name=_SELECTION_UID_LIST_FILE_NAME,
            record_selection_rule="pago_technical_prefilter.retained",
            record_selection_config_sha256="policy-sha",
            dataset_kind="annotation_enriched_proteome",
        )
        return payload, metadata_root, selection_root, derived_root

    def test_resolve_produces_honest_provenance_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            payload, _, _, _ = self._resolve(temporary_directory, ["12", "10"])

            manifest = payload["manifest"]
            self.assertEqual(manifest["artifact_type"], ARTIFACT_TYPE)
            self.assertNotEqual(manifest["artifact_type"], "ncbi_protein_fasta_snapshot")
            self.assertEqual(
                manifest["derived_from_artifact_type"], _SELECTION_ARTIFACT_TYPE
            )
            self.assertEqual(
                manifest["record_selection_rule"],
                "pago_technical_prefilter.retained",
            )
            self.assertEqual(manifest["dataset_kind"], "annotation_enriched_proteome")
            self.assertEqual(manifest["requested_uid_count"], 2)
            self.assertEqual(manifest["resolved_uid_count"], 2)
            self.assertEqual(manifest["fasta_record_count"], 2)
            self.assertEqual(manifest["record_order"], "as_selected")
            self.assertIn("source_record_ids_sha256", manifest)
            self.assertIn("fasta_file_sha256", manifest)
            self.assertEqual(
                manifest["fasta_file_name"], "protein_sequences.fasta"
            )
            self.assertEqual(
                manifest["search_query"],
                "(PIWI[All Fields] OR Argonaute[All Fields])",
            )

            fasta_uids = parse_protein_uids_from_fasta_deflines(
                fasta_file_path=str(payload["fasta_file_path"])
            )
            self.assertEqual(fasta_uids, ["12", "10"])

    def test_validation_detects_reordered_or_tampered_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            payload, _, _, derived_root = self._resolve(temporary_directory, ["10", "11", "12"])
            latest_directory = derived_root / "latest"
            manifest_payload = payload["manifest"]

            fasta_file_path = latest_directory / "protein_sequences.fasta"
            fasta_text = fasta_file_path.read_text(encoding="utf-8")
            records = [block for block in fasta_text.split(">") if block.strip()]
            reordered = ">" + ">".join(reversed(records))
            fasta_file_path.write_text(reordered, encoding="utf-8")

            with self.assertRaisesRegex(RuntimeError, "mismatch"):
                _validate_loaded_derived_protein_fasta_payload(
                    snapshot_directory=latest_directory,
                    manifest_payload=manifest_payload,
                )

            # With fasta_file_sha256 patched to the tampered value, the
            # record-id set check still catches the reordering.
            from src.pago_pipeline.storage import sha256_of_file as _sha

            patched_manifest = {
                **manifest_payload,
                "fasta_file_sha256": _sha(input_file_path=fasta_file_path),
            }
            with self.assertRaisesRegex(RuntimeError, "record-id set"):
                _validate_loaded_derived_protein_fasta_payload(
                    snapshot_directory=latest_directory,
                    manifest_payload=patched_manifest,
                )

    def test_latest_unavailable_when_selection_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            payload, metadata_root, selection_root, derived_root = self._resolve(
                temporary_directory, ["10", "11"]
            )
            self.assertTrue(
                latest_derived_protein_fasta_snapshot_is_available(
                    snapshot_root_directory=derived_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    source_selection_snapshot_root_directory=selection_root,
                    selection_artifact_type=_SELECTION_ARTIFACT_TYPE,
                    selection_uid_list_file_name=_SELECTION_UID_LIST_FILE_NAME,
                )
            )
            _write_selection_latest_snapshot(
                snapshot_root_directory=selection_root,
                retained_uids=["10", "11", "12"],
                snapshot_directory_name="selection_2",
            )
            self.assertFalse(
                latest_derived_protein_fasta_snapshot_is_available(
                    snapshot_root_directory=derived_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    source_selection_snapshot_root_directory=selection_root,
                    selection_artifact_type=_SELECTION_ARTIFACT_TYPE,
                    selection_uid_list_file_name=_SELECTION_UID_LIST_FILE_NAME,
                )
            )

    def test_reuse_latest_or_create_reuses_without_resaving(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            payload, metadata_root, selection_root, derived_root = self._resolve(
                temporary_directory, ["10", "11"]
            )
            with mock.patch(
                "src.pago_pipeline.derived_protein_fasta_snapshot"
                ".save_derived_protein_fasta_snapshot"
            ) as save_mock:
                resolve_derived_protein_fasta_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=derived_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    source_selection_snapshot_root_directory=selection_root,
                    selection_artifact_type=_SELECTION_ARTIFACT_TYPE,
                    selection_uid_list_file_name=_SELECTION_UID_LIST_FILE_NAME,
                    record_selection_rule="pago_technical_prefilter.retained",
                    record_selection_config_sha256="policy-sha",
                    dataset_kind="annotation_enriched_proteome",
                )
                save_mock.assert_not_called()


if __name__ == "__main__":
    unittest.main()
