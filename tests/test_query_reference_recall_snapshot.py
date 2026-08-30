from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.query_reference_recall_snapshot import (
    ARTIFACT_TYPE,
    DEFAULT_MANIFEST_FILE_NAME,
    latest_query_reference_recall_snapshot_is_available,
    resolve_query_reference_recall_snapshot,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic

_FIXTURE_CSV = (
    Path(__file__).resolve().parent / "fixtures" / "query_recall_reference_set.csv"
)


def _write_metadata_latest_snapshot(
    *,
    snapshot_root_directory: Path,
    accession_versions: list[str],
    snapshot_directory_name: str = "metadata_1",
) -> Path:
    directory = snapshot_root_directory / "snapshots" / snapshot_directory_name
    directory.mkdir(parents=True, exist_ok=True)
    csv_file_path = directory / "protein_metadata.csv"
    pd.DataFrame(
        {
            "protein_uid": [str(100 + i) for i in range(len(accession_versions))],
            "gbseq__accession_version": accession_versions,
        }
    ).to_csv(csv_file_path, index=False)
    write_json_atomic(payload={"checks": {}}, output_file_path=directory / "qc_report.json")
    write_json_atomic(
        payload={
            "artifact_type": "ncbi_protein_metadata_snapshot",
            "snapshot_format_version": "1.0",
            "csv_file_name": "protein_metadata.csv",
            "csv_file_sha256": sha256_of_file(input_file_path=csv_file_path),
            "qc_report_file_name": "qc_report.json",
            "immutable_snapshot_directory_name": snapshot_directory_name,
            "immutable_snapshot_relative_path": (
                f"snapshots/{snapshot_directory_name}"
            ),
            "row_count": len(accession_versions),
            "search_query": (
                "(PIWI[All Fields] OR Argonaute[All Fields]) AND "
                "(Bacteria[Organism] OR Archaea[Organism])"
            ),
            "translated_query": "translated",
        },
        output_file_path=directory / DEFAULT_MANIFEST_FILE_NAME,
    )
    latest = snapshot_root_directory / "latest"
    if latest.exists():
        shutil.rmtree(latest)
    shutil.copytree(directory, latest)
    return directory


class QueryReferenceRecallSnapshotTests(unittest.TestCase):
    def test_resolve_computes_recall_against_committed_reference_set(self) -> None:
        reference_accessions = list(
            pd.read_csv(_FIXTURE_CSV, dtype=str)["accession"]
        )
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            metadata_root = temporary_directory / "metadata"
            recall_root = temporary_directory / "recall"
            # Retrieve exactly the first two reference accessions, plus noise.
            _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                accession_versions=[
                    reference_accessions[0],
                    reference_accessions[1],
                    "WP_999999999.1",
                ],
            )

            payload = resolve_query_reference_recall_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                snapshot_root_directory=recall_root,
                source_metadata_snapshot_root_directory=metadata_root,
                query_recall_reference_set_csv_path=_FIXTURE_CSV,
            )
            manifest = payload["manifest"]
            self.assertEqual(manifest["artifact_type"], ARTIFACT_TYPE)
            self.assertEqual(manifest["recovered_count"], 2)
            self.assertEqual(
                manifest["reference_count"], len(reference_accessions)
            )
            self.assertIn("overall_reference_recall", manifest["stratum_recall"])
            self.assertIn("long_a_reference_recall", manifest["stratum_recall"])
            # The committed reference set has LONG_A / LONG_B / SHORT but no
            # PIWI_RE, so that stratum must be NOT_EVALUABLE (null), never 0.0.
            self.assertIsNone(manifest["stratum_recall"]["piwi_re_reference_recall"])
            self.assertEqual(
                manifest["stratum_recall_status"]["piwi_re_reference_recall"],
                "NOT_EVALUABLE",
            )
            for evaluable_metric in (
                "overall_reference_recall",
                "long_a_reference_recall",
                "long_b_reference_recall",
                "short_reference_recall",
            ):
                self.assertEqual(
                    manifest["stratum_recall_status"][evaluable_metric],
                    "EVALUABLE",
                )
            self.assertEqual(
                manifest["query_recall_reference_set_csv_sha256"],
                sha256_of_file(input_file_path=_FIXTURE_CSV),
            )
            self.assertEqual(int(payload["detail"]["recovered"].sum()), 2)

    def test_reuse_and_invalidation_on_reference_set_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            metadata_root = temporary_directory / "metadata"
            recall_root = temporary_directory / "recall"
            reference_copy = temporary_directory / "ref.csv"
            reference_copy.write_bytes(_FIXTURE_CSV.read_bytes())
            _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                accession_versions=["WP_000000000.1"],
            )
            resolve_query_reference_recall_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                snapshot_root_directory=recall_root,
                source_metadata_snapshot_root_directory=metadata_root,
                query_recall_reference_set_csv_path=reference_copy,
            )
            self.assertTrue(
                latest_query_reference_recall_snapshot_is_available(
                    snapshot_root_directory=recall_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    query_recall_reference_set_csv_path=reference_copy,
                )
            )
            with mock.patch(
                "src.pago_pipeline.query_reference_recall_snapshot"
                ".save_query_reference_recall_snapshot"
            ) as save_mock:
                resolve_query_reference_recall_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=recall_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    query_recall_reference_set_csv_path=reference_copy,
                )
                save_mock.assert_not_called()

            reference_copy.write_text(
                reference_copy.read_text(encoding="utf-8")
                + "WP_777.1,Zz,Zz sp,PAGO,PIWI_RE,,src,LITERATURE_PHYLOGENETIC,provisional,x\n",
                encoding="utf-8",
            )
            self.assertFalse(
                latest_query_reference_recall_snapshot_is_available(
                    snapshot_root_directory=recall_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    query_recall_reference_set_csv_path=reference_copy,
                )
            )


if __name__ == "__main__":
    unittest.main()
