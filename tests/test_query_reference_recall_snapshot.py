from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.query_reference_recall import (
    build_matching_strategy_sha256,
    protein_sequence_sha256,
)
from src.pago_pipeline.query_reference_recall_snapshot import (
    ARTIFACT_TYPE,
    DEFAULT_MANIFEST_FILE_NAME,
    SNAPSHOT_FORMAT_VERSION,
    latest_query_reference_recall_snapshot_is_available,
    resolve_query_reference_recall_snapshot,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic

_FIXTURE_CSV = (
    Path(__file__).resolve().parent / "fixtures" / "query_recall_reference_set.csv"
)

# A valid 12-column reference row (matches the committed fixture schema).
_EXTRA_REFERENCE_ROW = (
    "WP_7770001.1,Zz,Zz sp,PAGO,LONG_A,"
    + protein_sequence_sha256("MZZZZ")
    + ",5,,src,LITERATURE_PHYLOGENETIC,provisional,synthetic test row\n"
)


def _write_metadata_latest_snapshot(
    *,
    snapshot_root_directory: Path,
    accession_versions: list[str],
    sequences: list[str] | None = None,
    snapshot_directory_name: str = "metadata_1",
) -> Path:
    directory = snapshot_root_directory / "snapshots" / snapshot_directory_name
    directory.mkdir(parents=True, exist_ok=True)
    csv_file_path = directory / "protein_metadata.csv"
    frame = pd.DataFrame(
        {
            "protein_uid": [str(100 + i) for i in range(len(accession_versions))],
            "gbseq__accession_version": accession_versions,
        }
    )
    if sequences is not None:
        frame["gbseq__sequence"] = sequences
    frame.to_csv(csv_file_path, index=False)
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
        reference_frame = pd.read_csv(_FIXTURE_CSV, dtype=str)
        reference_accessions = list(reference_frame["accession"])
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
            self.assertEqual(
                manifest["snapshot_format_version"], SNAPSHOT_FORMAT_VERSION
            )
            self.assertEqual(
                manifest["matching_strategy_sha256"],
                build_matching_strategy_sha256(),
            )
            self.assertEqual(manifest["exact_recovered_count"], 2)
            self.assertEqual(manifest["equivalent_recovered_count"], 2)
            self.assertEqual(
                manifest["reference_count"], len(reference_accessions)
            )
            for reading in ("stratum_exact_recall", "stratum_equivalent_recall"):
                self.assertIn("overall_reference_recall", manifest[reading])
                self.assertIn("long_a_reference_recall", manifest[reading])
            # The two retrieved accessions are both LONG_A (TtAgo, MjAgo).
            self.assertAlmostEqual(
                manifest["stratum_exact_recall"]["long_a_reference_recall"],
                2 / 8,
            )
            # LONG_B / SHORT / PIWI_RE all present but unrecovered here: a
            # legitimate EVALUABLE 0.0, never NOT_EVALUABLE.
            for evaluable_metric in (
                "overall_reference_recall",
                "long_a_reference_recall",
                "long_b_reference_recall",
                "short_reference_recall",
                "piwi_re_reference_recall",
            ):
                self.assertEqual(
                    manifest["stratum_recall_status"][evaluable_metric],
                    "EVALUABLE",
                )
                self.assertEqual(
                    manifest["stratum_exact_recall"]["piwi_re_reference_recall"],
                    0.0,
                )
            self.assertEqual(
                manifest["query_recall_reference_set_csv_sha256"],
                sha256_of_file(input_file_path=_FIXTURE_CSV),
            )
            detail = payload["detail"]
            self.assertEqual(int(detail["recovered"].sum()), 2)
            self.assertIn("match_method", detail.columns)
            self.assertIn("sequence_sha256", detail.columns)

    def test_recovery_by_sequence_identity_under_alias_accession(self) -> None:
        # A reference protein retrieved under a different accession, recognised
        # purely by identical (whitespace-stripped, uppercased) protein sequence.
        alias_sequence = "M" + "AGIDENTITY" * 12
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            metadata_root = temporary_directory / "metadata"
            recall_root = temporary_directory / "recall"
            reference_csv = temporary_directory / "ref.csv"
            reference_csv.write_text(
                "accession,protein_short_name,organism,ago_family,clade,"
                "sequence_sha256,sequence_length,uniprot_accession,"
                "reference_label_source,reference_label_evidence,"
                "verification_status,notes\n"
                "ABP72561.1,RsAgo,Cereibacter sphaeroides,PAGO,LONG_B,"
                + protein_sequence_sha256(alias_sequence)
                + ",120,A4WYU7,src,EXPERIMENTAL,verified,alias test\n"
                "WP_010870838.1,MjAgo,Methanocaldococcus jannaschii,PAGO,LONG_A,"
                + protein_sequence_sha256("MJNOTRETRIEVED")
                + ",100,Q58717,src,EXPERIMENTAL,verified,not retrieved\n",
                encoding="utf-8",
            )
            _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                accession_versions=["A4WYU7.1", "WP_999999999.1"],
                sequences=[alias_sequence.lower(), "OTHERSEQUENCE"],
            )

            payload = resolve_query_reference_recall_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                snapshot_root_directory=recall_root,
                source_metadata_snapshot_root_directory=metadata_root,
                query_recall_reference_set_csv_path=reference_csv,
            )
            manifest = payload["manifest"]
            self.assertEqual(manifest["exact_recovered_count"], 0)
            self.assertEqual(manifest["equivalent_recovered_count"], 1)
            self.assertEqual(
                manifest["stratum_exact_recall"]["long_b_reference_recall"], 0.0
            )
            self.assertEqual(
                manifest["stratum_equivalent_recall"]["long_b_reference_recall"],
                1.0,
            )
            detail = payload["detail"].set_index("accession")
            rsago = detail.loc["ABP72561.1"]
            self.assertEqual(rsago["match_method"], "SEQUENCE_SHA256")
            self.assertEqual(rsago["matched_accession"], "A4WYU7.1")
            self.assertEqual(str(rsago["matched_protein_uid"]), "100")
            self.assertTrue(bool(rsago["recovered"]))
            self.assertFalse(bool(rsago["recovered_exact_accession"]))

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
                reference_copy.read_text(encoding="utf-8") + _EXTRA_REFERENCE_ROW,
                encoding="utf-8",
            )
            self.assertFalse(
                latest_query_reference_recall_snapshot_is_available(
                    snapshot_root_directory=recall_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    query_recall_reference_set_csv_path=reference_copy,
                )
            )

    def test_snapshot_written_by_earlier_methodology_is_not_reused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            metadata_root = temporary_directory / "metadata"
            recall_root = temporary_directory / "recall"
            _write_metadata_latest_snapshot(
                snapshot_root_directory=metadata_root,
                accession_versions=["WP_000000000.1"],
            )
            resolve_query_reference_recall_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest_or_create,
                snapshot_root_directory=recall_root,
                source_metadata_snapshot_root_directory=metadata_root,
                query_recall_reference_set_csv_path=_FIXTURE_CSV,
            )
            latest_manifest = recall_root / "latest" / DEFAULT_MANIFEST_FILE_NAME
            import json

            payload = json.loads(latest_manifest.read_text(encoding="utf-8"))
            payload["matching_strategy_sha256"] = "0" * 64
            latest_manifest.write_text(json.dumps(payload), encoding="utf-8")

            self.assertFalse(
                latest_query_reference_recall_snapshot_is_available(
                    snapshot_root_directory=recall_root,
                    source_metadata_snapshot_root_directory=metadata_root,
                    query_recall_reference_set_csv_path=_FIXTURE_CSV,
                )
            )


if __name__ == "__main__":
    unittest.main()
