from __future__ import annotations

import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.query_reference_recall import build_matching_strategy_sha256
from src.pago_pipeline.query_reference_recall_snapshot import (
    ARTIFACT_TYPE,
    DEFAULT_DETAIL_FILE_NAME,
    DEFAULT_MANIFEST_FILE_NAME,
    DEFAULT_SUMMARY_FILE_NAME,
    SNAPSHOT_FORMAT_VERSION,
    latest_query_reference_recall_snapshot_is_available,
    load_latest_query_reference_recall_snapshot,
    resolve_query_reference_recall_snapshot,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic

_TTAGO_ACCESSION = "WP_011174533.1"
_METADATA_SNAPSHOT_DIRECTORY_NAME = "metadata_1"

# Byte identity of the summary CSV for a fixed logical input, pinned so that
# the Linux and Windows jobs both prove the written bytes agree. Recomputed
# deliberately if the summary schema or the serialized representation changes.
_EXPECTED_SUMMARY_SHA256 = (
    "50c8ee73abc28ec100a2f171f9deaafdee77b41d23f9bd95b34c8d36dcfc9119"
)


def _write_metadata_snapshot(
    *,
    snapshot_root_directory: Path,
    rows: list[tuple[str, str, str]],
    declare_csv_sha256: bool = True,
    directory_name: str = _METADATA_SNAPSHOT_DIRECTORY_NAME,
) -> Path:
    """Build a complete synthetic metadata snapshot, offline.

    No network and no committed retrieval are needed: the upstream contract is
    a manifest plus a flattened CSV plus a QC report, all of which are written
    here.
    """
    directory = snapshot_root_directory / "snapshots" / directory_name
    directory.mkdir(parents=True, exist_ok=True)
    csv_file_path = directory / "protein_metadata.csv"
    with csv_file_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "protein_uid",
                "gbseq__accession_version",
                "gbseq__sequence",
            ],
        )
        writer.writeheader()
        for protein_uid, accession, sequence in rows:
            writer.writerow(
                {
                    "protein_uid": protein_uid,
                    "gbseq__accession_version": accession,
                    "gbseq__sequence": sequence,
                }
            )

    qc_report_file_path = directory / "qc_report.json"
    write_json_atomic(payload={"checks": {}}, output_file_path=qc_report_file_path)

    manifest_payload: dict[str, object] = {
        "artifact_type": "ncbi_protein_metadata_snapshot",
        "snapshot_format_version": "1.0",
        "csv_file_name": csv_file_path.name,
        "qc_report_file_name": qc_report_file_path.name,
        "immutable_snapshot_directory_name": directory_name,
        "immutable_snapshot_relative_path": f"snapshots/{directory_name}",
        "row_count": len(rows),
        "search_query": "(PIWI[All Fields] OR Argonaute[All Fields])",
        "translated_query": "translated",
    }
    if declare_csv_sha256:
        manifest_payload["csv_file_sha256"] = sha256_of_file(
            input_file_path=csv_file_path
        )
    write_json_atomic(
        payload=manifest_payload,
        output_file_path=directory / DEFAULT_MANIFEST_FILE_NAME,
    )

    latest_directory = snapshot_root_directory / "latest"
    if latest_directory.exists():
        shutil.rmtree(latest_directory)
    shutil.copytree(directory, latest_directory)
    return directory


def _rewrite_metadata_csv(
    *, snapshot_root_directory: Path, rows: list[tuple[str, str, str]]
) -> None:
    """Change the metadata CSV bytes while leaving both manifests untouched."""
    for directory in (
        snapshot_root_directory / "snapshots" / _METADATA_SNAPSHOT_DIRECTORY_NAME,
        snapshot_root_directory / "latest",
    ):
        csv_file_path = directory / "protein_metadata.csv"
        with csv_file_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "protein_uid",
                    "gbseq__accession_version",
                    "gbseq__sequence",
                ],
            )
            writer.writeheader()
            for protein_uid, accession, sequence in rows:
                writer.writerow(
                    {
                        "protein_uid": protein_uid,
                        "gbseq__accession_version": accession,
                        "gbseq__sequence": sequence,
                    }
                )


def _mutate_recall_manifest(recall_root_directory: Path, **overrides: object) -> None:
    manifest_file_path = (
        recall_root_directory / "latest" / DEFAULT_MANIFEST_FILE_NAME
    )
    payload = json.loads(manifest_file_path.read_text(encoding="utf-8"))
    payload.update(overrides)
    manifest_file_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )


class _SnapshotHarness(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory()
        root = Path(self._temporary_directory.name)
        self.metadata_root = root / "metadata"
        self.recall_root = root / "recall"
        self.metadata_root.mkdir(parents=True, exist_ok=True)

    def tearDown(self) -> None:
        self._temporary_directory.cleanup()

    def _resolve(self, snapshot_mode: SnapshotMode | str) -> dict[str, object]:
        return resolve_query_reference_recall_snapshot(
            snapshot_mode=snapshot_mode,
            snapshot_root_directory=self.recall_root,
            source_metadata_snapshot_root_directory=self.metadata_root,
        )


class SnapshotRoundTripTests(_SnapshotHarness):
    """T-R1 … T-R4, T-X1, T-X2."""

    def setUp(self) -> None:
        super().setUp()
        _write_metadata_snapshot(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MAAA")],
        )

    def test_create_new_then_reuse_latest_returns_the_same_result(self) -> None:
        created = self._resolve(SnapshotMode.create_new)
        reused = self._resolve(SnapshotMode.reuse_latest)
        self.assertEqual(created["summary_rows"], reused["summary_rows"])
        self.assertEqual(created["detail_rows"], reused["detail_rows"])
        self.assertEqual(created["manifest"], reused["manifest"])
        manifest = created["manifest"]
        self.assertEqual(manifest["artifact_type"], ARTIFACT_TYPE)
        self.assertEqual(
            manifest["snapshot_format_version"], SNAPSHOT_FORMAT_VERSION
        )
        self.assertEqual(manifest["reference_count"], 21)
        self.assertEqual(manifest["intended_use"], "NCBI_QUERY_COVERAGE_AUDIT")
        self.assertIs(manifest["detector_validation_eligible"], False)
        self.assertIs(manifest["sequence_evidence_complete"], True)
        self.assertIs(manifest["source_metadata_csv_sha256_declared_upstream"], True)

    def test_reuse_latest_or_create_reuses_a_compatible_result(self) -> None:
        created = self._resolve(SnapshotMode.create_new)
        reused = self._resolve(SnapshotMode.reuse_latest_or_create)
        self.assertEqual(
            created["snapshot_directory"].name,
            Path(str(reused["manifest"]["immutable_snapshot_directory_name"])).name,
        )

    def test_reuse_latest_without_a_stored_result_raises(self) -> None:
        with self.assertRaisesRegex(FileNotFoundError, "create_new"):
            self._resolve(SnapshotMode.reuse_latest)

    def test_two_independent_runs_agree_byte_for_byte(self) -> None:
        first = self._resolve(SnapshotMode.create_new)
        second = self._resolve(SnapshotMode.create_new)
        # Two runs inside the same second must each claim their own immutable
        # directory rather than colliding on the second-resolution name.
        self.assertNotEqual(
            first["snapshot_directory"], second["snapshot_directory"]
        )
        for key in ("summary_file_path", "detail_file_path"):
            self.assertEqual(
                Path(str(first[key])).read_bytes(),
                Path(str(second[key])).read_bytes(),
                key,
            )
        varying = {
            "snapshot_created_at_utc",
            "immutable_snapshot_directory_name",
            "immutable_snapshot_relative_path",
            "output_files",
        }
        first_manifest = {
            key: value
            for key, value in dict(first["manifest"]).items()
            if key not in varying
        }
        second_manifest = {
            key: value
            for key, value in dict(second["manifest"]).items()
            if key not in varying
        }
        self.assertEqual(first_manifest, second_manifest)
        # Only the recorded paths differ inside output_files; every hash agrees.
        for key in ("summary_file", "detail_file"):
            self.assertEqual(
                dict(first["manifest"])["output_files"][key]["sha256"],
                dict(second["manifest"])["output_files"][key]["sha256"],
                key,
            )

    def test_written_csv_bytes_are_platform_independent(self) -> None:
        created = self._resolve(SnapshotMode.create_new)
        summary_bytes = Path(str(created["summary_file_path"])).read_bytes()
        self.assertIn(b"\r\n", summary_bytes)
        self.assertEqual(summary_bytes.replace(b"\r\n", b"").count(b"\n"), 0)
        self.assertEqual(
            hashlib.sha256(summary_bytes).hexdigest(), _EXPECTED_SUMMARY_SHA256
        )

    def test_not_evaluable_strata_survive_the_csv_round_trip(self) -> None:
        created = self._resolve(SnapshotMode.create_new)
        # The committed panel fills every stratum, so assert the serialized
        # representation of the booleans and metrics that are present.
        for row in created["summary_rows"]:
            self.assertIn(
                row["recall_status"], {"EVALUABLE", "NOT_EVALUABLE"}
            )
            self.assertIn(
                row["retrieval_equivalent_recall_is_lower_bound"],
                {"true", "false"},
            )
            if row["recall_status"] == "NOT_EVALUABLE":
                self.assertEqual(row["retrieval_equivalent_recall"], "")
                self.assertEqual(row["retrieval_equivalent_recall_is_lower_bound"], "false")


class SnapshotInvalidationTests(_SnapshotHarness):
    """T-I1 … T-I10."""

    def setUp(self) -> None:
        super().setUp()
        _write_metadata_snapshot(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MAAA")],
        )
        self._resolve(SnapshotMode.create_new)

    def _available(self) -> bool:
        return latest_query_reference_recall_snapshot_is_available(
            snapshot_root_directory=self.recall_root,
            source_metadata_snapshot_root_directory=self.metadata_root,
            panel_directory=None,
        )

    def _available_with_panel(self) -> bool:
        from src.pago_pipeline.query_recall_reference_panel import (
            DEFAULT_PANEL_DIRECTORY,
        )

        return latest_query_reference_recall_snapshot_is_available(
            snapshot_root_directory=self.recall_root,
            source_metadata_snapshot_root_directory=self.metadata_root,
            panel_directory=DEFAULT_PANEL_DIRECTORY,
        )

    def test_a_compatible_stored_result_is_available(self) -> None:
        self.assertTrue(self._available())
        self.assertTrue(self._available_with_panel())

    def test_changed_panel_lock_hash_invalidates(self) -> None:
        _mutate_recall_manifest(self.recall_root, panel_lock_sha256="0" * 64)
        self.assertFalse(self._available_with_panel())

    def test_changed_matching_strategy_hash_invalidates(self) -> None:
        _mutate_recall_manifest(self.recall_root, matching_strategy_sha256="0" * 64)
        self.assertFalse(self._available())

    def test_changed_snapshot_format_version_invalidates(self) -> None:
        _mutate_recall_manifest(self.recall_root, snapshot_format_version="1.1")
        self.assertFalse(self._available())

    def test_changed_source_manifest_hash_invalidates(self) -> None:
        _mutate_recall_manifest(
            self.recall_root, source_metadata_manifest_sha256="0" * 64
        )
        self.assertFalse(self._available())

    def test_changed_metadata_csv_with_declared_upstream_hash_rejects_the_input(
        self,
    ) -> None:
        # The upstream manifest declares csv_file_sha256, so the upstream
        # loader rejects the input itself — in every mode, before computing.
        _rewrite_metadata_csv(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MZZZ")],
        )
        for snapshot_mode in (
            SnapshotMode.create_new,
            SnapshotMode.reuse_latest,
            SnapshotMode.reuse_latest_or_create,
        ):
            with self.subTest(snapshot_mode=snapshot_mode):
                with self.assertRaises(RuntimeError):
                    self._resolve(snapshot_mode)

    def test_changed_metadata_csv_without_declared_upstream_hash_blocks_reuse(
        self,
    ) -> None:
        # Rebuild the upstream snapshot without csv_file_sha256, so the
        # upstream loader accepts a changed CSV and P11's own comparison is the
        # only thing that notices.
        shutil.rmtree(self.metadata_root)
        self.metadata_root.mkdir(parents=True, exist_ok=True)
        _write_metadata_snapshot(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MAAA")],
            declare_csv_sha256=False,
        )
        shutil.rmtree(self.recall_root, ignore_errors=True)
        created = self._resolve(SnapshotMode.create_new)
        self.assertIs(
            dict(created["manifest"])["source_metadata_csv_sha256_declared_upstream"],
            False,
        )
        self.assertTrue(self._available())

        _rewrite_metadata_csv(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MZZZ")],
        )
        self.assertFalse(self._available())
        with self.assertRaisesRegex(FileNotFoundError, "create_new"):
            self._resolve(SnapshotMode.reuse_latest)
        recreated = self._resolve(SnapshotMode.reuse_latest_or_create)
        self.assertEqual(
            dict(recreated["manifest"])["source_metadata_csv_sha256"],
            sha256_of_file(
                input_file_path=self.metadata_root
                / "snapshots"
                / _METADATA_SNAPSHOT_DIRECTORY_NAME
                / "protein_metadata.csv"
            ),
        )

    def test_a_methodology_change_alone_invalidates_a_stored_result(self) -> None:
        # Same data inputs, different declared policy: the strategy hash is the
        # thing that must notice.
        self.assertTrue(self._available())
        _mutate_recall_manifest(
            self.recall_root,
            matching_strategy_sha256=hashlib.sha256(
                (build_matching_strategy_sha256() + "policy-changed").encode("utf-8")
            ).hexdigest(),
        )
        self.assertFalse(self._available())

    def test_a_missing_output_file_fails_integrity(self) -> None:
        (self.recall_root / "latest" / DEFAULT_DETAIL_FILE_NAME).unlink()
        self.assertFalse(self._available())
        with self.assertRaises(FileNotFoundError):
            self._resolve(SnapshotMode.reuse_latest)

    def test_a_tampered_output_file_fails_integrity(self) -> None:
        tampered = self.recall_root / "latest" / DEFAULT_SUMMARY_FILE_NAME
        tampered.write_text("stratum,metric_name\r\n", encoding="utf-8")
        self.assertFalse(self._available())

    def test_reuse_latest_does_not_fall_back_to_creating(self) -> None:
        _mutate_recall_manifest(self.recall_root, snapshot_format_version="1.1")
        saved_before = sorted(
            path.name for path in (self.recall_root / "snapshots").iterdir()
        )
        with self.assertRaisesRegex(FileNotFoundError, "create_new"):
            self._resolve(SnapshotMode.reuse_latest)
        saved_after = sorted(
            path.name for path in (self.recall_root / "snapshots").iterdir()
        )
        self.assertEqual(saved_before, saved_after)


class ImmutableSourceValidationTests(_SnapshotHarness):
    """latest/ selects the source; the immutable snapshot is what is validated."""

    def setUp(self) -> None:
        super().setUp()
        self.immutable_directory = _write_metadata_snapshot(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MAAA")],
        )

    def _rewrite_immutable_csv(self, rows: list[tuple[str, str, str]]) -> None:
        csv_file_path = self.immutable_directory / "protein_metadata.csv"
        with csv_file_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "protein_uid",
                    "gbseq__accession_version",
                    "gbseq__sequence",
                ],
            )
            writer.writeheader()
            for protein_uid, accession, sequence in rows:
                writer.writerow(
                    {
                        "protein_uid": protein_uid,
                        "gbseq__accession_version": accession,
                        "gbseq__sequence": sequence,
                    }
                )

    def _mutate_immutable_manifest(self, **overrides: object) -> None:
        manifest_file_path = self.immutable_directory / DEFAULT_MANIFEST_FILE_NAME
        payload = json.loads(manifest_file_path.read_text(encoding="utf-8"))
        payload.update(overrides)
        manifest_file_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
        )

    def test_a_corrupt_immutable_csv_is_rejected_even_when_latest_validates(
        self,
    ) -> None:
        # latest/ is untouched and still validates. Only the immutable CSV
        # changes, so its bytes no longer match the hash its own manifest
        # declares. The input must be rejected before any computation, in every
        # mode, rather than measured over the corrupt bytes.
        self._rewrite_immutable_csv([("9", "WP_999999.1", "MZZZ")])
        for snapshot_mode in (
            SnapshotMode.create_new,
            SnapshotMode.reuse_latest,
            SnapshotMode.reuse_latest_or_create,
        ):
            with self.subTest(snapshot_mode=snapshot_mode):
                with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                    self._resolve(snapshot_mode)

    def test_a_divergent_immutable_manifest_is_rejected(self) -> None:
        # Provenance must not mix one manifest's fields with another's hash.
        self._mutate_immutable_manifest(search_query="IMMUTABLE-ONLY-QUERY")
        with self.assertRaisesRegex(RuntimeError, "differ"):
            self._resolve(SnapshotMode.create_new)

    def test_provenance_is_copied_from_the_manifest_actually_used(self) -> None:
        created = self._resolve(SnapshotMode.create_new)
        manifest = dict(created["manifest"])
        immutable_manifest = json.loads(
            (self.immutable_directory / DEFAULT_MANIFEST_FILE_NAME).read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            manifest["search_query"], immutable_manifest["search_query"]
        )
        self.assertEqual(
            manifest["translated_query"], immutable_manifest["translated_query"]
        )
        self.assertEqual(
            manifest["source_metadata_manifest_sha256"],
            sha256_of_file(
                input_file_path=self.immutable_directory / DEFAULT_MANIFEST_FILE_NAME
            ),
        )
        self.assertEqual(
            manifest["source_metadata_csv_sha256"],
            sha256_of_file(
                input_file_path=self.immutable_directory / "protein_metadata.csv"
            ),
        )
        self.assertEqual(
            manifest["source_metadata_snapshot_relative_path"],
            immutable_manifest["immutable_snapshot_relative_path"],
        )

    def test_a_source_without_a_declared_csv_hash_still_resolves(self) -> None:
        shutil.rmtree(self.metadata_root)
        self.metadata_root.mkdir(parents=True, exist_ok=True)
        immutable_directory = _write_metadata_snapshot(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MAAA")],
            declare_csv_sha256=False,
        )
        created = self._resolve(SnapshotMode.create_new)
        manifest = dict(created["manifest"])
        self.assertIs(
            manifest["source_metadata_csv_sha256_declared_upstream"], False
        )
        self.assertEqual(
            manifest["source_metadata_csv_sha256"],
            sha256_of_file(
                input_file_path=immutable_directory / "protein_metadata.csv"
            ),
        )

    def test_a_missing_immutable_source_does_not_fall_back_to_latest(self) -> None:
        shutil.rmtree(self.immutable_directory)
        with self.assertRaisesRegex(
            (FileNotFoundError, RuntimeError), "immutable"
        ):
            self._resolve(SnapshotMode.create_new)

    def test_a_source_that_declares_no_immutable_path_is_rejected(self) -> None:
        for directory in (self.immutable_directory, self.metadata_root / "latest"):
            manifest_file_path = directory / DEFAULT_MANIFEST_FILE_NAME
            payload = json.loads(manifest_file_path.read_text(encoding="utf-8"))
            payload.pop("immutable_snapshot_relative_path", None)
            manifest_file_path.write_text(
                json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
            )
        with self.assertRaisesRegex(RuntimeError, "immutable_snapshot_relative_path"):
            self._resolve(SnapshotMode.create_new)


class RequiredOutputFileEntryTests(_SnapshotHarness):
    """A manifest declaring fewer than both outputs describes an incomplete result."""

    def setUp(self) -> None:
        super().setUp()
        _write_metadata_snapshot(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, "MAAA")],
        )
        self._resolve(SnapshotMode.create_new)
        self.latest_manifest_file_path = (
            self.recall_root / "latest" / DEFAULT_MANIFEST_FILE_NAME
        )

    def _latest_output_files(self) -> dict[str, object]:
        payload = json.loads(
            self.latest_manifest_file_path.read_text(encoding="utf-8")
        )
        return dict(payload["output_files"])

    def _available(self) -> bool:
        return latest_query_reference_recall_snapshot_is_available(
            snapshot_root_directory=self.recall_root,
            source_metadata_snapshot_root_directory=self.metadata_root,
        )

    def test_a_complete_manifest_loads_and_is_reusable(self) -> None:
        self.assertTrue(self._available())
        loaded = self._resolve(SnapshotMode.reuse_latest)
        self.assertEqual(len(loaded["summary_rows"]), 6)
        self.assertEqual(len(loaded["detail_rows"]), 21)

    def test_empty_output_files_is_unavailable_and_fails_to_load(self) -> None:
        _mutate_recall_manifest(self.recall_root, output_files={})
        self.assertFalse(self._available())
        with self.assertRaisesRegex(RuntimeError, "summary_file"):
            load_latest_query_reference_recall_snapshot(
                snapshot_root_directory=self.recall_root
            )

    def test_a_missing_summary_entry_is_rejected(self) -> None:
        output_files = self._latest_output_files()
        del output_files["summary_file"]
        _mutate_recall_manifest(self.recall_root, output_files=output_files)
        self.assertFalse(self._available())
        with self.assertRaisesRegex(RuntimeError, "summary_file"):
            load_latest_query_reference_recall_snapshot(
                snapshot_root_directory=self.recall_root
            )

    def test_a_missing_detail_entry_is_rejected(self) -> None:
        output_files = self._latest_output_files()
        del output_files["detail_file"]
        _mutate_recall_manifest(self.recall_root, output_files=output_files)
        self.assertFalse(self._available())
        with self.assertRaisesRegex(RuntimeError, "detail_file"):
            load_latest_query_reference_recall_snapshot(
                snapshot_root_directory=self.recall_root
            )

    def test_a_structurally_invalid_required_entry_is_rejected(self) -> None:
        output_files = self._latest_output_files()
        output_files["detail_file"] = "not-an-object"
        _mutate_recall_manifest(self.recall_root, output_files=output_files)
        self.assertFalse(self._available())
        with self.assertRaisesRegex(RuntimeError, "detail_file"):
            load_latest_query_reference_recall_snapshot(
                snapshot_root_directory=self.recall_root
            )

    def test_a_required_entry_without_a_hash_is_rejected(self) -> None:
        output_files = self._latest_output_files()
        entry = dict(output_files["summary_file"])
        del entry["sha256"]
        output_files["summary_file"] = entry
        _mutate_recall_manifest(self.recall_root, output_files=output_files)
        self.assertFalse(self._available())
        with self.assertRaisesRegex(RuntimeError, "summary_file"):
            load_latest_query_reference_recall_snapshot(
                snapshot_root_directory=self.recall_root
            )


class SnapshotEvidenceTests(_SnapshotHarness):
    """Completeness and the use contract survive persistence."""

    def test_incomplete_sequence_evidence_is_recorded_in_the_manifest(self) -> None:
        _write_metadata_snapshot(
            snapshot_root_directory=self.metadata_root,
            rows=[("1", _TTAGO_ACCESSION, ""), ("2", "WP_999999.1", "")],
        )
        created = self._resolve(SnapshotMode.create_new)
        manifest = dict(created["manifest"])
        self.assertIs(manifest["sequence_evidence_complete"], False)
        self.assertEqual(manifest["retrieved_row_count"], 2)
        self.assertEqual(manifest["retrieved_rows_with_sequence_count"], 0)
        lower_bound = manifest["stratum_equivalent_is_lower_bound"]
        # TtAgo is recovered by accession, so LONG_A is unsaturated only
        # because the panel holds nine LONG_A references.
        self.assertIs(lower_bound["long_a_reference_recall"], True)
        self.assertIs(lower_bound["overall_reference_recall"], True)
        for value in lower_bound.values():
            self.assertIsInstance(value, bool)
        for value in manifest["stratum_recall_status"].values():
            self.assertIsInstance(value, str)


if __name__ == "__main__":
    unittest.main()
