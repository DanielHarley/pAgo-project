from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.apaz_hmm_build_snapshot import (
    load_apaz_hmm_build_snapshot_by_directory,
    save_apaz_hmm_build_snapshot,
)
from src.pago_pipeline.storage import write_json_atomic
from tests.apaz_hmm_test_support import (
    deterministic_fake_hmm_builder,
    read_json,
    write_apaz_reference_fixture,
)


def _save_fixture_snapshot(directory: Path):
    reference_directory = directory / "apaz_reference"
    snapshot_root = directory / "snapshot_root"
    write_apaz_reference_fixture(reference_directory=reference_directory)
    return save_apaz_hmm_build_snapshot(
        snapshot_root_directory=snapshot_root,
        reference_directory=reference_directory,
        hmm_builder=deterministic_fake_hmm_builder,
        builder_identity="fixture-builder/1",
    )


class ApazHmmBuildSnapshotTests(unittest.TestCase):
    def test_snapshot_round_trips_and_records_structural_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            result = _save_fixture_snapshot(Path(name))
            loaded = load_apaz_hmm_build_snapshot_by_directory(
                snapshot_directory=result.snapshot_directory
            )
            manifest = loaded["manifest"]
            self.assertEqual(manifest["required_partition"], "BUILD")
            self.assertTrue(manifest["build_protocol"])
            self.assertEqual(manifest["model_count"], 6)
            for entry in manifest["models"]:
                self.assertGreaterEqual(entry["match_state_count"], 1)
                self.assertGreaterEqual(entry["hmm_file_size_bytes"], 1)
            self.assertEqual(
                set(loaded["hmm_file_path_by_model_id"]),
                {"global", "Ia", "Ib", "IIa", "IIb", "III"},
            )

    def test_rejects_policy_metadata_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            result = _save_fixture_snapshot(Path(name))
            manifest = read_json(result.manifest_file_path)
            manifest["random_seed"] = 43
            write_json_atomic(
                payload=manifest, output_file_path=result.manifest_file_path
            )
            with self.assertRaisesRegex(RuntimeError, "build_policy"):
                load_apaz_hmm_build_snapshot_by_directory(
                    snapshot_directory=result.snapshot_directory
                )

    def test_rejects_a_tampered_output_hash(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            result = _save_fixture_snapshot(Path(name))
            manifest = read_json(result.manifest_file_path)
            manifest["output_files"]["Ia_hmm"]["sha256"] = "0" * 64
            for entry in manifest["models"]:
                if entry["model_id"] == "Ia":
                    entry["hmm_sha256"] = "0" * 64
            write_json_atomic(
                payload=manifest, output_file_path=result.manifest_file_path
            )
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                load_apaz_hmm_build_snapshot_by_directory(
                    snapshot_directory=result.snapshot_directory
                )

    def test_rejects_a_missing_model_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            result = _save_fixture_snapshot(Path(name))
            (result.snapshot_directory / "apaz_III.hmm").unlink()
            with self.assertRaises((FileNotFoundError, RuntimeError)):
                load_apaz_hmm_build_snapshot_by_directory(
                    snapshot_directory=result.snapshot_directory
                )

    def test_rejects_a_stripped_provenance_field(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            result = _save_fixture_snapshot(Path(name))
            manifest = read_json(result.manifest_file_path)
            del manifest["build_protocol"]
            write_json_atomic(
                payload=manifest, output_file_path=result.manifest_file_path
            )
            with self.assertRaisesRegex(RuntimeError, "build_protocol"):
                load_apaz_hmm_build_snapshot_by_directory(
                    snapshot_directory=result.snapshot_directory
                )


if __name__ == "__main__":
    unittest.main()
