from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.apaz_hmm_build_snapshot import (
    latest_apaz_hmm_build_snapshot_is_available,
    load_apaz_hmm_build_snapshot_by_directory,
    save_apaz_hmm_build_snapshot,
)
from src.pago_pipeline.storage import write_json_atomic
from tests.apaz_hmm_test_support import (
    deterministic_fake_hmm_builder,
    read_json,
    refresh_apaz_locked_seed,
    refresh_apaz_partitions_lock_hash,
    write_apaz_reference_fixture,
)

_FIXTURE_BUILDER_IDENTITY = "fixture-builder/1"


def _save_fixture_snapshot(directory: Path):
    reference_directory = directory / "apaz_reference"
    snapshot_root = directory / "snapshot_root"
    write_apaz_reference_fixture(reference_directory=reference_directory)
    return save_apaz_hmm_build_snapshot(
        snapshot_root_directory=snapshot_root,
        reference_directory=reference_directory,
        hmm_builder=deterministic_fake_hmm_builder,
        builder_identity=_FIXTURE_BUILDER_IDENTITY,
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


class ApazHmmBuildSnapshotReuseInvalidationTests(unittest.TestCase):
    """A previously compatible ``latest`` snapshot must stop being reusable
    the moment any of its declared source hashes no longer matches the
    current reference directory: the seed lock itself, the partition table,
    or an individual seed alignment.
    """

    def _write_fixture_and_publish_latest(self, directory: Path) -> Path:
        reference_directory = directory / "apaz_reference"
        snapshot_root = directory / "snapshot_root"
        write_apaz_reference_fixture(reference_directory=reference_directory)
        save_apaz_hmm_build_snapshot(
            snapshot_root_directory=snapshot_root,
            reference_directory=reference_directory,
            hmm_builder=deterministic_fake_hmm_builder,
            builder_identity=_FIXTURE_BUILDER_IDENTITY,
        )
        self.assertTrue(
            latest_apaz_hmm_build_snapshot_is_available(
                snapshot_root_directory=snapshot_root,
                reference_directory=reference_directory,
                builder_identity=_FIXTURE_BUILDER_IDENTITY,
            ),
            "Precondition failed: the freshly published snapshot must start "
            "out compatible with its own reference directory.",
        )
        return snapshot_root

    def test_rejects_reuse_after_the_seed_lock_hash_changes(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            snapshot_root = self._write_fixture_and_publish_latest(directory)
            reference_directory = directory / "apaz_reference"
            lock_path = reference_directory / "seeds_lock.json"
            # Change only the lock file's own bytes (an additive, unchecked
            # field), so seed_lock_sha256 changes while every declared
            # partition/seed hash it protects remains internally consistent.
            payload = read_json(lock_path)
            payload["fixture_revision_marker"] = "changed"
            write_json_atomic(payload=payload, output_file_path=lock_path)

            self.assertFalse(
                latest_apaz_hmm_build_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    reference_directory=reference_directory,
                    builder_identity=_FIXTURE_BUILDER_IDENTITY,
                )
            )

    def test_rejects_reuse_after_the_partitions_hash_changes(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            snapshot_root = self._write_fixture_and_publish_latest(directory)
            reference_directory = directory / "apaz_reference"
            partitions_path = reference_directory / "apaz_partitions.csv"
            # Reorder the rows: partition membership is accession-keyed, so
            # this changes the file's bytes (and therefore its SHA-256)
            # without changing any accession's partition assignment.
            with partitions_path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            reordered_rows = [rows[-1], *rows[:-1]]
            with partitions_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(reordered_rows)
            refresh_apaz_partitions_lock_hash(reference_directory=reference_directory)

            self.assertFalse(
                latest_apaz_hmm_build_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    reference_directory=reference_directory,
                    builder_identity=_FIXTURE_BUILDER_IDENTITY,
                )
            )

    def test_rejects_reuse_after_an_individual_seed_hash_changes(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            snapshot_root = self._write_fixture_and_publish_latest(directory)
            reference_directory = directory / "apaz_reference"
            seed_path = reference_directory / "apaz_Ia.sto"
            # A trailing blank line is skipped unconditionally by the
            # Stockholm reader regardless of the // terminator, so this
            # changes apaz_Ia.hmm's declared seed_sha256 without altering
            # its parsed accessions, sequences, or alignment length.
            with seed_path.open("a", encoding="utf-8") as handle:
                handle.write("\n")
            refresh_apaz_locked_seed(
                lock_file_path=reference_directory / "seeds_lock.json",
                seed_file_path=seed_path,
            )

            self.assertFalse(
                latest_apaz_hmm_build_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    reference_directory=reference_directory,
                    builder_identity=_FIXTURE_BUILDER_IDENTITY,
                )
            )


if __name__ == "__main__":
    unittest.main()
