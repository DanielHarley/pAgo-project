from __future__ import annotations

import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

if "Bio" not in sys.modules:
    fake_bio_module = types.ModuleType("Bio")
    fake_bio_module.__version__ = "test"
    fake_entrez_module = types.ModuleType("Bio.Entrez")
    fake_entrez_module.urlopen = lambda *args, **kwargs: None
    fake_seqio_module = types.ModuleType("Bio.SeqIO")
    fake_bio_module.Entrez = fake_entrez_module
    fake_bio_module.SeqIO = fake_seqio_module
    sys.modules["Bio"] = fake_bio_module
    sys.modules["Bio.Entrez"] = fake_entrez_module
    sys.modules["Bio.SeqIO"] = fake_seqio_module

from src.pago_pipeline.ncbi_snapshot import SnapshotMode
from src.pago_pipeline.pfam_hmm_bundle_snapshot import (
    ARTIFACT_TYPE,
    SNAPSHOT_FORMAT_VERSION,
    latest_pfam_hmm_bundle_snapshot_is_available,
    load_pfam_hmm_bundle_snapshot_by_directory,
    resolve_pfam_hmm_bundle_snapshot,
    save_pfam_hmm_bundle_snapshot,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic
from tests.pfam_hmm_test_support import (
    build_fake_pyhmmer,
    read_registry,
    write_registry,
    write_test_reference_bundle,
)


class PfamHmmBundleSnapshotTests(unittest.TestCase):
    def _resolve(self, temporary_directory: Path):
        registry_path, hmm_directory = write_test_reference_bundle(
            reference_directory=temporary_directory / "reference"
        )
        snapshot_root = temporary_directory / "snapshots_root"
        fake_pyhmmer = build_fake_pyhmmer()
        payload = resolve_pfam_hmm_bundle_snapshot(
            snapshot_mode=SnapshotMode.reuse_latest_or_create,
            snapshot_root_directory=snapshot_root,
            registry_file_path=registry_path,
            hmm_directory=hmm_directory,
            pyhmmer_module=fake_pyhmmer,
        )
        return payload, registry_path, hmm_directory, snapshot_root, fake_pyhmmer

    def test_resolve_persists_complete_provenance_and_pressed_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            payload, registry_path, _, snapshot_root, _ = self._resolve(
                temporary_directory
            )
            manifest = payload["manifest"]

            self.assertEqual(manifest["artifact_type"], ARTIFACT_TYPE)
            self.assertEqual(
                manifest["snapshot_format_version"], SNAPSHOT_FORMAT_VERSION
            )
            self.assertEqual(
                manifest["source_registry_sha256"],
                sha256_of_file(input_file_path=registry_path),
            )
            self.assertEqual(manifest["source_release"], "test-release-1")
            self.assertEqual(manifest["source_model_count"], 2)
            self.assertEqual(manifest["bundle_model_count"], 2)
            self.assertEqual(manifest["pressed_output_count"], 4)
            self.assertEqual(
                manifest["canonical_model_accessions"],
                ["PF02171.22", "PF16487.10"],
            )
            self.assertEqual(len(manifest["derivation_identity_sha256"]), 64)
            self.assertEqual(len(payload["pressed_file_paths"]), 4)
            self.assertTrue((snapshot_root / "latest" / "manifest.json").is_file())
            self.assertEqual(
                set(manifest["output_files"]),
                {
                    "bundle_hmm",
                    "pressed_h3f",
                    "pressed_h3i",
                    "pressed_h3m",
                    "pressed_h3p",
                },
            )

    def test_reuse_latest_or_create_does_not_rebuild_a_matching_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            _, registry_path, hmm_directory, snapshot_root, fake_pyhmmer = self._resolve(
                temporary_directory
            )
            with mock.patch(
                "src.pago_pipeline.pfam_hmm_bundle_snapshot"
                ".save_pfam_hmm_bundle_snapshot"
            ) as save_mock:
                resolve_pfam_hmm_bundle_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=snapshot_root,
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=fake_pyhmmer,
                )
                save_mock.assert_not_called()

    def test_latest_is_invalidated_by_reference_or_pyhmmer_identity_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            _, registry_path, hmm_directory, snapshot_root, fake_pyhmmer = self._resolve(
                temporary_directory
            )
            self.assertTrue(
                latest_pfam_hmm_bundle_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=fake_pyhmmer,
                )
            )
            self.assertFalse(
                latest_pfam_hmm_bundle_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=build_fake_pyhmmer(version="changed-version"),
                )
            )

            registry = read_registry(registry_path)
            registry["source_release"] = "test-release-2"
            write_registry(registry_path, registry)
            self.assertFalse(
                latest_pfam_hmm_bundle_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=fake_pyhmmer,
                )
            )

    def test_loader_rejects_tampered_bundle_even_when_other_outputs_are_intact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            payload, _, _, _, fake_pyhmmer = self._resolve(temporary_directory)
            bundle_file_path = Path(payload["bundle_file_path"])
            bundle_file_path.write_bytes(bundle_file_path.read_bytes() + b"tamper")

            with self.assertRaisesRegex(RuntimeError, "size mismatch|hash mismatch"):
                load_pfam_hmm_bundle_snapshot_by_directory(
                    snapshot_directory=payload["snapshot_directory"],
                    pyhmmer_module=fake_pyhmmer,
                )

    def test_loader_rejects_inconsistent_derivation_or_runtime_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            payload, _, _, _, fake_pyhmmer = self._resolve(temporary_directory)
            manifest_file_path = Path(payload["manifest_file_path"])
            manifest = dict(payload["manifest"])
            manifest["source_release"] = "tampered-release"
            write_json_atomic(
                payload=manifest,
                output_file_path=manifest_file_path,
            )
            with self.assertRaisesRegex(RuntimeError, "derivation identity"):
                load_pfam_hmm_bundle_snapshot_by_directory(
                    snapshot_directory=payload["snapshot_directory"],
                    pyhmmer_module=fake_pyhmmer,
                )

            write_json_atomic(
                payload=payload["manifest"],
                output_file_path=manifest_file_path,
            )
            with self.assertRaisesRegex(RuntimeError, "PyHMMER version mismatch"):
                load_pfam_hmm_bundle_snapshot_by_directory(
                    snapshot_directory=payload["snapshot_directory"],
                    pyhmmer_module=build_fake_pyhmmer(version="other-version"),
                )

    def test_snapshot_creation_requires_the_press_step(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, hmm_directory = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference"
            )
            snapshot_root = temporary_directory / "snapshots_root"
            with self.assertRaisesRegex(RuntimeError, "does not expose hmmer.hmmpress"):
                save_pfam_hmm_bundle_snapshot(
                    snapshot_root_directory=snapshot_root,
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=build_fake_pyhmmer(hmmpress_function=None),
                )
            snapshots_directory = snapshot_root / "snapshots"
            self.assertEqual(
                list(snapshots_directory.iterdir()) if snapshots_directory.exists() else [],
                [],
            )

    def test_failed_press_cleans_incomplete_immutable_snapshot_and_preserves_latest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, hmm_directory = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference"
            )
            snapshot_root = temporary_directory / "snapshots_root"
            latest_directory = snapshot_root / "latest"
            latest_directory.mkdir(parents=True)
            sentinel_file_path = latest_directory / "sentinel.txt"
            sentinel_file_path.write_text("previous latest", encoding="utf-8")

            def failing_hmmpress(hmms, output) -> None:
                list(hmms)
                Path(str(output) + ".h3f").write_bytes(b"partial")
                raise RuntimeError("simulated press failure")

            with self.assertRaisesRegex(RuntimeError, "simulated press failure"):
                save_pfam_hmm_bundle_snapshot(
                    snapshot_root_directory=snapshot_root,
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=build_fake_pyhmmer(
                        hmmpress_function=failing_hmmpress
                    ),
                )

            snapshots_directory = snapshot_root / "snapshots"
            self.assertEqual(
                list(snapshots_directory.iterdir()) if snapshots_directory.exists() else [],
                [],
            )
            self.assertEqual(
                sentinel_file_path.read_text(encoding="utf-8"),
                "previous latest",
            )


if __name__ == "__main__":
    unittest.main()
