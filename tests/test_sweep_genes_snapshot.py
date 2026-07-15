import json
import importlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SWEEP_ROOT = PROJECT_ROOT / "scripts" / "sweep-2.1.3.0"

loaded_bio_module = sys.modules.get("Bio")
loaded_seqio_module = sys.modules.get("Bio.SeqIO")
if (
    loaded_bio_module is not None
    and not hasattr(loaded_bio_module, "__file__")
) or (
    loaded_seqio_module is not None
    and not hasattr(loaded_seqio_module, "parse")
):
    for module_name in ("Bio.SeqIO", "Bio.Entrez", "Bio"):
        sys.modules.pop(module_name, None)

real_seqio = importlib.import_module("Bio.SeqIO")

from src.pago_pipeline.storage import sha256_of_file
import src.pago_pipeline.sweep_genes_snapshot as sweep_genes_snapshot_module
from src.pago_pipeline.sweep_genes_snapshot import (
    load_latest_sweep_genes_snapshot,
    resolve_sweep_genes_snapshot,
    save_sweep_genes_snapshot,
)


def _write_minimal_fasta_latest_snapshot(directory: Path) -> Path:
    snapshot_directory = directory / "fasta_snapshot" / "latest"
    snapshot_directory.mkdir(parents=True, exist_ok=True)

    fasta_file_path = snapshot_directory / "protein_sequences.fasta"
    fasta_file_path.write_text(
        (
            ">protein_uid=1|accession=A|length=6|organism=Specimen_1 Example one\n"
            "ACDEFG\n"
            ">protein_uid=2|accession=B|length=6|organism=Specimen_2 Example two\n"
            "ACDFGH\n"
        ),
        encoding="utf-8",
    )

    manifest_payload = {
        "artifact_type": "ncbi_protein_fasta_snapshot",
        "fasta_file_name": fasta_file_path.name,
        "fasta_file_sha256": sha256_of_file(input_file_path=fasta_file_path),
        "immutable_snapshot_relative_path": "snapshots/example_fasta_snapshot",
        "immutable_snapshot_directory_name": "example_fasta_snapshot",
        "search_query": "pAgo[Protein Name]",
        "translated_query": "pAgo[Protein Name]",
    }
    manifest_file_path = snapshot_directory / "manifest.json"
    manifest_file_path.write_text(
        json.dumps(manifest_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    return snapshot_directory


def _write_minimal_sweep_genes_snapshot(
    *,
    snapshot_directory: Path,
    source_fasta_snapshot_directory: Path,
    scripts_sweep_root_directory: Path,
    source_fasta_manifest_override: dict[str, object] | None = None,
) -> Path:
    snapshot_directory.mkdir(parents=True, exist_ok=True)
    source_manifest_file_path = source_fasta_snapshot_directory / "manifest.json"
    source_manifest_payload = json.loads(
        source_manifest_file_path.read_text(encoding="utf-8")
    )
    if source_fasta_manifest_override is not None:
        source_manifest_payload.update(source_fasta_manifest_override)

    embeddings_file_path = snapshot_directory / "sweep_genes_embeddings_2D.npy"
    sequence_metadata_file_path = snapshot_directory / "sequence_metadata.csv"
    profiling_log_file_path = snapshot_directory / "profiling_log.csv"

    np.save(embeddings_file_path, np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32))
    sequence_metadata_file_path.write_text(
        "sequence_index,record_id,description,sequence_length\n"
        "0,protein_uid=1,record one,6\n"
        "1,protein_uid=2,record two,6\n",
        encoding="utf-8",
    )
    profiling_log_file_path.write_text(
        "work_unit_kind,elapsed_seconds\nsnapshot_total,0.1\n",
        encoding="utf-8",
    )

    manifest_payload = {
        "artifact_type": "sweep_genes_embedding_snapshot",
        "manifest_file_name": "manifest.json",
        "embeddings_file_name": embeddings_file_path.name,
        "sequence_metadata_file_name": sequence_metadata_file_path.name,
        "profiling_log_file_name": profiling_log_file_path.name,
        "sequence_count": 2,
        "embedding_dimension": 2,
        "projected_dimensions_per_mask": 2,
        "mask_count": 1,
        "mask_length": 3,
        "masks": [[1, 0, 1]],
        "composition": "binary",
        "projection": True,
        "random_seed": 7,
        "chunk_size": 1,
        "n_jobs": 1,
        "dtype": "float32",
        "immutable_snapshot_directory_name": snapshot_directory.name,
        "immutable_snapshot_relative_path": f"snapshots/{snapshot_directory.name}",
        "sweep_package_root_directory": str(scripts_sweep_root_directory),
        "source_fasta_snapshot_relative_path": source_manifest_payload.get(
            "immutable_snapshot_relative_path"
        ),
        "source_fasta_snapshot_directory_name": source_manifest_payload.get(
            "immutable_snapshot_directory_name"
        ),
        "source_fasta_file_name": source_manifest_payload.get("fasta_file_name"),
        "source_fasta_file_sha256": source_manifest_payload.get("fasta_file_sha256"),
        "source_fasta_manifest_sha256": sha256_of_file(
            input_file_path=source_manifest_file_path
        ),
        "embeddings_file_sha256": sha256_of_file(input_file_path=embeddings_file_path),
        "sequence_metadata_file_sha256": sha256_of_file(
            input_file_path=sequence_metadata_file_path
        ),
        "profiling_log_file_sha256": sha256_of_file(
            input_file_path=profiling_log_file_path
        ),
    }
    manifest_file_path = snapshot_directory / "manifest.json"
    manifest_file_path.write_text(
        json.dumps(manifest_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    return snapshot_directory


class ResolveSweepGenesSnapshotTests(unittest.TestCase):
    def test_latest_availability_requires_sweep_artifact_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_fasta_snapshot_directory = _write_minimal_fasta_latest_snapshot(
                temporary_directory
            )
            sweep_snapshot_root_directory = temporary_directory / "sweep_genes_snapshot"
            latest_directory = sweep_snapshot_root_directory / "latest"

            _write_minimal_sweep_genes_snapshot(
                snapshot_directory=latest_directory,
                source_fasta_snapshot_directory=source_fasta_snapshot_directory,
                scripts_sweep_root_directory=temporary_directory / "unused_sweep",
            )
            manifest_file_path = latest_directory / "manifest.json"
            manifest_payload = json.loads(manifest_file_path.read_text(encoding="utf-8"))
            manifest_payload.pop("embeddings_file_sha256")
            manifest_file_path.write_text(
                json.dumps(manifest_payload, indent=2) + "\n",
                encoding="utf-8",
            )

            self.assertFalse(
                sweep_genes_snapshot_module.latest_sweep_genes_snapshot_is_available(
                    snapshot_root_directory=sweep_snapshot_root_directory
                )
            )

    def test_manifest_compatibility_rejects_different_sweep_package_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_fasta_snapshot_directory = _write_minimal_fasta_latest_snapshot(
                temporary_directory
            )
            snapshot_directory = (
                temporary_directory / "sweep_genes_snapshot" / "snapshots" / "candidate"
            )
            saved_sweep_root_directory = temporary_directory / "saved_sweep"
            requested_sweep_root_directory = temporary_directory / "requested_sweep"

            _write_minimal_sweep_genes_snapshot(
                snapshot_directory=snapshot_directory,
                source_fasta_snapshot_directory=source_fasta_snapshot_directory,
                scripts_sweep_root_directory=saved_sweep_root_directory,
            )
            source_manifest_file_path = source_fasta_snapshot_directory / "manifest.json"
            manifest_payload = json.loads(
                (snapshot_directory / "manifest.json").read_text(encoding="utf-8")
            )
            source_fasta_snapshot_payload = {
                "manifest": json.loads(
                    source_manifest_file_path.read_text(encoding="utf-8")
                ),
                "manifest_file_path": source_manifest_file_path,
            }

            self.assertFalse(
                sweep_genes_snapshot_module._sweep_genes_snapshot_manifest_matches_request(
                    manifest_payload=manifest_payload,
                    source_fasta_snapshot_payload=source_fasta_snapshot_payload,
                    scripts_sweep_root_directory=requested_sweep_root_directory,
                    masks=[[1, 0, 1]],
                    projected_dimensions_per_mask=2,
                    composition="binary",
                    projection=True,
                    random_seed=7,
                    chunk_size=1,
                    n_jobs=1,
                )
            )

    def test_reuse_latest_rejects_changed_source_fasta_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_fasta_snapshot_directory = _write_minimal_fasta_latest_snapshot(
                temporary_directory
            )
            sweep_snapshot_root_directory = temporary_directory / "sweep_genes_snapshot"
            _write_minimal_sweep_genes_snapshot(
                snapshot_directory=sweep_snapshot_root_directory / "latest",
                source_fasta_snapshot_directory=source_fasta_snapshot_directory,
                scripts_sweep_root_directory=temporary_directory / "unused_sweep",
            )

            source_manifest_file_path = source_fasta_snapshot_directory / "manifest.json"
            source_manifest_payload = json.loads(
                source_manifest_file_path.read_text(encoding="utf-8")
            )
            source_manifest_payload["immutable_snapshot_directory_name"] = (
                "replacement_fasta_snapshot"
            )
            source_manifest_payload["immutable_snapshot_relative_path"] = (
                "snapshots/replacement_fasta_snapshot"
            )
            source_manifest_file_path.write_text(
                json.dumps(source_manifest_payload, indent=2) + "\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(FileNotFoundError, "current source FASTA"):
                resolve_sweep_genes_snapshot(
                    snapshot_mode="reuse_latest",
                    snapshot_root_directory=sweep_snapshot_root_directory,
                    source_fasta_snapshot_root_directory=(
                        temporary_directory / "fasta_snapshot"
                    ),
                    scripts_sweep_root_directory=temporary_directory / "unused_sweep",
                )

    def test_reuse_latest_or_create_republishes_matching_immutable_without_loading_stale_latest_array(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_fasta_snapshot_directory = _write_minimal_fasta_latest_snapshot(
                temporary_directory
            )
            sweep_snapshot_root_directory = temporary_directory / "sweep_genes_snapshot"
            stale_latest_directory = sweep_snapshot_root_directory / "latest"
            matching_snapshot_directory = (
                sweep_snapshot_root_directory / "snapshots" / "matching_snapshot"
            )

            _write_minimal_sweep_genes_snapshot(
                snapshot_directory=stale_latest_directory,
                source_fasta_snapshot_directory=source_fasta_snapshot_directory,
                scripts_sweep_root_directory=temporary_directory / "unused_sweep",
                source_fasta_manifest_override={"fasta_file_sha256": "stale"},
            )
            _write_minimal_sweep_genes_snapshot(
                snapshot_directory=matching_snapshot_directory,
                source_fasta_snapshot_directory=source_fasta_snapshot_directory,
                scripts_sweep_root_directory=temporary_directory / "unused_sweep",
            )

            real_numpy_load = np.load
            with patch(
                "src.pago_pipeline.sweep_genes_snapshot.np.load",
                wraps=real_numpy_load,
            ) as numpy_load_mock:
                resolved_snapshot = resolve_sweep_genes_snapshot(
                    snapshot_mode="reuse_latest_or_create",
                    snapshot_root_directory=sweep_snapshot_root_directory,
                    source_fasta_snapshot_root_directory=(
                        temporary_directory / "fasta_snapshot"
                    ),
                    scripts_sweep_root_directory=temporary_directory / "unused_sweep",
                    masks=[[1, 0, 1]],
                    projected_dimensions_per_mask=2,
                    random_seed=7,
                    chunk_size=1,
                    n_jobs=1,
                )

            try:
                self.assertEqual(
                    Path(resolved_snapshot["snapshot_directory"]),
                    matching_snapshot_directory,
                )
                self.assertEqual(resolved_snapshot["embeddings"].shape, (2, 2))
                loaded_paths = [
                    Path(call.args[0]).resolve()
                    for call in numpy_load_mock.call_args_list
                    if call.args
                ]
                self.assertNotIn(
                    (stale_latest_directory / "sweep_genes_embeddings_2D.npy").resolve(),
                    loaded_paths,
                )

                latest_manifest_payload = json.loads(
                    (stale_latest_directory / "manifest.json").read_text(
                        encoding="utf-8"
                    )
                )
                self.assertEqual(
                    latest_manifest_payload["immutable_snapshot_directory_name"],
                    "matching_snapshot",
                )
            finally:
                resolved_snapshot["embeddings"]._mmap.close()


@unittest.skipUnless(
    SWEEP_ROOT.exists(),
    "SWeeP 2.1.3.0 is an external local dependency and is not committed.",
)
class SaveSweepGenesSnapshotTests(unittest.TestCase):
    def _write_minimal_fasta_snapshot(self, directory: Path) -> Path:
        snapshot_directory = directory / "fasta_snapshot" / "latest"
        snapshot_directory.mkdir(parents=True, exist_ok=True)

        fasta_file_path = snapshot_directory / "protein_sequences.fasta"
        fasta_file_path.write_text(
            (
                ">protein_uid=1|accession=A|length=6|organism=Specimen_1 Example one\n"
                "ACDEFG\n"
                ">protein_uid=2|accession=B|length=6|organism=Specimen_2 Example two\n"
                "ACDFGH\n"
            ),
            encoding="utf-8",
        )

        manifest_payload = {
            "artifact_type": "ncbi_protein_fasta_snapshot",
            "fasta_file_name": fasta_file_path.name,
            "fasta_file_sha256": sha256_of_file(input_file_path=fasta_file_path),
            "immutable_snapshot_relative_path": "snapshots/example_fasta_snapshot",
            "immutable_snapshot_directory_name": "example_fasta_snapshot",
            "search_query": "pAgo[Protein Name]",
            "translated_query": "pAgo[Protein Name]",
        }
        manifest_file_path = snapshot_directory / "manifest.json"
        manifest_file_path.write_text(
            json.dumps(manifest_payload, indent=2) + "\n",
            encoding="utf-8",
        )
        return snapshot_directory

    def test_saves_and_loads_latest_sweep_genes_snapshot(self) -> None:
        sweep_genes_snapshot_module.SeqIO = real_seqio

        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_fasta_snapshot_directory = self._write_minimal_fasta_snapshot(
                temporary_directory
            )
            sweep_snapshot_root_directory = temporary_directory / "sweep_genes_snapshot"

            save_sweep_genes_snapshot(
                snapshot_root_directory=sweep_snapshot_root_directory,
                source_fasta_snapshot_directory=source_fasta_snapshot_directory,
                scripts_sweep_root_directory=SWEEP_ROOT,
                masks=[[1, 0, 1, 0, 1]],
                projected_dimensions_per_mask=4,
                random_seed=7,
                chunk_size=1,
                n_jobs=1,
            )

            loaded_snapshot = load_latest_sweep_genes_snapshot(
                snapshot_root_directory=sweep_snapshot_root_directory
            )

            self.assertEqual(loaded_snapshot["embeddings"].shape, (2, 4))
            self.assertEqual(len(loaded_snapshot["sequence_metadata"]), 2)
            self.assertTrue(loaded_snapshot["profiling_log_file_path"].exists())
            self.assertEqual(
                loaded_snapshot["manifest"]["artifact_type"],
                "sweep_genes_embedding_snapshot",
            )
            self.assertEqual(loaded_snapshot["manifest"]["embedding_dimension"], 4)
            self.assertEqual(loaded_snapshot["manifest"]["chunk_size"], 1)
            self.assertEqual(loaded_snapshot["manifest"]["n_jobs"], 1)
            self.assertEqual(
                loaded_snapshot["manifest"]["profiling_log_file_name"],
                "profiling_log.csv",
            )
            self.assertIsInstance(loaded_snapshot["embeddings"], np.memmap)
            self.assertTrue(np.isfinite(loaded_snapshot["embeddings"]).all())
            loaded_snapshot["embeddings"]._mmap.close()


if __name__ == "__main__":
    unittest.main()
