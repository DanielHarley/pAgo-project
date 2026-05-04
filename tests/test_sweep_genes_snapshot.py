import json
import importlib
import sys
import tempfile
import unittest
from pathlib import Path

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
    save_sweep_genes_snapshot,
)


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
