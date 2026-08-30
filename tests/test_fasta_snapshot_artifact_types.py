from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.ncbi_fasta_snapshot import (
    DEFAULT_FASTA_SNAPSHOT_ARTIFACT_TYPES,
    load_fasta_snapshot_by_directory,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic
from src.pago_pipeline.sweep_genes_snapshot import (
    SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES,
)


def _write_fasta_snapshot_dir(
    *,
    directory: Path,
    artifact_type: str,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    fasta_file_path = directory / "protein_sequences.fasta"
    fasta_file_path.write_text(">protein_uid=1\nMKT\n", encoding="utf-8")
    write_json_atomic(
        payload={
            "artifact_type": artifact_type,
            "fasta_file_name": "protein_sequences.fasta",
            "fasta_file_sha256": sha256_of_file(input_file_path=fasta_file_path),
        },
        output_file_path=directory / "manifest.json",
    )


class FastaSnapshotArtifactTypeTests(unittest.TestCase):
    def test_default_still_rejects_non_ncbi_fasta_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            directory = Path(temporary_directory_name) / "snap"
            _write_fasta_snapshot_dir(
                directory=directory,
                artifact_type="derived_protein_fasta_snapshot",
            )
            with self.assertRaisesRegex(RuntimeError, "artifact_type mismatch"):
                load_fasta_snapshot_by_directory(snapshot_directory=directory)

            self.assertEqual(
                DEFAULT_FASTA_SNAPSHOT_ARTIFACT_TYPES,
                ("ncbi_protein_fasta_snapshot",),
            )

    def test_sweep_supported_types_accept_derived_fasta_but_not_arbitrary(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            base = Path(temporary_directory_name)

            derived_directory = base / "derived"
            _write_fasta_snapshot_dir(
                directory=derived_directory,
                artifact_type="derived_protein_fasta_snapshot",
            )
            payload = load_fasta_snapshot_by_directory(
                snapshot_directory=derived_directory,
                allowed_artifact_types=SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES,
            )
            self.assertEqual(
                payload["manifest"]["artifact_type"],
                "derived_protein_fasta_snapshot",
            )

            ncbi_directory = base / "ncbi"
            _write_fasta_snapshot_dir(
                directory=ncbi_directory,
                artifact_type="ncbi_protein_fasta_snapshot",
            )
            self.assertEqual(
                load_fasta_snapshot_by_directory(
                    snapshot_directory=ncbi_directory,
                    allowed_artifact_types=SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES,
                )["manifest"]["artifact_type"],
                "ncbi_protein_fasta_snapshot",
            )

            bogus_directory = base / "bogus"
            _write_fasta_snapshot_dir(
                directory=bogus_directory,
                artifact_type="some_other_snapshot",
            )
            with self.assertRaisesRegex(RuntimeError, "artifact_type mismatch"):
                load_fasta_snapshot_by_directory(
                    snapshot_directory=bogus_directory,
                    allowed_artifact_types=SUPPORTED_SOURCE_FASTA_ARTIFACT_TYPES,
                )


if __name__ == "__main__":
    unittest.main()
