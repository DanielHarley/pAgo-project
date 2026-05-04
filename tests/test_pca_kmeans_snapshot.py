import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

if "Bio" not in sys.modules:
    fake_bio_module = types.ModuleType("Bio")
    fake_bio_module.__version__ = "test"
    fake_entrez_module = types.ModuleType("Bio.Entrez")
    fake_seqio_module = types.ModuleType("Bio.SeqIO")
    fake_bio_module.Entrez = fake_entrez_module
    fake_bio_module.SeqIO = fake_seqio_module
    sys.modules["Bio"] = fake_bio_module
    sys.modules["Bio.Entrez"] = fake_entrez_module
    sys.modules["Bio.SeqIO"] = fake_seqio_module

from src.pago_pipeline.pca_kmeans_snapshot import load_pca_kmeans_snapshot_by_directory
from src.pago_pipeline.storage import sha256_of_file


class LoadPcaKMeansSnapshotTests(unittest.TestCase):
    def _write_minimal_snapshot(self, snapshot_directory: Path) -> Path:
        snapshot_directory.mkdir(parents=True, exist_ok=True)

        pca_coordinates_file_path = snapshot_directory / "pca_coordinates_3D.npy"
        explained_variance_ratio_file_path = (
            snapshot_directory / "explained_variance_ratio_3D.npy"
        )
        cluster_assignments_file_path = snapshot_directory / "cluster_assignments.csv"
        stability_grid_file_path = snapshot_directory / "stability_grid.csv"
        profiling_log_file_path = snapshot_directory / "profiling_log.csv"
        alignment_report_file_path = snapshot_directory / "alignment_report.json"

        np.save(
            pca_coordinates_file_path,
            np.asarray([[1.0, 2.0, 3.0]], dtype=np.float32),
        )
        np.save(
            explained_variance_ratio_file_path,
            np.asarray([0.6, 0.3, 0.1], dtype=np.float32),
        )
        pd.DataFrame(
            {
                "cluster_label": [0, 1],
                "pc1": [1.0, 2.0],
                "pc2": [3.0, 4.0],
                "mixed_metadata": ["alpha", 42],
            }
        ).to_csv(cluster_assignments_file_path, index=False)
        pd.DataFrame(
            {
                "candidate_pca_component_count": [3],
                "candidate_cluster_count_k": [2],
                "sampled_silhouette_value": [0.5],
            }
        ).to_csv(stability_grid_file_path, index=False)
        pd.DataFrame(
            {
                "work_unit_kind": ["snapshot_total"],
                "elapsed_seconds": [0.1],
            }
        ).to_csv(profiling_log_file_path, index=False)
        alignment_report_file_path.write_text(
            json.dumps({"row_count": 2}, indent=2) + "\n",
            encoding="utf-8",
        )

        manifest_payload = {
            "artifact_type": "pca_kmeans_snapshot",
            "pca_coordinates_file_name": pca_coordinates_file_path.name,
            "explained_variance_ratio_file_name": explained_variance_ratio_file_path.name,
            "cluster_assignments_file_name": cluster_assignments_file_path.name,
            "stability_grid_file_name": stability_grid_file_path.name,
            "profiling_log_file_name": profiling_log_file_path.name,
            "alignment_report_file_name": alignment_report_file_path.name,
            "pca_coordinates_file_sha256": sha256_of_file(
                input_file_path=pca_coordinates_file_path
            ),
            "explained_variance_ratio_file_sha256": sha256_of_file(
                input_file_path=explained_variance_ratio_file_path
            ),
            "cluster_assignments_file_sha256": sha256_of_file(
                input_file_path=cluster_assignments_file_path
            ),
            "stability_grid_file_sha256": sha256_of_file(
                input_file_path=stability_grid_file_path
            ),
            "profiling_log_file_sha256": sha256_of_file(
                input_file_path=profiling_log_file_path
            ),
            "alignment_report_file_sha256": sha256_of_file(
                input_file_path=alignment_report_file_path
            ),
        }
        manifest_file_path = snapshot_directory / "manifest.json"
        manifest_file_path.write_text(
            json.dumps(manifest_payload, indent=2) + "\n",
            encoding="utf-8",
        )
        return cluster_assignments_file_path

    def test_load_snapshot_reads_cluster_assignments_with_low_memory_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name) / "snapshot"
            cluster_assignments_file_path = self._write_minimal_snapshot(
                snapshot_directory
            )

            with patch(
                "src.pago_pipeline.pca_kmeans_snapshot.pd.read_csv",
                wraps=pd.read_csv,
            ) as read_csv_mock:
                loaded_snapshot = load_pca_kmeans_snapshot_by_directory(
                    snapshot_directory=snapshot_directory
                )

            try:
                self.assertEqual(len(loaded_snapshot["cluster_assignments"]), 2)
                cluster_assignment_calls = [
                    call
                    for call in read_csv_mock.call_args_list
                    if call.args and Path(call.args[0]) == cluster_assignments_file_path
                ]
                self.assertEqual(len(cluster_assignment_calls), 1)
                self.assertEqual(
                    cluster_assignment_calls[0].kwargs.get("low_memory"),
                    False,
                )
            finally:
                loaded_snapshot["pca_coordinates"]._mmap.close()
                loaded_snapshot["explained_variance_ratio"]._mmap.close()


if __name__ == "__main__":
    unittest.main()
