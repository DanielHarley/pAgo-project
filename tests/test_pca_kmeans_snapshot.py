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

import src.pago_pipeline.pca_kmeans_snapshot as pca_kmeans_snapshot_module
from src.pago_pipeline.pca_kmeans_snapshot import (
    load_pca_kmeans_snapshot_by_directory,
    resolve_pca_kmeans_snapshot,
)
from src.pago_pipeline.storage import sha256_of_file


def _write_minimal_source_sweep_latest_snapshot(root_directory: Path) -> Path:
    snapshot_directory = root_directory / "latest"
    snapshot_directory.mkdir(parents=True, exist_ok=True)

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
        "immutable_snapshot_directory_name": "source_sweep_snapshot",
        "immutable_snapshot_relative_path": "snapshots/source_sweep_snapshot",
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


def _write_minimal_source_metadata_latest_snapshot(root_directory: Path) -> Path:
    snapshot_directory = root_directory / "latest"
    snapshot_directory.mkdir(parents=True, exist_ok=True)

    csv_file_path = snapshot_directory / "protein_metadata.csv"
    qc_report_file_path = snapshot_directory / "qc_report.json"
    csv_file_path.write_text("protein_uid\n1\n2\n", encoding="utf-8")
    qc_report_file_path.write_text("{}\n", encoding="utf-8")

    manifest_payload = {
        "artifact_type": "ncbi_protein_metadata_snapshot",
        "csv_file_name": csv_file_path.name,
        "qc_report_file_name": qc_report_file_path.name,
        "csv_file_sha256": sha256_of_file(input_file_path=csv_file_path),
        "qc_report_file_sha256": sha256_of_file(input_file_path=qc_report_file_path),
        "immutable_snapshot_directory_name": "source_metadata_snapshot",
        "immutable_snapshot_relative_path": "snapshots/source_metadata_snapshot",
    }
    manifest_file_path = snapshot_directory / "manifest.json"
    manifest_file_path.write_text(
        json.dumps(manifest_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    return snapshot_directory


def _write_minimal_pca_kmeans_snapshot(
    *,
    snapshot_directory: Path,
    source_sweep_snapshot_directory: Path,
    source_metadata_snapshot_directory: Path,
    source_sweep_manifest_override: dict[str, object] | None = None,
) -> Path:
    snapshot_directory.mkdir(parents=True, exist_ok=True)

    source_sweep_manifest_file_path = source_sweep_snapshot_directory / "manifest.json"
    source_sweep_manifest_payload = json.loads(
        source_sweep_manifest_file_path.read_text(encoding="utf-8")
    )
    if source_sweep_manifest_override is not None:
        source_sweep_manifest_payload.update(source_sweep_manifest_override)

    source_metadata_manifest_file_path = (
        source_metadata_snapshot_directory / "manifest.json"
    )
    source_metadata_manifest_payload = json.loads(
        source_metadata_manifest_file_path.read_text(encoding="utf-8")
    )

    pca_coordinates_file_path = snapshot_directory / "pca_coordinates_2D.npy"
    explained_variance_ratio_file_path = (
        snapshot_directory / "explained_variance_ratio_2D.npy"
    )
    cluster_assignments_file_path = snapshot_directory / "cluster_assignments.csv"
    stability_grid_file_path = snapshot_directory / "stability_grid.csv"
    profiling_log_file_path = snapshot_directory / "profiling_log.csv"
    alignment_report_file_path = snapshot_directory / "alignment_report.json"

    np.save(
        pca_coordinates_file_path,
        np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
    )
    np.save(
        explained_variance_ratio_file_path,
        np.asarray([0.6, 0.4], dtype=np.float64),
    )
    pd.DataFrame({"cluster_label": [0, 1], "pc1": [1.0, 3.0]}).to_csv(
        cluster_assignments_file_path,
        index=False,
    )
    pd.DataFrame({"pca_component_count": [2], "k": [2]}).to_csv(
        stability_grid_file_path,
        index=False,
    )
    pd.DataFrame({"work_unit_kind": ["snapshot_total"]}).to_csv(
        profiling_log_file_path,
        index=False,
    )
    alignment_report_file_path.write_text(
        json.dumps({"sequence_row_count": 2}, indent=2) + "\n",
        encoding="utf-8",
    )

    manifest_payload = {
        "artifact_type": "pca_kmeans_snapshot",
        "manifest_file_name": "manifest.json",
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
        "immutable_snapshot_directory_name": snapshot_directory.name,
        "immutable_snapshot_relative_path": f"snapshots/{snapshot_directory.name}",
        "source_sweep_snapshot_relative_path": source_sweep_manifest_payload.get(
            "immutable_snapshot_relative_path"
        ),
        "source_sweep_snapshot_directory_name": source_sweep_manifest_payload.get(
            "immutable_snapshot_directory_name"
        ),
        "source_sweep_manifest_sha256": sha256_of_file(
            input_file_path=source_sweep_manifest_file_path
        ),
        "source_metadata_snapshot_relative_path": source_metadata_manifest_payload.get(
            "immutable_snapshot_relative_path"
        ),
        "source_metadata_snapshot_directory_name": source_metadata_manifest_payload.get(
            "immutable_snapshot_directory_name"
        ),
        "source_metadata_manifest_sha256": sha256_of_file(
            input_file_path=source_metadata_manifest_file_path
        ),
        "parameter_search_configuration": {
            "pca_component_count_grid": [2],
            "kmeans_cluster_count_grid": [
                int(value)
                for value in pca_kmeans_snapshot_module.DEFAULT_KMEANS_CLUSTER_COUNT_GRID
            ],
            "kmeans_initialization_repeat_count": (
                pca_kmeans_snapshot_module.DEFAULT_KMEANS_INITIALIZATION_REPEAT_COUNT
            ),
            "subsample_repeat_count": (
                pca_kmeans_snapshot_module.DEFAULT_SUBSAMPLE_REPEAT_COUNT
            ),
            "subsample_fraction": pca_kmeans_snapshot_module.DEFAULT_SUBSAMPLE_FRACTION,
            "pca_svd_solver": pca_kmeans_snapshot_module.DEFAULT_PCA_SVD_SOLVER,
            "pca_random_state": pca_kmeans_snapshot_module.DEFAULT_PCA_RANDOM_STATE,
            "kmeans_n_init": pca_kmeans_snapshot_module.DEFAULT_KMEANS_N_INIT,
            "silhouette_sample_size": (
                pca_kmeans_snapshot_module.DEFAULT_SILHOUETTE_SAMPLE_SIZE
            ),
            "silhouette_random_state": (
                pca_kmeans_snapshot_module.DEFAULT_SILHOUETTE_RANDOM_STATE
            ),
            "subsample_random_state": (
                pca_kmeans_snapshot_module.DEFAULT_SUBSAMPLE_RANDOM_STATE
            ),
            "minimum_acceptable_init_ari_min": (
                pca_kmeans_snapshot_module.DEFAULT_MINIMUM_ACCEPTABLE_INIT_ARI_MIN
            ),
            "minimum_acceptable_subsample_ari_min": (
                pca_kmeans_snapshot_module.DEFAULT_MINIMUM_ACCEPTABLE_SUBSAMPLE_ARI_MIN
            ),
            "export_projection_component_count": (
                pca_kmeans_snapshot_module.DEFAULT_EXPORT_PROJECTION_COMPONENT_COUNT
            ),
        },
    }
    manifest_file_path = snapshot_directory / "manifest.json"
    manifest_file_path.write_text(
        json.dumps(manifest_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    return snapshot_directory


class ResolvePcaKMeansSnapshotTests(unittest.TestCase):
    def test_latest_availability_requires_pca_kmeans_artifact_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_sweep_snapshot_directory = _write_minimal_source_sweep_latest_snapshot(
                temporary_directory / "sweep_genes"
            )
            source_metadata_snapshot_directory = (
                _write_minimal_source_metadata_latest_snapshot(
                    temporary_directory / "metadata"
                )
            )
            pca_snapshot_root_directory = temporary_directory / "pca_kmeans"
            latest_directory = pca_snapshot_root_directory / "latest"

            _write_minimal_pca_kmeans_snapshot(
                snapshot_directory=latest_directory,
                source_sweep_snapshot_directory=source_sweep_snapshot_directory,
                source_metadata_snapshot_directory=source_metadata_snapshot_directory,
            )
            manifest_file_path = latest_directory / "manifest.json"
            manifest_payload = json.loads(manifest_file_path.read_text(encoding="utf-8"))
            manifest_payload.pop("pca_coordinates_file_sha256")
            manifest_file_path.write_text(
                json.dumps(manifest_payload, indent=2) + "\n",
                encoding="utf-8",
            )

            self.assertFalse(
                pca_kmeans_snapshot_module.latest_pca_kmeans_snapshot_is_available(
                    snapshot_root_directory=pca_snapshot_root_directory
                )
            )

    def test_reuse_latest_or_create_republishes_matching_immutable_without_loading_stale_latest_arrays(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            source_sweep_snapshot_root_directory = temporary_directory / "sweep_genes"
            source_metadata_snapshot_root_directory = temporary_directory / "metadata"
            source_sweep_snapshot_directory = (
                _write_minimal_source_sweep_latest_snapshot(
                    source_sweep_snapshot_root_directory
                )
            )
            source_metadata_snapshot_directory = (
                _write_minimal_source_metadata_latest_snapshot(
                    source_metadata_snapshot_root_directory
                )
            )

            pca_snapshot_root_directory = temporary_directory / "pca_kmeans"
            stale_latest_directory = pca_snapshot_root_directory / "latest"
            matching_snapshot_directory = (
                pca_snapshot_root_directory / "snapshots" / "matching_snapshot"
            )

            _write_minimal_pca_kmeans_snapshot(
                snapshot_directory=stale_latest_directory,
                source_sweep_snapshot_directory=source_sweep_snapshot_directory,
                source_metadata_snapshot_directory=source_metadata_snapshot_directory,
                source_sweep_manifest_override={
                    "immutable_snapshot_directory_name": "stale_source"
                },
            )
            _write_minimal_pca_kmeans_snapshot(
                snapshot_directory=matching_snapshot_directory,
                source_sweep_snapshot_directory=source_sweep_snapshot_directory,
                source_metadata_snapshot_directory=source_metadata_snapshot_directory,
            )

            real_numpy_load = np.load
            with patch(
                "src.pago_pipeline.pca_kmeans_snapshot.np.load",
                wraps=real_numpy_load,
            ) as numpy_load_mock:
                resolved_snapshot = resolve_pca_kmeans_snapshot(
                    snapshot_mode="reuse_latest_or_create",
                    snapshot_root_directory=pca_snapshot_root_directory,
                    source_sweep_snapshot_root_directory=(
                        source_sweep_snapshot_root_directory
                    ),
                    source_metadata_snapshot_root_directory=(
                        source_metadata_snapshot_root_directory
                    ),
                    pca_component_count_grid=(2,),
                )

            try:
                self.assertEqual(
                    Path(resolved_snapshot["snapshot_directory"]),
                    matching_snapshot_directory,
                )
                self.assertEqual(resolved_snapshot["pca_coordinates"].shape, (2, 2))
                loaded_paths = [
                    Path(call.args[0]).resolve()
                    for call in numpy_load_mock.call_args_list
                    if call.args
                ]
                self.assertNotIn(
                    (stale_latest_directory / "pca_coordinates_2D.npy").resolve(),
                    loaded_paths,
                )
                self.assertNotIn(
                    (
                        stale_latest_directory / "explained_variance_ratio_2D.npy"
                    ).resolve(),
                    loaded_paths,
                )
                self.assertNotIn(
                    (
                        source_sweep_snapshot_directory
                        / "sweep_genes_embeddings_2D.npy"
                    ).resolve(),
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
                resolved_snapshot["pca_coordinates"]._mmap.close()
                resolved_snapshot["explained_variance_ratio"]._mmap.close()


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
