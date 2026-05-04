import sys
import types
import unittest

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

import src.pago_pipeline.pca_plot.pca_plot_dataset as dataset_module
import src.pago_pipeline.pca_plot.pca_plot_render as render_module
import src.pago_pipeline.pca_plot.pca_plot_snapshot as snapshot_module


class BuildPcaPlotDataframeTests(unittest.TestCase):
    def test_auto_mode_builds_2d_plot_when_pc3_is_missing(self) -> None:
        cluster_assignments_dataframe = pd.DataFrame(
            {
                "cluster_label": [0, 1],
                "pc1": [1.0, 2.0],
                "pc2": [3.0, 4.0],
                "description": ["pAgo protein", "helicase"],
                "gbseq__organism": ["Species A", "Species B"],
            }
        )

        plot_dataframe = dataset_module.build_pca_plot_dataframe(
            cluster_assignments_dataframe=cluster_assignments_dataframe,
            plot_dimension_mode="auto",
        )

        self.assertEqual(plot_dataframe["plot_dimension_count"].unique().tolist(), [2])
        self.assertNotIn("plot_pc3", plot_dataframe.columns)
        self.assertEqual(plot_dataframe["plot_dimension_label"].unique().tolist(), ["2d"])
        self.assertEqual(plot_dataframe["pago_label"].tolist(), ["pAgo", "Other"])

    def test_auto_mode_builds_3d_plot_when_pc3_is_available(self) -> None:
        cluster_assignments_dataframe = pd.DataFrame(
            {
                "cluster_label": [0, 1],
                "pc1": [1.0, 2.0],
                "pc2": [3.0, 4.0],
                "pc3": [5.0, 6.0],
                "description": ["protein a", "protein b"],
                "gbseq__organism": ["Species A", "Species B"],
            }
        )

        plot_dataframe = dataset_module.build_pca_plot_dataframe(
            cluster_assignments_dataframe=cluster_assignments_dataframe,
            plot_dimension_mode="auto",
        )

        self.assertEqual(plot_dataframe["plot_dimension_count"].unique().tolist(), [3])
        self.assertIn("plot_pc3", plot_dataframe.columns)
        self.assertEqual(plot_dataframe["plot_dimension_label"].unique().tolist(), ["3d"])

    def test_explicit_2d_mode_ignores_available_pc3(self) -> None:
        cluster_assignments_dataframe = pd.DataFrame(
            {
                "cluster_label": [0],
                "pc1": [1.0],
                "pc2": [2.0],
                "pc3": [3.0],
            }
        )

        plot_dataframe = dataset_module.build_pca_plot_dataframe(
            cluster_assignments_dataframe=cluster_assignments_dataframe,
            plot_dimension_mode="2d",
        )

        self.assertEqual(plot_dataframe["plot_dimension_count"].unique().tolist(), [2])
        self.assertNotIn("plot_pc3", plot_dataframe.columns)


class BuildPcaPlotFigureTests(unittest.TestCase):
    def test_build_pca_plot_figure_returns_2d_trace_for_2d_dataframe(self) -> None:
        plot_dataframe = dataset_module.build_pca_plot_dataframe(
            cluster_assignments_dataframe=pd.DataFrame(
                {
                    "cluster_label": [0, 1],
                    "pc1": [1.0, 2.0],
                    "pc2": [3.0, 4.0],
                }
            ),
            plot_dimension_mode="2d",
        )

        figure_object = render_module.build_pca_plot_figure(
            plot_dataframe=plot_dataframe,
            figure_title="PCA 2D",
        )

        self.assertEqual(figure_object.data[0].type, "scatter")

    def test_build_pca_plot_figure_returns_3d_trace_for_3d_dataframe(self) -> None:
        plot_dataframe = dataset_module.build_pca_plot_dataframe(
            cluster_assignments_dataframe=pd.DataFrame(
                {
                    "cluster_label": [0, 1],
                    "pc1": [1.0, 2.0],
                    "pc2": [3.0, 4.0],
                    "pc3": [5.0, 6.0],
                }
            ),
            plot_dimension_mode="3d",
        )

        figure_object = render_module.build_pca_plot_figure(
            plot_dataframe=plot_dataframe,
            figure_title="PCA 3D",
        )

        self.assertEqual(figure_object.data[0].type, "scatter3d")


if __name__ == "__main__":
    unittest.main()
