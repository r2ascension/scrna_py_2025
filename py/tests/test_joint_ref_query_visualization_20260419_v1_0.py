import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd


MODULE_PATH = pathlib.Path("/home/h2048/script/py/joint_ref_query_visualization_20260419_v1_0.py")


def load_module():
    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Module under test does not exist yet: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "joint_ref_query_visualization_20260419_v1_0",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class JointRefQueryVisualizationHelpersTest(unittest.TestCase):
    def test_get_mode_config_returns_expected_paths(self):
        module = load_module()
        allcell_cfg = module.get_mode_config("allcell")
        bcell_cfg = module.get_mode_config("bcell")

        self.assertEqual(allcell_cfg["mode"], "allcell")
        self.assertIn("query_mapped_to_reference.h5ad", str(allcell_cfg["query_h5ad"]))
        self.assertIn("joint_ref_query_visualization/allcell", str(allcell_cfg["output_dir"]))

        self.assertEqual(bcell_cfg["mode"], "bcell")
        self.assertIn(
            "bcell_reference_plus_query_schpl_v1_2_scanvi_umap_20260413.h5ad",
            str(bcell_cfg["merged_h5ad"]),
        )
        self.assertIn("joint_ref_query_visualization/bcell", str(bcell_cfg["output_dir"]))

    def test_resolve_umap_key_prefers_requested_then_fallbacks(self):
        module = load_module()
        keys = ["X_umap_scanvi", "X_umap"]
        self.assertEqual(module.resolve_umap_key(keys, preferred="X_umap"), "X_umap")
        self.assertEqual(
            module.resolve_umap_key(["X_umap_scanvi", "X_scvi"], preferred="X_umap"),
            "X_umap_scanvi",
        )

        with self.assertRaises(KeyError):
            module.resolve_umap_key(["X_scvi", "X_pca"])

    def test_filter_available_markers_drops_missing_and_empty_groups(self):
        module = load_module()
        markers = {
            "A": ["G1", "G2"],
            "B": ["G3"],
            "C": ["G9"],
        }

        filtered, missing = module.filter_available_markers(markers, ["G1", "G3"])

        self.assertEqual(filtered, {"A": ["G1"], "B": ["G3"]})
        self.assertEqual(missing["A"], ["G2"])
        self.assertEqual(missing["C"], ["G9"])

    def test_compute_group_proportions_normalizes_within_sample_and_source(self):
        module = load_module()
        obs = pd.DataFrame(
            {
                "sample": ["s1", "s1", "s1", "s2", "s2", "s3"],
                "data_source": ["reference", "reference", "reference", "query", "query", "query"],
                "plot_label": ["A", "A", "B", "A", "B", "B"],
            }
        )

        proportions = module.compute_group_proportions(
            obs,
            sample_key="sample",
            label_key="plot_label",
            source_key="data_source",
        )

        s1 = proportions[(proportions["sample"] == "s1") & (proportions["data_source"] == "reference")]
        s2 = proportions[(proportions["sample"] == "s2") & (proportions["data_source"] == "query")]

        self.assertAlmostEqual(float(s1.loc[s1["plot_label"] == "A", "proportion"].iloc[0]), 2 / 3)
        self.assertAlmostEqual(float(s1.loc[s1["plot_label"] == "B", "proportion"].iloc[0]), 1 / 3)
        self.assertAlmostEqual(float(s2.loc[s2["plot_label"] == "A", "proportion"].iloc[0]), 1 / 2)
        self.assertAlmostEqual(float(s2.loc[s2["plot_label"] == "B", "proportion"].iloc[0]), 1 / 2)

    def test_build_plot_label_uses_reference_and_query_specific_columns(self):
        module = load_module()
        obs = pd.DataFrame(
            {
                "data_source": ["reference", "reference", "query", "query"],
                "ref_label": ["RefA", None, "ignore_me", "ignore_me_too"],
                "qry_label": ["ignore", "ignore", "QryA", None],
                "qry_pred": ["ignore", "ignore", "PredA", "PredB"],
            }
        )

        labels = module.build_plot_label(
            obs,
            source_key="data_source",
            ref_label_key="ref_label",
            query_label_key="qry_label",
            query_fallback_key="qry_pred",
            unknown_label="Unknown",
        )

        self.assertEqual(labels.tolist(), ["RefA", "Unknown", "QryA", "PredB"])

    def test_select_best_bcell_candidate_prefers_complete_fields(self):
        module = load_module()
        candidates = [
            {
                "path": "/tmp/older.h5ad",
                "has_data_source": True,
                "has_umap": True,
                "has_query_final": False,
                "marker_hits": 6,
            },
            {
                "path": "/tmp/newer.h5ad",
                "has_data_source": True,
                "has_umap": True,
                "has_query_final": True,
                "marker_hits": 4,
            },
        ]

        best = module.select_best_bcell_candidate(candidates)
        self.assertEqual(best["path"], "/tmp/newer.h5ad")

    def test_plotting_smoke_writes_expected_outputs(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir)
            obs = pd.DataFrame(
                {
                    "data_source": ["reference", "reference", "query", "query"],
                    "plot_label": ["A", "B", "A", "B"],
                    "source_label": [
                        "reference · A",
                        "reference · B",
                        "query · A",
                        "query · B",
                    ],
                    "sample": ["s1", "s1", "s2", "s2"],
                    "mapping_confidence": [np.nan, np.nan, 0.91, 0.67],
                },
                index=[f"cell_{i}" for i in range(4)],
            )
            var = pd.DataFrame(index=pd.Index(["G1", "G2", "G3"], name="gene"))
            adata = ad.AnnData(
                X=np.array(
                    [
                        [2.0, 0.0, 1.0],
                        [0.0, 3.0, 0.0],
                        [1.0, 1.0, 0.0],
                        [0.0, 2.0, 4.0],
                    ],
                    dtype=np.float32,
                ),
                obs=obs,
                var=var,
            )
            adata.obsm["X_umap"] = np.array(
                [[0.0, 0.0], [1.0, 0.2], [0.1, 1.0], [1.1, 1.2]], dtype=np.float32
            )

            marker_dict = {"A": ["G1", "G2"], "B": ["G3"]}

            module.plot_umap_overview(
                adata,
                output_dir=outdir,
                title="smoke overview",
                label_palette={"A": "#1f77b4", "B": "#ff7f0e", "Unknown": "#999999"},
            )
            module.plot_dotplot(
                adata,
                marker_dict=marker_dict,
                output_dir=outdir,
                title="smoke dotplot",
                group_order=["reference · A", "query · A", "reference · B", "query · B"],
            )
            module.plot_feature_grid(
                adata,
                genes=["G1", "G2", "G3"],
                output_dir=outdir,
                title="smoke featureplot",
            )
            module.plot_proportion_violin(
                adata.obs,
                output_dir=outdir,
                title="smoke violin",
                label_order=["A", "B"],
            )

            for expected in [
                "joint_umap_overview.pdf",
                "joint_dotplot.pdf",
                "joint_featureplot.pdf",
                "joint_proportion_violin.pdf",
            ]:
                self.assertTrue((outdir / expected).exists(), expected)

    def test_write_mode_summary_handles_dataframe_details(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir)
            obs = pd.DataFrame(
                {
                    "data_source": ["reference", "query"],
                    "plot_label": ["A", "A"],
                },
                index=["cell_1", "cell_2"],
            )
            adata = ad.AnnData(
                X=np.array([[1.0], [2.0]], dtype=np.float32),
                obs=obs,
                var=pd.DataFrame(index=pd.Index(["G1"], name="gene")),
            )
            availability = pd.DataFrame(
                {
                    "gene": ["G1"],
                    "in_reference": [True],
                    "in_query": [True],
                }
            )

            module.write_mode_summary(
                mode="bcell",
                output_dir=outdir,
                plot_adata=adata,
                availability=availability,
                details={
                    "selected_input": "/tmp/example.h5ad",
                    "marker_availability": availability,
                },
            )

            summary = json.loads((outdir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["mode"], "bcell")
            self.assertEqual(summary["selected_input"], "/tmp/example.h5ad")


if __name__ == "__main__":
    unittest.main()