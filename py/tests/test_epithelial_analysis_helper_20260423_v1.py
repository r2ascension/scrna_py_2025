import importlib.util
import pathlib
import sys
import unittest

import anndata as ad
import numpy as np
import pandas as pd

GENERIC_MODULE_PATH = pathlib.Path("/home/h2048/script/py/anndata_expression_viz_helper_20260423_v1.py")
MODULE_PATH = pathlib.Path("/home/h2048/script/py/epithelial_analysis_helper_20260423_v1.py")


def load_module():
    if not GENERIC_MODULE_PATH.exists():
        raise FileNotFoundError(f"Generic helper missing: {GENERIC_MODULE_PATH}")
    generic_spec = importlib.util.spec_from_file_location(
        "anndata_expression_viz_helper_20260423_v1",
        GENERIC_MODULE_PATH,
    )
    generic_module = importlib.util.module_from_spec(generic_spec)
    assert generic_spec.loader is not None
    sys.modules[generic_spec.name] = generic_module
    generic_spec.loader.exec_module(generic_module)

    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Module under test does not exist yet: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "epithelial_analysis_helper_20260423_v1",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class EpithelialAnalysisHelperTest(unittest.TestCase):
    def _make_subcluster_adata(self):
        obs = pd.DataFrame(
            {"subcluster": ["AT1 _0 ", "AT2_0", "MysteryCluster"]},
            index=["c1", "c2", "c3"],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2"], name="gene"))
        return ad.AnnData(
            X=np.array([[1.0, 0.0], [0.0, 1.0], [2.0, 2.0]], dtype=np.float32),
            obs=obs,
            var=var,
        )

    def _make_review_adata(self):
        obs = pd.DataFrame(
            {"subcluster": ["C1", "C1", "C2", "C2", "C3"]},
            index=["r1", "r2", "r3", "r4", "r5"],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2", "G3", "G4"], name="gene"))
        X = np.array(
            [
                [3.0, 0.0, 0.0, 0.0],
                [1.0, 2.0, 0.0, 0.0],
                [0.0, 0.0, 4.0, 0.0],
                [0.0, 0.0, 2.0, 5.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            dtype=np.float32,
        )
        return ad.AnnData(X=X, obs=obs, var=var)

    def test_default_schema_objects_are_nonempty(self):
        module = load_module()
        self.assertIn("AT1 _0", module.build_default_l3_mapping())
        self.assertIn("AT1_Canonical", module.build_default_marker_panels())
        self.assertIn("Alveolar", module.build_default_lineage_groups())
        self.assertIn("AT1_Canonical", module.build_default_l3_colors())
        self.assertGreater(len(module.build_misannotation_review_specs()), 0)

    def test_apply_l3_annotations_strips_whitespace_and_maps_clusters(self):
        module = load_module()
        adata = self._make_subcluster_adata()
        mapped = module.apply_l3_annotations(adata, module.build_default_l3_mapping())
        self.assertEqual(mapped.obs["cell_type_L3"].astype(str).tolist()[0], "AT1_Canonical")
        self.assertEqual(mapped.obs["cell_type_L3"].astype(str).tolist()[1], "AT2_Canonical")

    def test_apply_l3_annotations_fills_unmapped_with_original_cluster(self):
        module = load_module()
        adata = self._make_subcluster_adata()
        mapped = module.apply_l3_annotations(adata, module.build_default_l3_mapping())
        self.assertEqual(mapped.obs["cell_type_L3"].astype(str).tolist()[2], "MysteryCluster")
        self.assertEqual(mapped.uns["cell_type_L3_summary"]["n_unmapped"], 1)

    def test_build_misannotation_review_table_keeps_cluster_order(self):
        module = load_module()
        adata = self._make_review_adata()
        specs = [
            {
                "cluster": "C1",
                "wrong_label": "WrongA",
                "correct_label": "CorrectA",
                "wrong_genes": ["G1"],
                "correct_genes": ["G2"],
                "note": "first",
            },
            {
                "cluster": "C2",
                "wrong_label": "WrongB",
                "correct_label": "CorrectB",
                "wrong_genes": ["G3"],
                "correct_genes": ["G4"],
                "note": "second",
            },
        ]
        df = module.build_misannotation_review_table(adata, specs=specs, cluster_col="subcluster")
        self.assertEqual(df["review_cluster"].drop_duplicates().tolist(), ["C1", "C2"])

    def test_build_misannotation_review_table_contains_expected_columns(self):
        module = load_module()
        adata = self._make_review_adata()
        specs = [
            {
                "cluster": "C1",
                "wrong_label": "WrongA",
                "correct_label": "CorrectA",
                "wrong_genes": ["G1"],
                "correct_genes": ["G2"],
                "note": "first",
            }
        ]
        df = module.build_misannotation_review_table(adata, specs=specs, cluster_col="subcluster")
        self.assertTrue({"review_cluster", "gene", "pct_expr", "mean_expr", "plot_label", "x", "y"}.issubset(df.columns))


if __name__ == "__main__":
    unittest.main()
