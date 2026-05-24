import importlib.util
import pathlib
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

MODULE_PATH = pathlib.Path("/home/h2048/script/py/anndata_expression_viz_helper_20260423_v1.py")


def load_module():
    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Module under test does not exist yet: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "anndata_expression_viz_helper_20260423_v1",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class AnnDataExpressionVizHelperTest(unittest.TestCase):
    def _make_dense_adata(self, with_raw: bool = False):
        obs = pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["c1", "c2", "c3"])
        var = pd.DataFrame(index=pd.Index(["Gene1", "Gene2", "Gene3"], name="gene"))
        adata = ad.AnnData(
            X=np.array([[1.0, 0.0, 4.0], [3.0, 2.0, 0.0], [0.0, 5.0, 6.0]], dtype=np.float32),
            obs=obs,
            var=var,
        )
        if with_raw:
            raw = adata.copy()
            raw.var_names = pd.Index(["GeneA", "GeneB", "GeneC"], name="gene")
            adata.raw = raw
        return adata

    def _make_sparse_adata(self):
        obs = pd.DataFrame({"cell_type": ["A", "B", "B"]}, index=["s1", "s2", "s3"])
        var = pd.DataFrame(index=pd.Index(["Gene1", "Gene2"], name="gene"))
        return ad.AnnData(
            X=sparse.csr_matrix(np.array([[1.0, 0.0], [0.0, 4.0], [2.0, 0.0]], dtype=np.float32)),
            obs=obs,
            var=var,
        )

    def test_helper_exports_exist(self):
        module = load_module()
        self.assertTrue(callable(module.resolve_expression_source))
        self.assertTrue(callable(module.gene_available))
        self.assertTrue(callable(module.extract_gene_vector))
        self.assertTrue(callable(module.get_available_markers))

    def test_resolve_expression_source_prefers_raw_when_present(self):
        module = load_module()
        adata = self._make_dense_adata(with_raw=True)
        info = module.resolve_expression_source(adata)
        self.assertTrue(info["use_raw"])
        self.assertIn("GeneA", info["gene_universe"])
        self.assertEqual(info["source_name"], ".raw")

    def test_resolve_expression_source_falls_back_to_x_when_raw_missing(self):
        module = load_module()
        adata = self._make_dense_adata(with_raw=False)
        info = module.resolve_expression_source(adata)
        self.assertFalse(info["use_raw"])
        self.assertEqual(info["source_name"], ".X")
        self.assertIn("Gene1", info["gene_universe"])

    def test_get_available_markers_filters_missing_genes(self):
        module = load_module()
        available, missing = module.get_available_markers(
            {"A": ["Gene1", "Missing"]},
            ["Gene1", "Gene2"],
        )
        self.assertEqual(available, {"A": ["Gene1"]})
        self.assertIn("Missing", missing)

    def test_extract_gene_vector_handles_dense_and_sparse(self):
        module = load_module()
        sparse_adata = self._make_sparse_adata()
        dense_vec = module.extract_gene_vector(self._make_dense_adata(), "Gene1")
        sparse_vec = module.extract_gene_vector(sparse_adata, "Gene1")
        self.assertEqual(dense_vec.shape, (3,))
        self.assertEqual(sparse_vec.shape, (3,))
        self.assertTrue(np.allclose(sparse_vec, np.array([1.0, 0.0, 2.0], dtype=np.float32)))

    def test_build_mean_expression_matrix_returns_grouped_dataframe(self):
        module = load_module()
        adata = self._make_dense_adata()
        df = module.build_mean_expression_matrix(adata, ["Gene1", "Gene2"], groupby="cell_type")
        self.assertEqual(list(df.index), ["A", "B"])
        self.assertEqual(list(df.columns), ["Gene1", "Gene2"])
        self.assertAlmostEqual(df.loc["A", "Gene1"], 2.0)
        self.assertAlmostEqual(df.loc["B", "Gene2"], 5.0)

    def test_prepare_output_dirs_creates_expected_subdirs(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            dirs = module.prepare_output_dirs(tmpdir)
            self.assertTrue(dirs["figures"].is_dir())
            self.assertTrue(dirs["tables"].is_dir())
            self.assertTrue(dirs["artifacts"].is_dir())


if __name__ == "__main__":
    unittest.main()
