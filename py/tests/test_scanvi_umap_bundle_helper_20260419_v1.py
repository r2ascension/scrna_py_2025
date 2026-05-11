import importlib.util
import json
import pathlib
import sys
import tempfile
import types
import unittest

import anndata as ad
import numpy as np
import pandas as pd

MODULE_PATH = pathlib.Path("/home/h2048/script/py/scanvi_umap_bundle_helper_20260419_v1.py")


def load_module():
    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Module under test does not exist yet: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "scanvi_umap_bundle_helper_20260419_v1",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FakeUMAP:
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.embedding_ = None

    def fit_transform(self, coords):
        coords = np.asarray(coords, dtype=np.float32)
        self.embedding_ = np.column_stack([coords[:, 0], coords[:, 0] + 1.0]).astype(np.float32)
        return self.embedding_

    def transform(self, coords):
        coords = np.asarray(coords, dtype=np.float32)
        return np.column_stack([coords[:, 0], coords[:, 0] + 1.0]).astype(np.float32)


class ScanviUmapBundleHelperTest(unittest.TestCase):
    def _make_adata(self):
        obs = pd.DataFrame(
            {
                "data_source": ["reference", "reference", "query"],
                "sample": ["s1", "s1", "s2"],
            },
            index=["cell1", "cell2", "cell3"],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2"], name="gene"))
        adata = ad.AnnData(
            X=np.array([[1.0, 0.0], [0.0, 1.0], [3.0, 2.0]], dtype=np.float32),
            obs=obs,
            var=var,
        )
        adata.obsm["X_scanvi"] = np.array(
            [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]], dtype=np.float32
        )
        return adata

    def test_resolve_latent_key_prefers_requested_then_fallbacks(self):
        module = load_module()
        self.assertEqual(module.resolve_latent_key(["X_scanvi", "X_scvi"], preferred="X_scanvi"), "X_scanvi")
        self.assertEqual(module.resolve_latent_key(["X_scANVI", "X_scvi"], preferred="missing"), "X_scANVI")
        with self.assertRaises(KeyError):
            module.resolve_latent_key(["X_pca", "X_harmony"])

    def test_store_umap_coordinates_sets_aliases_and_default(self):
        module = load_module()
        adata = self._make_adata()
        coords = np.array([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]], dtype=np.float64)

        info = module.store_umap_coordinates(
            adata,
            coords,
            primary_umap_key="X_umap_scANVI",
            alias_keys=("X_umap_scanvi", "X_umap_scanvi_corrected"),
            set_default_x_umap=True,
        )

        self.assertEqual(info["primary_umap_key"], "X_umap_scANVI")
        for key in ["X_umap_scANVI", "X_umap_scanvi", "X_umap_scanvi_corrected", "X_umap"]:
            self.assertIn(key, adata.obsm)
            self.assertEqual(adata.obsm[key].dtype, np.float32)
        self.assertTrue(np.allclose(adata.obsm["X_umap"], coords.astype(np.float32)))

    def test_fit_bundle_writes_operator_and_manifest(self):
        module = load_module()
        adata = self._make_adata()
        dumped = {}

        module._load_umap_class = lambda: FakeUMAP
        module._load_joblib = lambda: types.SimpleNamespace(
            dump=lambda obj, path: (dumped.setdefault("operator", obj), pathlib.Path(path).write_text("ok", encoding="utf-8")),
            load=lambda path: dumped["operator"],
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            result = module.fit_bundle(
                adata,
                "X_scanvi",
                tmpdir,
                primary_umap_key="X_umap_scanvi",
                alias_keys=("X_umap_scanvi_corrected",),
                source_h5ad="/tmp/source.h5ad",
                extra_manifest={"lineage": "bcell"},
            )
            self.assertTrue(pathlib.Path(result["operator_path"]).exists())
            self.assertTrue(pathlib.Path(result["manifest_path"]).exists())
            self.assertIn("scanvi_umap_bundle", adata.uns)
            self.assertEqual(result["manifest"]["source_h5ad"], "/tmp/source.h5ad")
            self.assertEqual(result["manifest"]["lineage"], "bcell")
            self.assertEqual(result["manifest"]["data_source_counts"], {"reference": 2, "query": 1})
            self.assertIn("X_umap", adata.obsm)

            payload = json.loads(pathlib.Path(result["manifest_path"]).read_text(encoding="utf-8"))
            self.assertEqual(payload["primary_umap_key"], "X_umap_scanvi")

    def test_refresh_from_existing_latent_can_write_sidecar(self):
        module = load_module()
        adata = self._make_adata()
        dumped = {}

        module._load_umap_class = lambda: FakeUMAP
        module._load_joblib = lambda: types.SimpleNamespace(
            dump=lambda obj, path: (dumped.setdefault("operator", obj), pathlib.Path(path).write_text("ok", encoding="utf-8")),
            load=lambda path: dumped["operator"],
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            sidecar = pathlib.Path(tmpdir) / "refreshed.h5ad"
            result = module.refresh_from_existing_latent(
                adata,
                output_dir=tmpdir,
                primary_umap_key="X_umap_refined",
                alias_keys=("X_umap_scanvi",),
                source_h5ad="/tmp/original.h5ad",
                write_sidecar_h5ad=True,
                sidecar_h5ad_path=sidecar,
            )
            self.assertEqual(result["manifest"]["mode"], "refresh_from_existing_latent")
            self.assertTrue(sidecar.exists())
            reloaded = ad.read_h5ad(sidecar)
            self.assertIn("X_umap", reloaded.obsm)
            self.assertIn("X_umap_refined", reloaded.obsm)
            self.assertIn("scanvi_umap_refresh", reloaded.uns)

    def test_project_query_supports_loaded_or_direct_operator(self):
        module = load_module()
        operator = FakeUMAP()
        direct = module.project_query(np.array([[3.0, 4.0]], dtype=np.float32), operator=operator)
        self.assertTrue(np.allclose(direct, np.array([[3.0, 4.0]], dtype=np.float32)[:, [0, 0]] + np.array([[0.0, 1.0]], dtype=np.float32)))

        module._load_joblib = lambda: types.SimpleNamespace(load=lambda path: operator)
        loaded = module.project_query(np.array([[5.0, 6.0]], dtype=np.float32), operator_path="/tmp/fake.joblib")
        self.assertEqual(loaded.dtype, np.float32)
        self.assertTrue(np.allclose(loaded, np.array([[5.0, 6.0]], dtype=np.float32)[:, [0, 0]] + np.array([[0.0, 1.0]], dtype=np.float32)))


if __name__ == "__main__":
    unittest.main()
