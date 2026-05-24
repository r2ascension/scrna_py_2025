#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import importlib.util
import sys
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SCRIPT_PATH = Path("/home/h2048/script/py/epithelial_milopy_tissue_celltype_20260419_v1.py")


def load_module(testcase: unittest.TestCase):
    if not SCRIPT_PATH.exists():
        testcase.fail(f"Expected Milo module is missing: {SCRIPT_PATH}")
    spec = importlib.util.spec_from_file_location("epithelial_milo_20260419_v1", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        testcase.fail(f"Unable to load module spec from: {SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TestEpithelialMilopyHelpers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_module(unittest.TestCase("runTest"))

    def make_adata(self):
        obs_rows = []

        def add_cells(label_l2, label_l3, sample, tissue, n):
            for _ in range(n):
                obs_rows.append(
                    {
                        "sample": sample,
                        "tissue": tissue,
                        "cell_type_L2": label_l2,
                        "cell_type_L3": label_l3,
                    }
                )

        # Ready-for-analysis L2/L3 groups: enough cells, >=2 tissues, >=2 samples/tissue.
        add_cells("Secretory_Lineage", "Goblet", "nose_s1", "nose", 6)
        add_cells("Secretory_Lineage", "Goblet", "nose_s2", "nose", 6)
        add_cells("Secretory_Lineage", "Goblet", "lung_s1", "lung parenchyma", 6)
        add_cells("Secretory_Lineage", "Goblet", "lung_s2", "lung parenchyma", 6)

        # Only one tissue after filtering -> should be skipped.
        add_cells("Alveolar", "AT2", "lung_s3", "lung parenchyma", 7)
        add_cells("Alveolar", "AT2", "lung_s4", "lung parenchyma", 7)

        # Blank labels should be ignored by summary.
        add_cells("", "", "noise_s1", "nose", 4)

        obs = pd.DataFrame(obs_rows)
        obs.index = [f"cell_{i}" for i in range(obs.shape[0])]
        x = np.zeros((obs.shape[0], 1), dtype=np.float32)
        adata = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=["placeholder_feature"]))
        adata.obsm["X_scanvi"] = np.zeros((adata.n_obs, 5), dtype=np.float32)
        adata.obsm["X_umap_scanvi"] = np.zeros((adata.n_obs, 2), dtype=np.float32)
        return adata

    def test_choose_available_key_prefers_first_present_candidate(self):
        key = self.mod.choose_available_key(
            available_keys=["X_scanvi", "X_scvi", "X_umap_scanvi"],
            candidates=["X_scanvi_major", "X_scanvi", "X_scvi"],
            kind="latent embedding",
        )
        self.assertEqual(key, "X_scanvi")

    def test_summarize_level_groups_excludes_blank_labels(self):
        adata = self.make_adata()
        cfg = self.mod.MiloRunConfig(
            min_cells_per_group=10,
            min_cells_per_sample=5,
            min_samples_per_tissue=2,
        )
        summary = self.mod.summarize_level_groups(adata, "cell_type_L3", cfg)
        labels = summary["cell_type"].tolist()
        self.assertIn("Goblet", labels)
        self.assertIn("AT2", labels)
        self.assertNotIn("", labels)

    def test_prepare_level_subset_marks_ready_and_skipped_groups(self):
        adata = self.make_adata()
        cfg = self.mod.MiloRunConfig(
            min_cells_per_group=10,
            min_cells_per_sample=5,
            min_samples_per_tissue=2,
        )

        goblet_subset, goblet_info = self.mod.prepare_level_subset(
            adata,
            label_col="cell_type_L3",
            label_value="Goblet",
            cfg=cfg,
        )
        self.assertIsNotNone(goblet_subset)
        self.assertEqual(goblet_info["status"], "ready")
        self.assertEqual(goblet_subset.obs["tissue"].nunique(), 2)
        self.assertEqual(goblet_subset.obs["sample"].nunique(), 4)

        at2_subset, at2_info = self.mod.prepare_level_subset(
            adata,
            label_col="cell_type_L3",
            label_value="AT2",
            cfg=cfg,
        )
        self.assertIsNone(at2_subset)
        self.assertEqual(at2_info["status"], "skipped")
        self.assertIn("fewer than 2 tissues", at2_info["reason"])


if __name__ == "__main__":
    unittest.main()
