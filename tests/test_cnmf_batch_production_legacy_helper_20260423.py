import importlib.util
import pathlib
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd


MODULE_PATH = pathlib.Path("/home/h2048/script/py/cnmf_batch_production_20251227_v1_2.py")


def load_module():
    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Module under test does not exist yet: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "cnmf_batch_production_20251227_v1_2",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class CnmfLegacyBatchHelperIntegrationTest(unittest.TestCase):
    def test_get_lineage_config_returns_bcell_latest_h5ad(self):
        module = load_module()
        cfg = module.get_lineage_config("bcell")

        self.assertEqual(cfg["mode"], "bcell")
        self.assertIn(
            "bcell_tissue_comparison_final.h5ad",
            str(cfg["input_h5ad"]),
        )
        self.assertIn("bcell", str(cfg["output_root"]))

    def test_get_lineage_config_raises_for_unknown_mode(self):
        module = load_module()
        with self.assertRaises((KeyError, ValueError)):
            module.get_lineage_config("unknown_lineage")

    def test_build_run_settings_smoke_reduces_workload_and_keeps_batch_aware_only(self):
        module = load_module()
        settings = module.build_run_settings(mode="bcell", run_type="smoke")

        self.assertEqual(settings["mode"], "bcell")
        self.assertEqual(settings["run_type"], "smoke")
        self.assertEqual(settings["tracks"], ["batch_aware"])
        self.assertIsInstance(settings["subset_n_cells"], int)
        self.assertGreater(settings["subset_n_cells"], 0)
        self.assertEqual(settings["cnmf_config"]["n_workers"], 1)
        self.assertLess(settings["cnmf_config"]["n_iter"], module.CNMF_HELPER.DEFAULT_CNMF_CONFIG["n_iter"])
        self.assertLess(len(settings["k_range"]), len(module.CNMF_HELPER.K_RANGE_MEDIUM))
        self.assertIn("/bcell/smoke", str(settings["output_dir"]))

    def test_build_run_settings_full_defaults_to_batch_aware_only(self):
        module = load_module()
        settings = module.build_run_settings(mode="bcell", run_type="full")

        self.assertEqual(settings["mode"], "bcell")
        self.assertEqual(settings["run_type"], "full")
        self.assertEqual(settings["tracks"], ["batch_aware"])
        self.assertIsNone(settings["subset_n_cells"])
        self.assertEqual(settings["cnmf_config"]["n_iter"], module.CNMF_HELPER.DEFAULT_CNMF_CONFIG["n_iter"])
        self.assertIn("/bcell/full", str(settings["output_dir"]))

    def test_build_run_settings_keeps_explicit_dual_track_override(self):
        module = load_module()
        settings = module.build_run_settings(
            mode="bcell",
            run_type="full",
            tracks=["uncorrected", "batch_aware"],
        )

        self.assertEqual(settings["tracks"], ["uncorrected", "batch_aware"])

    def test_execute_helper_track_delegates_to_run_cnmf_full(self):
        module = load_module()
        settings = module.build_run_settings(mode="bcell", run_type="smoke")

        obs = pd.DataFrame(
            {
                "sample": ["s1", "s1", "s2", "s2"],
                "tissue": ["nose", "nose", "lung", "lung"],
                "cell_type_L3": ["Naive_B", "Memory_B", "Naive_B", "Memory_B"],
            },
            index=[f"cell_{i}" for i in range(4)],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2", "G3"], name="gene"))
        adata = ad.AnnData(X=np.ones((4, 3), dtype=np.float32), obs=obs, var=var)
        adata.layers["counts"] = adata.X.copy()

        called = {}
        original = module.CNMF_HELPER.run_cnmf_full

        def fake_run_cnmf_full(*, adata, output_dir, run_name, k_range, celltype_col, batch_col, use_batch_hvg, cnmf_config, viz_config):
            called["n_obs"] = adata.n_obs
            called["output_dir"] = pathlib.Path(output_dir)
            called["run_name"] = run_name
            called["k_range"] = list(k_range)
            called["celltype_col"] = celltype_col
            called["batch_col"] = batch_col
            called["use_batch_hvg"] = use_batch_hvg
            called["cnmf_config"] = dict(cnmf_config)
            called["viz_config"] = dict(viz_config)
            return {
                "success": True,
                "run_name": run_name,
                "k_range": list(k_range),
                "recommendation": {"recommended_k": k_range[0] if k_range else None},
            }

        module.CNMF_HELPER.run_cnmf_full = fake_run_cnmf_full
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                result = module.execute_helper_track(
                    adata=adata,
                    dataset_id="bcell_latest",
                    track_name="batch_aware",
                    run_settings=settings,
                    output_dir=pathlib.Path(tmpdir),
                    celltype_col="cell_type_L3",
                    batch_col="sample",
                )
        finally:
            module.CNMF_HELPER.run_cnmf_full = original

        self.assertTrue(result["success"])
        self.assertEqual(called["n_obs"], 4)
        self.assertTrue(called["use_batch_hvg"])
        self.assertEqual(called["batch_col"], "sample")
        self.assertEqual(called["celltype_col"], "cell_type_L3")
        self.assertEqual(called["cnmf_config"]["n_workers"], 1)
        self.assertEqual(called["run_name"], "bcell_latest_batch_aware_smoke")

    def test_run_lineage_pipeline_smoke_subsets_and_runs_batch_aware_once(self):
        module = load_module()

        obs = pd.DataFrame(
            {
                "sample": [f"s{i%4}" for i in range(12)],
                "tissue": ["nose" if i % 2 == 0 else "lung" for i in range(12)],
                "cell_type_L3": ["Naive_B" if i % 3 else "Memory_B" for i in range(12)],
            },
            index=[f"cell_{i}" for i in range(12)],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2", "G3"], name="gene"))
        adata = ad.AnnData(X=np.ones((12, 3), dtype=np.float32), obs=obs, var=var)
        adata.layers["counts"] = adata.X.copy()

        calls = []
        original_load = module.load_and_validate_dataset
        original_ensure = module.ensure_counts_layer
        original_execute = module.execute_helper_track

        def fake_load(_path):
            return adata.copy(), {
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "celltype_col": "cell_type_L3",
                "batch_col": "sample",
                "has_harmony": False,
            }

        def fake_execute(*, adata, dataset_id, track_name, run_settings, output_dir, celltype_col, batch_col):
            calls.append(
                {
                    "n_obs": adata.n_obs,
                    "dataset_id": dataset_id,
                    "track_name": track_name,
                    "output_dir": pathlib.Path(output_dir),
                    "celltype_col": celltype_col,
                    "batch_col": batch_col,
                    "run_type": run_settings["run_type"],
                }
            )
            return {"success": True, "run_name": f"{dataset_id}_{track_name}"}

        module.load_and_validate_dataset = fake_load
        module.ensure_counts_layer = lambda _adata: True
        module.execute_helper_track = fake_execute
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                result = module.run_lineage_pipeline(
                    mode="bcell",
                    run_type="smoke",
                    output_root=pathlib.Path(tmpdir),
                )

                self.assertTrue((pathlib.Path(tmpdir) / "smoke" / "processing_result.json").exists())
        finally:
            module.load_and_validate_dataset = original_load
            module.ensure_counts_layer = original_ensure
            module.execute_helper_track = original_execute

        self.assertTrue(result["success"])
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["track_name"], "batch_aware")
        self.assertEqual(calls[0]["run_type"], "smoke")
        self.assertEqual(calls[0]["celltype_col"], "cell_type_L3")
        self.assertEqual(calls[0]["batch_col"], "sample")
        self.assertLessEqual(calls[0]["n_obs"], adata.n_obs)

    def test_run_lineage_pipeline_prefers_specific_celltype_column_over_generic_metadata(self):
        module = load_module()

        obs = pd.DataFrame(
            {
                "sample": ["s1", "s1", "s2", "s2"],
                "cell_type": ["B", "B", "B", "B"],
                "cell_type_L3": ["Naive_B", "Memory_B", "Naive_B", "Memory_B"],
            },
            index=[f"cell_{i}" for i in range(4)],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2"], name="gene"))
        adata = ad.AnnData(X=np.ones((4, 2), dtype=np.float32), obs=obs, var=var)
        adata.layers["counts"] = adata.X.copy()

        seen = {}
        original_load = module.load_and_validate_dataset
        original_ensure = module.ensure_counts_layer
        original_execute = module.execute_helper_track

        def fake_load(_path):
            return adata.copy(), {
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "celltype_col": "cell_type",
                "batch_col": "sample",
            }

        def fake_execute(*, adata, dataset_id, track_name, run_settings, output_dir, celltype_col, batch_col):
            seen["celltype_col"] = celltype_col
            seen["batch_col"] = batch_col
            return {"success": True}

        module.load_and_validate_dataset = fake_load
        module.ensure_counts_layer = lambda _adata: True
        module.execute_helper_track = fake_execute
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                result = module.run_lineage_pipeline(
                    mode="bcell",
                    run_type="smoke",
                    output_root=pathlib.Path(tmpdir),
                )
        finally:
            module.load_and_validate_dataset = original_load
            module.ensure_counts_layer = original_ensure
            module.execute_helper_track = original_execute

        self.assertTrue(result["success"])
        self.assertEqual(seen["batch_col"], "sample")
        self.assertEqual(seen["celltype_col"], "cell_type_L3")

    def test_main_defaults_to_helper_orchestrator_not_legacy_discovery(self):
        module = load_module()

        called = {}
        original_helper = module.run_helper_orchestrator
        original_legacy = module.run_legacy_discovery_pipeline

        def fake_helper(*, mode, run_type, output_root, tracks):
            called["mode"] = mode
            called["run_type"] = run_type
            called["output_root"] = pathlib.Path(output_root)
            called["tracks"] = list(tracks)
            return {"success": True, "execution_mode": "helper_orchestrator"}

        def fail_legacy(*args, **kwargs):
            raise AssertionError("legacy discovery should not run on the default entrypoint")

        module.run_helper_orchestrator = fake_helper
        module.run_legacy_discovery_pipeline = fail_legacy
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                summary = module.main(["--output-root", tmpdir])
        finally:
            module.run_helper_orchestrator = original_helper
            module.run_legacy_discovery_pipeline = original_legacy

        self.assertTrue(summary["success"])
        self.assertEqual(called["mode"], "bcell")
        self.assertEqual(called["run_type"], "full")
        self.assertEqual(called["tracks"], ["batch_aware"])

    def test_main_legacy_discovery_requires_explicit_flag(self):
        module = load_module()

        called = {}
        original_helper = module.run_helper_orchestrator
        original_legacy = module.run_legacy_discovery_pipeline

        def fail_helper(*args, **kwargs):
            raise AssertionError("helper orchestrator should be bypassed when --legacy-discovery is set")

        def fake_legacy(*, output_dir=None):
            called["output_dir"] = pathlib.Path(output_dir)
            return {"success": True, "execution_mode": "legacy_discovery"}

        module.run_helper_orchestrator = fail_helper
        module.run_legacy_discovery_pipeline = fake_legacy
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                summary = module.main(["--legacy-discovery", "--output-root", tmpdir])
        finally:
            module.run_helper_orchestrator = original_helper
            module.run_legacy_discovery_pipeline = original_legacy

        self.assertTrue(summary["success"])
        self.assertEqual(summary["execution_mode"], "legacy_discovery")
        self.assertTrue(str(called["output_dir"]).endswith(tmpdir))


if __name__ == "__main__":
    unittest.main()
