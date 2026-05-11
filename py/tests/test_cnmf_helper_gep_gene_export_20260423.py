import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd


MODULE_PATH = pathlib.Path("/home/h2048/script/py/cnmf_helper_20260419_v1_1.py")


def load_module():
    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Module under test does not exist yet: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "cnmf_helper_20260419_v1_1",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class CnmfHelperGepGeneExportTest(unittest.TestCase):
    def test_extract_top_genes_with_scores_per_gep_returns_ranked_long_table(self):
        module = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            cnmf_run_dir = pathlib.Path(tmpdir)
            spectra = pd.DataFrame(
                {
                    "G1": [0.10, 0.80],
                    "G2": [0.70, 0.20],
                    "G3": [0.40, 0.90],
                },
                index=["Program_1", "Program_2"],
            )
            spectra.to_csv(
                cnmf_run_dir / "toy_run.gene_spectra_score.k_2.dt_0_1.txt",
                sep="\t",
            )

            long_df = module.extract_top_genes_with_scores_per_gep(
                cnmf_run_dir=cnmf_run_dir,
                name="toy_run",
                k=2,
                n_top=2,
            )

        self.assertIsNotNone(long_df)
        self.assertEqual(
            list(long_df.columns),
            ["k", "gep", "rank", "gene", "score"],
        )
        self.assertEqual(len(long_df), 4)

        gep1 = long_df[long_df["gep"] == "GEP_1"].sort_values("rank")
        self.assertEqual(gep1["gene"].tolist(), ["G2", "G3"])
        self.assertTrue(np.allclose(gep1["score"].tolist(), [0.70, 0.40]))

        gep2 = long_df[long_df["gep"] == "GEP_2"].sort_values("rank")
        self.assertEqual(gep2["gene"].tolist(), ["G3", "G1"])
        self.assertTrue(np.allclose(gep2["score"].tolist(), [0.90, 0.80]))

    def test_run_cnmf_full_exports_gep_gene_tables_and_records_paths(self):
        module = load_module()

        obs = pd.DataFrame(
            {
                "sample": ["s1", "s1", "s2", "s2"],
                "cell_type": ["B", "B", "B", "B"],
            },
            index=[f"cell_{i}" for i in range(4)],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2", "G3"], name="gene"))
        adata = ad.AnnData(X=np.ones((4, 3), dtype=np.float32), obs=obs, var=var)
        adata.layers["counts"] = adata.X.copy()

        original_cnmf_available = module.CNMF_AVAILABLE
        original_select_hvg_robust = module.select_hvg_robust
        original_prepare_cnmf_inputs = module.prepare_cnmf_inputs
        original_run_cnmf_pipeline = module.run_cnmf_pipeline
        original_calculate_k_stability_metrics = module.calculate_k_stability_metrics
        original_generate_all_visualizations = module.generate_all_visualizations

        module.CNMF_AVAILABLE = True

        def fake_select_hvg_robust(*args, **kwargs):
            return ["G1", "G2", "G3"], "non-batch-seurat_v3", {"cell_cycle": [], "stress_ieg": []}

        def fake_prepare_cnmf_inputs(_adata, _hvg_genes, output_dir):
            hvg_counts = output_dir / "cnmf_input_hvg_counts.h5ad"
            tp10k = output_dir / "cnmf_input_tp10k.h5ad"
            hvg_txt = output_dir / "cnmf_input_hvg_genes.txt"
            hvg_counts.write_text("dummy")
            tp10k.write_text("dummy")
            hvg_txt.write_text("G1\nG2\nG3\n")
            return hvg_counts, tp10k, hvg_txt

        def fake_run_cnmf_pipeline(*, hvg_counts_path, tp10k_path, hvg_txt_path, k_range, output_dir, name, cnmf_config):
            run_dir = output_dir / "cnmf_output" / name
            run_dir.mkdir(parents=True, exist_ok=True)
            spectra = pd.DataFrame(
                {
                    "G1": [0.10, 0.80],
                    "G2": [0.70, 0.20],
                    "G3": [0.40, 0.90],
                },
                index=["Program_1", "Program_2"],
            )
            spectra.to_csv(
                run_dir / f"{name}.gene_spectra_score.k_2.dt_0_1.txt",
                sep="\t",
            )
            return True

        def fake_calculate_k_stability_metrics(*, cnmf_output_dir, name, k_range, output_dir, viz_config):
            recommendation = {
                "recommended_k": 2,
                "confidence": "medium",
                "reasoning": ["test recommendation"],
                "alternative_k": [],
            }
            (output_dir / "k_selection_recommendation.json").write_text(json.dumps(recommendation))
            return {
                2: {
                    "mean_usage_entropy": 0.5,
                    "usage_sparsity_gini": 0.5,
                    "max_usage_concentration": 0.8,
                    "mean_pairwise_distance": 1.0,
                }
            }

        module.select_hvg_robust = fake_select_hvg_robust
        module.prepare_cnmf_inputs = fake_prepare_cnmf_inputs
        module.run_cnmf_pipeline = fake_run_cnmf_pipeline
        module.calculate_k_stability_metrics = fake_calculate_k_stability_metrics
        module.generate_all_visualizations = lambda **kwargs: {2: {"local_density": True, "clustergram": True, "usage_heatmap": True}}

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                result = module.run_cnmf_full(
                    adata=adata,
                    output_dir=pathlib.Path(tmpdir),
                    run_name="toy run",
                    k_range=[2],
                    celltype_col="cell_type",
                    batch_col="sample",
                    use_batch_hvg=True,
                )

                self.assertTrue(result["success"])
                self.assertIn("gep_gene_tables", result["paths"])
                self.assertEqual(
                    sorted(result["paths"]["gep_gene_tables"].keys()),
                    ["directory", "long_by_k", "wide_by_k"],
                )

                wide_path = pathlib.Path(result["paths"]["gep_gene_tables"]["wide_by_k"]["2"])
                long_path = pathlib.Path(result["paths"]["gep_gene_tables"]["long_by_k"]["2"])
                self.assertTrue(wide_path.exists())
                self.assertTrue(long_path.exists())

                long_df = pd.read_csv(long_path, sep="\t")
                self.assertEqual(
                    list(long_df.columns),
                    ["k", "gep", "rank", "gene", "score"],
                )
                self.assertIn("run_summary.json", {p.name for p in pathlib.Path(tmpdir).iterdir()})
        finally:
            module.CNMF_AVAILABLE = original_cnmf_available
            module.select_hvg_robust = original_select_hvg_robust
            module.prepare_cnmf_inputs = original_prepare_cnmf_inputs
            module.run_cnmf_pipeline = original_run_cnmf_pipeline
            module.calculate_k_stability_metrics = original_calculate_k_stability_metrics
            module.generate_all_visualizations = original_generate_all_visualizations


if __name__ == "__main__":
    unittest.main()
