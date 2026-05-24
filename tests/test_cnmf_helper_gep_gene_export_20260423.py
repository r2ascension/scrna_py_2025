import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd


matplotlib.use("Agg")


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
    def test_load_usage_matrix_falls_back_to_npz_when_consensus_txt_is_blank(self):
        module = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            cnmf_run_dir = pathlib.Path(tmpdir)
            (cnmf_run_dir / "cnmf_tmp").mkdir(parents=True, exist_ok=True)

            blank_usage = pd.DataFrame(index=["cell_1", "cell_2"], columns=[1, 2])
            blank_usage.to_csv(
                cnmf_run_dir / "toy_run.usages.k_2.dt_0_1.consensus.txt",
                sep="\t",
            )

            np.savez(
                cnmf_run_dir / "cnmf_tmp" / "toy_run.usages.k_2.dt_0_1.consensus.df.npz",
                data=np.array([[0.9, 0.1], [0.2, 0.8]], dtype=np.float32),
                index=np.array(["cell_1", "cell_2"], dtype=object),
                columns=np.array([1, 2], dtype=np.int64),
            )

            usage_df = module.load_usage_matrix(
                cnmf_run_dir=cnmf_run_dir,
                name="toy_run",
                k=2,
            )

        self.assertIsNotNone(usage_df)
        self.assertEqual(usage_df.columns.tolist(), ["GEP_1", "GEP_2"])
        self.assertTrue(np.allclose(usage_df.loc["cell_1"].tolist(), [0.9, 0.1]))
        self.assertTrue(np.allclose(usage_df.loc["cell_2"].tolist(), [0.2, 0.8]))

    def test_load_usage_matrix_reconstructs_from_tpm_when_txt_and_npz_are_blank(self):
        module = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            cnmf_run_dir = pathlib.Path(tmpdir)
            cnmf_tmp_dir = cnmf_run_dir / "cnmf_tmp"
            cnmf_tmp_dir.mkdir(parents=True, exist_ok=True)

            blank_usage = pd.DataFrame(index=["cell_1", "cell_2"], columns=[1, 2])
            blank_usage.to_csv(
                cnmf_run_dir / "toy_run.usages.k_2.dt_0_1.consensus.txt",
                sep="\t",
            )
            np.savez(
                cnmf_tmp_dir / "toy_run.usages.k_2.dt_0_1.consensus.df.npz",
                data=np.full((2, 2), np.nan, dtype=np.float32),
                index=np.array(["cell_1", "cell_2"], dtype=object),
                columns=np.array([1, 2], dtype=np.int64),
            )

            np.savez(
                cnmf_tmp_dir / "toy_run.gene_spectra_tpm.k_2.dt_0_1.df.npz",
                data=np.array(
                    [
                        [9000.0, 1000.0, 0.0],
                        [0.0, 1000.0, 9000.0],
                    ],
                    dtype=np.float32,
                ),
                index=np.array(["Program_1", "Program_2"], dtype=object),
                columns=np.array(["G1", "G2", "G3"], dtype=object),
            )

            obs = pd.DataFrame(index=pd.Index(["cell_1", "cell_2"], name="cell"))
            var = pd.DataFrame(index=pd.Index(["G1", "G2", "G3"], name="gene"))
            adata = ad.AnnData(
                X=np.array(
                    [
                        [9000.0, 1000.0, 0.0],
                        [0.0, 1000.0, 9000.0],
                    ],
                    dtype=np.float32,
                ),
                obs=obs,
                var=var,
            )
            adata.write_h5ad(cnmf_tmp_dir / "toy_run.tpm.h5ad")

            usage_df = module.load_usage_matrix(
                cnmf_run_dir=cnmf_run_dir,
                name="toy_run",
                k=2,
            )

            cache_path = cnmf_tmp_dir / "toy_run.usages.k_2.reconstructed.df.npz"
            self.assertIsNotNone(usage_df)
            self.assertTrue(cache_path.exists())
            self.assertEqual(usage_df.columns.tolist(), ["GEP_1", "GEP_2"])
            self.assertTrue(np.allclose(usage_df.sum(axis=1).values, [1.0, 1.0], atol=1e-5))
            self.assertGreater(usage_df.loc["cell_1", "GEP_1"], 0.9)
            self.assertGreater(usage_df.loc["cell_2", "GEP_2"], 0.9)

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

    def test_extract_top_genes_filters_clone_style_lncRNA_symbols_before_ranking(self):
        module = load_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            cnmf_run_dir = pathlib.Path(tmpdir)
            spectra = pd.DataFrame(
                {
                    "AL114490.2": [0.99, 0.98],
                    "AC124014.1": [0.97, 0.96],
                    "MS4A1": [0.80, 0.10],
                    "CD79A": [0.70, 0.20],
                    "BANK1": [0.60, 0.30],
                    "CD3D": [0.10, 0.90],
                    "IL7R": [0.20, 0.80],
                },
                index=["Program_1", "Program_2"],
            )
            spectra.to_csv(
                cnmf_run_dir / "toy_run.gene_spectra_score.k_2.dt_0_1.txt",
                sep="\t",
            )

            wide_df = module.extract_top_genes_per_gep(
                cnmf_run_dir=cnmf_run_dir,
                name="toy_run",
                k=2,
                n_top=2,
                lineage_context="B cell",
            )
            long_df = module.extract_top_genes_with_scores_per_gep(
                cnmf_run_dir=cnmf_run_dir,
                name="toy_run",
                k=2,
                n_top=2,
                lineage_context="B cell",
            )

        self.assertIsNotNone(wide_df)
        self.assertIsNotNone(long_df)
        self.assertEqual(wide_df["GEP_1"].tolist(), ["MS4A1", "CD79A"])
        self.assertEqual(wide_df["GEP_2"].tolist(), ["CD3D", "IL7R"])
        self.assertFalse(long_df["gene"].str.match(r"^(AC|AL)[0-9]{6}(\\.[0-9]+)?$").any())

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

    def test_plot_gep_usage_heatmap_uses_per_cell_matrix_and_exports_row_order(self):
        module = load_module()

        obs = pd.DataFrame(
            {
                "cell_type": ["L3_A", "L3_A", "L3_B", "L3_B"],
                "cell_type_L2": ["L2_A", "L2_A", "L2_B", "L2_B"],
                "cell_type_L3": ["L3_A", "L3_A", "L3_B", "L3_B"],
            },
            index=[f"cell_{i}" for i in range(4)],
        )
        var = pd.DataFrame(index=pd.Index(["G1", "G2"], name="gene"))
        adata = ad.AnnData(X=np.ones((4, 2), dtype=np.float32), obs=obs, var=var)
        usage_df = pd.DataFrame(
            {
                "GEP_1": [0.9, 0.7, 0.2, 0.1],
                "GEP_2": [0.1, 0.3, 0.8, 0.9],
            },
            index=obs.index,
        )

        captured = {}
        original_loader = module.load_usage_matrix
        from matplotlib.axes import Axes
        original_imshow = Axes.imshow

        def fake_load_usage_matrix(cnmf_run_dir, name, k):
            return usage_df.copy()

        def recording_imshow(self, arr, *args, **kwargs):
            captured["shape"] = tuple(arr.shape)
            return original_imshow(self, arr, *args, **kwargs)

        module.load_usage_matrix = fake_load_usage_matrix
        Axes.imshow = recording_imshow

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                viz_dir = pathlib.Path(tmpdir)
                success = module.plot_gep_usage_heatmap(
                    cnmf_output_dir=viz_dir,
                    name="toy_run",
                    k=2,
                    adata=adata,
                    celltype_col="cell_type",
                    viz_dir=viz_dir,
                    viz_config={"figure_format": "png", "dpi": 72},
                )

                self.assertTrue(success)
                self.assertEqual(captured.get("shape"), (4, 2))
                self.assertTrue((viz_dir / "gep_usage_heatmap_k2.png").exists())

                row_order_path = viz_dir / "gep_usage_heatmap_k2_cell_order.tsv"
                self.assertTrue(row_order_path.exists())
                row_order = pd.read_csv(row_order_path, sep="\t")
                self.assertEqual(row_order.shape[0], 4)
                self.assertEqual(
                    list(row_order.columns[:6]),
                    ["row_index", "cell", "cell_type", "dominant_gep", "dominant_usage", "cell_type_L2"],
                )
                self.assertEqual(
                    row_order.loc[:, ["cell_type_L2", "cell_type_L3"]].iloc[0].tolist(),
                    ["L2_A", "L3_A"],
                )
        finally:
            module.load_usage_matrix = original_loader
            Axes.imshow = original_imshow


if __name__ == "__main__":
    unittest.main()
