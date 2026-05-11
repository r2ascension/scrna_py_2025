#!/usr/bin/env Rscript

worker_path <- "/home/h2048/script/R/program_lineage_method_worker_20260507.R"
source(worker_path)

python_cmd <- "/home/h2048/miniconda3/envs/bbknn_env/bin/python"
if (!file.exists(python_cmd)) stop(sprintf("Missing Python executable: %s", python_cmd), call. = FALSE)

tmp_root <- tempfile("pa_harmony2_method_")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
h5ad_path <- file.path(tmp_root, "tiny_harmony2_input.h5ad")
create_py <- file.path(tmp_root, "create_tiny_harmony2_input.py")
writeLines(c(
  "import anndata as ad",
  "import numpy as np",
  "import pandas as pd",
  "from scipy import sparse",
  "rng = np.random.default_rng(20260511)",
  "n_cells = 24",
  "obs = pd.DataFrame({",
  "    'dataset': ['A'] * 12 + ['B'] * 12,",
  "    'cell_type_L3': ['Type1'] * 8 + ['Type2'] * 8 + ['Type3'] * 8,",
  "    'sample': ['S1'] * 6 + ['S2'] * 6 + ['S3'] * 6 + ['S4'] * 6,",
  "}, index=[f'cell_{i}' for i in range(n_cells)])",
  "X = sparse.csr_matrix(rng.poisson(1.0, size=(n_cells, 10)).astype(np.float32))",
  "adata = ad.AnnData(X=X, obs=obs)",
  "adata.obsm['X_pca'] = rng.normal(size=(n_cells, 6)).astype(np.float32)",
  paste0("adata.write_h5ad(r'", normalizePath(h5ad_path, winslash = "/", mustWork = FALSE), "')")
), create_py)
cmd_status <- system2(python_cmd, args = create_py, stdout = TRUE, stderr = TRUE)
if (!file.exists(h5ad_path)) {
  stop(sprintf("Failed to create test h5ad. Python output:\n%s", paste(cmd_status, collapse = "\n")), call. = FALSE)
}

cfg <- list(
  lineage = "tiny_lineage",
  method = "harmony2",
  lineage_dir = file.path(tmp_root, "tiny_lineage"),
  output_dir = file.path(tmp_root, "tiny_lineage", "harmony2_full"),
  status_tsv = file.path(tmp_root, "status.tsv"),
  run_stamp = "20260511_test",
  h5ad_path = h5ad_path,
  batch_col = "dataset",
  celltype_col = "cell_type_L3",
  sample_col = "sample",
  harmony2_python = python_cmd,
  harmony2_helper = "/home/h2048/script/py/harmony2_helper_20260511_v1.py",
  harmony2_n_pcs = 6L,
  harmony2_resolutions = c(0.4, 0.8),
  harmony2_default_resolution = 0.8,
  harmony2_metric_sample_size = 24L,
  harmony2_max_plot_cells = 24L,
  seed = 42L
)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
res <- run_compute_method(cfg, cfg$output_dir)
if (!identical(res$status, "ok")) {
  cat("[DEBUG] harmony2 result:\n")
  print(res)
  if (!is.null(res$command_result$output)) cat(paste(res$command_result$output, collapse = "\n"), "\n")
}
stopifnot(identical(res$status, "ok"))
stopifnot(file.exists(file.path(cfg$output_dir, "harmony2_summary.json")))
stopifnot(file.exists(file.path(cfg$output_dir, "harmony2_result.h5ad")))
stopifnot(file.exists(file.path(cfg$output_dir, "harmony2_metrics.tsv")))
stopifnot(file.exists(file.path(cfg$output_dir, "harmony2_umap_overview.png")))

llm_output_dir <- file.path(cfg$lineage_dir, "llm_parallel", "harmony2")
dir.create(llm_output_dir, recursive = TRUE, showWarnings = FALSE)
llm_cfg <- modifyList(cfg, list(method = "llm_harmony2", output_dir = llm_output_dir, llm_enable_live = FALSE))
llm_res <- run_llm_summary_worker(llm_cfg, llm_output_dir)
stopifnot(identical(llm_res$status, "ok"))
stopifnot(nrow(llm_res$index) > 0L)
stopifnot(any(basename(llm_res$index$path) == "harmony2_metrics.tsv"))
stopifnot(any(basename(llm_res$index$path) == "harmony2_umap_overview.png"))
stopifnot(file.exists(file.path(llm_output_dir, "harmony2_LLM_prompt.md")))

cat("[OK] harmony2 worker and LLM discovery smoke test passed.\n")
