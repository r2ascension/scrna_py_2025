#!/usr/bin/env Rscript

worker_path <- "/home/h2048/script/R/program_lineage_method_worker_20260507.R"
source(worker_path)

pa_source_coexpr_helper()
stopifnot(exists("hdwgcna_plot_summary", mode = "function", inherits = TRUE))
stopifnot(exists("covarnet_plot_all", mode = "function", inherits = TRUE))

tmp_root <- tempfile("pa_llm_discovery_")
lineage_dir <- file.path(tmp_root, "test_lineage")
source_dir <- file.path(lineage_dir, "covarnet_full")
output_dir <- file.path(lineage_dir, "llm_parallel", "covarnet")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(
  data.frame(celltype = "T_cell", n_nodes = 100L, n_edges = 250L, mean_degree = 5.0),
  file.path(source_dir, "covarnet_network_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(celltype = "T_cell", direction = "positive", n_edges = 180L),
  file.path(source_dir, "covarnet_edge_direction_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(celltype = "T_cell", gene = "CD3D", degree = 18L, hub_score = 0.9),
  file.path(source_dir, "covarnet_node_degree_all.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(gene = "MALAT1", reason = "lncRNA"),
  file.path(source_dir, "gene_exclusion_summary.csv"),
  row.names = FALSE
)
writeLines("# Gene-exclusion sidecar should not be primary evidence", file.path(source_dir, "gene_exclusion_LLM_prompt.md"))
writeLines("# Gene-exclusion interpretation should not be primary evidence", file.path(source_dir, "gene_exclusion_LLM_interpretation.md"))
utils::write.csv(
  data.frame(gene_a = "CD3D", gene_b = "CD247", r = 0.81, direction = "positive", celltype = "T_cell"),
  file.path(source_dir, "covarnet_T_cell_edges.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(gene = "CD3D", degree = 18L, hub_score = 0.9, celltype = "T_cell"),
  file.path(source_dir, "covarnet_T_cell_nodes.csv"),
  row.names = FALSE
)

idx <- discover_llm_source_files(source_dir)
stopifnot(nrow(idx) >= 4L)
stopifnot(identical(idx$role[[1]], "visualization_evidence"))
stopifnot(!any(basename(idx$path) %in% c("gene_exclusion_LLM_prompt.md", "gene_exclusion_LLM_interpretation.md")))
stopifnot(any(idx$role == "qc_context" & basename(idx$path) == "gene_exclusion_summary.csv"))

legacy_covarnet_dir <- file.path(tmp_root, "legacy_covarnet_full")
dir.create(legacy_covarnet_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(gene_a = "MS4A1", gene_b = "CD79A", r = 0.74, direction = "positive", celltype = "B_cell"),
  file.path(legacy_covarnet_dir, "covarnet_B_cell_edges.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(gene = "MS4A1", degree = 21L, hub_score = 0.95, celltype = "B_cell"),
  file.path(legacy_covarnet_dir, "covarnet_B_cell_nodes.csv"),
  row.names = FALSE
)
writeLines("# Legacy gene-exclusion sidecar", file.path(legacy_covarnet_dir, "gene_exclusion_LLM_interpretation.md"))
legacy_cov_idx <- discover_llm_source_files(legacy_covarnet_dir)
stopifnot(nrow(legacy_cov_idx) == 2L)
stopifnot(all(legacy_cov_idx$role == "program_evidence"))
stopifnot(all(basename(legacy_cov_idx$path) %in% c("covarnet_B_cell_edges.csv", "covarnet_B_cell_nodes.csv")))

legacy_hdwgcna_dir <- file.path(tmp_root, "legacy_hdwgcna_full", "hdwgcna_B_cell")
dir.create(legacy_hdwgcna_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(gene_name = "MS4A1", module = "blue", kME = 0.83),
  file.path(legacy_hdwgcna_dir, "hdwgcna_module_membership.csv"),
  row.names = FALSE
)
legacy_hdw_idx <- discover_llm_source_files(dirname(legacy_hdwgcna_dir))
stopifnot(nrow(legacy_hdw_idx) == 1L)
stopifnot(identical(legacy_hdw_idx$role[[1]], "program_evidence"))
stopifnot(identical(basename(legacy_hdw_idx$path[[1]]), "hdwgcna_module_membership.csv"))

cfg <- list(
  lineage = "test_lineage",
  method = "llm_covarnet",
  lineage_dir = lineage_dir,
  output_dir = output_dir,
  status_tsv = file.path(tmp_root, "status.tsv"),
  run_stamp = "20260509_test",
  llm_enable_live = FALSE
)
writeLines("# Prior CoVarNet interpretation\n\n上一轮解释应作为 continuity context。", file.path(output_dir, "covarnet_LLM_interpretation.md"))
writeLines("role\text\tpath\trel_path\tsize_bytes\nprogram_evidence\tcsv\told.csv\told.csv\t10", file.path(output_dir, "covarnet_llm_input_index.tsv"))

res <- run_llm_summary_worker(cfg, output_dir)
stopifnot(identical(res$status, "ok"))
stopifnot(any(res$index$role == "prior_llm_context"))
stopifnot(any(res$index$role == "program_evidence" & basename(res$index$path) == "covarnet_T_cell_edges.csv"))

prompt_path <- file.path(output_dir, "covarnet_LLM_prompt.md")
prompt <- readLines(prompt_path, warn = FALSE)
stopifnot(any(grepl("covarnet_network_summary.csv", prompt, fixed = TRUE)))
stopifnot(any(grepl("Prior CoVarNet interpretation", prompt, fixed = TRUE)))
stopifnot(any(grepl("visualization_evidence", prompt, fixed = TRUE)))
stopifnot(!any(grepl("gene_exclusion_LLM_prompt", prompt, fixed = TRUE)))

cat("[OK] LLM source discovery prioritizes visualization CSVs and excludes gene-exclusion LLM sidecars.\n")