#!/usr/bin/env Rscript

source("/home/h2048/script/R/launch_cnmf_l2_unit_llm_20260521.R")

tmp_root <- tempfile("cnmf_l2_combine_")
source_dir <- file.path(tmp_root, "bcell", "cnmf_by_celltype_l2", "Plasma", "cnmf_full")
dir.create(file.path(source_dir, "gep_gene_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_dir, "visualizations"), recursive = TRUE, showWarnings = FALSE)

writeLines(
  '{"status":"ok","lineage":"bcell","celltype_l2":"Plasma","safe_celltype":"Plasma"}',
  file.path(source_dir, "cnmf_l2_status.json")
)
utils::write.table(
  data.frame(
    k = c(3L, 3L, 3L, 3L, 5L, 5L),
    gep = c("GEP_1", "GEP_1", "GEP_2", "GEP_2", "GEP_1", "GEP_1"),
    rank = c(1L, 2L, 1L, 2L, 1L, 2L),
    gene = c("JCHAIN", "XBP1", "MKI67", "TOP2A", "IGHG1", "MZB1"),
    score = c(0.92, 0.81, 0.75, 0.63, 0.90, 0.78)
  ),
  file = file.path(source_dir, "gep_gene_tables", "gep_gene_scores_k3.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  data.frame(
    k = c(5L, 5L),
    gep = c("GEP_1", "GEP_1"),
    rank = c(1L, 2L),
    gene = c("IGHG1", "MZB1"),
    score = c(0.90, 0.78)
  ),
  file = file.path(source_dir, "gep_gene_tables", "gep_gene_scores_k5.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.csv(
  data.frame(cell = c("cell_1", "cell_2"), cell_type_L1 = c("Lymphoid", "Lymphoid"), cell_type_L2 = c("Plasma", "Plasma")),
  file.path(source_dir, "cell_metadata_for_l2_cnmf.csv"),
  row.names = FALSE
)
writeLines('{"recommended_k":3}', file.path(source_dir, "k_selection_recommendation.json"))
writeLines("png", file.path(source_dir, "visualizations", "gep_top_gene_barplots_k3.png"))
writeLines("png", file.path(source_dir, "visualizations", "gep_top_gene_barplots_k5.png"))

idx <- build_l2_cnmf_unit_rows(run_root = tmp_root, lineages = "bcell", all_cnmf_k = TRUE)
stopifnot(nrow(idx) == 3L)

for (i in seq_len(nrow(idx))) {
  task <- idx[i, , drop = FALSE]
  write_json_unit_llm(list(status = if (grepl("GEP_2$", task$unit_id[[1]])) "error" else "ok", error = if (grepl("GEP_2$", task$unit_id[[1]])) "mock failure" else NULL), task$status_json[[1]])
  write_lines_unit_llm(
    c(
      sprintf("# %s", task$unit_id[[1]]),
      "",
      sprintf("Interpretation for %s.", task$unit_label[[1]])
    ),
    task$interpretation_md[[1]]
  )
}

combined <- combine_cnmf_unit_llm_by_celltype(
  idx,
  run_root = tmp_root,
  output_subdir = "cnmf_by_celltype_l2",
  level_tag = "l2",
  celltype_field = "celltype_l2"
)

stopifnot(is.data.frame(combined$master_manifest_df))
stopifnot(nrow(combined$master_manifest_df) == 1L)
stopifnot(file.exists(combined$master_manifest_tsv))

combined_md <- file.path(tmp_root, "bcell", "cnmf_by_celltype_l2", "llm_parallel", "cnmf", "celltypes", "Plasma", "cnmf_l2_Plasma_LLM_combined.md")
combined_group_tsv <- file.path(tmp_root, "bcell", "cnmf_by_celltype_l2", "llm_parallel", "cnmf", "celltypes", "Plasma", "cnmf_l2_Plasma_LLM_combined_manifest.tsv")
stopifnot(file.exists(combined_md))
stopifnot(file.exists(combined_group_tsv))

combined_lines <- readLines(combined_md, warn = FALSE)
stopifnot(any(grepl("# bcell / L2 / Plasma cNMF LLM 合并汇总", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("## Plasma | k3 / GEP_1", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("## Plasma | k3 / GEP_2", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("## Plasma | k5 / GEP_1", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("Interpretation for Plasma | k3 / GEP_1.", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("mock failure", combined_lines, fixed = TRUE)))
stopifnot(!any(grepl("# cnmf_l2_Plasma_k3_GEP_1", combined_lines, fixed = TRUE)))

manifest_df <- utils::read.delim(combined_group_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(manifest_df) == 3L)
stopifnot(identical(as.character(manifest_df$unit_id), c("cnmf_l2_Plasma_k3_GEP_1", "cnmf_l2_Plasma_k3_GEP_2", "cnmf_l2_Plasma_k5_GEP_1")))

cat("[OK] L2 cNMF combined celltype LLM report bundles all unit interpretations for the same L2 cell type.\n")