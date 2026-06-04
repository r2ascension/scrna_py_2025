#!/usr/bin/env Rscript

source("/home/h2048/script/R/launch_cnmf_l3_unit_llm_20260519.R")

tmp_root <- tempfile("cnmf_l3_combine_")
source_dir <- file.path(tmp_root, "bcell", "cnmf_by_celltype", "Plasma_IgA", "cnmf_full")
dir.create(file.path(source_dir, "gep_gene_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_dir, "visualizations"), recursive = TRUE, showWarnings = FALSE)

writeLines(
  '{"status":"ok","lineage":"bcell","celltype_l2":"Plasma","celltype_l3":"Plasma_IgA","safe_celltype":"Plasma_IgA"}',
  file.path(source_dir, "cnmf_l3_status.json")
)
utils::write.table(
  data.frame(
    k = c(4L, 4L, 4L, 4L),
    gep = c("GEP_1", "GEP_1", "GEP_2", "GEP_2"),
    rank = c(1L, 2L, 1L, 2L),
    gene = c("IGHA1", "JCHAIN", "MKI67", "TOP2A"),
    score = c(0.94, 0.86, 0.72, 0.61)
  ),
  file = file.path(source_dir, "gep_gene_tables", "gep_gene_scores_k4.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.csv(
  data.frame(cell = c("cell_1", "cell_2"), cell_type_L1 = c("Lymphoid", "Lymphoid"), cell_type_L2 = c("Plasma", "Plasma"), cell_type_L3 = c("Plasma_IgA", "Plasma_IgA")),
  file.path(source_dir, "cell_metadata_for_l3_cnmf.csv"),
  row.names = FALSE
)
writeLines('{"recommended_k":4}', file.path(source_dir, "k_selection_recommendation.json"))
writeLines("png", file.path(source_dir, "visualizations", "gep_top_gene_barplots_k4.png"))

idx <- build_l3_cnmf_unit_rows(run_root = tmp_root, lineages = "bcell", all_cnmf_k = TRUE)
stopifnot(nrow(idx) == 2L)

for (i in seq_len(nrow(idx))) {
  task <- idx[i, , drop = FALSE]
  write_json_unit_llm(list(status = "ok"), task$status_json[[1]])
  write_lines_unit_llm(c(sprintf("# %s", task$unit_id[[1]]), "", sprintf("Interpretation for %s.", task$unit_label[[1]])), task$interpretation_md[[1]])
}

combined <- combine_cnmf_unit_llm_by_celltype(
  idx,
  run_root = tmp_root,
  output_subdir = "cnmf_by_celltype",
  level_tag = "l3",
  celltype_field = "celltype_l3"
)

stopifnot(is.data.frame(combined$master_manifest_df))
stopifnot(nrow(combined$master_manifest_df) == 1L)

combined_md <- file.path(tmp_root, "bcell", "cnmf_by_celltype", "llm_parallel", "cnmf", "celltypes", "Plasma_IgA", "cnmf_l3_Plasma_IgA_LLM_combined.md")
stopifnot(file.exists(combined_md))
combined_lines <- readLines(combined_md, warn = FALSE)
stopifnot(any(grepl("# bcell / L3 / Plasma_IgA cNMF LLM 合并汇总", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("celltype_l2: `Plasma`", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("celltype_l3: `Plasma_IgA`", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("## Plasma > Plasma_IgA | k4 / GEP_1", combined_lines, fixed = TRUE)))
stopifnot(any(grepl("## Plasma > Plasma_IgA | k4 / GEP_2", combined_lines, fixed = TRUE)))

cat("[OK] L3 cNMF combined celltype LLM report bundles all unit interpretations for the same L3 cell type.\n")