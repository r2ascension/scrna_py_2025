#!/usr/bin/env Rscript

source("/home/h2048/script/R/launch_cnmf_l3_unit_llm_20260519.R")

tmp_root <- tempfile("cnmf_l3_paths_")
lineage_dir <- file.path(tmp_root, "bcell", "cnmf_by_celltype")
dir.create(lineage_dir, recursive = TRUE, showWarnings = FALSE)

write_cnmf_fixture <- function(celltype_root, celltype_l3, safe_celltype, k, gep, status = "ok") {
  source_dir <- file.path(celltype_root, "cnmf_full")
  dir.create(file.path(source_dir, "gep_gene_tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(source_dir, "visualizations"), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data.frame(
      k = c(k, k),
      gep = c(gep, gep),
      rank = c(1L, 2L),
      gene = c("JCHAIN", "XBP1"),
      score = c(0.9, 0.8)
    ),
    file = file.path(source_dir, "gep_gene_tables", sprintf("gep_gene_scores_k%d.tsv", k)),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  utils::write.table(
    data.frame(cell = c("cell1", "cell2"), cell_type_L2 = c("Plasma", "Plasma"), cell_type_L3 = c(celltype_l3, celltype_l3)),
    file = file.path(source_dir, "cell_metadata_for_l3_cnmf.csv"),
    sep = ",",
    row.names = FALSE,
    quote = TRUE
  )
  writeLines(sprintf('{"status":"%s","lineage":"bcell","celltype_l3":"%s","safe_celltype":"%s","path_layout":"celltype_method"}', status, celltype_l3, safe_celltype), file.path(source_dir, "cnmf_l3_status.json"))
  writeLines(sprintf('{"recommended_k":%d}', k), file.path(source_dir, "k_selection_recommendation.json"))
  writeLines("png", file.path(source_dir, "visualizations", sprintf("gep_top_gene_barplots_k%d.png", k)))
  source_dir
}

canonical_root <- file.path(lineage_dir, "Memory_B")
legacy_root_same <- file.path(lineage_dir, "cnmf_Memory_B")
legacy_root_only <- file.path(lineage_dir, "cnmf_Plasma_IgA")

canonical_source <- write_cnmf_fixture(canonical_root, "Memory B", "Memory_B", 3L, "GEP_1")
legacy_same_source <- write_cnmf_fixture(legacy_root_same, "Memory B", "Memory_B", 9L, "GEP_9")
writeLines('{"status":"ok","lineage":"bcell","celltype_l3":"Memory B","safe_celltype":"Memory_B"}', file.path(legacy_same_source, "cnmf_l3_status.json"))
legacy_only_source <- write_cnmf_fixture(legacy_root_only, "Plasma IgA", "Plasma_IgA", 5L, "GEP_2")
writeLines('{"status":"ok","lineage":"bcell","celltype_l3":"Plasma IgA","safe_celltype":"Plasma_IgA"}', file.path(legacy_only_source, "cnmf_l3_status.json"))

idx <- build_l3_cnmf_unit_rows(run_root = tmp_root, lineages = "bcell", all_cnmf_k = TRUE)

stopifnot(is.data.frame(idx))
stopifnot(nrow(idx) == 2L)
stopifnot(any(idx$unit_id == "cnmf_l3_Memory_B_k3_GEP_1"))
stopifnot(any(idx$unit_id == "cnmf_l3_Plasma_IgA_k5_GEP_2"))
stopifnot(!any(idx$unit_id == "cnmf_l3_Memory_B_k9_GEP_9"))

memory_row <- idx[idx$unit_id == "cnmf_l3_Memory_B_k3_GEP_1", , drop = FALSE]
plasma_row <- idx[idx$unit_id == "cnmf_l3_Plasma_IgA_k5_GEP_2", , drop = FALSE]

stopifnot(nrow(memory_row) == 1L)
stopifnot(nrow(plasma_row) == 1L)
stopifnot(identical(normalizePath(memory_row$source_dir[[1]], winslash = "/", mustWork = FALSE), normalizePath(canonical_source, winslash = "/", mustWork = FALSE)))
stopifnot(identical(normalizePath(plasma_row$source_dir[[1]], winslash = "/", mustWork = FALSE), normalizePath(legacy_only_source, winslash = "/", mustWork = FALSE)))
stopifnot(grepl("Plasma > Plasma IgA", plasma_row$unit_label[[1]], fixed = TRUE))

cat("[OK] L3 cNMF unit LLM index prefers canonical celltype/method paths and falls back to legacy directories.\n")