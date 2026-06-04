#!/usr/bin/env Rscript

script_path <- "/home/h2048/script/R/launch_cnmf_l2_unit_llm_20260521.R"
source(script_path)

tmp_root <- tempfile("cnmf_l2_unit_llm_")
source_dir <- file.path(tmp_root, "bcell", "cnmf_by_celltype_l2", "Memory_B", "cnmf_full")
dir.create(file.path(source_dir, "gep_gene_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(source_dir, "visualizations"), recursive = TRUE, showWarnings = FALSE)

writeLines(
  '{"status":"ok","lineage":"bcell","celltype_l2":"Memory_B","safe_celltype":"Memory_B"}',
  file.path(source_dir, "cnmf_l2_status.json")
)
utils::write.table(
  data.frame(
    k = c(3L, 3L, 3L, 3L),
    gep = c("GEP_1", "GEP_1", "GEP_2", "GEP_2"),
    rank = c(1L, 2L, 1L, 2L),
    gene = c("MS4A1", "CD79A", "CD74", "HLA-DRA"),
    score = c(0.90, 0.80, 0.70, 0.60)
  ),
  file = file.path(source_dir, "gep_gene_tables", "gep_gene_scores_k3.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.csv(
  data.frame(
    cell = c("cell_1", "cell_2"),
    cell_type_L1 = c("Lymphoid", "Lymphoid"),
    cell_type_L2 = c("Memory_B", "Memory_B"),
    tissue = c("nasal", "polyp")
  ),
  file.path(source_dir, "cell_metadata_for_l2_cnmf.csv"),
  row.names = FALSE
)
writeLines('{"recommended_k":3}', file.path(source_dir, "k_selection_recommendation.json"))
writeLines("png", file.path(source_dir, "visualizations", "gep_top_gene_barplots_k3.png"))

idx <- build_l2_cnmf_unit_rows(run_root = tmp_root, lineages = "bcell", all_cnmf_k = TRUE)

stopifnot(is.data.frame(idx))
stopifnot(nrow(idx) == 2L)

unit <- idx[idx$unit_id == "cnmf_l2_Memory_B_k3_GEP_1", , drop = FALSE]
stopifnot(nrow(unit) == 1L)
stopifnot(identical(as.character(unit$unit_label)[1], "Memory_B | k3 / GEP_1"))

extra <- parse_task_extra_unit_llm(unit[1, , drop = FALSE])
stopifnot(identical(as.character(extra$celltype_l2)[1], "Memory_B"))
stopifnot(identical(as.character(extra$celltype_l1)[1], "Lymphoid"))

prompt <- build_unit_prompt(unit[1, , drop = FALSE])
stopifnot(any(grepl("celltype_hierarchy: `L2=Memory_B`", prompt, fixed = TRUE)))
stopifnot(any(grepl("celltype_l2: `Memory_B`", prompt, fixed = TRUE)))

cat("[OK] L2 cNMF unit index carries L2 hierarchy into per-GEP prompt metadata.\n")
