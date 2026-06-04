#!/usr/bin/env Rscript

unit_llm_path <- "/home/h2048/script/R/launch_program_unit_llm_20260513.R"
source(unit_llm_path)

tmp_root <- tempfile("program_unit_llm_")
lineage_dir <- file.path(tmp_root, "bcell")
dir.create(lineage_dir, recursive = TRUE, showWarnings = FALSE)

cnmf_dir <- file.path(lineage_dir, "cnmf_full")
dir.create(file.path(cnmf_dir, "gep_gene_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cnmf_dir, "visualizations"), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  data.frame(
    k = c(3L, 3L, 3L, 3L),
    gep = c("GEP_1", "GEP_1", "GEP_2", "GEP_2"),
    rank = c(1L, 2L, 1L, 2L),
    gene = c("MS4A1", "CD79A", "CD3D", "IL7R"),
    score = c(0.9, 0.8, 0.7, 0.6)
  ),
  file = file.path(cnmf_dir, "gep_gene_tables", "gep_gene_scores_k3.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  data.frame(
    k = c(5L, 5L, 5L, 5L, 5L, 5L),
    gep = c("GEP_1", "GEP_1", "GEP_3", "GEP_3", "GEP_3", "GEP_3"),
    rank = c(1L, 2L, 1L, 2L, 3L, 4L),
    gene = c("IGHM", "CD74", "AL114490.2", "RP11-1A.1", "JCHAIN", "XBP1"),
    score = c(0.95, 0.75, 0.99, 0.97, 0.85, 0.65)
  ),
  file = file.path(cnmf_dir, "gep_gene_tables", "gep_gene_scores_k5.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.csv(
  data.frame(row = 1L, GEP_1 = "IGHM", GEP_3 = "AL114490.2"),
  file = file.path(cnmf_dir, "gep_gene_tables", "gep_top_genes_k5.csv"),
  row.names = FALSE
)
writeLines('{"recommended_k":3}', file.path(cnmf_dir, "k_selection_recommendation.json"))
writeLines("png", file.path(cnmf_dir, "visualizations", "gep_top_gene_barplots_k3.png"))
writeLines("png", file.path(cnmf_dir, "visualizations", "gep_top_gene_barplots_k5.png"))

pycogaps_dir <- file.path(lineage_dir, "pycogaps_full")
dir.create(pycogaps_dir, recursive = TRUE, showWarnings = FALSE)
writeLines('{"Pattern1":["MS4A1","CD79A"],"Pattern2":["CD3D","IL7R"]}', file.path(pycogaps_dir, "pycogaps_top_genes_by_pattern.json"))
utils::write.table(
  data.frame(gene = c("MS4A1", "CD79A", "CD3D", "IL7R"), Pattern1 = c(1, 0.8, 0.1, 0.2), Pattern2 = c(0.1, 0.2, 1, 0.9)),
  file = file.path(pycogaps_dir, "pycogaps_gene_patterns.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

covarnet_dir <- file.path(lineage_dir, "covarnet_full")
dir.create(covarnet_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(gene = c("MS4A1", "CD79A"), degree = c(10, 8), hub_score = c(0.9, 0.8)), file.path(covarnet_dir, "covarnet_Naive_B_nodes.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene_a = "MS4A1", gene_b = "CD79A", r = 0.8), file.path(covarnet_dir, "covarnet_Naive_B_edges.csv"), row.names = FALSE)

hdwgcna_dir <- file.path(lineage_dir, "hdwgcna_full", "hdwgcna_Naive_B")
dir.create(hdwgcna_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(gene_name = c("MS4A1", "CD79A", "CD3D"), module = c("blue", "blue", "grey"), color = c("blue", "blue", "grey"), kME_blue = c(0.9, 0.8, 0.1)),
  file.path(hdwgcna_dir, "hdwgcna_module_membership.csv"),
  row.names = FALSE
)
writeLines("pdf", file.path(hdwgcna_dir, "hdwgcna_module_trait_condition.pdf"))
writeLines("pdf", file.path(hdwgcna_dir, "hdwgcna_ME_heatmap_by_condition.pdf"))
utils::write.csv(data.frame(module = "blue", Healthy = 1), file.path(hdwgcna_dir, "hdwgcna_ME_mean_by_condition.csv"), row.names = FALSE)
utils::write.csv(data.frame(module = "blue", nasal = 1), file.path(hdwgcna_dir, "hdwgcna_ME_mean_by_tissue.csv"), row.names = FALSE)
writeLines("pdf", file.path(hdwgcna_dir, "hdwgcna_ME_heatmap_by_tissue.pdf"))
writeLines("pdf", file.path(hdwgcna_dir, "hdwgcna_hub_dotplot.pdf"))
writeLines("pdf", file.path(hdwgcna_dir, "hdwgcna_hub_rank_barplot.pdf"))
writeLines('{"status":"ok"}', file.path(hdwgcna_dir, "hdwgcna_celltype_status.json"))

idx <- build_unit_llm_task_index(
  run_root = tmp_root,
  lineages = "bcell",
  source_methods = c("cnmf", "pycogaps", "covarnet", "hdwgcna")
)

stopifnot(is.data.frame(idx))
stopifnot(nrow(idx) == 8L)
stopifnot(sum(idx$source_method == "cnmf" & idx$unit_type == "gep") == 4L)
stopifnot(sum(idx$source_method == "pycogaps" & idx$unit_type == "pattern") == 2L)
stopifnot(sum(idx$source_method == "covarnet" & idx$unit_type == "celltype_network") == 1L)
stopifnot(sum(idx$source_method == "hdwgcna" & idx$unit_type == "hdwgcna_module") == 1L)
stopifnot(any(idx$unit_id == "cnmf_k3_GEP_1"))
stopifnot(any(idx$unit_id == "cnmf_k5_GEP_3"))
stopifnot(any(idx$source_method == "cnmf" & idx$unit_id == "cnmf_k3_GEP_1" & idx$unit_label == "k3 / GEP_1"))
stopifnot(any(idx$source_method == "cnmf" & idx$unit_id == "cnmf_k5_GEP_3" & idx$unit_label == "k5 / GEP_3"))
cnmf_k5_gep3 <- idx[idx$unit_id == "cnmf_k5_GEP_3", , drop = FALSE]
stopifnot(nrow(cnmf_k5_gep3) == 1L)
unit_evidence <- strsplit(cnmf_k5_gep3$evidence_paths, ";", fixed = TRUE)[[1]]
unit_gene_file <- unit_evidence[grepl("cnmf_k5_GEP_3_gene_scores[.]csv$", unit_evidence)]
stopifnot(length(unit_gene_file) == 1L)
stopifnot(file.exists(unit_gene_file))
unit_gene_tbl <- utils::read.csv(unit_gene_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(unique(unit_gene_tbl$k), 5L))
stopifnot(identical(unique(unit_gene_tbl$gep), "GEP_3"))
stopifnot(!any(unit_gene_tbl$gene %in% c("AL114490.2", "RP11-1A.1")))
stopifnot(all(unit_gene_tbl$gene %in% c("JCHAIN", "XBP1")))
stopifnot(!any(grepl("gep_gene_scores_k5[.]tsv$|gep_top_genes_k5[.]csv$", unit_evidence)))
cnmf_prompt <- build_unit_prompt(cnmf_k5_gep3[1, , drop = FALSE])
stopifnot(any(grepl("当前 unit 是一个明确的 `k × GEP` 单元", cnmf_prompt, fixed = TRUE)))
stopifnot(any(idx$unit_id == "pycogaps_Pattern1"))
stopifnot(any(idx$unit_id == "covarnet_Naive_B"))
hdwgcna_blue <- idx[idx$unit_id == "hdwgcna_Naive_B_blue", , drop = FALSE]
stopifnot(nrow(hdwgcna_blue) == 1L)
hdwgcna_evidence <- strsplit(hdwgcna_blue$evidence_paths, ";", fixed = TRUE)[[1]]
stopifnot(!any(grepl("condition|disease|CRSwNP|Healthy|hub_dotplot", basename(hdwgcna_evidence), ignore.case = TRUE)))
stopifnot(any(grepl("tissue", basename(hdwgcna_evidence), ignore.case = TRUE)))
hdwgcna_prompt <- build_unit_prompt(hdwgcna_blue[1, , drop = FALSE])
stopifnot(any(grepl("所有输入细胞均按健康样本/健康组织背景解释", hdwgcna_prompt, fixed = TRUE)))
stopifnot(any(grepl("输出正文不得出现这些字面字符串", hdwgcna_prompt, fixed = TRUE)))
stopifnot(all(grepl("llm_parallel", idx$output_dir, fixed = TRUE)))

cat("[OK] unit-level LLM task index splits cNMF GEPs, PyCoGAPS patterns, and CoVarNet celltype networks.\n")