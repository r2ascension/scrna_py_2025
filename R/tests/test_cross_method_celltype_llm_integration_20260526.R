#!/usr/bin/env Rscript

source("/home/h2048/script/R/launch_cross_method_celltype_llm_20260526.R")

tmp_root <- tempfile("cross_method_llm_")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

write_stub_report <- function(path, title, body_lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(sprintf("# %s", title), "", body_lines), path)
}

path_l2_md <- file.path(tmp_root, "bcell", "cnmf_by_celltype_l2", "llm_parallel", "cnmf", "celltypes", "Plasma", "cnmf_l2_Plasma_LLM_combined.md")
path_l2_tsv <- file.path(tmp_root, "bcell", "cnmf_by_celltype_l2", "llm_parallel", "cnmf", "celltypes", "Plasma", "cnmf_l2_Plasma_LLM_combined_manifest.tsv")
path_l2_json <- file.path(tmp_root, "bcell", "cnmf_by_celltype_l2", "llm_parallel", "cnmf", "celltypes", "Plasma", "cnmf_l2_Plasma_LLM_combined_manifest.json")
write_stub_report(path_l2_md, "bcell / L2 / Plasma cNMF LLM 合并汇总", c("L2 plasma context line."))
writeLines("colA\nvalue", path_l2_tsv)
writeLines("{}", path_l2_json)

path_l3_md <- file.path(tmp_root, "bcell", "cnmf_by_celltype", "llm_parallel", "cnmf", "celltypes", "Plasma_IgA", "cnmf_l3_Plasma_IgA_LLM_combined.md")
path_l3_tsv <- file.path(tmp_root, "bcell", "cnmf_by_celltype", "llm_parallel", "cnmf", "celltypes", "Plasma_IgA", "cnmf_l3_Plasma_IgA_LLM_combined_manifest.tsv")
path_l3_json <- file.path(tmp_root, "bcell", "cnmf_by_celltype", "llm_parallel", "cnmf", "celltypes", "Plasma_IgA", "cnmf_l3_Plasma_IgA_LLM_combined_manifest.json")
write_stub_report(path_l3_md, "bcell / L3 / Plasma_IgA cNMF LLM 合并汇总", c("L3 plasma IgA exact line."))
writeLines("colA\nvalue", path_l3_tsv)
writeLines("{}", path_l3_json)

cov_status <- file.path(tmp_root, "bcell", "llm_parallel", "covarnet", "units", "covarnet_Plasma_IgA", "covarnet_Plasma_IgA_LLM_status.json")
cov_md <- file.path(tmp_root, "bcell", "llm_parallel", "covarnet", "units", "covarnet_Plasma_IgA", "covarnet_Plasma_IgA_LLM_interpretation.md")
dir.create(dirname(cov_status), recursive = TRUE, showWarnings = FALSE)
writeLines('{"status":"ok"}', cov_status)
write_stub_report(cov_md, "covarnet_Plasma_IgA", c("CoVarNet exact plasma IgA line."))

hd_status <- file.path(tmp_root, "bcell", "llm_parallel", "hdwgcna", "units", "hdwgcna_Plasma_IgA_turquoise", "hdwgcna_Plasma_IgA_turquoise_LLM_status.json")
hd_md <- file.path(tmp_root, "bcell", "llm_parallel", "hdwgcna", "units", "hdwgcna_Plasma_IgA_turquoise", "hdwgcna_Plasma_IgA_turquoise_LLM_interpretation.md")
dir.create(dirname(hd_status), recursive = TRUE, showWarnings = FALSE)
writeLines('{"status":"ok"}', hd_status)
write_stub_report(hd_md, "hdwgcna_Plasma_IgA_turquoise", c("hdWGCNA turquoise module line."))

py_status <- file.path(tmp_root, "bcell", "llm_parallel", "pycogaps", "pycogaps_LLM_status.json")
py_md <- file.path(tmp_root, "bcell", "llm_parallel", "pycogaps", "pycogaps_LLM_interpretation.md")
py_idx <- file.path(tmp_root, "bcell", "llm_parallel", "pycogaps", "pycogaps_llm_input_index.tsv")
dir.create(dirname(py_status), recursive = TRUE, showWarnings = FALSE)
writeLines('{"status":"ok"}', py_status)
write_stub_report(py_md, "pycogaps global", c("Global pyCoGAPS context line."))
writeLines("idx\n1", py_idx)

l2_records <- data.frame(
  lineage = "bcell",
  level_tag = "l2",
  celltype_field = "celltype_l2",
  celltype_value = "Plasma",
  safe_celltype = "Plasma",
  celltype_l1 = "Lymphoid",
  celltype_l2 = "Plasma",
  celltype_l3 = NA_character_,
  total_units = 3L,
  ok_units = 3L,
  error_units = 0L,
  group_status = "ok",
  combined_md = path_l2_md,
  combined_manifest_tsv = path_l2_tsv,
  combined_manifest_json = path_l2_json,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

l3_records <- data.frame(
  lineage = "bcell",
  level_tag = "l3",
  celltype_field = "celltype_l3",
  celltype_value = "Plasma_IgA",
  safe_celltype = "Plasma_IgA",
  celltype_l1 = "Lymphoid",
  celltype_l2 = "Plasma",
  celltype_l3 = "Plasma_IgA",
  total_units = 2L,
  ok_units = 2L,
  error_units = 0L,
  group_status = "ok",
  combined_md = path_l3_md,
  combined_manifest_tsv = path_l3_tsv,
  combined_manifest_json = path_l3_json,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

root_unit_records <- data.frame(
  lineage = c("bcell", "bcell"),
  source_method = c("covarnet", "hdwgcna"),
  unit_type = c("celltype_network", "hdwgcna_module"),
  unit_id = c("covarnet_Plasma_IgA", "hdwgcna_Plasma_IgA_turquoise"),
  unit_label = c("Plasma_IgA", "Plasma_IgA / turquoise"),
  celltype = c("Plasma_IgA", "Plasma_IgA"),
  safe_celltype = c("Plasma_IgA", "Plasma_IgA"),
  module = c(NA_character_, "turquoise"),
  pattern = c(NA_character_, NA_character_),
  report_md = c(cov_md, hd_md),
  interpretation_md = c(cov_md, hd_md),
  status_json = c(cov_status, hd_status),
  status = c("ok", "ok"),
  error = c(NA_character_, NA_character_),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

global_reports <- data.frame(
  lineage = "bcell",
  source_method = "pycogaps",
  report_md = py_md,
  status_json = py_status,
  status = "ok",
  input_index_tsv = py_idx,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

combined <- write_cross_method_celltype_dossiers(
  run_root = tmp_root,
  l2_records = l2_records,
  l3_records = l3_records,
  root_unit_records = root_unit_records,
  global_reports = global_reports,
  output_subdir = "llm_cross_method"
)

stopifnot(is.data.frame(combined$master_manifest_df))
stopifnot(nrow(combined$master_manifest_df) == 2L)
stopifnot(file.exists(combined$master_manifest_tsv))
stopifnot(file.exists(combined$master_manifest_json))

plasma_iga_md <- file.path(tmp_root, "bcell", "llm_cross_method", "celltypes", "Plasma_IgA", "cross_method_Plasma_IgA_LLM_dossier.md")
plasma_md <- file.path(tmp_root, "bcell", "llm_cross_method", "celltypes", "Plasma", "cross_method_Plasma_LLM_dossier.md")
stopifnot(file.exists(plasma_iga_md))
stopifnot(file.exists(plasma_md))

plasma_iga_lines <- readLines(plasma_iga_md, warn = FALSE)
stopifnot(any(grepl("# bcell / Plasma_IgA 跨方法 LLM 总整合", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("## cNMF L3：Plasma_IgA", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("L3 plasma IgA exact line.", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("## cNMF L2 背景：Plasma", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("L2 plasma context line.", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("## CoVarNet：Plasma_IgA", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("CoVarNet exact plasma IgA line.", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("## hdWGCNA：Plasma_IgA / turquoise", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("hdWGCNA turquoise module line.", plasma_iga_lines, fixed = TRUE)))
stopifnot(any(grepl("pycogaps", plasma_iga_lines, fixed = TRUE)))
stopifnot(!any(grepl("Global pyCoGAPS context line.", plasma_iga_lines, fixed = TRUE)))

plasma_lines <- readLines(plasma_md, warn = FALSE)
stopifnot(any(grepl("# bcell / Plasma 跨方法 LLM 总整合", plasma_lines, fixed = TRUE)))
stopifnot(any(grepl("## 关联 L3 子类型", plasma_lines, fixed = TRUE)))
stopifnot(any(grepl("Plasma_IgA", plasma_lines, fixed = TRUE)))

task_df <- data.frame(
  lineage = "bcell",
  source_method = "hdwgcna",
  unit_type = "hdwgcna_module",
  unit_id = "hdwgcna_Memory_B_turquoise",
  unit_label = "Memory_B / turquoise",
  source_dir = "/tmp/source",
  output_dir = "/tmp/out",
  evidence_paths = "",
  payload_json = "/tmp/out/payload.json",
  prompt_md = "/tmp/out/prompt.md",
  interpretation_md = hd_md,
  status_json = hd_status,
  extra_json = '{"celltype":"Memory_B","module":"turquoise"}',
  stringsAsFactors = FALSE,
  check.names = FALSE
)
task_records <- unit_llm_index_records(task_df)
stopifnot(identical(task_records$celltype[[1]], "Memory_B"))
stopifnot(identical(task_records$safe_celltype[[1]], "Memory_B"))
stopifnot(identical(task_records$module[[1]], "turquoise"))

cat("[OK] cross-method celltype dossier bundles exact-method evidence, parent L2 context, and lineage-global references.\n")