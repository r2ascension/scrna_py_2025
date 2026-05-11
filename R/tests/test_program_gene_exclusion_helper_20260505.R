#!/usr/bin/env Rscript

source('/home/h2048/script/R/program_gene_exclusion_helper_20260505_v1.R')

suppressPackageStartupMessages({
  library(Seurat)
})

old_deepseek_key <- Sys.getenv('DEEPSEEK_API_KEY', unset = NA_character_)
on.exit({
  if (is.na(old_deepseek_key)) {
    Sys.unsetenv('DEEPSEEK_API_KEY')
  } else {
    Sys.setenv(DEEPSEEK_API_KEY = old_deepseek_key)
  }
}, add = TRUE)
Sys.setenv(DEEPSEEK_API_KEY = 'your_deepseek_api_key_here')

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

packet_non_b <- pa_build_gene_exclusion_packet(
  feature_names = c('IGHG1', 'MT-CO1', 'RPS12', 'LINC00152', 'MS4A1'),
  lineage_context = 'myeloid'
)
prompt_non_b <- pa_build_gene_exclusion_llm_prompt(packet_non_b)
assert_true(grepl('Category summary', prompt_non_b$user_prompt, fixed = TRUE), 'LLM prompt should include category summary')
assert_true(grepl('single-cell RNA-seq gene-program QC expert', prompt_non_b$system_prompt, fixed = TRUE), 'LLM system prompt should define QC expert role')
assert_true(packet_non_b$audit_table$should_exclude[packet_non_b$audit_table$gene_symbol == 'IGHG1'], 'IGHG1 should be excluded outside B/plasma lineage')
assert_true(packet_non_b$audit_table$should_exclude[packet_non_b$audit_table$gene_symbol == 'MT-CO1'], 'MT-CO1 should be excluded')
assert_true(packet_non_b$audit_table$should_exclude[packet_non_b$audit_table$gene_symbol == 'RPS12'], 'RPS12 should be excluded')
assert_true(packet_non_b$audit_table$should_exclude[packet_non_b$audit_table$gene_symbol == 'LINC00152'], 'LINC00152 should be excluded by lncRNA heuristic')
assert_true(!packet_non_b$audit_table$should_exclude[packet_non_b$audit_table$gene_symbol == 'MS4A1'], 'MS4A1 should be retained')

packet_b <- pa_build_gene_exclusion_packet(
  feature_names = c('IGHG1', 'MS4A1'),
  lineage_context = 'Bcell'
)
assert_true(!packet_b$audit_table$should_exclude[packet_b$audit_table$gene_symbol == 'IGHG1'], 'IGHG1 should be retained for B/plasma lineage context')

counts <- matrix(
  c(5, 3, 4, 6, 8,
    2, 1, 2, 1, 3,
    7, 8, 6, 7, 6),
  nrow = 5,
  dimnames = list(c('IGHG1', 'MT-CO1', 'RPS12', 'LINC00152', 'MS4A1'), paste0('cell_', 1:3))
)
seu <- CreateSeuratObject(counts = counts)
seu <- NormalizeData(seu, verbose = FALSE)
seu$cell_type_final_l3 <- c('Mono', 'Mono', 'Mono')

out_dir <- file.path(tempdir(), 'program_gene_exclusion_helper_test')
res <- pa_apply_gene_exclusion_to_seurat(
  seurat_obj = seu,
  output_dir = out_dir,
  gene_exclusion_config = list(
    lineage_context = 'myeloid',
    llm_config = list(enabled = TRUE, model = 'deepseek-chat')
  )
)

assert_true(nrow(res$seurat_obj) == 1L, 'Only one non-excluded feature should remain in myeloid test object')
assert_true(identical(rownames(res$seurat_obj), 'MS4A1'), 'Filtered Seurat object should retain only MS4A1')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_summary.csv')), 'Summary CSV should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_audit.csv')), 'Audit CSV should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_summary.png')), 'Summary PNG should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_dashboard.png')), 'Dashboard PNG should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_dashboard.pdf')), 'Dashboard PDF should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_visualization_manifest.json')), 'Visualization manifest should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_LLM_prompt.md')), 'LLM prompt markdown should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_LLM_interpretation.md')), 'LLM interpretation markdown should be written even without live key')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_LLM_status.json')), 'LLM status JSON should be written')
assert_true(file.exists(file.path(out_dir, 'gene_exclusion_manifest.json')), 'Manifest JSON should be written')
assert_true(identical(res$manifest$llm$status, 'skipped_no_live_key'), 'LLM should skip safely when key is missing or placeholder')
assert_true(file.exists(res$manifest$visualizations$dashboard_plot$png), 'Dashboard PNG path should be recorded in manifest')

llm_md <- paste(readLines(file.path(out_dir, 'gene_exclusion_LLM_interpretation.md'), warn = FALSE), collapse = '\n')
assert_true(grepl('LLM status', llm_md, fixed = TRUE), 'LLM interpretation should record status')

cat('program gene exclusion helper tests passed.\n')
