#!/usr/bin/env Rscript

bundle_path <- "/home/h2048/script/R/program_architecture_bundle_20260428_v1.R"
source(bundle_path)

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

assert_identical <- function(x, y, msg) {
  if (!identical(x, y)) stop(msg, call. = FALSE)
}

required_fns <- c(
  "pa_run_hdwgcna_runner",
  "pa_import_hdwgcna_results",
  "pa_run_cnmf_runner",
  "pa_import_cnmf_results",
  "pa_run_covarnet_runner",
  "pa_run_paga_runner",
  "pa_run_slingshot_runner",
  "pa_run_cytotrace2_runner",
  "pa_run_enrichment_runner",
  "pa_run_interpret_agent_runner",
  "pa_run_report_runner"
)

for (fn in required_fns) {
  assert_true(exists(fn, mode = "function"), sprintf("%s should be exported by the bundle", fn))
}

make_tiny_seurat <- function() {
  suppressPackageStartupMessages({
    library(Seurat)
  })
  set.seed(42)
  genes <- c("MS4A1","CD79A","CD79B","BANK1","HLA-DRA","CD74","MKI67","TOP2A")
  cells <- paste0("cell_", seq_len(12))
  counts <- matrix(rpois(length(genes) * length(cells), lambda = 5), nrow = length(genes), dimnames = list(genes, cells))
  counts[1:4, 1:6] <- counts[1:4, 1:6] + 5
  counts[5:8, 7:12] <- counts[5:8, 7:12] + 5
  seu <- CreateSeuratObject(counts = counts)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu$cell_type_final_l3 <- rep(c("Naive_B", "Memory_B"), each = 6)
  seu$sample <- rep(c("S1", "S2", "S3"), length.out = 12)
  seu$dataset <- rep(c("D1", "D2"), each = 6)
  seu$slingshot_cluster <- rep(c("C1", "C2", "C3"), each = 4)
  umap <- rbind(
    cbind(seq(0, 5), rep(0, 6)),
    cbind(seq(0, 5), rep(3, 6))
  )
  rownames(umap) <- colnames(seu)
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  seu[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = DefaultAssay(seu))
  seu
}

tmp_dir <- file.path(tempdir(), "pa_full_compute_runner")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# hdWGCNA runner: dry-run smoke
seu <- make_tiny_seurat()
hd_plan <- pa_run_hdwgcna_runner(
  seurat_obj = seu,
  output_dir = file.path(tmp_dir, "hdwgcna"),
  celltypes = "Naive_B",
  dry_run = TRUE
)
assert_identical(hd_plan$status, "dry_run", "hdWGCNA dry-run should return dry_run status")

# CoVarNet runner: real lightweight compute
covar_res <- pa_run_covarnet_runner(
  seurat_obj = seu,
  output_dir = file.path(tmp_dir, "covarnet"),
  celltypes = c("Naive_B", "Memory_B"),
  celltype_col = "cell_type_final_l3",
  cor_thr = 0,
  pval_thr = 1,
  min_cells = 1,
  n_genes = 6
)
assert_true(is.list(covar_res), "CoVarNet runner should return a list")
assert_true(is.data.frame(covar_res$program_tbl), "CoVarNet runner should return program_tbl")
assert_true(nrow(covar_res$program_tbl) >= 1, "CoVarNet runner should emit at least one program row")

# Slingshot runner: real lightweight compute
slingshot_res <- pa_run_slingshot_runner(
  seurat_obj = seu,
  output_dir = file.path(tmp_dir, "slingshot"),
  cluster_col = "slingshot_cluster",
  reduced_dim = "umap",
  start_cluster = "C1"
)
assert_identical(slingshot_res$trajectory_packet$engine, "Slingshot", "Slingshot runner should return Slingshot packet")
assert_true(file.exists(file.path(tmp_dir, "slingshot", "slingshot_pseudotime.tsv")), "Slingshot pseudotime file should be written")

# cNMF runner: dry-run + importer smoke
mock_cnmf_dir <- file.path(tmp_dir, "cnmf")
dir.create(file.path(mock_cnmf_dir, "gep_gene_tables"), recursive = TRUE, showWarnings = FALSE)
writeLines(c("{}"), file.path(mock_cnmf_dir, "cnmf_input_placeholder.h5ad"))
cnmf_plan <- pa_run_cnmf_runner(
  adata_h5ad_path = file.path(mock_cnmf_dir, "cnmf_input_placeholder.h5ad"),
  output_dir = mock_cnmf_dir,
  run_name = "demo_cnmf",
  dry_run = TRUE
)
assert_identical(cnmf_plan$status, "dry_run", "cNMF dry-run should return dry_run status")

jsonlite::write_json(
  list(recommended_k = 4, confidence = "medium", reasoning = c("mock")),
  file.path(mock_cnmf_dir, "k_selection_recommendation.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
jsonlite::write_json(
  list(success = TRUE, recommendation = list(recommended_k = 4)),
  file.path(mock_cnmf_dir, "run_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
utils::write.table(
  data.frame(k = 4, gep = c("GEP_1", "GEP_1", "GEP_2", "GEP_2"), rank = c(1,2,1,2), gene = c("MS4A1","CD79A","MKI67","TOP2A"), score = c(2,1.5,3,2.5)),
  file = file.path(mock_cnmf_dir, "gep_gene_tables", "gep_gene_scores_k4.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cnmf_import <- pa_import_cnmf_results(mock_cnmf_dir)
assert_true(is.data.frame(cnmf_import$program_tbl), "cNMF importer should return program_tbl")
assert_true(nrow(cnmf_import$program_tbl) == 2, "cNMF importer should reconstruct two GEP programs")

# PAGA runner: dry-run smoke
paga_plan <- pa_run_paga_runner(
  adata_h5ad_path = file.path(mock_cnmf_dir, "cnmf_input_placeholder.h5ad"),
  output_dir = file.path(tmp_dir, "paga"),
  groupby_key = "cell_type_final_l3",
  dry_run = TRUE
)
assert_identical(paga_plan$status, "dry_run", "PAGA dry-run should return dry_run status")

# CytoTRACE2 runner: importer mode
cyto_csv <- file.path(tmp_dir, "cytotrace2_scores.csv")
utils::write.csv(
  data.frame(cell_id = colnames(seu), cytotrace2_score = seq(0.9, 0.1, length.out = ncol(seu))),
  cyto_csv,
  row.names = FALSE,
  quote = FALSE
)
cyto_res <- pa_run_cytotrace2_runner(
  output_dir = file.path(tmp_dir, "cytotrace2"),
  score_csv = cyto_csv,
  pseudotime_table = slingshot_res$trajectory_packet$pseudotime_table
)
assert_identical(cyto_res$validation_packet$validator, "CytoTRACE2", "CytoTRACE2 runner should return validation packet")

# Enrichment runner
enrich_res <- pa_run_enrichment_runner(
  gene_vector = c("MKI67","TOP2A","CCNB1","CDK1","UBE2C","BUB1","BIRC5","MCM2"),
  output_dir = file.path(tmp_dir, "enrichment"),
  species = "human",
  sources = c("GO_BP")
)
assert_true(is.list(enrich_res$enrichment_packet), "Enrichment runner should return enrichment packet")

# interpret_agent runner: dry-run prompt generation
interpret_res <- pa_run_interpret_agent_runner(
  enrichment_packet = enrich_res$enrichment_packet,
  output_dir = file.path(tmp_dir, "interpret"),
  context_str = "Tiny smoke-test cell state comparison",
  dry_run = TRUE
)
assert_true(is.list(interpret_res$interpret_packet), "interpret_agent runner should return packet")
assert_true(file.exists(file.path(tmp_dir, "interpret", "interpret_agent_prompt.txt")), "interpret_agent dry-run should write prompt")

# REPORT runner
unit <- pa_new_analysis_unit(
  lineage = "Bcell",
  state_level = "L3",
  contrast_id = "smoke",
  condition_a = "A",
  condition_b = "B",
  trajectory_enabled = TRUE,
  output_dir = file.path(tmp_dir, "unit")
)
report_res <- pa_run_report_runner(
  unit = unit,
  output_dir = file.path(tmp_dir, "report"),
  program_registry = covar_res$registry,
  trajectory_branch_packet = pa_build_trajectory_branch_packet(
    topology_screen = list(screen_type = "PAGA"),
    primary_trajectory = slingshot_res$trajectory_packet,
    maturity_validation = cyto_res$validation_packet
  ),
  synthesis_packet = pa_build_synthesis_packet(
    enrichment_packet = enrich_res$enrichment_packet,
    interpret_agent_packet = interpret_res$interpret_packet
  )
)
assert_true(file.exists(report_res$report_packet$report_md_path), "REPORT runner should write REPORT.md")

cat("All full compute runner helper tests passed.\n")
