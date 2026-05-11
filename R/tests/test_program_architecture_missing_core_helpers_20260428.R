#!/usr/bin/env Rscript

bundle_path <- "/home/h2048/script/R/program_architecture_bundle_20260428_v1.R"
source(bundle_path)

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

assert_identical <- function(x, y, msg) {
  if (!identical(x, y)) {
    stop(sprintf("%s\nExpected: %s\nActual: %s", msg, paste(capture.output(str(y)), collapse = " "), paste(capture.output(str(x)), collapse = " ")), call. = FALSE)
  }
}

unit <- pa_new_analysis_unit(
  lineage = "Bcell",
  state_level = "L3",
  contrast_id = "nose_vs_lung",
  condition_a = "nose",
  condition_b = "lung",
  trajectory_enabled = TRUE,
  output_dir = tempdir()
)

cnmf_tbl <- data.frame(
  program_id = c("cnmf_1", "cnmf_2"),
  stringsAsFactors = FALSE
)
cnmf_tbl$gene_vector <- list(c("MS4A1", "CD79A"), c("JCHAIN", "XBP1"))
cnmf_tbl$usage_column <- c("GEP1", "GEP2")
cnmf_registry <- pa_register_cnmf_programs(unit$unit_id[[1]], cnmf_tbl, score_object_path = "cnmf_usages.rds")
assert_true(is.data.frame(cnmf_registry), "cNMF registry should be data.frame")
assert_identical(cnmf_registry$source_type[[1]], "cNMF", "cNMF source type should be recorded")

hd_tbl <- data.frame(
  module_id = c("blue", "turquoise"),
  stringsAsFactors = FALSE
)
hd_tbl$gene_vector <- list(c("CD74", "HLA-DRA"), c("BANK1", "CD79B"))
hd_registry <- pa_register_hdwgcna_programs(unit$unit_id[[1]], hd_tbl)
assert_true(is.data.frame(hd_registry), "hdWGCNA registry should be data.frame")
assert_identical(hd_registry$source_type[[1]], "hdWGCNA", "hdWGCNA source type should be recorded")

covar_tbl <- data.frame(
  program_id = c("covarnet_naive_b"),
  stringsAsFactors = FALSE
)
covar_tbl$gene_vector <- list(c("MS4A1", "HVCN1", "CD79A"))
covar_tbl$hub_genes <- list(c("MS4A1", "CD79A"))
covar_registry <- pa_register_covarnet_programs(unit$unit_id[[1]], covar_tbl)
assert_true(is.data.frame(covar_registry), "CoVarNet registry should be data.frame")
assert_identical(covar_registry$source_type[[1]], "CoVarNet", "CoVarNet source type should be recorded")

paga_csv <- file.path(tempdir(), "paga_connectivities.csv")
write.csv(
  data.frame(cell_a = c(0, 0.4), cell_b = c(0.4, 0), row.names = c("Naive_B", "Memory_B")),
  paga_csv,
  quote = FALSE
)
paga_import <- pa_import_paga_connectivities_csv(paga_csv)
assert_true(is.matrix(paga_import$connectivity_matrix), "PAGA import should return connectivity matrix")

paga_packet <- pa_build_paga_topology_packet(
  connectivity_matrix = paga_import$connectivity_matrix,
  group_labels = paga_import$group_labels,
  source_path = paga_csv,
  group_key = "cell_type_L3"
)
assert_identical(paga_packet$screen_type, "PAGA", "PAGA packet should mark screen type")

slingshot_packet <- pa_build_slingshot_trajectory_packet(
  pseudotime_table = data.frame(cell_id = c("c1", "c2"), lineage = c("L1", "L1"), pseudotime = c(0.1, 0.8), stringsAsFactors = FALSE),
  lineage_summary = data.frame(lineage = "L1", n_cells = 2L, root_state = "Naive_B", terminal_state = "Plasma", stringsAsFactors = FALSE),
  branch_summary = data.frame(branch_id = "B1", from = "Naive_B", to = "Plasma", stringsAsFactors = FALSE)
)
assert_identical(slingshot_packet$engine, "Slingshot", "Slingshot packet should mark engine")

cyto_packet <- pa_build_cytotrace2_validation_packet(
  score_table = data.frame(cell_id = c("c1", "c2"), cytotrace2_score = c(0.9, 0.2), stringsAsFactors = FALSE),
  maturity_summary = data.frame(lineage = "L1", rho_with_pseudotime = -0.8, stringsAsFactors = FALSE),
  direction_consistency = "consistent"
)
assert_identical(cyto_packet$validator, "CytoTRACE2", "CytoTRACE2 packet should mark validator")

traj_packet <- pa_build_trajectory_branch_packet(
  topology_screen = paga_packet,
  primary_trajectory = slingshot_packet,
  maturity_validation = cyto_packet,
  preferred_engine = "Slingshot"
)
assert_identical(traj_packet$preferred_engine, "Slingshot", "trajectory branch packet should keep preferred engine")

enrich_packet <- pa_build_enrichment_packet(
  ora_table = data.frame(term = "B_CELL_RECEPTOR_SIGNALING", padj = 0.01, stringsAsFactors = FALSE),
  source_dbs = c("Hallmark", "KEGG")
)
assert_true(is.list(enrich_packet), "enrichment packet should be list")

interpret_packet <- pa_build_interpret_agent_packet(
  structured_table = data.frame(celltype_label = "Naive_B", narrative = "测试", stringsAsFactors = FALSE),
  model = "deepseek-v4-flash"
)
assert_true(is.list(interpret_packet), "interpret packet should be list")

report_file <- file.path(tempdir(), "REPORT.md")
writeLines(c("# Demo report", "ok"), report_file)
report_packet <- pa_build_report_packet(
  report_md_path = report_file,
  sections = c("overview", "trajectory")
)
assert_true(is.list(report_packet), "report packet should be list")

synthesis_packet <- pa_build_synthesis_packet(
  enrichment_packet = enrich_packet,
  interpret_agent_packet = interpret_packet,
  report_packet = report_packet
)
assert_true(is.list(synthesis_packet), "synthesis packet should be list")

evidence_packet <- pa_build_evidence_packet(
  unit_id = unit$unit_id[[1]],
  aggregation_context = data.frame(dummy = 1),
  trajectory_branch_packet = traj_packet,
  synthesis_packet = synthesis_packet
)
assert_true(!is.null(evidence_packet$trajectory_branch_packet), "evidence packet should carry trajectory packet")
assert_true(!is.null(evidence_packet$synthesis_packet), "evidence packet should carry synthesis packet")

out_dir <- file.path(tempdir(), "pa_missing_helpers_run")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
run_result <- pa_run_program_architecture_unit(
  unit = unit,
  discrete_outputs = list(deg_table = data.frame(gene = "MS4A1", direction = "up", contrast_id = unit$contrast_id[[1]], stringsAsFactors = FALSE), enrichment_list = list()),
  program_inputs = list(program_registry = pa_bind_program_registries(cnmf_registry, hd_registry, covar_registry), trajectory_branch_packet = traj_packet, synthesis_packet = synthesis_packet),
  rewiring_inputs = list(network_summary = list(), conflict_summary = list()),
  outdir = out_dir
)
assert_true(file.exists(file.path(out_dir, "trajectory_branch_packet.rds")), "trajectory packet should be written")
assert_true(file.exists(file.path(out_dir, "synthesis_packet.rds")), "synthesis packet should be written")
assert_true(is.list(run_result$evidence_packet), "runner should return extended evidence packet")

cat("All missing core helper tests passed.\n")
