#!/usr/bin/env Rscript

bundle_path <- "/home/h2048/script/R/program_architecture_bundle_20260428_v1.R"
source(bundle_path)

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

assert_equal <- function(x, y, msg) {
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
  output_dir = tempdir()
)
assert_true(is.data.frame(unit), "analysis unit should be a data.frame")
assert_true("unit_id" %in% colnames(unit), "analysis unit should contain unit_id")
assert_equal(unit$unit_id[[1]], "Bcell__L3__nose_vs_lung__cell", "unit_id should follow expected convention")

agg <- pa_register_aggregation(
  unit_id = unit$unit_id[[1]],
  aggregation_type = "metacell",
  source_assay = "RNA",
  source_layer = "counts",
  lineage = "Bcell",
  group_keys = c("sample", "tissue", "cell_type_L3"),
  k = 30,
  min_cells = 20,
  n_units = 5L
)
assert_true(is.data.frame(agg), "aggregation registry row should be a data.frame")
assert_equal(agg$aggregation_type[[1]], "metacell", "aggregation type should be preserved")
assert_equal(agg$group_keys[[1]], c("sample", "tissue", "cell_type_L3"), "group keys should be stored as list-column")

program_tbl <- data.frame(
  program_id = c("cnmf_1", "wgcna_blue"),
  stringsAsFactors = FALSE
)
program_tbl$gene_vector <- list(c("MS4A1", "CD79A", "HLA-DRA"), c("CD74", "HLA-DPA1"))
program_tbl$optional_weight <- list(c(0.8, 0.7, 0.6), c(0.5, 0.4))

prog <- pa_register_programs(
  unit_id = unit$unit_id[[1]],
  source_type = "cNMF",
  source_subtype = "GEP",
  program_tbl = program_tbl,
  score_level = "cell"
)
assert_true(is.data.frame(prog), "program registry should be a data.frame")
assert_equal(nrow(prog), 2L, "program registry should contain two programs")
assert_equal(prog$program_id[[1]], "cnmf_1", "program id should be preserved")

validation_packet <- pa_build_program_validation_packet(
  dme_table = data.frame(module = "blue", logFC = 0.8, stringsAsFactors = FALSE),
  trait_correlation = data.frame(trait = "nose", rho = 0.5, stringsAsFactors = FALSE),
  preservation_stats = data.frame(module = "blue", z = 6.1, stringsAsFactors = FALSE),
  robustness_summary = list(stable = TRUE)
)
assert_true(is.list(validation_packet), "validation packet should be a list")
assert_true("dme_table" %in% names(validation_packet), "validation packet should contain dme_table")

deg_tbl <- data.frame(
  gene = c("MS4A1", "CD79A", "XBP1", "JCHAIN"),
  direction = c("up", "up", "down", "down"),
  contrast_id = c("nose_vs_lung", "nose_vs_lung", "nose_vs_lung", "nose_vs_lung"),
  stringsAsFactors = FALSE
)

crosswalk <- pa_build_deg_program_crosswalk(
  deg_tbl = deg_tbl,
  program_registry = prog,
  direction_col = "direction",
  gene_col = "gene"
)
assert_true(is.data.frame(crosswalk), "crosswalk should be a data.frame")
assert_true(nrow(crosswalk) >= 2L, "crosswalk should have rows")
assert_true(all(c("target_program_id", "overlap_n", "dominant_rank") %in% colnames(crosswalk)), "crosswalk should expose ranking columns")

network_summary <- list(edges_changed = 12L, gained_hubs = c("MS4A1"))
evidence_packet <- pa_build_evidence_packet(
  unit_id = unit$unit_id[[1]],
  aggregation_context = agg,
  de_table = deg_tbl,
  enrichment_list = list(Hallmark = data.frame(term = "INTERFERON_GAMMA_RESPONSE", stringsAsFactors = FALSE)),
  program_crosswalk = crosswalk,
  validation_summary = validation_packet,
  network_summary = network_summary,
  context_str = "demo context"
)
assert_true(is.list(evidence_packet), "evidence packet should be a list")
assert_true(identical(evidence_packet$unit_id, unit$unit_id[[1]]), "evidence packet should retain unit_id")

fig_reg <- pa_register_figure(
  figure_id = "Fig2",
  panel_id = "B",
  source_unit = unit$unit_id[[1]],
  input_tables = c("crosswalk.tsv", "validation.tsv"),
  claim_sentence = "程序层和离散层在该比较中一致指向 B 细胞活化。",
  support_packet = evidence_packet,
  output_path = file.path(tempdir(), "Fig2B.png")
)
assert_true(is.data.frame(fig_reg), "figure registry row should be a data.frame")

bundle_dir <- file.path(tempdir(), "reviewer_bundle_test")
out_bundle <- pa_write_reviewer_bundle(
  output_dir = bundle_dir,
  unit = unit,
  evidence_packet = evidence_packet,
  figure_registry = fig_reg,
  program_registry = prog,
  validation_packet = validation_packet
)
assert_true(file.exists(out_bundle$methods_snippet), "methods snippet should be written")
assert_true(file.exists(out_bundle$result_snippet), "result snippet should be written")
assert_true(file.exists(out_bundle$figure_registry_tsv), "figure registry TSV should be written")

run_outdir <- file.path(tempdir(), "program_architecture_run")
run_result <- pa_run_program_architecture_unit(
  unit = unit,
  discrete_outputs = list(deg_table = deg_tbl, enrichment_list = list()),
  program_inputs = list(
    program_registry = prog,
    dme_table = validation_packet$dme_table,
    trait_correlation = validation_packet$trait_correlation,
    preservation_stats = validation_packet$preservation_stats,
    robustness_summary = validation_packet$robustness_summary,
    score_summary = list(cell = "ok")
  ),
  rewiring_inputs = list(network_summary = network_summary, conflict_summary = list()),
  outdir = run_outdir
)
assert_true(file.exists(file.path(run_outdir, "analysis_unit.rds")), "analysis unit RDS should be written")
assert_true(is.list(run_result$evidence_packet), "runner should return evidence packet")

cat("All program architecture helper tests passed.\n")
