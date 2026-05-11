#!/usr/bin/env Rscript
# ==============================================================================
# Smoke tests for CD8 TCR STARTRAC panel runner helpers (2026-05-08 v1)
# ==============================================================================

options(warn = 1)

RUNNER_PATH <- "/home/h2048/script/R/cd8_tcr_startrac_panel_runner_20260508_v1.R"
source(RUNNER_PATH)

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

assert_equal <- function(x, y, msg) {
  if (!identical(x, y)) {
    stop(sprintf("%s\nExpected: %s\nActual: %s", msg, paste(y, collapse = ","), paste(x, collapse = ",")), call. = FALSE)
  }
}

assert_near <- function(x, y, tol = 1e-8, msg = "values differ") {
  if (is.na(x) || is.na(y) || abs(x - y) > tol) {
    stop(sprintf("%s\nExpected: %.10f\nActual: %.10f", msg, y, x), call. = FALSE)
  }
}

meta <- data.frame(
  cell_id = paste0("cell", 1:10),
  sample = c("S1", "S1", "S1", "S2", "S2", "S2", "S2", "S2", "S2", "S2"),
  group = c("NP", "NP", "NP", "CBL", "CBL", "CBL", "CBL", "CBL", "CBL", "CBL"),
  cell_type_L2 = rep("CD8 T cells", 10),
  cell_type_L3 = c("CD8_Tem", "CD8_Tem", "CD8_Trm", "CD8_Tem", "CD8_Tem", "CD8_Tem", "CD8_Trm", "CD8_Trm", "CD8_Trm", "CD8_Trm"),
  clone_id = c("clA", "clB", "clB", "clC", "clC", "clC", "clC", "clC", "clD", "clD"),
  gzmk_expr = c(0, 3, 4, 1, 2, 3, 4, 5, 0, 0),
  stringsAsFactors = FALSE
)

status <- cd8tcr_compute_clone_status(meta, clone_col = "clone_id")
assert_equal(
  status$clone_status[match(c("clA", "clB", "clC"), status$clone_id)],
  c("singleton", "small", "medium"),
  "clone status bins should follow observed clone sizes"
)

div <- cd8tcr_compute_diversity_metrics(c(clA = 1, clB = 1, clC = 2))
assert_near(div$shannon, -(1/4 * log(1/4) + 1/4 * log(1/4) + 2/4 * log(2/4)), msg = "Shannon diversity mismatch")
assert_near(div$inverse_simpson, 1 / ((1/4)^2 + (1/4)^2 + (2/4)^2), msg = "Inverse Simpson mismatch")
assert_near(div$chao1, 5, msg = "Chao1 should use singleton/doubleton correction")
assert_true(div$ace >= 2, "ACE should be at least observed richness for this simple vector")

sample_div <- cd8tcr_compute_diversity_by_sample(meta, sample_col = "sample", group_col = "group", clone_col = "clone_id")
assert_equal(sort(sample_div$sample), c("S1", "S2"), "diversity should be sample-level")
assert_true(all(c("shannon", "inverse_simpson", "inverse_pielou", "chao1", "ace") %in% colnames(sample_div)), "diversity output missing expected columns")

pseudo_meta <- data.frame(
  cell_id = paste0("cell", 1:4),
  clone_id = paste0("cell", 1:4),
  stringsAsFactors = FALSE
)
pseudo_check <- cd8tcr_validate_tcr_metadata(pseudo_meta, cell_col = "cell_id", clone_col = "clone_id")
assert_true(!pseudo_check$has_real_clones, "cell-name clone IDs must not pass real clone validation")
assert_true(any(grepl("cell_id", pseudo_check$messages)), "pseudo-clone validation should explain the cell_id problem")

cd8_prepared <- cd8tcr_prepare_cd8_metadata(
  meta,
  l2_col = "cell_type_L2",
  l3_col = "cell_type_L3",
  gzmk_expr_col = "gzmk_expr",
  gzmk_threshold = 2
)
assert_equal(nrow(cd8_prepared), 10L, "all input rows are CD8 cells in this fixture")
assert_equal(sum(cd8_prepared$gzmk_positive), 6L, "GZMK positivity should use >= threshold")

panel_status <- cd8tcr_build_panel_status(has_tcr = FALSE, has_expression = TRUE)
assert_equal(panel_status$status[panel_status$panel == "f"], "available", "Monocle3 expression panel should be available without TCR")
assert_equal(panel_status$status[panel_status$panel == "h"], "blocked_missing_tcr", "TCRmatch panel should be blocked without TCR")

tmp_root <- file.path(tempdir(), sprintf("cd8tcr_runner_test_%s", Sys.getpid()))
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
metadata_path <- file.path(tmp_root, "metadata.tsv")
utils::write.table(meta, metadata_path, sep = "\t", row.names = FALSE, quote = FALSE)

cfg <- cd8tcr_default_config()
cfg$metadata_table <- metadata_path
cfg$input_rds <- ""
cfg$input_h5ad <- ""
cfg$tcr_metadata_path <- ""
cfg$sample_metadata_path <- ""
cfg$output_dir <- file.path(tmp_root, "out")
cfg$group_col <- "group"
cfg$sample_col <- "sample"
cfg$l2_col <- "cell_type_L2"
cfg$l3_col <- "cell_type_L3"
cfg$gzmk_expr_col <- "gzmk_expr"
cfg$gzmk_threshold <- 2

runner_result <- cd8tcr_run_pipeline(cfg)
assert_true(runner_result$has_tcr, "pipeline should detect repeated real clone IDs from metadata")
assert_true(file.exists(file.path(cfg$output_dir, "CD8_TCR_STARTRAC_REPORT.md")), "pipeline should write report")
assert_true(file.exists(file.path(cfg$output_dir, "tables", "panel_a_clone_status_by_group.tsv")), "pipeline should write panel a source table")
assert_true(file.exists(file.path(cfg$output_dir, "tables", "panel_b_cd8_diversity_by_sample.tsv")), "pipeline should write panel b source table")
assert_true(file.exists(file.path(cfg$output_dir, "tables", "panel_d_gzmk_cd8_top100_clonotypes.tsv")), "pipeline should write panel d source table")

missing_group_cfg <- cfg
missing_group_cfg$output_dir <- file.path(tmp_root, "out_missing_disease_group")
missing_group_cfg$panel_ab_group_col <- "disease_group"
missing_group_result <- cd8tcr_run_pipeline(missing_group_cfg)
assert_true(missing_group_result$has_tcr, "fixture still has real TCR clones when disease group is missing")
assert_equal(
  missing_group_result$panel_status$status[missing_group_result$panel_status$panel == "a"],
  "blocked_missing_disease_group",
  "panel a should be blocked when explicit disease group metadata is missing"
)
assert_true(
  !file.exists(file.path(missing_group_cfg$output_dir, "tables", "panel_a_clone_status_by_group.tsv")),
  "panel a table must not be written without explicit disease group metadata"
)

sample_design <- data.frame(
  sample = c(paste0("CBL", 1:5), paste0("NPBL", 1:8), paste0("CIT", 1:4), paste0("NP", 1:8)),
  disease_group = c(rep("CBL", 5), rep("NP-BL", 8), rep("CIT", 4), rep("NP", 8)),
  stringsAsFactors = FALSE
)

panel_ab_rows <- list()
panel_ab_tcr <- list()
for (sample_index in seq_len(nrow(sample_design))) {
  sample_id <- sample_design$sample[[sample_index]]
  clone_sizes <- c(2 + sample_index %% 4, 1 + sample_index %% 3, 1)
  cell_counter <- 1L
  for (clone_index in seq_along(clone_sizes)) {
    clone_id <- sprintf("%s_clone_%d", sample_id, clone_index)
    for (cell_in_clone in seq_len(clone_sizes[[clone_index]])) {
      cell_id <- sprintf("%s_cell_%03d", sample_id, cell_counter)
      panel_ab_rows[[length(panel_ab_rows) + 1L]] <- data.frame(
        cell_id = cell_id,
        sample = sample_id,
        group = "respiratory_tissue",
        cell_type_L2 = "CD8 T cells",
        cell_type_L3 = if (cell_counter %% 2 == 0L) "GZMK+ CD8 T cells" else "CD8_Tem",
        gzmk_expr = if (cell_counter %% 2 == 0L) 3 else 0,
        stringsAsFactors = FALSE
      )
      panel_ab_tcr[[length(panel_ab_tcr) + 1L]] <- data.frame(
        cell_id = cell_id,
        clonotype_id = clone_id,
        stringsAsFactors = FALSE
      )
      cell_counter <- cell_counter + 1L
    }
  }
}

panel_ab_meta <- do.call(rbind, panel_ab_rows)
panel_ab_tcr <- do.call(rbind, panel_ab_tcr)
panel_ab_metadata_path <- file.path(tmp_root, "panel_ab_metadata.tsv")
panel_ab_tcr_path <- file.path(tmp_root, "panel_ab_tcr.tsv")
panel_ab_sample_path <- file.path(tmp_root, "panel_ab_sample_metadata.tsv")
utils::write.table(panel_ab_meta, panel_ab_metadata_path, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(panel_ab_tcr, panel_ab_tcr_path, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(sample_design, panel_ab_sample_path, sep = "\t", row.names = FALSE, quote = FALSE)

panel_ab_cfg <- cd8tcr_default_config()
panel_ab_cfg$metadata_table <- panel_ab_metadata_path
panel_ab_cfg$input_rds <- ""
panel_ab_cfg$input_h5ad <- ""
panel_ab_cfg$tcr_metadata_path <- panel_ab_tcr_path
panel_ab_cfg$sample_metadata_path <- panel_ab_sample_path
panel_ab_cfg$output_dir <- file.path(tmp_root, "out_panel_ab")
panel_ab_cfg$group_col <- "group"
panel_ab_cfg$panel_ab_group_col <- "disease_group"
panel_ab_cfg$sample_col <- "sample"
panel_ab_cfg$l2_col <- "cell_type_L2"
panel_ab_cfg$l3_col <- "cell_type_L3"
panel_ab_cfg$gzmk_expr_col <- "gzmk_expr"
panel_ab_cfg$gzmk_threshold <- 2
panel_ab_cfg$strict_four_group <- TRUE

panel_ab_result <- cd8tcr_run_pipeline(panel_ab_cfg)
assert_equal(
  panel_ab_result$panel_status$status[panel_ab_result$panel_status$panel == "a"],
  "available",
  "panel a should be available with real clonotypes and explicit four-group mapping"
)
assert_equal(
  panel_ab_result$panel_status$status[panel_ab_result$panel_status$panel == "b"],
  "available",
  "panel b should be available with real clonotypes and explicit four-group mapping"
)
assert_true(file.exists(file.path(panel_ab_cfg$output_dir, "tables", "four_group_design_validation.tsv")), "four-group validation table should be written")
four_group_validation <- cd8tcr_read_table(file.path(panel_ab_cfg$output_dir, "tables", "four_group_design_validation.tsv"))
assert_true(all(four_group_validation$ok), "fixture should match 5/8/4/8 sample design")
panel_a_table <- cd8tcr_read_table(file.path(panel_ab_cfg$output_dir, "tables", "panel_a_clone_status_by_group.tsv"))
assert_true("disease_group" %in% colnames(panel_a_table), "panel a should be grouped by disease_group")
diversity_stats <- cd8tcr_read_table(file.path(panel_ab_cfg$output_dir, "tables", "panel_b_cd8_diversity_stats.tsv"))
assert_true(any(diversity_stats$method == "dunn.test.BH"), "panel b should use Dunn multiple comparisons with BH adjustment")
assert_true(!any(grepl("wilcox", diversity_stats$method, ignore.case = TRUE)), "panel b stats must not silently fall back to Wilcoxon")
assert_true(file.exists(file.path(panel_ab_cfg$output_dir, "tables", "panel_b_cd8_diversity_plot_long.tsv")), "panel b should write a long-format plot source table")

bad_sample_design <- sample_design[sample_design$sample != "NP8", , drop = FALSE]
bad_panel_ab_meta <- panel_ab_meta[panel_ab_meta$sample %in% bad_sample_design$sample, , drop = FALSE]
bad_panel_ab_tcr <- panel_ab_tcr[panel_ab_tcr$cell_id %in% bad_panel_ab_meta$cell_id, , drop = FALSE]
bad_metadata_path <- file.path(tmp_root, "bad_panel_ab_metadata.tsv")
bad_tcr_path <- file.path(tmp_root, "bad_panel_ab_tcr.tsv")
bad_sample_path <- file.path(tmp_root, "bad_panel_ab_sample_metadata.tsv")
utils::write.table(bad_panel_ab_meta, bad_metadata_path, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(bad_panel_ab_tcr, bad_tcr_path, sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(bad_sample_design, bad_sample_path, sep = "\t", row.names = FALSE, quote = FALSE)

bad_panel_ab_cfg <- panel_ab_cfg
bad_panel_ab_cfg$metadata_table <- bad_metadata_path
bad_panel_ab_cfg$tcr_metadata_path <- bad_tcr_path
bad_panel_ab_cfg$sample_metadata_path <- bad_sample_path
bad_panel_ab_cfg$output_dir <- file.path(tmp_root, "out_bad_panel_ab")
bad_panel_ab_result <- cd8tcr_run_pipeline(bad_panel_ab_cfg)
assert_equal(
  bad_panel_ab_result$panel_status$status[bad_panel_ab_result$panel_status$panel == "a"],
  "blocked_invalid_four_group",
  "strict sample-count mismatch should block panel a"
)
assert_equal(
  bad_panel_ab_result$panel_status$status[bad_panel_ab_result$panel_status$panel == "b"],
  "blocked_invalid_four_group",
  "strict sample-count mismatch should block panel b"
)
bad_four_group_validation <- cd8tcr_read_table(file.path(bad_panel_ab_cfg$output_dir, "tables", "four_group_design_validation.tsv"))
assert_true(!all(bad_four_group_validation$ok), "bad fixture should record failed 5/8/4/8 sample validation")
assert_true(!file.exists(file.path(bad_panel_ab_cfg$output_dir, "tables", "panel_a_clone_status_by_group.tsv")), "strict invalid design should not write panel a table")
assert_true(!file.exists(file.path(bad_panel_ab_cfg$output_dir, "tables", "panel_b_cd8_diversity_by_sample.tsv")), "strict invalid design should not write panel b table")

cat("[TEST] cd8_tcr_startrac_panel_runner helper smoke tests passed\n")