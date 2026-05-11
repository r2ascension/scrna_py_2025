#!/usr/bin/env Rscript

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

helper_path <- "/home/h2048/script/R/stromal_smc_trajectory_interpret_helper_20260429.R"
sys.source(helper_path, envir = globalenv())

trend_summary <- data.frame(
  method = c(rep("ToyMethod", 6), rep("ToyMethod", 3)),
  trajectory_id = c(rep("Lineage1", 6), rep("Lineage2", 3)),
  gene = c(rep("ACTA2", 6), rep("TAGLN", 3)),
  bin_id = c(1:6, 1:3),
  bin_mid = c(seq(1/12, 11/12, length.out = 6), c(0.15, 0.45, 0.75)),
  mean_scaled_pseudotime = c(seq(0.08, 0.92, length.out = 6), c(0.16, 0.48, 0.80)),
  n_cells = c(rep(12L, 6), rep(8L, 3)),
  mean_pseudotime = c(seq(0.1, 0.9, length.out = 6), c(0.2, 0.5, 0.8)),
  mean_curve_weight = c(rep(0.9, 6), rep(0.75, 3)),
  mean_expression = c(2, 3, 5, 6, 8, 9, 1.5, 1.7, 1.6),
  median_expression = c(2, 3, 5, 6, 8, 9, 1.5, 1.7, 1.6),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

smoothed_summary <- smcti_add_smoothed_expression_to_summary(
  trend_summary,
  smoother_method = "smooth_spline",
  smoother_min_points = 5L
)

assert_true("smoothed_expression" %in% colnames(smoothed_summary), "Expected smoothed_expression column in trend summary")
assert_true("smoothing_applied" %in% colnames(smoothed_summary), "Expected smoothing_applied column in trend summary")
assert_true("smoother_method" %in% colnames(smoothed_summary), "Expected smoother_method column in trend summary")

acta2_rows <- smoothed_summary[smoothed_summary$gene == "ACTA2", , drop = FALSE]
tagln_rows <- smoothed_summary[smoothed_summary$gene == "TAGLN", , drop = FALSE]
assert_true(all(is.finite(acta2_rows$smoothed_expression)), "Expected finite smoothed values for adequately sized groups")
assert_true(all(acta2_rows$smoothing_applied), "Expected smoothing to be applied to ACTA2 group")
assert_true(!any(tagln_rows$smoothing_applied), "Expected short groups to fall back without smoothing")

plot_df <- smcti_build_gene_trend_plot_data(
  smoothed_summary,
  smoother_method = "smooth_spline",
  smoother_grid_n = 41L,
  smoother_min_points = 5L
)

assert_true(all(c("plot_layer", "plot_expression") %in% colnames(plot_df)), "Expected layered plot columns in plot data")
assert_true(any(plot_df$plot_layer == "observed_bin"), "Expected observed bin rows in plot data")
assert_true(any(plot_df$plot_layer == "smoothed_curve"), "Expected smoothed curve rows in plot data")
assert_true(all(c("mean_expression", "median_expression") %in% colnames(plot_df)), "Expected original summary statistics to survive in plot data")

curve_rows <- plot_df[plot_df$plot_layer == "smoothed_curve" & plot_df$gene == "ACTA2", , drop = FALSE]
assert_true(nrow(curve_rows) >= 20L, "Expected dense smoothed curve rows for sufficiently sized groups")
assert_true(all(is.finite(curve_rows$plot_expression)), "Expected finite smoothed plot expressions")

cat("smc_anno gene trend smoothing smoke test passed.\n")
