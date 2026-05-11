#!/usr/bin/env Rscript
# ==============================================================================
# Program Validation Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) {
  source(PA_CORE_HELPER_PATH_20260428_V1)
}

pa_build_program_validation_packet <- function(
  dme_table = NULL,
  trait_correlation = NULL,
  preservation_stats = NULL,
  projection_stats = NULL,
  robustness_summary = NULL,
  roe_table = NULL,
  roe_matrix = NULL,
  rogue_table = NULL,
  rogue_summary = NULL,
  tree_summary = NULL
) {
  list(
    dme_table = dme_table,
    trait_correlation = trait_correlation,
    preservation_stats = preservation_stats,
    projection_stats = projection_stats,
    robustness_summary = robustness_summary,
    roe_table = roe_table,
    roe_matrix = roe_matrix,
    rogue_table = rogue_table,
    rogue_summary = rogue_summary,
    tree_summary = tree_summary
  )
}

pa_validation_packet_summary_lines <- function(validation_packet) {
  if (is.null(validation_packet) || !is.list(validation_packet)) return("- 未提供 program validation packet。")
  c(
    sprintf("- DME 表：%s", if (is.data.frame(validation_packet$dme_table)) sprintf("%d 行", nrow(validation_packet$dme_table)) else "无"),
    sprintf("- trait correlation：%s", if (is.data.frame(validation_packet$trait_correlation)) sprintf("%d 行", nrow(validation_packet$trait_correlation)) else "无"),
    sprintf("- preservation：%s", if (is.data.frame(validation_packet$preservation_stats)) sprintf("%d 行", nrow(validation_packet$preservation_stats)) else "无"),
    sprintf("- robustness：%s", if (is.null(validation_packet$robustness_summary)) "无" else "已提供"),
    sprintf(
      "- Ro/e：%s",
      if (is.matrix(validation_packet$roe_matrix) || is.data.frame(validation_packet$roe_table)) {
        if (is.matrix(validation_packet$roe_matrix)) {
          sprintf("%d x %d", nrow(validation_packet$roe_matrix), ncol(validation_packet$roe_matrix))
        } else {
          sprintf("%d 行", nrow(validation_packet$roe_table))
        }
      } else {
        "无"
      }
    ),
    sprintf("- ROGUE：%s", if (is.data.frame(validation_packet$rogue_table)) sprintf("%d 行", nrow(validation_packet$rogue_table)) else "无"),
    sprintf("- ggtree：%s", if (is.data.frame(validation_packet$tree_summary) && nrow(validation_packet$tree_summary) > 0L) sprintf("%d 个树", nrow(validation_packet$tree_summary)) else "无")
  )
}

if (sys.nframe() == 0) {
  cat("Program Validation Helper (2026-04-28 v1) loaded.\n")
}