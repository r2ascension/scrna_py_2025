#!/usr/bin/env Rscript
# ==============================================================================
# Program Architecture Helper Bundle (2026-04-28 v1)
# ==============================================================================

pa_bundle_paths_20260428_v1 <- c(
  "/home/h2048/script/R/program_architecture_core_20260428_v1.R",
  "/home/h2048/script/R/aggregation_registry_helper_20260428_v1.R",
  "/home/h2048/script/R/program_registry_helper_20260428_v1.R",
  "/home/h2048/script/R/program_source_helper_20260428_v1.R",
  "/home/h2048/script/R/trajectory_branch_helper_20260428_v1.R",
  "/home/h2048/script/R/program_support_helper_20260429_v1.R",
  "/home/h2048/script/R/synthesis_bridge_helper_20260428_v1.R",
  "/home/h2048/script/R/program_validation_helper_20260428_v1.R",
  "/home/h2048/script/R/program_crosswalk_helper_20260428_v1.R",
  "/home/h2048/script/R/evidence_packet_helper_20260428_v1.R",
  "/home/h2048/script/R/reviewer_bundle_helper_20260428_v1.R"
)

missing_paths <- pa_bundle_paths_20260428_v1[!file.exists(pa_bundle_paths_20260428_v1)]
if (length(missing_paths) > 0L) {
  stop(sprintf(
    "Program architecture helper bundle is missing required files: %s",
    paste(missing_paths, collapse = ", ")
  ), call. = FALSE)
}

invisible(lapply(pa_bundle_paths_20260428_v1, source))

if (sys.nframe() == 0) {
  cat("Program Architecture Helper Bundle (2026-04-28 v1) loaded.\n")
}