#!/usr/bin/env Rscript
# ==============================================================================
# Aggregation Registry Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) {
  source(PA_CORE_HELPER_PATH_20260428_V1)
}

pa_register_aggregation <- function(
  unit_id,
  aggregation_type,
  source_assay = NA_character_,
  source_layer = NA_character_,
  lineage = NA_character_,
  group_keys = character(),
  k = NA_real_,
  max_shared = NA_real_,
  min_cells = NA_real_,
  replicate_col = NA_character_,
  n_units = NA_integer_,
  qc_summary_path = NA_character_
) {
  data.frame(
    aggregation_id = paste(pa_scalar_chr(unit_id, "unit_id"), pa_scalar_chr(aggregation_type, "aggregation_type"), sep = "__"),
    unit_id = pa_scalar_chr(unit_id, "unit_id"),
    aggregation_type = pa_scalar_chr(aggregation_type, "aggregation_type"),
    source_assay = pa_null_coalesce(source_assay, NA_character_),
    source_layer = pa_null_coalesce(source_layer, NA_character_),
    lineage = pa_null_coalesce(lineage, NA_character_),
    group_keys = I(list(pa_unique_chr(group_keys))),
    k = as.numeric(k),
    max_shared = as.numeric(max_shared),
    min_cells = as.numeric(min_cells),
    replicate_col = pa_null_coalesce(replicate_col, NA_character_),
    n_units = as.integer(n_units),
    qc_summary_path = pa_null_coalesce(qc_summary_path, NA_character_),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

pa_build_metacell_qc_packet <- function(
  grouped_cells = NULL,
  total_umis = NULL,
  outlier_fraction = NULL,
  dissolve_fraction = NULL,
  purity_score = NULL,
  inner_fold_summary = NULL,
  inner_stdev_summary = NULL,
  lateral_gene_summary = NULL,
  noisy_gene_summary = NULL
) {
  list(
    grouped_cells = grouped_cells,
    total_umis = total_umis,
    outlier_fraction = outlier_fraction,
    dissolve_fraction = dissolve_fraction,
    purity_score = purity_score,
    inner_fold_summary = inner_fold_summary,
    inner_stdev_summary = inner_stdev_summary,
    lateral_gene_summary = lateral_gene_summary,
    noisy_gene_summary = noisy_gene_summary
  )
}

if (sys.nframe() == 0) {
  cat("Aggregation Registry Helper (2026-04-28 v1) loaded.\n")
}