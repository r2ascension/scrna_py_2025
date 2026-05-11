#!/usr/bin/env Rscript
# ==============================================================================
# Program Registry Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) {
  source(PA_CORE_HELPER_PATH_20260428_V1)
}

pa_register_programs <- function(
  unit_id,
  source_type,
  source_subtype,
  program_tbl,
  score_level = NA_character_,
  score_object_path = NA_character_,
  validation_status = "raw"
) {
  pa_validate_required_columns(program_tbl, c("program_id", "gene_vector"), "program_tbl")
  unit_id <- pa_scalar_chr(unit_id, "unit_id")
  source_type <- pa_scalar_chr(source_type, "source_type")
  source_subtype <- pa_scalar_chr(source_subtype, "source_subtype")

  weight_col <- if ("optional_weight" %in% colnames(program_tbl)) program_tbl$optional_weight else vector("list", nrow(program_tbl))

  data.frame(
    program_id = pa_safe_trim(program_tbl$program_id),
    unit_id = rep(unit_id, nrow(program_tbl)),
    source_type = rep(source_type, nrow(program_tbl)),
    source_subtype = rep(source_subtype, nrow(program_tbl)),
    gene_vector = I(lapply(program_tbl$gene_vector, pa_normalize_gene_vector)),
    optional_weight = I(as.list(weight_col)),
    score_level = rep(pa_null_coalesce(score_level, NA_character_), nrow(program_tbl)),
    score_object_path = rep(pa_null_coalesce(score_object_path, NA_character_), nrow(program_tbl)),
    provenance_unit = rep(unit_id, nrow(program_tbl)),
    validation_status = rep(pa_null_coalesce(validation_status, "raw"), nrow(program_tbl)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

if (sys.nframe() == 0) {
  cat("Program Registry Helper (2026-04-28 v1) loaded.\n")
}