#!/usr/bin/env Rscript
# ==============================================================================
# Program Crosswalk Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
PA_PROGRAM_REGISTRY_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_registry_helper_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) {
  source(PA_CORE_HELPER_PATH_20260428_V1)
}
if (!exists("pa_register_programs", mode = "function")) {
  source(PA_PROGRAM_REGISTRY_HELPER_PATH_20260428_V1)
}

pa_empty_crosswalk <- function() {
  data.frame(
    contrast_id = character(),
    evidence_type = character(),
    direction = character(),
    target_program_id = character(),
    overlap_n = integer(),
    overlap_frac = numeric(),
    program_coverage = numeric(),
    enrichment_p = numeric(),
    dominant_rank = integer(),
    interpretation_tag = character(),
    overlap_genes = I(list()),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

pa_build_deg_program_crosswalk <- function(
  deg_tbl,
  program_registry,
  direction_col = "direction",
  gene_col = "gene",
  contrast_col = "contrast_id"
) {
  pa_validate_required_columns(deg_tbl, c(gene_col, direction_col), "deg_tbl")
  pa_validate_required_columns(program_registry, c("program_id", "gene_vector"), "program_registry")

  if (nrow(deg_tbl) == 0L || nrow(program_registry) == 0L) return(pa_empty_crosswalk())

  directions <- pa_unique_chr(deg_tbl[[direction_col]])
  if (length(directions) == 0L) return(pa_empty_crosswalk())

  contrast_id <- if (contrast_col %in% colnames(deg_tbl)) {
    first_contrast <- pa_safe_trim(deg_tbl[[contrast_col]])
    first_contrast <- first_contrast[nzchar(first_contrast)]
    if (length(first_contrast) == 0L) NA_character_ else first_contrast[[1]]
  } else {
    NA_character_
  }

  rows <- vector("list", length = 0L)
  for (cur_dir in directions) {
    deg_genes <- pa_normalize_gene_vector(deg_tbl[[gene_col]][pa_safe_trim(deg_tbl[[direction_col]]) == cur_dir])
    for (i in seq_len(nrow(program_registry))) {
      pg <- pa_normalize_gene_vector(program_registry$gene_vector[[i]])
      overlap <- intersect(deg_genes, pg)
      rows[[length(rows) + 1L]] <- data.frame(
        contrast_id = pa_null_coalesce(contrast_id, NA_character_),
        evidence_type = "DEG",
        direction = cur_dir,
        target_program_id = program_registry$program_id[[i]],
        overlap_n = as.integer(length(overlap)),
        overlap_frac = if (length(deg_genes) == 0L) 0 else length(overlap) / length(deg_genes),
        program_coverage = if (length(pg) == 0L) 0 else length(overlap) / length(pg),
        enrichment_p = NA_real_,
        dominant_rank = NA_integer_,
        interpretation_tag = if (length(overlap) == 0L) "none" else "candidate_program_overlap",
        overlap_genes = I(list(overlap)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  if (nrow(out) == 0L) return(pa_empty_crosswalk())

  split_idx <- split(seq_len(nrow(out)), out$direction)
  for (dir_name in names(split_idx)) {
    idx <- split_idx[[dir_name]]
    ord <- order(-out$overlap_n[idx], -out$overlap_frac[idx], out$target_program_id[idx])
    ranked_idx <- idx[ord]
    out$dominant_rank[ranked_idx] <- seq_along(ranked_idx)
    out[idx, ] <- out[ranked_idx, ]
  }
  rownames(out) <- NULL
  out
}

if (sys.nframe() == 0) {
  cat("Program Crosswalk Helper (2026-04-28 v1) loaded.\n")
}