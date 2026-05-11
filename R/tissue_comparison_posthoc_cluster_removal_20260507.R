#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

source("/home/h2048/script/R/tissue_comparison_llm_outlier_visual_batch_20260507.R")

OUTPUT_STEM <- "posthoc_cluster_removal_20260507"

REMOVAL_TARGETS <- list(
  list(
    lineage_tag = "BCELL",
    output_dir = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415",
    cluster_kind = "choir",
    remove_clusters = c(14),
    removal_basis = "CHOIR outlier c14"
  ),
  list(
    lineage_tag = "EPITHELIAL",
    output_dir = "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun",
    cluster_kind = "leiden",
    remove_clusters = c(14, 17),
    removal_basis = "OFA outlier c14/c17"
  ),
  list(
    lineage_tag = "STROMAL_ENDOTHELIAL",
    output_dir = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414",
    cluster_kind = "choir",
    remove_clusters = c(6, 52),
    removal_basis = "OFA outlier c6/c52"
  )
)

count_meta_col <- function(obj, colname) {
  if (!(colname %in% colnames(obj@meta.data))) return(data.table())
  vals <- sanitize_text(obj@meta.data[[colname]], default = "<missing>")
  dt <- data.table(level = vals)
  out <- dt[, .N, by = level][order(-N, level)]
  out[, column := colname]
  setcolorder(out, c("column", "level", "N"))
  out
}

count_cluster_vector <- function(cluster_vec) {
  dt <- data.table(cluster_id = as.character(cluster_vec))
  dt[!is.na(cluster_id) & nzchar(cluster_id), .N, by = cluster_id][order(as.integer(cluster_id))]
}

process_removal_target <- function(target) {
  object_path <- find_final_object_path(target$output_dir)
  base_stub <- tools::file_path_sans_ext(basename(object_path))
  removal_dir <- file.path(target$output_dir, OUTPUT_STEM)
  dir.create(removal_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n=== %s | %s ===\n", target$lineage_tag, target$output_dir))
  obj <- readRDS(object_path)
  obj <- obj[, sort(colnames(obj))]

  assign_df <- read_cluster_assignments(target$output_dir, target$cluster_kind)
  best_match <- match_cluster_vector(colnames(obj), assign_df)
  if (is.na(best_match$rate) || best_match$rate < 0.80) {
    stop(sprintf(
      "Cluster assignment match rate too low for %s: %.3f",
      target$lineage_tag,
      best_match$rate
    ))
  }

  cluster_vec <- as.integer(best_match$cluster_vec)
  remove_clusters_chr <- as.character(target$remove_clusters)
  remove_mask <- as.character(cluster_vec) %in% remove_clusters_chr
  removed_cells <- colnames(obj)[remove_mask]
  kept_cells <- colnames(obj)[!remove_mask]

  if (length(removed_cells) == 0L) {
    stop(sprintf("No cells matched requested removal clusters for %s", target$lineage_tag))
  }

  filtered_obj <- subset(obj, cells = kept_cells)

  removed_dt <- data.table(
    cell = removed_cells,
    cluster_id = as.character(cluster_vec[remove_mask]),
    cluster_kind = target$cluster_kind,
    lineage_tag = target$lineage_tag,
    removal_basis = target$removal_basis
  )[order(as.integer(cluster_id), cell)]

  cluster_before <- count_cluster_vector(cluster_vec)
  cluster_after <- count_cluster_vector(cluster_vec[!remove_mask])
  cluster_summary <- merge(
    cluster_before,
    cluster_after,
    by = "cluster_id",
    all = TRUE,
    suffixes = c("_before", "_after")
  )
  cluster_summary[is.na(N_before), N_before := 0L]
  cluster_summary[is.na(N_after), N_after := 0L]
  cluster_summary[, removed_n := N_before - N_after]
  cluster_summary[, cluster_id_int := suppressWarnings(as.integer(cluster_id))]
  setorder(cluster_summary, cluster_id_int)
  cluster_summary[, cluster_id_int := NULL]

  meta_counts <- rbindlist(list(
    count_meta_col(obj, "tissue")[, state := "before"],
    count_meta_col(filtered_obj, "tissue")[, state := "after"],
    count_meta_col(obj, "cell_type_L2")[, state := "before"],
    count_meta_col(filtered_obj, "cell_type_L2")[, state := "after"],
    count_meta_col(obj, "cell_type_L3")[, state := "before"],
    count_meta_col(filtered_obj, "cell_type_L3")[, state := "after"]
  ), fill = TRUE)

  removed_cells_path <- file.path(removal_dir, sprintf("%s_removed_cells.tsv", base_stub))
  cluster_summary_path <- file.path(removal_dir, sprintf("%s_cluster_counts_before_after.tsv", base_stub))
  meta_counts_path <- file.path(removal_dir, sprintf("%s_meta_counts_before_after.tsv", base_stub))
  filtered_rds_path <- file.path(
    removal_dir,
    sprintf("%s_rm_%s_%s.rds", base_stub, target$cluster_kind, paste(remove_clusters_chr, collapse = "_"))
  )
  summary_path <- file.path(removal_dir, "removal_summary.tsv")

  fwrite(removed_dt, removed_cells_path, sep = "\t")
  fwrite(cluster_summary, cluster_summary_path, sep = "\t")
  fwrite(meta_counts, meta_counts_path, sep = "\t")
  saveRDS(filtered_obj, filtered_rds_path)

  summary_dt <- data.table(
    lineage_tag = target$lineage_tag,
    output_dir = target$output_dir,
    cluster_kind = target$cluster_kind,
    remove_clusters = paste(remove_clusters_chr, collapse = ","),
    removal_basis = target$removal_basis,
    original_cells = ncol(obj),
    removed_cells = length(removed_cells),
    retained_cells = ncol(filtered_obj),
    match_mode = best_match$mode,
    match_rate = sprintf("%.4f", best_match$rate),
    filtered_rds = filtered_rds_path,
    removed_cells_tsv = removed_cells_path,
    cluster_summary_tsv = cluster_summary_path,
    meta_counts_tsv = meta_counts_path
  )
  fwrite(summary_dt, summary_path, sep = "\t")

  review_dir <- file.path(target$output_dir, "figures", OUTPUT_STEM)
  review_readme <- file.path(target$output_dir, "figures", "llm_outlier_review_20260507", "README.md")
  readme_lines <- c(
    sprintf("# %s posthoc cluster removal", target$lineage_tag),
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    sprintf("- Source output: `%s`", target$output_dir),
    sprintf("- Source final object: `%s`", basename(object_path)),
    sprintf("- Removal basis: %s", target$removal_basis),
    sprintf("- Cluster kind used for matching: `%s`", target$cluster_kind),
    sprintf("- Removed clusters: `%s`", paste(sprintf("c%s", remove_clusters_chr), collapse = ", ")),
    sprintf("- Cell matching mode: `%s` | rate=%s", best_match$mode, sprintf("%.4f", best_match$rate)),
    sprintf("- Cells removed: %d / %d", length(removed_cells), ncol(obj)),
    sprintf("- Filtered final RDS: `%s`", basename(filtered_rds_path)),
    "",
    "## Files",
    "",
    sprintf("- [`%s`](%s)", basename(filtered_rds_path), basename(filtered_rds_path)),
    sprintf("- [`%s`](%s)", basename(removed_cells_path), basename(removed_cells_path)),
    sprintf("- [`%s`](%s)", basename(cluster_summary_path), basename(cluster_summary_path)),
    sprintf("- [`%s`](%s)", basename(meta_counts_path), basename(meta_counts_path)),
    sprintf("- [`%s`](%s)", basename(summary_path), basename(summary_path)),
    ""
  )
  if (file.exists(review_readme)) {
    readme_lines <- c(
      readme_lines,
      "## Related outlier review",
      "",
      sprintf("- [`README.md`](../figures/llm_outlier_review_20260507/README.md)"),
      ""
    )
  }
  writeLines(readme_lines, file.path(removal_dir, "README.md"))

  summary_dt
}

cat("=== Tissue Comparison Posthoc Cluster Removal (2026-05-07) ===\n")
rows <- lapply(REMOVAL_TARGETS, function(target) {
  tryCatch(
    process_removal_target(target),
    error = function(e) {
      data.table(
        lineage_tag = target$lineage_tag,
        output_dir = target$output_dir,
        cluster_kind = target$cluster_kind,
        remove_clusters = paste(target$remove_clusters, collapse = ","),
        removal_basis = target$removal_basis,
        original_cells = NA_integer_,
        removed_cells = NA_integer_,
        retained_cells = NA_integer_,
        match_mode = NA_character_,
        match_rate = NA_character_,
        filtered_rds = NA_character_,
        removed_cells_tsv = NA_character_,
        cluster_summary_tsv = NA_character_,
        meta_counts_tsv = NA_character_,
        error = conditionMessage(e)
      )
    }
  )
})
summary_dt <- rbindlist(rows, fill = TRUE)
print(summary_dt)
if ("error" %in% colnames(summary_dt) && any(!is.na(summary_dt$error) & nzchar(summary_dt$error))) {
  stop("One or more posthoc removal targets failed")
}
cat("[DONE] Posthoc filtered final objects created successfully.\n")
