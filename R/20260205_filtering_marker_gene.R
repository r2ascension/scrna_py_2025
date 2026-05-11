#!/usr/bin/env Rscript
# ============================================================================
# Get top-N marker genes per cluster from a Seurat-style DE CSV
# Input CSV columns (required): cluster, gene, p_val_adj, avg_log2FC (or avg_logFC)
# Output: one CSV with top30 per cluster
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------ config ---------------------------------------
in_csv <- '/home/h2048/data/R/0128/stromal_interpret_v1_1/all_markers.csv'
out_csv <- if (length(commandArgs(trailingOnly = TRUE)) >= 2) {
  commandArgs(trailingOnly = TRUE)[2]
} else {
  "/home/h2048/data/R/0128/stromal_interpret_v1_1/top30_per_cluster.csv"
}
top_n <- 30
# -----------------------------------------------------------------------------
dt <- fread(in_csv)
# required
stopifnot(all(c("cluster", "gene") %in% names(dt)))

# choose columns
fc_col   <- if ("avg_log2FC" %in% names(dt)) "avg_log2FC" else if ("avg_logFC" %in% names(dt)) "avg_logFC" else NA_character_
padj_col <- if ("p_val_adj" %in% names(dt)) "p_val_adj" else NA_character_

if (is.na(fc_col)) stop("Need avg_log2FC or avg_logFC")

dt[, cluster := as.character(cluster)]
dt[, gene    := as.character(gene)]

# correct ordering: padj asc, FC desc
if (!is.na(padj_col)) {
  setorderv(dt, c("cluster", padj_col, fc_col), c(1, 1, -1))
} else {
  setorderv(dt, c("cluster", padj_col, fc_col), c(1, 1, -1))
}

top_dt <- dt[, head(.SD, top_n), by = cluster]

# optional: de-dup within cluster
top_dt <- top_dt[!duplicated(paste(cluster, gene))]

fwrite(top_dt, out_csv)
message("Wrote: ", out_csv, " (", nrow(top_dt), " rows)")