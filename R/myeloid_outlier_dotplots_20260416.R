#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(stringr)
})

OBJECT_PATH <- "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414/myeloid_tissue_comparison_final.rds"
CHOIR_STRUCTURED_PATH <- "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414/reports/choir_llm_structured.tsv"
OUTPUT_ROOT <- "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414"
FIGURE_DIR <- file.path(OUTPUT_ROOT, "figures", "outlier_dotplots_20260416")
SUMMARY_DIR <- file.path(OUTPUT_ROOT, "reports", "choir", "outlier_dotplots_20260416")
GENE_NOISE_REGEX <- "^(MT-|RPL|RPS|H3F3A$|H3F3B$|MIR|SNOR|SNORA|SNORD|RN7|RNU|CTD-|AC[0-9]|RP[0-9]|LLNLF|ENSG|RF[0-9]|TRIR$|SEM1$|FAU$|TMSB4X$|EEF1A1$|RPSA$|RPS3A$|RPS6$|RPLP0$|RPLP1$|RPLP2$)"

CURATED_FEATURES <- list(
  `2` = list(
    up = c("FCN1", "SERPINB2", "S100A12", "CD300E", "EREG", "VCAN", "IL1B"),
    down = c("APOC1", "FABP4", "MARCO", "MRC1", "GPNMB")
  ),
  `3` = list(
    up = c("ADORA3", "CCL2", "MERTK", "F13A1", "PLTP", "RNASE1", "FCGR2B"),
    down = c("S100A12", "MME", "SELL", "APOBEC3A", "RND3")
  ),
  `7` = list(
    up = c("CD207", "CD1A", "ALOX15", "CD1C", "FCER1A", "CLEC10A"),
    down = c("MARCO", "CCL18", "FABP4", "S100A9", "C5AR1")
  ),
  `8` = list(
    up = c("PLAU", "ADORA2A", "IRAK2", "PLPP3", "CXCL8", "CSF1", "VEGFA", "ICAM1"),
    down = c("VCAN", "AXL", "C5AR1", "PLIN2", "S100A8")
  )
)

save_plot_dual <- function(plot_obj, out_prefix, width = 12, height = 6) {
  ggsave(paste0(out_prefix, ".pdf"), plot = plot_obj, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), plot = plot_obj, width = width, height = height, dpi = 300, limitsize = FALSE)
}

ordered_unique <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  x[!duplicated(x)]
}

clean_gene_vector <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(character())
  vals <- unlist(strsplit(paste(x[!is.na(x)], collapse = ","), ",\\s*"))
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  ordered_unique(vals)
}

find_choir_col <- function(obj) {
  preferred <- "CHOIR_clusters_0.2"
  meta_cols <- colnames(obj@meta.data)
  if (preferred %in% meta_cols) return(preferred)
  fallback <- grep("^CHOIR_clusters", meta_cols, value = TRUE)
  if (length(fallback) == 0) stop("No CHOIR cluster column found in object metadata.")
  fallback[[1]]
}

read_marker_table <- function(cluster_id) {
  marker_path <- file.path(OUTPUT_ROOT, "reports", "choir", sprintf("ofa_%s", cluster_id), "markers.csv")
  if (!file.exists(marker_path)) stop(sprintf("Marker file not found for cluster %s: %s", cluster_id, marker_path))
  dt <- fread(marker_path)
  dt[, gene := as.character(gene)]
  dt[, padj_use := if ("padj" %in% names(dt)) as.numeric(padj) else as.numeric(p_val_adj)]
  dt[, pct_delta := pct.1 - pct.2]
  dt[, marker_path := marker_path]
  dt
}

read_structured_tsv <- function(path) {
  as.data.table(
    read.delim(
      file = path,
      sep = "\t",
      header = TRUE,
      quote = '"',
      fill = TRUE,
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

select_cluster_genes <- function(marker_dt, llm_key_drivers, curated_spec, n_up = 6, n_down = 5) {
  llm_genes <- clean_gene_vector(llm_key_drivers)
  curated_up <- clean_gene_vector(curated_spec$up)
  curated_down <- clean_gene_vector(curated_spec$down)
  preferred_order <- ordered_unique(c(llm_genes, curated_up, curated_down))

  present_preferred <- preferred_order[preferred_order %in% marker_dt$gene]
  preferred_dt <- marker_dt[match(present_preferred, gene, nomatch = 0)]

  up_genes <- preferred_order[preferred_order %in% preferred_dt[avg_log2FC > 0, gene]]
  down_genes <- preferred_order[preferred_order %in% preferred_dt[avg_log2FC < 0, gene]]

  fallback_up <- marker_dt[
    avg_log2FC > 0 &
      pct_delta >= 0.05 &
      !grepl(GENE_NOISE_REGEX, gene),
    .(gene, avg_log2FC, pct_delta, padj_use)
  ][order(-avg_log2FC, -pct_delta, padj_use), gene]

  fallback_down <- marker_dt[
    avg_log2FC < 0 &
      (pct.2 - pct.1) >= 0.05 &
      !grepl(GENE_NOISE_REGEX, gene),
    .(gene, avg_log2FC, pct_diff = pct.2 - pct.1, padj_use)
  ][order(avg_log2FC, -pct_diff, padj_use), gene]

  up_final <- ordered_unique(c(up_genes, fallback_up))[seq_len(min(n_up, length(ordered_unique(c(up_genes, fallback_up)))))]
  down_final <- ordered_unique(c(down_genes, fallback_down))[seq_len(min(n_down, length(ordered_unique(c(down_genes, fallback_down)))))]

  selection_dt <- marker_dt[gene %in% c(up_final, down_final), .(gene, avg_log2FC, pct.1, pct.2, padj = padj_use)]
  selection_dt[, direction := ifelse(avg_log2FC >= 0, "up_in_cluster", "up_in_rest")]
  selection_dt[, selected_reason := fifelse(
    gene %in% llm_genes, "llm_key_driver",
    fifelse(gene %in% c(curated_up, curated_down), "curated_fallback", "top_deg_fallback")
  )]
  selection_dt[, gene := factor(gene, levels = c(up_final, down_final))]
  selection_dt <- selection_dt[order(gene)]
  selection_dt[, gene := as.character(gene)]

  list(
    up = up_final,
    down = down_final,
    summary = selection_dt
  )
}

plot_one_cluster <- function(obj, choir_col, structured_row, marker_dt, genes_sel) {
  cluster_id_chr <- as.character(structured_row$cluster_id)
  cluster_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[choir_col]]) == cluster_id_chr]
  if (length(cluster_cells) == 0) stop(sprintf("No cells found for cluster %s", cluster_id_chr))

  plot_group <- ifelse(as.character(obj@meta.data[[choir_col]]) == cluster_id_chr,
                       sprintf("cluster_%s", cluster_id_chr),
                       "rest")
  obj$outlier_plot_group_20260416 <- factor(plot_group, levels = c(sprintf("cluster_%s", cluster_id_chr), "rest"))

  feature_groups <- list()
  if (length(genes_sel$up) > 0) feature_groups[["up in cluster"]] <- genes_sel$up
  if (length(genes_sel$down) > 0) feature_groups[["up in rest"]] <- genes_sel$down
  if (length(feature_groups) == 0) stop(sprintf("No plot genes available for cluster %s", cluster_id_chr))

  title_main <- sprintf("Myeloid outlier dotplot | CHOIR cluster %s vs rest", cluster_id_chr)
  subtitle_text <- sprintf(
    "%s | annotated L3: %s | n_cluster=%d",
    structured_row$discovery_assessment,
    structured_row$annotated_l3_dominant,
    length(cluster_cells)
  )

  p <- DotPlot(
    obj,
    features = feature_groups,
    group.by = "outlier_plot_group_20260416",
    assay = DefaultAssay(obj)
  ) +
    RotatedAxis() +
    labs(
      title = title_main,
      subtitle = subtitle_text,
      x = NULL,
      y = NULL,
      color = "avg.exp.scaled",
      size = "pct.exp"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )

  file_stub <- sprintf("myeloid_cluster_%s_vs_rest_dotplot", cluster_id_chr)
  out_prefix <- file.path(FIGURE_DIR, file_stub)
  save_plot_dual(p, out_prefix, width = 11, height = 6)

  summary_dt <- copy(genes_sel$summary)
  summary_dt[, cluster_id := cluster_id_chr]
  summary_dt[, discovery_assessment := structured_row$discovery_assessment]
  summary_dt[, annotated_l3_dominant := structured_row$annotated_l3_dominant]
  summary_dt[, annotated_l3_secondary := structured_row$annotated_l3_secondary]
  summary_dt[, key_drivers_raw := structured_row$key_drivers]
  summary_dt[, integrated_comment := structured_row$integrated_diagnostic_comment]
  fwrite(summary_dt, file.path(SUMMARY_DIR, sprintf("myeloid_cluster_%s_selected_deg_pct.csv", cluster_id_chr)))

  data.table(
    cluster_id = cluster_id_chr,
    discovery_assessment = structured_row$discovery_assessment,
    annotated_l3_dominant = structured_row$annotated_l3_dominant,
    annotated_l3_secondary = structured_row$annotated_l3_secondary,
    n_cluster = length(cluster_cells),
    plot_prefix = out_prefix,
    marker_path = unique(marker_dt$marker_path),
    up_genes = paste(genes_sel$up, collapse = ", "),
    down_genes = paste(genes_sel$down, collapse = ", ")
  )
}

message("[1/4] Loading myeloid object and structured CHOIR outlier calls ...")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(OBJECT_PATH)
if ("RNA" %in% names(obj@assays)) DefaultAssay(obj) <- "RNA"
choir_col <- find_choir_col(obj)
obj@meta.data[[choir_col]] <- as.character(obj@meta.data[[choir_col]])

structured_dt <- read_structured_tsv(CHOIR_STRUCTURED_PATH)
structured_dt[, cluster_id := str_match(comparison, "cluster_(\\d+)_vs_rest")[, 2]]
structured_dt <- structured_dt[
  !is.na(cluster_id) & cluster_id %in% names(CURATED_FEATURES),
  .(cluster_id, comparison, key_drivers, discovery_assessment,
    annotated_l3_dominant, annotated_l3_secondary, integrated_diagnostic_comment)
]
structured_dt[, cluster_id := as.character(cluster_id)]
structured_dt <- structured_dt[order(as.integer(cluster_id))]

message(sprintf("[OK] Target clusters: %s", paste(structured_dt$cluster_id, collapse = ", ")))

message("[2/4] Selecting DEG + pct genes for each target cluster ...")
manifest_list <- list()
for (i in seq_len(nrow(structured_dt))) {
  row_i <- structured_dt[i]
  cid <- row_i$cluster_id
  curated_spec <- CURATED_FEATURES[[cid]]
  if (is.null(curated_spec)) stop(sprintf("No curated feature fallback defined for cluster %s", cid))
  marker_dt <- read_marker_table(cid)
  genes_sel <- select_cluster_genes(marker_dt, row_i$key_drivers, curated_spec)
  manifest_list[[cid]] <- plot_one_cluster(obj, choir_col, row_i, marker_dt, genes_sel)
  message(sprintf("[OK] Cluster %s done | up=%d | down=%d", cid, length(genes_sel$up), length(genes_sel$down)))
}

message("[3/4] Writing manifest and quick README ...")
manifest_dt <- rbindlist(manifest_list, fill = TRUE)
fwrite(manifest_dt, file.path(SUMMARY_DIR, "myeloid_outlier_dotplot_manifest.tsv"), sep = "\t")

readme_lines <- c(
  "# Myeloid outlier dotplots (2026-04-16)",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Object: %s", basename(OBJECT_PATH)),
  sprintf("CHOIR cluster column: %s", choir_col),
  "",
  "## Target clusters",
  paste0("- cluster ", manifest_dt$cluster_id, ": ", manifest_dt$discovery_assessment,
         " | dominant=", manifest_dt$annotated_l3_dominant,
         ifelse(is.na(manifest_dt$annotated_l3_secondary) | !nzchar(manifest_dt$annotated_l3_secondary), "", paste0(" | secondary=", manifest_dt$annotated_l3_secondary))),
  "",
  "## Files",
  "- `myeloid_outlier_dotplot_manifest.tsv`: cluster-level manifest",
  "- `myeloid_cluster_<id>_selected_deg_pct.csv`: selected genes with log2FC / pct.1 / pct.2 / padj",
  "- `myeloid_cluster_<id>_vs_rest_dotplot.{png,pdf}`: dotplot output",
  ""
)
writeLines(readme_lines, file.path(SUMMARY_DIR, "README.md"))

message("[4/4] Done.")
message(sprintf("[OUT] Figures: %s", FIGURE_DIR))
message(sprintf("[OUT] Summaries: %s", SUMMARY_DIR))
