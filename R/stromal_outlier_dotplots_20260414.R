suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

save_plot_dual <- function(plot_obj, out_prefix, width = 11, height = 6) {
  ggsave(paste0(out_prefix, ".pdf"), plot = plot_obj, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), plot = plot_obj, width = width, height = height, dpi = 300, limitsize = FALSE)
}

find_choir_col <- function(obj) {
  preferred <- "CHOIR_clusters_0.2"
  meta_cols <- colnames(obj@meta.data)
  if (preferred %in% meta_cols) {
    return(preferred)
  }
  fallback <- grep("^CHOIR_clusters", meta_cols, value = TRUE)
  if (length(fallback) == 0) {
    stop("No CHOIR cluster column found in object metadata.")
  }
  fallback[[1]]
}

selected_marker_table <- function(marker_path, genes) {
  dt <- fread(marker_path)
  dt[, gene := as.character(gene)]
  dt <- dt[gene %in% genes, .(gene, avg_log2FC, pct.1, pct.2, padj)]
  dt[, direction := ifelse(avg_log2FC >= 0, "up_in_outlier", "down_in_outlier")]
  dt[, gene := factor(gene, levels = genes)]
  dt <- dt[order(gene)]
  dt[, gene := as.character(gene)]
  dt
}

prepare_feature_groups <- function(feature_groups, gene_universe) {
  out <- lapply(feature_groups, function(x) unique(x[x %in% gene_universe]))
  out[lengths(out) > 0]
}

plot_one_outlier <- function(obj, choir_col, cluster_id, feature_groups, marker_path,
                             lineage_key, lineage_label, out_dir, summary_dir) {
  cluster_id_chr <- as.character(cluster_id)
  meta <- obj@meta.data
  cells_in_cluster <- rownames(meta)[as.character(meta[[choir_col]]) == cluster_id_chr]
  if (length(cells_in_cluster) == 0) {
    stop(sprintf("No cells found for %s cluster %s", lineage_key, cluster_id_chr))
  }

  feature_groups_use <- prepare_feature_groups(feature_groups, rownames(obj))
  genes_use <- unique(unlist(feature_groups_use, use.names = FALSE))
  if (length(genes_use) == 0) {
    stop(sprintf("No selected genes are present in the object for %s cluster %s", lineage_key, cluster_id_chr))
  }

  plot_group <- ifelse(as.character(meta[[choir_col]]) == cluster_id_chr,
                       sprintf("cluster_%s", cluster_id_chr),
                       "rest")
  obj$outlier_plot_group_20260414 <- factor(
    plot_group,
    levels = c(sprintf("cluster_%s", cluster_id_chr), "rest")
  )

  summary_dt <- selected_marker_table(marker_path, genes_use)
  summary_path <- file.path(summary_dir, sprintf("%s_cluster_%s_selected_deg_pct.csv", lineage_key, cluster_id_chr))
  fwrite(summary_dt, summary_path)

  p <- DotPlot(
    obj,
    features = feature_groups_use,
    group.by = "outlier_plot_group_20260414",
    assay = DefaultAssay(obj)
  ) +
    RotatedAxis() +
    labs(
      title = sprintf("%s | CHOIR cluster %s vs rest", lineage_label, cluster_id_chr),
      subtitle = sprintf("group.by = %s | n_outlier = %d | genes selected from %s", choir_col, length(cells_in_cluster), basename(dirname(marker_path))),
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

  out_prefix <- file.path(out_dir, sprintf("%s_cluster_%s_vs_rest_dotplot", lineage_key, cluster_id_chr))
  save_plot_dual(p, out_prefix, width = 11, height = 6)

  message(sprintf("[OK] %s cluster %s -> %s", lineage_key, cluster_id_chr, out_prefix))
  invisible(list(plot = out_prefix, summary = summary_path, n_outlier = length(cells_in_cluster)))
}

lineage_specs <- list(
  endothelial = list(
    label = "Stromal Endothelial",
    rds = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414/stromal_endothelial_tissue_comparison_final.rds",
    root = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414",
    targets = list(
      `5` = list(
        marker_path = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414/reports/choir/ofa_5/markers.csv",
        feature_groups = list(
          `outlier-up` = c("FBLN1", "SFRP2", "COL1A2", "DCN", "LUM", "COL1A1"),
          `endothelial-low` = c("TM4SF1", "VWF", "PECAM1", "CLDN5")
        )
      ),
      `33` = list(
        marker_path = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414/reports/choir/ofa_33/markers.csv",
        feature_groups = list(
          `outlier-up` = c("ELTD1", "MSMB", "SCGB1A1", "SCGB3A1", "PRKCDBP", "GPR116"),
          `endothelial-context` = c("VWF", "PECAM1", "RNASE1", "TM4SF1")
        )
      )
    )
  ),
  fibroblast = list(
    label = "Stromal Fibroblast",
    rds = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_0_20260414/stromal_fibroblast_tissue_comparison_final.rds",
    root = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_0_20260414",
    targets = list(
      `23` = list(
        marker_path = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_0_20260414/reports/choir/ofa_23/markers.csv",
        feature_groups = list(
          `outlier-up` = c("SERPINA3", "SFTPC", "SFTPB", "ATP6V0C", "NME1-NME2", "MT1X"),
          `fibro-low` = c("COL1A1", "COL3A1", "COL1A2", "IGFBP5")
        )
      ),
      `59` = list(
        marker_path = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_0_20260414/reports/choir/ofa_59/markers.csv",
        feature_groups = list(
          `epithelial-like` = c("SERPINA3", "SFTPC", "SFTPA1", "SFTPB", "SCGB3A1", "DBNDD2"),
          `fibro-context` = c("C7", "LAMA2", "ELN", "DCN")
        )
      )
    )
  ),
  smc = list(
    label = "Stromal SMC/Pericyte",
    rds = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_0_20260414/stromal_smc_tissue_comparison_final.rds",
    root = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_0_20260414",
    targets = list(
      `6` = list(
        marker_path = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_0_20260414/reports/choir/ofa_6/markers.csv",
        feature_groups = list(
          `neuronal-like` = c("MUSTN1", "DGKB", "NMNAT2", "GRID2", "RGS6", "SDK1"),
          `smc-context` = c("RGS5", "CNN1", "DES", "PDE3A")
        )
      )
    )
  )
)

all_outputs <- list()
for (lineage_key in names(lineage_specs)) {
  spec <- lineage_specs[[lineage_key]]
  message(sprintf("=== %s ===", lineage_key))
  obj <- readRDS(spec$rds)
  if ("RNA" %in% names(obj@assays)) {
    DefaultAssay(obj) <- "RNA"
  }
  choir_col <- find_choir_col(obj)
  out_dir <- file.path(spec$root, "figures", "outlier_dotplots_20260414")
  summary_dir <- file.path(spec$root, "reports", "choir", "outlier_dotplots_20260414")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  lineage_outputs <- list()
  for (cluster_id in names(spec$targets)) {
    target <- spec$targets[[cluster_id]]
    lineage_outputs[[cluster_id]] <- plot_one_outlier(
      obj = obj,
      choir_col = choir_col,
      cluster_id = cluster_id,
      feature_groups = target$feature_groups,
      marker_path = target$marker_path,
      lineage_key = lineage_key,
      lineage_label = spec$label,
      out_dir = out_dir,
      summary_dir = summary_dir
    )
  }
  all_outputs[[lineage_key]] <- lineage_outputs
  rm(obj)
  invisible(gc())
}

manifest_rows <- rbindlist(
  lapply(names(all_outputs), function(lineage_key) {
    lineage_outputs <- all_outputs[[lineage_key]]
    rbindlist(lapply(names(lineage_outputs), function(cluster_id) {
      rec <- lineage_outputs[[cluster_id]]
      data.table(
        lineage = lineage_key,
        cluster_id = cluster_id,
        plot_prefix = rec$plot,
        summary_csv = rec$summary,
        n_outlier = rec$n_outlier
      )
    }))
  }),
  fill = TRUE
)
manifest_path <- "/home/h2048/data/R/0414/stromal_outlier_dotplots_20260414_manifest.tsv"
fwrite(manifest_rows, manifest_path, sep = "\t")
message(sprintf("[OK] manifest -> %s", manifest_path))
