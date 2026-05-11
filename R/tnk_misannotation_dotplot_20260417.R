#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

OBJECT_PATH <- "/home/h2048/data/R/0414/tnk_tissue_comparison_v2_6_3_20260414_relabel_helper/tnk_tissue_comparison_final.rds"
OUTPUT_DIR <- "/home/h2048/data/R/0414/tnk_tissue_comparison_v2_6_3_20260414_relabel_helper/figures/tnk_misannotation_dotplots_20260417"
CLUSTER_COL <- "CHOIR_clusters_0.2"
LABEL_COL <- "cell_type_L3"

suspect_tbl <- data.frame(
  cluster_id = c("23", "25", "28", "31"),
  dominant_wrong_label = c("CD8 Temra", "CD8 Trm", "CD4 Tcm / CD4 Naive", "CD8 Trm"),
  corrected_interpretation = c(
    "CD8 cytotoxic T cells with alveolar/secretory RNA contamination",
    "stress / heat-shock activated mixed T-cell state",
    "CD4 naive-Tcm-like T cells with goblet / secretory contamination",
    "CD8 T cells with severe secretory / epithelial RNA contamination"
  ),
  evidence_note = c(
    "LLM: mixed alignment; strong SFTPB/SFTPC-driven contamination signal overlaid on CD8 Temra program.",
    "LLM: mixed alignment; stress-activation genes dominate while canonical Trm / cytotoxic markers are suppressed.",
    "LLM: mixed alignment; CCR7/SELL/LEF1 coexist with MUC5AC/TFF3/KRT19 goblet-like program.",
    "LLM: low alignment; ZG16B/LTF/PZP/SFTPC-like secretory program conflicts with annotated CD8 Trm/Temra state."
  ),
  stringsAsFactors = FALSE
)

suspect_tbl$reference_labels <- list(
  c("CD8 Temra"),
  c("CD8 Trm", "CD4 Th1"),
  c("CD4 Tcm", "CD4 Naive"),
  c("CD8 Trm", "CD8 Temra")
)

suspect_tbl$wrong_genes <- list(
  c("CD8B", "NKG7", "PRF1", "GZMB", "GNLY", "FGFBP2", "CX3CR1", "FCGR3A"),
  c("ITGAE", "CXCR6", "ZNF683", "XCL1", "IFNG", "TBX21", "PRF1", "GZMB"),
  c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "LTB", "MAL", "MALAT1"),
  c("ITGAE", "CXCR6", "ZNF683", "FGFBP2", "PRF1", "GZMB", "NKG7", "GNLY", "FCER1G", "SH2D1B")
)

suspect_tbl$corrected_genes <- list(
  c("SFTPB", "SFTPC", "STATH", "LTF", "ZG16B", "BPIFB2", "TFF3", "PIGR"),
  c("HSPA1A", "HSPA1B", "DNAJB1", "JUN", "FOS", "EGR1", "TNF", "RHOB", "ATF3", "CDKN1A"),
  c("MUC5AC", "TFF3", "KRT19", "PIGR", "BPIFB1", "BPIFB2", "STATH", "SCGB3A1"),
  c("ZG16B", "LTF", "BPIFB2", "PZP", "SFTPC", "STATH", "TFF3", "SCGB3A1")
)

make_grouped_dotplot <- function(obj_sub, group_col, wrong_genes, corrected_genes, title, subtitle, file_stub) {
  wrong_present <- intersect(wrong_genes, rownames(obj_sub))
  corrected_present <- intersect(corrected_genes, rownames(obj_sub))

  if (length(wrong_present) < 2 || length(corrected_present) < 2) {
    message(sprintf("[SKIP] %s: insufficient marker genes present", file_stub))
    return(NULL)
  }

  features_present <- c(wrong_present, corrected_present)
  split_at <- length(wrong_present)

  p <- DotPlot(
    object = obj_sub,
    features = features_present,
    group.by = group_col,
    cols = c("#f3f4f6", "#2563eb"),
    dot.scale = 6,
    scale = FALSE
  ) +
    RotatedAxis() +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      color = "Avg.Exp",
      size = "Pct.Exp",
      caption = "左侧为原标注支持基因；右侧为修正解释 / 污染来源支持基因"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0)
    )

  if (split_at > 0 && split_at < length(features_present)) {
    p <- p + geom_vline(xintercept = split_at + 0.5, linetype = "dashed", color = "grey55")
  }

  width <- max(10, 0.45 * length(features_present) + 4)
  height <- 5.5
  ggsave(file.path(OUTPUT_DIR, paste0(file_stub, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(OUTPUT_DIR, paste0(file_stub, ".pdf")), p, width = width, height = height, device = cairo_pdf)

  data.table(
    file_stub = file_stub,
    n_wrong_genes_present = length(wrong_present),
    n_corrected_genes_present = length(corrected_present),
    wrong_genes_present = paste(wrong_present, collapse = ", "),
    corrected_genes_present = paste(corrected_present, collapse = ", ")
  )
}

message("[1/5] Loading TNK object ...")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(OBJECT_PATH)
if (!(CLUSTER_COL %in% colnames(obj@meta.data))) {
  stop(sprintf("Cluster column not found: %s", CLUSTER_COL))
}
if (!(LABEL_COL %in% colnames(obj@meta.data))) {
  stop(sprintf("Label column not found: %s", LABEL_COL))
}
if ("RNA" %in% names(obj@assays)) {
  DefaultAssay(obj) <- "RNA"
}
obj@meta.data[[CLUSTER_COL]] <- as.character(obj@meta.data[[CLUSTER_COL]])
obj@meta.data[[LABEL_COL]] <- as.character(obj@meta.data[[LABEL_COL]])

message("[2/5] Building suspect cluster summary ...")
cluster_counts <- table(obj@meta.data[[CLUSTER_COL]], obj@meta.data[[LABEL_COL]])

summary_list <- vector("list", length = nrow(suspect_tbl))
plot_logs <- vector("list", length = nrow(suspect_tbl) + 1L)

for (i in seq_len(nrow(suspect_tbl))) {
  cluster_id <- suspect_tbl$cluster_id[[i]]
  target_cells <- rownames(obj@meta.data)[obj@meta.data[[CLUSTER_COL]] == cluster_id]
  ref_labels <- suspect_tbl$reference_labels[[i]]
  reference_cells <- rownames(obj@meta.data)[
    obj@meta.data[[CLUSTER_COL]] != cluster_id & obj@meta.data[[LABEL_COL]] %in% ref_labels
  ]

  cluster_tab <- sort(cluster_counts[cluster_id, ], decreasing = TRUE)
  cluster_tab <- cluster_tab[cluster_tab > 0]

  summary_list[[i]] <- data.table(
    cluster_id = cluster_id,
    target_n_cells = length(target_cells),
    reference_n_cells = length(reference_cells),
    dominant_wrong_label = suspect_tbl$dominant_wrong_label[[i]],
    corrected_interpretation = suspect_tbl$corrected_interpretation[[i]],
    reference_labels = paste(ref_labels, collapse = "; "),
    top_L3_breakdown = paste(names(cluster_tab), cluster_tab, sep = "=", collapse = "; "),
    evidence_note = suspect_tbl$evidence_note[[i]]
  )

  if (length(target_cells) == 0 || length(reference_cells) == 0) {
    message(sprintf("[SKIP] cluster %s: target or reference cells missing", cluster_id))
    next
  }

  obj_sub <- subset(obj, cells = c(target_cells, reference_cells))
  obj_sub$misannotation_view <- ifelse(
    Cells(obj_sub) %in% target_cells,
    paste0("CHOIR_", cluster_id, " target"),
    paste0("reference: ", paste(ref_labels, collapse = " + "))
  )
  obj_sub$misannotation_view <- factor(
    obj_sub$misannotation_view,
    levels = c(
      paste0("CHOIR_", cluster_id, " target"),
      paste0("reference: ", paste(ref_labels, collapse = " + "))
    )
  )

  plot_logs[[i]] <- cbind(
    data.table(cluster_id = cluster_id),
    make_grouped_dotplot(
      obj_sub = obj_sub,
      group_col = "misannotation_view",
      wrong_genes = suspect_tbl$wrong_genes[[i]],
      corrected_genes = suspect_tbl$corrected_genes[[i]],
      title = sprintf("TNK misannotation check: CHOIR_%s", cluster_id),
      subtitle = sprintf(
        "reference label(s): %s | interpreted as: %s",
        paste(ref_labels, collapse = " + "),
        suspect_tbl$corrected_interpretation[[i]]
      ),
      file_stub = sprintf("dotplot_cluster_%s_wrong_vs_corrected", cluster_id)
    )
  )
}

summary_dt <- rbindlist(summary_list, fill = TRUE)
fwrite(summary_dt, file.path(OUTPUT_DIR, "misannotation_cluster_summary.tsv"), sep = "\t")

message("[3/5] Drawing overview dotplot across suspect clusters ...")
overview_cells <- rownames(obj@meta.data)[obj@meta.data[[CLUSTER_COL]] %in% suspect_tbl$cluster_id]
obj_overview <- subset(obj, cells = overview_cells)
cluster_label_map <- setNames(
  paste0(
    "CHOIR_", suspect_tbl$cluster_id,
    "\n", suspect_tbl$dominant_wrong_label
  ),
  suspect_tbl$cluster_id
)
obj_overview$misannotation_cluster <- unname(cluster_label_map[obj_overview@meta.data[[CLUSTER_COL]]])
obj_overview$misannotation_cluster <- factor(
  obj_overview$misannotation_cluster,
  levels = unname(cluster_label_map)
)

overview_genes <- unique(c(
  "CD8B", "PRF1", "GZMB", "FGFBP2",
  "ITGAE", "CXCR6", "ZNF683",
  "CCR7", "SELL", "TCF7", "LEF1",
  "HSPA1A", "HSPA1B", "DNAJB1", "JUN", "FOS",
  "SFTPB", "SFTPC", "MUC5AC", "TFF3", "KRT19",
  "ZG16B", "LTF", "BPIFB2", "PZP", "STATH"
))
overview_genes <- intersect(overview_genes, rownames(obj_overview))

if (length(overview_genes) >= 8) {
  p_overview <- DotPlot(
    object = obj_overview,
    features = overview_genes,
    group.by = "misannotation_cluster",
    cols = c("#f3f4f6", "#7c3aed"),
    dot.scale = 6,
    scale = FALSE
  ) +
    RotatedAxis() +
    labs(
      title = "TNK misannotation overview across suspect CHOIR clusters",
      subtitle = "一图总览：左半偏原标注 T-cell markers，右半偏修正解释 / 污染来源 markers",
      x = NULL,
      y = NULL,
      color = "Avg.Exp",
      size = "Pct.Exp"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )

  width_overview <- max(12, 0.42 * length(overview_genes) + 4)
  ggsave(file.path(OUTPUT_DIR, "dotplot_overview_suspect_clusters.png"), p_overview, width = width_overview, height = 6, dpi = 320)
  ggsave(file.path(OUTPUT_DIR, "dotplot_overview_suspect_clusters.pdf"), p_overview, width = width_overview, height = 6, device = cairo_pdf)
  plot_logs[[length(plot_logs)]] <- data.table(
    cluster_id = "overview",
    file_stub = "dotplot_overview_suspect_clusters",
    n_wrong_genes_present = NA_integer_,
    n_corrected_genes_present = NA_integer_,
    wrong_genes_present = paste(overview_genes, collapse = ", "),
    corrected_genes_present = NA_character_
  )
}

message("[4/5] Writing plot log and README ...")
plot_log_dt <- rbindlist(plot_logs, fill = TRUE)
if (nrow(plot_log_dt) > 0) {
  fwrite(plot_log_dt, file.path(OUTPUT_DIR, "dotplot_generation_summary.tsv"), sep = "\t")
}

readme_lines <- c(
  "# TNK misannotation dotplots",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Included suspect clusters",
  "- CHOIR_23: CD8 Temra-like cluster with alveolar / secretory RNA contamination signals (e.g. SFTPB/SFTPC).",
  "- CHOIR_25: stress / heat-shock dominated T-cell state, mismatching canonical CD8 Trm labeling.",
  "- CHOIR_28: CD4 Naive/Tcm-like cluster mixed with goblet / secretory markers (e.g. MUC5AC/TFF3/KRT19).",
  "- CHOIR_31: low-alignment CD8 cluster with strong secretory / epithelial contamination signatures (e.g. ZG16B/LTF/PZP/SFTPC).",
  "",
  "## Plot design",
  "- Each `dotplot_cluster_*` figure compares one suspect CHOIR cluster against a reference pool built from its dominant original label(s).",
  "- Gene order is always: original-label-support markers first, corrected-interpretation / contamination markers second.",
  "- `dotplot_overview_suspect_clusters.*` provides a cross-cluster overview using a compact shared marker panel.",
  "",
  "## Files",
  "- `misannotation_cluster_summary.tsv`: cluster-level summary, dominant labels, and evidence notes.",
  "- `dotplot_generation_summary.tsv`: genes actually found in the object and used for each figure.",
  "- `dotplot_cluster_*_wrong_vs_corrected.*`: per-cluster target-vs-reference dotplots.",
  "- `dotplot_overview_suspect_clusters.*`: combined overview across all suspect clusters.",
  ""
)
writeLines(readme_lines, file.path(OUTPUT_DIR, "README.md"))

message("[5/5] Done.")
message(sprintf("Output directory: %s", OUTPUT_DIR))