#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

OBJECT_PATH <- "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413/tnk_tissue_comparison_final.rds"
CHOIR_SCREEN_PATH <- "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413/reports/llm_choir_discovery_screen.tsv"
OFA_SCREEN_PATH <- "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413/reports/llm_ofa_discovery_screen.tsv"
OUTPUT_DIR <- "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413/figures/llm_focus_dotplots_20260413"
CLUSTER_COL <- "CHOIR_clusters_0.2"

DISCOVERY_LINEAGE_GENES <- c(
  "CCR7", "SELL", "IL7R", "CXCR5", "FOXP3", "IL2RA", "CTLA4",
  "CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "GZMK", "FGFBP2",
  "FCGR3A", "GNLY", "XCL1", "CXCR6", "ZNF683", "IL17A", "IL26",
  "TRDC", "TRGC2", "KLRB1", "ZBTB16", "KIT", "IL23R", "AHR"
)

DISCOVERY_STATE_GENES <- c(
  "HSPA1A", "HSPA1B", "HSPA6", "DNAJB1", "MT1X", "MT2A",
  "MKI67", "TYMS", "BIRC5", "UBE2C", "KIF18B", "CDC45",
  "CSF2", "TNFSF11", "RGS1", "LMNA", "PHLDA1", "DUSP4"
)

OUTLIER_CONTAM_GENES <- c(
  "SCGB1A1", "SCGB3A1", "BPIFB1", "BPIFA1", "TFF3", "MSMB",
  "KRT19", "MGP", "S100A8", "S100A9", "TPSAB1", "TPSB2",
  "IGKC", "IGLC1", "IGHG3", "HBB", "HBA2", "SLPI", "SERPINB3",
  "PZP", "KCNQ5"
)

OUTLIER_LINEAGE_QC_GENES <- c(
  "CD4", "IL7R", "CCR7", "SELL", "CXCR5", "FOXP3", "CD8B",
  "NKG7", "FGFBP2", "FCGR3A", "GNLY", "TRDC", "TRGC2",
  "HSPA1A", "HSPA6", "DNAJB1", "MKI67", "TYMS", "BIRC5"
)

CONFLICT_OVERVIEW_GENES <- c(
  "CD4", "IL7R", "CCR7", "SELL", "CD8B", "NKG7", "FGFBP2",
  "FCGR3A", "GNLY", "TRDC", "HSPA1A", "HSPA6", "DNAJB1",
  "SCGB1A1", "SCGB3A1", "BPIFB1", "MSMB", "S100A8", "TPSAB1",
  "IGKC", "HBB", "MKI67", "TYMS", "MT1X", "IL17A", "IL26"
)

read_llm_screen <- function(path, source_tag) {
  fread(path) %>%
    mutate(
      cluster_id = stringr::str_match(record_label, "cluster_(\\d+)_vs_rest")[, 2] %>% as.integer(),
      source_tag = source_tag
    ) %>%
    filter(!is.na(cluster_id)) %>%
    transmute(
      cluster_id,
      !!paste0(source_tag, "_annotation") := annotation_label,
      !!paste0(source_tag, "_class") := biological_signal_class,
      !!paste0(source_tag, "_confidence") := confidence,
      !!paste0(source_tag, "_short_call") := short_call,
      !!paste0(source_tag, "_evidence") := evidence_summary
    )
}

class_support_label <- function(primary_source, secondary_source, class_name) {
  if (isTRUE(primary_source) && isTRUE(secondary_source)) return("CHOIR+OFA")
  if (isTRUE(primary_source)) return("CHOIR")
  if (isTRUE(secondary_source)) return("OFA")
  NA_character_
}

plot_focus_dotplot <- function(obj, cluster_df, focus_name, genes, file_stub, subtitle) {
  if (nrow(cluster_df) == 0) {
    message(sprintf("[SKIP] %s: no clusters", focus_name))
    return(invisible(NULL))
  }

  keep_clusters <- as.character(cluster_df$cluster_id)
  keep_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[CLUSTER_COL]]) %in% keep_clusters]
  if (length(keep_cells) == 0) {
    message(sprintf("[SKIP] %s: no cells", focus_name))
    return(invisible(NULL))
  }

  obj_sub <- subset(obj, cells = keep_cells)
  label_map <- setNames(cluster_df$plot_label, as.character(cluster_df$cluster_id))
  obj_sub$llm_focus_cluster <- unname(label_map[as.character(obj_sub@meta.data[[CLUSTER_COL]])])
  obj_sub$llm_focus_cluster <- factor(
    obj_sub$llm_focus_cluster,
    levels = cluster_df$plot_label
  )

  genes_present <- intersect(genes, rownames(obj_sub))
  if (length(genes_present) < 3) {
    message(sprintf("[SKIP] %s: fewer than 3 genes present", focus_name))
    return(invisible(NULL))
  }

  plot_width <- max(12, 0.42 * length(genes_present) + 4)
  plot_height <- max(8, 0.38 * nrow(cluster_df) + 2.5)

  p <- DotPlot(
    object = obj_sub,
    features = genes_present,
    group.by = "llm_focus_cluster",
    cols = c("#f3f4f6", "#2563eb"),
    dot.scale = 6
  ) +
    RotatedAxis() +
    labs(
      title = sprintf("TNK LLM-focused dotplot: %s", focus_name),
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      color = "Avg.Exp",
      size = "Pct.Exp"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(OUTPUT_DIR, paste0(file_stub, ".png")), p, width = plot_width, height = plot_height, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, paste0(file_stub, ".pdf")), p, width = plot_width, height = plot_height, device = cairo_pdf)

  data.table(
    focus_name = focus_name,
    file_stub = file_stub,
    n_clusters = nrow(cluster_df),
    n_cells = length(keep_cells),
    n_genes_requested = length(genes),
    n_genes_present = length(genes_present),
    genes_present = paste(genes_present, collapse = ", ")
  )
}

message("[1/5] Loading object and LLM screens ...")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(OBJECT_PATH)
if (!(CLUSTER_COL %in% colnames(obj@meta.data))) {
  stop(sprintf("Cluster column not found: %s", CLUSTER_COL))
}
if ("RNA" %in% names(obj@assays)) {
  DefaultAssay(obj) <- "RNA"
}
obj@meta.data[[CLUSTER_COL]] <- as.character(obj@meta.data[[CLUSTER_COL]])

choir_screen <- read_llm_screen(CHOIR_SCREEN_PATH, "choir")
ofa_screen <- read_llm_screen(OFA_SCREEN_PATH, "ofa")

message("[2/5] Merging CHOIR and OFA focus calls ...")
focus_tbl <- full_join(choir_screen, ofa_screen, by = "cluster_id") %>%
  mutate(
    choir_discovery = choir_class == "potential_discovery",
    ofa_discovery = ofa_class == "potential_discovery",
    choir_outlier = choir_class == "likely_outlier",
    ofa_outlier = ofa_class == "likely_outlier",
    any_discovery = choir_discovery | ofa_discovery,
    any_outlier = choir_outlier | ofa_outlier,
    focus_group = case_when(
      any_discovery & any_outlier ~ "conflicted",
      any_discovery ~ "discovery",
      any_outlier ~ "outlier",
      TRUE ~ NA_character_
    ),
    discovery_support = vapply(
      seq_len(n()),
      function(i) class_support_label(choir_discovery[i], ofa_discovery[i], "potential_discovery"),
      character(1)
    ),
    outlier_support = vapply(
      seq_len(n()),
      function(i) class_support_label(choir_outlier[i], ofa_outlier[i], "likely_outlier"),
      character(1)
    ),
    dominant_annotation = dplyr::coalesce(choir_annotation, ofa_annotation),
    dominant_short_call = dplyr::coalesce(choir_short_call, ofa_short_call),
    plot_label = paste0("c", cluster_id)
  ) %>%
  filter(!is.na(focus_group)) %>%
  arrange(focus_group, cluster_id)

if (nrow(focus_tbl) == 0) {
  stop("No potential_discovery / likely_outlier clusters were found in the LLM screens.")
}

summary_tbl <- focus_tbl %>%
  mutate(
    cluster_size = as.integer(table(obj@meta.data[[CLUSTER_COL]])[as.character(cluster_id)]),
    cluster_size = ifelse(is.na(cluster_size), 0L, cluster_size)
  ) %>%
  select(
    cluster_id, cluster_size, focus_group, dominant_annotation, dominant_short_call,
    choir_class, ofa_class, discovery_support, outlier_support,
    choir_confidence, ofa_confidence, choir_short_call, ofa_short_call,
    choir_evidence, ofa_evidence
  )

fwrite(summary_tbl, file.path(OUTPUT_DIR, "llm_focus_cluster_summary.tsv"), sep = "\t")

message("[3/5] Building focused cluster groups ...")
discovery_tbl <- focus_tbl %>% filter(focus_group == "discovery") %>% arrange(cluster_id)
outlier_tbl <- focus_tbl %>% filter(focus_group == "outlier") %>% arrange(cluster_id)
conflict_tbl <- focus_tbl %>% filter(focus_group == "conflicted") %>% arrange(cluster_id)

message(sprintf("  discovery-only clusters: %d", nrow(discovery_tbl)))
message(sprintf("  outlier-only clusters:   %d", nrow(outlier_tbl)))
message(sprintf("  conflicted clusters:     %d", nrow(conflict_tbl)))

message("[4/5] Drawing dotplots ...")
plot_logs <- bind_rows(
  plot_focus_dotplot(
    obj, discovery_tbl, "Discovery-only clusters: lineage / identity",
    DISCOVERY_LINEAGE_GENES,
    "dotplot_discovery_lineage",
    "Potential discoveries flagged by CHOIR/OFA without any outlier call"
  ),
  plot_focus_dotplot(
    obj, discovery_tbl, "Discovery-only clusters: stress / proliferation / state",
    DISCOVERY_STATE_GENES,
    "dotplot_discovery_state",
    "Stress-response, metal-response, cytokine, and cell-cycle programs"
  ),
  plot_focus_dotplot(
    obj, outlier_tbl, "Outlier-only clusters: contamination / exclusion markers",
    OUTLIER_CONTAM_GENES,
    "dotplot_outlier_contamination",
    "Likely-removal clusters flagged without any competing discovery call"
  ),
  plot_focus_dotplot(
    obj, outlier_tbl, "Outlier-only clusters: lineage mismatch / QC markers",
    OUTLIER_LINEAGE_QC_GENES,
    "dotplot_outlier_lineage_qc",
    "Lineage-conflict, stress-response, and proliferation markers"
  ),
  plot_focus_dotplot(
    obj, conflict_tbl, "Conflicted clusters: discovery vs outlier disagreement",
    CONFLICT_OVERVIEW_GENES,
    "dotplot_conflicted_overview",
    "Clusters called discovery by one screen but outlier by the other"
  )
)

if (nrow(plot_logs) > 0) {
  fwrite(plot_logs, file.path(OUTPUT_DIR, "dotplot_generation_summary.tsv"), sep = "\t")
}

message("[5/5] Writing quick README ...")
readme_lines <- c(
  "# TNK LLM-focused dotplots",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Cluster grouping logic",
  "- discovery: flagged as `potential_discovery` by CHOIR and/or OFA, and **not** flagged as outlier by either screen.",
  "- outlier: flagged as `likely_outlier` by CHOIR and/or OFA, and **not** flagged as discovery by either screen.",
  "- conflicted: discovery in one screen but outlier in the other; inspect manually before filtering.",
  "",
  "## Files",
  "- `llm_focus_cluster_summary.tsv`: merged CHOIR/OFA focus calls per cluster",
  "- `dotplot_generation_summary.tsv`: gene panels actually used for each dotplot",
  "- `dotplot_discovery_lineage.*`: discovery-only clusters, lineage/state markers",
  "- `dotplot_discovery_state.*`: discovery-only clusters, stress/proliferation markers",
  "- `dotplot_outlier_contamination.*`: outlier-only clusters, contamination markers",
  "- `dotplot_outlier_lineage_qc.*`: outlier-only clusters, lineage mismatch / QC markers",
  "- `dotplot_conflicted_overview.*`: conflicted clusters (if any)",
  ""
)
writeLines(readme_lines, file.path(OUTPUT_DIR, "README.md"))

message("[DONE] LLM-focused TNK dotplots written to:")
message(sprintf("  %s", OUTPUT_DIR))
