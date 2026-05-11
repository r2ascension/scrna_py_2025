#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260414_v2.R")

OUTPUT_STEM <- "llm_outlier_review_20260507"
PLOT_CLUSTER_COL <- "llm_visual_cluster_id"
TOP_MARKERS_PER_CLUSTER <- 6L
FEATUREPLOT_GENES_PER_CLUSTER <- 4L
DOTPLOT_GENES_PER_CLUSTER <- 4L
LINEAGE_GENES_PER_PANEL <- 4L
CONTAMINATION_FEATUREPLOT_GENES_PER_CLUSTER <- 6L
CONTAMINATION_DOTPLOT_MAX_GENES <- 18L
TECHNICAL_GENE_REGEX <- "^(MT-|RPS|RPL|MALAT1$|NEAT1$|MTRNR)"

sanitize_text <- function(x, default = "") {
  x <- as.character(x)
  x[is.na(x)] <- default
  x <- gsub("[\r\n]+", " ", x)
  trimws(x)
}

normalize_yes <- function(x) {
  val <- tolower(trimws(as.character(x)))
  val %in% c("yes", "y", "true", "1")
}

first_existing_col <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0L) return(NULL)
  hits[[1]]
}

find_final_object_path <- function(output_dir) {
  hits <- list.files(output_dir, pattern = "_final\\.rds$", full.names = TRUE)
  if (length(hits) == 0L) {
    stop(sprintf("No *_final.rds file found under: %s", output_dir))
  }
  hits[[1]]
}

screen_label_from_key <- function(key) {
  switch(
    tolower(key),
    choir = "CHOIR",
    leiden = "Leiden",
    ofa = "OFA",
    toupper(key)
  )
}

build_default_targets <- function() {
  list(
    list(
      lineage_tag = "BCELL",
      output_dir = "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415",
      requested_removal_clusters = list(choir = c(14)),
      focus_genes_by_screen_cluster = list(
        choir = list(
          `14` = c("SFTPC", "SFTPA1", "SFTPA2", "SFTPB")
        )
      )
    ),
    list(
      lineage_tag = "MYELOID",
      output_dir = "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416"
    ),
    list(
      lineage_tag = "EPITHELIAL",
      output_dir = "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun",
      cluster_kind = "leiden",
      requested_removal_clusters = list(ofa = c(14, 17)),
      focus_genes_by_screen_cluster = list(
        ofa = list(
          `14` = c("HBB", "IGKC", "HBA1", "HBA2", "IGLC2"),
          `17` = c("HBB", "HBA1", "HBA2", "IGKC", "IGLC2")
        )
      )
    ),
    list(
      lineage_tag = "STROMAL_ENDOTHELIAL",
      output_dir = "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414",
      requested_removal_clusters = list(ofa = c(6, 52)),
      focus_genes_by_screen_cluster = list(
        ofa = list(
          `6` = c("IGHG3", "JCHAIN", "IGHA1", "IGHA2", "IGHM", "IGLC2"),
          `52` = c("AL445259.1", "PRDM16-DT", "LINC00841", "AC090004.2")
        )
      )
    ),
    list(
      lineage_tag = "STROMAL_FIBROBLAST",
      output_dir = "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414"
    ),
    list(
      lineage_tag = "STROMAL_SMC",
      output_dir = "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414"
    ),
    list(
      lineage_tag = "TNK",
      output_dir = "/home/h2048/data/R/0414/tnk_tissue_comparison_v2_6_3_20260414_relabel_helper",
      focus_genes_by_screen_cluster = list(
        choir = list(
          `23` = c("SFTPB", "SFTPC", "SFTPA1", "SFTPA2"),
          `28` = c("MUC5AC", "TFF3", "SCGB1A1", "SCGB3A1"),
          `31` = c("SCGB1A1", "SCGB3A1", "BPIFA1", "TFF3")
        ),
        ofa = list(
          `23` = c("SFTPB", "SFTPC", "SFTPA1", "SFTPA2"),
          `28` = c("MUC5AC", "TFF3", "SCGB1A1", "SCGB3A1")
        )
      )
    )
  )
}

default_contamination_panel_genes <- function(lineage_tag) {
  switch(
    toupper(lineage_tag),
    BCELL = c("SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SCGB1A1", "SCGB3A1", "BPIFA1", "BPIFB1", "MUC5AC", "MUC5B", "TFF3", "KRT17", "HBB", "HBA1", "HBA2"),
    MYELOID = c("HBB", "HBA1", "HBA2", "IGKC", "IGLC1", "IGLC2", "IGHG3", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SCGB1A1", "TPSAB1", "TPSB2"),
    EPITHELIAL = c("HBB", "HBA1", "HBA2", "IGKC", "IGLC1", "IGLC2", "IGHG1", "IGHG3", "IGHA1", "IGHA2", "IGHM", "JCHAIN", "TPSAB1", "TPSB2", "S100A8", "S100A9"),
    TNK = c("SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SCGB1A1", "SCGB3A1", "BPIFA1", "BPIFB1", "MUC5AC", "MUC5B", "TFF3", "HBB", "HBA1", "HBA2", "IGKC", "IGLC1", "IGLC2", "IGHG1", "IGHG3", "IGHA1", "IGHA2", "IGHM", "JCHAIN", "TPSAB1", "TPSB2"),
    STROMAL_ENDOTHELIAL = c("IGHG3", "IGHA1", "IGHA2", "IGHM", "IGKC", "IGLC1", "IGLC2", "JCHAIN", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SCGB1A1", "SCGB3A1", "HBB", "HBA1", "HBA2", "HSPA1A", "HSPA6"),
    STROMAL_FIBROBLAST = c("GNLY", "NKG7", "CCL4", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "IGKC", "IGLC1", "IGHA1", "IGHA2", "JCHAIN", "GRIA4", "CHRNA6", "HBB", "HBA1", "HBA2"),
    STROMAL_SMC = c("SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "IGKC", "IGLC1", "JCHAIN", "HBB", "HBA1", "HBA2", "GNLY", "NKG7"),
    character()
  )
}

get_requested_removal_clusters <- function(target, screen_key = NULL) {
  cfg <- target$requested_removal_clusters
  if (is.null(cfg)) return(character())
  if (is.null(screen_key)) return(character())
  vals <- cfg[[screen_key]]
  if (is.null(vals)) return(character())
  unique(as.character(vals))
}

get_focus_genes_from_config <- function(target, screen_key, cluster_id) {
  cfg <- target$focus_genes_by_screen_cluster
  if (is.null(cfg)) return(character())
  screen_cfg <- cfg[[screen_key]]
  if (is.null(screen_cfg)) return(character())
  genes <- screen_cfg[[as.character(cluster_id)]]
  if (is.null(genes)) return(character())
  unique(sanitize_text(genes))
}

detect_primary_cluster_kind <- function(output_dir, preferred_kind = NULL) {
  if (!is.null(preferred_kind) && preferred_kind %in% c("choir", "leiden")) {
    return(preferred_kind)
  }
  reports_dir <- file.path(output_dir, "reports")
  if (file.exists(file.path(reports_dir, "llm_choir_discovery_screen.tsv")) &&
      dir.exists(file.path(reports_dir, "choir"))) {
    return("choir")
  }
  if (file.exists(file.path(reports_dir, "llm_leiden_discovery_screen.tsv")) &&
      dir.exists(file.path(reports_dir, "leiden"))) {
    return("leiden")
  }
  NA_character_
}

screen_md_candidates <- function(output_dir, screen_key, cluster_kind = NULL) {
  candidates <- switch(
    tolower(screen_key),
    choir = c("LLM_CHOIR_DISCOVERY_REVIEW.md", "LLM_CHOIR_INTERPRETATION.md", "LLM_INTERPRETATION.md"),
    leiden = c("LLM_LEIDEN_DISCOVERY_REVIEW.md", "LLM_LEIDEN_INTERPRETATION.md", "LLM_INTERPRETATION.md"),
    ofa = c("LLM_OFA_DISCOVERY_REVIEW.md", "LLM_INTERPRETATION.md"),
    c("LLM_INTERPRETATION.md")
  )
  full <- file.path(output_dir, candidates)
  full[file.exists(full)]
}

build_screen_specs <- function(output_dir, preferred_kind = NULL) {
  reports_dir <- file.path(output_dir, "reports")
  primary_kind <- detect_primary_cluster_kind(output_dir, preferred_kind = preferred_kind)
  specs <- list()

  if (!is.na(primary_kind)) {
    primary_screen_path <- file.path(reports_dir, sprintf("llm_%s_discovery_screen.tsv", primary_kind))
    if (file.exists(primary_screen_path)) {
      specs[[length(specs) + 1L]] <- list(
        key = primary_kind,
        label = screen_label_from_key(primary_kind),
        screen_path = primary_screen_path,
        cluster_kind = primary_kind,
        cluster_dir = file.path(reports_dir, primary_kind),
        source_md = screen_md_candidates(output_dir, primary_kind, cluster_kind = primary_kind)
      )
    }
  }

  ofa_screen_path <- file.path(reports_dir, "llm_ofa_discovery_screen.tsv")
  if (file.exists(ofa_screen_path) && !is.na(primary_kind)) {
    specs[[length(specs) + 1L]] <- list(
      key = "ofa",
      label = "OFA",
      screen_path = ofa_screen_path,
      cluster_kind = primary_kind,
      cluster_dir = file.path(reports_dir, primary_kind),
      source_md = screen_md_candidates(output_dir, "ofa", cluster_kind = primary_kind)
    )
  }

  specs
}

read_cluster_assignments <- function(output_dir, cluster_kind) {
  reports_subdir <- file.path(output_dir, "reports", cluster_kind)
  csv_name <- switch(
    tolower(cluster_kind),
    choir = "choir_clusters.csv",
    leiden = "leiden_clusters.csv",
    stop(sprintf("Unsupported cluster kind: %s", cluster_kind))
  )
  csv_path <- file.path(reports_subdir, csv_name)
  if (!file.exists(csv_path)) {
    stop(sprintf("Cluster assignment CSV not found: %s", csv_path))
  }

  dt <- fread(csv_path)
  cell_col <- first_existing_col(dt, c("cell", "cell_id", "barcode"))
  cluster_col <- first_existing_col(dt, c("choir_cluster", "cluster_id", "cluster", "leiden_cluster"))
  if (is.null(cell_col) || is.null(cluster_col)) {
    stop(sprintf("Could not detect cell/cluster columns in: %s", csv_path))
  }

  out <- data.frame(
    cell = sanitize_text(dt[[cell_col]]),
    cluster_id = suppressWarnings(as.integer(as.character(dt[[cluster_col]]))),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$cluster_id) & nzchar(out$cell), , drop = FALSE]
  unique(out)
}

match_cluster_vector <- function(obj_cells, assign_df) {
  exact_idx <- match(obj_cells, assign_df$cell)
  best <- list(
    cluster_vec = assign_df$cluster_id[exact_idx],
    rate = mean(!is.na(exact_idx)),
    mode = "exact"
  )

  try_mode <- function(mode_name, lhs, rhs) {
    idx <- match(lhs, rhs)
    rate <- mean(!is.na(idx))
    list(cluster_vec = assign_df$cluster_id[idx], rate = rate, mode = mode_name)
  }

  if (any(grepl("::", assign_df$cell))) {
    cand <- try_mode("strip_assignment_prefix", obj_cells, sub("^.*::", "", assign_df$cell))
    if (!is.na(cand$rate) && cand$rate > best$rate) best <- cand
  }
  if (any(grepl("::", obj_cells))) {
    cand <- try_mode("strip_object_prefix", sub("^.*::", "", obj_cells), assign_df$cell)
    if (!is.na(cand$rate) && cand$rate > best$rate) best <- cand
  }
  if (any(grepl("::", assign_df$cell)) && any(grepl("::", obj_cells))) {
    cand <- try_mode(
      "strip_both_prefixes",
      sub("^.*::", "", obj_cells),
      sub("^.*::", "", assign_df$cell)
    )
    if (!is.na(cand$rate) && cand$rate > best$rate) best <- cand
  }

  best
}

ensure_expression_ready <- function(obj) {
  if ("RNA" %in% names(obj@assays)) {
    DefaultAssay(obj) <- "RNA"
  }

  has_data <- FALSE
  has_data <- tryCatch({
    dat <- suppressWarnings(GetAssayData(obj, assay = DefaultAssay(obj), slot = "data"))
    nrow(dat) > 0 && ncol(dat) > 0
  }, error = function(e) FALSE)

  if (!has_data) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  obj
}

add_plot_cluster_column <- function(obj, assign_df, plot_col = PLOT_CLUSTER_COL) {
  best_match <- match_cluster_vector(colnames(obj), assign_df)
  if (is.na(best_match$rate) || best_match$rate < 0.80) {
    stop(sprintf(
      "Cluster assignment match rate too low (%.3f) for plotting column %s",
      best_match$rate,
      plot_col
    ))
  }

  cluster_levels <- sort(unique(assign_df$cluster_id))
  obj[[plot_col]] <- factor(
    as.character(best_match$cluster_vec),
    levels = as.character(cluster_levels)
  )

  attr(obj[[plot_col]], "match_mode") <- best_match$mode
  attr(obj[[plot_col]], "match_rate") <- best_match$rate
  obj
}

pick_umap_reduction <- function(obj, lineage_tag) {
  available <- tryCatch(Reductions(obj), error = function(e) character())
  defaults <- tc_default_final_schema_lineage(lineage_tag)
  reduction <- tc_pick_existing_reduction_name(
    available_reductions = available,
    candidates = unique(c("umap", defaults$default_umap_reduction_candidates)),
    fallback_regex = "^umap"
  )
  if (is.na(reduction) || !nzchar(reduction)) {
    stop(sprintf(
      "No UMAP reduction found for %s. Available reductions: %s",
      lineage_tag,
      paste(available, collapse = ", ")
    ))
  }
  reduction
}

point_size_for_obj <- function(obj) {
  n_cells <- ncol(obj)
  if (n_cells > 150000) return(0.06)
  if (n_cells > 70000) return(0.12)
  if (n_cells > 30000) return(0.18)
  0.25
}

resolve_marker_panels <- function(lineage_tag) {
  switch(
    toupper(lineage_tag),
    BCELL = tc_default_bcell_marker_panels(),
    TNK = tc_default_tnk_marker_panels(),
    EPITHELIAL = tc_default_epithelial_marker_panels(),
    EPI = tc_default_epithelial_marker_panels(),
    MYELOID = tc_default_myeloid_marker_panels(),
    STROMAL_ENDOTHELIAL = tc_default_stromal_endothelial_marker_panels(),
    STROMAL_FIBROBLAST = tc_default_stromal_fibroblast_marker_panels(),
    STROMAL_SMC = tc_default_stromal_smc_marker_panels(),
    list()
  )
}

flatten_lineage_panel_genes <- function(obj, lineage_tag, n_per_panel = LINEAGE_GENES_PER_PANEL) {
  panels <- resolve_marker_panels(lineage_tag)
  normalized <- tc_normalize_marker_panels(panels, available_features = rownames(obj))
  genes <- unlist(lapply(normalized, function(panel) head(panel$genes, n_per_panel)), use.names = FALSE)
  unique(genes)
}

read_outlier_screen <- function(spec) {
  dt <- fread(spec$screen_path, sep = "\t")
  if (nrow(dt) == 0L) return(data.table())
  if (!"record_label" %in% colnames(dt)) {
    stop(sprintf("record_label column not found in %s", spec$screen_path))
  }

  dt[, cluster_id := suppressWarnings(as.integer(str_match(as.character(record_label), "cluster_(\\d+)_vs_rest")[, 2]))]
  dt[, outlier_flag_norm := normalize_yes(outlier_flag)]
  dt[, biological_signal_class_norm := tolower(trimws(as.character(biological_signal_class)))]

  out <- dt[!is.na(cluster_id) & (outlier_flag_norm | biological_signal_class_norm == "likely_outlier")]
  if (nrow(out) == 0L) return(data.table())

  out[, `:=`(
    screen_key = spec$key,
    screen_label = spec$label,
    annotation_label = sanitize_text(annotation_label),
    short_call = sanitize_text(short_call),
    evidence_summary = sanitize_text(evidence_summary),
    confidence = sanitize_text(confidence),
    followup = sanitize_text(followup)
  )]

  out <- unique(out, by = c("screen_key", "cluster_id"))
  setorder(out, cluster_id)
  out
}

read_marker_table <- function(marker_path) {
  if (!file.exists(marker_path)) return(data.frame())
  tryCatch(fread(marker_path), error = function(e) data.frame())
}

select_top_marker_genes <- function(marker_path,
                                    available_features,
                                    n = TOP_MARKERS_PER_CLUSTER) {
  markers <- read_marker_table(marker_path)
  if (!is.data.frame(markers) || nrow(markers) == 0L) return(character())

  gene_col <- first_existing_col(markers, c("gene", "Gene", "feature"))
  lfc_col <- first_existing_col(markers, c("avg_log2FC", "avg_logFC", "log2FoldChange"))
  padj_col <- first_existing_col(markers, c("p_val_adj", "padj", "p_adj"))

  if (is.null(gene_col) || is.null(lfc_col)) return(character())

  genes <- sanitize_text(markers[[gene_col]])
  lfc <- suppressWarnings(as.numeric(markers[[lfc_col]]))
  padj <- if (!is.null(padj_col)) suppressWarnings(as.numeric(markers[[padj_col]])) else rep(NA_real_, length(genes))

  build_candidate_genes <- function(use_padj = TRUE, exclude_technical = TRUE) {
    keep <- !is.na(lfc) & lfc > 0 & nzchar(genes)
    if (use_padj) {
      keep <- keep & (is.na(padj) | padj <= 0.05)
    }
    candidate_genes <- genes[keep]
    candidate_lfc <- lfc[keep]
    ord <- order(candidate_lfc, decreasing = TRUE, na.last = TRUE)
    candidate_genes <- candidate_genes[ord]
    if (exclude_technical) {
      candidate_genes <- candidate_genes[!grepl(TECHNICAL_GENE_REGEX, candidate_genes, ignore.case = TRUE)]
    }
    tc_match_features_to_available(candidate_genes, available_features)
  }

  candidates <- build_candidate_genes(use_padj = TRUE, exclude_technical = TRUE)
  if (length(candidates) < min(3L, n)) {
    candidates <- build_candidate_genes(use_padj = FALSE, exclude_technical = TRUE)
  }
  if (length(candidates) < min(3L, n)) {
    candidates <- build_candidate_genes(use_padj = FALSE, exclude_technical = FALSE)
  }

  unique(head(candidates, n))
}

select_ranked_candidate_genes <- function(marker_path,
                                          candidate_genes,
                                          available_features,
                                          n = CONTAMINATION_FEATUREPLOT_GENES_PER_CLUSTER,
                                          fallback_genes = character()) {
  candidate_genes <- unique(tc_match_features_to_available(candidate_genes, available_features))
  fallback_genes <- unique(tc_match_features_to_available(fallback_genes, available_features))
  if (length(candidate_genes) == 0L && length(fallback_genes) == 0L) return(character())

  markers <- read_marker_table(marker_path)
  ranked <- character()
  if (is.data.frame(markers) && nrow(markers) > 0L) {
    gene_col <- first_existing_col(markers, c("gene", "Gene", "feature"))
    lfc_col <- first_existing_col(markers, c("avg_log2FC", "avg_logFC", "log2FoldChange"))
    if (!is.null(gene_col) && !is.null(lfc_col)) {
      genes <- sanitize_text(markers[[gene_col]])
      lfc <- suppressWarnings(as.numeric(markers[[lfc_col]]))
      df <- data.frame(gene = genes, lfc = lfc, stringsAsFactors = FALSE)
      df <- df[df$gene %in% candidate_genes, , drop = FALSE]
      df <- df[!duplicated(df$gene), , drop = FALSE]
      if (nrow(df) > 0L) {
        df <- df[order(df$lfc, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
        ranked <- tc_match_features_to_available(df$gene, available_features)
      }
    }
  }

  out <- unique(c(ranked, candidate_genes, fallback_genes))
  head(out, n)
}

resolve_contamination_focus_genes <- function(target,
                                              spec_key,
                                              cluster_id,
                                              marker_path,
                                              available_features,
                                              top_markers = character()) {
  explicit_genes <- get_focus_genes_from_config(target, spec_key, cluster_id)
  default_genes <- default_contamination_panel_genes(target$lineage_tag)
  candidates <- unique(c(explicit_genes, default_genes))
  select_ranked_candidate_genes(
    marker_path = marker_path,
    candidate_genes = candidates,
    available_features = available_features,
    n = CONTAMINATION_FEATUREPLOT_GENES_PER_CLUSTER,
    fallback_genes = top_markers
  )
}

save_plot_png <- function(plot_obj, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
  invisible(path)
}

draw_outlier_overview_umap <- function(obj,
                                       cluster_ids,
                                       reduction,
                                       plot_col,
                                       output_path,
                                       title) {
  cluster_ids <- as.character(sort(unique(cluster_ids)))
  meta_cluster <- as.character(obj@meta.data[[plot_col]])
  overview_levels <- c("background", paste0("c", cluster_ids))
  obj$llm_outlier_overview <- factor(
    ifelse(meta_cluster %in% cluster_ids, paste0("c", meta_cluster), "background"),
    levels = overview_levels
  )

  colors <- c(
    background = "#d1d5db",
    stats::setNames(grDevices::hcl.colors(length(cluster_ids), palette = "Dark 3"), paste0("c", cluster_ids))
  )

  p <- DimPlot(
    obj,
    reduction = reduction,
    group.by = "llm_outlier_overview",
    cols = colors,
    pt.size = point_size_for_obj(obj),
    label = TRUE,
    repel = TRUE,
    shuffle = FALSE
  ) +
    ggtitle(title) +
    coord_equal() +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  save_plot_png(p, output_path, width = 10, height = 8)
}

draw_cluster_highlight_umap <- function(obj,
                                        cluster_id,
                                        reduction,
                                        plot_col,
                                        output_path,
                                        title,
                                        subtitle = "") {
  cluster_id <- as.character(cluster_id)
  meta_cluster <- as.character(obj@meta.data[[plot_col]])
  obj$llm_cluster_highlight <- factor(
    ifelse(meta_cluster == cluster_id, paste0("c", cluster_id), "other"),
    levels = c("other", paste0("c", cluster_id))
  )

  p <- DimPlot(
    obj,
    reduction = reduction,
    group.by = "llm_cluster_highlight",
    cols = c(other = "#e5e7eb", stats::setNames("#dc2626", paste0("c", cluster_id))),
    pt.size = point_size_for_obj(obj),
    shuffle = FALSE
  ) +
    labs(title = title, subtitle = subtitle) +
    coord_equal() +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "none"
    )

  save_plot_png(p, output_path, width = 8.5, height = 7)
}

draw_featureplot_grid <- function(obj,
                                  genes,
                                  reduction,
                                  output_path,
                                  title,
                                  subtitle = "") {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(obj)]
  if (length(genes) == 0L) return(invisible(NULL))

  ncol_plot <- min(2L, length(genes))
  nrow_plot <- ceiling(length(genes) / ncol_plot)

  p <- FeaturePlot(
    obj,
    features = genes,
    reduction = reduction,
    combine = TRUE,
    ncol = ncol_plot,
    order = TRUE,
    cols = c("#f3f4f6", "#b91c1c")
  ) &
    theme_classic(base_size = 11) &
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold")
    )

  p <- p + patchwork::plot_annotation(title = title, subtitle = subtitle)
  save_plot_png(p, output_path, width = max(8, 4.5 * ncol_plot), height = 3.9 * nrow_plot + 1.2)
}

build_plot_labels <- function(cluster_ids, annotation_labels) {
  lbl <- ifelse(
    nzchar(annotation_labels),
    sprintf("c%s\n%s", cluster_ids, str_trunc(annotation_labels, 28)),
    sprintf("c%s", cluster_ids)
  )
  stats::setNames(lbl, as.character(cluster_ids))
}

draw_outlier_dotplot <- function(obj,
                                 cluster_ids,
                                 plot_labels,
                                 genes,
                                 plot_col,
                                 output_path,
                                 title,
                                 subtitle = "") {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(obj)]
  cluster_ids <- as.character(sort(unique(cluster_ids)))
  if (length(cluster_ids) == 0L || length(genes) < 2L) return(invisible(NULL))

  keep_cells <- colnames(obj)[as.character(obj@meta.data[[plot_col]]) %in% cluster_ids]
  if (length(keep_cells) == 0L) return(invisible(NULL))

  obj_sub <- subset(obj, cells = keep_cells)
  obj_sub$llm_dotplot_group <- factor(
    unname(plot_labels[as.character(obj_sub@meta.data[[plot_col]])]),
    levels = unname(plot_labels[cluster_ids])
  )

  p <- DotPlot(
    object = obj_sub,
    features = genes,
    group.by = "llm_dotplot_group",
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
      size = "Pct.Exp"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )

  width <- max(10, 0.4 * length(genes) + 4)
  height <- max(6, 0.45 * length(cluster_ids) + 2.5)
  save_plot_png(p, output_path, width = width, height = height)
}

rel_from_review_dir <- function(review_dir, path) {
  review_dir <- normalizePath(review_dir, winslash = "/", mustWork = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  output_dir <- normalizePath(dirname(dirname(review_dir)), winslash = "/", mustWork = FALSE)

  if (startsWith(path, paste0(review_dir, "/"))) {
    return(sub(paste0("^", review_dir, "/"), "", path))
  }
  if (startsWith(path, paste0(output_dir, "/"))) {
    return(file.path("..", "..", sub(paste0("^", output_dir, "/"), "", path)))
  }
  basename(path)
}

build_cluster_section <- function(row) {
  lines <- c(
    sprintf("### c%s · %s", row$cluster_id, row$annotation_label),
    "",
    sprintf("- confidence: `%s`", row$confidence),
    sprintf("- class: `%s` | outlier_flag: `%s`", row$biological_signal_class, row$outlier_flag),
    if (row$user_remove_request == "yes") "- user removal request: **yes**" else NULL,
    if (nzchar(row$short_call)) sprintf("- short call: %s", row$short_call) else NULL,
    if (nzchar(row$evidence_summary)) sprintf("- evidence: %s", row$evidence_summary) else NULL,
    if (nzchar(row$top_marker_genes)) sprintf("- top markers used: `%s`", row$top_marker_genes) else NULL,
    if (nzchar(row$contamination_focus_genes)) sprintf("- contamination / wrong-lineage genes: `%s`", row$contamination_focus_genes) else NULL,
    ""
  )

  if (!is.na(row$highlight_umap_relpath) && nzchar(row$highlight_umap_relpath)) {
    lines <- c(lines, sprintf("![c%s highlight](%s)", row$cluster_id, row$highlight_umap_relpath), "")
  }
  if (!is.na(row$featureplot_relpath) && nzchar(row$featureplot_relpath)) {
    lines <- c(lines, sprintf("![c%s featureplot](%s)", row$cluster_id, row$featureplot_relpath), "")
  }
  if (!is.na(row$contamination_featureplot_relpath) && nzchar(row$contamination_featureplot_relpath)) {
    lines <- c(lines, sprintf("![c%s contamination featureplot](%s)", row$cluster_id, row$contamination_featureplot_relpath), "")
  }
  lines
}

build_screen_section <- function(spec,
                                 review_dir,
                                 summary_path,
                                 summary_dt,
                                 overview_path = NA_character_,
                                 top_dotplot_path = NA_character_,
                                 lineage_dotplot_path = NA_character_,
                                 contamination_dotplot_path = NA_character_) {
  if (nrow(summary_dt) == 0L) {
    return(c(
      sprintf("## %s outliers", spec$label),
      "",
      "- No clusters were flagged as outliers in the current screen.",
      ""
    ))
  }

  source_md_lines <- vapply(
    spec$source_md,
    function(path) sprintf("- [%s](%s)", basename(path), rel_from_review_dir(review_dir, path)),
    character(1)
  )
  source_screen_rel <- rel_from_review_dir(review_dir, spec$screen_path)
  requested_tbl <- summary_dt[summary_dt$user_remove_request == "yes", , drop = FALSE]

  lines <- c(
    sprintf("## %s outliers (%d)", spec$label, nrow(summary_dt)),
    "",
    sprintf("- Source screen: [`%s`](%s)", basename(spec$screen_path), source_screen_rel),
    if (length(source_md_lines) > 0L) "- Source LLM markdowns:" else NULL,
    if (length(source_md_lines) > 0L) source_md_lines else NULL,
    sprintf("- Summary table: [`%s`](%s)", basename(summary_path), rel_from_review_dir(review_dir, summary_path)),
    if (nrow(requested_tbl) > 0L) sprintf("- User-requested removal clusters: %s", paste(sprintf("c%s", requested_tbl$cluster_id), collapse = ", ")) else NULL,
    ""
  )

  if (!is.na(overview_path) && nzchar(overview_path)) {
    lines <- c(lines, sprintf("![%s outlier overview](%s)", spec$label, rel_from_review_dir(review_dir, overview_path)), "")
  }
  if (!is.na(top_dotplot_path) && nzchar(top_dotplot_path)) {
    lines <- c(lines, sprintf("![%s outlier top-marker dotplot](%s)", spec$label, rel_from_review_dir(review_dir, top_dotplot_path)), "")
  }
  if (!is.na(lineage_dotplot_path) && nzchar(lineage_dotplot_path)) {
    lines <- c(lines, sprintf("![%s outlier lineage-panel dotplot](%s)", spec$label, rel_from_review_dir(review_dir, lineage_dotplot_path)), "")
  }
  if (!is.na(contamination_dotplot_path) && nzchar(contamination_dotplot_path)) {
    lines <- c(lines, sprintf("![%s contamination dotplot](%s)", spec$label, rel_from_review_dir(review_dir, contamination_dotplot_path)), "")
  }

  for (i in seq_len(nrow(summary_dt))) {
    lines <- c(lines, build_cluster_section(summary_dt[i, ]))
  }

  lines
}

process_target <- function(target) {
  output_dir <- target$output_dir
  lineage_tag <- target$lineage_tag
  object_path <- if (!is.null(target$object_path)) target$object_path else find_final_object_path(output_dir)
  preferred_kind <- if (!is.null(target$cluster_kind)) target$cluster_kind else NULL

  cat(sprintf("\n=== %s | %s ===\n", lineage_tag, output_dir))
  review_dir <- file.path(output_dir, "figures", OUTPUT_STEM)
  tables_dir <- file.path(review_dir, "tables")
  umap_dir <- file.path(review_dir, "umap")
  feature_dir <- file.path(review_dir, "featureplots")
  dotplot_dir <- file.path(review_dir, "dotplots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dotplot_dir, recursive = TRUE, showWarnings = FALSE)

  screen_specs <- build_screen_specs(output_dir, preferred_kind = preferred_kind)
  if (length(screen_specs) == 0L) {
    warning(sprintf("No supported outlier screens found under: %s", output_dir))
    return(data.table(
      lineage_tag = lineage_tag,
      output_dir = output_dir,
      review_dir = review_dir,
      status = "skipped",
      n_screens = 0L,
      n_outlier_clusters = 0L,
      n_requested_removal_clusters = 0L,
      note = "no outlier discovery screens found"
    ))
  }

  primary_kind <- screen_specs[[1]]$cluster_kind
  obj <- readRDS(object_path)
  obj <- obj[, sort(colnames(obj))]
  obj <- ensure_expression_ready(obj)
  assignments <- read_cluster_assignments(output_dir, primary_kind)
  obj <- add_plot_cluster_column(obj, assignments, plot_col = PLOT_CLUSTER_COL)
  umap_reduction <- pick_umap_reduction(obj, lineage_tag)
  lineage_panel_genes <- flatten_lineage_panel_genes(obj, lineage_tag)

  readme_lines <- c(
    sprintf("# %s LLM outlier visual review", lineage_tag),
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    sprintf("- Output directory: `%s`", output_dir),
    sprintf("- Final object: `%s`", basename(object_path)),
    sprintf("- Plotting reduction: `%s`", umap_reduction),
    sprintf("- Plot cluster assignment: `%s`", primary_kind),
    ""
  )

  combined_rows <- list()
  screen_summary_rows <- list()

  for (spec in screen_specs) {
    cat(sprintf("[screen] %s\n", spec$label))
    screen_dt <- read_outlier_screen(spec)
    summary_path <- file.path(tables_dir, sprintf("%s_outlier_summary.tsv", spec$key))

    if (nrow(screen_dt) == 0L) {
      fwrite(data.table(), summary_path, sep = "\t")
      readme_lines <- c(readme_lines, build_screen_section(spec, review_dir, summary_path, data.frame()))
      screen_summary_rows[[length(screen_summary_rows) + 1L]] <- data.table(
        screen_key = spec$key,
        screen_label = spec$label,
        n_outlier_clusters = 0L,
        n_requested_removal_clusters = 0L
      )
      next
    }

    cluster_ids <- as.character(sort(unique(screen_dt$cluster_id)))
    plot_labels <- build_plot_labels(
      cluster_ids = screen_dt$cluster_id,
      annotation_labels = screen_dt$annotation_label
    )

    overview_path <- file.path(umap_dir, sprintf("%s_outlier_overview.png", spec$key))
    draw_outlier_overview_umap(
      obj = obj,
      cluster_ids = cluster_ids,
      reduction = umap_reduction,
      plot_col = PLOT_CLUSTER_COL,
      output_path = overview_path,
      title = sprintf("%s outlier overview", spec$label)
    )

    cluster_rows <- list()
    top_gene_union <- character()
    contamination_gene_union <- character()

    for (i in seq_len(nrow(screen_dt))) {
      row <- screen_dt[i]
      cluster_id <- as.integer(row$cluster_id)
      marker_path <- file.path(spec$cluster_dir, sprintf("ofa_%s", cluster_id), "markers.csv")
      top_markers <- select_top_marker_genes(
        marker_path = marker_path,
        available_features = rownames(obj),
        n = TOP_MARKERS_PER_CLUSTER
      )
      feature_genes <- head(top_markers, FEATUREPLOT_GENES_PER_CLUSTER)
      dotplot_genes <- head(top_markers, DOTPLOT_GENES_PER_CLUSTER)
      top_gene_union <- unique(c(top_gene_union, dotplot_genes))

      contamination_genes <- resolve_contamination_focus_genes(
        target = target,
        spec_key = spec$key,
        cluster_id = cluster_id,
        marker_path = marker_path,
        available_features = rownames(obj),
        top_markers = top_markers
      )
      contamination_gene_union <- unique(c(contamination_gene_union, contamination_genes))

      requested_remove_ids <- get_requested_removal_clusters(target, spec$key)
      user_remove_request <- if (as.character(cluster_id) %in% requested_remove_ids) "yes" else "no"

      cluster_stub <- sprintf("%s_cluster_%02d", spec$key, cluster_id)
      highlight_path <- file.path(umap_dir, sprintf("%s_highlight.png", cluster_stub))
      featureplot_path <- file.path(feature_dir, sprintf("%s_top_markers.png", cluster_stub))
      contamination_featureplot_path <- file.path(feature_dir, sprintf("%s_contamination_markers.png", cluster_stub))

      draw_cluster_highlight_umap(
        obj = obj,
        cluster_id = cluster_id,
        reduction = umap_reduction,
        plot_col = PLOT_CLUSTER_COL,
        output_path = highlight_path,
        title = sprintf("%s outlier c%s", spec$label, cluster_id),
        subtitle = str_trunc(row$short_call, 120)
      )

      if (length(feature_genes) > 0L) {
        draw_featureplot_grid(
          obj = obj,
          genes = feature_genes,
          reduction = umap_reduction,
          output_path = featureplot_path,
          title = sprintf("%s c%s top-marker FeaturePlot", spec$label, cluster_id),
          subtitle = paste(feature_genes, collapse = ", ")
        )
      }

      if (length(contamination_genes) > 0L) {
        draw_featureplot_grid(
          obj = obj,
          genes = contamination_genes,
          reduction = umap_reduction,
          output_path = contamination_featureplot_path,
          title = sprintf("%s c%s contamination / wrong-lineage FeaturePlot", spec$label, cluster_id),
          subtitle = paste(contamination_genes, collapse = ", ")
        )
      }

      cluster_rows[[length(cluster_rows) + 1L]] <- data.table(
        screen_key = spec$key,
        screen_label = spec$label,
        cluster_id = cluster_id,
        record_id = sanitize_text(row$record_id),
        record_label = sanitize_text(row$record_label),
        annotation_label = sanitize_text(row$annotation_label),
        biological_signal_class = sanitize_text(row$biological_signal_class),
        confidence = sanitize_text(row$confidence),
        short_call = sanitize_text(row$short_call),
        outlier_flag = sanitize_text(row$outlier_flag),
        evidence_summary = sanitize_text(row$evidence_summary),
        followup = sanitize_text(row$followup),
        marker_path = marker_path,
        top_marker_genes = paste(top_markers, collapse = ", "),
        contamination_focus_genes = paste(contamination_genes, collapse = ", "),
        user_remove_request = user_remove_request,
        highlight_umap_relpath = rel_from_review_dir(review_dir, highlight_path),
        featureplot_relpath = if (file.exists(featureplot_path)) rel_from_review_dir(review_dir, featureplot_path) else "",
        contamination_featureplot_relpath = if (file.exists(contamination_featureplot_path)) rel_from_review_dir(review_dir, contamination_featureplot_path) else ""
      )
    }

    summary_dt <- rbindlist(cluster_rows, fill = TRUE)
    fwrite(summary_dt, summary_path, sep = "\t")
    combined_rows[[length(combined_rows) + 1L]] <- summary_dt

    top_dotplot_path <- file.path(dotplot_dir, sprintf("%s_outlier_dotplot_top_markers.png", spec$key))
    draw_outlier_dotplot(
      obj = obj,
      cluster_ids = cluster_ids,
      plot_labels = plot_labels,
      genes = unique(top_gene_union),
      plot_col = PLOT_CLUSTER_COL,
      output_path = top_dotplot_path,
      title = sprintf("%s outlier dotplot: top markers", spec$label),
      subtitle = "Top positive DE genes per outlier cluster"
    )

    lineage_dotplot_path <- file.path(dotplot_dir, sprintf("%s_outlier_dotplot_lineage_panels.png", spec$key))
    draw_outlier_dotplot(
      obj = obj,
      cluster_ids = cluster_ids,
      plot_labels = plot_labels,
      genes = lineage_panel_genes,
      plot_col = PLOT_CLUSTER_COL,
      output_path = lineage_dotplot_path,
      title = sprintf("%s outlier dotplot: lineage panels", spec$label),
      subtitle = "Default lineage marker panels for mismatch / contamination review"
    )

    contamination_dotplot_path <- file.path(dotplot_dir, sprintf("%s_outlier_dotplot_contamination_markers.png", spec$key))
    draw_outlier_dotplot(
      obj = obj,
      cluster_ids = cluster_ids,
      plot_labels = plot_labels,
      genes = unique(head(contamination_gene_union, CONTAMINATION_DOTPLOT_MAX_GENES)),
      plot_col = PLOT_CLUSTER_COL,
      output_path = contamination_dotplot_path,
      title = sprintf("%s outlier dotplot: contamination / wrong-lineage markers", spec$label),
      subtitle = "Contamination genes or wrong-lineage markers highlighted in the outlier review"
    )

    readme_lines <- c(
      readme_lines,
      build_screen_section(
        spec = spec,
        review_dir = review_dir,
        summary_path = summary_path,
        summary_dt = summary_dt,
        overview_path = if (file.exists(overview_path)) overview_path else NA_character_,
        top_dotplot_path = if (file.exists(top_dotplot_path)) top_dotplot_path else NA_character_,
        lineage_dotplot_path = if (file.exists(lineage_dotplot_path)) lineage_dotplot_path else NA_character_,
        contamination_dotplot_path = if (file.exists(contamination_dotplot_path)) contamination_dotplot_path else NA_character_
      )
    )

    screen_summary_rows[[length(screen_summary_rows) + 1L]] <- data.table(
      screen_key = spec$key,
      screen_label = spec$label,
      n_outlier_clusters = nrow(summary_dt),
      n_requested_removal_clusters = sum(summary_dt$user_remove_request == "yes")
    )
  }

  combined_dt <- if (length(combined_rows) > 0L) rbindlist(combined_rows, fill = TRUE) else data.table()
  combined_path <- file.path(tables_dir, "all_outlier_clusters.tsv")
  fwrite(combined_dt, combined_path, sep = "\t")

  screen_summary_dt <- if (length(screen_summary_rows) > 0L) rbindlist(screen_summary_rows, fill = TRUE) else data.table()
  screen_summary_path <- file.path(tables_dir, "screen_outlier_counts.tsv")
  fwrite(screen_summary_dt, screen_summary_path, sep = "\t")

  readme_lines <- c(
    readme_lines,
    "## Files",
    "",
    sprintf("- [`%s`](%s): combined outlier rows across all processed screens", basename(combined_path), rel_from_review_dir(review_dir, combined_path)),
    sprintf("- [`%s`](%s): outlier cluster counts per screen", basename(screen_summary_path), rel_from_review_dir(review_dir, screen_summary_path)),
    ""
  )

  writeLines(readme_lines, file.path(review_dir, "README.md"))

  data.table(
    lineage_tag = lineage_tag,
    output_dir = output_dir,
    review_dir = review_dir,
    status = "ok",
    n_screens = nrow(screen_summary_dt),
    n_outlier_clusters = nrow(combined_dt),
    n_requested_removal_clusters = sum(combined_dt$user_remove_request == "yes"),
    note = sprintf("umap=%s; plot_cluster=%s", umap_reduction, primary_kind)
  )
}

run_outlier_visual_batch <- function(targets = build_default_targets()) {
  cat("=== Tissue Comparison LLM Outlier Visual Batch (2026-05-07) ===\n")
  batch_rows <- list()

  for (target in targets) {
    res <- tryCatch(
      process_target(target),
      error = function(e) {
        data.table(
          lineage_tag = target$lineage_tag,
          output_dir = target$output_dir,
          review_dir = file.path(target$output_dir, "figures", OUTPUT_STEM),
          status = "error",
          n_screens = NA_integer_,
          n_outlier_clusters = NA_integer_,
          n_requested_removal_clusters = NA_integer_,
          note = conditionMessage(e)
        )
      }
    )
    batch_rows[[length(batch_rows) + 1L]] <- res
    cat(sprintf(
      "[result] %s | %s | outliers=%s | requested_remove=%s\n",
      res$lineage_tag,
      res$status,
      res$n_outlier_clusters,
      res$n_requested_removal_clusters
    ))
    gc(verbose = FALSE)
  }

  batch_summary <- rbindlist(batch_rows, fill = TRUE)
  print(batch_summary)
  batch_summary
}

if (sys.nframe() == 0) {
  batch_summary <- run_outlier_visual_batch()
  if (any(batch_summary$status == "error")) {
    stop("One or more targets failed during LLM outlier visualization batch generation")
  }
  cat("[DONE] All target outlier review figures generated successfully.\n")
}
