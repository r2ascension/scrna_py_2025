#!/usr/bin/env Rscript
# ==============================================================================
# BayesPrism Downstream Analysis Pipeline
# ==============================================================================
# Version:  v1.1 (2026-03-08)
# Status:   Production-ready
# Changes from v1.0:
#   [P0] Z/theta cell type alignment by name+order (not just count)
#   [P0] Enrichment: specificity-based top genes (log2 ct vs others mean)
#   [P0] PCA/correlation: remove zero-variance cell types before computation
#   [P0] Dotplot: floor p.adjust to .Machine$double.xmin before log10 scale
#   [P0] LLM parser: line-by-line state machine (no fragile regex on multi-line)
#   [P1] theta_cv: validate and reorder sample rows to match theta
#   [P1] Stacked bar: implement real sort by dominant cell type (remove dead code)
#   [P1] Remove Z_var dead code
#   [P1] GMT loader: toupper() genes for consistent case
#   [P1] DeepSeek API: encode="json" (cleaner than raw+toJSON)
#   [P1] LLM condition: nzchar() instead of != ""
#   [P1] build_celltype_prompt: sort by specificity column
#
# Input  (from bayesprism_deconvolution_v3.0.R):
#   - theta_final.csv   : cell type fractions (samples x cell types)
#   - theta_cv.csv      : coefficient of variation (uncertainty)
#   - Z_matrix.rds      : deconvolved expression (genes x cell types)
#
# Output:
#   - Extended visualizations (stacked bar, jitter, PCA, correlation)
#   - GMT-based functional enrichment per cell type
#   - DeepSeek LLM interpretation of cell type biology
#   - DOWNSTREAM_REPORT.md
#
# Notes:
#   - n=5 samples -> descriptive/exploratory only; no formal statistical tests
#   - All labels/comments in English (publication standard)
# ==============================================================================

# ==============================================================================
# 0. Configuration  (MODIFY THESE)
# ==============================================================================

# ----- Input paths -----
BAYESPRISM_OUTPUT_DIR <- "/home/h2048/data/R/bulk0303/bayesprism_out"
THETA_CSV <- file.path(BAYESPRISM_OUTPUT_DIR, "theta_final.csv")
THETA_CV_CSV <- file.path(BAYESPRISM_OUTPUT_DIR, "theta_cv.csv")
Z_RDS <- file.path(BAYESPRISM_OUTPUT_DIR, "Z_matrix.rds")

# ----- Output -----
OUTPUT_DIR <- file.path(BAYESPRISM_OUTPUT_DIR, "downstream_v1_1")

# ----- Reference databases (MSigDB GMT) -----
GMT_HALLMARK <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"

# ----- Analysis parameters -----
N_TOP_GENES_ENRICHMENT <- 100 # Top specific genes per cell type for enrichment
N_TOP_TERMS_PLOT <- 15 # Top enriched terms to display per cell type
MIN_GENE_SET_SIZE <- 10
MAX_GENE_SET_SIZE <- 500
ENRICH_PVAL_CUTOFF <- 0.05
ENRICH_QVAL_CUTOFF <- 0.20

# ----- LLM (DeepSeek API) -----
RUN_LLM_INTERPRETATION <- TRUE
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
DEEPSEEK_MODEL <- "deepseek-chat"
LLM_MAX_TOKENS <- 800
LLM_TEMPERATURE <- 0.3
# Study context injected into LLM prompt
STUDY_CONTEXT <- paste(
  "BayesPrism deconvolution of bulk RNA-seq from Eustachian tube epithelium.",
  "Reference: LungMap single-cell atlas. Disease context: upper airway disease",
  "including chronic rhinosinusitis (CRSwNP) and type 2 inflammation.",
  "Focus on cell type functional roles in airway biology and disease pathogenesis."
)

# ----- Visualization -----
FIGURE_WIDTH_BASE <- 10
FIGURE_HEIGHT_BASE <- 7

# ----- Reproducibility -----
set.seed(42)

# ==============================================================================
# Thread limiting (prevent BLAS over-subscription)
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS = "4",
  MKL_NUM_THREADS = "4",
  OPENBLAS_NUM_THREADS = "4",
  NUMEXPR_NUM_THREADS = "4"
)

# ==============================================================================
# 1. Libraries
# ==============================================================================

cat("\n=== Loading Libraries ===\n")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(corrplot)
  library(clusterProfiler)
  library(dplyr)
  library(tidyr)
  library(httr)
  library(jsonlite)
})

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)

# Null-coalescing operator
`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) == 1 && !is.na(a)) a else b
}

cat("Output directory:", OUTPUT_DIR, "\n")

# ==============================================================================
# 2. Load BayesPrism Results  [v1.1: full name+order alignment]
# ==============================================================================

cat("\n=== Step 1/7: Loading BayesPrism Outputs ===\n")

# ---- theta: samples x cell types ----
theta_raw <- fread(THETA_CSV, header = TRUE)
stopifnot("theta_final.csv must have >= 2 columns" = ncol(theta_raw) >= 2)

sample_names <- theta_raw[[1]]
theta <- data.matrix(theta_raw[, -1, with = FALSE])
rownames(theta) <- sample_names
cat(sprintf(
  "  theta:    %d samples x %d cell types\n",
  nrow(theta),
  ncol(theta)
))

# ---- theta_cv: uncertainty ----
# [v1.1 P1] Validate both sample rows and cell type columns, reorder to match theta
theta_cv_raw <- fread(THETA_CV_CSV, header = TRUE)
stopifnot("theta_cv.csv must have >= 2 columns" = ncol(theta_cv_raw) >= 2)

theta_cv <- data.matrix(theta_cv_raw[, -1, with = FALSE])
rownames(theta_cv) <- theta_cv_raw[[1]]
cat(sprintf(
  "  theta_cv: %d samples x %d cell types\n",
  nrow(theta_cv),
  ncol(theta_cv)
))

if (!setequal(rownames(theta), rownames(theta_cv))) {
  stop("theta and theta_cv sample names do not match")
}
if (!setequal(colnames(theta), colnames(theta_cv))) {
  stop("theta and theta_cv cell type names do not match")
}
theta_cv <- theta_cv[rownames(theta), colnames(theta), drop = FALSE]

# ---- Z matrix: genes x cell types ----
# [v1.1 P0] Alignment by name+order, not just column count
Z <- readRDS(Z_RDS)
if (!is.matrix(Z)) {
  Z <- as.matrix(Z)
}

if (!is.null(colnames(Z)) && setequal(colnames(theta), colnames(Z))) {
  # Correct orientation, reorder columns to match theta
  Z <- Z[, colnames(theta), drop = FALSE]
} else if (!is.null(rownames(Z)) && setequal(colnames(theta), rownames(Z))) {
  cat("  [INFO] Z matched theta cell types by rownames -> transposing\n")
  Z <- t(Z)
  Z <- Z[, colnames(theta), drop = FALSE]
} else if (nrow(Z) == ncol(theta) && ncol(Z) > nrow(Z)) {
  # Heuristic: more columns than rows -> likely celltypes x genes
  cat("  [INFO] Z dimensions suggest cell types x genes -> transposing\n")
  Z <- t(Z)
} else if (ncol(Z) != ncol(theta)) {
  stop(sprintf(
    "Z/theta mismatch: theta has %d cell types, Z has %d columns. Check input files.",
    ncol(theta),
    ncol(Z)
  ))
}

# Final name validation after orientation fix
if (!is.null(colnames(Z)) && !setequal(colnames(theta), colnames(Z))) {
  stop("theta and Z cell type names do not match after orientation correction")
}
Z <- Z[, colnames(theta), drop = FALSE] # enforce identical order

storage.mode(theta) <- "double"
storage.mode(theta_cv) <- "double"
storage.mode(Z) <- "double"

cat(sprintf("  Z matrix: %d genes x %d cell types\n", nrow(Z), ncol(Z)))

n_samples <- nrow(theta)
n_cell_types <- ncol(theta)
cell_types <- colnames(theta)
cat(sprintf("  Samples: %s\n", paste(rownames(theta), collapse = ", ")))
cat(sprintf(
  "  Cell types (%d): %s\n",
  n_cell_types,
  paste(cell_types, collapse = ", ")
))

# Color palette (color-blind tolerant)
if (n_cell_types <= 8) {
  ct_colors <- brewer.pal(max(n_cell_types, 3), "Set2")
} else if (n_cell_types <= 12) {
  ct_colors <- brewer.pal(12, "Paired")
} else {
  ct_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_cell_types)
}
names(ct_colors) <- cell_types

# [v1.1 P3] Pre-compute non-constant theta for PCA/correlation
# Remove zero-variance cell types to avoid NaN/unstable results
theta_sd <- apply(theta, 2, sd, na.rm = TRUE)
theta_nonconst <- theta[, theta_sd > 0, drop = FALSE]
n_removed_cv <- ncol(theta) - ncol(theta_nonconst)
if (n_removed_cv > 0) {
  cat(sprintf(
    "  [INFO] %d zero-variance cell type(s) excluded from PCA/correlation: %s\n",
    n_removed_cv,
    paste(setdiff(cell_types, colnames(theta_nonconst)), collapse = ", ")
  ))
}

# ==============================================================================
# 3. Basic Composition Visualizations
# ==============================================================================

cat("\n=== Step 2/7: Composition Visualizations ===\n")

theta_dt <- as.data.table(theta, keep.rownames = "sample")
theta_long <- melt(
  theta_dt,
  id.vars = "sample",
  variable.name = "cell_type",
  value.name = "fraction"
)
theta_long[, cell_type := factor(cell_type, levels = cell_types)]

# [v1.1 P1] Stacked bar: real sort by dominant cell type, then dominant fraction
sample_order_dt <- theta_long[, .SD[which.max(fraction)], by = sample]
sample_order <- sample_order_dt[order(cell_type, -fraction), sample]
theta_long[, sample := factor(sample, levels = sample_order)]

# ---- 3a. Stacked bar ----
tryCatch(
  {
    p_stack <- ggplot(
      theta_long,
      aes(x = sample, y = fraction, fill = cell_type)
    ) +
      geom_col(width = 0.72, color = "white", linewidth = 0.25) +
      scale_fill_manual(values = ct_colors) +
      scale_y_continuous(expand = c(0, 0), labels = scales::percent_format()) +
      labs(
        title = "BayesPrism: Cell Type Composition per Sample",
        x = "Bulk Sample",
        y = "Estimated Fraction",
        fill = "Cell Type"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right"
      )
    ggsave(
      file.path(OUTPUT_DIR, "figures", "01_stacked_bar.pdf"),
      p_stack,
      width = FIGURE_WIDTH_BASE,
      height = FIGURE_HEIGHT_BASE
    )
    cat("  Saved: 01_stacked_bar.pdf\n")
  },
  error = function(e) warning("Stacked bar: ", e$message)
)

# ---- 3b. Jitter + mean bar ----
mean_dt <- theta_long[, .(mean_frac = mean(fraction)), by = cell_type]

tryCatch(
  {
    p_jitter <- ggplot(
      theta_long,
      aes(x = cell_type, y = fraction, color = cell_type)
    ) +
      geom_col(
        data = mean_dt,
        aes(x = cell_type, y = mean_frac, fill = cell_type),
        inherit.aes = FALSE,
        width = 0.55,
        alpha = 0.35,
        color = NA
      ) +
      geom_jitter(width = 0.15, size = 3, alpha = 0.85) +
      scale_color_manual(values = ct_colors) +
      scale_fill_manual(values = ct_colors) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(
        title = sprintf(
          "Cell Type Fraction Distribution (n = %d samples)",
          n_samples
        ),
        x = "Cell Type",
        y = "Estimated Fraction",
        caption = "Bars = mean; dots = individual samples"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank()
      )
    ggsave(
      file.path(OUTPUT_DIR, "figures", "02_fraction_jitter.pdf"),
      p_jitter,
      width = FIGURE_WIDTH_BASE,
      height = FIGURE_HEIGHT_BASE
    )
    cat("  Saved: 02_fraction_jitter.pdf\n")
  },
  error = function(e) warning("Jitter plot: ", e$message)
)

# ---- 3c. Theta heatmap ----
tryCatch(
  {
    pdf(
      file.path(OUTPUT_DIR, "figures", "03_theta_heatmap.pdf"),
      width = 10,
      height = max(4, n_samples * 0.8 + 2)
    )
    pheatmap(
      theta,
      cluster_rows = n_samples > 2,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "#4393c3", "#08306b"))(100),
      border_color = "grey80",
      main = "Cell Type Fractions (theta)",
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = 9,
      fontsize_col = 10,
      fontsize_row = 10,
      angle_col = 45
    )
    dev.off()
    cat("  Saved: 03_theta_heatmap.pdf\n")
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("Theta heatmap: ", e$message)
  }
)

# ---- 3d. Uncertainty (CV) heatmap ----
tryCatch(
  {
    pdf(
      file.path(OUTPUT_DIR, "figures", "04_theta_cv_heatmap.pdf"),
      width = 10,
      height = max(4, n_samples * 0.8 + 2)
    )
    pheatmap(
      theta_cv,
      cluster_rows = n_samples > 2,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "#f4a582", "#b2182b"))(100),
      border_color = "grey80",
      main = "Estimation Uncertainty (Coefficient of Variation)",
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = 9,
      fontsize_col = 10,
      fontsize_row = 10,
      angle_col = 45
    )
    dev.off()
    cat("  Saved: 04_theta_cv_heatmap.pdf\n")
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("CV heatmap: ", e$message)
  }
)

# ---- 3e. Sample PCA (theta-space) ----
# [v1.1 P0] Use theta_nonconst; unified tryCatch with internal ggrepel fallback
if (n_samples >= 3 && ncol(theta_nonconst) >= 2) {
  tryCatch(
    {
      theta_scaled <- scale(theta_nonconst)
      pca_res <- prcomp(theta_scaled, center = FALSE, scale. = FALSE)
      var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
      pca_dt <- data.table(
        sample = rownames(theta_nonconst),
        PC1 = pca_res$x[, 1],
        PC2 = pca_res$x[, 2]
      )

      # Prefer ggrepel for label placement; fall back to geom_text
      label_layer <- tryCatch(
        ggrepel::geom_text_repel(
          aes(label = sample),
          size = 3.5,
          max.overlaps = 20
        ),
        error = function(e) {
          geom_text(aes(label = sample), vjust = -0.8, size = 3.5)
        }
      )

      p_pca <- ggplot(pca_dt, aes(x = PC1, y = PC2)) +
        geom_point(size = 5, color = "#2166ac") +
        label_layer +
        labs(
          title = "Sample PCA in Cell Type Composition Space",
          x = sprintf("PC1 (%.1f%%)", var_exp[1]),
          y = sprintf("PC2 (%.1f%%)", var_exp[2]),
          caption = sprintf(
            "PCA on %d variable cell types",
            ncol(theta_nonconst)
          )
        ) +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor = element_blank())
      ggsave(
        file.path(OUTPUT_DIR, "figures", "05_theta_PCA.pdf"),
        p_pca,
        width = 8,
        height = 7
      )
      cat("  Saved: 05_theta_PCA.pdf\n")
    },
    error = function(e) warning("PCA plot: ", e$message)
  )
} else {
  cat(sprintf(
    "  Skipping PCA: need >= 3 samples and >= 2 variable cell types (have %d and %d)\n",
    n_samples,
    ncol(theta_nonconst)
  ))
}

# ---- 3f. Cell type variability across samples ----
tryCatch(
  {
    cv_across <- apply(theta, 2, function(x) sd(x) / (mean(x) + 1e-10))
    cv_df <- data.frame(
      cell_type = names(cv_across),
      cv = cv_across,
      mean_frac = colMeans(theta)
    )
    cv_df$cell_type <- factor(
      cv_df$cell_type,
      levels = cv_df$cell_type[order(cv_df$cv)]
    )

    p_cv <- ggplot(cv_df, aes(x = cell_type, y = cv, fill = mean_frac)) +
      geom_col() +
      scale_fill_gradientn(
        colors = c("#deebf7", "#084594"),
        name = "Mean\nFraction"
      ) +
      coord_flip() +
      labs(
        title = "Cell Type Variability Across Samples",
        x = "Cell Type",
        y = "Coefficient of Variation (between samples)",
        caption = "Higher CV = more variable across samples"
      ) +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank())
    ggsave(
      file.path(OUTPUT_DIR, "figures", "06_celltype_CV.pdf"),
      p_cv,
      width = 9,
      height = max(5, n_cell_types * 0.4 + 2)
    )
    cat("  Saved: 06_celltype_CV.pdf\n")
  },
  error = function(e) warning("CV plot: ", e$message)
)

# ==============================================================================
# 4. Cell Type Correlation  [v1.1 P0: use theta_nonconst]
# ==============================================================================

cat("\n=== Step 3/7: Cell Type Correlation Analysis ===\n")

if (n_samples >= 3 && ncol(theta_nonconst) >= 2) {
  tryCatch(
    {
      cor_mat <- cor(
        theta_nonconst,
        method = "pearson",
        use = "pairwise.complete.obs"
      )

      fwrite(
        as.data.table(cor_mat, keep.rownames = "cell_type"),
        file.path(OUTPUT_DIR, "tables", "celltype_correlation_matrix.csv")
      )

      pdf(
        file.path(OUTPUT_DIR, "figures", "07_celltype_correlation.pdf"),
        width = max(8, ncol(cor_mat) * 0.8 + 2),
        height = max(8, ncol(cor_mat) * 0.8 + 2)
      )
      corrplot(
        cor_mat,
        method = "color",
        type = "upper",
        tl.col = "black",
        tl.srt = 45,
        addCoef.col = "black",
        number.cex = 0.75,
        col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(200),
        title = "Cell Type Fraction Correlation (Pearson)",
        mar = c(0, 0, 2, 0)
      )
      dev.off()
      cat("  Saved: 07_celltype_correlation.pdf\n")

      cor_pairs_idx <- which(
        abs(cor_mat) > 0.7 & lower.tri(cor_mat),
        arr.ind = TRUE
      )
      if (nrow(cor_pairs_idx) > 0) {
        cor_pairs_dt <- data.table(
          CellType1 = rownames(cor_mat)[cor_pairs_idx[, 1]],
          CellType2 = colnames(cor_mat)[cor_pairs_idx[, 2]],
          Pearson_r = cor_mat[cor_pairs_idx]
        )
        cor_pairs_dt <- cor_pairs_dt[order(-abs(Pearson_r))]
        fwrite(
          cor_pairs_dt,
          file.path(OUTPUT_DIR, "tables", "celltype_strong_correlations.csv")
        )
        cat(sprintf(
          "  Strong correlations (|r|>0.7): %d pairs\n",
          nrow(cor_pairs_dt)
        ))
      }
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      warning("Correlation: ", e$message)
    }
  )
} else {
  cat(sprintf(
    "  Skipping correlation: need >= 3 samples and >= 2 variable cell types\n"
  ))
}

# ==============================================================================
# 5. Z Matrix Analysis
# ==============================================================================

cat("\n=== Step 4/7: Z Matrix Analysis ===\n")

# Technical gene filter pattern (shared across all gene-selection steps)
technical_pattern <- "^MT-|^RPS|^RPL|^MRPS|^MRPL|^HBA|^HBB|^HBG|^HSPA|^HSPB|^MALAT1|^NEAT1|^XIST"

# ---- 5a. Marker gene validation heatmap ----
marker_genes <- list(
  Epithelial = c(
    "EPCAM",
    "KRT5",
    "TP63",
    "KRT14",
    "MUC5AC",
    "MUC5B",
    "FOXJ1",
    "RSPH1",
    "SCGB1A1",
    "SCGB3A1"
  ),
  T_NK = c(
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    "CD8B",
    "GNLY",
    "NKG7",
    "GZMA",
    "GZMB",
    "KLRD1"
  ),
  Myeloid = c(
    "CD14",
    "FCGR3A",
    "CD68",
    "CD163",
    "CD1C",
    "CLEC9A",
    "FCER1A",
    "ITGAX",
    "MRC1",
    "CCL18"
  ),
  B_Plasma = c(
    "MS4A1",
    "CD79A",
    "CD79B",
    "IGHG1",
    "IGHM",
    "IGHA1",
    "MZB1",
    "XBP1",
    "PRDM1"
  ),
  Stromal = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "ACTA2",
    "PDGFRA",
    "PDGFRB",
    "FAP"
  ),
  Endothelial = c(
    "PECAM1",
    "VWF",
    "CDH5",
    "CLDN5",
    "FLT1",
    "KDR",
    "ACKR1",
    "RAMP2"
  ),
  Mast = c("KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2", "FCER1G"),
  ILC_Eosinophil = c(
    "IL5RA",
    "SIGLEC8",
    "CCR3",
    "PRG2",
    "EPX",
    "GATA2",
    "PTGDR2"
  )
)

all_markers <- unique(unlist(marker_genes))
present_mask <- all_markers %in% rownames(Z)
markers_in_Z <- all_markers[present_mask]
cat(sprintf(
  "  Marker genes in Z: %d / %d (%.0f%%)\n",
  sum(present_mask),
  length(all_markers),
  100 * mean(present_mask)
))

if (length(markers_in_Z) >= 5) {
  Z_mark <- Z[markers_in_Z, , drop = FALSE]
  ann_row <- data.frame(
    Lineage = rep(names(marker_genes), sapply(marker_genes, length))
  )
  rownames(ann_row) <- unlist(marker_genes)
  ann_row <- ann_row[rownames(Z_mark), , drop = FALSE]

  lineage_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set1"))(length(marker_genes)),
    names(marker_genes)
  )

  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "08_Z_marker_heatmap.pdf"),
        width = max(10, n_cell_types * 1.0 + 4),
        height = max(8, length(markers_in_Z) * 0.28 + 3)
      )
      pheatmap(
        log1p(Z_mark),
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_row = ann_row,
        annotation_colors = list(Lineage = lineage_colors),
        color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
        border_color = NA,
        main = "Deconvolved Marker Gene Expression\n(log1p, row-scaled Z matrix)",
        fontsize_row = 7,
        fontsize_col = 10,
        angle_col = 45
      )
      dev.off()
      cat("  Saved: 08_Z_marker_heatmap.pdf\n")
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      warning("Marker heatmap: ", e$message)
    }
  )

  fwrite(
    as.data.table(Z_mark, keep.rownames = "gene"),
    file.path(OUTPUT_DIR, "tables", "Z_marker_expression.csv")
  )
}

# ---- 5b. Top cell-type-specific genes (for enrichment) ----
# [v1.1 P0] Use specificity score (log2 ratio vs other cell types) instead of
#            absolute expression, to avoid housekeeping gene enrichment bias.
# [v1.1 P1] Removed dead Z_var line from v1.0.

get_top_specific_genes <- function(Z, cell_types, n_top = 100, tech_pat) {
  res <- lapply(cell_types, function(ct) {
    ct_expr <- Z[, ct]
    other_cols <- setdiff(cell_types, ct)

    if (length(other_cols) > 0) {
      other_mean <- rowMeans(Z[, other_cols, drop = FALSE])
    } else {
      other_mean <- rep(0, nrow(Z))
    }

    specificity <- log2((ct_expr + 1e-8) / (other_mean + 1e-8))

    dt <- data.table(
      gene = rownames(Z),
      cell_type = ct,
      expression = as.numeric(ct_expr),
      specificity = as.numeric(specificity)
    )

    dt <- dt[is.finite(specificity)]
    dt <- dt[!grepl(tech_pat, gene)]
    dt[, gene := toupper(gene)]
    dt <- dt[!duplicated(gene)]
    setorder(dt, -specificity, -expression)
    dt[seq_len(min(n_top, .N))]
  })
  rbindlist(res, use.names = TRUE)
}

top_genes_dt <- get_top_specific_genes(
  Z = Z,
  cell_types = cell_types,
  n_top = N_TOP_GENES_ENRICHMENT,
  tech_pat = technical_pattern
)

fwrite(
  top_genes_dt,
  file.path(OUTPUT_DIR, "tables", "top_specific_genes_per_celltype.csv")
)
cat(sprintf(
  "  Top specific genes per cell type: %d rows\n",
  nrow(top_genes_dt)
))

# ==============================================================================
# 6. GMT-Based Functional Enrichment
# ==============================================================================

cat("\n=== Step 5/7: GMT-Based Functional Enrichment ===\n")

# Helper: load GMT as TERM2GENE data.table
# [v1.1 P1] toupper() genes for case-consistent matching with top_genes_dt
load_gmt_as_term2gene <- function(gmt_path, name = "GMT") {
  if (!file.exists(gmt_path)) {
    cat(sprintf("  [WARN] GMT not found: %s\n", gmt_path))
    return(NULL)
  }
  lines <- readLines(gmt_path, warn = FALSE)
  term2gene_list <- lapply(lines, function(ln) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) {
      return(NULL)
    }
    term_name <- parts[1]
    genes <- toupper(parts[-(1:2)]) # [v1.1] force uppercase
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) == 0) {
      return(NULL)
    }
    data.table(term = term_name, gene = genes)
  })
  term2gene_list <- term2gene_list[!vapply(term2gene_list, is.null, logical(1))]
  dt <- rbindlist(term2gene_list)
  cat(sprintf(
    "  Loaded %s: %d terms, %d genes\n",
    name,
    length(unique(dt$term)),
    length(unique(dt$gene))
  ))
  dt
}

# Helper: run compareCluster for one GMT
# [v1.1 P1+P0] Added universe_genes for proper background; used by all enrichment calls
run_compare_cluster <- function(
  gene_celltype_df,
  term2gene_dt,
  db_name,
  universe_genes = NULL,
  pval = ENRICH_PVAL_CUTOFF,
  qval = ENRICH_QVAL_CUTOFF,
  minGS = MIN_GENE_SET_SIZE,
  maxGS = MAX_GENE_SET_SIZE
) {
  if (is.null(term2gene_dt) || nrow(term2gene_dt) == 0) {
    return(NULL)
  }
  tryCatch(
    {
      t0 <- proc.time()["elapsed"]
      res <- compareCluster(
        gene ~ cell_type,
        data = gene_celltype_df,
        fun = enricher,
        TERM2GENE = as.data.frame(term2gene_dt),
        universe = universe_genes,
        pvalueCutoff = pval,
        qvalueCutoff = qval,
        minGSSize = minGS,
        maxGSSize = maxGS
      )
      elapsed <- round(proc.time()["elapsed"] - t0, 1)
      n_sig <- if (!is.null(res)) nrow(res@compareClusterResult) else 0
      cat(sprintf(
        "    %s: %d significant terms (%.1f s)\n",
        db_name,
        n_sig,
        elapsed
      ))
      res
    },
    error = function(e) {
      cat(sprintf("    [WARN] %s failed: %s\n", db_name, e$message))
      NULL
    }
  )
}

# Load GMT databases
cat("  Loading GMT databases...\n")
term2gene_hallmark <- load_gmt_as_term2gene(GMT_HALLMARK, "Hallmark")
term2gene_go_all <- load_gmt_as_term2gene(GMT_GO_ALL, "GO-all (BP+MF+CC)")

# Subset GO sub-ontologies from combined file
if (!is.null(term2gene_go_all)) {
  term2gene_go_bp <- term2gene_go_all[grepl("^GOBP_", term)]
  term2gene_go_mf <- term2gene_go_all[grepl("^GOMF_", term)]
  term2gene_go_cc <- term2gene_go_all[grepl("^GOCC_", term)]
  cat(sprintf(
    "    GO BP: %d | GO MF: %d | GO CC: %d terms\n",
    length(unique(term2gene_go_bp$term)),
    length(unique(term2gene_go_mf$term)),
    length(unique(term2gene_go_cc$term))
  ))
} else {
  term2gene_go_bp <- term2gene_go_mf <- term2gene_go_cc <- NULL
}

# Background universe: all non-technical genes detected in Z
# [v1.1 P0] Proper ORA background; prevents inflated p-values
gene_universe <- unique(toupper(rownames(Z)[
  !grepl(technical_pattern, rownames(Z))
]))
cat(sprintf("  Background universe: %d genes\n", length(gene_universe)))

# Input gene-celltype data frame
gene_ct_df <- unique(as.data.frame(top_genes_dt[, .(gene, cell_type)]))

enrichment_results <- list()
cat("  Running enrichment (sequential)...\n")
enrichment_results[["Hallmark"]] <- run_compare_cluster(
  gene_ct_df,
  term2gene_hallmark,
  "Hallmark",
  gene_universe
)
enrichment_results[["GO_BP"]] <- run_compare_cluster(
  gene_ct_df,
  term2gene_go_bp,
  "GO_BP",
  gene_universe
)
enrichment_results[["GO_MF"]] <- run_compare_cluster(
  gene_ct_df,
  term2gene_go_mf,
  "GO_MF",
  gene_universe
)
enrichment_results[["GO_CC"]] <- run_compare_cluster(
  gene_ct_df,
  term2gene_go_cc,
  "GO_CC",
  gene_universe
)

# Save enrichment objects
saveRDS(
  enrichment_results,
  file.path(OUTPUT_DIR, "reports", "enrichment_results.rds")
)
cat("  Saved: enrichment_results.rds\n")

# Export CSVs + dotplots
for (db_name in names(enrichment_results)) {
  res <- enrichment_results[[db_name]]
  if (is.null(res) || nrow(res@compareClusterResult) == 0) {
    next
  }

  fwrite(
    as.data.table(res@compareClusterResult),
    file.path(OUTPUT_DIR, "tables", sprintf("enrichment_%s.csv", db_name))
  )

  # [v1.1 P0] Floor p.adjust before log10 color scale to prevent Inf/-Inf crash
  tryCatch(
    {
      plot_data <- as.data.table(res@compareClusterResult)

      plot_data <- plot_data %>%
        group_by(Cluster) %>%
        slice_min(
          order_by = p.adjust,
          n = N_TOP_TERMS_PLOT,
          with_ties = FALSE
        ) %>%
        ungroup()
      plot_data <- as.data.table(plot_data)

      plot_data[,
        p_adjust_plot := pmax(as.numeric(p.adjust), .Machine$double.xmin)
      ]

      n_terms <- length(unique(plot_data$Description))
      n_clusters <- length(unique(plot_data$Cluster))
      pw <- max(10, min(30, n_clusters * 1.2 + 4))
      ph <- max(6, min(28, n_terms * 0.32 + 3))

      p_dot <- ggplot(
        plot_data,
        aes(x = Cluster, y = Description, size = Count, color = p_adjust_plot)
      ) +
        geom_point() +
        scale_color_gradientn(
          colors = c("#b2182b", "#ef8a62", "#fddbc7", "#d1e5f0", "#2166ac"),
          name = "adj. p-value",
          trans = "log10"
        ) +
        scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
        labs(
          title = sprintf(
            "Functional Enrichment: %s\n(Top %d terms per cell type)",
            db_name,
            N_TOP_TERMS_PLOT
          ),
          x = "Cell Type",
          y = NULL
        ) +
        theme_bw(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 8),
          panel.grid.minor = element_blank()
        )
      ggsave(
        file.path(
          OUTPUT_DIR,
          "figures",
          sprintf("09_enrichment_%s_dotplot.pdf", db_name)
        ),
        p_dot,
        width = pw,
        height = ph,
        limitsize = FALSE
      )
      cat(sprintf("  Saved: 09_enrichment_%s_dotplot.pdf\n", db_name))
    },
    error = function(e) warning(sprintf("Dotplot %s: %s", db_name, e$message))
  )
}

# ==============================================================================
# 7. DeepSeek LLM Interpretation
# ==============================================================================

cat("\n=== Step 6/7: LLM Interpretation (DeepSeek) ===\n")

# Helper: call DeepSeek API
# [v1.1 P1] encode="json" instead of raw+toJSON (cleaner, more robust)
call_deepseek <- function(
  prompt,
  api_key,
  model = DEEPSEEK_MODEL,
  max_tokens = LLM_MAX_TOKENS,
  temperature = LLM_TEMPERATURE
) {
  if (!nzchar(api_key)) {
    stop("DEEPSEEK_API_KEY is not set")
  }

  body <- list(
    model = model,
    messages = list(list(role = "user", content = prompt)),
    max_tokens = max_tokens,
    temperature = temperature
  )

  resp <- POST(
    url = "https://api.deepseek.com/v1/chat/completions",
    add_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ),
    body = body,
    encode = "json", # [v1.1] httr handles serialization
    timeout(120)
  )

  if (http_error(resp)) {
    stop(sprintf(
      "API error %d: %s",
      status_code(resp),
      content(resp, "text", encoding = "UTF-8")
    ))
  }

  content(resp, "parsed", simplifyVector = FALSE)$choices[[1]]$message$content
}

# Helper: parse LLM response
# [v1.1 P0] Line-by-line state machine: handles multi-line Reasoning correctly
parse_llm_response <- function(text) {
  fields <- c(
    "Cell Type Identity",
    "Confidence",
    "Key Functions",
    "Disease Relevance",
    "Reasoning"
  )
  result <- setNames(as.list(rep(NA_character_, length(fields))), fields)
  lines <- unlist(strsplit(text, "\n", fixed = TRUE))
  current_field <- NULL

  for (ln in lines) {
    ln_trim <- trimws(ln)
    if (!nzchar(ln_trim)) {
      next
    }

    matched <- FALSE
    for (f in fields) {
      prefix <- paste0(f, ":")
      if (startsWith(tolower(ln_trim), tolower(prefix))) {
        value <- trimws(substr(ln_trim, nchar(prefix) + 1, nchar(ln_trim)))
        result[[f]] <- value
        current_field <- f
        matched <- TRUE
        break
      }
    }
    # Continuation: append to current field (handles multi-line Reasoning)
    if (
      !matched && !is.null(current_field) && !is.na(result[[current_field]])
    ) {
      result[[current_field]] <- trimws(paste(result[[current_field]], ln_trim))
    }
  }
  result
}

# Helper: build prompt for one cell type
# [v1.1 P1] Use specificity column (not expression) for gene ordering in prompt
build_celltype_prompt <- function(
  ct,
  top_genes,
  enrichment_results,
  theta_row,
  study_ctx
) {
  mean_frac <- round(mean(theta_row), 4)
  cv_val <- round(sd(theta_row) / (mean(theta_row) + 1e-10), 3)

  ct_sub <- top_genes[cell_type == ct][order(-specificity, -expression)]
  ct_genes <- ct_sub[seq_len(min(30, .N)), gene]
  genes_str <- paste(ct_genes, collapse = ", ")

  enrich_lines <- c()
  for (db_name in names(enrichment_results)) {
    res <- enrichment_results[[db_name]]
    if (is.null(res)) {
      next
    }
    ccr <- as.data.table(res@compareClusterResult)
    ct_terms <- ccr[Cluster == ct][order(p.adjust)][seq_len(min(5, .N))]
    if (nrow(ct_terms) == 0) {
      next
    }
    enrich_lines <- c(
      enrich_lines,
      sprintf(
        "[%s] %s",
        db_name,
        paste(head(ct_terms$Description, 5), collapse = "; ")
      )
    )
  }
  enrich_str <- if (length(enrich_lines) > 0) {
    paste(enrich_lines, collapse = "\n")
  } else {
    "No significant enrichment detected."
  }

  sprintf(
    "You are an expert in respiratory biology and single-cell genomics.

Study context:
%s

Analyze the deconvolved cell type '%s' from BayesPrism bulk RNA-seq deconvolution.

Quantitative summary:
- Estimated mean fraction across samples: %.1f%%
- Sample-to-sample CV: %.3f (higher = more variable)

Top cell-type-specific genes (ranked by specificity score):
%s

Functional enrichment (top terms per database):
%s

Provide a concise biological interpretation using EXACTLY this format:

Cell Type Identity: [specific cell subtype name or state]
Confidence: [High / Medium / Low]
Key Functions: [3-5 key biological functions, semicolon-separated]
Disease Relevance: [relevance to airway inflammation / CRSwNP / Type 2 inflammation]
Reasoning: [2-3 sentences explaining the evidence and interpretation]

IMPORTANT: Include all five fields. Keep total response under 250 words.",
    study_ctx,
    ct,
    mean_frac * 100,
    cv_val,
    genes_str,
    enrich_str
  )
}

# ---- Run LLM for all cell types ----
llm_results <- list()
llm_raw_text <- list()

# [v1.1 P1] Use nzchar() for robust empty-string check
if (RUN_LLM_INTERPRETATION && nzchar(DEEPSEEK_API_KEY)) {
  cat(sprintf(
    "  Running LLM for %d cell types (sequential)...\n",
    n_cell_types
  ))

  for (ct in cell_types) {
    cat(sprintf("    Processing: %s ... ", ct))

    success <- FALSE
    attempts <- 0
    wait_sec <- 5

    while (!success && attempts < 4) {
      attempts <- attempts + 1
      result <- tryCatch(
        {
          prompt <- build_celltype_prompt(
            ct = ct,
            top_genes = top_genes_dt,
            enrichment_results = enrichment_results,
            theta_row = theta[, ct],
            study_ctx = STUDY_CONTEXT
          )
          raw_text <- call_deepseek(prompt, api_key = DEEPSEEK_API_KEY)
          list(raw = raw_text, parsed = parse_llm_response(raw_text), ok = TRUE)
        },
        error = function(e) {
          list(raw = NULL, parsed = NULL, ok = FALSE, err = e$message)
        }
      )

      if (result$ok) {
        llm_raw_text[[ct]] <- result$raw
        llm_results[[ct]] <- result$parsed
        success <- TRUE
        cat("OK\n")
      } else {
        cat(sprintf("[attempt %d: %s] ", attempts, result$err))
        Sys.sleep(wait_sec)
        wait_sec <- wait_sec * 2
      }
    }

    if (!success) {
      cat("FAILED\n")
      llm_results[[ct]] <- setNames(
        as.list(rep(NA_character_, 5)),
        c(
          "Cell Type Identity",
          "Confidence",
          "Key Functions",
          "Disease Relevance",
          "Reasoning"
        )
      )
    }
    Sys.sleep(1.5) # polite rate limiting
  }

  llm_dt <- rbindlist(lapply(cell_types, function(ct) {
    r <- llm_results[[ct]]
    data.table(
      Cell_Type = ct,
      Identity = r[["Cell Type Identity"]] %||% NA_character_,
      Confidence = r[["Confidence"]] %||% NA_character_,
      Key_Functions = r[["Key Functions"]] %||% NA_character_,
      Disease_Relevance = r[["Disease Relevance"]] %||% NA_character_,
      Reasoning = r[["Reasoning"]] %||% NA_character_
    )
  }))

  fwrite(
    llm_dt,
    file.path(OUTPUT_DIR, "tables", "LLM_celltype_interpretation.csv")
  )
  saveRDS(
    list(results = llm_results, raw = llm_raw_text),
    file.path(OUTPUT_DIR, "reports", "LLM_results_full.rds")
  )
  cat(sprintf(
    "  LLM interpretation: %d / %d successful\n",
    sum(!is.na(llm_dt$Identity)),
    n_cell_types
  ))
} else {
  cat(
    "  [SKIPPED] RUN_LLM_INTERPRETATION = FALSE or DEEPSEEK_API_KEY not set\n"
  )
  llm_dt <- NULL
}

# ==============================================================================
# 8. DOWNSTREAM_REPORT.md
# ==============================================================================

cat("\n=== Step 7/7: Generating DOWNSTREAM_REPORT.md ===\n")

report_path <- file.path(OUTPUT_DIR, "DOWNSTREAM_REPORT.md")

tryCatch(
  {
    sink(report_path)

    cat("# BayesPrism Downstream Analysis Report\n\n")
    cat(sprintf("**Date:** %s  \n", format(Sys.Date(), "%Y-%m-%d")))
    cat(sprintf("**Version:** v1.1  \n"))
    cat(sprintf("**Input dir:** %s  \n", BAYESPRISM_OUTPUT_DIR))
    cat(sprintf("**Output dir:** %s  \n\n", OUTPUT_DIR))
    cat("---\n\n")

    cat("## 1. Dataset Summary\n\n")
    cat(sprintf(
      "- **Samples:** %d (`%s`)  \n",
      n_samples,
      paste(rownames(theta), collapse = "`, `")
    ))
    cat(sprintf("- **Cell types:** %d  \n", n_cell_types))
    cat(sprintf("- **Z matrix genes:** %d  \n\n", nrow(Z)))
    cat(sprintf(
      "> Note: n = %d -> all results are descriptive/exploratory.\n\n",
      n_samples
    ))

    cat("## 2. Cell Type Fractions (Mean across samples)\n\n")
    cat("| Cell Type | Mean Fraction | SD | Between-sample CV |\n")
    cat("|-----------|:-------------:|:--:|:-----------------:|\n")
    for (ct in cell_types) {
      m <- mean(theta[, ct])
      s <- sd(theta[, ct])
      cv <- s / (m + 1e-10)
      cat(sprintf("| %s | %.3f | %.3f | %.3f |\n", ct, m, s, cv))
    }
    cat("\n")

    cat("## 3. Functional Enrichment Summary\n\n")
    for (db_name in names(enrichment_results)) {
      res <- enrichment_results[[db_name]]
      if (is.null(res)) {
        cat(sprintf(
          "### %s\nNo results (GMT not found or no significant hits).\n\n",
          db_name
        ))
        next
      }
      ccr <- as.data.table(res@compareClusterResult)
      cat(sprintf("### %s (%d significant terms)\n\n", db_name, nrow(ccr)))
      for (ct in cell_types) {
        ct_top <- ccr[Cluster == ct][order(p.adjust)][seq_len(min(3, .N))]
        if (nrow(ct_top) == 0) {
          next
        }
        cat(sprintf(
          "**%s:** %s  \n",
          ct,
          paste(ct_top$Description, collapse = "; ")
        ))
      }
      cat("\n")
    }

    cat("## 4. LLM Cell Type Interpretations\n\n")
    if (!is.null(llm_dt) && nrow(llm_dt) > 0) {
      for (i in seq_len(nrow(llm_dt))) {
        r <- llm_dt[i]
        cat(sprintf("### %s\n\n", r$Cell_Type))
        cat(sprintf("- **Identity:** %s  \n", r$Identity %||% "N/A"))
        cat(sprintf("- **Confidence:** %s  \n", r$Confidence %||% "N/A"))
        cat(sprintf("- **Key Functions:** %s  \n", r$Key_Functions %||% "N/A"))
        cat(sprintf(
          "- **Disease Relevance:** %s  \n",
          r$Disease_Relevance %||% "N/A"
        ))
        cat(sprintf("- **Reasoning:** %s  \n\n", r$Reasoning %||% "N/A"))
      }
    } else {
      cat("LLM interpretation was not run or returned no results.\n\n")
    }

    cat("## 5. Output Files\n\n")
    for (sub in c("figures", "tables", "reports")) {
      cat(sprintf("### %s\n", sub))
      fs <- list.files(file.path(OUTPUT_DIR, sub), full.names = FALSE)
      for (f in sort(fs)) {
        cat(sprintf("- `%s/%s`\n", sub, f))
      }
      cat("\n")
    }

    cat("---\n\n## 6. Session Info\n\n```\n")
    print(sessionInfo())
    cat("```\n")

    sink()
    cat("  Saved: DOWNSTREAM_REPORT.md\n")
  },
  error = function(e) {
    while (sink.number() > 0) {
      sink()
    }
    warning("Report generation: ", e$message)
  },
  finally = {
    while (sink.number() > 0) {
      sink()
    }
  }
)

# ==============================================================================
# Final Summary
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n", sep = "")
cat("BayesPrism Downstream Analysis v1.1 COMPLETE\n")
cat(strrep("=", 80), "\n", sep = "")
cat(sprintf("Output: %s\n\n", OUTPUT_DIR))
cat("Key outputs:\n")
cat("  figures/01_stacked_bar.pdf              <- composition overview\n")
cat(
  "  figures/02_fraction_jitter.pdf          <- per-cell-type distributions\n"
)
cat("  figures/03_theta_heatmap.pdf            <- fraction heatmap\n")
cat("  figures/04_theta_cv_heatmap.pdf         <- estimation uncertainty\n")
cat("  figures/05_theta_PCA.pdf                <- sample PCA\n")
cat("  figures/06_celltype_CV.pdf              <- variability ranking\n")
cat("  figures/07_celltype_correlation.pdf     <- inter-type correlation\n")
cat("  figures/08_Z_marker_heatmap.pdf         <- marker gene validation\n")
cat("  figures/09_enrichment_*_dotplot.pdf     <- functional enrichment\n")
cat(
  "  tables/top_specific_genes_per_celltype.csv <- specificity-ranked genes\n"
)
cat("  tables/LLM_celltype_interpretation.csv  <- DeepSeek interpretation\n")
cat("  DOWNSTREAM_REPORT.md                    <- integrated report\n")
cat(strrep("=", 80), "\n", sep = "")


# ==============================================================================
# Extended Marker Heatmap: Ionocyte + Tuft Cell
# ==============================================================================
# Appended block — shares Z and all variables already in memory.
# Regenerates the marker heatmap with two rare respiratory epithelial cell types
# added: Ionocyte (FOXI1+/CFTR+) and Tuft_Cell (POU2F3+/TRPM5+).
# Output: figures/08b_Z_marker_heatmap_extended.pdf
# ==============================================================================

cat("\n=== Extended Marker Heatmap: Ionocyte + Tuft Cell ===\n")

marker_genes_extended <- list(
  # ---- Original 8 groups (identical to Step 4) ----
  Epithelial = c(
    "EPCAM",
    "KRT5",
    "TP63",
    "KRT14",
    "MUC5AC",
    "MUC5B",
    "FOXJ1",
    "RSPH1",
    "SCGB1A1",
    "SCGB3A1"
  ),
  T_NK = c(
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    "CD8B",
    "GNLY",
    "NKG7",
    "GZMA",
    "GZMB",
    "KLRD1"
  ),
  Myeloid = c(
    "CD14",
    "FCGR3A",
    "CD68",
    "CD163",
    "CD1C",
    "CLEC9A",
    "FCER1A",
    "ITGAX",
    "MRC1",
    "CCL18"
  ),
  B_Plasma = c(
    "MS4A1",
    "CD79A",
    "CD79B",
    "IGHG1",
    "IGHM",
    "IGHA1",
    "MZB1",
    "XBP1",
    "PRDM1"
  ),
  Stromal = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "ACTA2",
    "PDGFRA",
    "PDGFRB",
    "FAP"
  ),
  Endothelial = c(
    "PECAM1",
    "VWF",
    "CDH5",
    "CLDN5",
    "FLT1",
    "KDR",
    "ACKR1",
    "RAMP2"
  ),
  Mast = c(
    "KIT",
    "TPSAB1",
    "TPSB2",
    "CPA3",
    "MS4A2",
    "FCER1G"
  ),
  ILC_Eosinophil = c(
    "IL5RA",
    "SIGLEC8",
    "CCR3",
    "PRG2",
    "EPX",
    "GATA2",
    "PTGDR2"
  ),
  # ---- Two new groups ----
  # Ionocyte: rare airway cells for ion/pH regulation (CFTR-expressing)
  Ionocyte = c(
    "FOXI1", # master TF, most specific ionocyte marker
    "CFTR", # chloride channel, clinically important
    "ATP6V1G3", # vacuolar H+-ATPase subunit
    "ATP6V0D2", # vacuolar H+-ATPase subunit
    "ASCL3", # transcription factor
    "BSND", # barttin (Cl- channel subunit)
    "CLCNKB" # kidney/airway chloride channel
  ),
  # Tuft_Cell (Brush cell): chemosensory; POU2F3 master TF
  Tuft_Cell = c(
    "POU2F3", # master TF, definitive tuft marker
    "TRPM5", # taste signal transduction channel
    "DCLK1", # tuft cell kinase marker
    "CHAT", # choline acetyltransferase
    "SH2D6", # tuft-specific signaling adaptor
    "AVIL", # actin-binding, brush border
    "LRMP" # lymphoid-restricted membrane protein
  )
)

all_markers_ext <- unique(unlist(marker_genes_extended))
present_mask_ext <- all_markers_ext %in% rownames(Z)
markers_in_Z_ext <- all_markers_ext[present_mask_ext]

# Report per-group detection
for (grp in names(marker_genes_extended)) {
  grp_markers <- marker_genes_extended[[grp]]
  n_found <- sum(grp_markers %in% rownames(Z))
  cat(sprintf(
    "  %-20s %d / %d markers detected\n",
    grp,
    n_found,
    length(grp_markers)
  ))
}

if (length(markers_in_Z_ext) < 5) {
  warning("Extended heatmap: fewer than 5 markers found in Z — skipping")
} else {
  Z_mark_ext <- Z[markers_in_Z_ext, , drop = FALSE]

  # Row annotation
  ann_row_ext <- data.frame(
    Lineage = rep(
      names(marker_genes_extended),
      sapply(marker_genes_extended, length)
    )
  )
  rownames(ann_row_ext) <- unlist(marker_genes_extended)
  ann_row_ext <- ann_row_ext[rownames(Z_mark_ext), , drop = FALSE]

  # Color palette: 10 groups — use a 10-color qualitative palette
  n_groups_ext <- length(marker_genes_extended)
  lineage_pal_ext <- c(
    brewer.pal(8, "Set1"),
    brewer.pal(max(3, n_groups_ext - 8), "Dark2")
  )[seq_len(n_groups_ext)]
  lineage_colors_ext <- setNames(lineage_pal_ext, names(marker_genes_extended))

  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "08b_Z_marker_heatmap_extended.pdf"),
        width = max(10, n_cell_types * 1.0 + 4),
        height = max(8, length(markers_in_Z_ext) * 0.28 + 3)
      )
      pheatmap(
        log1p(Z_mark_ext),
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_row = ann_row_ext,
        annotation_colors = list(Lineage = lineage_colors_ext),
        color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
        border_color = NA,
        main = paste0(
          "Deconvolved Marker Gene Expression (Extended)\n",
          "log1p, row-scaled Z matrix | + Ionocyte & Tuft Cell"
        ),
        fontsize_row = 7,
        fontsize_col = 10,
        angle_col = 45
      )
      dev.off()
      cat("  Saved: figures/08b_Z_marker_heatmap_extended.pdf\n")
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      warning("Extended marker heatmap: ", e$message)
    }
  )

  fwrite(
    as.data.table(Z_mark_ext, keep.rownames = "gene"),
    file.path(OUTPUT_DIR, "tables", "Z_marker_expression_extended.csv")
  )
  cat("  Saved: tables/Z_marker_expression_extended.csv\n")
  cat(sprintf(
    "  Total markers in plot: %d / %d (%.0f%%)\n",
    length(markers_in_Z_ext),
    length(all_markers_ext),
    100 * mean(present_mask_ext)
  ))
}
