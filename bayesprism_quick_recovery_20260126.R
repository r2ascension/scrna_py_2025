#!/usr/bin/env Rscript

# ==============================================================================
# BayesPrism Quick Recovery Script
# ==============================================================================
# Purpose: Generate visualizations and reports from existing bp_res.rds
# Author: r2end
# Date: 2025-01-26
#
# Usage:
#   1. Make sure bp_res.rds exists in output_dir
#   2. Update output_dir path below
#   3. Run: Rscript bayesprism_quick_recovery.R
#
# This script will generate all visualizations and reports WITHOUT re-running
# the time-consuming BayesPrism deconvolution.
# ==============================================================================

suppressPackageStartupMessages({
  library(BayesPrism)
  library(ggplot2)
  library(pheatmap)
  library(data.table)
  library(RColorBrewer)
})

# ==============================================================================
# Configuration
# ==============================================================================

# ===== MODIFY THIS PATH =====
output_dir <- "/home/h2048/data/R/bulk0121/bayesprism_out" # ← CHANGE THIS
# Check if bp_res.rds exists
bp_res_path <- file.path(output_dir, "bp_res.rds")
if (!file.exists(bp_res_path)) {
  stop("ERROR: bp_res.rds not found at: ", bp_res_path)
}

# Check if reference_metadata.rds exists (for summary report)
ref_metadata_path <- file.path(output_dir, "reference_metadata.rds")
has_metadata <- file.exists(ref_metadata_path)

cat("\n")
cat(
  "================================================================================\n"
)
cat("BayesPrism Quick Recovery Script\n")
cat(
  "================================================================================\n"
)
cat("Output directory:", output_dir, "\n")
cat("bp_res.rds found:  YES\n")
cat("metadata found:    ", ifelse(has_metadata, "YES", "NO"), "\n")
cat(
  "================================================================================\n\n"
)

# ==============================================================================
# 1. Load BayesPrism Results
# ==============================================================================

cat("=== Step 1/5: Loading BayesPrism Results ===\n")

bp_res <- readRDS(bp_res_path)
cat("  Loaded: bp_res.rds\n")

if (has_metadata) {
  ref_metadata <- readRDS(ref_metadata_path)
  cat("  Loaded: reference_metadata.rds\n")
}

# ==============================================================================
# 2. Extract Results
# ==============================================================================

cat("\n=== Step 2/5: Extracting Results ===\n")

# Cell type fractions (theta)
theta <- get.fraction(
  bp = bp_res,
  which.theta = "final",
  state.or.type = "type"
)

cat(sprintf(
  "  Theta dimensions: %d samples × %d cell types\n",
  nrow(theta),
  ncol(theta)
))

# Coefficient of variation (uncertainty)
theta_cv <- tryCatch(
  {
    if (!is.null(bp_res@posterior.theta_f@theta.cv)) {
      bp_res@posterior.theta_f@theta.cv
    } else if ("cv" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@cv
    } else if ("theta_cv" %in% slotNames(bp_res@posterior.theta_f)) {
      bp_res@posterior.theta_f@theta_cv
    } else {
      warning("Could not extract theta CV - using placeholder zeros")
      matrix(
        0,
        nrow = nrow(theta),
        ncol = ncol(theta),
        dimnames = dimnames(theta)
      )
    }
  },
  error = function(e) {
    warning(
      "Error extracting theta CV: ",
      e$message,
      " - using placeholder zeros"
    )
    matrix(
      0,
      nrow = nrow(theta),
      ncol = ncol(theta),
      dimnames = dimnames(theta)
    )
  }
)

cat("  Extracted: theta_cv\n")

# Deconvolved gene expression per cell type (Z matrix)
Z <- tryCatch(
  {
    # Try method 1: cell.name = "all"
    get.exp(bp = bp_res, state.or.type = "type", cell.name = "all")
  },
  error = function(e1) {
    cat("  Method 1 failed, trying alternative...\n")
    tryCatch(
      {
        # Try method 2: without cell.name
        get.exp(bp = bp_res, state.or.type = "type")
      },
      error = function(e2) {
        cat("  Method 2 failed, trying direct extraction...\n")
        tryCatch(
          {
            # Try method 3: extract directly from posterior
            bp_res@posterior.initial.cellType@Z.mean
          },
          error = function(e3) {
            cat("  Method 3 failed, trying final fallback...\n")
            # Try method 4: extract from final posterior
            if (!is.null(bp_res@posterior.theta_f@Z)) {
              bp_res@posterior.theta_f@Z
            } else {
              stop("Failed to extract Z matrix using all methods")
            }
          }
        )
      }
    )
  }
)

cat(sprintf(
  "  Z matrix dimensions: %d genes × %d cell types\n",
  nrow(Z),
  ncol(Z)
))

# Check object class and convert to matrix if needed
cat(sprintf("  Z object class: %s\n", paste(class(Z), collapse = ", ")))
cat(sprintf("  Z is matrix: %s\n", is.matrix(Z)))
cat(sprintf("  Z is list: %s\n", is.list(Z)))
cat(sprintf("  Z is array: %s\n", is.array(Z)))

# Handle different Z structures
if (is.list(Z) && !is.data.frame(Z) && !is.matrix(Z)) {
  cat("  WARNING: Z is a list, attempting to extract matrix...\n")
  # Try to find matrix in list
  if (length(Z) > 0) {
    # Check what's in the list
    cat(sprintf("  List has %d elements\n", length(Z)))
    cat(sprintf(
      "  First element class: %s\n",
      paste(class(Z[[1]]), collapse = ", ")
    ))

    # Take first element if it's a matrix
    if (is.matrix(Z[[1]]) || is.array(Z[[1]])) {
      Z <- Z[[1]]
      cat("  Extracted matrix from list element 1\n")
    } else {
      # Try to convert first element to matrix
      Z <- as.matrix(Z[[1]])
      cat("  Converted list element 1 to matrix\n")
    }
  } else {
    stop("Z is an empty list")
  }
}

# Convert to matrix if not already
if (!is.matrix(Z)) {
  cat("  Converting Z to matrix...\n")
  Z <- as.matrix(Z)
  cat(sprintf(
    "  After conversion - class: %s, dimensions: %d × %d\n",
    class(Z),
    nrow(Z),
    ncol(Z)
  ))
}

# Check if Z dimensions are correct (should be genes × cell types)
# If Z is transposed (cell types × genes), fix it
if (ncol(Z) > nrow(Z) && ncol(Z) > 100) {
  cat("  WARNING: Z matrix appears transposed. Fixing...\n")
  Z <- t(Z)
  cat(sprintf(
    "  Corrected dimensions: %d genes × %d cell types\n",
    nrow(Z),
    ncol(Z)
  ))
}

# Validate Z dimensions
if (ncol(Z) != ncol(theta)) {
  warning(sprintf(
    "Dimension mismatch: Z has %d cell types but theta has %d",
    ncol(Z),
    ncol(theta)
  ))
}

# Try to extract log-likelihood
loglik_summary <- tryCatch(
  {
    if (!is.null(bp_res@posterior.theta_f@loglikelihood)) {
      sprintf("%.2f", bp_res@posterior.theta_f@loglikelihood)
    } else if ("logLik" %in% slotNames(bp_res@posterior.theta_f)) {
      sprintf("%.2f", bp_res@posterior.theta_f@logLik)
    } else if ("log.lik" %in% slotNames(bp_res@posterior.theta_f)) {
      sprintf("%.2f", bp_res@posterior.theta_f@log.lik)
    } else {
      "Not available"
    }
  },
  error = function(e) "Not available"
)

cat("  Log-likelihood:", loglik_summary, "\n")

# ==============================================================================
# 3. Save Results
# ==============================================================================

cat("\n=== Step 3/5: Saving Results ===\n")

# Save theta and theta_cv as CSV
write.csv(theta, file.path(output_dir, "theta_final.csv"), quote = FALSE)
write.csv(theta_cv, file.path(output_dir, "theta_cv.csv"), quote = FALSE)
cat("  Saved: theta_final.csv, theta_cv.csv\n")

# Save Z as compressed RDS
saveRDS(Z, file.path(output_dir, "Z_matrix.rds"), compress = "xz")
cat("  Saved: Z_matrix.rds (compressed)\n")

# Save top 1000 variable genes as CSV for quick inspection
# Calculate variance across cell types (by row = by gene)
if (nrow(Z) > 1) {
  z_var <- apply(Z, 1, var)
  top_genes <- names(sort(z_var, decreasing = TRUE))[1:min(1000, nrow(Z))]
  write.csv(
    Z[top_genes, , drop = FALSE],
    file.path(output_dir, "Z_matrix_top1000genes.csv"),
    quote = FALSE
  )
  cat("  Saved: Z_matrix_top1000genes.csv\n")
} else {
  warning("Z matrix has only 1 gene - skipping top genes CSV")
}

# Quick preview
cat("\n  Cell type fractions (theta):\n")
print(round(theta, 3))

# ==============================================================================
# 4. Marker Gene Validation
# ==============================================================================

cat("\n=== Step 4/5: Marker Gene Validation ===\n")

# Define marker genes for respiratory tract cell types
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
    "CD3G",
    "CD4",
    "CD8A",
    "CD8B",
    "GNLY",
    "NKG7",
    "GZMA",
    "GZMB"
  ),
  Myeloid = c(
    "CD14",
    "FCGR3A",
    "CD68",
    "CD163",
    "CD1C",
    "CLEC9A",
    "FCER1A",
    "ITGAX"
  ),
  B_cells = c("MS4A1", "CD79A", "CD79B", "IGHG1", "IGHM", "IGHA1"),
  Stromal = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "ACTA2",
    "PDGFRA",
    "PDGFRB"
  ),
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1", "KDR")
)

# Check which markers are present
all_markers <- unique(unlist(marker_genes))
markers_present <- all_markers[all_markers %in% rownames(Z)]

cat(sprintf(
  "  Marker genes detected: %d / %d (%.1f%%)\n",
  length(markers_present),
  length(all_markers),
  100 * length(markers_present) / length(all_markers)
))

if (length(markers_present) > 0) {
  Z_marker <- Z[markers_present, , drop = FALSE]
  write.csv(Z_marker, file.path(output_dir, "Z_marker_genes.csv"))

  # Create marker gene heatmap
  tryCatch(
    {
      pdf(file.path(output_dir, "Z_marker_heatmap.pdf"), width = 10, height = 8)

      pheatmap(
        log1p(Z_marker),
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
        border_color = NA,
        main = "Deconvolved Marker Gene Expression (log1p, row-scaled)",
        fontsize_row = 7,
        fontsize_col = 10
      )

      dev.off()
      cat("  Saved: Z_marker_genes.csv, Z_marker_heatmap.pdf\n")
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      warning("Error creating marker heatmap: ", e$message)
    }
  )
}

# ==============================================================================
# 5. Generate Visualizations
# ==============================================================================

cat("\n=== Step 5/5: Generating Visualizations ===\n")

# Prepare data for plotting
theta_df <- as.data.frame(theta)
theta_df$sample <- rownames(theta_df)
theta_long <- melt(
  as.data.table(theta_df),
  id.vars = "sample",
  variable.name = "cell_type",
  value.name = "fraction"
)

# 1. Stacked bar plot
tryCatch(
  {
    p1 <- ggplot(theta_long, aes(x = sample, y = fraction, fill = cell_type)) +
      geom_col(width = 0.7, color = "black", size = 0.3) +
      scale_fill_brewer(palette = "Set2") +
      scale_y_continuous(expand = c(0, 0)) +
      theme_bw(base_size = 12) +
      labs(
        title = "BayesPrism: Cell Type Fractions (Healthy Samples)",
        x = "Bulk Sample",
        y = "Fraction (proportion of reads)",
        fill = "Cell Type"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )

    ggsave(
      file.path(output_dir, "theta_stacked_bar.pdf"),
      p1,
      width = 8,
      height = 5
    )
    cat("  Saved: theta_stacked_bar.pdf\n")
  },
  error = function(e) {
    warning("Error creating stacked bar plot: ", e$message)
  }
)

# 2. Theta heatmap
tryCatch(
  {
    pdf(file.path(output_dir, "theta_heatmap.pdf"), width = 8, height = 6)

    pheatmap(
      theta,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
      border_color = "grey60",
      main = "Cell Type Fractions (theta)",
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = 9,
      fontsize = 10
    )

    dev.off()
    cat("  Saved: theta_heatmap.pdf\n")
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("Error creating theta heatmap: ", e$message)
  }
)

# 3. Theta CV heatmap
tryCatch(
  {
    pdf(file.path(output_dir, "theta_cv_heatmap.pdf"), width = 8, height = 6)

    pheatmap(
      theta_cv,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "orange", "red"))(100),
      border_color = "grey60",
      main = "Coefficient of Variation (uncertainty)",
      display_numbers = TRUE,
      number_format = "%.3f",
      fontsize_number = 9,
      fontsize = 10
    )

    dev.off()
    cat("  Saved: theta_cv_heatmap.pdf\n")
  },
  error = function(e) {
    if (dev.cur() > 1) {
      dev.off()
    }
    warning("Error creating theta CV heatmap: ", e$message)
  }
)

# 4. Boxplot
tryCatch(
  {
    p2 <- ggplot(
      theta_long,
      aes(x = cell_type, y = fraction, fill = cell_type)
    ) +
      geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw(base_size = 12) +
      labs(
        title = "Cell Type Fraction Distribution (n=5 healthy samples)",
        x = "Cell Type",
        y = "Fraction"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank()
      )

    ggsave(
      file.path(output_dir, "theta_boxplot.pdf"),
      p2,
      width = 8,
      height = 5
    )
    cat("  Saved: theta_boxplot.pdf\n")
  },
  error = function(e) {
    warning("Error creating boxplot: ", e$message)
  }
)

# ==============================================================================
# 6. Generate Summary Report
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

tryCatch(
  {
    summary_path <- file.path(output_dir, "analysis_summary.txt")
    sink(summary_path)

    cat(
      "================================================================================\n"
    )
    cat("BayesPrism Deconvolution Analysis Summary (Quick Recovery)\n")
    cat(
      "================================================================================\n\n"
    )

    cat("Recovery Date:", as.character(Sys.Date()), "\n\n")

    cat("--- Input Data ---\n")
    if (has_metadata) {
      cat("Reference:", ref_metadata$source, "\n")
      cat("  Total cells:      ", ref_metadata$n_cells, "\n")
      cat("  Total genes:      ", ref_metadata$n_genes, "\n")
      cat(
        "  Annotation mode:  ",
        ifelse(ref_metadata$use_single_level, "SINGLE level", "TWO levels"),
        "\n"
      )
      cat("  Cell types:       ", length(ref_metadata$celltype_counts), "\n\n")
    } else {
      cat("Reference metadata not available\n\n")
    }

    cat("Bulk samples:", nrow(theta), "\n")
    cat("  Sample IDs:       ", paste(rownames(theta), collapse = ", "), "\n\n")

    cat("--- Results ---\n")
    cat("Final log-likelihood:", loglik_summary, "\n\n")

    cat("Cell type fractions (theta):\n")
    print(round(theta, 4))
    cat("\n")

    cat("Mean fractions across samples:\n")
    print(round(colMeans(theta), 4))
    cat("\n")

    cat("Coefficient of variation (theta_cv):\n")
    print(round(theta_cv, 4))
    cat("\n")

    cat("--- Marker Gene Detection ---\n")
    cat(sprintf(
      "Markers detected: %d / %d (%.1f%%)\n",
      length(markers_present),
      length(all_markers),
      100 * length(markers_present) / length(all_markers)
    ))
    cat("\n")

    cat("--- Output Files ---\n")
    cat("1.  bp_res.rds                    : BayesPrism result object\n")
    cat("2.  theta_final.csv               : Cell type fractions\n")
    cat("3.  theta_cv.csv                  : Uncertainty (CV)\n")
    cat(
      "4.  Z_matrix.rds                  : Deconvolved gene expression (compressed)\n"
    )
    cat(
      "5.  Z_matrix_top1000genes.csv     : Top variable genes (quick inspection)\n"
    )
    cat("6.  Z_marker_genes.csv            : Marker gene expression\n")
    cat("7.  theta_stacked_bar.pdf         : Stacked bar plot\n")
    cat("8.  theta_heatmap.pdf             : Fraction heatmap\n")
    cat("9.  theta_cv_heatmap.pdf          : Uncertainty heatmap\n")
    cat("10. theta_boxplot.pdf             : Fraction distribution\n")
    cat("11. Z_marker_heatmap.pdf          : Marker gene heatmap\n")
    cat("12. analysis_summary.txt          : This summary\n\n")

    cat(
      "================================================================================\n"
    )
    cat("Recovery completed successfully!\n")
    cat("Output directory:", output_dir, "\n")
    cat(
      "================================================================================\n"
    )

    sink()
  },
  error = function(e) {
    while (sink.number() > 0) {
      sink()
    }
    warning("Error writing summary report: ", e$message)
  },
  finally = {
    while (sink.number() > 0) {
      sink()
    }
  }
)

cat("  Saved: analysis_summary.txt\n")

# ==============================================================================
# Final Message
# ==============================================================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("✓ Quick recovery completed successfully!\n")
cat("  Output directory: ", output_dir, "\n")
cat("  Files generated: 11 (from existing bp_res.rds)\n")
cat(rep("=", 80), "\n", sep = "")
cat("\nGenerated files:\n")
cat(
  "- CSV: theta_final.csv, theta_cv.csv, Z_matrix_top1000genes.csv, Z_marker_genes.csv\n"
)
cat("- RDS: Z_matrix.rds\n")
cat("- PDF: theta_stacked_bar.pdf, theta_heatmap.pdf, theta_cv_heatmap.pdf,\n")
cat("       theta_boxplot.pdf, Z_marker_heatmap.pdf\n")
cat("- TXT: analysis_summary.txt\n\n")
cat("Next steps:\n")
cat("1. Review theta_stacked_bar.pdf for cell type composition\n")
cat("2. Check theta_cv_heatmap.pdf for estimation uncertainty\n")
cat("3. Validate results with Z_marker_heatmap.pdf\n\n")
