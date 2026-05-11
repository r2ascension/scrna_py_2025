#!/usr/bin/env Rscript

# ============================================================
# GEO scRNA-seq Complete Processing Pipeline (with EmptyDrops)
# Author: r2end
# Date: 2024-12-11
#
# Features:
# - Dynamic GSE/GSM/HRA/HRR detection
# - Auto-detect data format (H5/10x/matrix)
# - Handle non-standard 10x file naming
# - Smart matrix header detection (no regex tricks)
# - Filter out ADT data
# - EmptyDrops filtering per sample
# - QC metrics (no additional filtering)
# - Memory-efficient processing
# - Individual Seurat objects per sample
# ============================================================

# ===== Configuration Section =====
DATA_ROOT <- "/home/h2048/data/source/1210"
METADATA_PATH <- "/home/h2048/data/source/reference/metadata_full_atlas_20251210.csv"
OUTPUT_DIR <- "/home/h2048/data/R/1210/Polyp"
QC_PLOTS_DIR <- "/home/h2048/data/R/1210/Polyp/qc_plots"

# EmptyDrops parameters
RUN_EMPTYDROPS <- TRUE
EMPTYDROPS_LOWER <- 100
EMPTYDROPS_FDR_THRESHOLD <- 0.01

# QC plot settings
MAX_PLOT_POINTS <- 5000

# QC thresholds (for visualization reference lines only)
QC_THRESHOLDS <- list(
  nFeature_min = 200,
  nFeature_max = 6000,
  nCount_min = 500,
  nCount_max = 50000,
  percent_mt_max = 20,
  percent_rb_min = 5
)

# ===== Load Libraries =====
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(hdf5r)
  library(ggplot2)
  library(patchwork)
  library(DropletUtils)
})

cat("============================================================\n")
cat("GEO/HRA scRNA-seq Complete Processing Pipeline\n")
cat("============================================================\n")
cat("Data root:  ", DATA_ROOT, "\n")
cat("Metadata:   ", METADATA_PATH, "\n")
cat("Output:     ", OUTPUT_DIR, "\n")
cat("QC plots:   ", QC_PLOTS_DIR, "\n")
cat("EmptyDrops: ", ifelse(RUN_EMPTYDROPS, "ENABLED", "DISABLED"), "\n")
if (RUN_EMPTYDROPS) {
  cat("  Lower bound:", EMPTYDROPS_LOWER, "UMIs\n")
  cat("  FDR threshold:", EMPTYDROPS_FDR_THRESHOLD, "\n")
}
cat("============================================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ===== Load Metadata =====
cat("Loading metadata...\n")

if (!file.exists(METADATA_PATH)) {
  stop("Metadata file not found: ", METADATA_PATH)
}

metadata <- fread(METADATA_PATH, header = TRUE)

cat("  Original columns:", paste(names(metadata), collapse = ", "), "\n")

# ===== Standardize Metadata Column Names =====
cat("Standardizing metadata columns...\n")

cols <- names(metadata)

# Fix typo: filtered_or_Falset -> filtered_or_raw
if ("filtered_or_Falset" %in% cols && !"filtered_or_raw" %in% cols) {
  setnames(metadata, "filtered_or_Falset", "filtered_or_raw")
  cat("  Renamed 'filtered_or_Falset' to 'filtered_or_raw'\n")
}

# Fix typo: reference geFalseme -> reference_genome
ref_cols <- grep("^reference", cols, value = TRUE, ignore.case = TRUE)
if (length(ref_cols) == 1 && !"reference_genome" %in% cols) {
  setnames(metadata, ref_cols, "reference_genome")
  cat("  Renamed '", ref_cols, "' to 'reference_genome'\n", sep = "")
}

# ===== Process Metadata Columns =====
cat("Processing metadata columns...\n")

# Create sample_id column from 'sample' for indexing
if ("sample" %in% names(metadata)) {
  metadata[, sample_id := sample]
  cat("  Created 'sample_id' column from 'sample' (sample column preserved)\n")
} else {
  stop("'sample' column not found in metadata")
}

# Clean sample_id column (remove spaces and tabs)
metadata[, sample_id := gsub("[\\s\\t]+", "", sample_id)]

# Check for duplicates
if (any(duplicated(metadata$sample_id))) {
  warning(
    "Duplicated sample_ids found in metadata; only first occurrence will be used."
  )
}

# Set key for fast lookup
setkey(metadata, sample_id)

cat("  Total metadata records:", nrow(metadata), "\n")
cat("  Unique sample_ids:", uniqueN(metadata$sample_id), "\n")
cat("  Final columns:", paste(names(metadata), collapse = ", "), "\n\n")

# ===== Dynamic GSE/HRA Detection =====
cat("Scanning for GSE/HRA directories...\n")

all_dirs <- list.dirs(DATA_ROOT, recursive = FALSE, full.names = FALSE)
DATASET_IDS <- all_dirs[grepl("^(GSE|HRA)\\d+$", all_dirs)]

if (length(DATASET_IDS) == 0) {
  stop("No GSE/HRA directories found in: ", DATA_ROOT)
}

cat("  Found", length(DATASET_IDS), "datasets:\n")
cat("  ", paste(DATASET_IDS, collapse = ", "), "\n\n")

# ===== Helper Functions =====

# Read 10x files manually (handles non-standard naming)
read_10x_manual <- function(gsm_dir, sample_prefix = NULL) {
  files <- list.files(gsm_dir, full.names = FALSE)

  # Find barcodes file
  barcode_pattern <- if (is.null(sample_prefix)) {
    "barcodes\\.tsv"
  } else {
    paste0(sample_prefix, ".*barcodes\\.tsv")
  }
  barcode_file <- grep(
    barcode_pattern,
    files,
    value = TRUE,
    ignore.case = TRUE
  )[1]

  # Find features/genes file
  feature_pattern <- if (is.null(sample_prefix)) {
    "(features|genes)\\.tsv"
  } else {
    paste0(sample_prefix, ".*(features|genes)\\.tsv")
  }
  feature_file <- grep(
    feature_pattern,
    files,
    value = TRUE,
    ignore.case = TRUE
  )[1]

  # Find matrix file
  matrix_pattern <- if (is.null(sample_prefix)) {
    "matrix\\.mtx"
  } else {
    paste0(sample_prefix, ".*matrix\\.mtx")
  }
  matrix_file <- grep(matrix_pattern, files, value = TRUE, ignore.case = TRUE)[
    1
  ]

  if (is.na(barcode_file) || is.na(feature_file) || is.na(matrix_file)) {
    return(NULL)
  }

  barcode_path <- file.path(gsm_dir, barcode_file)
  feature_path <- file.path(gsm_dir, feature_file)
  matrix_path <- file.path(gsm_dir, matrix_file)

  # Read files
  barcodes <- readLines(barcode_path)
  features <- read.table(
    feature_path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  mat <- readMM(matrix_path)

  # Set dimensions
  colnames(mat) <- barcodes
  rownames(mat) <- features[, 1] # Use first column as gene names

  return(mat)
}

# Detect data format in sample directory
detect_data_format <- function(sample_dir, sample_id) {
  files <- list.files(sample_dir, full.names = FALSE)

  # Check for H5 format
  h5_files <- grep("\\.(h5|hdf5)$", files, value = TRUE, ignore.case = TRUE)
  if (length(h5_files) > 0) {
    return(list(
      format = "h5",
      path = file.path(sample_dir, h5_files[1])
    ))
  }

  # Check for 10x format (barcodes + features/genes + matrix)
  has_barcodes <- any(grepl("barcodes\\.tsv", files, ignore.case = TRUE))
  has_features <- any(grepl(
    "(features|genes)\\.tsv",
    files,
    ignore.case = TRUE
  ))
  has_matrix <- any(grepl("matrix\\.mtx", files, ignore.case = TRUE))

  if (has_barcodes && has_features && has_matrix) {
    # Check if files have standard names
    standard_barcodes <- any(files %chin% c("barcodes.tsv.gz", "barcodes.tsv"))

    return(list(
      format = "10x",
      path = sample_dir,
      has_prefix = !standard_barcodes,
      sample_prefix = if (!standard_barcodes) sample_id else NULL
    ))
  }

  # Check for matrix format
  matrix_files <- grep(
    "\\.(txt|tsv|csv|gz)$",
    files,
    value = TRUE,
    ignore.case = TRUE
  )
  matrix_files <- matrix_files[
    !grepl(
      "features|barcodes|genes|meta|anno",
      matrix_files,
      ignore.case = TRUE
    )
  ]

  if (length(matrix_files) > 0) {
    prio_idx <- grep(
      "expr|expression|count|matrix",
      matrix_files,
      ignore.case = TRUE
    )
    if (length(prio_idx) > 0) {
      chosen <- matrix_files[prio_idx[1]]
    } else {
      chosen <- matrix_files[1]
    }

    return(list(
      format = "matrix",
      path = file.path(sample_dir, chosen)
    ))
  }

  return(list(format = "unknown", path = NULL))
}
#' Load counts data and filter ADT
load_counts_data <- function(format_info, sample_id) {
  format_type <- format_info$format
  path <- format_info$path

  tryCatch(
    {
      if (format_type == "h5") {
        cat("    Loading H5 format...\n")
        counts_list <- Read10X_h5(path, use.names = TRUE)

        if (is.list(counts_list) && !inherits(counts_list, "Matrix")) {
          cat("    Multiple assays detected:", names(counts_list), "\n")
          if ("Gene Expression" %in% names(counts_list)) {
            counts <- counts_list[["Gene Expression"]]
            cat("    Using 'Gene Expression' assay\n")
          } else {
            counts <- counts_list[[1]]
            cat("    Using first assay:", names(counts_list)[1], "\n")
          }
        } else {
          counts <- counts_list
        }
      } else if (format_type == "10x") {
        cat("    Loading 10x format...\n")

        # Try standard Read10X first
        if (!format_info$has_prefix) {
          counts_list <- Read10X(path, gene.column = 1)

          if (is.list(counts_list) && !inherits(counts_list, "Matrix")) {
            cat("    Multiple assays detected:", names(counts_list), "\n")
            if ("Gene Expression" %in% names(counts_list)) {
              counts <- counts_list[["Gene Expression"]]
              cat("    Using 'Gene Expression' assay (ADT filtered)\n")
            } else {
              counts <- counts_list[[1]]
              cat("    Using first assay:", names(counts_list)[1], "\n")
            }
          } else {
            counts <- counts_list
          }
        } else {
          # Manual reading for non-standard naming
          cat("    Detected non-standard file naming, reading manually...\n")
          counts <- read_10x_manual(path, format_info$sample_prefix)

          if (is.null(counts)) {
            stop("Failed to read 10x files manually")
          }
        }
      } else if (format_type == "matrix") {
        cat("    Loading matrix format...\n")

        # ---- 1. 先用首行判断分隔符和初步 header ----
        if (grepl("\\.gz$", path)) {
          con <- gzfile(path, "rt")
        } else {
          con <- file(path, "rt")
        }
        first_line <- readLines(con, n = 1)
        close(con)

        delim <- if (grepl("\t", first_line)) "\t" else ","

        tokens <- strsplit(first_line, split = delim, fixed = TRUE)[[1]]
        use_header <- FALSE
        if (length(tokens) > 1) {
          suppressWarnings({
            num_flags <- !is.na(as.numeric(tokens[-1]))
          })
          # 如果除了第一列之外全是数字，就当作「没有 header」
          use_header <- !(all(num_flags))
        }
        cat("    Header detection: use_header =", use_header, "\n")

        mat <- fread(path, sep = delim, header = use_header, data.table = FALSE)

        # ---- 2. 二次检查：如果用了 header，但列名还是数值，说明刚才误判了 ----
        if (use_header && ncol(mat) > 1) {
          cn <- colnames(mat)
          suppressWarnings({
            cn_numeric <- !is.na(as.numeric(cn[-1]))
          })
          if (all(cn_numeric)) {
            cat(
              "    Detected numeric column names after header read; re-reading as headerless matrix...\n"
            )
            mat <- fread(path, sep = delim, header = FALSE, data.table = FALSE)
            use_header <- FALSE
          }
        }

        # DEBUG: Print raw matrix info
        cat(
          "    Raw matrix dim (including first column):",
          nrow(mat),
          "x",
          ncol(mat),
          "\n"
        )
        cat("    Preview first column (expect gene names):\n")
        print(head(mat[[1]], 5))
        cat("    Preview column names (expect cell IDs or sample IDs):\n")
        print(head(colnames(mat), 5))

        # 第一列作为基因名
        genes <- mat[, 1]
        mat <- mat[, -1, drop = FALSE]

        # Sanity check: must have count columns
        if (ncol(mat) == 0) {
          stop(
            "Expression matrix has only one column (gene names). No count columns detected. Check delimiter/header in: ",
            basename(path)
          )
        }

        # 转成普通 matrix，然后加行名
        mat_numeric <- as.matrix(mat)
        rownames(mat_numeric) <- genes

        # 如果最终判定为「无 header」，给细胞生成合成 ID
        if (!use_header) {
          colnames(mat_numeric) <- paste0(
            sample_id,
            "_Cell",
            seq_len(ncol(mat_numeric))
          )
          cat(
            "    No valid cell header detected; generated synthetic cell IDs like:",
            colnames(mat_numeric)[1],
            ",",
            colnames(mat_numeric)[2],
            "...\n"
          )
        }

        cat("    After setting rownames:\n")
        cat(
          "      Dimensions:",
          nrow(mat_numeric),
          "rows x",
          ncol(mat_numeric),
          "cols\n"
        )
        cat(
          "      Example rownames:",
          paste(head(rownames(mat_numeric), 3), collapse = ", "),
          "\n"
        )
        cat(
          "      Example colnames:",
          paste(head(colnames(mat_numeric), 3), collapse = ", "),
          "\n"
        )

        # 转 sparse matrix
        counts <- as(mat_numeric, "sparseMatrix")

        cat("    Final counts matrix:\n")
        cat(
          "      Dimensions:",
          nrow(counts),
          "genes x",
          ncol(counts),
          "cells\n"
        )
        cat(
          "      Example gene names (rownames):",
          paste(head(rownames(counts), 3), collapse = ", "),
          "\n"
        )
        cat(
          "      Example cell IDs (colnames):",
          paste(head(colnames(counts), 3), collapse = ", "),
          "\n"
        )
      } else {
        stop("Unknown format")
      }

      if (!inherits(counts, "sparseMatrix")) {
        counts <- as(counts, "sparseMatrix")
      }

      cat("    Loaded:", nrow(counts), "genes x", ncol(counts), "cells\n")

      return(counts)
    },
    error = function(e) {
      cat("    ❌ Error loading", sample_id, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}

# Run EmptyDrops filtering
run_emptydrops_filter <- function(counts, sample_id) {
  cells_before <- ncol(counts)

  # Skip if too few cells
  if (cells_before < 100) {
    cat("    Skipping EmptyDrops (too few cells:", cells_before, ")\n")
    return(list(
      counts = counts,
      cells_before = cells_before,
      cells_after = cells_before,
      cells_removed = 0,
      method = "skipped_fewcells"
    ))
  }

  tryCatch(
    {
      cat("    Running EmptyDrops...\n")

      set.seed(42)
      e.out <- emptyDrops(counts, lower = EMPTYDROPS_LOWER)

      is.cell <- e.out$FDR <= EMPTYDROPS_FDR_THRESHOLD
      is.cell[is.na(is.cell)] <- FALSE

      high.count <- colSums(counts) > EMPTYDROPS_LOWER
      is.cell <- is.cell | (is.na(e.out$FDR) & high.count)

      counts_filtered <- counts[, is.cell, drop = FALSE]
      cells_after <- ncol(counts_filtered)
      cells_removed <- cells_before - cells_after

      cat("    EmptyDrops Results:\n")
      cat("      Cells before:", format(cells_before, big.mark = ","), "\n")
      cat("      Cells after:", format(cells_after, big.mark = ","), "\n")
      cat(
        "      Cells removed:",
        format(cells_removed, big.mark = ","),
        "(",
        round(cells_removed / cells_before * 100, 1),
        "%)\n"
      )

      return(list(
        counts = counts_filtered,
        cells_before = cells_before,
        cells_after = cells_after,
        cells_removed = cells_removed,
        method = "emptydrops"
      ))
    },
    error = function(e) {
      cat("    ⚠️  EmptyDrops failed:", conditionMessage(e), "\n")
      cat("    Proceeding without EmptyDrops filtering\n")
      return(list(
        counts = counts,
        cells_before = cells_before,
        cells_after = cells_before,
        cells_removed = 0,
        method = "failed"
      ))
    }
  )
}

# Extract orig.ident from path
extract_orig_ident <- function(path) {
  if (is.null(path) || length(path) == 0) {
    return("Unknown")
  }

  filename <- basename(path)
  orig_ident <- gsub("_features\\.tsv.*", "", filename, ignore.case = TRUE)
  orig_ident <- gsub("_genes\\.tsv.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_matrix\\.mtx.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_barcodes\\.tsv.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_filtered.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub("_raw.*", "", orig_ident, ignore.case = TRUE)
  orig_ident <- gsub(
    "\\.(h5|hdf5|txt|tsv|csv|gz)$",
    "",
    orig_ident,
    ignore.case = TRUE
  )

  return(orig_ident)
}

# Calculate QC metrics (human-specific)
calculate_qc_metrics <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-"
  )
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^RP[SL]"
  )
  return(seurat_obj)
}

# Generate QC plot
generate_qc_plot <- function(
  seurat_obj,
  sample_id,
  output_path,
  max_points = MAX_PLOT_POINTS
) {
  qc_data <- as.data.table(seurat_obj@meta.data)

  if (nrow(qc_data) > max_points) {
    sample_idx <- sample.int(nrow(qc_data), max_points)
    qc_data_sample <- qc_data[sample_idx]
  } else {
    qc_data_sample <- qc_data
  }

  p1 <- ggplot(qc_data, aes(x = "Sample", y = nFeature_RNA)) +
    geom_violin(fill = "#3498db", alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = nFeature_RNA),
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nFeature_min,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nFeature_max,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    labs(title = "Features per Cell", x = "", y = "nFeature_RNA") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p2 <- ggplot(qc_data, aes(x = "Sample", y = nCount_RNA)) +
    geom_violin(fill = "#2ecc71", alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = nCount_RNA),
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nCount_min,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$nCount_max,
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    scale_y_log10() +
    labs(title = "Counts per Cell", x = "", y = "nCount_RNA (log10)") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p3 <- ggplot(qc_data, aes(x = "Sample", y = percent.mt)) +
    geom_violin(fill = "#e74c3c", alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = percent.mt),
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$percent_mt_max,
      linetype = "dashed",
      color = "darkred",
      linewidth = 0.5
    ) +
    labs(title = "Mitochondrial %", x = "", y = "% MT genes") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p4 <- ggplot(qc_data, aes(x = "Sample", y = percent.rb)) +
    geom_violin(fill = "#9b59b6", alpha = 0.7) +
    geom_jitter(
      data = qc_data_sample,
      aes(x = "Sample", y = percent.rb),
      width = 0.2,
      size = 0.5,
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = QC_THRESHOLDS$percent_rb_min,
      linetype = "dashed",
      color = "darkviolet",
      linewidth = 0.5
    ) +
    labs(title = "Ribosomal %", x = "", y = "% RB genes") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  p5 <- ggplot(qc_data_sample, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", color = "red", linewidth = 0.8) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Feature vs Count",
      x = "nCount_RNA (log10)",
      y = "nFeature_RNA (log10)"
    ) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p6 <- ggplot(qc_data_sample, aes(x = nCount_RNA, y = percent.mt)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_x_log10() +
    labs(title = "MT% vs Count", x = "nCount_RNA (log10)", y = "% MT genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  layout <- "ABC\nDEF"

  combined <- p1 +
    p2 +
    p3 +
    p4 +
    p5 +
    p6 +
    plot_layout(design = layout) +
    plot_annotation(
      title = paste0("Quality Control Metrics: ", sample_id),
      subtitle = paste0(
        "Cells: ",
        format(ncol(seurat_obj), big.mark = ","),
        " | ",
        "Genes: ",
        format(nrow(seurat_obj), big.mark = ","),
        " | ",
        "Median nFeature: ",
        format(round(median(qc_data$nFeature_RNA)), big.mark = ","),
        " | ",
        "Median nCount: ",
        format(round(median(qc_data$nCount_RNA)), big.mark = ","),
        if (nrow(qc_data) > max_points) {
          paste0(
            " | Sampled ",
            format(max_points, big.mark = ","),
            " cells for scatter plots"
          )
        } else {
          ""
        }
      ),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    )

  ggsave(
    filename = output_path,
    plot = combined,
    width = 12,
    height = 8,
    dpi = 300
  )
  cat("    QC plot saved:", basename(output_path), "\n")
}

# ===== Main Processing Loop =====

processing_summary <- data.table(
  Dataset = character(),
  Sample_ID = character(),
  sample = character(),
  Format = character(),
  Cells_Raw = integer(),
  Cells_Filtered = integer(),
  Cells_Removed = integer(),
  EmptyDrops_Method = character(),
  Genes = integer(),
  Median_nFeature = numeric(),
  Median_nCount = numeric(),
  Median_MT = numeric(),
  Median_RB = numeric(),
  Status = character()
)

total_processed <- 0
total_failed <- 0
total_cells_raw <- 0
total_cells_filtered <- 0

for (dataset_id in DATASET_IDS) {
  dataset_dir <- file.path(DATA_ROOT, dataset_id)

  if (!dir.exists(dataset_dir)) {
    cat("⚠️  Skipping", dataset_id, "(directory not found)\n\n")
    next
  }

  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Processing Dataset:", dataset_id, "\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

  # Find all sample directories (GSM* or HRR*)
  sample_dirs <- list.dirs(dataset_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[grepl("^(GSM|HRR)\\d+$", basename(sample_dirs))]

  if (length(sample_dirs) == 0) {
    cat("⚠️  No sample directories found in", dataset_id, "\n\n")
    next
  }

  cat("Found", length(sample_dirs), "samples\n\n")

  for (sample_dir in sample_dirs) {
    sample_id <- basename(sample_dir)
    cat("📦", sample_id, "\n")

    # Detect data format
    format_info <- detect_data_format(sample_dir, sample_id)
    cat("  Format:", format_info$format, "\n")

    if (format_info$format == "unknown") {
      cat("  ⚠️  Skipping (unknown format)\n\n")
      processing_summary <- rbind(
        processing_summary,
        data.table(
          Dataset = dataset_id,
          Sample_ID = sample_id,
          sample = NA_character_,
          Format = "unknown",
          Cells_Raw = NA_integer_,
          Cells_Filtered = NA_integer_,
          Cells_Removed = NA_integer_,
          EmptyDrops_Method = "skipped",
          Genes = NA_integer_,
          Median_nFeature = NA_real_,
          Median_nCount = NA_real_,
          Median_MT = NA_real_,
          Median_RB = NA_real_,
          Status = "skipped"
        )
      )
      total_failed <- total_failed + 1
      next
    }

    # Load counts
    counts <- load_counts_data(format_info, sample_id)

    if (is.null(counts)) {
      processing_summary <- rbind(
        processing_summary,
        data.table(
          Dataset = dataset_id,
          Sample_ID = sample_id,
          sample = NA_character_,
          Format = format_info$format,
          Cells_Raw = NA_integer_,
          Cells_Filtered = NA_integer_,
          Cells_Removed = NA_integer_,
          EmptyDrops_Method = "failed",
          Genes = NA_integer_,
          Median_nFeature = NA_real_,
          Median_nCount = NA_real_,
          Median_MT = NA_real_,
          Median_RB = NA_real_,
          Status = "failed"
        )
      )
      total_failed <- total_failed + 1
      next
    }

    # Run EmptyDrops filtering
    if (RUN_EMPTYDROPS) {
      ed_result <- run_emptydrops_filter(counts, sample_id)
      counts <- ed_result$counts
      cells_raw <- ed_result$cells_before
      cells_filtered <- ed_result$cells_after
      cells_removed <- ed_result$cells_removed
      ed_method <- ed_result$method
    } else {
      cells_raw <- ncol(counts)
      cells_filtered <- ncol(counts)
      cells_removed <- 0
      ed_method <- "disabled"
    }

    # Check for zero cells or genes after EmptyDrops
    if (ncol(counts) == 0 || nrow(counts) == 0) {
      cat(
        "  ⚠️  No valid cells or genes detected (",
        nrow(counts),
        " genes x ",
        ncol(counts),
        " cells). Skipping sample.\n\n",
        sep = ""
      )

      sample_val <- if (sample_id %in% metadata$sample_id) {
        metadata[sample_id]$sample
      } else {
        NA_character_
      }

      processing_summary <- rbind(
        processing_summary,
        data.table(
          Dataset = dataset_id,
          Sample_ID = sample_id,
          sample = sample_val,
          Format = format_info$format,
          Cells_Raw = cells_raw,
          Cells_Filtered = 0,
          Cells_Removed = cells_raw,
          EmptyDrops_Method = if (RUN_EMPTYDROPS) {
            "failed_nocells"
          } else {
            "disabled_nocells"
          },
          Genes = nrow(counts),
          Median_nFeature = NA_real_,
          Median_nCount = NA_real_,
          Median_MT = NA_real_,
          Median_RB = NA_real_,
          Status = "failed_nocells"
        )
      )
      total_failed <- total_failed + 1
      next
    }

    # Create Seurat object
    orig_ident <- extract_orig_ident(format_info$path)

    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = orig_ident,
      min.cells = 0,
      min.features = 0
    )

    cat(
      "  Created Seurat object:",
      format(ncol(seurat_obj), big.mark = ","),
      "cells x",
      format(nrow(seurat_obj), big.mark = ","),
      "genes\n"
    )

    # Add basic metadata
    seurat_obj$sample_id <- sample_id
    seurat_obj$dataset <- dataset_id
    seurat_obj$orig.ident <- orig_ident
    seurat_obj$emptydrops_filtered <- RUN_EMPTYDROPS &&
      ed_method == "emptydrops"

    # Add metadata from CSV
    if (sample_id %in% metadata$sample_id) {
      meta_row <- metadata[sample_id]

      seurat_obj$sample <- meta_row$sample
      seurat_obj$GEO <- meta_row$GEO
      seurat_obj$tissue <- meta_row$tissue
      seurat_obj$tissue_level_2 <- meta_row$tissue_level_2
      seurat_obj$disease_level_1 <- meta_row$disease_level_1
      seurat_obj$disease_level_2 <- meta_row$disease_level_2
      seurat_obj$condition <- meta_row$condition
      seurat_obj$tissue_sampling_method <- meta_row$tissue_sampling_method
      seurat_obj$filtered_or_raw <- meta_row$filtered_or_raw
      seurat_obj$target_cell <- meta_row$target_cell
      seurat_obj$frozen_or_fresh <- meta_row$frozen_or_fresh
      seurat_obj$sex <- meta_row$sex
      seurat_obj$age <- meta_row$age
      seurat_obj$platform <- meta_row$platform
      seurat_obj$reference_genome <- meta_row$reference_genome
      seurat_obj$assay <- meta_row$assay
      seurat_obj$tissue_dissociation_protocol <- meta_row$tissue_dissociation_protocol

      cat("  ✓ Metadata added\n")
    } else {
      cat("  ⚠️  No metadata found for", sample_id, "\n")
    }

    # Calculate QC metrics
    cat("  Calculating QC metrics...\n")
    seurat_obj <- calculate_qc_metrics(seurat_obj)

    # Get QC stats
    qc_stats <- seurat_obj@meta.data
    median_nfeature <- median(qc_stats$nFeature_RNA)
    median_ncount <- median(qc_stats$nCount_RNA)
    median_mt <- median(qc_stats$percent.mt)
    median_rb <- median(qc_stats$percent.rb)

    cat("  QC Summary:\n")
    cat(
      "    Median nFeature:",
      format(round(median_nfeature), big.mark = ","),
      "\n"
    )
    cat(
      "    Median nCount:",
      format(round(median_ncount), big.mark = ","),
      "\n"
    )
    cat("    Median MT%:", round(median_mt, 2), "\n")
    cat("    Median RB%:", round(median_rb, 2), "\n")

    # Generate QC plot
    qc_plot_path <- file.path(QC_PLOTS_DIR, paste0(sample_id, "_qc.png"))
    generate_qc_plot(seurat_obj, sample_id, qc_plot_path)

    # Save Seurat object
    output_file <- file.path(OUTPUT_DIR, paste0(sample_id, "_seurat.rds"))
    saveRDS(seurat_obj, output_file, compress = "gzip")
    cat("  ✓ Seurat object saved:", basename(output_file), "\n")

    # Get sample value for summary
    sample_val <- if (sample_id %in% metadata$sample_id) {
      metadata[sample_id]$sample
    } else {
      NA_character_
    }

    # Update summary
    processing_summary <- rbind(
      processing_summary,
      data.table(
        Dataset = dataset_id,
        Sample_ID = sample_id,
        sample = sample_val,
        Format = format_info$format,
        Cells_Raw = cells_raw,
        Cells_Filtered = cells_filtered,
        Cells_Removed = cells_removed,
        EmptyDrops_Method = ed_method,
        Genes = nrow(seurat_obj),
        Median_nFeature = median_nfeature,
        Median_nCount = median_ncount,
        Median_MT = median_mt,
        Median_RB = median_rb,
        Status = "success"
      )
    )

    total_processed <- total_processed + 1
    total_cells_raw <- total_cells_raw + cells_raw
    total_cells_filtered <- total_cells_filtered + cells_filtered

    # Memory cleanup
    rm(counts, seurat_obj, qc_stats)
    gc(verbose = FALSE)

    cat("\n")
  }
}

# ===== Final Summary =====

cat("============================================================\n")
cat("Processing Complete\n")
cat("============================================================\n\n")

cat("Summary Statistics:\n")
cat("  Total samples processed:", total_processed, "\n")
cat("  Total samples failed:", total_failed, "\n")
cat("  Total cells (raw):", format(total_cells_raw, big.mark = ","), "\n")
cat(
  "  Total cells (filtered):",
  format(total_cells_filtered, big.mark = ","),
  "\n"
)
if (RUN_EMPTYDROPS) {
  cells_removed_total <- total_cells_raw - total_cells_filtered
  cat(
    "  Cells removed by EmptyDrops:",
    format(cells_removed_total, big.mark = ","),
    "(",
    round(cells_removed_total / total_cells_raw * 100, 1),
    "%)\n"
  )
}

if (all(is.na(processing_summary$Genes))) {
  max_genes <- NA_integer_
} else {
  max_genes <- max(processing_summary$Genes, na.rm = TRUE)
}
cat("  Total genes (max):", format(max_genes, big.mark = ","), "\n\n")

cat("Format Distribution:\n")
print(table(processing_summary$Format))
cat("\n")

cat("Status Distribution:\n")
print(table(processing_summary$Status))
cat("\n")

if (RUN_EMPTYDROPS) {
  cat("EmptyDrops Method Distribution:\n")
  print(table(processing_summary$EmptyDrops_Method))
  cat("\n")
}

# Save processing summary
summary_file <- file.path(OUTPUT_DIR, "processing_summary.csv")
fwrite(processing_summary, summary_file)
cat("Processing summary saved:", summary_file, "\n\n")

# Generate summary statistics by Dataset
cat("Per-Dataset Statistics:\n")
dataset_summary <- processing_summary[
  Status == "success",
  .(
    Samples = .N,
    Total_Cells_Raw = sum(Cells_Raw),
    Total_Cells_Filtered = sum(Cells_Filtered),
    Pct_Removed = round(
      (sum(Cells_Raw) - sum(Cells_Filtered)) / sum(Cells_Raw) * 100,
      1
    ),
    Mean_Cells = round(mean(Cells_Filtered)),
    Median_nFeature = round(median(Median_nFeature)),
    Median_nCount = round(median(Median_nCount)),
    Median_MT = round(median(Median_MT), 2),
    Median_RB = round(median(Median_RB), 2)
  ),
  by = Dataset
]

print(dataset_summary)

dataset_summary_file <- file.path(OUTPUT_DIR, "dataset_summary.csv")
fwrite(dataset_summary, dataset_summary_file)
cat("\nDataset summary saved:", dataset_summary_file, "\n\n")

cat("============================================================\n")
cat("✨ Pipeline Complete!\n")
cat("============================================================\n")
cat("Output locations:\n")
cat("  Seurat objects:", OUTPUT_DIR, "\n")
cat("  QC plots:", QC_PLOTS_DIR, "\n")
cat("  Summary files:", OUTPUT_DIR, "\n")
cat("============================================================\n")

# Print file counts
rds_count <- length(list.files(OUTPUT_DIR, pattern = "_seurat\\.rds$"))
qc_count <- length(list.files(QC_PLOTS_DIR, pattern = "_qc\\.png$"))

cat("\nGenerated files:\n")
cat("  RDS files:", rds_count, "\n")
cat("  QC plots:", qc_count, "\n")

if (rds_count != qc_count || rds_count != total_processed) {
  cat("\n⚠️  Warning: File count mismatch detected\n")
  cat("  Expected:", total_processed, "files per type\n")
  cat("  Check processing_summary.csv for details\n")
}

cat("\n============================================================\n")
