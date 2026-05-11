#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Subcluster LLM Interpretation Analysis v1.3.1-HOTFIX
# ==============================================================================
#
# Author: r2end
# Date: 2026-01-26
# Version: v1.3.1-HOTFIX (Production-grade fixes applied)
#
# CRITICAL FIXES IN v1.3.1-HOTFIX:
#   ✅ P0-1: Fixed per-celltype subcluster selection (uses L2→L3 mapping, not string prefix)
#   ✅ P0-2: Initialize LLM result objects to prevent undefined errors
#   ✅ P0-3: Auto-select future plan based on OS + limit BLAS threads
#   ✅ P0-4: Auto-detect logFC column name (avg_log2FC vs avg_logFC)
#   ✅ P1-2: Optimized parallel marker discovery memory usage
#   ✅ P1-3: Clean CellMarker marker field (split comma-separated genes)
#   ✅ P1-4: Filter invalid gene symbols before GO enrichment
#   ✅ P2-1: Save sessionInfo and QC checkpoint
#   ✅ P2-2: Force factor conversion for clustering columns
#   ✅ P2-3: Handle empty compareClusterResult safely
#
# SKIPPED (per user request):
#   - P1-5: LLM input limiting (controlled by package)
#   - P1-1: Normalize/Scale logic review (analysis-specific)
#
# Key Features:
#   - Uses clusterProfiler package's interpret() functions
#   - Parallel FindAllMarkers (10-20x faster)
#   - Multi-database enrichment (CellMarker + GO)
#   - Three LLM tasks: annotation, phenotyping, interpretation
#   - Memory-efficient processing
#   - Production-grade error handling
#   - Cross-platform compatible (Windows/Linux/Mac)
#
# Requirements:
#   - Custom clusterProfiler package with interpret() functions
#   - fanyi package for LLM API calls
#   - All standard Seurat/enrichment dependencies
#
# Runtime: ~20-30 minutes
# Memory: < 40GB RAM
#
# ==============================================================================

# Input:
#   - adata_bcell_subclustered_FINAL_v2_20260119.h5ad (from Python pipeline)
#     Required fields:
#       * cell_type_L2: scANVI predictions (Memory B, Naive B, Plasma, etc.)
#       * cell_type_L3: Subclusters (Memory_B_c0, Memory_B_c1, etc.)
#       * subcluster_id: Cluster IDs (0, 1, 2, ...)
#     Optional fields:
#       * percent.mt: Mitochondrial percentage (if present, will be used for QC)
#
# Output:
#   - Comprehensive LLM interpretation report
#   - Cell subtype annotations with confidence scores
#   - Functional phenotype characterizations
#   - Publication-ready mechanistic narratives
#   - High-quality visualizations
#   - sessionInfo.txt for reproducibility
#   - seurat_post_qc.rds checkpoint
#
# ==============================================================================
#
# Key Features:
#   - Uses GetSeurat() for seamless Python-R integration
#   - Parallel FindAllMarkers (10-20x faster)
#   - Multi-database enrichment (CellMarker + GO)
#   - Three LLM tasks: annotation, phenotyping, interpretation
#   - Memory-efficient processing
#   - Production-grade error handling
#   - Uses clusterProfiler package's built-in interpret() functions
#   - Secure API key management via environment variables
#
# Changes in v1.3-PACKAGE:
#   - Removed custom interpret() implementation
#   - Uses clusterProfiler::interpret() from package
#   - Uses clusterProfiler::interpret_hierarchical()
#   - Uses clusterProfiler::interpret_agent()
#   - Relies on package's call_llm_fanyi() internal function
#   - Simplified code by delegating to package
#
# Requirements:
#   - Custom clusterProfiler package with interpret() functions installed
#   - fanyi package for LLM API calls
#   - All other standard dependencies
#
# Runtime: ~20-30 minutes (depends on API speed)
# Memory: < 40GB RAM
#
# ==============================================================================

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# ----- Input/Output Paths -----
H5AD_PATH <- "/home/h2048/data/py/0119/bcell_analysis/results/subcluster_v2_20260119/adata_bcell_subclustered_FINAL_v2_20260119.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0126/bcell_interpret_v1"

# ----- Cell Type/Subcluster Columns -----
CELLTYPE_L2_COL <- "cell_type_L2" # scANVI level (Memory B, Naive B, etc.)
CELLTYPE_L3_COL <- "cell_type_L3" # With subclusters (Memory_B_c0, etc.)
SUBCLUSTER_ID_COL <- "subcluster_id" # Pure cluster ID (0, 1, 2, ...)

# ----- Analysis Parameters -----
MIN_CELLS_FOR_ANALYSIS <- 50 # Minimum cells for reliable enrichment
TOP_N_MARKERS <- 50 # Top N markers per subcluster for enrichment
MIN_CELLS_FOR_INDIVIDUAL_ANALYSIS <- 200 # Min cells to analyze celltype separately

# ----- Marker Gene Settings -----
MARKER_METHOD <- "wilcox" # "wilcox" (fast, Seurat v5), "t" (fast), "MAST" (slow)
MARKER_LOGFC_THRESHOLD <- 0.5
MARKER_PADJ_THRESHOLD <- 0.05
MIN_PCT <- 0.25

# ----- Parallel Processing -----
N_CORES <- 8 # Adjust based on your system
FUTURE_MAX_SIZE <- 10 * 1024^3 # 10GB per worker

# ----- LLM API Settings -----
USE_DEEPSEEK <- TRUE
API_KEY_ENV <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05" # Environment variable name (NOT the key itself!)
DEEPSEEK_MODEL <- "deepseek-chat" # Model name for fanyi::chat_request

# ----- Reference Databases -----
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.xlsx"

# ----- Gene Filtering -----
MIN_CELLS_PER_GENE <- 3
MIN_GENES_PER_CELL <- 200

# ----- Output Options -----
SAVE_INDIVIDUAL_REPORTS <- TRUE
GENERATE_COMBINED_REPORT <- TRUE
SAVE_ENRICHMENT_OBJECTS <- TRUE
GENERATE_DOTPLOTS <- TRUE

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("B Cell Subcluster LLM Interpretation Analysis v1.0-OPTIMIZED\n")
cat(
  "================================================================================\n\n"
)

cat("[INFO] Loading required libraries...\n")

suppressPackageStartupMessages({
  library(reticulate) # Python interface
  library(SCNT) # GetSeurat function
  library(Seurat)
  library(future) # Parallel processing backend
  library(future.apply) # future_lapply function
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rio) # For reading Excel/various file formats
  library(fanyi) # For LLM API calls (chat_request)
})

cat("[OK] Libraries loaded\n\n")

# ==============================================================================
# VERIFY SEURAT VERSION AND TEST METHOD COMPATIBILITY
# ==============================================================================

cat("[INFO] Checking Seurat version compatibility...\n")

seurat_version <- packageVersion("Seurat")
cat(sprintf("  Seurat version: %s\n", seurat_version))

# Auto-correct test method name for Seurat v5
if (MARKER_METHOD == "wilcoxon" && seurat_version >= "5.0.0") {
  cat("[WARN] Detected 'wilcoxon' test with Seurat v5\n")
  cat("  Auto-correcting to 'wilcox' for compatibility\n")
  MARKER_METHOD <- "wilcox"
} else if (MARKER_METHOD == "wilcox" && seurat_version < "5.0.0") {
  cat("[WARN] Detected 'wilcox' test with Seurat v4\n")
  cat("  Auto-correcting to 'wilcoxon' for compatibility\n")
  MARKER_METHOD <- "wilcoxon"
}

cat(sprintf("  Final test method: %s\n", MARKER_METHOD))
cat("[OK] Version compatibility verified\n\n")

# ==============================================================================
# SETUP PARALLEL PROCESSING
# ==============================================================================

cat("[INFO] Setting up parallel processing...\n")
cat(sprintf("  Using %d cores\n", N_CORES))

# P0-3 FIX: Auto-select plan based on OS
# Windows doesn't support multicore (fork-based)
if (.Platform$OS.type == "windows") {
  cat("  [INFO] Windows detected - using multisession\n")
  plan("multisession", workers = N_CORES)
} else {
  cat("  [INFO] Unix-like OS - using multicore\n")
  plan("multicore", workers = N_CORES)
}

options(future.globals.maxSize = FUTURE_MAX_SIZE)

# P0-3 FIX: Limit BLAS/OpenMP threads to prevent CPU oversubscription
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)
cat("  [INFO] BLAS threads limited to 1 per worker\n")

cat("[OK] Parallel processing configured\n\n")

# ==============================================================================
# SETUP PYTHON ENVIRONMENT
# ==============================================================================

cat("[INFO] Configuring Python environment for GetSeurat...\n")

use_condaenv("bbknn_env", required = TRUE)
py_config()

cat("[OK] Python environment configured\n\n")

# ==============================================================================
# SETUP OUTPUT DIRECTORY
# ==============================================================================

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("[INFO] Created output directory: %s\n", OUTPUT_DIR))
}

# Create subdirectories
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "enrichment_objects"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "marker_genes"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "per_celltype"), showWarnings = FALSE)

cat("[OK] Output directory structure created\n\n")

# ==============================================================================
# SETUP LLM API
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("LLM API Configuration\n")
cat(
  "================================================================================\n\n"
)

if (USE_DEEPSEEK) {
  api_key <- Sys.getenv(API_KEY_ENV)

  if (api_key == "") {
    cat("[WARN] No API key found in environment variable:", API_KEY_ENV, "\n")
    cat("[WARN] LLM interpretation will be skipped\n")
    cat("[INFO] To enable LLM features, set your API key:\n")
    cat(sprintf("  Sys.setenv(%s = 'sk-your-actual-key')\n", API_KEY_ENV))
    cat("[INFO] Or add to your .Renviron file:\n")
    cat(sprintf("  echo '%s=sk-your-actual-key' >> ~/.Renviron\n", API_KEY_ENV))
    USE_LLM <- FALSE
  } else {
    cat("[OK] Found API key in environment\n")
    cat(sprintf("  Model: %s\n", DEEPSEEK_MODEL))

    # Configure fanyi with DeepSeek API
    tryCatch(
      {
        set_translate_option(
          key = api_key,
          source = "deepseek"
        )

        # Test API with a simple request
        test_response <- chat_request(
          x = "test",
          model = DEEPSEEK_MODEL
        )

        cat("[OK] DeepSeek API configured via fanyi\n")
        cat("[OK] API connection verified\n")
        USE_LLM <- TRUE
      },
      error = function(e) {
        cat(
          "[WARN] Failed to configure DeepSeek API:",
          conditionMessage(e),
          "\n"
        )
        cat("[WARN] LLM interpretation will be skipped\n")
        USE_LLM <- FALSE
      }
    )
  }
} else {
  cat("[INFO] LLM interpretation disabled\n")
  USE_LLM <- FALSE
}

cat("\n")

# ==============================================================================
# P0-2 FIX: Initialize LLM result objects
# ==============================================================================

# Initialize to NULL to prevent "object not found" errors when USE_LLM = FALSE
annotation_results <- NULL
phenotype_results <- NULL
celltype_interpretations <- list()

cat("[INFO] LLM result objects initialized\n\n")

# ==============================================================================
# PACKAGE FUNCTIONS USAGE NOTE
# ==============================================================================

# This script uses interpret() functions from the custom clusterProfiler package
# The package provides:
#   - interpret(enrichment_results, context, task, model)
#   - interpret_hierarchical(minor_data, major_data, mapping, task, model)
#   - interpret_agent(enrichment_results, model)
#
# These functions are already implemented in the package and will use
# the API key configured via fanyi::set_translate_option()
#
# No custom interpret() function is defined here - we rely on the package

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Fast gene filtering with compiled regex
filter_genes_optimized <- function(seurat_obj, min_cells = 3) {
  cat("[INFO] Filtering low-quality genes (optimized)...\n")

  DefaultAssay(seurat_obj) <- "RNA"
  all_genes <- rownames(seurat_obj)

  # Combined pattern for efficiency
  combined_pattern <- paste(
    c(
      "^MT-", # Mitochondrial
      "^RPS|^RPL|^MRPS|^MRPL", # Ribosomal
      "^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$", # Ribosomal pseudogenes
      "^(AC|AL|AP|BX|Z)[0-9]+\\.", # Unannotated
      "^RP[0-9]+-", # Ribosomal protein genes
      "^CTD-|^CTB-|^CTC-", # Chromosome transferred
      "-OT[0-9]+$", # Overlapping transcripts
      "^LOC[0-9]+" # LOC genes
    ),
    collapse = "|"
  )

  genes_to_remove <- grep(combined_pattern, all_genes, value = TRUE)
  cat(sprintf("  Pattern-based removal: %d genes\n", length(genes_to_remove)))

  genes_to_keep <- setdiff(all_genes, genes_to_remove)
  seurat_obj <- subset(seurat_obj, features = genes_to_keep)

  # Min cells filter using sparse matrix operations
  counts_mat <- tryCatch(
    LayerData(seurat_obj, layer = "counts"),
    error = function(e) GetAssayData(seurat_obj, slot = "counts")
  )

  gene_ncells <- Matrix::rowSums(counts_mat > 0)
  keep_genes <- names(gene_ncells[gene_ncells >= min_cells])

  before_n <- nrow(seurat_obj)
  seurat_obj <- subset(seurat_obj, features = keep_genes)
  after_n <- nrow(seurat_obj)

  cat(sprintf(
    "  Final: %d -> %d genes (removed %d)\n",
    before_n,
    after_n,
    before_n - after_n
  ))

  return(seurat_obj)
}

# Auto-detect correct test name for Seurat version
get_wilcox_test_name <- function() {
  seurat_version <- packageVersion("Seurat")
  if (seurat_version >= "5.0.0") {
    return("wilcox")
  } else {
    return("wilcoxon")
  }
}

# Parallel FindAllMarkers
find_markers_parallel <- function(
  seurat_obj,
  group_by,
  test_use = "wilcoxon",
  only_pos = TRUE,
  min_pct = 0.25,
  logfc_threshold = 0.25
) {
  cat("[INFO] Finding markers in parallel...\n")

  Idents(seurat_obj) <- group_by
  clusters <- levels(Idents(seurat_obj))
  cat(sprintf(
    "  Processing %d clusters with %d cores\n",
    length(clusters),
    N_CORES
  ))

  # Parallel processing by cluster
  marker_list <- future_lapply(
    clusters,
    function(cluster_id) {
      tryCatch(
        {
          FindMarkers(
            seurat_obj,
            ident.1 = cluster_id,
            only.pos = only_pos,
            min.pct = min_pct,
            logfc.threshold = logfc_threshold,
            test.use = test_use,
            verbose = FALSE
          )
        },
        error = function(e) {
          cat(sprintf(
            "  [WARN] Error in cluster %s: %s\n",
            cluster_id,
            e$message
          ))
          return(NULL)
        }
      )
    },
    future.seed = TRUE
  )

  names(marker_list) <- clusters
  marker_list <- marker_list[!sapply(marker_list, is.null)]

  # P1-2 FIX: Combine results with optimized column selection
  # Keep only essential columns to reduce memory footprint
  markers_df <- bind_rows(lapply(names(marker_list), function(cid) {
    df <- marker_list[[cid]]
    df$cluster <- cid
    df$gene <- rownames(df)

    # Keep only essential columns to reduce memory
    essential_cols <- c("gene", "cluster", "p_val", "p_val_adj")
    optional_cols <- c("avg_log2FC", "avg_logFC", "pct.1", "pct.2")
    available_optional <- optional_cols[optional_cols %in% colnames(df)]

    df[, c(essential_cols, available_optional)]
  }))

  cat(sprintf("  [OK] Found %d total markers\n", nrow(markers_df)))

  return(markers_df)
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Loading Data\n")
cat(
  "================================================================================\n\n"
)

cat("[INFO] Loading h5ad using GetSeurat()...\n")
cat("  Input:", H5AD_PATH, "\n")

start_time <- Sys.time()

seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

load_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\n[OK] Data loaded successfully (%.1f sec)\n", load_time))
cat(sprintf("  Cells: %d\n", ncol(seurat_obj)))
cat(sprintf("  Genes: %d\n", nrow(seurat_obj)))
cat(sprintf("  Assays: %s\n", paste(names(seurat_obj@assays), collapse = ", ")))

# Check required columns
required_cols <- c(CELLTYPE_L2_COL, CELLTYPE_L3_COL, SUBCLUSTER_ID_COL)
missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))

if (length(missing_cols) > 0) {
  stop(
    "[ERROR] Missing required metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("\n[INFO] Cell type distribution (Level 2):\n")
celltype_l2_counts <- table(seurat_obj@meta.data[[CELLTYPE_L2_COL]])
for (ct in names(celltype_l2_counts)) {
  cat(sprintf("  %-30s: %6d cells\n", ct, celltype_l2_counts[ct]))
}

cat("\n[INFO] Subcluster distribution (Level 3):\n")
celltype_l3_counts <- table(seurat_obj@meta.data[[CELLTYPE_L3_COL]])
cat(sprintf("  Total subclusters: %d\n", length(celltype_l3_counts)))
cat(sprintf("  Median size: %d cells\n", median(celltype_l3_counts)))
cat(sprintf(
  "  Range: %d - %d cells\n",
  min(celltype_l3_counts),
  max(celltype_l3_counts)
))

# ==============================================================================
# DATA PREPROCESSING
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Data Preprocessing\n")
cat(
  "================================================================================\n\n"
)

# Pre-filter low-quality cells
cat("[INFO] Pre-filtering low-quality cells...\n")
before_cells <- ncol(seurat_obj)

# Check if percent.mt exists
has_percent_mt <- "percent.mt" %in% colnames(seurat_obj@meta.data)

if (has_percent_mt) {
  cat("  Using percent.mt for quality filtering (threshold: <20%)\n")
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= MIN_GENES_PER_CELL &
      percent.mt < 20 &
      nCount_RNA > 0
  )
} else {
  cat("  [INFO] percent.mt not found - skipping mitochondrial filtering\n")
  cat("  [INFO] This is expected if MT genes were already filtered in Python\n")
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= MIN_GENES_PER_CELL &
      nCount_RNA > 0
  )
}

after_cells <- ncol(seurat_obj)
cat(sprintf(
  "  Cells: %d -> %d (removed %d)\n",
  before_cells,
  after_cells,
  before_cells - after_cells
))

# Filter genes
seurat_obj <- filter_genes_optimized(seurat_obj, min_cells = MIN_CELLS_PER_GENE)

# Normalize
cat("\n[INFO] Normalizing data...\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  nfeatures = 2000,
  verbose = FALSE
)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
cat("[OK] Normalization complete\n")

# ==============================================================================
# LOAD CELLMARKER DATABASE
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Loading CellMarker Database\n")
cat(
  "================================================================================\n\n"
)

cat("[INFO] Reading CellMarker database from:", CELLMARKER_PATH, "\n")

cellmarker_db <- import(CELLMARKER_PATH)

cat(sprintf("[OK] Loaded %d records\n", nrow(cellmarker_db)))

# Filter for B cell-relevant tissues
b_cell_keywords <- c(
  "B cell",
  "Plasma",
  "Lymph",
  "Blood",
  "Immune",
  "Spleen",
  "Bone marrow",
  "Lymph node"
)

if ("tissue_type" %in% colnames(cellmarker_db)) {
  cellmarker_filtered <- cellmarker_db %>%
    filter(grepl(
      paste(b_cell_keywords, collapse = "|"),
      tissue_type,
      ignore.case = TRUE
    ))

  if (nrow(cellmarker_filtered) > 0) {
    cat(sprintf(
      "[INFO] Filtered to %d B cell-relevant records\n",
      nrow(cellmarker_filtered)
    ))
    cellmarker_db <- cellmarker_filtered
  }
}

# P1-3 FIX: Clean and split CellMarker marker field
# Many CellMarker entries have comma/semicolon-separated gene lists
# Without splitting, enrichment will treat "CD79A,MS4A1" as a single gene
cat("[INFO] Preparing TERM2GENE with marker field cleaning...\n")

cellmarker_term2gene <- cellmarker_db %>%
  select(cell_name, marker) %>%
  # Split marker field by common delimiters
  tidyr::separate_rows(marker, sep = "[,;/\\s]+") %>%
  # Standardize to uppercase for consistency
  dplyr::mutate(marker = toupper(trimws(marker))) %>%
  # Remove empty entries
  dplyr::filter(marker != "" & !is.na(marker)) %>%
  # Remove duplicates
  distinct() %>%
  # Rename to standard enricher format
  dplyr::rename(term = cell_name, gene = marker)

cat(sprintf(
  "[OK] TERM2GENE prepared with %d unique cell type-gene pairs (after splitting)\n",
  nrow(cellmarker_term2gene)
))

# Show sample
sample_terms <- head(unique(cellmarker_term2gene$term), 3)
cat("[INFO] Sample terms:\n")
for (term in sample_terms) {
  genes <- cellmarker_term2gene$gene[cellmarker_term2gene$term == term]
  cat(sprintf(
    "  %s: %d genes (e.g., %s)\n",
    term,
    length(genes),
    paste(head(genes, 3), collapse = ", ")
  ))
}

# ==============================================================================
# COMPUTE MARKER GENES (PARALLEL)
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Computing Marker Genes (Parallel)\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("[INFO] Using method: %s\n", MARKER_METHOD))
cat(sprintf("  Grouping by: %s\n", CELLTYPE_L3_COL))
cat(sprintf(
  "  Thresholds: logFC >= %.2f, padj <= %.3f\n",
  MARKER_LOGFC_THRESHOLD,
  MARKER_PADJ_THRESHOLD
))

# P2-2 FIX: Force factor conversion for clustering columns
# This ensures proper Idents/levels behavior in FindMarkers
cat("[INFO] Converting clustering columns to factor...\n")
seurat_obj@meta.data[[CELLTYPE_L2_COL]] <- as.factor(seurat_obj@meta.data[[
  CELLTYPE_L2_COL
]])
seurat_obj@meta.data[[CELLTYPE_L3_COL]] <- as.factor(seurat_obj@meta.data[[
  CELLTYPE_L3_COL
]])
cat(sprintf(
  "  %s: %d levels\n",
  CELLTYPE_L2_COL,
  nlevels(seurat_obj@meta.data[[CELLTYPE_L2_COL]])
))
cat(sprintf(
  "  %s: %d levels\n",
  CELLTYPE_L3_COL,
  nlevels(seurat_obj@meta.data[[CELLTYPE_L3_COL]])
))

marker_start <- Sys.time()

all_markers <- find_markers_parallel(
  seurat_obj,
  group_by = CELLTYPE_L3_COL,
  test_use = MARKER_METHOD,
  only_pos = TRUE,
  min_pct = MIN_PCT,
  logfc_threshold = MARKER_LOGFC_THRESHOLD
)

marker_time <- as.numeric(difftime(Sys.time(), marker_start, units = "secs"))
cat(sprintf(
  "\n[OK] Marker finding complete (%.1f sec, %.1f markers/sec)\n",
  marker_time,
  nrow(all_markers) / marker_time
))

# Filter by p-value
all_markers <- all_markers %>%
  filter(p_val_adj <= MARKER_PADJ_THRESHOLD)

cat(sprintf("[OK] After p-value filter: %d markers\n", nrow(all_markers)))

# Save markers
marker_output <- file.path(
  OUTPUT_DIR,
  "marker_genes",
  "all_subcluster_markers.csv"
)
write.csv(all_markers, marker_output, row.names = FALSE)
cat(sprintf("[OK] Saved markers to: %s\n", basename(marker_output)))

# P0-4 FIX: Auto-detect logFC column name for Seurat version compatibility
logfc_col <- if ("avg_log2FC" %in% colnames(all_markers)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(all_markers)) {
  "avg_logFC"
} else {
  stop("Cannot find logFC column (expected avg_log2FC or avg_logFC)")
}
cat(sprintf("[INFO] Using logFC column: %s\n", logfc_col))

# Prepare top markers per cluster
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(.data[[logfc_col]])) %>%
  slice_head(n = TOP_N_MARKERS) %>%
  ungroup()

cat(sprintf(
  "[INFO] Selected top %d markers per cluster for enrichment\n",
  TOP_N_MARKERS
))

# ==============================================================================
# ENRICHMENT ANALYSIS
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Enrichment Analysis (CellMarker + GO)\n")
cat(
  "================================================================================\n\n"
)

# Get valid subclusters
subcluster_sizes <- table(seurat_obj@meta.data[[CELLTYPE_L3_COL]])
valid_subclusters <- names(subcluster_sizes[
  subcluster_sizes >= MIN_CELLS_FOR_ANALYSIS
])

cat(sprintf(
  "[INFO] Analyzing %d subclusters with >= %d cells\n",
  length(valid_subclusters),
  MIN_CELLS_FOR_ANALYSIS
))

# Prepare data
top_markers_filtered <- top_markers %>%
  filter(cluster %in% valid_subclusters) %>%
  select(gene, cluster)

# CellMarker enrichment
cat("\n[INFO] Running CellMarker enrichment...\n")

cellmarker_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers_filtered,
  fun = enricher,
  TERM2GENE = cellmarker_term2gene,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# P2-3 FIX: Safe handling of empty compareClusterResult
cat(sprintf("[OK] CellMarker enrichment complete\n"))
ccr_cm <- cellmarker_enrich@compareClusterResult
sig_n_cm <- if (is.null(ccr_cm) || nrow(ccr_cm) == 0) {
  0
} else {
  sum(ccr_cm$p.adjust < 0.05, na.rm = TRUE)
}
cat(sprintf("  Total results: %d\n", if (is.null(ccr_cm)) 0 else nrow(ccr_cm)))
cat(sprintf("  Significant terms (p.adjust < 0.05): %d\n", sig_n_cm))

if (SAVE_ENRICHMENT_OBJECTS) {
  saveRDS(
    cellmarker_enrich,
    file.path(OUTPUT_DIR, "enrichment_objects", "cellmarker_enrichment.rds")
  )
}

# GO enrichment
cat("\n[INFO] Running GO Biological Process enrichment...\n")

# P1-4 FIX: Filter for valid gene symbols before GO enrichment
# Many gene names may not be recognized by org.Hs.eg.db
cat("[INFO] Validating gene symbols against org.Hs.eg.db...\n")

valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
top_markers_go <- top_markers_filtered %>%
  dplyr::mutate(gene = toupper(gene)) %>%
  dplyr::filter(gene %in% valid_symbols)

# Report filtering stats per cluster
genes_per_cluster_before <- top_markers_filtered %>%
  group_by(cluster) %>%
  summarise(n_before = n(), .groups = "drop")
genes_per_cluster_after <- top_markers_go %>%
  group_by(cluster) %>%
  summarise(n_after = n(), .groups = "drop")

filter_stats <- left_join(
  genes_per_cluster_before,
  genes_per_cluster_after,
  by = "cluster"
) %>%
  mutate(n_after = ifelse(is.na(n_after), 0, n_after))

cat(sprintf(
  "  Total genes: %d -> %d (%.1f%% valid)\n",
  nrow(top_markers_filtered),
  nrow(top_markers_go),
  100 * nrow(top_markers_go) / nrow(top_markers_filtered)
))

# Warn about clusters with very few valid genes
low_gene_clusters <- filter_stats %>%
  filter(n_after < 5)

if (nrow(low_gene_clusters) > 0) {
  cat("[WARN] Some clusters have < 5 valid genes:\n")
  for (i in 1:nrow(low_gene_clusters)) {
    cat(sprintf(
      "  %s: %d genes\n",
      low_gene_clusters$cluster[i],
      low_gene_clusters$n_after[i]
    ))
  }
}

go_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers_go, # Use validated gene list
  fun = enrichGO,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# P2-3 FIX: Safe handling of empty compareClusterResult
cat(sprintf("[OK] GO enrichment complete\n"))
ccr_go <- go_enrich@compareClusterResult
sig_n_go <- if (is.null(ccr_go) || nrow(ccr_go) == 0) {
  0
} else {
  sum(ccr_go$p.adjust < 0.05, na.rm = TRUE)
}
cat(sprintf("  Total results: %d\n", if (is.null(ccr_go)) 0 else nrow(ccr_go)))
cat(sprintf("  Significant terms (p.adjust < 0.05): %d\n", sig_n_go))

if (SAVE_ENRICHMENT_OBJECTS) {
  saveRDS(
    go_enrich,
    file.path(OUTPUT_DIR, "enrichment_objects", "go_enrichment.rds")
  )
}

# ==============================================================================
# VISUALIZATION
# ==============================================================================

if (GENERATE_DOTPLOTS) {
  cat("\n")
  cat(
    "================================================================================\n"
  )
  cat("Generating Enrichment Visualizations\n")
  cat(
    "================================================================================\n\n"
  )

  # CellMarker dotplot
  cat("[INFO] Creating CellMarker dotplot...\n")

  pdf(
    file.path(OUTPUT_DIR, "figures", "cellmarker_enrichment_dotplot.pdf"),
    width = 16,
    height = 12
  )
  tryCatch(
    {
      print(
        dotplot(cellmarker_enrich, showCategory = 10, font.size = 7) +
          ggtitle("CellMarker Enrichment by B Cell Subcluster") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      cat("[OK] Saved: cellmarker_enrichment_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] CellMarker plot failed:", conditionMessage(e), "\n")
    }
  )
  dev.off()

  # GO dotplot
  cat("[INFO] Creating GO dotplot...\n")

  pdf(
    file.path(OUTPUT_DIR, "figures", "go_enrichment_dotplot.pdf"),
    width = 32,
    height = 28
  )
  tryCatch(
    {
      print(
        dotplot(go_enrich, showCategory = 15, font.size = 6) +
          ggtitle("GO Biological Process Enrichment by B Cell Subcluster") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      cat("[OK] Saved: go_enrichment_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] GO plot failed:", conditionMessage(e), "\n")
    }
  )
  dev.off()
}

# ==============================================================================
# LLM INTERPRETATION - ALL SUBCLUSTERS
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("LLM-Powered Interpretation (All Subclusters)\n")
cat(
  "================================================================================\n\n"
)

if (!USE_LLM) {
  cat("[WARN] LLM interpretation skipped (API not configured)\n")
  cat("[INFO] You can still manually inspect enrichment results\n")
} else {
  # Task 1: Cell subtype annotation
  cat("[INFO] Task 1: Cell subtype annotation (all subclusters)...\n")
  cat("  This may take several minutes...\n\n")

  annotation_results <- tryCatch(
    {
      # Prepare context for B cell analysis
      context_text <- paste(
        "B cells from chronic rhinosinusitis with nasal polyps (CRSwNP) study.",
        "These are subclusters of B cell populations derived from scANVI annotation.",
        "Major types include Memory B cells, Naive B cells, Plasma cells, and",
        "Plasmablasts. The disease is characterized by Type 2 inflammation with",
        "IL-4/IL-13 signaling, eosinophilia, and local antibody production.",
        "Looking for evidence of class-switching, activation status, antibody",
        "production, and tissue residency markers."
      )

      # Call package's interpret() with enrichment results
      # The package function expects compareClusterResult or data.frame
      interpret(
        go_enrich, # Use GO enrichment (primary source)
        task = "annotation",
        model = DEEPSEEK_MODEL,
        context = context_text
      )
    },
    error = function(e) {
      cat("[ERROR] Annotation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(annotation_results)) {
    cat("\n[OK] Annotation complete!\n\n")

    # Print summary
    cat(
      "================================================================================\n"
    )
    cat("ANNOTATION SUMMARY (ALL SUBCLUSTERS)\n")
    cat(
      "================================================================================\n\n"
    )

    for (i in seq_along(annotation_results)) {
      result <- annotation_results[[i]]
      cat(sprintf(
        "%-35s: %s (%s)\n",
        names(annotation_results)[i],
        result$cell_type,
        result$confidence
      ))
    }

    # Save results
    annotation_df <- data.frame(
      subcluster = names(annotation_results),
      cell_type = sapply(annotation_results, function(x) x$cell_type),
      confidence = sapply(annotation_results, function(x) x$confidence),
      stringsAsFactors = FALSE
    )

    write.csv(
      annotation_df,
      file.path(OUTPUT_DIR, "subcluster_annotations_llm.csv"),
      row.names = FALSE
    )

    saveRDS(
      annotation_results,
      file.path(OUTPUT_DIR, "reports", "annotation_results_full.rds")
    )

    cat("\n[OK] Saved annotation results\n")
  }

  # Task 2: Functional phenotyping
  cat("\n[INFO] Task 2: Functional phenotyping (all subclusters)...\n")

  phenotype_results <- tryCatch(
    {
      context_text <- paste(
        "B cells from CRSwNP patients showing Type 2 inflammation.",
        "Looking for functional states like:",
        "- Activation status (resting vs activated)",
        "- Antibody production (IgE vs IgG vs IgA)",
        "- Proliferation and differentiation",
        "- Inflammatory or regulatory phenotypes",
        "- Tissue residency vs circulating"
      )

      interpret(
        go_enrich,
        task = "phenotype", # Note: test file uses "phenotype" not "phenotyping"
        model = DEEPSEEK_MODEL,
        context = context_text
      )
    },
    error = function(e) {
      cat("[ERROR] Phenotyping failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(phenotype_results)) {
    cat("\n[OK] Phenotyping complete!\n\n")

    # Save results
    phenotype_df <- data.frame(
      subcluster = names(phenotype_results),
      phenotype = sapply(phenotype_results, function(x) x$phenotype),
      confidence = sapply(phenotype_results, function(x) x$confidence),
      stringsAsFactors = FALSE
    )

    write.csv(
      phenotype_df,
      file.path(OUTPUT_DIR, "subcluster_phenotypes_llm.csv"),
      row.names = FALSE
    )

    saveRDS(
      phenotype_results,
      file.path(OUTPUT_DIR, "reports", "phenotype_results_full.rds")
    )

    cat("[OK] Saved phenotype results\n")
  }
}

# ==============================================================================
# PER-CELLTYPE DETAILED ANALYSIS
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Per-Celltype Detailed Analysis\n")
cat(
  "================================================================================\n\n"
)

if (USE_LLM) {
  # Get unique L2 celltypes
  celltypes_l2 <- unique(seurat_obj@meta.data[[CELLTYPE_L2_COL]])
  celltypes_l2 <- celltypes_l2[!is.na(celltypes_l2)]

  cat(sprintf(
    "[INFO] Analyzing %d cell types individually\n",
    length(celltypes_l2)
  ))

  celltype_interpretations <- list()

  for (celltype in celltypes_l2) {
    cat(sprintf("\n--- %s ---\n", celltype))

    # P0-1 FIX: Robust subcluster selection using L2-to-L3 mapping
    # Don't rely on string prefix assumptions (space vs underscore issues)
    # Instead, directly filter by L2 and get the corresponding L3 values
    celltype_mask <- seurat_obj@meta.data[[CELLTYPE_L2_COL]] == celltype
    celltype_subclusters <- unique(seurat_obj@meta.data[[CELLTYPE_L3_COL]][
      celltype_mask
    ])
    celltype_subclusters <- celltype_subclusters[!is.na(celltype_subclusters)]

    n_cells <- sum(celltype_mask)
    n_subclusters <- length(celltype_subclusters)

    cat(sprintf("  Cells: %d, Subclusters: %d\n", n_cells, n_subclusters))

    # Skip if too few cells or only 1 subcluster
    if (n_cells < MIN_CELLS_FOR_INDIVIDUAL_ANALYSIS || n_subclusters <= 1) {
      cat("  Skipping (too few cells or subclusters)\n")
      next
    }

    # Filter enrichment for this celltype
    cm_filtered <- cellmarker_enrich@compareClusterResult %>%
      filter(Cluster %in% celltype_subclusters)

    go_filtered <- go_enrich@compareClusterResult %>%
      filter(Cluster %in% celltype_subclusters)

    if (nrow(cm_filtered) == 0 && nrow(go_filtered) == 0) {
      cat("  Skipping (no enrichment results)\n")
      next
    }

    # Create subset enrichment objects
    cm_subset <- cellmarker_enrich
    cm_subset@compareClusterResult <- cm_filtered
    go_subset <- go_enrich
    go_subset@compareClusterResult <- go_filtered

    # Generate mechanistic interpretation
    cat("  Generating mechanistic interpretation...\n")

    interpretation <- tryCatch(
      {
        context_text <- paste(
          celltype,
          "cells from CRSwNP patients.",
          "These subclusters represent functional states within the",
          celltype,
          "population. Looking to understand:",
          "- Activation and differentiation states",
          "- Antibody class-switching patterns",
          "- Response to Type 2 cytokines (IL-4/IL-13)",
          "- Tissue retention vs trafficking signatures",
          "- Interactions with other immune cells"
        )

        interpret(
          go_subset,
          task = "interpretation",
          model = DEEPSEEK_MODEL,
          context = context_text
        )
      },
      error = function(e) {
        cat("  [WARN] Interpretation failed:", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    if (!is.null(interpretation)) {
      celltype_interpretations[[celltype]] <- interpretation

      # Save individual celltype report
      if (SAVE_INDIVIDUAL_REPORTS) {
        report_file <- file.path(
          OUTPUT_DIR,
          "per_celltype",
          paste0(gsub(" ", "_", celltype), "_interpretation.md")
        )

        sink(report_file)
        cat("# ", celltype, " Interpretation Report\n\n", sep = "")
        cat(
          "**Generated:** ",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          "\n\n",
          sep = ""
        )
        cat("---\n\n")

        cat("## Overview\n\n")
        if (!is.null(interpretation$overview)) {
          cat(interpretation$overview, "\n\n")
        }

        cat("## Key Mechanisms\n\n")
        if (!is.null(interpretation$key_mechanisms)) {
          cat(interpretation$key_mechanisms, "\n\n")
        }

        cat("## Pathway Crosstalk\n\n")
        if (!is.null(interpretation$crosstalk)) {
          cat(interpretation$crosstalk, "\n\n")
        }

        cat("## Testable Hypotheses\n\n")
        if (!is.null(interpretation$hypothesis)) {
          cat(interpretation$hypothesis, "\n\n")
        }

        cat("## Publication-Ready Narrative\n\n")
        if (!is.null(interpretation$narrative)) {
          cat(interpretation$narrative, "\n\n")
        }

        sink()

        cat("  [OK] Saved:", basename(report_file), "\n")
      }
    }
  }

  # Save all celltype interpretations
  saveRDS(
    celltype_interpretations,
    file.path(OUTPUT_DIR, "reports", "celltype_interpretations_full.rds")
  )
}

# ==============================================================================
# GENERATE COMBINED REPORT
# ==============================================================================

if (GENERATE_COMBINED_REPORT && USE_LLM) {
  cat("\n")
  cat(
    "================================================================================\n"
  )
  cat("Generating Combined Report\n")
  cat(
    "================================================================================\n\n"
  )

  report_file <- file.path(OUTPUT_DIR, "BCELL_INTERPRETATION_REPORT.md")

  sink(report_file)

  cat("# B Cell Subcluster Interpretation Report\n\n")
  cat(
    "**Generated:** ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "\n\n",
    sep = ""
  )
  cat("---\n\n")

  cat("## Overview\n\n")
  cat(
    "LLM-powered interpretations of B cell subclusters from CRSwNP analysis.\n\n"
  )

  cat("### Dataset Summary\n\n")
  cat(sprintf("- Total B cells: %d\n", ncol(seurat_obj)))
  cat(sprintf(
    "- Cell types (Level 2): %d\n",
    length(unique(seurat_obj@meta.data[[CELLTYPE_L2_COL]]))
  ))
  cat(sprintf("- Subclusters analyzed: %d\n", length(valid_subclusters)))
  cat(sprintf("- Minimum cells per cluster: %d\n", MIN_CELLS_FOR_ANALYSIS))
  cat("\n")

  if (!is.null(annotation_results)) {
    cat("---\n\n")
    cat("## Cell Subtype Annotations\n\n")

    for (i in seq_along(annotation_results)) {
      result <- annotation_results[[i]]
      cat(sprintf("### %s\n\n", names(annotation_results)[i]))
      cat(sprintf("**Cell Type:** %s  \n", result$cell_type))
      cat(sprintf("**Confidence:** %s  \n\n", result$confidence))

      if (!is.null(result$reasoning)) {
        cat("**Reasoning:**\n\n")
        cat(result$reasoning, "\n\n")
      }

      if (!is.null(result$supporting_markers)) {
        cat("**Supporting Markers:**\n\n")
        cat(result$supporting_markers, "\n\n")
      }

      cat("---\n\n")
    }
  }

  if (!is.null(phenotype_results)) {
    cat("---\n\n")
    cat("## Functional Phenotypes\n\n")

    for (i in seq_along(phenotype_results)) {
      result <- phenotype_results[[i]]
      cat(sprintf("### %s\n\n", names(phenotype_results)[i]))
      cat(sprintf("**Phenotype:** %s  \n", result$phenotype))
      cat(sprintf("**Confidence:** %s  \n\n", result$confidence))

      if (!is.null(result$reasoning)) {
        cat("**Reasoning:**\n\n")
        cat(result$reasoning, "\n\n")
      }

      cat("---\n\n")
    }
  }

  if (length(celltype_interpretations) > 0) {
    cat("---\n\n")
    cat("## Per-Celltype Mechanistic Interpretations\n\n")

    for (celltype in names(celltype_interpretations)) {
      interpretation <- celltype_interpretations[[celltype]]

      cat(sprintf("### %s\n\n", celltype))

      if (!is.null(interpretation$narrative)) {
        cat(interpretation$narrative, "\n\n")
      }

      cat("---\n\n")
    }
  }

  sink()

  cat("[OK] Report saved:", basename(report_file), "\n")
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

cat("\n")
cat(
  "================================================================================\n"
)
cat("ANALYSIS COMPLETE\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Total time: %.1f minutes\n", total_time))
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat("  - BCELL_INTERPRETATION_REPORT.md\n")
cat("  - subcluster_annotations_llm.csv\n")
cat("  - subcluster_phenotypes_llm.csv\n")
cat("  - figures/cellmarker_enrichment_dotplot.pdf\n")
cat("  - figures/go_enrichment_dotplot.pdf\n")
cat("  - marker_genes/all_subcluster_markers.csv\n")
if (SAVE_INDIVIDUAL_REPORTS) {
  cat("  - per_celltype/[celltype]_interpretation.md\n")
}
cat("\n")

if (USE_LLM) {
  cat("LLM interpretation: SUCCESS\n")
  cat("  Load results:\n")
  cat("    readRDS('reports/annotation_results_full.rds')\n")
  cat("    readRDS('reports/phenotype_results_full.rds')\n")
  cat("    readRDS('reports/celltype_interpretations_full.rds')\n\n")
} else {
  cat("LLM interpretation: SKIPPED (configure API key to enable)\n\n")
}

cat("Confidence score distribution:\n")
if (!is.null(annotation_results)) {
  conf_table <- table(sapply(annotation_results, function(x) x$confidence))
  for (conf in names(conf_table)) {
    cat(sprintf("  %s: %d subclusters\n", conf, conf_table[conf]))
  }
}

cat("\n")
cat(
  "================================================================================\n"
)
cat("Saving Session Info and Checkpoints\n")
cat(
  "================================================================================\n\n"
)

# P2-1 FIX: Save sessionInfo for reproducibility
cat("[INFO] Saving sessionInfo...\n")
sessioninfo_path <- file.path(OUTPUT_DIR, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), sessioninfo_path)
cat(sprintf("[OK] Saved: %s\n", basename(sessioninfo_path)))

# P2-1 FIX: Save QC-filtered Seurat object as checkpoint
cat("[INFO] Saving post-QC Seurat object checkpoint...\n")
checkpoint_path <- file.path(OUTPUT_DIR, "seurat_post_qc.rds")
tryCatch(
  {
    saveRDS(seurat_obj, checkpoint_path)
    cat(sprintf(
      "[OK] Saved: %s (%.1f MB)\n",
      basename(checkpoint_path),
      file.size(checkpoint_path) / 1024^2
    ))
  },
  error = function(e) {
    cat(sprintf("[WARN] Failed to save checkpoint: %s\n", conditionMessage(e)))
  }
)

cat("\n")
cat(
  "================================================================================\n"
)
cat("END OF SCRIPT\n")
cat(
  "================================================================================\n"
)


API_KEY_ENV <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05"
MODEL <- "deepseek-chat"

# 运行前在 shell 里设置（推荐），或在 R 里临时设置：
Sys.setenv(DEEPSEEK_API_KEY = "sk-ed1879cf6fa14b04aac9cb6c078a3d05")

fanyi::set_translate_option(
  key = Sys.getenv(API_KEY_ENV),
  source = "deepseek"
)

ctx <- paste(
  "Human normal tissue B cells (non-diseased baseline).",
  "Goal: annotate B-cell subclusters and summarize baseline functional states.",
  "Focus on: naive vs memory vs GC-related programs, class-switching, plasma cell differentiation,",
  "antigen presentation, proliferation/cell cycle, interferon-like signatures (if present), and tissue residency/trafficking."
)

# # 设置 API 密钥
# Sys.setenv(DEEPSEEK_API_KEY = "sk-ed1879cf6fa14b04aac9cb6c078a3d05")
# fanyi::set_translate_option(
#   key = Sys.getenv("DEEPSEEK_API_KEY"),
#   source = "deepseek"
# )
# test_response <- tryCatch(
#   {
#     fanyi::chat_request("test", model = "deepseek-chat")
#   },
#   error = function(e) {
#     cat("[ERROR] Failed to connect to DeepSeek API:", conditionMessage(e), "\n")
#     return(NULL)
#   }
# )

# if (!is.null(test_response)) {
#   cat("[OK] API connection successful\n")
# }

MODEL <- "deepseek-chat"

y_anno <- clusterProfiler::interpret(
  list(go_bp_enrich),
  task = "annotation",
  context = ctx,
  model = MODEL
)

y_pheno <- clusterProfiler::interpret(
  list(go_bp_enrich),
  task = "phenotype", # 若你的包要求 "phenotyping" 就改这个字符串
  context = ctx,
  model = MODEL
)

y_mech <- clusterProfiler::interpret(
  go_bp_enrich,
  task = "interpretation",
  context = ctx,
  model = MODEL
)
