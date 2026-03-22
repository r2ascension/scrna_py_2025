#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Subcluster LLM Interpretation Analysis v1.3.1-HOTFIX
# ==============================================================================

# Author: r2end
# Date: 2026-01-26
# Version: v1.3.1-HOTFIX (Production-grade fixes applied)

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
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv" # PanglaoDB 数据路径

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
  library(fanyi) # For LLM API calls (chat_request)
})

cat("[OK] Libraries loaded\n\n")

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
# LOAD PANGLAODB DATA
# ==============================================================================

cat("[INFO] Loading PanglaoDB database from:", PANGLAODB_PATH, "\n")

# Load PanglaoDB data
panglaodb_data <- import(PANGLAODB_PATH)

# Format PanglaoDB for clusterProfiler (TERM2GENE format)
panglaodb_term2gene <- panglaodb_data %>%
  select(cell_type, marker) %>%
  separate_rows(marker, sep = ",") %>%
  mutate(marker = toupper(trimws(marker))) %>%
  filter(marker != "" & !is.na(marker)) %>%
  distinct() %>%
  rename(term = cell_type, gene = marker)

cat(
  "[INFO] PanglaoDB TERM2GENE prepared with",
  nrow(panglaodb_term2gene),
  "unique cell type-gene pairs\n"
)

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

# Parallel marker discovery
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

# ==============================================================================
# ENRICHMENT ANALYSIS (PanglaoDB)
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("Enrichment Analysis (PanglaoDB + GO)\n")
cat(
  "================================================================================\n\n"
)

# Filter markers by p-value
all_markers <- all_markers %>%
  filter(p_val_adj <= MARKER_PADJ_THRESHOLD)

# Prepare data for enrichment analysis
top_markers_filtered <- all_markers %>%
  select(gene, cluster)

# PanglaoDB enrichment
cat("\n[INFO] Running PanglaoDB enrichment...\n")

panglaodb_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers_filtered,
  fun = enricher,
  TERM2GENE = panglaodb_term2gene,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# ==============================================================================
# GO ENRICHMENT ANALYSIS
# ==============================================================================

cat("\n[INFO] Running GO Biological Process enrichment...\n")

go_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers_filtered,
  fun = enrichGO,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

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

  # PanglaoDB DotPlot
  cat("[INFO] Creating PanglaoDB enrichment dotplot...\n")

  pdf(
    file.path(OUTPUT_DIR, "figures", "panglaodb_enrichment_dotplot.pdf"),
    width = 16,
    height = 12
  )
  tryCatch(
    {
      print(
        dotplot(panglaodb_enrich, showCategory = 10, font.size = 7) +
          ggtitle("PanglaoDB Enrichment by B Cell Subcluster") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      cat("[OK] Saved: panglaodb_enrichment_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] PanglaoDB plot failed:", conditionMessage(e), "\n")
    }
  )
  dev.off()

  # GO DotPlot
  cat("[INFO] Creating GO enrichment dotplot...\n")

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
} else {
  # Task 1: Cell subtype annotation
  annotation_results <- tryCatch(
    {
      interpret(
        panglaodb_enrich,
        task = "annotation",
        context = "Normal human tissue B cells (non-diseased baseline). Focus on naive vs memory vs GC-related programs.",
        model = DEEPSEEK_MODEL
      )
    },
    error = function(e) {
      cat("[ERROR] Annotation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(annotation_results)) {
    cat("\n[OK] Annotation complete!\n")
    # Save results
    write.csv(
      annotation_results,
      file.path(OUTPUT_DIR, "subcluster_annotations_llm.csv")
    )
  }
}

cat("\n")
cat(
  "================================================================================\n"
)
cat("END OF SCRIPT\n")
cat(
  "================================================================================\n"
)
