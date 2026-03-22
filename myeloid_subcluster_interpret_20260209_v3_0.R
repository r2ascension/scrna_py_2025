#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Cell Subcluster Interpretation Pipeline - v3.1-MYELOID
# ==============================================================================
#
# Version: v3.1-MYELOID (2026-02-09)
# Status: Production-ready Myeloid specialized version
# Base: v3.0-HOTFIX + Multi-ssGSEA + Load existing enrichment
#
# Key Features:
#   ✅ Myeloid expert context (M1/M2, monocyte states, DC maturation)
#   ✅ ALL ssGSEA methods (hallmark/GO BP/MF/CC/KEGG)
#   ✅ Load existing enrichment RDS (optional)
#   ✅ No cell type filtering (analyze all)
#   ✅ All v1.1 bug fixes integrated
#
# New in v3.1:
#   • Run ALL ssGSEA methods in single script
#   • Optional loading of pre-computed enrichment
#   • No cell type filtering by default
#
# ==============================================================================

# ==============================================================================
# Configuration Parameters
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0128/myeloid_analysis_unified/results/subcluster_unified_v2_20260128/adata_myeloid_subclustered_FINAL_v2_20260128.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0131/myeloid_interpret_v2_6"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

# MSigDB GO GMT Files
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"


# Optional: Custom myeloid marker database
CUSTOM_MARKERS_CSV <- NULL  # e.g., "/path/to/myeloid_markers.csv"

# Analysis Scope
# ⭐ UPDATED: No cell type filtering - analyze ALL subclusters
CELLTYPE_SPECIFIC <- FALSE  # Changed to FALSE - analyze all cells
# If you want to filter to specific celltype, set to TRUE and specify below:
# CELLTYPE_TO_ANALYZE <- "Myeloid"  # Or "Myeloid cells"
CELLTYPE_COLUMN <- "cell_type_L2"
SUBCLUSTER_COLUMN <- "cell_type_L3"

# Data Processing
FORCE_NORMALIZE <- FALSE    # Set TRUE to force LogNormalize even if data appears normalized

# Analysis Methods
RUN_SSGSEA <- TRUE
RUN_ENRICHMENT <- TRUE
RUN_INTERPRET_AGENT <- TRUE

# ⭐ NEW: Load pre-computed enrichment results (skip computation if available)
LOAD_EXISTING_ENRICHMENT <- FALSE  # Set TRUE to load from reports/*.rds
ENRICHMENT_DIR <- NULL  # Path to directory with *.rds files (default: OUTPUT_DIR/reports)

# ⭐ UPDATED: ssGSEA Configuration - Run ALL methods
SSGSEA_METHODS <- c("hallmark", "go_bp", "go_mf", "go_cc", "kegg")  # Run all
# Available: "hallmark", "go_bp", "go_mf", "go_cc", "kegg", "custom"
SSGSEA_CUSTOM_GMT <- NULL  # Path if using "custom"
SSGSEA_N_TOP <- 20  # Top pathways per method for LLM context

# LLM Configuration
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
INTERPRET_AGENT_MODEL <- "deepseek-reasoner"
INTERPRET_AGENT_N_PATHWAYS <- 50
INTERPRET_AGENT_ADD_PPI <- TRUE  # Inflammatory cascades and signaling important

# Computational Parameters
N_CORES <- 4
TOP_N_MARKERS <- 50
MEMORY_LIMIT_GB <- 48

# ==============================================================================
# Thread Limiting (v1.1 bug fix)
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ==============================================================================
# Load Libraries
# ==============================================================================

cat("\n=== Loading Libraries ===\n")

suppressPackageStartupMessages({
  library(reticulate)
  library(SCNT)
  library(Seurat)
  library(clusterProfiler)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(future)
  library(future.apply)
  library(httr)
  library(jsonlite)
  library(rlang)
  library(GSVA)
  library(BiocParallel)
})

# v1.1 fix: Use multisession (spawn) instead of multicore (fork)
plan("multisession", workers = N_CORES)
options(future.globals.maxSize = MEMORY_LIMIT_GB * 1024^3)

# Setup Python
use_condaenv("bbknn_env", required = TRUE)

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports", "interpret_agent_txt"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports", "llm_raw"),
           recursive = TRUE, showWarnings = FALSE)

cat("[OK] Libraries loaded and directories created\n\n")

# ==============================================================================
# Validate API Key
# ==============================================================================

cat("=== Configuring DeepSeek API ===\n")

if (nchar(DEEPSEEK_API_KEY) < 10) {
  stop(
    "ERROR: DEEPSEEK_API_KEY environment variable not set or invalid\n",
    "Please set it: export DEEPSEEK_API_KEY='your-key-here'"
  )
}

# FIX 0-D: Remove fanyi dependency, use httr directly for API testing
test_response <- tryCatch(
  {
    payload <- list(
      model = "deepseek-reasoner",
      messages = list(list(role = "user", content = "test")),
      max_tokens = 10
    )
    
    response <- httr::POST(
      url = "https://api.deepseek.com/v1/chat/completions",
      httr::add_headers(
        "Authorization" = paste("Bearer", DEEPSEEK_API_KEY),
        "Content-Type" = "application/json"
      ),
      body = payload,
      encode = "json",
      httr::timeout(30)
    )
    
    if (httr::status_code(response) == 200) {
      cat("[OK] DeepSeek API connection successful\n\n")
      TRUE
    } else {
      cat(sprintf("[ERROR] API test failed with status %d\n", httr::status_code(response)))
      content <- httr::content(response, as = "text", encoding = "UTF-8")
      cat("Response:", content, "\n")
      FALSE
    }
  },
  error = function(e) {
    cat("[ERROR] Failed to connect to DeepSeek API:", conditionMessage(e), "\n")
    FALSE
  }
)

if (!test_response) {
  stop("Failed to connect to DeepSeek API. Please check your API key and network.")
}

# ==============================================================================
# Load Seurat Data
# ==============================================================================

cat("=== Loading Seurat Data ===\n")

seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

cat(sprintf("\nLoaded: %d cells x %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))

# v1.1 fix: Fail-fast validation for required columns
need_cols <- c(CELLTYPE_COLUMN, SUBCLUSTER_COLUMN)
missing <- setdiff(need_cols, colnames(seurat_obj@meta.data))
if (length(missing) > 0) {
  stop("ERROR: Missing required metadata columns: ", paste(missing, collapse = ", "))
}

cat("[OK] Required metadata columns present\n")

# Filter to myeloid cells if analyzing specific celltype
if (CELLTYPE_SPECIFIC && !is.null(CELLTYPE_TO_ANALYZE)) {
  n_before <- ncol(seurat_obj)
  
  # FIX 0-F: Use safer subset method instead of get() with NSE
  # Handle both exact match and partial match
  celltype_values <- seurat_obj@meta.data[[CELLTYPE_COLUMN]]
  
  if (CELLTYPE_TO_ANALYZE %in% celltype_values) {
    # Exact match
    keep_cells <- rownames(seurat_obj@meta.data)[celltype_values == CELLTYPE_TO_ANALYZE]
    seurat_obj <- subset(seurat_obj, cells = keep_cells)
  } else {
    # Try partial match
    matching_celltypes <- grep(CELLTYPE_TO_ANALYZE, 
                               unique(celltype_values),
                               value = TRUE, ignore.case = TRUE)
    if (length(matching_celltypes) > 0) {
      cat(sprintf("[INFO] Partial match found: %s\n", 
                  paste(matching_celltypes, collapse = ", ")))
      keep_cells <- rownames(seurat_obj@meta.data)[celltype_values %in% matching_celltypes]
      seurat_obj <- subset(seurat_obj, cells = keep_cells)
    } else {
      stop(sprintf("No cells found for celltype '%s'!", CELLTYPE_TO_ANALYZE))
    }
  }
  
  n_after <- ncol(seurat_obj)
  cat(sprintf("\n[INFO] Filtered to myeloid cells: %d → %d cells\n", n_before, n_after))
  
  if (n_after == 0) {
    stop("No cells remaining after filtering!")
  }
}

# FIX 0-E: Safe data slot access (avoid @x on non-sparse matrices)
# FIX 1.1: Use explicit FORCE_NORMALIZE parameter
data_slot <- GetAssayData(seurat_obj, slot = "data")

if (FORCE_NORMALIZE) {
  cat("[INFO] FORCE_NORMALIZE=TRUE, normalizing data\n")
  skip_normalize <- FALSE
} else {
  # Check if data appears to be log-normalized
  data_max <- if (inherits(data_slot, "dgCMatrix")) {
    if (length(data_slot@x) > 0) max(data_slot@x) else 0
  } else {
    max(data_slot)
  }
  
  if (data_max < 20 && data_max > 0) {
    cat(sprintf("[INFO] Data appears log-normalized (max=%.2f), skipping NormalizeData\n", data_max))
    cat("[INFO] Set FORCE_NORMALIZE=TRUE to override this behavior\n")
    skip_normalize <- TRUE
  } else {
    skip_normalize <- FALSE
  }
}

if (!skip_normalize) {
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 1e4,
    verbose = FALSE
  )
  cat("[OK] Data normalized\n")
}

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = FALSE
)

cat("[OK] Data preparation complete\n\n")
required_cols <- c("cell_type_L2", "cell_type_L3")
available_cols <- colnames(seurat_obj@meta.data)

for (col in required_cols) {
  if (!col %in% available_cols) {
    stop(sprintf("Required column '%s' not found in metadata!", col))
  }
}

cat(sprintf(
  "✓ cell_type_L2: %d unique types\n",
  length(unique(seurat_obj$cell_type_L2))
))
cat(sprintf(
  "✓ cell_type_L3: %d unique subclusters (raw)\n",
  length(unique(seurat_obj$cell_type_L3))
))

# Check for NAs
na_l2 <- sum(is.na(seurat_obj$cell_type_L2))
na_l3 <- sum(is.na(seurat_obj$cell_type_L3))

if (na_l2 > 0) {
  cat(sprintf("⚠️  Warning: %d cells have NA in cell_type_L2\n", na_l2))
}
if (na_l3 > 0) {
  cat(sprintf("⚠️  Warning: %d cells have NA in cell_type_L3\n", na_l3))
}

# Print example L3 labels (before standardization)
cat("\nExample L3 labels (raw format, will be standardized later):\n")
print(head(sort(unique(seurat_obj$cell_type_L3)), 10))
cat("\n")
cat(
  "Note: L3 labels will be converted to hierarchical format {cell_type_L2}_c{id}\n"
)
cat("      Example: 0 → 'Alveolar macrophages_c0'\n\n")

cat("\n=== Standardizing L3 Labels ===\n")

# Current L3 format: numeric IDs (0, 1, 2, 3, 4)
# Target format: {cell_type_L2}_c{subcluster_id}
# Example: "Alveolar macrophages_c0", "Classical monocytes_c1"

# Show current structure
cat("\nCurrent L3 labels (before standardization):\n")
print(table(seurat_obj$cell_type_L3, useNA = "ifany"))

# Step 1: Preserve original numeric subcluster IDs
seurat_obj$subcluster_id <- as.integer(as.character(seurat_obj$cell_type_L3))

cat(sprintf(
  "\n✓ Preserved original subcluster IDs: %d unique values\n",
  length(unique(seurat_obj$subcluster_id[!is.na(seurat_obj$subcluster_id)]))
))

# Step 2: Create hierarchical L3 labels
seurat_obj$cell_type_L3 <- ifelse(
  is.na(seurat_obj$cell_type_L2) | is.na(seurat_obj$subcluster_id),
  NA_character_,
  paste0(seurat_obj$cell_type_L2, "_c", seurat_obj$subcluster_id)
)

# Step 3: Convert to ordered factor (L2 alphabetical, then subcluster_id numeric)
tmp_df <- seurat_obj@meta.data[, c(
  "cell_type_L2",
  "subcluster_id",
  "cell_type_L3"
)]
tmp_df <- tmp_df[!is.na(tmp_df$cell_type_L3), ]
ordered_levels <- unique(
  tmp_df[order(tmp_df$cell_type_L2, tmp_df$subcluster_id), "cell_type_L3"]
)
seurat_obj$cell_type_L3 <- factor(
  seurat_obj$cell_type_L3,
  levels = ordered_levels
)

# Validation
cat("\n=== Validation Results ===\n")
cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))
cat(sprintf(
  "Cells with valid L3 labels: %d\n",
  sum(!is.na(seurat_obj$cell_type_L3))
))
cat(sprintf(
  "Unique L3 subclusters: %d\n",
  length(levels(seurat_obj$cell_type_L3))
))

# Show final structure
cat("\nFinal L3 labels (hierarchical format):\n")
l3_table <- table(seurat_obj$cell_type_L3, useNA = "ifany")
print(l3_table)

# Per-celltype summary
cat("\n=== Per-Celltype Subcluster Summary ===\n")
l2_summary <- seurat_obj@meta.data %>%
  filter(!is.na(cell_type_L2) & !is.na(cell_type_L3)) %>%
  group_by(cell_type_L2) %>%
  summarise(
    n_cells = n(),
    n_subclusters = n_distinct(subcluster_id),
    subcluster_range = sprintf(
      "c%d-c%d",
      min(subcluster_id),
      max(subcluster_id)
    ),
    .groups = "drop"
  ) %>%
  arrange(cell_type_L2)

print(l2_summary)

cat("\n")
# ==============================================================================
# Compute Marker Genes
# ==============================================================================

cat("=== Computing Marker Genes ===\n")

Idents(seurat_obj) <- SUBCLUSTER_COLUMN
clusters <- levels(Idents(seurat_obj))

cat(sprintf("Finding markers for %d myeloid cell subclusters...\n", length(clusters)))

# v1.1: Report cluster sizes
cluster_sizes <- table(Idents(seurat_obj))
cat("\nmyeloid cell subcluster sizes:\n")
print(cluster_sizes)
cat("\n")

marker_list <- future_lapply(
  clusters,
  function(cluster_id) {
    tryCatch(
      {
        # FIX 1.2: Explicitly specify slot="data"
        FindMarkers(
          seurat_obj,
          ident.1 = cluster_id,
          only.pos = TRUE,
          min.pct = 0.25,
          logfc.threshold = 0.5,
          test.use = "wilcox",
          slot = "data",
          verbose = FALSE
        )
      },
      error = function(e) {
        cat(sprintf("[WARN] FindMarkers failed for %s: %s\n",
                    cluster_id, conditionMessage(e)))
        NULL
      }
    )
  },
  future.seed = TRUE
)

names(marker_list) <- clusters
marker_list <- marker_list[!sapply(marker_list, is.null)]

# v1.1: Report skipped clusters
skipped <- setdiff(clusters, names(marker_list))
if (length(skipped) > 0) {
  cat(sprintf("[WARN] Skipped %d clusters (FindMarkers failed):\n", length(skipped)))
  cat(paste("  -", skipped, collapse = "\n"), "\n\n")
}

all_markers <- bind_rows(lapply(names(marker_list), function(cid) {
  df <- marker_list[[cid]]
  df$cluster <- cid
  df$gene <- rownames(df)
  df
}))

all_markers <- all_markers %>% filter(p_val_adj < 0.05)

cat(sprintf("[OK] Found %d significant markers across %d clusters\n",
            nrow(all_markers), length(unique(all_markers$cluster))))

write.csv(all_markers, file.path(OUTPUT_DIR, "all_markers.csv"), row.names = FALSE)

# ==============================================================================
# Filter Confounding Genes
# ==============================================================================

cat("\n=== Filtering Confounding Genes ===\n")

genes_to_filter <- c(
  grep("^MT-", rownames(seurat_obj), value = TRUE),
  grep("^RP[SL]", rownames(seurat_obj), value = TRUE),
  # Stress/IEG
  "FOS", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "ZFP36",
  "DUSP1", "DUSP2", "IER2", "IER3", "ATF3", "BTG2", "FOSB",
  "NR4A1", "NR4A2", "NR4A3",
  # Heat shock
  "HSP90AA1", "HSPA1A", "HSPA1B", "DNAJB1", "HSPA6"
)

cat(sprintf("Filtering %d potentially confounding genes\n", length(genes_to_filter)))

# v1.1 fix: Auto-detect logFC column name
lfc_col <- if ("avg_log2FC" %in% colnames(all_markers)) {
  "avg_log2FC"
} else if ("avg_logFC" %in% colnames(all_markers)) {
  "avg_logFC"
} else {
  stop("ERROR: Cannot find logFC column (tried 'avg_log2FC' and 'avg_logFC')")
}

cat(sprintf("[INFO] Using logFC column: %s\n", lfc_col))

top_markers <- all_markers %>%
  filter(!gene %in% genes_to_filter) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(.data[[lfc_col]])) %>%
  slice_head(n = TOP_N_MARKERS) %>%
  ungroup() %>%
  mutate(
    gene = toupper(gene),
    cluster = as.character(cluster)
  ) %>%
  dplyr::select(gene, cluster)

universe_genes <- setdiff(toupper(rownames(seurat_obj)), toupper(genes_to_filter))

cat(sprintf("Universe: %d genes (filtered %d)\n", 
            length(universe_genes), length(genes_to_filter)))

# v1.1: Report markers per cluster
markers_per_cluster <- top_markers %>%
  group_by(cluster) %>%
  summarise(n_markers = n())
cat("\nMarkers per cluster:\n")
print(markers_per_cluster)
cat("\n")

write.csv(top_markers, file.path(OUTPUT_DIR, "top_markers_filtered.csv"),
          row.names = FALSE)

# ==============================================================================
# Load Reference Databases
# ==============================================================================

cat("\n=== Loading Reference Databases ===\n")

# --- CellMarker ---
cellmarker_term2gene <- NULL
tryCatch({
  db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)
  db <- db %>% filter(grepl("Human", species, ignore.case = TRUE))
  
  term2gene_list <- list()
  for (i in 1:nrow(db)) {
    row <- db[i, ]
    cell_type <- row$cell_name
    markers_raw <- row$marker
    
    if (!is.na(markers_raw) && markers_raw != "") {
      markers <- unlist(strsplit(markers_raw, "[,;\\s]+"))
      markers <- gsub('["\r\n\\[\\]]', '', markers)
      markers <- unique(toupper(trimws(markers[markers != "" & !is.na(markers)])))
      
      if (length(markers) > 0) {
        term2gene_list[[length(term2gene_list) + 1]] <- data.frame(
          term = cell_type, gene = markers, stringsAsFactors = FALSE
        )
      }
    }
  }
  
  cellmarker_term2gene <- bind_rows(term2gene_list)
  cat(sprintf("[OK] CellMarker: %d pairs, %d cell types\n",
              nrow(cellmarker_term2gene), 
              length(unique(cellmarker_term2gene$term))))
}, error = function(e) {
  cat("[WARN] Failed to load CellMarker:", conditionMessage(e), "\n")
})

# --- PanglaoDB ---
panglaodb_term2gene <- NULL
tryCatch({
  db <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)
  setnames(db, old = c("official gene symbol", "cell type"),
           new = c("gene_symbol", "cell_type"), skip_absent = TRUE)
  db <- db %>% filter(grepl("Hs", species, fixed = TRUE))
  
  panglaodb_term2gene <- db %>%
    dplyr::select(cell_type, gene_symbol) %>%
    mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct() %>%
    dplyr::rename(term = cell_type, gene = gene_symbol)
  
  cat(sprintf("[OK] PanglaoDB: %d pairs, %d cell types\n",
              nrow(panglaodb_term2gene),
              length(unique(panglaodb_term2gene$term))))
}, error = function(e) {
  cat("[WARN] Failed to load PanglaoDB:", conditionMessage(e), "\n")
})

# --- Custom Markers (if provided) ---
custom_term2gene <- NULL
if (!is.null(CUSTOM_MARKERS_CSV) && file.exists(CUSTOM_MARKERS_CSV)) {
  tryCatch({
    custom_db <- read.csv(CUSTOM_MARKERS_CSV, stringsAsFactors = FALSE)
    if (all(c("term", "gene") %in% colnames(custom_db))) {
      custom_term2gene <- custom_db %>%
        mutate(
          term = as.character(term),
          gene = toupper(trimws(as.character(gene)))
        ) %>%
        filter(gene != "" & !is.na(gene)) %>%
        distinct()
      
      cat(sprintf("[OK] Custom markers: %d pairs, %d terms\n",
                  nrow(custom_term2gene),
                  length(unique(custom_term2gene$term))))
    }
  }, error = function(e) {
    cat("[WARN] Failed to load custom markers:", conditionMessage(e), "\n")
  })
}

# --- MSigDB ---
gmt_all <- NULL
tryCatch({
  gmt_all <- read.gmt(MSIGDB_GMT_PATH)
  gmt_all <- gmt_all %>% mutate(gene = toupper(gene))
  cat(sprintf("[OK] MSigDB: %d gene sets\n", length(unique(gmt_all$term))))
}, error = function(e) {
  cat("[WARN] Failed to load MSigDB GMT:", conditionMessage(e), "\n")
})

# Extract Hallmark
hallmark_term2gene <- NULL
if (!is.null(gmt_all)) {
  hallmark_term2gene <- gmt_all %>% filter(grepl("^HALLMARK_", term))
  if (nrow(hallmark_term2gene) > 0) {
    cat(sprintf("[OK] Hallmark: %d pathways\n", 
                length(unique(hallmark_term2gene$term))))
  }
}

# Extract KEGG
msigdb_kegg_term2gene <- NULL
if (!is.null(gmt_all)) {
  msigdb_kegg_term2gene <- gmt_all %>% filter(grepl("^KEGG_", term))
  if (nrow(msigdb_kegg_term2gene) > 0) {
    cat(sprintf("[OK] MSigDB KEGG: %d pathways\n",
                length(unique(msigdb_kegg_term2gene$term))))
  }
}

# --- GO Terms (GMT-based, zero gene loss) ---
go_bp_term2gene <- NULL
go_mf_term2gene <- NULL
go_cc_term2gene <- NULL

tryCatch({
  gmt_go <- read.gmt(GMT_GO_ALL)
  gmt_go <- gmt_go %>% mutate(gene = toupper(gene))
  
  go_bp_term2gene <- gmt_go %>% filter(grepl("^GOBP_", term))
  go_mf_term2gene <- gmt_go %>% filter(grepl("^GOMF_", term))
  go_cc_term2gene <- gmt_go %>% filter(grepl("^GOCC_", term))
  
  cat(sprintf("[OK] GO terms: BP=%d, MF=%d, CC=%d\n",
              length(unique(go_bp_term2gene$term)),
              length(unique(go_mf_term2gene$term)),
              length(unique(go_cc_term2gene$term))))
}, error = function(e) {
  cat("[WARN] Failed to load GO GMT:", conditionMessage(e), "\n")
})

# ==============================================================================
# Load Existing Enrichment Results (Optional)
# ==============================================================================

if (LOAD_EXISTING_ENRICHMENT) {
  cat("\n=== Loading Existing Enrichment Results ===\n")
  
  enrich_dir <- if (!is.null(ENRICHMENT_DIR)) {
    ENRICHMENT_DIR
  } else {
    file.path(OUTPUT_DIR, "reports")
  }
  
  if (!dir.exists(enrich_dir)) {
    cat(sprintf("[WARN] Enrichment directory not found: %s\n", enrich_dir))
    cat("[INFO] Will run enrichment from scratch\n")
    LOAD_EXISTING_ENRICHMENT <- FALSE
  } else {
    # List of expected enrichment files
    enrich_files <- list(
      custom = "custom_enrich.rds",
      cellmarker = "cellmarker_enrich.rds",
      panglaodb = "panglaodb_enrich.rds",
      go_bp = "go_bp_enrich.rds",
      go_mf = "go_mf_enrich.rds",
      go_cc = "go_cc_enrich.rds",
      hallmark = "hallmark_enrich.rds",
      kegg = "kegg_enrich.rds"
    )
    
    enrichment_results <- list()
    
    for (db_name in names(enrich_files)) {
      file_path <- file.path(enrich_dir, enrich_files[[db_name]])
      
      if (file.exists(file_path)) {
        tryCatch({
          loaded_obj <- readRDS(file_path)
          enrichment_results[[db_name]] <- loaded_obj
          cat(sprintf("[OK] Loaded %s\n", enrich_files[[db_name]]))
        }, error = function(e) {
          cat(sprintf("[WARN] Failed to load %s: %s\n", 
                     enrich_files[[db_name]], conditionMessage(e)))
        })
      } else {
        cat(sprintf("[INFO] Not found: %s\n", enrich_files[[db_name]]))
      }
    }
    
    n_loaded <- length(enrichment_results)
    cat(sprintf("\n[SUMMARY] Loaded %d enrichment result(s)\n", n_loaded))
    
    if (n_loaded == 0) {
      cat("[INFO] No enrichment results loaded, will run from scratch\n")
      LOAD_EXISTING_ENRICHMENT <- FALSE
      enrichment_results <- list()
    } else {
      cat("[INFO] Using loaded enrichment results, skipping computation\n")
      RUN_ENRICHMENT <- FALSE  # Skip running enrichment
    }
  }
}

# ==============================================================================
# Traditional Enrichment Analysis (compareCluster)
# ==============================================================================

if (!exists("enrichment_results")) {
  enrichment_results <- list()
}

if (RUN_ENRICHMENT) {
  cat("\n=== Running Traditional Enrichment Analysis (compareCluster) ===\n")
  
  run_compareCluster_enrichment <- function(term2gene_df, db_name) {
    if (is.null(term2gene_df) || nrow(term2gene_df) == 0) {
      return(NULL)
    }
    
    cat(sprintf("\n--- Database: %s ---\n", db_name))
    
    # FIX 1.4: Use smaller minGSSize for cell type marker databases
    # CellMarker/PanglaoDB often have small but highly specific gene sets
    min_size <- if (db_name %in% c("CellMarker", "PanglaoDB", "Custom")) 5 else 10
    max_size <- 500
    
    result <- tryCatch({
      compareCluster(
        gene ~ cluster,
        data = top_markers,
        fun = enricher,
        TERM2GENE = term2gene_df,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2,
        minGSSize = min_size,
        maxGSSize = max_size
      )
    }, error = function(e) {
      cat("[WARN]", db_name, "enrichment failed:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(result)) {
      ccr <- result@compareClusterResult
      n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
      cat(sprintf("[OK] %s: Found %d significant terms (minGSSize=%d)\n", 
                  db_name, n_sig, min_size))
      
      saveRDS(result, file.path(OUTPUT_DIR, "reports", 
                                paste0(tolower(gsub(" ", "_", db_name)), "_enrich.rds")))
    }
    
    return(result)
  }
  
  # Run enrichment for each database
  enrichment_results$custom <- run_compareCluster_enrichment(custom_term2gene, "Custom")
  enrichment_results$cellmarker <- run_compareCluster_enrichment(cellmarker_term2gene, "CellMarker")
  enrichment_results$panglaodb <- run_compareCluster_enrichment(panglaodb_term2gene, "PanglaoDB")
  enrichment_results$go_bp <- run_compareCluster_enrichment(go_bp_term2gene, "GO BP")
  enrichment_results$go_mf <- run_compareCluster_enrichment(go_mf_term2gene, "GO MF")
  enrichment_results$go_cc <- run_compareCluster_enrichment(go_cc_term2gene, "GO CC")
  enrichment_results$hallmark <- run_compareCluster_enrichment(hallmark_term2gene, "Hallmark")
  enrichment_results$kegg <- run_compareCluster_enrichment(msigdb_kegg_term2gene, "KEGG")
  
  # Remove NULL entries
  enrichment_results <- enrichment_results[!sapply(enrichment_results, is.null)]
  
  cat("\n[OK] Traditional enrichment complete\n")
}

# ==============================================================================
# Pseudobulk ssGSEA (Primary Evidence) - ALL METHODS
# ==============================================================================

ssgsea_results_all <- list()  # Store all ssGSEA results by method

if (RUN_SSGSEA) {
  cat("\n=== Running Pseudobulk ssGSEA (ALL METHODS) ===\n")
  
  # Helper functions
  term2gene_to_list <- function(term2gene_df) {
    split(toupper(term2gene_df$gene), term2gene_df$term)
  }
  
  map_gene_sets_to_features <- function(gene_sets, features) {
    feats_upper <- toupper(features)
    keep <- !duplicated(feats_upper)
    m <- setNames(features[keep], feats_upper[keep])
    gs2 <- lapply(gene_sets, function(g) unique(na.omit(m[toupper(g)])))
    gs2[lengths(gs2) > 0]
  }
  
  filter_gs_size <- function(gene_sets, minSize = 10, maxSize = 500) {
    gene_sets[lengths(gene_sets) >= minSize & lengths(gene_sets) <= maxSize]
  }
  
  # ⭐ Loop through ALL ssGSEA methods
  for (method in SSGSEA_METHODS) {
    cat(sprintf("\n--- ssGSEA Method: %s ---\n", method))
    
    # Select gene sets based on method
    selected_term2gene <- switch(
      method,
      "hallmark" = hallmark_term2gene,
      "go_bp" = go_bp_term2gene,
      "go_mf" = go_mf_term2gene,
      "go_cc" = go_cc_term2gene,
      "kegg" = msigdb_kegg_term2gene,
      "custom" = if (!is.null(SSGSEA_CUSTOM_GMT) && file.exists(SSGSEA_CUSTOM_GMT)) {
        read.gmt(SSGSEA_CUSTOM_GMT) %>% mutate(gene = toupper(gene))
      } else {
        NULL
      },
      NULL
    )
    
    if (is.null(selected_term2gene) || nrow(selected_term2gene) == 0) {
      cat(sprintf("[WARN] Method '%s': No gene sets available, skipping\n", method))
      next
    }
    
    # Prepare gene sets
    features <- rownames(seurat_obj[["RNA"]])
    ssgsea_gene_sets <- term2gene_to_list(selected_term2gene)
    ssgsea_gene_sets <- map_gene_sets_to_features(ssgsea_gene_sets, features)
    ssgsea_gene_sets <- filter_gs_size(ssgsea_gene_sets, minSize = 10, maxSize = 500)
    
    if (length(ssgsea_gene_sets) == 0) {
      cat(sprintf("[WARN] Method '%s': No valid gene sets after filtering, skipping\n", method))
      next
    }
    
    cat(sprintf("[OK] Prepared %d gene sets\n", length(ssgsea_gene_sets)))
    
    # Compute pseudobulk expression (only for genes needed by this method)
    needed_genes <- unique(unlist(ssgsea_gene_sets, use.names = FALSE))
    
    cat(sprintf("[INFO] Computing pseudobulk for %d genes...\n", length(needed_genes)))
    avg_expr <- tryCatch({
      Seurat::AverageExpression(
        seurat_obj,
        assays = "RNA",
        slot = "data",
        group.by = SUBCLUSTER_COLUMN,
        features = needed_genes,
        verbose = FALSE
      )[["RNA"]]
    }, error = function(e) {
      cat(sprintf("[ERROR] Pseudobulk failed: %s\n", conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(avg_expr)) {
      cat(sprintf("[WARN] Method '%s': Pseudobulk failed, skipping\n", method))
      next
    }
    
    # Run ssGSEA
    cat("[INFO] Running ssGSEA...\n")
    bp <- BiocParallel::SnowParam(workers = N_CORES, type = "SOCK", progressbar = FALSE)
    
    ssgsea_scores <- tryCatch({
      if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
        param <- GSVA::ssgseaParam(
          exprData = as.matrix(avg_expr),
          geneSets = ssgsea_gene_sets,
          alpha = 0.25,
          normalize = TRUE,
          minSize = 10,
          maxSize = 500
        )
        GSVA::gsva(param, BPPARAM = bp, verbose = FALSE)
      } else {
        GSVA::gsva(
          as.matrix(avg_expr),
          ssgsea_gene_sets,
          method = "ssgsea",
          ssgsea.norm = TRUE,
          verbose = FALSE
        )
      }
    }, error = function(e) {
      cat(sprintf("[ERROR] ssGSEA failed: %s\n", conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(ssgsea_scores)) {
      cat(sprintf("[WARN] Method '%s': ssGSEA computation failed, skipping\n", method))
      next
    }
    
    # Compute Z-scores
    ssgsea_z <- t(scale(t(ssgsea_scores)))
    
    cat(sprintf("[OK] ssGSEA complete: %d pathways x %d clusters\n",
                nrow(ssgsea_scores), ncol(ssgsea_scores)))
    
    # Store results
    ssgsea_results_all[[method]] <- list(
      scores = ssgsea_scores,
      z_scores = ssgsea_z,
      gene_sets = ssgsea_gene_sets
    )
    
    # Save to disk
    saveRDS(ssgsea_scores, 
            file.path(OUTPUT_DIR, "reports", sprintf("ssgsea_scores_%s.rds", method)))
    saveRDS(ssgsea_z,
            file.path(OUTPUT_DIR, "reports", sprintf("ssgsea_z_%s.rds", method)))
    
    # Export top pathways per cluster
    ssgsea_top_per_cluster <- lapply(colnames(ssgsea_scores), function(cid) {
      sc <- ssgsea_scores[, cid]
      z <- ssgsea_z[, cid]
      ord <- order(z, sc, decreasing = TRUE, na.last = TRUE)
      top_idx <- ord[1:min(SSGSEA_N_TOP * 1.5, length(ord))]
      data.frame(
        method = method,
        cluster = cid,
        pathway = names(sc)[top_idx],
        score = sc[top_idx],
        z_score = z[top_idx],
        rank = 1:length(top_idx),
        stringsAsFactors = FALSE
      )
    })
    
    ssgsea_top_df <- bind_rows(ssgsea_top_per_cluster)
    write.csv(ssgsea_top_df,
              file.path(OUTPUT_DIR, "reports", sprintf("ssgsea_top_pathways_%s.csv", method)),
              row.names = FALSE)
    
    cat(sprintf("[OK] Method '%s' results saved\n", method))
  }
  
  # Summary
  cat(sprintf("\n[SUMMARY] Completed %d ssGSEA methods\n", length(ssgsea_results_all)))
  for (method_name in names(ssgsea_results_all)) {
    n_pathways <- nrow(ssgsea_results_all[[method_name]]$scores)
    cat(sprintf("  - %s: %d pathways\n", method_name, n_pathways))
  }
  
  # ⭐ Use primary method for downstream (first successful method, or hallmark if available)
  primary_method <- if ("hallmark" %in% names(ssgsea_results_all)) {
    "hallmark"
  } else {
    names(ssgsea_results_all)[1]
  }
  
  if (!is.null(primary_method)) {
    cat(sprintf("\n[INFO] Using '%s' as primary method for LLM context\n", primary_method))
    ssgsea_scores <- ssgsea_results_all[[primary_method]]$scores
    ssgsea_z <- ssgsea_results_all[[primary_method]]$z_scores
    ssgsea_gene_sets <- ssgsea_results_all[[primary_method]]$gene_sets
  } else {
    cat("[WARN] No ssGSEA methods succeeded\n")
    ssgsea_scores <- NULL
    ssgsea_z <- NULL
    ssgsea_gene_sets <- NULL
  }
}

# ==============================================================================
# Myeloid Cell Expert Context Generator
# ==============================================================================

generate_myeloid_context <- function() {
  paste(
    "CRITICAL CONTEXT FOR MYELOID CELL ANNOTATION:",
    "",
    "These are myeloid cell subclusters from respiratory tissues (normal/disease study).",
    "",
    "=== LINEAGE GATE VERIFICATION (Important) ===",
    "Before finalizing annotation, verify lineage markers:",
    "• Myeloid: LYZ+, CD68+, CD14+ (monocyte/macrophage)",
    "• DC lineage: FCER1A+, CD1C+ (cDC2), XCR1+/CLEC9A+ (cDC1), LILRA4+ (pDC)",
    "• Contamination check:",
    "  - T cells: CD3D/CD3E (if dominant → T cell contamination)",
    "  - B cells: CD79A/MS4A1 (if dominant → B cell contamination)",
    "  - Epithelial: EPCAM/KRT (if dominant → epithelial contamination)",
    "",
    "⚠️ State vs Identity Warning:",
    "If markers are mainly OXPHOS (ATP5*, COX*, ND*) or stress (HSP*, FOS, JUN),",
    "note this is likely a METABOLIC/ACTIVATION STATE, not cell identity.",
    "",
    "Key myeloid biology principles:",
    "",
    "1. MONOCYTE STATES:",
    "   - Classical monocytes: CD14++, CD16-, S100A8+, S100A9+, FCN1+",
    "     * Circulating, pro-inflammatory potential",
    "     * Can differentiate into macrophages/DCs in tissues",
    "   - Intermediate monocytes: CD14+, CD16+, HLA-DR+",
    "     * Inflammatory, antigen presentation",
    "   - Non-classical monocytes: CD14+, CD16++, FCGR3A+",
    "     * Patrolling, vascular surveillance",
    "     * CX3CR1+, SELL low",
    "",
    "2. MACROPHAGE POLARIZATION (CRITICAL):",
    "   - M1 (Pro-inflammatory): CXCL9+, CXCL10+, IL1B+, TNF+, NOS2+",
    "     * Classical activation, pathogen killing",
    "     * High MHC-II, co-stimulatory molecules",
    "   - M2 (Anti-inflammatory/Tissue repair): CD163+, MRC1+, MSR1+",
    "     * Alternative activation subtypes:",
    "       - M2a: IL4/IL13 induced, CD163+, CD209+",
    "       - M2b: Immune complex induced, IL10+",
    "       - M2c: IL10/TGFb induced, tissue remodeling",
    "   - Mixed/Intermediate: Co-expression of M1/M2 markers",
    "   - Tissue-resident macrophages (TRMs):",
    "     * Alveolar macrophages: FABP4+, MARCO+",
    "     * Interstitial macrophages: LYVE1+, FOLR2+",
    "",
    "3. DENDRITIC CELLS:",
    "   - cDC1 (Type 1 conventional DC): XCR1+, CLEC9A+, IRF8+",
    "     * Cross-presentation, anti-viral immunity",
    "     * BATF3-dependent development",
    "   - cDC2 (Type 2 conventional DC): CD1C+, FCER1A+, IRF4+",
    "     * Th2/Th17 priming, bacterial response",
    "   - pDC (Plasmacytoid DC): LILRA4+, IRF7+, IL3RA+",
    "     * Type I interferon production",
    "     * Anti-viral immunity",
    "   - Maturation states:",
    "     * Immature: Low CD80/CD86, high endocytosis",
    "     * Mature: High CD80/CD86/CD83, CCR7+, IL12+",
    "     * Tolerogenic: PDCD1LG1+, IDO1+",
    "",
    "4. MAST CELLS:",
    "   - KIT+, TPSAB1+, CPA3+, HPGDS+",
    "   - Granule-associated: TRYPTASE, CHYMASE, HISTAMINE",
    "   - IgE receptor: FCER1A+",
    "   - Activated: CPA3+, IL1RL1+, IL33 responsive",
    "",
    "5. FUNCTIONAL STATES:",
    "   - Resting: Low activation markers, tissue homeostasis",
    "   - Activated: CD86+, CD80+, IL1B+, TNF+",
    "   - Phagocytic: High CD68, MSR1, MARCO",
    "   - Antigen-presenting: High HLA-DR/DP/DQ, CD86+",
    "   - Proliferating: MKI67+, TOP2A+, STMN1+",
    "   - Inflammatory: CXCL8+, CCL2+, IL1B+, NFKB pathway",
    "",
    "6. TISSUE CONTEXT (Respiratory tissues):",
    "   - Alveolar macrophages dominate normal lung",
    "   - Interstitial macrophages in connective tissue",
    "   - Monocyte-derived cells increase in inflammation",
    "   - DCs enriched at mucosal surfaces",
    "",
    "ANNOTATION INSTRUCTIONS:",
    "For EACH cluster, specify:",
    "  a) Cell type: Monocyte / Macrophage / DC / Mast cell",
    "  b) Subtype:",
    "     - Monocytes: Classical / Intermediate / Non-classical",
    "     - Macrophages: M1 / M2 / Mixed / Alveolar / Interstitial",
    "     - DCs: cDC1 / cDC2 / pDC / Mature / Immature",
    "  c) Functional state: Resting / Activated / Inflammatory / Phagocytic",
    "  d) Tissue context: if applicable (e.g., alveolar, interstitial)",
    "",
    "Example good annotations:",
    "  - 'Classical monocytes (CD14++, S100A8+, FCN1+)'",
    "  - 'M2-polarized alveolar macrophages (CD163+, FABP4+, tissue-resident)'",
    "  - 'Mature cDC2 dendritic cells (CD1C+, CD86+, antigen-presenting)'",
    "  - 'Inflammatory M1 macrophages (CXCL10+, IL1B+, TNF+)'",
    "",
    "Always base reasoning on:",
    "  1. Specific marker gene expression",
    "  2. Enriched pathways (Hallmark/KEGG/GO for activation states)",
    "  3. ssGSEA evidence (if available)",
    "  4. Biological plausibility in respiratory tissue",
    sep = "\n"
  )
}

myeloid_context <- generate_myeloid_context()

# ==============================================================================
# LLM Interpretation Functions
# ==============================================================================

# Standard interpret() wrapper for T cells
interpret_tcell <- function(cluster_id, marker_genes, enrichment_data = NULL,
                           ssgsea_data = NULL, task = "annotation",
                           api_key = DEEPSEEK_API_KEY) {
  cat(sprintf("\n--- Myeloid Cell Interpretation: %s (Task: %s) ---\n",
              cluster_id, task))
  
  markers_text <- paste(head(marker_genes, 30), collapse = ", ")
  
  # Enrichment summary
  enrichment_text <- ""
  if (!is.null(enrichment_data) && is.list(enrichment_data)) {
    enrichment_summaries <- c()
    for (db_name in names(enrichment_data)) {
      db_results <- enrichment_data[[db_name]]
      if (!is.null(db_results) && nrow(db_results) > 0) {
        top_terms <- head(db_results$Description, 5)
        enrichment_summaries <- c(
          enrichment_summaries,
          sprintf("%s: %s", db_name, paste(top_terms, collapse = "; "))
        )
      }
    }
    if (length(enrichment_summaries) > 0) {
      enrichment_text <- paste("\n\nEnrichment Results:",
                               paste(enrichment_summaries, collapse = "\n"),
                               sep = "\n")
    }
  }
  
  # ssGSEA summary
  ssgsea_text <- ""
  if (!is.null(ssgsea_data) && nrow(ssgsea_data) > 0) {
    # FIX 1.3: Use SSGSEA_N_TOP parameter
    top_pathways <- ssgsea_data %>%
      arrange(desc(z_score), desc(score)) %>%
      head(SSGSEA_N_TOP)
    
    pathway_lines <- sprintf("  %s (Z=%.2f, Score=%.3f)",
                            top_pathways$pathway,
                            top_pathways$z_score,
                            top_pathways$score)
    
    ssgsea_text <- paste("\n\nssGSEA Top Pathways (by specificity):",
                        paste(pathway_lines, collapse = "\n"),
                        sep = "\n")
  }
  
  # Construct prompt
  if (task == "annotation") {
    prompt <- sprintf(
      "%s

Based on marker genes, enrichment, and ssGSEA evidence, provide:

Cell Type: [specific myeloid cell subtype]
Confidence: [High/Medium/Low]
Regulatory Drivers: [driver1; driver2; driver3]
Markers: [marker1; marker2; marker3]
Reasoning: [brief explanation]

Marker Genes: %s%s%s

CRITICAL: Include ALL FIVE fields with exact format.",
      myeloid_context, markers_text, enrichment_text, ssgsea_text
    )
  } else if (task == "phenotype") {
    prompt <- sprintf(
      "%s

Based on marker genes, enrichment, and ssGSEA evidence, provide:

Functional Phenotype: [dominant functional state]
Confidence: [High/Medium/Low]
Regulatory Drivers: [driver1; driver2; driver3]
Key Processes: [process1; process2; process3]
Reasoning: [brief explanation]
Network Evidence: [network observations]

Marker Genes: %s%s%s

CRITICAL: Include ALL SIX fields with exact format.",
      myeloid_context, markers_text, enrichment_text, ssgsea_text
    )
  } else {
    stop("Invalid task type. Must be 'annotation' or 'phenotype'")
  }
  
  # Call API
  result <- tryCatch({
    # FIX 0-A: Use body as list with encode="json" to avoid double JSON encoding
    payload <- list(
      model = "deepseek-reasoner",
      messages = list(list(role = "user", content = prompt)),
      temperature = 0.3
    )
    
    response <- httr::POST(
      url = "https://api.deepseek.com/v1/chat/completions",
      httr::add_headers(
        "Authorization" = paste("Bearer", api_key),
        "Content-Type" = "application/json"
      ),
      body = payload,
      encode = "json",
      httr::timeout(120)
    )
    
    if (httr::status_code(response) == 200) {
      content <- httr::content(response, as = "parsed")
      interpretation <- content$choices[[1]]$message$content
      cat(sprintf("[OK] Received interpretation (%d chars)\n", nchar(interpretation)))
      
      # Save raw prompt and response for debugging
      raw_dir <- file.path(dirname(dirname(getwd())), "reports", "llm_raw")
      if (!dir.exists(raw_dir)) {
        raw_dir <- file.path(getwd(), "reports", "llm_raw")
      }
      if (dir.exists(raw_dir)) {
        raw_file <- file.path(raw_dir, sprintf("%s_%s.txt", cluster_id, task))
        tryCatch({
          writeLines(c(
            "=== PROMPT ===",
            prompt,
            "",
            "=== RESPONSE ===",
            interpretation
          ), raw_file)
        }, error = function(e) NULL)
      }
      
      interpretation
    } else {
      cat(sprintf("[ERROR] API returned status %d\n", httr::status_code(response)))
      NULL
    }
  }, error = function(e) {
    cat(sprintf("[ERROR] API call failed: %s\n", conditionMessage(e)))
    NULL
  })
  
  return(result)
}

# Parse interpretation
parse_interpretation <- function(text, task = "annotation") {
  if (is.null(text) || nchar(text) == 0 || is.na(text)) {
    return(list(raw = NA_character_))
  }
  
  # FIX 0-C: Use DOTALL mode (?s) to capture multi-line fields
  extract_field <- function(text, field_name) {
    pattern <- sprintf("(?s)%s:\\s*(.+?)(?=\\n[A-Za-z][A-Za-z ]*:|$)", field_name)
    m <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
    if (length(m) >= 2) {
      trimws(m[2])
    } else {
      NA_character_
    }
  }
  
  parsed <- list(raw = text)
  
  if (task == "annotation") {
    parsed$cell_type <- extract_field(text, "Cell Type")
    parsed$confidence <- extract_field(text, "Confidence")
    parsed$regulatory_drivers <- extract_field(text, "Regulatory Drivers")
    parsed$markers <- extract_field(text, "Markers")
    parsed$reasoning <- extract_field(text, "Reasoning")
  } else if (task == "phenotype") {
    parsed$functional_phenotype <- extract_field(text, "Functional Phenotype")
    parsed$confidence <- extract_field(text, "Confidence")
    parsed$regulatory_drivers <- extract_field(text, "Regulatory Drivers")
    parsed$key_processes <- extract_field(text, "Key Processes")
    parsed$reasoning <- extract_field(text, "Reasoning")
    parsed$network_evidence <- extract_field(text, "Network Evidence")
  }
  
  return(parsed)
}

# ==============================================================================
# Run Standard LLM Interpretations
# ==============================================================================

cat("\n=== Running Standard LLM Interpretations ===\n")

old_plan <- future::plan()
future::plan("sequential")  # Avoid rate limiting

clusters_for_interpretation <- unique(top_markers$cluster)

# Annotation task
cat("\n--- Annotation Task ---\n")
annotation_results <- future_lapply(
  clusters_for_interpretation,
  function(cluster_id) {
    markers <- top_markers %>% filter(cluster == cluster_id) %>% pull(gene)
    
    # Get enrichment for this cluster
    cluster_enrichment <- lapply(enrichment_results, function(db_result) {
      if (!is.null(db_result) && "compareClusterResult" %in% class(db_result)) {
        ccr <- db_result@compareClusterResult
        ccr %>% filter(Cluster == cluster_id)
      } else {
        NULL
      }
    })
    cluster_enrichment <- cluster_enrichment[sapply(cluster_enrichment, 
                                                    function(x) !is.null(x) && nrow(x) > 0)]
    
    # Get ssGSEA for this cluster
    cluster_ssgsea <- NULL
    if (!is.null(ssgsea_scores) && cluster_id %in% colnames(ssgsea_scores)) {
      ssgsea_df <- data.frame(
        pathway = rownames(ssgsea_scores),
        score = ssgsea_scores[, cluster_id],
        z_score = ssgsea_z[, cluster_id],
        stringsAsFactors = FALSE
      )
      cluster_ssgsea <- ssgsea_df %>% arrange(desc(z_score), desc(score)) %>% head(20)
    }
    
    interpretation <- interpret_tcell(
      cluster_id = cluster_id,
      marker_genes = markers,
      enrichment_data = cluster_enrichment,
      ssgsea_data = cluster_ssgsea,
      task = "annotation"
    )
    
    parsed <- parse_interpretation(interpretation, task = "annotation")
    parsed$cluster <- cluster_id
    return(parsed)
  },
  future.seed = TRUE
)
names(annotation_results) <- clusters_for_interpretation

annotation_df <- bind_rows(lapply(annotation_results, function(x) {
  data.frame(
    Cluster = x$cluster,
    Cell_Type = x$cell_type %||% NA_character_,
    Confidence = x$confidence %||% NA_character_,
    Regulatory_Drivers = x$regulatory_drivers %||% NA_character_,
    Markers = x$markers %||% NA_character_,
    Reasoning = x$reasoning %||% NA_character_,
    stringsAsFactors = FALSE
  )
}))

write.csv(annotation_df, file.path(OUTPUT_DIR, "annotation_results.csv"),
          row.names = FALSE)
saveRDS(annotation_results, file.path(OUTPUT_DIR, "reports", "annotation_results.rds"))

cat(sprintf("[OK] Annotation complete: %d/%d successful\n",
            sum(!is.na(annotation_df$Cell_Type)), nrow(annotation_df)))

# Phenotype task
cat("\n--- Phenotype Task ---\n")
phenotype_results <- future_lapply(
  clusters_for_interpretation,
  function(cluster_id) {
    markers <- top_markers %>% filter(cluster == cluster_id) %>% pull(gene)
    
    cluster_enrichment <- lapply(enrichment_results, function(db_result) {
      if (!is.null(db_result) && "compareClusterResult" %in% class(db_result)) {
        ccr <- db_result@compareClusterResult
        ccr %>% filter(Cluster == cluster_id)
      } else {
        NULL
      }
    })
    cluster_enrichment <- cluster_enrichment[sapply(cluster_enrichment,
                                                    function(x) !is.null(x) && nrow(x) > 0)]
    
    cluster_ssgsea <- NULL
    if (!is.null(ssgsea_scores) && cluster_id %in% colnames(ssgsea_scores)) {
      ssgsea_df <- data.frame(
        pathway = rownames(ssgsea_scores),
        score = ssgsea_scores[, cluster_id],
        z_score = ssgsea_z[, cluster_id],
        stringsAsFactors = FALSE
      )
      cluster_ssgsea <- ssgsea_df %>% arrange(desc(z_score), desc(score)) %>% head(20)
    }
    
    interpretation <- interpret_tcell(
      cluster_id = cluster_id,
      marker_genes = markers,
      enrichment_data = cluster_enrichment,
      ssgsea_data = cluster_ssgsea,
      task = "phenotype"
    )
    
    parsed <- parse_interpretation(interpretation, task = "phenotype")
    parsed$cluster <- cluster_id
    return(parsed)
  },
  future.seed = TRUE
)
names(phenotype_results) <- clusters_for_interpretation

phenotype_df <- bind_rows(lapply(phenotype_results, function(x) {
  data.frame(
    Cluster = x$cluster,
    Functional_Phenotype = x$functional_phenotype %||% NA_character_,
    Confidence = x$confidence %||% NA_character_,
    Regulatory_Drivers = x$regulatory_drivers %||% NA_character_,
    Key_Processes = x$key_processes %||% NA_character_,
    Reasoning = x$reasoning %||% NA_character_,
    Network_Evidence = x$network_evidence %||% NA_character_,
    stringsAsFactors = FALSE
  )
}))

write.csv(phenotype_df, file.path(OUTPUT_DIR, "phenotype_results.csv"),
          row.names = FALSE)
saveRDS(phenotype_results, file.path(OUTPUT_DIR, "reports", "phenotype_results.rds"))

cat(sprintf("[OK] Phenotype complete: %d/%d successful\n",
            sum(!is.na(phenotype_df$Functional_Phenotype)), nrow(phenotype_df)))

future::plan(old_plan)

# ==============================================================================
# Configure DeepSeek API
# ==============================================================================

cat("\n=== Configuring DeepSeek API ===\n")

# 检查 API key
if (is.null(DEEPSEEK_API_KEY) || DEEPSEEK_API_KEY == "") {
  stop("DEEPSEEK_API_KEY not set! Please set: export DEEPSEEK_API_KEY='your-key'")
}

# ⭐ 关键：正确设置 fanyi 的 translate option
fanyi::set_translate_option(
  key = DEEPSEEK_API_KEY,
  source = "deepseek"
)

cat("[OK] DeepSeek API configured\n\n")


source('/home/h2048/script/R/interpret.R')
source('/home/h2048/script/R/interpret_agent_hotfix.R')
# ==============================================================================
# interpret_agent Integration (Deep Mode)
# ==============================================================================

if (RUN_INTERPRET_AGENT && length(enrichment_results) > 0) {
  cat("\n=== Deep Mode: interpret_agent ===\n")
  
  # Build global gene fold change vector
  fc_tbl <- all_markers %>%
    filter(!gene %in% genes_to_filter) %>%
    mutate(gene = toupper(gene)) %>%
    group_by(gene) %>%
    summarise(fc = .data[[lfc_col]][which.max(abs(.data[[lfc_col]]))], .groups = "drop")
  
  gene_fc <- fc_tbl$fc
  names(gene_fc) <- fc_tbl$gene
  
  # Combine TERM2GENE with DB prefixes
  term2gene_blocks <- list()
  
  if (!is.null(go_bp_term2gene) && nrow(go_bp_term2gene) > 0) {
    term2gene_blocks$GO_BP <- go_bp_term2gene %>%
      mutate(term = paste0("GO_BP|", term))
  }
  if (!is.null(go_mf_term2gene) && nrow(go_mf_term2gene) > 0) {
    term2gene_blocks$GO_MF <- go_mf_term2gene %>%
      mutate(term = paste0("GO_MF|", term))
  }
  if (!is.null(go_cc_term2gene) && nrow(go_cc_term2gene) > 0) {
    term2gene_blocks$GO_CC <- go_cc_term2gene %>%
      mutate(term = paste0("GO_CC|", term))
  }
  if (!is.null(hallmark_term2gene) && nrow(hallmark_term2gene) > 0) {
    term2gene_blocks$HALLMARK <- hallmark_term2gene %>%
      mutate(term = paste0("HALLMARK|", term))
  }
  if (!is.null(msigdb_kegg_term2gene) && nrow(msigdb_kegg_term2gene) > 0) {
    term2gene_blocks$KEGG <- msigdb_kegg_term2gene %>%
      mutate(term = paste0("KEGG|", term))
  }
  if (!is.null(cellmarker_term2gene) && nrow(cellmarker_term2gene) > 0) {
    term2gene_blocks$CELLMARKER <- cellmarker_term2gene %>%
      mutate(term = paste0("CELLMARKER|", term))
  }
  if (!is.null(panglaodb_term2gene) && nrow(panglaodb_term2gene) > 0) {
    term2gene_blocks$PANGLAODB <- panglaodb_term2gene %>%
      mutate(term = paste0("PANGLAODB|", term))
  }
  
  combined_term2gene <- bind_rows(term2gene_blocks) %>%
    mutate(term = as.character(term), gene = toupper(gene)) %>%
    distinct()
  
  if (nrow(combined_term2gene) == 0) {
    cat("[WARN] Combined TERM2GENE empty, skip interpret_agent\n")
  } else {
    # FIX 1.5: Filter to genes present in expression matrix to avoid mismatches
    available_genes <- toupper(rownames(seurat_obj))
    combined_term2gene <- combined_term2gene %>% filter(gene %in% available_genes)
    combined_universe <- unique(combined_term2gene$gene)
    
    cat(sprintf("[OK] Combined TERM2GENE: %d pairs, %d genes (filtered to available)\n",
                nrow(combined_term2gene), length(combined_universe)))
    
    agent_results <- list()
    agent_rows <- list()
    
    # Helper function
    to_scalar <- function(x) {
      if (is.null(x) || length(x) == 0) return("")
      if (is.character(x) && length(x) == 1) return(x)
      if (is.character(x)) return(paste(x, collapse = "; "))
      if (is.list(x)) return(paste(capture.output(str(x, max.level = 2)), collapse = "\n"))
      as.character(x)
    }
    
    clusters_for_agent <- sort(unique(top_markers$cluster))
    
    for (cid in clusters_for_agent) {
      cat(sprintf("\n[interpret_agent] Cluster: %s\n", cid))
      
      genes_c <- top_markers %>%
        filter(cluster == cid) %>% pull(gene) %>% unique()
      
      if (length(genes_c) < 10) {
        cat(sprintf("[WARN] Too few genes (%d), skip\n", length(genes_c)))
        next
      }
      
      # Per-cluster enrichment
      er <- tryCatch({
        enricher(
          gene = genes_c,
          TERM2GENE = combined_term2gene,
          universe = combined_universe,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          minGSSize = 10,
          maxGSSize = 500
        )
      }, error = function(e) {
        cat("[WARN] enricher failed:", conditionMessage(e), "\n")
        NULL
      })
      
      if (is.null(er) || is.null(er@result) || nrow(er@result) == 0) {
        cat("[WARN] No enriched terms, skip interpret_agent\n")
        next
      }
      
      # Run interpret_agent
      ia <- tryCatch({
        interpret_agent(
          x = er,
          context = myeloid_context,
          n_pathways = INTERPRET_AGENT_N_PATHWAYS,
          model = INTERPRET_AGENT_MODEL,
          api_key = DEEPSEEK_API_KEY,
          add_ppi = INTERPRET_AGENT_ADD_PPI,
          gene_fold_change = gene_fc
        )
      }, error = function(e) {
        cat("[ERROR] interpret_agent failed:", conditionMessage(e), "\n")
        NULL
      })
      
      if (is.null(ia)) next
      
      # Normalize result
      res_one <- if (inherits(ia, "interpretation_list") && length(ia) >= 1) {
        ia[[1]]
      } else {
        ia
      }
      if (is.list(res_one)) res_one$cluster <- cid
      
      agent_results[[cid]] <- res_one
      
      # Save per-cluster txt
      txt <- capture.output(print(res_one))
      writeLines(txt, file.path(OUTPUT_DIR, "reports", "interpret_agent_txt",
                               paste0(cid, ".txt")))
      
      # Build CSV row
      agent_rows[[length(agent_rows) + 1]] <- data.frame(
        Cluster = cid,
        Overview = to_scalar(res_one$overview),
        Key_Mechanisms = to_scalar(res_one$key_mechanisms),
        Regulatory_Drivers = to_scalar(res_one$regulatory_drivers),
        Narrative = to_scalar(res_one$narrative),
        stringsAsFactors = FALSE
      )
      
      cat("[OK] interpret_agent done\n")
    }
    
    # Save all
    saveRDS(agent_results, file.path(OUTPUT_DIR, "reports", "interpret_agent_results.rds"))
    
    if (length(agent_rows) > 0) {
      agent_df <- bind_rows(agent_rows)
      write.csv(agent_df, file.path(OUTPUT_DIR, "interpret_agent_results.csv"),
                row.names = FALSE)
      cat(sprintf("\n[OK] interpret_agent CSV: %d clusters\n", nrow(agent_df)))
    }
  }
} else {
  cat("\n[INFO] Skipping interpret_agent (disabled or no enrichment)\n")
}

# ==============================================================================
# Generate Report
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

report_file <- file.path(OUTPUT_DIR, "REPORT.md")

tryCatch({
  sink(report_file)
  
  cat("# Myeloid Cell Subcluster Interpretation Report v3.0-MYELOID\n\n")
  cat("**Generated:** ", format(Sys.time()), "\n\n")
  cat("**Pipeline:** Myeloid cell specialized v3.0 (ssGSEA + LLM + interpret_agent)\n\n")
  cat("---\n\n")
  
  cat("## Dataset Summary\n\n")
  cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("- Genes: %d\n", nrow(seurat_obj)))
  cat(sprintf("- Subclusters analyzed: %d\n", length(unique(top_markers$cluster))))
  cat("\n")
  
  if (RUN_SSGSEA && !is.null(ssgsea_scores)) {
    cat(sprintf("- ssGSEA method: %s\n", SSGSEA_METHOD))
    cat(sprintf("- ssGSEA pathways: %d\n", nrow(ssgsea_scores)))
  }
  cat("\n---\n\n")
  
  # Annotation results
  if (file.exists(file.path(OUTPUT_DIR, "annotation_results.csv"))) {
    cat("## Cell Type Annotations\n\n")
    ann_df <- read.csv(file.path(OUTPUT_DIR, "annotation_results.csv"),
                      stringsAsFactors = FALSE)
    
    for (i in 1:nrow(ann_df)) {
      row <- ann_df[i, ]
      cat(sprintf("### %s\n\n", row$Cluster))
      cat(sprintf("**Cell Type:** %s  \n", row$Cell_Type))
      cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))
      
      if (!is.na(row$Regulatory_Drivers) && row$Regulatory_Drivers != "") {
        cat("**Regulatory Drivers:**  \n")
        cat(gsub("; ", "  \n- ", paste0("- ", row$Regulatory_Drivers)), "\n\n")
      }
      
      if (!is.na(row$Reasoning) && row$Reasoning != "") {
        cat("**Reasoning:**  \n", row$Reasoning, "\n\n")
      }
      
      cat("---\n\n")
    }
  }
  
  # Phenotype results
  if (file.exists(file.path(OUTPUT_DIR, "phenotype_results.csv"))) {
    cat("## Functional Phenotypes\n\n")
    phen_df <- read.csv(file.path(OUTPUT_DIR, "phenotype_results.csv"),
                       stringsAsFactors = FALSE)
    
    for (i in 1:nrow(phen_df)) {
      row <- phen_df[i, ]
      cat(sprintf("### %s\n\n", row$Cluster))
      cat(sprintf("**Phenotype:** %s  \n", row$Functional_Phenotype))
      cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))
      
      if (!is.na(row$Key_Processes) && row$Key_Processes != "") {
        cat("**Key Processes:**  \n")
        processes <- strsplit(row$Key_Processes, "; ")[[1]]
        cat(paste0("- ", head(processes, 5), collapse = "\n"), "\n\n")
      }
      
      cat("---\n\n")
    }
  }
  
  # interpret_agent results
  if (file.exists(file.path(OUTPUT_DIR, "interpret_agent_results.csv"))) {
    cat("## Deep Mode Interpretations\n\n")
    agent_df <- read.csv(file.path(OUTPUT_DIR, "interpret_agent_results.csv"),
                        stringsAsFactors = FALSE)
    
    for (i in 1:nrow(agent_df)) {
      row <- agent_df[i, ]
      cat(sprintf("### %s\n\n", row$Cluster))
      
      if (!is.na(row$Overview) && row$Overview != "") {
        cat("**Overview:**  \n", row$Overview, "\n\n")
      }
      
      if (!is.na(row$Regulatory_Drivers) && row$Regulatory_Drivers != "") {
        cat("**Regulatory Drivers:**  \n", row$Regulatory_Drivers, "\n\n")
      }
      
      if (!is.na(row$Narrative) && row$Narrative != "") {
        cat("**Narrative:**  \n", row$Narrative, "\n\n")
      }
      
      cat("---\n\n")
    }
  }
  
  sink()
  cat("[OK] Report saved: REPORT.md\n")
}, error = function(e) {
  tryCatch(sink(), error = function(e) NULL)
  cat("[WARN] Failed to generate report:", conditionMessage(e), "\n")
})

# ==============================================================================
# Complete
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("MYELOID ANALYSIS COMPLETE - v3.0-MYELOID\n")
cat("================================================================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat("  - REPORT.md ⭐⭐⭐                    Comprehensive summary\n")
cat("  - annotation_results.csv ⭐⭐        Cell Type annotations\n")
cat("  - phenotype_results.csv ⭐⭐         Functional phenotypes\n")
if (RUN_INTERPRET_AGENT) {
  cat("  - interpret_agent_results.csv ⭐⭐⭐  Deep mode interpretations\n")
}
if (RUN_SSGSEA) {
  cat("  - reports/ssgsea_scores_pseudobulk.rds  ssGSEA evidence\n")
  cat("  - reports/ssgsea_top_pathways.csv   Top pathways per cluster\n")
}
cat("\n")

cat("Key Features (v3.0-MYELOID):\n")
cat("  ✅ Myeloid cell expert context (M1/M2/monocyte states/DC maturation)\n")
cat("  ✅ ssGSEA primary evidence (Hallmark optimized for polarization)\n")
cat("  ✅ All v1.1 bug fixes integrated\n")
cat("  ✅ compareCluster enrichment (per-cluster)\n")
cat("  ✅ Dual-task LLM (annotation + phenotype)\n")
cat("  ✅ interpret_agent with regulatory networks\n")
cat("  ✅ Multisession parallel (fork-safe)\n")
cat("\n")

cat("Improvements over standard pipelines:\n")
cat("  • ssGSEA self-ranking evidence (no p-value dependency)\n")
cat("  • Z-score pathway specificity (auto-filter housekeeping)\n")
cat("  • PPI network expansion (enrichment + ssGSEA genes)\n")
cat("  • Enhanced LLM context with ssGSEA pathways\n")
cat("  • Memory-optimized pseudobulk computation\n")
cat("\n")

cat("================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")