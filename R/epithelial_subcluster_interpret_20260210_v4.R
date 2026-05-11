#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Cell Subcluster Interpretation Pipeline - v3.2.1-EPITHELIAL-PRODUCTION
# ==============================================================================
#
# Version: v3.2.1-EPITHELIAL-PRODUCTION (2026-02-10)
# Status: Production-ready with all P0/P1 fixes
# Base: v3.2 + Review fixes
#
# Key Features:
#   ✅ Epithelial expert context (Basal/Goblet/Ciliated differentiation)
#   ✅ Custom epithelial marker database (13 subtypes)
#   ✅ ALL ssGSEA methods (hallmark/GO BP/MF/CC/KEGG)
#   ✅ Dual-task LLM: Annotation + Phenotype
#   ✅ interpret_agent with DeepSeek-Reasoner
#   ✅ Load existing enrichment with TERM2GENE rebuild ⭐ FIXED
#   ✅ GSVA v1.40+ compatibility
#   ✅ API retry with exponential backoff ⭐ NEW
#   ✅ Enhanced regex parsing ⭐ NEW
#
# Fixes in v3.2.1:
#   • P0-1: LOAD_EXISTING_ENRICHMENT完整修复（TERM2GENE重建）
#   • P0-2: KEGG文件名统一
#   • P0-3: pseudobulk单细胞cluster类型统一
#   • P0-4: ssGSEA基因大小写对齐（内部转换）
#   • P0-5: compareCluster添加universe参数
#   • P0-6: fanyi依赖检查
#   • P1-3: API retry机制（指数退避）
#   • P1-4: parse_interpretation增强regex
#   • P2: genes_to_filter去重 + ignore.case
#   • P2: REPORT.md sink()保护
#
# ==============================================================================

# ==============================================================================
# Configuration Parameters
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0210/epithelial_interpret_v3_2_1"

# Reference Database Paths
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"

# Analysis Scope (Epithelial specific)
CLUSTER_COLUMN <- "subcluster"  # Main subcluster identities
CELLTYPE_COLUMN <- "celltypist_pred"  # Major cell type labels (optional)

# Analysis Methods
RUN_SSGSEA <- TRUE
RUN_ENRICHMENT <- TRUE
RUN_LLM_ANNOTATION <- TRUE
RUN_LLM_PHENOTYPE <- TRUE
RUN_INTERPRET_AGENT <- TRUE

# ⭐ FIXED: Load pre-computed enrichment (now rebuilds TERM2GENE)
LOAD_EXISTING_ENRICHMENT <- FALSE
ENRICHMENT_DIR <- NULL  # Path to directory with *.rds files

# ssGSEA Configuration - Run ALL methods
SSGSEA_METHOD <- c("hallmark", "go_bp", "go_mf", "go_cc", "kegg")
SSGSEA_CUSTOM_GMT <- NULL
SSGSEA_N_TOP <- 20

# LLM Configuration
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
LLM_MODEL <- "deepseek-reasoner"
LLM_MAX_RETRIES <- 3  # ⭐ NEW: API retry attempts
LLM_TIMEOUT_SECONDS <- 120  # API timeout
INTERPRET_AGENT_N_PATHWAYS <- 50
INTERPRET_AGENT_ADD_PPI <- TRUE  # ⚠️ Epithelial - PPI disabled

# Computational Parameters
N_CORES <- 1
TOP_N_MARKERS <- 50
MEMORY_LIMIT_GB <- 48

# ==============================================================================
# Thread Limiting
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

# ⭐ P0-6 FIX: Check fanyi dependency for interpret_agent
if (RUN_INTERPRET_AGENT) {
  if (!requireNamespace("fanyi", quietly = TRUE)) {
    stop(
      "ERROR: 'fanyi' package required for interpret_agent but not installed.\n",
      "Install with: install.packages('fanyi')\n",
      "Or set RUN_INTERPRET_AGENT=FALSE to skip deep mode analysis."
    )
  }
  library(fanyi)
  cat("[OK] fanyi package loaded for interpret_agent\n")
}

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

# ⭐ P1-3 FIX: API test with retry
test_api_connection <- function(api_key, model, max_retries = 3) {
  for (attempt in 1:max_retries) {
    result <- tryCatch(
      {
        payload <- list(
          model = model,
          messages = list(list(role = "user", content = "test")),
          max_tokens = 10
        )
        
        response <- httr::POST(
          url = "https://api.deepseek.com/v1/chat/completions",
          httr::add_headers(
            "Authorization" = paste("Bearer", api_key),
            "Content-Type" = "application/json"
          ),
          body = payload,
          encode = "json",
          httr::timeout(30)
        )
        
        status <- httr::status_code(response)
        
        if (status == 200) {
          return(list(success = TRUE, message = "API connection successful"))
        } else if (status %in% c(429, 500, 502, 503)) {
          # Retryable errors
          wait_time <- 2^attempt
          cat(sprintf("[WARN] API test failed with status %d, retry %d/%d after %ds\n",
                     status, attempt, max_retries, wait_time))
          Sys.sleep(wait_time)
          next
        } else {
          # Non-retryable error
          return(list(success = FALSE, 
                     message = sprintf("API test failed with status %d", status)))
        }
      },
      error = function(e) {
        if (attempt < max_retries) {
          wait_time <- 2^attempt
          cat(sprintf("[WARN] API test error: %s, retry %d/%d after %ds\n",
                     conditionMessage(e), attempt, max_retries, wait_time))
          Sys.sleep(wait_time)
          return(NULL)  # Continue to next attempt
        } else {
          return(list(success = FALSE, 
                     message = sprintf("API test error: %s", conditionMessage(e))))
        }
      }
    )
    
    if (!is.null(result)) {
      if (result$success) {
        cat(sprintf("[OK] %s\n\n", result$message))
        return(TRUE)
      } else {
        cat(sprintf("[ERROR] %s\n", result$message))
        return(FALSE)
      }
    }
  }
  
  cat("[ERROR] API test failed after all retries\n")
  return(FALSE)
}

# ==============================================================================
# Load Seurat Data
# ==============================================================================

cat("=== Loading Seurat Data ===\n")

seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

cat(sprintf("\nLoaded: %d cells x %d genes\n", ncol(seurat_obj), nrow(seurat_obj)))

# Validate required columns
need_cols <- c(CLUSTER_COLUMN)
missing <- setdiff(need_cols, colnames(seurat_obj@meta.data))
if (length(missing) > 0) {
  stop("ERROR: Missing required metadata column: ", paste(missing, collapse = ", "))
}

cat("[OK] Required metadata columns present\n")

# Check for cell type column (optional)
has_celltype_column <- CELLTYPE_COLUMN %in% colnames(seurat_obj@meta.data)
if (has_celltype_column) {
  cat(sprintf("[INFO] Cell type column '%s' found\n", CELLTYPE_COLUMN))
  celltype_dist <- table(seurat_obj@meta.data[[CELLTYPE_COLUMN]])
  cat("\nCell type distribution:\n")
  print(celltype_dist)
  cat("\n")
}

cat(sprintf("\nSubcluster distribution ('%s'):\n", CLUSTER_COLUMN))
cluster_counts <- table(seurat_obj@meta.data[[CLUSTER_COLUMN]], useNA = "ifany")
print(cluster_counts)
cat("\n")

# ==============================================================================
# Normalization
# ==============================================================================

cat("=== Normalization ===\n")

# ⭐ P1-2 IMPROVEMENT: Check if counts slot exists
has_counts <- tryCatch({
  counts_data <- GetAssayData(seurat_obj, slot = "counts")
  !is.null(counts_data) && sum(counts_data) > 0
}, error = function(e) FALSE)

if (has_counts) {
  cat("[INFO] Counts slot detected, performing normalization\n")
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
} else {
  cat("[WARN] No counts slot found, assuming data is pre-normalized\n")
  cat("[INFO] Skipping NormalizeData step\n")
}

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                  nfeatures = 4000, verbose = FALSE)

cat("[OK] Data preparation complete\n")

# ==============================================================================
# Compute Marker Genes (Parallel)
# ==============================================================================

cat("\n=== Computing Marker Genes (Parallel) ===\n")

# ⭐ P1-1 NOTE: Using future_lapply with multisession
# For very large Seurat objects (>10GB), consider:
#   - Using FindAllMarkers() instead
#   - Or plan("sequential") for marker computation
#   - Or plan("multicore") on Linux if fork is allowed
cat("[INFO] Using multisession parallelization (N_CORES=%d)\n", N_CORES)
cat("[INFO] This may use significant memory for large objects\n\n")

Idents(seurat_obj) <- CLUSTER_COLUMN
clusters <- levels(Idents(seurat_obj))

cat(sprintf("Finding markers for %d subclusters...\n", length(clusters)))

marker_list <- future_lapply(
  clusters,
  function(cluster_id) {
    tryCatch(
      {
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
        cat(sprintf("[WARN] Failed to find markers for %s: %s\n",
                   cluster_id, conditionMessage(e)))
        NULL
      }
    )
  },
  future.seed = TRUE
)

names(marker_list) <- clusters
marker_list <- marker_list[!sapply(marker_list, is.null)]

# Report skipped clusters
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
# Prepare Top Markers for Enrichment (Epithelial specific filtering)
# ==============================================================================

cat("\n=== Preparing Top Markers ===\n")

# ⭐ P2-1 FIX: genes_to_filter with unique() and ignore.case
genes_to_filter <- unique(c(
  grep("^MT-", rownames(seurat_obj), value = TRUE, ignore.case = TRUE),
  grep("^RP[SL]", rownames(seurat_obj), value = TRUE, ignore.case = TRUE),
  # Immediate early genes
  "FOS", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
  "ZFP36", "DUSP1", "DUSP2", "IER2", "IER3", "ATF3", "BTG2",
  "FOSB", "NR4A1", "NR4A2", "NR4A3",
  # Heat shock proteins
  "HSP90AA1", "HSPA1A", "HSPA1B", "DNAJB1", "HSPA6",
  # Cell cycle (G2/M and S phase)
  "MKI67", "TOP2A", "UBE2C", "CENPF", "PCNA",
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7"
))

cat(sprintf("Filtering %d potentially confounding genes (deduplicated)\n", 
           length(genes_to_filter)))

# Auto-detect logFC column name
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
  mutate(gene = toupper(gene)) %>%
  dplyr::select(gene, cluster)

cat(sprintf("Selected top %d clean markers per cluster\n", TOP_N_MARKERS))
cat(sprintf("Total markers for enrichment: %d\n", nrow(top_markers)))

write.csv(top_markers, file.path(OUTPUT_DIR, "top_markers_filtered.csv"),
         row.names = FALSE)

# Detailed filtering info
filtered_info <- data.frame(
  category = c(
    "MT genes", "Ribosomal genes", "Stress/IEG genes",
    "Cell cycle genes", "HSP genes", "Total filtered"
  ),
  count = c(
    sum(grepl("^MT-", genes_to_filter, ignore.case = TRUE)),
    sum(grepl("^RP[SL]", genes_to_filter, ignore.case = TRUE)),
    sum(genes_to_filter %in% c(
      "FOS", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
      "ZFP36", "DUSP1", "DUSP2", "IER2", "IER3", "ATF3", "BTG2",
      "FOSB", "NR4A1", "NR4A2", "NR4A3"
    )),
    sum(genes_to_filter %in% c(
      "MKI67", "TOP2A", "UBE2C", "CENPF", "PCNA",
      "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7"
    )),
    sum(genes_to_filter %in% c(
      "HSP90AA1", "HSPA1A", "HSPA1B", "DNAJB1", "HSPA6"
    )),
    length(genes_to_filter)
  )
)

write.csv(filtered_info, file.path(OUTPUT_DIR, "filtered_genes_info.csv"),
         row.names = FALSE)

cat("[OK] Saved filtered_genes_info.csv\n")

# ==============================================================================
# Custom Epithelial Marker Database
# ==============================================================================

cat("\n=== Creating Custom Epithelial Marker Database ===\n")

# ⭐ P2-2 FIX: Removed CC10 (non-standard SYMBOL)
epithelial_markers_db <- data.frame(
  subtype = c(
    "Basal_Proliferative", "Basal_Quiescent", "Basal_Squamous",
    "Goblet", "Club", "Serous",
    "Ciliated_Mature", "Deuterosomal",
    "Ionocyte", "Tuft", "Neuroendocrine",
    "Suprabasal", "Intermediate"
  ),
  markers = c(
    "TP63,KRT5,KRT14,ITGA6,ITGB4,COL17A1,MKI67,TOP2A",
    "TP63,KRT5,KRT14,ITGA6,ITGB4,COL17A1",
    "TP63,KRT5,KRT14,KRT6A,KRT16,S100A2,SPRR3",
    "MUC5AC,MUC5B,TFF3,SPDEF,AGR2,WFDC2",
    "SCGB1A1,SCGB3A1,CYP2F1,BPIFA1",  # Removed CC10
    "LTF,LCN2,MUC7,PRB3,PRB4",
    "FOXJ1,DNAI1,DNAI2,RSPH1,TUBA1A,TUBB4B,PIFO",
    "DEUP1,CCNO,CDC20B,FOXN4,PLK4",
    "FOXI1,CFTR,ATP6V1B1,ATP6V0D2,ASCL3",
    "POU2F3,TRPM5,GFI1B,AVIL,IL25",
    "CHGA,CHGB,SYP,INSM1,ASCL1,CALCA",
    "KRT13,KRT4,CLDN4,DSP,PKP1",
    "KRT8,KRT18,KRT19,MUC1,EPCAM"
  ),
  stringsAsFactors = FALSE
)

epithelial_markers_term2gene <- epithelial_markers_db %>%
  tidyr::separate_rows(markers, sep = ",") %>%
  dplyr::mutate(markers = trimws(toupper(markers))) %>%
  dplyr::filter(markers != "") %>%
  dplyr::select(term = subtype, gene = markers) %>%
  distinct()

cat(sprintf(
  "[OK] Created custom epithelial markers: %d subtypes, %d markers\n",
  length(unique(epithelial_markers_term2gene$term)),
  nrow(epithelial_markers_term2gene)
))

# ⭐ P0-5 FIX: Define universe for enrichment
universe_genes <- setdiff(
  toupper(rownames(seurat_obj)),
  toupper(genes_to_filter)
)

cat(sprintf(
  "\nDefined universe: %d genes (total: %d, filtered: %d)\n",
  length(universe_genes),
  nrow(seurat_obj),
  length(genes_to_filter)
))

# ==============================================================================
# ⭐ P0-1 FIX: Build TERM2GENE (Always, even in LOAD mode)
# ==============================================================================

cat("\n=== Building TERM2GENE Databases (Always) ===\n")

# This section now runs regardless of LOAD_EXISTING_ENRICHMENT
# Only compareCluster() computation is skipped in LOAD mode

# --- CellMarker ---
cellmarker_term2gene <- NULL
tryCatch({
  db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)
  db <- db %>% filter(grepl("Human", species, ignore.case = TRUE))
  
  cat(sprintf("[OK] Loaded %d CellMarker entries\n", nrow(db)))
  
  # Epithelial-related filtering
  epithelial_related <- db %>%
    filter(
      grepl("Epithelial|Basal|Goblet|Ciliated|Club|Secretory|Airway|Respiratory|Lung",
            cell_name, ignore.case = TRUE) |
      grepl("Lung|Airway|Respiratory|Nasal|Sinus|Trachea|Bronchus|Bronchi",
            tissue_type, ignore.case = TRUE)
    )
  
  cat(sprintf("[OK] Found %d epithelial-related entries\n", nrow(epithelial_related)))
  
  if (nrow(epithelial_related) >= 30) {
    db <- epithelial_related
    cat("[INFO] Using epithelial-specific subset for enrichment\n")
  }
  
  # Parse markers
  term2gene_list <- list()
  
  for (i in 1:nrow(db)) {
    row <- db[i, ]
    cell_type <- row$cell_name
    markers_raw <- row$marker
    
    if (!is.na(markers_raw) && markers_raw != "") {
      markers <- unlist(strsplit(markers_raw, "[,;\\s]+"))
      markers <- gsub('["\r\n\\[\\]]', '', markers)
      markers <- trimws(markers)
      markers <- toupper(markers)
      markers <- unique(markers[markers != "" & !is.na(markers)])
      
      if (length(markers) > 0) {
        term2gene_list[[length(term2gene_list) + 1]] <- data.frame(
          term = cell_type,
          gene = markers,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  cellmarker_term2gene <- bind_rows(term2gene_list)
  cat(sprintf("[OK] Prepared CellMarker TERM2GENE: %d pairs\n", nrow(cellmarker_term2gene)))
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
  cat(sprintf("[OK] Loaded %d PanglaoDB human markers\n", nrow(db)))
  
  panglaodb_term2gene <- db %>%
    dplyr::select(cell_type, gene_symbol) %>%
    mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct() %>%
    dplyr::rename(term = cell_type, gene = gene_symbol)
  
  cat(sprintf("[OK] Prepared PanglaoDB TERM2GENE: %d pairs\n", nrow(panglaodb_term2gene)))
}, error = function(e) {
  cat("[WARN] Failed to load PanglaoDB:", conditionMessage(e), "\n")
})

# --- GO Terms (GMT-based) ---
go_bp_gmt <- NULL
go_mf_gmt <- NULL
go_cc_gmt <- NULL

tryCatch({
  if (!file.exists(GMT_GO_ALL)) {
    stop("GMT file not found at: ", GMT_GO_ALL)
  }
  
  cat("Loading c5.all GMT file...\n")
  go_all_gmt <- read.gmt(GMT_GO_ALL)
  
  go_bp_gmt <- go_all_gmt[grep("^GOBP_", go_all_gmt$term), ]
  go_mf_gmt <- go_all_gmt[grep("^GOMF_", go_all_gmt$term), ]
  go_cc_gmt <- go_all_gmt[grep("^GOCC_", go_all_gmt$term), ]
  
  cat(sprintf("✓ GO BP: %d gene sets loaded\n", length(unique(go_bp_gmt$term))))
  cat(sprintf("✓ GO MF: %d gene sets loaded\n", length(unique(go_mf_gmt$term))))
  cat(sprintf("✓ GO CC: %d gene sets loaded\n", length(unique(go_cc_gmt$term))))
  
  # Gene coverage check
  all_genes_in_gmt <- unique(go_all_gmt$gene)
  markers_in_gmt <- sum(unique(top_markers$gene) %in% all_genes_in_gmt)
  total_markers <- length(unique(top_markers$gene))
  cat(sprintf("✓ Gene coverage: %d/%d markers (%.1f%%) found in GMT\n",
             markers_in_gmt, total_markers, markers_in_gmt / total_markers * 100))
}, error = function(e) {
  cat("[ERROR] Failed to load GO GMT:", conditionMessage(e), "\n")
})

# --- MSigDB (Hallmark + KEGG) ---
hallmark_term2gene <- NULL
kegg_term2gene <- NULL

tryCatch({
  read_gmt <- function(gmt_file) {
    lines <- readLines(gmt_file)
    gene_sets_list <- lapply(lines, function(line) {
      parts <- strsplit(line, "\t")[[1]]
      list(name = parts[1], description = parts[2], genes = parts[-(1:2)])
    })
    term2gene_list <- lapply(gene_sets_list, function(gs) {
      if (length(gs$genes) > 0) {
        data.frame(term = rep(gs$name, length(gs$genes)),
                  gene = gs$genes, stringsAsFactors = FALSE)
      }
    })
    do.call(rbind, term2gene_list)
  }
  
  all_genesets <- read_gmt(MSIGDB_GMT_PATH)
  cat(sprintf("[OK] Loaded %d MSigDB gene set entries\n", nrow(all_genesets)))
  
  hallmark_term2gene <- all_genesets %>%
    filter(grepl("^HALLMARK_", term)) %>%
    dplyr::select(term, gene)
  
  kegg_term2gene <- all_genesets %>%
    filter(grepl("KEGG_", term)) %>%
    dplyr::select(term, gene)
  
  if (nrow(hallmark_term2gene) > 0) {
    cat(sprintf("✓ Hallmark: %d pathways\n", length(unique(hallmark_term2gene$term))))
  }
  if (nrow(kegg_term2gene) > 0) {
    cat(sprintf("✓ KEGG: %d pathways\n", length(unique(kegg_term2gene$term))))
  }
}, error = function(e) {
  cat("[ERROR] Failed to load MSigDB GMT:", conditionMessage(e), "\n")
})

cat("\n[OK] All TERM2GENE databases built\n")

# ==============================================================================
# Load or Run Enrichment Analysis
# ==============================================================================

# Initialize enrichment objects
epithelial_markers_enrich <- NULL
cellmarker_enrich <- NULL
panglaodb_enrich <- NULL
go_bp_enrich <- NULL
go_mf_enrich <- NULL
go_cc_enrich <- NULL
hallmark_enrich <- NULL
kegg_enrich <- NULL  # ⭐ P0-2 FIX: Renamed from msigdb_kegg_enrich

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
    # ⭐ P0-2 FIX: Updated file list with consistent naming
    rds_files <- c(
      "epithelial_markers_enrich.rds", "cellmarker_enrich.rds",
      "panglaodb_enrich.rds", "go_bp_enrich.rds", "go_mf_enrich.rds",
      "go_cc_enrich.rds", "hallmark_enrich.rds", "kegg_enrich.rds"
    )
    
    for (rds_file in rds_files) {
      rds_path <- file.path(enrich_dir, rds_file)
      if (file.exists(rds_path)) {
        obj_name <- gsub(".rds$", "", rds_file)
        assign(obj_name, readRDS(rds_path))
        cat(sprintf("  [OK] Loaded %s\n", rds_file))
      }
    }
    
    cat("\n[INFO] Loaded enrichment results, skipping compareCluster computation\n")
    cat("[INFO] TERM2GENE databases were rebuilt for interpret_agent\n")
  }
}

if (!LOAD_EXISTING_ENRICHMENT && RUN_ENRICHMENT) {
  cat("\n=== Running Enrichment Analysis (compareCluster) ===\n")
  
  # Helper function
  run_compareCluster_enrichment <- function(term2gene_df, db_name, 
                                           output_filename = NULL) {
    if (is.null(term2gene_df) || nrow(term2gene_df) == 0) {
      return(NULL)
    }
    
    cat(sprintf("\n--- Database: %s ---\n", db_name))
    
    # Smaller minGSSize for cell type marker databases
    min_size <- if (db_name %in% c("Epithelial Markers", "CellMarker", "PanglaoDB")) 5 else 10
    max_size <- 500
    
    result <- tryCatch({
      # ⭐ P0-5 FIX: Added universe parameter
      compareCluster(
        gene ~ cluster,
        data = top_markers,
        fun = enricher,
        TERM2GENE = term2gene_df,
        universe = universe_genes,  # ⭐ NEW
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
      
      # Determine output filename
      if (is.null(output_filename)) {
        output_filename <- paste0(tolower(gsub(" ", "_", db_name)), "_enrich.rds")
      }
      
      saveRDS(result, file.path(OUTPUT_DIR, "reports", output_filename))
    }
    
    return(result)
  }
  
  # Run enrichment for each database
  epithelial_markers_enrich <- run_compareCluster_enrichment(
    epithelial_markers_term2gene, "Epithelial Markers"
  )
  cellmarker_enrich <- run_compareCluster_enrichment(cellmarker_term2gene, "CellMarker")
  panglaodb_enrich <- run_compareCluster_enrichment(panglaodb_term2gene, "PanglaoDB")
  go_bp_enrich <- run_compareCluster_enrichment(go_bp_gmt, "GO BP")
  go_mf_enrich <- run_compareCluster_enrichment(go_mf_gmt, "GO MF")
  go_cc_enrich <- run_compareCluster_enrichment(go_cc_gmt, "GO CC")
  hallmark_enrich <- run_compareCluster_enrichment(hallmark_term2gene, "Hallmark")
  
  # ⭐ P0-2 FIX: KEGG with explicit filename
  kegg_enrich <- run_compareCluster_enrichment(
    kegg_term2gene, "KEGG", output_filename = "kegg_enrich.rds"
  )
  
  cat("\n[OK] Enrichment analysis complete\n")
}

# Collect all enrichment results
enrichment_results <- list(
  epithelial_markers = epithelial_markers_enrich,
  cellmarker = cellmarker_enrich,
  panglaodb = panglaodb_enrich,
  go_bp = go_bp_enrich,
  go_mf = go_mf_enrich,
  go_cc = go_cc_enrich,
  hallmark = hallmark_enrich,
  kegg = kegg_enrich
)

# Remove NULL entries
enrichment_results <- enrichment_results[!sapply(enrichment_results, is.null)]

# ==============================================================================
# ssGSEA Analysis (ALL Methods) - FIXED for GSVA v1.40+
# ==============================================================================

ssgsea_results_all <- list()

if (RUN_SSGSEA) {
  cat("\n=== ssGSEA Analysis (ALL Methods) ===\n")
  
  # ⭐ P0-4 FIX: Helper function with internal gene symbol alignment
  run_ssgsea_method <- function(method_name, term2gene, sample_avg) {
    cat(sprintf("\n--- Running ssGSEA: %s ---\n", method_name))
    
    # Filter gene sets by size (10-500 genes)
    geneset_list <- split(term2gene$gene, term2gene$term)
    geneset_list <- geneset_list[sapply(geneset_list, length) >= 10 &
                                  sapply(geneset_list, length) <= 500]
    
    cat(sprintf("  Gene sets: %d (after size filtering)\n", length(geneset_list)))
    
    if (length(geneset_list) == 0) {
      cat("  [SKIP] No valid gene sets\n")
      return(NULL)
    }
    
    # ⭐ P0-4 FIX: Internal gene symbol case alignment (Method B)
    # Create a working copy with uppercased rownames
    sample_avg_upper <- sample_avg
    original_rownames <- rownames(sample_avg)
    rownames(sample_avg_upper) <- toupper(original_rownames)
    
    # Handle potential duplicates from uppercasing
    if (any(duplicated(rownames(sample_avg_upper)))) {
      cat("  [INFO] Detected duplicate gene symbols after uppercasing, using make.unique\n")
      rownames(sample_avg_upper) <- make.unique(rownames(sample_avg_upper))
    }
    
    # Ensure gene sets are also uppercase
    geneset_list <- lapply(geneset_list, toupper)
    
    # Find common genes
    common_genes <- intersect(rownames(sample_avg_upper), unique(unlist(geneset_list)))
    cat(sprintf("  Common genes: %d\n", length(common_genes)))
    
    if (length(common_genes) < 100) {
      cat("  [SKIP] Too few common genes (<100)\n")
      return(NULL)
    }
    
    # Subset matrix
    mat_subset <- sample_avg_upper[common_genes, , drop = FALSE]
    
    # ⭐ P0-3 FIX: GSVA v1.40+ compatibility - use param instead of BPPARAM
    bpparam <- BiocParallel::SnowParam(workers = N_CORES, progressbar = FALSE)
    
    # Run GSVA with ssgsea method
    scores <- tryCatch({
      gsva(
        expr = as.matrix(mat_subset),
        gset.idx.list = geneset_list,
        method = "ssgsea",
        kcdf = "Gaussian",
        abs.ranking = FALSE,
        min.sz = 10,
        max.sz = 500,
        verbose = FALSE,
        param = bpparam  # ⭐ FIXED: Use param instead of BPPARAM
      )
    }, error = function(e) {
      cat("  [ERROR]", conditionMessage(e), "\n")
      NULL
    })
    
    if (is.null(scores)) return(NULL)
    
    # Z-score normalization (for relative comparison across clusters)
    scores_z <- t(scale(t(scores)))
    
    cat(sprintf("  [OK] Computed %d pathway scores\n", nrow(scores_z)))
    
    # Save scores
    saveRDS(scores_z, file.path(OUTPUT_DIR, "reports",
                               paste0("ssgsea_scores_", method_name, ".rds")))
    
    # Extract top N pathways per cluster
    # ⭐ P2-3 NOTE: Sorted by z-score (relative activation across clusters)
    top_list <- lapply(colnames(scores_z), function(cid) {
      vals <- scores_z[, cid]
      top_idx <- order(vals, decreasing = TRUE)[1:min(SSGSEA_N_TOP, length(vals))]
      data.frame(
        cluster = cid,
        pathway = rownames(scores_z)[top_idx],
        zscore = vals[top_idx],
        method = method_name,
        stringsAsFactors = FALSE
      )
    })
    
    top_df <- bind_rows(top_list)
    write.csv(top_df,
             file.path(OUTPUT_DIR, "reports",
                      paste0("ssgsea_top_pathways_", method_name, ".csv")),
             row.names = FALSE)
    
    cat(sprintf("  [OK] Saved top %d pathways per cluster\n", SSGSEA_N_TOP))
    
    return(scores_z)
  }
  
  # Compute Pseudobulk Profiles
  cat("\nComputing pseudobulk profiles...\n")
  
  meta_sub <- seurat_obj@meta.data[, c(CLUSTER_COLUMN), drop = FALSE]
  colnames(meta_sub) <- "cluster"
  meta_sub$cluster <- as.character(meta_sub$cluster)
  
  expr_mat <- GetAssayData(seurat_obj, slot = "data")
  
  # ⭐ P0-3 FIX: Memory-optimized pseudobulk with consistent return type
  sample_avg <- sapply(unique(meta_sub$cluster), function(cid) {
    cells_in_cluster <- rownames(meta_sub)[meta_sub$cluster == cid]
    # Always use rowMeans, even for single-cell clusters (returns numeric vector)
    Matrix::rowMeans(expr_mat[, cells_in_cluster, drop = FALSE])
  })
  
  cat(sprintf("[OK] Pseudobulk: %d genes x %d clusters\n",
             nrow(sample_avg), ncol(sample_avg)))
  
  # Load Gene Sets and Run ssGSEA for Each Method
  
  # Hallmark
  if ("hallmark" %in% SSGSEA_METHOD) {
    cat("\n[INFO] Loading Hallmark gene sets...\n")
    if (!is.null(hallmark_term2gene) && nrow(hallmark_term2gene) > 0) {
      ssgsea_results_all[["hallmark"]] <- run_ssgsea_method(
        "hallmark", hallmark_term2gene, sample_avg
      )
    }
  }
  
  # GO BP/MF/CC
  if (any(c("go_bp", "go_mf", "go_cc") %in% SSGSEA_METHOD)) {
    cat("\n[INFO] Loading GO gene sets...\n")
    
    if ("go_bp" %in% SSGSEA_METHOD && !is.null(go_bp_gmt) && nrow(go_bp_gmt) > 0) {
      ssgsea_results_all[["go_bp"]] <- run_ssgsea_method("go_bp", go_bp_gmt, sample_avg)
    }
    
    if ("go_mf" %in% SSGSEA_METHOD && !is.null(go_mf_gmt) && nrow(go_mf_gmt) > 0) {
      ssgsea_results_all[["go_mf"]] <- run_ssgsea_method("go_mf", go_mf_gmt, sample_avg)
    }
    
    if ("go_cc" %in% SSGSEA_METHOD && !is.null(go_cc_gmt) && nrow(go_cc_gmt) > 0) {
      ssgsea_results_all[["go_cc"]] <- run_ssgsea_method("go_cc", go_cc_gmt, sample_avg)
    }
  }
  
  # KEGG
  if ("kegg" %in% SSGSEA_METHOD) {
    cat("\n[INFO] Loading KEGG gene sets...\n")
    if (!is.null(kegg_term2gene) && nrow(kegg_term2gene) > 0) {
      ssgsea_results_all[["kegg"]] <- run_ssgsea_method("kegg", kegg_term2gene, sample_avg)
    }
  }
  
  # Custom GMT (if provided)
  if ("custom" %in% SSGSEA_METHOD && !is.null(SSGSEA_CUSTOM_GMT)) {
    if (file.exists(SSGSEA_CUSTOM_GMT)) {
      cat("\n[INFO] Loading custom GMT...\n")
      custom_gmt <- read.gmt(SSGSEA_CUSTOM_GMT)
      if (nrow(custom_gmt) > 0) {
        ssgsea_results_all[["custom"]] <- run_ssgsea_method("custom", custom_gmt, sample_avg)
      }
    } else {
      cat("[WARN] Custom GMT file not found:", SSGSEA_CUSTOM_GMT, "\n")
    }
  }
  
  cat("\n[OK] All ssGSEA methods complete\n")
  cat(sprintf("  Computed %d method(s)\n", length(ssgsea_results_all)))
  
  # Use primary method for downstream
  primary_method <- if ("hallmark" %in% names(ssgsea_results_all)) {
    "hallmark"
  } else if (length(ssgsea_results_all) > 0) {
    names(ssgsea_results_all)[1]
  } else {
    NULL
  }
  
  if (!is.null(primary_method)) {
    cat(sprintf("\n[INFO] Using '%s' as primary method for LLM context\n", primary_method))
    ssgsea_scores <- ssgsea_results_all[[primary_method]]
  } else {
    ssgsea_scores <- NULL
  }
  
} else {
  cat("\n[INFO] Skipping ssGSEA (RUN_SSGSEA=FALSE)\n")
  ssgsea_scores <- NULL
}

# ==============================================================================
# Epithelial Cell Expert Context Generator
# ==============================================================================

generate_epithelial_context <- function() {
  paste(
    "CRITICAL CONTEXT FOR EPITHELIAL CELL ANNOTATION:",
    "",
    "These are respiratory epithelial cells from nasal cavity, sinus, bronchi, and lung.",
    "Data source: CRSwNP (Chronic Rhinosinusitis with Nasal Polyps) study.",
    "",
    "=== LINEAGE GATE VERIFICATION (Important) ===",
    "Before finalizing annotation, verify lineage markers:",
    "• Epithelial: EPCAM+, CDH1+ (if absent → potential contamination)",
    "• Basal cells: TP63+, KRT5+, KRT14+ (stem/progenitor)",
    "• Secretory: MUC5AC/MUC5B/SCGB1A1",
    "• Ciliated: FOXJ1+, DNAI1+, DNAI2+",
    "• Immune: CD45/PTPRC (if dominant → immune contamination)",
    "",
    "⚠️ State vs Identity Warning:",
    "If markers are mainly OXPHOS (ATP5*, COX*, ND*) or stress (HSP*, FOS, JUN),",
    "note this is likely a METABOLIC/ACTIVATION STATE, not cell identity.",
    "",
    "Key epithelial biology principles:",
    "",
    "1. EPITHELIAL CELL TYPES:",
    "   - Basal cells: TP63+KRT5+KRT14+ stem/progenitor cells",
    "     * Quiescent basal: Reserve stem cells in homeostasis",
    "     * Proliferative basal: MKI67+ actively dividing progenitors",
    "     * Squamous basal: KRT6A+KRT16+ squamous metaplasia precursors",
    "",
    "   - Secretory cells: Mucus and antimicrobial production",
    "     * Goblet cells: MUC5AC+MUC5B+ mucin-producing",
    "     * Club cells: SCGB1A1+ secretory, Clara cell-like",
    "     * Serous cells: LTF+LCN2+ antimicrobial protein secretion",
    "",
    "   - Ciliated cells: FOXJ1+ motile cilia for mucociliary clearance",
    "     * Mature ciliated: DNAI1+DNAI2+ fully differentiated",
    "     * Deuterosomal cells: DEUP1+CCNO+ intermediate differentiating state",
    "",
    "   - Specialized cells: Rare chemosensory and regulatory functions",
    "     * Ionocytes: FOXI1+CFTR+ ion transport regulation",
    "     * Tuft cells: POU2F3+TRPM5+ chemosensory, immune regulation",
    "     * Neuroendocrine: CHGA+CHGB+ hormonal signaling",
    "",
    "   - Transitional states:",
    "     * Suprabasal: KRT13+KRT4+ intermediate between basal and differentiated",
    "     * Intermediate: KRT8+KRT18+KRT19+ multi-lineage potential",
    "",
    "2. DIFFERENTIATION TRAJECTORIES:",
    "   - Basal → Secretory lineage:",
    "     Basal → Suprabasal → Secretory (Goblet/Club/Serous)",
    "     Markers: TP63 loss, secretory gene upregulation (MUC5AC, SCGB1A1)",
    "",
    "   - Basal → Ciliated lineage:",
    "     Basal → Deuterosomal → Mature Ciliated",
    "     Markers: TP63 loss, FOXJ1 gain, centriole/basal body assembly",
    "",
    "   - Basal → Specialized lineage:",
    "     Basal → Specialized (Ionocyte/Tuft/NE)",
    "     Markers: Cell-type specific TF programs (FOXI1, POU2F3, ASCL1)",
    "",
    "3. FUNCTIONAL STATES:",
    "   - Proliferation: Quiescent (G0) vs cycling (G1/S/G2/M)",
    "     MKI67+, TOP2A+, PCNA+ indicate active proliferation",
    "",
    "   - Differentiation stage: Undifferentiated basal vs intermediate vs terminal",
    "     Loss of TP63, gain of lineage-specific TFs",
    "",
    "   - Secretory activity: Mucus production level, antimicrobial output",
    "     High MUC5AC/MUC5B in active goblet cells",
    "",
    "   - Ciliary function: Ciliary assembly, motility, mucociliary clearance",
    "     FOXJ1+ with high DNAI1/DNAI2 for functional cilia",
    "",
    "   - Barrier integrity: Tight junction formation, epithelial polarization",
    "     CLDN, OCLN, ZO-1 markers",
    "",
    "   - Stress response: Oxidative stress, inflammation sensing",
    "     HSP genes, inflammatory cytokines (if present)",
    "",
    "4. TISSUE CONTEXT (Respiratory tissues):",
    "   - Basal cells maintain tissue homeostasis",
    "   - Goblet cell hyperplasia in CRSwNP",
    "   - Ciliated cell dysfunction in disease",
    "   - Mucus hypersecretion pathology",
    "",
    "ANNOTATION INSTRUCTIONS:",
    "For EACH cluster, specify:",
    "  a) Lineage: Basal / Secretory / Ciliated / Specialized / Transitional",
    "  b) Differentiation stage: Quiescent / Proliferative / Intermediate / Terminal",
    "  c) Functional state: Resting / Activated / Stressed / Proliferating",
    "  d) Cell subtype: Specific marker-based identity",
    "",
    "Example good annotations:",
    "  - 'Basal cells (quiescent, TP63+KRT5+, reserve stem)'",
    "  - 'Goblet cells (MUC5AC+, mucus hypersecretion, activated)'",
    "  - 'Mature ciliated cells (FOXJ1+DNAI1+, functional cilia)'",
    "  - 'Deuterosomal cells (DEUP1+CCNO+, ciliated intermediate)'",
    "  - 'Proliferative basal cells (TP63+MKI67+, actively dividing)'",
    "",
    "Always base reasoning on:",
    "  1. Specific marker gene expression",
    "  2. Enriched pathways (GO/Hallmark/KEGG)",
    "  3. ssGSEA evidence (if available)",
    "  4. Biological plausibility in respiratory epithelium",
    sep = "\n"
  )
}

epithelial_context <- generate_epithelial_context()

# ==============================================================================
# LLM Interpretation Functions (Annotation + Phenotype)
# ==============================================================================

# ⭐ P1-3 FIX: API call with retry and exponential backoff
call_deepseek_api <- function(prompt, api_key, model, max_retries = 3, 
                              timeout_seconds = 120) {
  for (attempt in 1:max_retries) {
    result <- tryCatch(
      {
        payload <- list(
          model = model,
          messages = list(list(role = "user", content = prompt)),
          temperature = 0.3,
          max_tokens = 2000  # Limit output length
        )
        
        response <- httr::POST(
          url = "https://api.deepseek.com/v1/chat/completions",
          httr::add_headers(
            "Authorization" = paste("Bearer", api_key),
            "Content-Type" = "application/json"
          ),
          body = payload,
          encode = "json",
          httr::timeout(timeout_seconds)
        )
        
        status <- httr::status_code(response)
        
        if (status == 200) {
          content <- httr::content(response, as = "parsed")
          interpretation <- content$choices[[1]]$message$content
          return(list(success = TRUE, content = interpretation))
        } else if (status %in% c(429, 500, 502, 503)) {
          # Retryable errors - rate limit or server errors
          wait_time <- 2^attempt  # Exponential backoff: 2s, 4s, 8s
          cat(sprintf("  [WARN] API returned status %d, retry %d/%d after %ds\n",
                     status, attempt, max_retries, wait_time))
          Sys.sleep(wait_time)
          next
        } else {
          # Non-retryable error
          return(list(success = FALSE, 
                     message = sprintf("API returned status %d", status)))
        }
      },
      error = function(e) {
        if (attempt < max_retries) {
          wait_time <- 2^attempt
          cat(sprintf("  [WARN] API error: %s, retry %d/%d after %ds\n",
                     conditionMessage(e), attempt, max_retries, wait_time))
          Sys.sleep(wait_time)
          return(NULL)  # Continue to next attempt
        } else {
          return(list(success = FALSE, 
                     message = sprintf("API error: %s", conditionMessage(e))))
        }
      }
    )
    
    if (!is.null(result)) {
      return(result)
    }
  }
  
  return(list(success = FALSE, message = "Max retries exceeded"))
}

interpret_epithelial <- function(cluster_id, marker_genes, enrichment_data = NULL,
                                 ssgsea_data = NULL, task = "annotation",
                                 api_key = DEEPSEEK_API_KEY) {
  cat(sprintf("\n--- Epithelial Interpretation: %s (Task: %s) ---\n",
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
    top_pathways <- ssgsea_data %>%
      arrange(desc(zscore)) %>%
      head(SSGSEA_N_TOP)
    
    pathway_lines <- sprintf("  %s (Z=%.2f)",
                            top_pathways$pathway,
                            top_pathways$zscore)
    
    ssgsea_text <- paste("\n\nssGSEA Top Pathways (by relative activation):",
                        paste(pathway_lines, collapse = "\n"),
                        sep = "\n")
  }
  
  # Construct prompt
  if (task == "annotation") {
    prompt <- sprintf(
      "%s

Based on marker genes, enrichment, and ssGSEA evidence, provide:

Cell Type: [specific epithelial subtype]
Confidence: [High/Medium/Low]
Regulatory Drivers: [driver1; driver2; driver3]
Markers: [marker1; marker2; marker3]
Reasoning: [brief explanation]

Marker Genes: %s%s%s

CRITICAL: Include ALL FIVE fields with exact format.",
      epithelial_context, markers_text, enrichment_text, ssgsea_text
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
      epithelial_context, markers_text, enrichment_text, ssgsea_text
    )
  } else {
    stop("Invalid task type. Must be 'annotation' or 'phenotype'")
  }
  
  # Call API with retry
  result <- call_deepseek_api(
    prompt = prompt,
    api_key = api_key,
    model = LLM_MODEL,
    max_retries = LLM_MAX_RETRIES,
    timeout_seconds = LLM_TIMEOUT_SECONDS
  )
  
  if (result$success) {
    interpretation <- result$content
    cat(sprintf("[OK] Received interpretation (%d chars)\n", nchar(interpretation)))
    
    # Save raw prompt and response for debugging
    raw_dir <- file.path(OUTPUT_DIR, "reports", "llm_raw")
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
    
    return(interpretation)
  } else {
    cat(sprintf("[ERROR] API call failed: %s\n", result$message))
    return(NULL)
  }
}

# ⭐ P1-4 FIX: Enhanced regex parsing (Method B)
parse_interpretation <- function(text, task = "annotation") {
  if (is.null(text) || nchar(text) == 0 || is.na(text)) {
    return(list(raw = NA_character_))
  }
  
  # Enhanced regex with more flexible field name matching
  extract_field <- function(text, field_name) {
    # Create flexible pattern that allows:
    # - Multiple spaces/tabs
    # - Case variations (Cell Type / Cell type / cell type)
    # - With or without colon spacing
    field_pattern <- gsub(" ", "\\\\s+", field_name)  # Allow multiple spaces
    pattern <- sprintf(
      "(?si)%s\\s*:\\s*(.+?)(?=\\n\\s*[A-Z][A-Za-z][A-Za-z ]*\\s*:|$)",
      field_pattern
    )
    
    m <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
    if (length(m) >= 2) {
      trimws(m[2])
    } else {
      # Fallback: try without case sensitivity and looser spacing
      pattern_loose <- sprintf("(?i)%s[:\\s]+(.+?)(?=\\n|$)", 
                              gsub(" ", "\\s*", field_name))
      m2 <- regmatches(text, regexec(pattern_loose, text, perl = TRUE))[[1]]
      if (length(m2) >= 2) {
        # Take first line only
        first_line <- strsplit(trimws(m2[2]), "\n")[[1]][1]
        return(trimws(first_line))
      }
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
# Run Standard LLM Interpretations (Annotation + Phenotype)
# ==============================================================================

if (RUN_LLM_ANNOTATION || RUN_LLM_PHENOTYPE) {
  cat("\n=== Running Standard LLM Interpretations ===\n")
  
  old_plan <- future::plan()
  future::plan("sequential")  # Avoid rate limiting
  
  clusters_for_interpretation <- unique(top_markers$cluster)
  
  # Annotation task
  if (RUN_LLM_ANNOTATION) {
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
            zscore = ssgsea_scores[, cluster_id],
            stringsAsFactors = FALSE
          )
          cluster_ssgsea <- ssgsea_df %>% arrange(desc(zscore)) %>% head(20)
        }
        
        interpretation <- interpret_epithelial(
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
  }
  
  # Phenotype task
  if (RUN_LLM_PHENOTYPE) {
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
            zscore = ssgsea_scores[, cluster_id],
            stringsAsFactors = FALSE
          )
          cluster_ssgsea <- ssgsea_df %>% arrange(desc(zscore)) %>% head(20)
        }
        
        interpretation <- interpret_epithelial(
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
  }
  
  future::plan(old_plan)
}

# ==============================================================================
# Configure DeepSeek API for interpret_agent
# ==============================================================================

if (RUN_INTERPRET_AGENT) {
  cat("\n=== Configuring DeepSeek API for interpret_agent ===\n")
  
  if (is.null(DEEPSEEK_API_KEY) || DEEPSEEK_API_KEY == "") {
    stop("DEEPSEEK_API_KEY not set! Please set: export DEEPSEEK_API_KEY='your-key'")
  }
  
  fanyi::set_translate_option(
    key = DEEPSEEK_API_KEY,
    source = "deepseek"
  )
  
  cat("[OK] DeepSeek API configured\n\n")
  source('/home/h2048/script/R/interpret.R')
  source('/home/h2048/script/R/interpret_agent_hotfix.R')
}

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
  
  if (!is.null(go_bp_gmt) && nrow(go_bp_gmt) > 0) {
    term2gene_blocks$GO_BP <- go_bp_gmt %>%
      mutate(term = paste0("GO_BP|", term))
  }
  if (!is.null(go_mf_gmt) && nrow(go_mf_gmt) > 0) {
    term2gene_blocks$GO_MF <- go_mf_gmt %>%
      mutate(term = paste0("GO_MF|", term))
  }
  if (!is.null(go_cc_gmt) && nrow(go_cc_gmt) > 0) {
    term2gene_blocks$GO_CC <- go_cc_gmt %>%
      mutate(term = paste0("GO_CC|", term))
  }
  
  if (!is.null(epithelial_markers_term2gene) && nrow(epithelial_markers_term2gene) > 0) {
    term2gene_blocks$EPITHELIAL <- epithelial_markers_term2gene %>%
      mutate(term = paste0("EPITHELIAL|", term))
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
    # Filter to genes present in expression matrix
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
          context = epithelial_context,
          n_pathways = INTERPRET_AGENT_N_PATHWAYS,
          model = LLM_MODEL,
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
# Generate Summary Report
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

report_file <- file.path(OUTPUT_DIR, "REPORT.md")

# ⭐ P2-4 FIX: on.exit() to ensure sink is always restored
on.exit({
  tryCatch(sink(), error = function(e) NULL)
}, add = TRUE)

tryCatch({
  sink(report_file)
  
  cat("# Epithelial Cell Subcluster Interpretation Report v3.2.1-EPITHELIAL-PRODUCTION\n\n")
  cat("**Generated:** ", format(Sys.time()), "\n\n")
  cat("**Pipeline:** Epithelial cell specialized v3.2.1 (All P0/P1 fixes applied)\n\n")
  cat("---\n\n")
  
  # Dataset Summary
  cat("## Dataset Summary\n\n")
  cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("- Genes: %d\n", nrow(seurat_obj)))
  cat(sprintf("- Subclusters analyzed: %d\n", length(unique(top_markers$cluster))))
  cat("\n")
  
  if (has_celltype_column) {
    cat("### Cell Type Distribution\n\n")
    celltype_counts <- table(seurat_obj@meta.data[[CELLTYPE_COLUMN]])
    for (ct in names(celltype_counts)) {
      cat(sprintf("- %s: %d cells\n", ct, celltype_counts[ct]))
    }
    cat("\n")
  }
  
  if (RUN_SSGSEA && length(ssgsea_results_all) > 0) {
    cat(sprintf("- ssGSEA methods: %s\n", paste(names(ssgsea_results_all), collapse = ", ")))
  }
  cat("\n---\n\n")
  
  # Cell Type Annotations
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
  
  # Functional Phenotypes
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
  
  # interpret_agent Results
  if (file.exists(file.path(OUTPUT_DIR, "interpret_agent_results.csv"))) {
    cat("## Deep Mode Interpretations (DeepSeek-Reasoner)\n\n")
    agent_df <- read.csv(file.path(OUTPUT_DIR, "interpret_agent_results.csv"),
                        stringsAsFactors = FALSE)
    
    for (i in 1:nrow(agent_df)) {
      row <- agent_df[i, ]
      cat(sprintf("### %s\n\n", row$Cluster))
      
      if (!is.na(row$Overview) && row$Overview != "") {
        cat("**Overview:**  \n", row$Overview, "\n\n")
      }
      
      if (!is.na(row$Key_Mechanisms) && row$Key_Mechanisms != "" && 
          row$Key_Mechanisms != "See overview") {
        cat("**Key Mechanisms:**  \n", row$Key_Mechanisms, "\n\n")
      }
      
      if (!is.na(row$Regulatory_Drivers) && row$Regulatory_Drivers != "" &&
          row$Regulatory_Drivers != "See overview") {
        cat("**Regulatory Drivers:**  \n", row$Regulatory_Drivers, "\n\n")
      }
      
      if (!is.na(row$Narrative) && row$Narrative != "" &&
          row$Narrative != row$Overview) {
        cat("**Biological Narrative:**  \n", row$Narrative, "\n\n")
      }
      
      cat("---\n\n")
    }
  }
  
  # ssGSEA Summary
  if (RUN_SSGSEA && length(ssgsea_results_all) > 0) {
    cat("## ssGSEA Pathway Analysis Summary\n\n")
    
    for (method_name in names(ssgsea_results_all)) {
      cat(sprintf("### %s Pathways\n\n", toupper(method_name)))
      
      top_file <- file.path(OUTPUT_DIR, "reports",
                           paste0("ssgsea_top_pathways_", method_name, ".csv"))
      
      if (file.exists(top_file)) {
        top_df <- read.csv(top_file, stringsAsFactors = FALSE)
        
        # Show top 3 pathways for first 3 clusters as example
        example_clusters <- unique(top_df$cluster)[1:min(3, length(unique(top_df$cluster)))]
        
        for (cid in example_clusters) {
          cluster_top <- top_df %>%
            filter(cluster == cid) %>%
            head(3)
          
          if (nrow(cluster_top) > 0) {
            cat(sprintf("**%s:**  \n", cid))
            for (j in 1:nrow(cluster_top)) {
              cat(sprintf("- %s (z-score: %.2f)\n",
                         cluster_top$pathway[j],
                         cluster_top$zscore[j]))
            }
            cat("\n")
          }
        }
        
        cat(sprintf("*See full results in: %s*\n\n", basename(top_file)))
      }
      
      cat("---\n\n")
    }
  }
  
  # Pipeline Information
  cat("## Pipeline Details\n\n")
  cat("### Configuration\n\n")
  cat(sprintf("- H5AD Path: %s\n", basename(H5AD_PATH)))
  cat(sprintf("- Output Directory: %s\n", OUTPUT_DIR))
  cat(sprintf("- Cluster Column: %s\n", CLUSTER_COLUMN))
  cat(sprintf("- Top Markers per Cluster: %d\n", TOP_N_MARKERS))
  cat("\n")
  
  cat("### Methods Enabled\n\n")
  cat(sprintf("- Enrichment Analysis: %s\n", ifelse(RUN_ENRICHMENT || LOAD_EXISTING_ENRICHMENT, "Yes", "No")))
  cat(sprintf("- ssGSEA Analysis: %s\n", ifelse(RUN_SSGSEA, "Yes", "No")))
  cat(sprintf("- LLM Annotation: %s\n", ifelse(RUN_LLM_ANNOTATION, "Yes", "No")))
  cat(sprintf("- LLM Phenotype: %s\n", ifelse(RUN_LLM_PHENOTYPE, "Yes", "No")))
  cat(sprintf("- interpret_agent: %s\n", ifelse(RUN_INTERPRET_AGENT, "Yes", "No")))
  if (RUN_INTERPRET_AGENT) {
    cat(sprintf("  • Model: %s\n", LLM_MODEL))
    cat(sprintf("  • Max Retries: %d\n", LLM_MAX_RETRIES))
    cat(sprintf("  • PPI Network: %s\n", ifelse(INTERPRET_AGENT_ADD_PPI, "Enabled", "Disabled")))
  }
  cat("\n")
  
  cat("### Computational Resources\n\n")
  cat(sprintf("- Cores: %d\n", N_CORES))
  cat(sprintf("- Memory Limit: %d GB\n", MEMORY_LIMIT_GB))
  cat("\n")
  
  cat("---\n\n")
  
  cat("## Output Files\n\n")
  cat("### Main Results\n\n")
  cat("- `all_markers.csv` - All significant markers\n")
  cat("- `top_markers_filtered.csv` - Filtered markers for enrichment\n")
  cat("- `filtered_genes_info.csv` - Gene filtering statistics\n")
  if (RUN_LLM_ANNOTATION) {
    cat("- `annotation_results.csv` ⭐⭐⭐ - Cell type annotations\n")
  }
  if (RUN_LLM_PHENOTYPE) {
    cat("- `phenotype_results.csv` ⭐⭐⭐ - Functional phenotypes\n")
  }
  if (RUN_INTERPRET_AGENT) {
    cat("- `interpret_agent_results.csv` ⭐⭐⭐ - Deep mode interpretations\n")
  }
  cat("\n")
  
  cat("### Enrichment Results\n\n")
  cat("- `reports/epithelial_markers_enrich.rds` - Custom epithelial markers\n")
  cat("- `reports/cellmarker_enrich.rds` - CellMarker database\n")
  cat("- `reports/panglaodb_enrich.rds` - PanglaoDB markers\n")
  cat("- `reports/go_bp_enrich.rds` - GO Biological Process\n")
  cat("- `reports/go_mf_enrich.rds` - GO Molecular Function\n")
  cat("- `reports/go_cc_enrich.rds` - GO Cellular Component\n")
  cat("- `reports/hallmark_enrich.rds` - Hallmark pathways\n")
  cat("- `reports/kegg_enrich.rds` - KEGG pathways\n")
  cat("\n")
  
  if (RUN_SSGSEA && length(ssgsea_results_all) > 0) {
    cat("### ssGSEA Results\n\n")
    for (method_name in names(ssgsea_results_all)) {
      cat(sprintf("- `reports/ssgsea_scores_%s.rds` - Full scores\n", method_name))
      cat(sprintf("- `reports/ssgsea_top_pathways_%s.csv` - Top pathways\n", method_name))
    }
    cat("\n")
  }
  
  cat("### Additional Outputs\n\n")
  cat("- `reports/annotation_results.rds` - Full annotation results\n")
  cat("- `reports/phenotype_results.rds` - Full phenotype results\n")
  cat("- `reports/interpret_agent_txt/` - Per-cluster detailed interpretations\n")
  cat("- `reports/llm_raw/` - Raw LLM responses (for debugging)\n")
  cat("\n")
  
  sink()
  cat("[OK] Report saved: REPORT.md\n")
}, error = function(e) {
  cat("[ERROR] Failed to generate report:", conditionMessage(e), "\n")
})

# ==============================================================================
# Complete
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("EPITHELIAL ANALYSIS COMPLETE - v3.2.1-EPITHELIAL-PRODUCTION\n")
cat("================================================================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat("  - REPORT.md ⭐⭐⭐                    Comprehensive summary\n")
cat("  - all_markers.csv                   All significant markers\n")
cat("  - top_markers_filtered.csv          Filtered markers for enrichment\n")
if (RUN_LLM_ANNOTATION) {
  cat("  - annotation_results.csv ⭐⭐⭐       Cell type annotations\n")
}
if (RUN_LLM_PHENOTYPE) {
  cat("  - phenotype_results.csv ⭐⭐⭐        Functional phenotypes\n")
}
if (RUN_INTERPRET_AGENT) {
  cat("  - interpret_agent_results.csv ⭐⭐⭐  Deep mode interpretations\n")
}
if (RUN_SSGSEA && length(ssgsea_results_all) > 0) {
  cat("  - reports/ssgsea_scores_[method].rds  ALL ssGSEA results\n")
  cat("  - reports/ssgsea_top_pathways_[method].csv  Top pathways per method\n")
}
cat("\n")

cat("Fixes Applied in v3.2.1:\n")
cat("  ✅ P0-1: LOAD_EXISTING_ENRICHMENT complete fix (TERM2GENE always rebuilt)\n")
cat("  ✅ P0-2: KEGG filename consistency (kegg_enrich.rds)\n")
cat("  ✅ P0-3: Pseudobulk type consistency (rowMeans with drop=FALSE)\n")
cat("  ✅ P0-4: ssGSEA gene symbol case alignment (internal conversion)\n")
cat("  ✅ P0-5: compareCluster universe parameter added\n")
cat("  ✅ P0-6: fanyi dependency check with requireNamespace\n")
cat("  ✅ P1-2: Counts slot detection before normalization\n")
cat("  ✅ P1-3: API retry with exponential backoff (3 attempts)\n")
cat("  ✅ P1-4: Enhanced regex parsing for field extraction\n")
cat("  ✅ P2-1: genes_to_filter deduplication and ignore.case\n")
cat("  ✅ P2-2: Removed CC10 from epithelial marker database\n")
cat("  ✅ P2-4: REPORT.md sink() protection with on.exit\n")
cat("\n")

cat("Key Features:\n")
cat("  ✅ Epithelial expert context (Basal/Goblet/Ciliated differentiation)\n")
cat("  ✅ Custom epithelial marker database (13 subtypes)\n")
cat("  ✅ ALL ssGSEA methods (hallmark/GO BP/MF/CC/KEGG)\n")
cat("  ✅ Dual-task LLM: Annotation + Phenotype\n")
cat("  ✅ interpret_agent with DeepSeek-Reasoner\n")
cat("  ✅ Load existing enrichment with TERM2GENE rebuild\n")
cat("  ✅ GSVA v1.40+ compatibility\n")
cat("  ✅ Production-grade error handling and retry logic\n")
cat("\n")

cat("Analysis modes:\n")
cat("  1. Annotation mode: Cell type identification with marker genes\n")
cat("  2. Phenotype mode: Functional state characterization\n")
cat("  3. Deep mode (interpret_agent): Network-based interpretation\n")
cat("\n")

cat("================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")