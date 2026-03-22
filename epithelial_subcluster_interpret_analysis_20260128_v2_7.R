#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Cell Subcluster LLM Interpretation - PRODUCTION VERSION v2.7
# ==============================================================================
#
# Version: v2.7-FIXED (2026-01-29)
# Status: Production-ready with GMT-based GO enrichment + FIXED LLM parsing
# Adapted from: bcell_subcluster_interpret_analysis_20260127_v2_1.R (WORKING)
#
# Key Fixes from v2.6 → v2.7:
#   ✅ FIXED parse_interpretation() function (proper regex pattern matching)
#   ✅ FIXED field extraction based on task type (annotation vs phenotype)
#   ✅ FIXED LLM prompt format for better structured responses
#   ✅ Added debug output to verify field extraction
#   ✅ Improved error handling and fallback mechanisms
#
# Key Features:
#   ✅ GMT-based GO enrichment (NO gene ID conversion loss!)
#   ✅ 8 database support (Epithelial, CellMarker, PanglaoDB, GO BP/MF/CC, Hallmark, KEGG)
#   ✅ Complete interpret() field extraction
#   ✅ Annotation + Phenotype dual tasks
#   ✅ Comprehensive CSV outputs and REPORT.md
#
# Epithelial-specific:
#   - 14 major cell types with 48 total subclusters
#   - Uses 'subcluster' column for clustering
#   - Adapted marker databases for respiratory epithelium
#
# ==============================================================================

# ==============================================================================
# Configuration Parameters
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0204/epithelial_interpret_v2_7_FIXED"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

# MSigDB GO GMT Files (v2.7 - for zero gene loss enrichment)
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"

# DeepSeek API Key (from environment variable)
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
if (DEEPSEEK_API_KEY == "") {
  stop(
    "ERROR: Missing environment variable 'DEEPSEEK_API_KEY'. Please set it before running."
  )
}

# Analysis Parameters
N_CORES <- 4
TOP_N_MARKERS <- 50

# ======================================================================
# Deep Mode (interpret_agent) - Dual-track
# ======================================================================
RUN_INTERPRET_AGENT <- TRUE
INTERPRET_AGENT_MODEL <- "deepseek-reasoner"
INTERPRET_AGENT_N_PATHWAYS <- 50
INTERPRET_AGENT_ADD_PPI <- TRUE   # 上皮先别开，PPI 会引入额外依赖/波动


# Column names in Seurat object
CLUSTER_COLUMN <- "subcluster" # Main subcluster identities
CELLTYPE_COLUMN <- "celltypist_pred" # Major cell type labels

# ==============================================================================
# Load Libraries
# ==============================================================================

cat("\n=== Loading Libraries ===\n")

library(reticulate)
library(SCNT)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(future)
library(future.apply)
library(httr)
library(jsonlite)
library(rlang)

# Setup parallel
plan("sequential")
options(future.globals.maxSize = 30 * 1024^3)

# Setup Python
use_condaenv("bbknn_env", required = TRUE)

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)

cat("[OK] Libraries loaded and directories created\n\n")

cat("\n=== Configuring DeepSeek API ===\n")
# DeepSeek API Key (from environment variable)
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
if (nchar(DEEPSEEK_API_KEY) < 10) {
  stop(
    "ERROR: DEEPSEEK_API_KEY environment variable not set or invalid\n",
    "Please set it: export DEEPSEEK_API_KEY='your-key-here'"
  )
}
Sys.setenv(DEEPSEEK_API_KEY = DEEPSEEK_API_KEY)

fanyi::set_translate_option(
  key = DEEPSEEK_API_KEY,
  source = "deepseek"
)

test_response <- tryCatch(
  {
    fanyi::chat_request("test", model = "deepseek-reasoner")
  },
  error = function(e) {
    cat("[ERROR] Failed to connect to DeepSeek API:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(test_response)) {
  cat("[OK] DeepSeek API connection successful\n")
} else {
  stop("Failed to connect to DeepSeek API. Please check your API key.")
}


# ==============================================================================
# Load Seurat Data
# ==============================================================================

cat("=== Loading Seurat Data ===\n")

seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

cat(sprintf(
  "\nLoaded: %d cells x %d genes\n",
  ncol(seurat_obj),
  nrow(seurat_obj)
))

# Check available columns
cat("\nAvailable metadata columns:\n")
print(head(colnames(seurat_obj@meta.data), 20))

# Verify cluster column exists
if (!CLUSTER_COLUMN %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Column '%s' not found in metadata!", CLUSTER_COLUMN))
}

cat(sprintf("\nUsing cluster column: %s\n", CLUSTER_COLUMN))

# Verify celltype column exists (with fallback)
has_celltype_column <- CELLTYPE_COLUMN %in% colnames(seurat_obj@meta.data)
if (!has_celltype_column) {
  cat(sprintf("\n[WARN] Column '%s' not found in metadata!\n", CELLTYPE_COLUMN))
  cat(sprintf(
    "[INFO] Will use '%s' as celltype column for per-celltype analysis\n",
    CLUSTER_COLUMN
  ))
  CELLTYPE_COLUMN_ACTUAL <- CLUSTER_COLUMN
} else {
  cat(sprintf("Using celltype column: %s\n", CELLTYPE_COLUMN))
  CELLTYPE_COLUMN_ACTUAL <- CELLTYPE_COLUMN
}

# Normalize data
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = FALSE
)

# ==============================================================================
# Compute Marker Genes (Parallel)
# ==============================================================================

cat("\n=== Computing Marker Genes (Parallel) ===\n")

Idents(seurat_obj) <- CLUSTER_COLUMN
clusters <- levels(Idents(seurat_obj))

cat(sprintf("Finding markers for %d subclusters...\n", length(clusters)))

# Display cell type distribution
if (has_celltype_column) {
  celltype_dist <- table(seurat_obj@meta.data[[CELLTYPE_COLUMN_ACTUAL]])
  cat("\nCell type distribution:\n")
  print(celltype_dist)
  cat("\n")
}

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
          verbose = FALSE
        )
      },
      error = function(e) {
        cat(sprintf(
          "[WARN] Failed to find markers for %s: %s\n",
          cluster_id,
          conditionMessage(e)
        ))
        NULL
      }
    )
  },
  future.seed = TRUE
)

names(marker_list) <- clusters
marker_list <- marker_list[!sapply(marker_list, is.null)]

all_markers <- bind_rows(lapply(names(marker_list), function(cid) {
  df <- marker_list[[cid]]
  df$cluster <- cid
  df$gene <- rownames(df)
  df
}))

all_markers <- all_markers %>% filter(p_val_adj < 0.05)

cat(sprintf(
  "[OK] Found %d significant markers across %d clusters\n",
  nrow(all_markers),
  length(unique(all_markers$cluster))
))

write.csv(
  all_markers,
  file.path(OUTPUT_DIR, "all_markers.csv"),
  row.names = FALSE
)

# ==============================================================================
# Prepare Top Markers for Enrichment
# ==============================================================================

cat("\n=== Preparing Top Markers ===\n")

# Genes to filter (MT, ribosomal, stress response, cell cycle)
genes_to_filter <- c(
  grep("^MT-", rownames(seurat_obj), value = TRUE),
  grep("^RP[SL]", rownames(seurat_obj), value = TRUE),
  # Immediate early genes
  "FOS",
  "JUN",
  "JUNB",
  "JUND",
  "EGR1",
  "EGR2",
  "EGR3",
  "ZFP36",
  "DUSP1",
  "DUSP2",
  "IER2",
  "IER3",
  "ATF3",
  "BTG2",
  "FOSB",
  "NR4A1",
  "NR4A2",
  "NR4A3",
  # Heat shock proteins
  "HSP90AA1",
  "HSPA1A",
  "HSPA1B",
  "DNAJB1",
  "HSPA6",
  # Cell cycle (G2/M and S phase)
  "MKI67",
  "TOP2A",
  "UBE2C",
  "CENPF",
  "PCNA",
  "MCM2",
  "MCM3",
  "MCM4",
  "MCM5",
  "MCM6",
  "MCM7"
)

cat(sprintf(
  "Filtering %d potentially confounding genes\n",
  length(genes_to_filter)
))

top_markers <- all_markers %>%
  filter(!gene %in% genes_to_filter) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = TOP_N_MARKERS) %>%
  ungroup() %>%
  mutate(gene = toupper(gene)) %>%
  dplyr::select(gene, cluster)

cat(sprintf("Selected top %d clean markers per cluster\n", TOP_N_MARKERS))
cat(sprintf("Total markers for enrichment: %d\n", nrow(top_markers)))

write.csv(
  top_markers,
  file.path(OUTPUT_DIR, "top_markers_filtered.csv"),
  row.names = FALSE
)

# Define universe for enrichment (all detectable genes minus filtered genes)
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

filtered_info <- data.frame(
  category = c(
    "MT genes",
    "Ribosomal genes",
    "Stress/IEG genes",
    "Cell cycle genes",
    "Total filtered"
  ),
  count = c(
    sum(grepl("^MT-", genes_to_filter)),
    sum(grepl("^RP[SL]", genes_to_filter)),
    sum(
      genes_to_filter %in%
        c(
          "FOS",
          "JUN",
          "JUNB",
          "JUND",
          "EGR1",
          "EGR2",
          "EGR3",
          "ZFP36",
          "DUSP1",
          "DUSP2",
          "IER2",
          "IER3",
          "ATF3",
          "BTG2",
          "FOSB",
          "NR4A1",
          "NR4A2",
          "NR4A3"
        )
    ),
    sum(
      genes_to_filter %in%
        c(
          "MKI67",
          "TOP2A",
          "UBE2C",
          "CENPF",
          "PCNA",
          "MCM2",
          "MCM3",
          "MCM4",
          "MCM5",
          "MCM6",
          "MCM7"
        )
    ),
    length(genes_to_filter)
  )
)

write.csv(
  filtered_info,
  file.path(OUTPUT_DIR, "filtered_genes_info.csv"),
  row.names = FALSE
)

# ==============================================================================
# Load CellMarker Database
# ==============================================================================

cat("\n=== Loading CellMarker Database ===\n")

cellmarker_db <- NULL
cellmarker_term2gene <- NULL

cellmarker_db <- tryCatch(
  {
    db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)
    cat(sprintf("[OK] Loaded %d CellMarker entries\n", nrow(db)))
    db
  },
  error = function(e) {
    cat("[WARN] Failed to load CellMarker:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(cellmarker_db)) {
  # Filter for human entries
  cellmarker_db <- cellmarker_db %>%
    filter(grepl("Human", species, ignore.case = TRUE))

  cat(sprintf("[OK] Filtered to %d human entries\n", nrow(cellmarker_db)))

  # Filter for lung and airway related tissues
  tissue_keywords <- c(
    "Lung",
    "Airway",
    "Bronchus",
    "Bronchi",
    "Alveolar",
    "Respiratory",
    "Nasal",
    "Sinus",
    "Trachea",
    "Epithelial",
    "Epithelium"
  )

  tissue_pattern <- paste(tissue_keywords, collapse = "|")
  epithelial_related <- cellmarker_db %>%
    filter(
      grepl(tissue_pattern, tissue_type, ignore.case = TRUE) |
        grepl(tissue_pattern, cancer_type, ignore.case = TRUE) |
        grepl(
          "Epithelial|Basal|Goblet|Ciliated|Club|Secretory",
          cell_name,
          ignore.case = TRUE
        )
    )

  cat(sprintf(
    "[OK] Found %d epithelial/respiratory entries\n",
    nrow(epithelial_related)
  ))

  # Use epithelial subset if sufficient, otherwise use full database
  if (nrow(epithelial_related) >= 50) {
    cellmarker_db <- epithelial_related
    cat("[INFO] Using epithelial-specific subset for enrichment\n")
  } else {
    cat("[INFO] Using full human database for broader coverage\n")
  }

  # Create term2gene mapping (row by row to handle marker field properly)
  term2gene_list <- list()

  for (i in 1:nrow(cellmarker_db)) {
    row <- cellmarker_db[i, ]
    cell_type <- row$cell_name
    markers_raw <- row$marker

    if (!is.na(markers_raw) && markers_raw != "") {
      # Split markers by common delimiters
      markers <- unlist(strsplit(markers_raw, "[,;\\s]+"))
      # Clean up markers
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

  # Validate the data
  if (nrow(cellmarker_term2gene) == 0) {
    cat("[ERROR] CellMarker term2gene is empty after processing!\n")
    cellmarker_term2gene <- NULL
  } else {
    cat(sprintf(
      "[OK] Prepared CellMarker TERM2GENE: %d pairs\n",
      nrow(cellmarker_term2gene)
    ))
    cat(sprintf(
      "    Cell types: %d\n",
      length(unique(cellmarker_term2gene$term))
    ))
    cat(sprintf("    Genes: %d\n", length(unique(cellmarker_term2gene$gene))))
  }
}

# ==============================================================================
# Load PanglaoDB Database
# ==============================================================================

cat("\n=== Loading PanglaoDB Database ===\n")

panglaodb_db <- NULL
panglaodb_term2gene <- NULL

panglaodb_db <- tryCatch(
  {
    db <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)

    # Rename columns (PanglaoDB has spaces in column names)
    setnames(
      db,
      old = c("official gene symbol", "cell type"),
      new = c("gene_symbol", "cell_type"),
      skip_absent = TRUE
    )

    # Filter for human (Hs)
    db <- db %>% filter(grepl("Hs", species, fixed = TRUE))

    cat(sprintf("[OK] Loaded %d human markers from PanglaoDB\n", nrow(db)))
    db
  },
  error = function(e) {
    cat("[WARN] Failed to load PanglaoDB:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(panglaodb_db)) {
  # Filter for epithelial and lung-related cell types
  celltype_keywords <- c(
    "Epithelial",
    "Basal",
    "Goblet",
    "Ciliated",
    "Club",
    "Secretory",
    "AT1",
    "AT2",
    "Alveolar",
    "Airway",
    "Ionocyte",
    "Tuft",
    "Brush",
    "SMG",
    "Gland",
    "Duct",
    "Mucous",
    "Serous",
    "Suprabasal"
  )

  celltype_pattern <- paste(celltype_keywords, collapse = "|")
  db_filtered <- panglaodb_db %>%
    filter(grepl(celltype_pattern, cell_type, ignore.case = TRUE))

  cat(sprintf(
    "[OK] Filtered to %d epithelial entries (from %d total)\n",
    nrow(db_filtered),
    nrow(panglaodb_db)
  ))

  # Use filtered subset if sufficient
  if (nrow(db_filtered) >= 100) {
    panglaodb_db <- db_filtered
    cat("[INFO] Using epithelial-specific subset\n")
  } else {
    cat("[INFO] Using full database for broader coverage\n")
  }

  # Create term2gene mapping
  panglaodb_term2gene <- panglaodb_db %>%
    dplyr::select(cell_type, gene_symbol) %>%
    mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct() %>%
    dplyr::rename(term = cell_type, gene = gene_symbol)

  # Validate the data
  if (nrow(panglaodb_term2gene) == 0) {
    cat("[ERROR] PanglaoDB term2gene is empty after processing!\n")
    panglaodb_term2gene <- NULL
  } else {
    cat(sprintf(
      "[OK] Prepared PanglaoDB TERM2GENE: %d pairs\n",
      nrow(panglaodb_term2gene)
    ))
    cat(sprintf(
      "    Cell types: %d\n",
      length(unique(panglaodb_term2gene$term))
    ))
    cat(sprintf("    Genes: %d\n", length(unique(panglaodb_term2gene$gene))))
  }
}

# ==============================================================================
# Custom Epithelial Markers Database
# ==============================================================================

cat("\n=== Creating Custom Epithelial Markers Database ===\n")

epithelial_markers_db <- data.frame(
  subtype = c(
    # Basal cells
    rep("Basal", 6),
    # Suprabasal
    rep("Suprabasal", 4),
    # Dividing Basal
    rep("Dividing_Basal", 5),
    # Secretory - Goblet
    rep("Secretory_Goblet", 7),
    # Secretory - Club
    rep("Secretory_Club", 5),
    # Ciliated
    rep("Ciliated", 7),
    # Deuterosomal
    rep("Deuterosomal", 4),
    # Ionocyte
    rep("Ionocyte", 4),
    # Brush/Tuft
    rep("Brush_Tuft", 5),
    # AT1
    rep("AT1", 5),
    # AT2
    rep("AT2", 6),
    # SMG Basal
    rep("SMG_Basal", 4),
    # SMG Duct
    rep("SMG_Duct", 5),
    # SMG Mucous
    rep("SMG_Mucous", 4),
    # SMG Serous
    rep("SMG_Serous", 5)
  ),
  markers = c(
    # Basal
    "TP63",
    "KRT5",
    "KRT14",
    "KRT15",
    "NGFR",
    "DLK2",
    # Suprabasal
    "KRT13",
    "KRT4",
    "SPRR1A",
    "SPRR1B",
    # Dividing Basal
    "MKI67",
    "TOP2A",
    "UBE2C",
    "CENPF",
    "STMN1",
    # Secretory - Goblet
    "MUC5AC",
    "MUC5B",
    "TFF3",
    "SPDEF",
    "FOXA3",
    "AGR2",
    "FCGBP",
    # Secretory - Club
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "CYP2F1",
    "LYPD2",
    # Ciliated
    "FOXJ1",
    "RSPH1",
    "DNAI1",
    "DNAH5",
    "TUBA1A",
    "TUBB4B",
    "CAPS",
    # Deuterosomal
    "DEUP1",
    "CCNO",
    "CDC20B",
    "CEP78",
    # Ionocyte
    "FOXI1",
    "CFTR",
    "ATP6V0D2",
    "ATP6V1G3",
    # Brush/Tuft
    "POU2F3",
    "TRPM5",
    "GFI1B",
    "AVIL",
    "ALOX5",
    # AT1
    "AGER",
    "PDPN",
    "HOPX",
    "CAV1",
    "EMP2",
    # AT2
    "SFTPC",
    "SFTPB",
    "SFTPA1",
    "LAMP3",
    "ABCA3",
    "LPCAT1",
    # SMG Basal
    "KRT5",
    "TP63",
    "NGFR",
    "ACTA2",
    # SMG Duct
    "KRT7",
    "KRT19",
    "CFTR",
    "AQP5",
    "SCNN1A",
    # SMG Mucous
    "MUC5B",
    "BPIFA2",
    "PRR4",
    "STATH",
    # SMG Serous
    "LTF",
    "LYZ",
    "DMBT1",
    "AZGP1",
    "PRB3"
  ),
  stringsAsFactors = FALSE
)

epithelial_markers_term2gene <- epithelial_markers_db %>%
  mutate(
    term = subtype,
    gene = toupper(trimws(markers))
  ) %>%
  dplyr::select(term, gene) %>%
  distinct()

cat(sprintf(
  "Created custom epithelial markers: %d unique subtypes, %d markers\n",
  length(unique(epithelial_markers_term2gene$term)),
  nrow(epithelial_markers_term2gene)
))

# ==============================================================================
# Load MSigDB GMT Files
# ==============================================================================

cat("\n=== Loading MSigDB GMT Files ===\n")

# Read main MSigDB GMT once
gmt_all <- NULL
tryCatch(
  {
    gmt_all <- read.gmt(MSIGDB_GMT_PATH)
    gmt_all <- gmt_all %>% mutate(gene = toupper(gene))
    cat(sprintf(
      "[OK] Loaded MSigDB GMT: %d gene sets\n",
      length(unique(gmt_all$term))
    ))
  },
  error = function(e) {
    cat("[WARN] Failed to load MSigDB GMT:", conditionMessage(e), "\n")
  }
)

# Extract Hallmark pathways
hallmark_term2gene <- NULL
if (!is.null(gmt_all)) {
  hallmark_term2gene <- gmt_all %>%
    filter(grepl("^HALLMARK_", term))

  if (nrow(hallmark_term2gene) > 0) {
    cat(sprintf(
      "[OK] Extracted Hallmark: %d pathways\n",
      length(unique(hallmark_term2gene$term))
    ))
  } else {
    cat("[WARN] No Hallmark pathways found in GMT\n")
    hallmark_term2gene <- NULL
  }
}

# Extract KEGG pathways from MSigDB
msigdb_kegg_term2gene <- NULL
if (!is.null(gmt_all)) {
  msigdb_kegg_term2gene <- gmt_all %>%
    filter(grepl("^KEGG_", term))

  if (nrow(msigdb_kegg_term2gene) > 0) {
    cat(sprintf(
      "[OK] Extracted MSigDB KEGG: %d pathways\n",
      length(unique(msigdb_kegg_term2gene$term))
    ))
  } else {
    cat("[WARN] No KEGG pathways found in GMT\n")
    msigdb_kegg_term2gene <- NULL
  }
}

# Load GO terms (v2.7 - GMT-based, zero gene loss)
go_bp_term2gene <- NULL
go_mf_term2gene <- NULL
go_cc_term2gene <- NULL

tryCatch(
  {
    gmt_go <- read.gmt(GMT_GO_ALL)
    gmt_go <- gmt_go %>% mutate(gene = toupper(gene))

    # Split by GO category
    go_bp_term2gene <<- gmt_go %>%
      filter(grepl("^GOBP_", term))
    go_mf_term2gene <<- gmt_go %>%
      filter(grepl("^GOMF_", term))
    go_cc_term2gene <<- gmt_go %>%
      filter(grepl("^GOCC_", term))

    cat(sprintf(
      "[OK] Loaded GO terms from GMT:\n"
    ))
    cat(sprintf(
      "    - GO BP: %d terms\n",
      length(unique(go_bp_term2gene$term))
    ))
    cat(sprintf(
      "    - GO MF: %d terms\n",
      length(unique(go_mf_term2gene$term))
    ))
    cat(sprintf(
      "    - GO CC: %d terms\n",
      length(unique(go_cc_term2gene$term))
    ))
  },
  error = function(e) {
    cat("[WARN] Failed to load GO GMT:", conditionMessage(e), "\n")
  }
)

# ==============================================================================
# Enrichment Analysis Function
# ==============================================================================

run_enrichment <- function(
  gene_list,
  term2gene_df,
  analysis_name = "enrichment",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  universe = NULL,
  minGSSize = 5,
  maxGSSize = 500
) {
  if (is.null(term2gene_df) || nrow(term2gene_df) == 0) {
    cat(sprintf("[SKIP] %s: No term2gene mapping available\n", analysis_name))
    return(NULL)
  }

  cat(sprintf("\n--- Running %s ---\n", analysis_name))

  result <- tryCatch(
    {
      enricher(
        gene = unique(gene_list),
        universe = universe,
        TERM2GENE = term2gene_df,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        pAdjustMethod = "BH",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      )
    },
    error = function(e) {
      cat(sprintf(
        "[ERROR] %s failed: %s\n",
        analysis_name,
        conditionMessage(e)
      ))
      NULL
    }
  )

  if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
    cat(sprintf(
      "[OK] %s: %d significant terms\n",
      analysis_name,
      nrow(as.data.frame(result))
    ))
  } else {
    cat(sprintf("[WARN] %s: No significant terms found\n", analysis_name))
  }

  return(result)
}

# ==============================================================================
# Run Enrichment Analyses (Parallel by cluster)
# ==============================================================================

cat("\n=== Running Enrichment Analyses ===\n")

clusters_for_enrichment <- unique(top_markers$cluster)
cat(sprintf("Processing %d clusters...\n", length(clusters_for_enrichment)))

# Create database list
database_configs <- list(
  list(
    name = "epithelial_markers",
    term2gene = epithelial_markers_term2gene,
    filename = "epithelial_markers"
  ),
  list(
    name = "cellmarker",
    term2gene = cellmarker_term2gene,
    filename = "cellmarker"
  ),
  list(
    name = "panglaodb",
    term2gene = panglaodb_term2gene,
    filename = "panglaodb"
  ),
  list(
    name = "go_bp",
    term2gene = go_bp_term2gene,
    filename = "go_bp"
  ),
  list(
    name = "go_mf",
    term2gene = go_mf_term2gene,
    filename = "go_mf"
  ),
  list(
    name = "go_cc",
    term2gene = go_cc_term2gene,
    filename = "go_cc"
  ),
  list(
    name = "hallmark",
    term2gene = hallmark_term2gene,
    filename = "hallmark"
  ),
  list(
    name = "msigdb_kegg",
    term2gene = msigdb_kegg_term2gene,
    filename = "msigdb_kegg"
  )
)

# Initialize results storage
enrichment_results <- list()

for (db_config in database_configs) {
  db_name <- db_config$name
  cat(sprintf("\n=== Processing Database: %s ===\n", db_name))

  if (is.null(db_config$term2gene)) {
    cat(sprintf("[SKIP] %s: Not available\n", db_name))
    next
  }

  # Run enrichment for each cluster
  cluster_results <- future_lapply(
    clusters_for_enrichment,
    function(cluster_id) {
      genes <- top_markers %>%
        filter(cluster == cluster_id) %>%
        pull(gene) %>%
        unique()

      result <- run_enrichment(
        gene_list = genes,
        term2gene_df = db_config$term2gene,
        analysis_name = sprintf("%s (Cluster %s)", db_name, cluster_id),
        universe = universe_genes
      )

      if (!is.null(result)) {
        result@result$cluster <- cluster_id
      }

      return(result)
    },
    future.seed = TRUE
  )

  names(cluster_results) <- clusters_for_enrichment
  cluster_results <- cluster_results[!sapply(cluster_results, is.null)]

  # Combine results
  if (length(cluster_results) > 0) {
    combined_df <- bind_rows(lapply(cluster_results, function(x) {
      if (!is.null(x) && nrow(as.data.frame(x)) > 0) {
        as.data.frame(x)
      } else {
        NULL
      }
    }))

    if (nrow(combined_df) > 0) {
      enrichment_results[[db_name]] <- combined_df

      # Save results
      write.csv(
        combined_df,
        file.path(
          OUTPUT_DIR,
          "reports",
          paste0(db_config$filename, "_enrich.csv")
        ),
        row.names = FALSE
      )

      # Save RDS
      saveRDS(
        cluster_results,
        file.path(
          OUTPUT_DIR,
          "reports",
          paste0(db_config$filename, "_enrich.rds")
        )
      )

      cat(sprintf(
        "[OK] %s: Saved %d results across %d clusters\n",
        db_name,
        nrow(combined_df),
        length(cluster_results)
      ))
    }
  }
}

# ==============================================================================
# Generate Dotplots
# ==============================================================================

cat("\n=== Generating Dotplots ===\n")

for (db_name in names(enrichment_results)) {
  cat(sprintf("Creating dotplot for %s...\n", db_name))

  tryCatch(
    {
      enrich_df <- enrichment_results[[db_name]]

      # Select top 5 terms per cluster
      plot_data <- enrich_df %>%
        group_by(cluster) %>%
        arrange(p.adjust) %>%
        slice_head(n = 5) %>%
        ungroup()

      if (nrow(plot_data) == 0) {
        cat(sprintf("[SKIP] %s: No data to plot\n", db_name))
        next
      }

      # Create dotplot
      p <- ggplot(plot_data, aes(x = cluster, y = Description)) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = "red", high = "blue") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)
        ) +
        labs(
          title = sprintf("%s Enrichment", db_name),
          x = "Cluster",
          y = "Term",
          color = "Adjusted p-value",
          size = "Gene Count"
        )

      # Calculate dimensions with size limits to prevent oversized PDFs
      n_clusters <- length(unique(plot_data$cluster))
      n_terms <- length(unique(plot_data$Description))

      plot_width <- max(12, min(30, n_clusters * 0.8))
      plot_height <- max(8, min(24, n_terms * 0.3))

      cat(sprintf(
        "[INFO] Dotplot size: %.1f x %.1f inches (%d clusters, %d terms)\n",
        plot_width,
        plot_height,
        n_clusters,
        n_terms
      ))

      ggsave(
        file.path(OUTPUT_DIR, "figures", paste0(db_name, "_dotplot.pdf")),
        p,
        width = plot_width,
        height = plot_height,
        limitsize = FALSE
      )

      cat(sprintf("[OK] Saved dotplot: %s_dotplot.pdf\n", db_name))
    },
    error = function(e) {
      cat(sprintf(
        "[WARN] Failed to create dotplot for %s: %s\n",
        db_name,
        conditionMessage(e)
      ))
    }
  )
}

# ==============================================================================
# LLM Interpretation Function (FIXED v2.7)
# ==============================================================================

cat("\n=== Setting up LLM Interpretation (FIXED v2.7) ===\n")

# Function to call DeepSeek API
interpret <- function(
  cluster_id,
  marker_genes,
  enrichment_data = NULL,
  task = "annotation",
  api_key = DEEPSEEK_API_KEY
) {
  cat(sprintf(
    "\n--- LLM Interpretation: Cluster %s (Task: %s) ---\n",
    cluster_id,
    task
  ))

  # Prepare marker genes text
  markers_text <- paste(head(marker_genes, 30), collapse = ", ")

  # Prepare enrichment summary
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
      enrichment_text <- paste(
        "\n\nEnrichment Analysis Results:",
        paste(enrichment_summaries, collapse = "\n"),
        sep = "\n"
      )
    }
  }

  # Construct prompt based on task
  if (task == "annotation") {
    prompt <- sprintf(
      "You are an expert in respiratory epithelial cell biology. Based on the following marker genes and enrichment results for a cell cluster, provide EXACTLY these fields in this format:

Cell Type: [specific cell subtype, e.g., 'Basal cycling', 'Secretory goblet mucin-high']
Confidence: [High/Medium/Low]
Regulatory Drivers: [driver1; driver2; driver3]
Markers: [marker1; marker2; marker3]
Reasoning: [brief explanation]

Marker Genes: %s%s

CRITICAL: You must include ALL FIVE fields (Cell Type, Confidence, Regulatory Drivers, Markers, Reasoning) in your response with the exact format shown above.",
      markers_text,
      enrichment_text
    )
  } else if (task == "phenotype") {
    prompt <- sprintf(
      "You are an expert in respiratory epithelial cell biology. Based on the following marker genes and enrichment results for a cell cluster, provide EXACTLY these fields in this format:

Functional Phenotype: [dominant functional state, e.g., 'Mucus hypersecretion', 'Active ciliogenesis']
Confidence: [High/Medium/Low]
Regulatory Drivers: [driver1; driver2; driver3]
Key Processes: [process1; process2; process3]
Reasoning: [brief explanation]
Network Evidence: [network observations]

Marker Genes: %s%s

CRITICAL: You must include ALL SIX fields (Functional Phenotype, Confidence, Regulatory Drivers, Key Processes, Reasoning, Network Evidence) in your response with the exact format shown above.",
      markers_text,
      enrichment_text
    )
  } else {
    stop("Invalid task type. Must be 'annotation' or 'phenotype'")
  }

  # Call DeepSeek API
  result <- tryCatch(
    {
      response <- httr::POST(
        url = "https://api.deepseek.com/v1/chat/completions",
        httr::add_headers(
          "Authorization" = paste("Bearer", api_key),
          "Content-Type" = "application/json"
        ),
        body = jsonlite::toJSON(
          list(
            model = "deepseek-reasoner",
            messages = list(
              list(role = "user", content = prompt)
            ),
            temperature = 0.3
          ),
          auto_unbox = TRUE
        ),
        encode = "json"
      )

      if (httr::status_code(response) == 200) {
        content <- httr::content(response, as = "parsed")
        interpretation <- content$choices[[1]]$message$content

        # Display truncated version (first 500 chars)
        cat(sprintf(
          "[OK] Received interpretation (%d chars)\n",
          nchar(interpretation)
        ))
        cat("Preview:\n")
        cat(substr(interpretation, 1, 500))
        if (nchar(interpretation) > 500) {
          cat("...[truncated]")
        }
        cat("\n")

        interpretation
      } else {
        cat(sprintf(
          "[ERROR] API returned status %d\n",
          httr::status_code(response)
        ))
        NULL
      }
    },
    error = function(e) {
      cat(sprintf("[ERROR] API call failed: %s\n", conditionMessage(e)))
      NULL
    }
  )

  return(result)
}

# Function to parse interpretation (FIXED v2.7)
parse_interpretation <- function(text, task = "annotation") {
  if (is.null(text) || nchar(text) == 0 || is.na(text)) {
    cat("[WARN] Empty or NA text provided to parse_interpretation\n")
    return(list(raw = NA_character_))
  }

  cat(sprintf(
    "\n[DEBUG] Parsing %s task response (%d chars)\n",
    task,
    nchar(text)
  ))

  # Extract fields using improved regex
  extract_field <- function(text, field_name, debug = FALSE) {
    # Build pattern: field name followed by colon, then capture everything until next field or end
    pattern <- sprintf(
      "%s:\\s*(.+?)(?=\n[A-Z][a-z]+(\\s+[A-Z][a-z]+)*:|$)",
      field_name
    )

    if (debug) {
      cat(sprintf("[DEBUG] Extracting field: %s\n", field_name))
      cat(sprintf("[DEBUG] Pattern: %s\n", pattern))
    }

    match <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]

    if (length(match) >= 2) {
      # Remove the field name prefix and trim
      value <- trimws(gsub(sprintf("^%s:\\s*", field_name), "", match[1]))
      if (debug) {
        cat(sprintf("[DEBUG] Extracted: %s\n", substr(value, 1, 100)))
      }
      return(value)
    } else {
      if (debug) {
        cat(sprintf("[DEBUG] No match found for field: %s\n", field_name))
      }
      return(NA_character_)
    }
  }

  parsed <- list(raw = text)

  if (task == "annotation") {
    parsed$cell_type <- extract_field(text, "Cell Type", debug = TRUE)
    parsed$confidence <- extract_field(text, "Confidence", debug = TRUE)
    parsed$regulatory_drivers <- extract_field(
      text,
      "Regulatory Drivers",
      debug = TRUE
    )
    parsed$markers <- extract_field(text, "Markers", debug = TRUE)
    parsed$reasoning <- extract_field(text, "Reasoning", debug = TRUE)

    # Debug output
    cat(sprintf("[DEBUG] Parsed fields:\n"))
    cat(sprintf(
      "  Cell Type: %s\n",
      ifelse(is.na(parsed$cell_type), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Confidence: %s\n",
      ifelse(is.na(parsed$confidence), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Regulatory Drivers: %s\n",
      ifelse(is.na(parsed$regulatory_drivers), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Markers: %s\n",
      ifelse(is.na(parsed$markers), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Reasoning: %s\n",
      ifelse(is.na(parsed$reasoning), "NA", "FOUND")
    ))
  } else if (task == "phenotype") {
    parsed$functional_phenotype <- extract_field(
      text,
      "Functional Phenotype",
      debug = TRUE
    )
    parsed$confidence <- extract_field(text, "Confidence", debug = TRUE)
    parsed$regulatory_drivers <- extract_field(
      text,
      "Regulatory Drivers",
      debug = TRUE
    )
    parsed$key_processes <- extract_field(text, "Key Processes", debug = TRUE)
    parsed$reasoning <- extract_field(text, "Reasoning", debug = TRUE)
    parsed$network_evidence <- extract_field(
      text,
      "Network Evidence",
      debug = TRUE
    )

    # Debug output
    cat(sprintf("[DEBUG] Parsed fields:\n"))
    cat(sprintf(
      "  Functional Phenotype: %s\n",
      ifelse(is.na(parsed$functional_phenotype), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Confidence: %s\n",
      ifelse(is.na(parsed$confidence), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Regulatory Drivers: %s\n",
      ifelse(is.na(parsed$regulatory_drivers), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Key Processes: %s\n",
      ifelse(is.na(parsed$key_processes), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Reasoning: %s\n",
      ifelse(is.na(parsed$reasoning), "NA", "FOUND")
    ))
    cat(sprintf(
      "  Network Evidence: %s\n",
      ifelse(is.na(parsed$network_evidence), "NA", "FOUND")
    ))
  }

  return(parsed)
}

# ==============================================================================
# Run LLM Interpretations (Annotation Task)
# ==============================================================================

cat("\n=== Running LLM Interpretations (Annotation Task) ===\n")

# Force sequential execution for LLM calls to avoid rate limiting
old_plan <- future::plan()
future::plan("sequential")
cat("[INFO] Using sequential execution for LLM calls to avoid rate limiting\n")

annotation_results <- future_lapply(
  clusters_for_enrichment,
  function(cluster_id) {
    # Get marker genes
    markers <- top_markers %>%
      filter(cluster == cluster_id) %>%
      pull(gene)

    # Get enrichment data for this cluster
    cluster_enrichment <- lapply(enrichment_results, function(db_results) {
      db_results %>% filter(cluster == cluster_id)
    })
    cluster_enrichment <- cluster_enrichment[
      sapply(cluster_enrichment, nrow) > 0
    ]

    # Call LLM
    interpretation <- interpret(
      cluster_id = cluster_id,
      marker_genes = markers,
      enrichment_data = cluster_enrichment,
      task = "annotation"
    )

    # Parse result
    parsed <- parse_interpretation(interpretation, task = "annotation")
    parsed$cluster <- cluster_id

    return(parsed)
  },
  future.seed = TRUE
)

names(annotation_results) <- clusters_for_enrichment

# Restore original plan
future::plan(old_plan)

# Convert to data frame
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

write.csv(
  annotation_df,
  file.path(OUTPUT_DIR, "annotation_results.csv"),
  row.names = FALSE
)

saveRDS(
  annotation_results,
  file.path(OUTPUT_DIR, "reports", "annotation_results.rds")
)

cat(sprintf(
  "\n[OK] Saved annotation results for %d clusters\n",
  nrow(annotation_df)
))

# Display summary
cat("\nAnnotation Summary:\n")
cat(sprintf("  Total clusters: %d\n", nrow(annotation_df)))
cat(sprintf("  Successful parses: %d\n", sum(!is.na(annotation_df$Cell_Type))))
cat(sprintf("  Failed parses: %d\n", sum(is.na(annotation_df$Cell_Type))))

# ==============================================================================
# Run LLM Interpretations (Phenotype Task)
# ==============================================================================

cat("\n=== Running LLM Interpretations (Phenotype Task) ===\n")

# Force sequential execution for LLM calls
old_plan <- future::plan()
future::plan("sequential")
cat("[INFO] Using sequential execution for LLM calls to avoid rate limiting\n")

phenotype_results <- future_lapply(
  clusters_for_enrichment,
  function(cluster_id) {
    # Get marker genes
    markers <- top_markers %>%
      filter(cluster == cluster_id) %>%
      pull(gene)

    # Get enrichment data for this cluster
    cluster_enrichment <- lapply(enrichment_results, function(db_results) {
      db_results %>% filter(cluster == cluster_id)
    })
    cluster_enrichment <- cluster_enrichment[
      sapply(cluster_enrichment, nrow) > 0
    ]

    # Call LLM
    interpretation <- interpret(
      cluster_id = cluster_id,
      marker_genes = markers,
      enrichment_data = cluster_enrichment,
      task = "phenotype"
    )

    # Parse result
    parsed <- parse_interpretation(interpretation, task = "phenotype")
    parsed$cluster <- cluster_id

    return(parsed)
  },
  future.seed = TRUE
)

names(phenotype_results) <- clusters_for_enrichment

# Restore original plan
future::plan(old_plan)

# Convert to data frame
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

write.csv(
  phenotype_df,
  file.path(OUTPUT_DIR, "phenotype_results.csv"),
  row.names = FALSE
)

saveRDS(
  phenotype_results,
  file.path(OUTPUT_DIR, "reports", "phenotype_results.rds")
)

cat(sprintf(
  "\n[OK] Saved phenotype results for %d clusters\n",
  nrow(phenotype_df)
))

# Display summary
cat("\nPhenotype Summary:\n")
cat(sprintf("  Total clusters: %d\n", nrow(phenotype_df)))
cat(sprintf(
  "  Successful parses: %d\n",
  sum(!is.na(phenotype_df$Functional_Phenotype))
))
cat(sprintf(
  "  Failed parses: %d\n",
  sum(is.na(phenotype_df$Functional_Phenotype))
))
# ======================================================================
# Deep Mode: clusterProfiler::interpret_agent (Dual-track)
# ======================================================================

cat("\n=== Deep Mode Interpretation: interpret_agent (Dual-track) ===\n")

if (isTRUE(RUN_INTERPRET_AGENT)) {
  
  if (!requireNamespace("fanyi", quietly = TRUE)) {
    stop("Package 'fanyi' is required by clusterProfiler::interpret_agent(). Please install it.")
  }
  
  # ---- fold-change vector (gene -> logFC)
  lfc_col <- dplyr::case_when(
    "avg_log2FC" %in% colnames(all_markers) ~ "avg_log2FC",
    "avg_logFC"  %in% colnames(all_markers) ~ "avg_logFC",
    TRUE ~ NA_character_
  )
  if (is.na(lfc_col)) {
    cat("[WARN] No avg_log2FC/avg_logFC found in all_markers -> gene_fold_change will be NULL\n")
    gene_fc <- NULL
  } else {
    fc_tbl <- all_markers %>%
      dplyr::mutate(gene = toupper(gene)) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        fc = .data[[lfc_col]][which.max(abs(.data[[lfc_col]]))],
        .groups = "drop"
      )
    gene_fc <- fc_tbl$fc
    names(gene_fc) <- fc_tbl$gene
  }
  
  # ---- combine TERM2GENE with DB prefix
  term2gene_blocks <- list(
    EPITHELIAL = if (!is.null(epithelial_markers_term2gene) && nrow(epithelial_markers_term2gene) > 0)
      epithelial_markers_term2gene %>% dplyr::mutate(term = paste0("EPITHELIAL|", term)) else NULL,
    
    CELLMARKER = if (!is.null(cellmarker_term2gene) && nrow(cellmarker_term2gene) > 0)
      cellmarker_term2gene %>% dplyr::mutate(term = paste0("CELLMARKER|", term)) else NULL,
    
    PANGLAODB = if (!is.null(panglaodb_term2gene) && nrow(panglaodb_term2gene) > 0)
      panglaodb_term2gene %>% dplyr::mutate(term = paste0("PANGLAODB|", term)) else NULL,
    
    GO_BP = if (!is.null(go_bp_term2gene) && nrow(go_bp_term2gene) > 0)
      go_bp_term2gene %>% dplyr::mutate(term = paste0("GO_BP|", term)) else NULL,
    
    GO_MF = if (!is.null(go_mf_term2gene) && nrow(go_mf_term2gene) > 0)
      go_mf_term2gene %>% dplyr::mutate(term = paste0("GO_MF|", term)) else NULL,
    
    GO_CC = if (!is.null(go_cc_term2gene) && nrow(go_cc_term2gene) > 0)
      go_cc_term2gene %>% dplyr::mutate(term = paste0("GO_CC|", term)) else NULL,
    
    HALLMARK = if (!is.null(hallmark_term2gene) && nrow(hallmark_term2gene) > 0)
      hallmark_term2gene %>% dplyr::mutate(term = paste0("HALLMARK|", term)) else NULL,
    
    KEGG = if (!is.null(msigdb_kegg_term2gene) && nrow(msigdb_kegg_term2gene) > 0)
      msigdb_kegg_term2gene %>% dplyr::mutate(term = paste0("KEGG|", term)) else NULL
  )
  
  combined_term2gene <- dplyr::bind_rows(term2gene_blocks) %>%
    dplyr::mutate(
      term = as.character(term),
      gene = toupper(as.character(gene))
    ) %>%
    dplyr::distinct()
  
  if (nrow(combined_term2gene) == 0) {
    cat("[WARN] combined TERM2GENE is empty -> skip interpret_agent\n")
  } else {
    
    combined_universe <- unique(combined_term2gene$gene)
    
    cat(sprintf("[OK] interpret_agent TERM2GENE combined: %d pairs, %d genes\n",
                nrow(combined_term2gene), length(combined_universe)))
    
    # ---- context: 上皮专用（你可以再加疾病/部位信息）
    deep_context <- paste(
      "Respiratory epithelium subclusters (human).",
      "Focus on epithelial subtype identity (basal/suprabasal/secretory/ciliated/ionocyte/tuft/AT1/AT2/SMG) and state (cycling, stress, differentiation, mucin program, ciliogenesis).",
      "Use marker evidence and enriched terms; avoid overinterpreting housekeeping pathways.",
      sep = "\n"
    )
    
    # ---- output dirs
    dir.create(file.path(OUTPUT_DIR, "reports", "interpret_agent_txt"),
               recursive = TRUE, showWarnings = FALSE)
    
    to_scalar <- function(x) {
      if (is.null(x) || length(x) == 0) return("")
      if (is.character(x) && length(x) == 1) return(x)
      if (is.character(x)) return(paste(x, collapse = "; "))
      if (is.list(x)) {
        if (requireNamespace("jsonlite", quietly = TRUE)) {
          return(jsonlite::toJSON(x, auto_unbox = TRUE))
        }
        return(paste(capture.output(str(x, max.level = 2)), collapse = "\n"))
      }
      as.character(x)
    }
    
    agent_results <- list()
    agent_rows <- list()
    
    clusters_for_agent <- sort(unique(top_markers$cluster))
    
    for (cid in clusters_for_agent) {
      cat(sprintf("\n[interpret_agent] Cluster: %s\n", cid))
      
      genes_c <- top_markers %>% dplyr::filter(cluster == cid) %>% dplyr::pull(gene) %>% unique()
      if (length(genes_c) < 10) {
        cat(sprintf("[WARN] Too few genes (%d), skip %s\n", length(genes_c), cid))
        next
      }
      
      er <- tryCatch(
        {
          clusterProfiler::enricher(
            gene = genes_c,
            TERM2GENE = combined_term2gene,
            universe = universe_genes,   # 你前面定义的“可检测基因-过滤基因”
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500
          )
        },
        error = function(e) {
          cat("[WARN] enricher failed:", conditionMessage(e), "\n")
          NULL
        }
      )
      
      if (is.null(er) || nrow(as.data.frame(er)) == 0) {
        cat(sprintf("[WARN] No enriched terms, skip interpret_agent: %s\n", cid))
        next
      }
      
      ia <- tryCatch(
        {
          clusterProfiler::interpret_agent(
            x = er,
            context = deep_context,
            n_pathways = INTERPRET_AGENT_N_PATHWAYS,
            model = INTERPRET_AGENT_MODEL,
            api_key = DEEPSEEK_API_KEY,
            add_ppi = INTERPRET_AGENT_ADD_PPI,
            gene_fold_change = gene_fc
          )
        },
        error = function(e) {
          cat("[ERROR] interpret_agent failed:", conditionMessage(e), "\n")
          NULL
        }
      )
      if (is.null(ia)) next
      
      # normalize output
      res_one <- ia
      if (inherits(ia, "interpretation_list") && length(ia) >= 1) res_one <- ia[[1]]
      if (is.list(res_one)) res_one$cluster <- cid
      
      agent_results[[cid]] <- res_one
      
      # per cluster txt
      txt <- capture.output(print(res_one))
      writeLines(txt, file.path(OUTPUT_DIR, "reports", "interpret_agent_txt", paste0(cid, ".txt")))
      
      # csv row
      agent_rows[[length(agent_rows) + 1]] <- data.frame(
        Cluster = cid,
        Overview = to_scalar(res_one$overview),
        Key_Mechanisms = to_scalar(res_one$key_mechanisms),
        Hypothesis = to_scalar(res_one$hypothesis),
        Narrative = to_scalar(res_one$narrative),
        Regulatory_Drivers = to_scalar(res_one$regulatory_drivers),
        Network_Evidence = to_scalar(res_one$network_evidence),
        Data_Source = to_scalar(res_one$data_source),
        stringsAsFactors = FALSE
      )
      
      cat("[OK] interpret_agent done:", cid, "\n")
    }
    
    saveRDS(agent_results, file.path(OUTPUT_DIR, "reports", "interpret_agent_results.rds"))
    
    if (length(agent_rows) > 0) {
      agent_df <- dplyr::bind_rows(agent_rows)
      
      # ---- join annotation confidence/cell_type if available
      annot_path <- file.path(OUTPUT_DIR, "annotation_results.csv")
      if (file.exists(annot_path)) {
        annot_csv <- read.csv(annot_path, stringsAsFactors = FALSE) %>%
          dplyr::select(Cluster, Cell_Type, Confidence) %>%
          dplyr::rename(Annotation_Cell_Type = Cell_Type, Annotation_Confidence = Confidence)
        
        agent_df <- agent_df %>%
          dplyr::left_join(annot_csv, by = "Cluster")
      }
      
      write.csv(agent_df, file.path(OUTPUT_DIR, "interpret_agent_results.csv"), row.names = FALSE)
      cat(sprintf("\n[OK] interpret_agent CSV exported: %d clusters\n", nrow(agent_df)))
    } else {
      cat("\n[WARN] interpret_agent produced no cluster results.\n")
    }
  }
  
} else {
  cat("[INFO] RUN_INTERPRET_AGENT=FALSE, skip interpret_agent\n")
}

# ==============================================================================
# Per-Celltype Detailed Interpretation
# ==============================================================================

cat("\n=== Generating Per-Celltype Interpretations ===\n")

if (has_celltype_column) {
  celltypes <- unique(seurat_obj@meta.data[[CELLTYPE_COLUMN_ACTUAL]])

  celltype_interpretations <- list()

  for (celltype in celltypes) {
    cat(sprintf("\nProcessing celltype: %s\n", celltype))

    # Get subclusters for this celltype
    celltype_cells <- seurat_obj@meta.data[[CELLTYPE_COLUMN_ACTUAL]] == celltype
    subclusters <- unique(seurat_obj@meta.data[celltype_cells, CLUSTER_COLUMN])

    if (length(subclusters) == 0) {
      cat(sprintf("[SKIP] No subclusters found for %s\n", celltype))
      next
    }

    # Compile summary
    summary_text <- sprintf("# %s - Detailed Subcluster Analysis\n\n", celltype)
    summary_text <- paste0(
      summary_text,
      sprintf("Total subclusters: %d\n\n", length(subclusters))
    )

    for (subcluster in subclusters) {
      # Get annotation
      annot <- annotation_results[[subcluster]]
      pheno <- phenotype_results[[subcluster]]

      summary_text <- paste0(
        summary_text,
        sprintf("\n## Subcluster: %s\n\n", subcluster)
      )

      if (!is.null(annot)) {
        summary_text <- paste0(
          summary_text,
          sprintf("**Cell Type:** %s\n", annot$cell_type %||% "NA")
        )
        summary_text <- paste0(
          summary_text,
          sprintf("**Confidence:** %s\n\n", annot$confidence %||% "NA")
        )
        summary_text <- paste0(
          summary_text,
          sprintf("**Markers:** %s\n\n", annot$markers %||% "NA")
        )
      }

      if (!is.null(pheno)) {
        summary_text <- paste0(
          summary_text,
          sprintf(
            "**Functional Phenotype:** %s\n\n",
            pheno$functional_phenotype %||% "NA"
          )
        )
        summary_text <- paste0(
          summary_text,
          sprintf("**Key Processes:** %s\n\n", pheno$key_processes %||% "NA")
        )
      }

      summary_text <- paste0(summary_text, "---\n")
    }

    # Save to file
    output_file <- file.path(
      OUTPUT_DIR,
      "reports",
      paste0(gsub(" ", "_", celltype), "_interpretation.txt")
    )
    writeLines(summary_text, output_file)

    celltype_interpretations[[celltype]] <- summary_text
    cat(sprintf("[OK] Saved interpretation for %s\n", celltype))
  }

  # Save combined RDS
  saveRDS(
    celltype_interpretations,
    file.path(OUTPUT_DIR, "reports", "celltype_interpretations.rds")
  )

  cat(sprintf(
    "\n[OK] Completed interpretations for %d celltypes\n",
    length(celltype_interpretations)
  ))
}

# ==============================================================================
# Generate Summary Report
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

report_file <- file.path(OUTPUT_DIR, "REPORT.md")

tryCatch(
  {
    sink(report_file)

    cat("# Epithelial Cell Subcluster Interpretation Report v2.7-FIXED\n\n")
    cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
    cat("**Status:** FIXED LLM parsing with improved field extraction\n\n")
    cat(
      "**Databases:** Epithelial Markers, CellMarker, PanglaoDB, GO BP/MF/CC, Hallmark, MSigDB KEGG\n\n"
    )
    cat("**LLM Model:** DeepSeek Chat\n\n")
    cat("---\n\n")

    cat("## Dataset Summary\n\n")
    cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
    cat(sprintf(
      "- Total subclusters: %d\n",
      length(unique(seurat_obj@meta.data[[CLUSTER_COLUMN]]))
    ))
    if (has_celltype_column) {
      cat(sprintf(
        "- Major cell types: %d\n",
        length(unique(seurat_obj@meta.data[[CELLTYPE_COLUMN_ACTUAL]]))
      ))
    }
    cat("\n")

    # Cell type distribution
    if (has_celltype_column) {
      cat("### Cell Type Distribution\n\n")
      celltype_counts <- table(seurat_obj@meta.data[[CELLTYPE_COLUMN_ACTUAL]])
      for (ct in names(sort(celltype_counts, decreasing = TRUE))) {
        cat(sprintf("- %s: %d cells\n", ct, celltype_counts[ct]))
      }
      cat("\n---\n\n")
    }

    # Annotation results
    if (file.exists(file.path(OUTPUT_DIR, "annotation_results.csv"))) {
      cat("## Cell Subtype Annotations (LLM-Generated)\n\n")

      annotation_csv <- read.csv(
        file.path(OUTPUT_DIR, "annotation_results.csv"),
        stringsAsFactors = FALSE
      )

      cat(sprintf(
        "**Success Rate:** %d/%d (%.1f%%)\n\n",
        sum(!is.na(annotation_csv$Cell_Type)),
        nrow(annotation_csv),
        100 * sum(!is.na(annotation_csv$Cell_Type)) / nrow(annotation_csv)
      ))

      for (i in 1:nrow(annotation_csv)) {
        row <- annotation_csv[i, ]

        # Skip if failed to parse
        if (is.na(row$Cell_Type)) {
          cat(sprintf("### %s ❌ PARSING FAILED\n\n", row$Cluster))
          next
        }

        cat(sprintf("### %s\n\n", row$Cluster))
        cat(sprintf("**Cell Type:** %s  \n", row$Cell_Type))
        cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))

        if (!is.na(row$Regulatory_Drivers) && row$Regulatory_Drivers != "NA") {
          cat("**Regulatory Drivers:**  \n")
          cat(
            gsub("; ", "  \n- ", paste0("- ", row$Regulatory_Drivers)),
            "\n\n"
          )
        }

        if (!is.na(row$Markers) && row$Markers != "NA") {
          cat("**Key Markers:**  \n")
          markers <- strsplit(row$Markers, "; ")[[1]]
          top_markers <- markers[1:min(8, length(markers))]
          cat(paste0("- ", top_markers, collapse = "\n"), "\n")
          if (length(markers) > 8) {
            cat(sprintf("  *(and %d more...)*\n", length(markers) - 8))
          }
          cat("\n")
        }

        cat("**Reasoning:**  \n")
        cat(row$Reasoning, "\n\n")
        cat("---\n\n")
      }
    }

    # Phenotype results
    if (file.exists(file.path(OUTPUT_DIR, "phenotype_results.csv"))) {
      cat("## Functional Phenotypes (LLM-Generated)\n\n")

      phenotype_csv <- read.csv(
        file.path(OUTPUT_DIR, "phenotype_results.csv"),
        stringsAsFactors = FALSE
      )

      cat(sprintf(
        "**Success Rate:** %d/%d (%.1f%%)\n\n",
        sum(!is.na(phenotype_csv$Functional_Phenotype)),
        nrow(phenotype_csv),
        100 *
          sum(!is.na(phenotype_csv$Functional_Phenotype)) /
          nrow(phenotype_csv)
      ))

      for (i in 1:nrow(phenotype_csv)) {
        row <- phenotype_csv[i, ]

        # Skip if failed to parse
        if (is.na(row$Functional_Phenotype)) {
          cat(sprintf("### %s ❌ PARSING FAILED\n\n", row$Cluster))
          next
        }

        cat(sprintf("### %s\n\n", row$Cluster))
        cat(sprintf(
          "**Functional Phenotype:** %s  \n",
          row$Functional_Phenotype
        ))
        cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))

        if (!is.na(row$Key_Processes) && row$Key_Processes != "NA") {
          cat("**Key Processes:**  \n")
          processes <- strsplit(row$Key_Processes, "; ")[[1]]
          cat(paste0("- ", processes, collapse = "\n"), "\n\n")
        }

        cat("**Reasoning:**  \n")
        cat(row$Reasoning, "\n\n")
        cat("---\n\n")
      }
    }
    agent_csv_path <- file.path(OUTPUT_DIR, "interpret_agent_results.csv")
    if (file.exists(agent_csv_path)) {
      cat("## Deep Mode Interpretations (interpret_agent)\n\n")
      df <- read.csv(agent_csv_path, stringsAsFactors = FALSE)
      
      for (i in 1:nrow(df)) {
        row <- df[i, ]
        cat(sprintf("### %s\n\n", row$Cluster))
        if (!is.na(row$Annotation_Cell_Type)) {
          cat(sprintf("**Annotation:** %s (%s)\n\n",
                      row$Annotation_Cell_Type, row$Annotation_Confidence))
        }
        cat("**Overview:**  \n"); cat(row$Overview, "\n\n")
        if (!is.na(row$Regulatory_Drivers) && nchar(row$Regulatory_Drivers) > 0) {
          cat("**Regulatory Drivers:**  \n"); cat(row$Regulatory_Drivers, "\n\n")
        }
        cat("**Narrative:**  \n"); cat(row$Narrative, "\n\n")
        cat("---\n\n")
      }
    }
    sink()
    cat("[OK] Report saved: REPORT.md\n")
  },
  error = function(e) {
    tryCatch(sink(), error = function(e) NULL)
    cat("[WARN] Failed to generate report:", conditionMessage(e), "\n")
  }
)

# ==============================================================================
# Complete
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("ANALYSIS COMPLETE - Epithelial Subcluster Interpretation v2.7-FIXED\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat(
  "  - REPORT.md                          Comprehensive interpretation report ⭐ NEW\n"
)
cat(
  "  - annotation_results.csv             Cell_Type + Markers + Success Rate\n"
)
cat(
  "  - phenotype_results.csv              Functional_Phenotype + Key_Processes\n"
)
cat("  - all_markers.csv                    All significant marker genes\n")
cat("  - top_markers_filtered.csv           Filtered markers for enrichment\n")
cat("\n")

cat("Enrichment results:\n")
for (db_name in names(enrichment_results)) {
  cat(sprintf("  - reports/%s_enrich.csv\n", db_name))
}
cat("\n")

cat("Figures:\n")
for (db_name in names(enrichment_results)) {
  cat(sprintf("  - figures/%s_dotplot.pdf\n", db_name))
}
cat("\n")

cat("v2.7-FIXED Improvements:\n")
cat("  ✅ Fixed parse_interpretation() with robust regex patterns\n")
cat("  ✅ Added debug output to track field extraction\n")
cat("  ✅ Improved LLM prompt format for structured responses\n")
cat("  ✅ Enhanced error handling and logging\n")
cat("  ✅ Success rate reporting in REPORT.md\n")
cat("\n")

cat(
  "================================================================================\n"
)
cat("DONE\n")
cat(
  "================================================================================\n"
)
