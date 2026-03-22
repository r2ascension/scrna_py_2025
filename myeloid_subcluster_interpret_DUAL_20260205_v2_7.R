#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Cell Subcluster Dual LLM Interpretation - v2.7 ENHANCED
# ==============================================================================
#
# Version: v2.7 Enhanced (2026-02-05)
# Status: Production + Deep Analysis Mode
# Features:
#   ✅ Standard interpret() for ALL clusters (fast)
#   ✅ interpret_agent() for ALL celltypes (deep, multi-agent)
#   ✅ PPI network integration in agent mode
#   ✅ 7 databases: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB
#
# New in v2.7 Enhanced:
#   - Dual interpretation strategy (standard + deep)
#   - Agent-based analysis with Cleaner→Detective→Synthesizer pipeline
#   - PPI network evidence integration for ALL celltypes
#   - Comprehensive mechanistic insights
#
# ==============================================================================

# ==============================================================================
# Configuration Parameters
# ==============================================================================
# Input/Output Paths
H5AD_PATH <- "/path/to/your/myeloid_data.h5ad"
OUTPUT_DIR <- "/path/to/output"

# Reference Database Paths
CELLMARKER_PATH <- "/path/to/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/path/to/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/path/to/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/path/to/c5.all.v2025.1.Hs.symbols.gmt"

# DeepSeek API Key
DEEPSEEK_API_KEY <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05"

# Analysis Parameters
N_CORES <- 8
TOP_N_MARKERS <- 50

# ⭐ NEW: Agent Mode Configuration
ENABLE_AGENT_MODE <- TRUE  # Set to FALSE to skip agent analysis
# Agent analysis will be applied to ALL celltypes (no filtering)

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
library(fanyi)
library(future)
library(future.apply)

# Setup parallel
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 10 * 1024^3)

# Setup Python
use_condaenv("bbknn_env", required = TRUE)

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports_agent"), showWarnings = FALSE)

cat("[OK] Libraries loaded and directories created\n\n")

# ==============================================================================
# Configure DeepSeek API
# ==============================================================================

cat("\n=== Configuring DeepSeek API ===\n")

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
# Load Data (same as v2.6)
# ==============================================================================

cat("=== Loading Seurat Data ===\n")

seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

cat(sprintf(
  "\nLoaded: %d cells x %d genes\n",
  ncol(seurat_obj),
  nrow(seurat_obj)
))

# Validate data structure
cat("\n=== Validating Data Structure ===\n")

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
# Standardize L3 Labels (same as v2.6)
# ==============================================================================

cat("\n=== Standardizing L3 Labels ===\n")

cat("\nCurrent L3 labels (before standardization):\n")
print(table(seurat_obj$cell_type_L3, useNA = "ifany"))

seurat_obj$subcluster_id <- as.integer(as.character(seurat_obj$cell_type_L3))

cat(sprintf(
  "\n✓ Preserved original subcluster IDs: %d unique values\n",
  length(unique(seurat_obj$subcluster_id[!is.na(seurat_obj$subcluster_id)]))
))

seurat_obj$cell_type_L3 <- ifelse(
  is.na(seurat_obj$cell_type_L2) | is.na(seurat_obj$subcluster_id),
  NA_character_,
  paste0(seurat_obj$cell_type_L2, "_c", seurat_obj$subcluster_id)
)

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

cat("\nFinal L3 labels (hierarchical format):\n")
l3_table <- table(seurat_obj$cell_type_L3, useNA = "ifany")
print(l3_table)

cat("\n")

# ==============================================================================
# Compute Marker Genes (Parallel) - same as v2.6
# ==============================================================================

cat("\n=== Computing Marker Genes (Parallel) ===\n")

Idents(seurat_obj) <- "cell_type_L3"
clusters <- levels(Idents(seurat_obj))

cat(sprintf("Finding markers for %d clusters...\n", length(clusters)))

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
      error = function(e) NULL
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

cat(sprintf("[OK] Found %d significant markers\n", nrow(all_markers)))

write.csv(
  all_markers,
  file.path(OUTPUT_DIR, "all_markers.csv"),
  row.names = FALSE
)

# ==============================================================================
# Prepare Top Markers (same as v2.6)
# ==============================================================================

cat("\n=== Preparing Top Markers ===\n")

genes_to_filter <- c(
  grep("^MT-", rownames(seurat_obj), value = TRUE),
  grep("^RP[SL]", rownames(seurat_obj), value = TRUE),
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
  "HSP90AA1",
  "HSPA1A",
  "HSPA1B",
  "DNAJB1"
)

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

# ==============================================================================
# Load Databases and Run Enrichment (same as v2.6)
# ==============================================================================
# [KEEPING ALL DATABASE LOADING AND ENRICHMENT CODE FROM v2.6]
# This section includes:
# - CellMarker database
# - PanglaoDB database
# - GO GMT files (BP/MF/CC)
# - MSigDB (Hallmark + KEGG)
# - All enrichment analyses

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
  cellmarker_db <- cellmarker_db %>%
    filter(grepl("Human", species, ignore.case = TRUE))

  myeloid_related <- cellmarker_db %>%
    filter(
      grepl(
        "Monocyte|Macrophage|Dendritic|DC|Myeloid|Mast cell",
        cell_name,
        ignore.case = TRUE
      ) |
        grepl(
          "Blood|Lymph|Spleen|Bone marrow|Lung|Intestine|Liver",
          tissue_type,
          ignore.case = TRUE
        )
    )

  if (nrow(myeloid_related) >= 50) {
    cellmarker_db <- myeloid_related
  }

  term2gene_list <- list()

  for (i in 1:nrow(cellmarker_db)) {
    row <- cellmarker_db[i, ]
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

  cat(sprintf(
    "[OK] Prepared CellMarker TERM2GENE: %d pairs\n",
    nrow(cellmarker_term2gene)
  ))
}

cat("\n=== Loading PanglaoDB Database ===\n")

panglaodb_db <- NULL
panglaodb_term2gene <- NULL

panglaodb_db <- tryCatch(
  {
    db <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)
    setnames(
      db,
      old = c("official gene symbol", "cell type"),
      new = c("gene_symbol", "cell_type"),
      skip_absent = TRUE
    )
    db <- db %>% filter(grepl("Hs", species, fixed = TRUE))
    cat(sprintf("[OK] Loaded %d human markers\n", nrow(db)))
    db
  },
  error = function(e) {
    cat("[WARN] Failed to load PanglaoDB:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(panglaodb_db)) {
  panglaodb_term2gene <- panglaodb_db %>%
    dplyr::select(cell_type, gene_symbol) %>%
    mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct() %>%
    dplyr::rename(term = cell_type, gene = gene_symbol)

  cat(sprintf(
    "[OK] Prepared PanglaoDB TERM2GENE: %d pairs\n",
    nrow(panglaodb_term2gene)
  ))
}

cat("\n=== Loading GO Gene Sets from GMT ===\n")

if (!file.exists(GMT_GO_ALL)) {
  stop("GMT file not found: ", GMT_GO_ALL)
}

cat("Loading c5.all GMT file...\n")
go_all_gmt <- read.gmt(GMT_GO_ALL)

go_bp_gmt <- go_all_gmt[grep("^GOBP_", go_all_gmt$term), ]
go_mf_gmt <- go_all_gmt[grep("^GOMF_", go_all_gmt$term), ]
go_cc_gmt <- go_all_gmt[grep("^GOCC_", go_all_gmt$term), ]

cat(sprintf("✓ GO BP: %d gene sets loaded\n", length(unique(go_bp_gmt$term))))
cat(sprintf("✓ GO MF: %d gene sets loaded\n", length(unique(go_mf_gmt$term))))
cat(sprintf("✓ GO CC: %d gene sets loaded\n", length(unique(go_cc_gmt$term))))

cat("\n=== GO Enrichment (GMT-based) ===\n")

go_bp_enrich <- tryCatch(
  {
    cat("Running GO BP enrichment...\n")
    compareCluster(
      gene ~ cluster,
      data = top_markers,
      fun = enricher,
      TERM2GENE = go_bp_gmt,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
  },
  error = function(e) {
    cat("[WARN] GO BP enrichment failed:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(go_bp_enrich)) {
  cat(sprintf(
    "[OK] GO BP: Found %d significant terms\n",
    sum(go_bp_enrich@compareClusterResult$p.adjust < 0.05)
  ))
  saveRDS(go_bp_enrich, file.path(OUTPUT_DIR, "reports", "go_bp_enrich.rds"))
}

go_mf_enrich <- tryCatch(
  {
    cat("Running GO MF enrichment...\n")
    compareCluster(
      gene ~ cluster,
      data = top_markers,
      fun = enricher,
      TERM2GENE = go_mf_gmt,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
  },
  error = function(e) NULL
)

if (!is.null(go_mf_enrich)) {
  saveRDS(go_mf_enrich, file.path(OUTPUT_DIR, "reports", "go_mf_enrich.rds"))
}

go_cc_enrich <- tryCatch(
  {
    cat("Running GO CC enrichment...\n")
    compareCluster(
      gene ~ cluster,
      data = top_markers,
      fun = enricher,
      TERM2GENE = go_cc_gmt,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
  },
  error = function(e) NULL
)

if (!is.null(go_cc_enrich)) {
  saveRDS(go_cc_enrich, file.path(OUTPUT_DIR, "reports", "go_cc_enrich.rds"))
}

cat("\n=== MSigDB Enrichment (Hallmark + KEGG) ===\n")

hallmark_enrich <- NULL
msigdb_kegg_enrich <- NULL

read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  gene_sets_list <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    list(
      name = parts[1],
      description = parts[2],
      genes = parts[-(1:2)]
    )
  })
  term2gene_list <- lapply(gene_sets_list, function(gs) {
    if (length(gs$genes) > 0) {
      data.frame(
        term = rep(gs$name, length(gs$genes)),
        gene = gs$genes,
        stringsAsFactors = FALSE
      )
    }
  })
  do.call(rbind, term2gene_list)
}

tryCatch(
  {
    all_genesets <- read_gmt(MSIGDB_GMT_PATH)
    
    hallmark_term2gene <- all_genesets %>%
      filter(grepl("^HALLMARK_", term)) %>%
      dplyr::select(term, gene)

    kegg_term2gene <- all_genesets %>%
      filter(grepl("KEGG_", term)) %>%
      dplyr::select(term, gene)

    if (nrow(hallmark_term2gene) > 0) {
      hallmark_enrich <- tryCatch(
        {
          compareCluster(
            gene ~ cluster,
            data = top_markers,
            fun = enricher,
            TERM2GENE = hallmark_term2gene,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500
          )
        },
        error = function(e) NULL
      )

      if (!is.null(hallmark_enrich)) {
        saveRDS(
          hallmark_enrich,
          file.path(OUTPUT_DIR, "reports", "hallmark_enrich.rds")
        )
      }
    }

    if (nrow(kegg_term2gene) > 0) {
      msigdb_kegg_enrich <- tryCatch(
        {
          compareCluster(
            gene ~ cluster,
            data = top_markers,
            fun = enricher,
            TERM2GENE = kegg_term2gene,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500
          )
        },
        error = function(e) NULL
      )

      if (!is.null(msigdb_kegg_enrich)) {
        saveRDS(
          msigdb_kegg_enrich,
          file.path(OUTPUT_DIR, "reports", "msigdb_kegg_enrich.rds")
        )
      }
    }
  },
  error = function(e) {
    cat("[ERROR] Failed to load GMT file\n")
  }
)

cellmarker_enrich <- NULL
if (!is.null(cellmarker_term2gene)) {
  cellmarker_enrich <- tryCatch(
    {
      compareCluster(
        gene ~ cluster,
        data = top_markers,
        fun = enricher,
        TERM2GENE = cellmarker_term2gene,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2
      )
    },
    error = function(e) NULL
  )

  if (!is.null(cellmarker_enrich)) {
    saveRDS(
      cellmarker_enrich,
      file.path(OUTPUT_DIR, "reports", "cellmarker_enrich.rds")
    )
  }
}

panglaodb_enrich <- NULL
if (!is.null(panglaodb_term2gene)) {
  panglaodb_enrich <- tryCatch(
    {
      compareCluster(
        gene ~ cluster,
        data = top_markers,
        fun = enricher,
        TERM2GENE = panglaodb_term2gene,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2
      )
    },
    error = function(e) NULL
  )

  if (!is.null(panglaodb_enrich)) {
    saveRDS(
      panglaodb_enrich,
      file.path(OUTPUT_DIR, "reports", "panglaodb_enrich.rds")
    )
  }
}

# ==============================================================================
# ⭐ NEW: PART 1 - Standard interpret() Analysis (ALL Clusters)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("PART 1: STANDARD interpret() ANALYSIS (Fast Mode)\n")
cat("================================================================================\n\n")

# Prepare enrichment list
enrichment_list <- list()
if (!is.null(go_bp_enrich)) enrichment_list[[length(enrichment_list) + 1]] <- go_bp_enrich
if (!is.null(go_mf_enrich)) enrichment_list[[length(enrichment_list) + 1]] <- go_mf_enrich
if (!is.null(go_cc_enrich)) enrichment_list[[length(enrichment_list) + 1]] <- go_cc_enrich
if (!is.null(hallmark_enrich)) enrichment_list[[length(enrichment_list) + 1]] <- hallmark_enrich
if (!is.null(msigdb_kegg_enrich)) enrichment_list[[length(enrichment_list) + 1]] <- msigdb_kegg_enrich
if (!is.null(cellmarker_enrich)) enrichment_list[[length(enrichment_list) + 1]] <- cellmarker_enrich
if (!is.null(panglaodb_enrich)) enrichment_list[[length(enrichment_list) + 1]] <- panglaodb_enrich

cat(sprintf("Using %d database(s) for standard interpretation\n\n", length(enrichment_list)))

# --- Task 1: Annotation ---
cat("=== Task 1: Cell Type Annotation (Standard) ===\n")

annotation_results_standard <- NULL

if (length(enrichment_list) > 0) {
  annotation_results_standard <- tryCatch(
    {
      interpret(
        x = enrichment_list,
        context = paste(
          "Myeloid cells from normal respiratory tract tissues.",
          "These are subclusters of 9 major myeloid cell populations.",
          "Focus on identifying functional states (M1/M2 polarization, activation states)."
        ),
        task = "annotation",
        n_pathways = 15,
        model = "deepseek-chat",
        api_key = DEEPSEEK_API_KEY
      )
    },
    error = function(e) {
      cat("[WARN] Standard annotation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}

if (!is.null(annotation_results_standard)) {
  saveRDS(
    annotation_results_standard,
    file.path(OUTPUT_DIR, "reports", "annotation_results_standard.rds")
  )
  
  # Save CSV
  if (is.list(annotation_results_standard) && length(annotation_results_standard) > 0) {
    annotation_table <- data.frame(
      Cluster = character(),
      Cell_Type = character(),
      Confidence = character(),
      Reasoning = character(),
      Markers = character(),
      stringsAsFactors = FALSE
    )
    
    for (cluster_name in names(annotation_results_standard)) {
      cluster_result <- annotation_results_standard[[cluster_name]]
      
      annotation_table <- rbind(
        annotation_table,
        data.frame(
          Cluster = cluster_name,
          Cell_Type = if (!is.null(cluster_result$cell_type)) cluster_result$cell_type else "Unknown",
          Confidence = if (!is.null(cluster_result$confidence)) cluster_result$confidence else "NA",
          Reasoning = if (!is.null(cluster_result$reasoning)) cluster_result$reasoning else "",
          Markers = if (!is.null(cluster_result$markers)) paste(cluster_result$markers, collapse = "; ") else "NA",
          stringsAsFactors = FALSE
        )
      )
    }
    
    write.csv(
      annotation_table,
      file.path(OUTPUT_DIR, "annotation_results_standard.csv"),
      row.names = FALSE
    )
  }
  
  cat("[OK] Standard annotation complete\n\n")
}

# --- Task 2: Phenotyping ---
cat("=== Task 2: Functional Phenotyping (Standard) ===\n")

phenotype_results_standard <- NULL

if (!is.null(go_bp_enrich)) {
  phenotype_results_standard <- tryCatch(
    {
      interpret(
        x = go_bp_enrich,
        context = paste(
          "Myeloid cells from normal respiratory tract.",
          "Identify baseline functional states and activation signatures."
        ),
        task = "phenotyping",
        n_pathways = 30,
        model = "deepseek-chat",
        api_key = DEEPSEEK_API_KEY
      )
    },
    error = function(e) {
      cat("[WARN] Standard phenotyping failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}

if (!is.null(phenotype_results_standard)) {
  saveRDS(
    phenotype_results_standard,
    file.path(OUTPUT_DIR, "reports", "phenotype_results_standard.rds")
  )
  
  # Save CSV
  if (is.list(phenotype_results_standard) && length(phenotype_results_standard) > 0) {
    phenotype_table <- data.frame(
      Cluster = character(),
      Functional_Phenotype = character(),
      Confidence = character(),
      Reasoning = character(),
      Key_Processes = character(),
      stringsAsFactors = FALSE
    )
    
    for (cluster_name in names(phenotype_results_standard)) {
      cluster_result <- phenotype_results_standard[[cluster_name]]
      
      phenotype_table <- rbind(
        phenotype_table,
        data.frame(
          Cluster = cluster_name,
          Functional_Phenotype = if (!is.null(cluster_result$phenotype)) cluster_result$phenotype else "Unknown",
          Confidence = if (!is.null(cluster_result$confidence)) cluster_result$confidence else "NA",
          Reasoning = if (!is.null(cluster_result$reasoning)) cluster_result$reasoning else "",
          Key_Processes = if (!is.null(cluster_result$key_processes)) paste(cluster_result$key_processes, collapse = "; ") else "NA",
          stringsAsFactors = FALSE
        )
      )
    }
    
    write.csv(
      phenotype_table,
      file.path(OUTPUT_DIR, "phenotype_results_standard.csv"),
      row.names = FALSE
    )
  }
  
  cat("[OK] Standard phenotyping complete\n\n")
}

# ==============================================================================
# ⭐ NEW: PART 2 - Deep interpret_agent() Analysis (Selected Celltypes)
# ==============================================================================

if (ENABLE_AGENT_MODE) {
  cat("\n")
  cat("================================================================================\n")
  cat("PART 2: DEEP interpret_agent() ANALYSIS (Multi-Agent Mode)\n")
  cat("================================================================================\n\n")
  
  # Get all L2 celltypes
  celltypes_l2 <- unique(seurat_obj@meta.data$cell_type_L2)
  celltypes_l2 <- celltypes_l2[!is.na(celltypes_l2)]
  
  cat(sprintf("Analyzing ALL %d celltypes with agent mode:\n", length(celltypes_l2)))
  cat(paste0("  - ", celltypes_l2, collapse = "\n"), "\n\n")
  
  agent_interpretations <- list()
  
  for (celltype in celltypes_l2) {
    cat(sprintf("\n--- Deep Analysis: %s ---\n", celltype))
    
    # Get subclusters for this celltype
    celltype_subclusters <- unique(
      seurat_obj@meta.data$cell_type_L3[
        seurat_obj@meta.data$cell_type_L2 == celltype
      ]
    )
    celltype_subclusters <- celltype_subclusters[!is.na(celltype_subclusters)]
    
    n_subclusters <- length(celltype_subclusters)
    cat(sprintf("  Subclusters: %d\n", n_subclusters))
    
    if (n_subclusters <= 1) {
      cat("  [SKIP] Only 1 subcluster\n")
      next
    }
    
    # Filter GO BP enrichment for this celltype
    celltype_go_bp <- NULL
    if (!is.null(go_bp_enrich)) {
      go_filtered <- go_bp_enrich@compareClusterResult %>%
        filter(Cluster %in% celltype_subclusters)
      
      if (nrow(go_filtered) > 0) {
        go_subset <- go_bp_enrich
        go_subset@compareClusterResult <- go_filtered
        celltype_go_bp <- go_subset
        cat(sprintf("  GO BP terms: %d\n", nrow(go_filtered)))
      }
    }
    
    if (is.null(celltype_go_bp)) {
      cat("  [SKIP] No GO BP enrichment results\n")
      next
    }
    
    # ⭐ Call interpret_agent() with PPI network integration
    cat("  Calling interpret_agent() (Cleaner→Detective→Synthesizer)...\n")
    
    agent_result <- tryCatch(
      {
        interpret_agent(
          x = celltype_go_bp,
          context = paste(
            celltype,
            "cells from normal respiratory tract.",
            "Focus on functional heterogeneity among subclusters.",
            "Key questions:",
            "- What functional states distinguish the subclusters?",
            "- Are there proliferative vs quiescent populations?",
            "- What polarization/activation markers define each subcluster?",
            "- For macrophages: M1 vs M2 polarization?",
            "- For monocytes: classical vs non-classical states?",
            "- For DCs: maturation level and antigen presentation?"
          ),
          n_pathways = 50,  # More pathways for agent filtering
          model = "deepseek-chat",
          api_key = DEEPSEEK_API_KEY,
          add_ppi = TRUE  # ⭐ Enable PPI network integration
        )
      },
      error = function(e) {
        cat("  [WARN] Agent analysis failed:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    
    if (!is.null(agent_result)) {
      agent_interpretations[[celltype]] <- agent_result
      cat("  [OK] Agent analysis complete\n")
      
      # Save individual report
      report_file <- file.path(
        OUTPUT_DIR,
        "reports_agent",
        paste0(gsub(" ", "_", celltype), "_agent_interpretation.txt")
      )
      
      tryCatch(
        {
          report_lines <- c()
          
          if (is.list(agent_result)) {
            report_lines <- c(
              sprintf("# %s Deep Agent Interpretation\n", celltype),
              sprintf("Generated: %s\n", format(Sys.time())),
              sprintf("Total Subclusters: %d\n\n", length(agent_result)),
              paste(rep("=", 80), collapse = ""),
              "\n\n"
            )
            
            for (cluster_name in names(agent_result)) {
              if (cluster_name %in% celltype_subclusters ||
                  grepl(celltype, cluster_name, fixed = TRUE)) {
                cluster_result <- agent_result[[cluster_name]]
                
                report_lines <- c(
                  report_lines,
                  sprintf("## %s\n\n", cluster_name)
                )
                
                if (!is.null(cluster_result$overview)) {
                  report_lines <- c(
                    report_lines,
                    "### Overview (Agent Synthesizer)\n",
                    cluster_result$overview,
                    "\n\n"
                  )
                }
                
                if (!is.null(cluster_result$regulatory_drivers)) {
                  report_lines <- c(
                    report_lines,
                    "### Regulatory Drivers (Agent Detective)\n",
                    paste(
                      "-",
                      cluster_result$regulatory_drivers,
                      collapse = "\n"
                    ),
                    "\n\n"
                  )
                }
                
                if (!is.null(cluster_result$key_mechanisms)) {
                  report_lines <- c(
                    report_lines,
                    "### Key Mechanisms\n",
                    cluster_result$key_mechanisms,
                    "\n\n"
                  )
                }
                
                if (!is.null(cluster_result$network_evidence)) {
                  report_lines <- c(
                    report_lines,
                    "### PPI Network Evidence (Agent Detective)\n",
                    cluster_result$network_evidence,
                    "\n\n"
                  )
                }
                
                if (!is.null(cluster_result$hypothesis)) {
                  report_lines <- c(
                    report_lines,
                    "### Hypothesis\n",
                    cluster_result$hypothesis,
                    "\n\n"
                  )
                }
                
                if (!is.null(cluster_result$narrative)) {
                  report_lines <- c(
                    report_lines,
                    "### Narrative (Agent Synthesizer)\n",
                    cluster_result$narrative,
                    "\n\n"
                  )
                }
                
                report_lines <- c(report_lines, "---\n\n")
              }
            }
            
            report_text <- paste(report_lines, collapse = "")
            writeLines(report_text, report_file)
            cat("  [OK] Saved agent report\n")
          }
        },
        error = function(e) {
          cat("  [WARN] Failed to save agent report\n")
        }
      )
    }
  }
  
  if (length(agent_interpretations) > 0) {
    saveRDS(
      agent_interpretations,
      file.path(OUTPUT_DIR, "reports_agent", "agent_interpretations_ALL.rds")
    )
    cat(sprintf(
      "\n[OK] Completed agent analysis for %d celltypes\n",
      length(agent_interpretations)
    ))
  }
}

# ==============================================================================
# Final Summary Report
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS COMPLETE - DUAL MODE v2.7\n")
cat("================================================================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key Outputs:\n\n")

cat("STANDARD interpret() Results:\n")
cat("  - annotation_results_standard.csv    Cell type annotations (fast)\n")
cat("  - phenotype_results_standard.csv     Functional phenotypes (fast)\n")
cat("  - reports/annotation_results_standard.rds\n")
cat("  - reports/phenotype_results_standard.rds\n\n")

if (ENABLE_AGENT_MODE) {
  cat("DEEP interpret_agent() Results (ALL celltypes):\n")
  cat("  - reports_agent/[celltype]_agent_interpretation.txt\n")
  cat("  - reports_agent/agent_interpretations_ALL.rds\n\n")
}

cat("Enrichment Data (7 databases):\n")
cat("  - reports/go_bp_enrich.rds\n")
cat("  - reports/go_mf_enrich.rds\n")
cat("  - reports/go_cc_enrich.rds\n")
cat("  - reports/hallmark_enrich.rds\n")
cat("  - reports/msigdb_kegg_enrich.rds\n")
cat("  - reports/cellmarker_enrich.rds\n")
cat("  - reports/panglaodb_enrich.rds\n\n")

cat("Key Features of This Analysis:\n")
cat("  ✅ Standard interpret() for ALL clusters (fast, comprehensive)\n")
cat("  ✅ Deep interpret_agent() for ALL celltypes (slow, mechanistic)\n")
cat("  ✅ Multi-agent pipeline: Cleaner → Detective → Synthesizer\n")
cat("  ✅ PPI network integration in agent mode\n")
cat("  ✅ 7 databases: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB\n\n")

cat("Analysis Summary:\n")
celltypes_l2 <- unique(seurat_obj@meta.data$cell_type_L2)
celltypes_l2 <- celltypes_l2[!is.na(celltypes_l2)]
cat(sprintf("  Total celltypes analyzed: %d\n", length(celltypes_l2)))
cat(sprintf("  Total subclusters: %d\n", length(levels(seurat_obj$cell_type_L3))))
cat("\n")

cat("================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")
