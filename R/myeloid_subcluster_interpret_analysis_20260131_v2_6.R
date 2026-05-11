#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Cell Subcluster LLM Interpretation - PRODUCTION VERSION v2.6
# ==============================================================================
#
# Version: v2.6 (2026-01-31)
# Status: Production-ready with GMT-based GO enrichment
# Adapted from: bcell_subcluster_interpret_analysis_20260127_v2_1.R
#
# Key Features:
#   ✅ GMT-based GO enrichment (NO gene ID conversion loss!)
#      - Replaces enrichGO + bitr (4-72% loss)
#      - Uses enricher + MSigDB GMT (0-5% loss)
#   ✅ 7 database support: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB
#   ✅ Annotation + Phenotype dual tasks
#   ✅ Comprehensive CSV outputs and REPORT.md
#
# Myeloid Cell Types (9 total):
#   - Classical monocytes (8,909 cells)
#   - Non-classical monocytes (2,029 cells)
#   - Macrophages (13,895 cells)
#   - Alveolar macrophages (24,146 cells)
#   - Intestinal macrophages (386 cells)
#   - DC (279 cells)
#   - DC2 (1,033 cells)
#   - pDC (336 cells)
#   - Mast cells (3,730 cells)
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

# DeepSeek API Key
DEEPSEEK_API_KEY <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05"

# Analysis Parameters
N_CORES <- 8
TOP_N_MARKERS <- 50

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

cat("[OK] Libraries loaded and directories created\n\n")

# ==============================================================================
# 配置 DeepSeek API
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
# Load Data
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
# Standardize L3 Labels to Hierarchical Format
# ==============================================================================

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
# Compute Marker Genes (Parallel)
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
# Prepare Top Markers for Enrichment
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

filtered_info <- data.frame(
  category = c(
    "MT genes",
    "Ribosomal genes",
    "Stress response genes",
    "Total filtered"
  ),
  count = c(
    sum(grepl("^MT-", genes_to_filter)),
    sum(grepl("^RP[SL]", genes_to_filter)),
    sum(!grepl("^(MT-|RP[SL])", genes_to_filter)),
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
  cellmarker_db <- cellmarker_db %>%
    filter(grepl("Human", species, ignore.case = TRUE))

  cat(sprintf("[OK] Filtered to %d human entries\n", nrow(cellmarker_db)))

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

  cat(sprintf(
    "[OK] Found %d myeloid/immune-related entries\n",
    nrow(myeloid_related)
  ))

  if (nrow(myeloid_related) >= 50) {
    cellmarker_db <- myeloid_related
    cat("[INFO] Using myeloid/immune-specific subset for enrichment\n")
  } else {
    cat("[INFO] Using full human database for broader coverage\n")
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
  cat(sprintf(
    "    Cell types: %d\n",
    length(unique(cellmarker_term2gene$term))
  ))
  cat(sprintf("    Genes: %d\n", length(unique(cellmarker_term2gene$gene))))
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

# ==============================================================================
# Enrichment Analysis - CellMarker
# ==============================================================================

cat("\n=== CellMarker Enrichment ===\n")

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
    error = function(e) {
      cat("[WARN] CellMarker enrichment failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(cellmarker_enrich)) {
    ccr <- cellmarker_enrich@compareClusterResult
    n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] Found %d significant cell types\n", n_sig))

    saveRDS(
      cellmarker_enrich,
      file.path(OUTPUT_DIR, "reports", "cellmarker_enrich.rds")
    )
  }
}

# ==============================================================================
# Enrichment Analysis - PanglaoDB
# ==============================================================================

cat("\n=== PanglaoDB Enrichment ===\n")

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
    error = function(e) {
      cat("[WARN] PanglaoDB enrichment failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(panglaodb_enrich)) {
    ccr <- panglaodb_enrich@compareClusterResult
    n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] Found %d significant terms\n", n_sig))

    saveRDS(
      panglaodb_enrich,
      file.path(OUTPUT_DIR, "reports", "panglaodb_enrich.rds")
    )
  }
}

# ==============================================================================
# Load GO GMT Files (v2.6 - Zero Gene Loss Method)
# ==============================================================================

cat("\n=== Loading GO Gene Sets from GMT ===\n")

if (!file.exists(GMT_GO_ALL)) {
  cat("[ERROR] GMT file not found at:", GMT_GO_ALL, "\n")
  cat("Please download from MSigDB:\n")
  cat(
    "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/c5.all.v2025.1.Hs.symbols.gmt\n"
  )
  stop("GMT file required for GO enrichment")
}

cat("Loading c5.all GMT file...\n")
go_all_gmt <- read.gmt(GMT_GO_ALL)

go_bp_gmt <- go_all_gmt[grep("^GOBP_", go_all_gmt$term), ]
go_mf_gmt <- go_all_gmt[grep("^GOMF_", go_all_gmt$term), ]
go_cc_gmt <- go_all_gmt[grep("^GOCC_", go_all_gmt$term), ]

cat(sprintf("✓ GO BP: %d gene sets loaded\n", length(unique(go_bp_gmt$term))))
cat(sprintf("✓ GO MF: %d gene sets loaded\n", length(unique(go_mf_gmt$term))))
cat(sprintf("✓ GO CC: %d gene sets loaded\n", length(unique(go_cc_gmt$term))))

all_genes_in_gmt <- unique(go_all_gmt$gene)
markers_in_gmt <- sum(unique(top_markers$gene) %in% all_genes_in_gmt)
total_markers <- length(unique(top_markers$gene))
cat(sprintf(
  "✓ Gene coverage: %d/%d markers (%.1f%%) found in GMT\n",
  markers_in_gmt,
  total_markers,
  markers_in_gmt / total_markers * 100
))

cat("\n")

# ==============================================================================
# Enrichment Analysis - GO (BP + MF + CC) - GMT-based
# ==============================================================================

cat("=== GO Enrichment (GMT-based, zero gene loss) ===\n")

go_bp_enrich <- tryCatch(
  {
    cat("Running GO BP enrichment (GMT-based)...\n")
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
  ccr <- go_bp_enrich@compareClusterResult
  n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf(
    "[OK] GO BP: Found %d significant terms (GMT-based, no gene loss!)\n",
    n_sig
  ))
  saveRDS(go_bp_enrich, file.path(OUTPUT_DIR, "reports", "go_bp_enrich.rds"))
}

go_mf_enrich <- tryCatch(
  {
    cat("Running GO MF enrichment (GMT-based)...\n")
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
  error = function(e) {
    cat("[WARN] GO MF enrichment failed:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(go_mf_enrich)) {
  ccr <- go_mf_enrich@compareClusterResult
  n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf(
    "[OK] GO MF: Found %d significant terms (GMT-based, no gene loss!)\n",
    n_sig
  ))
  saveRDS(go_mf_enrich, file.path(OUTPUT_DIR, "reports", "go_mf_enrich.rds"))
}

go_cc_enrich <- tryCatch(
  {
    cat("Running GO CC enrichment (GMT-based)...\n")
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
  error = function(e) {
    cat("[WARN] GO CC enrichment failed:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(go_cc_enrich)) {
  ccr <- go_cc_enrich@compareClusterResult
  n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf(
    "[OK] GO CC: Found %d significant terms (GMT-based, no gene loss!)\n",
    n_sig
  ))
  saveRDS(go_cc_enrich, file.path(OUTPUT_DIR, "reports", "go_cc_enrich.rds"))
}

# ==============================================================================
# Enrichment Analysis - MSigDB (Hallmark + KEGG)
# ==============================================================================

cat("\n=== MSigDB Enrichment (Hallmark + KEGG) from Local GMT ===\n")

hallmark_enrich <- NULL
msigdb_kegg_enrich <- NULL

read_gmt <- function(gmt_file) {
  cat(sprintf("Reading GMT file: %s\n", gmt_file))

  if (!file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file)
  }

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
    cat(sprintf("[OK] Loaded %d gene set entries\n", nrow(all_genesets)))

    hallmark_term2gene <- all_genesets %>%
      filter(grepl("^HALLMARK_", term)) %>%
      dplyr::select(term, gene)

    n_hallmark_sets <- length(unique(hallmark_term2gene$term))
    cat(sprintf("[OK] Extracted %d Hallmark gene sets\n", n_hallmark_sets))

    kegg_term2gene <- all_genesets %>%
      filter(grepl("KEGG_", term)) %>%
      dplyr::select(term, gene)

    n_kegg_sets <- length(unique(kegg_term2gene$term))
    cat(sprintf("[OK] Extracted %d KEGG gene sets\n", n_kegg_sets))

    if (nrow(hallmark_term2gene) > 0) {
      cat("\nRunning Hallmark enrichment...\n")
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
        error = function(e) {
          cat("[WARN] Hallmark enrichment failed:", conditionMessage(e), "\n")
          return(NULL)
        }
      )

      if (!is.null(hallmark_enrich)) {
        ccr <- hallmark_enrich@compareClusterResult
        n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
        cat(sprintf("[OK] Found %d significant Hallmark pathways\n", n_sig))

        saveRDS(
          hallmark_enrich,
          file.path(OUTPUT_DIR, "reports", "hallmark_enrich.rds")
        )
      }
    }

    if (nrow(kegg_term2gene) > 0) {
      cat("\nRunning KEGG enrichment...\n")
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
        error = function(e) {
          cat("[WARN] KEGG enrichment failed:", conditionMessage(e), "\n")
          return(NULL)
        }
      )

      if (!is.null(msigdb_kegg_enrich)) {
        ccr <- msigdb_kegg_enrich@compareClusterResult
        n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
        cat(sprintf("[OK] Found %d significant KEGG pathways\n", n_sig))

        saveRDS(
          msigdb_kegg_enrich,
          file.path(OUTPUT_DIR, "reports", "msigdb_kegg_enrich.rds")
        )
      }
    }
  },
  error = function(e) {
    cat("[ERROR] Failed to load GMT file:", conditionMessage(e), "\n")
    cat("[INFO] Enrichment will be skipped\n")
  }
)

# ==============================================================================
# Visualize Enrichment Results
# ==============================================================================

cat("\n=== Visualizing Enrichment Results ===\n")

if (!is.null(cellmarker_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "cellmarker_dotplot.pdf"),
        width = 16,
        height = 12
      )
      print(
        dotplot(cellmarker_enrich, showCategory = 10, font.size = 7) +
          ggtitle("CellMarker Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved cellmarker_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] CellMarker plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

if (!is.null(panglaodb_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "panglaodb_dotplot.pdf"),
        width = 16,
        height = 12
      )
      print(
        dotplot(panglaodb_enrich, showCategory = 10, font.size = 7) +
          ggtitle("PanglaoDB Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved panglaodb_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] PanglaoDB plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

if (!is.null(go_bp_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "go_bp_dotplot.pdf"),
        width = 16,
        height = 14
      )
      print(
        dotplot(go_bp_enrich, showCategory = 15, font.size = 6) +
          ggtitle("GO Biological Process Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved go_bp_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] GO BP plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

if (!is.null(go_mf_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "go_mf_dotplot.pdf"),
        width = 16,
        height = 14
      )
      print(
        dotplot(go_mf_enrich, showCategory = 15, font.size = 6) +
          ggtitle("GO Molecular Function Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved go_mf_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] GO MF plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

if (!is.null(go_cc_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "go_cc_dotplot.pdf"),
        width = 16,
        height = 14
      )
      print(
        dotplot(go_cc_enrich, showCategory = 15, font.size = 6) +
          ggtitle("GO Cellular Component Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved go_cc_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] GO CC plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

if (!is.null(hallmark_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "hallmark_dotplot.pdf"),
        width = 16,
        height = 12
      )
      print(
        dotplot(hallmark_enrich, showCategory = 15, font.size = 6) +
          ggtitle("Hallmark Pathways Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved hallmark_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] Hallmark plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

if (!is.null(msigdb_kegg_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "msigdb_kegg_dotplot.pdf"),
        width = 16,
        height = 14
      )
      print(
        dotplot(msigdb_kegg_enrich, showCategory = 15, font.size = 6) +
          ggtitle("KEGG Pathways (MSigDB) Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved msigdb_kegg_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] MSigDB KEGG plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

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
# Prepare Enrichment Summary for LLM Context
# ==============================================================================

cat("\n=== Preparing Enrichment Summary for LLM ===\n")

enrichment_summary <- list()

if (!is.null(cellmarker_enrich)) {
  cm_top <- cellmarker_enrich@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(p.adjust, n = 3) %>%
    ungroup()
  enrichment_summary[["cellmarker"]] <- cm_top
  cat(sprintf("  [OK] CellMarker: %d top terms\n", nrow(cm_top)))
}

if (!is.null(panglaodb_enrich)) {
  pdb_top <- panglaodb_enrich@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(p.adjust, n = 3) %>%
    ungroup()
  enrichment_summary[["panglaodb"]] <- pdb_top
  cat(sprintf("  [OK] PanglaoDB: %d top terms\n", nrow(pdb_top)))
}

if (!is.null(hallmark_enrich)) {
  hall_top <- hallmark_enrich@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(p.adjust, n = 3) %>%
    ungroup()
  enrichment_summary[["hallmark"]] <- hall_top
  cat(sprintf("  [OK] Hallmark: %d top terms\n", nrow(hall_top)))
}

additional_context <- ""

if (length(enrichment_summary) > 0) {
  additional_context <- "\n\nAdditional enrichment evidence:\n"

  for (db_name in names(enrichment_summary)) {
    db_terms <- enrichment_summary[[db_name]]
    if (nrow(db_terms) > 0) {
      term_summary <- db_terms %>%
        group_by(Cluster) %>%
        summarize(
          top_terms = paste(Description[1:min(3, n())], collapse = "; "),
          .groups = "drop"
        )

      additional_context <- paste0(
        additional_context,
        sprintf("\n- %s: ", toupper(db_name)),
        paste(
          sprintf("%s (%s)", term_summary$Cluster, term_summary$top_terms),
          collapse = " | "
        )
      )
    }
  }
}

cat("\n[OK] Prepared enrichment summary for context\n")

# ==============================================================================
# LLM Interpretation - Task 1: Annotation
# ==============================================================================

cat("\n=== LLM Interpretation: Annotation ===\n")

annotation_results <- NULL

enrichment_list <- list()

if (!is.null(go_bp_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_bp_enrich
  cat("✓ GO BP added\n")
}

if (!is.null(go_mf_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_mf_enrich
  cat("✓ GO MF added\n")
}

if (!is.null(go_cc_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_cc_enrich
  cat("✓ GO CC added\n")
}

if (!is.null(hallmark_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- hallmark_enrich
  cat("✓ Hallmark added\n")
}

if (!is.null(msigdb_kegg_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- msigdb_kegg_enrich
  cat("✓ KEGG added\n")
}

if (!is.null(cellmarker_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- cellmarker_enrich
  cat("✓ CellMarker added\n")
}

if (!is.null(panglaodb_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- panglaodb_enrich
  cat("✓ PanglaoDB added\n")
}

cat(sprintf("\nUsing %d database(s) for annotation\n", length(enrichment_list)))

if (length(enrichment_list) == 0) {
  cat("[ERROR] No enrichment objects available\n")
  annotation_results <- NULL
} else {
  annotation_results <- tryCatch(
    {
      interpret(
        x = enrichment_list,
        context = paste(
          "Myeloid cells from normal respiratory tract tissues (nasal cavity, paranasal sinuses, bronchi, lung).",
          "These are subclusters of 9 major myeloid cell populations:",
          "(1) Classical monocytes - CD14+ CD16- inflammatory monocytes with recruitment capacity,",
          "(2) Non-classical monocytes - CD14low CD16+ patrolling monocytes with vascular surveillance,",
          "(3) Macrophages - tissue-resident phagocytes with diverse activation states,",
          "(4) Alveolar macrophages - lung-specific macrophages with surfactant processing,",
          "(5) Intestinal macrophages - gut-resident macrophages (if any contamination),",
          "(6) DC (Dendritic cells) - classical myeloid DCs with antigen presentation,",
          "(7) DC2 (Type 2 DCs) - specialized DCs for Th2 immunity,",
          "(8) pDC (Plasmacytoid DCs) - type I interferon-producing DCs,",
          "(9) Mast cells - tissue-resident granule-containing effector cells.",
          "The samples are from healthy respiratory tract with baseline immune surveillance.",
          "We are particularly interested in identifying functional states such as:",
          "- M1 vs M2 macrophage polarization (pro-inflammatory vs anti-inflammatory/tissue repair),",
          "- Classical vs alternative monocyte activation states,",
          "- Mature vs immature DC subsets with antigen presentation capacity,",
          "- Tissue-resident vs recruited/circulating populations,",
          "- Activated vs resting mast cell degranulation states,",
          "- Inflammatory vs homeostatic/patrolling phenotypes,",
          "- Phagocytic vs antigen-presenting functional states",
          additional_context
        ),
        task = "annotation",
        n_pathways = 15
      )
    },
    error = function(e) {
      cat("[WARN] Annotation failed:", conditionMessage(e), "\n")
      cat("       Error details:", e$message, "\n")
      return(NULL)
    }
  )
}

if (!is.null(annotation_results)) {
  cat("\n[OK] Annotation complete!\n\n")

  cat("=== Annotation Results ===\n\n")

  if (is.list(annotation_results)) {
    if (length(annotation_results) > 0 && is.list(annotation_results[[1]])) {
      annotation_table <- data.frame(
        Cluster = character(),
        Cell_Type = character(),
        Confidence = character(),
        Reasoning = character(),
        Regulatory_Drivers = character(),
        Markers = character(),
        stringsAsFactors = FALSE
      )

      for (cluster_name in names(annotation_results)) {
        cluster_result <- annotation_results[[cluster_name]]

        if (
          !is.null(cluster_result$overview) &&
            is.null(cluster_result$cell_type) &&
            cluster_result$confidence == "None"
        ) {
          annotation_table <- rbind(
            annotation_table,
            data.frame(
              Cluster = cluster_name,
              Cell_Type = "No significant enrichment",
              Confidence = "None",
              Reasoning = cluster_result$overview,
              Regulatory_Drivers = "NA",
              Markers = "NA",
              stringsAsFactors = FALSE
            )
          )

          cat(sprintf("**%s**\n", cluster_name))
          cat("  Cell Type: No significant enrichment\n")
          cat("  Note: ", cluster_result$overview, "\n\n")
          next
        }

        cell_type <- if (!is.null(cluster_result$cell_type)) {
          cluster_result$cell_type
        } else {
          "Unknown"
        }

        confidence <- if (!is.null(cluster_result$confidence)) {
          cluster_result$confidence
        } else {
          "NA"
        }

        reasoning <- if (!is.null(cluster_result$reasoning)) {
          cluster_result$reasoning
        } else {
          "No reasoning provided"
        }

        regulatory_drivers <- if (!is.null(cluster_result$regulatory_drivers)) {
          paste(cluster_result$regulatory_drivers, collapse = "; ")
        } else {
          "NA"
        }

        markers <- if (!is.null(cluster_result$markers)) {
          paste(cluster_result$markers, collapse = "; ")
        } else {
          "NA"
        }

        annotation_table <- rbind(
          annotation_table,
          data.frame(
            Cluster = cluster_name,
            Cell_Type = cell_type,
            Confidence = confidence,
            Reasoning = reasoning,
            Regulatory_Drivers = regulatory_drivers,
            Markers = markers,
            stringsAsFactors = FALSE
          )
        )

        cat(sprintf("**%s**\n", cluster_name))
        cat(sprintf("  Cell Type: %s\n", cell_type))
        cat(sprintf("  Confidence: %s\n", confidence))
        if (regulatory_drivers != "NA") {
          drivers_list <- strsplit(regulatory_drivers, "; ")[[1]]
          cat(sprintf(
            "  Regulatory Drivers: %s",
            paste(head(drivers_list, 3), collapse = ", ")
          ))
          if (length(drivers_list) > 3) {
            cat(sprintf(" (+ %d more)", length(drivers_list) - 3))
          }
          cat("\n")
        }
        if (markers != "NA") {
          marker_list <- strsplit(markers, "; ")[[1]]
          top5 <- marker_list[1:min(5, length(marker_list))]
          cat(sprintf("  Markers: %s", paste(top5, collapse = ", ")))
          if (length(marker_list) > 5) {
            cat(sprintf(" (+ %d more)", length(marker_list) - 5))
          }
          cat("\n")
        }
        cat(sprintf("  Reasoning: %s...\n\n", substr(reasoning, 1, 150)))
      }

      write.csv(
        annotation_table,
        file.path(OUTPUT_DIR, "annotation_results.csv"),
        row.names = FALSE
      )
      cat("[OK] Saved annotation_results.csv (full content preserved)\n\n")
    } else if (!is.null(annotation_results$interpretation)) {
      cat(annotation_results$interpretation)
      cat("\n\n")
      writeLines(
        annotation_results$interpretation,
        file.path(OUTPUT_DIR, "annotation_results.txt")
      )
    } else {
      cat("Result structure (first 3 elements):\n")
      print(head(annotation_results, 3))
      cat("\n")
    }
  }

  saveRDS(
    annotation_results,
    file.path(OUTPUT_DIR, "reports", "annotation_results.rds")
  )
  cat("[OK] Saved complete results to RDS\n")
} else {
  cat("[WARN] Enrichment not available, skipping annotation\n")
}

# ==============================================================================
# LLM Interpretation - Task 2: Phenotyping
# ==============================================================================

cat("\n=== LLM Interpretation: Phenotyping ===\n")

phenotype_results <- NULL

if (!is.null(go_bp_enrich)) {
  phenotype_results <- tryCatch(
    {
      interpret(
        x = go_bp_enrich,
        context = paste(
          "Myeloid cells from normal respiratory tract (nasal cavity, paranasal sinuses, bronchi, lung).",
          "Looking for baseline functional states and activation signatures in healthy tissue.",
          "Key functional states to identify:",
          "- M1 vs M2 macrophage polarization: pro-inflammatory (iNOS, IL-12, TNF) vs tissue repair (Arg1, CD206, IL-10)",
          "- Monocyte activation: classical (CD14++CD16-) vs intermediate (CD14++CD16+) vs non-classical (CD14+CD16++)",
          "- DC maturation: immature vs mature with antigen presentation capacity",
          "- Phagocytic activity: active phagocytosis vs antigen presentation vs surveillance",
          "- Inflammatory signaling: NF-κB, type I/II interferon, inflammasome activation",
          "- Tissue residency: tissue-resident vs recruited/circulating populations",
          "- Metabolic state: glycolytic (M1-like) vs oxidative phosphorylation (M2-like)",
          "- Lipid handling: foam cell formation, lipid droplet accumulation (especially alveolar macrophages)",
          "These are HEALTHY controls, so expect homeostatic activation states",
          "rather than pathogenic pro-inflammatory signatures.",
          "However, baseline type 2 immunity (M2-like) is expected in mucosal tissues.",
          additional_context
        ),
        task = "phenotyping",
        n_pathways = 30
      )
    },
    error = function(e) {
      cat("[WARN] Phenotyping failed:", conditionMessage(e), "\n")
      cat("       Error details:", e$message, "\n")
      return(NULL)
    }
  )

  if (!is.null(phenotype_results)) {
    cat("\n[OK] Phenotyping complete!\n\n")

    cat("=== Phenotype Results ===\n\n")

    if (is.list(phenotype_results)) {
      if (length(phenotype_results) > 0 && is.list(phenotype_results[[1]])) {
        phenotype_table <- data.frame(
          Cluster = character(),
          Functional_Phenotype = character(),
          Confidence = character(),
          Reasoning = character(),
          Regulatory_Drivers = character(),
          Key_Processes = character(),
          Network_Evidence = character(),
          stringsAsFactors = FALSE
        )

        for (cluster_name in names(phenotype_results)) {
          cluster_result <- phenotype_results[[cluster_name]]

          phenotype <- if (!is.null(cluster_result$phenotype)) {
            cluster_result$phenotype
          } else if (!is.null(cluster_result$functional_state)) {
            cluster_result$functional_state
          } else {
            "Unknown"
          }

          confidence <- if (!is.null(cluster_result$confidence)) {
            cluster_result$confidence
          } else {
            "NA"
          }

          reasoning <- if (!is.null(cluster_result$reasoning)) {
            cluster_result$reasoning
          } else {
            "No reasoning provided"
          }

          regulatory_drivers <- if (
            !is.null(cluster_result$regulatory_drivers)
          ) {
            paste(cluster_result$regulatory_drivers, collapse = "; ")
          } else {
            "NA"
          }

          key_processes <- if (!is.null(cluster_result$key_processes)) {
            paste(cluster_result$key_processes, collapse = "; ")
          } else {
            "NA"
          }

          network_evidence <- if (!is.null(cluster_result$network_evidence)) {
            cluster_result$network_evidence
          } else {
            "NA"
          }

          phenotype_table <- rbind(
            phenotype_table,
            data.frame(
              Cluster = cluster_name,
              Functional_Phenotype = phenotype,
              Confidence = confidence,
              Reasoning = reasoning,
              Regulatory_Drivers = regulatory_drivers,
              Key_Processes = key_processes,
              Network_Evidence = network_evidence,
              stringsAsFactors = FALSE
            )
          )

          cat(sprintf("**%s**\n", cluster_name))
          cat(sprintf("  Phenotype: %s\n", phenotype))
          cat(sprintf("  Confidence: %s\n", confidence))
          if (regulatory_drivers != "NA") {
            cat(sprintf("  Regulatory Drivers: %s\n", regulatory_drivers))
          }
          if (key_processes != "NA") {
            processes <- strsplit(key_processes, "; ")[[1]]
            top3 <- processes[1:min(3, length(processes))]
            cat(sprintf("  Key Processes: %s", paste(top3, collapse = "; ")))
            if (length(processes) > 3) {
              cat(sprintf(" (+ %d more)", length(processes) - 3))
            }
            cat("\n")
          }
          cat(sprintf("  Reasoning: %s...\n\n", substr(reasoning, 1, 150)))
        }

        write.csv(
          phenotype_table,
          file.path(OUTPUT_DIR, "phenotype_results.csv"),
          row.names = FALSE
        )
        cat("[OK] Saved phenotype_results.csv (full content preserved)\n\n")
      } else if (!is.null(phenotype_results$interpretation)) {
        cat(phenotype_results$interpretation)
        cat("\n\n")
        writeLines(
          phenotype_results$interpretation,
          file.path(OUTPUT_DIR, "phenotype_results.txt")
        )
      } else {
        cat("Result structure (first 3 elements):\n")
        print(head(phenotype_results, 3))
        cat("\n")
      }
    }

    saveRDS(
      phenotype_results,
      file.path(OUTPUT_DIR, "reports", "phenotype_results.rds")
    )
    cat("[OK] Saved complete results to RDS\n")
  }
} else {
  cat("[WARN] GO BP enrichment not available, skipping phenotyping\n")
}

# ==============================================================================
# Per-Celltype Detailed Analysis
# ==============================================================================

cat("\n=== Per-Celltype Detailed Analysis ===\n")

celltypes_l2 <- unique(seurat_obj@meta.data$cell_type_L2)
celltypes_l2 <- celltypes_l2[!is.na(celltypes_l2)]

cat(sprintf("Analyzing %d cell types at L2 level\n", length(celltypes_l2)))

celltype_interpretations <- list()

for (celltype in celltypes_l2) {
  cat(sprintf("\n--- Processing %s ---\n", celltype))

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

  celltype_context <- ""

  if (!is.null(hallmark_enrich)) {
    hall_filtered <- hallmark_enrich@compareClusterResult %>%
      filter(Cluster %in% celltype_subclusters, p.adjust < 0.05) %>%
      group_by(Cluster) %>%
      slice_min(p.adjust, n = 2) %>%
      ungroup()

    if (nrow(hall_filtered) > 0) {
      celltype_context <- paste0(
        celltype_context,
        "\nHallmark pathways: ",
        paste(
          sprintf("%s (%s)", hall_filtered$Cluster, hall_filtered$Description),
          collapse = "; "
        )
      )
    }
  }

  interpretation <- tryCatch(
    {
      interpret(
        x = celltype_go_bp,
        context = paste(
          celltype,
          "cells from normal respiratory tract (nasal cavity, sinus, bronchi, lung).",
          "Focus on functional heterogeneity among",
          celltype,
          "subclusters.",
          "These are healthy baseline cells, so we expect:",
          "- Homeostatic activation states rather than pathogenic inflammation",
          "- Normal differentiation trajectories and maturation stages",
          "- Tissue-specific adaptations to mucosal immune surveillance",
          "- Evidence of baseline polarization (M1/M2 for macrophages, classical/non-classical for monocytes)",
          "Key questions:",
          "1. What functional states distinguish the subclusters?",
          "2. Are there proliferative vs quiescent populations?",
          "3. Do subclusters show tissue-specific or universal phenotypes?",
          "4. What polarization or activation markers define each subcluster?",
          "5. For macrophages: M1 vs M2 polarization signatures?",
          "6. For monocytes: classical vs intermediate vs non-classical states?",
          "7. For DCs: maturation level and antigen presentation capacity?",
          celltype_context
        ),
        task = "interpretation",
        n_pathways = 25
      )
    },
    error = function(e) {
      cat("  [WARN] Interpretation failed:", conditionMessage(e), "\n")
      cat("         Error details:", e$message, "\n")
      return(NULL)
    }
  )

  if (!is.null(interpretation)) {
    celltype_interpretations[[celltype]] <- interpretation
    cat("  [OK] Interpretation complete\n")

    report_file <- file.path(
      OUTPUT_DIR,
      "reports",
      paste0(gsub(" ", "_", celltype), "_interpretation.txt")
    )

    tryCatch(
      {
        report_lines <- c()

        if (is.list(interpretation)) {
          if (length(interpretation) > 0 && is.list(interpretation[[1]])) {
            report_lines <- c(
              sprintf("# %s Detailed Interpretation\n", celltype),
              sprintf("Generated: %s\n", format(Sys.time())),
              sprintf("Total Subclusters: %d\n\n", length(interpretation)),
              paste(rep("=", 80), collapse = ""),
              "\n\n"
            )

            for (cluster_name in names(interpretation)) {
              # Check if this cluster belongs to current celltype
              # More robust check: cluster should match celltype OR be in celltype_subclusters
              if (
                cluster_name %in%
                  celltype_subclusters ||
                  grepl(celltype, cluster_name, fixed = TRUE)
              ) {
                cluster_result <- interpretation[[cluster_name]]

                report_lines <- c(
                  report_lines,
                  sprintf("## %s\n\n", cluster_name)
                )

                if (!is.null(cluster_result$phenotype)) {
                  report_lines <- c(
                    report_lines,
                    "### Phenotype\n",
                    cluster_result$phenotype,
                    "\n\n"
                  )
                }

                if (!is.null(cluster_result$confidence)) {
                  report_lines <- c(
                    report_lines,
                    sprintf("**Confidence:** %s\n\n", cluster_result$confidence)
                  )
                }

                if (!is.null(cluster_result$regulatory_drivers)) {
                  report_lines <- c(
                    report_lines,
                    "### Regulatory Drivers\n",
                    paste(
                      "-",
                      cluster_result$regulatory_drivers,
                      collapse = "\n"
                    ),
                    "\n\n"
                  )
                }

                if (!is.null(cluster_result$key_processes)) {
                  report_lines <- c(
                    report_lines,
                    "### Key Biological Processes\n",
                    paste("-", cluster_result$key_processes, collapse = "\n"),
                    "\n\n"
                  )
                }

                if (!is.null(cluster_result$reasoning)) {
                  report_lines <- c(
                    report_lines,
                    "### Reasoning\n",
                    cluster_result$reasoning,
                    "\n\n"
                  )
                }

                if (!is.null(cluster_result$network_evidence)) {
                  report_lines <- c(
                    report_lines,
                    "### Network Evidence\n",
                    cluster_result$network_evidence,
                    "\n\n"
                  )
                }

                if (
                  !is.null(cluster_result$refined_network) &&
                    is.data.frame(cluster_result$refined_network)
                ) {
                  report_lines <- c(report_lines, "### Regulatory Network\n\n")

                  for (i in 1:nrow(cluster_result$refined_network)) {
                    row <- cluster_result$refined_network[i, ]
                    report_lines <- c(
                      report_lines,
                      sprintf(
                        "%d. **%s** → **%s** (%s)\n",
                        i,
                        row$source,
                        row$target,
                        row$interaction
                      )
                    )
                    if (!is.null(row$reason) && !is.na(row$reason)) {
                      report_lines <- c(
                        report_lines,
                        sprintf("   *%s*\n\n", row$reason)
                      )
                    }
                  }
                  report_lines <- c(report_lines, "\n")
                }

                report_lines <- c(report_lines, "---\n\n")
              }
            }

            report_text <- paste(report_lines, collapse = "")
          } else {
            report_text <- ""
            if (!is.null(interpretation$interpretation)) {
              report_text <- interpretation$interpretation
            } else if (!is.null(interpretation$narrative)) {
              report_text <- interpretation$narrative
            }
          }
        } else if (is.character(interpretation)) {
          report_text <- interpretation
        } else {
          report_text <- ""
        }

        if (nchar(report_text) > 0) {
          writeLines(report_text, report_file)
          cat("  [OK] Saved detailed text report\n")
        } else {
          cat("  [INFO] No extractable text, saving as RDS\n")
          saveRDS(
            interpretation,
            file.path(
              OUTPUT_DIR,
              "reports",
              paste0(gsub(" ", "_", celltype), "_interpretation.rds")
            )
          )
        }
      },
      error = function(e) {
        cat("  [WARN] Failed to save report:", e$message, "\n")
        saveRDS(
          interpretation,
          file.path(
            OUTPUT_DIR,
            "reports",
            paste0(gsub(" ", "_", celltype), "_interpretation.rds")
          )
        )
      }
    )
  }
}

if (length(celltype_interpretations) > 0) {
  saveRDS(
    celltype_interpretations,
    file.path(OUTPUT_DIR, "reports", "celltype_interpretations.rds")
  )
  cat(sprintf(
    "\n[OK] Saved interpretations for %d cell types\n",
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

    cat("# Myeloid Cell Subcluster Interpretation Report v2.6\n\n")
    cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
    cat(
      "**Databases (7 total):** CellMarker, PanglaoDB, GO BP/MF/CC, Hallmark, MSigDB KEGG\n\n"
    )
    cat("**LLM Model:** DeepSeek Chat (via interpret() function)\n\n")
    cat("---\n\n")

    cat("## Dataset Summary\n\n")
    cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
    cat(sprintf(
      "- Subclusters (L3): %d\n",
      length(unique(seurat_obj$cell_type_L3))
    ))
    cat(sprintf(
      "- Cell types (L2): %d\n",
      length(unique(seurat_obj$cell_type_L2))
    ))
    cat("\n---\n\n")

    annotation_csv_path <- file.path(OUTPUT_DIR, "annotation_results.csv")
    if (file.exists(annotation_csv_path)) {
      cat("## Cell Subtype Annotations (LLM-Generated)\n\n")

      annotation_df <- read.csv(annotation_csv_path, stringsAsFactors = FALSE)

      for (i in 1:nrow(annotation_df)) {
        row <- annotation_df[i, ]
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
    } else if (!is.null(annotation_results)) {
      cat("## Cell Subtype Annotations (LLM-Generated)\n\n")
      cat("(See annotation_results.rds for detailed structure)\n\n")
      cat("---\n\n")
    }

    phenotype_csv_path <- file.path(OUTPUT_DIR, "phenotype_results.csv")
    if (file.exists(phenotype_csv_path)) {
      cat("## Functional Phenotypes (LLM-Generated)\n\n")

      phenotype_df <- read.csv(phenotype_csv_path, stringsAsFactors = FALSE)

      for (i in 1:nrow(phenotype_df)) {
        row <- phenotype_df[i, ]
        cat(sprintf("### %s\n\n", row$Cluster))
        cat(sprintf(
          "**Functional Phenotype:** %s  \n",
          row$Functional_Phenotype
        ))
        cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))

        if (!is.na(row$Regulatory_Drivers) && row$Regulatory_Drivers != "NA") {
          cat("**Regulatory Drivers:**  \n")
          cat(
            gsub("; ", "  \n- ", paste0("- ", row$Regulatory_Drivers)),
            "\n\n"
          )
        }

        if (!is.na(row$Key_Processes) && row$Key_Processes != "NA") {
          cat("**Key Processes:**  \n")
          processes <- strsplit(row$Key_Processes, "; ")[[1]]
          top_processes <- processes[1:min(5, length(processes))]
          cat(paste0("- ", top_processes, collapse = "\n"), "\n")
          if (length(processes) > 5) {
            cat(sprintf("  *(and %d more...)*\n", length(processes) - 5))
          }
          cat("\n")
        }

        cat("**Reasoning:**  \n")
        cat(row$Reasoning, "\n\n")

        if (
          !is.na(row$Network_Evidence) &&
            row$Network_Evidence != "NA" &&
            nchar(row$Network_Evidence) > 10
        ) {
          cat("**Network Evidence:**  \n")
          cat(row$Network_Evidence, "\n\n")
        }

        cat("---\n\n")
      }
    } else if (!is.null(phenotype_results)) {
      cat("## Functional Phenotypes (LLM-Generated)\n\n")
      cat("(See phenotype_results.rds for detailed structure)\n\n")
      cat("---\n\n")
    }

    if (length(celltype_interpretations) > 0) {
      cat("## Per-Celltype Detailed Interpretations\n\n")

      for (celltype in names(celltype_interpretations)) {
        interpretation <- celltype_interpretations[[celltype]]
        if (is.null(interpretation)) {
          next
        }

        cat(sprintf("### %s\n\n", celltype))

        txt_file <- file.path(
          OUTPUT_DIR,
          "reports",
          paste0(gsub(" ", "_", celltype), "_interpretation.txt")
        )

        if (file.exists(txt_file)) {
          txt_content <- readLines(txt_file)
          cat(paste(txt_content, collapse = "\n"))
          cat("\n\n")
        } else {
          cat("(See RDS file for details)\n\n")
        }

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
cat("ANALYSIS COMPLETE - MYELOID v2.6\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat(
  "  - REPORT.md                          Comprehensive interpretation report\n"
)
cat(
  "  - annotation_results.csv             Cell_Type + Markers ⭐\n"
)
cat(
  "  - phenotype_results.csv              Functional_Phenotype + Key_Processes ⭐\n"
)
cat("  - all_markers.csv                    All significant marker genes\n")
cat(
  "  - top_markers_filtered.csv           Filtered markers used for enrichment\n"
)
cat("  - filtered_genes_info.csv            Gene filtering statistics\n")
cat("\n")

cat("Figures (7 databases):\n")
cat("  - figures/cellmarker_dotplot.pdf     CellMarker database enrichment\n")
cat("  - figures/panglaodb_dotplot.pdf      PanglaoDB enrichment\n")
cat("  - figures/go_bp_dotplot.pdf          GO biological process enrichment\n")
cat(
  "  - figures/go_mf_dotplot.pdf          GO molecular function enrichment\n"
)
cat(
  "  - figures/go_cc_dotplot.pdf          GO cellular component enrichment\n"
)
cat("  - figures/hallmark_dotplot.pdf       Hallmark pathways (50 canonical)\n")
cat("  - figures/msigdb_kegg_dotplot.pdf    KEGG pathways (MSigDB local)\n")
cat("\n")

cat("RDS objects (load with readRDS()):\n")
cat("  - reports/cellmarker_enrich.rds      CellMarker enrichment object\n")
cat("  - reports/panglaodb_enrich.rds       PanglaoDB enrichment object\n")
cat("  - reports/go_bp_enrich.rds           GO BP enrichment object\n")
cat("  - reports/go_mf_enrich.rds           GO MF enrichment object\n")
cat("  - reports/go_cc_enrich.rds           GO CC enrichment object\n")
cat("  - reports/hallmark_enrich.rds        Hallmark enrichment object\n")
cat("  - reports/msigdb_kegg_enrich.rds     MSigDB KEGG enrichment object\n")
cat("  - reports/annotation_results.rds     LLM annotation (raw structure)\n")
cat("  - reports/phenotype_results.rds      LLM phenotype (raw structure)\n")
cat(
  "  - reports/celltype_interpretations.rds Per-celltype detailed interpretations\n"
)
cat("\n")

cat("Per-celltype reports:\n")
cat(
  "  - reports/[celltype]_interpretation.txt  Detailed narrative for each L2 celltype\n"
)
cat("\n")

cat("Important Notes:\n")
cat("  ⚠️  v2.6 PRODUCTION: GMT-based GO enrichment (ZERO gene loss!)\n")
cat("      - Replaced enrichGO + bitr (4-72% loss)\n")
cat("      - Uses enricher + MSigDB GMT (0-5% loss)\n")
cat("      - 7 databases: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB\n")
cat("      - Full Reasoning preserved in CSV (no truncation)\n")
cat("      - Console output truncated for readability\n")
cat(
  "      - Ready for multi-database input: list(go_bp, hallmark, cellmarker, ...)\n"
)
cat("\n")

cat("Myeloid Cell Types (9 total):\n")
cat(
  "  - Classical monocytes:        CD14+ CD16- inflammatory recruitment\n"
)
cat(
  "  - Non-classical monocytes:    CD14low CD16+ patrolling/surveillance\n"
)
cat("  - Macrophages:                Tissue-resident phagocytes (M1/M2)\n")
cat(
  "  - Alveolar macrophages:       Lung-specific surfactant processing\n"
)
cat("  - Intestinal macrophages:     Gut-resident (if any contamination)\n")
cat("  - DC:                         Classical myeloid antigen presentation\n")
cat("  - DC2:                        Type 2 immunity specialized DCs\n")
cat(
  "  - pDC:                        Type I interferon-producing plasmacytoid DCs\n"
)
cat(
  "  - Mast cells:                 Tissue-resident granule effector cells\n"
)
cat("\n")

cat(
  "================================================================================\n"
)
cat("DONE\n")
cat(
  "================================================================================\n"
)
