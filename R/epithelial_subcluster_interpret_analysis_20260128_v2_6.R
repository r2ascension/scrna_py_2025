#!/usr/bin/env Rscript
# ==============================================================================
# Epithelial Cell Subcluster LLM Interpretation - PRODUCTION VERSION v2.6
# ==============================================================================
#
# Version: v2.6 (2026-01-28)
# Status: Production-ready with GMT-based GO enrichment
# Adapted from: bcell_subcluster_interpret_analysis_20260127_v2_1.R
#
# Key Features:
#   ✅ GMT-based GO enrichment (NO gene ID conversion loss!)
#      - Replaces enrichGO + bitr (4-72% loss)
#      - Uses enricher + MSigDB GMT (0-5% loss)
#   ✅ 7-8 database support
#   ✅ Complete interpret() field extraction
#   ✅ Multi-database support: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB
#   ✅ Annotation + Phenotype dual tasks
#   ✅ Comprehensive CSV outputs and REPORT.md
#
# Epithelial-specific:
#   - 14 major cell types with 48 total subclusters
#   - Uses 'subcluster' column for clustering
#   - Adapted marker databases for respiratory epithelium
#
# ==============================================================================
# =======================
# Force single-thread everywhere
# =======================
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# data.table 也有自己的线程池
suppressWarnings({
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(1)
  }
})

# 可选：如果装了 RhpcBLASctl，就更稳
suppressWarnings({
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
})

# ==============================================================================
# Configuration Parameters
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0128/epithelial_interpret_v2_6"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

# MSigDB GO GMT Files (v2.6 - for zero gene loss enrichment)
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"

# DeepSeek API Key
DEEPSEEK_API_KEY <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05"

# Analysis Parameters
N_CORES <- 8
TOP_N_MARKERS <- 50

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
library(fanyi)
library(future)
library(future.apply)

# Setup parallel
plan("sequential") # ✅ 不再 fork/并行
options(future.globals.maxSize = 30 * 1024^3)

# Setup Python
use_condaenv("bbknn_env", required = TRUE)

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)

cat("[OK] Libraries loaded and directories created\n\n")

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
cat(sprintf("Using celltype column: %s\n", CELLTYPE_COLUMN))

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
if (CELLTYPE_COLUMN %in% colnames(seurat_obj@meta.data)) {
  celltype_dist <- table(seurat_obj@meta.data[[CELLTYPE_COLUMN]])
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

cellmarker_term2gene <- NULL

cellmarker_db <- tryCatch(
  {
    db <- fread(CELLMARKER_PATH, header = TRUE, stringsAsFactors = FALSE)

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
    db_filtered <- db %>%
      filter(
        grepl(tissue_pattern, tissueType, ignore.case = TRUE) |
          grepl(tissue_pattern, cancerType, ignore.case = TRUE)
      )

    cat(sprintf(
      "[OK] Loaded %d epithelial/respiratory entries from CellMarker (from %d total)\n",
      nrow(db_filtered),
      nrow(db)
    ))

    # Create term2gene mapping
    cellmarker_term2gene <<- db_filtered %>%
      mutate(
        geneSymbol = strsplit(as.character(geneSymbol), ", ")
      ) %>%
      unnest(geneSymbol) %>%
      mutate(
        term = paste(cellName, tissueType, sep = " | "),
        gene = toupper(trimws(geneSymbol))
      ) %>%
      filter(nchar(gene) > 0) %>%
      dplyr::select(term, gene) %>%
      distinct()

    cat(sprintf(
      "Created CellMarker term2gene: %d unique terms\n",
      length(unique(cellmarker_term2gene$term))
    ))

    db_filtered
  },
  error = function(e) {
    cat("[WARN] Failed to load CellMarker database:", conditionMessage(e), "\n")
    NULL
  }
)

# ==============================================================================
# Load PanglaoDB Database
# ==============================================================================

cat("\n=== Loading PanglaoDB Database ===\n")

panglaodb_term2gene <- NULL

panglaodb_db <- tryCatch(
  {
    db <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)

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
    db_filtered <- db %>%
      filter(grepl(celltype_pattern, cell.type, ignore.case = TRUE))

    cat(sprintf(
      "[OK] Loaded %d epithelial entries from PanglaoDB (from %d total)\n",
      nrow(db_filtered),
      nrow(db)
    ))

    # Create term2gene mapping
    panglaodb_term2gene <<- db_filtered %>%
      mutate(
        term = cell.type,
        gene = toupper(trimws(official.gene.symbol))
      ) %>%
      filter(nchar(gene) > 0) %>%
      dplyr::select(term, gene) %>%
      distinct()

    cat(sprintf(
      "Created PanglaoDB term2gene: %d unique terms\n",
      length(unique(panglaodb_term2gene$term))
    ))

    db_filtered
  },
  error = function(e) {
    cat("[WARN] Failed to load PanglaoDB database:", conditionMessage(e), "\n")
    NULL
  }
)

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

# Load Hallmark
hallmark_term2gene <- NULL
tryCatch(
  {
    gmt <- read.gmt(MSIGDB_GMT_PATH)
    hallmark_term2gene <<- gmt %>%
      filter(grepl("^HALLMARK_", term)) %>%
      mutate(gene = toupper(gene))

    cat(sprintf(
      "[OK] Loaded Hallmark: %d pathways\n",
      length(unique(hallmark_term2gene$term))
    ))
  },
  error = function(e) {
    cat("[WARN] Failed to load Hallmark:", conditionMessage(e), "\n")
  }
)

# Load KEGG from MSigDB
msigdb_kegg_term2gene <- NULL
tryCatch(
  {
    gmt <- read.gmt(MSIGDB_GMT_PATH)
    msigdb_kegg_term2gene <<- gmt %>%
      filter(grepl("^KEGG_", term)) %>%
      mutate(gene = toupper(gene))

    cat(sprintf(
      "[OK] Loaded MSigDB KEGG: %d pathways\n",
      length(unique(msigdb_kegg_term2gene$term))
    ))
  },
  error = function(e) {
    cat("[WARN] Failed to load MSigDB KEGG:", conditionMessage(e), "\n")
  }
)

# Load GO terms (v2.6 NEW - GMT-based, zero gene loss)
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
  qvalueCutoff = 0.2
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
        TERM2GENE = term2gene_df,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        pAdjustMethod = "BH"
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
        analysis_name = sprintf("%s (Cluster %s)", db_name, cluster_id)
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

      ggsave(
        file.path(OUTPUT_DIR, "figures", paste0(db_name, "_dotplot.pdf")),
        p,
        width = max(12, length(unique(plot_data$cluster)) * 0.8),
        height = max(8, length(unique(plot_data$Description)) * 0.3),
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
# LLM Interpretation Function
# ==============================================================================

cat("\n=== Setting up LLM Interpretation ===\n")

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
      "You are an expert in respiratory epithelial cell biology. Based on the following marker genes and enrichment results for a cell cluster, provide:\n\n1. **Cell Type**: The most specific epithelial cell subtype (e.g., 'Basal cycling', 'Secretory goblet mucin-high', 'Ciliated mature', 'AT2 activated', 'SMG serous')\n2. **Confidence**: High/Medium/Low\n3. **Regulatory Drivers**: Key transcription factors or signaling pathways (semicolon-separated)\n4. **Markers**: Top marker genes that define this cell type (semicolon-separated)\n5. **Reasoning**: Brief explanation of your annotation\n\nMarker Genes: %s%s\n\nRespond ONLY in this exact format:\nCell Type: [cell type]\nConfidence: [High/Medium/Low]\nRegulatory Drivers: [driver1; driver2; driver3]\nMarkers: [marker1; marker2; marker3]\nReasoning: [your reasoning]",
      markers_text,
      enrichment_text
    )
  } else if (task == "phenotype") {
    prompt <- sprintf(
      "You are an expert in respiratory epithelial cell biology. Based on the following marker genes and enrichment results for a cell cluster, provide:\n\n1. **Functional Phenotype**: The dominant functional state (e.g., 'Mucus hypersecretion', 'Active ciliogenesis', 'Basal stem-like', 'Surfactant production', 'Glandular differentiation')\n2. **Confidence**: High/Medium/Low\n3. **Regulatory Drivers**: Key transcription factors or pathways (semicolon-separated)\n4. **Key Processes**: Major biological processes active in this cluster (semicolon-separated)\n5. **Reasoning**: Brief explanation\n6. **Network Evidence**: Any gene regulatory networks or pathway interactions you observe\n\nMarker Genes: %s%s\n\nRespond ONLY in this exact format:\nFunctional Phenotype: [phenotype]\nConfidence: [High/Medium/Low]\nRegulatory Drivers: [driver1; driver2; driver3]\nKey Processes: [process1; process2; process3]\nReasoning: [your reasoning]\nNetwork Evidence: [network observations]",
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
            model = "deepseek-chat",
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

# Function to parse interpretation
parse_interpretation <- function(text, task = "annotation") {
  if (is.null(text) || nchar(text) == 0) {
    return(list(raw = NA_character_))
  }

  # Extract fields using regex
  extract_field <- function(text, field_name) {
    pattern <- sprintf(
      "%s:\\s*(.+?)(?=\n[A-Z][a-z]+ [A-Z][a-z]+:|$)",
      field_name
    )
    match <- regmatches(text, regexpr(pattern, text, perl = TRUE))
    if (length(match) > 0) {
      trimws(gsub(sprintf("^%s:\\s*", field_name), "", match))
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
# Run LLM Interpretations (Parallel)
# ==============================================================================

cat("\n=== Running LLM Interpretations (Annotation Task) ===\n")

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
  "[OK] Saved annotation results for %d clusters\n",
  nrow(annotation_df)
))

# ==============================================================================
# Run LLM Interpretations (Phenotype Task)
# ==============================================================================

cat("\n=== Running LLM Interpretations (Phenotype Task) ===\n")

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
  "[OK] Saved phenotype results for %d clusters\n",
  nrow(phenotype_df)
))

# ==============================================================================
# Per-Celltype Detailed Interpretation
# ==============================================================================

cat("\n=== Generating Per-Celltype Interpretations ===\n")

if (CELLTYPE_COLUMN %in% colnames(seurat_obj@meta.data)) {
  celltypes <- unique(seurat_obj@meta.data[[CELLTYPE_COLUMN]])

  celltype_interpretations <- list()

  for (celltype in celltypes) {
    cat(sprintf("\nProcessing celltype: %s\n", celltype))

    # Get subclusters for this celltype
    celltype_cells <- seurat_obj@meta.data[[CELLTYPE_COLUMN]] == celltype
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

    cat("# Epithelial Cell Subcluster Interpretation Report v2.6\n\n")
    cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
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
    cat(sprintf(
      "- Major cell types: %d\n",
      length(unique(seurat_obj@meta.data[[CELLTYPE_COLUMN]]))
    ))
    cat("\n")

    # Cell type distribution
    cat("### Cell Type Distribution\n\n")
    celltype_counts <- table(seurat_obj@meta.data[[CELLTYPE_COLUMN]])
    for (ct in names(sort(celltype_counts, decreasing = TRUE))) {
      cat(sprintf("- %s: %d cells\n", ct, celltype_counts[ct]))
    }
    cat("\n---\n\n")

    # Annotation results
    if (file.exists(file.path(OUTPUT_DIR, "annotation_results.csv"))) {
      cat("## Cell Subtype Annotations (LLM-Generated)\n\n")

      annotation_csv <- read.csv(
        file.path(OUTPUT_DIR, "annotation_results.csv"),
        stringsAsFactors = FALSE
      )

      for (i in 1:nrow(annotation_csv)) {
        row <- annotation_csv[i, ]
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

      for (i in 1:nrow(phenotype_csv)) {
        row <- phenotype_csv[i, ]
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
cat("ANALYSIS COMPLETE - Epithelial Subcluster Interpretation v2.6\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat(
  "  - REPORT.md                          Comprehensive interpretation report\n"
)
cat("  - annotation_results.csv             Cell_Type + Markers\n")
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

cat(
  "================================================================================\n"
)
cat("DONE\n")
cat(
  "================================================================================\n"
)
