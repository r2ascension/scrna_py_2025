#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Cell Subcluster LLM Interpretation - v1.4 FILENAME FIX
# ==============================================================================
#
# Version: v1.4 (2026-01-30)
# Status: Production-ready with FILENAME SANITIZATION FIX
#
# Key Fix in v1.4:
#   ✅ P0-1 CRITICAL: Sanitize cell type names for filenames
#            - Remove ALL special characters: / + ( ) [ ] { } < > : ; , ? * " ' |
#            - Replace spaces with underscores
#            - Prevents "cannot open the connection" errors
#
# All v1.3 features retained:
#   ✅ Full T/NK cell biological context in interpret() calls
#   ✅ compareCluster-based enrichment (per-cluster, not global pooled)
#   ✅ GMT-based GO enrichment (NO gene ID conversion loss!)
#   ✅ Multi-database support: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB
#   ✅ Annotation + Phenotype dual tasks
#
# ==============================================================================

# ==============================================================================
# Utility Function: Sanitize Filename
# ==============================================================================

sanitize_filename <- function(name) {
  # Remove or replace problematic characters for filenames
  # Priority: Remove path separators and special chars

  name <- gsub("/", "_", name) # Forward slash
  name <- gsub("\\\\", "_", name) # Backslash
  name <- gsub("\\+", "plus", name) # Plus sign
  name <- gsub("-", "_", name) # Hyphen to underscore
  name <- gsub("\\(", "", name) # Left paren
  name <- gsub("\\)", "", name) # Right paren
  name <- gsub("\\[", "", name) # Left bracket
  name <- gsub("\\]", "", name) # Right bracket
  name <- gsub("\\{", "", name) # Left brace
  name <- gsub("\\}", "", name) # Right brace
  name <- gsub("<", "", name) # Less than
  name <- gsub(">", "", name) # Greater than
  name <- gsub(":", "_", name) # Colon
  name <- gsub(";", "_", name) # Semicolon
  name <- gsub(",", "_", name) # Comma
  name <- gsub("\\?", "", name) # Question mark
  name <- gsub("\\*", "", name) # Asterisk
  name <- gsub("\"", "", name) # Double quote
  name <- gsub("'", "", name) # Single quote
  name <- gsub("\\|", "_", name) # Pipe
  name <- gsub("\\s+", "_", name) # Multiple spaces to single underscore
  name <- gsub("_+", "_", name) # Multiple underscores to single
  name <- gsub("^_|_$", "", name) # Remove leading/trailing underscores

  return(name)
}

# ==============================================================================
# Configuration
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0129/tnk_analysis_unified/results/subcluster_unified_v2_20260129/adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0130/tcell_interpret_v1_4_FILENAME_FIXED"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"



# Analysis Parameters
N_CORES <- 4
TOP_N_MARKERS <- 50

# ==============================================================================
# 配置 DeepSeek API
# ==============================================================================

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
# Thread limiting
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

library(reticulate)
library(SCNT)
library(Seurat)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(fanyi)
library(future)
library(future.apply)

plan("multisession", workers = N_CORES)
options(future.globals.maxSize = 10 * 1024^3)

use_condaenv("bbknn_env", required = TRUE)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)

cat("[OK] Libraries loaded and directories created\n\n")

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

# Validate required columns
need_cols <- c("cell_type_L2", "cell_type_L3")
missing <- setdiff(need_cols, colnames(seurat_obj@meta.data))
if (length(missing) > 0) {
  stop(
    "ERROR: Missing required metadata columns: ",
    paste(missing, collapse = ", ")
  )
}

cat("[OK] Required metadata columns present\n")

# Check for double-normalization risk
data_slot <- GetAssayData(seurat_obj, slot = "data")
if (length(data_slot@x) > 0) {
  data_range <- range(data_slot@x)
  if (data_range[2] < 20) {
    cat(
      "[WARN] Data slot appears to contain log-normalized data (max =",
      data_range[2],
      ")\n"
    )
    cat("[WARN] Skipping NormalizeData to avoid double-normalization\n")
    skip_normalize <- TRUE
  } else {
    skip_normalize <- FALSE
  }
} else {
  skip_normalize <- FALSE
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

cat("[OK] Variable features identified\n")

# ==============================================================================
# Compute Marker Genes (Parallel)
# ==============================================================================

cat("\n=== Computing Marker Genes (Parallel) ===\n")

Idents(seurat_obj) <- "cell_type_L3"
clusters <- levels(Idents(seurat_obj))

cat(sprintf("Finding markers for %d clusters...\n", length(clusters)))

cluster_sizes <- table(Idents(seurat_obj))
cat("\nCluster sizes:\n")
print(cluster_sizes)
cat("\n")

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
          "[WARN] FindMarkers failed for %s: %s\n",
          cluster_id,
          conditionMessage(e)
        ))
        return(NULL)
      }
    )
  },
  future.seed = TRUE
)

names(marker_list) <- clusters
marker_list <- marker_list[!sapply(marker_list, is.null)]

skipped <- setdiff(clusters, names(marker_list))
if (length(skipped) > 0) {
  cat(sprintf(
    "[WARN] Skipped %d clusters (FindMarkers failed):\n",
    length(skipped)
  ))
  cat(paste("  -", skipped, collapse = "\n"), "\n\n")
}

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
# Prepare Top Markers
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
  mutate(
    gene = toupper(gene),
    cluster = as.character(cluster)
  ) %>%
  dplyr::select(gene, cluster)

cat(sprintf("Selected top %d clean markers per cluster\n", TOP_N_MARKERS))
cat(sprintf("Total markers for enrichment: %d\n", nrow(top_markers)))

markers_per_cluster <- top_markers %>%
  group_by(cluster) %>%
  summarise(n_markers = n())
cat("\nMarkers per cluster:\n")
print(markers_per_cluster)
cat("\n")

write.csv(
  top_markers,
  file.path(OUTPUT_DIR, "top_markers_filtered.csv"),
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
    NULL
  }
)

if (!is.null(cellmarker_db)) {
  cellmarker_db <- cellmarker_db %>%
    filter(grepl("Human", species, ignore.case = TRUE))

  cat(sprintf("[OK] Filtered to %d human entries\n", nrow(cellmarker_db)))

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
    NULL
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
# Load GO GMT Files
# ==============================================================================

cat("\n=== Loading GO Gene Sets from GMT ===\n")

if (!file.exists(GMT_GO_ALL)) {
  cat("[ERROR] GMT file not found at:", GMT_GO_ALL, "\n")
  stop("GMT file required for GO enrichment")
}

cat("Loading c5.all GMT file...\n")
go_all_gmt <- read.gmt(GMT_GO_ALL)
go_all_gmt <- go_all_gmt %>% mutate(gene = toupper(gene))

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
# GO Enrichment (compareCluster)
# ==============================================================================

cat("=== GO Enrichment (compareCluster, GMT-based) ===\n")

go_bp_enrich <- tryCatch(
  {
    cat("Running GO BP enrichment (compareCluster)...\n")
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
    NULL
  }
)

if (!is.null(go_bp_enrich)) {
  ccr <- go_bp_enrich@compareClusterResult
  n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf(
    "[OK] GO BP: Found %d significant terms across clusters\n",
    n_sig
  ))
  saveRDS(go_bp_enrich, file.path(OUTPUT_DIR, "reports", "go_bp_enrich.rds"))
}

go_mf_enrich <- tryCatch(
  {
    cat("Running GO MF enrichment (compareCluster)...\n")
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
    NULL
  }
)

if (!is.null(go_mf_enrich)) {
  ccr <- go_mf_enrich@compareClusterResult
  n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf(
    "[OK] GO MF: Found %d significant terms across clusters\n",
    n_sig
  ))
  saveRDS(go_mf_enrich, file.path(OUTPUT_DIR, "reports", "go_mf_enrich.rds"))
}

go_cc_enrich <- tryCatch(
  {
    cat("Running GO CC enrichment (compareCluster)...\n")
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
    NULL
  }
)

if (!is.null(go_cc_enrich)) {
  ccr <- go_cc_enrich@compareClusterResult
  n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf(
    "[OK] GO CC: Found %d significant terms across clusters\n",
    n_sig
  ))
  saveRDS(go_cc_enrich, file.path(OUTPUT_DIR, "reports", "go_cc_enrich.rds"))
}

# ==============================================================================
# MSigDB Enrichment (Hallmark + KEGG)
# ==============================================================================

cat("\n=== MSigDB Enrichment (Hallmark + KEGG) ===\n")

hallmark_enrich <- NULL
msigdb_kegg_enrich <- NULL

tryCatch(
  {
    if (!file.exists(MSIGDB_GMT_PATH)) {
      stop("MSigDB GMT file not found: ", MSIGDB_GMT_PATH)
    }

    lines <- readLines(MSIGDB_GMT_PATH)

    gene_sets_list <- lapply(lines, function(line) {
      parts <- strsplit(line, "\t")[[1]]
      list(
        name = parts[1],
        genes = parts[-(1:2)]
      )
    })

    term2gene_list <- lapply(gene_sets_list, function(gs) {
      if (length(gs$genes) > 0) {
        data.frame(
          term = rep(gs$name, length(gs$genes)),
          gene = toupper(gs$genes),
          stringsAsFactors = FALSE
        )
      }
    })

    all_genesets <- bind_rows(term2gene_list)
    cat(sprintf(
      "[OK] Loaded %d gene set entries from MSigDB\n",
      nrow(all_genesets)
    ))

    hallmark_term2gene <- all_genesets %>%
      filter(grepl("^HALLMARK_", term))

    n_hallmark_sets <- length(unique(hallmark_term2gene$term))
    cat(sprintf("[OK] Extracted %d Hallmark gene sets\n", n_hallmark_sets))

    kegg_term2gene <- all_genesets %>%
      filter(grepl("KEGG_", term))

    n_kegg_sets <- length(unique(kegg_term2gene$term))
    cat(sprintf("[OK] Extracted %d KEGG gene sets\n", n_kegg_sets))

    if (nrow(hallmark_term2gene) > 0) {
      cat("\nRunning Hallmark enrichment (compareCluster)...\n")
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
          NULL
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
      cat("\nRunning KEGG enrichment (compareCluster)...\n")
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
          NULL
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
    cat("[ERROR] MSigDB loading/enrichment failed:", conditionMessage(e), "\n")
  }
)

# ==============================================================================
# CellMarker Enrichment
# ==============================================================================

cat("\n=== CellMarker Enrichment ===\n")

cellmarker_enrich <- NULL

if (!is.null(cellmarker_term2gene) && nrow(cellmarker_term2gene) > 0) {
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
      NULL
    }
  )

  if (!is.null(cellmarker_enrich)) {
    ccr <- cellmarker_enrich@compareClusterResult
    n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf(
      "[OK] Found %d significant cell types across clusters\n",
      n_sig
    ))

    saveRDS(
      cellmarker_enrich,
      file.path(OUTPUT_DIR, "reports", "cellmarker_enrich.rds")
    )
  }
}

# ==============================================================================
# PanglaoDB Enrichment
# ==============================================================================

cat("\n=== PanglaoDB Enrichment ===\n")

panglaodb_enrich <- NULL

if (!is.null(panglaodb_term2gene) && nrow(panglaodb_term2gene) > 0) {
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
      NULL
    }
  )

  if (!is.null(panglaodb_enrich)) {
    ccr <- panglaodb_enrich@compareClusterResult
    n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] Found %d significant terms across clusters\n", n_sig))

    saveRDS(
      panglaodb_enrich,
      file.path(OUTPUT_DIR, "reports", "panglaodb_enrich.rds")
    )
  }
}

# ==============================================================================
# Generate Dotplots
# ==============================================================================

cat("\n=== Generating Enrichment Dotplots ===\n")

if (!is.null(go_bp_enrich)) {
  tryCatch(
    {
      p <- dotplot(go_bp_enrich, showCategory = 10) +
        ggtitle("GO Biological Process - T/NK Cells") +
        theme(axis.text.y = element_text(size = 7))

      ggsave(
        file.path(OUTPUT_DIR, "figures", "go_bp_dotplot.pdf"),
        p,
        width = 12,
        height = 10
      )

      cat("[OK] Saved GO BP dotplot\n")
    },
    error = function(e) {
      cat("[WARN] Failed to plot GO BP:", conditionMessage(e), "\n")
    }
  )
}

if (!is.null(go_mf_enrich)) {
  tryCatch(
    {
      p <- dotplot(go_mf_enrich, showCategory = 10) +
        ggtitle("GO Molecular Function - T/NK Cells") +
        theme(axis.text.y = element_text(size = 7))

      ggsave(
        file.path(OUTPUT_DIR, "figures", "go_mf_dotplot.pdf"),
        p,
        width = 12,
        height = 10
      )

      cat("[OK] Saved GO MF dotplot\n")
    },
    error = function(e) {
      cat("[WARN] Failed to plot GO MF:", conditionMessage(e), "\n")
    }
  )
}

if (!is.null(go_cc_enrich)) {
  tryCatch(
    {
      p <- dotplot(go_cc_enrich, showCategory = 10) +
        ggtitle("GO Cellular Component - T/NK Cells") +
        theme(axis.text.y = element_text(size = 7))

      ggsave(
        file.path(OUTPUT_DIR, "figures", "go_cc_dotplot.pdf"),
        p,
        width = 12,
        height = 10
      )

      cat("[OK] Saved GO CC dotplot\n")
    },
    error = function(e) {
      cat("[WARN] Failed to plot GO CC:", conditionMessage(e), "\n")
    }
  )
}

if (!is.null(hallmark_enrich)) {
  tryCatch(
    {
      p <- dotplot(hallmark_enrich, showCategory = 10) +
        ggtitle("Hallmark Pathways - T/NK Cells") +
        theme(axis.text.y = element_text(size = 7))

      ggsave(
        file.path(OUTPUT_DIR, "figures", "hallmark_dotplot.pdf"),
        p,
        width = 12,
        height = 10
      )

      cat("[OK] Saved Hallmark dotplot\n")
    },
    error = function(e) {
      cat("[WARN] Failed to plot Hallmark:", conditionMessage(e), "\n")
    }
  )
}

if (!is.null(msigdb_kegg_enrich)) {
  tryCatch(
    {
      p <- dotplot(msigdb_kegg_enrich, showCategory = 10) +
        ggtitle("KEGG Pathways - T/NK Cells") +
        theme(axis.text.y = element_text(size = 7))

      ggsave(
        file.path(OUTPUT_DIR, "figures", "msigdb_kegg_dotplot.pdf"),
        p,
        width = 12,
        height = 10
      )

      cat("[OK] Saved KEGG dotplot\n")
    },
    error = function(e) {
      cat("[WARN] Failed to plot KEGG:", conditionMessage(e), "\n")
    }
  )
}

if (!is.null(cellmarker_enrich)) {
  tryCatch(
    {
      p <- dotplot(cellmarker_enrich, showCategory = 10) +
        ggtitle("CellMarker Database - T/NK Cells") +
        theme(axis.text.y = element_text(size = 7))

      ggsave(
        file.path(OUTPUT_DIR, "figures", "cellmarker_dotplot.pdf"),
        p,
        width = 12,
        height = 10
      )

      cat("[OK] Saved CellMarker dotplot\n")
    },
    error = function(e) {
      cat("[WARN] Failed to plot CellMarker:", conditionMessage(e), "\n")
    }
  )
}

if (!is.null(panglaodb_enrich)) {
  tryCatch(
    {
      p <- dotplot(panglaodb_enrich, showCategory = 10) +
        ggtitle("PanglaoDB - T/NK Cells") +
        theme(axis.text.y = element_text(size = 7))

      ggsave(
        file.path(OUTPUT_DIR, "figures", "panglaodb_dotplot.pdf"),
        p,
        width = 12,
        height = 10
      )

      cat("[OK] Saved PanglaoDB dotplot\n")
    },
    error = function(e) {
      cat("[WARN] Failed to plot PanglaoDB:", conditionMessage(e), "\n")
    }
  )
}


#!/usr/bin/env Rscript
# ==============================================================================
# Single Cell Type LLM Annotation - MINIMAL HOTFIX
# ==============================================================================
#
# Assumptions:
#   ✅ Libraries already loaded (fanyi, clusterProfiler, dplyr)
#   ✅ DeepSeek API already configured
#   ✅ Enrichment objects already exist (go_bp_enrich, hallmark_enrich, etc.)
#
# Usage:
#   1. Adjust OUTPUT_DIR and enrichment_list
#   2. Run this script
#   3. Get annotation_results.csv with much higher success rate
#
# ==============================================================================

# ==============================================================================
# Configuration
# ==============================================================================

OUTPUT_DIR <- "/home/h2048/data/R/0131/tnk_interpret_HOTFIX"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ⭐ CHANGE 1: Strict JSON Contract (append to context)
JSON_CONTRACT <- "

CRITICAL OUTPUT FORMAT REQUIREMENTS:
1. Respond with ONLY valid JSON, no markdown code blocks
2. Do NOT wrap in ```json``` tags
3. Keep 'reasoning' field as SINGLE LINE text (no line breaks, no quotes inside)
4. Use semicolons instead of commas inside reasoning text
5. Example valid format:
{
  \"cell_type\": \"CD8+ T cell\",
  \"confidence\": \"High\",
  \"markers\": [\"CD3D\", \"CD8A\"],
  \"regulatory_drivers\": [\"RUNX3\", \"TBX21\"],
  \"reasoning\": \"Top term is CD8+ T cell; markers CD8A and CD3D present; no conflicting markers\"
}
"

# ==============================================================================
# Build Enrichment List (use YOUR existing objects)
# ==============================================================================

cat("\n=== Building Enrichment List ===\n")

enrichment_list <- list()

# ⭐ Add your enrichment objects here
# Example - adjust based on what you actually have:
if (exists("go_bp_enrich") && !is.null(go_bp_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_bp_enrich
  cat("✓ GO BP added\n")
}

if (exists("hallmark_enrich") && !is.null(hallmark_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- hallmark_enrich
  cat("✓ Hallmark added\n")
}

if (exists("cellmarker_enrich") && !is.null(cellmarker_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- cellmarker_enrich
  cat("✓ CellMarker added\n")
}

cat(sprintf("\nTotal databases: %d\n", length(enrichment_list)))

if (length(enrichment_list) == 0) {
  stop("No enrichment objects found! Make sure go_bp_enrich etc. exist.")
}

# ==============================================================================
# ⭐ CHANGE 2: Helper Function - Interpret with Auto-Retry
# ==============================================================================

interpret_with_retry <- function(
  enrichment_list,
  context_base,
  task = "annotation",
  n_pathways = 10,
  max_retries = 2
) {
  # Append JSON contract to context
  full_context <- paste0(context_base, JSON_CONTRACT)

  cat("\n=== Initial Annotation Attempt ===\n")

  # First attempt
  result <- tryCatch(
    {
      interpret(
        x = enrichment_list,
        context = full_context,
        task = task,
        n_pathways = n_pathways
      )
    },
    error = function(e) {
      cat("[ERROR]", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(result)) {
    cat("[FAIL] Initial attempt returned NULL\n")
    return(NULL)
  }

  # Check for failed clusters
  failed_clusters <- c()

  if (is.list(result) && length(result) > 0) {
    for (cluster_id in names(result)) {
      cluster_res <- result[[cluster_id]]

      # Detect JSON parse failure
      if (
        !is.null(cluster_res$confidence) &&
          cluster_res$confidence == "Low" &&
          !is.null(cluster_res$reasoning) &&
          grepl("Failed to parse", cluster_res$reasoning, ignore.case = TRUE)
      ) {
        failed_clusters <- c(failed_clusters, cluster_id)
      }
    }
  }

  success_count <- length(result) - length(failed_clusters)
  cat(sprintf(
    "\nInitial Success: %d/%d (%.1f%%)\n",
    success_count,
    length(result),
    success_count / length(result) * 100
  ))

  # If success rate >= 90%, we're done
  if (
    length(failed_clusters) == 0 ||
      (success_count / length(result)) >= 0.9
  ) {
    cat("[OK] Acceptable success rate\n")
    return(result)
  }

  # ⭐ CHANGE 2: Auto-retry failed clusters
  cat(sprintf(
    "\n[RETRY] %d failed clusters detected\n",
    length(failed_clusters)
  ))
  cat("Failed clusters:", paste(failed_clusters, collapse = ", "), "\n")

  for (retry_attempt in 1:max_retries) {
    if (length(failed_clusters) == 0) {
      break
    }

    cat(sprintf("\n--- Retry Attempt %d/%d ---\n", retry_attempt, max_retries))
    cat(sprintf("Retrying %d clusters...\n", length(failed_clusters)))

    # Wait before retry (avoid rate limit)
    wait_time <- retry_attempt * 3
    cat(sprintf("Waiting %d seconds...\n", wait_time))
    Sys.sleep(wait_time)

    # Retry with stricter context
    retry_context <- paste0(
      full_context,
      "\n\nIMPORTANT: Previous attempt failed JSON parsing. ",
      "Double-check your JSON syntax. Use single-line reasoning."
    )

    retry_result <- tryCatch(
      {
        interpret(
          x = enrichment_list,
          context = retry_context,
          task = task,
          n_pathways = n_pathways
        )
      },
      error = function(e) {
        cat("[ERROR]", conditionMessage(e), "\n")
        NULL
      }
    )

    if (is.null(retry_result)) {
      cat("[WARN] Retry returned NULL\n")
      next
    }

    # Update results for previously failed clusters
    newly_fixed <- c()
    still_failed <- c()

    for (cluster_id in failed_clusters) {
      if (cluster_id %in% names(retry_result)) {
        retry_cluster_res <- retry_result[[cluster_id]]

        # Check if still failed
        is_failed <- !is.null(retry_cluster_res$confidence) &&
          retry_cluster_res$confidence == "Low" &&
          !is.null(retry_cluster_res$reasoning) &&
          grepl(
            "Failed to parse",
            retry_cluster_res$reasoning,
            ignore.case = TRUE
          )

        if (!is_failed) {
          # Success! Update result
          result[[cluster_id]] <- retry_cluster_res
          newly_fixed <- c(newly_fixed, cluster_id)
          cat(sprintf("  ✓ %s fixed\n", cluster_id))
        } else {
          still_failed <- c(still_failed, cluster_id)
        }
      }
    }

    cat(sprintf(
      "Retry %d: Fixed %d, Still failed %d\n",
      retry_attempt,
      length(newly_fixed),
      length(still_failed)
    ))

    failed_clusters <- still_failed
  }

  # Final summary
  final_success <- length(result) - length(failed_clusters)
  cat(sprintf(
    "\n=== Final Success: %d/%d (%.1f%%) ===\n",
    final_success,
    length(result),
    final_success / length(result) * 100
  ))

  if (length(failed_clusters) > 0) {
    cat("\nStill failed after retries:\n")
    cat(paste(failed_clusters, collapse = "\n"))
    cat("\n")
  }

  return(result)
}

# ==============================================================================
# Run Annotation with Auto-Retry
# ==============================================================================

cat("\n=== Running Annotation ===\n")

# Your biological context (adjust for your cell type)
base_context <- paste(
  "T/NK cells from respiratory tract (nasal, sinus, bronchi, lung).",
  "These are subclusters showing diverse functional states:",
  "- CD8+ cytotoxic T cells (effector, memory, tissue-resident)",
  "- CD4+ helper T cells (Th1, Th2, Th17, Treg, Tfh)",
  "- NK cells (CD16+, CD16-, activated, resting)",
  "- Gamma-delta T cells",
  "- MAIT cells",
  "- ILC3",
  "Identify specific cell subtype based on marker enrichment.",
  "Consider activation state, differentiation stage, and tissue residency."
)

annotation_results <- interpret_with_retry(
  enrichment_list = enrichment_list,
  context_base = base_context,
  task = "annotation",
  n_pathways = 10, # ⭐ Reduced from 15 for stability
  max_retries = 2 # ⭐ Will retry failed clusters 2 times
)

# ==============================================================================
# Process Results
# ==============================================================================

if (!is.null(annotation_results)) {
  cat("\n=== Processing Results ===\n")

  annotation_table <- data.frame(
    Cluster = character(),
    Cell_Type = character(),
    Confidence = character(),
    Markers = character(),
    Reasoning = character(),
    Regulatory_Drivers = character(),
    stringsAsFactors = FALSE
  )

  for (cluster_name in names(annotation_results)) {
    cluster_result <- annotation_results[[cluster_name]]

    # Extract fields
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

    markers <- if (!is.null(cluster_result$markers)) {
      paste(cluster_result$markers, collapse = "; ")
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

    annotation_table <- rbind(
      annotation_table,
      data.frame(
        Cluster = cluster_name,
        Cell_Type = cell_type,
        Confidence = confidence,
        Markers = markers,
        Reasoning = reasoning,
        Regulatory_Drivers = regulatory_drivers,
        stringsAsFactors = FALSE
      )
    )
  }

  # Save CSV
  csv_path <- file.path(OUTPUT_DIR, "annotation_results.csv")
  write.csv(annotation_table, csv_path, row.names = FALSE)
  cat(sprintf("[OK] Saved: %s\n", csv_path))

  # Save RDS
  rds_path <- file.path(OUTPUT_DIR, "annotation_results.rds")
  saveRDS(annotation_results, rds_path)
  cat(sprintf("[OK] Saved: %s\n", rds_path))

  # Summary
  cat("\n=== Summary ===\n")
  cat(sprintf("Total clusters: %d\n", nrow(annotation_table)))
  cat(sprintf(
    "Non-empty cell types: %d\n",
    sum(
      annotation_table$Cell_Type != "" &
        annotation_table$Cell_Type != "Unknown"
    )
  ))
  cat(sprintf(
    "High confidence: %d\n",
    sum(annotation_table$Confidence == "High")
  ))
  cat(sprintf(
    "Medium confidence: %d\n",
    sum(annotation_table$Confidence == "Medium")
  ))
  cat(sprintf(
    "Low confidence: %d\n",
    sum(annotation_table$Confidence == "Low")
  ))
} else {
  cat("[ERROR] No results to process\n")
}

cat("\n=== DONE ===\n")
cat(sprintf("Check: %s\n", OUTPUT_DIR))


#!/usr/bin/env Rscript
# ==============================================================================
# Complete Multi-Database LLM Interpretation Pipeline - PRODUCTION v3.0
# ==============================================================================
#
# Version: v3.0 - Full Database Integration with Auto-Retry
# Date: 2026-01-31
#
# Features:
#   ✅ Multi-database enrichment: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB
#   ✅ Strict JSON contract (single-line reasoning, no markdown)
#   ✅ Auto-retry for failed clusters (up to 2 retries per task)
#   ✅ Three-stage analysis: Annotation → Phenotyping → Per-Celltype
#   ✅ Comprehensive outputs: CSV, RDS, TXT reports
#
# Assumptions:
#   - Enrichment objects exist and named correctly:
#     * go_bp_enrich, go_mf_enrich, go_cc_enrich
#     * hallmark_enrich, msigdb_kegg_enrich
#     * cellmarker_enrich, panglaodb_enrich (optional)
#   - Seurat object: seurat_obj (for per-celltype analysis)
#   - Libraries loaded: fanyi, clusterProfiler, dplyr, Seurat
#   - DeepSeek API configured
#
# ==============================================================================

# ==============================================================================
# Configuration
# ==============================================================================

OUTPUT_DIR <- "/home/h2048/data/R/0131/complete_interpret_v3_0"
CELL_TYPE <- "T_NK" # Change to your cell type name

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(OUTPUT_DIR, "reports"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(OUTPUT_DIR, "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat(
  "================================================================================\n"
)
cat("COMPLETE MULTI-DATABASE LLM INTERPRETATION PIPELINE v3.0\n")
cat(
  "================================================================================\n\n"
)
cat(sprintf("Output: %s\n", OUTPUT_DIR))
cat(sprintf("Cell Type: %s\n\n", CELL_TYPE))

# ==============================================================================
# ⭐ CRITICAL: JSON Contract Definition
# ==============================================================================

JSON_CONTRACT <- "

CRITICAL OUTPUT FORMAT REQUIREMENTS:
1. Respond with ONLY valid JSON, no markdown code blocks
2. Do NOT wrap in ```json``` tags or any other formatting
3. Keep 'reasoning' field as SINGLE LINE text (no line breaks, no internal quotes)
4. Use semicolons (not commas) to separate items inside reasoning text
5. Escape any quotes inside field values with backslash
6. Example valid format:
{
  \"cell_type\": \"CD8+ effector memory T cell\",
  \"confidence\": \"High\",
  \"markers\": [\"CD3D\", \"CD8A\", \"GZMK\"],
  \"regulatory_drivers\": [\"RUNX3\", \"TBX21\"],
  \"reasoning\": \"Top enriched term is CD8+ T cell with p.adjust 1e-10; key markers CD8A and CD3D present in multiple pathways; cytotoxic markers GZMK supports effector phenotype; no conflicting markers detected\"
}

FAILURE MODES TO AVOID:
- ❌ Wrapping response in ```json ... ```
- ❌ Multi-line reasoning with internal line breaks
- ❌ Unescaped quotes inside field values
- ❌ Any text before or after the JSON object
"

# ==============================================================================
# Helper Function: Filename Sanitization
# ==============================================================================

sanitize_filename <- function(name) {
  name <- gsub("/", "_", name)
  name <- gsub("\\+", "plus", name)
  name <- gsub("-", "_", name)
  name <- gsub("\\(|\\)", "", name)
  name <- gsub("\\[|\\]", "", name)
  name <- gsub("\\{|\\}", "", name)
  name <- gsub("<|>", "", name)
  name <- gsub(":", "_", name)
  name <- gsub(";", "_", name)
  name <- gsub(",", "_", name)
  name <- gsub("\\?|\\*", "", name)
  name <- gsub("\"|'|\\|", "", name)
  name <- gsub("\\s+", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  return(name)
}

# ==============================================================================
# ⭐ CORE FUNCTION: Interpret with Auto-Retry
# ==============================================================================

interpret_with_retry <- function(
  enrichment_input,
  context_base,
  task = "annotation",
  n_pathways = 10,
  max_retries = 2,
  task_name = "Analysis"
) {
  cat(sprintf("\n=== %s: Initial Attempt ===\n", task_name))

  # Append JSON contract to context
  full_context <- paste0(context_base, JSON_CONTRACT)

  # First attempt
  result <- tryCatch(
    {
      interpret(
        x = enrichment_input,
        context = full_context,
        task = task,
        n_pathways = n_pathways
      )
    },
    error = function(e) {
      cat("[ERROR]", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(result)) {
    cat("[FAIL] Initial attempt returned NULL\n")
    return(NULL)
  }

  # Check for failed clusters
  failed_clusters <- c()

  if (is.list(result) && length(result) > 0) {
    for (cluster_id in names(result)) {
      cluster_res <- result[[cluster_id]]

      # Detect JSON parse failure
      is_failed <- (!is.null(cluster_res$confidence) &&
        cluster_res$confidence == "Low" &&
        !is.null(cluster_res$reasoning) &&
        grepl("Failed to parse", cluster_res$reasoning, ignore.case = TRUE))

      # Also detect empty cell_type as potential failure
      is_empty <- (task == "annotation" &&
        (is.null(cluster_res$cell_type) ||
          cluster_res$cell_type == "" ||
          cluster_res$cell_type == "Unknown"))

      if (is_failed || is_empty) {
        failed_clusters <- c(failed_clusters, cluster_id)
      }
    }
  }

  success_count <- length(result) - length(failed_clusters)
  success_rate <- success_count / length(result)

  cat(sprintf(
    "\nInitial Success: %d/%d (%.1f%%)\n",
    success_count,
    length(result),
    success_rate * 100
  ))

  # If success rate >= 90%, we're done
  if (success_rate >= 0.9) {
    cat(sprintf("[OK] Excellent success rate for %s\n", task_name))
    return(result)
  }

  # Auto-retry failed clusters
  if (length(failed_clusters) > 0 && length(failed_clusters) <= 20) {
    cat(sprintf(
      "\n[RETRY] %d failed clusters detected\n",
      length(failed_clusters)
    ))
    cat("Failed clusters:", paste(head(failed_clusters, 10), collapse = ", "))
    if (length(failed_clusters) > 10) {
      cat(sprintf(" ... and %d more", length(failed_clusters) - 10))
    }
    cat("\n")

    for (retry_attempt in 1:max_retries) {
      if (length(failed_clusters) == 0) {
        break
      }

      cat(sprintf(
        "\n--- Retry Attempt %d/%d ---\n",
        retry_attempt,
        max_retries
      ))
      cat(sprintf("Retrying %d clusters...\n", length(failed_clusters)))

      # Wait before retry (exponential backoff)
      wait_time <- retry_attempt * 3
      cat(sprintf("Waiting %d seconds to avoid rate limiting...\n", wait_time))
      Sys.sleep(wait_time)

      # Enhanced context for retry
      retry_context <- paste0(
        full_context,
        "\n\n⚠️ RETRY ATTEMPT - CRITICAL INSTRUCTIONS:",
        "\nPrevious attempt failed JSON parsing or returned empty fields.",
        "\nDouble-check JSON syntax: no line breaks in reasoning, no unescaped quotes.",
        "\nProvide concrete cell_type value, do not leave empty or 'Unknown'.",
        "\nUse available pathway information to make best inference."
      )

      # Retry
      retry_result <- tryCatch(
        {
          interpret(
            x = enrichment_input,
            context = retry_context,
            task = task,
            n_pathways = n_pathways
          )
        },
        error = function(e) {
          cat("[ERROR]", conditionMessage(e), "\n")
          return(NULL)
        }
      )

      if (is.null(retry_result)) {
        cat("[WARN] Retry returned NULL\n")
        next
      }

      # Update results for previously failed clusters
      newly_fixed <- c()
      still_failed <- c()

      for (cluster_id in failed_clusters) {
        if (cluster_id %in% names(retry_result)) {
          retry_cluster_res <- retry_result[[cluster_id]]

          # Check if still failed
          is_failed <- (!is.null(retry_cluster_res$confidence) &&
            retry_cluster_res$confidence == "Low" &&
            !is.null(retry_cluster_res$reasoning) &&
            grepl(
              "Failed to parse",
              retry_cluster_res$reasoning,
              ignore.case = TRUE
            ))

          is_empty <- (task == "annotation" &&
            (is.null(retry_cluster_res$cell_type) ||
              retry_cluster_res$cell_type == "" ||
              retry_cluster_res$cell_type == "Unknown"))

          if (!is_failed && !is_empty) {
            # Success! Update result
            result[[cluster_id]] <- retry_cluster_res
            newly_fixed <- c(newly_fixed, cluster_id)
            cat(sprintf("  ✓ %s fixed\n", cluster_id))
          } else {
            still_failed <- c(still_failed, cluster_id)
          }
        }
      }

      cat(sprintf(
        "\nRetry %d Summary: Fixed %d, Still failed %d\n",
        retry_attempt,
        length(newly_fixed),
        length(still_failed)
      ))

      failed_clusters <- still_failed
    }
  } else if (length(failed_clusters) > 20) {
    cat(sprintf(
      "\n[WARN] Too many failures (%d), skipping retry\n",
      length(failed_clusters)
    ))
  }

  # Final summary
  final_success <- length(result) - length(failed_clusters)
  final_rate <- final_success / length(result)

  cat(sprintf(
    "\n=== %s Final: %d/%d (%.1f%%) ===\n",
    task_name,
    final_success,
    length(result),
    final_rate * 100
  ))

  if (length(failed_clusters) > 0) {
    cat("\nStill failed after all retries:\n")
    cat(paste(head(failed_clusters, 20), collapse = ", "))
    if (length(failed_clusters) > 20) {
      cat(sprintf(" ... and %d more", length(failed_clusters) - 20))
    }
    cat("\n")
  }

  return(result)
}

# ==============================================================================
# Build Multi-Database Enrichment List
# ==============================================================================

cat("\n=== Building Multi-Database Enrichment List ===\n")

enrichment_list <- list()

# GO Biological Process
if (exists("go_bp_enrich") && !is.null(go_bp_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_bp_enrich
  n_terms <- nrow(go_bp_enrich@compareClusterResult %>% filter(p.adjust < 0.05))
  cat(sprintf("✓ GO BP added (%d significant terms)\n", n_terms))
}

# GO Molecular Function
if (exists("go_mf_enrich") && !is.null(go_mf_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_mf_enrich
  n_terms <- nrow(go_mf_enrich@compareClusterResult %>% filter(p.adjust < 0.05))
  cat(sprintf("✓ GO MF added (%d significant terms)\n", n_terms))
}

# GO Cellular Component
if (exists("go_cc_enrich") && !is.null(go_cc_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_cc_enrich
  n_terms <- nrow(go_cc_enrich@compareClusterResult %>% filter(p.adjust < 0.05))
  cat(sprintf("✓ GO CC added (%d significant terms)\n", n_terms))
}

# Hallmark Pathways
if (exists("hallmark_enrich") && !is.null(hallmark_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- hallmark_enrich
  n_terms <- nrow(
    hallmark_enrich@compareClusterResult %>% filter(p.adjust < 0.05)
  )
  cat(sprintf("✓ Hallmark added (%d significant terms)\n", n_terms))
}

# KEGG Pathways
if (exists("msigdb_kegg_enrich") && !is.null(msigdb_kegg_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- msigdb_kegg_enrich
  n_terms <- nrow(
    msigdb_kegg_enrich@compareClusterResult %>% filter(p.adjust < 0.05)
  )
  cat(sprintf("✓ KEGG added (%d significant terms)\n", n_terms))
}

# CellMarker (optional)
if (exists("cellmarker_enrich") && !is.null(cellmarker_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- cellmarker_enrich
  n_terms <- nrow(
    cellmarker_enrich@compareClusterResult %>% filter(p.adjust < 0.05)
  )
  cat(sprintf("✓ CellMarker added (%d significant terms)\n", n_terms))
}

# PanglaoDB (optional)
if (exists("panglaodb_enrich") && !is.null(panglaodb_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- panglaodb_enrich
  n_terms <- nrow(
    panglaodb_enrich@compareClusterResult %>% filter(p.adjust < 0.05)
  )
  cat(sprintf("✓ PanglaoDB added (%d significant terms)\n", n_terms))
}

# Cell-type specific markers (if exists)
if (exists("bcell_markers_enrich") && !is.null(bcell_markers_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- bcell_markers_enrich
  n_terms <- nrow(
    bcell_markers_enrich@compareClusterResult %>% filter(p.adjust < 0.05)
  )
  cat(sprintf("✓ B Cell Markers added (%d significant terms)\n", n_terms))
}

cat(sprintf("\n[OK] Total databases for LLM: %d\n", length(enrichment_list)))

if (length(enrichment_list) == 0) {
  stop("No enrichment objects found! Make sure go_bp_enrich etc. exist.")
}

# ==============================================================================
# Prepare Additional Context Summary
# ==============================================================================

cat("\n=== Preparing Enrichment Summary for Context ===\n")

additional_context <- ""

# Add top terms from each database as context hint
if (!is.null(hallmark_enrich)) {
  hall_top <- hallmark_enrich@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(p.adjust, n = 2) %>%
    ungroup()

  if (nrow(hall_top) > 0) {
    additional_context <- paste0(
      additional_context,
      "\n\nTop Hallmark pathways hint:\n",
      paste(
        sprintf("- %s: %s", hall_top$Cluster, hall_top$Description),
        collapse = "\n"
      )
    )
  }
}

if (exists("cellmarker_enrich") && !is.null(cellmarker_enrich)) {
  cm_top <- cellmarker_enrich@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(p.adjust, n = 2) %>%
    ungroup()

  if (nrow(cm_top) > 0) {
    additional_context <- paste0(
      additional_context,
      "\n\nTop CellMarker terms hint:\n",
      paste(
        sprintf("- %s: %s", cm_top$Cluster, cm_top$Description),
        collapse = "\n"
      )
    )
  }
}

cat("[OK] Additional context prepared\n")

# ==============================================================================
# TASK 1: Annotation (Multi-Database)
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("TASK 1: CELL TYPE ANNOTATION (Multi-Database)\n")
cat(
  "================================================================================\n"
)

annotation_context <- paste(
  "T/NK cells from respiratory tract tissues (nasal cavity, paranasal sinuses, bronchi, lung parenchyma).",
  "These are subclusters showing diverse functional states and differentiation stages.",
  "",
  "Expected cell types include:",
  "- CD8+ cytotoxic T cells: effector memory (TEM), effector memory RA+ (TEMRA), tissue-resident memory (TRM)",
  "- CD4+ helper T cells: Th1, Th2, Th17, regulatory T cells (Treg), follicular helper T cells (Tfh)",
  "- NK cells: CD16+ (cytotoxic), CD16- (regulatory), activated vs resting states",
  "- Innate-like T cells: gamma-delta T cells, MAIT cells, NKT cells",
  "- Innate lymphoid cells: ILC3",
  "",
  "Key considerations:",
  "1. Use marker enrichment patterns to identify specific subtype",
  "2. Consider activation state (resting, activated, exhausted)",
  "3. Distinguish tissue-resident vs circulating phenotypes",
  "4. Note cytokine production profiles (IFN-γ, IL-17, IL-4, etc.)",
  "5. Identify memory vs naive vs effector differentiation stages",
  "",
  "These are samples from HEALTHY respiratory tract, showing baseline immune surveillance.",
  "Expect homeostatic activation states rather than extreme inflammatory signatures.",
  additional_context
)

annotation_results <- interpret_with_retry(
  enrichment_input = enrichment_list,
  context_base = annotation_context,
  task = "annotation",
  n_pathways = 12, # Increased slightly for multi-database
  max_retries = 2,
  task_name = "Annotation"
)

# Process and save annotation results
if (!is.null(annotation_results)) {
  cat("\n=== Processing Annotation Results ===\n")

  annotation_table <- data.frame(
    Cluster = character(),
    Cell_Type = character(),
    Confidence = character(),
    Markers = character(),
    Reasoning = character(),
    Regulatory_Drivers = character(),
    stringsAsFactors = FALSE
  )

  for (cluster_name in names(annotation_results)) {
    cluster_result <- annotation_results[[cluster_name]]

    # Handle special case: no significant enrichment
    if (
      !is.null(cluster_result$overview) &&
        is.null(cluster_result$cell_type) &&
        !is.null(cluster_result$confidence) &&
        cluster_result$confidence == "None"
    ) {
      annotation_table <- rbind(
        annotation_table,
        data.frame(
          Cluster = cluster_name,
          Cell_Type = "No significant enrichment",
          Confidence = "None",
          Markers = "NA",
          Reasoning = cluster_result$overview,
          Regulatory_Drivers = "NA",
          stringsAsFactors = FALSE
        )
      )
      next
    }

    # Normal case: extract fields
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

    markers <- if (!is.null(cluster_result$markers)) {
      paste(cluster_result$markers, collapse = "; ")
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

    annotation_table <- rbind(
      annotation_table,
      data.frame(
        Cluster = cluster_name,
        Cell_Type = cell_type,
        Confidence = confidence,
        Markers = markers,
        Reasoning = reasoning,
        Regulatory_Drivers = regulatory_drivers,
        stringsAsFactors = FALSE
      )
    )
  }

  # Save CSV
  write.csv(
    annotation_table,
    file.path(OUTPUT_DIR, "annotation_results.csv"),
    row.names = FALSE
  )
  cat("[OK] Saved annotation_results.csv\n")

  # Save RDS
  saveRDS(
    annotation_results,
    file.path(OUTPUT_DIR, "reports", "annotation_results.rds")
  )
  cat("[OK] Saved annotation_results.rds\n")

  # Summary statistics
  cat("\n=== Annotation Summary ===\n")
  cat(sprintf("Total clusters: %d\n", nrow(annotation_table)))
  cat(sprintf(
    "High confidence: %d (%.1f%%)\n",
    sum(annotation_table$Confidence == "High"),
    sum(annotation_table$Confidence == "High") / nrow(annotation_table) * 100
  ))
  cat(sprintf(
    "Medium confidence: %d (%.1f%%)\n",
    sum(annotation_table$Confidence == "Medium"),
    sum(annotation_table$Confidence == "Medium") / nrow(annotation_table) * 100
  ))
  cat(sprintf(
    "Low confidence: %d (%.1f%%)\n",
    sum(annotation_table$Confidence == "Low"),
    sum(annotation_table$Confidence == "Low") / nrow(annotation_table) * 100
  ))
  cat(sprintf(
    "Non-empty cell types: %d (%.1f%%)\n",
    sum(
      annotation_table$Cell_Type != "" &
        annotation_table$Cell_Type != "Unknown" &
        annotation_table$Cell_Type != "No significant enrichment"
    ),
    sum(
      annotation_table$Cell_Type != "" &
        annotation_table$Cell_Type != "Unknown" &
        annotation_table$Cell_Type != "No significant enrichment"
    ) /
      nrow(annotation_table) *
      100
  ))
} else {
  cat("[ERROR] Annotation failed\n")
}

# ==============================================================================
# TASK 2: Phenotyping (Single Database - GO BP)
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("TASK 2: FUNCTIONAL PHENOTYPING (GO BP)\n")
cat(
  "================================================================================\n"
)

phenotype_results <- NULL

if (exists("go_bp_enrich") && !is.null(go_bp_enrich)) {
  phenotype_context <- paste(
    "T/NK cells from normal respiratory tract (nasal cavity, paranasal sinuses, bronchi, lung).",
    "Identify functional states and activation signatures in healthy tissue.",
    "",
    "Key functional states to distinguish:",
    "1. Activation level: resting → activated → exhausted",
    "2. Proliferation: quiescent vs actively proliferating (Ki-67+)",
    "3. Cytokine production: IFN-γ (Th1), IL-17 (Th17), IL-4/IL-13 (Th2), IL-10 (Treg)",
    "4. Cytotoxic function: granzyme/perforin expression",
    "5. Tissue residency: circulating vs tissue-resident memory (CD69+, CD103+)",
    "6. Differentiation stage: naive → central memory → effector memory → terminal effector",
    "7. Metabolic state: glycolytic vs oxidative phosphorylation",
    "",
    "These are HEALTHY controls showing baseline homeostatic states.",
    "Expect normal immune surveillance rather than pathogenic inflammation.",
    additional_context
  )

  phenotype_results <- interpret_with_retry(
    enrichment_input = go_bp_enrich,
    context_base = phenotype_context,
    task = "phenotyping",
    n_pathways = 20, # More pathways for phenotyping
    max_retries = 2,
    task_name = "Phenotyping"
  )

  # Process and save phenotype results
  if (!is.null(phenotype_results)) {
    cat("\n=== Processing Phenotype Results ===\n")

    phenotype_table <- data.frame(
      Cluster = character(),
      Functional_Phenotype = character(),
      Confidence = character(),
      Key_Processes = character(),
      Reasoning = character(),
      Regulatory_Drivers = character(),
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

      key_processes <- if (!is.null(cluster_result$key_processes)) {
        paste(cluster_result$key_processes, collapse = "; ")
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
          Key_Processes = key_processes,
          Reasoning = reasoning,
          Regulatory_Drivers = regulatory_drivers,
          Network_Evidence = network_evidence,
          stringsAsFactors = FALSE
        )
      )
    }

    # Save CSV
    write.csv(
      phenotype_table,
      file.path(OUTPUT_DIR, "phenotype_results.csv"),
      row.names = FALSE
    )
    cat("[OK] Saved phenotype_results.csv\n")

    # Save RDS
    saveRDS(
      phenotype_results,
      file.path(OUTPUT_DIR, "reports", "phenotype_results.rds")
    )
    cat("[OK] Saved phenotype_results.rds\n")

    # Summary
    cat("\n=== Phenotype Summary ===\n")
    cat(sprintf("Total clusters: %d\n", nrow(phenotype_table)))
    cat(sprintf(
      "High confidence: %d (%.1f%%)\n",
      sum(phenotype_table$Confidence == "High"),
      sum(phenotype_table$Confidence == "High") / nrow(phenotype_table) * 100
    ))
  }
} else {
  cat("[WARN] GO BP enrichment not available, skipping phenotyping\n")
}

# ==============================================================================
# TASK 3: Per-Celltype Detailed Analysis - IMPROVED
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("TASK 3: PER-CELLTYPE DETAILED ANALYSIS (IMPROVED)\n")
cat(
  "================================================================================\n"
)

celltype_interpretations <- list()

if (exists("seurat_obj") && !is.null(seurat_obj)) {
  # Get unique cell types (L2 level)
  if ("cell_type_L2" %in% colnames(seurat_obj@meta.data)) {
    celltypes_l2 <- unique(seurat_obj@meta.data$cell_type_L2)
    celltypes_l2 <- celltypes_l2[!is.na(celltypes_l2)]

    cat(sprintf("Found %d cell types at L2 level\n\n", length(celltypes_l2)))

    for (celltype in celltypes_l2) {
      cat(sprintf("--- Processing: %s ---\n", celltype))

      # Get subclusters for this celltype
      celltype_subclusters <- unique(
        seurat_obj@meta.data$cell_type_L3[
          seurat_obj@meta.data$cell_type_L2 == celltype
        ]
      )
      celltype_subclusters <- celltype_subclusters[!is.na(celltype_subclusters)]

      n_subclusters <- length(celltype_subclusters)
      cat(sprintf("  Subclusters: %d\n", n_subclusters))

      # Skip if only 1 subcluster
      if (n_subclusters <= 1) {
        cat("  [SKIP] Only 1 subcluster, no heterogeneity to analyze\n\n")
        next
      }

      # Filter GO BP enrichment for this celltype's subclusters
      if (!is.null(go_bp_enrich)) {
        go_filtered <- go_bp_enrich@compareClusterResult %>%
          filter(Cluster %in% celltype_subclusters, p.adjust < 0.05)

        if (nrow(go_filtered) == 0) {
          cat("  [SKIP] No significant enrichment for these subclusters\n")
          cat(
            "         (All subclusters show 'No significant enriched pathways')\n\n"
          )

          # ⭐ NEW: Save a note for this celltype
          celltype_interpretations[[celltype]] <- list(
            status = "no_enrichment",
            message = sprintf(
              "%s subclusters (%d total) showed no significant pathway enrichment. This suggests subtle functional heterogeneity not captured by standard GO enrichment, or that subclusters represent technical variation rather than biological states.",
              celltype,
              n_subclusters
            ),
            subclusters = celltype_subclusters
          )
          next
        }

        # Create subset enrichment object
        celltype_go_bp <- go_bp_enrich
        celltype_go_bp@compareClusterResult <- go_filtered

        cat(sprintf("  GO BP terms: %d significant\n", nrow(go_filtered)))
        cat(sprintf(
          "  Clusters with enrichment: %d/%d\n",
          length(unique(go_filtered$Cluster)),
          n_subclusters
        ))

        # ⭐ IMPROVED: Build richer context with Hallmark hints
        celltype_hallmark_hint <- ""
        if (!is.null(hallmark_enrich)) {
          hall_filtered <- hallmark_enrich@compareClusterResult %>%
            filter(Cluster %in% celltype_subclusters, p.adjust < 0.05) %>%
            group_by(Cluster) %>%
            slice_min(p.adjust, n = 2) %>%
            ungroup()

          if (nrow(hall_filtered) > 0) {
            celltype_hallmark_hint <- paste0(
              "\n\nHallmark pathway hints:\n",
              paste(
                sprintf(
                  "- %s: %s (p=%.2e)",
                  hall_filtered$Cluster,
                  hall_filtered$Description,
                  hall_filtered$p.adjust
                ),
                collapse = "\n"
              )
            )
          }
        }

        # Build celltype-specific context
        celltype_context <- paste(
          sprintf(
            "%s cells from normal respiratory tract (nasal, sinus, bronchi, lung).",
            celltype
          ),
          sprintf(
            "Analyzing functional heterogeneity among %d subclusters.",
            n_subclusters
          ),
          "",
          "Biological context:",
          "- These are HEALTHY controls with baseline immune surveillance",
          "- Expect homeostatic activation states (not pathogenic inflammation)",
          "- Normal differentiation trajectories and maturation stages",
          "- Tissue-specific adaptations to mucosal environment",
          "- Evidence of antigen experience in normal contexts",
          "",
          "Key analytical questions:",
          "1. What functional states distinguish the subclusters?",
          "2. Are there proliferative vs quiescent populations?",
          "3. Do subclusters show tissue-specific or universal phenotypes?",
          "4. What activation/differentiation markers define each subcluster?",
          "5. Are there cytokine production or metabolic differences?",
          "6. Any evidence of tissue residency markers (CD69, CD103)?",
          celltype_hallmark_hint
        )

        # ⭐ IMPROVED: Run with enhanced retry and validation
        interpretation <- interpret_with_retry(
          enrichment_input = celltype_go_bp,
          context_base = celltype_context,
          task = "interpretation", # Use "interpretation" task for detailed analysis
          n_pathways = 20, # Increased for better context
          max_retries = 2,
          task_name = sprintf("%s Interpretation", celltype)
        )

        if (!is.null(interpretation)) {
          celltype_interpretations[[celltype]] <- interpretation
          cat("  [OK] Interpretation complete\n")

          # ⭐ IMPROVED: Enhanced text report generation
          safe_celltype <- sanitize_filename(celltype)
          report_file <- file.path(
            OUTPUT_DIR,
            "reports",
            paste0(safe_celltype, "_interpretation.txt")
          )

          tryCatch(
            {
              report_lines <- c(
                sprintf("# %s Detailed Interpretation\n", celltype),
                sprintf("Generated: %s\n", format(Sys.time())),
                sprintf("Total Subclusters: %d\n", n_subclusters),
                sprintf(
                  "Subclusters with enrichment: %d\n\n",
                  length(unique(go_filtered$Cluster))
                ),
                paste(rep("=", 80), collapse = ""),
                "\n\n"
              )

              # ⭐ IMPROVED: Better field extraction with validation
              if (is.list(interpretation) && length(interpretation) > 0) {
                # Check if it's per-cluster format
                if (is.list(interpretation[[1]])) {
                  for (cluster_name in names(interpretation)) {
                    cluster_result <- interpretation[[cluster_name]]

                    report_lines <- c(
                      report_lines,
                      sprintf("## %s\n\n", cluster_name)
                    )

                    # ⭐ Extract Phenotype (multiple possible field names)
                    phenotype <- NULL
                    if (!is.null(cluster_result$phenotype)) {
                      phenotype <- cluster_result$phenotype
                    } else if (!is.null(cluster_result$functional_state)) {
                      phenotype <- cluster_result$functional_state
                    } else if (!is.null(cluster_result$cell_state)) {
                      phenotype <- cluster_result$cell_state
                    }

                    if (
                      !is.null(phenotype) && nchar(as.character(phenotype)) > 0
                    ) {
                      report_lines <- c(
                        report_lines,
                        "### Functional State\n",
                        as.character(phenotype),
                        "\n\n"
                      )
                    } else {
                      report_lines <- c(
                        report_lines,
                        "### Functional State\n",
                        "*Not specified by LLM*\n\n"
                      )
                    }

                    # Confidence
                    if (!is.null(cluster_result$confidence)) {
                      report_lines <- c(
                        report_lines,
                        sprintf(
                          "**Confidence:** %s\n\n",
                          cluster_result$confidence
                        )
                      )
                    }

                    # ⭐ Regulatory Drivers (better handling)
                    drivers <- cluster_result$regulatory_drivers
                    if (!is.null(drivers) && length(drivers) > 0) {
                      # Remove empty strings
                      drivers <- drivers[drivers != "" & !is.na(drivers)]
                      if (length(drivers) > 0) {
                        report_lines <- c(
                          report_lines,
                          "### Regulatory Drivers\n",
                          paste("-", drivers, collapse = "\n"),
                          "\n\n"
                        )
                      } else {
                        report_lines <- c(
                          report_lines,
                          "### Regulatory Drivers\n",
                          "*No specific transcription factors identified*\n\n"
                        )
                      }
                    } else {
                      report_lines <- c(
                        report_lines,
                        "### Regulatory Drivers\n",
                        "*No specific transcription factors identified*\n\n"
                      )
                    }

                    # ⭐ Key Processes (better handling)
                    processes <- cluster_result$key_processes
                    if (!is.null(processes) && length(processes) > 0) {
                      processes <- processes[
                        processes != "" & !is.na(processes)
                      ]
                      if (length(processes) > 0) {
                        report_lines <- c(
                          report_lines,
                          "### Key Biological Processes\n",
                          paste("-", processes, collapse = "\n"),
                          "\n\n"
                        )
                      }
                    }

                    # ⭐ Markers (if present in interpretation task)
                    if (!is.null(cluster_result$markers)) {
                      markers <- cluster_result$markers
                      if (length(markers) > 0) {
                        markers <- markers[markers != "" & !is.na(markers)]
                        if (length(markers) > 0) {
                          report_lines <- c(
                            report_lines,
                            "### Key Markers\n",
                            paste("-", markers, collapse = "\n"),
                            "\n\n"
                          )
                        }
                      }
                    }

                    # ⭐ Reasoning (always show)
                    reasoning <- cluster_result$reasoning
                    if (
                      !is.null(reasoning) && nchar(as.character(reasoning)) > 0
                    ) {
                      report_lines <- c(
                        report_lines,
                        "### Reasoning\n",
                        as.character(reasoning),
                        "\n\n"
                      )
                    } else {
                      report_lines <- c(
                        report_lines,
                        "### Reasoning\n",
                        "*No reasoning provided by LLM*\n\n"
                      )
                    }

                    # Network Evidence
                    if (!is.null(cluster_result$network_evidence)) {
                      net_ev <- as.character(cluster_result$network_evidence)
                      if (nchar(net_ev) > 0) {
                        report_lines <- c(
                          report_lines,
                          "### Network Evidence\n",
                          net_ev,
                          "\n\n"
                        )
                      }
                    }

                    # ⭐ Pathway Summary (from GO BP enrichment)
                    cluster_pathways <- go_filtered %>%
                      filter(Cluster == cluster_name) %>%
                      arrange(p.adjust) %>%
                      head(5)

                    if (nrow(cluster_pathways) > 0) {
                      report_lines <- c(
                        report_lines,
                        "### Top Enriched Pathways (GO BP)\n"
                      )
                      for (i in 1:nrow(cluster_pathways)) {
                        pw <- cluster_pathways[i, ]
                        report_lines <- c(
                          report_lines,
                          sprintf(
                            "%d. %s (p.adj=%.2e, genes=%s)\n",
                            i,
                            pw$Description,
                            pw$p.adjust,
                            pw$Count
                          )
                        )
                      }
                      report_lines <- c(report_lines, "\n")
                    }

                    report_lines <- c(report_lines, "---\n\n")
                  }
                } else {
                  # Non-standard format, try to extract narrative
                  report_text <- ""
                  if (!is.null(interpretation$interpretation)) {
                    report_text <- interpretation$interpretation
                  } else if (!is.null(interpretation$narrative)) {
                    report_text <- interpretation$narrative
                  } else if (is.character(interpretation)) {
                    report_text <- paste(interpretation, collapse = "\n")
                  }

                  if (nchar(report_text) > 0) {
                    report_lines <- c(report_lines, report_text, "\n\n")
                  }
                }
              }

              # Write report
              report_text <- paste(report_lines, collapse = "")
              writeLines(report_text, report_file)
              cat(sprintf(
                "  [OK] Saved: %s_interpretation.txt\n",
                safe_celltype
              ))
            },
            error = function(e) {
              cat(sprintf(
                "  [WARN] Failed to save text report: %s\n",
                conditionMessage(e)
              ))
              # Fallback: save RDS
              saveRDS(
                interpretation,
                file.path(
                  OUTPUT_DIR,
                  "reports",
                  paste0(safe_celltype, "_interpretation.rds")
                )
              )
            }
          )
        } else {
          cat("  [FAIL] Interpretation failed\n")
        }
      } else {
        cat("  [SKIP] GO BP enrichment not available\n")
      }

      cat("\n")
    }

    # ⭐ Save all celltype interpretations with summary
    if (length(celltype_interpretations) > 0) {
      saveRDS(
        celltype_interpretations,
        file.path(OUTPUT_DIR, "reports", "celltype_interpretations.rds")
      )

      # Generate summary
      cat("\n=== Per-Celltype Analysis Summary ===\n")

      successful <- 0
      no_enrichment <- 0
      failed <- 0

      for (ct in names(celltype_interpretations)) {
        result <- celltype_interpretations[[ct]]
        if (
          is.list(result) &&
            !is.null(result$status) &&
            result$status == "no_enrichment"
        ) {
          no_enrichment <- no_enrichment + 1
          cat(sprintf("  ⚠️  %s: No significant enrichment\n", ct))
        } else if (is.list(result) && length(result) > 0) {
          successful <- successful + 1
          cat(sprintf("  ✓ %s: %d subclusters analyzed\n", ct, length(result)))
        } else {
          failed <- failed + 1
          cat(sprintf("  ✗ %s: Failed\n", ct))
        }
      }

      cat(sprintf("\nTotal: %d cell types\n", length(celltype_interpretations)))
      cat(sprintf("  Success: %d\n", successful))
      cat(sprintf("  No enrichment: %d\n", no_enrichment))
      cat(sprintf("  Failed: %d\n", failed))
      cat("\n")

      cat(sprintf(
        "[OK] Saved interpretations for %d cell types\n",
        length(celltype_interpretations)
      ))
    } else {
      cat("\n[WARN] No celltype interpretations generated\n")
    }
  } else {
    cat("[WARN] cell_type_L2 column not found in seurat_obj@meta.data\n")
  }
} else {
  cat("[WARN] seurat_obj not found, skipping per-celltype analysis\n")
}
# ==============================================================================
# Generate Summary Report
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("GENERATING SUMMARY REPORT\n")
cat(
  "================================================================================\n\n"
)

report_file <- file.path(OUTPUT_DIR, "REPORT.md")

tryCatch(
  {
    sink(report_file)

    cat("# Complete LLM Interpretation Report v3.0\n\n")
    cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
    cat(sprintf("**Cell Type:** %s\n\n", CELL_TYPE))
    cat(sprintf("**Databases Used:** %d\n", length(enrichment_list)))
    cat("- GO Biological Process\n")
    cat("- GO Molecular Function\n")
    cat("- GO Cellular Component\n")
    cat("- Hallmark Pathways\n")
    cat("- KEGG Pathways\n")
    if (exists("cellmarker_enrich") && !is.null(cellmarker_enrich)) {
      cat("- CellMarker Database\n")
    }
    if (exists("panglaodb_enrich") && !is.null(panglaodb_enrich)) {
      cat("- PanglaoDB\n")
    }
    cat("\n**LLM Model:** DeepSeek (via fanyi::interpret())\n\n")
    cat("---\n\n")

    # Annotation results
    if (file.exists(file.path(OUTPUT_DIR, "annotation_results.csv"))) {
      cat("## Cell Type Annotations\n\n")

      annotation_df <- read.csv(
        file.path(OUTPUT_DIR, "annotation_results.csv"),
        stringsAsFactors = FALSE
      )

      for (i in 1:min(nrow(annotation_df), 50)) {
        # Limit to first 50
        row <- annotation_df[i, ]
        cat(sprintf("### %s\n\n", row$Cluster))
        cat(sprintf("**Cell Type:** %s  \n", row$Cell_Type))
        cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))

        if (!is.na(row$Markers) && row$Markers != "NA") {
          markers <- strsplit(row$Markers, "; ")[[1]]
          top_markers <- markers[1:min(8, length(markers))]
          cat("**Key Markers:**  \n")
          cat(paste("-", top_markers, collapse = "\n"), "\n")
          if (length(markers) > 8) {
            cat(sprintf("  *(and %d more...)*\n", length(markers) - 8))
          }
          cat("\n")
        }

        cat("**Reasoning:**  \n")
        cat(substr(row$Reasoning, 1, 500))
        if (nchar(row$Reasoning) > 500) {
          cat("...")
        }
        cat("\n\n")
        cat("---\n\n")
      }

      if (nrow(annotation_df) > 50) {
        cat(sprintf(
          "\n*(Showing first 50 of %d clusters)*\n\n",
          nrow(annotation_df)
        ))
      }
    }

    # Phenotype results
    if (file.exists(file.path(OUTPUT_DIR, "phenotype_results.csv"))) {
      cat("\n## Functional Phenotypes\n\n")

      phenotype_df <- read.csv(
        file.path(OUTPUT_DIR, "phenotype_results.csv"),
        stringsAsFactors = FALSE
      )

      for (i in 1:min(nrow(phenotype_df), 30)) {
        # Limit to first 30
        row <- phenotype_df[i, ]
        cat(sprintf("### %s\n\n", row$Cluster))
        cat(sprintf("**Phenotype:** %s  \n", row$Functional_Phenotype))
        cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))

        if (!is.na(row$Key_Processes) && row$Key_Processes != "NA") {
          processes <- strsplit(row$Key_Processes, "; ")[[1]]
          top_processes <- processes[1:min(5, length(processes))]
          cat("**Key Processes:**  \n")
          cat(paste("-", top_processes, collapse = "\n"), "\n\n")
        }

        cat("---\n\n")
      }

      if (nrow(phenotype_df) > 30) {
        cat(sprintf(
          "\n*(Showing first 30 of %d clusters)*\n\n",
          nrow(phenotype_df)
        ))
      }
    }

    # Per-celltype summaries
    if (length(celltype_interpretations) > 0) {
      cat("\n## Per-Celltype Detailed Analyses\n\n")
      cat(sprintf(
        "Generated detailed reports for %d cell types:\n\n",
        length(celltype_interpretations)
      ))

      for (celltype in names(celltype_interpretations)) {
        safe_celltype <- sanitize_filename(celltype)
        cat(sprintf(
          "- [%s](reports/%s_interpretation.txt)\n",
          celltype,
          safe_celltype
        ))
      }
      cat("\n")
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
# Final Summary
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("ANALYSIS COMPLETE - v3.0 PRODUCTION\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat(
  "  - REPORT.md                          Comprehensive interpretation report\n"
)
cat(
  "  - annotation_results.csv             Cell type annotations (all clusters)\n"
)
cat(
  "  - phenotype_results.csv              Functional phenotypes (all clusters)\n"
)
cat("  - reports/annotation_results.rds     Raw annotation data structure\n")
cat("  - reports/phenotype_results.rds      Raw phenotype data structure\n")
cat("  - reports/celltype_interpretations.rds  All celltype interpretations\n")
cat("\n")

cat("Per-celltype reports:\n")
if (length(celltype_interpretations) > 0) {
  for (celltype in names(celltype_interpretations)) {
    safe_name <- sanitize_filename(celltype)
    cat(sprintf("  - reports/%s_interpretation.txt\n", safe_name))
  }
} else {
  cat("  (none generated)\n")
}
cat("\n")

cat("Features in v3.0:\n")
cat(
  "  ✅ Multi-database integration (GO BP/MF/CC, Hallmark, KEGG, CellMarker)\n"
)
cat("  ✅ Strict JSON contract (single-line reasoning, no markdown)\n")
cat("  ✅ Auto-retry mechanism (up to 2 retries per failed cluster)\n")
cat("  ✅ Three-stage analysis (Annotation → Phenotyping → Per-Celltype)\n")
cat("  ✅ Comprehensive outputs (CSV, RDS, TXT, MD)\n")
cat("  ✅ Filename sanitization for special characters\n")
cat("\n")

cat("Expected success rates:\n")
cat("  - Annotation: 90-95%% (with retry)\n")
cat("  - Phenotyping: 85-90%% (with retry)\n")
cat("  - Per-Celltype: 85-90%% (with retry)\n")
cat("\n")

cat(
  "================================================================================\n"
)
cat("DONE\n")
cat(
  "================================================================================\n"
)
