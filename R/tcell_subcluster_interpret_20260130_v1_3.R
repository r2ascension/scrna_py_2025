#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Cell Subcluster LLM Interpretation - CORRECT VERSION v1.3
# ==============================================================================
#
# Version: v1.3 (2026-01-30)
# Status: Production-ready with CORRECT interpret() usage
#
# Key Fix in v1.3:
#   ✅ P0-6: ADDED BACK context parameter (IT IS SUPPORTED!)
#            - interpret() DOES accept 'context' parameter
#            - interpret() DOES accept 'task' parameter
#            - v1.2 was WRONG to remove these
#
# All v1.1 features retained:
#   ✅ compareCluster-based enrichment (per-cluster, not global pooled)
#   ✅ GMT-based GO enrichment (NO gene ID conversion loss!)
#   ✅ Multi-database support: GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB
#   ✅ Annotation + Phenotype dual tasks
#
# CRITICAL NOTE:
#   fanyi::interpret() official parameters:
#   - x (enrichment object) ✅
#   - context (biological context) ✅ SUPPORTED!
#   - n_pathways (number of pathways) ✅
#   - model (LLM model) ✅
#   - api_key (API key) ✅
#   - task (task type) ✅ SUPPORTED!
#   - prior, add_ppi, gene_fold_change ✅
#
#   Unsupported (v1.1 errors):
#   - save_history ❌
#   - history_dir ❌
#   - additional_context ❌ (use 'context' instead)
#
# ==============================================================================

# ==============================================================================
# Configuration
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0129/tnk_analysis_unified/results/subcluster_unified_v2_20260129/adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0130/tcell_interpret_v1_3_CONTEXT_FIXED"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"

# DeepSeek API Key (from environment variable)
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
if (nchar(DEEPSEEK_API_KEY) < 10) {
  stop(
    "ERROR: DEEPSEEK_API_KEY environment variable not set or invalid\n",
    "Please set it: export DEEPSEEK_API_KEY='your-key-here'"
  )
}

# Analysis Parameters
N_CORES <- 4
TOP_N_MARKERS <- 50

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

# ==============================================================================
# ⭐ FIX: LLM Interpretation - Annotation Task
# ==============================================================================

cat("\n=== LLM Interpretation: Annotation Task ===\n")

enrich_databases <- list()

if (!is.null(go_bp_enrich)) {
  enrich_databases$go_bp <- go_bp_enrich
}
if (!is.null(go_mf_enrich)) {
  enrich_databases$go_mf <- go_mf_enrich
}
if (!is.null(go_cc_enrich)) {
  enrich_databases$go_cc <- go_cc_enrich
}
if (!is.null(hallmark_enrich)) {
  enrich_databases$hallmark <- hallmark_enrich
}
if (!is.null(msigdb_kegg_enrich)) {
  enrich_databases$kegg <- msigdb_kegg_enrich
}
if (!is.null(cellmarker_enrich)) {
  enrich_databases$cellmarker <- cellmarker_enrich
}
if (!is.null(panglaodb_enrich)) {
  enrich_databases$panglaodb <- panglaodb_enrich
}

cat(sprintf(
  "Using %d enrichment databases for interpretation\n",
  length(enrich_databases)
))

annotation_results <- NULL

if (length(enrich_databases) > 0) {
  cat("\n[INFO] Calling interpret() for annotation task...\n")
  cat("[NOTE] v1.3 FIX: Adding biological context back!\n\n")

  annotation_results <- tryCatch(
    {
      # ⭐ v1.3 FIX: Added back context parameter (IT IS SUPPORTED!)
      interpret(
        x = enrich_databases,
        context = paste(
          "CRITICAL CONTEXT FOR T/NK CELL ANNOTATION:",
          "",
          "These are T/NK cell subclusters from respiratory tissues (CRSwNP study).",
          "",
          "Key T/NK cell biology principles:",
          "",
          "1. CD4+ vs CD8+ DISTINCTION (CRITICAL):",
          "   - CD4+ T cells: Helper T cell lineage",
          "     * Th1 (TBX21+, IFNG+, CXCR3+): Anti-viral, cellular immunity",
          "     * Th2 (GATA3+, IL4/IL5/IL13+, CCR4+): Allergy, Type 2 inflammation",
          "     * Th17 (RORC+, IL17A+, CCR6+): Neutrophil recruitment, barrier immunity",
          "     * Tfh (BCL6+, CXCR5+, PDCD1+): B cell help, germinal center",
          "     * Treg (FOXP3+, IL2RA+, CTLA4+): Immunosuppression, tolerance",
          "   - CD8+ T cells: Cytotoxic T cell lineage",
          "     * Cytotoxic effectors (GZMB+, PRF1+, NKG7+): Direct cell killing",
          "     * TRM (CD69+, ITGAE/CD103+, S1PR1low): Tissue-resident memory",
          "     * TEMRA (KLRG1+, CD45RA+): Terminal effector memory",
          "",
          "2. NK CELLS (CD3- cytotoxic lymphocytes):",
          "   - CD16+ NK (FCGR3A+): Mature, potent cytotoxicity",
          "   - CD16- NK (FCGR3A-): Immature or regulatory",
          "   - Key markers: NCAM1, NKG7, GNLY, KLRD1, KLRC1",
          "",
          "3. MEMORY STATES:",
          "   - Naïve: CCR7+, SELL+, CD45RA+, TCF7+ (antigen-inexperienced)",
          "   - TCM (Central memory): CCR7+, CD45RO+ (lymphoid trafficking)",
          "   - TEM (Effector memory): CCR7-, CD45RO+ (peripheral circulation)",
          "   - TEMRA (Terminal effector): CD45RA+, KLRG1+ (mainly CD8+)",
          "   - TRM (Tissue-resident): CD69+, ITGAE+, S1PR1low (tissue-locked)",
          "",
          "4. FUNCTIONAL STATES:",
          "   - Resting: Low activation markers",
          "   - Activated: CD69+, CD25+ (IL2RA), HLA-DR+",
          "   - Exhausted: PDCD1+, LAG3+, TIGIT+, HAVCR2+, TOX+",
          "   - Proliferating: MKI67+, TOP2A+, STMN1+",
          "   - Cytotoxic program: GZMB, PRF1, NKG7, GNLY, GZMK, GZMA",
          "",
          "ANNOTATION INSTRUCTIONS:",
          "For EACH cluster, specify:",
          "  a) Lineage: CD4+ / CD8+ / NK / Innate-like T",
          "  b) Memory state: Naïve / TCM / TEM / TEMRA / TRM",
          "  c) Functional state: Resting / Activated / Exhausted / Proliferating",
          "  d) Effector phenotype: Cytotoxic / Th1 / Th2 / Th17 / Tfh / Treg / NK",
          sep = "\n"
        ),
        task = "annotation",
        model = "deepseek-reasoner",
        api_key = DEEPSEEK_API_KEY,
        n_pathways = 15
      )
    },
    error = function(e) {
      cat(
        "[ERROR] Annotation interpretation failed:",
        conditionMessage(e),
        "\n"
      )
      NULL
    }
  )

  if (!is.null(annotation_results)) {
    saveRDS(
      annotation_results,
      file.path(OUTPUT_DIR, "reports", "annotation_results.rds")
    )

    cat("[OK] Annotation results saved\n")

    # Extract structured data
    annotation_df_list <- lapply(names(annotation_results), function(cid) {
      result <- annotation_results[[cid]]

      data.frame(
        Cluster = cid,
        Cell_Type = ifelse(!is.null(result$cell_type), result$cell_type, ""),
        Confidence = ifelse(!is.null(result$confidence), result$confidence, ""),
        Markers = ifelse(
          !is.null(result$markers),
          paste(result$markers, collapse = "; "),
          ""
        ),
        Regulatory_Drivers = ifelse(
          !is.null(result$regulatory_drivers),
          paste(result$regulatory_drivers, collapse = "; "),
          ""
        ),
        Reasoning = ifelse(!is.null(result$reasoning), result$reasoning, ""),
        stringsAsFactors = FALSE
      )
    })

    annotation_df <- bind_rows(annotation_df_list)

    write.csv(
      annotation_df,
      file.path(OUTPUT_DIR, "annotation_results.csv"),
      row.names = FALSE
    )

    cat("[OK] Annotation CSV exported\n")
  }
} else {
  cat("[WARN] No enrichment databases available for annotation\n")
}

# ==============================================================================
# ⭐ FIX: LLM Interpretation - Phenotype Task
# ==============================================================================

cat("\n=== LLM Interpretation: Phenotype Task ===\n")

phenotype_results <- NULL

if (length(enrich_databases) > 0) {
  cat("\n[INFO] Calling interpret() for phenotyping task...\n")

  phenotype_results <- tryCatch(
    {
      # ⭐ v1.3 FIX: Added back context parameter
      interpret(
        x = enrich_databases,
        context = paste(
          "CRITICAL CONTEXT FOR T/NK CELL FUNCTIONAL PHENOTYPING:",
          "",
          "These are T/NK cell functional states in respiratory inflammation (CRSwNP).",
          "",
          "Key functional phenotypes to consider:",
          "",
          "1. CD4+ T CELL FUNCTIONS:",
          "   - Th1 effectors: IFN-γ production, macrophage activation, anti-viral",
          "   - Th2 effectors: IL-4/5/13 production, eosinophil recruitment, IgE promotion",
          "   - Th17 effectors: IL-17A/F production, neutrophil recruitment, barrier defense",
          "   - Tfh: B cell help, antibody class switching, germinal center formation",
          "   - Treg: Suppression (IL-10, TGF-β), maintaining tolerance",
          "",
          "2. CD8+ T CELL FUNCTIONS:",
          "   - Cytotoxic killing: Granzyme/perforin-mediated target cell lysis",
          "   - IFN-γ production: Antiviral and inflammatory responses",
          "   - Tissue surveillance: TRM-mediated rapid local immunity",
          "   - Memory recall: Fast reactivation upon reinfection",
          "",
          "3. NK CELL FUNCTIONS:",
          "   - Innate cytotoxicity: Direct killing without prior sensitization",
          "   - ADCC (CD16+): Antibody-dependent cellular cytotoxicity",
          "   - Cytokine production: IFN-γ, TNF for immune amplification",
          "   - Immunoregulation: Shaping adaptive immune responses",
          "",
          "4. ACTIVATION VS EXHAUSTION:",
          "   - Activated state: Recent TCR engagement, proliferation, effector functions",
          "   - Exhausted state: Chronic stimulation, dampened effector functions, high PD-1/LAG3/TIGIT",
          "",
          "PHENOTYPING INSTRUCTIONS:",
          "For EACH cluster, describe:",
          "  a) Primary effector function (cytotoxic/helper/regulatory/innate)",
          "  b) Key biological processes",
          "  c) Activation/exhaustion state",
          "  d) Tissue localization strategy (resident vs circulating)",
          "  e) Relevance to respiratory disease pathology",
          sep = "\n"
        ),
        task = "phenotyping",
        model = "deepseek-chat",
        api_key = DEEPSEEK_API_KEY,
        n_pathways = 30
      )
    },
    error = function(e) {
      cat("[ERROR] Phenotype interpretation failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(phenotype_results)) {
    saveRDS(
      phenotype_results,
      file.path(OUTPUT_DIR, "reports", "phenotype_results.rds")
    )

    cat("[OK] Phenotype results saved\n")

    # Extract structured data
    phenotype_df_list <- lapply(names(phenotype_results), function(cid) {
      result <- phenotype_results[[cid]]

      data.frame(
        Cluster = cid,
        Functional_Phenotype = ifelse(
          !is.null(result$phenotype),
          result$phenotype,
          ""
        ),
        Confidence = ifelse(!is.null(result$confidence), result$confidence, ""),
        Key_Processes = ifelse(
          !is.null(result$key_processes),
          paste(result$key_processes, collapse = "; "),
          ""
        ),
        Regulatory_Drivers = ifelse(
          !is.null(result$regulatory_drivers),
          paste(result$regulatory_drivers, collapse = "; "),
          ""
        ),
        Reasoning = ifelse(!is.null(result$reasoning), result$reasoning, ""),
        Network_Evidence = ifelse(
          !is.null(result$network_evidence),
          result$network_evidence,
          ""
        ),
        stringsAsFactors = FALSE
      )
    })

    phenotype_df <- bind_rows(phenotype_df_list)

    write.csv(
      phenotype_df,
      file.path(OUTPUT_DIR, "phenotype_results.csv"),
      row.names = FALSE
    )

    cat("[OK] Phenotype CSV exported\n")
  }
} else {
  cat("[WARN] No enrichment databases available for phenotyping\n")
}

# ==============================================================================
# ⭐ FIX: Per-Celltype Detailed Interpretation
# ==============================================================================

cat("\n=== Per-Celltype Detailed Interpretation ===\n")

celltype_interpretations <- list()

if (!is.null(seurat_obj$cell_type_L2)) {
  unique_celltypes <- unique(seurat_obj$cell_type_L2)

  cat(sprintf(
    "Interpreting %d unique cell types (L2)...\n",
    length(unique_celltypes)
  ))

  for (celltype in unique_celltypes) {
    cat(sprintf("\n--- Interpreting: %s ---\n", celltype))

    # Get subclusters of this celltype
    subclusters <- seurat_obj@meta.data %>%
      filter(cell_type_L2 == celltype) %>%
      pull(cell_type_L3) %>%
      unique()

    if (length(subclusters) == 0) {
      cat("[WARN] No subclusters found, skipping\n")
      next
    }

    cat(sprintf(
      "Found %d subclusters: %s\n",
      length(subclusters),
      paste(subclusters, collapse = ", ")
    ))

    # Filter enrichment results to this celltype's subclusters
    celltype_enrich <- list()

    for (db_name in names(enrich_databases)) {
      db_result <- enrich_databases[[db_name]]

      if (!is.null(db_result) && "compareClusterResult" %in% class(db_result)) {
        ccr <- db_result@compareClusterResult
        ccr_filtered <- ccr %>% filter(Cluster %in% subclusters)

        if (nrow(ccr_filtered) > 0) {
          filtered_result <- db_result
          filtered_result@compareClusterResult <- ccr_filtered
          celltype_enrich[[db_name]] <- filtered_result
        }
      }
    }

    if (length(celltype_enrich) == 0) {
      cat("[WARN] No enrichment results for this celltype, skipping\n")
      next
    }

    cat(sprintf(
      "Using %d enrichment databases for this celltype\n",
      length(celltype_enrich)
    ))

    # Interpret this celltype
    interpretation <- tryCatch(
      {
        # ⭐ v1.3 FIX: Added back context parameter
        interpret(
          x = celltype_enrich,
          context = paste(
            sprintf("DETAILED INTERPRETATION FOR: %s", celltype),
            sprintf(
              "SUBCLUSTERS (%d): %s",
              length(subclusters),
              paste(subclusters, collapse = ", ")
            ),
            "",
            "Provide a comprehensive biological interpretation:",
            "1. Cell identity and lineage (CD4+/CD8+/NK/innate-like T)",
            "2. Functional state and activation status",
            "3. Memory differentiation stage",
            "4. Biological role in respiratory tissue immunity",
            "5. Subcluster heterogeneity and what drives it",
            "6. Disease relevance in CRSwNP pathogenesis",
            "7. Comparison to related T/NK cell subtypes",
            "",
            "Connect findings across all available databases (GO, Hallmark, CellMarker, etc.)",
            sep = "\n"
          ),
          task = "annotation",
          model = "deepseek-chat",
          api_key = DEEPSEEK_API_KEY,
          n_pathways = 25
        )
      },
      error = function(e) {
        cat("[ERROR] Interpretation failed:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (!is.null(interpretation)) {
      celltype_interpretations[[celltype]] <- interpretation

      # Save as text file
      txt_file <- file.path(
        OUTPUT_DIR,
        "reports",
        paste0(gsub(" ", "_", celltype), "_interpretation.txt")
      )

      tryCatch(
        {
          sink(txt_file)
          cat(sprintf("# %s - Detailed Interpretation\n\n", celltype))
          cat(sprintf(
            "**Subclusters:** %s\n\n",
            paste(subclusters, collapse = ", ")
          ))

          for (subcluster in names(interpretation)) {
            cat(sprintf("## %s\n\n", subcluster))
            result <- interpretation[[subcluster]]

            if (!is.null(result$cell_type)) {
              cat(sprintf("**Cell Type:** %s\n", result$cell_type))
            }
            if (!is.null(result$confidence)) {
              cat(sprintf("**Confidence:** %s\n\n", result$confidence))
            }
            if (!is.null(result$reasoning)) {
              cat("**Reasoning:**\n")
              cat(result$reasoning, "\n\n")
            }

            cat("---\n\n")
          }

          sink()
          cat(sprintf("[OK] Saved interpretation: %s\n", basename(txt_file)))
        },
        error = function(e) {
          tryCatch(sink(), error = function(e) NULL)
          cat("[WARN] Failed to save text file:", conditionMessage(e), "\n")
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
      "\n[OK] Saved %d per-celltype interpretations\n",
      length(celltype_interpretations)
    ))
  }
} else {
  cat(
    "[WARN] No cell_type_L2 column found, skipping per-celltype interpretation\n"
  )
}

# ==============================================================================
# Generate Summary Report
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

report_file <- file.path(OUTPUT_DIR, "REPORT.md")

tryCatch(
  {
    sink(report_file)

    cat("# T/NK Cell Subcluster Interpretation Report v1.3\n\n")
    cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
    cat("**Databases:** CellMarker, PanglaoDB, GO BP/MF/CC, Hallmark, KEGG\n\n")
    cat("**LLM Model:** DeepSeek (via interpret() with biological context)\n\n")
    cat("**Analysis Method:** compareCluster (per-cluster enrichment)\n\n")
    cat("**✨ v1.3 NEW:** Full T/NK cell biological context restored!\n\n")
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

    # Annotation results
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
    }

    # Phenotype results
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
    }

    # Per-celltype interpretations
    if (length(celltype_interpretations) > 0) {
      cat("## Per-Celltype Detailed Interpretations\n\n")

      for (celltype in names(celltype_interpretations)) {
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
# Completion
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("T/NK CELL ANALYSIS COMPLETE - v1.3 (CONTEXT FIXED)\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat(
  "  - REPORT.md                          Comprehensive interpretation report\n"
)
cat("  - annotation_results.csv             Cell_Type + Markers ⭐\n")
cat(
  "  - phenotype_results.csv              Functional_Phenotype + Key_Processes ⭐\n"
)
cat("  - all_markers.csv                    All significant marker genes\n")
cat(
  "  - top_markers_filtered.csv           Filtered markers used for enrichment\n"
)
cat("\n")

cat("Figures:\n")
cat("  - figures/cellmarker_dotplot.pdf     CellMarker database enrichment\n")
cat("  - figures/panglaodb_dotplot.pdf      PanglaoDB enrichment\n")
cat("  - figures/go_bp_dotplot.pdf          GO biological process\n")
cat("  - figures/go_mf_dotplot.pdf          GO molecular function\n")
cat("  - figures/go_cc_dotplot.pdf          GO cellular component\n")
cat("  - figures/hallmark_dotplot.pdf       Hallmark pathways\n")
cat("  - figures/msigdb_kegg_dotplot.pdf    KEGG pathways\n")
cat("\n")

cat("Version history:\n")
cat("  v1.3 (2026-01-30): ✅ CORRECT! Added back context parameter (P0-6)\n")
cat("                     - interpret() DOES support 'context' and 'task'!\n")
cat("                     - Now with full T/NK cell biological context\n")
cat("  v1.2 (2026-01-30): ❌ Over-correction (removed ALL parameters)\n")
cat("  v1.1 (2026-01-29): ❌ Used fake parameters (save_history, etc.)\n")
cat("\n")

cat("Key improvements in v1.3:\n")
cat("  ✓ Biological context restored for all 3 interpret() calls\n")
cat("  ✓ CD4+/CD8+ distinction explicitly described\n")
cat("  ✓ Memory states (Naïve/TCM/TEM/TEMRA/TRM) defined\n")
cat("  ✓ Functional states (activated/exhausted) specified\n")
cat("  ✓ T/NK cell-specific effector phenotypes (Th1/Th2/Th17/Tfh/Treg/NK)\n")
cat("  ✓ Disease context (CRSwNP respiratory inflammation) included\n")
cat("\n")

cat(
  "================================================================================\n"
)
cat("DONE\n")
cat(
  "================================================================================\n"
)
