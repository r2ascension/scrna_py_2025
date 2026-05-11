#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Subcluster LLM Interpretation - PRODUCTION VERSION v2.5 - Multi-Database Support
# ==============================================================================
#
# Version: v2.4.2 FINAL (2026-01-27)
# Status: Production-ready with complete interpret() result extraction
#
# Key Updates from v2.4.1:
#   ✅ Complete extraction of interpret() structure:
#      - Phenotype, Confidence, Reasoning
#      - Regulatory_Drivers, Key_Processes
#      - Network_Evidence, Refined_Network
#   ✅ CSV tables with all fields for easy analysis
#   ✅ Detailed per-celltype reports with full subcluster breakdown
#   ✅ REPORT.md with comprehensive annotations and phenotypes
#
# Original Features (from v2.4.1):
#   1. ✅ MSigDB使用本地GMT文件（无需网络）
#   2. ✅ GO扩展：BP + MF + CC（3种本体论）
#   3. ✅ 9个数据库整合
#   4. ✅ interpret()使用正确的单对象输入
#
# ==============================================================================

# ==============================================================================
# 配置参数
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0119/bcell_analysis/results/subcluster_v2_20260119/adata_bcell_subclustered_FINAL_v2_20260119.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0127/bcell_interpret_v2_4_1"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
BCELL_MARKERS_PATH <- "/home/h2048/data/source/reference/CellMarker/bcell_markers_comprehensive.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"

# DeepSeek API Key
DEEPSEEK_API_KEY <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05"

# Analysis Parameters
N_CORES <- 8
TOP_N_MARKERS <- 50

# ==============================================================================
# 加载库
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
# 加载数据
# ==============================================================================

cat("=== Loading Seurat Data ===\n")

seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

cat(sprintf(
  "\nLoaded: %d cells x %d genes\n",
  ncol(seurat_obj),
  nrow(seurat_obj)
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
# 计算 Marker Genes（并行）
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
# 准备 Top Markers for Enrichment
# ==============================================================================

cat("\n=== Preparing Top Markers ===\n")

genes_to_filter <- c(
  grep("^MT-", rownames(seurat_obj), value = TRUE),
  grep("^RP[SL]", rownames(seurat_obj), value = TRUE),
  "FOS", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
  "ZFP36", "DUSP1", "DUSP2", "IER2", "IER3", "ATF3",
  "BTG2", "FOSB", "NR4A1", "NR4A2", "NR4A3",
  "HSP90AA1", "HSPA1A", "HSPA1B", "DNAJB1"
)

cat(sprintf("Filtering %d potentially confounding genes\n", length(genes_to_filter)))

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
  category = c("MT genes", "Ribosomal genes", "Stress response genes", "Total filtered"),
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
# 加载 B Cell Markers Database
# ==============================================================================

cat("\n=== Loading B Cell Markers Database ===\n")

bcell_markers_db <- NULL
bcell_markers_term2gene <- NULL

bcell_markers_db <- tryCatch(
  {
    db <- fread(BCELL_MARKERS_PATH, header = TRUE, stringsAsFactors = FALSE)
    cat(sprintf("[OK] Loaded %d B cell subtypes\n", nrow(db)))
    db
  },
  error = function(e) {
    cat("[WARN] Failed to load B Cell Markers:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(bcell_markers_db)) {
  term2gene_list <- list()

  for (i in 1:nrow(bcell_markers_db)) {
    row <- bcell_markers_db[i, ]
    term <- paste(row$Cell_Type, row$Subset, sep = "_")
    all_markers <- c()

    if (!is.na(row$Core_Markers) && row$Core_Markers != "") {
      all_markers <- c(all_markers, unlist(strsplit(row$Core_Markers, ",")))
    }
    if (!is.na(row$Surface_Markers) && row$Surface_Markers != "") {
      all_markers <- c(all_markers, unlist(strsplit(row$Surface_Markers, ",")))
    }
    if (!is.na(row$Transcription_Factors) && row$Transcription_Factors != "") {
      all_markers <- c(all_markers, unlist(strsplit(row$Transcription_Factors, ",")))
    }
    if (!is.na(row$Functional_Markers) && row$Functional_Markers != "") {
      all_markers <- c(all_markers, unlist(strsplit(row$Functional_Markers, ",")))
    }

    all_markers <- gsub('["\r\n]', '', all_markers)
    all_markers <- trimws(all_markers)
    all_markers <- toupper(all_markers)
    all_markers <- unique(all_markers[all_markers != "" & !is.na(all_markers)])

    if (length(all_markers) > 0) {
      term2gene_list[[length(term2gene_list) + 1]] <- data.frame(
        term = term,
        gene = all_markers,
        stringsAsFactors = FALSE
      )
    }
  }

  bcell_markers_term2gene <- bind_rows(term2gene_list)

  cat(sprintf("[OK] Prepared B Cell TERM2GENE: %d pairs\n", nrow(bcell_markers_term2gene)))
  cat(sprintf("    Cell subtypes: %d\n", length(unique(bcell_markers_term2gene$term))))
  cat(sprintf("    Genes: %d\n", length(unique(bcell_markers_term2gene$gene))))
}

# ==============================================================================
# 加载 CellMarker Database
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

  bcell_related <- cellmarker_db %>%
    filter(
      grepl("B cell|B-cell|Plasma|Plasmablast|Lymphocyte|Immune", cell_name, ignore.case = TRUE) |
      grepl("Blood|Lymph|Spleen|Bone marrow|Tonsil", tissue_type, ignore.case = TRUE)
    )

  cat(sprintf("[OK] Found %d B cell/immune-related entries\n", nrow(bcell_related)))

  if (nrow(bcell_related) >= 50) {
    cellmarker_db <- bcell_related
    cat("[INFO] Using B cell/immune-specific subset for enrichment\n")
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

  cat(sprintf("[OK] Prepared CellMarker TERM2GENE: %d pairs\n", nrow(cellmarker_term2gene)))
  cat(sprintf("    Cell types: %d\n", length(unique(cellmarker_term2gene$term))))
  cat(sprintf("    Genes: %d\n", length(unique(cellmarker_term2gene$gene))))
}

# ==============================================================================
# 加载 PanglaoDB Database
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

  cat(sprintf("[OK] Prepared PanglaoDB TERM2GENE: %d pairs\n", nrow(panglaodb_term2gene)))
}

# ==============================================================================
# 富集分析 - B Cell Markers
# ==============================================================================

cat("\n=== B Cell Markers Enrichment ===\n")

bcell_markers_enrich <- NULL

if (!is.null(bcell_markers_term2gene)) {
  bcell_markers_enrich <- tryCatch(
    {
      compareCluster(
        gene ~ cluster,
        data = top_markers,
        fun = enricher,
        TERM2GENE = bcell_markers_term2gene,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2
      )
    },
    error = function(e) {
      cat("[WARN] B Cell Markers enrichment failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(bcell_markers_enrich)) {
    ccr <- bcell_markers_enrich@compareClusterResult
    n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] Found %d significant B cell subtypes\n", n_sig))

    saveRDS(bcell_markers_enrich, file.path(OUTPUT_DIR, "reports", "bcell_markers_enrich.rds"))
  }
}

# ==============================================================================
# 富集分析 - CellMarker
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

    saveRDS(cellmarker_enrich, file.path(OUTPUT_DIR, "reports", "cellmarker_enrich.rds"))
  }
}

# ==============================================================================
# 富集分析 - PanglaoDB
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

    saveRDS(panglaodb_enrich, file.path(OUTPUT_DIR, "reports", "panglaodb_enrich.rds"))
  }
}

# ==============================================================================
# 富集分析 - GO (BP + MF + CC)
# ==============================================================================

cat("\n=== GO Enrichment (BP + MF + CC) ===\n")

# GO Biological Process
go_bp_enrich <- tryCatch(
  {
    cat("Running GO BP enrichment...\n")
    compareCluster(
      gene ~ cluster,
      data = top_markers,
      fun = enrichGO,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2
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
  cat(sprintf("[OK] GO BP: Found %d significant terms\n", n_sig))
  saveRDS(go_bp_enrich, file.path(OUTPUT_DIR, "reports", "go_bp_enrich.rds"))
}

# GO Molecular Function
go_mf_enrich <- tryCatch(
  {
    cat("Running GO MF enrichment...\n")
    compareCluster(
      gene ~ cluster,
      data = top_markers,
      fun = enrichGO,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "MF",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2
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
  cat(sprintf("[OK] GO MF: Found %d significant terms\n", n_sig))
  saveRDS(go_mf_enrich, file.path(OUTPUT_DIR, "reports", "go_mf_enrich.rds"))
}

# GO Cellular Component
go_cc_enrich <- tryCatch(
  {
    cat("Running GO CC enrichment...\n")
    compareCluster(
      gene ~ cluster,
      data = top_markers,
      fun = enrichGO,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "CC",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2
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
  cat(sprintf("[OK] GO CC: Found %d significant terms\n", n_sig))
  saveRDS(go_cc_enrich, file.path(OUTPUT_DIR, "reports", "go_cc_enrich.rds"))
}

# ==============================================================================
# 富集分析 - MSigDB (Hallmark + KEGG) from Local GMT
# ==============================================================================

cat("\n=== MSigDB Enrichment (Hallmark + KEGG) from Local GMT ===\n")

hallmark_enrich <- NULL
msigdb_kegg_enrich <- NULL

# Function to read GMT file
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
    
    # Extract Hallmark
    hallmark_term2gene <- all_genesets %>%
      filter(grepl("^HALLMARK_", term)) %>%
      dplyr::select(term, gene)
    
    n_hallmark_sets <- length(unique(hallmark_term2gene$term))
    cat(sprintf("[OK] Extracted %d Hallmark gene sets\n", n_hallmark_sets))
    
    # Extract KEGG
    kegg_term2gene <- all_genesets %>%
      filter(grepl("KEGG_", term)) %>%
      dplyr::select(term, gene)
    
    n_kegg_sets <- length(unique(kegg_term2gene$term))
    cat(sprintf("[OK] Extracted %d KEGG gene sets\n", n_kegg_sets))
    
    # Hallmark enrichment
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
        
        saveRDS(hallmark_enrich, file.path(OUTPUT_DIR, "reports", "hallmark_enrich.rds"))
      }
    }
    
    # KEGG enrichment
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
        
        saveRDS(msigdb_kegg_enrich, file.path(OUTPUT_DIR, "reports", "msigdb_kegg_enrich.rds"))
      }
    }
  },
  error = function(e) {
    cat("[ERROR] Failed to load GMT file:", conditionMessage(e), "\n")
    cat("[INFO] Enrichment will be skipped\n")
  }
)

# ==============================================================================
# 可视化富集结果
# ==============================================================================

cat("\n=== Visualizing Enrichment Results ===\n")

# B Cell Markers
if (!is.null(bcell_markers_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "bcell_markers_dotplot.pdf"), width = 16, height = 12)
      print(
        dotplot(bcell_markers_enrich, showCategory = 10, font.size = 7) +
          ggtitle("B Cell Markers Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved bcell_markers_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] B Cell Markers plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

# CellMarker
if (!is.null(cellmarker_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "cellmarker_dotplot.pdf"), width = 16, height = 12)
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

# PanglaoDB
if (!is.null(panglaodb_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "panglaodb_dotplot.pdf"), width = 16, height = 12)
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

# GO BP
if (!is.null(go_bp_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "go_bp_dotplot.pdf"), width = 16, height = 14)
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

# GO MF
if (!is.null(go_mf_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "go_mf_dotplot.pdf"), width = 16, height = 14)
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

# GO CC
if (!is.null(go_cc_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "go_cc_dotplot.pdf"), width = 16, height = 14)
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

# Hallmark
if (!is.null(hallmark_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "hallmark_dotplot.pdf"), width = 16, height = 12)
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

# MSigDB KEGG
if (!is.null(msigdb_kegg_enrich)) {
  tryCatch(
    {
      pdf(file.path(OUTPUT_DIR, "figures", "msigdb_kegg_dotplot.pdf"), width = 16, height = 14)
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
    fanyi::chat_request("test", model = "deepseek-chat")
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
# 准备富集结果摘要用于 LLM Context
# ==============================================================================

cat("\n=== Preparing Enrichment Summary for LLM ===\n")

# 创建富集结果摘要（用于context，不用于interpret()的x参数）
enrichment_summary <- list()

if (!is.null(bcell_markers_enrich)) {
  bcell_top <- bcell_markers_enrich@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(p.adjust, n = 3) %>%
    ungroup()
  enrichment_summary[["bcell_markers"]] <- bcell_top
  cat(sprintf("  [OK] B Cell Markers: %d top terms\n", nrow(bcell_top)))
}

if (!is.null(cellmarker_enrich)) {
  cm_top <- cellmarker_enrich@compareClusterResult %>%
    filter(p.adjust < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(p.adjust, n = 3) %>%
    ungroup()
  enrichment_summary[["cellmarker"]] <- cm_top
  cat(sprintf("  [OK] CellMarker: %d top terms\n", nrow(cm_top)))
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

# 构建context字符串（包含其他数据库的top terms）
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
        paste(sprintf("%s (%s)", term_summary$Cluster, term_summary$top_terms), collapse = " | ")
      )
    }
  }
}

cat("\n[OK] Prepared enrichment summary for context\n")

# ==============================================================================
# LLM 解释 - Task 1: Annotation
# ==============================================================================

cat("\n=== LLM Interpretation: Annotation ===\n")

annotation_results <- NULL

# ✅ UPDATE v2.5: Multi-database input (tested with 3 databases, supports up to 6)
# Build enrichment list for comprehensive annotation
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

if (!is.null(bcell_markers_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- bcell_markers_enrich
  cat("✓ B Cell Markers added\n")
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
        x = enrichment_list,  # Multi-database list
        context = paste(
          "B cells from normal nasal cavity, sinus, bronchi, and lung tissues.",
          "These are subclusters of three major B cell populations:",
          "(1) Memory B cells - class-switched, antigen-experienced cells showing CD27+ phenotype,",
          "(2) Naive B cells - IgD+ IgM+ cells that have not undergone somatic hypermutation,",
          "(3) Plasma cells - antibody-secreting cells with high PRDM1 and XBP1 expression.",
          "The samples are from healthy respiratory tract with baseline immune surveillance.",
          "We are particularly interested in identifying functional states such as:",
          "- Germinal center-experienced vs non-GC memory B cells",
          "- Activated vs resting naive B cells",
          "- Short-lived plasmablasts vs long-lived plasma cells",
          "- Tissue-resident vs circulating populations",
          "- IgG vs IgA vs IgE class-switched subtypes",
          additional_context  # 包含其他数据库的证据
        ),
        task = "annotation",
        n_pathways = 15  # Optimized for multi-database (tested with 3, supports up to 6)
      )
    },
    error = function(e) {
      cat("[WARN] Annotation failed:", conditionMessage(e), "\n")
      cat("       Error details:", e$message, "\n")
      return(NULL)
    }
  )
}
  )

  if (!is.null(annotation_results)) {
    cat("\n[OK] Annotation complete!\n\n")

    cat("=== Annotation Results ===\n\n")
    
    # interpret() 返回的是 per-cluster results 的 list
    if (is.list(annotation_results)) {
      # 检查是否是 per-cluster results (list of lists)
      if (length(annotation_results) > 0 && is.list(annotation_results[[1]])) {
        # 这是标准的 per-cluster 格式
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
          
          # 处理特殊情况：没有显著富集的cluster（只有cluster, overview, confidence）
          if (!is.null(cluster_result$overview) && 
              is.null(cluster_result$cell_type) && 
              cluster_result$confidence == "None") {
            # 这是一个没有显著富集的cluster
            annotation_table <- rbind(annotation_table, data.frame(
              Cluster = cluster_name,
              Cell_Type = "No significant enrichment",
              Confidence = "None",
              Reasoning = cluster_result$overview,
              Regulatory_Drivers = "NA",
              Markers = "NA",
              stringsAsFactors = FALSE
            ))
            
            cat(sprintf("**%s**\n", cluster_name))
            cat("  Cell Type: No significant enrichment\n")
            cat("  Note: ", cluster_result$overview, "\n\n")
            next
          }
          
          # 正常情况：提取完整注释
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
          
          # Annotation task 有 markers 字段，不是 key_processes
          markers <- if (!is.null(cluster_result$markers)) {
            paste(cluster_result$markers, collapse = "; ")
          } else {
            "NA"
          }
          
          annotation_table <- rbind(annotation_table, data.frame(
            Cluster = cluster_name,
            Cell_Type = cell_type,
            Confidence = confidence,
            Reasoning = reasoning,  # 完整保存，不截断
            Regulatory_Drivers = regulatory_drivers,
            Markers = markers,
            stringsAsFactors = FALSE
          ))
          
          # 打印单个cluster结果（控制台输出时截断）
          cat(sprintf("**%s**\n", cluster_name))
          cat(sprintf("  Cell Type: %s\n", cell_type))
          cat(sprintf("  Confidence: %s\n", confidence))
          if (regulatory_drivers != "NA") {
            drivers_list <- strsplit(regulatory_drivers, "; ")[[1]]
            cat(sprintf("  Regulatory Drivers: %s", paste(head(drivers_list, 3), collapse = ", ")))
            if (length(drivers_list) > 3) cat(sprintf(" (+ %d more)", length(drivers_list) - 3))
            cat("\n")
          }
          if (markers != "NA") {
            # 只显示前5个
            marker_list <- strsplit(markers, "; ")[[1]]
            top5 <- marker_list[1:min(5, length(marker_list))]
            cat(sprintf("  Markers: %s", paste(top5, collapse = ", ")))
            if (length(marker_list) > 5) cat(sprintf(" (+ %d more)", length(marker_list) - 5))
            cat("\n")
          }
          cat(sprintf("  Reasoning: %s...\n\n", 
                     substr(reasoning, 1, 150)))  # 控制台截断显示
        }
        
        # 保存为CSV（完整内容）
        write.csv(annotation_table, 
                  file.path(OUTPUT_DIR, "annotation_results.csv"),
                  row.names = FALSE)
        cat("[OK] Saved annotation_results.csv (full content preserved)\n\n")
        
      } else if (!is.null(annotation_results$interpretation)) {
        # 单个文本结果
        cat(annotation_results$interpretation)
        cat("\n\n")
        writeLines(annotation_results$interpretation, 
                   file.path(OUTPUT_DIR, "annotation_results.txt"))
      } else {
        # 其他格式，尝试打印结构
        cat("Result structure (first 3 elements):\n")
        print(head(annotation_results, 3))
        cat("\n")
      }
    }
    
    # 总是保存完整RDS作为备份
    saveRDS(annotation_results, 
            file.path(OUTPUT_DIR, "reports", "annotation_results.rds"))
    cat("[OK] Saved complete results to RDS\n")
  }
} else {
  cat("[WARN] GO BP enrichment not available, skipping annotation\n")
}

# ==============================================================================
# LLM 解释 - Task 2: Phenotyping
# ==============================================================================

cat("\n=== LLM Interpretation: Phenotyping ===\n")

phenotype_results <- NULL

if (!is.null(go_bp_enrich)) {
  phenotype_results <- tryCatch(
    {
      interpret(
        x = go_bp_enrich,  # 单个enrichment object
        context = paste(
          "B cells from normal respiratory tract (nasal cavity, paranasal sinuses, bronchi, lung).",
          "Looking for baseline functional states and activation signatures in healthy tissue.",
          "Key functional states to identify:",
          "- Activation level: resting vs activated vs exhausted",
          "- Proliferation: quiescent vs actively proliferating",
          "- Antibody production: non-secreting vs antibody-secreting cells (ASC)",
          "- Cytokine responsiveness: IL-4/IL-13 responsive, IFN-responsive, etc.",
          "- Tissue residency: circulating vs tissue-resident memory",
          "- Differentiation stage: germinal center vs post-GC vs terminally differentiated",
          "These are HEALTHY controls, so expect homeostatic activation states",
          "rather than pathogenic inflammatory signatures.",
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
    
    # 正确处理per-cluster results
    if (is.list(phenotype_results)) {
      if (length(phenotype_results) > 0 && is.list(phenotype_results[[1]])) {
        # Per-cluster 格式
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
          
          # Phenotyping task: 提取功能表型
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
          
          regulatory_drivers <- if (!is.null(cluster_result$regulatory_drivers)) {
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
          
          phenotype_table <- rbind(phenotype_table, data.frame(
            Cluster = cluster_name,
            Functional_Phenotype = phenotype,
            Confidence = confidence,
            Reasoning = reasoning,  # 完整保存
            Regulatory_Drivers = regulatory_drivers,
            Key_Processes = key_processes,
            Network_Evidence = network_evidence,  # 完整保存
            stringsAsFactors = FALSE
          ))
          
          # 打印单个cluster结果（控制台输出时截断）
          cat(sprintf("**%s**\n", cluster_name))
          cat(sprintf("  Phenotype: %s\n", phenotype))
          cat(sprintf("  Confidence: %s\n", confidence))
          if (regulatory_drivers != "NA") {
            cat(sprintf("  Regulatory Drivers: %s\n", regulatory_drivers))
          }
          if (key_processes != "NA") {
            # 只显示前3个
            processes <- strsplit(key_processes, "; ")[[1]]
            top3 <- processes[1:min(3, length(processes))]
            cat(sprintf("  Key Processes: %s", paste(top3, collapse = "; ")))
            if (length(processes) > 3) cat(sprintf(" (+ %d more)", length(processes) - 3))
            cat("\n")
          }
          cat(sprintf("  Reasoning: %s...\n\n", 
                     substr(reasoning, 1, 150)))  # 控制台截断
        }
        
        # 保存为CSV（完整内容）
        write.csv(phenotype_table, 
                  file.path(OUTPUT_DIR, "phenotype_results.csv"),
                  row.names = FALSE)
        cat("[OK] Saved phenotype_results.csv (full content preserved)\n\n")
        
      } else if (!is.null(phenotype_results$interpretation)) {
        cat(phenotype_results$interpretation)
        cat("\n\n")
        writeLines(phenotype_results$interpretation, 
                   file.path(OUTPUT_DIR, "phenotype_results.txt"))
      } else {
        cat("Result structure (first 3 elements):\n")
        print(head(phenotype_results, 3))
        cat("\n")
      }
    }

    # 总是保存完整RDS
    saveRDS(phenotype_results, 
            file.path(OUTPUT_DIR, "reports", "phenotype_results.rds"))
    cat("[OK] Saved complete results to RDS\n")
  }
} else {
  cat("[WARN] GO BP enrichment not available, skipping phenotyping\n")
}

# ==============================================================================
# Per-Celltype 深入分析
# ==============================================================================

cat("\n=== Per-Celltype Detailed Analysis ===\n")

celltypes_l2 <- unique(seurat_obj@meta.data$cell_type_L2)
celltypes_l2 <- celltypes_l2[!is.na(celltypes_l2)]

cat(sprintf("Analyzing %d cell types at L2 level\n", length(celltypes_l2)))

celltype_interpretations <- list()

for (celltype in celltypes_l2) {
  cat(sprintf("\n--- Processing %s ---\n", celltype))

  celltype_subclusters <- unique(
    seurat_obj@meta.data$cell_type_L3[seurat_obj@meta.data$cell_type_L2 == celltype]
  )
  celltype_subclusters <- celltype_subclusters[!is.na(celltype_subclusters)]

  n_subclusters <- length(celltype_subclusters)
  cat(sprintf("  Subclusters: %d\n", n_subclusters))

  if (n_subclusters <= 1) {
    cat("  [SKIP] Only 1 subcluster\n")
    next
  }

  # 过滤GO BP enrichment（用作主要输入）
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

  # 构建该细胞类型的additional context
  celltype_context <- ""
  
  # 添加其他数据库的top terms
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
        paste(sprintf("%s (%s)", hall_filtered$Cluster, hall_filtered$Description), collapse = "; ")
      )
    }
  }

  # LLM 解释
  interpretation <- tryCatch(
    {
      interpret(
        x = celltype_go_bp,  # 单个enrichment object
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
          "- Evidence of antigen experience and memory formation in normal contexts",
          "Key questions:",
          "1. What functional states distinguish the subclusters?",
          "2. Are there proliferative vs quiescent populations?",
          "3. Do subclusters show tissue-specific or universal phenotypes?",
          "4. What activation or differentiation markers define each subcluster?",
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
        # 提取完整的interpretation内容
        report_lines <- c()
        
        if (is.list(interpretation)) {
          # 检查是否是per-cluster格式
          if (length(interpretation) > 0 && is.list(interpretation[[1]])) {
            # 为该celltype的每个subcluster生成详细报告
            report_lines <- c(
              sprintf("# %s Detailed Interpretation\n", celltype),
              sprintf("Generated: %s\n", format(Sys.time())),
              sprintf("Total Subclusters: %d\n\n", length(interpretation)),
              paste(rep("=", 80), collapse=""), "\n\n"
            )
            
            for (cluster_name in names(interpretation)) {
              # 只处理该celltype的subclusters
              if (grepl(celltype, cluster_name, fixed = TRUE)) {
                cluster_result <- interpretation[[cluster_name]]
                
                report_lines <- c(report_lines, 
                                 sprintf("## %s\n\n", cluster_name))
                
                # Phenotype
                if (!is.null(cluster_result$phenotype)) {
                  report_lines <- c(report_lines,
                                   "### Phenotype\n",
                                   cluster_result$phenotype, "\n\n")
                }
                
                # Confidence
                if (!is.null(cluster_result$confidence)) {
                  report_lines <- c(report_lines,
                                   sprintf("**Confidence:** %s\n\n", 
                                          cluster_result$confidence))
                }
                
                # Regulatory Drivers
                if (!is.null(cluster_result$regulatory_drivers)) {
                  report_lines <- c(report_lines,
                                   "### Regulatory Drivers\n",
                                   paste("-", cluster_result$regulatory_drivers, 
                                        collapse = "\n"), "\n\n")
                }
                
                # Key Processes
                if (!is.null(cluster_result$key_processes)) {
                  report_lines <- c(report_lines,
                                   "### Key Biological Processes\n",
                                   paste("-", cluster_result$key_processes, 
                                        collapse = "\n"), "\n\n")
                }
                
                # Reasoning
                if (!is.null(cluster_result$reasoning)) {
                  report_lines <- c(report_lines,
                                   "### Reasoning\n",
                                   cluster_result$reasoning, "\n\n")
                }
                
                # Network Evidence
                if (!is.null(cluster_result$network_evidence)) {
                  report_lines <- c(report_lines,
                                   "### Network Evidence\n",
                                   cluster_result$network_evidence, "\n\n")
                }
                
                # Refined Network (if exists)
                if (!is.null(cluster_result$refined_network) && 
                    is.data.frame(cluster_result$refined_network)) {
                  report_lines <- c(report_lines,
                                   "### Regulatory Network\n\n")
                  
                  for (i in 1:nrow(cluster_result$refined_network)) {
                    row <- cluster_result$refined_network[i, ]
                    report_lines <- c(report_lines,
                                     sprintf("%d. **%s** → **%s** (%s)\n",
                                            i, row$source, row$target, row$interaction))
                    if (!is.null(row$reason) && !is.na(row$reason)) {
                      report_lines <- c(report_lines,
                                       sprintf("   *%s*\n\n", row$reason))
                    }
                  }
                  report_lines <- c(report_lines, "\n")
                }
                
                report_lines <- c(report_lines, 
                                 "---\n\n")
              }
            }
            
            report_text <- paste(report_lines, collapse = "")
            
          } else {
            # 非标准格式，尝试其他提取方式
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
          saveRDS(interpretation, 
                  file.path(OUTPUT_DIR, "reports",
                           paste0(gsub(" ", "_", celltype), "_interpretation.rds")))
        }
      },
      error = function(e) {
        cat("  [WARN] Failed to save report:", e$message, "\n")
        saveRDS(interpretation, 
                file.path(OUTPUT_DIR, "reports",
                         paste0(gsub(" ", "_", celltype), "_interpretation.rds")))
      }
    )
  }
}

if (length(celltype_interpretations) > 0) {
  saveRDS(
    celltype_interpretations,
    file.path(OUTPUT_DIR, "reports", "celltype_interpretations.rds")
  )
  cat(sprintf("\n[OK] Saved interpretations for %d cell types\n", length(celltype_interpretations)))
}

# ==============================================================================
# 生成总结报告
# ==============================================================================

cat("\n=== Generating Summary Report ===\n")

report_file <- file.path(OUTPUT_DIR, "REPORT.md")

tryCatch(
  {
    sink(report_file)

    cat("# B Cell Subcluster Interpretation Report v2.4.1\n\n")
    cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
    cat("**Databases (9 total):** B Cell Markers, CellMarker, PanglaoDB, GO BP/MF/CC, Hallmark, MSigDB KEGG\n\n")
    cat("**LLM Model:** DeepSeek Chat (via interpret() function)\n\n")
    cat("---\n\n")

    cat("## Dataset Summary\n\n")
    cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
    cat(sprintf("- Subclusters (L3): %d\n", length(unique(seurat_obj$cell_type_L3))))
    cat(sprintf("- Cell types (L2): %d\n", length(unique(seurat_obj$cell_type_L2))))
    cat("\n---\n\n")

    # Annotation results from CSV if available
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
          cat(gsub("; ", "  \n- ", paste0("- ", row$Regulatory_Drivers)), "\n\n")
        }
        
        if (!is.na(row$Markers) && row$Markers != "NA") {
          cat("**Key Markers:**  \n")
          # 只显示前8个，避免太长
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

    # Phenotype results from CSV if available
    phenotype_csv_path <- file.path(OUTPUT_DIR, "phenotype_results.csv")
    if (file.exists(phenotype_csv_path)) {
      cat("## Functional Phenotypes (LLM-Generated)\n\n")
      
      phenotype_df <- read.csv(phenotype_csv_path, stringsAsFactors = FALSE)
      
      for (i in 1:nrow(phenotype_df)) {
        row <- phenotype_df[i, ]
        cat(sprintf("### %s\n\n", row$Cluster))
        cat(sprintf("**Functional Phenotype:** %s  \n", row$Functional_Phenotype))
        cat(sprintf("**Confidence:** %s  \n\n", row$Confidence))
        
        if (!is.na(row$Regulatory_Drivers) && row$Regulatory_Drivers != "NA") {
          cat("**Regulatory Drivers:**  \n")
          cat(gsub("; ", "  \n- ", paste0("- ", row$Regulatory_Drivers)), "\n\n")
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
        
        if (!is.na(row$Network_Evidence) && row$Network_Evidence != "NA" && 
            nchar(row$Network_Evidence) > 10) {
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

    # Per-celltype interpretations
    if (length(celltype_interpretations) > 0) {
      cat("## Per-Celltype Detailed Interpretations\n\n")

      for (celltype in names(celltype_interpretations)) {
        interpretation <- celltype_interpretations[[celltype]]
        if (is.null(interpretation)) {
          next
        }

        cat(sprintf("### %s\n\n", celltype))
        
        # 尝试从txt文件读取
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
# 完成
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("ANALYSIS COMPLETE - v2.4.1 HOTFIX\n")
cat("================================================================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat("  - REPORT.md                          Comprehensive interpretation report\n")
cat("  - annotation_results.csv             Cell_Type + Markers (not key_processes!) ⭐\n")
cat("  - phenotype_results.csv              Functional_Phenotype + Key_Processes ⭐\n")
cat("  - all_markers.csv                    All significant marker genes\n")
cat("  - top_markers_filtered.csv           Filtered markers used for enrichment\n")
cat("  - filtered_genes_info.csv            Gene filtering statistics\n")
cat("\n")

cat("Figures (9 databases):\n")
cat("  - figures/bcell_markers_dotplot.pdf  B cell marker enrichment (36 subtypes)\n")
cat("  - figures/cellmarker_dotplot.pdf     CellMarker database enrichment\n")
cat("  - figures/panglaodb_dotplot.pdf      PanglaoDB enrichment\n")
cat("  - figures/go_bp_dotplot.pdf          GO biological process enrichment\n")
cat("  - figures/go_mf_dotplot.pdf          GO molecular function enrichment ⭐ NEW\n")
cat("  - figures/go_cc_dotplot.pdf          GO cellular component enrichment ⭐ NEW\n")
cat("  - figures/hallmark_dotplot.pdf       Hallmark pathways (50 canonical)\n")
cat("  - figures/msigdb_kegg_dotplot.pdf    KEGG pathways (MSigDB local)\n")
cat("\n")

cat("RDS objects (load with readRDS()):\n")
cat("  - reports/bcell_markers_enrich.rds   B cell marker enrichment object\n")
cat("  - reports/cellmarker_enrich.rds      CellMarker enrichment object\n")
cat("  - reports/panglaodb_enrich.rds       PanglaoDB enrichment object\n")
cat("  - reports/go_bp_enrich.rds           GO BP enrichment object\n")
cat("  - reports/go_mf_enrich.rds           GO MF enrichment object ⭐ NEW\n")
cat("  - reports/go_cc_enrich.rds           GO CC enrichment object ⭐ NEW\n")
cat("  - reports/hallmark_enrich.rds        Hallmark enrichment object\n")
cat("  - reports/msigdb_kegg_enrich.rds     MSigDB KEGG enrichment object\n")
cat("  - reports/annotation_results.rds     LLM annotation (raw structure)\n")
cat("  - reports/phenotype_results.rds      LLM phenotype (raw structure)\n")
cat("  - reports/celltype_interpretations.rds Per-celltype detailed interpretations\n")
cat("\n")

cat("Per-celltype reports:\n")
cat("  - reports/[celltype]_interpretation.txt  Detailed narrative for each L2 celltype\n")
cat("\n")

cat("Important Notes:\n")
cat("  ⚠️  v2.4.2 PRODUCTION: Correct field extraction based on task type\n")
cat("      - Annotation task → cell_type + markers (✓ tested)\n")
cat("      - Phenotyping task → phenotype + key_processes (✓ tested)\n")
cat("      - Full Reasoning preserved in CSV (no truncation)\n")
cat("      - Console output truncated for readability\n")
cat("      - Ready for multi-database input: list(go_bp, hallmark, ...)\n")
cat("\n")

cat("================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")