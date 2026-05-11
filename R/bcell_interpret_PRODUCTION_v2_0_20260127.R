#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Subcluster LLM Interpretation - PRODUCTION VERSION
# ==============================================================================
#
# Version: v2.0 (2026-01-27)
# Status: Production-ready with QUICK_REFERENCE_MEMORY v2.13 compliance
#
# 关键修正点（基于QUICK_REFERENCE_MEMORY v2.13）：
#   1. ✅ interpret()使用list输入而非合并的compareClusterResult
#   2. ✅ 详细的context描述（200+词）提高注释准确性
#   3. ✅ 过滤MT/ribo/stress基因避免技术偏差
#   4. ✅ 整合B Cell Markers Comprehensive数据库（36个亚型）
#   5. ✅ 统一API设置（仅一次，位置优化）
#   6. ✅ 并行处理（future）加速marker计算
#   7. ✅ 多数据库证据整合（B Cell + PanglaoDB + GO）
#
# 运行方式：
#   - 在 RStudio 中逐段运行（Ctrl+Enter）推荐
#   - 或直接 source("bcell_interpret_FIXED_v2_20260127.R")
#
# ==============================================================================

# ==============================================================================
# 配置参数
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0119/bcell_analysis/results/subcluster_v2_20260119/adata_bcell_subclustered_FINAL_v2_20260119.h5ad"
OUTPUT_DIR <- "/home/h2048/data/r/0127/bcell_interpret_fixed"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
BCELL_MARKERS_PATH <- "/home/h2048/data/source/reference/CellMarker/bcell_markers_comprehensive.csv"

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
          test.use = "wilcoxon",
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

# 合并结果
all_markers <- bind_rows(lapply(names(marker_list), function(cid) {
  df <- marker_list[[cid]]
  df$cluster <- cid
  df$gene <- rownames(df)
  df
}))

all_markers <- all_markers %>% filter(p_val_adj < 0.05)

cat(sprintf("[OK] Found %d significant markers\n", nrow(all_markers)))

# 保存
write.csv(
  all_markers,
  file.path(OUTPUT_DIR, "all_markers.csv"),
  row.names = FALSE
)

# ==============================================================================
# 准备 Top Markers for Enrichment
# ==============================================================================

cat("\n=== Preparing Top Markers ===\n")

# 定义需要过滤的基因（MT/ribo/stress响应基因）
genes_to_filter <- c(
  # 线粒体基因
  grep("^MT-", rownames(seurat_obj), value = TRUE),
  # 核糖体基因
  grep("^RP[SL]", rownames(seurat_obj), value = TRUE),
  # 立即早期应激响应基因
  "FOS", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "ZFP36",
  "DUSP1", "DUSP2", "IER2", "IER3", "ATF3", "BTG2", "FOSB", "NR4A1",
  "NR4A2", "NR4A3", "HSP90AA1", "HSPA1A", "HSPA1B", "DNAJB1"
)

cat(sprintf("Filtering %d potentially confounding genes\n", length(genes_to_filter)))

top_markers <- all_markers %>%
  filter(!gene %in% genes_to_filter) %>%  # 过滤干扰基因
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = TOP_N_MARKERS) %>%
  ungroup() %>%
  mutate(gene = toupper(gene)) %>%
  select(gene, cluster)

cat(sprintf("Selected top %d clean markers per cluster\n", TOP_N_MARKERS))
cat(sprintf("Total markers for enrichment: %d\n", nrow(top_markers)))

# 保存过滤后的marker列表（用于富集）
write.csv(
  top_markers,
  file.path(OUTPUT_DIR, "top_markers_filtered.csv"),
  row.names = FALSE
)

# 保存过滤信息
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
# 加载 B Cell Markers Comprehensive Database
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
  # 构建 TERM2GENE（展开所有marker列）
  term2gene_list <- list()
  
  # 处理每一行
  for (i in 1:nrow(bcell_markers_db)) {
    row <- bcell_markers_db[i, ]
    
    # 组合 Cell_Type 和 Subset 作为 term
    term <- paste(row$Cell_Type, row$Subset, sep = "_")
    
    # 提取所有marker类型的基因
    all_markers <- c()
    
    # Core markers
    if (!is.na(row$Core_Markers) && row$Core_Markers != "") {
      core <- unlist(strsplit(row$Core_Markers, ","))
      all_markers <- c(all_markers, core)
    }
    
    # Surface markers
    if (!is.na(row$Surface_Markers) && row$Surface_Markers != "") {
      surface <- unlist(strsplit(row$Surface_Markers, ","))
      all_markers <- c(all_markers, surface)
    }
    
    # Transcription factors
    if (!is.na(row$Transcription_Factors) && row$Transcription_Factors != "") {
      tfs <- unlist(strsplit(row$Transcription_Factors, ","))
      all_markers <- c(all_markers, tfs)
    }
    
    # Functional markers
    if (!is.na(row$Functional_Markers) && row$Functional_Markers != "") {
      func <- unlist(strsplit(row$Functional_Markers, ","))
      all_markers <- c(all_markers, func)
    }
    
    # 清理基因名（去除引号、空格、特殊字符）
    all_markers <- gsub('["\r\n]', '', all_markers)
    all_markers <- trimws(all_markers)
    all_markers <- toupper(all_markers)
    all_markers <- unique(all_markers[all_markers != "" & !is.na(all_markers)])
    
    # 创建 term-gene 对
    if (length(all_markers) > 0) {
      term2gene_list[[length(term2gene_list) + 1]] <- data.frame(
        term = term,
        gene = all_markers,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # 合并所有 term-gene 对
  bcell_markers_term2gene <- bind_rows(term2gene_list)
  
  cat(sprintf("[OK] Prepared B Cell TERM2GENE: %d pairs\n", nrow(bcell_markers_term2gene)))
  cat(sprintf("    Cell subtypes: %d\n", length(unique(bcell_markers_term2gene$term))))
  cat(sprintf("    Genes: %d\n", length(unique(bcell_markers_term2gene$gene))))
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
    select(cell_type, gene_symbol) %>%
    mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct() %>%
    rename(term = cell_type, gene = gene_symbol)
  
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
    
    saveRDS(
      bcell_markers_enrich,
      file.path(OUTPUT_DIR, "reports", "bcell_markers_enrich.rds")
    )
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
    
    saveRDS(
      panglaodb_enrich,
      file.path(OUTPUT_DIR, "reports", "panglaodb_enrich.rds")
    )
  }
}

# ==============================================================================
# 富集分析 - GO Biological Process
# ==============================================================================

cat("\n=== GO Enrichment ===\n")

go_enrich <- tryCatch(
  {
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
    cat("[WARN] GO enrichment failed:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(go_enrich)) {
  ccr <- go_enrich@compareClusterResult
  n_sig <- sum(ccr$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf("[OK] Found %d significant GO terms\n", n_sig))
  
  saveRDS(go_enrich, file.path(OUTPUT_DIR, "reports", "go_enrich.rds"))
}

# ==============================================================================
# 可视化富集结果
# ==============================================================================

cat("\n=== Visualizing Enrichment Results ===\n")

# B Cell Markers
if (!is.null(bcell_markers_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "bcell_markers_dotplot.pdf"),
        width = 16,
        height = 12
      )
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

# PanglaoDB
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

# GO
if (!is.null(go_enrich)) {
  tryCatch(
    {
      pdf(
        file.path(OUTPUT_DIR, "figures", "go_dotplot.pdf"),
        width = 16,
        height = 14
      )
      print(
        dotplot(go_enrich, showCategory = 15, font.size = 6) +
          ggtitle("GO BP Enrichment") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      dev.off()
      cat("[OK] Saved go_dotplot.pdf\n")
    },
    error = function(e) {
      cat("[WARN] GO plot failed\n")
      tryCatch(dev.off(), error = function(e) NULL)
    }
  )
}

# ==============================================================================
# 配置 DeepSeek API（统一位置，仅一次）
# ==============================================================================

cat("\n=== Configuring DeepSeek API ===\n")

# 设置环境变量
Sys.setenv(DEEPSEEK_API_KEY = DEEPSEEK_API_KEY)

# 配置 fanyi 包
fanyi::set_translate_option(
  key = DEEPSEEK_API_KEY,
  source = "deepseek"
)

# 测试连接
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
# 准备富集结果列表用于 LLM 解释
# ==============================================================================

cat("\n=== Preparing Enrichment Results for LLM ===\n")

# 使用list格式（推荐方式，符合clusterProfiler::interpret()文档）
enrich_list <- list()

if (!is.null(bcell_markers_enrich)) {
  enrich_list[["B_Cell_Markers"]] <- bcell_markers_enrich
  cat(sprintf("  [OK] B Cell Markers: %d terms\n", 
              nrow(bcell_markers_enrich@compareClusterResult)))
}

if (!is.null(panglaodb_enrich)) {
  enrich_list[["PanglaoDB"]] <- panglaodb_enrich
  cat(sprintf("  [OK] PanglaoDB: %d terms\n", 
              nrow(panglaodb_enrich@compareClusterResult)))
}

if (!is.null(go_enrich)) {
  enrich_list[["GO_BP"]] <- go_enrich
  cat(sprintf("  [OK] GO BP: %d terms\n", 
              nrow(go_enrich@compareClusterResult)))
}

if (length(enrich_list) == 0) {
  cat("[WARN] No enrichment results available for LLM interpretation\n")
} else {
  cat(sprintf("\n[OK] Prepared %d enrichment sources for LLM\n", length(enrich_list)))
}

# ==============================================================================
# LLM 解释 - Task 1: Annotation
# ==============================================================================

cat("\n=== LLM Interpretation: Annotation ===\n")

annotation_results <- NULL

if (length(enrich_list) > 0) {
  annotation_results <- tryCatch(
    {
      interpret(
        enrich_list,
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
          "- IgG vs IgA vs IgE class-switched subtypes"
        ),
        task = "annotation"
      )
    },
    error = function(e) {
      cat("[WARN] Annotation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (!is.null(annotation_results)) {
    cat("\n[OK] Annotation complete!\n\n")
    
    # 打印摘要
    for (nm in names(annotation_results)) {
      res <- annotation_results[[nm]]
      ct_print <- if (is.null(res$cell_type)) "Unknown" else res$cell_type
      conf_print <- if (is.null(res$confidence)) "None" else res$confidence
      
      cat(sprintf("%-40s: %s (%s)\n", nm, ct_print, conf_print))
    }
    
    # 清洗并保存结果
    clean_data_list <- lapply(names(annotation_results), function(sub_name) {
      res <- annotation_results[[sub_name]]
      
      ct <- if (is.null(res$cell_type)) "Unknown" else res$cell_type
      
      if (ct == "Unknown") {
        return(NULL)
      }
      
      conf <- if (is.null(res$confidence)) "NA" else res$confidence
      
      return(data.frame(
        subcluster = sub_name,
        cell_type = ct,
        confidence = conf,
        stringsAsFactors = FALSE
      ))
    })
    
    annotation_df <- dplyr::bind_rows(clean_data_list)
    
    if (nrow(annotation_df) > 0) {
      write.csv(
        annotation_df,
        file.path(OUTPUT_DIR, "annotations.csv"),
        row.names = FALSE
      )
      cat("\n[OK] Saved annotations.csv\n")
    }
    
    saveRDS(
      annotation_results,
      file.path(OUTPUT_DIR, "reports", "annotation_results.rds")
    )
    
    # 验证置信度分布
    cat("\nConfidence Distribution:\n")
    conf_counts <- table(sapply(annotation_results, function(x) {
      if (is.null(x$confidence)) "NA" else x$confidence
    }))
    print(conf_counts)
    
    # 警告：低置信度占比过高
    low_conf_count <- sum(conf_counts[c("Low", "NA")], na.rm = TRUE)
    total_count <- sum(conf_counts)
    if (low_conf_count / total_count > 0.3) {
      cat("\n⚠️  WARNING: >30% clusters have low/NA confidence!\n")
      cat("   Recommendations:\n")
      cat("   1. Check if enrichment databases cover your cell types\n")
      cat("   2. Consider lowering pvalueCutoff in enrichment (e.g., 0.1)\n")
      cat("   3. Increase TOP_N_MARKERS (e.g., 100 instead of 50)\n")
      cat("   4. Manually verify low-confidence annotations\n\n")
    }
  }
}

# ==============================================================================
# LLM 解释 - Task 2: Phenotyping
# ==============================================================================

cat("\n=== LLM Interpretation: Phenotyping ===\n")

phenotype_results <- NULL

if (length(enrich_list) > 0) {
  phenotype_results <- tryCatch(
    {
      interpret(
        enrich_list,
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
          "rather than pathogenic inflammatory signatures."
        ),
        task = "phenotyping"
      )
    },
    error = function(e) {
      cat("[WARN] Phenotyping failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (!is.null(phenotype_results)) {
    cat("\n[OK] Phenotyping complete!\n\n")
    
    # 打印摘要
    for (nm in names(phenotype_results)) {
      res <- phenotype_results[[nm]]
      ph_print <- if (is.null(res$phenotype)) "Unknown" else res$phenotype
      conf_print <- if (is.null(res$confidence)) "None" else res$confidence
      
      cat(sprintf("%-40s: %s (%s)\n", nm, ph_print, conf_print))
    }
    
    # 清洗并保存结果
    clean_pheno_list <- lapply(names(phenotype_results), function(sub_name) {
      res <- phenotype_results[[sub_name]]
      
      if (is.null(res)) {
        return(NULL)
      }
      
      ph <- if (is.null(res$phenotype)) "Unknown" else res$phenotype
      
      if (ph == "Unknown") {
        return(NULL)
      }
      
      conf <- if (is.null(res$confidence)) "NA" else res$confidence
      
      return(data.frame(
        subcluster = sub_name,
        phenotype = ph,
        confidence = conf,
        stringsAsFactors = FALSE
      ))
    })
    
    phenotype_df <- dplyr::bind_rows(clean_pheno_list)
    
    if (nrow(phenotype_df) > 0) {
      write.csv(
        phenotype_df,
        file.path(OUTPUT_DIR, "phenotypes.csv"),
        row.names = FALSE
      )
      cat("\n[OK] Saved phenotypes.csv\n")
    }
    
    saveRDS(
      phenotype_results,
      file.path(OUTPUT_DIR, "reports", "phenotype_results.rds")
    )
  }
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
  
  # 获取该细胞类型的所有L3亚群
  celltype_subclusters <- unique(
    seurat_obj@meta.data$cell_type_L3[seurat_obj@meta.data$cell_type_L2 == celltype]
  )
  celltype_subclusters <- celltype_subclusters[!is.na(celltype_subclusters)]
  
  n_subclusters <- length(celltype_subclusters)
  cat(sprintf("  Subclusters: %d\n", n_subclusters))
  
  # 跳过只有1个亚群的
  if (n_subclusters <= 1) {
    cat("  [SKIP] Only 1 subcluster\n")
    next
  }
  
  # 构建该细胞类型的enrichment list
  celltype_enrich_list <- list()
  
  # 过滤B Cell Markers enrichment
  if (!is.null(bcell_markers_enrich)) {
    bcell_filtered <- bcell_markers_enrich@compareClusterResult %>%
      filter(Cluster %in% celltype_subclusters)
    
    if (nrow(bcell_filtered) > 0) {
      bcell_subset <- bcell_markers_enrich
      bcell_subset@compareClusterResult <- bcell_filtered
      celltype_enrich_list[["B_Cell_Markers"]] <- bcell_subset
    }
  }
  
  # 过滤PanglaoDB enrichment
  if (!is.null(panglaodb_enrich)) {
    pang_filtered <- panglaodb_enrich@compareClusterResult %>%
      filter(Cluster %in% celltype_subclusters)
    
    if (nrow(pang_filtered) > 0) {
      pang_subset <- panglaodb_enrich
      pang_subset@compareClusterResult <- pang_filtered
      celltype_enrich_list[["PanglaoDB"]] <- pang_subset
    }
  }
  
  # 过滤GO enrichment
  if (!is.null(go_enrich)) {
    go_filtered <- go_enrich@compareClusterResult %>%
      filter(Cluster %in% celltype_subclusters)
    
    if (nrow(go_filtered) > 0) {
      go_subset <- go_enrich
      go_subset@compareClusterResult <- go_filtered
      celltype_enrich_list[["GO_BP"]] <- go_subset
    }
  }
  
  # 检查是否有富集结果
  if (length(celltype_enrich_list) == 0) {
    cat("  [SKIP] No enrichment results\n")
    next
  }
  
  cat(sprintf("  Enrichment sources: %d\n", length(celltype_enrich_list)))
  
  # LLM 解释
  interpretation <- tryCatch(
    {
      interpret(
        celltype_enrich_list,
        context = paste(
          celltype, "cells from normal respiratory tract (nasal cavity, sinus, bronchi, lung).",
          "Focus on functional heterogeneity among", celltype, "subclusters.",
          "These are healthy baseline cells, so we expect:",
          "- Homeostatic activation states rather than pathogenic inflammation",
          "- Normal differentiation trajectories and maturation stages",
          "- Tissue-specific adaptations to mucosal immune surveillance",
          "- Evidence of antigen experience and memory formation in normal contexts",
          "Key questions:",
          "1. What functional states distinguish the subclusters?",
          "2. Are there proliferative vs quiescent populations?",
          "3. Do subclusters show tissue-specific or universal phenotypes?",
          "4. What activation or differentiation markers define each subcluster?"
        ),
        task = "interpretation"
      )
    },
    error = function(e) {
      cat("  [WARN] Interpretation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (!is.null(interpretation)) {
    celltype_interpretations[[celltype]] <- interpretation
    cat("  [OK] Interpretation complete\n")
    
    # 保存个别报告
    report_file <- file.path(
      OUTPUT_DIR,
      "reports",
      paste0(gsub(" ", "_", celltype), "_interpretation.txt")
    )
    
    tryCatch(
      {
        sink(report_file)
        cat("# ", celltype, " Interpretation\n\n", sep = "")
        cat("Generated: ", format(Sys.time()), "\n\n", sep = "")
        
        if (!is.null(interpretation$overview)) {
          cat("## Overview\n\n", interpretation$overview, "\n\n", sep = "")
        }
        
        if (!is.null(interpretation$key_mechanisms)) {
          cat("## Key Mechanisms\n\n", interpretation$key_mechanisms, "\n\n", sep = "")
        }
        
        if (!is.null(interpretation$narrative)) {
          cat("## Narrative\n\n", interpretation$narrative, "\n\n", sep = "")
        }
        
        sink()
        cat("  [OK] Saved report\n")
      },
      error = function(e) {
        tryCatch(sink(), error = function(e) NULL)
        cat("  [WARN] Failed to save report\n")
      }
    )
  }
}

# 保存所有细胞类型解释
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
    
    cat("# B Cell Subcluster Interpretation Report\n\n")
    cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
    cat("**Databases:** B Cell Markers Comprehensive, PanglaoDB, GO BP\n\n")
    cat("---\n\n")
    
    cat("## Dataset Summary\n\n")
    cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
    cat(sprintf("- Subclusters (L3): %d\n", length(unique(seurat_obj$cell_type_L3))))
    cat(sprintf("- Cell types (L2): %d\n", length(unique(seurat_obj$cell_type_L2))))
    cat("\n---\n\n")
    
    # Annotation Section
    if (!is.null(annotation_results)) {
      cat("## Cell Subtype Annotations\n\n")
      
      for (nm in names(annotation_results)) {
        result <- annotation_results[[nm]]
        if (is.null(result)) next
        
        ct <- if (is.null(result$cell_type)) "Unknown" else result$cell_type
        conf <- if (is.null(result$confidence)) "NA" else result$confidence
        reason <- if (is.null(result$reasoning)) "No reasoning provided." else result$reasoning
        
        cat(sprintf("### %s\n\n", nm))
        cat(sprintf("**Cell Type:** %s  \n", ct))
        cat(sprintf("**Confidence:** %s  \n\n", conf))
        cat("**Reasoning:**\n\n", reason, "\n\n", sep = "")
        cat("---\n\n")
      }
    }
    
    # Phenotype Section
    if (!is.null(phenotype_results)) {
      cat("## Functional Phenotypes\n\n")
      
      for (nm in names(phenotype_results)) {
        result <- phenotype_results[[nm]]
        if (is.null(result)) next
        
        ph <- if (is.null(result$phenotype)) "Unknown" else result$phenotype
        conf <- if (is.null(result$confidence)) "NA" else result$confidence
        
        cat(sprintf("### %s\n\n", nm))
        cat(sprintf("**Phenotype:** %s  \n", ph))
        cat(sprintf("**Confidence:** %s  \n\n", conf))
        cat("---\n\n")
      }
    }
    
    # Per-Celltype Interpretations
    if (length(celltype_interpretations) > 0) {
      cat("## Per-Celltype Detailed Interpretations\n\n")
      
      for (celltype in names(celltype_interpretations)) {
        interpretation <- celltype_interpretations[[celltype]]
        if (is.null(interpretation)) next
        
        cat(sprintf("### %s\n\n", celltype))
        
        content <- if (!is.null(interpretation$narrative)) {
          interpretation$narrative
        } else if (!is.null(interpretation$overview)) {
          interpretation$overview
        } else {
          "No detailed interpretation available."
        }
        
        cat(content, "\n\n", sep = "")
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
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat("  - REPORT.md                          Comprehensive interpretation report\n")
cat("  - annotations.csv                    Cell type annotations with confidence\n")
cat("  - phenotypes.csv                     Functional phenotypes\n")
cat("  - all_markers.csv                    All significant marker genes\n")
cat("  - top_markers_filtered.csv           Filtered markers used for enrichment\n")
cat("  - filtered_genes_info.csv            Gene filtering statistics\n")
cat("\n")

cat("Figures:\n")
cat("  - figures/bcell_markers_dotplot.pdf  B cell marker enrichment (36 subtypes)\n")
cat("  - figures/panglaodb_dotplot.pdf      PanglaoDB enrichment\n")
cat("  - figures/go_dotplot.pdf             GO biological process enrichment\n")
cat("\n")

cat("RDS objects (load with readRDS()):\n")
cat("  - reports/bcell_markers_enrich.rds   B cell marker enrichment object\n")
cat("  - reports/panglaodb_enrich.rds       PanglaoDB enrichment object\n")
cat("  - reports/go_enrich.rds              GO enrichment object\n")
cat("  - reports/annotation_results.rds     LLM annotation results (full)\n")
cat("  - reports/phenotype_results.rds      LLM phenotype results (full)\n")
cat("  - reports/celltype_interpretations.rds Per-celltype detailed interpretations\n")
cat("\n")

cat("Per-celltype reports:\n")
cat("  - reports/[celltype]_interpretation.txt  Detailed narrative for each L2 celltype\n")
cat("\n")

if (!is.null(annotation_results)) {
  cat("Confidence distribution (Annotation):\n")
  conf_table <- table(sapply(annotation_results, function(x) {
    if (is.null(x$confidence)) "NA" else x$confidence
  }))
  for (conf in names(conf_table)) {
    cat(sprintf("  %-10s: %d\n", conf, conf_table[conf]))
  }
}

cat("\n")
cat("================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")
