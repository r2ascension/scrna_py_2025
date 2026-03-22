#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Subcluster LLM Interpretation - SIMPLE PANGLAODB VERSION
# ==============================================================================
#
# 特点：
#   - 无函数封装，直接线性执行
#   - 无限制条件，遇错跳过
#   - PanglaoDB数据库（更全面）
#   - 直接调用 interpret()
#   - 可视化结构，逐步运行
#
# 运行方式：
#   - 在 RStudio 中逐段运行（Ctrl+Enter）
#   - 或直接 source()
#
# ==============================================================================

# ==============================================================================
# 配置参数（修改这里）
# ==============================================================================

H5AD_PATH <- "/home/h2048/data/py/0119/bcell_analysis/results/subcluster_v2_20260119/adata_bcell_subclustered_FINAL_v2_20260119.h5ad"
OUTPUT_DIR <- "/home/h2048/data/r/0126/bcell_interpret_panglaodb"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"

# API Key
Sys.setenv(DEEPSEEK_API_KEY = "sk-ed1879cf6fa14b04aac9cb6c078a3d05")

# 参数
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

# Setup API
set_translate_api("deepseek", api_key = Sys.getenv("DEEPSEEK_API_KEY"))

# Create output dirs
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reports"), showWarnings = FALSE)

cat("[OK] Setup complete\n\n")

# ==============================================================================
# 加载数据
# ==============================================================================

cat("=== Loading Data ===\n")

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
DefaultAssay(seurat_obj) <- "RNA"
# 并行 FindMarkers
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

cat(sprintf("[OK] Found %d markers\n", nrow(all_markers)))

# 保存
write.csv(
  all_markers,
  file.path(OUTPUT_DIR, "all_markers.csv"),
  row.names = FALSE
)

# ==============================================================================
# 准备 Top Markers
# ==============================================================================

cat("\n=== Preparing Top Markers ===\n")

top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = TOP_N_MARKERS) %>%
  ungroup() %>%
  mutate(gene = toupper(gene)) %>% # 确保大写
  select(gene, cluster)

cat(sprintf("Selected top %d markers per cluster\n", TOP_N_MARKERS))

# ==============================================================================
# 加载 PanglaoDB 数据库
# ==============================================================================

cat("\n=== Loading PanglaoDB Database ===\n")

panglaodb_db <- NULL

panglaodb_db <- tryCatch(
  {
    db <- fread(PANGLAODB_PATH, header = TRUE, stringsAsFactors = FALSE)

    # 重命名列（移除空格）
    setnames(
      db,
      old = c("official gene symbol", "cell type"),
      new = c("gene_symbol", "cell_type"),
      skip_absent = TRUE
    )

    # 过滤人类markers
    db <- db %>% filter(grepl("Hs", species, fixed = TRUE))

    cat(sprintf("[OK] Loaded %d human markers\n", nrow(db)))

    db
  },
  error = function(e) {
    cat("[WARN] Failed to load PanglaoDB:", conditionMessage(e), "\n")
    cat("[INFO] Continuing without PanglaoDB\n")
    return(NULL)
  }
)

# 准备 TERM2GENE
panglaodb_term2gene <- NULL

if (!is.null(panglaodb_db)) {
  panglaodb_term2gene <- panglaodb_db %>%
    select(cell_type, gene_symbol) %>%
    mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct() %>%
    rename(term = cell_type, gene = gene_symbol)

  cat(sprintf("[OK] Prepared TERM2GENE: %d pairs\n", nrow(panglaodb_term2gene)))
  cat(sprintf("    Cell types: %d\n", length(unique(panglaodb_term2gene$term))))
  cat(sprintf("    Genes: %d\n", length(unique(panglaodb_term2gene$gene))))
}

# ==============================================================================
# PanglaoDB 富集
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
# GO 富集
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
  cat(sprintf("[OK] Found %d significant terms\n", n_sig))

  saveRDS(go_enrich, file.path(OUTPUT_DIR, "reports", "go_enrich.rds"))
}

# ==============================================================================
# 可视化 - PanglaoDB
# ==============================================================================

cat("\n=== Visualizing PanglaoDB ===\n")

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

# ==============================================================================
# 可视化 - GO
# ==============================================================================

cat("\n=== Visualizing GO ===\n")

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
# 设置 API 密钥
Sys.setenv(DEEPSEEK_API_KEY = "sk-ed1879cf6fa14b04aac9cb6c078a3d05")
fanyi::set_translate_option(
  key = Sys.getenv("DEEPSEEK_API_KEY"),
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
  cat("[OK] API connection successful\n")
}
# ==============================================================================
# LLM 解释 - Task 1: Annotation
# ==============================================================================

cat("\n=== LLM Interpretation: Annotation ===\n")

annotation_results <- NULL

# 准备富集对象列表
enrich_list <- list()
if (!is.null(panglaodb_enrich)) {
  enrich_list[[length(enrich_list) + 1]] <- panglaodb_enrich
}
if (!is.null(go_enrich)) {
  enrich_list[[length(enrich_list) + 1]] <- go_enrich
}

if (length(enrich_list) > 0) {
  annotation_results <- tryCatch(
    {
      interpret(
        enrich_list,
        context = "B cells from normal nasal cavity, sinus, bronchi, and lung. Subclusters of Memory B, Naive B, Plasma cells.",
        task = "annotation"
      )
    },
    error = function(e) {
      cat("[WARN] Annotation failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(annotation_results)) {
    cat("\n[OK] Annotation complete! Processing results...\n\n")

    # --- 步骤 1: 打印摘要 (仅用于显示，不影响后续保存) ---
    for (nm in names(annotation_results)) {
      res <- annotation_results[[nm]]
      ct_print <- if (is.null(res$cell_type)) "Unknown" else res$cell_type
      conf_print <- if (is.null(res$confidence)) "None" else res$confidence

      cat(sprintf("%-35s: %s (%s)\n", nm, ct_print, conf_print))
    }

    # --- 步骤 2: 清洗数据并构建 DataFrame (核心修改) ---
    # 我们使用 lapply 手动提取，而不是让 bind_rows 去猜结构
    # 这样可以完美解决 "Can't recycle" 的长度不一致报错

    clean_data_list <- lapply(names(annotation_results), function(sub_name) {
      res <- annotation_results[[sub_name]]

      # 1. 安全获取 cell_type
      ct <- if (is.null(res$cell_type)) "Unknown" else res$cell_type

      # 2. 【核心逻辑】如果是 Unknown，直接返回 NULL (bind_rows 会自动忽略 NULL)
      if (ct == "Unknown") {
        return(NULL)
      }

      # 3. 安全获取 confidence
      conf <- if (is.null(res$confidence)) "NA" else res$confidence

      # 4. 返回一个标准的单行 dataframe
      return(data.frame(
        subcluster = sub_name,
        cell_type = ct,
        confidence = conf,
        stringsAsFactors = FALSE
      ))
    })

    # --- 步骤 3: 合并结果 ---
    # bind_rows 会自动忽略上面返回 NULL 的项
    annotation_df <- dplyr::bind_rows(clean_data_list)

    # --- 步骤 4: 保存 ---
    if (nrow(annotation_df) > 0) {
      write.csv(
        annotation_df,
        file.path(OUTPUT_DIR, "annotations.csv"),
        row.names = FALSE
      )
      cat("\n[OK] Saved annotations.csv (Filtered 'Unknown' clusters)\n")
    } else {
      cat(
        "\n[WARN] No valid annotations found after filtering 'Unknown'. CSV not saved.\n"
      )
    }

    saveRDS(
      annotation_results,
      file.path(OUTPUT_DIR, "reports", "annotation_results.rds")
    )
  }
}
library(dplyr)

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
        context = "B cells from normal nasal cavity, sinus, bronchi, and lung. Looking for baseline activation, antibody production, and proliferation states.",
        task = "phenotyping"
      )
    },
    error = function(e) {
      cat("[WARN] Phenotyping failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(phenotype_results)) {
    cat("\n[OK] Phenotyping complete! Processing results...\n\n")

    # --- 步骤 1: 打印摘要 ---
    for (nm in names(phenotype_results)) {
      res <- phenotype_results[[nm]]
      ph_print <- if (is.null(res$phenotype)) "Unknown" else res$phenotype
      conf_print <- if (is.null(res$confidence)) "None" else res$confidence

      cat(sprintf("%-35s: %s (%s)\n", nm, ph_print, conf_print))
    }

    # --- 步骤 2: 清洗数据并构建 DataFrame (防报错版) ---
    clean_pheno_list <- lapply(names(phenotype_results), function(sub_name) {
      res <- phenotype_results[[sub_name]]

      # 如果结果本身为空，跳过
      if (is.null(res)) {
        return(NULL)
      }

      # 安全获取字段
      ph <- if (is.null(res$phenotype)) "Unknown" else res$phenotype

      # 如果表型是 Unknown，根据需要决定是否跳过 (这里选择跳过以保持整洁)
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

    # --- 步骤 3: 合并结果 ---
    phenotype_df <- dplyr::bind_rows(clean_pheno_list)

    # --- 步骤 4: 保存 ---
    if (nrow(phenotype_df) > 0) {
      write.csv(
        phenotype_df,
        file.path(OUTPUT_DIR, "phenotypes.csv"),
        row.names = FALSE
      )
      cat("\n[OK] Saved phenotypes.csv (Filtered 'Unknown')\n")
    } else {
      cat(
        "\n[WARN] No valid phenotypes found after filtering. CSV not saved.\n"
      )
    }

    saveRDS(
      phenotype_results,
      file.path(OUTPUT_DIR, "reports", "phenotype_results.rds")
    )
  }
}

# ==============================================================================
# 按细胞类型分组分析
# ==============================================================================

cat("\n=== Per-Celltype Analysis ===\n")

# 获取 L2 细胞类型
celltypes_l2 <- unique(seurat_obj@meta.data$cell_type_L2)
celltypes_l2 <- celltypes_l2[!is.na(celltypes_l2)]

cat(sprintf("Analyzing %d cell types\n", length(celltypes_l2)))

celltype_interpretations <- list()

for (celltype in celltypes_l2) {
  cat(sprintf("\n--- %s ---\n", celltype))

  # 获取该细胞类型的所有L3亚群
  celltype_mask <- seurat_obj@meta.data$cell_type_L2 == celltype
  celltype_subclusters <- unique(seurat_obj@meta.data$cell_type_L3[
    celltype_mask
  ])
  celltype_subclusters <- celltype_subclusters[!is.na(celltype_subclusters)]

  n_subclusters <- length(celltype_subclusters)

  cat(sprintf("  Subclusters: %d\n", n_subclusters))

  # 跳过只有 1 个亚群的
  if (n_subclusters <= 1) {
    cat("  Skipping (only 1 subcluster)\n")
    next
  }

  # 过滤富集结果
  panglaodb_filtered <- NULL
  go_filtered <- NULL

  if (!is.null(panglaodb_enrich)) {
    # 确保 compareClusterResult 存在且是 data.frame
    if (!is.null(panglaodb_enrich@compareClusterResult)) {
      panglaodb_filtered <- panglaodb_enrich@compareClusterResult %>%
        filter(Cluster %in% celltype_subclusters)
    }
  }

  if (!is.null(go_enrich)) {
    if (!is.null(go_enrich@compareClusterResult)) {
      go_filtered <- go_enrich@compareClusterResult %>%
        filter(Cluster %in% celltype_subclusters)
    }
  }

  # 检查是否有结果
  has_results <- FALSE
  if (!is.null(panglaodb_filtered) && nrow(panglaodb_filtered) > 0) {
    has_results <- TRUE
  }
  if (!is.null(go_filtered) && nrow(go_filtered) > 0) {
    has_results <- TRUE
  }

  if (!has_results) {
    cat("  Skipping (no enrichment results)\n")
    next
  }

  # 创建子集富集对象
  enrich_subset <- list()

  if (!is.null(panglaodb_filtered) && nrow(panglaodb_filtered) > 0) {
    panglaodb_subset <- panglaodb_enrich
    panglaodb_subset@compareClusterResult <- panglaodb_filtered
    enrich_subset[[length(enrich_subset) + 1]] <- panglaodb_subset
  }

  if (!is.null(go_filtered) && nrow(go_filtered) > 0) {
    go_subset <- go_enrich
    go_subset@compareClusterResult <- go_filtered
    enrich_subset[[length(enrich_subset) + 1]] <- go_subset
  }

  # LLM 解释
  interpretation <- tryCatch(
    {
      interpret(
        enrich_subset,
        context = paste(
          celltype,
          "cells from normal nasal cavity, sinus, bronchi, and lung. Looking for functional states, homeostatic activation, and differentiation."
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
        cat("Generated:", format(Sys.time()), "\n\n")

        # 安全打印各个部分
        if (!is.null(interpretation$overview)) {
          cat("## Overview\n\n", interpretation$overview, "\n\n")
        }

        if (!is.null(interpretation$key_mechanisms)) {
          cat("## Key Mechanisms\n\n", interpretation$key_mechanisms, "\n\n")
        }

        if (!is.null(interpretation$narrative)) {
          cat("## Narrative\n\n", interpretation$narrative, "\n\n")
        }

        sink()
        cat("  [OK] Saved:", basename(report_file), "\n")
      },
      error = function(e) {
        tryCatch(sink(), error = function(e) NULL)
        cat("  [WARN] Failed to save report file\n")
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
    cat("**Database:** PanglaoDB & GO\n\n")
    cat("---\n\n")

    cat("## Dataset Summary\n\n")
    cat(sprintf("- Total cells: %d\n", ncol(seurat_obj)))
    cat(sprintf("- Subclusters: %d\n", length(unique(seurat_obj$cell_type_L3))))
    cat("\n---\n\n")

    # --- Annotation Section ---
    if (!is.null(annotation_results)) {
      cat("## Cell Subtype Annotations\n\n")

      for (nm in names(annotation_results)) {
        result <- annotation_results[[nm]]
        if (is.null(result)) {
          next
        }

        # 安全提取
        ct <- if (is.null(result$cell_type)) "Unknown" else result$cell_type
        conf <- if (is.null(result$confidence)) "NA" else result$confidence
        reason <- if (is.null(result$reasoning)) {
          "No reasoning provided."
        } else {
          result$reasoning
        }

        cat(sprintf("### %s\n\n", nm))
        cat(sprintf("**Cell Type:** %s  \n", ct))
        cat(sprintf("**Confidence:** %s  \n\n", conf))
        cat("**Reasoning:**\n\n", reason, "\n\n")
        cat("---\n\n")
      }
    }

    # --- Phenotype Section ---
    if (!is.null(phenotype_results)) {
      cat("## Functional Phenotypes\n\n")

      for (nm in names(phenotype_results)) {
        result <- phenotype_results[[nm]]
        if (is.null(result)) {
          next
        }

        ph <- if (is.null(result$phenotype)) "Unknown" else result$phenotype
        conf <- if (is.null(result$confidence)) "NA" else result$confidence

        cat(sprintf("### %s\n\n", nm))
        cat(sprintf("**Phenotype:** %s  \n", ph))
        cat(sprintf("**Confidence:** %s  \n\n", conf))
        cat("---\n\n")
      }
    }

    # --- Per-Celltype Interpretation Section ---
    if (length(celltype_interpretations) > 0) {
      cat("## Per-Celltype Interpretations\n\n")

      for (celltype in names(celltype_interpretations)) {
        interpretation <- celltype_interpretations[[celltype]]
        if (is.null(interpretation)) {
          next
        }

        cat(sprintf("### %s\n\n", celltype))

        # 优先打印 Narrative，如果没有则打印 Overview
        content <- if (!is.null(interpretation$narrative)) {
          interpretation$narrative
        } else if (!is.null(interpretation$overview)) {
          interpretation$overview
        } else {
          "No detailed interpretation available."
        }

        cat(content, "\n\n")
        cat("---\n\n")
      }
    }

    sink()
    cat("[OK] Report saved:", basename(report_file), "\n")
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
cat(
  "================================================================================\n"
)
cat("ANALYSIS COMPLETE\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
cat("  - REPORT.md\n")
cat("  - annotations.csv\n")
cat("  - phenotypes.csv\n")
cat("  - figures/panglaodb_dotplot.pdf\n")
cat("  - figures/go_dotplot.pdf\n")
cat("  - reports/[celltype]_interpretation.txt\n")
cat("\n")

cat("Load results:\n")
cat("  readRDS('reports/annotation_results.rds')\n")
cat("  readRDS('reports/phenotype_results.rds')\n")
cat("  readRDS('reports/celltype_interpretations.rds')\n")
cat("\n")

if (!is.null(annotation_results)) {
  cat("Confidence distribution:\n")
  conf_table <- table(sapply(annotation_results, function(x) x$confidence))
  for (conf in names(conf_table)) {
    cat(sprintf("  %s: %d\n", conf, conf_table[conf]))
  }
}

cat("\n")
cat(
  "================================================================================\n"
)
cat("DONE\n")
cat(
  "================================================================================\n"
)
