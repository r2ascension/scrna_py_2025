################################################################################
# 高效的基因标准化与可用性追踪系统（批次级别）
# Efficient Gene Standardization & Availability Tracking System (Batch-level)
#
# 核心改进：
# 1. try-catch兼容不同Seurat版本
# 2. 批次级别可用性矩阵（基因×批次，而非样本×基因）
# 3. 每个批次单独保存CSV
# 4. Checkpoint机制
# 5. 修复DietSeurat默认assay错误
################################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)
library(org.Hs.eg.db)
# library(org.Mm.eg.db)  # 小鼠数据
library(AnnotationDbi)

# ═══════════════════════════════════════════════════════════════════════════
# 辅助函数：兼容性封装（try-catch策略）
# ═══════════════════════════════════════════════════════════════════════════

#' 获取表达矩阵（兼容所有Seurat版本）
get_counts_matrix <- function(seurat_obj, assay = "RNA") {
  # 尝试不同的API调用方式
  counts <- tryCatch(
    {
      # 尝试1: Seurat V5 LayerData
      LayerData(seurat_obj, assay = assay, layer = "counts")
    },
    error = function(e1) {
      tryCatch(
        {
          # 尝试2: Seurat V5 GetAssayData with layer
          GetAssayData(seurat_obj, assay = assay, layer = "counts")
        },
        error = function(e2) {
          tryCatch(
            {
              # 尝试3: Seurat V4 GetAssayData with slot
              GetAssayData(seurat_obj, assay = assay, slot = "counts")
            },
            error = function(e3) {
              # 尝试4: 直接访问
              seurat_obj[[assay]]@counts
            }
          )
        }
      )
    }
  )

  return(counts)
}

#' 创建Assay对象（兼容所有Seurat版本）
create_assay_object <- function(counts_matrix) {
  assay <- tryCatch(
    {
      # 尝试1: Seurat V5
      CreateAssay5Object(counts = counts_matrix)
    },
    error = function(e) {
      # 尝试2: Seurat V4
      CreateAssayObject(counts = counts_matrix)
    }
  )

  return(assay)
}


# ═══════════════════════════════════════════════════════════════════════════
# 第一部分：基因名映射
# ═══════════════════════════════════════════════════════════════════════════

#' 将基因名映射到官方符号
map_gene_names <- function(genes, gene_db = org.Hs.eg.db) {
  gene_mapping <- data.frame(
    original_name = genes,
    official_symbol = NA_character_,
    mapping_source = NA_character_,
    stringsAsFactors = FALSE
  )

  # Layer 1: SYMBOL
  symbol_match <- tryCatch(
    {
      mapIds(
        gene_db,
        keys = genes,
        column = "SYMBOL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
    },
    error = function(e) rep(NA, length(genes))
  )

  matched <- !is.na(symbol_match)
  gene_mapping$official_symbol[matched] <- symbol_match[matched]
  gene_mapping$mapping_source[matched] <- "SYMBOL"

  # Layer 2: ALIAS
  unmatched <- is.na(gene_mapping$official_symbol)
  if (sum(unmatched) > 0) {
    alias_match <- tryCatch(
      {
        mapIds(
          gene_db,
          keys = genes[unmatched],
          column = "SYMBOL",
          keytype = "ALIAS",
          multiVals = "first"
        )
      },
      error = function(e) rep(NA, sum(unmatched))
    )

    alias_matched <- !is.na(alias_match)
    gene_mapping$official_symbol[unmatched][alias_matched] <- alias_match[
      alias_matched
    ]
    gene_mapping$mapping_source[unmatched][alias_matched] <- "ALIAS"
  }

  # Layer 3: ENSEMBL
  unmatched <- is.na(gene_mapping$official_symbol)
  if (sum(unmatched) > 0) {
    ensembl_match <- tryCatch(
      {
        mapIds(
          gene_db,
          keys = genes[unmatched],
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        )
      },
      error = function(e) rep(NA, sum(unmatched))
    )

    ensembl_matched <- !is.na(ensembl_match)
    gene_mapping$official_symbol[unmatched][ensembl_matched] <- ensembl_match[
      ensembl_matched
    ]
    gene_mapping$mapping_source[unmatched][ensembl_matched] <- "ENSEMBL"
  }

  # Layer 4: 保持原名
  unmatched <- is.na(gene_mapping$official_symbol)
  gene_mapping$official_symbol[unmatched] <- genes[unmatched]
  gene_mapping$mapping_source[unmatched] <- "UNMAPPED"

  return(gene_mapping)
}


# ═══════════════════════════════════════════════════════════════════════════
# 第二部分：处理同义基因（删除重复）
# ═══════════════════════════════════════════════════════════════════════════

#' 处理同义基因：只保留第一个，删除其他
remove_synonym_duplicates <- function(seurat_obj, gene_mapping, assay = "RNA") {
  # 识别同义基因组
  symbol_counts <- table(gene_mapping$official_symbol)
  duplicate_symbols <- names(symbol_counts)[symbol_counts > 1]

  if (length(duplicate_symbols) == 0) {
    # 无同义基因，只需重命名
    counts_matrix <- get_counts_matrix(seurat_obj, assay)
    rownames(counts_matrix) <- gene_mapping$official_symbol
    seurat_obj[[assay]] <- create_assay_object(counts_matrix)

    return(list(
      cleaned_obj = seurat_obj,
      synonym_report = NULL
    ))
  }

  cat("  发现", length(duplicate_symbols), "组同义基因\n")

  # 记录哪些基因被删除
  synonym_report <- list()
  genes_to_keep <- logical(nrow(seurat_obj))
  names(genes_to_keep) <- rownames(seurat_obj)
  genes_to_keep[] <- TRUE

  for (symbol in duplicate_symbols) {
    original_names <- gene_mapping$original_name[
      gene_mapping$official_symbol == symbol
    ]
    genes_to_keep[original_names[-1]] <- FALSE

    synonym_report[[symbol]] <- data.frame(
      official_symbol = symbol,
      n_synonyms = length(original_names),
      kept = original_names[1],
      removed = paste(original_names[-1], collapse = " | "),
      stringsAsFactors = FALSE
    )
  }

  # 过滤
  seurat_obj <- seurat_obj[genes_to_keep, ]

  # 重命名
  gene_mapping_kept <- gene_mapping[
    gene_mapping$original_name %in% rownames(seurat_obj),
  ]
  gene_mapping_kept <- gene_mapping_kept[
    match(rownames(seurat_obj), gene_mapping_kept$original_name),
  ]

  counts_matrix <- get_counts_matrix(seurat_obj, assay)
  rownames(counts_matrix) <- gene_mapping_kept$official_symbol
  seurat_obj[[assay]] <- create_assay_object(counts_matrix)

  synonym_report_df <- if (length(synonym_report) > 0) {
    do.call(rbind, synonym_report)
  } else {
    NULL
  }

  cat("  ✓ 删除", sum(!genes_to_keep), "个重复基因\n")

  return(list(
    cleaned_obj = seurat_obj,
    synonym_report = synonym_report_df
  ))
}


# ═══════════════════════════════════════════════════════════════════════════
# 第三部分：标准化单个RDS
# ═══════════════════════════════════════════════════════════════════════════

#' 标准化单个RDS文件
standardize_single_rds <- function(
  rds_path,
  gene_db = org.Hs.eg.db,
  output_dir = "standardized_output",
  assay = "RNA"
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  file_base <- tools::file_path_sans_ext(basename(rds_path))

  cat("\n========================================\n")
  cat("处理:", file_base, "\n")
  cat("========================================\n")

  # 读取并清理
  cat("[1/4] 读取并清理RDS...\n")
  seurat_obj <- readRDS(rds_path)
  cat("  原始大小:", format(object.size(seurat_obj), units = "MB"), "\n")

  # 检查目标assay是否存在 (兼容不同Seurat版本)
  available_assays <- tryCatch(
    {
      # 尝试1: 直接获取names
      names(seurat_obj@assays)
    },
    error = function(e) {
      tryCatch(
        {
          # 尝试2: 使用Assays函数
          as.character(Assays(seurat_obj))
        },
        error = function(e2) {
          # 尝试3: 使用names(Assays())
          names(Assays(seurat_obj))
        }
      )
    }
  )
  cat("  可用Assays:", paste(available_assays, collapse = ", "), "\n")

  if (!assay %in% available_assays) {
    stop(
      "指定的assay '",
      assay,
      "' 不存在。请从以下选择: ",
      paste(available_assays, collapse = ", ")
    )
  }

  # 关键修复: 设置默认assay为目标assay（避免DietSeurat错误）
  cat("  当前默认Assay:", DefaultAssay(seurat_obj), "\n")
  cat("  设置默认Assay为:", assay, "\n")
  DefaultAssay(seurat_obj) <- assay

  seurat_obj <- DietSeurat(
    seurat_obj,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL,
    assays = assay
  )

  cat("  清理后:", format(object.size(seurat_obj), units = "MB"), "\n")

  original_genes <- rownames(seurat_obj[[assay]])
  n_cells <- ncol(seurat_obj)
  cat("  细胞:", n_cells, ", 基因:", length(original_genes), "\n")

  # 基因映射
  cat("\n[2/4] 基因名映射...\n")
  gene_mapping <- map_gene_names(original_genes, gene_db)
  n_mapped <- sum(gene_mapping$mapping_source != "UNMAPPED")
  cat(
    "  映射:",
    n_mapped,
    "/",
    length(original_genes),
    sprintf("(%.1f%%)\n", 100 * n_mapped / length(original_genes))
  )

  # 处理同义基因
  cat("\n[3/4] 处理同义基因...\n")
  result <- remove_synonym_duplicates(seurat_obj, gene_mapping, assay)
  seurat_obj <- result$cleaned_obj

  cat("  原始:", length(original_genes), "→ 标准化:", nrow(seurat_obj), "\n")

  # 保存
  cat("\n[4/4] 保存结果...\n")

  cleaned_rds_dir <- file.path(output_dir, "cleaned_rds")
  dir.create(cleaned_rds_dir, showWarnings = FALSE, recursive = TRUE)
  cleaned_rds_path <- file.path(
    cleaned_rds_dir,
    paste0(file_base, "_cleaned.rds")
  )
  saveRDS(seurat_obj, cleaned_rds_path)
  cat("  ✓ RDS:", basename(cleaned_rds_path), "\n")

  # 保存映射和报告
  gene_mapping$file_source <- file_base
  mapping_dir <- file.path(output_dir, "gene_mappings")
  dir.create(mapping_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(
    gene_mapping,
    file.path(mapping_dir, paste0(file_base, "_mapping.csv"))
  )

  if (!is.null(result$synonym_report)) {
    synonym_dir <- file.path(output_dir, "synonym_reports")
    dir.create(synonym_dir, showWarnings = FALSE, recursive = TRUE)
    fwrite(
      result$synonym_report,
      file.path(synonym_dir, paste0(file_base, "_synonyms.csv"))
    )
  }

  final_result <- list(
    file_name = file_base,
    cleaned_rds = cleaned_rds_path,
    final_genes = rownames(seurat_obj),
    n_cells = n_cells
  )

  rm(seurat_obj, result)
  gc(verbose = FALSE)

  return(final_result)
}


# ═══════════════════════════════════════════════════════════════════════════
# 第四部分：构建批次级别基因可用性矩阵
# ═══════════════════════════════════════════════════════════════════════════

#' 构建批次级别基因可用性矩阵（每个批次一个CSV）
#'
#' @param cleaned_rds_dir 清理后RDS目录
#' @param batch_column metadata批次列名
#' @param output_dir 输出目录
#' @return 批次信息列表
#'
build_batch_gene_availability <- function(
  cleaned_rds_dir,
  batch_column = "orig.ident",
  output_dir = "batch_gene_availability"
) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║  构建批次级别基因可用性矩阵                           ║\n")
  cat("║  (每个批次单独保存CSV)                                ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  rds_files <- list.files(
    cleaned_rds_dir,
    pattern = "_cleaned\\.rds$",
    full.names = TRUE
  )
  cat("找到", length(rds_files), "个RDS文件\n\n")

  # 收集所有批次的基因信息
  all_batches <- list()
  all_genes <- c()

  cat("扫描批次和基因...\n")
  for (i in seq_along(rds_files)) {
    cat(sprintf("  [%d/%d] %s\n", i, length(rds_files), basename(rds_files[i])))

    seurat_obj <- readRDS(rds_files[i])

    # 确保使用正确的assay (兼容不同Seurat版本)
    available_assays <- tryCatch(
      {
        names(seurat_obj@assays)
      },
      error = function(e) {
        tryCatch(
          {
            as.character(Assays(seurat_obj))
          },
          error = function(e2) {
            names(Assays(seurat_obj))
          }
        )
      }
    )

    if (!"RNA" %in% available_assays) {
      cat("    ⚠️  警告: 没有RNA assay，使用默认assay\n")
    } else {
      DefaultAssay(seurat_obj) <- "RNA"
    }

    seurat_obj <- DietSeurat(
      seurat_obj,
      counts = TRUE,
      data = FALSE,
      scale.data = FALSE,
      dimreducs = NULL,
      graphs = NULL
    )

    genes <- rownames(seurat_obj)
    all_genes <- unique(c(all_genes, genes))

    # 获取批次信息
    if (
      !is.null(batch_column) && batch_column %in% colnames(seurat_obj@meta.data)
    ) {
      batch_ids <- unique(seurat_obj@meta.data[[batch_column]])

      # 对每个批次记录其基因列表
      for (batch_id in batch_ids) {
        if (!batch_id %in% names(all_batches)) {
          all_batches[[batch_id]] <- list(
            genes = genes,
            source_file = basename(rds_files[i])
          )
        }
      }
    } else {
      # 使用文件名作为批次
      batch_id <- gsub("_cleaned\\.rds$", "", basename(rds_files[i]))
      all_batches[[batch_id]] <- list(
        genes = genes,
        source_file = basename(rds_files[i])
      )
    }

    rm(seurat_obj)
    gc(verbose = FALSE)
  }

  cat("\n汇总:\n")
  cat("  总批次数:", length(all_batches), "\n")
  cat("  总基因数:", length(all_genes), "\n\n")

  # 为每个批次生成可用性向量并保存
  cat("生成批次级别可用性文件...\n")

  batch_summary <- data.frame(
    batch_id = character(),
    n_genes = integer(),
    source_file = character(),
    csv_file = character(),
    stringsAsFactors = FALSE
  )

  for (batch_id in names(all_batches)) {
    cat(sprintf("  处理批次: %s\n", batch_id))

    batch_genes <- all_batches[[batch_id]]$genes

    # 创建可用性向量：该批次有的基因=1，没有的=0
    availability <- ifelse(all_genes %in% batch_genes, 1, 0)

    # 保存为CSV
    availability_df <- data.frame(
      gene = all_genes,
      available = availability,
      stringsAsFactors = FALSE
    )

    safe_batch_id <- gsub("[^a-zA-Z0-9_-]", "_", batch_id)
    csv_file <- file.path(output_dir, paste0("batch_", safe_batch_id, ".csv"))
    fwrite(availability_df, csv_file)

    # 记录
    batch_summary <- rbind(
      batch_summary,
      data.frame(
        batch_id = batch_id,
        n_genes = sum(availability),
        source_file = all_batches[[batch_id]]$source_file,
        csv_file = basename(csv_file),
        stringsAsFactors = FALSE
      )
    )
  }

  # 保存批次汇总
  summary_file <- file.path(output_dir, "batch_summary.csv")
  fwrite(batch_summary, summary_file)

  # 保存完整基因列表
  gene_list_file <- file.path(output_dir, "all_genes.txt")
  writeLines(all_genes, gene_list_file)

  cat("\n✓ 完成!\n")
  cat("  批次数:", length(all_batches), "\n")
  cat("  汇总文件:", basename(summary_file), "\n")
  cat("  基因列表:", basename(gene_list_file), "\n\n")

  return(invisible(list(
    batches = all_batches,
    all_genes = all_genes,
    summary = batch_summary,
    output_dir = output_dir
  )))
}


# ═══════════════════════════════════════════════════════════════════════════
# 第五部分：智能merge
# ═══════════════════════════════════════════════════════════════════════════

#' 智能merge
smart_merge <- function(
  cleaned_rds_dir,
  add_cell_ids = TRUE,
  output_file = "merged_seurat.rds"
) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║  智能Merge                                             ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")

  rds_files <- list.files(
    cleaned_rds_dir,
    pattern = "_cleaned\\.rds$",
    full.names = TRUE
  )
  cat("找到", length(rds_files), "个RDS文件\n\n")

  # 读取对象
  cat("读取RDS文件...\n")
  obj_list <- lapply(seq_along(rds_files), function(i) {
    cat(sprintf("  [%d/%d] %s\n", i, length(rds_files), basename(rds_files[i])))
    obj <- readRDS(rds_files[i])
    # 确保默认assay正确 (兼容不同Seurat版本)
    available_assays <- tryCatch(
      {
        names(obj@assays)
      },
      error = function(e) {
        tryCatch(
          {
            as.character(Assays(obj))
          },
          error = function(e2) {
            names(Assays(obj))
          }
        )
      }
    )

    if ("RNA" %in% available_assays) {
      DefaultAssay(obj) <- "RNA"
    }
    obj
  })

  file_names <- gsub("_cleaned\\.rds$", "", basename(rds_files))

  # Merge
  cat("\n执行merge...\n")
  if (add_cell_ids) {
    merged <- merge(
      x = obj_list[[1]],
      y = obj_list[-1],
      add.cell.ids = file_names
    )
  } else {
    merged <- merge(x = obj_list[[1]], y = obj_list[-1])
  }

  assay_use <- DefaultAssay(merged)

  # Join all layers
  cat("  Join layers...\n")
  print(Layers(merged[[assay_use]]))
  merged[[assay_use]] <- JoinLayers(merged[[assay_use]])
  print(Layers(merged[[assay_use]]))

  # Get counts matrix
  cnt <- LayerData(
    object = merged,
    assay = assay_use,
    layer = "counts"
  )

  # Filter genes (min.cells = 3)
  cat("  过滤基因 (min.cells = 3)...\n")
  gene_ncells <- Matrix::rowSums(cnt > 0)
  keep_genes <- names(gene_ncells[gene_ncells >= 3])
  cat("  保留基因数:", length(keep_genes), "\n")

  merged <- subset(merged, features = keep_genes)

  cat("✓ merge完成\n")
  cat("  细胞数:", ncol(merged), "\n")
  cat("  基因数:", nrow(merged), "\n\n")

  # 保存
  cat("保存:", output_file, "\n")
  saveRDS(merged, output_file)
  cat("✓ 完成\n\n")

  return(merged)
}


# ═══════════════════════════════════════════════════════════════════════════
# 第六部分：完整工作流程（带Checkpoint）
# ═══════════════════════════════════════════════════════════════════════════

#' 完整工作流程（带Checkpoint）
#'
#' @param original_rds_dir 原始RDS目录
#' @param output_dir 输出目录
#' @param species 物种
#' @param batch_column metadata批次列名
#' @param assay Assay名称
#' @param start_from_checkpoint 从哪个checkpoint开始 ("none", "standardization", "availability")
#' @return 结果列表
#'
run_complete_workflow <- function(
  original_rds_dir,
  output_dir = "integration_output",
  species = "human",
  batch_column = "orig.ident",
  assay = "RNA",
  start_from_checkpoint = "none"
) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║  完整的基因标准化与智能整合工作流程                      ║\n")
  cat("║  - 批次级别可用性矩阵                                    ║\n")
  cat("║  - Checkpoint机制                                        ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  start_time <- Sys.time()
  cat("开始时间:", format(start_time), "\n")
  cat("物种:", species, "\n")
  cat("输入:", original_rds_dir, "\n")
  cat("输出:", output_dir, "\n")
  cat("Checkpoint:", start_from_checkpoint, "\n\n")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 选择数据库
  gene_db <- if (species == "human") {
    org.Hs.eg.db
  } else if (species == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("请安装: BiocManager::install('org.Mm.eg.db')")
    }
    library(org.Mm.eg.db)
    org.Mm.eg.db
  } else {
    stop("不支持的物种")
  }

  standardization_dir <- file.path(output_dir, "01_standardization")
  cleaned_rds_dir <- file.path(standardization_dir, "cleaned_rds")
  availability_dir <- file.path(output_dir, "02_batch_gene_availability")

  # ═══ Checkpoint 1: 标准化 ═══
  if (start_from_checkpoint == "none") {
    cat("═══════════════════════════════════════════════════════════\n")
    cat("阶段1: 基因名标准化\n")
    cat("═══════════════════════════════════════════════════════════\n")

    rds_files <- list.files(
      original_rds_dir,
      pattern = "\\.rds$",
      full.names = TRUE
    )
    cat("找到", length(rds_files), "个RDS文件\n\n")

    standardization_results <- lapply(seq_along(rds_files), function(i) {
      cat("\n进度:", i, "/", length(rds_files), "\n")
      result <- standardize_single_rds(
        rds_path = rds_files[i],
        gene_db = gene_db,
        output_dir = standardization_dir,
        assay = assay
      )
      gc(verbose = FALSE)
      result
    })

    # 保存checkpoint
    checkpoint_file <- file.path(output_dir, "checkpoint_standardization.rds")
    saveRDS(standardization_results, checkpoint_file)
    cat("\n✓ Checkpoint已保存:", basename(checkpoint_file), "\n")
  } else if (start_from_checkpoint == "standardization") {
    cat("从Checkpoint恢复: 跳过标准化阶段\n")
    checkpoint_file <- file.path(output_dir, "checkpoint_standardization.rds")
    if (file.exists(checkpoint_file)) {
      standardization_results <- readRDS(checkpoint_file)
      cat("✓ 已加载checkpoint\n")
    } else {
      stop("Checkpoint文件不存在:", checkpoint_file)
    }
  } else {
    standardization_results <- NULL
  }

  # ═══ Checkpoint 2: 批次可用性矩阵 ═══
  if (start_from_checkpoint %in% c("none", "standardization")) {
    cat("\n\n═══════════════════════════════════════════════════════════\n")
    cat("阶段2: 构建批次级别基因可用性矩阵\n")
    cat("═══════════════════════════════════════════════════════════\n")

    availability_result <- build_batch_gene_availability(
      cleaned_rds_dir = cleaned_rds_dir,
      batch_column = batch_column,
      output_dir = availability_dir
    )

    # 保存checkpoint
    checkpoint_file <- file.path(output_dir, "checkpoint_availability.rds")
    saveRDS(availability_result, checkpoint_file)
    cat("✓ Checkpoint已保存:", basename(checkpoint_file), "\n")
  } else if (start_from_checkpoint == "availability") {
    cat("从Checkpoint恢复: 跳过可用性矩阵构建\n")
    checkpoint_file <- file.path(output_dir, "checkpoint_availability.rds")
    if (file.exists(checkpoint_file)) {
      availability_result <- readRDS(checkpoint_file)
      cat("✓ 已加载checkpoint\n")
    } else {
      stop("Checkpoint文件不存在:", checkpoint_file)
    }
  } else {
    availability_result <- NULL
  }

  # ═══ 阶段3: 智能merge ═══
  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("阶段3: 智能Merge\n")
  cat("═══════════════════════════════════════════════════════════\n")

  merged_file <- file.path(output_dir, "03_merged_seurat.rds")

  merged_obj <- smart_merge(
    cleaned_rds_dir = cleaned_rds_dir,
    add_cell_ids = TRUE,
    output_file = merged_file
  )

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

  # ═══ 生成报告 ═══
  report_file <- file.path(output_dir, "REPORT.txt")
  sink(report_file)

  cat("═══════════════════════════════════════════════════════════\n")
  cat("基因标准化与智能整合工作流程报告\n")
  cat("═══════════════════════════════════════════════════════════\n\n")
  cat("完成时间:", format(end_time), "\n")
  cat("总耗时:", sprintf("%.2f", elapsed), "分钟\n")
  cat("物种:", species, "\n\n")

  cat("───────────────────────────────────────────────────────────\n")
  cat("核心改进:\n")
  cat("───────────────────────────────────────────────────────────\n")
  cat("1. try-catch兼容不同Seurat版本\n")
  cat("2. 批次级别可用性矩阵（而非样本级别）\n")
  cat("3. 每个批次单独保存CSV（节省内存）\n")
  cat("4. Checkpoint机制（可中断恢复）\n")
  cat("5. 修复DietSeurat默认assay错误\n\n")

  if (!is.null(standardization_results)) {
    cat("───────────────────────────────────────────────────────────\n")
    cat("阶段1: 基因标准化\n")
    cat("───────────────────────────────────────────────────────────\n")
    for (res in standardization_results) {
      cat(sprintf(
        "%s: %d 基因, %d 细胞\n",
        res$file_name,
        length(res$final_genes),
        res$n_cells
      ))
    }
  }

  if (!is.null(availability_result)) {
    cat("\n───────────────────────────────────────────────────────────\n")
    cat("阶段2: 批次级别基因可用性\n")
    cat("───────────────────────────────────────────────────────────\n")
    cat("批次数:", length(availability_result$batches), "\n")
    cat("总基因数:", length(availability_result$all_genes), "\n")
    cat("输出目录:", basename(availability_result$output_dir), "\n")
  }

  cat("\n───────────────────────────────────────────────────────────\n")
  cat("阶段3: Merge\n")
  cat("───────────────────────────────────────────────────────────\n")
  cat("文件:", basename(merged_file), "\n")
  cat("总细胞:", ncol(merged_obj), "\n")
  cat("总基因:", nrow(merged_obj), "\n")

  cat("\n═══════════════════════════════════════════════════════════\n")
  cat("输出文件:\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("01_standardization/\n")
  cat("  ├── cleaned_rds/             # 标准化后的RDS\n")
  cat("  ├── gene_mappings/           # 基因映射表\n")
  cat("  └── synonym_reports/         # 同义基因删除报告\n")
  cat("02_batch_gene_availability/    # 批次级别可用性矩阵\n")
  cat("  ├── batch_*.csv              # 每个批次一个CSV\n")
  cat("  ├── batch_summary.csv        # 批次汇总\n")
  cat("  └── all_genes.txt            # 完整基因列表\n")
  cat("03_merged_seurat.rds           # merge后的对象\n")
  cat("checkpoint_*.rds               # Checkpoint文件\n")
  cat("REPORT.txt                     # 本报告\n")

  sink()

  cat("\n报告已保存:", report_file, "\n")
  cat("\n╔══════════════════════════════════════════════════════════╗\n")
  cat(
    "║  完成！总耗时:",
    sprintf("%.2f", elapsed),
    "分钟",
    rep(" ", max(0, 28 - nchar(sprintf("%.2f", elapsed)))),
    "║\n"
  )
  cat("╚══════════════════════════════════════════════════════════╝\n\n")

  return(list(
    standardization = standardization_results,
    availability = availability_result,
    merged_object = merged_obj,
    output_dir = output_dir,
    elapsed_time = elapsed
  ))
}


# ═══════════════════════════════════════════════════════════════════════════
# 批次NA值诊断工具
# Diagnosis Tool for NA Batch Values
# ═══════════════════════════════════════════════════════════════════════════

#' 诊断所有RDS文件中的批次信息
#' @param cleaned_rds_dir 清理后RDS目录
#' @param batch_column metadata批次列名
diagnose_batch_issues <- function(cleaned_rds_dir, batch_column = "dataset") {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║  批次信息诊断                                          ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")

  rds_files <- list.files(
    cleaned_rds_dir,
    pattern = "_cleaned\\.rds$",
    full.names = TRUE
  )

  cat("找到", length(rds_files), "个RDS文件\n\n")

  diagnosis_results <- list()

  for (i in seq_along(rds_files)) {
    file_name <- basename(rds_files[i])
    cat("═══════════════════════════════════════════════════════\n")
    cat(sprintf("[%d/%d] 检查: %s\n", i, length(rds_files), file_name))
    cat("═══════════════════════════════════════════════════════\n")

    # 读取对象
    seurat_obj <- readRDS(rds_files[i])

    # 确保默认assay正确 (兼容不同Seurat版本)
    available_assays <- tryCatch(
      {
        names(seurat_obj@assays)
      },
      error = function(e) {
        tryCatch(
          {
            as.character(Assays(seurat_obj))
          },
          error = function(e2) {
            names(Assays(seurat_obj))
          }
        )
      }
    )

    if ("RNA" %in% available_assays) {
      DefaultAssay(seurat_obj) <- "RNA"
    }

    seurat_obj <- DietSeurat(
      seurat_obj,
      counts = TRUE,
      data = FALSE,
      scale.data = FALSE,
      dimreducs = NULL,
      graphs = NULL
    )

    # 基本信息
    cat("细胞数:", ncol(seurat_obj), "\n")
    cat("基因数:", nrow(seurat_obj), "\n")

    # 检查metadata列
    cat("\nMetadata列:\n")
    print(colnames(seurat_obj@meta.data))

    # 检查批次列是否存在
    if (!batch_column %in% colnames(seurat_obj@meta.data)) {
      cat("\n⚠️  警告: 批次列", batch_column, "不存在!\n")
      cat(
        "可用的列:",
        paste(colnames(seurat_obj@meta.data), collapse = ", "),
        "\n"
      )

      diagnosis_results[[file_name]] <- list(
        status = "missing_column",
        available_columns = colnames(seurat_obj@meta.data)
      )
    } else {
      # 检查批次值
      batch_values <- seurat_obj@meta.data[[batch_column]]
      unique_batches <- unique(batch_values)

      cat("\n批次列:", batch_column, "\n")
      cat("唯一批次数:", length(unique_batches), "\n")

      # 检查NA值
      na_count <- sum(is.na(batch_values))
      if (na_count > 0) {
        cat("\n⚠️  发现NA值:\n")
        cat("  NA细胞数:", na_count, "\n")
        cat(
          "  NA比例:",
          sprintf("%.2f%%", 100 * na_count / length(batch_values)),
          "\n"
        )

        # 显示前几个NA细胞的信息
        na_cells <- which(is.na(batch_values))
        cat("\n前10个NA细胞示例:\n")
        print(head(seurat_obj@meta.data[na_cells, ], 10))
      }

      # 显示所有批次值（包括NA）
      cat("\n批次值统计:\n")
      batch_table <- table(batch_values, useNA = "ifany")
      print(batch_table)

      diagnosis_results[[file_name]] <- list(
        status = "ok",
        batch_column = batch_column,
        unique_batches = unique_batches,
        na_count = na_count,
        na_percentage = 100 * na_count / length(batch_values),
        batch_distribution = as.data.frame(batch_table)
      )
    }

    cat("\n")
    rm(seurat_obj)
    gc(verbose = FALSE)
  }

  # 汇总报告
  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║  诊断汇总                                              ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")

  files_with_missing_column <- sum(sapply(diagnosis_results, function(x) {
    x$status == "missing_column"
  }))
  files_with_na <- sum(sapply(diagnosis_results, function(x) {
    x$status == "ok" && x$na_count > 0
  }))

  cat("文件总数:", length(diagnosis_results), "\n")
  cat("缺少批次列的文件:", files_with_missing_column, "\n")
  cat("含有NA批次值的文件:", files_with_na, "\n\n")

  if (files_with_na > 0) {
    cat("含有NA的文件详情:\n")
    for (fname in names(diagnosis_results)) {
      res <- diagnosis_results[[fname]]
      if (res$status == "ok" && res$na_count > 0) {
        cat(sprintf(
          "  - %s: %d个NA细胞 (%.2f%%)\n",
          fname,
          res$na_count,
          res$na_percentage
        ))
      }
    }
  }

  return(invisible(diagnosis_results))
}


#' 修复NA批次值
#' @param seurat_obj Seurat对象
#' @param batch_column 批次列名
#' @param file_name 文件名(用于生成替代批次ID)
#' @param na_strategy 处理NA的策略: "use_filename", "use_default", "remove_cells"
fix_na_batches <- function(
  seurat_obj,
  batch_column = "dataset",
  file_name = NULL,
  na_strategy = "use_filename"
) {
  if (!batch_column %in% colnames(seurat_obj@meta.data)) {
    stop("批次列不存在:", batch_column)
  }

  batch_values <- seurat_obj@meta.data[[batch_column]]
  na_indices <- is.na(batch_values)
  na_count <- sum(na_indices)

  if (na_count == 0) {
    cat("✓ 无NA值,无需修复\n")
    return(seurat_obj)
  }

  cat("发现", na_count, "个NA批次值\n")

  if (na_strategy == "use_filename") {
    # 使用文件名替代NA
    if (is.null(file_name)) {
      stop("na_strategy='use_filename'时必须提供file_name参数")
    }

    replacement_value <- gsub("_cleaned\\.rds$", "", file_name)
    seurat_obj@meta.data[[batch_column]][na_indices] <- replacement_value
    cat("✓ 已将NA替换为:", replacement_value, "\n")
  } else if (na_strategy == "use_default") {
    # 使用固定默认值
    seurat_obj@meta.data[[batch_column]][na_indices] <- "Unknown_Batch"
    cat("✓ 已将NA替换为: Unknown_Batch\n")
  } else if (na_strategy == "remove_cells") {
    # 删除NA细胞
    seurat_obj <- seurat_obj[, !na_indices]
    cat("✓ 已删除", na_count, "个NA细胞\n")
    cat("剩余细胞:", ncol(seurat_obj), "\n")
  } else {
    stop("不支持的na_strategy:", na_strategy)
  }

  return(seurat_obj)
}


#' 批量修复所有RDS文件中的NA批次值
batch_fix_na_issues <- function(
  cleaned_rds_dir,
  batch_column = "dataset",
  na_strategy = "use_filename",
  output_dir = NULL
) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║  批量修复NA批次值                                      ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")

  if (is.null(output_dir)) {
    output_dir <- cleaned_rds_dir
    cat("⚠️  将覆盖原文件\n\n")
  } else {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    cat("输出目录:", output_dir, "\n\n")
  }

  rds_files <- list.files(
    cleaned_rds_dir,
    pattern = "_cleaned\\.rds$",
    full.names = TRUE
  )

  cat("找到", length(rds_files), "个RDS文件\n\n")

  fix_summary <- data.frame(
    file_name = character(),
    original_na_count = integer(),
    action_taken = character(),
    final_cells = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(rds_files)) {
    file_name <- basename(rds_files[i])
    cat(sprintf("[%d/%d] 处理: %s\n", i, length(rds_files), file_name))

    seurat_obj <- readRDS(rds_files[i])

    # 确保默认assay正确 (兼容不同Seurat版本)
    available_assays <- tryCatch(
      {
        names(seurat_obj@assays)
      },
      error = function(e) {
        tryCatch(
          {
            as.character(Assays(seurat_obj))
          },
          error = function(e2) {
            names(Assays(seurat_obj))
          }
        )
      }
    )

    if ("RNA" %in% available_assays) {
      DefaultAssay(seurat_obj) <- "RNA"
    }

    original_cells <- ncol(seurat_obj)

    # 检查是否有NA
    if (batch_column %in% colnames(seurat_obj@meta.data)) {
      na_count <- sum(is.na(seurat_obj@meta.data[[batch_column]]))

      if (na_count > 0) {
        cat("  发现", na_count, "个NA值,执行修复...\n")
        seurat_obj <- fix_na_batches(
          seurat_obj,
          batch_column,
          file_name,
          na_strategy
        )
        action <- paste0("fixed_", na_strategy)
      } else {
        cat("  ✓ 无NA值\n")
        action <- "no_action_needed"
      }
    } else {
      cat("  ⚠️  批次列不存在\n")
      na_count <- NA
      action <- "column_missing"
    }

    # 保存
    output_path <- file.path(output_dir, file_name)
    saveRDS(seurat_obj, output_path)

    fix_summary <- rbind(
      fix_summary,
      data.frame(
        file_name = file_name,
        original_na_count = na_count,
        action_taken = action,
        final_cells = ncol(seurat_obj),
        stringsAsFactors = FALSE
      )
    )

    rm(seurat_obj)
    gc(verbose = FALSE)
    cat("\n")
  }

  # 保存修复报告
  summary_file <- file.path(output_dir, "na_fix_summary.csv")
  write.csv(fix_summary, summary_file, row.names = FALSE)

  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║  修复完成                                              ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")
  cat("修复报告:", summary_file, "\n\n")

  return(invisible(fix_summary))
}


# ═══════════════════════════════════════════════════════════════════════════
# 使用示例
# ═══════════════════════════════════════════════════════════════════════════

# # 设置工作目录
setwd("/home/h2048/data/R/1124")

# 设置Seurat版本
options(Seurat.object.assay.version = "v4")

# 完整流程（从头开始）
results <- run_complete_workflow(
  original_rds_dir = "/home/h2048/data/source/final",
  output_dir = "/home/h2048/data/R/1124/merge2",
  species = "human",
  batch_column = "study",
  assay = "RNA",
  start_from_checkpoint = "none"
)
#
# # 从checkpoint恢复（跳过标准化）
# results <- run_complete_workflow(
#   original_rds_dir = "/home/h2048/data/source/final",
#   output_dir = "/home/h2048/data/R/1124/merge2",
#   species = "human",
#   batch_column = "study",
#   assay = "RNA",
#   start_from_checkpoint = "standardization"
# )
#
# # 批次诊断
# diagnosis <- diagnose_batch_issues(
#   cleaned_rds_dir = "/home/h2048/data/R/1124/merge2/01_standardization/cleaned_rds",
#   batch_column = "study"
# )
#
# # 批量修复NA值
# fix_summary <- batch_fix_na_issues(
#   cleaned_rds_dir = "/home/h2048/data/R/1124/merge2/01_standardization/cleaned_rds",
#   batch_column = "study",
#   na_strategy = "use_filename",
#   output_dir = "/home/h2048/data/R/1124/merge2/01_standardization/cleaned_rds_fixed"
# )

merged_object <- results$merged_object
head(merged_object)
merged_object$percent_sarscov2 <- NULL
# 假设数据在Seurat对象中
merged_object@meta.data$ann_level_2[is.na(
  merged_object@meta.data$ann_level_2
)] <- merged_object@meta.data$CellType[is.na(
  merged_object@meta.data$ann_level_2
)]
merged_object@meta.data$batch[is.na(
  merged_object@meta.data$batch
)] <- merged_object@meta.data$dataset[is.na(merged_object@meta.data$batch)]

# 检查结果
cat(sprintf("study空值数: %d\n", sum(is.na(merged_object@meta.data$study))))
cat(sprintf("batch空值数: %d\n", sum(is.na(merged_object@meta.data$batch))))
table(merged_object@meta.data$entropy_original_ann_level_2_clean_leiden_3)
table(merged_object$ann_level_1, merged_object$cellType)
table(merged_object$CellType, merged_object$dataset)
getwd()
library(Seurat)
library(SCNT)
library(reticulate)

# 1）确认 Python 环境里有 anndata
py_config()
# 或者指定环境：
# use_python("/path/to/python", required = TRUE)
use_condaenv("bbknn_env", required = TRUE)

saveRDS(merged_object, 'merged_object_sc.rds')

# 在导出前，检查并转换所有metadata列为兼容类型
meta <- merged_object@meta.data

# 找出所有非标准类型的列
for (col in colnames(meta)) {
  col_class <- class(meta[[col]])[1]

  # 转换所有非character/numeric/logical的列为character
  if (!col_class %in% c("character", "numeric", "integer", "logical")) {
    cat(sprintf("Converting %s (%s) to character\n", col, col_class))
    meta[[col]] <- as.character(meta[[col]])
  }

  # 特别处理NA值（确保不是R的NA类型）
  if (any(is.na(meta[[col]]))) {
    if (is.character(meta[[col]])) {
      meta[[col]][is.na(meta[[col]])] <- "Unknown"
    }
  }
}

# 写回Seurat对象
merged_object@meta.data <- meta

# 2）如果是由多个对象 merge 出来的 Seurat v5 对象，先 JoinLayers（很重要）
assay_use <- DefaultAssay(merged_object)
merged_object[[assay_use]] <- JoinLayers(merged_object[[assay_use]])

# 3）直接导出为 h5ad（单细胞模式）
GetH5ad(
  merged_object,
  output_path = "merged_object_sc.h5ad",
  mode = "sc", # 单细胞 / 非空间
  assay = assay_use # 比如 "RNA"
)
