#!/usr/bin/env Rscript
# ==============================================================================
# 表型→基因关联 路线D：模块评分 → 加权 HPO 词条映射 → Phen2Gene 融合 → 基因排序
# Phenotype→Gene Pipeline D: Module Score → Weighted HPO → Phen2Gene Fusion
# ==============================================================================
#
# 工作流摘要 (Workflow Summary):
#   1. 将 HPO gene sets（或 MSigDB HPO gene sets）作为"表型程序"
#   2. 用 AddModuleScore（Seurat）/ UCell / AUCell 对每个细胞计算表型程序活性分数
#   3. 按细胞类型汇总分数 → 得到每个细胞类型对每个 HPO term 的程序活性强度
#   4. 取各细胞类型中活性最强的 top-K HPO terms → 加权 HPO term 集合
#   5. 选项 D5a: 直接将 top HPO IDs 作为查询，调用 Phen2Gene REST API → 基因排序
#   6. 选项 D5b: 用 HPO→gene 反向投影（与路线A一致）→ 基因排序
#   7. 融合：将 scRNA 模块分数作为额外证据，对步骤 5 基因排序做二次加权
#
# 输入 (Inputs):
#   seurat_obj    : Seurat 对象（含 raw/normalized counts、细胞类型注释）
#   hpo_gene_path : HPO phenotype_to_genes.txt 本地路径
#   target_hpo_ids: 目标 HPO term ID 列表（可选；None=自动从模块分数推断）
#
# 输出 (Outputs):
#   <output_dir>/module_scores/         : 模块分数矩阵（细胞×HPO term）
#   <output_dir>/celltype_hpo_scores/   : 细胞类型×HPO term 汇总分数
#   <output_dir>/gene_ranking_D.csv     : 融合排序结果
#
# 依赖包 (Required Packages):
#   Seurat, UCell（可选）, dplyr, ggplot2
#   httr, jsonlite（Phen2Gene API 可选）
#
# 参考 (References):
#   - Seurat AddModuleScore:
#     https://satijalab.org/seurat/reference/addmodulescore
#   - UCell (Andreatta & Carmona, 2021):
#     https://github.com/carmonalab/UCell
#   - AUCell (Aibar et al., 2017):
#     https://bioconductor.org/packages/AUCell
#   - Phen2Gene (Zhao et al., 2020):
#     https://github.com/WGLab/Phen2Gene
#
# 版本 (Version): v1.0  2026-03-27
# ==============================================================================

# ==============================================================================
# 0. 配置参数 Configuration
# ==============================================================================

SEURAT_RDS_PATH  <- NULL   # 替换为你的 RDS 路径
HPO_GENE_PATH    <- NULL   # HPO phenotype_to_genes.txt 路径（NULL=下载）
OUTPUT_DIR       <- "./results_pipeline_D"

CELL_TYPE_COL    <- "Annotation_2"
SAMPLE_COL       <- "sample"
CONDITION_COL    <- "disease"

# 目标 HPO terms（NULL=从模块分数自动推断 top-K terms）
TARGET_HPO_IDS   <- NULL

# 模块评分参数
SCORE_METHOD     <- "AddModuleScore"  # "AddModuleScore" | "UCell" | "AUCell"
MAX_HPO_TERMS    <- 200   # 计算模块分数的 HPO terms 最大数量（太多会很慢）
TOP_HPO_PER_CELLTYPE <- 20  # 每个细胞类型保留 top HPO terms 数

# Phen2Gene API
USE_PHEN2GENE_API <- FALSE
PHEN2GENE_API_URL <- "https://phen2gene.wglab.org/api/query"
PHEN2GENE_TOP_K   <- 200

# 基因排序参数
FINAL_TOP_K <- 200  # 最终返回 top-k 基因数

# ==============================================================================
# 1. 加载依赖包
# ==============================================================================

.load_packages_D <- function() {
  pkgs <- c("Seurat", "dplyr", "ggplot2")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("请安装 %s", pkg))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  if (SCORE_METHOD == "UCell") {
    if (!requireNamespace("UCell", quietly = TRUE)) {
      message("[WARN] UCell 未安装，回退到 AddModuleScore")
      SCORE_METHOD <<- "AddModuleScore"
    } else {
      suppressPackageStartupMessages(library(UCell))
    }
  }
  if (SCORE_METHOD == "AUCell") {
    if (!requireNamespace("AUCell", quietly = TRUE)) {
      message("[WARN] AUCell 未安装，回退到 AddModuleScore")
      SCORE_METHOD <<- "AddModuleScore"
    } else {
      suppressPackageStartupMessages(library(AUCell))
    }
  }
  if (USE_PHEN2GENE_API) {
    for (pkg in c("httr", "jsonlite")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("[WARN] %s 未安装，跳过 Phen2Gene API", pkg))
        USE_PHEN2GENE_API <<- FALSE
      }
    }
  }
}

# ==============================================================================
# 2. HPO 基因集加载
# ==============================================================================

#' 加载 HPO phenotype_to_genes.txt
#'
#' @param local_path 本地文件路径（NULL=下载）
#' @param max_terms  最多加载的 HPO term 数量（避免计算太多模块分数）
#' @return list(hpo_df, hpo_gene_sets)
load_hpo_for_module_score <- function(local_path = NULL,
                                       max_terms = MAX_HPO_TERMS) {
  message("[Step 2] 加载 HPO 基因集...")

  if (!is.null(local_path) && file.exists(local_path)) {
    hpo_raw <- read.table(local_path, sep = "\t", header = FALSE,
                          quote = "", comment.char = "#", fill = TRUE)
    colnames(hpo_raw)[c(1, 4)] <- c("hpo_id", "gene_symbol")
    hpo_df <- hpo_raw[, c("hpo_id", "gene_symbol")]
  } else {
    message("  尝试下载 HPO phenotype_to_genes.txt...")
    url <- paste0("https://purl.obolibrary.org/obo/hp/hpoa/",
                  "phenotype_to_genes.txt")
    tryCatch({
      tmp <- tempfile(fileext = ".txt")
      download.file(url, destfile = tmp, quiet = TRUE, method = "libcurl")
      hpo_raw <- read.table(tmp, sep = "\t", header = FALSE,
                             quote = "", fill = TRUE)
      colnames(hpo_raw)[c(1, 4)] <- c("hpo_id", "gene_symbol")
      hpo_df <- hpo_raw[, c("hpo_id", "gene_symbol")]
    }, error = function(e) {
      message("  HPO 下载失败，使用 Toy HPO")
      hpo_df <<- .make_toy_hpo_df()
    })
  }

  hpo_df$gene_symbol <- toupper(hpo_df$gene_symbol)
  hpo_df <- hpo_df[!is.na(hpo_df$hpo_id) & !is.na(hpo_df$gene_symbol), ]

  # 按 HPO term 拆分为基因集列表（格式：list(hpo_id = c(gene1, gene2, ...))）
  hpo_gene_sets <- split(hpo_df$gene_symbol, hpo_df$hpo_id)

  # 过滤：基因集大小在 5–500 之间
  gs_sizes <- lengths(hpo_gene_sets)
  hpo_gene_sets <- hpo_gene_sets[gs_sizes >= 5 & gs_sizes <= 500]

  # 如果 term 太多，只取前 max_terms 个（可按基因集大小排序）
  if (length(hpo_gene_sets) > max_terms) {
    message(sprintf("  HPO terms 过多（%d），截取前 %d 个",
                    length(hpo_gene_sets), max_terms))
    hpo_gene_sets <- hpo_gene_sets[seq_len(max_terms)]
  }

  message(sprintf("  保留 %d 个 HPO terms 用于模块评分", length(hpo_gene_sets)))
  list(hpo_df = hpo_df, hpo_gene_sets = hpo_gene_sets)
}

# ==============================================================================
# 3. 模块评分：AddModuleScore / UCell / AUCell
# ==============================================================================

#' 计算 HPO 表型程序模块分数
#'
#' @param seurat_obj Seurat 对象（需已完成 NormalizeData）
#' @param hpo_gene_sets list，HPO term → gene 向量
#' @param method      "AddModuleScore" | "UCell" | "AUCell"
#' @param output_dir  输出目录
#' @return seurat_obj（meta.data 新增模块分数列）
compute_module_scores <- function(seurat_obj,
                                   hpo_gene_sets,
                                   method = SCORE_METHOD,
                                   output_dir = file.path(OUTPUT_DIR, "module_scores")) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message(sprintf("[Step 3] 用 %s 计算 %d 个 HPO 表型程序的模块分数...",
                  method, length(hpo_gene_sets)))

  # 确保数据已归一化
  if (!("data" %in% names(seurat_obj@assays$RNA@layers)) &&
      !all(is.finite(seurat_obj@assays$RNA@data@x))) {
    message("  数据未归一化，执行 NormalizeData...")
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  }

  # 只保留在 Seurat 对象中存在的基因
  all_seurat_genes <- toupper(rownames(seurat_obj))
  hpo_gene_sets <- lapply(hpo_gene_sets, function(gs) {
    gs[gs %in% all_seurat_genes]
  })
  hpo_gene_sets <- hpo_gene_sets[lengths(hpo_gene_sets) >= 3]

  if (method == "UCell") {
    # UCell：基于 Mann-Whitney U 统计量，对细胞单独评分
    seurat_obj <- UCell::AddModuleScore_UCell(
      seurat_obj,
      features = hpo_gene_sets,
      name     = "",    # 不添加后缀
      verbose  = FALSE
    )
  } else if (method == "AUCell") {
    # AUCell：为每个细胞计算基因集的 AUC 活性分数
    expr_mat <- tryCatch(
      GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
    )
    cells_rankings <- AUCell::AUCell_buildRankings(
      as.matrix(expr_mat), plotStats = FALSE, verbose = FALSE
    )
    cells_AUC <- AUCell::AUCell_calcAUC(
      hpo_gene_sets, cells_rankings, verbose = FALSE
    )
    # 将 AUC 写入 meta.data
    auc_mat <- t(as.matrix(AUCell::getAUC(cells_AUC)))
    colnames(auc_mat) <- paste0(colnames(auc_mat), "_AUC")
    seurat_obj <- AddMetaData(seurat_obj, metadata = as.data.frame(auc_mat))
  } else {
    # AddModuleScore（Seurat 内置，Tirosh et al. 方法）
    # 由于 AddModuleScore 一次只能添加一个 list，批量循环
    batch_size <- 50  # 每批计算 50 个 terms，避免 Seurat 名称截断
    hpo_names  <- names(hpo_gene_sets)
    for (batch_start in seq(1, length(hpo_gene_sets), by = batch_size)) {
      batch_end  <- min(batch_start + batch_size - 1, length(hpo_gene_sets))
      batch_sets <- hpo_gene_sets[batch_start:batch_end]
      # Seurat AddModuleScore 需要 feature list，并用 name 参数作为前缀
      score_names <- paste0("HPO_", gsub(":", "_", names(batch_sets)))
      seurat_obj <- tryCatch(
        AddModuleScore(
          seurat_obj,
          features = batch_sets,
          name     = "HPO_score_",
          ctrl     = 100,
          seed     = 42
        ),
        error = function(e) {
          message(sprintf("  AddModuleScore 批次 %d-%d 失败: %s",
                          batch_start, batch_end, e$message))
          seurat_obj  # 返回未修改的对象
        }
      )
    }
  }

  message(sprintf("  模块分数计算完成，meta.data 列数: %d",
                  ncol(seurat_obj@meta.data)))

  # 保存每个细胞的模块分数
  score_cols <- grep("HPO_score_", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(score_cols) > 0) {
    score_df <- seurat_obj@meta.data[, score_cols, drop = FALSE]
    score_df$cell_id   <- rownames(score_df)
    score_df$cell_type <- seurat_obj@meta.data[[CELL_TYPE_COL]]
    write.csv(score_df, file.path(output_dir, "cell_module_scores.csv"),
              row.names = FALSE)
    message(sprintf("  细胞模块分数保存至: %s/cell_module_scores.csv", output_dir))
  }

  seurat_obj
}

# ==============================================================================
# 4. 按细胞类型汇总模块分数 → 加权 HPO term 集合
# ==============================================================================

#' 按细胞类型汇总模块分数，得到每个细胞类型的 HPO term 活性排名
#'
#' @param seurat_obj 已添加模块分数的 Seurat 对象
#' @param cell_type_col 细胞类型列名
#' @param hpo_ids   HPO term ID 列表（用于映射 score 列名）
#' @param top_n     每个细胞类型保留 top_n HPO terms
#' @param output_dir 输出目录
#' @return data.frame：cell_type, hpo_id, mean_score, rank
aggregate_celltype_hpo_scores <- function(seurat_obj,
                                           cell_type_col = CELL_TYPE_COL,
                                           hpo_ids,
                                           top_n = TOP_HPO_PER_CELLTYPE,
                                           output_dir = file.path(OUTPUT_DIR, "celltype_hpo_scores")) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("[Step 4] 按细胞类型汇总 HPO 模块分数...")

  score_cols <- grep("HPO_score_", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(score_cols) == 0) {
    message("  [WARN] 未找到模块分数列（HPO_score_*），跳过汇总")
    return(data.frame())
  }

  meta <- seurat_obj@meta.data
  cell_types <- unique(meta[[cell_type_col]])

  all_scores <- lapply(cell_types, function(ct) {
    ct_mask  <- meta[[cell_type_col]] == ct
    ct_scores <- colMeans(meta[ct_mask, score_cols, drop = FALSE], na.rm = TRUE)
    df <- data.frame(
      cell_type  = ct,
      score_col  = names(ct_scores),
      mean_score = as.numeric(ct_scores),
      stringsAsFactors = FALSE
    )
    df <- df[order(df$mean_score, decreasing = TRUE), ]
    # 映射 score_col → hpo_id
    df$hpo_id <- gsub("^HPO_score_", "", df$score_col)
    df$hpo_id <- gsub("_", ":", df$hpo_id)
    utils::head(df, top_n)
  })

  ct_hpo_df <- do.call(rbind, all_scores)
  write.csv(ct_hpo_df, file.path(output_dir, "celltype_hpo_scores.csv"),
            row.names = FALSE)
  message(sprintf("  细胞类型×HPO 汇总保存至: %s/celltype_hpo_scores.csv", output_dir))
  ct_hpo_df
}

#' 从细胞类型 HPO 分数汇总得到全局加权 HPO term 集合
#'
#' @param ct_hpo_df aggregate_celltype_hpo_scores() 的返回值
#' @return data.frame：hpo_id, total_weight, source_cell_types
build_weighted_hpo_from_scores <- function(ct_hpo_df) {
  message("[Step 4b] 构建全局加权 HPO term 集合...")
  if (nrow(ct_hpo_df) == 0) return(data.frame())

  result <- ct_hpo_df %>%
    dplyr::group_by(hpo_id) %>%
    dplyr::summarise(
      total_weight      = sum(mean_score, na.rm = TRUE),
      source_cell_types = paste(unique(cell_type), collapse = "|"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_weight))

  message(sprintf("  全局 HPO term 数: %d", nrow(result)))
  result
}

# ==============================================================================
# 5. Phen2Gene API 调用（可选）
# ==============================================================================

#' 调用 Phen2Gene REST API
#'
#' @param hpo_ids  HPO ID 向量
#' @param top_k    返回 top-k 基因
#' @return data.frame（Gene, Score, Rank）或空 data.frame
query_phen2gene_D <- function(hpo_ids, top_k = PHEN2GENE_TOP_K) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    message("[SKIP] httr/jsonlite 未安装，跳过 Phen2Gene API")
    return(data.frame())
  }
  message(sprintf("[Step 5] 调用 Phen2Gene API，%d 个 HPO terms...",
                  length(hpo_ids)))

  payload <- jsonlite::toJSON(
    list(HPO_list = paste(hpo_ids, collapse = ";")),
    auto_unbox = TRUE
  )
  resp <- tryCatch(
    httr::POST(PHEN2GENE_API_URL,
               httr::content_type_json(),
               body = payload),
    error = function(e) {
      message(sprintf("  Phen2Gene API 调用失败: %s", e$message))
      NULL
    }
  )
  if (is.null(resp)) return(data.frame())

  result <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                       simplifyDataFrame = TRUE),
    error = function(e) NULL
  )
  if (is.null(result) || !is.data.frame(result)) {
    message("  Phen2Gene 返回格式异常")
    return(data.frame())
  }

  result <- utils::head(result[order(result$Score, decreasing = TRUE), ], top_k)
  message(sprintf("  Phen2Gene 返回 %d 个候选基因", nrow(result)))
  result
}

# ==============================================================================
# 6. HPO 反向投影基因排序（与路线A步骤6相同的实现）
# ==============================================================================

rank_genes_from_weighted_hpo <- function(weighted_hpo, hpo_df, top_k = FINAL_TOP_K) {
  if (nrow(weighted_hpo) == 0) {
    message("[WARN] 加权 HPO 集合为空，无法排序")
    return(data.frame())
  }
  message("[Step 6] HPO 反向投影基因排序...")
  `%||%` <- function(x, y) if (!is.null(x)) x else y

  hpo_gene_map <- split(toupper(hpo_df$gene_symbol), hpo_df$hpo_id)

  gene_scores <- list()
  for (i in seq_len(nrow(weighted_hpo))) {
    term <- weighted_hpo$hpo_id[i]
    wt   <- weighted_hpo$total_weight[i]
    gs   <- hpo_gene_map[[term]]
    if (is.null(gs)) next
    for (g in gs) {
      gene_scores[[g]] <- (gene_scores[[g]] %||% 0) + wt
    }
  }
  if (length(gene_scores) == 0) return(data.frame())

  gene_df <- data.frame(
    gene            = names(gene_scores),
    module_score_sum = unlist(gene_scores),
    stringsAsFactors = FALSE
  )
  gene_df <- gene_df[order(gene_df$module_score_sum, decreasing = TRUE), ]
  utils::head(gene_df, top_k)
}

# ==============================================================================
# 7. 融合 scRNA 证据：对基因排序做二次加权
# ==============================================================================

#' 融合 Phen2Gene 排序与 scRNA 模块分数证据
#'
#' 策略：
#   final_score = phen2gene_rank_score × (1 + module_score_bonus)
#   其中 module_score_bonus = gene_module_score_sum / max(gene_module_score_sum)
#
#' @param phen2gene_df  query_phen2gene_D() 的返回值
#' @param module_df     rank_genes_from_weighted_hpo() 的返回值
#' @param top_k         返回 top-k 基因
#' @return data.frame：gene, phen2gene_score, module_score_bonus, final_score, rank
fuse_rankings <- function(phen2gene_df, module_df, top_k = FINAL_TOP_K) {
  message("[Step 7] 融合 Phen2Gene 排序与模块分数证据...")

  # 若 Phen2Gene 结果为空，直接返回模块分数排序
  if (nrow(phen2gene_df) == 0 && nrow(module_df) > 0) {
    module_df$final_score <- module_df$module_score_sum
    module_df$rank <- seq_len(nrow(module_df))
    return(utils::head(module_df, top_k))
  }
  # 若模块分数为空，直接返回 Phen2Gene 排序
  if (nrow(module_df) == 0 && nrow(phen2gene_df) > 0) {
    phen2gene_df$final_score <- phen2gene_df$Score
    phen2gene_df$rank <- seq_len(nrow(phen2gene_df))
    return(utils::head(phen2gene_df, top_k))
  }
  if (nrow(phen2gene_df) == 0 && nrow(module_df) == 0) {
    return(data.frame())
  }

  # 标准化 Phen2Gene 得分
  pg2 <- phen2gene_df
  pg2$gene            <- toupper(pg2$Gene)
  pg2$phen2gene_score <- pg2$Score / max(pg2$Score, na.rm = TRUE)

  # 标准化模块分数
  mod <- module_df
  mod$gene <- toupper(mod$gene)
  mod$module_bonus <- mod$module_score_sum /
    max(mod$module_score_sum, na.rm = TRUE)

  # 合并
  fused <- merge(pg2[, c("gene", "phen2gene_score")],
                 mod[, c("gene", "module_bonus")],
                 by = "gene", all.x = TRUE)
  fused$module_bonus[is.na(fused$module_bonus)] <- 0
  fused$final_score <- fused$phen2gene_score * (1 + fused$module_bonus)
  fused <- fused[order(fused$final_score, decreasing = TRUE), ]
  fused$rank <- seq_len(nrow(fused))

  utils::head(fused, top_k)
}

# ==============================================================================
# 8. 主流程
# ==============================================================================

#' 运行路线 D 完整流程
#'
#' @param seurat_obj  Seurat 对象（NULL=使用 Toy Demo）
#' @param hpo_path    HPO 文件路径
#' @param target_hpo_ids 目标 HPO IDs（NULL=自动推断）
#' @param output_dir  输出目录
#' @return list(ct_hpo_scores, weighted_hpo, phen2gene_result, gene_ranking)
run_pipeline_D <- function(seurat_obj     = NULL,
                            hpo_path       = HPO_GENE_PATH,
                            target_hpo_ids = TARGET_HPO_IDS,
                            output_dir     = OUTPUT_DIR) {
  .load_packages_D()
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Toy Demo ---
  if (is.null(seurat_obj)) {
    message("=== [Toy Demo] 未提供 Seurat 对象，生成模拟数据 ===")
    seurat_obj <- .make_toy_seurat_D()
  }

  # Step 2: 加载 HPO 基因集
  hpo_data      <- load_hpo_for_module_score(hpo_path)
  hpo_df        <- hpo_data$hpo_df
  hpo_gene_sets <- hpo_data$hpo_gene_sets

  # 若指定了目标 HPO，仅计算这些 term 的模块分数
  if (!is.null(target_hpo_ids)) {
    hpo_gene_sets <- hpo_gene_sets[
      names(hpo_gene_sets) %in% target_hpo_ids
    ]
    if (length(hpo_gene_sets) == 0) {
      message("[WARN] 指定的 target_hpo_ids 在 HPO 数据中找不到，使用全部 terms")
      hpo_gene_sets <- hpo_data$hpo_gene_sets
    }
  }

  # Step 3: 模块评分
  seurat_scored <- tryCatch(
    compute_module_scores(seurat_obj, hpo_gene_sets,
                          method = SCORE_METHOD,
                          output_dir = file.path(output_dir, "module_scores")),
    error = function(e) {
      message(sprintf("[WARN] 模块评分失败: %s", e$message))
      seurat_obj
    }
  )

  # Step 4: 汇总细胞类型×HPO 分数
  ct_hpo_scores <- aggregate_celltype_hpo_scores(
    seurat_scored,
    cell_type_col = CELL_TYPE_COL,
    hpo_ids       = names(hpo_gene_sets),
    top_n         = TOP_HPO_PER_CELLTYPE,
    output_dir    = file.path(output_dir, "celltype_hpo_scores")
  )

  # Step 4b: 全局加权 HPO 集合
  weighted_hpo <- build_weighted_hpo_from_scores(ct_hpo_scores)
  if (nrow(weighted_hpo) > 0) {
    write.csv(weighted_hpo,
              file.path(output_dir, "weighted_hpo_terms_D.csv"),
              row.names = FALSE)
  }

  # Step 5: 可选 Phen2Gene API
  phen2gene_result <- data.frame()
  if (USE_PHEN2GENE_API && nrow(weighted_hpo) > 0) {
    top_hpo_ids <- utils::head(weighted_hpo$hpo_id, 10)
    phen2gene_result <- query_phen2gene_D(top_hpo_ids)
    if (nrow(phen2gene_result) > 0) {
      write.csv(phen2gene_result,
                file.path(output_dir, "phen2gene_result_D.csv"),
                row.names = FALSE)
    }
  }

  # Step 6: 模块分数反向投影基因排序
  module_gene_ranking <- rank_genes_from_weighted_hpo(weighted_hpo, hpo_df)

  # Step 7: 融合排序
  gene_ranking <- fuse_rankings(phen2gene_result, module_gene_ranking)

  # 若融合结果为空，回退到模块分数排序
  if (nrow(gene_ranking) == 0) {
    gene_ranking <- module_gene_ranking
  }

  if (nrow(gene_ranking) > 0) {
    out_file <- file.path(output_dir, "gene_ranking_D.csv")
    write.csv(gene_ranking, out_file, row.names = FALSE)
    message(sprintf("[Done] 基因排序结果保存至: %s", out_file))
    message("\n=== Top-20 候选基因（模块分数×HPO 融合排序） ===")
    print(utils::head(gene_ranking, 20))
  }

  invisible(list(
    ct_hpo_scores    = ct_hpo_scores,
    weighted_hpo     = weighted_hpo,
    phen2gene_result = phen2gene_result,
    gene_ranking     = gene_ranking
  ))
}

# ==============================================================================
# Toy 数据生成
# ==============================================================================

.make_toy_seurat_D <- function(n_cells = 600, n_genes = 500, seed = 42) {
  set.seed(seed)
  counts <- matrix(
    stats::rnbinom(n_cells * n_genes, mu = 2, size = 1),
    nrow = n_genes, ncol = n_cells
  )
  gene_names <- paste0("GENE", seq_len(n_genes))
  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Cell", seq_len(n_cells))

  meta <- data.frame(
    sample_id  = rep(paste0("S", 1:6), each = n_cells / 6),
    disease    = rep(c("CRSwNP", "Control"), each = n_cells / 2),
    Annotation_2 = rep(c("CD4_T", "CD8_T", "Macrophage"), each = n_cells / 3),
    row.names  = colnames(counts),
    stringsAsFactors = FALSE
  )

  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  message(sprintf("  Toy Seurat: %d 细胞 × %d 基因", ncol(obj), nrow(obj)))
  obj
}

.make_toy_hpo_df <- function() {
  genes <- paste0("GENE", seq_len(500))
  hpo_ids <- paste0("HP:000", sprintf("%04d", seq_len(50)))
  set.seed(1)
  rows <- do.call(rbind, lapply(hpo_ids, function(h) {
    data.frame(hpo_id = h,
               gene_symbol = sample(genes, 10, replace = FALSE),
               stringsAsFactors = FALSE)
  }))
  rows
}

# ==============================================================================
# 脚本入口
# ==============================================================================

if (!interactive()) {
  seurat_input <- if (!is.null(SEURAT_RDS_PATH) && file.exists(SEURAT_RDS_PATH)) {
    message(sprintf("加载 Seurat RDS: %s", SEURAT_RDS_PATH))
    readRDS(SEURAT_RDS_PATH)
  } else {
    NULL
  }

  result <- run_pipeline_D(
    seurat_obj = seurat_input,
    output_dir = OUTPUT_DIR
  )

  message("\n=== 路线 D 流程完成 ===")
  message(sprintf("结果保存至: %s", OUTPUT_DIR))
}
