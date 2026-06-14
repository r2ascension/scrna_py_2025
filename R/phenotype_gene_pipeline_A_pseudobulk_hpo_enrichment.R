#!/usr/bin/env Rscript
# ==============================================================================
# 表型→基因关联 路线A：Pseudobulk DE → HPO富集 → 语义相似度 → 基因排序
# Phenotype→Gene Pipeline A: Pseudobulk DE → HPO Enrichment → Gene Ranking
# ==============================================================================
#
# 工作流摘要 (Workflow Summary):
#   1. 从 Seurat 对象按 (sample × cell_type) 聚合计数 → pseudobulk 矩阵
#   2. 用 DESeq2 / edgeR 对每种细胞类型做差异表达分析
#   3. 将每种细胞类型的上调/差异基因投影到 HPO 基因集 → ORA 富集
#   4. 整合各细胞类型的显著 HPO term → 加权 HPO term 集合
#   5. 选项 A4a: 用 HPO term 语义相似度矩阵 + 基因–表型注释反推基因排序
#   6. 选项 A4b: 调用 Phen2Gene REST API（需网络），输出候选基因排序
#
# 输入 (Inputs):
#   seurat_obj     : Seurat 对象，包含 raw counts、细胞类型注释、样本 ID、分组
#   condition_col  : 比较条件列名（如 "disease"，取值 "CRSwNP" vs "Control"）
#   cell_type_col  : 细胞类型列名（如 "Annotation_2"）
#   sample_col     : 样本 ID 列名（如 "sample"）
#
# 输出 (Outputs):
#   <output_dir>/pseudobulk_de/    : 各细胞类型 DESeq2/edgeR 结果 CSV
#   <output_dir>/hpo_enrichment/   : 各细胞类型 HPO ORA 结果 CSV + 图
#   <output_dir>/gene_ranking/     : 汇总基因排序表 CSV
#
# 依赖包 (Required Packages):
#   Seurat, DESeq2, edgeR (任选其一), clusterProfiler, org.Hs.eg.db
#   dplyr, ggplot2, httr (Phen2Gene REST 可选)
#
# 参考 (References):
#   - Crowell et al. (2020) muscat: multi-sample, multi-group scRNA DS
#     Nature Communications. https://doi.org/10.1038/s41467-020-19894-4
#   - Seurat AggregateExpression pseudobulk:
#     https://satijalab.org/seurat/articles/de_vignette
#   - HPO gene sets / phenotype_to_genes.txt:
#     https://hpo.jax.org/app/download/annotation
#   - Phen2Gene REST API: https://github.com/WGLab/Phen2Gene
#
# 版本 (Version): v1.0  2026-03-27
# ==============================================================================

# ==============================================================================
# 0. 配置参数 Configuration
# ==============================================================================

# --- 真实路径替换区 (Replace paths for real analysis) ---
SEURAT_RDS_PATH   <- NULL   # 替换为你的 RDS 路径，如 "/data/T_object.rds"
HPO_ANNO_PATH     <- NULL   # HPO phenotype_to_genes.txt 本地路径（NULL=下载）
MSIGDB_HPO_GMT    <- NULL   # MSigDB HPO GMT 本地路径（NULL=用 clusterProfiler msigdbr）
OUTPUT_DIR        <- "./results_pipeline_A"

# --- 分析参数 ---
CONDITION_COL     <- "disease"        # 比较条件列
CONDITION_TEST    <- "CRSwNP"         # 实验组
CONDITION_CTRL    <- "Control"        # 对照组
CELL_TYPE_COL     <- "Annotation_2"   # 细胞类型列
SAMPLE_COL        <- "sample"         # 样本 ID 列
MIN_CELL_PER_SAMPLE <- 3              # 每个 pseudobulk 样本最小细胞数
MIN_SAMPLE_PER_GROUP <- 3            # 每组最小样本数
DE_PVAL_CUTOFF    <- 0.05             # DE p.adj 阈值
DE_LOGFC_CUTOFF   <- 0.5             # |log2FC| 阈值
ORA_PVAL_CUTOFF   <- 0.05            # HPO ORA p.adj 阈值
TOP_HPO_PER_CELLTYPE <- 20           # 每种细胞类型保留 top HPO terms 数
USE_PHEN2GENE_API <- FALSE           # 是否调用 Phen2Gene REST API
PHEN2GENE_API_URL <- "https://phen2gene.wglab.org/api/query"  # Phen2Gene API

# ==============================================================================
# 1. 加载依赖包 Load Packages
# ==============================================================================

.load_packages <- function() {
  pkgs <- c("Seurat", "dplyr", "ggplot2", "DESeq2")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("请安装 %s 包: BiocManager::install('%s')", pkg, pkg))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # 富集包（可选）
  for (pkg in c("clusterProfiler", "org.Hs.eg.db", "msigdbr")) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    } else {
      message(sprintf("[WARN] %s 未安装，部分功能将跳过", pkg))
    }
  }
  if (USE_PHEN2GENE_API && !requireNamespace("httr", quietly = TRUE)) {
    message("[WARN] httr 未安装，将跳过 Phen2Gene API 调用")
  }
}

# ==============================================================================
# 2. Pseudobulk 聚合 Pseudobulk Aggregation
# ==============================================================================

#' 按 (sample × cell_type) 聚合 raw counts
#'
#' @param seurat_obj Seurat 对象，包含 raw counts
#' @param cell_type_col 细胞类型列名
#' @param sample_col 样本 ID 列名
#' @param condition_col 条件列名
#' @param min_cell 每个 pseudobulk 样本最小细胞数
#' @return list(counts_list, meta_list)，按细胞类型分组
aggregate_pseudobulk <- function(seurat_obj,
                                  cell_type_col = "Annotation_2",
                                  sample_col = "sample",
                                  condition_col = "disease",
                                  min_cell = 3) {
  message("[Step 2] 聚合 pseudobulk 计数矩阵...")

  # 兼容 Seurat v4/v5
  meta <- seurat_obj@meta.data

  # 检查必要列
  req_cols <- c(cell_type_col, sample_col, condition_col)
  missing <- setdiff(req_cols, colnames(meta))
  if (length(missing) > 0) {
    stop(sprintf("元数据缺少列: %s", paste(missing, collapse = ", ")))
  }

  cell_types <- unique(meta[[cell_type_col]])
  message(sprintf("  发现 %d 种细胞类型", length(cell_types)))

  counts_list <- list()
  meta_list   <- list()

  for (ct in cell_types) {
    ct_cells  <- rownames(meta)[meta[[cell_type_col]] == ct]
    ct_seurat <- seurat_obj[, ct_cells]

    # 用 AggregateExpression（Seurat v5 原生 pseudobulk 函数）
    # 若为 Seurat v4，回退到手动 tapply 聚合
    pb <- tryCatch(
      AggregateExpression(
        ct_seurat,
        group.by      = c(sample_col, condition_col),
        assays        = "RNA",
        slot          = "counts",
        return.seurat = FALSE,
        verbose       = FALSE
      )$RNA,
      error = function(e) {
        # Seurat v4 回退：手动按 sample 汇总
        message(sprintf("    AggregateExpression 失败，使用手动聚合: %s", e$message))
        .manual_aggregate(ct_seurat, sample_col)
      }
    )

    # 提取 pseudobulk 元数据（样本 × 条件）
    pb_meta <- .build_pb_meta(ct_seurat, colnames(pb), sample_col, condition_col)

    # 过滤细胞数不足的 pseudobulk 样本
    cell_counts <- table(ct_seurat@meta.data[[sample_col]])
    keep_samples <- names(cell_counts)[cell_counts >= min_cell]
    keep_idx <- pb_meta[[sample_col]] %in% keep_samples
    pb      <- pb[, keep_idx, drop = FALSE]
    pb_meta <- pb_meta[keep_idx, , drop = FALSE]

    if (ncol(pb) < 2) {
      message(sprintf("    [SKIP] %s: 有效 pseudobulk 样本 < 2，跳过", ct))
      next
    }

    counts_list[[ct]] <- as.matrix(pb)
    meta_list[[ct]]   <- pb_meta
    message(sprintf("    %s: %d pseudobulk 样本，%d 基因",
                    ct, ncol(pb), nrow(pb)))
  }
  list(counts = counts_list, meta = meta_list)
}

.manual_aggregate <- function(seurat_sub, sample_col) {
  # 手动按 sample_id 汇总（Seurat v4 兼容）
  meta <- seurat_sub@meta.data
  counts_mat <- tryCatch(
    LayerData(seurat_sub, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(seurat_sub, assay = "RNA", slot = "counts")
  )
  samples <- unique(meta[[sample_col]])
  pb_list <- lapply(samples, function(s) {
    s_cells <- rownames(meta)[meta[[sample_col]] == s]
    Matrix::rowSums(counts_mat[, s_cells, drop = FALSE])
  })
  do.call(cbind, stats::setNames(pb_list, samples))
}

.build_pb_meta <- function(seurat_sub, pb_colnames, sample_col, condition_col) {
  meta <- seurat_sub@meta.data
  sample_condition <- unique(meta[, c(sample_col, condition_col)])
  rownames(sample_condition) <- NULL
  sample_condition
}

# ==============================================================================
# 3. 差异表达分析 Differential Expression (DESeq2)
# ==============================================================================

#' 对每种细胞类型运行 DESeq2 pseudobulk DE
#'
#' @param pb_data aggregate_pseudobulk() 的返回值
#' @param condition_col 条件列名
#' @param test 实验组值
#' @param ctrl 对照组值
#' @param pval_cutoff p.adj 阈值
#' @param logfc_cutoff |log2FC| 阈值
#' @param output_dir DE 结果输出目录
#' @return list（按细胞类型），每个元素为 data.frame DE 结果
run_deseq2_pseudobulk <- function(pb_data,
                                   condition_col = "disease",
                                   test = "CRSwNP",
                                   ctrl = "Control",
                                   pval_cutoff  = 0.05,
                                   logfc_cutoff = 0.5,
                                   min_samples  = 3,
                                   output_dir   = "./pseudobulk_de") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("[Step 3] 运行 DESeq2 pseudobulk 差异表达...")

  de_results <- list()

  for (ct in names(pb_data$counts)) {
    pb_mat  <- pb_data$counts[[ct]]
    pb_meta <- pb_data$meta[[ct]]

    # 过滤：每组至少 min_samples 个样本
    grp_counts <- table(pb_meta[[condition_col]])
    if (!all(c(test, ctrl) %in% names(grp_counts)) ||
        any(grp_counts[c(test, ctrl)] < min_samples)) {
      message(sprintf("  [SKIP] %s: 组内样本数不足 %d", ct, min_samples))
      next
    }

    # 对齐列名
    common_samples <- intersect(colnames(pb_mat), pb_meta[[SAMPLE_COL]])
    if (length(common_samples) == 0) {
      # 当 pb_mat 列名格式为 "sample_condition"，尝试用 meta 行顺序对应
      pb_meta_ord <- pb_meta
    } else {
      pb_meta_ord <- pb_meta[match(colnames(pb_mat), pb_meta[[SAMPLE_COL]]), ]
    }
    pb_meta_ord <- pb_meta_ord[!is.na(pb_meta_ord[[condition_col]]), ]
    pb_mat      <- pb_mat[, rownames(pb_meta_ord) %in% seq_len(ncol(pb_mat))]

    # 构建 DESeq2 对象
    col_data <- data.frame(
      condition = factor(pb_meta[[condition_col]], levels = c(ctrl, test)),
      row.names = colnames(pb_mat)
    )
    col_data <- col_data[!is.na(col_data$condition), , drop = FALSE]
    pb_mat   <- pb_mat[, rownames(col_data), drop = FALSE]

    # 过滤低表达基因（至少在一半样本中 count ≥ 1）
    keep_genes <- rowSums(pb_mat >= 1) >= (ncol(pb_mat) / 2)
    pb_mat     <- pb_mat[keep_genes, , drop = FALSE]

    message(sprintf("  %s: %d 样本 × %d 基因", ct, ncol(pb_mat), nrow(pb_mat)))

    dds <- tryCatch({
      DESeq2::DESeqDataSetFromMatrix(
        countData = round(pb_mat),
        colData   = col_data,
        design    = ~condition
      )
    }, error = function(e) {
      message(sprintf("    DESeq2 初始化失败: %s", e$message))
      return(NULL)
    })
    if (is.null(dds)) next

    dds <- tryCatch(
      DESeq2::DESeq(dds, quiet = TRUE),
      error = function(e) {
        message(sprintf("    DESeq2 运行失败: %s", e$message))
        NULL
      }
    )
    if (is.null(dds)) next

    res <- DESeq2::results(dds, contrast = c("condition", test, ctrl),
                           independentFiltering = TRUE)
    res_df <- as.data.frame(res)
    res_df$gene        <- rownames(res_df)
    res_df$cell_type   <- ct

    # 标注显著性
    res_df$significant <- !is.na(res_df$padj) &
      res_df$padj < pval_cutoff &
      abs(res_df$log2FoldChange) > logfc_cutoff

    # 保存全部结果
    out_file <- file.path(output_dir,
                          sprintf("de_%s.csv", gsub("[^A-Za-z0-9_]", "_", ct)))
    write.csv(res_df, out_file, row.names = FALSE)
    message(sprintf("    显著 DE 基因: %d  → 保存至 %s",
                    sum(res_df$significant, na.rm = TRUE), out_file))

    de_results[[ct]] <- res_df
  }

  de_results
}

# ==============================================================================
# 4. HPO 基因集富集分析 HPO ORA Enrichment
# ==============================================================================

#' 下载并解析 HPO phenotype_to_genes.txt
#'
#' @param local_path 本地文件路径（NULL=在线下载）
#' @return data.frame：hpo_id, hpo_name, gene_symbol
load_hpo_gene_sets <- function(local_path = NULL) {
  if (!is.null(local_path) && file.exists(local_path)) {
    message("[Step 4] 从本地加载 HPO 基因集: ", local_path)
    hpo_raw <- read.table(local_path, sep = "\t", header = TRUE,
                          quote = "", comment.char = "#")
  } else {
    message("[Step 4] 尝试下载 HPO phenotype_to_genes.txt...")
    tmp <- tempfile(fileext = ".txt")
    url <- paste0("https://purl.obolibrary.org/obo/hp/hpoa/",
                  "phenotype_to_genes.txt")
    tryCatch(
      download.file(url, destfile = tmp, quiet = TRUE, method = "libcurl"),
      error = function(e) stop("HPO 下载失败，请手动下载并设置 HPO_ANNO_PATH")
    )
    hpo_raw <- read.table(tmp, sep = "\t", header = FALSE, quote = "",
                          skip = 1, fill = TRUE)
    colnames(hpo_raw) <- c("hpo_id", "hpo_name", "entrez_id", "gene_symbol",
                            "additional", "source", "disease_id")[seq_len(ncol(hpo_raw))]
  }
  # 标准化列名
  if (!"gene_symbol" %in% colnames(hpo_raw)) {
    # phenotype_to_genes.txt 新版格式尝试列 4
    colnames(hpo_raw)[4] <- "gene_symbol"
  }
  if (!"hpo_id" %in% colnames(hpo_raw)) {
    colnames(hpo_raw)[1] <- "hpo_id"
  }
  hpo_raw[, c("hpo_id", "gene_symbol")]
}

#' 用 clusterProfiler ORA 对每种细胞类型的 DE 基因做 HPO 富集
#'
#' @param de_results run_deseq2_pseudobulk() 的返回值
#' @param hpo_df    load_hpo_gene_sets() 的返回值
#' @param universe  背景基因集（默认=所有被测基因的并集）
#' @param pval_cutoff ORA p.adj 阈值
#' @param top_n 保留每种细胞类型 top_n HPO terms
#' @param output_dir 输出目录
#' @return list（按细胞类型），每个元素为 clusterProfiler ORA 结果
run_hpo_ora <- function(de_results,
                         hpo_df,
                         universe       = NULL,
                         pval_cutoff    = 0.05,
                         top_n          = 20,
                         output_dir     = "./hpo_enrichment") {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("请安装 clusterProfiler: BiocManager::install('clusterProfiler')")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("[Step 4] 运行 HPO ORA 富集分析...")

  # 构建 TERM2GENE data.frame（clusterProfiler 格式）
  term2gene <- data.frame(
    term   = hpo_df$hpo_id,
    gene   = toupper(hpo_df$gene_symbol),
    stringsAsFactors = FALSE
  )

  # 背景基因（全被测基因）
  if (is.null(universe)) {
    universe <- toupper(unique(unlist(
      lapply(de_results, function(x) x$gene)
    )))
  }

  ora_results <- list()

  for (ct in names(de_results)) {
    de_df <- de_results[[ct]]

    # 取上调显著 DE 基因
    sig_up <- toupper(de_df$gene[
      !is.na(de_df$padj) & de_df$padj < DE_PVAL_CUTOFF &
        de_df$log2FoldChange > DE_LOGFC_CUTOFF
    ])

    if (length(sig_up) < 3) {
      message(sprintf("  [SKIP] %s: 上调 DE 基因 < 3，跳过 HPO ORA", ct))
      next
    }

    ora <- tryCatch(
      clusterProfiler::enricher(
        gene         = sig_up,
        TERM2GENE    = term2gene,
        universe     = universe,
        pvalueCutoff = pval_cutoff,
        pAdjustMethod = "BH",
        minGSSize    = 5,
        maxGSSize    = 500
      ),
      error = function(e) {
        message(sprintf("  [WARN] %s ORA 失败: %s", ct, e$message))
        NULL
      }
    )

    if (is.null(ora) || nrow(as.data.frame(ora)) == 0) {
      message(sprintf("  %s: 无显著 HPO terms", ct))
      next
    }

    ora_df <- as.data.frame(ora)
    ora_df <- ora_df[order(ora_df$p.adjust), ]
    ora_df <- utils::head(ora_df, top_n)

    # 保存
    out_file <- file.path(output_dir,
                          sprintf("hpo_ora_%s.csv", gsub("[^A-Za-z0-9_]", "_", ct)))
    write.csv(ora_df, out_file, row.names = FALSE)
    message(sprintf("  %s: %d 显著 HPO terms → %s", ct, nrow(ora_df), out_file))

    # 绘制 barplot
    tryCatch({
      p <- ggplot2::ggplot(
        utils::head(ora_df, 15),
        ggplot2::aes(x = -log10(p.adjust),
                     y = stats::reorder(ID, -log10(p.adjust)))
      ) +
        ggplot2::geom_bar(stat = "identity", fill = "#4E9FD4") +
        ggplot2::labs(
          title = sprintf("HPO ORA: %s", ct),
          x = "-log10(p.adj)", y = "HPO Term"
        ) +
        ggplot2::theme_minimal(base_size = 10)
      ggplot2::ggsave(
        file.path(output_dir, sprintf("hpo_barplot_%s.pdf",
                                      gsub("[^A-Za-z0-9_]", "_", ct))),
        p, width = 8, height = 5
      )
    }, error = function(e) NULL)

    ora_results[[ct]] <- ora_df
  }

  ora_results
}

# ==============================================================================
# 5. 构建加权 HPO 词条集合 Build Weighted HPO Term Set
# ==============================================================================

#' 整合各细胞类型 ORA 结果，生成全局加权 HPO term 集合
#'
#' @param ora_results run_hpo_ora() 的返回值
#' @param weight_by  用哪一列做权重（"Count"=基因数, "-log10pAdj"=显著性）
#' @return data.frame：hpo_id, weight, source_cell_types
build_weighted_hpo_set <- function(ora_results,
                                    weight_by = "-log10pAdj") {
  message("[Step 5] 构建加权 HPO term 集合...")

  all_rows <- lapply(names(ora_results), function(ct) {
    df <- ora_results[[ct]]
    df$cell_type <- ct
    df$weight <- if (weight_by == "Count") {
      df$Count
    } else {
      -log10(pmax(df$p.adjust, 1e-300))
    }
    df[, c("ID", "weight", "cell_type")]
  })
  all_df <- do.call(rbind, all_rows)

  # 按 HPO term 聚合（多个细胞类型对同一 term 的权重求和）
  hpo_agg <- all_df %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      total_weight      = sum(weight),
      source_cell_types = paste(unique(cell_type), collapse = "|"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_weight))

  message(sprintf("  汇总得到 %d 个 HPO terms", nrow(hpo_agg)))
  hpo_agg
}

# ==============================================================================
# 6. 基于 HPO→基因 映射的基因排序 Gene Ranking via HPO→Gene Back-Projection
# ==============================================================================

#' 将加权 HPO term 集合反向投影到基因，得到候选基因排序
#'
#' @param weighted_hpo build_weighted_hpo_set() 的返回值
#' @param hpo_df       load_hpo_gene_sets() 的返回值（HPO→gene 映射）
#' @param de_results   run_deseq2_pseudobulk() 的返回值（用于过滤与加权）
#' @param top_k        返回 top-k 基因
#' @return data.frame：gene, phenotype_score, hpo_terms, supported_by_de
rank_genes_from_hpo <- function(weighted_hpo, hpo_df, de_results, top_k = 200) {
  message("[Step 6] 通过 HPO→基因 映射排序候选基因...")

  # 构建 HPO→gene 查找表
  hpo_gene_map <- split(toupper(hpo_df$gene_symbol), hpo_df$hpo_id)

  # 对每个 HPO term，用其权重给关联基因累积得分
  gene_scores <- list()
  for (i in seq_len(nrow(weighted_hpo))) {
    term    <- weighted_hpo$ID[i]
    wt      <- weighted_hpo$total_weight[i]
    genes_i <- hpo_gene_map[[term]]
    if (is.null(genes_i)) next
    for (g in genes_i) {
      gene_scores[[g]] <- (gene_scores[[g]] %||% 0) + wt
    }
  }

  if (length(gene_scores) == 0) {
    message("  [WARN] 未能从 HPO 映射得到任何基因")
    return(data.frame())
  }

  gene_df <- data.frame(
    gene             = names(gene_scores),
    phenotype_score  = unlist(gene_scores),
    stringsAsFactors = FALSE
  )
  gene_df <- gene_df[order(gene_df$phenotype_score, decreasing = TRUE), ]

  # 整合 DE 证据：标注在哪些细胞类型中为显著上调
  de_upregulated <- lapply(de_results, function(df) {
    df$gene[!is.na(df$padj) & df$padj < DE_PVAL_CUTOFF &
              df$log2FoldChange > DE_LOGFC_CUTOFF]
  })
  gene_df$supported_by_de <- sapply(toupper(gene_df$gene), function(g) {
    cts_supporting <- names(de_upregulated)[sapply(de_upregulated, function(v) g %in% toupper(v))]
    if (length(cts_supporting) == 0) "" else paste(cts_supporting, collapse = "|")
  })

  # 添加 HPO terms 注释
  gene_df$hpo_terms <- sapply(toupper(gene_df$gene), function(g) {
    terms_for_gene <- weighted_hpo$ID[
      sapply(weighted_hpo$ID, function(t) g %in% (hpo_gene_map[[t]] %||% character(0)))
    ]
    if (length(terms_for_gene) == 0) "" else paste(utils::head(terms_for_gene, 5), collapse = "|")
  })

  gene_df <- utils::head(gene_df, top_k)
  message(sprintf("  Top-%d 候选基因排序完成", nrow(gene_df)))
  gene_df
}

# NULL 合并运算符（R < 4.4 兼容）
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ==============================================================================
# 7. Phen2Gene REST API 调用（可选）Phen2Gene API (Optional)
# ==============================================================================

#' 调用 Phen2Gene REST API，输入 HPO ID 列表，返回候选基因排序
#'
#' @param hpo_ids HPO ID 向量（如 c("HP:0000118", "HP:0001250")）
#' @param top_k   返回 top-k 基因
#' @return data.frame：Gene, Score, Rank（或空 data.frame 如果失败）
query_phen2gene <- function(hpo_ids, top_k = 100) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    message("[SKIP] httr 未安装，跳过 Phen2Gene API 调用")
    return(data.frame())
  }
  message(sprintf("[Step 7] 调用 Phen2Gene API，HPO terms: %d 个...",
                  length(hpo_ids)))

  payload <- list(HPO_list = paste(hpo_ids, collapse = ";"))
  resp <- tryCatch(
    httr::POST(PHEN2GENE_API_URL,
               httr::content_type_json(),
               body = jsonlite::toJSON(payload, auto_unbox = TRUE)),
    error = function(e) {
      message(sprintf("  Phen2Gene API 调用失败: %s", e$message))
      NULL
    }
  )
  if (is.null(resp)) return(data.frame())

  content_text <- httr::content(resp, as = "text", encoding = "UTF-8")
  result <- tryCatch(
    jsonlite::fromJSON(content_text, simplifyDataFrame = TRUE),
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
# 8. 主流程 Main Pipeline
# ==============================================================================

#' 运行路线 A 完整流程
#'
#' @param seurat_obj   Seurat 对象（NULL=使用内置 Toy 数据演示）
#' @param condition_col 条件列名
#' @param test          实验组值
#' @param ctrl          对照组值
#' @param cell_type_col 细胞类型列名
#' @param sample_col    样本列名
#' @param hpo_path      HPO 文件本地路径
#' @param output_dir    输出目录
#' @return list(de_results, ora_results, weighted_hpo, gene_ranking, phen2gene_result)
run_pipeline_A <- function(seurat_obj      = NULL,
                            condition_col   = CONDITION_COL,
                            test            = CONDITION_TEST,
                            ctrl            = CONDITION_CTRL,
                            cell_type_col   = CELL_TYPE_COL,
                            sample_col      = SAMPLE_COL,
                            hpo_path        = HPO_ANNO_PATH,
                            output_dir      = OUTPUT_DIR) {
  .load_packages()
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # --- 真实数据 or Toy Demo ---
  if (is.null(seurat_obj)) {
    message("=== [Toy Demo] 未提供 Seurat 对象，生成模拟数据 ===")
    seurat_obj <- .make_toy_seurat()
    condition_col <- "condition"
    test          <- "case"
    ctrl          <- "control"
    cell_type_col <- "cell_type"
    sample_col    <- "sample_id"
  }

  # Step 2: Pseudobulk 聚合
  pb_data <- aggregate_pseudobulk(
    seurat_obj,
    cell_type_col = cell_type_col,
    sample_col    = sample_col,
    condition_col = condition_col,
    min_cell      = MIN_CELL_PER_SAMPLE
  )

  # Step 3: DESeq2 DE
  de_dir <- file.path(output_dir, "pseudobulk_de")
  de_results <- run_deseq2_pseudobulk(
    pb_data,
    condition_col = condition_col,
    test          = test,
    ctrl          = ctrl,
    pval_cutoff   = DE_PVAL_CUTOFF,
    logfc_cutoff  = DE_LOGFC_CUTOFF,
    min_samples   = MIN_SAMPLE_PER_GROUP,
    output_dir    = de_dir
  )

  if (length(de_results) == 0) {
    message("[WARN] 无任何细胞类型通过 DE 分析，流程终止")
    return(invisible(NULL))
  }

  # Step 4: HPO 基因集加载
  hpo_df <- tryCatch(
    load_hpo_gene_sets(hpo_path),
    error = function(e) {
      message(sprintf("[WARN] HPO 加载失败: %s\n  使用内置 Toy HPO", e$message))
      .make_toy_hpo()
    }
  )

  # Step 4: HPO ORA
  ora_dir <- file.path(output_dir, "hpo_enrichment")
  ora_results <- run_hpo_ora(
    de_results,
    hpo_df,
    pval_cutoff = ORA_PVAL_CUTOFF,
    top_n       = TOP_HPO_PER_CELLTYPE,
    output_dir  = ora_dir
  )

  if (length(ora_results) == 0) {
    message("[WARN] 无显著 HPO ORA 结果，跳过后续步骤")
    return(invisible(list(de_results = de_results)))
  }

  # Step 5: 加权 HPO term 集合
  weighted_hpo <- build_weighted_hpo_set(ora_results)
  write.csv(weighted_hpo, file.path(output_dir, "weighted_hpo_terms.csv"),
            row.names = FALSE)

  # Step 6: 基因排序
  gene_ranking <- rank_genes_from_hpo(weighted_hpo, hpo_df, de_results)
  gene_rank_file <- file.path(output_dir, "gene_ranking", "gene_ranking.csv")
  dir.create(dirname(gene_rank_file), showWarnings = FALSE, recursive = TRUE)
  write.csv(gene_ranking, gene_rank_file, row.names = FALSE)
  message(sprintf("[Done] 基因排序结果保存至 %s", gene_rank_file))

  # Step 7: Phen2Gene（可选）
  phen2gene_result <- data.frame()
  if (USE_PHEN2GENE_API && nrow(weighted_hpo) > 0) {
    top_hpo_ids <- utils::head(weighted_hpo$ID, 10)
    phen2gene_result <- query_phen2gene(top_hpo_ids)
    if (nrow(phen2gene_result) > 0) {
      write.csv(phen2gene_result,
                file.path(output_dir, "gene_ranking", "phen2gene_result.csv"),
                row.names = FALSE)
    }
  }

  # 打印 top-20 结果
  if (nrow(gene_ranking) > 0) {
    message("\n=== Top-20 候选基因 (HPO 反向投影得分) ===")
    print(utils::head(gene_ranking, 20))
  }

  invisible(list(
    de_results       = de_results,
    ora_results      = ora_results,
    weighted_hpo     = weighted_hpo,
    gene_ranking     = gene_ranking,
    phen2gene_result = phen2gene_result
  ))
}

# ==============================================================================
# 辅助：Toy 数据生成 (Toy Data Helpers)
# ==============================================================================

.make_toy_seurat <- function(n_cells = 600, n_genes = 500, seed = 42) {
  set.seed(seed)
  counts <- matrix(
    stats::rnbinom(n_cells * n_genes, mu = 2, size = 1),
    nrow = n_genes, ncol = n_cells
  )
  gene_names <- paste0("GENE", seq_len(n_genes))
  # 注入少量差异基因（前 30 个基因在 case 中高表达）
  case_cells <- seq(1, n_cells, by = 2)
  counts[1:30, case_cells] <- counts[1:30, case_cells] + 8
  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Cell", seq_len(n_cells))

  meta <- data.frame(
    sample_id  = rep(paste0("S", 1:6), each = n_cells / 6),
    condition  = rep(c("case", "control"), each = n_cells / 2),
    cell_type  = rep(c("CD4_T", "CD8_T", "Macrophage"), each = n_cells / 3),
    row.names  = colnames(counts),
    stringsAsFactors = FALSE
  )

  obj <- Seurat::CreateSeuratObject(
    counts = counts, meta.data = meta
  )
  message(sprintf("  Toy Seurat: %d 细胞 × %d 基因", ncol(obj), nrow(obj)))
  obj
}

.make_toy_hpo <- function() {
  genes <- paste0("GENE", seq_len(500))
  hpo_ids <- paste0("HP:000", sprintf("%04d", seq_len(50)))
  # 每个 HPO term 随机关联 10 个基因
  set.seed(1)
  rows <- do.call(rbind, lapply(hpo_ids, function(h) {
    data.frame(hpo_id = h,
               gene_symbol = sample(genes, 10, replace = FALSE),
               stringsAsFactors = FALSE)
  }))
  rows
}

# ==============================================================================
# 脚本入口 Script Entry Point
# ==============================================================================

if (!interactive()) {
  # 命令行运行：Rscript R/phenotype_gene_pipeline_A_pseudobulk_hpo_enrichment.R
  # 可通过修改顶部配置参数（SEURAT_RDS_PATH 等）切换为真实数据
  seurat_input <- if (!is.null(SEURAT_RDS_PATH) && file.exists(SEURAT_RDS_PATH)) {
    message(sprintf("加载 Seurat RDS: %s", SEURAT_RDS_PATH))
    readRDS(SEURAT_RDS_PATH)
  } else {
    NULL  # 使用 Toy Demo
  }

  result <- run_pipeline_A(
    seurat_obj    = seurat_input,
    output_dir    = OUTPUT_DIR
  )

  message("\n=== 路线 A 流程完成 ===")
  message(sprintf("结果保存至: %s", normalizePath(OUTPUT_DIR, mustWork = FALSE)))
}
