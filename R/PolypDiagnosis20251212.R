suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

peek_seurat_rds <- function(
  rds_path,
  mt_pattern = "^MT-",
  rb_pattern = "^(RPL|RPS)"
) {
  cat("============================================================\n")
  cat("RDS:", rds_path, "\n")
  obj <- readRDS(rds_path)

  cat("Seurat object:\n")
  print(obj)

  assay <- DefaultAssay(obj)
  cat("\nDefault assay:", assay, "\n")
  counts <- GetAssayData(obj, assay = assay, slot = "counts")

  cat("\nCounts class:", paste(class(counts), collapse = ","), "\n")
  cat("Counts dim  :", paste(dim(counts), collapse = " x "), "\n")
  cat(
    "nnzero      :",
    if (inherits(counts, "dgCMatrix")) length(counts@x) else NA,
    "\n"
  )
  cat("Any NA in rownames? ", anyNA(rownames(counts)), "\n")
  cat("Any duplicated genes? ", any(duplicated(rownames(counts))), "\n")

  genes <- rownames(counts)
  mt_hits <- grep(mt_pattern, genes)
  rb_hits <- grep(rb_pattern, genes)

  cat("\nMT genes matched:", length(mt_hits), "  (pattern:", mt_pattern, ")\n")
  cat("RB genes matched:", length(rb_hits), "  (pattern:", rb_pattern, ")\n")
  if (length(mt_hits) > 0) {
    cat("MT examples:", paste(head(genes[mt_hits], 10), collapse = ", "), "\n")
  }
  if (length(rb_hits) > 0) {
    cat("RB examples:", paste(head(genes[rb_hits], 10), collapse = ", "), "\n")
  }

  # 重新计算一份（不依赖你脚本里当时写入的 percent.mt / percent.rb）
  obj[["percent.mt_recalc"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent.rb_recalc"]] <- PercentageFeatureSet(obj, pattern = rb_pattern)

  cat("\nQC (recalculated) summary:\n")
  print(summary(obj$percent.mt_recalc))
  print(summary(obj$percent.rb_recalc))

  # 看看 meta 里原本有没有/是否全 NA
  qc_cols <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"),
    colnames(obj@meta.data)
  )
  cat("\nQC columns present in meta:", paste(qc_cols, collapse = ", "), "\n")
  if (length(qc_cols) > 0) {
    print(head(obj@meta.data[, qc_cols, drop = FALSE], 5))
  }

  invisible(obj)
}
setwd('/home/h2048/data/R/1210/Polyp/')
# 示例：先看最异常的 GSM8499586
peek_seurat_rds("GSM8499586_seurat.rds")


library(Seurat)
library(Matrix)

obj <- readRDS("GSM8499586_seurat.rds")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

nCount <- Matrix::colSums(counts)
nFeature <- Matrix::colSums(counts > 0)

bad <- which(
  is.na(obj$nCount_RNA) |
    is.na(obj$nFeature_RNA) |
    is.na(obj$percent.mt) |
    is.na(obj$percent.rb) |
    is.na(nCount) |
    is.na(nFeature) |
    nCount == 0 |
    nFeature == 0
)

length(bad)
head(colnames(obj)[bad], 20)
cbind(
  nCount = nCount[bad],
  nFeature = nFeature[bad],
  obj@meta.data[
    bad,
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
  ]
)


library(Seurat)
library(Matrix)

repair_qc_and_drop_bad <- function(
  obj,
  assay = "RNA",
  mt_pattern = "^MT-",
  rb_pattern = "^(RPL|RPS)"
) {
  counts <- GetAssayData(obj, assay = assay, layer = "counts")

  nCount <- Matrix::colSums(counts)
  nFeature <- Matrix::colSums(counts > 0)

  # 写回（覆盖旧的/可能错误的 meta）
  obj$nCount_RNA <- as.numeric(nCount)
  obj$nFeature_RNA <- as.numeric(nFeature)

  mt_genes <- grep(mt_pattern, rownames(counts), value = TRUE)
  rb_genes <- grep(rb_pattern, rownames(counts), value = TRUE)

  mt_sum <- if (length(mt_genes)) {
    Matrix::colSums(counts[mt_genes, , drop = FALSE])
  } else {
    rep(0, ncol(counts))
  }
  rb_sum <- if (length(rb_genes)) {
    Matrix::colSums(counts[rb_genes, , drop = FALSE])
  } else {
    rep(0, ncol(counts))
  }

  denom <- pmax(nCount, 1) # 防止除零导致 NA
  obj$percent.mt <- as.numeric(100 * mt_sum / denom)
  obj$percent.rb <- as.numeric(100 * rb_sum / denom)

  # 关键：先剔除所有 NA 或 0-count（否则后面任何逻辑筛选都可能产生 NA 下标）
  drop <- is.na(obj$nCount_RNA) |
    is.na(obj$nFeature_RNA) |
    is.na(obj$percent.mt) |
    is.na(obj$percent.rb) |
    obj$nCount_RNA <= 0 |
    obj$nFeature_RNA <= 0

  if (any(drop)) {
    obj <- subset(obj, cells = colnames(obj)[!drop])
  }
  obj
}

obj <- readRDS("GSM8499586_seurat.rds")
obj <- repair_qc_and_drop_bad(obj)


library(Matrix)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

# 用原始对象读（剔除前）更好：先 readRDS，再取 counts，再检查
obj0 <- readRDS("GSM8499586_seurat.rds")
counts0 <- GetAssayData(obj0, assay = "RNA", layer = "counts")

# 1) 检查矩阵是否有 NA（如果这里 TRUE，说明读入/构建就异常）
anyNA(counts0)

# 2) 定点检查这两列
cells <- c("CELL670_N1", "CELL959_N1")
Matrix::colSums(counts0[, cells, drop = FALSE])
Matrix::colSums(counts0[, cells, drop = FALSE] > 0)

# 3) 看这两列有哪些非零条目（正常“空细胞”应该是 0 个非零；NA 情况会更诡异）
lapply(cells, \(cc) {
  v <- counts0[, cc]
  list(
    nnz = length(v@x),
    anyNA_x = anyNA(v@x),
    head_x = head(v@x, 10)
  )
})
