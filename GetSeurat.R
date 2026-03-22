GetSeurat <- function(
  h5ad_path,
  assay = "RNA",
  prefer_raw = TRUE,
  prefer_layer_counts = TRUE,
  validate_counts = TRUE,
  debug = FALSE
) {
  library(reticulate)
  library(Seurat)
  library(Matrix)

  if (debug) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("GetSeurat(): reading h5ad -> Seurat (", h5ad_path, ")\n", sep = "")
    cat(strrep("=", 70), "\n\n", sep = "")
  }

  anndata <- reticulate::import("anndata", convert = FALSE)
  adata <- anndata$read_h5ad(h5ad_path)

  # ---------------------------------------------------------------------------
  # 1. 选择使用的表达矩阵：raw / layers['counts'] / X
  # ---------------------------------------------------------------------------
  use_raw <- FALSE
  use_layer_counts <- FALSE
  counts_py <- NULL
  var_names <- NULL

  ## 1.1 优先 adata$raw$X
  if (prefer_raw && !reticulate::py_is_null_xptr(adata$raw)) {
    if (debug) {
      cat("  Checking adata$raw$X ...\n")
    }
    tryCatch(
      {
        if (!is.null(adata$raw$X)) {
          counts_py <- adata$raw$X
          # 优先从 raw.var 拿基因名，退回 var.index
          if (!is.null(adata$raw$var)) {
            var_names <- as.character(
              reticulate::py_to_r(adata$raw$var$index$to_list())
            )
          } else {
            var_names <- as.character(
              reticulate::py_to_r(adata$var$index$to_list())
            )
          }
          use_raw <- TRUE
          message("Using raw counts matrix from adata$raw$X")
        } else {
          if (debug) cat("  adata$raw$X is NULL, skip.\n")
        }
      },
      error = function(e) {
        if (debug) {
          cat(
            "  Failed to use adata$raw$X: ",
            conditionMessage(e),
            "\n",
            sep = ""
          )
        }
        use_raw <<- FALSE
      }
    )
  }

  ## 1.2 其次 adata.layers['counts'] 或 adata.layers['raw_counts']
  if (!use_raw && prefer_layer_counts) {
    if (debug) {
      cat("  Checking adata$layers['counts'] / ['raw_counts'] ...\n")
    }
    layer_keys <- as.character(reticulate::py_to_r(adata$layers$keys()))

    if ("counts" %in% layer_keys) {
      # 更稳的写法是用 [[ ]] 访问，而不是 get()
      counts_py <- adata$layers[["counts"]]
      var_names <- as.character(
        reticulate::py_to_r(adata$var$index$to_list())
      )
      use_layer_counts <- TRUE
      message("Using counts matrix from adata$layers['counts']")
    } else if ("raw_counts" %in% layer_keys) {
      counts_py <- adata$layers[["raw_counts"]]
      var_names <- as.character(
        reticulate::py_to_r(adata$var$index$to_list())
      )
      use_layer_counts <- TRUE
      message("Using counts matrix from adata$layers['raw_counts']")
    } else {
      if (debug) {
        cat(
          "  No 'counts' or 'raw_counts' in adata$layers. Keys: ",
          paste(layer_keys, collapse = ", "),
          "\n",
          sep = ""
        )
      }
    }
  }

  ## 1.3 最后强制 fallback 到 adata.X（以 counts_py 是否为 NULL 为准）
  if (is.null(counts_py)) {
    message("Using matrix from adata$X (likely normalized)")
    counts_py <- adata$X
    var_names <- as.character(
      reticulate::py_to_r(adata$var$index$to_list())
    )
  }

  ## 1.4 兜底：如果真的是 NULL，就直接报错（防止后面 as.matrix(NULL)）
  if (is.null(counts_py)) {
    stop(
      "counts_py is NULL: no usable matrix found in adata$raw$X, ",
      "adata$layers['counts'/'raw_counts'], or adata$X.\n",
      "请在 Python 侧检查 AnnData 对象的 raw / layers / X 是否存在。"
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Python -> R：确保得到 dgCMatrix
  # ---------------------------------------------------------------------------
  counts_r <- reticulate::py_to_r(counts_py)

  if (!(inherits(counts_r, "dgCMatrix") || inherits(counts_r, "dgRMatrix"))) {
    counts_r <- as(as.matrix(counts_r), "dgCMatrix")
  } else if (inherits(counts_r, "dgRMatrix")) {
    counts_r <- as(counts_r, "dgCMatrix")
  }

  # 转置：AnnData 是 cells x genes → Seurat 要 genes x cells
  counts_r <- Matrix::t(counts_r)

  # obs / var names
  # （和你原代码保持一致，用 obs$index / var$index）
  obs_names <- as.character(
    reticulate::py_to_r(adata$obs$index$to_list())
  )

  if (is.null(var_names)) {
    var_names <- as.character(
      reticulate::py_to_r(adata$var$index$to_list())
    )
  }

  if (length(var_names) != nrow(counts_r)) {
    stop(
      "Number of gene names (",
      length(var_names),
      ") does not match number of rows in counts (",
      nrow(counts_r),
      ")"
    )
  }
  if (length(obs_names) != ncol(counts_r)) {
    stop(
      "Number of cell names (",
      length(obs_names),
      ") does not match number of columns in counts (",
      ncol(counts_r),
      ")"
    )
  }

  rownames(counts_r) <- var_names
  colnames(counts_r) <- obs_names

  if (debug) {
    cat(sprintf(
      "  Counts matrix: %d genes x %d cells\n",
      nrow(counts_r),
      ncol(counts_r)
    ))
    cat(
      "  Example genes: ",
      paste(head(rownames(counts_r), 5), collapse = ", "),
      "\n"
    )
    cat(
      "  Example cells: ",
      paste(head(colnames(counts_r), 5), collapse = ", "),
      "\n\n"
    )
  }

  # ---------------------------------------------------------------------------
  # 3. 可选：验证是否像 raw counts（和 GetH5ad 对齐）
  # ---------------------------------------------------------------------------
  if (validate_counts) {
    if (debug) {
      cat(
        strrep("-", 70),
        "\nRaw counts validation\n",
        strrep("-", 70),
        "\n",
        sep = ""
      )
    }

    nc <- ncol(counts_r)
    ng <- nrow(counts_r)
    sample_size_cells <- min(3000, nc)
    sample_size_genes <- min(3000, ng)

    sample_counts <- counts_r[
      1:sample_size_genes,
      1:sample_size_cells,
      drop = FALSE
    ]
    sample_counts <- as.matrix(sample_counts)

    max_val <- max(sample_counts)
    mean_val <- mean(sample_counts)
    is_integer_like <- all(sample_counts == round(sample_counts))

    if (debug) {
      cat(sprintf(
        "  Sample: %d genes x %d cells\n",
        sample_size_genes,
        sample_size_cells
      ))
      cat(sprintf("  Max value  : %.2f\n", max_val))
      cat(sprintf("  Mean value : %.4f\n", mean_val))
      cat(sprintf("  Integer-like: %s\n", is_integer_like))
    }

    if (max_val < 20 && mean_val < 2 && !is_integer_like) {
      stop(
        "\n❌ VALIDATION FAILED: matrix looks log-transformed instead of raw counts.\n",
        "   Max value  : ",
        round(max_val, 2),
        "\n",
        "   Mean value : ",
        round(mean_val, 4),
        "\n",
        "   Integer-like: ",
        is_integer_like,
        "\n\n",
        "下游 Monocle3/DE 等分析通常要求 raw counts。\n",
        "请在 Python 侧确认 adata.raw 或 layers['counts'] 存放的是原始 UMI。\n",
        "如已确认是 raw counts，可调用 GetSeurat(..., validate_counts = FALSE)。\n"
      )
    }

    if (is_integer_like && max_val < 10) {
      warning(
        "\n⚠️  WARNING: max count value is very low (",
        round(max_val, 2),
        ").\n",
        "   可能原因：\n",
        "   1) 测序深度极低；\n",
        "   2) 部分下采样 / 过滤后矩阵；\n",
        "   3) 不是原始 UMI 计数。\n",
        "如这是预期，请忽略此警告。\n"
      )
    }

    if (is_integer_like && max_val >= 10 && debug) {
      cat("\n✅ VALIDATION PASSED: counts look like raw counts.\n\n")
    }
  }

  # ---------------------------------------------------------------------------
  # 4. meta.data
  # ---------------------------------------------------------------------------
  meta_data <- reticulate::py_to_r(adata$obs)
  # 确保 meta_data 行名和细胞名一致
  rownames(meta_data) <- obs_names

  if (debug) {
    cat("  meta.data columns:", ncol(meta_data), "\n")
  }

  # ---------------------------------------------------------------------------
  # 5. 创建 Seurat 对象
  # ---------------------------------------------------------------------------
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts_r,
    meta.data = meta_data,
    assay = assay
  )

  # ---------------------------------------------------------------------------
  # 6. obsm -> Seurat reductions
  # ---------------------------------------------------------------------------
  if (debug) {
    cat("\nImporting obsm as reductions ...\n")
  }

  obsm_dict <- reticulate::py_to_r(
    reticulate::py_call(adata$obsm$as_dict)
  )

  if (length(obsm_dict) > 0) {
    for (key in names(obsm_dict)) {
      embedding <- obsm_dict[[key]]
      embedding <- reticulate::py_to_r(embedding)
      embedding <- as.matrix(embedding)

      # 检查行数是否匹配细胞数
      if (nrow(embedding) != ncol(seurat_obj)) {
        if (debug) {
          cat(
            "  Skip obsm key ",
            key,
            " (nrow != ncol(seurat_obj))\n",
            sep = ""
          )
        }
        next
      }

      rownames(embedding) <- colnames(seurat_obj)
      colnames(embedding) <- paste0(key, "_", seq_len(ncol(embedding)))

      fixed_key <- gsub("^X_", "", key) # X_umap -> umap
      # Seurat 的 DimReduc key 最后一般加 "_"，但这里你之前是直接用名字，就保持不动
      seurat_obj[[fixed_key]] <- Seurat::CreateDimReducObject(
        embeddings = embedding,
        key = fixed_key,
        assay = assay
      )

      if (debug) {
        cat(
          "  Added reduction: ",
          fixed_key,
          " (from obsm key ",
          key,
          ", dims = ",
          ncol(embedding),
          ")\n",
          sep = ""
        )
      }
    }
  } else if (debug) {
    cat("  No obsm entries found.\n")
  }

  if (debug) {
    cat("\nGetSeurat() finished.\n")
    cat(
      "  Seurat object: ",
      nrow(seurat_obj),
      " genes x ",
      ncol(seurat_obj),
      " cells\n",
      sep = ""
    )
  }

  return(seurat_obj)
}

# 示例用法：
# scnew <- GetSeurat("sce.h5ad")
# scnew <- GetSeurat("sce.h5ad", debug = TRUE)
# scnew <- GetSeurat("sce.h5ad", validate_counts = FALSE)  # 已确认是 raw，可跳过检查
