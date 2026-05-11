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

  fail <- function(...) {
    stop(..., call. = FALSE)
  }

  assert_scalar_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      fail("`", name, "` must be a single TRUE/FALSE value.")
    }
  }

  assert_nonempty_string <- function(x, name) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
      fail("`", name, "` must be a non-empty string.")
    }
  }

  preview_values <- function(x, n = 5L) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    paste(utils::head(x, n), collapse = ", ")
  }

  get_py_shape <- function(py_obj) {
    tryCatch(
      as.integer(reticulate::py_to_r(py_obj$shape)),
      error = function(e) NULL
    )
  }

  py_list_keys <- function(mapping) {
    tryCatch(
      as.character(reticulate::py_to_r(py_builtins$list(mapping$keys()))),
      error = function(e) character(0)
    )
  }

  assert_valid_dimnames <- function(x, label) {
    if (anyNA(x) || any(!nzchar(trimws(x)))) {
      fail(label, " contains NA/empty values.")
    }
    dup_n <- sum(duplicated(x))
    if (dup_n > 0L) {
      dup_examples <- preview_values(x[duplicated(x)])
      fail(
        label,
        " contains duplicated values (",
        dup_n,
        " duplicates). Examples: ",
        dup_examples
      )
    }
  }

  normalize_reduction_name <- function(key) {
    fixed_key <- gsub("^(X_)+", "", key)
    fixed_key <- trimws(fixed_key)
    fixed_key <- gsub("[[:space:]]+", "_", fixed_key)
    fixed_key <- gsub("[^A-Za-z0-9_]", "_", fixed_key)
    fixed_key <- gsub("_+", "_", fixed_key)
    fixed_key <- gsub("^_+|_+$", "", fixed_key)

    if (!nzchar(fixed_key)) {
      fixed_key <- "reduction"
    }
    if (!grepl("^[A-Za-z]", fixed_key)) {
      fixed_key <- paste0("reduction_", fixed_key)
    }

    fixed_key
  }

  normalize_dimred_key <- function(name) {
    dimred_key <- gsub("[^A-Za-z0-9]", "", name)
    if (!nzchar(dimred_key)) {
      dimred_key <- "reduction"
    }
    if (!grepl("^[A-Za-z]", dimred_key)) {
      dimred_key <- paste0("Reduction", dimred_key)
    }
    paste0(dimred_key, "_")
  }

  assert_nonempty_string(h5ad_path, "h5ad_path")
  assert_nonempty_string(assay, "assay")
  assert_scalar_flag(prefer_raw, "prefer_raw")
  assert_scalar_flag(prefer_layer_counts, "prefer_layer_counts")
  assert_scalar_flag(validate_counts, "validate_counts")
  assert_scalar_flag(debug, "debug")

  if (!file.exists(h5ad_path)) {
    fail("H5AD file does not exist: ", h5ad_path)
  }
  h5ad_path <- normalizePath(h5ad_path, winslash = "/", mustWork = TRUE)

  if (debug) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("GetSeurat(): reading h5ad -> Seurat (", h5ad_path, ")\n", sep = "")
    cat(strrep("=", 70), "\n\n", sep = "")
  }

  anndata <- reticulate::import("anndata", convert = FALSE)
  py_builtins <- reticulate::import_builtins(convert = FALSE)
  adata <- anndata$read_h5ad(h5ad_path)

  # ---------------------------------------------------------------------------
  # 1. 选择使用的表达矩阵：raw / layers['counts'] / X
  # ---------------------------------------------------------------------------
  use_raw <- FALSE
  use_layer_counts <- FALSE
  counts_py <- NULL
  counts_source <- NULL
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
          counts_source <- "adata$raw$X"
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

    layer_keys <- py_list_keys(adata$layers)

    counts_layer <- tryCatch(adata$layers[["counts"]], error = function(e) NULL)
    raw_counts_layer <- tryCatch(adata$layers[["raw_counts"]], error = function(e) NULL)

    if (!is.null(counts_layer)) {
      counts_py <- counts_layer
      counts_source <- "adata$layers['counts']"
      var_names <- as.character(
        reticulate::py_to_r(adata$var$index$to_list())
      )
      use_layer_counts <- TRUE
      message("Using counts matrix from adata$layers['counts']")
    } else if (!is.null(raw_counts_layer)) {
      counts_py <- raw_counts_layer
      counts_source <- "adata$layers['raw_counts']"
      var_names <- as.character(
        reticulate::py_to_r(adata$var$index$to_list())
      )
      use_layer_counts <- TRUE
      message("Using counts matrix from adata$layers['raw_counts']")
    } else if ("counts" %in% layer_keys) {
      counts_py <- adata$layers[["counts"]]
      counts_source <- "adata$layers['counts']"
      var_names <- as.character(
        reticulate::py_to_r(adata$var$index$to_list())
      )
      use_layer_counts <- TRUE
      message("Using counts matrix from adata$layers['counts']")
    } else if ("raw_counts" %in% layer_keys) {
      counts_py <- adata$layers[["raw_counts"]]
      counts_source <- "adata$layers['raw_counts']"
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
    counts_source <- "adata$X"
    var_names <- as.character(
      reticulate::py_to_r(adata$var$index$to_list())
    )
  }

  ## 1.4 兜底：如果真的是 NULL，就直接报错（防止后面 as.matrix(NULL)）
  if (is.null(counts_py)) {
    fail(
      "counts_py is NULL: no usable matrix found in adata$raw$X, ",
      "adata$layers['counts'/'raw_counts'], or adata$X.\n",
      "请在 Python 侧检查 AnnData 对象的 raw / layers / X 是否存在。"
    )
  }

  counts_shape <- get_py_shape(counts_py)
  if (!is.null(counts_shape)) {
    if (length(counts_shape) != 2L) {
      fail(
        "Selected matrix from ",
        counts_source,
        " is not 2-dimensional. shape=",
        paste(counts_shape, collapse = " x ")
      )
    }
    if (any(counts_shape <= 0L)) {
      fail(
        "Selected matrix from ",
        counts_source,
        " has non-positive shape: ",
        paste(counts_shape, collapse = " x ")
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 2. Python -> R：确保得到 dgCMatrix
  # ---------------------------------------------------------------------------
  counts_r <- reticulate::py_to_r(counts_py)

  if (is.null(counts_r)) {
    fail("Failed to convert ", counts_source, " from Python to R.")
  }

  if (inherits(counts_r, "Matrix")) {
    if (!inherits(counts_r, "dgCMatrix")) {
      counts_r <- methods::as(counts_r, "TsparseMatrix")
      counts_r <- methods::as(counts_r, "CsparseMatrix")
    }
  } else {
    counts_r <- methods::as(as.matrix(counts_r), "dgCMatrix")
  }

  # 转置：AnnData 是 cells x genes → Seurat 要 genes x cells
  counts_r <- Matrix::t(counts_r)
  if (!inherits(counts_r, "dgCMatrix")) {
    counts_r <- methods::as(counts_r, "TsparseMatrix")
    counts_r <- methods::as(counts_r, "CsparseMatrix")
  }
  if (!inherits(counts_r, "dgCMatrix")) {
    fail("Internal error: counts matrix could not be normalized to dgCMatrix.")
  }
  if (nrow(counts_r) == 0L || ncol(counts_r) == 0L) {
    fail("Converted matrix from ", counts_source, " is empty after transpose.")
  }
  if (length(counts_r@x) > 0L && (anyNA(counts_r@x) || any(!is.finite(counts_r@x)))) {
    fail("Detected NA/NaN/Inf values in matrix from ", counts_source, ".")
  }
  if ((use_raw || use_layer_counts) && length(counts_r@x) > 0L && any(counts_r@x < 0)) {
    fail(
      "Detected negative values in matrix from ",
      counts_source,
      ", which should contain raw counts."
    )
  }

  # obs / var names
  # （和你原代码保持一致，用 obs$index / var$index）
  obs_names <- as.character(
    reticulate::py_to_r(adata$obs$index$to_list())
  )
  assert_valid_dimnames(obs_names, "Cell names (adata$obs index)")

  if (is.null(var_names)) {
    var_names <- as.character(
      reticulate::py_to_r(adata$var$index$to_list())
    )
  }
  assert_valid_dimnames(var_names, "Gene names")

  if (!is.null(counts_shape) && length(var_names) != counts_shape[[2]]) {
    fail(
      "Selected matrix from ",
      counts_source,
      " has ",
      counts_shape[[2]],
      " features, but extracted ",
      length(var_names),
      " gene names."
    )
  }
  if (!is.null(counts_shape) && length(obs_names) != counts_shape[[1]]) {
    fail(
      "Selected matrix from ",
      counts_source,
      " has ",
      counts_shape[[1]],
      " cells, but extracted ",
      length(obs_names),
      " obs names."
    )
  }

  if (length(var_names) != nrow(counts_r)) {
    fail(
      "Number of gene names (",
      length(var_names),
      ") does not match number of rows in counts (",
      nrow(counts_r),
      ")"
    )
  }
  if (length(obs_names) != ncol(counts_r)) {
    fail(
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
    cat("  Selected matrix source: ", counts_source, "\n", sep = "")
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
      fail(
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
  meta_data <- as.data.frame(meta_data, stringsAsFactors = FALSE)
  if (nrow(meta_data) != length(obs_names)) {
    fail(
      "meta.data row count (",
      nrow(meta_data),
      ") does not match number of cells (",
      length(obs_names),
      ")."
    )
  }
  if (anyDuplicated(colnames(meta_data)) > 0L) {
    fail(
      "meta.data contains duplicated column names. Examples: ",
      preview_values(colnames(meta_data)[duplicated(colnames(meta_data))])
    )
  }
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

  obsm_keys <- py_list_keys(adata$obsm)
  if (debug && length(obsm_keys) == 0L) {
    cat("  No readable obsm keys found via Python list(keys()).\n")
  }
  obsm_keys <- unique(obsm_keys[!is.na(obsm_keys) & nzchar(obsm_keys)])

  if (length(obsm_keys) > 0) {
    for (key in obsm_keys) {
      embedding <- tryCatch(
        reticulate::py_to_r(adata$obsm[[key]]),
        error = function(e) {
          if (debug) {
            cat(
              "  Skip obsm key ",
              key,
              " (read failed: ",
              conditionMessage(e),
              ")\n",
              sep = ""
            )
          }
          NULL
        }
      )
      if (is.null(embedding)) {
        next
      }

      if (is.vector(embedding) && !is.list(embedding)) {
        embedding <- matrix(embedding, ncol = 1L)
      } else {
        embedding <- tryCatch(
          as.matrix(embedding),
          error = function(e) NULL
        )
      }

      if (is.null(embedding) || length(dim(embedding)) != 2L) {
        if (debug) {
          cat(
            "  Skip obsm key ",
            key,
            " (not a 2D matrix after conversion)\n",
            sep = ""
          )
        }
        next
      }

      if (nrow(embedding) != ncol(seurat_obj) && ncol(embedding) == ncol(seurat_obj)) {
        if (debug) {
          cat("  Transposing obsm key ", key, " to align cells x dims\n", sep = "")
        }
        embedding <- t(embedding)
      }

      if (nrow(embedding) != ncol(seurat_obj)) {
        if (debug) {
          cat(
            "  Skip obsm key ",
            key,
            " (nrow != ncol(seurat_obj) after alignment)\n",
            sep = ""
          )
        }
        next
      }
      if (ncol(embedding) < 1L) {
        if (debug) {
          cat("  Skip obsm key ", key, " (no embedding dimensions)\n", sep = "")
        }
        next
      }

      embedding_numeric <- suppressWarnings(
        matrix(
          as.numeric(embedding),
          nrow = nrow(embedding),
          ncol = ncol(embedding),
          dimnames = dimnames(embedding)
        )
      )
      if (anyNA(embedding_numeric) || any(!is.finite(embedding_numeric))) {
        if (debug) {
          cat(
            "  Skip obsm key ",
            key,
            " (contains NA/NaN/Inf after numeric coercion)\n",
            sep = ""
          )
        }
        next
      }
      embedding <- embedding_numeric

      fixed_key <- normalize_reduction_name(key)
      if (fixed_key %in% names(seurat_obj@reductions)) {
        if (debug) {
          cat(
            "  Skip obsm key ",
            key,
            " (name collision after normalization: ",
            fixed_key,
            ")\n",
            sep = ""
          )
        }
        next
      }

      dimred_key <- normalize_dimred_key(fixed_key)
      rownames(embedding) <- colnames(seurat_obj)
      colnames(embedding) <- paste0(dimred_key, seq_len(ncol(embedding)))

      seurat_obj[[fixed_key]] <- tryCatch(
        Seurat::CreateDimReducObject(
          embeddings = embedding,
          key = dimred_key,
          assay = assay
        ),
        error = function(e) {
          fail(
            "Failed to create reduction '",
            fixed_key,
            "' from obsm key '",
            key,
            "': ",
            conditionMessage(e)
          )
        }
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
