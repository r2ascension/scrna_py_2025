#!/usr/bin/env Rscript
# ==============================================================================
# hdWGCNA + CoVarNet Helper Functions v1.1
# ==============================================================================
#
# Purpose:
#   Production-ready reusable wrappers for gene co-expression network analysis
#   in single-cell data (Seurat v5 + hdWGCNA + CoVarNet / igraph).
#
# Two workflows:
#   [A] hdWGCNA   -- metacell-based WGCNA on Seurat objects
#   [B] CoVarNet  -- cell-type-stratified gene co-variation networks
#
# Architecture (v1.1 changes from v1.0):
#   P0  Counts guard now uses Seurat v5 Layers() API -- same pattern as template
#   P0  Dead params n_hvg / target_use wired up correctly or removed
#   P0  NormalizeMetacells() now receives target.use (dot notation)
#   P0  ident.group no longer hardcoded; derived from group_by[1]
#   P0  cor = +-1 NaN guard: clamp before t-stat computation
#   P1  One-shot runner split into global setup (once) + per-celltype loop
#   P1  Usage example uses only defined variables, directly runnable
#   P2  Sparse-aware row variance in covarnet_compute_cor()
#   P2  All file-writing helpers call dir.create() internally
#
# Dependencies:
#   devtools::install_github("smorabit/hdWGCNA", ref="dev")
#   BiocManager::install("WGCNA")
#   install.packages(c("igraph", "ggraph", "tidygraph", "pheatmap",
#                      "tidyverse", "cowplot", "patchwork", "viridis"))
#
# Author:  r2end
# Date:    2026-04-23
# Version: v1.1
# ==============================================================================


# ==============================================================================
# 0. Thread control
# ==============================================================================

Sys.setenv(
  OMP_NUM_THREADS      = "4",
  MKL_NUM_THREADS      = "4",
  OPENBLAS_NUM_THREADS = "4"
)

options(stringsAsFactors = FALSE)

PA_GENE_EXCLUSION_HELPER_PATH_20260505_V1 <- "/home/h2048/script/R/program_gene_exclusion_helper_20260505_v1.R"
if (!exists("pa_apply_gene_exclusion_to_seurat", mode = "function")) {
  source(PA_GENE_EXCLUSION_HELPER_PATH_20260505_V1)
}


# ##############################################################################
# >>>>>>>>>>  SECTION 1: GLOBAL CONFIG (EDIT THIS)  <<<<<<<<<<<<<<<<<<<<<<<<<<<
# ##############################################################################

# ----- Input / output -----
H5AD_PATH  <- "/home/h2048/data/py/MMDD/your_adata.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0423/coexpression_v1"

# ----- Seurat metadata column names -----
CELLTYPE_COL  <- "cell_type_L2"   # column used to split cell types
SAMPLE_COL    <- "sample_id"      # donor / sample ID
CONDITION_COL <- "disease_status" # CRSwNP / Healthy / etc.
TISSUE_COL    <- "tissue"         # nasal / sinus / lung / etc.

# ----- hdWGCNA global -----
WGCNA_NAME      <- "hdWGCNA"
# Gene selection mode: "fraction" uses GENE_FRACTION; "n_top" uses GENE_N_TOP
GENE_SELECT_MODE <- "fraction"
GENE_FRACTION   <- 0.05     # keep genes detected in >= X fraction of cells
GENE_N_TOP      <- 3000     # used only when GENE_SELECT_MODE == "n_top"
N_METACELL_K    <- 25       # k for KNN metacell aggregation
METACELL_TARGET <- 5e4      # library-size target for NormalizeMetacells()
SOFT_POWER      <- NULL     # NULL = auto-select (R^2 >= 0.80)
NETWORK_TYPE    <- "signed hybrid"
TOMTYPE         <- "unsigned"
N_HUB_GENES    <- 15
N_CORES_WGCNA  <- 8

# ----- CoVarNet -----
COVAR_CELLTYPES  <- NULL    # NULL = all cell types; or c("Naive_B", "Memory_B")
COVAR_MIN_CELLS  <- 50      # skip cell type if n_cells < threshold
COVAR_N_GENES    <- 500     # top variable genes per cell type
COVAR_COR_METHOD <- "pearson"
COVAR_COR_THR    <- 0.30    # |r| edge inclusion threshold
COVAR_PVAL_THR   <- 0.05    # FDR threshold for edges

# ----- Plotting -----
FIG_WIDTH  <- 10
FIG_HEIGHT <- 8
DPI        <- 300


# ==============================================================================
# 1. Library loader
# ==============================================================================

load_coexpr_libs <- function() {
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat("Loading co-expression libraries\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")

  pkgs_required <- c("Seurat", "WGCNA", "igraph",
                     "ggraph", "tidygraph", "tidyverse",
                     "cowplot", "patchwork", "viridis", "pheatmap")
  pkgs_optional <- c("hdWGCNA", "Matrix")

  for (pkg in pkgs_required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("[ERROR] Required package not found: %s", pkg))
    }
    library(pkg, character.only = TRUE, quietly = TRUE)
    cat(sprintf("[OK]   %-18s loaded\n", pkg))
  }

  for (pkg in pkgs_optional) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE, quietly = TRUE)
      cat(sprintf("[OK]   %-18s loaded\n", pkg))
    } else {
      cat(sprintf("[WARN] %-18s not found -- some steps may fail\n", pkg))
    }
  }

  suppressMessages(enableWGCNAThreads(nThreads = N_CORES_WGCNA))
  cat(sprintf("[OK] WGCNA threads = %d\n\n", N_CORES_WGCNA))
  invisible(NULL)
}


# ==============================================================================
# 2. Input validation helpers
# ==============================================================================

# P0 fix: use Seurat v5 Layers() API, consistent with project template
check_counts_layer <- function(seurat_obj) {
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("[ERROR] 'RNA' assay not found in Seurat object.")
  }
  if (!"counts" %in% Layers(seurat_obj[["RNA"]])) {
    stop(paste(
      "[ERROR] RNA assay missing 'counts' layer.",
      "Ensure GetSeurat() or your loader preserves raw counts.",
      "Check: Layers(obj[[\"RNA\"]])"
    ))
  }
  invisible(TRUE)
}

hdwgcna_ensure_scaled_data <- function(seurat_obj,
                                       features   = NULL,
                                       assay      = "RNA",
                                       wgcna_name = WGCNA_NAME) {
  if (!assay %in% names(seurat_obj@assays)) {
    stop(sprintf("[ERROR] Assay '%s' not found in Seurat object.", assay))
  }

  if (is.null(features)) {
    features <- tryCatch(
      GetWGCNAGenes(seurat_obj, wgcna_name = wgcna_name),
      error = function(e) NULL
    )
  }

  features <- unique(stats::na.omit(as.character(features)))
  features <- intersect(features, rownames(seurat_obj))
  if (length(features) == 0L) features <- rownames(seurat_obj)

  layers <- Layers(seurat_obj[[assay]])
  needs_scaling <- TRUE

  if ("scale.data" %in% layers) {
    scaled_features <- tryCatch(
      rownames(LayerData(seurat_obj[[assay]], layer = "scale.data")),
      error = function(e) character()
    )
    if (length(scaled_features) > 0L && all(features %in% scaled_features)) {
      needs_scaling <- FALSE
    }
  }

  if (!needs_scaling) {
    return(seurat_obj)
  }

  cat(sprintf("[INFO] Running ScaleData on %d features before ModuleEigengenes\n",
              length(features)))

  tryCatch(
    ScaleData(seurat_obj, features = features, verbose = FALSE),
    error = function(e) {
      cat(sprintf("[ERROR] ScaleData: %s\n", e$message))
      stop(e)
    }
  )
}

# Ensure output directory exists (called internally by all file-writing helpers)
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

hdwgcna_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

hdwgcna_write_json <- function(x, path) {
  ensure_dir(dirname(path))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  }
  saveRDS(x, file = sub("\\.json$", ".rds", path))
  invisible(path)
}

hdwgcna_read_json <- function(path, default = NULL) {
  if (file.exists(path) && requireNamespace("jsonlite", quietly = TRUE)) {
    return(tryCatch(
      jsonlite::read_json(path, simplifyVector = TRUE),
      error = function(e) default
    ))
  }
  rds_path <- sub("\\.json$", ".rds", path)
  if (file.exists(rds_path)) {
    return(tryCatch(readRDS(rds_path), error = function(e) default))
  }
  default
}

hdwgcna_celltype_status_path <- function(ct_dir) {
  file.path(ct_dir, "hdwgcna_celltype_status.json")
}

hdwgcna_write_celltype_status <- function(record, ct_dir) {
  record$updated_at <- hdwgcna_timestamp()
  if (is.null(record$output_dir)) record$output_dir <- ct_dir
  hdwgcna_write_json(record, hdwgcna_celltype_status_path(ct_dir))
  invisible(record)
}

hdwgcna_read_celltype_status <- function(ct_dir) {
  record <- hdwgcna_read_json(hdwgcna_celltype_status_path(ct_dir), default = NULL)
  if (!is.list(record) || is.null(record$status)) return(NULL)
  record$status <- as.character(record$status)[1]
  record
}

hdwgcna_existing_celltype_record <- function(ct, ct_dir, resume_skip_statuses = c("ok", "no_modules")) {
  record <- hdwgcna_read_celltype_status(ct_dir)
  if (!is.null(record) && record$status %in% resume_skip_statuses) {
    record$celltype <- if (is.null(record$celltype)) ct else as.character(record$celltype)[1]
    record$output_dir <- ct_dir
    record$skipped_existing <- TRUE
    return(record)
  }

  membership_csv <- file.path(ct_dir, "hdwgcna_module_membership.csv")
  if ("ok" %in% resume_skip_statuses && file.exists(membership_csv)) {
    mods <- tryCatch(utils::read.csv(membership_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    module_ids <- if (!is.null(mods) && "module" %in% colnames(mods)) {
      setdiff(unique(as.character(mods$module)), "grey")
    } else {
      character()
    }
    record <- list(
      celltype = ct,
      status = "ok",
      output_dir = ct_dir,
      module_ids = module_ids,
      n_modules = length(module_ids),
      membership_csv = membership_csv,
      inferred_from = "hdwgcna_module_membership.csv",
      skipped_existing = TRUE
    )
    hdwgcna_write_celltype_status(record, ct_dir)
    return(record)
  }

  NULL
}

hdwgcna_write_partial_results <- function(results, out_dir) {
  if (length(results) == 0L) return(invisible(NULL))
  saveRDS(results, file = file.path(out_dir, "hdwgcna_celltype_results_partial.rds"))
  invisible(NULL)
}

hdwgcna_is_timeout_error <- function(e) {
  inherits(e, "TimeoutException") || grepl("timeout|time limit", conditionMessage(e), ignore.case = TRUE)
}

hdwgcna_with_timeout <- function(expr, timeout_sec = NULL) {
  timeout_sec <- suppressWarnings(as.numeric(timeout_sec)[1])
  if (length(timeout_sec) == 0L || is.na(timeout_sec) || timeout_sec <= 0) {
    return(force(expr))
  }
  if (requireNamespace("R.utils", quietly = TRUE)) {
    return(R.utils::withTimeout(force(expr), timeout = timeout_sec, onTimeout = "error"))
  }
  old_limits <- setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  force(expr)
}


# ==============================================================================
# 3. hdWGCNA -- Setup (run ONCE per Seurat object)
# ==============================================================================

# ----- 3.1 SetupForWGCNA -----
# P0 fix: both gene_select modes ("fraction" / "n_top") are wired up correctly
hdwgcna_setup <- function(seurat_obj,
                           wgcna_name       = WGCNA_NAME,
                           gene_select_mode = GENE_SELECT_MODE,
                           gene_fraction    = GENE_FRACTION,
                           gene_n_top       = GENE_N_TOP) {
  cat("\n--- hdWGCNA: SetupForWGCNA ---\n")
  check_counts_layer(seurat_obj)

  if (!gene_select_mode %in% c("fraction", "n_top")) {
    stop(sprintf("[ERROR] gene_select_mode must be 'fraction' or 'n_top', got '%s'",
                 gene_select_mode))
  }

  cat(sprintf("[INFO] Mode: %s", gene_select_mode))
  if (gene_select_mode == "fraction") {
    cat(sprintf(", fraction = %.2f\n", gene_fraction))
    obj <- tryCatch(
      SetupForWGCNA(
        seurat_obj,
        gene_select = "fraction",
        fraction    = gene_fraction,
        wgcna_name  = wgcna_name
      ),
      error = function(e) {
        cat(sprintf("[ERROR] SetupForWGCNA (fraction): %s\n", e$message))
        stop(e)
      }
    )
  } else {
    # n_top mode: select top variable genes by expression variability
    cat(sprintf(", n_top = %d\n", gene_n_top))
    obj <- tryCatch(
      SetupForWGCNA(
        seurat_obj,
        gene_select = "variable",
        nfeatures   = gene_n_top,
        wgcna_name  = wgcna_name
      ),
      error = function(e) {
        cat(sprintf("[ERROR] SetupForWGCNA (n_top): %s\n", e$message))
        stop(e)
      }
    )
  }

  n_sel <- length(GetWGCNAGenes(obj, wgcna_name))
  cat(sprintf("[OK] SetupForWGCNA done -- genes selected: %d\n", n_sel))

  if (n_sel < 200) {
    cat(sprintf(
      "[WARN] Only %d genes selected. Consider lowering gene_fraction or",
      n_sel
    ))
    cat(" increasing gene_n_top.\n")
  }

  obj
}


# ----- 3.2 MetacellsByGroups + NormalizeMetacells (run ONCE) -----
# P0 fix: target.use now passed to NormalizeMetacells()
# P0 fix: ident.group derived from group_by[1] instead of hardcoded CELLTYPE_COL
hdwgcna_metacells <- function(seurat_obj,
                               group_by        = c(CELLTYPE_COL, SAMPLE_COL),
                               metacell_k      = N_METACELL_K,
                               metacell_target = METACELL_TARGET,
                               wgcna_name      = WGCNA_NAME) {
  cat("\n--- hdWGCNA: MetacellsByGroups ---\n")
  cat(sprintf("[INFO] k = %d, group_by = %s, ident.group = '%s'\n",
              metacell_k,
              paste(group_by, collapse = " + "),
              group_by[1]))

  obj <- tryCatch(
    MetacellsByGroups(
      seurat_obj,
      group.by    = group_by,
      k           = metacell_k,
      ident.group = group_by[1],   # P0: was hardcoded CELLTYPE_COL
      wgcna_name  = wgcna_name
    ),
    error = function(e) {
      cat(sprintf("[ERROR] MetacellsByGroups: %s\n", e$message))
      stop(e)
    }
  )

  cat(sprintf("[OK] Metacells created\n"))
  cat(sprintf("[INFO] Normalizing metacells (target.use = %.0e)\n",
              metacell_target))

  obj <- tryCatch(
    NormalizeMetacells(
      obj,
      target.use = metacell_target,   # P0: was missing
      wgcna_name = wgcna_name
    ),
    error = function(e) {
      cat(sprintf("[ERROR] NormalizeMetacells: %s\n", e$message))
      stop(e)
    }
  )

  cat(sprintf("[OK] Metacells normalized\n"))
  obj
}


# ==============================================================================
# 4. hdWGCNA -- Per-cell-type steps (loop these)
# ==============================================================================

# ----- 4.1 SetDatExpr -----
hdwgcna_set_datexpr <- function(seurat_obj,
                                 group_name,
                                 group_by   = CELLTYPE_COL,
                                 wgcna_name = WGCNA_NAME) {
  cat(sprintf("\n--- hdWGCNA: SetDatExpr [%s] ---\n", group_name))

  obj <- tryCatch(
    SetDatExpr(
      seurat_obj,
      group.by   = group_by,
      group_name = group_name,
      wgcna_name = wgcna_name
    ),
    error = function(e) {
      cat(sprintf("[ERROR] SetDatExpr: %s\n", e$message))
      stop(e)
    }
  )

  dat <- GetDatExpr(obj, wgcna_name = wgcna_name)
  cat(sprintf("[OK] datExpr: %d metacells x %d genes\n", nrow(dat), ncol(dat)))

  if (nrow(dat) < 30) {
    cat(sprintf(
      "[WARN] Only %d metacells for '%s'. Consider smaller k or merging subtypes.\n",
      nrow(dat), group_name
    ))
  }

  obj
}


# ----- 4.2 Soft power selection -----
hdwgcna_test_soft_power <- function(seurat_obj,
                                     out_dir    = OUTPUT_DIR,
                                     wgcna_name = WGCNA_NAME) {
  cat("\n--- hdWGCNA: TestSoftPowers ---\n")
  ensure_dir(out_dir)   # P2: always ensure dir before writing

  obj <- tryCatch(
    TestSoftPowers(seurat_obj, wgcna_name = wgcna_name),
    error = function(e) {
      cat(sprintf("[ERROR] TestSoftPowers: %s\n", e$message))
      stop(e)
    }
  )

  pdf_path <- file.path(out_dir, "hdwgcna_soft_power.pdf")
  pdf(pdf_path, width = FIG_WIDTH, height = FIG_HEIGHT)
  print(PlotSoftPowers(obj, wgcna_name = wgcna_name))
  dev.off()
  cat(sprintf("[OK] Soft power plot: %s\n", pdf_path))

  power_table  <- GetPowerTable(obj, wgcna_name = wgcna_name)
  first_good   <- which(power_table$SFT.R.sq >= 0.80)[1]
  recommended  <- if (!is.na(first_good)) power_table$Power[first_good] else {
    cat(sprintf(
      "[WARN] No power achieves R^2 >= 0.80. Using conservative fallback = 12.\n"
    ))
    12L
  }

  cat(sprintf("[OK] Recommended soft_power = %d (R^2 = %.3f)\n",
              recommended,
              power_table$SFT.R.sq[power_table$Power == recommended]))

  list(obj              = obj,
       recommended_power = recommended,
       power_table       = power_table)
}


# ----- 4.3 Construct network -----
hdwgcna_construct_network <- function(seurat_obj,
                                       soft_power = SOFT_POWER,
                                       net_type   = NETWORK_TYPE,
                                       tom_type   = TOMTYPE,
                                       wgcna_name = WGCNA_NAME) {
  if (is.null(soft_power)) {
    stop("[ERROR] soft_power is NULL. Run hdwgcna_test_soft_power() first.")
  }

  cat(sprintf("\n--- hdWGCNA: ConstructNetwork (power=%d, type='%s') ---\n",
              soft_power, net_type))

  obj <- tryCatch(
    ConstructNetwork(
      seurat_obj,
      soft_power    = soft_power,
      networkType   = net_type,
      TOMType       = tom_type,
      wgcna_name    = wgcna_name,
      setDatExpr    = FALSE,
      overwrite_tom = TRUE
    ),
    error = function(e) {
      cat(sprintf("[ERROR] ConstructNetwork: %s\n", e$message))
      stop(e)
    }
  )

  mods   <- GetModules(obj, wgcna_name = wgcna_name)
  n_mods <- length(setdiff(unique(mods$module), "grey"))
  cat(sprintf("[OK] Network built -- %d modules (+ grey unassigned)\n", n_mods))
  obj
}


# ==============================================================================
# 5. hdWGCNA -- Downstream analysis
# ==============================================================================

# ----- 5.1 Module eigengenes -----
hdwgcna_module_eigengenes <- function(seurat_obj,
                                       group_by   = CELLTYPE_COL,
                                       wgcna_name = WGCNA_NAME) {
  cat("\n--- hdWGCNA: ModuleEigengenes ---\n")

  mods <- tryCatch(GetModules(seurat_obj, wgcna_name = wgcna_name), error = function(e) NULL)
  gene_col <- c("gene_name", "gene")[c("gene_name", "gene") %in% colnames(mods)][1]
  scale_features <- if (!is.null(mods) && !is.na(gene_col)) {
    unique(mods[mods$module != "grey", gene_col])
  } else {
    NULL
  }

  seurat_obj <- hdwgcna_ensure_scaled_data(
    seurat_obj,
    features = scale_features,
    wgcna_name = wgcna_name
  )

  obj <- tryCatch(
    ModuleEigengenes(seurat_obj, group.by = group_by, wgcna_name = wgcna_name),
    error = function(e) {
      cat(sprintf("[ERROR] ModuleEigengenes: %s\n", e$message))
      stop(e)
    }
  )

  me <- GetMEs(obj, wgcna_name = wgcna_name)
  cat(sprintf("[OK] MEs: %d cells x %d modules\n", nrow(me), ncol(me)))
  obj
}


# ----- 5.2 Hub genes -----
hdwgcna_hub_genes <- function(seurat_obj,
                               n_hubs     = N_HUB_GENES,
                               wgcna_name = WGCNA_NAME) {
  cat(sprintf("\n--- hdWGCNA: ModuleConnectivity + Hub Genes (n=%d) ---\n",
              n_hubs))

  obj <- tryCatch(
    ModuleConnectivity(seurat_obj, wgcna_name = wgcna_name),
    error = function(e) {
      cat(sprintf("[ERROR] ModuleConnectivity: %s\n", e$message))
      stop(e)
    }
  )

  hubs <- GetHubGenes(obj, n_hubs = n_hubs, wgcna_name = wgcna_name)
  cat(sprintf("[OK] Hub genes: %d modules\n", length(unique(hubs$module))))
  list(obj = obj, hub_genes = hubs)
}


# ----- 5.3 Module-trait correlation -----
hdwgcna_trait_correlation <- function(seurat_obj,
                                       trait_col  = CONDITION_COL,
                                       out_dir    = OUTPUT_DIR,
                                       wgcna_name = WGCNA_NAME) {
  cat(sprintf("\n--- hdWGCNA: Trait Correlation [%s] ---\n", trait_col))
  ensure_dir(out_dir)   # P2

  meta <- seurat_obj@meta.data
  if (!trait_col %in% colnames(meta)) {
    cat(sprintf("[WARN] Column '%s' not found -- skipping trait correlation.\n",
                trait_col))
    return(invisible(NULL))
  }

  # One-hot encode factor levels into numeric trait matrix
  trait_levels <- sort(unique(na.omit(meta[[trait_col]])))
  trait_mat <- sapply(trait_levels, function(lv) {
    as.integer(meta[[trait_col]] == lv)
  })
  rownames(trait_mat) <- rownames(meta)
  colnames(trait_mat) <- make.names(trait_levels)

  me          <- GetMEs(seurat_obj, wgcna_name = wgcna_name)
  shared_cells <- intersect(rownames(me), rownames(trait_mat))

  if (length(shared_cells) < 20) {
    cat(sprintf("[WARN] Only %d shared cells -- skipping trait correlation.\n",
                length(shared_cells)))
    return(invisible(NULL))
  }

  me_sub    <- me[shared_cells, , drop = FALSE]
  trait_sub <- trait_mat[shared_cells, , drop = FALSE]
  n         <- length(shared_cells)

  cor_mat  <- cor(me_sub, trait_sub, use = "pairwise.complete.obs")

  # P0 fix: clamp before t-stat to avoid NaN at r = +/-1
  cor_clamped <- pmin(pmax(cor_mat, -0.9999), 0.9999)
  t_stat      <- cor_clamped * sqrt((n - 2) / (1 - cor_clamped^2))
  pval_mat    <- 2 * pt(-abs(t_stat), df = n - 2)
  padj_mat    <- apply(pval_mat, 2, function(col) p.adjust(col, method = "BH"))

  # Build tidy data frame for ggplot
  cor_df  <- as.data.frame(cor_mat) %>%
    tibble::rownames_to_column("module") %>%
    tidyr::pivot_longer(-module, names_to = "trait", values_to = "r")

  padj_df <- as.data.frame(padj_mat) %>%
    tibble::rownames_to_column("module") %>%
    tidyr::pivot_longer(-module, names_to = "trait", values_to = "padj")

  plot_df <- left_join(cor_df, padj_df, by = c("module", "trait")) %>%
    mutate(
      sig_label = dplyr::case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ ""
      )
    )

  p <- ggplot(plot_df, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig_label), size = 3.5) +
    scale_fill_gradient2(
      low  = "#2166AC", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-1, 1), name = "Pearson r"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank()
    ) +
    labs(
      title = sprintf("Module-Trait Correlation (%s)", trait_col),
      x = "Trait", y = "Module"
    )

  pdf_path <- file.path(out_dir,
                         sprintf("hdwgcna_module_trait_%s.pdf",
                                 make.names(trait_col)))
  pdf(pdf_path, width = FIG_WIDTH, height = FIG_HEIGHT)
  print(p)
  dev.off()
  cat(sprintf("[OK] Trait correlation plot: %s\n", pdf_path))

  invisible(list(cor_mat  = cor_mat,
                 pval_mat = pval_mat,
                 padj_mat = padj_mat,
                 plot     = p,
                 plot_df  = plot_df))
}


# ----- 5.4 Module UMAP -----
hdwgcna_plot_module_umap <- function(seurat_obj,
                                      out_dir    = OUTPUT_DIR,
                                      wgcna_name = WGCNA_NAME) {
  cat("\n--- hdWGCNA: Module UMAP ---\n")
  ensure_dir(out_dir)   # P2

  mods      <- GetModules(seurat_obj, wgcna_name = wgcna_name)
  mod_names <- setdiff(unique(mods$module), "grey")

  plots <- lapply(mod_names, function(mod) {
    tryCatch(
      ModuleFeaturePlot(seurat_obj, features = mod, wgcna_name = wgcna_name),
      error = function(e) {
        cat(sprintf("  [WARN] UMAP for '%s': %s\n", mod, e$message))
        NULL
      }
    )
  })
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0) {
    cat("[WARN] No module UMAP plots produced.\n")
    return(invisible(NULL))
  }

  n_col    <- min(4L, ceiling(sqrt(length(plots))))
  pdf_path <- file.path(out_dir, "hdwgcna_module_umap.pdf")
  pdf(pdf_path, width = FIG_WIDTH * 2, height = FIG_HEIGHT * 2)
  print(wrap_plots(plots, ncol = n_col))
  dev.off()

  cat(sprintf("[OK] Module UMAP (%d modules): %s\n", length(plots), pdf_path))
  invisible(plots)
}


# ----- 5.5 Hub gene dot plot (condition) -----
hdwgcna_hub_dotplot <- function(seurat_obj,
                                 hub_genes_df,
                                 group_col  = CONDITION_COL,
                                 out_dir    = OUTPUT_DIR,
                                 wgcna_name = WGCNA_NAME,
                                 n_modules  = 6L) {
  cat("\n--- hdWGCNA: Hub Gene Dot Plot ---\n")
  ensure_dir(out_dir)   # P2

  mods     <- GetModules(seurat_obj, wgcna_name = wgcna_name)
  mod_sizes <- sort(
    table(mods$module[mods$module != "grey"]),
    decreasing = TRUE
  )
  top_mods <- names(head(mod_sizes, n_modules))

  # gene_name column differs by hdWGCNA version; handle both
  gene_col  <- if ("gene_name" %in% colnames(hub_genes_df)) "gene_name" else "gene"
  plot_genes <- hub_genes_df %>%
    filter(module %in% top_mods) %>%
    arrange(module) %>%
    pull(!!sym(gene_col)) %>%
    unique() %>%
    intersect(rownames(seurat_obj))

  if (length(plot_genes) == 0) {
    cat("[WARN] No hub genes found in Seurat features -- skipping dotplot.\n")
    return(invisible(NULL))
  }

  Idents(seurat_obj) <- group_col

  p <- DotPlot(seurat_obj, features = plot_genes, dot.scale = 6) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    ) +
    scale_color_viridis_c(option = "plasma", name = "Avg Expr") +
    labs(
      title = sprintf("hdWGCNA Hub Genes by %s", group_col),
      x     = "Gene", y = group_col
    )

  pdf_path <- file.path(out_dir, "hdwgcna_hub_dotplot.pdf")
  pdf(pdf_path,
      width  = FIG_WIDTH,
      height = max(8, length(plot_genes) * 0.25))
  print(p)
  dev.off()

  cat(sprintf("[OK] Hub dotplot: %s\n", pdf_path))
  invisible(p)
}


# ----- 5.6 Export module membership CSV -----
hdwgcna_export_modules <- function(seurat_obj,
                                    out_dir    = OUTPUT_DIR,
                                    wgcna_name = WGCNA_NAME) {
  cat("\n--- hdWGCNA: Export Module Membership ---\n")
  ensure_dir(out_dir)   # P2

  mods <- GetModules(seurat_obj, wgcna_name = wgcna_name)

  size_tbl <- mods %>%
    filter(module != "grey") %>%
    dplyr::count(module, name = "n_genes") %>%
    arrange(desc(n_genes))

  cat(sprintf("[INFO] %d non-grey modules:\n", nrow(size_tbl)))
  for (i in seq_len(nrow(size_tbl))) {
    cat(sprintf("  %-30s %d genes\n",
                size_tbl$module[i], size_tbl$n_genes[i]))
  }

  csv_path <- file.path(out_dir, "hdwgcna_module_membership.csv")
  write.csv(mods, csv_path, row.names = FALSE)
  cat(sprintf("[OK] Module membership saved: %s\n", csv_path))

  invisible(list(modules = mods, sizes = size_tbl))
}


# ==============================================================================
# 6. hdWGCNA -- One-shot runner (P1 refactored)
# ==============================================================================
# Architecture change from v1.0:
#   Global step (setup + metacells) is done ONCE outside the cell-type loop.
#   Per-celltype step (datexpr + soft power + network + downstream) is looped.
#   This avoids repeating expensive setup O(n_celltypes) times.
#
# Usage pattern:
#   obj <- hdwgcna_global_setup(seurat_obj)
#   results <- hdwgcna_run_celltypes(obj, celltypes = c("Naive_B", "Memory_B"))

hdwgcna_global_setup <- function(seurat_obj,
                                  out_dir         = OUTPUT_DIR,
                                  group_by        = c(CELLTYPE_COL, SAMPLE_COL),
                                  gene_select_mode = GENE_SELECT_MODE,
                                  gene_fraction    = GENE_FRACTION,
                                  gene_n_top       = GENE_N_TOP,
                                  metacell_k      = N_METACELL_K,
                                  metacell_target = METACELL_TARGET,
                                  gene_exclusion_config = list(),
                                  wgcna_name      = WGCNA_NAME) {
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat("hdWGCNA: Global setup (run once)\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")

  out_dir <- ensure_dir(out_dir)
  exclusion_res <- pa_apply_gene_exclusion_to_seurat(
    seurat_obj = seurat_obj,
    output_dir = out_dir,
    gene_exclusion_config = gene_exclusion_config,
    prefix = "gene_exclusion"
  )
  seurat_obj <- exclusion_res$seurat_obj
  cat(sprintf(
    "[INFO] Gene exclusion kept %d/%d features (removed %d)\n",
    exclusion_res$manifest$n_kept_features,
    exclusion_res$manifest$n_total_features,
    exclusion_res$manifest$n_excluded_features
  ))

  obj <- hdwgcna_setup(
    seurat_obj,
    wgcna_name       = wgcna_name,
    gene_select_mode = gene_select_mode,
    gene_fraction    = gene_fraction,
    gene_n_top       = gene_n_top
  )

  obj <- hdwgcna_metacells(
    obj,
    group_by        = group_by,
    metacell_k      = metacell_k,
    metacell_target = metacell_target,
    wgcna_name      = wgcna_name
  )

  cat("\n[OK] Global setup complete. Proceed with hdwgcna_run_celltypes().\n")
  invisible(obj)
}


hdwgcna_run_celltypes <- function(seurat_obj,
                                   celltypes  = NULL,
                                   soft_power = SOFT_POWER,
                                   out_dir    = OUTPUT_DIR,
                                   group_by   = CELLTYPE_COL,
                                   wgcna_name = WGCNA_NAME,
                                   resume_celltypes = FALSE,
                                   resume_skip_statuses = c("ok", "no_modules"),
                                   celltype_timeout_sec = NULL) {
  if (is.null(celltypes)) {
    celltypes <- sort(unique(as.character(seurat_obj@meta.data[[CELLTYPE_COL]])))
  }

  resume_celltypes <- isTRUE(resume_celltypes)
  resume_skip_statuses <- unique(as.character(resume_skip_statuses))
  celltype_timeout_sec <- suppressWarnings(as.numeric(celltype_timeout_sec)[1])
  if (length(celltype_timeout_sec) == 0L || is.na(celltype_timeout_sec)) celltype_timeout_sec <- 0

  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("hdWGCNA: Per-celltype loop (%d cell types)\n", length(celltypes)))
  cat(sprintf("[INFO] resume_celltypes = %s; celltype_timeout_sec = %s\n",
              ifelse(resume_celltypes, "TRUE", "FALSE"),
              ifelse(celltype_timeout_sec > 0, as.character(celltype_timeout_sec), "disabled")))
  cat(paste0(rep("=", 80), collapse = ""), "\n")

  ensure_dir(out_dir)
  results <- list()

  for (ct in celltypes) {
    ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
    ct_dir  <- file.path(out_dir, sprintf("hdwgcna_%s", ct_safe))
    ensure_dir(ct_dir)

    cat(sprintf("\n%s\n", paste(rep("-", 60), collapse = "")))
    cat(sprintf("[INFO] Cell type: %s\n", ct))

    if (resume_celltypes) {
      existing_record <- hdwgcna_existing_celltype_record(ct, ct_dir, resume_skip_statuses = resume_skip_statuses)
      if (!is.null(existing_record)) {
        results[[ct]] <- existing_record
        cat(sprintf("[SKIP] %s existing status=%s -- %s\n", ct, existing_record$status, ct_dir))
        hdwgcna_write_partial_results(results, out_dir)
        next
      }
    }

    started_at <- hdwgcna_timestamp()
    hdwgcna_write_celltype_status(
      list(
        celltype = ct,
        status = "started",
        output_dir = ct_dir,
        started_at = started_at,
        timeout_sec = ifelse(celltype_timeout_sec > 0, celltype_timeout_sec, NA_real_)
      ),
      ct_dir
    )

    # Per-celltype steps
    obj <- tryCatch(hdwgcna_with_timeout({
      o <- hdwgcna_set_datexpr(seurat_obj,
                                group_name = ct,
                                group_by   = group_by,
                                wgcna_name = wgcna_name)

      sp_res <- hdwgcna_test_soft_power(o, out_dir = ct_dir,
                                         wgcna_name = wgcna_name)
      o       <- sp_res$obj
      sp_used <- if (!is.null(soft_power)) soft_power else sp_res$recommended_power
      cat(sprintf("[INFO] soft_power used = %d\n", sp_used))

      o <- hdwgcna_construct_network(o, soft_power = sp_used,
                                      wgcna_name = wgcna_name)

      mods <- tryCatch(GetModules(o, wgcna_name = wgcna_name), error = function(e) NULL)
      module_ids <- if (!is.null(mods) && "module" %in% colnames(mods)) {
        setdiff(unique(mods$module), "grey")
      } else {
        character()
      }

      if (length(module_ids) == 0L) {
        cat(sprintf("[WARN] %s: no non-grey modules detected -- skipping downstream steps.\n", ct))
        results[[ct]] <- list(
          celltype    = ct,
          status      = "no_modules",
          hub_genes   = NULL,
          module_ids  = character(),
          output_dir  = ct_dir,
          started_at  = started_at
        )
        hdwgcna_write_celltype_status(results[[ct]][setdiff(names(results[[ct]]), "hub_genes")], ct_dir)
        hdwgcna_write_partial_results(results, out_dir)
        cat(sprintf("[OK] %s complete -- %s (no non-grey modules)\n", ct, ct_dir))
      } else {
        o <- hdwgcna_module_eigengenes(o, wgcna_name = wgcna_name)

        hub_res <- hdwgcna_hub_genes(o, wgcna_name = wgcna_name)
        o        <- hub_res$obj

        hdwgcna_trait_correlation(o, out_dir = ct_dir, wgcna_name = wgcna_name)
        hdwgcna_plot_module_umap(o,  out_dir = ct_dir, wgcna_name = wgcna_name)
        hdwgcna_export_modules(o,    out_dir = ct_dir, wgcna_name = wgcna_name)
        hdwgcna_hub_dotplot(o,
                             hub_genes_df = hub_res$hub_genes,
                             out_dir      = ct_dir,
                             wgcna_name   = wgcna_name)

        if (exists("hdwgcna_plot_summary", mode = "function", inherits = TRUE)) {
          tryCatch(
            hdwgcna_plot_summary(
              seurat_obj = o,
              hub_genes_df = hub_res$hub_genes,
              group_cols = c(CONDITION_COL, TISSUE_COL),
              out_dir = ct_dir,
              wgcna_name = wgcna_name,
              make_umap = FALSE,
              make_violin = TRUE
            ),
            error = function(e) cat(sprintf("[WARN] Enhanced hdWGCNA visualizations failed for %s: %s\n", ct, e$message))
          )
        }

        results[[ct]] <- list(
          celltype   = ct,
          status     = "ok",
          hub_genes  = hub_res$hub_genes,
          module_ids = module_ids,
          output_dir = ct_dir,
          started_at = started_at
        )
        hdwgcna_write_celltype_status(results[[ct]][setdiff(names(results[[ct]]), "hub_genes")], ct_dir)
        hdwgcna_write_partial_results(results, out_dir)

        cat(sprintf("[OK] %s complete -- %s\n", ct, ct_dir))
      }

      o
    }, timeout_sec = celltype_timeout_sec),
    error = function(e) {
      status <- if (hdwgcna_is_timeout_error(e)) "timeout" else "error"
      results[[ct]] <<- list(
        celltype = ct,
        status = status,
        error = conditionMessage(e),
        output_dir = ct_dir,
        started_at = started_at,
        timeout_sec = ifelse(celltype_timeout_sec > 0, celltype_timeout_sec, NA_real_)
      )
      hdwgcna_write_celltype_status(results[[ct]], ct_dir)
      hdwgcna_write_partial_results(results, out_dir)
      if (identical(status, "timeout")) {
        cat(sprintf("[TIMEOUT] %s exceeded %.0f sec: %s -- skipping.\n", ct, celltype_timeout_sec, conditionMessage(e)))
      } else {
        cat(sprintf("[ERROR] %s failed: %s -- skipping.\n", ct, conditionMessage(e)))
      }
      NULL
    })
  }

  status_vec <- if (length(results) == 0L) {
    character()
  } else {
    vapply(results, function(x) {
      if (is.null(x$status)) "ok" else as.character(x$status)[1]
    }, character(1))
  }

  n_ok <- sum(status_vec == "ok")
  n_no_modules <- sum(status_vec == "no_modules")
  n_timeout <- sum(status_vec == "timeout")
  n_error <- sum(status_vec == "error")
  n_skipped <- if (length(results) == 0L) 0L else sum(vapply(results, function(x) isTRUE(x$skipped_existing), logical(1)))

  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("[OK] hdWGCNA loop done. %d/%d cell types succeeded; %d had no non-grey modules; %d timed out; %d errored; %d skipped existing.\n",
              n_ok, length(celltypes), n_no_modules, n_timeout, n_error, n_skipped))

  invisible(results)
}


# ==============================================================================
# 7. CoVarNet -- Cell-type-stratified co-variation network
# ==============================================================================

# ----- 7.1 Compute gene-gene correlation per cell type -----
# P0 fix: cor_mat clamped before t-stat to avoid NaN at r = +/-1
# P2 fix: sparse-aware row variance via Matrix::rowMeans on squared deviations
covarnet_compute_cor <- function(seurat_obj,
                                  celltype_name,
                                  n_genes     = COVAR_N_GENES,
                                  cor_method  = COVAR_COR_METHOD,
                                  min_cells   = COVAR_MIN_CELLS) {
  cat(sprintf("\n--- CoVarNet: Correlation [%s] ---\n", celltype_name))

  cells <- colnames(seurat_obj)[
    seurat_obj@meta.data[[CELLTYPE_COL]] == celltype_name
  ]

  if (length(cells) < min_cells) {
    cat(sprintf("[WARN] %s: %d cells < min_cells=%d -- skipping.\n",
                celltype_name, length(cells), min_cells))
    return(NULL)
  }

  sub_obj <- subset(seurat_obj, cells = cells)
  expr    <- GetAssayData(sub_obj, layer = "data")   # log-norm, likely sparse

  # P2: sparse-aware row variance
  # var(x) = E[x^2] - E[x]^2, computed without densifying full matrix
  row_means  <- Matrix::rowMeans(expr)
  row_means2 <- Matrix::rowMeans(expr^2)
  row_vars   <- row_means2 - row_means^2

  if (length(row_vars) < n_genes) {
    cat(sprintf("[WARN] %s: only %d genes available, requested %d\n",
                celltype_name, length(row_vars), n_genes))
  }

  top_genes <- names(sort(row_vars, decreasing = TRUE))[
    seq_len(min(n_genes, length(row_vars)))
  ]

  # Dense conversion only on the small HVG subset
  expr_sub <- t(as.matrix(expr[top_genes, , drop = FALSE]))  # cells x genes

  cat(sprintf("[INFO] %s: %d cells x %d genes, %s correlation\n",
              celltype_name, nrow(expr_sub), ncol(expr_sub), cor_method))

  cor_mat <- cor(expr_sub, method = cor_method, use = "pairwise.complete.obs")

  # P0 fix: clamp to (-0.9999, 0.9999) before t-stat
  n           <- nrow(expr_sub)
  cor_clamped <- pmin(pmax(cor_mat, -0.9999), 0.9999)
  t_stat      <- cor_clamped * sqrt((n - 2) / (1 - cor_clamped^2))
  pval_mat    <- 2 * pt(-abs(t_stat), df = n - 2)
  diag(pval_mat) <- 1   # self-pair -> p = 1

  cat(sprintf("[OK] %s: correlation matrix %d x %d\n",
              celltype_name, nrow(cor_mat), ncol(cor_mat)))

  list(
    cor_mat  = cor_mat,
    pval_mat = pval_mat,
    n_cells  = length(cells),
    genes    = top_genes,
    celltype = celltype_name
  )
}


# ----- 7.2 Build edge list with significance filter -----
covarnet_build_edges <- function(cor_result,
                                  cor_thr  = COVAR_COR_THR,
                                  pval_thr = COVAR_PVAL_THR) {
  if (is.null(cor_result)) return(NULL)

  ct    <- cor_result$celltype
  cr    <- cor_result$cor_mat
  pr    <- cor_result$pval_mat
  genes <- rownames(cr)

  idx <- which(upper.tri(cr), arr.ind = TRUE)
  edge_df <- data.frame(
    gene_a = genes[idx[, 1]],
    gene_b = genes[idx[, 2]],
    r      = cr[idx],
    pval   = pr[idx],
    stringsAsFactors = FALSE
  )

  edge_df$padj <- p.adjust(edge_df$pval, method = "BH")

  edge_df <- edge_df %>%
    filter(abs(r) >= cor_thr, padj < pval_thr) %>%
    mutate(
      direction = ifelse(r > 0, "positive", "negative"),
      weight    = abs(r),
      celltype  = ct
    ) %>%
    arrange(desc(weight))

  if (nrow(edge_df) == 0) {
    cat(sprintf("[WARN] %s: 0 edges pass thresholds (|r|>=%.2f, FDR<%.2f)\n",
                ct, cor_thr, pval_thr))
    return(NULL)
  }

  cat(sprintf("[OK] %s: %d edges (|r|>=%.2f, FDR<%.2f)\n",
              ct, nrow(edge_df), cor_thr, pval_thr))
  edge_df
}


# ----- 7.3 Build igraph from edge list -----
covarnet_build_graph <- function(edge_df) {
  if (is.null(edge_df) || nrow(edge_df) == 0) return(NULL)

  g <- igraph::graph_from_data_frame(
    edge_df[, c("gene_a", "gene_b", "r", "weight", "direction", "padj")],
    directed = FALSE
  )

  V(g)$degree      <- igraph::degree(g)
  V(g)$betweenness <- igraph::betweenness(g, normalized = TRUE)
  V(g)$hub_score   <- igraph::hub_score(g)$vector

  comm             <- igraph::cluster_louvain(g, weights = E(g)$weight)
  V(g)$community   <- comm$membership

  cat(sprintf("[OK] igraph: %d nodes, %d edges, %d communities\n",
              vcount(g), ecount(g), max(comm$membership)))
  g
}


# ----- 7.4 Extract hub genes from CoVarNet graph -----
covarnet_hub_genes <- function(graph, n_top = N_HUB_GENES) {
  if (is.null(graph)) return(NULL)

  node_df <- data.frame(
    gene        = V(graph)$name,
    degree      = V(graph)$degree,
    betweenness = V(graph)$betweenness,
    hub_score   = V(graph)$hub_score,
    community   = V(graph)$community,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(hub_score))

  top_hubs <- head(node_df, n_top)
  cat(sprintf("[INFO] Top hubs: %s\n",
              paste(head(top_hubs$gene, 5), collapse = ", ")))

  list(all_nodes = node_df, hub_genes = top_hubs)
}


# ----- 7.5 Full CoVarNet pipeline for all / selected cell types -----
covarnet_run_all <- function(seurat_obj,
                              celltypes = COVAR_CELLTYPES,
                              out_dir   = OUTPUT_DIR,
                              gene_exclusion_config = list(),
                              ...) {
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat("CoVarNet: Full pipeline\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  ensure_dir(out_dir)

  exclusion_res <- pa_apply_gene_exclusion_to_seurat(
    seurat_obj = seurat_obj,
    output_dir = out_dir,
    gene_exclusion_config = gene_exclusion_config,
    prefix = "gene_exclusion"
  )
  seurat_obj <- exclusion_res$seurat_obj
  cat(sprintf(
    "[INFO] Gene exclusion kept %d/%d features (removed %d)\n",
    exclusion_res$manifest$n_kept_features,
    exclusion_res$manifest$n_total_features,
    exclusion_res$manifest$n_excluded_features
  ))

  if (is.null(celltypes)) {
    celltypes <- sort(unique(seurat_obj@meta.data[[CELLTYPE_COL]]))
  }
  cat(sprintf("[INFO] %d cell types to process\n", length(celltypes)))

  results <- list()

  for (ct in celltypes) {
    ct_safe  <- gsub("[^A-Za-z0-9_]", "_", ct)

    cor_res  <- covarnet_compute_cor(seurat_obj, ct, ...)
    edge_df  <- covarnet_build_edges(cor_res)
    graph    <- covarnet_build_graph(edge_df)
    hub_res  <- covarnet_hub_genes(graph)

    if (!is.null(edge_df)) {
      write.csv(
        edge_df,
        file.path(out_dir, sprintf("covarnet_%s_edges.csv", ct_safe)),
        row.names = FALSE
      )
    }
    if (!is.null(hub_res)) {
      write.csv(
        hub_res$all_nodes,
        file.path(out_dir, sprintf("covarnet_%s_nodes.csv", ct_safe)),
        row.names = FALSE
      )
    }

    results[[ct]] <- list(
      cor_result = cor_res,
      edges      = edge_df,
      graph      = graph,
      hubs       = hub_res
    )
  }

  n_ok <- sum(sapply(results, function(x) !is.null(x$graph)))
  attr(results, "gene_exclusion_manifest") <- exclusion_res$manifest
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("[OK] CoVarNet done. Graphs produced for %d/%d cell types.\n",
              n_ok, length(celltypes)))
  invisible(results)
}


# ==============================================================================
# 8. Visualization helpers
# ==============================================================================

# ----- 8.1 Network plot (ggraph) for one cell type -----
plot_covarnet_graph <- function(graph,
                                 celltype_name = "Unknown",
                                 n_label_top   = 20L,
                                 out_dir       = OUTPUT_DIR) {
  if (is.null(graph) || vcount(graph) == 0) return(invisible(NULL))
  ensure_dir(out_dir)   # P2

  # Subset to top-degree nodes for readability
  deg       <- igraph::degree(graph)
  top_nodes <- names(sort(deg, decreasing = TRUE))[seq_len(
    min(200L, vcount(graph))
  )]
  g_sub <- igraph::induced_subgraph(graph, top_nodes)
  tg    <- tidygraph::as_tbl_graph(g_sub)

  p <- ggraph(tg, layout = "fr") +
    geom_edge_link(
      aes(alpha = weight, color = direction),
      show.legend = TRUE
    ) +
    geom_node_point(
      aes(size = degree, color = as.factor(community)),
      alpha = 0.8
    ) +
    geom_node_label(
      aes(label = ifelse(rank(-degree) <= n_label_top, name, "")),
      repel = TRUE, size = 2.5, max.overlaps = 20
    ) +
    scale_edge_color_manual(
      values = c("positive" = "#E41A1C", "negative" = "#377EB8"),
      name   = "Direction"
    ) +
    scale_edge_alpha(range = c(0.1, 0.6), guide = "none") +
    scale_size_continuous(name = "Degree", range = c(1, 6)) +
    scale_color_discrete(name = "Community") +
    theme_graph(base_family = "sans") +
    labs(
      title    = sprintf("CoVarNet: %s", celltype_name),
      subtitle = sprintf("%d nodes | %d edges", vcount(g_sub), ecount(g_sub))
    )

  ct_safe  <- gsub("[^A-Za-z0-9_]", "_", celltype_name)
  pdf_path <- file.path(out_dir,
                         sprintf("covarnet_%s_network.pdf", ct_safe))
  pdf(pdf_path, width = FIG_WIDTH, height = FIG_HEIGHT)
  print(p)
  dev.off()
  cat(sprintf("[OK] Network plot: %s\n", pdf_path))
  invisible(p)
}


# ----- 8.2 Cross-cell-type hub gene overlap heatmap -----
plot_hub_overlap <- function(covarnet_results,
                              n_top   = N_HUB_GENES,
                              out_dir = OUTPUT_DIR) {
  cat("\n--- Visualization: Hub Gene Overlap ---\n")
  ensure_dir(out_dir)   # P2

  hub_lists <- lapply(covarnet_results, function(res) {
    if (is.null(res$hubs)) return(character(0))
    head(res$hubs$hub_genes$gene, n_top)
  })
  hub_lists <- Filter(function(x) length(x) > 0, hub_lists)

  if (length(hub_lists) < 2) {
    cat("[WARN] Need >= 2 cell types with hub genes for overlap plot.\n")
    return(invisible(NULL))
  }

  all_genes <- unique(unlist(hub_lists))
  mat <- sapply(hub_lists, function(h) as.integer(all_genes %in% h))
  rownames(mat) <- all_genes

  shared_df <- data.frame(
    gene    = all_genes,
    n_types = rowSums(mat),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(n_types)) %>%
    filter(n_types > 1)

  cat(sprintf("[INFO] %d genes are hub in > 1 cell type\n", nrow(shared_df)))
  if (nrow(shared_df) > 0) {
    cat(sprintf("[INFO] Top shared hubs: %s\n",
                paste(head(shared_df$gene, 10), collapse = ", ")))
  }

  plot_mat <- mat[rownames(mat) %in% head(shared_df$gene, 50), , drop = FALSE]

  if (nrow(plot_mat) == 0) {
    cat("[WARN] No shared hub genes to plot.\n")
    return(invisible(NULL))
  }

  pdf_path <- file.path(out_dir, "covarnet_hub_overlap_heatmap.pdf")
  pdf(pdf_path,
      width  = FIG_WIDTH,
      height = max(6, nrow(plot_mat) * 0.2 + 2))

  pheatmap::pheatmap(
    plot_mat,
    color         = c("white", "#D6604D"),
    breaks        = c(-0.5, 0.5, 1.5),
    legend_breaks = c(0, 1),
    legend_labels = c("absent", "hub"),
    cluster_cols  = TRUE,
    cluster_rows  = TRUE,
    fontsize_row  = 7,
    fontsize_col  = 9,
    main = sprintf("Hub Gene Overlap (top %d per cell type)", n_top)
  )

  dev.off()
  cat(sprintf("[OK] Hub overlap heatmap: %s\n", pdf_path))
  invisible(list(shared = shared_df, matrix = mat))
}


# ==============================================================================
# 9. Usage example  (P1: all variables defined, directly runnable)
# ==============================================================================

# load_coexpr_libs()
#
# # --- Load Seurat object (uses project-standard GetSeurat) ---
# if (!file.exists(H5AD_PATH)) stop(sprintf("File not found: %s", H5AD_PATH))
# seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
#
# ensure_dir(OUTPUT_DIR)
#
# # ===========================================================================
# # Workflow A: hdWGCNA
# # ===========================================================================
#
# # Step 1: Global setup -- run ONCE regardless of how many cell types
# seurat_obj <- hdwgcna_global_setup(seurat_obj)
#
# # Step 2: Per-cell-type loop
# #   Pass NULL to process all cell types in CELLTYPE_COL,
# #   or specify a character vector.
# wgcna_results <- hdwgcna_run_celltypes(
#   seurat_obj = seurat_obj,
#   celltypes  = c("Naive_B", "Memory_B", "GC_B"),
#   out_dir    = file.path(OUTPUT_DIR, "hdwgcna")
# )
#
# # ===========================================================================
# # Workflow B: CoVarNet
# # ===========================================================================
#
# covar_dir <- file.path(OUTPUT_DIR, "covarnet")
#
# covar_results <- covarnet_run_all(
#   seurat_obj = seurat_obj,
#   celltypes  = c("Naive_B", "Memory_B", "GC_B", "Plasma_IgA"),
#   out_dir    = covar_dir
# )
#
# # Network plot per cell type
# for (ct in names(covar_results)) {
#   plot_covarnet_graph(
#     graph         = covar_results[[ct]]$graph,
#     celltype_name = ct,
#     out_dir       = covar_dir
#   )
# }
#
# # Cross-cell-type hub overlap
# plot_hub_overlap(covar_results, out_dir = covar_dir)
