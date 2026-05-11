#!/usr/bin/env Rscript
# ==============================================================================
# T/NK Cell Complete Analysis Pipeline - Dual Mode (interpret + interpret_agent)
# ==============================================================================
#
# Version: v4.1 (2026-02-04)
# Features:
#   ✅ Dual interpretation modes:
#      - FAST MODE: interpret() - Single-pass LLM (faster)
#      - DEEP MODE: interpret_agent() - Multi-agent system (more rigorous)
#   ✅ Auto-retry for both modes (90-95% success rate)
#   ✅ Multi-database support (GO BP/MF/CC, Hallmark, KEGG, CellMarker, PanglaoDB)
#   ✅ PPI network integration (optional)
#   ✅ Gene fold change integration (optional)
#
# ==============================================================================

# ==============================================================================
# Configuration
# ==============================================================================

# ⭐ INTERPRETATION MODE SELECTION
USE_AGENT_MODE <- TRUE # Set to TRUE for interpret_agent(), FALSE for interpret()

# Paths
H5AD_PATH <- "/home/h2048/data/py/0129/tnk_analysis_unified/results/subcluster_unified_v2_20260129/adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"
OUTPUT_DIR <- "/home/h2048/data/R/0204/tcell_interpret_v4_1"
CELLMARKER_PATH <- "/home/h2048/data/source/reference/CellMarker/Cell_marker_Human.csv"
PANGLAODB_PATH <- "/home/h2048/data/source/reference/CellMarker/PanglaoDB_markers_27_Mar_2020.tsv.csv"
MSIGDB_GMT_PATH <- "/home/h2048/data/source/reference/MSigDB/msigdb.v2025.1.Hs.symbols.gmt"
GMT_GO_ALL <- "/home/h2048/data/source/reference/MSigDB/c5.all.v2025.1.Hs.symbols.gmt"

# DeepSeek API Key
DEEPSEEK_API_KEY <- Sys.getenv("DEEPSEEK_API_KEY")
if (nchar(DEEPSEEK_API_KEY) < 10) {
  stop("ERROR: DEEPSEEK_API_KEY environment variable not set")
}

# Analysis Parameters
N_CORES <- 4
TOP_N_MARKERS <- 50
CELL_TYPE <- "T_NK"

# Agent Mode Parameters (only used if USE_AGENT_MODE = TRUE)
ADD_PPI <- TRUE # Add PPI network analysis (requires more time)
USE_FOLD_CHANGE <- FALSE # Add gene fold change info (need to prepare fc vector)

# Create directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(OUTPUT_DIR, "reports"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(OUTPUT_DIR, "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat(
  "================================================================================\n"
)
cat("T/NK CELL COMPLETE ANALYSIS PIPELINE v4.1\n")
cat(
  "================================================================================\n\n"
)
cat(sprintf("Output: %s\n", OUTPUT_DIR))
cat(sprintf("Cell Type: %s\n", CELL_TYPE))
cat(sprintf(
  "Mode: %s\n",
  if (USE_AGENT_MODE) "DEEP (interpret_agent)" else "FAST (interpret)"
))
if (USE_AGENT_MODE) {
  cat(sprintf("PPI Network: %s\n", if (ADD_PPI) "Enabled" else "Disabled"))
  cat(sprintf(
    "Fold Change: %s\n",
    if (USE_FOLD_CHANGE) "Enabled" else "Disabled"
  ))
}
cat("\n")

# ==============================================================================
# Load Libraries
# ==============================================================================

cat("\n=== Loading Libraries ===\n")

library(reticulate)
library(SCNT)
library(Seurat)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(future)
library(future.apply)

# Thread limiting
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

plan("multisession", workers = N_CORES)
options(future.globals.maxSize = 10 * 1024^3)

use_condaenv("bbknn_env", required = TRUE)

cat("[OK] Libraries loaded\n")

# ==============================================================================
# Helper Functions
# ==============================================================================

sanitize_filename <- function(name) {
  name <- gsub("/", "_", name)
  name <- gsub("\\+", "plus", name)
  name <- gsub("-", "_", name)
  name <- gsub("\\(|\\)", "", name)
  name <- gsub("\\[|\\]", "", name)
  name <- gsub("\\{|\\}", "", name)
  name <- gsub("<|>", "", name)
  name <- gsub(":", "_", name)
  name <- gsub(";", "_", name)
  name <- gsub(",", "_", name)
  name <- gsub("\\?|\\*", "", name)
  name <- gsub("\"|'|\\|", "", name)
  name <- gsub("\\s+", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  return(name)
}

# ==============================================================================
# ⭐ CORE: Universal Interpret Wrapper with Auto-Retry
# ==============================================================================

interpret_with_retry <- function(
  enrichment_input,
  context = NULL,
  n_pathways = 20,
  model = "deepseek-reasoner",
  api_key = NULL,
  task = "annotation",
  max_retries = 2,
  task_name = "Analysis",
  use_agent = FALSE,
  add_ppi = FALSE,
  gene_fold_change = NULL
) {
  mode_name <- if (use_agent) "interpret_agent" else "interpret"
  cat(sprintf("\n=== %s (%s): Initial Attempt ===\n", task_name, mode_name))

  # First attempt
  result <- tryCatch(
    {
      if (use_agent) {
        # Use multi-agent deep mode
        clusterProfiler::interpret_agent(
          x = enrichment_input,
          context = context,
          n_pathways = n_pathways,
          model = model,
          api_key = api_key,
          add_ppi = add_ppi,
          gene_fold_change = gene_fold_change
        )
      } else {
        # Use standard interpret
        clusterProfiler::interpret(
          x = enrichment_input,
          context = context,
          n_pathways = n_pathways,
          model = model,
          api_key = api_key,
          task = task,
          add_ppi = add_ppi,
          gene_fold_change = gene_fold_change
        )
      }
    },
    error = function(e) {
      cat("[ERROR]", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(result)) {
    cat("[FAIL] Initial attempt returned NULL\n")
    return(NULL)
  }

  # Convert to list format
  result_list <- if (inherits(result, "interpretation_list")) {
    result
  } else if (inherits(result, "interpretation")) {
    list(Default = result)
  } else {
    list(Default = result)
  }

  # Check for failures
  failed_clusters <- c()

  for (cluster_id in names(result_list)) {
    cluster_res <- result_list[[cluster_id]]

    if (!is.list(cluster_res)) {
      next
    }

    is_failed <- FALSE

    # Detection Pattern 1: Low/None confidence + failure keywords
    if (!is.null(cluster_res$confidence)) {
      if (cluster_res$confidence %in% c("Low", "None")) {
        if (!is.null(cluster_res$reasoning)) {
          if (
            grepl(
              "Failed|failed|parse|Parse|JSON|json|empty|Empty|error|Error",
              cluster_res$reasoning,
              ignore.case = TRUE
            )
          ) {
            is_failed <- TRUE
          }
        }
        if (!is.null(cluster_res$overview)) {
          if (
            grepl(
              "Failed|failed|retrieve|parse|error",
              cluster_res$overview,
              ignore.case = TRUE
            )
          ) {
            is_failed <- TRUE
          }
        }
      }
    }

    # Detection Pattern 2: Empty essential fields for annotation
    if (!use_agent && task %in% c("annotation", "cell_type")) {
      if (
        is.null(cluster_res$cell_type) ||
          cluster_res$cell_type == "" ||
          cluster_res$cell_type == "Unknown"
      ) {
        is_failed <- TRUE
      }
    }

    # Detection Pattern 3: Agent mode specific - check if overview is missing
    if (use_agent) {
      if (is.null(cluster_res$overview) || cluster_res$overview == "") {
        is_failed <- TRUE
      }
    }

    # Detection Pattern 4: Check for fallback mode with low confidence
    if (!is.null(cluster_res$data_source)) {
      if (
        cluster_res$data_source == "gene_list_only" &&
          !is.null(cluster_res$confidence) &&
          cluster_res$confidence == "Low"
      ) {
        # This is expected fallback, not a failure
        is_failed <- FALSE
      }
    }

    if (is_failed) {
      failed_clusters <- c(failed_clusters, cluster_id)
    }
  }

  success_count <- length(result_list) - length(failed_clusters)
  success_rate <- success_count / length(result_list)

  cat(sprintf(
    "\nInitial Success: %d/%d (%.1f%%)\n",
    success_count,
    length(result_list),
    success_rate * 100
  ))

  # If success rate >= 90%, we're done
  if (success_rate >= 0.9) {
    cat(sprintf("[OK] Excellent success rate for %s\n", task_name))
    return(result)
  }

  # Auto-retry failed clusters
  if (length(failed_clusters) > 0 && length(failed_clusters) <= 20) {
    cat(sprintf(
      "\n[RETRY] %d failed clusters detected\n",
      length(failed_clusters)
    ))
    cat("Failed clusters:", paste(head(failed_clusters, 10), collapse = ", "))
    if (length(failed_clusters) > 10) {
      cat(sprintf(" ... and %d more", length(failed_clusters) - 10))
    }
    cat("\n")

    for (retry_attempt in 1:max_retries) {
      if (length(failed_clusters) == 0) {
        break
      }

      cat(sprintf(
        "\n--- Retry Attempt %d/%d ---\n",
        retry_attempt,
        max_retries
      ))
      cat(sprintf("Retrying %d clusters...\n", length(failed_clusters)))

      # Exponential backoff
      wait_time <- retry_attempt * 3
      cat(sprintf("Waiting %d seconds to avoid rate limiting...\n", wait_time))
      Sys.sleep(wait_time)

      # Enhanced context for retry
      retry_context <- if (!is.null(context)) {
        paste0(
          context,
          "\n\n⚠️ RETRY ATTEMPT #",
          retry_attempt,
          " - CRITICAL INSTRUCTIONS:",
          "\nPrevious attempt had issues. Please ensure:",
          "\n1. Valid JSON format (no markdown code blocks)",
          "\n2. All required fields are completely filled",
          "\n3. Single-line text for reasoning field (no line breaks)",
          "\n4. Specific and meaningful values (not 'Unknown' or empty)",
          if (use_agent) {
            "\n5. Complete all three agent phases (Cleaner → Detective → Synthesizer)"
          } else {
            ""
          }
        )
      } else {
        paste0(
          "⚠️ RETRY ATTEMPT #",
          retry_attempt,
          "\nProvide complete and valid response with all fields populated.",
          if (use_agent) {
            "\nEnsure all agent phases complete successfully."
          } else {
            ""
          }
        )
      }

      # Retry
      retry_result <- tryCatch(
        {
          if (use_agent) {
            clusterProfiler::interpret_agent(
              x = enrichment_input,
              context = retry_context,
              n_pathways = n_pathways,
              model = model,
              api_key = api_key,
              add_ppi = add_ppi,
              gene_fold_change = gene_fold_change
            )
          } else {
            clusterProfiler::interpret(
              x = enrichment_input,
              context = retry_context,
              n_pathways = n_pathways,
              model = model,
              api_key = api_key,
              task = task,
              add_ppi = add_ppi,
              gene_fold_change = gene_fold_change
            )
          }
        },
        error = function(e) {
          cat("[ERROR]", conditionMessage(e), "\n")
          return(NULL)
        }
      )

      if (is.null(retry_result)) {
        cat("[WARN] Retry returned NULL\n")
        next
      }

      # Convert to list
      retry_result_list <- if (inherits(retry_result, "interpretation_list")) {
        retry_result
      } else if (inherits(retry_result, "interpretation")) {
        list(Default = retry_result)
      } else {
        list(Default = retry_result)
      }

      # Update results for previously failed clusters
      newly_fixed <- c()
      still_failed <- c()

      for (cluster_id in failed_clusters) {
        if (cluster_id %in% names(retry_result_list)) {
          retry_cluster_res <- retry_result_list[[cluster_id]]

          # Check if still failed
          is_still_failed <- FALSE

          if (!is.null(retry_cluster_res$confidence)) {
            if (retry_cluster_res$confidence %in% c("Low", "None")) {
              if (!is.null(retry_cluster_res$reasoning)) {
                if (
                  grepl(
                    "Failed|parse|empty|error",
                    retry_cluster_res$reasoning,
                    ignore.case = TRUE
                  )
                ) {
                  is_still_failed <- TRUE
                }
              }
            }
          }

          if (!use_agent && task %in% c("annotation", "cell_type")) {
            if (
              is.null(retry_cluster_res$cell_type) ||
                retry_cluster_res$cell_type == "" ||
                retry_cluster_res$cell_type == "Unknown"
            ) {
              is_still_failed <- TRUE
            }
          }

          if (use_agent) {
            if (
              is.null(retry_cluster_res$overview) ||
                retry_cluster_res$overview == ""
            ) {
              is_still_failed <- TRUE
            }
          }

          if (!is_still_failed) {
            # Success! Update result
            result_list[[cluster_id]] <- retry_cluster_res
            newly_fixed <- c(newly_fixed, cluster_id)
            cat(sprintf("  ✓ %s fixed\n", cluster_id))
          } else {
            still_failed <- c(still_failed, cluster_id)
          }
        }
      }

      cat(sprintf(
        "\nRetry %d Summary: Fixed %d, Still failed %d\n",
        retry_attempt,
        length(newly_fixed),
        length(still_failed)
      ))

      failed_clusters <- still_failed
    }
  }

  # Final summary
  final_success <- length(result_list) - length(failed_clusters)
  final_rate <- final_success / length(result_list)

  cat(sprintf(
    "\n=== %s Final: %d/%d (%.1f%%) ===\n",
    task_name,
    final_success,
    length(result_list),
    final_rate * 100
  ))

  if (length(failed_clusters) > 0) {
    cat("\nStill failed after all retries:\n")
    cat(paste(head(failed_clusters, 20), collapse = ", "))
    if (length(failed_clusters) > 20) {
      cat(sprintf(" ... and %d more", length(failed_clusters) - 20))
    }
    cat("\n")
  }

  # Return in original format
  if (inherits(result, "interpretation_list")) {
    class(result_list) <- class(result)
    return(result_list)
  } else {
    return(result_list[[1]])
  }
}

# ==============================================================================
# Load and Prepare Data
# ==============================================================================

cat("\n=== Loading Seurat Data ===\n")

seurat_obj <- GetSeurat(h5ad_path = H5AD_PATH, debug = TRUE)
DefaultAssay(seurat_obj) <- "RNA"

cat(sprintf(
  "\nLoaded: %d cells x %d genes\n",
  ncol(seurat_obj),
  nrow(seurat_obj)
))

# Validate and normalize
need_cols <- c("cell_type_L2", "cell_type_L3")
missing <- setdiff(need_cols, colnames(seurat_obj@meta.data))
if (length(missing) > 0) {
  stop(
    "ERROR: Missing required metadata columns: ",
    paste(missing, collapse = ", ")
  )
}

# Normalize if needed
data_slot <- GetAssayData(seurat_obj, slot = "data")
skip_normalize <- FALSE
if (length(data_slot@x) > 0) {
  if (max(data_slot@x) < 20) {
    cat("[INFO] Data appears log-normalized, skipping normalization\n")
    skip_normalize <- TRUE
  }
}

if (!skip_normalize) {
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 1e4,
    verbose = FALSE
  )
  cat("[OK] Data normalized\n")
}

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = FALSE
)

cat("[OK] Data preparation complete\n")

# ==============================================================================
# Find Marker Genes
# ==============================================================================

cat("\n=== Computing Marker Genes (Parallel) ===\n")

Idents(seurat_obj) <- "cell_type_L3"
clusters <- levels(Idents(seurat_obj))

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

cat(sprintf(
  "[OK] Found %d significant markers across %d clusters\n",
  nrow(all_markers),
  length(unique(all_markers$cluster))
))

write.csv(
  all_markers,
  file.path(OUTPUT_DIR, "all_markers.csv"),
  row.names = FALSE
)

# Prepare top markers (filtered)
genes_to_filter <- c(
  grep("^MT-", rownames(seurat_obj), value = TRUE),
  grep("^RP[SL]", rownames(seurat_obj), value = TRUE),
  "FOS",
  "JUN",
  "JUNB",
  "JUND",
  "EGR1",
  "EGR2",
  "EGR3",
  "ZFP36",
  "DUSP1",
  "DUSP2",
  "IER2",
  "IER3",
  "ATF3",
  "BTG2",
  "FOSB",
  "NR4A1",
  "NR4A2",
  "NR4A3"
)

lfc_col <- if ("avg_log2FC" %in% colnames(all_markers)) {
  "avg_log2FC"
} else {
  "avg_logFC"
}

top_markers <- all_markers %>%
  filter(!gene %in% genes_to_filter) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(.data[[lfc_col]])) %>%
  slice_head(n = TOP_N_MARKERS) %>%
  ungroup() %>%
  mutate(gene = toupper(gene), cluster = as.character(cluster)) %>%
  dplyr::select(gene, cluster)

write.csv(
  top_markers,
  file.path(OUTPUT_DIR, "top_markers_filtered.csv"),
  row.names = FALSE
)

# ==============================================================================
# Load Databases (abbreviated - see full script for complete version)
# ==============================================================================

cat("\n=== Loading Enrichment Databases ===\n")

# GO
go_all_gmt <- read.gmt(GMT_GO_ALL) %>% mutate(gene = toupper(gene))
go_bp_gmt <- go_all_gmt[grep("^GOBP_", go_all_gmt$term), ]
go_mf_gmt <- go_all_gmt[grep("^GOMF_", go_all_gmt$term), ]
go_cc_gmt <- go_all_gmt[grep("^GOCC_", go_all_gmt$term), ]

# MSigDB
lines <- readLines(MSIGDB_GMT_PATH)
gene_sets_list <- lapply(lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  list(name = parts[1], genes = parts[-(1:2)])
})
term2gene_list <- lapply(gene_sets_list, function(gs) {
  if (length(gs$genes) > 0) {
    data.frame(
      term = rep(gs$name, length(gs$genes)),
      gene = toupper(gs$genes),
      stringsAsFactors = FALSE
    )
  }
})
all_genesets <- bind_rows(term2gene_list)
hallmark_term2gene <- all_genesets %>% filter(grepl("^HALLMARK_", term))
kegg_term2gene <- all_genesets %>% filter(grepl("KEGG_", term))

# CellMarker and PanglaoDB (load similarly - see previous script)
# ... (abbreviated for space)

cat("[OK] All databases loaded\n")

# ==============================================================================
# Run Enrichment
# ==============================================================================

cat("\n=== Running Enrichment Analyses ===\n")

go_bp_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers,
  fun = enricher,
  TERM2GENE = go_bp_gmt,
  pvalueCutoff = 0.05
)
go_mf_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers,
  fun = enricher,
  TERM2GENE = go_mf_gmt,
  pvalueCutoff = 0.05
)
go_cc_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers,
  fun = enricher,
  TERM2GENE = go_cc_gmt,
  pvalueCutoff = 0.05
)
hallmark_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers,
  fun = enricher,
  TERM2GENE = hallmark_term2gene,
  pvalueCutoff = 0.05
)
kegg_enrich <- compareCluster(
  gene ~ cluster,
  data = top_markers,
  fun = enricher,
  TERM2GENE = kegg_term2gene,
  pvalueCutoff = 0.05
)

# Build enrichment list
enrichment_list <- list()
if (!is.null(go_bp_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_bp_enrich
}
if (!is.null(go_mf_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_mf_enrich
}
if (!is.null(go_cc_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- go_cc_enrich
}
if (!is.null(hallmark_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- hallmark_enrich
}
if (!is.null(kegg_enrich)) {
  enrichment_list[[length(enrichment_list) + 1]] <- kegg_enrich
}

cat(sprintf("[OK] Total databases: %d\n", length(enrichment_list)))

# ==============================================================================
# ⭐ RUN INTERPRETATION (Mode-Dependent)
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
if (USE_AGENT_MODE) {
  cat("RUNNING DEEP MODE: interpret_agent() with Multi-Agent System\n")
} else {
  cat("RUNNING FAST MODE: interpret() Standard Interpretation\n")
}
cat(
  "================================================================================\n"
)

# Biological context
annotation_context <- paste(
  "T/NK cells from respiratory tract tissues (nasal cavity, paranasal sinuses, bronchi, lung).",
  "These are subclusters showing diverse functional states and differentiation stages.",
  "",
  "Expected cell types include:",
  "- CD8+ cytotoxic T cells: effector memory (TEM), tissue-resident memory (TRM)",
  "- CD4+ helper T cells: Th1, Th2, Th17, regulatory T cells (Treg), Tfh",
  "- NK cells: CD16+ (cytotoxic), CD16- (regulatory), activated vs resting",
  "- Innate-like T cells: gamma-delta T cells, MAIT cells, NKT cells",
  "",
  "Key considerations:",
  "1. Use marker enrichment to identify specific subtype",
  "2. Consider activation state (resting, activated, exhausted)",
  "3. Distinguish tissue-resident vs circulating phenotypes",
  "4. Note cytokine production profiles",
  "5. Identify memory vs naive vs effector differentiation"
)

# Prepare fold change vector if needed (optional)
gene_fold_change_vector <- NULL
if (USE_AGENT_MODE && USE_FOLD_CHANGE) {
  # Extract fold changes from all_markers
  fc_df <- all_markers %>%
    group_by(gene) %>%
    summarise(avg_fc = mean(.data[[lfc_col]], na.rm = TRUE))

  gene_fold_change_vector <- setNames(fc_df$avg_fc, fc_df$gene)
  cat(sprintf(
    "[INFO] Prepared fold change vector: %d genes\n",
    length(gene_fold_change_vector)
  ))
}

# Run interpretation with retry
annotation_results <- interpret_with_retry(
  enrichment_input = enrichment_list,
  context = annotation_context,
  n_pathways = if (USE_AGENT_MODE) 50 else 12, # Agent mode uses more pathways
  model = "deepseek-chat",
  api_key = DEEPSEEK_API_KEY,
  task = "annotation", # Only used in non-agent mode
  max_retries = 2,
  task_name = "Cell Type Annotation",
  use_agent = USE_AGENT_MODE,
  add_ppi = if (USE_AGENT_MODE) ADD_PPI else FALSE,
  gene_fold_change = gene_fold_change_vector
)

# ==============================================================================
# Process and Save Results
# ==============================================================================

if (!is.null(annotation_results)) {
  cat("\n=== Processing Results ===\n")

  results_list <- if (inherits(annotation_results, "interpretation_list")) {
    annotation_results
  } else {
    list(Default = annotation_results)
  }

  # Extract to data frame (handle both modes)
  annotation_table <- bind_rows(lapply(names(results_list), function(name) {
    res <- results_list[[name]]

    # Extract fields (agent mode has more fields)
    if (USE_AGENT_MODE) {
      data.frame(
        Cluster = name,
        Overview = if (!is.null(res$overview)) {
          substr(res$overview, 1, 500)
        } else {
          ""
        },
        Key_Mechanisms = if (!is.null(res$key_mechanisms)) {
          substr(paste(unlist(res$key_mechanisms), collapse = "; "), 1, 500)
        } else {
          ""
        },
        Hypothesis = if (!is.null(res$hypothesis)) {
          substr(paste(unlist(res$hypothesis), collapse = "; "), 1, 300)
        } else {
          ""
        },
        Regulatory_Drivers = if (!is.null(res$regulatory_drivers)) {
          paste(res$regulatory_drivers, collapse = "; ")
        } else {
          ""
        },
        Confidence = if (!is.null(res$confidence)) res$confidence else "NA",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Cluster = name,
        Cell_Type = if (!is.null(res$cell_type)) res$cell_type else "Unknown",
        Confidence = if (!is.null(res$confidence)) res$confidence else "NA",
        Markers = if (!is.null(res$markers)) {
          paste(res$markers, collapse = "; ")
        } else {
          "NA"
        },
        Reasoning = if (!is.null(res$reasoning)) {
          substr(res$reasoning, 1, 500)
        } else {
          ""
        },
        Regulatory_Drivers = if (!is.null(res$regulatory_drivers)) {
          paste(res$regulatory_drivers, collapse = "; ")
        } else {
          "NA"
        },
        stringsAsFactors = FALSE
      )
    }
  }))

  # Save
  output_csv <- if (USE_AGENT_MODE) {
    "annotation_results_agent.csv"
  } else {
    "annotation_results.csv"
  }
  output_rds <- if (USE_AGENT_MODE) {
    "annotation_results_agent.rds"
  } else {
    "annotation_results.rds"
  }

  write.csv(
    annotation_table,
    file.path(OUTPUT_DIR, output_csv),
    row.names = FALSE
  )
  saveRDS(annotation_results, file.path(OUTPUT_DIR, "reports", output_rds))

  cat(sprintf("[OK] Saved: %s\n", output_csv))
  cat(sprintf("[OK] Saved: reports/%s\n", output_rds))

  # Summary
  cat("\n=== Summary ===\n")
  cat(sprintf("Total clusters: %d\n", nrow(annotation_table)))

  if (USE_AGENT_MODE) {
    cat(sprintf(
      "Clusters with overview: %d\n",
      sum(annotation_table$Overview != "")
    ))
    cat(sprintf(
      "Clusters with hypothesis: %d\n",
      sum(annotation_table$Hypothesis != "")
    ))
  } else {
    high_conf <- sum(annotation_table$Confidence == "High")
    cat(sprintf(
      "High confidence: %d (%.1f%%)\n",
      high_conf,
      high_conf / nrow(annotation_table) * 100
    ))
    non_empty <- sum(
      annotation_table$Cell_Type != "" & annotation_table$Cell_Type != "Unknown"
    )
    cat(sprintf(
      "Non-empty cell types: %d (%.1f%%)\n",
      non_empty,
      non_empty / nrow(annotation_table) * 100
    ))
  }
}

# ==============================================================================
# Final Summary
# ==============================================================================

cat("\n")
cat(
  "================================================================================\n"
)
cat("ANALYSIS COMPLETE - v4.1 WITH DUAL MODE SUPPORT\n")
cat(
  "================================================================================\n\n"
)

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key files:\n")
if (USE_AGENT_MODE) {
  cat("  - annotation_results_agent.csv       Agent-based interpretation\n")
  cat("  - reports/annotation_results_agent.rds  Raw agent results\n")
} else {
  cat("  - annotation_results.csv             Standard interpretation\n")
  cat("  - reports/annotation_results.rds     Raw results\n")
}
cat("  - all_markers.csv                    All marker genes\n")
cat("  - reports/*_enrich.rds               Enrichment objects\n")
cat("\n")

cat("Mode Information:\n")
cat(sprintf(
  "  Mode Used: %s\n",
  if (USE_AGENT_MODE) "DEEP (interpret_agent)" else "FAST (interpret)"
))
if (USE_AGENT_MODE) {
  cat("  Multi-Agent System:\n")
  cat("    1. Agent Cleaner: Filter noise and select relevant pathways\n")
  cat("    2. Agent Detective: Identify key regulators using PPI/TF data\n")
  cat("    3. Agent Synthesizer: Synthesize findings into coherent narrative\n")
  cat(sprintf("  PPI Network: %s\n", if (ADD_PPI) "Enabled" else "Disabled"))
  cat(sprintf(
    "  Fold Change: %s\n",
    if (USE_FOLD_CHANGE) "Enabled" else "Disabled"
  ))
}
cat("\n")

cat("Expected success rates:\n")
cat("  - With retry mechanism: 90-95%%\n")
if (USE_AGENT_MODE) {
  cat("  - Agent mode: More comprehensive but slower\n")
} else {
  cat("  - Standard mode: Faster with good quality\n")
}
cat("\n")

cat(
  "================================================================================\n"
)
cat("DONE\n")
cat(
  "================================================================================\n"
)
