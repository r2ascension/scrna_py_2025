#!/usr/bin/env Rscript
# ==============================================================================
# B Cell Tissue Comparison v2.6.0 - LLM Resume / Rebuild Runner
# ==============================================================================
#
# Purpose:
#   Reuse the validated 2026-04-08 B-cell downstream outputs, but rebuild the
#   LLM-relevant evidence under stricter 2026-04-10 rules:
#     - IG / MT / RP / ENSG / LINC-like genes are excluded from enrichment input
#       and LLM top-gene evidence.
#     - Energy / respiration pathways are not counted as top pathways.
#     - Comparison-style LLM tasks receive integrated up+down evidence together.
#     - ssGSEA evidence is aggregated across databases/methods instead of being
#       interpreted database-by-database.
#     - Output files are written back to the same locations as the previous run.
#
# Scope:
#   - reuse: final Seurat object, pseudobulk DE RDS, previous CHOIR cluster DE RDS
#   - rebuild: pairwise enrichment + integrated LLM, grouped ssGSEA + LLM,
#              CHOIR ssGSEA + OFA enrichment + integrated LLM
#   - preserve: visualizations / DE tables / final object paths
#
# Run:
#   Rscript /home/h2048/script/R/bcell_tissue_comparison_v2_6_20260410_llm_resume.R
#
# Date: 2026-04-10
# ==============================================================================

LINEAGE_TAG           <- "BCELL"
LINEAGE_DISPLAY       <- "B Cell"
LINEAGE_CONTEXT_LABEL <- "B cells and plasma cells"
LINEAGE_CONTEXT_LOWER <- "B-cell"
PIPELINE_VERSION_LABEL <- "v2.6.0-LLM-Resume-20260410"
GENERATED_BY_LABEL     <- "bcell_tissue_comparison_v2_6_20260410_llm_resume.R"
FINAL_FILE_PREFIX      <- "bcell_tissue_comparison_final"
H5AD_PATH              <- "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
OUTPUT_DIR             <- "/home/h2048/data/R/0408/bcell_tissue_comparison_v2_6_20260408"
PREVIOUS_OUTPUT_DIR    <- OUTPUT_DIR

INITIALIZE_ONLY <- TRUE
PIPELINE_TEST_MODE <- FALSE
USE_EXISTING_L2 <- FALSE
FILTER_IG_GENES_FOR_LLM_AND_ENRICHMENT <- TRUE
FILTER_TECHNICAL_GENES_FOR_LLM_AND_ENRICHMENT <- TRUE
DEPRIORITIZE_ENERGY_PATHWAYS <- TRUE
APPEND_ENERGY_PATHWAYS_AFTER_TOP <- FALSE
LLM_TOP_DEG_N <- 12L
INTERPRET_MULTI_DB_TERMS_PER_DB <- 5L
INTERPRET_MULTI_DB_MAX_DBS <- 6L
INTERPRET_SSGSEA_TERMS_PER_DB <- 10L
LLM_SSGSEA_TERMS_PER_DIRECTION <- 8L
LLM_EXTRA_RULES <- c(
  "Ignore MT-, IG, RPS, RPL, MRPS, MRPL, RP-, ENSG, and LINC genes for pathway computation and top-gene evidence.",
  "For all comparison-style tasks, integrate both directions into a single judgment instead of writing separate up/down interpretations.",
  "Energy / respiration pathways are not considered top pathways; continue selecting the next non-energy pathways."
)

source("/home/h2048/script/R/tissue_comparison_advanced_helper_20260408.R")

tryCatch(
  source("/home/h2048/script/R/bcell_tissue_comparison_v2_6_20260406.R"),
  bcell_pipeline_init_only = function(e) invisible(NULL)
)

if (!ENABLE_LLM) {
  stop("DEEPSEEK_API_KEY not available in the current environment; LLM rerun cannot proceed.")
}

cat("\n=== B-cell LLM Resume / Rebuild Runner ===\n")
cat(sprintf("Output dir: %s\n", OUTPUT_DIR))
cat(sprintf("Filtering IG genes: %s\n", tc_should_filter_ig_genes()))
cat(sprintf("Filtering technical genes: %s\n", tc_should_filter_technical_genes()))
cat(sprintf("Energy pathways excluded from top selection: %s\n\n", !isTRUE(APPEND_ENERGY_PATHWAYS_AFTER_TOP)))

previous_run <- tc_read_previous_run(PREVIOUS_OUTPUT_DIR, include_rds = TRUE, include_markdown = FALSE)
if (!isTRUE(previous_run$exists)) stop(sprintf("Previous output dir not found: %s", PREVIOUS_OUTPUT_DIR))

tc_require_rds <- function(name, rel_path = NULL) {
  x <- previous_run$rds[[name]]
  if (!is.null(x)) return(x)
  if (!is.null(rel_path)) {
    full_path <- file.path(PREVIOUS_OUTPUT_DIR, rel_path)
    if (!file.exists(full_path)) stop(sprintf("Missing previous RDS file: %s", full_path))
    x <- tryCatch(readRDS(full_path), error = function(e) NULL)
    if (!is.null(x)) return(x)
  }
  stop(sprintf("Missing previous RDS object: %s", name))
  x
}

obj_path <- file.path(OUTPUT_DIR, sprintf("%s.rds", FINAL_FILE_PREFIX))
if (!file.exists(obj_path)) stop(sprintf("Final Seurat object not found: %s", obj_path))
obj <- readRDS(obj_path)
if (!inherits(obj, "Seurat")) stop("Loaded final object is not a Seurat object.")

fig_dir   <- file.path(OUTPUT_DIR, "figures")
rpt_dir   <- file.path(OUTPUT_DIR, "reports")
de_dir    <- file.path(OUTPUT_DIR, "pseudobulk_de")
wx_dir    <- file.path(OUTPUT_DIR, "wilcox_exploratory")
de_l3_dir <- file.path(OUTPUT_DIR, "pseudobulk_de_L3")
wx_l3_dir <- file.path(OUTPUT_DIR, "wilcox_exploratory_L3")
for (d in c(fig_dir, rpt_dir, de_dir, wx_dir, de_l3_dir, wx_l3_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

meta <- obj@meta.data
tissues <- sort(unique(na.omit(as.character(meta[[TISSUE_COL]]))))
l2_types <- sort(unique(na.omit(as.character(meta[[CELLTYPE_L2_COL]]))))
l3_types <- sort(unique(na.omit(as.character(meta[[CELLTYPE_L3_COL]]))))
l3_l2_mapping_summary <- as.data.frame(table(
  L3 = obj@meta.data[[CELLTYPE_L3_COL]],
  L2 = obj@meta.data[[CELLTYPE_L2_COL]],
  useNA = "no"
), stringsAsFactors = FALSE) %>%
  dplyr::filter(Freq > 0) %>%
  dplyr::arrange(L3, dplyr::desc(Freq), L2)
umap_reduction <- pick_reduction(obj, UMAP_REDUCTION_PREFERRED)
choir_col <- grep("^CHOIR_clusters", colnames(obj@meta.data), value = TRUE)
choir_col <- if (length(choir_col) > 0) choir_col[[1]] else NULL

pb_de_all_prev        <- tc_require_rds("pseudobulk_de_all", "reports/pseudobulk_de_all.rds")
pb_de_l3_all_prev     <- tc_require_rds("pseudobulk_de_L3_all", "reports/pseudobulk_de_L3_all.rds")
wilcox_all_prev       <- tc_require_rds("wilcox_exploratory_all", "reports/wilcox_exploratory_all.rds")
wilcox_l3_all_prev    <- tc_require_rds("wilcox_exploratory_L3_all", "reports/wilcox_exploratory_L3_all.rds")
choir_ofa_prev        <- tc_require_rds("choir_ofa_all", "reports/choir/choir_ofa_all.rds")

rebuild_enrichment_from_de <- function(res, comp_dir, top_n = TOP_N_DEG_ENRICHMENT) {
  enr_by_direction <- list()
  for (dir_name in c("up", "down")) {
    genes <- if (dir_name == "up") {
      res$de_table %>%
        dplyr::filter(sig == "sig", log2FoldChange > 0) %>%
        dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
        utils::head(top_n) %>%
        dplyr::pull(gene)
    } else {
      res$de_table %>%
        dplyr::filter(sig == "sig", log2FoldChange < 0) %>%
        dplyr::arrange(log2FoldChange) %>%
        utils::head(top_n) %>%
        dplyr::pull(gene)
    }
    genes <- tc_filter_gene_symbols(genes)
    if (length(genes) < 5) {
      enr_by_direction[[dir_name]] <- NULL
      next
    }
    tested <- tc_filter_gene_symbols(res$tested_genes)
    enr_list <- list(
      GO_BP = run_gmt_enrichment(genes, go_bp_t2g, "GO_BP", tested),
      GO_MF = run_gmt_enrichment(genes, go_mf_t2g, "GO_MF", tested),
      GO_CC = run_gmt_enrichment(genes, go_cc_t2g, "GO_CC", tested),
      KEGG = run_gmt_enrichment(genes, kegg_t2g, "KEGG", tested),
      Hallmark = run_gmt_enrichment(genes, hallmark_t2g, "Hallmark", tested),
      CellMarker = run_gmt_enrichment(genes, cellmarker_t2g, "CellMarker", tested),
      PanglaoDB = run_gmt_enrichment(genes, panglaodb_t2g, "PanglaoDB", tested)
    )
    enr_list[[CUSTOM_DB_NAME]] <- run_gmt_enrichment(genes, custom_t2g, CUSTOM_DB_NAME, tested)
    enr_list <- enr_list[!vapply(enr_list, is.null, logical(1))]
    enr_by_direction[[dir_name]] <- enr_list

    enr_out <- file.path(comp_dir, paste0("enrichment_", dir_name))
    dir.create(enr_out, recursive = TRUE, showWarnings = FALSE)
    for (db in names(enr_list)) {
      er <- enr_list[[db]]
      er_df <- tryCatch(as.data.frame(er), error = function(e) data.frame())
      if (nrow(er_df) == 0) next
      data.table::fwrite(er_df, file.path(enr_out, paste0(db, ".csv")))
    }
    saveRDS(enr_list, file.path(enr_out, "all_enrichment.rds"))
  }
  enr_by_direction
}

rebuild_pairwise_level <- function(pb_de_level, group_col, level_name, out_dir) {
  enrich_level <- list()
  agent_level <- list()
  agent_structured_level <- list()
  for (group_name in names(pb_de_level)) {
    comps <- pb_de_level[[group_name]]
    if (is.null(comps) || length(comps) == 0) next
    group_l2_label <- resolve_group_l2_label(obj, group_col, group_name)
    enrich_level[[group_name]] <- list()
    agent_level[[group_name]] <- list()
    agent_structured_level[[group_name]] <- list()
    cat(sprintf("[PAIRWISE][%s] %s\n", level_name, group_name))
    for (comp_name in names(comps)) {
      res <- comps[[comp_name]]
      if (is.null(res) || is.null(res$de_table) || !is.data.frame(res$de_table) || nrow(res$de_table) == 0) next
      comp_dir <- file.path(out_dir, safe_name(group_name), safe_name(comp_name))
      dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(res$de_table, file.path(comp_dir, "DESeq2_results.csv"))
      enrich_level[[group_name]][[comp_name]] <- rebuild_enrichment_from_de(res, comp_dir)
      agent_level[[group_name]][[comp_name]] <- list()
      agent_structured_level[[group_name]][[comp_name]] <- list()

      parts <- strsplit(comp_name, "_vs_", fixed = TRUE)[[1]]
      if (length(parts) != 2) {
        cat(sprintf("  [WARN] unexpected comparison name: %s\n", comp_name))
        next
      }
      t2 <- parts[1]
      t1 <- parts[2]
      pairwise_llm_rec <- tryCatch(
        run_pairwise_integrated_llm(
          group_name = group_name,
          level_name = level_name,
          group_l2_label = group_l2_label,
          t1 = t1,
          t2 = t2,
          de_df = res$de_table,
          enrich_by_direction = enrich_level[[group_name]][[comp_name]],
          comp_dir = comp_dir
        ),
        error = function(e) {
          cat(sprintf("  [WARN] pairwise integrated LLM failed: %s\n", e$message))
          NULL
        }
      )
      if (!is.null(pairwise_llm_rec)) {
        agent_level[[group_name]][[comp_name]][["integrated_up_down"]] <- pairwise_llm_rec$raw_result
        agent_structured_level[[group_name]][[comp_name]][["integrated_up_down"]] <- pairwise_llm_rec
        cat(sprintf("  [OK] %s -> %s\n", comp_name, pairwise_llm_rec$status))
      }
      Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
    }
  }
  list(enrich = enrich_level, agent = agent_level, agent_structured = agent_structured_level)
}

run_grouped_ssgsea <- function(obj, group_col, level_name, methods, n_top, n_heatmap) {
  results_all <- list()
  group_key   <- paste0("ssgsea_group_", tolower(level_name))
  file_stem   <- level_file_stem("ssgsea", level_name)
  obj@meta.data[[group_key]] <- paste0(obj@meta.data[[TISSUE_COL]], "__", obj@meta.data[[group_col]])
  group_l2_map <- obj@meta.data %>%
    dplyr::transmute(
      group_id = .data[[group_key]],
      celltype_label = as.character(.data[[group_col]]),
      celltype_l2 = as.character(.data[[CELLTYPE_L2_COL]])
    ) %>%
    dplyr::group_by(group_id, celltype_label) %>%
    dplyr::summarise(
      celltype_l2 = {
        vals <- unique(celltype_l2[!is.na(celltype_l2) & nzchar(trimws(celltype_l2))])
        if (length(vals) == 1) vals else NA_character_
      },
      .groups = "drop"
    )
  ssgsea_groups <- sort(unique(obj@meta.data[[group_key]]))
  cat(sprintf("\n=== Grouped ssGSEA (average-expression; tissue x %s) ===\n", level_name))
  cat(sprintf("Methods: %s\n\n", paste(methods, collapse = ", ")))
  cat(sprintf("[INFO] ssGSEA groups (%d): %s\n\n", length(ssgsea_groups), paste(ssgsea_groups, collapse = ", ")))

  for (method in methods) {
    cat(sprintf("--- ssGSEA method [%s]: %s ---\n", level_name, method))
    selected_t2g <- switch(
      method,
      hallmark = hallmark_t2g,
      go_bp = go_bp_t2g,
      go_mf = go_mf_t2g,
      go_cc = go_cc_t2g,
      kegg = kegg_t2g,
      custom = if (!is.null(SSGSEA_CUSTOM_GMT) && file.exists(SSGSEA_CUSTOM_GMT)) {
        tryCatch(read.gmt(SSGSEA_CUSTOM_GMT) %>% dplyr::mutate(gene = toupper(gene)),
                 error = function(e) { cat("[WARN] custom GMT failed\n"); NULL })
      } else NULL,
      NULL
    )
    if (is.null(selected_t2g) || nrow(selected_t2g) == 0) {
      cat(sprintf("[WARN] Method '%s': no gene sets, skipping\n", method))
      next
    }
    features <- rownames(obj)
    gs_list  <- term2gene_to_list(selected_t2g)
    gs_list  <- map_gene_sets_to_features(gs_list, features)
    gs_list  <- filter_gs_size(gs_list)
    if (length(gs_list) == 0) {
      cat(sprintf("[WARN] Method '%s': no valid gene sets after filtering\n", method))
      next
    }
    cat(sprintf("[OK] %d gene sets prepared\n", length(gs_list)))
    needed_genes <- unique(unlist(gs_list, use.names = FALSE))
    cat(sprintf("[INFO] Computing group-average expression for %d genes...\n", length(needed_genes)))
    avg_expr <- tryCatch(
      AverageExpression(obj, assays = "RNA", slot = "data",
                        group.by = group_key, features = needed_genes, verbose = FALSE)[["RNA"]],
      error = function(e) { cat(sprintf("[ERROR] AverageExpression failed: %s\n", e$message)); NULL }
    )
    if (is.null(avg_expr)) {
      cat("[WARN] Skipping this method\n")
      next
    }
    cat("[INFO] Running ssGSEA...\n")
    bp <- BiocParallel::SnowParam(workers = N_CORES, type = "SOCK", progressbar = FALSE)
    ssgsea_scores <- tryCatch({
      if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
        param <- GSVA::ssgseaParam(exprData = as.matrix(avg_expr), geneSets = gs_list,
                                   alpha = 0.25, normalize = TRUE, minSize = 10, maxSize = 500)
        GSVA::gsva(param, BPPARAM = bp, verbose = FALSE)
      } else {
        GSVA::gsva(as.matrix(avg_expr), gs_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
      }
    }, error = function(e) { cat(sprintf("[ERROR] ssGSEA failed: %s\n", e$message)); NULL })
    if (is.null(ssgsea_scores)) {
      cat("[WARN] Skipping this method\n")
      next
    }
    ssgsea_z <- t(scale(t(ssgsea_scores)))
    cat(sprintf("[OK] ssGSEA done: %d pathways x %d groups\n", nrow(ssgsea_scores), ncol(ssgsea_scores)))

    results_all[[method]] <- list(scores = ssgsea_scores, z_scores = ssgsea_z, gene_sets = gs_list)
    saveRDS(ssgsea_scores, file.path(rpt_dir, sprintf("%s_scores_%s.rds", file_stem, method)))
    saveRDS(ssgsea_z,      file.path(rpt_dir, sprintf("%s_z_%s.rds", file_stem, method)))

    top_df <- extract_directional_ssgsea_top_rows(
      ssgsea_scores, ssgsea_z,
      top_n = n_top,
      group_field = "group",
      level_name = level_name,
      method = method
    )
    data.table::fwrite(top_df, file.path(rpt_dir, sprintf("%s_top_pathways_%s.csv", file_stem, method)))

    mean_abs_z <- rowMeans(abs(ssgsea_z), na.rm = TRUE)
    top_paths  <- names(sort(mean_abs_z, decreasing = TRUE))[seq_len(min(n_heatmap, length(mean_abs_z)))]
    heat_mat   <- ssgsea_z[top_paths, , drop = FALSE]
    heat_mat[is.nan(heat_mat) | is.infinite(heat_mat)] <- 0
    if (nrow(heat_mat) >= 3 && ncol(heat_mat) >= 2) {
      tryCatch({
        heat_path <- file.path(fig_dir, sprintf("%s_heatmap_%s", file_stem, method))
        ht <- pheatmap::pheatmap(
          heat_mat,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
          breaks = seq(-3, 3, length.out = 101),
          main = sprintf("%s ssGSEA Z-score (%s) — top %d pathways by tissue x %s",
                         LINEAGE_DISPLAY, method, nrow(heat_mat), level_name),
          fontsize_row = 7,
          fontsize_col = 9,
          cellwidth = 22,
          cellheight = 10,
          filename = paste0(heat_path, ".pdf"),
          width = max(10, ncol(heat_mat) * 1.8),
          height = max(8, nrow(heat_mat) * 0.35 + 3)
        )
        grDevices::png(paste0(heat_path, ".png"),
                       width = max(10, ncol(heat_mat) * 1.8),
                       height = max(8, nrow(heat_mat) * 0.35 + 3),
                       units = "in", res = 300)
        grid::grid.newpage(); grid::grid.draw(ht$gtable)
        grDevices::dev.off()
      }, error = function(e) cat(sprintf("[WARN] Heatmap failed for %s: %s\n", method, e$message)))
    }
    results_all[[method]][["top_df"]] <- top_df
    results_all[[method]][["llm"]] <- list()
    cat(sprintf("[OK] Method '%s' complete\n\n", method))
  }

  ssgsea_llm_combined <- list()
  ssgsea_method_names <- get_ssgsea_method_names(results_all)
  if (ENABLE_LLM && length(ssgsea_method_names) > 0) {
    cat(sprintf("[LLM-4] Interpreting ssGSEA groups with combined multi-database evidence (level=%s)...\n", level_name))
    ssgsea_llm_out <- file.path(rpt_dir, sprintf("ssgsea_llm_%s", tolower(level_name)))
    dir.create(ssgsea_llm_out, recursive = TRUE, showWarnings = FALSE)
    top_df_all <- dplyr::bind_rows(lapply(ssgsea_method_names, function(method) {
      results_all[[method]][["top_df"]]
    }))
    all_group_ids <- unique(top_df_all$group)
    for (gid in all_group_ids) {
      gid_df <- top_df_all[top_df_all$group == gid, , drop = FALSE]
      gid_meta <- group_l2_map[group_l2_map$group_id == gid, , drop = FALSE]
      gid_methods <- unique(as.character(gid_df$method))
      rec_gid <- tryCatch(
        run_ssgsea_group_llm(
          gid, gid_df, gid_methods, level_name, ssgsea_llm_out,
          celltype_l2 = if (nrow(gid_meta) >= 1) gid_meta$celltype_l2[1] else NA_character_
        ),
        error = function(e) {
          cat(sprintf("    [WARN] ssGSEA multi-db LLM failed for %s: %s\n", gid, e$message))
          NULL
        }
      )
      if (!is.null(rec_gid)) {
        ssgsea_llm_combined[[gid]] <- rec_gid
        cat(sprintf("    [OK] %s | dbs=%s | status=%s\n",
                    gid, paste(gid_methods, collapse = ", "), rec_gid$status))
      }
      Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
    }
    cat(sprintf("[OK] Combined ssGSEA LLM interpretations: %d / %d groups\n",
                length(ssgsea_llm_combined), length(all_group_ids)))
  }
  results_all[["llm_combined"]] <- ssgsea_llm_combined

  cat(sprintf("[SUMMARY] Completed %d ssGSEA methods for %s\n", length(ssgsea_method_names), level_name))
  for (mn in ssgsea_method_names) {
    cat(sprintf("  - %s: %d pathways x %d groups\n",
                mn, nrow(results_all[[mn]]$scores), ncol(results_all[[mn]]$scores)))
  }
  if (ENABLE_LLM) {
    cat(sprintf("  - combined multi-db ssGSEA LLM: %d records\n", length(ssgsea_llm_combined)))
  }
  cat("\n")
  results_all
}

run_choir_ssgsea_resume <- function(obj, choir_col, choir_dir) {
  ssgsea_choir_all <- list()
  if (is.null(choir_col) || !choir_col %in% colnames(obj@meta.data)) return(ssgsea_choir_all)
  cat("\n[CHOIR] Rebuilding ssGSEA per cluster\n")
  for (method in SSGSEA_CHOIR_METHODS) {
    cat(sprintf("  method: %s\n", method))
    selected_t2g_choir <- switch(
      method,
      hallmark = hallmark_t2g,
      go_bp = go_bp_t2g,
      go_mf = go_mf_t2g,
      go_cc = go_cc_t2g,
      kegg = kegg_t2g,
      custom = if (!is.null(SSGSEA_CUSTOM_GMT) && file.exists(SSGSEA_CUSTOM_GMT)) {
        tryCatch(read.gmt(SSGSEA_CUSTOM_GMT) %>% dplyr::mutate(gene = toupper(gene)), error = function(e) NULL)
      } else NULL,
      NULL
    )
    if (is.null(selected_t2g_choir) || nrow(selected_t2g_choir) == 0) next
    gs_choir <- term2gene_to_list(selected_t2g_choir)
    gs_choir <- map_gene_sets_to_features(gs_choir, rownames(obj))
    gs_choir <- filter_gs_size(gs_choir)
    if (length(gs_choir) == 0) next
    needed_genes_choir <- unique(unlist(gs_choir, use.names = FALSE))
    avg_choir <- tryCatch(
      AverageExpression(obj, assays = "RNA", slot = "data",
                        group.by = choir_col, features = needed_genes_choir, verbose = FALSE)[["RNA"]],
      error = function(e) {
        cat(sprintf("  [ERROR] AverageExpression %s: %s\n", method, e$message))
        NULL
      }
    )
    if (is.null(avg_choir)) next
    bp_choir <- BiocParallel::SnowParam(workers = N_CORES, type = "SOCK", progressbar = FALSE)
    ss_choir <- tryCatch({
      if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
        GSVA::gsva(
          GSVA::ssgseaParam(
            exprData = as.matrix(avg_choir),
            geneSets = gs_choir,
            alpha = 0.25,
            normalize = TRUE,
            minSize = 10,
            maxSize = 500
          ),
          BPPARAM = bp_choir,
          verbose = FALSE
        )
      } else {
        GSVA::gsva(as.matrix(avg_choir), gs_choir, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
      }
    }, error = function(e) {
      cat(sprintf("  [ERROR] ssGSEA %s: %s\n", method, e$message))
      NULL
    })
    if (is.null(ss_choir)) next
    sz_choir <- t(scale(t(ss_choir)))
    top_choir_df <- extract_directional_ssgsea_top_rows(
      ss_choir,
      sz_choir,
      top_n = SSGSEA_CHOIR_N_TOP,
      group_field = "cluster",
      method = method
    )
    ssgsea_choir_all[[method]] <- list(scores = ss_choir, z_scores = sz_choir, top_df = top_choir_df)
    saveRDS(ss_choir, file.path(choir_dir, sprintf("ssgsea_scores_%s.rds", method)))
    saveRDS(sz_choir, file.path(choir_dir, sprintf("ssgsea_z_%s.rds", method)))
    data.table::fwrite(top_choir_df, file.path(choir_dir, sprintf("ssgsea_top_pathways_%s.csv", method)))
  }
  ssgsea_choir_all
}

rebuild_choir_outputs <- function(obj, choir_ofa_prev, choir_col) {
  if (is.null(choir_col) || !choir_col %in% colnames(obj@meta.data)) {
    stop("CHOIR cluster column not found in final object.")
  }
  choir_dir <- file.path(rpt_dir, "choir")
  dir.create(choir_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.frame(cell = colnames(obj), choir_cluster = obj@meta.data[[choir_col]], stringsAsFactors = FALSE),
    file.path(choir_dir, "choir_clusters.csv")
  )
  ssgsea_choir_all <- run_choir_ssgsea_resume(obj, choir_col, choir_dir)
  choir_ofa_all <- list()
  cat("\n[CHOIR] Rebuilding OFA enrichment + integrated LLM\n")
  for (cl_name in names(choir_ofa_prev)) {
    prev_rec <- choir_ofa_prev[[cl_name]]
    de_ofa <- prev_rec$de_table
    if (is.null(de_ofa) || !is.data.frame(de_ofa) || nrow(de_ofa) == 0) next
    ofa_dir <- file.path(choir_dir, sprintf("ofa_%s", safe_name(cl_name)))
    dir.create(ofa_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(de_ofa, file.path(ofa_dir, "markers.csv"))
    all_tested_ofa <- tc_filter_gene_symbols(de_ofa$gene)
    ofa_enrich <- list()
    for (dir_name in c("up", "down")) {
      ofa_genes <- if (dir_name == "up") {
        de_ofa %>%
          dplyr::filter(padj < OFA_PADJ_THR, avg_log2FC > 0) %>%
          dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
          utils::head(OFA_TOP_N) %>%
          dplyr::pull(gene)
      } else {
        de_ofa %>%
          dplyr::filter(padj < OFA_PADJ_THR, avg_log2FC < 0) %>%
          dplyr::arrange(avg_log2FC) %>%
          utils::head(OFA_TOP_N) %>%
          dplyr::pull(gene)
      }
      ofa_genes <- tc_filter_gene_symbols(ofa_genes)
      if (length(ofa_genes) < 5) next
      enr_ofa <- list(
        GO_BP = run_gmt_enrichment(ofa_genes, go_bp_t2g, "GO_BP", all_tested_ofa),
        GO_MF = run_gmt_enrichment(ofa_genes, go_mf_t2g, "GO_MF", all_tested_ofa),
        GO_CC = run_gmt_enrichment(ofa_genes, go_cc_t2g, "GO_CC", all_tested_ofa),
        KEGG = run_gmt_enrichment(ofa_genes, kegg_t2g, "KEGG", all_tested_ofa),
        Hallmark = run_gmt_enrichment(ofa_genes, hallmark_t2g, "Hallmark", all_tested_ofa),
        CellMarker = run_gmt_enrichment(ofa_genes, cellmarker_t2g, "CellMarker", all_tested_ofa),
        PanglaoDB = run_gmt_enrichment(ofa_genes, panglaodb_t2g, "PanglaoDB", all_tested_ofa)
      )
      enr_ofa[[CUSTOM_DB_NAME]] <- run_gmt_enrichment(ofa_genes, custom_t2g, CUSTOM_DB_NAME, all_tested_ofa)
      enr_ofa <- enr_ofa[!vapply(enr_ofa, is.null, logical(1))]
      ofa_enrich[[dir_name]] <- enr_ofa
      enr_out_ofa <- file.path(ofa_dir, paste0("enrichment_", dir_name))
      dir.create(enr_out_ofa, recursive = TRUE, showWarnings = FALSE)
      for (db in names(enr_ofa)) {
        er <- enr_ofa[[db]]
        er_df <- tryCatch(as.data.frame(er), error = function(e) data.frame())
        if (nrow(er_df) == 0) next
        data.table::fwrite(er_df, file.path(enr_out_ofa, paste0(db, ".csv")))
      }
      saveRDS(enr_ofa, file.path(enr_out_ofa, "all_enrichment.rds"))
    }

    focal_cells <- colnames(obj)[obj@meta.data[[choir_col]] == cl_name]
    choir_cluster_l2 <- resolve_dominant_label(obj@meta.data[focal_cells, CELLTYPE_L2_COL])
    choir_llm_rec <- tryCatch(
      run_choir_cluster_llm(
        cluster_id = cl_name,
        de_df = de_ofa,
        ofa_enrich = ofa_enrich,
        ssgsea_choir_all = ssgsea_choir_all,
        choir_cluster_l2 = choir_cluster_l2
      ),
      error = function(e) {
        cat(sprintf("  [WARN] CHOIR %s LLM failed: %s\n", cl_name, e$message))
        NULL
      }
    )
    choir_llm <- list()
    if (!is.null(choir_llm_rec)) {
      choir_llm[["integrated"]] <- choir_llm_rec
      saveRDS(choir_llm_rec, file.path(ofa_dir, "llm_integrated_structured.rds"))
      if (nzchar(choir_llm_rec$raw_text)) {
        writeLines(choir_llm_rec$raw_text, file.path(ofa_dir, "llm_integrated_raw.txt"))
      }
      cat(sprintf("  [OK] cluster %s -> %s\n", cl_name, choir_llm_rec$status))
    }
    choir_ofa_all[[cl_name]] <- list(
      de_table = de_ofa,
      enrich = ofa_enrich,
      llm = choir_llm,
      n_focal = if (!is.null(prev_rec$n_focal)) prev_rec$n_focal else length(focal_cells),
      n_rest = if (!is.null(prev_rec$n_rest)) prev_rec$n_rest else (ncol(obj) - length(focal_cells))
    )
    Sys.sleep(STANDARDIZE_LLM_RETRY_SLEEP_SEC)
  }
  choir_results_all <- list(ssgsea = ssgsea_choir_all, ofa = choir_ofa_all)
  saveRDS(choir_results_all, file.path(choir_dir, "choir_results_all.rds"))
  saveRDS(choir_ofa_all, file.path(choir_dir, "choir_ofa_all.rds"))
  list(choir_results_all = choir_results_all, choir_ofa_all = choir_ofa_all)
}

pairwise_l2 <- rebuild_pairwise_level(pb_de_all_prev, CELLTYPE_L2_COL, "L2", de_dir)
pairwise_l3 <- rebuild_pairwise_level(pb_de_l3_all_prev, CELLTYPE_L3_COL, "L3", de_l3_dir)

enrich_all <- pairwise_l2$enrich
enrich_l3_all <- pairwise_l3$enrich
agent_all <- pairwise_l2$agent
agent_l3_all <- pairwise_l3$agent
agent_structured_all <- pairwise_l2$agent_structured
agent_structured_l3_all <- pairwise_l3$agent_structured

ssgsea_results_all <- run_grouped_ssgsea(
  obj = obj,
  group_col = CELLTYPE_L2_COL,
  level_name = "L2",
  methods = SSGSEA_METHODS,
  n_top = SSGSEA_N_TOP,
  n_heatmap = SSGSEA_N_HEATMAP
)
ssgsea_results_l3_all <- run_grouped_ssgsea(
  obj = obj,
  group_col = CELLTYPE_L3_COL,
  level_name = "L3",
  methods = SSGSEA_METHODS,
  n_top = SSGSEA_N_TOP,
  n_heatmap = SSGSEA_N_HEATMAP
)
saveRDS(ssgsea_results_all, file.path(rpt_dir, "ssgsea_results_all.rds"))
saveRDS(ssgsea_results_l3_all, file.path(rpt_dir, "ssgsea_results_L3_all.rds"))

choir_bundle <- rebuild_choir_outputs(obj, choir_ofa_prev, choir_col)
choir_results_all <- choir_bundle$choir_results_all
choir_ofa_all <- choir_bundle$choir_ofa_all

saveRDS(pb_de_all_prev, file.path(rpt_dir, "pseudobulk_de_all.rds"))
saveRDS(pb_de_l3_all_prev, file.path(rpt_dir, "pseudobulk_de_L3_all.rds"))
saveRDS(wilcox_all_prev, file.path(rpt_dir, "wilcox_exploratory_all.rds"))
saveRDS(wilcox_l3_all_prev, file.path(rpt_dir, "wilcox_exploratory_L3_all.rds"))
saveRDS(enrich_all, file.path(rpt_dir, "enrichment_all.rds"))
saveRDS(enrich_l3_all, file.path(rpt_dir, "enrichment_L3_all.rds"))
saveRDS(agent_all, file.path(rpt_dir, "interpret_agent_all.rds"))
saveRDS(agent_l3_all, file.path(rpt_dir, "interpret_agent_L3_all.rds"))
saveRDS(agent_structured_all, file.path(rpt_dir, "interpret_agent_structured_all.rds"))
saveRDS(agent_structured_l3_all, file.path(rpt_dir, "interpret_agent_structured_L3_all.rds"))
data.table::fwrite(l3_l2_mapping_summary, file.path(rpt_dir, "L3_to_L2_mapping_summary.csv"))
saveRDS(L3_TO_L2_REMAP, file.path(rpt_dir, "L3_to_L2_remap.rds"))

agent_structured_df_l2 <- flatten_interpretation_records(agent_structured_all)
agent_structured_df_l3 <- flatten_interpretation_records(agent_structured_l3_all)
agent_structured_df <- dplyr::bind_rows(agent_structured_df_l2, agent_structured_df_l3)
if (nrow(agent_structured_df) > 0) {
  data.table::fwrite(agent_structured_df, file.path(rpt_dir, "interpret_agent_structured.tsv"), sep = "\t")
}
if (nrow(agent_structured_df_l3) > 0) {
  data.table::fwrite(agent_structured_df_l3, file.path(rpt_dir, "interpret_agent_structured_L3.tsv"), sep = "\t")
}
write_interpretation_markdown(
  agent_structured_df,
  file.path(OUTPUT_DIR, "LLM_INTERPRETATION.md"),
  "# LLM Interpretation Summary (B Cell tissue comparison; resumed 2026-04-10)",
  include_raw = FALSE
)
write_interpretation_markdown(
  agent_structured_df,
  file.path(rpt_dir, "LLM_INTERPRETATION_FOR_LLM.md"),
  "# LLM Interpretation Structured Input",
  include_raw = TRUE
)

ssgsea_llm_df_l2 <- collect_ssgsea_llm_records(ssgsea_results_all)
ssgsea_llm_df_l3 <- collect_ssgsea_llm_records(ssgsea_results_l3_all)
ssgsea_llm_df <- dplyr::bind_rows(ssgsea_llm_df_l2, ssgsea_llm_df_l3)
saveRDS(ssgsea_llm_df_l2, file.path(rpt_dir, "ssgsea_llm_structured_L2_all.rds"))
saveRDS(ssgsea_llm_df_l3, file.path(rpt_dir, "ssgsea_llm_structured_L3_all.rds"))
saveRDS(ssgsea_llm_df, file.path(rpt_dir, "ssgsea_llm_structured_all.rds"))
if (nrow(ssgsea_llm_df) > 0) {
  data.table::fwrite(ssgsea_llm_df, file.path(rpt_dir, "ssgsea_llm_structured.tsv"), sep = "\t")
  write_interpretation_markdown(
    ssgsea_llm_df,
    file.path(OUTPUT_DIR, "LLM_SSGSEA_INTERPRETATION.md"),
    "# LLM ssGSEA Interpretation (B Cell tissue comparison; resumed 2026-04-10)",
    include_raw = FALSE
  )
}

choir_llm_df <- collect_choir_llm_records(choir_ofa_all)
saveRDS(choir_llm_df, file.path(rpt_dir, "choir_llm_structured_all.rds"))
if (nrow(choir_llm_df) > 0) {
  data.table::fwrite(choir_llm_df, file.path(rpt_dir, "choir_llm_structured.tsv"), sep = "\t")
  write_interpretation_markdown(
    choir_llm_df,
    file.path(OUTPUT_DIR, "LLM_CHOIR_INTERPRETATION.md"),
    "# LLM CHOIR Cluster Interpretation (B Cell tissue comparison; resumed 2026-04-10)",
    include_raw = FALSE
  )
}

rerun_report <- c(
  "# B Cell LLM Resume / Rebuild Report (2026-04-10)",
  "",
  sprintf("- Output directory: `%s`", OUTPUT_DIR),
  sprintf("- Final object reused: `%s`", basename(obj_path)),
  sprintf("- Pairwise LLM records: %d", nrow(agent_structured_df)),
  sprintf("- ssGSEA LLM records: %d", nrow(ssgsea_llm_df)),
  sprintf("- CHOIR LLM records: %d", nrow(choir_llm_df)),
  "",
  "## Rule changes applied",
  "",
  "- IG / MT / RP / ENSG / LINC-like genes are excluded from enrichment input and LLM top-gene evidence.",
  "- Energy / respiration pathways are skipped during top-pathway selection.",
  "- Pairwise and CHOIR comparison evidence is integrated across both directions before LLM interpretation.",
  "- ssGSEA evidence is aggregated across databases / methods instead of being interpreted database-by-database.",
  "",
  "## Updated outputs",
  "",
  "- `LLM_INTERPRETATION.md`",
  "- `LLM_SSGSEA_INTERPRETATION.md`",
  "- `LLM_CHOIR_INTERPRETATION.md`",
  "- `reports/interpret_agent_structured.tsv`",
  "- `reports/ssgsea_llm_structured.tsv`",
  "- `reports/choir_llm_structured.tsv`",
  "- `reports/enrichment_all.rds` / `reports/enrichment_L3_all.rds`",
  "- `reports/ssgsea_results_all.rds` / `reports/ssgsea_results_L3_all.rds`",
  "- `reports/choir/choir_results_all.rds` / `reports/choir/choir_ofa_all.rds`"
)
writeLines(rerun_report, file.path(OUTPUT_DIR, "REPORT_LLM_RERUN_20260410.md"))

cat("\n=== Resume / rebuild complete ===\n")
cat(sprintf("Pairwise LLM rows: %d\n", nrow(agent_structured_df)))
cat(sprintf("ssGSEA LLM rows: %d\n", nrow(ssgsea_llm_df)))
cat(sprintf("CHOIR LLM rows: %d\n", nrow(choir_llm_df)))
cat(sprintf("Updated report: %s\n", file.path(OUTPUT_DIR, "REPORT_LLM_RERUN_20260410.md")))
