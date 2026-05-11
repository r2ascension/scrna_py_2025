#!/usr/bin/env Rscript
# ==============================================================================
# Myeloid Tissue Comparison v1.2.2 - Resume Finalizer
# ==============================================================================
#
# Purpose:
#   - complete the interrupted v1.2.2 myeloid run in-place
#   - reuse finished DE / enrichment / ssGSEA / CHOIR / OFA artifacts
#   - rerun only the final discovery-screen LLM layer, rebuild REPORT.md,
#     and export the final Seurat / h5ad objects
#
# Run:
#   Rscript /home/h2048/script/R/myeloid_tissue_comparison_v1_2_2_resume_20260415.R
#
# Date: 2026-04-15
# ==============================================================================

source("/home/h2048/script/R/tissue_comparison_generic_wrapper_20260414_v3.R")

PIPELINE_CONFIG <- tc_build_generic_tissue_comparison_config(
  lineage = "MYELOID",
  overrides = list(
    OUTPUT_DIR = "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414",
    PREVIOUS_OUTPUT_DIR = "/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414",
    ALLOW_IN_PLACE_RESUME = TRUE,
    REQUIRE_PREVIOUS_OUTPUT_DIR = TRUE,
    REUSE_PREVIOUS_FINAL_OBJECT = FALSE,
    REUSE_PREVIOUS_OUTPUT_SUMMARY = TRUE
  )
)

loaded_env_files <- tc_load_env_candidates(PIPELINE_CONFIG$ENV_FILE_CANDIDATES)
llm_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
has_llm_key <- nchar(llm_key) >= 10 && !tc_wrapper_is_placeholder_secret(llm_key)
if (!has_llm_key) {
  stop(paste(
    "DEEPSEEK_API_KEY not found or is still a placeholder value in loaded .env files.",
    "Resume finalizer requires live LLM access because the interrupted stage is discovery-screen LLM output."
  ))
}
if (length(loaded_env_files) > 0) {
  cat(sprintf("[OK] Loaded environment file(s): %s\n", paste(loaded_env_files, collapse = ", ")))
}

tc_apply_named_list(PIPELINE_CONFIG, envir = environment())
invisible(tc_apply_advanced_shared_overrides(envir = environment(), overrides = SHARED_OVERRIDES))

GENERATED_BY_LABEL <- "myeloid_tissue_comparison_v1_2_2_resume_20260415.R"
PIPELINE_SUBTITLE <- paste(PIPELINE_SUBTITLE, "[resume finalize 2026-04-15]")
PIPELINE_CHANGELOG_TITLE <- sprintf("%s Resume Finalizer:", PIPELINE_VERSION_LABEL)
PIPELINE_CHANGELOG_LINES <- c(
  PIPELINE_CHANGELOG_LINES,
  "  [RESUME-1] Reused completed v1.2.2 artifacts to rerun only final discovery-screen LLM batches.",
  "  [RESUME-2] Rebuilt REPORT.md and exported final RDS / h5ad without recomputing pseudobulk DE, ssGSEA, or CHOIR."
)

INITIALIZE_ONLY <- TRUE
tryCatch(
  source(SHARED_ENGINE_PATH, local = environment()),
  bcell_pipeline_init_only = function(e) invisible(NULL)
)
INITIALIZE_ONLY <- FALSE

cat("\n=== Myeloid v1.2.2 Resume Finalizer ===\n")
cat(sprintf("Output dir: %s\n", OUTPUT_DIR))
cat(sprintf("Shared engine: %s\n", SHARED_ENGINE_PATH))
cat(sprintf("LLM enabled: %s\n\n", ENABLE_LLM))

resume_require_file <- function(path) {
  if (!file.exists(path)) stop(sprintf("Required file not found: %s", path))
  path
}

resume_require_rds <- function(rel_path) {
  full_path <- resume_require_file(file.path(OUTPUT_DIR, rel_path))
  readRDS(full_path)
}

resume_require_table <- function(rel_path, sep = "\t") {
  full_path <- resume_require_file(file.path(OUTPUT_DIR, rel_path))
  data.table::fread(full_path, sep = sep, data.table = FALSE)
}

obj <- NULL
final_rds_existing <- file.path(OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".rds"))
if (file.exists(final_rds_existing)) {
  obj <- safe_read_rds(final_rds_existing)
  if (!is.null(obj)) cat(sprintf("[OK] Loaded existing final object: %s\n", final_rds_existing))
}
if (is.null(obj)) {
  if (!file.exists(H5AD_PATH)) stop(sprintf("File not found: %s", H5AD_PATH))
  obj <- GetSeurat(
    h5ad_path = H5AD_PATH,
    prefer_raw = FALSE,
    prefer_layer_counts = TRUE,
    validate_counts = TRUE,
    debug = TRUE
  )
}
cat(sprintf("[OK] %d cells x %d genes\n", ncol(obj), nrow(obj)))
if (!"counts" %in% Layers(obj[["RNA"]])) stop("RNA assay missing 'counts' layer.")
cat("[OK] counts layer verified\n\n")

fig_dir   <- file.path(OUTPUT_DIR, "figures")
rpt_dir   <- file.path(OUTPUT_DIR, "reports")
de_dir    <- file.path(OUTPUT_DIR, "pseudobulk_de")
wx_dir    <- file.path(OUTPUT_DIR, "wilcox_exploratory")
de_l3_dir <- file.path(OUTPUT_DIR, "pseudobulk_de_L3")
wx_l3_dir <- file.path(OUTPUT_DIR, "wilcox_exploratory_L3")
for (d in c(fig_dir, rpt_dir, de_dir, wx_dir, de_l3_dir, wx_l3_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("=== Validating Metadata & L2/L3 Annotation ===\n")
meta <- obj@meta.data
required_cols <- c(TISSUE_COL, SAMPLE_COL, L3_SOURCE_COL)
if (USE_EXISTING_L2) required_cols <- c(required_cols, L2_SOURCE_COL)
for (col in required_cols) {
  if (!col %in% colnames(meta)) stop(sprintf("Missing required column: %s", col))
}
l3_vals <- as.character(meta[[L3_SOURCE_COL]])
if (USE_EXISTING_L2) {
  l2_vals <- as.character(meta[[L2_SOURCE_COL]])
  if (all(is.na(l2_vals)) || all(trimws(l2_vals) == "")) {
    stop(sprintf("Existing L2 column '%s' is empty.", L2_SOURCE_COL))
  }
  cat(sprintf("[OK] Using existing L2 from '%s'\n", L2_SOURCE_COL))
} else {
  l2_vals <- unname(L3_TO_L2_REMAP[l3_vals])
  unmapped <- unique(l3_vals[is.na(l2_vals)])
  if (length(unmapped) > 0) {
    stop(sprintf(
      "Unmapped L3 values in '%s': %s\nUpdate L3_TO_L2_REMAP.",
      L3_SOURCE_COL,
      paste(unmapped, collapse = ", ")
    ))
  }
  original_l2_backup_col <- paste0(CELLTYPE_L2_COL, "_input")
  if (CELLTYPE_L2_COL %in% colnames(obj@meta.data) && !original_l2_backup_col %in% colnames(obj@meta.data)) {
    obj@meta.data[[original_l2_backup_col]] <- obj@meta.data[[CELLTYPE_L2_COL]]
    cat(sprintf("[OK] Preserved original '%s' as '%s'\n", CELLTYPE_L2_COL, original_l2_backup_col))
  }
}
obj@meta.data[[CELLTYPE_L2_COL]] <- l2_vals
obj@meta.data[[CELLTYPE_L3_COL]] <- l3_vals
l3_l2_mapping_summary <- as.data.frame(table(
  L3 = obj@meta.data[[CELLTYPE_L3_COL]],
  L2 = obj@meta.data[[CELLTYPE_L2_COL]],
  useNA = "no"
), stringsAsFactors = FALSE) %>%
  dplyr::filter(Freq > 0) %>%
  dplyr::arrange(L3, dplyr::desc(Freq), L2)

bad_idx <- is.na(meta[[TISSUE_COL]]) | trimws(as.character(meta[[TISSUE_COL]])) == "" |
  is.na(meta[[SAMPLE_COL]]) | trimws(as.character(meta[[SAMPLE_COL]])) == "" |
  is.na(obj@meta.data[[CELLTYPE_L2_COL]]) | is.na(obj@meta.data[[CELLTYPE_L3_COL]]) |
  trimws(as.character(obj@meta.data[[CELLTYPE_L3_COL]])) == ""
if (sum(bad_idx) > 0) {
  cat(sprintf("[INFO] Dropping %d cells with NA tissue/sample/L2/L3\n", sum(bad_idx)))
  obj <- subset(obj, cells = colnames(obj)[!bad_idx])
}
meta <- obj@meta.data
tissues  <- sort(unique(na.omit(meta[[TISSUE_COL]])))
l2_types <- sort(unique(na.omit(meta[[CELLTYPE_L2_COL]])))
l3_types <- sort(unique(na.omit(meta[[CELLTYPE_L3_COL]])))
cat(sprintf("[OK] Cells=%d | Tissues=%s\n", ncol(obj), paste(tissues, collapse = ", ")))
cat(sprintf("[OK] L2 types (%d): %s\n", length(l2_types), paste(l2_types, collapse = ", ")))
cat(sprintf("[OK] L3 types (%d): %s\n", length(l3_types), paste(l3_types, collapse = ", ")))
input_had_data_layer <- "data" %in% Layers(obj[["RNA"]])
obj <- NormalizeData(obj, verbose = FALSE)
obj@misc$normalization_for_report <- list(
  rerun = TRUE,
  method = "LogNormalize",
  input_had_data_layer = input_had_data_layer,
  purpose = c("visualization", "marker_analysis"),
  date = as.character(Sys.time())
)
cat("[OK] data layer ready\n\n")

cat("=== Loading Saved Summary Artifacts ===\n")
pb_de_all            <- resume_require_rds("reports/pseudobulk_de_all.rds")
pb_de_l3_all         <- resume_require_rds("reports/pseudobulk_de_L3_all.rds")
wilcox_all           <- resume_require_rds("reports/wilcox_exploratory_all.rds")
wilcox_l3_all        <- resume_require_rds("reports/wilcox_exploratory_L3_all.rds")
enrich_all           <- resume_require_rds("reports/enrichment_all.rds")
enrich_l3_all        <- resume_require_rds("reports/enrichment_L3_all.rds")
agent_all            <- resume_require_rds("reports/interpret_agent_all.rds")
agent_l3_all         <- resume_require_rds("reports/interpret_agent_L3_all.rds")
agent_structured_all <- resume_require_rds("reports/interpret_agent_structured_all.rds")
agent_structured_l3_all <- resume_require_rds("reports/interpret_agent_structured_L3_all.rds")
ssgsea_results_all   <- resume_require_rds("reports/ssgsea_results_all.rds")
ssgsea_results_l3_all <- resume_require_rds("reports/ssgsea_results_L3_all.rds")
ssgsea_llm_df_l2     <- resume_require_rds("reports/ssgsea_llm_structured_L2_all.rds")
ssgsea_llm_df_l3     <- resume_require_rds("reports/ssgsea_llm_structured_L3_all.rds")
ssgsea_llm_df        <- resume_require_rds("reports/ssgsea_llm_structured_all.rds")
choir_ofa_all        <- resume_require_rds("reports/choir/choir_ofa_all.rds")
choir_llm_df         <- resume_require_rds("reports/choir_llm_structured_all.rds")
agent_structured_df_l2 <- flatten_interpretation_records(agent_structured_all)
agent_structured_df_l3 <- flatten_interpretation_records(agent_structured_l3_all)
agent_structured_df    <- dplyr::bind_rows(agent_structured_df_l2, agent_structured_df_l3)

marker_panel_viz <- NULL
marker_panel_summary_path <- file.path(rpt_dir, "marker_panel_visualization_summary.csv")
if (file.exists(marker_panel_summary_path)) {
  marker_panel_viz <- list(summary = data.table::fread(marker_panel_summary_path, data.table = FALSE))
}

choir_results_all <- list(ssgsea = list(), ofa = choir_ofa_all)
for (mn in unique(c(SSGSEA_CHOIR_METHODS, SSGSEA_METHODS, "hallmark", "go_bp"))) {
  score_path <- file.path(rpt_dir, "choir", sprintf("ssgsea_scores_%s.rds", mn))
  z_path <- file.path(rpt_dir, "choir", sprintf("ssgsea_z_%s.rds", mn))
  if (file.exists(score_path) && file.exists(z_path)) {
    choir_results_all[["ssgsea"]][[mn]] <- list(
      scores = readRDS(score_path),
      z = readRDS(z_path)
    )
  }
}
choir_reduction <- pick_reduction(obj, CHOIR_REDUCTION_CANDIDATES)
cat("[OK] Saved artifacts loaded\n\n")

write_interpretation_markdown(
  agent_structured_df,
  file.path(OUTPUT_DIR, "LLM_INTERPRETATION.md"),
  "# LLM Interpretation Summary (Normal Tissue Comparison)",
  include_raw = FALSE
)
write_interpretation_markdown(
  agent_structured_df,
  file.path(rpt_dir, "LLM_INTERPRETATION_FOR_LLM.md"),
  "# LLM Interpretation Structured Input",
  include_raw = TRUE
)

cat("=== Resume Discovery-Screen LLM ===\n")
ssgsea_discovery_screen_df <- tc_empty_discovery_screen_df()
choir_discovery_screen_df  <- tc_empty_discovery_screen_df()
ofa_discovery_screen_df    <- tc_empty_discovery_screen_df()
if (ENABLE_LLM) {
  if (nrow(ssgsea_llm_df) > 0) {
    ssgsea_screen_input <- ssgsea_llm_df %>%
      dplyr::mutate(
        record_id = sprintf("ssgsea_%04d", dplyr::row_number()),
        record_label = paste(celltype_level, comparison, direction, sep = " | "),
        annotation_label = ifelse(nzchar(celltype_label), celltype_label, celltype_l2),
        primary_text = ifelse(nzchar(cell_type_judgment), cell_type_judgment, overview),
        supporting_text = paste(
          ifelse(nzchar(annotated_l3_correspondence), annotated_l3_correspondence, ""),
          ifelse(nzchar(discovery_assessment), discovery_assessment, ""),
          ifelse(nzchar(outlier_assessment), outlier_assessment, ""),
          ifelse(nzchar(integrated_diagnostic_comment), integrated_diagnostic_comment, evidence),
          sep = "\n"
        )
      ) %>%
      dplyr::select(record_id, record_label, annotation_label, primary_text, supporting_text)
    ssgsea_discovery_screen_df <- run_and_write_discovery_screen(
      records_df = ssgsea_screen_input,
      family_label = "ssGSEA grouped interpretations",
      output_dir = OUTPUT_DIR,
      prefix = "llm_ssgsea",
      title = "# LLM ssGSEA Discovery / Outlier Review",
      batch_size = 12L
    )
    cat(sprintf("[OK] ssGSEA discovery screen: %d rows\n", nrow(ssgsea_discovery_screen_df)))
  }

  if (nrow(choir_llm_df) > 0) {
    choir_screen_input <- choir_llm_df %>%
      dplyr::mutate(
        record_id = sprintf("choir_%04d", dplyr::row_number()),
        record_label = paste(celltype_level, comparison, direction, sep = " | "),
        annotation_label = ifelse(nzchar(annotated_l3_dominant), annotated_l3_dominant, celltype_l2),
        primary_text = ifelse(nzchar(cell_type_judgment), cell_type_judgment, overview),
        supporting_text = paste(
          ifelse(nzchar(annotated_l3_correspondence), annotated_l3_correspondence, ""),
          ifelse(nzchar(discovery_assessment), discovery_assessment, ""),
          ifelse(nzchar(outlier_assessment), outlier_assessment, ""),
          ifelse(nzchar(integrated_diagnostic_comment), integrated_diagnostic_comment, evidence),
          sep = "\n"
        )
      ) %>%
      dplyr::select(record_id, record_label, annotation_label, primary_text, supporting_text)
    choir_discovery_screen_df <- run_and_write_discovery_screen(
      records_df = choir_screen_input,
      family_label = "CHOIR cluster interpretations",
      output_dir = OUTPUT_DIR,
      prefix = "llm_choir",
      title = "# LLM CHOIR Discovery / Outlier Review",
      batch_size = 10L
    )
    cat(sprintf("[OK] CHOIR discovery screen: %d rows\n", nrow(choir_discovery_screen_df)))
  }

  if (length(choir_ofa_all) > 0) {
    ofa_screen_input <- build_choir_ofa_screening_table(choir_ofa_all)
    if (nrow(ofa_screen_input) > 0) {
      ofa_discovery_screen_df <- run_and_write_discovery_screen(
        records_df = ofa_screen_input,
        family_label = "CHOIR OFA one-vs-rest summaries",
        output_dir = OUTPUT_DIR,
        prefix = "llm_ofa",
        title = "# LLM OFA Discovery / Outlier Review",
        batch_size = 10L
      )
      cat(sprintf("[OK] OFA discovery screen: %d rows\n", nrow(ofa_discovery_screen_df)))
    }
  }
} else {
  stop("ENABLE_LLM is FALSE after initialization; resume finalizer cannot complete LLM discovery screens.")
}
cat("\n")

cat("=== Generating REPORT.md ===\n")
md <- character()
add <- function(...) md <<- c(md, paste0(...))

add(REPORT_TITLE); add("")
add("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M")); add("")
add("**Pipeline:** ", PIPELINE_SUBTITLE); add("")
add("**Note:** Cross-site anatomical comparison of NORMAL tissues, NOT disease vs healthy."); add("")
add("---"); add("")

add("## 1. Data Overview"); add("")
add(sprintf("- **Input:** `%s`", basename(H5AD_PATH)))
add(sprintf("- **Total cells:** %s", format(ncol(obj), big.mark = ",")))
add(sprintf("- **Tissues:** %s", paste(tissues, collapse = ", ")))
add(sprintf("- **L2 subtypes (%d):** %s", length(l2_types), paste(l2_types, collapse = ", ")))
add(sprintf("- **L3 subtypes (%d):** %s", length(l3_types), paste(l3_types, collapse = ", ")))
add(sprintf("- **L3 source column:** `%s` -> standardized `%s`", L3_SOURCE_COL, CELLTYPE_L3_COL))
if (isTRUE(USE_EXISTING_L2)) {
  add(sprintf("- **L2 source column:** `%s` (reused from input metadata)", L2_SOURCE_COL))
  add(""); add("### Observed L3 / L2 Mapping"); add("")
  add(sprintf("| %s | Existing `%s` | Cells |", L3_TO_L2_TABLE_HEADER_LEFT, L2_SOURCE_COL)); add("|---|---|---|")
  for (i in seq_len(nrow(l3_l2_mapping_summary))) {
    add(sprintf("| %s | %s | %d |",
      l3_l2_mapping_summary$L3[i],
      l3_l2_mapping_summary$L2[i],
      l3_l2_mapping_summary$Freq[i]
    ))
  }
} else {
  add(""); add("### L3 -> L2 Remapping"); add("")
  add(sprintf("| %s | L2 (merged) |", L3_TO_L2_TABLE_HEADER_LEFT)); add("|---|---|")
  for (i in seq_along(L3_TO_L2_REMAP)) add(sprintf("| %s | %s |", names(L3_TO_L2_REMAP)[i], L3_TO_L2_REMAP[i]))
}
add("")

add("## 2. Visualization"); add("")
add("![UMAP tissue](figures/umap_tissue.png)"); add("")
add("![UMAP L2](figures/umap_celltype_L2.png)"); add("")
add("![UMAP L3](figures/umap_celltype_L3.png)"); add("")
add("![UMAP split](figures/umap_L2_split_tissue.png)"); add("")
add("![Dotplot](figures/dotplot_markers.png)"); add("")
add("![Composition](figures/composition_tissue_L2.png)"); add("")
add("![Sample composition](figures/composition_sample_level.png)"); add("")
add("![Heatmap](figures/heatmap_top_markers.png)"); add("")
if (!is.null(marker_panel_viz) && is.data.frame(marker_panel_viz$summary) && nrow(marker_panel_viz$summary) > 0 &&
    exists("tc_marker_panel_report_lines", mode = "function")) {
  for (line in tc_marker_panel_report_lines(marker_panel_viz$summary, figure_dir_rel = sprintf("figures/%s", MARKER_PANEL_FIG_SUBDIR))) {
    add(line)
  }
}

add("## 3. Pseudobulk DESeq2 (Primary Inference)"); add("")
add(sprintf("padj < %s, |log2FC| > %s", PADJ_THR, LFC_THR)); add("")
add("### 3.1 L2"); add("")
add("| L2 Subtype | Comparison | Samples (ref/case) | Up | Down | Total |"); add("|---|---|---|---|---|---|")
for (ct in names(pb_de_all)) for (comp in names(pb_de_all[[ct]])) {
  r <- pb_de_all[[ct]][[comp]]; if (is.null(r)) next
  add(sprintf("| %s | %s | %d / %d | %d | %d | %d |", ct, comp, r$n_samples_1, r$n_samples_2, r$n_up, r$n_down, r$n_up + r$n_down))
}
add(""); add("### 3.2 L3"); add("")
add("| L3 Subtype | Comparison | Samples (ref/case) | Up | Down | Total |"); add("|---|---|---|---|---|---|")
for (ct in names(pb_de_l3_all)) for (comp in names(pb_de_l3_all[[ct]])) {
  r <- pb_de_l3_all[[ct]][[comp]]; if (is.null(r)) next
  add(sprintf("| %s | %s | %d / %d | %d | %d | %d |", ct, comp, r$n_samples_1, r$n_samples_2, r$n_up, r$n_down, r$n_up + r$n_down))
}
add("")

add("## 4. Multi-Database Enrichment"); add("")
add(sprintf("Databases: GO BP/MF/CC, KEGG, Hallmark, CellMarker, PanglaoDB, %s", CUSTOM_DB_LABEL)); add("")
for (ct in names(enrich_all)) for (comp in names(enrich_all[[ct]])) for (dir_name in names(enrich_all[[ct]][[comp]])) {
  enr_l <- enrich_all[[ct]][[comp]][[dir_name]]
  if (is.null(enr_l) || length(enr_l) == 0) next
  add(sprintf("### %s | %s | %s", ct, comp, dir_name)); add("")
  for (db in names(enr_l)) {
    er <- enr_l[[db]]
    er_df <- tryCatch(as.data.frame(er), error = function(e) data.frame())
    if (is.null(er) || nrow(er_df) == 0) next
    top5 <- utils::head(er_df, 5)
    add(sprintf("**%s (top 5):**", db)); add("")
    add("| Term | p.adjust | Count |"); add("|---|---|---|")
    for (j in seq_len(nrow(top5))) add(sprintf("| %s | %.2e | %s |", top5$Description[j], top5$p.adjust[j], top5$Count[j]))
    add("")
  }
}

add("## 5. Grouped ssGSEA (Average-Expression)"); add("")
if (!RUN_SSGSEA || (length(ssgsea_results_all) == 0 && length(ssgsea_results_l3_all) == 0)) {
  add("ssGSEA not run or no methods succeeded."); add("")
} else {
  ssgsea_methods_report <- unique(c(get_ssgsea_method_names(ssgsea_results_all), get_ssgsea_method_names(ssgsea_results_l3_all)))
  add(sprintf("Methods: %s", paste(ssgsea_methods_report, collapse = ", "))); add("")
  add("### 5.1 Tissue x L2"); add("")
  add("| Method | Pathways | Groups | Output |"); add("|---|---|---|---|")
  for (mn in get_ssgsea_method_names(ssgsea_results_all)) {
    add(sprintf("| %s | %d | %d | `reports/ssgsea_scores_%s.rds` |",
      mn, nrow(ssgsea_results_all[[mn]]$scores), ncol(ssgsea_results_all[[mn]]$scores), mn))
  }
  add("")
  for (mn in get_ssgsea_method_names(ssgsea_results_all)) {
    add(sprintf("![ssGSEA heatmap %s](figures/ssgsea_heatmap_%s.png)", mn, mn)); add("")
  }
  add("### 5.2 Tissue x L3"); add("")
  add("| Method | Pathways | Groups | Output |"); add("|---|---|---|---|")
  for (mn in get_ssgsea_method_names(ssgsea_results_l3_all)) {
    add(sprintf("| %s | %d | %d | `reports/ssgsea_l3_scores_%s.rds` |",
      mn, nrow(ssgsea_results_l3_all[[mn]]$scores), ncol(ssgsea_results_l3_all[[mn]]$scores), mn))
  }
  add("")
  for (mn in get_ssgsea_method_names(ssgsea_results_l3_all)) {
    add(sprintf("![L3 ssGSEA heatmap %s](figures/ssgsea_l3_heatmap_%s.png)", mn, mn)); add("")
  }
  if (exists("ssgsea_llm_df") && nrow(ssgsea_llm_df) > 0) {
    add(sprintf("### 5.3 ssGSEA LLM Summary (%d records)", nrow(ssgsea_llm_df))); add("")
    add("Structured outputs: `reports/ssgsea_llm_structured.tsv` | `LLM_SSGSEA_INTERPRETATION.md`")
    add("- LLM is run once per tissue x cell type group using combined ssGSEA evidence across multiple databases/methods, rather than once per database."); add("")
  }
}

add("## 6. CHOIR Clustering + Per-cluster ssGSEA + OFA"); add("")
if (!RUN_CHOIR || length(choir_ofa_all) == 0) {
  add("CHOIR not run or no clusters processed."); add("")
} else {
  choir_col_report <- paste0("CHOIR_clusters_", CHOIR_ALPHA)
  if (!choir_col_report %in% colnames(obj@meta.data)) {
    fb <- grep("^CHOIR_clusters", colnames(obj@meta.data), value = TRUE)
    choir_col_report <- if (length(fb) > 0) fb[1] else NULL
  }
  if (!is.null(choir_col_report) && choir_col_report %in% colnames(obj@meta.data)) {
    n_clusters_report <- length(unique(na.omit(obj@meta.data[[choir_col_report]])))
  } else {
    choir_cluster_table <- resume_require_table("reports/choir/choir_clusters.csv", sep = ",")
    mapped <- rep(NA_character_, ncol(obj))
    names(mapped) <- colnames(obj)
    overlap <- intersect(choir_cluster_table$cell, colnames(obj))
    mapped[overlap] <- as.character(choir_cluster_table$choir_cluster[match(overlap, choir_cluster_table$cell)])
    obj@meta.data[[paste0("CHOIR_clusters_", CHOIR_ALPHA)]] <- mapped[colnames(obj)]
    choir_col_report <- paste0("CHOIR_clusters_", CHOIR_ALPHA)
    n_clusters_report <- length(unique(na.omit(obj@meta.data[[choir_col_report]])))
  }
  add(sprintf("- **CHOIR alpha:** %.3f | **Clusters:** %s | **Reduction:** `%s`",
    CHOIR_ALPHA,
    ifelse(is.na(n_clusters_report), "unknown", n_clusters_report),
    if (!is.null(choir_reduction)) choir_reduction else "N/A"
  )); add("")
  add("### 6.1 CHOIR UMAP"); add("")
  add("![CHOIR UMAP](figures/choir_umap.png)"); add("")
  add("![CHOIR UMAP by Tissue](figures/choir_umap_split_tissue.png)"); add("")
  if (length(choir_results_all[["ssgsea"]]) > 0) {
    add("### 6.2 ssGSEA per CHOIR cluster"); add("")
    for (mn in names(choir_results_all[["ssgsea"]])) {
      add(sprintf("![CHOIR ssGSEA %s](figures/ssgsea_choir_heatmap_%s.png)", mn, mn)); add("")
    }
  }
  add("CHOIR cluster cell-count tables: `reports/choir/choir_cluster_L2_counts.csv` | `reports/choir/choir_cluster_L3_counts.csv` | `reports/choir/choir_cluster_L2_L3_counts_long.csv`")
  add("")
  add("### 6.3 OFA One-vs-Rest Markers"); add("")
  add("| Cluster | n_focal | n_rest | Sig markers | Integrated LLM |"); add("|---|---|---|---|---|")
  for (cl_name in names(choir_ofa_all)) {
    r <- choir_ofa_all[[cl_name]]
    n_sig <- if (!is.null(r$de_table)) sum(r$de_table$padj < OFA_PADJ_THR, na.rm = TRUE) else 0
    has_integrated <- if (!is.null(r$llm$integrated)) r$llm$integrated$status else "skipped"
    add(sprintf("| %s | %d | %d | %d | %s |", cl_name, r$n_focal, r$n_rest, n_sig, has_integrated))
  }
  add("")
  if (exists("choir_llm_df") && nrow(choir_llm_df) > 0) {
    add(sprintf("### 6.4 CHOIR Cluster LLM Summary (%d records)", nrow(choir_llm_df))); add("")
    add("Structured outputs: `reports/choir_llm_structured.tsv` | `LLM_CHOIR_INTERPRETATION.md`")
    add("- Language: Chinese (Simplified). key_drivers: English gene symbols.")
    add("- Each CHOIR record now includes user-annotated L3 correspondence, annotation match degree, outlier assessment, discovery assessment, and an integrated diagnostic comment based on annotation counts + complete DEG + ssGSEA + OFA."); add("")
  }
}

add("## 7. LLM Interpretation (interpret_agent)"); add("")
if (!ENABLE_LLM) {
  add("**SKIPPED:** DEEPSEEK_API_KEY not set. Re-run with API key to enable."); add("")
} else {
  add("### 7.1 Tissue-Pair DESeq2 LLM"); add("")
  if (nrow(agent_structured_df) == 0) {
    add("No interpret_agent results generated."); add("")
  } else {
    add("Structured outputs: `reports/interpret_agent_structured.tsv` | `LLM_INTERPRETATION.md`")
    add("- Language: Chinese (Simplified). Validated gene->pathway map used.")
    add("- Pairwise comparisons are interpreted once per comparison using integrated up/down evidence across multiple databases."); add("")
    for (i in seq_len(nrow(agent_structured_df))) {
      rec <- agent_structured_df[i, , drop = FALSE]
      add(sprintf("#### %s | %s | %s | %s", rec$celltype_level, rec$celltype_label, rec$comparison, rec$direction)); add("")
      add(sprintf("- **Status:** %s | **Source DB:** %s", rec$status, ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
      if (nzchar(rec$overview)) add(sprintf("- **Overview:** %s", rec$overview))
      if (nzchar(rec$hypothesis)) add(sprintf("- **Hypothesis:** %s", rec$hypothesis))
      add("")
    }
  }
  add("### 7.2 ssGSEA Group LLM (Tissue x Cell Type)"); add("")
  if (!exists("ssgsea_llm_df") || nrow(ssgsea_llm_df) == 0) {
    add("No ssGSEA LLM results generated."); add("")
  } else {
    add(sprintf("**%d records** across %d unique tissue x cell type groups using combined ssGSEA database evidence.",
      nrow(ssgsea_llm_df), length(unique(ssgsea_llm_df$comparison))))
    add("Full details: `LLM_SSGSEA_INTERPRETATION.md`"); add("")
    for (i in seq_len(nrow(ssgsea_llm_df))) {
      rec <- ssgsea_llm_df[i, , drop = FALSE]
      add(sprintf("#### %s | %s | %s", rec$celltype_level, rec$comparison, rec$direction)); add("")
      add(sprintf("- **Status:** %s | **DB:** %s", rec$status, ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
      if (nzchar(rec$overview)) add(sprintf("- **Overview:** %s", rec$overview))
      if (nzchar(rec$hypothesis)) add(sprintf("- **Hypothesis:** %s", rec$hypothesis))
      add("")
    }
  }
  add("### 7.3 CHOIR Cluster LLM (One-vs-Rest)"); add("")
  if (!exists("choir_llm_df") || nrow(choir_llm_df) == 0) {
    add(if (!RUN_CHOIR) "CHOIR was not run; cluster LLM skipped." else "No CHOIR cluster LLM results generated."); add("")
  } else {
    add(sprintf("**%d records** (%d CHOIR clusters; integrated up/down evidence per cluster).",
      nrow(choir_llm_df), length(unique(choir_llm_df$comparison))))
    add("Full details: `LLM_CHOIR_INTERPRETATION.md`"); add("")
    for (i in seq_len(nrow(choir_llm_df))) {
      rec <- choir_llm_df[i, , drop = FALSE]
      add(sprintf("#### %s | %s | %s", rec$celltype_level, rec$comparison, rec$direction)); add("")
      add(sprintf("- **Status:** %s | **DB:** %s", rec$status, ifelse(nzchar(rec$source_db), rec$source_db, "NA")))
      if (nzchar(rec$overview)) add(sprintf("- **Overview:** %s", rec$overview))
      if (nzchar(rec$key_mechanisms)) add(sprintf("- **Key Mechanisms:** %s", rec$key_mechanisms))
      add("")
    }
  }
}

add("## 8. Methods"); add("")
add("- **DE:** Pseudobulk DESeq2 (Squair et al. 2021 Nat Commun)")
add("- **Exploratory:** Cell-level Wilcoxon (marker discovery ONLY, NOT inference)")
add(sprintf("- **Enrichment:** clusterProfiler::enricher() + 8 databases (GO BP/MF/CC, KEGG, Hallmark, CellMarker, PanglaoDB, %s)", CUSTOM_DB_LABEL))
add(sprintf("- **ssGSEA:** GSVA::ssgseaParam on tissue x L2/L3 grouped average expression (not sample-level pseudobulk); methods: %s", paste(SSGSEA_METHODS, collapse = ", ")))
add("- **ssGSEA LLM [LLM-4]:** standardize_result_with_llm once per tissue x cell-type group using combined top pathways across multiple ssGSEA databases")
add("- **Pairwise comparison LLM:** integrated up/down DEG plus combined multi-database enrichment are provided together in one record per comparison")
add("- **CHOIR cluster LLM [LLM-5]:** standardize_result_with_llm once per CHOIR cluster using integrated one-vs-rest DEG, cluster-level ssGSEA, and multi-database enrichment")
if (ENABLE_LLM) {
  add("- **LLM:** DeepSeek-based structured interpretation with integrated evidence prompts")
  add("- **LLM language:** Chinese (Simplified); gene symbols retained in English")
  add("- **Gene emphasis:** top-|log2FC| genes explicitly linked to enriched pathways")
} else {
  add("- **LLM:** SKIPPED (no API key)")
}
add("")
add("### Output Objects"); add("")
add(sprintf("- `%s.rds` -- Seurat object (cell_type_L2 + cell_type_L3)", FINAL_FILE_PREFIX))
add(sprintf("- `%s.h5ad` -- AnnData object (cell_type_L2 + cell_type_L3)", FINAL_FILE_PREFIX))
add("")
add("---")
add(sprintf("*Generated by %s*", GENERATED_BY_LABEL))

writeLines(md, file.path(OUTPUT_DIR, "REPORT.md"))
cat("[OK] REPORT.md written\n")

cat("\n=== Saving Final Object (RDS + h5ad) ===\n")
choir_col_for_export <- paste0("CHOIR_clusters_", CHOIR_ALPHA)
if (!choir_col_for_export %in% colnames(obj@meta.data)) {
  choir_cluster_csv <- resume_require_table("reports/choir/choir_clusters.csv", sep = ",")
  mapped <- rep(NA_character_, ncol(obj))
  names(mapped) <- colnames(obj)
  overlap <- intersect(choir_cluster_csv$cell, colnames(obj))
  mapped[overlap] <- as.character(choir_cluster_csv$choir_cluster[match(overlap, choir_cluster_csv$cell)])
  obj@meta.data[[choir_col_for_export]] <- mapped[colnames(obj)]
  cat(sprintf("[OK] Restored %s from reports/choir/choir_clusters.csv\n", choir_col_for_export))
}
if (!"cell_type_L3" %in% colnames(obj@meta.data)) {
  obj@meta.data[["cell_type_L3"]] <- as.character(obj@meta.data[[L3_SOURCE_COL]])
  cat(sprintf("[OK] Created cell_type_L3 from '%s'\n", L3_SOURCE_COL))
}
stopifnot("cell_type_L2" %in% colnames(obj@meta.data))
stopifnot("cell_type_L3" %in% colnames(obj@meta.data))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L2"]])))
stopifnot(all(!is.na(obj@meta.data[["cell_type_L3"]])))
cat(sprintf("[OK] cell_type_L2: %d unique\n", length(unique(obj@meta.data[["cell_type_L2"]]))))
cat(sprintf("[OK] cell_type_L3: %d unique\n", length(unique(obj@meta.data[["cell_type_L3"]]))))

rds_path <- file.path(OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".rds"))
saveRDS(obj, rds_path)
cat(sprintf("[OK] RDS saved: %s (%.1f MB)\n", basename(rds_path), file.size(rds_path) / 1e6))

h5ad_out_path <- file.path(OUTPUT_DIR, paste0(FINAL_FILE_PREFIX, ".h5ad"))
tryCatch({
  anndata <- reticulate::import("anndata", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  builtins <- reticulate::import("builtins", convert = FALSE)
  counts_mat <- GetAssayData(obj, layer = "counts")
  norm_layer_name <- if ("data" %in% Layers(obj[["RNA"]])) "data" else "counts"
  expr_mat <- GetAssayData(obj, layer = norm_layer_name)
  counts_scipy <- matrix_to_scipy_csr(counts_mat, scipy_sparse, np)
  expr_scipy <- if (identical(norm_layer_name, "counts")) counts_scipy else matrix_to_scipy_csr(expr_mat, scipy_sparse, np)
  meta_export <- sanitize_obs_for_h5ad(obj@meta.data)
  obs_df <- reticulate::r_to_py(meta_export)
  var_df <- data.frame(gene_symbol = rownames(obj), row.names = rownames(obj), stringsAsFactors = FALSE)
  var_py <- reticulate::r_to_py(var_df)
  adata <- anndata$AnnData(X = expr_scipy, obs = obs_df, var = var_py)
  adata$layers$`__setitem__`("counts", counts_scipy)
  adata$uns$`__setitem__`("X_layer", norm_layer_name)
  for (red_name in Reductions(obj)) {
    emb <- Embeddings(obj, reduction = red_name)
    adata$obsm$`__setitem__`(paste0("X_", red_name), np$array(emb, dtype = np$float32))
  }
  adata$write_h5ad(h5ad_out_path, compression = "gzip")
  cat(sprintf("[OK] h5ad saved: %s (%.1f MB)\n", basename(h5ad_out_path), file.size(h5ad_out_path) / 1e6))
  adata_check <- anndata$read_h5ad(h5ad_out_path)
  obs_cols <- reticulate::py_to_r(builtins$list(adata_check$obs$columns))
  layer_keys <- reticulate::py_to_r(builtins$list(adata_check$layers$keys()))
  stopifnot("cell_type_L2" %in% obs_cols, "cell_type_L3" %in% obs_cols, "counts" %in% layer_keys)
  cat(sprintf("[OK] h5ad verified: %d cells x %d genes, cell_type_L2 + L3 + counts layer present\n",
    reticulate::py_to_r(adata_check$n_obs), reticulate::py_to_r(adata_check$n_vars)))
  rm(adata, adata_check, builtins, counts_mat, expr_mat, counts_scipy, expr_scipy, meta_export, var_df)
  gc()
}, error = function(e) {
  cat(sprintf("[ERROR] h5ad export failed: %s\n", e$message))
  cat("[INFO] RDS saved successfully.\n")
})

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("MYELOID TISSUE COMPARISON RESUME FINALIZER COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")
cat(sprintf("Output: %s\n\n", OUTPUT_DIR))
cat("Key outputs:\n")
cat("  REPORT.md\n")
cat("  LLM_INTERPRETATION.md\n")
cat("  LLM_SSGSEA_INTERPRETATION.md\n")
cat("  LLM_CHOIR_INTERPRETATION.md\n")
cat("  LLM_SSGSEA_DISCOVERY_REVIEW.md\n")
cat("  LLM_CHOIR_DISCOVERY_REVIEW.md\n")
cat("  LLM_OFA_DISCOVERY_REVIEW.md\n")
cat(sprintf("  %s.rds / .h5ad\n", FINAL_FILE_PREFIX))
cat("\n")
cat(sprintf("%s\n", PIPELINE_CHANGELOG_TITLE))
for (line in PIPELINE_CHANGELOG_LINES) cat(sprintf("%s\n", line))
cat("\n")
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
cat("[OK] Done\n")
