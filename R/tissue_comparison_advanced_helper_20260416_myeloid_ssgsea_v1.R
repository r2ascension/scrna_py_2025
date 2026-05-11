#!/usr/bin/env Rscript
# ==============================================================================
# Tissue Comparison Advanced Helper Extensions (2026-04-16 Myeloid ssGSEA fix)
# ==============================================================================
#
# Purpose:
#   - keep older helper versions immutable
#   - preserve myeloid marker-panel defaults from 2026-04-14 v2
#   - make ssGSEA-only review and discovery screening more conservative so that
#     non-canonical pathway profiles are downgraded to mixed/uncertain unless
#     there is explicit evidence for contamination / misannotation / artifact
#
# Date: 2026-04-16
# ==============================================================================

BASE_ADVANCED_HELPER_PATH_20260416 <- "/home/h2048/script/R/tissue_comparison_advanced_helper_20260414_v2.R"
if (!file.exists(BASE_ADVANCED_HELPER_PATH_20260416)) {
  stop(sprintf("Base advanced helper not found: %s", BASE_ADVANCED_HELPER_PATH_20260416))
}
source(BASE_ADVANCED_HELPER_PATH_20260416)

tc_run_ssgsea_group_review_llm <- function(group_id,
                                           ss_ctx,
                                           evidence_text,
                                           source_db_label,
                                           annotated_label,
                                           annotated_level,
                                           warnings = character(),
                                           error_message = NULL,
                                           caller_env = parent.frame()) {
  enable_llm <- isTRUE(tc_lookup_in_caller("ENABLE_LLM", FALSE, caller_env))
  if (!enable_llm) {
    return(list(result = NULL, warnings = character(), error = "LLM disabled"))
  }

  lineage_context_lower <- tc_lookup_in_caller("LINEAGE_CONTEXT_LOWER", "tissue", caller_env)
  standardize_model <- tc_lookup_in_caller("STANDARDIZE_LLM_MODEL", tc_deepseek_default_chat_model(), caller_env)
  deepseek_api_key <- tc_lookup_in_caller("DEEPSEEK_API_KEY", Sys.getenv("DEEPSEEK_API_KEY", unset = ""), caller_env)
  standardize_retries <- as.integer(tc_lookup_in_caller("STANDARDIZE_LLM_MAX_RETRIES", 3L, caller_env))
  standardize_retry_sleep <- tc_lookup_in_caller("STANDARDIZE_LLM_RETRY_SLEEP_SEC", 2, caller_env)
  max_input_chars <- max(as.integer(tc_lookup_in_caller("STANDARDIZE_LLM_MAX_INPUT_CHARS", 12000L, caller_env)), 18000L)
  parse_json <- tc_lookup_in_caller("parse_standardized_json", tc_parse_json, caller_env)

  prompt <- paste(
    sprintf("You are reviewing a grouped ssGSEA profile in %s tissue comparison.", lineage_context_lower),
    "Return valid JSON only. No markdown, no code fences, no commentary.",
    "Use exactly these keys:",
    "cell_type_judgment, confidence, annotation_match_degree, annotated_l3_correspondence, outlier_assessment, discovery_assessment, integrated_diagnostic_comment, overview, key_mechanisms, hypothesis, narrative, key_drivers, evidence, limitations.",
    "",
    "Rules:",
    "1. All narrative fields must be Simplified Chinese strings; key_drivers must be an array of English gene symbols or marker names if needed.",
    "2. confidence must be one of: high, medium, low.",
    "3. annotation_match_degree must be one of: high, moderate, mixed, low.",
    "4. annotated_l3_correspondence must explicitly compare the ssGSEA-inferred state with the provided user annotation label.",
    "5. outlier_assessment must only call likely misannotation, severe contamination, or technical artifact when the pathways provide direct evidence for that conclusion (for example: coherent alternative-lineage programs, repeated mixed-lineage signals, or strong artifact/stress dominance). If the main issue is simply that canonical lineage pathways are weak, absent, or non-specific, explicitly say that current evidence is insufficient to call contamination.",
    "6. discovery_assessment must state whether the ssGSEA pattern suggests a plausible biological substate / tissue-adapted program / activation state worth follow-up; if not, say so explicitly.",
    "7. integrated_diagnostic_comment must integrate positive and negative ssGSEA signals together and explain whether this group is annotation-aligned, mixed/transition-like, likely outlier, or plausible biological discovery.",
    "8. Use ssGSEA as the main evidence; do not overstate mechanism beyond the pathways provided.",
    "9. overview must start with a concise cell-type/state judgment.",
    "10. hypothesis should be brief; use 'Not applicable.' if no strong hypothesis is justified.",
    sprintf("11. %s", tc_noninformative_gene_rule_text()),
    "12. Be conservative: broad stress/metabolic signatures alone do not automatically imply discovery.",
    "13. Lack of expected macrophage / monocyte / dendritic pathways alone is NOT enough to conclude contamination. For ssGSEA-only evidence, prefer mixed_or_uncertain style language unless there is explicit alternative-lineage or artifact evidence.",
    "",
    sprintf("Group: %s", group_id),
    sprintf("Annotated label: %s (%s)", annotated_label, annotated_level),
    sprintf("Source DB summary: %s", source_db_label),
    "Context:",
    tc_truncate_text(ss_ctx, max_chars = max_input_chars),
    "Evidence:",
    tc_truncate_text(evidence_text, max_chars = max_input_chars),
    "Warnings and errors:",
    tc_truncate_text(paste(c(warnings, error_message), collapse = "\n"), max_chars = 3000L),
    sep = "\n"
  )

  collected_warnings <- character()
  last_error <- NULL
  for (attempt in seq_len(max(1L, standardize_retries))) {
    if (attempt > 1) cat(sprintf("    [INFO] ssGSEA review LLM retry %d/%d\n", attempt, standardize_retries))
    response_text <- tryCatch(
      tc_safe_trim(tc_deepseek_chat_request(prompt, model = standardize_model, api_key = deepseek_api_key)),
      error = function(e) {
        last_error <<- conditionMessage(e)
        ""
      }
    )
    parsed <- parse_json(response_text)
    if (!is.null(parsed) && is.list(parsed)) {
      return(list(result = parsed, warnings = unique(collected_warnings), error = NULL))
    }
    if (nzchar(response_text)) {
      collected_warnings <- c(collected_warnings, sprintf("ssgsea review attempt %d returned non-JSON text", attempt))
    }
    if (attempt < standardize_retries) Sys.sleep(standardize_retry_sleep)
  }

  list(
    result = NULL,
    warnings = unique(collected_warnings),
    error = ifelse(is.null(last_error) || !nzchar(last_error), "ssGSEA review LLM failed", last_error)
  )
}

tc_build_record_screening_prompt <- function(batch_df, family_label) {
  required_cols <- c("record_id", "record_label", "annotation_label", "primary_text", "supporting_text")
  for (nm in setdiff(required_cols, names(batch_df))) batch_df[[nm]] <- ""
  payload <- jsonlite::toJSON(
    batch_df[, required_cols, drop = FALSE],
    dataframe = "rows",
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )

  is_ssgsea_family <- grepl("ssgsea", family_label %||% "", ignore.case = TRUE)
  discovery_rule <- if (is_ssgsea_family) {
    "5. potential_discovery means a coherent lineage-consistent state/substate or tissue-adapted program that is not better explained by contamination or obvious misannotation."
  } else {
    "5. potential_discovery means a coherent state/substate or tissue-adapted program that is not better explained by contamination or obvious misannotation."
  }
  outlier_rule <- if (is_ssgsea_family) {
    "6. likely_outlier means the record contains explicit, coherent evidence for misannotation, contamination, doublet, or technical/stress artifact that is stronger than any lineage-consistent interpretation. Lack of canonical pathways, weak pathway support, or broad non-specific programs alone is NOT sufficient."
  } else {
    "6. likely_outlier means evidence favors misannotation, contamination, doublet, or technical/stress artifact over a coherent lineage-consistent state."
  }
  mixed_rule <- if (is_ssgsea_family) {
    "8. mixed_or_uncertain means evidence is insufficient, indirect, or conflicting. For ssGSEA-only records, if the pattern is non-canonical but contamination is not directly supported, prefer mixed_or_uncertain over likely_outlier."
  } else {
    "8. mixed_or_uncertain means evidence is insufficient or conflicting."
  }
  extra_rule <- if (is_ssgsea_family) {
    "10. For ssGSEA-only screening, do not escalate to likely_outlier unless the text explicitly supports contamination / misannotation / artifact; otherwise keep the call conservative."
  } else {
    NULL
  }

  paste(
    sprintf("You are screening %s records to identify likely biological discoveries and likely outliers.", family_label),
    "Return valid JSON only. No markdown, no commentary.",
    "Return a JSON array with the same number of items and the same record_id values as the input.",
    "Each item must contain exactly these keys:",
    "record_id, biological_signal_class, confidence, short_call, discovery_flag, outlier_flag, evidence_summary, followup.",
    "Rules:",
    "1. biological_signal_class must be one of: aligned_state, potential_discovery, likely_outlier, mixed_or_uncertain.",
    "2. confidence must be one of: high, medium, low.",
    "3. discovery_flag and outlier_flag must be yes or no.",
    "4. short_call, evidence_summary, and followup must be concise Simplified Chinese.",
    discovery_rule,
    outlier_rule,
    "7. aligned_state means mostly matches annotation with no strong discovery or outlier signal.",
    mixed_rule,
    sprintf("9. %s", tc_noninformative_gene_rule_text()),
    extra_rule,
    "Input records:",
    payload,
    sep = "\n"
  )
}

if (sys.nframe() == 0) {
  cat("Tissue Comparison Advanced Helper extensions (2026-04-16 myeloid ssGSEA fix) loaded.\n")
  cat(sprintf("Base helper: %s\n", BASE_ADVANCED_HELPER_PATH_20260416))
}
