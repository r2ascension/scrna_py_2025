#!/usr/bin/env Rscript
# ==============================================================================
# Reviewer Bundle Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
PA_VALIDATION_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_validation_helper_20260428_v1.R"
PA_TRAJECTORY_BRANCH_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/trajectory_branch_helper_20260428_v1.R"
PA_SYNTHESIS_BRIDGE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/synthesis_bridge_helper_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) source(PA_CORE_HELPER_PATH_20260428_V1)
if (!exists("pa_validation_packet_summary_lines", mode = "function")) source(PA_VALIDATION_HELPER_PATH_20260428_V1)
if (!exists("pa_trajectory_packet_summary_lines", mode = "function")) source(PA_TRAJECTORY_BRANCH_HELPER_PATH_20260428_V1)
if (!exists("pa_synthesis_packet_summary_lines", mode = "function")) source(PA_SYNTHESIS_BRIDGE_HELPER_PATH_20260428_V1)

pa_register_figure <- function(
  figure_id,
  panel_id,
  source_unit,
  input_tables,
  claim_sentence,
  support_packet,
  caveat_sentence = NA_character_,
  output_path
) {
  data.frame(
    figure_id = pa_scalar_chr(figure_id, "figure_id"),
    panel_id = pa_scalar_chr(panel_id, "panel_id"),
    source_unit = pa_scalar_chr(source_unit, "source_unit"),
    input_tables = I(list(pa_unique_chr(input_tables))),
    claim_sentence = pa_scalar_chr(claim_sentence, "claim_sentence"),
    support_packet = I(list(support_packet)),
    caveat_sentence = pa_null_coalesce(caveat_sentence, NA_character_),
    output_path = pa_scalar_chr(output_path, "output_path"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

pa_write_figure_registry <- function(figure_registry, output_path) {
  pa_validate_required_columns(figure_registry, c("figure_id", "panel_id", "source_unit", "claim_sentence"), "figure_registry")
  pa_write_tsv(figure_registry, output_path)
}

pa_build_methods_snippet <- function(unit, program_registry = NULL, validation_packet = NULL, evidence_packet = NULL) {
  pa_validate_required_columns(unit, c("lineage", "state_level", "contrast_id", "aggregation_type"), "unit")
  program_sources <- if (is.data.frame(program_registry) && nrow(program_registry) > 0L) unique(program_registry$source_type) else unit$program_sources[[1]]
  c(
    "## Methods snippet",
    "",
    sprintf(
      "我们在 `%s` lineage 的 `%s` 层级上，对 `%s` 比较构建了统一的四层分析对象；基础聚合粒度设为 `%s`。",
      unit$lineage[[1]], unit$state_level[[1]], unit$contrast_id[[1]], unit$aggregation_type[[1]]
    ),
    sprintf(
      "程序层来源统一注册为：%s；所有结果通过 analysis unit、aggregation registry、program registry、validation packet 与 evidence packet 串联。",
      paste(program_sources, collapse = ", ")
    ),
    "验证轨输出统一记录 DME、trait correlation、preservation 与 robustness 等统计摘要，并为后续 enrichment、interpret_agent 与 reviewer appendix 提供固定上游 contract。",
    if (isTRUE(unit$trajectory_enabled[[1]])) {
      "轨迹层默认采用 PAGA 作为 topology screen、Slingshot 作为主 trajectory engine，并用 CytoTRACE2 作为方向/成熟度验证；Monocle3 作为可选兼容引擎保留。"
    } else {
      "轨迹层在该 unit 中未启用。"
    },
    "",
    "### Validation summary",
    pa_validation_packet_summary_lines(validation_packet),
    "",
    "### Trajectory / synthesis summary",
    pa_trajectory_packet_summary_lines(if (is.list(evidence_packet)) evidence_packet$trajectory_branch_packet else NULL),
    pa_synthesis_packet_summary_lines(if (is.list(evidence_packet)) evidence_packet$synthesis_packet else NULL)
  )
}

pa_build_result_snippet <- function(unit, evidence_packet, figure_registry = NULL, validation_packet = NULL) {
  pa_validate_required_columns(unit, c("contrast_id", "lineage", "state_level"), "unit")
  crosswalk_n <- if (!is.null(evidence_packet$program_crosswalk) && is.data.frame(evidence_packet$program_crosswalk)) nrow(evidence_packet$program_crosswalk) else 0L
  figure_n <- if (is.data.frame(figure_registry)) nrow(figure_registry) else 0L
  c(
    "## Results snippet",
    "",
    sprintf(
      "在 `%s` lineage 的 `%s` 层级比较 `%s` 时，我们将离散层差异结果与程序层 program registry 统一映射，得到 %d 条 crosswalk 记录。",
      unit$lineage[[1]], unit$state_level[[1]], unit$contrast_id[[1]], crosswalk_n
    ),
    sprintf(
      "这些结果随后被打包进 evidence packet，并可直接驱动 enrichment 汇总、LLM 解释以及 reviewer 输出；当前 figure registry 中登记了 %d 个图面板。",
      figure_n
    ),
    "如果 validation packet 提供了 preservation 或 trait-correlation 证据，则这些统计会与离散层 DE/DA 结果一起被保留，避免 reviewer 阶段重新手工拼接。",
    "",
    "### Validation summary",
    pa_validation_packet_summary_lines(validation_packet),
    "",
    "### Trajectory / synthesis summary",
    pa_trajectory_packet_summary_lines(evidence_packet$trajectory_branch_packet),
    pa_synthesis_packet_summary_lines(evidence_packet$synthesis_packet)
  )
}

pa_build_figure_legend <- function(figure_registry) {
  if (!is.data.frame(figure_registry) || nrow(figure_registry) == 0L) {
    return(c("## Figure legend", "", "当前没有登记 figure registry。"))
  }
  lines <- c("## Figure legend", "")
  for (i in seq_len(nrow(figure_registry))) {
    lines <- c(
      lines,
      sprintf("- %s%s：%s", figure_registry$figure_id[[i]], figure_registry$panel_id[[i]], figure_registry$claim_sentence[[i]])
    )
  }
  lines
}

pa_build_sensitivity_appendix <- function(unit, validation_packet = NULL) {
  pa_validate_required_columns(unit, c("contrast_id", "aggregation_type"), "unit")
  c(
    "## Sensitivity appendix",
    "",
    sprintf("- Analysis unit: `%s`", unit$unit_id[[1]]),
    sprintf("- Contrast: `%s`", unit$contrast_id[[1]]),
    sprintf("- Aggregation type: `%s`", unit$aggregation_type[[1]]),
    "- 建议在 reviewer response 中引用同名 validation packet，以说明程序层结果并非仅来自单一打分来源。",
    pa_validation_packet_summary_lines(validation_packet)
  )
}

pa_write_reviewer_bundle <- function(
  output_dir,
  unit,
  evidence_packet,
  figure_registry,
  program_registry = NULL,
  validation_packet = NULL
) {
  output_dir <- pa_prepare_output_dir(output_dir)
  methods_path <- file.path(output_dir, "methods_snippet.md")
  result_path <- file.path(output_dir, "result_snippet.md")
  legend_path <- file.path(output_dir, "figure_legend.md")
  appendix_path <- file.path(output_dir, "sensitivity_appendix.md")
  figure_tsv_path <- file.path(output_dir, "figure_registry.tsv")

  pa_write_markdown(pa_build_methods_snippet(unit, program_registry, validation_packet, evidence_packet), methods_path)
  pa_write_markdown(pa_build_result_snippet(unit, evidence_packet, figure_registry, validation_packet), result_path)
  pa_write_markdown(pa_build_figure_legend(figure_registry), legend_path)
  pa_write_markdown(pa_build_sensitivity_appendix(unit, validation_packet), appendix_path)
  if (is.data.frame(figure_registry) && nrow(figure_registry) > 0L) {
    pa_write_figure_registry(figure_registry, figure_tsv_path)
  } else {
    pa_write_tsv(data.frame(), figure_tsv_path)
  }

  invisible(list(
    methods_snippet = methods_path,
    result_snippet = result_path,
    figure_legend = legend_path,
    sensitivity_appendix = appendix_path,
    figure_registry_tsv = figure_tsv_path
  ))
}

if (sys.nframe() == 0) {
  cat("Reviewer Bundle Helper (2026-04-28 v1) loaded.\n")
}