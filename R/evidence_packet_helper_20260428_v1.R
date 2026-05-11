#!/usr/bin/env Rscript
# ==============================================================================
# Evidence Packet Helper (2026-04-28 v1)
# ==============================================================================

PA_CORE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_architecture_core_20260428_v1.R"
PA_AGGREGATION_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/aggregation_registry_helper_20260428_v1.R"
PA_VALIDATION_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_validation_helper_20260428_v1.R"
PA_CROSSWALK_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/program_crosswalk_helper_20260428_v1.R"
PA_TRAJECTORY_BRANCH_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/trajectory_branch_helper_20260428_v1.R"
PA_SYNTHESIS_BRIDGE_HELPER_PATH_20260428_V1 <- "/home/h2048/script/R/synthesis_bridge_helper_20260428_v1.R"
if (!exists("pa_new_analysis_unit", mode = "function")) source(PA_CORE_HELPER_PATH_20260428_V1)
if (!exists("pa_register_aggregation", mode = "function")) source(PA_AGGREGATION_HELPER_PATH_20260428_V1)
if (!exists("pa_build_program_validation_packet", mode = "function")) source(PA_VALIDATION_HELPER_PATH_20260428_V1)
if (!exists("pa_build_deg_program_crosswalk", mode = "function")) source(PA_CROSSWALK_HELPER_PATH_20260428_V1)
if (!exists("pa_build_trajectory_branch_packet", mode = "function")) source(PA_TRAJECTORY_BRANCH_HELPER_PATH_20260428_V1)
if (!exists("pa_build_synthesis_packet", mode = "function")) source(PA_SYNTHESIS_BRIDGE_HELPER_PATH_20260428_V1)

pa_build_evidence_packet <- function(
  unit_id,
  aggregation_context,
  de_table = NULL,
  enrichment_list = NULL,
  program_crosswalk = NULL,
  score_summary = NULL,
  network_summary = NULL,
  validation_summary = NULL,
  conflict_summary = NULL,
  context_str = NULL,
  trajectory_branch_packet = NULL,
  enrichment_packet = NULL,
  interpret_agent_packet = NULL,
  report_packet = NULL,
  synthesis_packet = NULL
) {
  list(
    unit_id = pa_scalar_chr(unit_id, "unit_id"),
    aggregation_context = aggregation_context,
    de_table = de_table,
    enrichment_list = enrichment_list,
    program_crosswalk = program_crosswalk,
    score_summary = score_summary,
    network_summary = network_summary,
    validation_summary = validation_summary,
    conflict_summary = conflict_summary,
    context_str = context_str,
    trajectory_branch_packet = trajectory_branch_packet,
    enrichment_packet = enrichment_packet,
    interpret_agent_packet = interpret_agent_packet,
    report_packet = report_packet,
    synthesis_packet = synthesis_packet,
    helper_version = PA_HELPER_VERSION_20260428_V1
  )
}

pa_run_program_architecture_unit <- function(
  unit,
  discrete_outputs,
  program_inputs,
  rewiring_inputs,
  outdir
) {
  pa_validate_required_columns(unit, c("unit_id", "lineage", "aggregation_type"), "unit")
  outdir <- pa_prepare_output_dir(outdir)
  unit_id <- unit$unit_id[[1]]

  discrete_outputs <- pa_null_coalesce(discrete_outputs, list())
  program_inputs <- pa_null_coalesce(program_inputs, list())
  rewiring_inputs <- pa_null_coalesce(rewiring_inputs, list())

  aggregation_registry <- pa_null_coalesce(
    program_inputs$aggregation_registry,
    pa_register_aggregation(
      unit_id = unit_id,
      aggregation_type = unit$aggregation_type[[1]],
      source_assay = "RNA",
      source_layer = unit$matrix_type[[1]],
      lineage = unit$lineage[[1]],
      group_keys = c("sample", unit$state_level[[1]])
    )
  )

  program_registry <- pa_null_coalesce(program_inputs$program_registry, data.frame())
  trajectory_branch_packet <- pa_null_coalesce(program_inputs$trajectory_branch_packet, rewiring_inputs$trajectory_branch_packet)
  enrichment_packet <- pa_null_coalesce(discrete_outputs$enrichment_packet, program_inputs$enrichment_packet)
  interpret_agent_packet <- pa_null_coalesce(program_inputs$interpret_agent_packet, rewiring_inputs$interpret_agent_packet)
  report_packet <- pa_null_coalesce(program_inputs$report_packet, rewiring_inputs$report_packet)
  synthesis_packet <- pa_null_coalesce(
    program_inputs$synthesis_packet,
    if (!is.null(enrichment_packet) || !is.null(interpret_agent_packet) || !is.null(report_packet)) {
      pa_build_synthesis_packet(
        enrichment_packet = enrichment_packet,
        interpret_agent_packet = interpret_agent_packet,
        report_packet = report_packet
      )
    } else {
      NULL
    }
  )
  validation_packet <- pa_null_coalesce(
    program_inputs$validation_packet,
    pa_build_program_validation_packet(
      dme_table = program_inputs$dme_table,
      trait_correlation = program_inputs$trait_correlation,
      preservation_stats = program_inputs$preservation_stats,
      projection_stats = program_inputs$projection_stats,
      robustness_summary = program_inputs$robustness_summary
    )
  )

  crosswalk <- pa_null_coalesce(
    program_inputs$program_crosswalk,
    if (!is.null(discrete_outputs$deg_table) && is.data.frame(program_registry) && nrow(program_registry) > 0L) {
      pa_build_deg_program_crosswalk(discrete_outputs$deg_table, program_registry)
    } else {
      pa_empty_crosswalk()
    }
  )

  evidence_packet <- pa_build_evidence_packet(
    unit_id = unit_id,
    aggregation_context = aggregation_registry,
    de_table = discrete_outputs$deg_table,
    enrichment_list = discrete_outputs$enrichment_list,
    program_crosswalk = crosswalk,
    score_summary = program_inputs$score_summary,
    network_summary = rewiring_inputs$network_summary,
    validation_summary = validation_packet,
    conflict_summary = rewiring_inputs$conflict_summary,
    context_str = rewiring_inputs$context_str,
    trajectory_branch_packet = trajectory_branch_packet,
    enrichment_packet = enrichment_packet,
    interpret_agent_packet = interpret_agent_packet,
    report_packet = report_packet,
    synthesis_packet = synthesis_packet
  )

  pa_write_rds(unit, file.path(outdir, "analysis_unit.rds"))
  pa_write_tsv(unit, file.path(outdir, "analysis_unit.tsv"))
  if (is.data.frame(aggregation_registry) && nrow(aggregation_registry) > 0L) {
    pa_write_rds(aggregation_registry, file.path(outdir, "aggregation_registry.rds"))
    pa_write_tsv(aggregation_registry, file.path(outdir, "aggregation_registry.tsv"))
  }
  if (is.data.frame(program_registry) && nrow(program_registry) > 0L) {
    pa_write_rds(program_registry, file.path(outdir, "program_registry.rds"))
    pa_write_tsv(program_registry, file.path(outdir, "program_registry.tsv"))
  }
  if (is.data.frame(crosswalk) && nrow(crosswalk) > 0L) {
    pa_write_rds(crosswalk, file.path(outdir, "program_crosswalk.rds"))
    pa_write_tsv(crosswalk, file.path(outdir, "program_crosswalk.tsv"))
  }
  pa_write_rds(validation_packet, file.path(outdir, "program_validation_packet.rds"))
  if (is.list(trajectory_branch_packet)) {
    pa_write_rds(trajectory_branch_packet, file.path(outdir, "trajectory_branch_packet.rds"))
  }
  if (is.list(synthesis_packet)) {
    pa_write_rds(synthesis_packet, file.path(outdir, "synthesis_packet.rds"))
  }
  pa_write_rds(evidence_packet, file.path(outdir, "evidence_packet.rds"))
  pa_write_json(
    list(
      unit_id = unit_id,
      contrast_id = unit$contrast_id[[1]],
      lineage = unit$lineage[[1]],
      aggregation_type = unit$aggregation_type[[1]],
      helper_version = PA_HELPER_VERSION_20260428_V1
    ),
    file.path(outdir, "program_architecture_manifest.json")
  )

  invisible(list(
    unit = unit,
    aggregation_registry = aggregation_registry,
    program_registry = program_registry,
    validation_packet = validation_packet,
    crosswalk = crosswalk,
    evidence_packet = evidence_packet
  ))
}

if (sys.nframe() == 0) {
  cat("Evidence Packet Helper (2026-04-28 v1) loaded.\n")
}