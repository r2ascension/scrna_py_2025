#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 2026-05-08 modified-lineage full reruns
# ------------------------------------------------------------------------------
# Flow:
#   1) Prepare filtered upstream h5ad inputs for epithelial/endothelial.
#   2) Retrain scVI/scANVI h5ad references for B, epithelial, endothelial, TNK.
#      Python training is intentionally sequential by default to avoid GPU OOM.
#   3) Launch R tissue-comparison wrappers in parallel; each wrapper uses the
#      20260508 batch-parallel LLM overlay so one lineage's LLM interpretation can
#      run while other lineages are still in downstream analysis.
#   4) Add pathway enrichment barplots after each lineage completes.
# ==============================================================================

ROOT="/home/h2048"
PY="${ROOT}/miniconda3/envs/scvi_env/bin/python"
RSCRIPT="/usr/bin/Rscript"
LOG_ROOT="${ROOT}/logs/20260508/modified_lineage_full_rerun"
mkdir -p "${LOG_ROOT}"

export PYTHONNOUSERSITE=1
export LLM_SCREEN_PARALLEL_WORKERS="${LLM_SCREEN_PARALLEL_WORKERS:-3}"
export LLM_SCREEN_PARALLEL_STAGGER_SEC="${LLM_SCREEN_PARALLEL_STAGGER_SEC:-0.5}"
export STANDARDIZE_LLM_RETRY_SLEEP_SEC="${STANDARDIZE_LLM_RETRY_SLEEP_SEC:-1}"
SCVI_ENV_LIB="${ROOT}/miniconda3/envs/scvi_env/lib"
SCVI_ENV_CUDA_LIB="${ROOT}/miniconda3/envs/scvi_env/targets/x86_64-linux/lib"
R_LD_LIBRARY_PATH="${SCVI_ENV_LIB}:${SCVI_ENV_CUDA_LIB}"

run_logged() {
  local name="$1"; shift
  local log_file="${LOG_ROOT}/${name}.log"
  echo "[$(date '+%F %T')] START ${name}" | tee "${log_file}"
  env -u LD_LIBRARY_PATH -u PYTHONPATH "$@" 2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} status=${status}" | tee -a "${log_file}"
  return "${status}"
}

filter_h5ad_inputs() {
  local epi_out="${ROOT}/data/py/0508/epithelial_input_rm_leiden14_17_20260508/epithelial_with_subclusters_rm_leiden14_17_20260508.h5ad"
  local epi_summary="${ROOT}/data/py/0508/epithelial_input_rm_leiden14_17_20260508/filter_summary.json"
  if [[ -s "${epi_out}" && -s "${epi_summary}" ]]; then
    echo "[SKIP] epithelial_filter_input already exists: ${epi_out}"
  else
    run_logged epithelial_filter_input \
      "${PY}" "${ROOT}/script/py/filter_h5ad_by_removed_cells_20260508.py" \
        --input-h5ad "${ROOT}/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad" \
        --removed-cells "${ROOT}/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/posthoc_cluster_removal_20260507/epithelial_tissue_comparison_final_removed_cells.tsv" \
        --output-h5ad "${epi_out}" \
        --summary-json "${epi_summary}" \
        --matched-cells-tsv "${ROOT}/data/py/0508/epithelial_input_rm_leiden14_17_20260508/matched_removed_cells.tsv" \
        --compression none
  fi

  local endo_out="${ROOT}/data/py/0508/stromal_input_rm_endothelial6_52_20260508/adata_stromal_subclustered_rm_endothelial6_52_20260508.h5ad"
  local endo_summary="${ROOT}/data/py/0508/stromal_input_rm_endothelial6_52_20260508/filter_summary.json"
  if [[ -s "${endo_out}" && -s "${endo_summary}" ]]; then
    echo "[SKIP] endothelial_filter_input already exists: ${endo_out}"
  else
    run_logged endothelial_filter_input \
      "${PY}" "${ROOT}/script/py/filter_h5ad_by_removed_cells_20260508.py" \
        --input-h5ad "${ROOT}/data/py/0120/stromal_analysis_unified/results/subcluster_unified_v2_20260128/adata_stromal_subclustered_FINAL_v2_20260128.h5ad" \
        --removed-cells "${ROOT}/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414/posthoc_cluster_removal_20260507/stromal_endothelial_tissue_comparison_final_removed_cells.tsv" \
        --output-h5ad "${endo_out}" \
        --summary-json "${endo_summary}" \
        --matched-cells-tsv "${ROOT}/data/py/0508/stromal_input_rm_endothelial6_52_20260508/matched_removed_cells.tsv" \
        --compression none
  fi
}

run_python_references() {
  local bcell_h5ad="${ROOT}/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508/bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad"
  if [[ -s "${bcell_h5ad}" ]]; then
    echo "[SKIP] bcell_scvi_scanvi_20260508 already exists: ${bcell_h5ad}"
  else
    run_logged bcell_scvi_scanvi_20260508 \
      "${PY}" "${ROOT}/script/py/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508.py"
  fi

  local epithelial_h5ad="${ROOT}/data/py/0508/epithelial_scanvi_rm_leiden14_17_20260508/epithelial_scanvi_rm_leiden14_17_SELF_for_R.h5ad"
  if [[ -s "${epithelial_h5ad}" ]]; then
    echo "[SKIP] epithelial_scvi_scanvi_20260508 already exists: ${epithelial_h5ad}"
  else
    run_logged epithelial_scvi_scanvi_20260508 \
      env EPI_SCANVI_INPUT_H5AD="${ROOT}/data/py/0508/epithelial_input_rm_leiden14_17_20260508/epithelial_with_subclusters_rm_leiden14_17_20260508.h5ad" \
          EPI_SCANVI_OUTPUT_DIR="${ROOT}/data/py/0508/epithelial_scanvi_rm_leiden14_17_20260508" \
          EPI_SCANVI_FINAL_H5AD="${epithelial_h5ad}" \
          "${PY}" "${ROOT}/script/py/epithelial_scvi_scanvi_20260417_v2_8_gpu_for_r.py"
  fi

  local endothelial_h5ad="${ROOT}/data/py/0508/stromal_branch_rerun_rm_endothelial6_52_20260508/endothelial/adata_endothelial_reference_v1_5_branchwise.h5ad"
  if [[ -s "${endothelial_h5ad}" ]]; then
    echo "[SKIP] endothelial_branch_scvi_scanvi_20260508 already exists: ${endothelial_h5ad}"
  else
    run_logged endothelial_branch_scvi_scanvi_20260508 \
      "${PY}" "${ROOT}/script/py/stromal_reintegration_branchwise_scvi_scanvi_20260407_v1_5.py" \
        --input-h5ad "${ROOT}/data/py/0508/stromal_input_rm_endothelial6_52_20260508/adata_stromal_subclustered_rm_endothelial6_52_20260508.h5ad" \
        --output-dir "${ROOT}/data/py/0508/stromal_branch_rerun_rm_endothelial6_52_20260508"
  fi

  local tnk_h5ad="${ROOT}/data/py/0508/tnk_scvi_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508/adata_tnk_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508.h5ad"
  if [[ -s "${tnk_h5ad}" ]]; then
    echo "[SKIP] tnk_scvi_scanvi_20260508 already exists: ${tnk_h5ad}"
  else
    run_logged tnk_scvi_scanvi_20260508 \
      "${PY}" "${ROOT}/script/py/tnk_scvi_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508.py"
  fi
}

run_r_lineage() {
  local name="$1"
  local wrapper="$2"
  local output_dir="$3"
  local log_file="${LOG_ROOT}/${name}_r_pipeline.log"
  echo "[$(date '+%F %T')] START ${name} R pipeline" | tee "${log_file}"
  env PYTHONNOUSERSITE=1 PYTHONPATH= RETICULATE_PYTHON="${PY}" LD_LIBRARY_PATH="${R_LD_LIBRARY_PATH}" "${RSCRIPT}" "${wrapper}" 2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} R pipeline status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi
  echo "[$(date '+%F %T')] START ${name} pathway barplots" | tee -a "${log_file}"
  env PYTHONNOUSERSITE=1 PYTHONPATH= RETICULATE_PYTHON="${PY}" LD_LIBRARY_PATH="${R_LD_LIBRARY_PATH}" "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_enrichment_barplots_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} pathway barplots status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi
  echo "[$(date '+%F %T')] START ${name} visual/LLM companions" | tee -a "${log_file}"
  env PYTHONNOUSERSITE=1 PYTHONPATH= RETICULATE_PYTHON="${PY}" LD_LIBRARY_PATH="${R_LD_LIBRARY_PATH}" "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_visual_llm_companion_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} visual/LLM companions status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi
  echo "[$(date '+%F %T')] START ${name} companion LLM batch" | tee -a "${log_file}"
  env PYTHONNOUSERSITE=1 PYTHONPATH= RETICULATE_PYTHON="${PY}" LD_LIBRARY_PATH="${R_LD_LIBRARY_PATH}" "${RSCRIPT}" "${ROOT}/script/R/smc_anno_figure_llm_batch_20260507.R" --output-dir "${output_dir}" --statuses queued_for_batch,error --timeout-sec 240 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} R pipeline+barplots+companions+LLM status=${status}" | tee -a "${log_file}"
  return "${status}"
}

run_r_pipelines_parallel() {
  run_r_lineage bcell \
    "${ROOT}/script/R/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508.R" \
    "${ROOT}/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508" &
  pid_b=$!

  run_r_lineage epithelial \
    "${ROOT}/script/R/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508.R" \
    "${ROOT}/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508" &
  pid_e=$!

  run_r_lineage endothelial \
    "${ROOT}/script/R/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508.R" \
    "${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508" &
  pid_endo=$!

  run_r_lineage tnk \
    "${ROOT}/script/R/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508.R" \
    "${ROOT}/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508" &
  pid_t=$!

  wait "${pid_b}"
  wait "${pid_e}"
  wait "${pid_endo}"
  wait "${pid_t}"
}

echo "=== 20260508 modified-lineage full reruns ==="
echo "Logs: ${LOG_ROOT}"
echo "LLM_SCREEN_PARALLEL_WORKERS=${LLM_SCREEN_PARALLEL_WORKERS}"

filter_h5ad_inputs
run_python_references
run_r_pipelines_parallel

echo "=== All modified-lineage reruns completed successfully ==="
