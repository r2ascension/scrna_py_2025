#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 2026-05-12 parallel front-load for modified-lineage recovery
# ------------------------------------------------------------------------------
# Keeps the existing 20260511 sequential launcher running B cell CHOIR, while
# starting independent work for other output directories:
#   - epithelial: full Leiden recovery (no CHOIR)
#   - endothelial: no-CHOIR front-load only; full CHOIR remains for later
# ===============================================================================

ROOT="/home/h2048"
PY="${ROOT}/miniconda3/envs/scvi_env/bin/python"
RSCRIPT="/usr/bin/Rscript"
LOG_ROOT="${ROOT}/logs/20260512/modified_lineage_frontload_parallel"
mkdir -p "${LOG_ROOT}"

export PYTHONNOUSERSITE=1
export LLM_SCREEN_PARALLEL_WORKERS="${LLM_SCREEN_PARALLEL_WORKERS:-2}"
export LLM_SCREEN_PARALLEL_STAGGER_SEC="${LLM_SCREEN_PARALLEL_STAGGER_SEC:-0.5}"
export STANDARDIZE_LLM_RETRY_SLEEP_SEC="${STANDARDIZE_LLM_RETRY_SLEEP_SEC:-1}"
SCVI_ENV_LIB="${ROOT}/miniconda3/envs/scvi_env/lib"
SCVI_ENV_CUDA_LIB="${ROOT}/miniconda3/envs/scvi_env/targets/x86_64-linux/lib"
R_LD_LIBRARY_PATH="${SCVI_ENV_LIB}:${SCVI_ENV_CUDA_LIB}"

run_r_env() {
  env PYTHONNOUSERSITE=1 \
      PYTHONPATH= \
      RETICULATE_PYTHON="${PY}" \
      LD_LIBRARY_PATH="${R_LD_LIBRARY_PATH}" \
      "$@"
}

run_epithelial_full_leiden() {
  local name="epithelial"
  local output_dir="${ROOT}/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508"
  local final_rds="${output_dir}/epithelial_tissue_comparison_final.rds"
  local final_h5ad="${output_dir}/epithelial_tissue_comparison_final.h5ad"
  local log_file="${LOG_ROOT}/${name}_full_leiden_frontload.log"

  echo "[$(date '+%F %T')] START ${name} full Leiden front-load" | tee "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/epithelial_tissue_comparison_v1_3_4_rm_leiden14_17_leiden_seq_resume_20260511.R" 2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} wrapper status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi
  if [[ ! -s "${final_rds}" || ! -s "${final_h5ad}" ]]; then
    echo "[$(date '+%F %T')] FAIL ${name}: final outputs missing" | tee -a "${log_file}"
    return 2
  fi

  echo "[$(date '+%F %T')] START ${name} pathway barplots" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_enrichment_barplots_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  [[ "${status}" -eq 0 ]] || return "${status}"

  echo "[$(date '+%F %T')] START ${name} visual/LLM companions" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_visual_llm_companion_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  [[ "${status}" -eq 0 ]] || return "${status}"

  echo "[$(date '+%F %T')] START ${name} companion LLM batch" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/smc_anno_figure_llm_batch_20260507.R" \
    --output-dir "${output_dir}" \
    --statuses queued_for_batch,error \
    --timeout-sec 240 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} full Leiden front-load status=${status}" | tee -a "${log_file}"
  return "${status}"
}

run_endothelial_pre_nochoir() {
  local name="endothelial"
  local output_dir="${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508"
  local final_rds="${output_dir}/stromal_endothelial_tissue_comparison_final.rds"
  local final_h5ad="${output_dir}/stromal_endothelial_tissue_comparison_final.h5ad"
  local log_file="${LOG_ROOT}/${name}_pre_nochoir_frontload.log"

  echo "[$(date '+%F %T')] START ${name} no-CHOIR front-load" | tee "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/stromal_endothelial_tissue_comparison_v1_1_4_rm_choir6_52_frontload_nochoir_20260512.R" 2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} no-CHOIR wrapper status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi
  if [[ ! -s "${final_rds}" || ! -s "${final_h5ad}" ]]; then
    echo "[$(date '+%F %T')] WARN ${name}: front-load final outputs missing; CHOIR wrapper can still regenerate later" | tee -a "${log_file}"
  fi

  echo "[$(date '+%F %T')] START ${name} pathway barplots from cached pre-stage artifacts" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_enrichment_barplots_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} no-CHOIR front-load status=${status}" | tee -a "${log_file}"
  return "${status}"
}

echo "=== 20260512 parallel front-load for modified-lineage recovery ==="
echo "Logs: ${LOG_ROOT}"
echo "LLM_SCREEN_PARALLEL_WORKERS=${LLM_SCREEN_PARALLEL_WORKERS}"
echo "R_LD_LIBRARY_PATH=${R_LD_LIBRARY_PATH}"

run_epithelial_full_leiden &
pid_epi=$!
run_endothelial_pre_nochoir &
pid_endo=$!

status=0
wait "${pid_epi}" || status=$?
wait "${pid_endo}" || status=$?

if [[ "${status}" -eq 0 ]]; then
  echo "=== Parallel front-load completed successfully ==="
else
  echo "=== Parallel front-load finished with status ${status} ==="
fi
exit "${status}"
