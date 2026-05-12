#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 2026-05-13 lightweight CHOIR restart for modified-lineage recovery
# ------------------------------------------------------------------------------
# Assumes stale heavy CHOIR processes have been stopped. Runs B cell and
# endothelial lightweight CHOIR wrappers in parallel, then regenerates pathway
# barplots + visual/LLM companions for each completed lineage.
# ===============================================================================

ROOT="/home/h2048"
PY="${ROOT}/miniconda3/envs/scvi_env/bin/python"
RSCRIPT="/usr/bin/Rscript"
LOG_ROOT="${ROOT}/logs/20260513/modified_lineage_choir_light_restart"
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

run_lineage() {
  local name="$1"
  local wrapper="$2"
  local output_dir="$3"
  local final_rds="$4"
  local final_h5ad="$5"
  local log_file="${LOG_ROOT}/${name}_choir_light_restart.log"

  echo "[$(date '+%F %T')] START ${name} lightweight CHOIR restart" | tee "${log_file}"
  echo "wrapper=${wrapper}" | tee -a "${log_file}"
  echo "output_dir=${output_dir}" | tee -a "${log_file}"

  run_r_env "${RSCRIPT}" "${wrapper}" 2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} wrapper status=${status}" | tee -a "${log_file}"
  if [[ "${status}" -ne 0 ]]; then
    return "${status}"
  fi
  if [[ ! -s "${final_rds}" || ! -s "${final_h5ad}" ]]; then
    echo "[$(date '+%F %T')] FAIL ${name}: final outputs missing" | tee -a "${log_file}"
    return 2
  fi

  echo "[$(date '+%F %T')] START ${name} pathway barplots" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_enrichment_barplots_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} pathway barplots status=${status}" | tee -a "${log_file}"
  if [[ "${status}" -ne 0 ]]; then
    return "${status}"
  fi

  echo "[$(date '+%F %T')] START ${name} visual/LLM companions" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_visual_llm_companion_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} visual/LLM companions status=${status}" | tee -a "${log_file}"
  if [[ "${status}" -ne 0 ]]; then
    return "${status}"
  fi

  echo "[$(date '+%F %T')] START ${name} companion LLM batch" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/smc_anno_figure_llm_batch_20260507.R" \
    --output-dir "${output_dir}" \
    --statuses queued_for_batch,error \
    --timeout-sec 240 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} lightweight CHOIR restart status=${status}" | tee -a "${log_file}"
  return "${status}"
}

echo "=== 20260513 lightweight CHOIR restart ==="
echo "Logs: ${LOG_ROOT}"
echo "LLM_SCREEN_PARALLEL_WORKERS=${LLM_SCREEN_PARALLEL_WORKERS}"
echo "R_LD_LIBRARY_PATH=${R_LD_LIBRARY_PATH}"

run_lineage bcell \
  "${ROOT}/script/R/bcell_tissue_comparison_v2_6_10_c22_c13_c25_c14drop_choir_light_20260513.R" \
  "${ROOT}/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508" \
  "${ROOT}/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/bcell_tissue_comparison_final.rds" \
  "${ROOT}/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/bcell_tissue_comparison_final.h5ad" &
pid_bcell=$!

run_lineage endothelial \
  "${ROOT}/script/R/stromal_endothelial_tissue_comparison_v1_1_5_rm_choir6_52_choir_light_20260513.R" \
  "${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508" \
  "${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508/stromal_endothelial_tissue_comparison_final.rds" \
  "${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508/stromal_endothelial_tissue_comparison_final.h5ad" &
pid_endo=$!

status=0
wait "${pid_bcell}" || status=$?
wait "${pid_endo}" || status=$?

if [[ "${status}" -eq 0 ]]; then
  echo "=== Lightweight CHOIR restart completed successfully ==="
else
  echo "=== Lightweight CHOIR restart finished with status ${status} ==="
fi
exit "${status}"
