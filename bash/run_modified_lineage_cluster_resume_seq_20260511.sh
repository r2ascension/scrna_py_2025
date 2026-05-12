#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 2026-05-11 sequential recovery for unfinished 2026-05-08 modified-lineage reruns
# ------------------------------------------------------------------------------
# B cell       : rerun downstream R stages + CHOIR sequentially.
# Epithelial   : rerun downstream R stages with Leiden backend, reusing
#                h5ad obs column `leiden_Epithelial_res0.8` (no CHOIR).
# Endothelial  : rerun downstream R stages + CHOIR sequentially, loading cached
#                pairwise/ssGSEA artifacts from the current 0508 output dir.
# ==============================================================================

ROOT="/home/h2048"
PY="${ROOT}/miniconda3/envs/scvi_env/bin/python"
RSCRIPT="/usr/bin/Rscript"
LOG_ROOT="${ROOT}/logs/20260511/modified_lineage_cluster_resume_seq"
mkdir -p "${LOG_ROOT}"

export PYTHONNOUSERSITE=1
export LLM_SCREEN_PARALLEL_WORKERS="${LLM_SCREEN_PARALLEL_WORKERS:-3}"
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

run_r_lineage() {
  local name="$1"
  local wrapper="$2"
  local output_dir="$3"
  local final_rds="$4"
  local final_h5ad="$5"
  local log_file="${LOG_ROOT}/${name}_r_pipeline_resume.log"

  echo "[$(date '+%F %T')] START ${name} sequential recovery" | tee "${log_file}"
  echo "wrapper=${wrapper}" | tee -a "${log_file}"
  echo "output_dir=${output_dir}" | tee -a "${log_file}"

  run_r_env "${RSCRIPT}" "${wrapper}" 2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} R wrapper status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi

  if [[ ! -s "${final_rds}" || ! -s "${final_h5ad}" ]]; then
    echo "[$(date '+%F %T')] FAIL ${name}: final outputs missing after wrapper" | tee -a "${log_file}"
    echo "  final_rds=${final_rds}" | tee -a "${log_file}"
    echo "  final_h5ad=${final_h5ad}" | tee -a "${log_file}"
    return 2
  fi

  echo "[$(date '+%F %T')] START ${name} pathway barplots" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_enrichment_barplots_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} pathway barplots status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi

  echo "[$(date '+%F %T')] START ${name} visual/LLM companions" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/tissue_comparison_visual_llm_companion_20260508.R" "${output_dir}" 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  if [[ "${status}" -ne 0 ]]; then
    echo "[$(date '+%F %T')] FAIL ${name} visual/LLM companions status=${status}" | tee -a "${log_file}"
    return "${status}"
  fi

  echo "[$(date '+%F %T')] START ${name} companion LLM batch" | tee -a "${log_file}"
  run_r_env "${RSCRIPT}" "${ROOT}/script/R/smc_anno_figure_llm_batch_20260507.R" \
    --output-dir "${output_dir}" \
    --statuses queued_for_batch,error \
    --timeout-sec 240 2>&1 | tee -a "${log_file}"
  status=${PIPESTATUS[0]}
  echo "[$(date '+%F %T')] END ${name} sequential recovery status=${status}" | tee -a "${log_file}"
  return "${status}"
}

echo "=== 20260511 sequential recovery for unfinished modified-lineage reruns ==="
echo "Logs: ${LOG_ROOT}"
echo "LLM_SCREEN_PARALLEL_WORKERS=${LLM_SCREEN_PARALLEL_WORKERS}"
echo "R_LD_LIBRARY_PATH=${R_LD_LIBRARY_PATH}"

run_r_lineage bcell \
  "${ROOT}/script/R/bcell_tissue_comparison_v2_6_9_c22_c13_c25_c14drop_choir_seq_resume_20260511.R" \
  "${ROOT}/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508" \
  "${ROOT}/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/bcell_tissue_comparison_final.rds" \
  "${ROOT}/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/bcell_tissue_comparison_final.h5ad"

run_r_lineage epithelial \
  "${ROOT}/script/R/epithelial_tissue_comparison_v1_3_4_rm_leiden14_17_leiden_seq_resume_20260511.R" \
  "${ROOT}/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508" \
  "${ROOT}/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/epithelial_tissue_comparison_final.rds" \
  "${ROOT}/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/epithelial_tissue_comparison_final.h5ad"

run_r_lineage endothelial \
  "${ROOT}/script/R/stromal_endothelial_tissue_comparison_v1_1_3_rm_choir6_52_choir_seq_resume_20260511.R" \
  "${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508" \
  "${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508/stromal_endothelial_tissue_comparison_final.rds" \
  "${ROOT}/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508/stromal_endothelial_tissue_comparison_final.h5ad"

echo "=== Sequential recovery completed successfully ==="
