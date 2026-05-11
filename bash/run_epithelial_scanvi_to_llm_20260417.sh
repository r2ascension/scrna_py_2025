#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/h2048"
PYTHON_BIN="/home/h2048/miniconda3/envs/scvi_env/bin/python"
PY_SCRIPT="/home/h2048/script/py/epithelial_scvi_scanvi_20260417_v2_8_gpu_for_r.py"
MILO_PYTHON_BIN="/home/h2048/miniconda3/envs/scarches_stable_pertpy/bin/python"
MILO_SCRIPT="/home/h2048/script/py/epithelial_milopy_tissue_celltype_20260419_v1.py"
R_BIN="/usr/bin/Rscript"
R_SCRIPT="/home/h2048/script/R/epithelial_tissue_comparison_v1_3_2_scanvi_rerun_20260417.R"
RUN_DATE="${RUN_DATE:-$(date '+%Y%m%d')}"
LOG_DIR="/home/h2048/logs/${RUN_DATE}"
LOG_FILE="${LOG_DIR}/epithelial_scanvi_to_llm_20260417.log"
PID_FILE="${LOG_DIR}/epithelial_scanvi_to_llm_20260417.pid"

mkdir -p "${LOG_DIR}"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${LOG_FILE}"
}

run_clean_env() {
  env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 "$@" 2>&1 | tee -a "${LOG_FILE}"
}

cd "${ROOT_DIR}"
printf '%s\n' "$$" > "${PID_FILE}"

log "[START] epithelial scanvi -> Milo -> R tissue comparison rerun"
log "[INFO] scANVI Python: ${PYTHON_BIN}"
log "[INFO] scANVI script: ${PY_SCRIPT}"
log "[INFO] Milo Python: ${MILO_PYTHON_BIN}"
log "[INFO] Milo script: ${MILO_SCRIPT}"
log "[INFO] R script: ${R_SCRIPT}"
log "[INFO] CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"

log "[STEP 1/3] GPU scANVI rerun"
run_clean_env "${PYTHON_BIN}" "${PY_SCRIPT}"

log "[STEP 2/3] pertpy Milo tissue DA for epithelial L2/L3"
run_clean_env "${MILO_PYTHON_BIN}" "${MILO_SCRIPT}"

log "[STEP 3/3] R tissue comparison rerun through LLM outputs"
run_clean_env "${R_BIN}" "${R_SCRIPT}"

log "[DONE] epithelial scanvi -> Milo -> R tissue comparison rerun completed"
