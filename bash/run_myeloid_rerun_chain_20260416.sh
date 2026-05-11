#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/h2048"
PYTHON_BIN="/home/h2048/miniconda3/envs/scvi_env/bin/python"
SCANVI_SCRIPT="/home/h2048/script/py/myeloid_scvi_scanvi_v2_4_cpu_rerun_20260416.py"
PATCH_SCRIPT="/home/h2048/script/py/myeloid_L3refined_tissueaware_patch_20260416_v1.py"
R_SCRIPT="/home/h2048/script/R/myeloid_tissue_comparison_v1_2_3_20260416.R"
VALIDATE_SCRIPT="/home/h2048/script/py/validate_myeloid_rerun_outputs_20260416.py"
LOG_DIR="/home/h2048/logs/20260417"
LOG_FILE="$LOG_DIR/myeloid_rerun_chain_20260416.log"

mkdir -p "$LOG_DIR"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG_FILE"
}

run_clean_env() {
  env -u LD_LIBRARY_PATH -u PYTHONPATH "$@" 2>&1 | tee -a "$LOG_FILE"
}

cd "$ROOT_DIR"

log "Starting myeloid rerun chain"
log "Step 1/4: scanvi CPU rerun"
run_clean_env "$PYTHON_BIN" "$SCANVI_SCRIPT"

log "Step 2/4: tissue-aware patch"
run_clean_env "$PYTHON_BIN" "$PATCH_SCRIPT"

log "Step 3/4: R tissue comparison rerun"
run_clean_env /usr/bin/Rscript "$R_SCRIPT"

log "Step 4/4: validate outputs"
run_clean_env "$PYTHON_BIN" "$VALIDATE_SCRIPT"

log "Myeloid rerun chain completed successfully"
