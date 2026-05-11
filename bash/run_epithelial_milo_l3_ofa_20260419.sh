#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/home/h2048"
MILO_PYTHON_BIN="/home/h2048/miniconda3/envs/scarches_stable_pertpy/bin/python"
MILO_SCRIPT="/home/h2048/script/py/epithelial_milopy_tissue_celltype_20260419_v1.py"
R_BIN="/usr/bin/Rscript"
L3_OFA_SCRIPT="/home/h2048/script/R/epithelial_l3_ofa_only_20260419.R"
MERGE_SCRIPT="/home/h2048/temp/merge_epithelial_l3_ofa_batches_20260419.R"
MILO_TEST_SCRIPT="/home/h2048/temp/test_epithelial_pertpy_milo_outputs_20260419.R"
L3_OFA_TEST_SCRIPT="/home/h2048/temp/test_epithelial_l3_ofa_outputs_20260419.R"
RUN_DATE="${RUN_DATE:-$(date '+%Y%m%d')}"
LOG_DIR="/home/h2048/logs/${RUN_DATE}"
LOG_FILE="${LOG_DIR}/run_epithelial_milo_l3_ofa_20260419.log"
PID_FILE="${LOG_DIR}/run_epithelial_milo_l3_ofa_20260419.pid"
REPORTS_DIR="/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/reports"
L3_OFA_REPORT_DIR="${REPORTS_DIR}/l3_ofa"
BATCH_DIR="/home/h2048/temp/epithelial_l3_ofa_batches_20260419"

mkdir -p "${LOG_DIR}" "${BATCH_DIR}"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${LOG_FILE}"
}

run_clean_env() {
  env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 "$@" 2>&1 | tee -a "${LOG_FILE}"
}

snapshot_l3_ofa_batch() {
  local batch_id="$1"
  cp "${L3_OFA_REPORT_DIR}/l3_ofa_vs_rest_summary.tsv" "${BATCH_DIR}/${batch_id}_l3_ofa_vs_rest_summary.tsv"
  cp "${L3_OFA_REPORT_DIR}/l3_ofa_inter_tissue_summary.tsv" "${BATCH_DIR}/${batch_id}_l3_ofa_inter_tissue_summary.tsv"
  cp "${L3_OFA_REPORT_DIR}/l3_ofa_vs_rest_all.rds" "${BATCH_DIR}/${batch_id}_l3_ofa_vs_rest_all.rds"
  cp "${L3_OFA_REPORT_DIR}/l3_ofa_inter_tissue_all.rds" "${BATCH_DIR}/${batch_id}_l3_ofa_inter_tissue_all.rds"
}

run_l3_batch() {
  local batch_id="$1"
  shift
  local append_args=()
  if [[ "$batch_id" != "batch1" ]]; then
    append_args=(--append)
  fi
  log "[L3-OFA] ${batch_id}: $*"
  run_clean_env "${R_BIN}" "${L3_OFA_SCRIPT}" "${append_args[@]}" "$@"
  snapshot_l3_ofa_batch "${batch_id}"
}

cd "${ROOT_DIR}"
printf '%s\n' "$$" > "${PID_FILE}"

log "[START] epithelial Milo + batched L3 OFA rerun"
log "[INFO] Milo Python: ${MILO_PYTHON_BIN}"
log "[INFO] Milo script: ${MILO_SCRIPT}"
log "[INFO] L3 OFA script: ${L3_OFA_SCRIPT}"

log "[STEP 1/4] Run epithelial pertpy Milo (L2 + L3)"
run_clean_env "${MILO_PYTHON_BIN}" "${MILO_SCRIPT}" --no-write-milo-h5ad

log "[STEP 2/4] Run batched epithelial L3 OFA"
run_l3_batch batch1 AT1_Canonical AT1_MatrixRemodeling AT2 Basal_Cycling Basal_EMT_ECM
run_l3_batch batch2 Basal_Inflammatory Basal_Progenitor Ciliated_Cycling_Immature Ciliated_Mature Ciliogenesis_Deuterosomal
run_l3_batch batch3a Club Goblet
run_l3_batch batch3b Goblet_Defense_DUOX2
run_l3_batch batch3c Ionocyte_Brush SMG_Duct_Secretory_Defense
run_l3_batch batch4a SMG_Mucous SMG_Serous Squamous_Metaplasia
run_l3_batch batch4b Suprabasal_Cycling Suprabasal_Progenitor

log "[STEP 3/4] Merge batched L3 OFA summaries/RDS"
run_clean_env "${R_BIN}" "${MERGE_SCRIPT}"

log "[STEP 4/4] Validate Milo and L3 OFA outputs"
run_clean_env "${R_BIN}" "${MILO_TEST_SCRIPT}"
run_clean_env "${R_BIN}" "${L3_OFA_TEST_SCRIPT}"

log "[DONE] epithelial Milo + batched L3 OFA rerun completed"
