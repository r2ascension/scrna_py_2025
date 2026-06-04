#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/home/h2048"
BASE_ENV="${BASE_ENV:-/home/h2048/miniconda3/envs/scarches_stable}"
OVERLAY="${OVERLAY:-/home/h2048/temp/cellcomm_overlay_20260526/site-packages}"
WRAPBIN="${WRAPBIN:-/home/h2048/temp/cellcomm_overlay_20260526/bin}"
PYTHON_BIN="${PYTHON_BIN:-${BASE_ENV}/bin/python}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/usr/bin/Rscript}"

PREPARED_MANIFEST="${PREPARED_MANIFEST:-/home/h2048/data/py/0525/non_unified_airway/communication_lineage_native_export_v2/prepared_inputs_manifest.tsv}"
PY_COMM_DIR="${PY_COMM_DIR:-/home/h2048/data/py/20260527/non_unified_airway/communication_lineage_native_export_v3}"
R_COMM_DIR="${R_COMM_DIR:-/home/h2048/data/R/20260527/non_unified_airway/communication}"
LOG_DIR="${LOG_DIR:-/home/h2048/logs/20260527}"

CELLCHAT_DIR="${CELLCHAT_DIR:-/home/h2048/data/R/0525/non_unified_airway/communication/cellchat_lineage_native_export_v2}"
LIANA_DIR="${LIANA_DIR:-${PY_COMM_DIR}/liana_20260527_v3}"
CELLPHONEDB_DIR="${CELLPHONEDB_DIR:-${PY_COMM_DIR}/cellphonedb_20260527_v3}"
CONSENSUS_DIR="${CONSENSUS_DIR:-${PY_COMM_DIR}/consensus_20260527_v3}"
SUMMARY_DIR="${SUMMARY_DIR:-${R_COMM_DIR}/summary_lineage_native_export_v3_20260527}"
MASTER_LOG="${MASTER_LOG:-${LOG_DIR}/non_unified_airway_multimethod_20260527_v3.log}"

TOP_N="${TOP_N:-25}"
LLM_ENABLED="${LLM_ENABLED:-true}"
RUN_CELLCHAT="${RUN_CELLCHAT:-auto}"
CELLCHAT_NBOOT="${CELLCHAT_NBOOT:-100}"
CELLCHAT_WORKERS="${CELLCHAT_WORKERS:-4}"
CELLPHONEDB_THREADS="${CELLPHONEDB_THREADS:-8}"
CELLPHONEDB_ITERATIONS="${CELLPHONEDB_ITERATIONS:-1000}"

FORCE_FLAG=""
if [[ "${FORCE:-false}" == "true" ]]; then
  FORCE_FLAG="--force"
fi

mkdir -p "${PY_COMM_DIR}" "${R_COMM_DIR}" "${LOG_DIR}"

export LD_LIBRARY_PATH="${BASE_ENV}/lib/python3.10/site-packages/torch/lib:${BASE_ENV}/lib"
export _SCARCHES_STABLE_LD_FIXED=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${OVERLAY}:/home/h2048/script/py"
export PATH="${WRAPBIN}:${PATH}"

cd "${REPO_ROOT}"

log() {
  local ts
  ts="$(date -Iseconds)"
  echo "[$ts] $*" | tee -a "${MASTER_LOG}"
}

run_step() {
  local step_name="$1"
  shift
  log "${step_name} :: START"
  "$@" 2>&1 | tee -a "${MASTER_LOG}"
  log "${step_name} :: DONE"
}

has_ok_summary() {
  local path="$1"
  [[ -f "${path}" ]] && grep -q '"status"[[:space:]]*:[[:space:]]*"ok"' "${path}"
}

should_run_cellchat() {
  case "${RUN_CELLCHAT}" in
    true)
      return 0
      ;;
    false)
      return 1
      ;;
    auto)
      if has_ok_summary "${CELLCHAT_DIR}/run_summary.json"; then
        return 1
      fi
      return 0
      ;;
    *)
      echo "Unsupported RUN_CELLCHAT=${RUN_CELLCHAT}; expected auto|true|false" >&2
      return 2
      ;;
  esac
}

log "[multimethod_cci_v3] repo_root=${REPO_ROOT}"
log "[multimethod_cci_v3] prepared_manifest=${PREPARED_MANIFEST}"
log "[multimethod_cci_v3] liana_dir=${LIANA_DIR}"
log "[multimethod_cci_v3] cellphonedb_dir=${CELLPHONEDB_DIR}"
log "[multimethod_cci_v3] consensus_dir=${CONSENSUS_DIR}"
log "[multimethod_cci_v3] summary_dir=${SUMMARY_DIR}"
log "[multimethod_cci_v3] cellchat_dir=${CELLCHAT_DIR}"
log "[multimethod_cci_v3] run_cellchat=${RUN_CELLCHAT} llm_enabled=${LLM_ENABLED} force=${FORCE:-false}"

if should_run_cellchat; then
  mkdir -p "${CELLCHAT_DIR}"
  run_step "Step 0/4 CellChat pairwise" \
    env -u LD_LIBRARY_PATH -u PYTHONPATH "${RSCRIPT_BIN}" \
      /home/h2048/script/R/non_unified_airway/06a_cellchat_pairwise_L2_L3.R \
      --prepared-manifest "${PREPARED_MANIFEST}" \
      --output-dir "${CELLCHAT_DIR}" \
      --nboot "${CELLCHAT_NBOOT}" \
      --workers "${CELLCHAT_WORKERS}" \
      ${FORCE_FLAG}
else
  log "[multimethod_cci_v3] Reusing existing CellChat outputs from ${CELLCHAT_DIR}"
fi

run_step "Step 1/4 LIANA pairwise" \
  "${PYTHON_BIN}" \
    /home/h2048/script/py/non_unified_airway/05b_run_liana_pairwise_L2_L3.py \
    --prepared-manifest "${PREPARED_MANIFEST}" \
    --out-dir "${LIANA_DIR}" \
    ${FORCE_FLAG}

run_step "Step 2/4 CellPhoneDB pairwise" \
  "${PYTHON_BIN}" \
    /home/h2048/script/py/non_unified_airway/05c_run_cellphonedb_pairwise_L2_L3_20260526_v2.py \
    --prepared-manifest "${PREPARED_MANIFEST}" \
    --out-dir "${CELLPHONEDB_DIR}" \
    --threads "${CELLPHONEDB_THREADS}" \
    --iterations "${CELLPHONEDB_ITERATIONS}" \
    ${FORCE_FLAG}

run_step "Step 3/4 multimethod consensus collector" \
  "${PYTHON_BIN}" \
    /home/h2048/script/py/non_unified_airway/05d_collect_cell_communication_consensus_20260525_v2.py \
    --liana-dir "${LIANA_DIR}" \
    --cellphonedb-dir "${CELLPHONEDB_DIR}" \
    --cellchat-dir "${CELLCHAT_DIR}" \
    --out-dir "${CONSENSUS_DIR}" \
    --top-n "${TOP_N}"

run_step "Step 4/4 visualization and LLM packaging" \
  env -u LD_LIBRARY_PATH -u PYTHONPATH "${RSCRIPT_BIN}" \
    /home/h2048/script/R/non_unified_airway/06b_cell_communication_visualization_summary_20260525_v2.R \
    --cellchat-dir "${CELLCHAT_DIR}" \
    --liana-dir "${LIANA_DIR}" \
    --cellphonedb-dir "${CELLPHONEDB_DIR}" \
    --consensus-dir "${CONSENSUS_DIR}" \
    --output-dir "${SUMMARY_DIR}" \
    --top-n "${TOP_N}" \
    --llm-enabled "${LLM_ENABLED}"

log "[multimethod_cci_v3] COMPLETE"
log "  LIANA_DIR=${LIANA_DIR}"
log "  CELLPHONEDB_DIR=${CELLPHONEDB_DIR}"
log "  CONSENSUS_DIR=${CONSENSUS_DIR}"
log "  SUMMARY_DIR=${SUMMARY_DIR}"