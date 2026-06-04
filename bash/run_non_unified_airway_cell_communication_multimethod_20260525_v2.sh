#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/home/h2048"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/usr/bin/Rscript}"
TOP_N="${TOP_N:-25}"
LLM_ENABLED="${LLM_ENABLED:-true}"
export PYTHONPATH="/home/h2048/script/py${PYTHONPATH:+:${PYTHONPATH}}"
FORCE_FLAG=""
if [[ "${FORCE:-false}" == "true" ]]; then
  FORCE_FLAG="--force"
fi

PREPARED_MANIFEST="${PREPARED_MANIFEST:-/home/h2048/data/py/0525/non_unified_airway/communication_lineage_native_export_v2/prepared_inputs_manifest.tsv}"
PY_COMM_DIR="${PY_COMM_DIR:-/home/h2048/data/py/0525/non_unified_airway/communication_lineage_native_export_v2}"
R_COMM_DIR="${R_COMM_DIR:-/home/h2048/data/R/0525/non_unified_airway/communication}"

LIANA_DIR="${LIANA_DIR:-${PY_COMM_DIR}/liana_20260525_v2}"
CELLPHONEDB_DIR="${CELLPHONEDB_DIR:-${PY_COMM_DIR}/cellphonedb_20260525_v2}"
CONSENSUS_DIR="${CONSENSUS_DIR:-${PY_COMM_DIR}/consensus_20260525_v2}"
CELLCHAT_DIR="${CELLCHAT_DIR:-${R_COMM_DIR}/cellchat_lineage_native_export_v2}"
SUMMARY_DIR="${SUMMARY_DIR:-${R_COMM_DIR}/summary_lineage_native_export_v2_20260525_v2}"

cd "${REPO_ROOT}"

echo "[multimethod_cci_v2] Step 1/4: LIANA pairwise"
"${PYTHON_BIN}" /home/h2048/script/py/non_unified_airway/05b_run_liana_pairwise_L2_L3.py \
  --prepared-manifest "${PREPARED_MANIFEST}" \
  --out-dir "${LIANA_DIR}" \
  ${FORCE_FLAG}

echo "[multimethod_cci_v2] Step 2/4: CellPhoneDB pairwise"
"${PYTHON_BIN}" /home/h2048/script/py/non_unified_airway/05c_run_cellphonedb_pairwise_L2_L3.py \
  --prepared-manifest "${PREPARED_MANIFEST}" \
  --out-dir "${CELLPHONEDB_DIR}" \
  ${FORCE_FLAG}

echo "[multimethod_cci_v2] Step 3/4: rich multimethod consensus collector"
"${PYTHON_BIN}" /home/h2048/script/py/non_unified_airway/05d_collect_cell_communication_consensus_20260525_v2.py \
  --liana-dir "${LIANA_DIR}" \
  --cellphonedb-dir "${CELLPHONEDB_DIR}" \
  --cellchat-dir "${CELLCHAT_DIR}" \
  --out-dir "${CONSENSUS_DIR}" \
  --top-n "${TOP_N}"

echo "[multimethod_cci_v2] Step 4/4: visualization + LLM packaging"
env -u LD_LIBRARY_PATH -u PYTHONPATH "${RSCRIPT_BIN}" /home/h2048/script/R/non_unified_airway/06b_cell_communication_visualization_summary_20260525_v2.R \
  --cellchat-dir "${CELLCHAT_DIR}" \
  --liana-dir "${LIANA_DIR}" \
  --cellphonedb-dir "${CELLPHONEDB_DIR}" \
  --consensus-dir "${CONSENSUS_DIR}" \
  --output-dir "${SUMMARY_DIR}" \
  --top-n "${TOP_N}" \
  --llm-enabled "${LLM_ENABLED}"

echo "[multimethod_cci_v2] Done"
echo "  LIANA_DIR=${LIANA_DIR}"
echo "  CELLPHONEDB_DIR=${CELLPHONEDB_DIR}"
echo "  CONSENSUS_DIR=${CONSENSUS_DIR}"
echo "  SUMMARY_DIR=${SUMMARY_DIR}"
