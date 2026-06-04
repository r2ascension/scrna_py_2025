#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Query CHOIR/Leiden Clustering + OFA + LLM Cleanup — 7-Lineage Orchestrator
# ==============================================================================
# Date: 2026-06-04
#
# For each lineage that has a processed query h5ad file, this script:
#   1. Runs a minimal R tissue-comparison wrapper (clustering + OFA + LLM only)
#   2. Runs the Python cleanup script to remove LLM-flagged outlier clusters
#
# Each lineage runs as an independent job; failures in one lineage do not
# block the others.  Re-running the script skips already-completed lineages.
# ==============================================================================

ROOT="/home/h2048"
PY="${ROOT}/miniconda3/envs/bbknn_env/bin/python"     # for direct Python execution (cleanup script)
RETICULATE_PY="${ROOT}/miniconda3/envs/scvi_env/bin/python"  # for R reticulate (scVI/anndata)
RSCRIPT="/usr/bin/Rscript"
LOG_ROOT="${ROOT}/logs/20260604/query_clustering_cleanup"
R_OUTPUT_BASE="${ROOT}/data/R/20260604"
PY_OUTPUT_BASE="${ROOT}/data/py/20260604"
SCRIPT_DIR="${ROOT}/script"

mkdir -p "${LOG_ROOT}" "${R_OUTPUT_BASE}" "${PY_OUTPUT_BASE}"

# ── Environment ──────────────────────────────────────────────────────────────

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
      RETICULATE_PYTHON="${RETICULATE_PY}" \
      LD_LIBRARY_PATH="${R_LD_LIBRARY_PATH}" \
      "$@"
}

# ── Lineage configuration ───────────────────────────────────────────────────
# Format: name|category|h5ad_path|wrapper_script|shared_engine
# Category A: canonical no-COVID merged files (ready for clustering)
# Category B: best-available merged files

LINEAGES=(
  # Category A — canonical no-COVID merged+scHPL
  "bcell|A|${ROOT}/data/py/0523/full_nocovid_20260523/bcell_merge_schpl_nocovid_20260523/bcell_reference_plus_query_schpl_v1_0_nocovid_20260523.h5ad|BCELL|${ROOT}/script/R/bcell_tissue_comparison_v2_6_2_20260414.R"
  "epithelial|A|${ROOT}/data/py/0523/full_nocovid_20260523/epithelial_merge_schpl_v1_0_nocovid_20260523/epithelial_reference_plus_query_schpl_v1_0_nocovid_20260523.h5ad|EPITHELIAL|${ROOT}/script/R/epithelial_tissue_comparison_v1_3_2_20260413.R"

  # Category B — best-available merged files
  "tnk|B|${ROOT}/data/py/20260318/tcell_only_merged_pipeline_v2_3_filtered/tcell_merged_v2_3_results.h5ad|TNK|${ROOT}/script/R/tnk_tissue_comparison_v2_6_2_20260413_rm_choir.R"
  "myeloid|B|${ROOT}/data/py/20260308/myeloid_only_merged_pipeline_v2_1_filtered/myeloid_merged_v2_1_results.h5ad|MYELOID|${ROOT}/script/R/myeloid_tissue_comparison_v1_2_3_20260416.R"
  "endothelial|B|${ROOT}/data/py/0406/stromal_schpl_v1_1_branchwise/branch_h5ad/adata_stromal_query_schpl_endothelial_v1_1_branchwise.h5ad|ENDOTHELIAL|${ROOT}/script/R/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414.R"
  "fibroblast|B|${ROOT}/data/py/0406/stromal_schpl_v1_1_branchwise/branch_h5ad/adata_stromal_query_schpl_fibroblast_v1_1_branchwise.h5ad|FIBROBLAST|${ROOT}/script/R/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414.R"
  "smc|B|${ROOT}/data/py/0406/stromal_schpl_v1_1_branchwise/branch_h5ad/adata_stromal_query_schpl_smc_v1_1_branchwise.h5ad|SMC|${ROOT}/script/R/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414.R"
)

PY_CLEANUP_SCRIPT="${ROOT}/script/py/core/integration/curate_query_cleanup_20260604.py"
CLUSTER_COLUMN="CHOIR_clusters_0.2"

# ── Resolve shared engines — use the latest available if specified is missing ─

declare -A ENGINE_FALLBACKS=(
  ["BCELL"]="${ROOT}/script/R/bcell_tissue_comparison_v2_6_2_20260414.R"
  ["EPITHELIAL"]="${ROOT}/script/R/epithelial_tissue_comparison_v1_3_2_20260413.R"
  ["TNK"]="${ROOT}/script/R/tnk_tissue_comparison_v2_6_2_20260413_rm_choir.R"
  ["MYELOID"]="${ROOT}/script/R/myeloid_tissue_comparison_v1_2_3_20260416.R"
  ["ENDOTHELIAL"]="${ROOT}/script/R/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414.R"
  ["FIBROBLAST"]="${ROOT}/script/R/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414.R"
  ["SMC"]="${ROOT}/script/R/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414.R"
)

# ── Helper: unified lineage R wrapper ────────────────────────────────────────

UNIFIED_WRAPPER="${ROOT}/script/R/query_choir_cleanup_20260604.R"
WRAPPER_PID_FILE="${LOG_ROOT}/.wrapper_pids"

> "${WRAPPER_PID_FILE}"

run_lineage_r() {
  local name="$1" category="$2" h5ad_path="$3" lineage_tag="$4" shared_engine="$5"
  local output_dir="${R_OUTPUT_BASE}/${name}_query_choir_cleanup_20260604"
  local final_rds="${output_dir}/${name}_tissue_comparison_final.rds"
  local final_h5ad="${output_dir}/${name}_tissue_comparison_final.h5ad"
  local log_file="${LOG_ROOT}/${name}_query_choir_cleanup.log"

  echo "[$(date '+%F %T')] START ${name} (category=${category})" | tee "${log_file}"
  echo "  h5ad=${h5ad_path}" | tee -a "${log_file}"
  echo "  engine=${shared_engine}" | tee -a "${log_file}"
  echo "  output=${output_dir}" | tee -a "${log_file}"

  # Check input h5ad exists
  if [[ ! -f "${h5ad_path}" ]]; then
    echo "[$(date '+%F %T')] SKIP ${name}: h5ad not found at ${h5ad_path}" | tee -a "${log_file}"
    return 2
  fi

  # Check shared engine exists
  if [[ ! -f "${shared_engine}" ]]; then
    echo "[$(date '+%F %T')] WARN ${name}: shared engine not found; trying fallback" | tee -a "${log_file}"
    shared_engine="${ENGINE_FALLBACKS[${lineage_tag}]}"
    if [[ ! -f "${shared_engine}" ]]; then
      echo "[$(date '+%F %T')] FAIL ${name}: no shared engine available" | tee -a "${log_file}"
      return 3
    fi
  fi

  # Skip if already completed (resume-safe)
  if [[ -s "${final_rds}" ]] && [[ -s "${final_h5ad}" ]] && \
     [[ -s "${output_dir}/reports/llm_choir_discovery_screen.tsv" ]]; then
    echo "[$(date '+%F %T')] SKIP ${name}: R outputs already exist (resume-safe)" | tee -a "${log_file}"
    return 0
  fi

  mkdir -p "${output_dir}"

  run_r_env "${RSCRIPT}" "${UNIFIED_WRAPPER}" \
    --lineage-tag "${lineage_tag}" \
    --h5ad-path "${h5ad_path}" \
    --output-dir "${output_dir}" \
    --shared-engine "${shared_engine}" \
    2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}

  echo "[$(date '+%F %T')] END ${name} R status=${status}" | tee -a "${log_file}"

  if [[ "${status}" -ne 0 ]]; then
    return "${status}"
  fi

  # Verify outputs
  if [[ ! -s "${final_rds}" ]] || [[ ! -s "${final_h5ad}" ]]; then
    echo "[$(date '+%F %T')] FAIL ${name}: final outputs missing after R run" | tee -a "${log_file}"
    return 4
  fi

  return 0
}

# ── Python cleanup stage ─────────────────────────────────────────────────────

run_lineage_py_cleanup() {
  local name="$1" h5ad_path="$2"
  local r_output_dir="${R_OUTPUT_BASE}/${name}_query_choir_cleanup_20260604"
  local cluster_csv="${r_output_dir}/reports/choir/choir_clusters.csv"
  local llm_screen_tsv="${r_output_dir}/reports/llm_choir_discovery_screen.tsv"
  local py_output_dir="${PY_OUTPUT_BASE}/${name}_query_cleaned_20260604"
  local output_h5ad="${py_output_dir}/${name}_query_cleaned.h5ad"
  local summary_json="${py_output_dir}/${name}_query_cleaned_summary.json"
  local removed_tsv="${py_output_dir}/${name}_query_cleaned_removed_cells.tsv"
  local log_file="${LOG_ROOT}/${name}_py_cleanup.log"

  echo "[$(date '+%F %T')] START ${name} Python cleanup" | tee "${log_file}"

  # Skip if already done
  if [[ -s "${output_h5ad}" ]] && [[ -s "${summary_json}" ]]; then
    echo "[$(date '+%F %T')] SKIP ${name}: Python cleanup already done" | tee -a "${log_file}"
    return 0
  fi

  # Check prerequisites
  if [[ ! -f "${cluster_csv}" ]]; then
    echo "[$(date '+%F %T')] FAIL ${name}: cluster CSV not found at ${cluster_csv}" | tee -a "${log_file}"
    return 5
  fi
  if [[ ! -f "${llm_screen_tsv}" ]]; then
    echo "[$(date '+%F %T')] FAIL ${name}: LLM screen TSV not found at ${llm_screen_tsv}" | tee -a "${log_file}"
    return 6
  fi

  mkdir -p "${py_output_dir}"

  "${PY}" "${PY_CLEANUP_SCRIPT}" \
    --input-h5ad "${h5ad_path}" \
    --cluster-csv "${cluster_csv}" \
    --llm-screen-tsv "${llm_screen_tsv}" \
    --output-h5ad "${output_h5ad}" \
    --summary-json "${summary_json}" \
    --removed-cells-tsv "${removed_tsv}" \
    --lineage-name "${name}" \
    --cluster-column "${CLUSTER_COLUMN}" \
    2>&1 | tee -a "${log_file}"
  local status=${PIPESTATUS[0]}

  echo "[$(date '+%F %T')] END ${name} Python cleanup status=${status}" | tee -a "${log_file}"
  return "${status}"
}

# ── Main execution ───────────────────────────────────────────────────────────

echo "=============================================="
echo "Query Clustering + OFA + LLM Cleanup Pipeline"
echo "Started: $(date '+%F %T')"
echo "=============================================="

declare -A JOB_PIDS
declare -A JOB_NAMES
FAILED_LINEAGES=()

# Phase 1: Run R clustering for all lineages (max 3 in parallel due to memory)
MAX_PARALLEL_R=3
RUNNING_COUNT=0
LINEAGE_INDEX=0

for lineage_entry in "${LINEAGES[@]}"; do
  IFS='|' read -r name category h5ad_path lineage_tag shared_engine <<< "${lineage_entry}"

  # Wait if we've reached parallel limit
  while [[ "${RUNNING_COUNT}" -ge "${MAX_PARALLEL_R}" ]]; do
    # Wait for any one to finish
    for pid in "${!JOB_PIDS[@]}"; do
      if ! kill -0 "${pid}" 2>/dev/null; then
        wait "${pid}" || true
        unset JOB_PIDS["${pid}"]
        RUNNING_COUNT=$((RUNNING_COUNT - 1))
        break
      fi
    done
    sleep 10
  done

  run_lineage_r "${name}" "${category}" "${h5ad_path}" "${lineage_tag}" "${shared_engine}" &
  local pid=$!
  JOB_PIDS["${pid}"]="${name}"
  JOB_NAMES["${pid}"]="${name}"
  RUNNING_COUNT=$((RUNNING_COUNT + 1))
  LINEAGE_INDEX=$((LINEAGE_INDEX + 1))
done

# Wait for remaining R jobs
echo "[$(date '+%F %T')] Waiting for remaining R jobs to finish..."
for pid in "${!JOB_PIDS[@]}"; do
  name="${JOB_NAMES[${pid}]}"
  if wait "${pid}"; then
    echo "[$(date '+%F %T')] R SUCCESS: ${name}"
  else
    local rc=$?
    echo "[$(date '+%F %T')] R FAILED (rc=${rc}): ${name} — will not run Python cleanup"
    FAILED_LINEAGES+=("${name}")
  fi
  unset JOB_PIDS["${pid}"]
done

echo ""
echo "=============================================="
echo "Phase 1 (R clustering+OFA+LLM) complete"
echo "Failed lineages: ${#FAILED_LINEAGES[@]}"
for f in "${FAILED_LINEAGES[@]}"; do echo "  - ${f}"; done
echo "=============================================="

# Phase 2: Python cleanup (sequential — fast, no parallel needed)
echo ""
echo "[$(date '+%F %T')] Starting Phase 2: Python cleanup"

for lineage_entry in "${LINEAGES[@]}"; do
  IFS='|' read -r name category h5ad_path lineage_tag shared_engine <<< "${lineage_entry}"

  # Skip if R phase failed for this lineage
  if printf '%s\n' "${FAILED_LINEAGES[@]}" | grep -qx "${name}"; then
    echo "[$(date '+%F %T')] SKIP ${name} Python cleanup: R phase failed"
    continue
  fi

  run_lineage_py_cleanup "${name}" "${h5ad_path}" || {
    echo "[$(date '+%F %T')] WARN ${name}: Python cleanup failed"
    FAILED_LINEAGES+=("${name}_py")
  }
done

echo ""
echo "=============================================="
echo "Pipeline complete: $(date '+%F %T')"
echo "Total lineages: ${#LINEAGES[@]}"
echo "Successful: $((${#LINEAGES[@]} - ${#FAILED_LINEAGES[@]}))"
echo "Failed: ${#FAILED_LINEAGES[@]}"
for f in "${FAILED_LINEAGES[@]}"; do echo "  - ${f}"; done
echo "=============================================="
echo ""
echo "Output locations:"
echo "  R outputs: ${R_OUTPUT_BASE}/"
echo "  Python outputs: ${PY_OUTPUT_BASE}/"
echo "  Logs: ${LOG_ROOT}/"

if [[ "${#FAILED_LINEAGES[@]}" -gt 0 ]]; then
  exit 1
fi
