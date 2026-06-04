#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/h2048"
DATE_TAG="${DATE_TAG:-20260523}"
BASE_OUT="${BASE_OUT:-/home/h2048/data/py/0523}"
ENV_PREFIX="${ENV_PREFIX:-/home/h2048/miniconda3/envs/scarches_pertpy_scmodels_20260511}"
PYTHON_BIN="${PYTHON_BIN:-python}"

SUBSET_DIR="${SUBSET_DIR:-$BASE_OUT/scarches_mapping_FIXED_v1_2_subsets_nocovid_${DATE_TAG}}"

EPITHELIAL_SCARCHES_OUT="${EPITHELIAL_SCARCHES_OUT:-$BASE_OUT/epithelial_scarches_query_v1_2_nocovid_${DATE_TAG}}"
EPITHELIAL_MERGE_OUT="${EPITHELIAL_MERGE_OUT:-$BASE_OUT/epithelial_merge_input_v1_0_nocovid_${DATE_TAG}}"
EPITHELIAL_SCHPL_OUT="${EPITHELIAL_SCHPL_OUT:-$BASE_OUT/epithelial_merge_schpl_v1_0_nocovid_${DATE_TAG}}"
EPITHELIAL_EVAL_OUT="${EPITHELIAL_EVAL_OUT:-$BASE_OUT/epithelial_nocovid_eval_${DATE_TAG}}"

BCELL_MAP_OUT="${BCELL_MAP_OUT:-$BASE_OUT/bcell_scarches_query_nocovid_${DATE_TAG}}"
BCELL_SCHPL_OUT="${BCELL_SCHPL_OUT:-$BASE_OUT/bcell_merge_schpl_nocovid_${DATE_TAG}}"

TCELL_OUT="${TCELL_OUT:-$BASE_OUT/tcell_only_merged_nocovid_${DATE_TAG}}"
TCELL_SCHPL_OUT="${TCELL_SCHPL_OUT:-$BASE_OUT/tcell_merge_schpl_nocovid_${DATE_TAG}}"

MYELOID_OUT="${MYELOID_OUT:-$BASE_OUT/myeloid_only_merged_nocovid_${DATE_TAG}}"
MYELOID_SCHPL_OUT="${MYELOID_SCHPL_OUT:-$BASE_OUT/myeloid_merge_schpl_nocovid_${DATE_TAG}}"

STROMAL_QUERY_OUT="${STROMAL_QUERY_OUT:-$BASE_OUT/stromal_scarches_query_branchwise_nocovid_${DATE_TAG}}"
STROMAL_SCHPL_OUT="${STROMAL_SCHPL_OUT:-$BASE_OUT/stromal_merge_schpl_nocovid_${DATE_TAG}}"
ALL_LINEAGE_EVAL_OUT="${ALL_LINEAGE_EVAL_OUT:-$BASE_OUT/all_lineage_eval_${DATE_TAG}}"

RUN_SUBSET="${RUN_SUBSET:-1}"
RUN_EPITHELIAL="${RUN_EPITHELIAL:-1}"
RUN_BCELL="${RUN_BCELL:-1}"
RUN_TCELL="${RUN_TCELL:-1}"
RUN_MYELOID="${RUN_MYELOID:-1}"
RUN_STROMAL="${RUN_STROMAL:-1}"
RUN_EVAL_AT_END="${RUN_EVAL_AT_END:-1}"
ALLOW_OVERWRITE="${ALLOW_OVERWRITE:-0}"

export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1
unset PYTHONPATH || true
export MPLBACKEND=Agg

section() {
  printf '\n%s\n' "================================================================================"
  printf '%s\n' "$1"
  printf '%s\n' "================================================================================"
}

require_file() {
  local file_path="$1"
  local description="$2"
  if [[ ! -f "$file_path" ]]; then
    printf '[ERROR] Missing %s: %s\n' "$description" "$file_path" >&2
    exit 1
  fi
}

if [[ -n "${ENV_PREFIX}" ]]; then
  source "$PROJECT_ROOT/miniconda3/etc/profile.d/conda.sh"
  conda activate "$ENV_PREFIX"
  PYTHON_BIN="${PYTHON_BIN:-python}"
fi

overwrite_flag=()
if [[ "$ALLOW_OVERWRITE" == "1" ]]; then
  overwrite_flag+=(--allow-overwrite)
fi

if [[ "$RUN_SUBSET" == "1" ]]; then
  section "[1/6] Build no-COVID subsets"
  CELLTYPE_SUBSET_INPUT_H5AD="$PROJECT_ROOT/data/py/0127/scarches_mapping_FIXED_v1_2/query_mapped_to_reference.h5ad" \
  CELLTYPE_SUBSET_OUTPUT_DIR="$SUBSET_DIR" \
  CELLTYPE_SUBSET_FILTER_COVID_RELATED=1 \
  CELLTYPE_SUBSET_COVID_PRIMARY_COLS="disease" \
  CELLTYPE_SUBSET_COVID_SECONDARY_COLS="COVID_status" \
  CELLTYPE_SUBSET_COVID_DATASET_COLS="dataset" \
  CELLTYPE_SUBSET_ALLOW_DATASET_FALLBACK=0 \
  "$PYTHON_BIN" "$PROJECT_ROOT/script/py/subset_by_celltype_20260127_v1_0.py"
fi

if [[ "$RUN_EPITHELIAL" == "1" ]]; then
  section "[2/6] Run epithelial no-COVID scArches + merge + scHPL + evaluation"
  QUERY_H5AD="$SUBSET_DIR/epithelial_cells.h5ad" \
  SCARCHES_OUT_DIR="$EPITHELIAL_SCARCHES_OUT" \
  MERGE_OUT_DIR="$EPITHELIAL_MERGE_OUT" \
  SCHPL_OUT_DIR="$EPITHELIAL_SCHPL_OUT" \
  EVAL_OUT_DIR="$EPITHELIAL_EVAL_OUT" \
  "$PROJECT_ROOT/script/bash/run_epithelial_nocovid_scarches_schpl_20260523.sh"
fi

if [[ "$RUN_BCELL" == "1" ]]; then
  section "[3/6] Run B-cell no-COVID query mapping + scHPL"
  BCELL_QUERY_H5AD="$SUBSET_DIR/b_cells.h5ad" \
  BCELL_OUTPUT_DIR="$BCELL_MAP_OUT" \
  "$PYTHON_BIN" "$PROJECT_ROOT/script/py/step2_map_query_20260204_v2_5_4.py"

  "$PYTHON_BIN" "$PROJECT_ROOT/script/py/lineage_merge_schpl_cli_20260523.py" \
    --lineage-name "B cell" \
    --lineage-slug "bcell" \
    --version "v1_0_nocovid_${DATE_TAG}" \
    --input-h5ad "$BCELL_MAP_OUT/reference_plus_query_merged_L2.h5ad" \
    --output-dir "$BCELL_SCHPL_OUT" \
    --reference-label-key "Cell_Type_L2" \
    --query-pred-key "Cell_Type_L2_pred" \
    --query-final-key "Cell_Type_L2_final" \
    --query-conf-key "mapping_confidence" \
    --latent-key "X_scANVI_L2" \
    --umap-key "X_umap" \
    --marker-layer-preferred "log1p" \
    --marker-genes CD19 MS4A1 CD27 IGHD IGHM MZB1 SDC1 JCHAIN \
    "${overwrite_flag[@]}"
fi

if [[ "$RUN_TCELL" == "1" ]]; then
  section "[4/6] Run T/NK no-COVID merge + scHPL"
  mkdir -p "$TCELL_OUT" "$TCELL_SCHPL_OUT"
  tcell_output_prefix="tcell_only_merged_nocovid"
  tcell_merge_log="$TCELL_OUT/step4_tnk_merge.log"
  {
    TCELL_QUERY_H5AD="$SUBSET_DIR/t_cells.h5ad" \
    TCELL_OUTPUT_DIR="$TCELL_OUT" \
    TCELL_OUTPUT_PREFIX="$tcell_output_prefix" \
    "$PYTHON_BIN" "$PROJECT_ROOT/script/py/tcell_only_ref_query_merge_pipeline_20260225_v1_1.py"
  } 2>&1 | tee "$tcell_merge_log"

  tcell_input_h5ad="$TCELL_OUT/${tcell_output_prefix}_results.h5ad"
  require_file "$tcell_input_h5ad" "T/NK merge output H5AD"

  tcell_schpl_log="$TCELL_SCHPL_OUT/step4_tnk_schpl.log"
  {
    "$PYTHON_BIN" "$PROJECT_ROOT/script/py/lineage_merge_schpl_cli_20260523.py" \
      --lineage-name "T/NK" \
      --lineage-slug "tnk" \
      --version "v1_0_nocovid_${DATE_TAG}" \
      --input-h5ad "$tcell_input_h5ad" \
      --output-dir "$TCELL_SCHPL_OUT" \
      --reference-label-key "cell_type_fine" \
      --query-pred-key "scanvi_pred" \
      --query-final-key "scanvi_pred" \
      --query-conf-key "scanvi_confidence" \
      --latent-key "X_scANVI" \
      --umap-key "X_umap" \
      --marker-layer-preferred "log1p" \
      --marker-genes CD3D TRBC2 CD4 IL7R CCR7 LTB CD8A NKG7 \
      "${overwrite_flag[@]}"
  } 2>&1 | tee "$tcell_schpl_log"
fi

if [[ "$RUN_MYELOID" == "1" ]]; then
  section "[5/6] Run myeloid no-COVID merge + scHPL"
  mkdir -p "$MYELOID_OUT" "$MYELOID_SCHPL_OUT"
  myeloid_output_prefix="myeloid_only_merged_nocovid"
  myeloid_merge_log="$MYELOID_OUT/step5_myeloid_merge.log"
  {
    MYELOID_QUERY_H5AD="$SUBSET_DIR/myeloid_cells.h5ad" \
    MYELOID_OUTPUT_DIR="$MYELOID_OUT" \
    MYELOID_OUTPUT_PREFIX="$myeloid_output_prefix" \
    "$PYTHON_BIN" "$PROJECT_ROOT/script/py/myeloid_only_ref_query_merge_pipeline_20260225_v1_1.py"
  } 2>&1 | tee "$myeloid_merge_log"

  myeloid_input_h5ad="$MYELOID_OUT/${myeloid_output_prefix}_results.h5ad"
  require_file "$myeloid_input_h5ad" "Myeloid merge output H5AD"

  myeloid_schpl_log="$MYELOID_SCHPL_OUT/step5_myeloid_schpl.log"
  {
    "$PYTHON_BIN" "$PROJECT_ROOT/script/py/lineage_merge_schpl_cli_20260523.py" \
      --lineage-name "Myeloid" \
      --lineage-slug "myeloid" \
      --version "v1_0_nocovid_${DATE_TAG}" \
      --input-h5ad "$myeloid_input_h5ad" \
      --output-dir "$MYELOID_SCHPL_OUT" \
      --reference-label-key "cell_type_fine" \
      --query-pred-key "scanvi_pred" \
      --query-final-key "scanvi_pred" \
      --query-conf-key "scanvi_confidence" \
      --latent-key "X_scANVI" \
      --umap-key "X_umap" \
      --marker-layer-preferred "log1p" \
      --marker-genes LYZ FCER1G CD14 FCGR3A CLEC10A XCR1 FCGR3B TPSB2 \
      "${overwrite_flag[@]}"
  } 2>&1 | tee "$myeloid_schpl_log"
fi

if [[ "$RUN_STROMAL" == "1" ]]; then
  section "[6/6] Run stromal branchwise no-COVID query mapping + scHPL"
  mkdir -p "$STROMAL_QUERY_OUT" "$STROMAL_SCHPL_OUT"
  stromal_query_log="$STROMAL_QUERY_OUT/step6_stromal_query.log"
  {
    "$PYTHON_BIN" "$PROJECT_ROOT/script/py/stromal_scarches_query_branchwise_20260407_v1_0.py" \
      --query-dir "$SUBSET_DIR" \
      --output-dir "$STROMAL_QUERY_OUT"
  } 2>&1 | tee "$stromal_query_log"

  stromal_manifest="$STROMAL_QUERY_OUT/branch_query_manifest.json"
  require_file "$stromal_manifest" "stromal branch query manifest"

  stromal_schpl_log="$STROMAL_SCHPL_OUT/step6_stromal_schpl.log"
  {
    "$PYTHON_BIN" "$PROJECT_ROOT/script/py/stromal_branchwise_merge_schpl_20260523_v1.py" \
      --query-manifest "$stromal_manifest" \
      --output-dir "$STROMAL_SCHPL_OUT" \
      "${overwrite_flag[@]}"
  } 2>&1 | tee "$stromal_schpl_log"
fi

if [[ "$RUN_EVAL_AT_END" == "1" ]]; then
  section "[7/7] Build all-lineage evaluation summary"
  mkdir -p "$ALL_LINEAGE_EVAL_OUT"
  eval_log="$ALL_LINEAGE_EVAL_OUT/step7_all_lineage_eval.log"
  {
    "$PYTHON_BIN" "$PROJECT_ROOT/script/py/evaluate_all_lineages_nocovid_scarches_schpl_20260523.py" \
      --base-out "$BASE_OUT" \
      --output-dir "$ALL_LINEAGE_EVAL_OUT"
  } 2>&1 | tee "$eval_log"
fi

section "Run summary"
printf 'Subsets          : %s\n' "$SUBSET_DIR"
printf 'Epithelial eval  : %s\n' "$EPITHELIAL_EVAL_OUT"
printf 'B-cell schpl     : %s\n' "$BCELL_SCHPL_OUT"
printf 'T/NK schpl       : %s\n' "$TCELL_SCHPL_OUT"
printf 'Myeloid schpl    : %s\n' "$MYELOID_SCHPL_OUT"
printf 'Stromal schpl    : %s\n' "$STROMAL_SCHPL_OUT"
printf 'All-lineage eval : %s\n' "$ALL_LINEAGE_EVAL_OUT"
printf '\n[OK] all-lineage no-COVID pipeline script finished\n'
