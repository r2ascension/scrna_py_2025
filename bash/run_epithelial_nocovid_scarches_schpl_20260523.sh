#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/home/h2048"
ENV_PREFIX="${ENV_PREFIX:-/home/h2048/miniconda3/envs/scarches_pertpy_scmodels_20260511}"
QUERY_H5AD="${QUERY_H5AD:-/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/epithelial_cells.h5ad}"
SCARCHES_OUT_DIR="${SCARCHES_OUT_DIR:-/home/h2048/data/py/0523/epithelial_scarches_query_v1_2_nocovid_20260523}"
MERGE_OUT_DIR="${MERGE_OUT_DIR:-/home/h2048/data/py/0523/epithelial_merge_input_v1_0_nocovid_20260523}"
SCHPL_OUT_DIR="${SCHPL_OUT_DIR:-/home/h2048/data/py/0523/epithelial_merge_schpl_v1_0_nocovid_20260523}"
EVAL_OUT_DIR="${EVAL_OUT_DIR:-/home/h2048/data/py/0523/epithelial_nocovid_eval_20260523}"

source "$PROJECT_ROOT/miniconda3/etc/profile.d/conda.sh"
conda activate "$ENV_PREFIX"

export EPITHELIAL_QUERY_H5AD="$QUERY_H5AD"
export EPITHELIAL_SCARCHES_OUT_DIR="$SCARCHES_OUT_DIR"
export EPITHELIAL_FILTER_COVID_RELATED=1
export EPITHELIAL_FILTER_PRIMARY_COLS="disease"
export EPITHELIAL_FILTER_SECONDARY_COLS="COVID_status"
export EPITHELIAL_FILTER_INCLUDE_DATASET_FALLBACK=0

python "$PROJECT_ROOT/script/py/epithelial_scarches_query_mapping_20260315_v1_1.py"

python "$PROJECT_ROOT/script/py/epithelial_merge_ref_query_20260523_v1.py" \
  --query-h5ad "$SCARCHES_OUT_DIR/epithelial_query_mapped_v1_2.h5ad" \
  --output-dir "$MERGE_OUT_DIR"

python "$PROJECT_ROOT/script/py/epithelial_merge_schpl_20260523_v1.py" \
  --input-h5ad "$MERGE_OUT_DIR/epithelial_reference_plus_query_merged_nocovid_20260523.h5ad" \
  --output-dir "$SCHPL_OUT_DIR"

python "$PROJECT_ROOT/script/py/evaluate_epithelial_nocovid_scarches_schpl_20260523.py" \
  --mapped-query-h5ad "$SCARCHES_OUT_DIR/epithelial_query_mapped_v1_2.h5ad" \
  --schpl-h5ad "$SCHPL_OUT_DIR/epithelial_reference_plus_query_schpl_v1_0_nocovid_20260523.h5ad" \
  --output-dir "$EVAL_OUT_DIR"

printf '\n[OK] epithelial no-COVID pipeline finished\n'
