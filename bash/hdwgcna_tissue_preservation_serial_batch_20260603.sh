#!/usr/bin/env bash
# Serial (one-at-a-time) hdWGCNA tissue preservation for all lineages
# Generated 2026-06-03 — avoids WGCNA parallel-worker memory contention
set -euo pipefail

RUNNER="/home/h2048/script/R/modules/hdwgcna/programs/hdwgcna_tissue_preservation_runner_20260531.R"
OUT_BASE="/home/h2048/data/R/20260603/hdwgcna_tissue_preservation_serial"
LOG_BASE="/home/h2048/logs/20260603"

mkdir -p "$OUT_BASE" "$LOG_BASE"

# ---- per-lineage config: rds|celltype_col|extra_args ----
declare -A LINEAGES
LINEAGES=(
  ["stromal_fibroblast"]="/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.rds|cell_type_L3|"
  ["bcell"]="/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/bcell_tissue_comparison_final.rds|cell_type_L3|"
  ["epithelial"]="/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/epithelial_tissue_comparison_final.rds|cell_type_L3_curated|"
  ["stromal_endothelial"]="/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508/stromal_endothelial_tissue_comparison_final.rds|cell_type_L3|"
  ["tnk"]="/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508/tnk_tissue_comparison_final.rds|cell_type_L3|--metacell-reduction scvi"
  ["myeloid"]="/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414/myeloid_tissue_comparison_final.rds|cell_type_L3|"
  ["stromal_smc"]="/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.rds|cell_type_L3|"
)

ORDER=(
  stromal_smc          # smallest, quick validation
  bcell
  stromal_endothelial
  tnk
  myeloid
  stromal_fibroblast
  epithelial           # largest, last
)

echo "=== hdWGCNA tissue preservation SERIAL batch ==="
echo "Started: $(date --iso-8601=seconds)"
echo "Runner:  $RUNNER"
echo ""

TOTAL=${#ORDER[@]}
IDX=0

for lineage in "${ORDER[@]}"; do
  IDX=$((IDX + 1))
  IFS='|' read -r rds_path ct_col extra_args <<< "${LINEAGES[$lineage]}"

  if [ ! -f "$rds_path" ]; then
    echo "[$IDX/$TOTAL SKIP] $lineage: RDS not found"
    continue
  fi

  out_dir="${OUT_BASE}/${lineage}"
  log_file="${LOG_BASE}/hdwgcna_tissue_preservation_serial_${lineage}_20260603.log"

  echo "[$IDX/$TOTAL] $lineage  @ $(date +%H:%M:%S)"
  echo "  RDS: $rds_path"
  echo "  CT:  $ct_col"
  echo "  EXTRA: ${extra_args:-none}"
  echo "  OUT: $out_dir"
  echo "  LOG: $log_file"

  start_ts=$(date +%s)

  /usr/bin/Rscript "$RUNNER" \
    --rds-path "$rds_path" \
    --lineage-name "$lineage" \
    --out-dir "$out_dir" \
    --celltype-col "$ct_col" \
    --tissue-col tissue \
    --sample-col sample \
    --max-shared 12 \
    --target-metacells 250 \
    --metacell-min-cells 25 \
    --soft-power-r2-cutoff 0.85 \
    --network-type 'signed hybrid' \
    --tom-type signed \
    --cor-type pearson \
    --deep-split 4 \
    --min-module-size 20 \
    --merge-cut-height 0.25 \
    --preservation-permutations 25 \
    ${extra_args:-} \
    >> "$log_file" 2>&1

  rc=$?
  end_ts=$(date +%s)
  elapsed=$((end_ts - start_ts))

  # quick summary
  summary="${out_dir}/hdwgcna_tissue_preservation_summary.json"
  if [ -f "$summary" ]; then
    py_out=$(python3 -c "
import json
d=json.load(open('$summary'))
print(f'networks={d[\"n_celltype_tissue_networks\"]} ok={d[\"n_ok_networks\"]} no_mod={d[\"n_no_module_networks\"]} skip={d[\"n_skipped_networks\"]} err={d[\"n_error_networks\"]}')
" 2>/dev/null)
    echo "  [$IDX/$TOTAL DONE] $lineage exit=$rc elapsed=${elapsed}s $py_out"
  else
    echo "  [$IDX/$TOTAL DONE] $lineage exit=$rc elapsed=${elapsed}s (no summary)"
  fi
  echo ""
done

echo "=== All lineages done @ $(date --iso-8601=seconds) ==="
echo "Summaries:"
for lineage in "${ORDER[@]}"; do
  summary="${OUT_BASE}/${lineage}/hdwgcna_tissue_preservation_summary.json"
  if [ -f "$summary" ]; then
    py_out=$(python3 -c "
import json
d=json.load(open('$summary'))
print(f'  {d[\"lineage_name\"]:25s} networks={d[\"n_celltype_tissue_networks\"]:3d} ok={d[\"n_ok_networks\"]:3d} skip={d[\"n_skipped_networks\"]:3d} err={d[\"n_error_networks\"]:3d}')
" 2>/dev/null)
    echo "$py_out"
  else
    printf "  %-25s (no summary)\n" "$lineage"
  fi
done
