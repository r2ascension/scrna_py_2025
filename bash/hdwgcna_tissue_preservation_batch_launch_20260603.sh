#!/usr/bin/env bash
# Batch launch hdWGCNA tissue preservation for all lineages
# Generated 2026-06-03
set -euo pipefail

RUNNER="/home/h2048/script/R/modules/hdwgcna/programs/hdwgcna_tissue_preservation_runner_20260531.R"
OUT_BASE="/home/h2048/data/R/20260603/hdwgcna_tissue_preservation"
LOG_BASE="/home/h2048/logs/20260603"

COMMON_ARGS_STR="--max-shared 12 --target-metacells 250 --metacell-min-cells 25 --soft-power-r2-cutoff 0.85 --network-type 'signed hybrid' --tom-type signed --cor-type pearson --deep-split 4 --min-module-size 20 --merge-cut-height 0.25 --preservation-permutations 25"

declare -A LINEAGES
LINEAGES=(
  ["stromal_fibroblast"]="/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.rds|cell_type_L3|tissue|sample"
  ["bcell"]="/home/h2048/data/R/0508/bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508/bcell_tissue_comparison_final.rds|cell_type_L3|tissue|sample"
  ["epithelial"]="/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/epithelial_tissue_comparison_final.rds|cell_type_L3_curated|tissue|sample"
  ["stromal_endothelial"]="/home/h2048/data/R/0508/stromal_endothelial_tissue_comparison_v1_1_2_rm_choir6_52_20260508/stromal_endothelial_tissue_comparison_final.rds|cell_type_L3|tissue|sample"
  ["tnk"]="/home/h2048/data/R/0508/tnk_tissue_comparison_v2_6_4_rm_choir23_28_31_41_ofa41_66_20260508/tnk_tissue_comparison_final.rds|cell_type_L3|tissue|sample"
  ["myeloid"]="/home/h2048/data/R/0414/myeloid_tissue_comparison_v1_2_2_20260414/myeloid_tissue_comparison_final.rds|cell_type_L3|tissue|sample"
  ["stromal_smc"]="/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.rds|cell_type_L3|tissue|sample"
)

echo "=== hdWGCNA tissue preservation batch launch ==="
echo "Started: $(date --iso-8601=seconds)"
echo "Runner:  $RUNNER"
echo ""

for lineage in "${!LINEAGES[@]}"; do
  IFS='|' read -r rds_path ct_col ts_col sm_col <<< "${LINEAGES[$lineage]}"

  if [ ! -f "$rds_path" ]; then
    echo "[SKIP] $lineage: RDS not found at $rds_path"
    continue
  fi

  out_dir="${OUT_BASE}/${lineage}"
  log_file="${LOG_BASE}/hdwgcna_tissue_preservation_${lineage}_20260603.log"

  echo "[LAUNCH] $lineage"
  echo "  RDS:     $rds_path"
  echo "  CT col:  $ct_col"
  echo "  Out:     $out_dir"
  echo "  Log:     $log_file"

  nohup bash -lc "
    echo '[$(date --iso-8601=seconds)] START $lineage'
    echo \"[env] host=\$(hostname) pid=\$\$\"
    echo \"[env] rds=$rds_path\"
    echo \"[env] ct_col=$ct_col ts_col=$ts_col sm_col=$sm_col\"
    echo \"[env] out_dir=$out_dir\"
    /usr/bin/Rscript '$RUNNER' \
      --rds-path '$rds_path' \
      --lineage-name '$lineage' \
      --out-dir '$out_dir' \
      --celltype-col '$ct_col' \
      --tissue-col '$ts_col' \
      --sample-col '$sm_col' \
      ${COMMON_ARGS_STR}
    rc=\$?
    echo '[$(date --iso-8601=seconds)] DONE $lineage exit=\$rc'
  " > "$log_file" 2>&1 < /dev/null &

  echo "  PID: $!"
  echo ""
done

echo "=== All lineages launched ==="
echo "Check status:"
for lineage in "${!LINEAGES[@]}"; do
  echo "  tail -5 ${LOG_BASE}/hdwgcna_tissue_preservation_${lineage}_20260603.log"
done
