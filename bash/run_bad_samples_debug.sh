#!/usr/bin/env bash
set -euo pipefail

BAD_TSV="${1:?need bad_samples.tsv}"
DATA_ROOT="${2:-/home/h2048/data/source/1210}"
DEBUG_R="${3:-/home/h2048/script/R/gsm_one_by_one_debug.R}"
OUT_BASE="${4:-/home/h2048/data/R/debug/bad_samples}"
LOG_TXT="${5:-${OUT_BASE}/all_bad_samples_debug.$(date +%Y%m%d_%H%M%S).txt}"

mkdir -p "$OUT_BASE"

{
  echo "============================================================"
  echo "Batch GSM debug run"
  echo "BAD_TSV:  $BAD_TSV"
  echo "DATA_ROOT:$DATA_ROOT"
  echo "OUT_BASE: $OUT_BASE"
  echo "LOG_TXT:  $LOG_TXT"
  echo "Time:     $(date)"
  echo "============================================================"
  echo
} > "$LOG_TXT"

# 跳过表头
tail -n +2 "$BAD_TSV" | while IFS=$'\t' read -r dataset_id sample_id format reason; do
  outdir="${OUT_BASE}/${dataset_id}_${sample_id}"
  mkdir -p "$outdir"

  {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "DATASET: $dataset_id  SAMPLE: $sample_id  FORMAT: $format  REASON: $reason"
    echo "OUTDIR : $outdir"
    echo "TIME   : $(date)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  } >> "$LOG_TXT"

  Rscript "$DEBUG_R" \
    --data_root "$DATA_ROOT" \
    --dataset_id "$dataset_id" \
    --sample_id "$sample_id" \
    --outdir "$outdir" \
    --run_emptydrops TRUE \
    --emptydrops_lower 100 \
    --emptydrops_fdr 0.01 \
    >> "$LOG_TXT" 2>&1 || {
      echo "[WARN] failed: ${dataset_id} ${sample_id}" >> "$LOG_TXT"
    }

  echo >> "$LOG_TXT"
done

echo "DONE. Combined log: $LOG_TXT"
