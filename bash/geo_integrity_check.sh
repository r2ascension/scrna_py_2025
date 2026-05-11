#!/usr/bin/env bash
set -u
set -o pipefail

# ============================================================
# GEO/HRA local integrity checker
# - gzip CRC test
# - MTX nnz vs actual entries (detect truncation)
# - 10x consistency (barcodes/features vs mtx dims)
# - H5 basic open & structure check (python+h5py)
# ============================================================

DATA_ROOT="/home/h2048/data/source/1210"
OUTDIR=""
LEVEL="full"   # full|quick  (quick skips MTX full scan; still gzip -t)

# -------------- arg parse --------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data_root) DATA_ROOT="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --level) LEVEL="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

ts="$(date +%Y%m%d_%H%M%S)"
if [[ -z "$OUTDIR" ]]; then
  OUTDIR="${DATA_ROOT%/}/_integrity_${ts}"
fi
mkdir -p "$OUTDIR"

REPORT_CSV="$OUTDIR/integrity_report.csv"
BAD_TSV="$OUTDIR/bad_samples.tsv"
SUMMARY_TXT="$OUTDIR/summary.txt"

echo "dataset_id,sample_id,format,file,check,status,message" > "$REPORT_CSV"
echo -e "dataset_id\tsample_id\tformat\treason" > "$BAD_TSV"
: > "$SUMMARY_TXT"

# -------------- helpers --------------
csv_escape() {
  # escape double quotes for CSV
  local s="$1"
  s="${s//\"/\"\"}"
  printf '%s' "$s"
}

log_csv() {
  local dataset_id="$1" sample_id="$2" format="$3" file="$4" check="$5" status="$6" msg="$7"
  printf "%s,%s,%s,%s,%s,%s,\"%s\"\n" \
    "$(csv_escape "$dataset_id")" \
    "$(csv_escape "$sample_id")" \
    "$(csv_escape "$format")" \
    "$(csv_escape "$file")" \
    "$(csv_escape "$check")" \
    "$(csv_escape "$status")" \
    "$(csv_escape "$msg")" >> "$REPORT_CSV"
}

mark_bad() {
  local dataset_id="$1" sample_id="$2" format="$3" reason="$4"
  echo -e "${dataset_id}\t${sample_id}\t${format}\t${reason}" >> "$BAD_TSV"
}

is_gz() {
  [[ "$1" =~ \.gz$ ]]
}

gzip_test() {
  local f="$1" dataset_id="$2" sample_id="$3" format="$4"
  if [[ ! -f "$f" ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "gzip_t" "FAIL" "file_missing"
    return 1
  fi
  if is_gz "$f"; then
    if gzip -t "$f" >/dev/null 2>&1; then
      log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "gzip_t" "OK" "crc_ok"
      return 0
    else
      log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "gzip_t" "FAIL" "crc_or_trailer_bad"
      return 1
    fi
  else
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "gzip_t" "SKIP" "not_gz"
    return 0
  fi
}

mtx_dims_and_missing_full() {
  # prints: nrow ncol nnz actual missing
  local mtx="$1"
  local cmd="cat"
  if is_gz "$mtx"; then cmd="zcat"; fi

  # full scan (accurate)
  $cmd "$mtx" | awk '
    BEGIN{seen=0; nnz=0; actual=0; nrow=0; ncol=0}
    /^%/ {next}
    {
      seen++
      if (seen==1) { nrow=$1; ncol=$2; nnz=$3; next }
      actual++
    }
    END{
      missing = nnz-actual
      printf "%d %d %d %d %d\n", nrow, ncol, nnz, actual, missing
    }'
}

mtx_dims_quick() {
  # quick: only read dims line, do not scan entries
  local mtx="$1"
  local cmd="cat"
  if is_gz "$mtx"; then cmd="zcat"; fi

  $cmd "$mtx" | awk '
    BEGIN{seen=0}
    /^%/ {next}
    { print $1, $2, $3; exit }'
}

count_lines_gz_or_plain() {
  local f="$1"
  if is_gz "$f"; then
    zcat "$f" | wc -l
  else
    wc -l < "$f"
  fi
}

check_mtx_file() {
  local dataset_id="$1" sample_id="$2" format="$3" mtx="$4"

  gzip_test "$mtx" "$dataset_id" "$sample_id" "$format" || {
    mark_bad "$dataset_id" "$sample_id" "$format" "mtx_gz_crc_fail"
    return 1
  }

  if [[ "$LEVEL" == "quick" ]]; then
    local dims
    dims="$(mtx_dims_quick "$mtx" 2>/dev/null || true)"
    if [[ -z "$dims" ]]; then
      log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$mtx")" "mtx_dims" "FAIL" "cannot_read_dims"
      mark_bad "$dataset_id" "$sample_id" "$format" "mtx_dims_unreadable"
      return 1
    fi
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$mtx")" "mtx_dims" "OK" "dims=${dims} (quick)"
    return 0
  fi

  local out
  out="$(mtx_dims_and_missing_full "$mtx" 2>/dev/null || true)"
  if [[ -z "$out" ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$mtx")" "mtx_scan" "FAIL" "cannot_scan_mtx"
    mark_bad "$dataset_id" "$sample_id" "$format" "mtx_scan_failed"
    return 1
  fi

  local nrow ncol nnz actual missing
  read -r nrow ncol nnz actual missing <<< "$out"

  if [[ "$missing" -ne 0 ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$mtx")" "mtx_nnz_vs_entries" "FAIL" \
      "dims=${nrow}x${ncol} nnz=${nnz} actual=${actual} missing=${missing} (TRUNCATED_OR_CORRUPT)"
    mark_bad "$dataset_id" "$sample_id" "$format" "mtx_truncated_missing_${missing}"
    return 1
  else
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$mtx")" "mtx_nnz_vs_entries" "OK" \
      "dims=${nrow}x${ncol} nnz=${nnz} actual=${actual} missing=0"
    return 0
  fi
}

check_10x_dir() {
  local dataset_id="$1" sample_id="$2" sample_dir="$3"
  local format="10x"

  # locate files (prefer .gz)
  local bar features mtx
  bar="$(ls -1 "$sample_dir"/*barcodes.tsv* 2>/dev/null | head -n 1 || true)"
  features="$(ls -1 "$sample_dir"/*features.tsv* "$sample_dir"/*genes.tsv* 2>/dev/null | head -n 1 || true)"
  mtx="$(ls -1 "$sample_dir"/*matrix.mtx* 2>/dev/null | head -n 1 || true)"

  if [[ -z "$bar" || -z "$features" || -z "$mtx" ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$sample_dir")" "10x_detect" "FAIL" \
      "missing_one_of(barcodes/features/matrix)"
    mark_bad "$dataset_id" "$sample_id" "$format" "10x_missing_files"
    return 1
  fi

  gzip_test "$bar" "$dataset_id" "$sample_id" "$format" || mark_bad "$dataset_id" "$sample_id" "$format" "barcodes_gz_crc_fail"
  gzip_test "$features" "$dataset_id" "$sample_id" "$format" || mark_bad "$dataset_id" "$sample_id" "$format" "features_gz_crc_fail"

  # MTX integrity
  check_mtx_file "$dataset_id" "$sample_id" "$format" "$mtx" || true

  # get dims (quick read ok even in full mode; we already scanned in full)
  local dims
  dims="$(mtx_dims_quick "$mtx" 2>/dev/null || true)"
  if [[ -z "$dims" ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$mtx")" "mtx_dims" "FAIL" "cannot_read_dims"
    mark_bad "$dataset_id" "$sample_id" "$format" "mtx_dims_unreadable"
    return 1
  fi
  local nrow ncol nnz
  read -r nrow ncol nnz <<< "$dims"

  # line counts
  local nbar nfeat
  nbar="$(count_lines_gz_or_plain "$bar" 2>/dev/null || echo -1)"
  nfeat="$(count_lines_gz_or_plain "$features" 2>/dev/null || echo -1)"

  if [[ "$nbar" -ne "$ncol" ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$bar")" "barcodes_n_vs_mtx_ncol" "FAIL" \
      "barcodes_lines=${nbar} mtx_ncol=${ncol}"
    mark_bad "$dataset_id" "$sample_id" "$format" "barcodes_count_mismatch"
  else
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$bar")" "barcodes_n_vs_mtx_ncol" "OK" \
      "barcodes_lines=${nbar} mtx_ncol=${ncol}"
  fi

  if [[ "$nfeat" -ne "$nrow" ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$features")" "features_n_vs_mtx_nrow" "FAIL" \
      "features_lines=${nfeat} mtx_nrow=${nrow}"
    mark_bad "$dataset_id" "$sample_id" "$format" "features_count_mismatch"
  else
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$features")" "features_n_vs_mtx_nrow" "OK" \
      "features_lines=${nfeat} mtx_nrow=${nrow}"
  fi

  return 0
}

check_text_matrix_file() {
  # for .csv/.tsv/.txt (optionally .gz)
  local dataset_id="$1" sample_id="$2" f="$3"
  local format="matrix"

  gzip_test "$f" "$dataset_id" "$sample_id" "$format" || {
    mark_bad "$dataset_id" "$sample_id" "$format" "matrixfile_gz_crc_fail"
    return 1
  }

  # sample a few lines for delimiter + column consistency
  local headn=200
  local tmp="$OUTDIR/.tmp_${sample_id}_head.txt"
  if is_gz "$f"; then
    zcat "$f" | head -n "$headn" > "$tmp" 2>/dev/null || true
  else
    head -n "$headn" "$f" > "$tmp" 2>/dev/null || true
  fi

  if [[ ! -s "$tmp" ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "matrix_sample_read" "FAIL" "cannot_read_head"
    mark_bad "$dataset_id" "$sample_id" "$format" "matrix_head_unreadable"
    rm -f "$tmp"
    return 1
  fi

  # detect delimiter by first line
  local first
  first="$(head -n 1 "$tmp")"
  local delim=","
  if [[ "$first" == *$'\t'* ]]; then delim=$'\t'; fi

  # column count consistency (head only; gzip -t already ensures full stream is intact)
  local nf0
  nf0="$(awk -F"$delim" 'NR==1{print NF; exit}' "$tmp")"
  local bad
  bad="$(awk -F"$delim" -v nf0="$nf0" 'NR>1 && NF!=nf0 {c++} END{print c+0}' "$tmp")"

  if [[ "$bad" -gt 0 ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "matrix_col_consistency(head)" "FAIL" \
      "head_nf0=${nf0} bad_lines_in_head=${bad}"
    mark_bad "$dataset_id" "$sample_id" "$format" "matrix_cols_inconsistent_head"
  else
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "matrix_col_consistency(head)" "OK" \
      "head_nf0=${nf0} bad_lines_in_head=0"
  fi

  rm -f "$tmp"
  return 0
}

check_h5_file() {
  local dataset_id="$1" sample_id="$2" f="$3"
  local format="h5"

  gzip_test "$f" "$dataset_id" "$sample_id" "$format" >/dev/null 2>&1 || true  # usually not gz

  # python+h5py checks
  python3 - <<PY 2>/dev/null
import sys
import os
f = r"""$f"""
try:
    import h5py
except Exception as e:
    print("NO_H5PY")
    sys.exit(3)

try:
    with h5py.File(f, "r") as h:
        # try 10x v3 structure first
        if "matrix" in h:
            g = h["matrix"]
            # expected arrays
            required = ["data", "indices", "indptr", "shape"]
            miss = [x for x in required if x not in g]
            if miss:
                print("OPEN_OK_BUT_MISSING_KEYS:" + ",".join(miss))
                sys.exit(2)
            data = g["data"]
            indices = g["indices"]
            indptr = g["indptr"]
            shape = tuple(g["shape"][()])
            # consistency
            nnz = int(data.shape[0])
            ok = (nnz == int(indices.shape[0])) and (nnz == int(indptr[-1]))
            if not ok:
                print(f"INCONSISTENT: nnz(data)={nnz} nnz(indices)={int(indices.shape[0])} indptr_last={int(indptr[-1])} shape={shape}")
                sys.exit(2)
            # check NaN/Inf
            import numpy as np
            arr = data[()]
            if np.isnan(arr).any() or np.isinf(arr).any():
                print("DATA_HAS_NAN_OR_INF")
                sys.exit(2)
            print(f"OK: shape={shape} nnz={nnz}")
            sys.exit(0)
        else:
            # not 10x v3; just list top-level keys
            keys = list(h.keys())
            print("OPEN_OK_KEYS:" + ",".join(keys[:50]))
            sys.exit(0)
except Exception as e:
    print("OPEN_FAIL:" + str(e))
    sys.exit(1)
PY

  local rc=$?
  if [[ $rc -eq 0 ]]; then
    local msg
    msg="$(python3 - <<PY 2>/dev/null
import h5py
f=r"""$f"""
with h5py.File(f,"r") as h:
    if "matrix" in h and all(k in h["matrix"] for k in ["shape","data","indices","indptr"]):
        shape=tuple(h["matrix"]["shape"][()])
        nnz=int(h["matrix"]["data"].shape[0])
        print(f"shape={shape} nnz={nnz}")
    else:
        print("opened_ok")
PY
)"
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "h5_open_check" "OK" "$msg"
    return 0
  elif [[ $rc -eq 3 ]]; then
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "h5_open_check" "SKIP" "python_h5py_not_available"
    return 0
  else
    local msg
    msg="$(python3 - <<PY 2>/dev/null
print("see_stderr_or_python_result")
PY
)"
    log_csv "$dataset_id" "$sample_id" "$format" "$(basename "$f")" "h5_open_check" "FAIL" "open_or_structure_fail"
    mark_bad "$dataset_id" "$sample_id" "$format" "h5_open_or_structure_fail"
    return 1
  fi
}

detect_format_and_check() {
  local dataset_id="$1" sample_id="$2" sample_dir="$3"

  # h5
  local h5
  h5="$(ls -1 "$sample_dir"/*.h5 "$sample_dir"/*.hdf5 2>/dev/null | head -n 1 || true)"
  if [[ -n "$h5" ]]; then
    check_h5_file "$dataset_id" "$sample_id" "$h5" || true
    return 0
  fi

  # 10x
  if ls "$sample_dir"/*matrix.mtx* >/dev/null 2>&1 && \
     ls "$sample_dir"/*barcodes.tsv* >/dev/null 2>&1 && \
     (ls "$sample_dir"/*features.tsv* >/dev/null 2>&1 || ls "$sample_dir"/*genes.tsv* >/dev/null 2>&1); then
    check_10x_dir "$dataset_id" "$sample_id" "$sample_dir" || true
    return 0
  fi

  # matrix text/csv
  local mf
  mf="$(ls -1 "$sample_dir"/*.csv "$sample_dir"/*.csv.gz "$sample_dir"/*.tsv "$sample_dir"/*.tsv.gz "$sample_dir"/*.txt "$sample_dir"/*.txt.gz 2>/dev/null | head -n 1 || true)"
  if [[ -n "$mf" ]]; then
    check_text_matrix_file "$dataset_id" "$sample_id" "$mf" || true
    return 0
  fi

  log_csv "$dataset_id" "$sample_id" "unknown" "$(basename "$sample_dir")" "detect_format" "SKIP" "no_supported_files_found"
  mark_bad "$dataset_id" "$sample_id" "unknown" "unknown_format"
  return 0
}

# -------------- main loop --------------
echo "DATA_ROOT=$DATA_ROOT" >> "$SUMMARY_TXT"
echo "OUTDIR=$OUTDIR" >> "$SUMMARY_TXT"
echo "LEVEL=$LEVEL" >> "$SUMMARY_TXT"
echo "" >> "$SUMMARY_TXT"

datasets=()
for d in "$DATA_ROOT"/GSE* "$DATA_ROOT"/HRA*; do
  [[ -d "$d" ]] || continue
  datasets+=("$d")
done

echo "Found datasets: ${#datasets[@]}" >> "$SUMMARY_TXT"

for dataset_dir in "${datasets[@]}"; do
  dataset_id="$(basename "$dataset_dir")"

  # samples
  for sample_dir in "$dataset_dir"/GSM* "$dataset_dir"/HRR*; do
    [[ -d "$sample_dir" ]] || continue
    sample_id="$(basename "$sample_dir")"
    detect_format_and_check "$dataset_id" "$sample_id" "$sample_dir"
  done
done

# -------------- summary --------------
total_rows=$(( $(wc -l < "$REPORT_CSV") - 1 ))
bad_n=$(( $(wc -l < "$BAD_TSV") - 1 ))

echo "" >> "$SUMMARY_TXT"
echo "Report rows (excluding header): $total_rows" >> "$SUMMARY_TXT"
echo "Bad samples: $bad_n" >> "$SUMMARY_TXT"

echo "Done."
echo "Report: $REPORT_CSV"
echo "Bad samples: $BAD_TSV"
echo "Summary: $SUMMARY_TXT"
