#!/bin/bash
# Fixed BBKNN Analysis Command
# Original issue: backslash had space after it, causing argument parsing error

cd /home/h2048/script/py
LOGDIR="/home/h2048/logs/$(date +%Y%m%d)"
mkdir -p "$LOGDIR"

# ⭐ FIXED: Removed space after backslash, added proper arguments
nohup python scanvi_bbknn_analysis_20251217_v1_6.py \
    --output_dir /home/h2048/data/py/1217 \
    --data_dir /home/h2048/data/core_data \
    --cell_types tcell bcell myeloid stromal_vascular \
    > "$LOGDIR"/scanvi_bbknn_analysis_20251217_v1_6.py_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Check process
sleep 1
ps aux | grep scanvi_bbknn_analysis_20251217_v1_6.py | grep -v grep

