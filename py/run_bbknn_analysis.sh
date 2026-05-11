#!/bin/bash
# BBKNN Analysis Runner Script
# Usage: 
#   ./run_bbknn_analysis.sh                    # Run all cell types
#   ./run_bbknn_analysis.sh tcell              # Run single cell type
#   ./run_bbknn_analysis.sh tcell bcell        # Run multiple cell types

cd /home/h2048/script/py || exit 1

# Create log directory
LOGDIR="/home/h2048/logs/$(date +%Y%m%d)"
mkdir -p "$LOGDIR"

# Default cell types if not provided
if [ $# -eq 0 ]; then
    CELL_TYPES=("tcell" "bcell" "myeloid" "stromal_vascular")
else
    CELL_TYPES=("$@")
fi

# Log file
LOG_FILE="$LOGDIR/scanvi_bbknn_analysis_20251217_v1_6.py_$(date +%Y%m%d_%H%M%S).log"

# Build command array (safer than string concatenation)
CMD_ARGS=(
    "scanvi_bbknn_analysis_20251217_v1_6.py"
    "--output_dir" "/home/h2048/data/py/1217"
    "--data_dir" "/home/h2048/data/core_data"
    "--cell_types"
)

# Add cell types to command
CMD_ARGS+=("${CELL_TYPES[@]}")

echo "=========================================="
echo "BBKNN Analysis Runner"
echo "=========================================="
echo "Cell types: ${CELL_TYPES[*]}"
echo "Output dir: /home/h2048/data/py/1217"
echo "Data dir: /home/h2048/data/core_data"
echo "Log file: $LOG_FILE"
echo "=========================================="

# Run in background
nohup python "${CMD_ARGS[@]}" > "$LOG_FILE" 2>&1 &
PID=$!

echo "Process started with PID: $PID"
echo "Monitor with: tail -f $LOG_FILE"
echo ""

# Check if process is running
sleep 2
if ps -p $PID > /dev/null 2>&1; then
    echo "✓ Process is running (PID: $PID)"
    echo ""
    ps aux | grep scanvi_bbknn_analysis_20251217_v1_6.py | grep -v grep
else
    echo "✗ Process failed to start. Check log: $LOG_FILE"
    echo ""
    tail -20 "$LOG_FILE"
    exit 1
fi

