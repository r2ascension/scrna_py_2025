#!/bin/bash
# scANVI Downstream Pipeline Runner Script
# Usage: ./run_downstream_pipeline.sh

cd /home/h2048/script/py
LOGDIR="/home/h2048/logs/$(date +%Y%m%d)"
mkdir -p "$LOGDIR"

# Run v1.6.2 pipeline
nohup python -u scanvi_downstream_pipeline_20251218_v1_6_2.py > "$LOGDIR"/scanvi_downstream_pipeline_20251218_v1_6_2.py_$(date +%Y%m%d_%H%M%S).log 2>&1 &
PID=$!

echo "=========================================="
echo "scANVI Downstream Pipeline v1.6.2"
echo "=========================================="
echo "Process started with PID: $PID"
echo "Log directory: $LOGDIR"
echo "Monitor with: tail -f $LOGDIR/scanvi_downstream_pipeline_20251218_v1_6_2.py_*.log"
echo "=========================================="

# Check if process is running
sleep 2
if ps -p $PID > /dev/null 2>&1; then
    echo "✓ Process is running"
    echo ""
    ps aux | grep scanvi_downstream_pipeline_20251218_v1_6_2.py | grep -v grep
else
    echo "✗ Process failed to start"
    echo "Check latest log:"
    ls -t "$LOGDIR"/scanvi_downstream_pipeline_20251218_v1_6_2.py_*.log 2>/dev/null | head -1 | xargs tail -20
fi

