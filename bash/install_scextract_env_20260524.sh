#!/usr/bin/env bash
set -euo pipefail

WORKSPACE_ROOT="/home/h2048"
CONDA_ROOT="$WORKSPACE_ROOT/miniconda3"
ENV_PATH="$CONDA_ROOT/envs/scextract_env"
REPO_DIR="$WORKSPACE_ROOT/tools/scExtract"
REQ_FILE="$REPO_DIR/requirements.txt"
PYTHON_BIN="$ENV_PATH/bin/python"

if [[ ! -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  echo "[ERROR] conda.sh not found under $CONDA_ROOT" >&2
  exit 1
fi

if [[ ! -d "$REPO_DIR" ]]; then
  echo "[ERROR] scExtract repo not found: $REPO_DIR" >&2
  exit 1
fi

source "$CONDA_ROOT/etc/profile.d/conda.sh"

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[INFO] Creating environment at $ENV_PATH"
  conda create -y -p "$ENV_PATH" python=3.11
else
  echo "[INFO] Reusing existing environment at $ENV_PATH"
fi

echo "[INFO] Upgrading packaging tools"
env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 "$PYTHON_BIN" -m pip install --upgrade pip setuptools wheel

if [[ -f "$REQ_FILE" ]]; then
  echo "[INFO] Installing scExtract upstream requirements"
  env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 "$PYTHON_BIN" -m pip install -r "$REQ_FILE"
fi

echo "[INFO] Installing workflow extras"
env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 "$PYTHON_BIN" -m pip install \
  anndata \
  scanpy \
  pandas \
  numpy \
  scipy \
  pyyaml \
  openpyxl \
  scikit-learn \
  pypdf \
  openai \
  pyfiglet \
  colorama \
  termcolor \
  celltypist \
  python-igraph \
  leidenalg \
  louvain

echo "[INFO] Verifying core imports"
env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 "$PYTHON_BIN" - <<'PY'
import sys
mods = ["anndata", "scanpy", "yaml", "openpyxl", "celltypist"]
for name in mods:
    __import__(name)
print("SCEXTRACT_ENV_OK", sys.executable)
PY

echo "[INFO] scExtract environment bootstrap complete"
