import argparse
import shutil
import subprocess
import sys
from pathlib import Path

#!/usr/bin/env python3
"""
seurat_to_ann.py

Convert a Seurat object (.rds or .h5Seurat) to an AnnData H5AD file using R (Seurat + SeuratDisk).

Usage:
    python seurat_to_ann.py input.rds [output.h5ad]
    python seurat_to_ann.py input.h5Seurat output.h5ad
"""


def main():
    p = argparse.ArgumentParser(description="Convert Seurat -> AnnData (.h5ad) using R (SeuratDisk).")
    p.add_argument("input", help="Path to Seurat object (.rds or .h5Seurat)")
    p.add_argument("output", nargs="?", help="Path for output .h5ad (optional). If omitted, same basename with .h5ad is used.")
    args = p.parse_args()

    rscript = shutil.which("Rscript")
    if not rscript:
        print("Rscript not found in PATH. Install R and ensure Rscript is available.", file=sys.stderr)
        sys.exit(2)

    inp = Path(args.input).expanduser().resolve()
    if not inp.exists():
        print(f"Input file does not exist: {inp}", file=sys.stderr)
        sys.exit(2)

    outp = Path(args.output) if args.output else inp.with_suffix(".h5ad")
    outp = outp.expanduser().resolve()

    # R code: uses commandArgs to receive input and output paths
    r_code = r'''
args <- commandArgs(trailingOnly=TRUE)
infile <- normalizePath(args[1])
outfile <- args[2]

# Load required package
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  stop("R package 'SeuratDisk' is required. Install it in R: install.packages('SeuratDisk') or Bioc/Github as appropriate.")
}

is_h5 <- grepl("\\.h5Seurat$", infile, ignore.case=TRUE)
if (is_h5) {
  SeuratDisk::Convert(infile, dest = "h5ad", overwrite = TRUE)
  base <- sub("\\.h5Seurat$", "", infile, ignore.case=TRUE)
} else {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("R package 'Seurat' is required to read .rds Seurat objects. Install it in R: install.packages('Seurat').")
  }
  # create temporary base in tempdir()
  tmpbase <- tempfile(pattern = "seurat_tmp_", tmpdir = tempdir())
  tmp_h5 <- paste0(tmpbase, ".h5Seurat")
  obj <- readRDS(infile)
  SeuratDisk::SaveH5Seurat(obj, filename = tmp_h5, overwrite = TRUE)
  SeuratDisk::Convert(tmp_h5, dest = "h5ad", overwrite = TRUE)
  base <- tmpbase
}
h5ad_file <- paste0(base, ".h5ad")
if (!file.exists(h5ad_file)) {
  stop("Expected .h5ad file not created: ", h5ad_file)
}
file.copy(h5ad_file, outfile, overwrite = TRUE)
cat("WROTE:", outfile, "\n")
'''

    try:
        proc = subprocess.run(
            [rscript, "-e", r_code, str(inp), str(outp)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print(proc.stdout.strip())
    except subprocess.CalledProcessError as e:
        msg = e.stderr.strip() or e.stdout.strip()
        print("R conversion failed:", file=sys.stderr)
        print(msg, file=sys.stderr)
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()