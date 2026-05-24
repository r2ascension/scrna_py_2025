#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
from pathlib import Path

import pyarrow.feather as feather


def convert_file(src: Path, dst: Path) -> dict:
    table = feather.read_table(src)
    columns = list(table.column_names)

    if "features" in columns:
        out = table
        mode = "already_compatible"
    elif "motifs" in columns:
        renamed = ["features" if c == "motifs" else c for c in columns]
        out = table.rename_columns(renamed)
        mode = "renamed_motifs_to_features"
    else:
        raise ValueError(f"No 'features' or 'motifs' column found in {src}")

    dst.parent.mkdir(parents=True, exist_ok=True)
    feather.write_feather(out, dst)

    return {
        "src": str(src),
        "dst": str(dst),
        "mode": mode,
        "nrows": out.num_rows,
        "ncols": len(out.column_names),
        "has_features": "features" in out.column_names,
        "has_motifs": "motifs" in out.column_names,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare R-SCENIC-compatible cisTarget feather databases.")
    parser.add_argument("--input-dir", default="/home/h2048/data/index_genome/cisTarget_databases")
    parser.add_argument("--output-dir", default="/home/h2048/data/index_genome/cisTarget_databases_rscenic")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    patterns = [
        "*10kb*.feather",
        "*500bp*.feather",
    ]

    files = []
    for pat in patterns:
        files.extend(sorted(input_dir.glob(pat)))
    files = sorted({p.resolve() for p in files})

    if not files:
        raise FileNotFoundError(f"No feather files found in {input_dir}")

    results = []
    for src in files:
        dst = output_dir / src.name
        results.append(convert_file(src, dst))

    summary = {
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "files": results,
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
