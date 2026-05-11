#!/usr/bin/env python3
"""
Build a corrected reference CSV from the GEO-aligned metadata table.

Design goals:
- Keep exactly the same schema as the original reference-style table
- Use current GEO-derived metadata as the authoritative source
- Produce a stable, review-friendly, sorted CSV
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print("Please install required packages: pip install pandas", file=sys.stderr)
    sys.exit(1)


REFERENCE_COLUMNS = [
    "GSE",
    "GSE_title",
    "GSE_type",
    "GSM",
    "title",
    "source_name_ch1",
    "organism_ch1",
    "characteristics_ch1",
    "platform_id",
    "instrument_model",
    "library_strategy",
    "library_source",
    "supplementary_file",
    "data_completeness",
]


def normalize_cell(value: object) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).replace("\r", " ").replace("\n", " ").split())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a corrected reference CSV from aligned GEO metadata")
    parser.add_argument("--input", required=True, help="Input aligned GEO metadata CSV")
    parser.add_argument("--output", required=True, help="Output corrected reference CSV")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()

    if not input_path.exists():
        print(f"Input file not found: {input_path}", file=sys.stderr)
        return 1

    df = pd.read_csv(input_path, dtype=str).fillna("")

    missing_columns = [col for col in REFERENCE_COLUMNS if col not in df.columns]
    if missing_columns:
        print(f"Missing required columns: {missing_columns}", file=sys.stderr)
        return 1

    corrected = df.loc[:, REFERENCE_COLUMNS].copy()
    for col in REFERENCE_COLUMNS:
        corrected[col] = corrected[col].map(normalize_cell)

    corrected = corrected.sort_values(["GSE", "GSM"], kind="stable", na_position="last")
    corrected = corrected.drop_duplicates(subset=["GSE", "GSM"], keep="first")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    corrected.to_csv(output_path, index=False, encoding="utf-8")

    print(f"Wrote corrected reference CSV: {output_path}")
    print(f"Rows: {len(corrected)}")
    print(f"Unique GSE: {corrected['GSE'].nunique()}")
    print(f"Unique GSM: {corrected['GSM'].nunique()}")
    print("Rows per GSE:")
    counts = corrected.groupby("GSE", dropna=False).size().sort_index()
    for gse_id, n_rows in counts.items():
        print(f"  {gse_id}: {n_rows}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
