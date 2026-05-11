#!/usr/bin/env python3
"""
Compare fetched GEO metadata with a reference CSV.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

KEY_COLUMNS = ["GSE", "GSM"]
COMPARE_COLUMNS = [
    "GSE_title",
    "GSE_type",
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



def normalize_value(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return " ".join(text.split())



def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in KEY_COLUMNS + COMPARE_COLUMNS:
        if col not in out.columns:
            out[col] = ""
        out[col] = out[col].map(normalize_value)
    return out.loc[:, KEY_COLUMNS + COMPARE_COLUMNS]



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare GEO fetched metadata to a reference CSV")
    parser.add_argument("--reference", required=True, help="Reference CSV path")
    parser.add_argument("--fetched", required=True, help="Fetched/aligned CSV path")
    parser.add_argument("--output-dir", required=True, help="Directory for comparison outputs")
    return parser.parse_args()



def main() -> int:
    args = parse_args()
    reference_path = Path(args.reference).resolve()
    fetched_path = Path(args.fetched).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    reference = normalize_df(pd.read_csv(reference_path, dtype=str, keep_default_na=False))
    fetched = normalize_df(pd.read_csv(fetched_path, dtype=str, keep_default_na=False))

    ref_indexed = reference.set_index(KEY_COLUMNS, drop=False)
    fet_indexed = fetched.set_index(KEY_COLUMNS, drop=False)

    ref_keys = set(ref_indexed.index)
    fet_keys = set(fet_indexed.index)

    missing_in_fetched = sorted(ref_keys - fet_keys)
    extra_in_fetched = sorted(fet_keys - ref_keys)
    common_keys = sorted(ref_keys & fet_keys)

    mismatch_rows: list[dict[str, str]] = []
    for key in common_keys:
        ref_row = ref_indexed.loc[key]
        fet_row = fet_indexed.loc[key]
        for col in COMPARE_COLUMNS:
            ref_val = normalize_value(ref_row[col])
            fet_val = normalize_value(fet_row[col])
            if ref_val != fet_val:
                mismatch_rows.append(
                    {
                        "GSE": key[0],
                        "GSM": key[1],
                        "column": col,
                        "reference": ref_val,
                        "fetched": fet_val,
                    }
                )

    missing_df = pd.DataFrame(missing_in_fetched, columns=KEY_COLUMNS)
    extra_df = pd.DataFrame(extra_in_fetched, columns=KEY_COLUMNS)
    mismatch_df = pd.DataFrame(mismatch_rows)

    missing_csv = output_dir / "missing_in_fetched.csv"
    extra_csv = output_dir / "extra_in_fetched.csv"
    mismatch_csv = output_dir / "field_mismatches.csv"
    summary_json = output_dir / "comparison_summary.json"

    missing_df.to_csv(missing_csv, index=False, encoding="utf-8")
    extra_df.to_csv(extra_csv, index=False, encoding="utf-8")
    mismatch_df.to_csv(mismatch_csv, index=False, encoding="utf-8")

    summary = {
        "reference_rows": int(len(reference)),
        "fetched_rows": int(len(fetched)),
        "common_keys": int(len(common_keys)),
        "missing_in_fetched": int(len(missing_in_fetched)),
        "extra_in_fetched": int(len(extra_in_fetched)),
        "field_mismatch_count": int(len(mismatch_rows)),
        "columns_compared": COMPARE_COLUMNS,
        "outputs": {
            "missing_in_fetched_csv": str(missing_csv),
            "extra_in_fetched_csv": str(extra_csv),
            "field_mismatches_csv": str(mismatch_csv),
        },
    }
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
