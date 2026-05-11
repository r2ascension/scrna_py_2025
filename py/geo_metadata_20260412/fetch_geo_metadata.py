#!/usr/bin/env python3
"""
Fetch GSM metadata from GEO for a list of GSE accessions.

Outputs:
  - GSE_GSM_metadata_full.csv: richer metadata table
  - GSE_GSM_metadata_aligned.csv: columns aligned to the user-provided reference CSV
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import Iterable

try:
    import GEOparse
    import pandas as pd
except ImportError:
    print("Please install required packages: pip install GEOparse pandas", file=sys.stderr)
    sys.exit(1)


DEFAULT_GSE_IDS = [
    "GSE136825",
    "GSE158277",
    "GSE189690",
    "GSE198950",
    "GSE172305",
    "GSE250580",
    "GSE313186",
    "GSE288682",
    "GSE72713",
    "GSE179265",
]

ALIGNED_COLUMNS = [
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

FULL_COLUMNS = [
    "GSE",
    "GSE_title",
    "GSE_summary",
    "GSE_type",
    "GSM",
    "title",
    "source_name_ch1",
    "organism_ch1",
    "characteristics_ch1",
    "molecule_ch1",
    "extract_protocol_ch1",
    "library_strategy",
    "library_source",
    "platform_id",
    "instrument_model",
    "description",
    "data_processing",
    "submission_date",
    "last_update_date",
    "contact_name",
    "geo_accession",
    "status",
    "type",
    "channel_count",
    "series_id",
    "supplementary_file",
    "data_completeness",
    "fetch_status",
    "fetch_error",
]


def first_value(values: Iterable[str] | None, default: str = "") -> str:
    if not values:
        return default
    if isinstance(values, str):
        return values
    values = list(values)
    return values[0] if values else default



def joined(values: Iterable[str] | None, sep: str = " | ") -> str:
    if not values:
        return ""
    if isinstance(values, str):
        return values
    return sep.join(str(v) for v in values if v is not None)



def normalize_text(value: str) -> str:
    return " ".join(str(value).replace("\n", " ").replace("\r", " ").split())



def fetch_gse_metadata(gse_id: str, cache_dir: Path, sleep_seconds: float = 1.0) -> list[dict[str, str]]:
    print(f"\n{'=' * 72}")
    print(f"Fetching {gse_id} ...")

    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(cache_dir), silent=True)
    except Exception as exc:  # pragma: no cover - depends on remote server
        print(f"  ERROR: Failed to fetch {gse_id}: {exc}")
        return [
            {
                "GSE": gse_id,
                "GSE_title": "",
                "GSE_summary": "",
                "GSE_type": "",
                "GSM": "",
                "title": "",
                "source_name_ch1": "",
                "organism_ch1": "",
                "characteristics_ch1": "",
                "molecule_ch1": "",
                "extract_protocol_ch1": "",
                "library_strategy": "",
                "library_source": "",
                "platform_id": "",
                "instrument_model": "",
                "description": "",
                "data_processing": "",
                "submission_date": "",
                "last_update_date": "",
                "contact_name": "",
                "geo_accession": "",
                "status": "",
                "type": "",
                "channel_count": "",
                "series_id": "",
                "supplementary_file": "",
                "data_completeness": "fetch_failed",
                "fetch_status": "error",
                "fetch_error": str(exc),
            }
        ]

    series_title = first_value(gse.metadata.get("title", [""]))
    series_summary = first_value(gse.metadata.get("summary", [""]))
    series_type = joined(gse.metadata.get("type", [""]), sep="; ")
    series_supplementary = joined(gse.metadata.get("supplementary_file", [""]))

    print(f"  Series title: {series_title}")
    print(f"  Samples found: {len(gse.gsms)}")

    rows: list[dict[str, str]] = []
    for gsm_name, gsm in sorted(gse.gsms.items()):
        supplementary = joined(gsm.metadata.get("supplementary_file", [""]))
        if not supplementary:
            supplementary = series_supplementary
        row = {
            "GSE": gse_id,
            "GSE_title": normalize_text(series_title),
            "GSE_summary": normalize_text(series_summary),
            "GSE_type": normalize_text(series_type),
            "GSM": gsm_name,
            "title": normalize_text(first_value(gsm.metadata.get("title", [""]))),
            "source_name_ch1": normalize_text(first_value(gsm.metadata.get("source_name_ch1", [""]))),
            "organism_ch1": normalize_text(first_value(gsm.metadata.get("organism_ch1", [""]))),
            "characteristics_ch1": normalize_text(joined(gsm.metadata.get("characteristics_ch1", [""]))),
            "molecule_ch1": normalize_text(first_value(gsm.metadata.get("molecule_ch1", [""]))),
            "extract_protocol_ch1": normalize_text(first_value(gsm.metadata.get("extract_protocol_ch1", [""]))),
            "library_strategy": normalize_text(first_value(gsm.metadata.get("library_strategy", [""]))),
            "library_source": normalize_text(first_value(gsm.metadata.get("library_source", [""]))),
            "platform_id": normalize_text(first_value(gsm.metadata.get("platform_id", [""]))),
            "instrument_model": normalize_text(first_value(gsm.metadata.get("instrument_model", [""]))),
            "description": normalize_text(joined(gsm.metadata.get("description", [""]))),
            "data_processing": normalize_text(first_value(gsm.metadata.get("data_processing", [""]))),
            "submission_date": normalize_text(first_value(gsm.metadata.get("submission_date", [""]))),
            "last_update_date": normalize_text(first_value(gsm.metadata.get("last_update_date", [""]))),
            "contact_name": normalize_text(first_value(gsm.metadata.get("contact_name", [""]))),
            "geo_accession": normalize_text(first_value(gsm.metadata.get("geo_accession", [""]))),
            "status": normalize_text(first_value(gsm.metadata.get("status", [""]))),
            "type": normalize_text(first_value(gsm.metadata.get("type", [""]))),
            "channel_count": normalize_text(first_value(gsm.metadata.get("channel_count", [""]))),
            "series_id": normalize_text(joined(gsm.metadata.get("series_id", [""]))),
            "supplementary_file": normalize_text(supplementary),
            "data_completeness": "complete",
            "fetch_status": "ok",
            "fetch_error": "",
        }
        rows.append(row)

    time.sleep(sleep_seconds)
    return rows



def build_aligned_dataframe(full_df: pd.DataFrame) -> pd.DataFrame:
    aligned = full_df.copy()
    for col in ALIGNED_COLUMNS:
        if col not in aligned.columns:
            aligned[col] = ""
    return aligned.loc[:, ALIGNED_COLUMNS].sort_values(["GSE", "GSM"], na_position="last")



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch GEO GSM metadata for a set of GSE accessions")
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent / "output"),
        help="Directory where output CSV files will be written",
    )
    parser.add_argument(
        "--cache-dir",
        default="/tmp/geo_cache_20260412",
        help="Directory for GEOparse cache files",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=1.0,
        help="Pause between GSE requests",
    )
    parser.add_argument(
        "--gse",
        nargs="*",
        default=DEFAULT_GSE_IDS,
        help="Override default GSE accession list",
    )
    return parser.parse_args()



def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    cache_dir = Path(args.cache_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, str]] = []
    for gse_id in args.gse:
        all_rows.extend(fetch_gse_metadata(gse_id=gse_id, cache_dir=cache_dir, sleep_seconds=args.sleep_seconds))

    if not all_rows:
        print("No rows fetched.", file=sys.stderr)
        return 1

    full_df = pd.DataFrame(all_rows)
    for col in FULL_COLUMNS:
        if col not in full_df.columns:
            full_df[col] = ""
    full_df = full_df.loc[:, FULL_COLUMNS].sort_values(["GSE", "GSM"], na_position="last")
    aligned_df = build_aligned_dataframe(full_df)

    full_csv = output_dir / "GSE_GSM_metadata_full.csv"
    aligned_csv = output_dir / "GSE_GSM_metadata_aligned.csv"
    full_df.to_csv(full_csv, index=False, encoding="utf-8")
    aligned_df.to_csv(aligned_csv, index=False, encoding="utf-8")

    ok_rows = (full_df["fetch_status"] == "ok").sum()
    error_rows = (full_df["fetch_status"] == "error").sum()

    print(f"\n{'=' * 72}")
    print("DONE")
    print(f"Output directory: {output_dir}")
    print(f"Full metadata CSV: {full_csv}")
    print(f"Aligned metadata CSV: {aligned_csv}")
    print(f"Requested GSE series: {len(args.gse)}")
    print(f"Rows written: {len(full_df)}")
    print(f"Successful GSM rows: {ok_rows}")
    print(f"Error placeholder rows: {error_rows}")
    print("\nRows per GSE:")
    counts = full_df.groupby("GSE", dropna=False).size().sort_index()
    for gse_id, n_rows in counts.items():
        print(f"  {gse_id}: {n_rows}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
