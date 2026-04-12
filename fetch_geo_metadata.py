#!/usr/bin/env python3
"""
fetch_geo_metadata.py - Fetch GSM metadata from GEO for a list of GSE accessions.

Usage:
    python fetch_geo_metadata.py

This script extracts GSM metadata from GEO for the following GSE datasets:
  GSE136825, GSE158277, GSE189690, GSE198950, GSE172305,
  GSE250580, GSE313186, GSE288682, GSE72713, GSE179265

Requirements:
    pip install GEOparse pandas

Output:
    GSE_GSM_metadata.csv - Comprehensive CSV with all sample metadata
"""

import os
import sys
import time

try:
    import GEOparse
    import pandas as pd
except ImportError:
    print("Please install required packages: pip install GEOparse pandas")
    sys.exit(1)

# GSE IDs extracted from the original file paths
GSE_IDS = [
    "GSE136825",  # GSE136825_genecounts_20190903.txt.gz
    "GSE158277",  # GSE158277_Supplementary_data_RNA.xlsx
    "GSE189690",  # GSE189690_gene_counts.txt.gz
    "GSE198950",  # GSE198950_CRS_GeneCount.txt.gz
    "GSE172305",  # GSE172305_RAW.tar
    "GSE250580",  # GSE250580_Combined_counts.csv.gz
    "GSE313186",  # GSE313186_Clinical_nasal_BPC_raw_featurecounts_matrix.txt.gz
    "GSE288682",  # GSE288682_count_matrix.txt.gz
    "GSE72713",   # GSE72713_RNA-seq_processed_data.txt.gz
    "GSE179265",  # GSE179265_readcount_genename.txt.gz
]

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "GSE_GSM_metadata.csv")


def fetch_gse_metadata(gse_id):
    """Fetch and parse metadata for all GSMs in a GSE."""
    print(f"\n{'='*60}")
    print(f"Fetching {gse_id}...")
    
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir="/tmp/geo_cache", silent=True)
    except Exception as e:
        print(f"  ERROR: Failed to fetch {gse_id}: {e}")
        return []
    
    # Get series-level metadata
    series_title = gse.metadata.get("title", [""])[0]
    series_summary = gse.metadata.get("summary", [""])[0]
    series_type = "; ".join(gse.metadata.get("type", [""]))
    platform_id = "; ".join(gse.metadata.get("platform_id", [""]))
    
    print(f"  Series title: {series_title}")
    print(f"  Number of samples: {len(gse.gsms)}")
    
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        row = {
            "GSE": gse_id,
            "GSE_title": series_title,
            "GSE_summary": series_summary[:200] + "..." if len(series_summary) > 200 else series_summary,
            "GSE_type": series_type,
            "GSM": gsm_name,
            "title": gsm.metadata.get("title", [""])[0],
            "source_name_ch1": gsm.metadata.get("source_name_ch1", [""])[0],
            "organism_ch1": gsm.metadata.get("organism_ch1", [""])[0],
            "characteristics_ch1": " | ".join(gsm.metadata.get("characteristics_ch1", [""])),
            "molecule_ch1": gsm.metadata.get("molecule_ch1", [""])[0],
            "extract_protocol_ch1": gsm.metadata.get("extract_protocol_ch1", [""])[0][:200],
            "library_strategy": gsm.metadata.get("library_strategy", [""])[0],
            "library_source": gsm.metadata.get("library_source", [""])[0],
            "platform_id": gsm.metadata.get("platform_id", [""])[0],
            "instrument_model": gsm.metadata.get("instrument_model", [""])[0],
            "description": " | ".join(gsm.metadata.get("description", [""])),
            "data_processing": gsm.metadata.get("data_processing", [""])[0][:200],
            "submission_date": gsm.metadata.get("submission_date", [""])[0],
            "last_update_date": gsm.metadata.get("last_update_date", [""])[0],
            "contact_name": gsm.metadata.get("contact_name", [""])[0],
            "geo_accession": gsm.metadata.get("geo_accession", [""])[0],
            "status": gsm.metadata.get("status", [""])[0],
            "type": gsm.metadata.get("type", [""])[0],
            "channel_count": gsm.metadata.get("channel_count", [""])[0],
            "series_id": " | ".join(gsm.metadata.get("series_id", [""])),
            "supplementary_file": " | ".join(gsm.metadata.get("supplementary_file", [""])),
        }
        rows.append(row)
    
    return rows


def main():
    os.makedirs("/tmp/geo_cache", exist_ok=True)
    
    all_rows = []
    for gse_id in GSE_IDS:
        rows = fetch_gse_metadata(gse_id)
        all_rows.extend(rows)
        time.sleep(2)  # Be polite to NCBI servers
    
    if not all_rows:
        print("\nNo samples found!")
        sys.exit(1)
    
    df = pd.DataFrame(all_rows)
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")
    
    print(f"\n{'='*60}")
    print(f"DONE!")
    print(f"Total GSE series: {len(GSE_IDS)}")
    print(f"Total GSM samples: {len(all_rows)}")
    print(f"Output: {OUTPUT_CSV}")
    print(f"Columns: {list(df.columns)}")
    
    # Print summary per GSE
    print(f"\nSummary per GSE:")
    for gse_id in GSE_IDS:
        n = len([r for r in all_rows if r["GSE"] == gse_id])
        print(f"  {gse_id}: {n} samples")


if __name__ == "__main__":
    main()
