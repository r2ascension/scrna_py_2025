#!/usr/bin/env python3
"""Prepare sample-level and cell-type-level site manifests from audited metadata."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Iterable

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PY_ROOT = SCRIPT_DIR.parent
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from non_unified_airway.common import DEFAULT_REPO_ROOT, build_sample_manifest, write_json, write_tsv

DEFAULT_METADATA_TSV = DEFAULT_REPO_ROOT / "data/py/0524/non_unified_airway/metadata_audit/metadata_canonical.tsv.gz"
DEFAULT_OUT_DIR = DEFAULT_REPO_ROOT / "data/py/0524/non_unified_airway/site_manifests"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-tsv", type=Path, default=DEFAULT_METADATA_TSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def summarize_celltypes(metadata: pd.DataFrame, celltype_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if celltype_col not in metadata.columns:
        empty = pd.DataFrame(columns=[celltype_col, "sample", "site_group", "disease_group", "n_cells"])
        return empty, empty
    sample_by_celltype = (
        metadata.groupby([celltype_col, "sample", "site_group", "disease_group"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values([celltype_col, "site_group", "sample"])
        .reset_index(drop=True)
    )
    celltype_site_summary = (
        metadata.groupby([celltype_col, "site_group", "analysis_compartment", "disease_group"], dropna=False)
        .agg(n_cells=(celltype_col, "size"), n_samples=("sample", "nunique"))
        .reset_index()
        .sort_values([celltype_col, "site_group", "disease_group"])
        .reset_index(drop=True)
    )
    return sample_by_celltype, celltype_site_summary


def subset_sample_manifest(sample_manifest: pd.DataFrame, *, disease_group: str | None = None, compartments: Iterable[str] | None = None) -> pd.DataFrame:
    out = sample_manifest.copy()
    if disease_group is not None:
        out = out.loc[out["disease_group"] == disease_group].copy()
    if compartments is not None:
        compartments = set(compartments)
        out = out.loc[out["analysis_compartment"].isin(compartments)].copy()
    return out.reset_index(drop=True)


def run_prepare_site_groups(metadata_tsv: Path, out_dir: Path) -> dict[str, int]:
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(metadata_tsv, sep="\t")
    sample_manifest = build_sample_manifest(metadata)

    write_tsv(out_dir / "allcells_sample_manifest.tsv", sample_manifest)
    write_tsv(
        out_dir / "healthy_mainline_sample_manifest.tsv",
        subset_sample_manifest(sample_manifest, disease_group="healthy"),
    )
    write_tsv(
        out_dir / "upper_airway_sample_manifest.tsv",
        subset_sample_manifest(sample_manifest, compartments=["upper_airway"]),
    )
    write_tsv(
        out_dir / "lower_airway_sample_manifest.tsv",
        subset_sample_manifest(sample_manifest, compartments=["lower_airway", "lung_parenchyma"]),
    )

    sample_by_l2, celltype_site_l2 = summarize_celltypes(metadata, "cell_type_L2")
    sample_by_l3, celltype_site_l3 = summarize_celltypes(metadata, "cell_type_L3")
    write_tsv(out_dir / "sample_by_celltype_L2.tsv", sample_by_l2)
    write_tsv(out_dir / "sample_by_celltype_L3.tsv", sample_by_l3)
    write_tsv(out_dir / "celltype_site_summary_L2.tsv", celltype_site_l2)
    write_tsv(out_dir / "celltype_site_summary_L3.tsv", celltype_site_l3)

    confidence_summary = (
        metadata.groupby(["confidence_group", "site_group", "disease_group"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["confidence_group", "site_group", "disease_group"])
        .reset_index(drop=True)
    )
    write_tsv(out_dir / "confidence_group_summary.tsv", confidence_summary)

    summary = {
        "n_cells": int(metadata.shape[0]),
        "n_sample_manifest_rows": int(sample_manifest.shape[0]),
        "n_healthy_mainline_rows": int(subset_sample_manifest(sample_manifest, disease_group="healthy").shape[0]),
        "n_upper_airway_rows": int(subset_sample_manifest(sample_manifest, compartments=["upper_airway"]).shape[0]),
        "n_lower_airway_rows": int(subset_sample_manifest(sample_manifest, compartments=["lower_airway", "lung_parenchyma"]).shape[0]),
        "n_celltype_l2_rows": int(sample_by_l2.shape[0]),
        "n_celltype_l3_rows": int(sample_by_l3.shape[0]),
    }
    write_json(out_dir / "site_manifest_summary.json", summary)
    return summary


def main() -> int:
    args = parse_args()
    summary = run_prepare_site_groups(args.metadata_tsv, args.out_dir)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
