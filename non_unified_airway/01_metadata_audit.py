#!/usr/bin/env python3
"""Run metadata audit for the non-unified airway workflow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PY_ROOT = SCRIPT_DIR.parent
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from non_unified_airway.common import (
    AUDIT_OBS_COLUMNS,
    DEFAULT_REPO_ROOT,
    add_canonical_fields,
    build_recommended_contrasts,
    build_sample_manifest,
    choose_first_existing_path,
    inspect_h5ad_light,
    read_obs_columns_h5py,
    recommend_model_formulas,
    write_json,
    write_tsv,
)

DEFAULT_CONFIG_MANIFEST = DEFAULT_REPO_ROOT / "data/py/0524/non_unified_airway/config/config_manifest.json"
DEFAULT_OUT_DIR = DEFAULT_REPO_ROOT / "data/py/0524/non_unified_airway/metadata_audit"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-manifest", type=Path, default=DEFAULT_CONFIG_MANIFEST)
    parser.add_argument("--h5ad", type=Path, default=None, help="Override the primary all-cells h5ad.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def load_config_manifest(config_manifest: Path) -> dict[str, Any]:
    if not config_manifest.exists():
        raise FileNotFoundError(f"Config manifest not found: {config_manifest}")
    return json.loads(config_manifest.read_text(encoding="utf-8"))


def pick_primary_h5ad(config: dict[str, Any], override_h5ad: Path | None) -> Path:
    if override_h5ad is not None:
        return override_h5ad
    candidates = config.get("allcells_reference_candidates", [])
    selected = choose_first_existing_path(candidates)
    if selected is None:
        raise FileNotFoundError(
            "No existing all-cells reference candidate found in config manifest. "
            "Run 00_config_non_unified_airway.py or pass --h5ad explicitly."
        )
    return selected


def build_crosstab(sample_manifest: pd.DataFrame, row_col: str, col_col: str, value_col: str = "sample") -> pd.DataFrame:
    if row_col not in sample_manifest.columns or col_col not in sample_manifest.columns:
        return pd.DataFrame(columns=[row_col, col_col, "n_samples"])
    table = (
        sample_manifest.groupby([row_col, col_col], dropna=False)[value_col]
        .nunique()
        .reset_index(name="n_samples")
        .sort_values([row_col, col_col])
        .reset_index(drop=True)
    )
    return table


def run_metadata_audit(h5ad_path: Path, out_dir: Path, config_manifest: Path | None = None) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)

    h5ad_info = inspect_h5ad_light(h5ad_path)
    metadata = read_obs_columns_h5py(h5ad_path, AUDIT_OBS_COLUMNS)
    annotated = add_canonical_fields(metadata)
    sample_manifest = build_sample_manifest(annotated)
    site_by_dataset = build_crosstab(sample_manifest, "site_group", "dataset")
    site_by_disease = build_crosstab(sample_manifest, "site_group", "disease_group")
    compartment_by_dataset = build_crosstab(sample_manifest, "analysis_compartment", "dataset")
    recommended_formulas = recommend_model_formulas(sample_manifest)
    recommended_contrasts = build_recommended_contrasts(sample_manifest)

    present_obs = pd.DataFrame({"obs_column": sorted(metadata.columns.tolist())})
    write_tsv(out_dir / "obs_columns_present.tsv", present_obs)
    write_tsv(out_dir / "metadata_canonical.tsv.gz", annotated.reset_index().rename(columns={"index": "cell_id"}), compression="gzip")
    write_tsv(out_dir / "sample_manifest.tsv", sample_manifest)
    write_tsv(out_dir / "site_by_dataset.tsv", site_by_dataset)
    write_tsv(out_dir / "site_by_disease.tsv", site_by_disease)
    write_tsv(out_dir / "analysis_compartment_by_dataset.tsv", compartment_by_dataset)
    write_tsv(out_dir / "recommended_contrasts.tsv", recommended_contrasts)
    write_json(out_dir / "recommended_model_formulas.json", recommended_formulas)

    summary = {
        "h5ad": h5ad_info,
        "config_manifest": str(config_manifest) if config_manifest is not None else None,
        "n_cells": int(annotated.shape[0]),
        "n_samples": int(sample_manifest["sample"].nunique()) if "sample" in sample_manifest.columns else 0,
        "n_site_groups": int(sample_manifest["site_group"].nunique()) if "site_group" in sample_manifest.columns else 0,
        "n_disease_groups": int(sample_manifest["disease_group"].nunique()) if "disease_group" in sample_manifest.columns else 0,
        "n_datasets": int(sample_manifest["dataset"].nunique()) if "dataset" in sample_manifest.columns else 0,
        "recommended_model_formulas": recommended_formulas,
    }
    write_json(out_dir / "preflight_summary.json", summary)
    return summary


def main() -> int:
    args = parse_args()
    config = load_config_manifest(args.config_manifest)
    h5ad_path = pick_primary_h5ad(config, args.h5ad)
    summary = run_metadata_audit(
        h5ad_path=h5ad_path,
        out_dir=args.out_dir,
        config_manifest=args.config_manifest,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
