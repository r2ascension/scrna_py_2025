#!/usr/bin/env python3
"""Audit healthy/site metadata contract for normal-airway ML inputs."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from normal_airway_ml_common_20260527 import (
    DEFAULT_OUTPUT_ROOT,
    build_sample_metadata,
    deep_get,
    ensure_dir,
    load_backed_obs_summary,
    load_config,
    prepare_obs_contract,
    summarize_top_counts,
    timestamp_slug,
    write_json,
)


def run_audit(input_h5ad: Path | str, cfg: dict[str, Any], output_dir: Path | str) -> dict[str, Any]:
    output_dir = ensure_dir(output_dir)
    inventory = load_backed_obs_summary(input_h5ad)
    backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        obs = backed.obs.copy()
        prepared, contract = prepare_obs_contract(obs, list(backed.obsm.keys()), cfg)
    finally:
        try:
            backed.file.close()
        except Exception:
            pass

    analysis = prepared.loc[contract["analysis_mask"]].copy()
    sample_meta = build_sample_metadata(prepared, contract["analysis_mask"]) if not analysis.empty else pd.DataFrame()

    summary = {
        "input_h5ad": str(Path(input_h5ad)),
        "shape": inventory["shape"],
        "obs_ncols": inventory["obs_ncols"],
        "layers": inventory["layers"],
        "obsm": inventory["obsm"],
        "resolved_columns": {
            k: contract.get(k)
            for k in ["sample_col", "sample_unit_col", "dataset_col", "study_col", "batch_col", "tissue_col", "condition_col", "cell_type_col", "latent_key"]
        },
        "n_analysis_cells": int(contract["n_analysis_cells"]),
        "n_analysis_samples": int(sample_meta.shape[0]) if not sample_meta.empty else 0,
        "site_counts_all": summarize_top_counts(prepared["site_label"], top_n=10),
        "healthy_site_counts": summarize_top_counts(analysis["site_label"], top_n=10) if not analysis.empty else {},
        "healthy_dataset_counts": summarize_top_counts(analysis["dataset_resolved"], top_n=10) if not analysis.empty else {},
        "healthy_celltype_counts": summarize_top_counts(analysis["cell_type_resolved"], top_n=10) if not analysis.empty else {},
        "healthy_samples_per_site": (
            analysis.groupby("site_label", observed=True)["sample_resolved"].nunique().sort_values(ascending=False).to_dict()
            if not analysis.empty
            else {}
        ),
        "site_by_dataset_shape": [0, 0],
        "healthy_condition_counts": summarize_top_counts(analysis["condition_resolved"], top_n=10) if not analysis.empty else {},
    }

    site_dataset = pd.DataFrame()
    site_sample = pd.DataFrame()
    site_celltype = pd.DataFrame()
    if not analysis.empty:
        site_dataset = pd.crosstab(analysis["site_label"], analysis["dataset_resolved"])
        site_sample = pd.crosstab(analysis["site_label"], analysis["sample_resolved"])
        site_celltype = pd.crosstab(analysis["site_label"], analysis["cell_type_resolved"])
        summary["site_by_dataset_shape"] = [int(site_dataset.shape[0]), int(site_dataset.shape[1])]

    prepared[[
        "sample_resolved",
        "sample_unit_resolved",
        "dataset_resolved",
        "study_resolved",
        "batch_resolved",
        "condition_resolved",
        "raw_tissue_label",
        "site_label",
        "site_axis",
        "is_healthy",
        "cell_type_resolved",
    ]].to_csv(output_dir / "resolved_obs_contract.tsv.gz", sep="\t", index_label="obs_name", compression="gzip")

    if not sample_meta.empty:
        sample_meta.to_csv(output_dir / "sample_metadata.tsv", sep="\t", index=False)
    if not site_dataset.empty:
        site_dataset.to_csv(output_dir / "site_by_dataset.tsv", sep="\t")
    if not site_sample.empty:
        site_sample.to_csv(output_dir / "site_by_sample.tsv", sep="\t")
    if not site_celltype.empty:
        site_celltype.to_csv(output_dir / "site_by_celltype.tsv", sep="\t")

    write_json(summary, output_dir / "audit_summary.json")
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Audit metadata contract for normal-airway ML inputs")
    parser.add_argument("--input-h5ad", default=None, help="Input h5ad path")
    parser.add_argument("--config-yaml", default=None, help="Config YAML path")
    parser.add_argument("--output-dir", default=None, help="Explicit output directory")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    cfg = load_config(args.config_yaml)
    input_h5ad = Path(args.input_h5ad or deep_get(cfg, "data_contract", "input_h5ad", default=""))
    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_h5ad}")
    output_dir = Path(args.output_dir) if args.output_dir else ensure_dir(Path(deep_get(cfg, "run", "output_root", default=str(DEFAULT_OUTPUT_ROOT))) / f"normal_airway_ml_audit_{timestamp_slug()}")
    summary = run_audit(input_h5ad=input_h5ad, cfg=cfg, output_dir=output_dir)
    print(pd.Series(summary).to_json(force_ascii=False, indent=2))


if __name__ == "__main__":
    main()
