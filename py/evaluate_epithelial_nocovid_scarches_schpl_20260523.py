#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Summarize the epithelial no-COVID scArches + scHPL rerun outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

from query_covid_filter_helper_20260523 import (
    build_covid_filter_debug_frame,
    build_removed_kept_value_counts,
    summarize_covid_filter_debug_frame,
)

DEFAULT_ORIGINAL_QUERY = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/epithelial_cells.h5ad"
DEFAULT_MAPPED_QUERY = "/home/h2048/data/py/0523/epithelial_scarches_query_v1_2_nocovid_20260523/epithelial_query_mapped_v1_2.h5ad"
DEFAULT_SCHPL_H5AD = "/home/h2048/data/py/0523/epithelial_merge_schpl_v1_0_nocovid_20260523/epithelial_reference_plus_query_schpl_v1_0_nocovid_20260523.h5ad"
DEFAULT_OUTPUT_DIR = "/home/h2048/data/py/0523/epithelial_nocovid_eval_20260523"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate epithelial no-COVID scArches + scHPL outputs")
    parser.add_argument("--original-query-h5ad", default=DEFAULT_ORIGINAL_QUERY)
    parser.add_argument("--mapped-query-h5ad", default=DEFAULT_MAPPED_QUERY)
    parser.add_argument("--schpl-h5ad", default=DEFAULT_SCHPL_H5AD)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    return parser


def read_obs(path: Path, columns: list[str]) -> pd.DataFrame:
    adata = sc.read_h5ad(path, backed="r")
    try:
        present = [col for col in columns if col in adata.obs.columns]
        obs = adata.obs[present].copy()
        obs.index = adata.obs_names.astype(str)
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()
    return obs


def value_counts_table(series: pd.Series, name: str) -> pd.DataFrame:
    table = series.astype("string").fillna("<NA>").value_counts(dropna=False).rename("n_cells").reset_index()
    table.columns = [name, "n_cells"]
    table["fraction"] = table["n_cells"] / max(int(table["n_cells"].sum()), 1)
    return table


def main() -> None:
    args = build_arg_parser().parse_args()
    original_path = Path(args.original_query_h5ad).expanduser().resolve()
    mapped_path = Path(args.mapped_query_h5ad).expanduser().resolve()
    schpl_path = Path(args.schpl_h5ad).expanduser().resolve()
    out_dir = Path(args.output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    original_obs = read_obs(original_path, ["disease", "COVID_status", "dataset", "condition", "sample", "tissue"])
    mapped_obs = read_obs(
        mapped_path,
        [
            "disease",
            "COVID_status",
            "dataset",
            "condition",
            "sample",
            "tissue",
            "scanvi_fine_pred",
            "scanvi_fine_pred_filtered",
            "scanvi_fine_conf",
            "mapping_confidence",
        ],
    )
    schpl_obs = read_obs(
        schpl_path,
        [
            "data_source",
            "scanvi_fine_pred",
            "scanvi_fine_pred_filtered",
            "scanvi_fine_conf",
            "schpl_pred",
            "schpl_prob",
            "schpl_rejected",
            "schpl_reject_type",
            "schpl_novel_candidate",
            "dataset",
            "sample",
            "tissue",
        ],
    )

    filter_debug = build_covid_filter_debug_frame(original_obs)
    filter_summary = summarize_covid_filter_debug_frame(filter_debug)
    filter_summary.to_csv(out_dir / "covid_filter_summary.tsv", sep="\t", index=False)
    filter_debug.to_csv(out_dir / "covid_filter_audit.tsv.gz", sep="\t", compression="gzip")
    for col, table in build_removed_kept_value_counts(
        original_obs,
        filter_debug["remove_covid_related"],
        columns=["disease", "COVID_status", "dataset", "condition", "tissue"],
    ).items():
        table.to_csv(out_dir / f"covid_filter_{col}_counts.tsv", sep="\t")

    mapped_conf = pd.to_numeric(mapped_obs.get("scanvi_fine_conf"), errors="coerce")
    mapping_summary = pd.DataFrame(
        [
            {"metric": "mapped_query_cells", "value": int(len(mapped_obs))},
            {"metric": "expected_cells_after_filter", "value": int((~filter_debug["remove_covid_related"]).sum())},
            {"metric": "median_scanvi_fine_conf", "value": float(mapped_conf.median()) if mapped_conf.notna().any() else None},
            {"metric": "mean_scanvi_fine_conf", "value": float(mapped_conf.mean()) if mapped_conf.notna().any() else None},
            {"metric": "low_conf_cells_lt_0.5", "value": int((mapped_conf < 0.5).sum()) if mapped_conf.notna().any() else None},
            {"metric": "low_conf_fraction_lt_0.5", "value": float((mapped_conf < 0.5).mean()) if mapped_conf.notna().any() else None},
        ]
    )
    mapping_summary.to_csv(out_dir / "mapping_summary.tsv", sep="\t", index=False)
    if "scanvi_fine_pred" in mapped_obs.columns:
        value_counts_table(mapped_obs["scanvi_fine_pred"], "scanvi_fine_pred").to_csv(
            out_dir / "scanvi_fine_pred_counts.tsv", sep="\t", index=False
        )
    if "scanvi_fine_pred_filtered" in mapped_obs.columns:
        value_counts_table(mapped_obs["scanvi_fine_pred_filtered"], "scanvi_fine_pred_filtered").to_csv(
            out_dir / "scanvi_fine_pred_filtered_counts.tsv", sep="\t", index=False
        )

    schpl_query = schpl_obs.copy()
    if "data_source" in schpl_query.columns:
        schpl_query = schpl_query[schpl_query["data_source"].astype(str) == "query"].copy()
    schpl_rejected = schpl_query.get("schpl_rejected", pd.Series([False] * len(schpl_query), index=schpl_query.index)).astype(bool)
    schpl_novel = schpl_query.get("schpl_novel_candidate", pd.Series([False] * len(schpl_query), index=schpl_query.index)).astype(bool)
    schpl_prob = pd.to_numeric(schpl_query.get("schpl_prob"), errors="coerce") if "schpl_prob" in schpl_query.columns else pd.Series(dtype=float)
    schpl_summary = pd.DataFrame(
        [
            {"metric": "schpl_query_cells", "value": int(len(schpl_query))},
            {"metric": "schpl_rejected_cells", "value": int(schpl_rejected.sum())},
            {"metric": "schpl_rejected_fraction", "value": float(schpl_rejected.mean()) if len(schpl_query) else None},
            {"metric": "schpl_novel_candidate_cells", "value": int(schpl_novel.sum())},
            {"metric": "median_schpl_prob", "value": float(schpl_prob.median()) if schpl_prob.notna().any() else None},
        ]
    )
    schpl_summary.to_csv(out_dir / "schpl_summary.tsv", sep="\t", index=False)
    if "schpl_reject_type" in schpl_query.columns:
        value_counts_table(schpl_query["schpl_reject_type"], "schpl_reject_type").to_csv(
            out_dir / "schpl_reject_type_counts.tsv", sep="\t", index=False
        )
    if {"scanvi_fine_pred_filtered", "schpl_pred"}.issubset(schpl_query.columns):
        accepted = schpl_query.loc[~schpl_rejected, ["scanvi_fine_pred_filtered", "schpl_pred"]].copy()
        if not accepted.empty:
            pd.crosstab(accepted["scanvi_fine_pred_filtered"], accepted["schpl_pred"]).to_csv(
                out_dir / "accepted_label_crosstab.tsv", sep="\t"
            )

    removed_cells = int(filter_debug["remove_covid_related"].sum())
    kept_cells = int((~filter_debug["remove_covid_related"]).sum())
    report = f"""
# Epithelial no-COVID rerun evaluation

## Input paths

- Original query: `{original_path}`
- Mapped query: `{mapped_path}`
- scHPL merged output: `{schpl_path}`

## COVID filter summary

- Original query cells: {len(original_obs):,}
- Removed as COVID-related: {removed_cells:,}
- Retained for rerun: {kept_cells:,}
- Mapped query cells observed: {len(mapped_obs):,}

## scArches / scANVI summary

- Median `scanvi_fine_conf`: {mapping_summary.loc[mapping_summary['metric'] == 'median_scanvi_fine_conf', 'value'].iloc[0]}
- Mean `scanvi_fine_conf`: {mapping_summary.loc[mapping_summary['metric'] == 'mean_scanvi_fine_conf', 'value'].iloc[0]}
- Low-confidence cells (<0.5): {mapping_summary.loc[mapping_summary['metric'] == 'low_conf_cells_lt_0.5', 'value'].iloc[0]}

## scHPL summary

- Query cells in scHPL output: {schpl_summary.loc[schpl_summary['metric'] == 'schpl_query_cells', 'value'].iloc[0]}
- Rejected cells: {schpl_summary.loc[schpl_summary['metric'] == 'schpl_rejected_cells', 'value'].iloc[0]}
- Novel candidate cells: {schpl_summary.loc[schpl_summary['metric'] == 'schpl_novel_candidate_cells', 'value'].iloc[0]}

## Output tables

- `covid_filter_summary.tsv`
- `covid_filter_audit.tsv.gz`
- `mapping_summary.tsv`
- `scanvi_fine_pred_counts.tsv`
- `scanvi_fine_pred_filtered_counts.tsv`
- `schpl_summary.tsv`
- `schpl_reject_type_counts.tsv` (if available)
- `accepted_label_crosstab.tsv` (if available)
""".strip() + "\n"
    report_path = out_dir / "epithelial_nocovid_evaluation_report.md"
    report_path.write_text(report, encoding="utf-8")

    print("=" * 80)
    print("Evaluation complete")
    print("=" * 80)
    print(f"Output dir: {out_dir}")
    print(f"Report    : {report_path}")


if __name__ == "__main__":
    main()
