#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc

DATE_TAG = "20260523"
DEFAULT_BASE_OUT = Path(f"/home/h2048/data/py/0523/full_nocovid_{DATE_TAG}")
DEFAULT_OUTPUT_DIR = DEFAULT_BASE_OUT / f"all_lineage_eval_{DATE_TAG}"
STROMAL_VERSION = "v1_0_branchwise_nocovid_20260523"


def _read_h5ad_metrics(
    path: Path,
    *,
    pred_key: str | None = None,
    final_key: str | None = None,
    conf_key: str | None = None,
    source_key: str = "data_source",
    query_value: str = "query",
) -> dict[str, Any]:
    if not path.exists():
        return {"status": "missing", "path": str(path)}

    adata = None
    try:
        adata = sc.read_h5ad(path, backed="r")
        metrics: dict[str, Any] = {
            "status": "ok",
            "path": str(path),
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        }
        obs = adata.obs
        if source_key in obs.columns:
            src = obs[source_key].astype(str)
            qry_mask = src.eq(query_value)
            ref_mask = src.eq("reference")
            metrics["n_query"] = int(qry_mask.sum())
            metrics["n_reference"] = int(ref_mask.sum())
        else:
            qry_mask = pd.Series(True, index=obs.index)
            metrics["n_query"] = int(adata.n_obs)
            metrics["n_reference"] = 0

        if pred_key and pred_key in obs.columns:
            pred = obs.loc[qry_mask, pred_key].astype(str)
            metrics["pred_top3"] = pred.value_counts().head(3).to_dict()
        if final_key and final_key in obs.columns:
            final = obs.loc[qry_mask, final_key].astype(str)
            metrics["final_top3"] = final.value_counts().head(3).to_dict()
            metrics["n_final_unknown"] = int((final == "Unknown").sum())
        if conf_key and conf_key in obs.columns:
            conf = pd.to_numeric(obs.loc[qry_mask, conf_key], errors="coerce")
            if int(conf.notna().sum()) > 0:
                metrics["confidence_median"] = float(conf.median())
                metrics["confidence_mean"] = float(conf.mean())

        for key in ["schpl_rejected", "schpl_novel_candidate", "schpl_pred"]:
            if key in obs.columns:
                vals = obs.loc[qry_mask, key]
                if key.endswith("_rejected") or key.endswith("_candidate"):
                    vals = vals.astype(bool)
                    metrics[f"n_{key}"] = int(vals.sum())
                else:
                    metrics[f"{key}_top3"] = vals.astype(str).value_counts().head(3).to_dict()
        return metrics
    finally:
        if adata is not None and getattr(adata, "file", None) is not None:
            adata.file.close()


def _lineage_specs(base_out: Path) -> list[dict[str, Any]]:
    return [
        {
            "lineage": "epithelial_query",
            "stage": "mapped_query",
            "path": base_out / f"epithelial_scarches_query_v1_2_nocovid_{DATE_TAG}" / "epithelial_query_mapped_v1_2.h5ad",
            "pred_key": "cell_type_scarches_pred",
            "final_key": "cell_type_scarches_final",
            "conf_key": "scarches_confidence",
            "source_key": "data_source",
            "query_value": "query",
        },
        {
            "lineage": "epithelial",
            "stage": "schpl",
            "path": base_out / f"epithelial_merge_schpl_v1_0_nocovid_{DATE_TAG}" / f"epithelial_reference_plus_query_schpl_v1_0_nocovid_{DATE_TAG}.h5ad",
            "pred_key": "cell_type_scarches_pred",
            "final_key": "cell_type_scarches_final",
            "conf_key": "scarches_confidence",
        },
        {
            "lineage": "bcell",
            "stage": "merged",
            "path": base_out / f"bcell_scarches_query_nocovid_{DATE_TAG}" / "reference_plus_query_merged_L2.h5ad",
            "pred_key": "Cell_Type_L2_pred",
            "final_key": "Cell_Type_L2_final",
            "conf_key": "mapping_confidence",
        },
        {
            "lineage": "bcell",
            "stage": "schpl",
            "path": base_out / f"bcell_merge_schpl_nocovid_{DATE_TAG}" / f"bcell_reference_plus_query_schpl_v1_0_nocovid_{DATE_TAG}.h5ad",
            "pred_key": "Cell_Type_L2_pred",
            "final_key": "Cell_Type_L2_final",
            "conf_key": "mapping_confidence",
        },
        {
            "lineage": "tnk",
            "stage": "merged",
            "path": base_out / f"tcell_only_merged_nocovid_{DATE_TAG}" / "tcell_only_merged_nocovid_results.h5ad",
            "pred_key": "scanvi_pred",
            "final_key": "scanvi_pred",
            "conf_key": "scanvi_confidence",
        },
        {
            "lineage": "tnk",
            "stage": "schpl",
            "path": base_out / f"tcell_merge_schpl_nocovid_{DATE_TAG}" / f"tnk_reference_plus_query_schpl_v1_0_nocovid_{DATE_TAG}.h5ad",
            "pred_key": "scanvi_pred",
            "final_key": "scanvi_pred",
            "conf_key": "scanvi_confidence",
        },
        {
            "lineage": "myeloid",
            "stage": "merged",
            "path": base_out / f"myeloid_only_merged_nocovid_{DATE_TAG}" / "myeloid_only_merged_nocovid_results.h5ad",
            "pred_key": "scanvi_pred",
            "final_key": "scanvi_pred",
            "conf_key": "scanvi_confidence",
        },
        {
            "lineage": "myeloid",
            "stage": "schpl",
            "path": base_out / f"myeloid_merge_schpl_nocovid_{DATE_TAG}" / f"myeloid_reference_plus_query_schpl_v1_0_nocovid_{DATE_TAG}.h5ad",
            "pred_key": "scanvi_pred",
            "final_key": "scanvi_pred",
            "conf_key": "scanvi_confidence",
        },
    ]


def _stromal_specs(base_out: Path) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    query_manifest = base_out / f"stromal_scarches_query_branchwise_nocovid_{DATE_TAG}" / "branch_query_manifest.json"
    if query_manifest.exists():
        with open(query_manifest) as f:
            qry = json.load(f)
        for branch, info in qry.get("branches", {}).items():
            specs.append(
                {
                    "lineage": f"stromal_{branch}",
                    "stage": "query_branch",
                    "path": Path(info["query_h5ad"]),
                    "pred_key": info.get("label_key", "cell_type_scarches_pred"),
                    "final_key": info.get("final_label_key", "cell_type_scarches_final"),
                    "conf_key": info.get("confidence_key", "scarches_confidence"),
                    "source_key": "data_source",
                    "query_value": "query",
                }
            )
    for branch in ["endothelial", "fibroblast", "smc"]:
        specs.append(
            {
                "lineage": f"stromal_{branch}",
                "stage": "schpl",
                "path": base_out / f"stromal_merge_schpl_nocovid_{DATE_TAG}" / branch / f"stromal_{branch}_reference_plus_query_schpl_{STROMAL_VERSION}.h5ad",
                "pred_key": "cell_type_scarches_pred",
                "final_key": "cell_type_scarches_final",
                "conf_key": "scarches_confidence",
            }
        )
    return specs


def build_summary(base_out: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for spec in _lineage_specs(base_out) + _stromal_specs(base_out):
        metrics = _read_h5ad_metrics(
            spec["path"],
            pred_key=spec.get("pred_key"),
            final_key=spec.get("final_key"),
            conf_key=spec.get("conf_key"),
            source_key=spec.get("source_key", "data_source"),
            query_value=spec.get("query_value", "query"),
        )
        row = {"lineage": spec["lineage"], "stage": spec["stage"]}
        row.update(metrics)
        rows.append(row)
    return pd.DataFrame(rows)


def write_markdown(summary_df: pd.DataFrame, output_path: Path) -> None:
    lines = [f"# All-lineage no-COVID evaluation ({DATE_TAG})", ""]
    if summary_df.empty:
        lines.append("No outputs found.")
    else:
        display_cols = [
            col for col in [
                "lineage",
                "stage",
                "status",
                "n_reference",
                "n_query",
                "n_final_unknown",
                "n_schpl_rejected",
                "n_schpl_novel_candidate",
                "confidence_median",
                "path",
            ] if col in summary_df.columns
        ]
        lines.append(summary_df[display_cols].to_markdown(index=False))
    output_path.write_text("\n".join(lines))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize all-lineage no-COVID outputs")
    parser.add_argument("--base-out", type=Path, default=DEFAULT_BASE_OUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_df = build_summary(args.base_out)
    csv_path = args.output_dir / "all_lineage_summary.csv"
    md_path = args.output_dir / "all_lineage_summary.md"
    summary_df.to_csv(csv_path, index=False)
    write_markdown(summary_df, md_path)

    print("=" * 80)
    print("All-lineage no-COVID evaluation summary")
    print("=" * 80)
    print(summary_df.to_string(index=False) if not summary_df.empty else "No outputs found")
    print(f"\nCSV: {csv_path}")
    print(f"MD : {md_path}")


if __name__ == "__main__":
    main()
