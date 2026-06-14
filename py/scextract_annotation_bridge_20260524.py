#!/usr/bin/env python3
"""Bridge selected scExtract-side annotations into downstream label contracts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd


SOURCE_COLUMN_MAP = {
    "scExtract": "cell_type_scextract",
    "celltypist_annotation": "cell_type_celltypist_scextract",
    "sctype_annotation": "cell_type_sctype_scextract",
    "mllmcelltype_annotation": "cell_type_mllmcelltype_scextract",
}


SOURCE_CHOICES = [
    "compare_only",
    "scExtract",
    "celltypist_annotation",
    "sctype_annotation",
    "mllmcelltype_annotation",
]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(payload: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def normalize_labels(series: pd.Series, unknown_label: str) -> pd.Series:
    values = series.astype("string").fillna(unknown_label).astype(str).str.strip()
    values = values.replace({"": unknown_label, "nan": unknown_label, "None": unknown_label, "<NA>": unknown_label})
    return values


def merge_rare_types(series: pd.Series, min_cells_per_type: int, unknown_label: str) -> tuple[pd.Series, int]:
    out = series.astype(str).copy()
    counts = out[out != unknown_label].value_counts()
    rare = counts[counts < int(min_cells_per_type)].index
    if len(rare) > 0:
        reassigned = int(out.isin(rare).sum())
        out.loc[out.isin(rare)] = unknown_label
        return out, reassigned
    return out, 0


def apply_bridge(
    adata: ad.AnnData,
    selected_label_source: str,
    source_column: str | None = None,
    unknown_label: str = "Unknown",
    min_cells_per_type: int = 10,
    confidence_column: str | None = None,
    min_confidence: float | None = None,
) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "selected_label_source": selected_label_source,
        "source_column": None,
        "status": "compare_only",
        "n_cells": int(adata.n_obs),
        "n_labeled": 0,
        "n_unknown": int(adata.n_obs),
        "n_unique_types": 0,
        "low_confidence_reassigned": 0,
        "rare_type_reassigned": 0,
    }

    adata.obs["annotation_source_selected"] = pd.Series(
        [selected_label_source] * adata.n_obs,
        index=adata.obs_names,
        dtype="string",
    )
    adata.obs["annotation_label_selected"] = pd.Series(
        [pd.NA] * adata.n_obs,
        index=adata.obs_names,
        dtype="string",
    )

    if selected_label_source == "compare_only":
        return summary

    resolved_column = source_column or SOURCE_COLUMN_MAP.get(selected_label_source, selected_label_source)
    if resolved_column not in adata.obs.columns:
        raise KeyError(f"Bridge source column not found: {resolved_column}")

    labels = normalize_labels(adata.obs[resolved_column], unknown_label)

    low_confidence_reassigned = 0
    if confidence_column and min_confidence is not None and confidence_column in adata.obs.columns:
        confidence = pd.to_numeric(adata.obs[confidence_column], errors="coerce")
        low_conf_mask = confidence.notna() & (confidence < float(min_confidence))
        low_confidence_reassigned = int(low_conf_mask.sum())
        if low_confidence_reassigned > 0:
            labels.loc[low_conf_mask] = unknown_label

    labels, rare_type_reassigned = merge_rare_types(labels, min_cells_per_type=min_cells_per_type, unknown_label=unknown_label)
    labels = labels.astype(str)

    categories = sorted(set(labels.tolist()) | {unknown_label})
    labels_cat = pd.Series(pd.Categorical(labels, categories=categories), index=adata.obs_names)

    adata.obs["labels_for_scanvi"] = labels_cat
    adata.obs["annotation_source_selected"] = pd.Series(
        [selected_label_source] * adata.n_obs,
        index=adata.obs_names,
        dtype="string",
    )
    adata.obs["annotation_label_selected"] = pd.Series(
        [resolved_column] * adata.n_obs,
        index=adata.obs_names,
        dtype="string",
    )
    adata.obs["scanvi_label_source"] = pd.Series(
        [selected_label_source] * adata.n_obs,
        index=adata.obs_names,
        dtype="string",
    )

    n_unknown = int((labels_cat.astype(str) == unknown_label).sum())
    n_labeled = int(adata.n_obs - n_unknown)
    unique_types = sorted({x for x in labels.astype(str).unique().tolist() if x != unknown_label})

    summary.update(
        {
            "status": "bridged",
            "source_column": resolved_column,
            "n_labeled": n_labeled,
            "n_unknown": n_unknown,
            "n_unique_types": len(unique_types),
            "unique_types_preview": unique_types[:20],
            "low_confidence_reassigned": low_confidence_reassigned,
            "rare_type_reassigned": rare_type_reassigned,
        }
    )
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Bridge selected scExtract annotations to downstream labels_for_scanvi.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--selected-label-source", required=True, choices=SOURCE_CHOICES)
    parser.add_argument("--source-column", default=None)
    parser.add_argument("--unknown-label", default="Unknown")
    parser.add_argument("--min-cells-per-type", type=int, default=10)
    parser.add_argument("--confidence-column", default=None)
    parser.add_argument("--min-confidence", type=float, default=None)
    parser.add_argument("--summary-json", default=None)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    adata = ad.read_h5ad(args.input_h5ad)
    summary = apply_bridge(
        adata,
        selected_label_source=args.selected_label_source,
        source_column=args.source_column,
        unknown_label=args.unknown_label,
        min_cells_per_type=args.min_cells_per_type,
        confidence_column=args.confidence_column,
        min_confidence=args.min_confidence,
    )
    output_path = Path(args.output_h5ad)
    ensure_dir(output_path.parent)
    adata.write_h5ad(output_path, compression="gzip")

    summary_path = Path(args.summary_json) if args.summary_json else output_path.with_suffix(".bridge_summary.json")
    write_json(summary, summary_path)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
