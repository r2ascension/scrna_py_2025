#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

DEFAULT_INPUT_H5AD = Path(
    "/home/h2048/data/py/20260531/main_lineage_h5ad_cleanup_20260531/epithelial/"
    "epithelial_tissue_comparison_final_fullgene_cleaned_20260531.h5ad"
)
DEFAULT_OUTPUT_H5AD = Path(
    "/home/h2048/data/py/20260531/main_lineage_h5ad_cleanup_20260531/epithelial/"
    "epithelial_tissue_comparison_final_fullgene_cleaned_l3refined_20260531.h5ad"
)
DEFAULT_SUMMARY_JSON = Path(
    "/home/h2048/data/py/20260531/main_lineage_h5ad_cleanup_20260531/epithelial/"
    "epithelial_l3_refine_summary_20260531.json"
)
DEFAULT_CHANGED_CELLS_TSV = Path(
    "/home/h2048/data/py/20260531/main_lineage_h5ad_cleanup_20260531/epithelial/"
    "epithelial_l3_refine_changed_cells_20260531.tsv"
)
DEFAULT_TARGET_COLUMNS = (
    "cell_type_scanvi_pred",
    "cell_type_L3",
    "cell_type_L3_curated",
    "scanvi_fine_pred",
    "scanvi_labels_fine",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Refine epithelial L3 labels in a copied h5ad without loading the full matrix into memory.")
    parser.add_argument("--input-h5ad", default=str(DEFAULT_INPUT_H5AD))
    parser.add_argument("--output-h5ad", default=str(DEFAULT_OUTPUT_H5AD))
    parser.add_argument("--summary-json", default=str(DEFAULT_SUMMARY_JSON))
    parser.add_argument("--changed-cells-tsv", default=str(DEFAULT_CHANGED_CELLS_TSV))
    parser.add_argument("--neuro-source-column", default="celltypist_pred_original")
    parser.add_argument("--neuro-source-label", default="Neuroendocrine")
    parser.add_argument("--neuro-target-label", default="Neuroendocrine")
    parser.add_argument(
        "--target-columns",
        default=",".join(DEFAULT_TARGET_COLUMNS),
        help="Comma-separated L3-like columns that should receive the Neuroendocrine relabel and Ionocyte_Brush rename.",
    )
    parser.add_argument("--rename-old", default="Ionocyte_Brush")
    parser.add_argument("--rename-new", default="Ionocyte")
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def normalize_string_series(series: pd.Series) -> pd.Series:
    normalized = series.astype("string")
    return normalized.str.strip()


def build_categories(original: pd.Series, updated: pd.Series, rename_old: str, rename_new: str, neuro_target_label: str) -> list[str]:
    original_categories: list[str] = []
    if isinstance(original.dtype, CategoricalDtype):
        original_categories = [str(value) for value in original.cat.categories.tolist()]
    else:
        original_categories = [str(value) for value in pd.unique(normalize_string_series(original).dropna())]

    ordered: list[str] = []
    for value in original_categories:
        mapped = rename_new if value == rename_old else value
        if mapped not in ordered:
            ordered.append(mapped)

    used_values = [str(value) for value in pd.unique(updated.dropna())]
    if neuro_target_label in used_values and neuro_target_label not in ordered:
        if rename_new in ordered:
            ordered.insert(ordered.index(rename_new) + 1, neuro_target_label)
        else:
            ordered.append(neuro_target_label)
    for value in used_values:
        if value not in ordered:
            ordered.append(value)

    used_set = set(used_values)
    return [value for value in ordered if value in used_set]


def build_updated_columns(
    obs: pd.DataFrame,
    target_columns: Iterable[str],
    neuro_source_column: str,
    neuro_source_label: str,
    neuro_target_label: str,
    rename_old: str,
    rename_new: str,
) -> tuple[dict[str, pd.Series], pd.Series, pd.DataFrame, dict[str, object]]:
    if neuro_source_column not in obs.columns:
        raise KeyError(f"Missing neuro source column: {neuro_source_column}")

    neuro_source = normalize_string_series(obs[neuro_source_column])
    neuro_mask = neuro_source.fillna("") == neuro_source_label

    updated_columns: dict[str, pd.Series] = {}
    summary_columns: dict[str, object] = {}
    changed_df = pd.DataFrame(index=obs.index)
    changed_df["neuro_source_value"] = neuro_source
    changed_df["neuro_source_match"] = neuro_mask.to_numpy()

    for column in target_columns:
        if column not in obs.columns:
            continue
        original = normalize_string_series(obs[column])
        updated = original.mask(original.fillna("") == rename_old, rename_new)
        updated = updated.mask(neuro_mask, neuro_target_label)
        updated_columns[column] = updated

        before_counts = original.value_counts(dropna=False)
        after_counts = updated.value_counts(dropna=False)
        changed_mask = (~original.fillna("<NA>").eq(updated.fillna("<NA>"))).to_numpy()
        changed_df[f"{column}__before"] = original
        changed_df[f"{column}__after"] = updated
        changed_df[f"{column}__changed"] = changed_mask

        summary_columns[column] = {
            "n_changed": int(changed_mask.sum()),
            "n_neuro_target": int((updated.fillna("") == neuro_target_label).sum()),
            "n_rename_target": int((updated.fillna("") == rename_new).sum()),
            "before": {str(k): int(v) for k, v in before_counts.items() if int(v) > 0},
            "after": {str(k): int(v) for k, v in after_counts.items() if int(v) > 0},
            "categories": build_categories(original, updated, rename_old=rename_old, rename_new=rename_new, neuro_target_label=neuro_target_label),
        }

    row_change_cols = [col for col in changed_df.columns if col.endswith("__changed")]
    changed_rows = changed_df.loc[changed_df[row_change_cols].any(axis=1)].copy() if row_change_cols else changed_df.iloc[0:0].copy()
    return updated_columns, neuro_mask, changed_rows, summary_columns


def choose_codes_dtype(n_categories: int) -> np.dtype:
    if n_categories < 127:
        return np.int8
    if n_categories < 32767:
        return np.int16
    return np.int32


def write_categorical_column(file_handle: h5py.File, column: str, updated: pd.Series, categories: list[str]) -> None:
    group = file_handle["obs"][column]
    group_attrs = dict(group.attrs)
    cat_attrs = dict(group["categories"].attrs)
    code_attrs = dict(group["codes"].attrs)

    categorical = pd.Categorical(updated, categories=categories, ordered=bool(group_attrs.get("ordered", False)))
    codes = categorical.codes.astype(choose_codes_dtype(len(categories)), copy=False)
    categories_array = np.asarray(categories, dtype=object)

    del group["categories"]
    del group["codes"]

    cat_ds = group.create_dataset("categories", data=categories_array, dtype=h5py.string_dtype(encoding="utf-8"))
    for key, value in cat_attrs.items():
        cat_ds.attrs[key] = value

    codes_ds = group.create_dataset("codes", data=codes, dtype=codes.dtype)
    for key, value in code_attrs.items():
        codes_ds.attrs[key] = value

    for key, value in group_attrs.items():
        group.attrs[key] = value


def main() -> None:
    args = parse_args()
    input_h5ad = Path(args.input_h5ad)
    output_h5ad = Path(args.output_h5ad)
    summary_json = Path(args.summary_json)
    changed_cells_tsv = Path(args.changed_cells_tsv)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    changed_cells_tsv.parent.mkdir(parents=True, exist_ok=True)

    target_columns = split_csv(args.target_columns)

    backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        obs_columns = [args.neuro_source_column, *target_columns]
        obs_columns = [column for column in obs_columns if column in backed.obs.columns]
        obs = backed.obs[obs_columns].copy()
        n_obs = int(backed.n_obs)
        n_vars = int(backed.n_vars)
    finally:
        file_manager = getattr(backed, "file", None)
        if file_manager is not None:
            file_manager.close()

    updated_columns, neuro_mask, changed_rows, summary_columns = build_updated_columns(
        obs=obs,
        target_columns=target_columns,
        neuro_source_column=args.neuro_source_column,
        neuro_source_label=args.neuro_source_label,
        neuro_target_label=args.neuro_target_label,
        rename_old=args.rename_old,
        rename_new=args.rename_new,
    )

    shutil.copy2(input_h5ad, output_h5ad)
    with h5py.File(output_h5ad, "r+") as handle:
        for column, updated in updated_columns.items():
            categories = list(summary_columns[column]["categories"])
            write_categorical_column(handle, column=column, updated=updated, categories=categories)

    changed_rows = changed_rows.reset_index().rename(columns={changed_rows.index.name or "index": "obs_name"})
    changed_rows.to_csv(changed_cells_tsv, sep="\t", index=False)

    summary = {
        "input_h5ad": str(input_h5ad),
        "output_h5ad": str(output_h5ad),
        "summary_json": str(summary_json),
        "changed_cells_tsv": str(changed_cells_tsv),
        "n_obs": n_obs,
        "n_vars": n_vars,
        "neuro_source_column": args.neuro_source_column,
        "neuro_source_label": args.neuro_source_label,
        "neuro_target_label": args.neuro_target_label,
        "rename_old": args.rename_old,
        "rename_new": args.rename_new,
        "n_neuro_source_cells": int(neuro_mask.sum()),
        "n_changed_rows": int(changed_rows.shape[0]),
        "target_columns": target_columns,
        "columns": summary_columns,
    }
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[refine] wrote changed cells: {changed_cells_tsv}")
    print(f"[refine] wrote summary: {summary_json}")


if __name__ == "__main__":
    main()
