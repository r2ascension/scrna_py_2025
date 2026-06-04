#!/usr/bin/env python3
"""Transfer merged/scHPL B-cell subtype labels back to the full-gene query object.

Repository-validated join contract for this pilot:
- full-gene query primary join key: `obs_names`
- merged/scHPL query primary join key: `obs_names` stripped from `qry::qry::` prefix
- `CellID` is *not* the correct bridge key here
- `cell_id` is mostly `Unknown` in the current query/merged assets
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

PY_ROOT = Path(__file__).resolve().parents[1]
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from normal_airway_ml_common_20260527 import deep_get, ensure_dir, load_config, write_json  # noqa: E402

CONFIG_KEY = "bcell_internal_site_core_validation"
DEFAULT_CONFIG_PATH = "/home/h2048/script/config/normal_airway_ml_bcell_internal_site_core_validation_20260531.yaml"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Transfer B-cell subtype labels from merged/scHPL query cells back to the full-gene query h5ad.")
    parser.add_argument("--config", type=Path, default=Path(DEFAULT_CONFIG_PATH), help="YAML config path.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional explicit run root override.")
    return parser.parse_args()


def resolve_run_root(cfg: dict[str, Any], output_dir: Path | None) -> Path:
    if output_dir is not None:
        return ensure_dir(output_dir)
    output_root = Path(str(deep_get(cfg, "run", "output_root", default="/home/h2048/data/py/20260531")))
    prefix = str(deep_get(cfg, "run", "run_name_prefix", default="bcell_internal_site_core_validation_20260531"))
    return ensure_dir(output_root / prefix)


def pick_first_existing(columns: list[str], candidates: list[str]) -> str:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise KeyError(f"None of the candidate columns exist: {candidates}")


def coerce_bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    text = values.astype("string").fillna("False").astype(str).str.strip().str.lower()
    return text.isin({"true", "1", "yes", "y", "t"})


def build_qc_table(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    total = int(df.shape[0])
    rows.extend(
        [
            {"section": "overall", "metric": "n_query_cells", "subtype_label": "", "rejected": "", "keep_for_validation": "", "n_cells": total, "fraction": 1.0},
            {"section": "overall", "metric": "n_joined_cells", "subtype_label": "", "rejected": "", "keep_for_validation": "", "n_cells": int(df["label_join_found"].sum()), "fraction": float(df["label_join_found"].mean())},
            {"section": "overall", "metric": "n_rejected_cells", "subtype_label": "", "rejected": "", "keep_for_validation": "", "n_cells": int(df["schpl_rejected"].sum()), "fraction": float(df["schpl_rejected"].mean())},
            {"section": "overall", "metric": "n_kept_for_validation", "subtype_label": "", "rejected": "", "keep_for_validation": "", "n_cells": int(df["keep_for_validation"].sum()), "fraction": float(df["keep_for_validation"].mean())},
        ]
    )
    subtype_counts = (
        df.groupby(["transferred_subtype_label", "schpl_rejected", "keep_for_validation"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    for row in subtype_counts.itertuples(index=False):
        rows.append(
            {
                "section": "subtype_counts",
                "metric": "cells_by_label",
                "subtype_label": str(row.transferred_subtype_label),
                "rejected": str(bool(row.schpl_rejected)),
                "keep_for_validation": str(bool(row.keep_for_validation)),
                "n_cells": int(row.n_cells),
                "fraction": float(row.n_cells / total) if total else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    run_root = resolve_run_root(cfg, args.output_dir)
    label_dir = ensure_dir(run_root / "labels")

    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    labels_cfg = deep_get(block, "labels", default={}) or {}
    join_cfg = deep_get(block, "query_join", default={}) or {}

    query_path = Path(str(deep_get(block, "disease_query_h5ad", default="")))
    merged_path = Path(str(deep_get(block, "merged_query_h5ad", default="")))
    merged_prefix = str(deep_get(join_cfg, "merged_query_prefix", default="qry::qry::"))

    subtype_col_candidates = [str(x) for x in deep_get(labels_cfg, "subtype_column_candidates", default=["viz__qry_label_final_only", "Cell_Type_L2_final", "Cell_Type_L2_pred"])]
    rejection_col_candidates = [str(x) for x in deep_get(labels_cfg, "rejection_column_candidates", default=["schpl_rejected"])]
    subtype_whitelist = {str(x) for x in deep_get(labels_cfg, "subtype_whitelist", default=["Naive_B", "Memory_B", "Plasma"])}
    exclude_labels = {str(x) for x in deep_get(labels_cfg, "exclude_labels", default=["Unknown", "Atypical_Memory_B", "GC_B"])}

    qa = ad.read_h5ad(query_path, backed="r")
    ma = ad.read_h5ad(merged_path, backed="r")
    try:
        qobs = qa.obs.copy()
        mobs = ma.obs.copy()
    finally:
        qa.file.close()
        ma.file.close()

    q_df = qobs.copy()
    q_df.insert(0, "obs_name", q_df.index.astype(str))
    q_df["join_id"] = q_df["obs_name"].astype(str)

    m_df = mobs.copy()
    m_df.insert(0, "merged_obs_name", m_df.index.astype(str))
    query_mask = ~m_df["merged_obs_name"].str.startswith("ref::")
    m_query = m_df.loc[query_mask].copy()
    if merged_prefix:
        m_query["join_id"] = m_query["merged_obs_name"].str.replace("^" + re.escape(merged_prefix), "", regex=True)
    else:
        m_query["join_id"] = m_query["merged_obs_name"].astype(str)

    subtype_col = pick_first_existing(list(m_query.columns), subtype_col_candidates)
    rejection_col = pick_first_existing(list(m_query.columns), rejection_col_candidates)

    if m_query["join_id"].duplicated().any():
        dupes = m_query.loc[m_query["join_id"].duplicated(keep=False), "join_id"].astype(str).head(10).tolist()
        raise ValueError(f"Merged query join_id is not unique; example duplicates: {dupes}")

    bridge = m_query[["join_id", subtype_col, rejection_col, "merged_obs_name"]].copy()
    bridge = bridge.rename(columns={subtype_col: "transferred_subtype_label", rejection_col: "schpl_rejected_raw"})
    bridge["transferred_subtype_label"] = bridge["transferred_subtype_label"].astype("string").fillna("Unknown").astype(str)
    bridge["schpl_rejected"] = coerce_bool(bridge["schpl_rejected_raw"])

    merged = q_df.merge(bridge[["join_id", "transferred_subtype_label", "schpl_rejected", "merged_obs_name"]], on="join_id", how="left")
    merged["label_join_found"] = merged["transferred_subtype_label"].notna()
    merged["transferred_subtype_label"] = merged["transferred_subtype_label"].fillna("Unknown").astype(str)
    merged["schpl_rejected"] = merged["schpl_rejected"].fillna(False).astype(bool)
    merged["subtype_whitelisted"] = merged["transferred_subtype_label"].isin(subtype_whitelist)
    merged["subtype_excluded"] = merged["transferred_subtype_label"].isin(exclude_labels)
    merged["keep_for_validation"] = merged["label_join_found"] & (~merged["schpl_rejected"]) & merged["subtype_whitelisted"] & (~merged["subtype_excluded"])
    merged["exclusion_reason"] = np.select(
        [
            ~merged["label_join_found"],
            merged["schpl_rejected"],
            ~merged["subtype_whitelisted"],
            merged["subtype_excluded"],
        ],
        [
            "missing_label_join",
            "schpl_rejected",
            "subtype_not_whitelisted",
            "subtype_explicitly_excluded",
        ],
        default="kept",
    )

    keep_cols_front = [
        "obs_name",
        "join_id",
        "merged_obs_name",
        "transferred_subtype_label",
        "schpl_rejected",
        "label_join_found",
        "subtype_whitelisted",
        "subtype_excluded",
        "keep_for_validation",
        "exclusion_reason",
    ]
    other_cols = [c for c in merged.columns if c not in keep_cols_front]
    merged = merged[keep_cols_front + other_cols].copy()

    transfer_path = label_dir / "query_label_transfer.tsv"
    qc_path = label_dir / "query_label_transfer_qc.tsv"
    merged.to_csv(transfer_path, sep="\t", index=False)
    qc_df = build_qc_table(merged)
    qc_df.to_csv(qc_path, sep="\t", index=False)

    manifest = {
        "config_path": str(args.config),
        "run_root": str(run_root),
        "query_path": str(query_path),
        "merged_path": str(merged_path),
        "selected_columns": {
            "subtype_column": subtype_col,
            "rejection_column": rejection_col,
            "query_join_key": "obs_names",
            "merged_join_key": "obs_names_without_prefix",
        },
        "files": {
            "query_label_transfer": str(transfer_path),
            "query_label_transfer_qc": str(qc_path),
        },
        "counts": {
            "n_query_cells": int(merged.shape[0]),
            "n_joined_cells": int(merged["label_join_found"].sum()),
            "n_kept_for_validation": int(merged["keep_for_validation"].sum()),
        },
    }
    write_json(manifest, label_dir / "label_transfer_manifest.json")

    print(json.dumps(manifest["counts"], ensure_ascii=False, indent=2))
    print(f"[ok] query label transfer -> {transfer_path}")
    print(f"[ok] label qc           -> {qc_path}")


if __name__ == "__main__":
    main()
