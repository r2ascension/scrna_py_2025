#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path

import numpy as np
import scanpy as sc

CURRENT_MERGED = Path(
    "/home/h2048/data/py/0330/bcell_merge_schpl_v3_0/reference_plus_query_merged_v3_0.h5ad"
)
SCANVI_SOURCE = Path(
    "/home/h2048/data/py/0208/merged_scanvi_L2_prod_v1/merged_scanvi_L2_prod.h5ad"
)
OUTPUT = Path(
    "/home/h2048/data/py/0413/bcell_merge_schpl_scanvi_umap_input_v1_0/reference_plus_query_merged_scanvi_umap_v1_0.h5ad"
)


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    print(f"[load] current merged: {CURRENT_MERGED}")
    cur = sc.read_h5ad(CURRENT_MERGED)
    print(f"[load] scanvi source : {SCANVI_SOURCE}")
    src = sc.read_h5ad(SCANVI_SOURCE)

    if "barcode" not in cur.obs.columns or "barcode" not in src.obs.columns:
        raise KeyError("Both input h5ad files must contain obs['barcode'] for alignment")
    if "X_umap_scanvi" not in src.obsm:
        raise KeyError("Historical scanvi source is missing obsm['X_umap_scanvi']")

    cur_barcodes = cur.obs["barcode"].astype(str).tolist()
    src_barcodes = src.obs["barcode"].astype(str).tolist()

    if len(set(cur_barcodes)) != len(cur_barcodes):
        raise ValueError("Current merged input contains duplicated barcodes")
    if len(set(src_barcodes)) != len(src_barcodes):
        raise ValueError("Historical scanvi source contains duplicated barcodes")

    src_index = {barcode: idx for idx, barcode in enumerate(src_barcodes)}
    missing = [barcode for barcode in cur_barcodes if barcode not in src_index]
    if missing:
        raise ValueError(
            f"Missing {len(missing)} barcodes in scanvi source; first few: {missing[:10]}"
        )

    src_positions = np.array([src_index[barcode] for barcode in cur_barcodes], dtype=int)
    scanvi_umap = np.asarray(src.obsm["X_umap_scanvi"])[src_positions]

    if scanvi_umap.shape[0] != cur.n_obs:
        raise AssertionError("Aligned scanvi UMAP row count does not match current merged input")

    cur.uns["umap_replaced_from"] = str(SCANVI_SOURCE)
    cur.uns["umap_replaced_key"] = "X_umap_scanvi"
    cur.uns["umap_replaced_note"] = (
        "For 2026-04-13 B-cell schpl rerun, X_umap was replaced with historical scanvi UMAP "
        "so all downstream plots consistently use scanvi coordinates."
    )

    if "X_umap" in cur.obsm:
        cur.obsm["X_umap_legacy_pre_scanvi_swap"] = np.asarray(cur.obsm["X_umap"]).copy()
    cur.obsm["X_umap_scanvi"] = scanvi_umap.copy()
    cur.obsm["X_umap"] = scanvi_umap.copy()

    print(f"[write] {OUTPUT}")
    cur.write_h5ad(OUTPUT, compression="gzip")
    print("[done] scanvi UMAP injected successfully")
    print(f"[summary] n_obs={cur.n_obs:,}, n_vars={cur.n_vars:,}")


if __name__ == "__main__":
    main()
