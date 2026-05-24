#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lineage_merge_schpl_core_20260402 import LineageMergeSchPLConfig, run_pipeline

PIPELINE_VERSION = "v1_0_branchwise_nocovid_20260523"
PIPELINE_TAG = "v1_0_branchwise_nocovid_20260523"
DEFAULT_REFERENCE_MANIFEST = Path(
    "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/branch_reference_manifest.json"
)
DEFAULT_QUERY_MANIFEST = Path(
    "/home/h2048/data/py/0407/stromal_scarches_query_branchwise_v1_0/branch_query_manifest.json"
)
DEFAULT_OUTPUT_DIR = Path(
    f"/home/h2048/data/py/0523/stromal_merge_schpl_{PIPELINE_TAG}"
)
BRANCH_MARKERS: dict[str, list[str]] = {
    "endothelial": ["CLDN5", "PECAM1", "VWF", "EMCN", "PROX1"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
    "smc": ["ACTA2", "TAGLN", "MYH11", "RGS5", "CSPG4"],
}


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    with open(path) as f:
        return json.load(f)


def _sanitize_object_cols(df: pd.DataFrame) -> None:
    for col in df.columns:
        s = df[col]
        if s.dtype != object:
            continue
        non_na = s.dropna()
        if len(non_na) == 0:
            df[col] = s.fillna("").astype(str)
            continue
        if non_na.map(lambda x: isinstance(x, (bool, np.bool_))).all():
            df[col] = s.fillna(False).astype(np.int8)
            continue
        if non_na.map(lambda x: not isinstance(x, str)).any():
            df[col] = s.map(lambda x: "" if pd.isna(x) else str(x))


def _ensure_counts_layer(adata: sc.AnnData) -> None:
    if "counts" not in adata.layers:
        adata.layers["counts"] = sparse.csr_matrix(adata.X)
    else:
        adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"])


def _build_log1p_layer(adata: sc.AnnData) -> None:
    tmp = sc.AnnData(X=adata.layers["counts"].copy(), obs=adata.obs.copy(), var=adata.var.copy())
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    adata.layers["log1p"] = sparse.csr_matrix(tmp.X)
    adata.X = adata.layers["log1p"]
    del tmp
    gc.collect()


def _merge_branch_reference_query(
    branch_name: str,
    ref_info: dict[str, Any],
    qry_info: dict[str, Any],
    branch_out_dir: Path,
) -> tuple[Path, dict[str, Any]]:
    ref_path = Path(ref_info["reference_h5ad"])
    qry_path = Path(qry_info["query_h5ad"])
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference h5ad not found for {branch_name}: {ref_path}")
    if not qry_path.exists():
        raise FileNotFoundError(f"Query h5ad not found for {branch_name}: {qry_path}")

    ad_ref = None
    ad_qry = None
    ad_ref = sc.read_h5ad(ref_path)
    ad_qry = sc.read_h5ad(qry_path)
    try:
        ad_ref.var_names_make_unique()
        ad_qry.var_names_make_unique()
        ad_ref.obs_names_make_unique()
        ad_qry.obs_names_make_unique()

        common_genes = [g for g in ad_ref.var_names.astype(str) if g in set(ad_qry.var_names.astype(str))]
        if len(common_genes) == 0:
            raise ValueError(f"0 shared genes for branch {branch_name}")

        ad_ref = ad_ref[:, common_genes].copy()
        ad_qry = ad_qry[:, common_genes].copy()
        _ensure_counts_layer(ad_ref)
        _ensure_counts_layer(ad_qry)

        ref_label_key = str(ref_info.get("label_key", "scanvi_label"))
        ref_latent_key = str(ref_info.get("latent_key", "X_scanvi"))
        qry_pred_key = str(qry_info.get("label_key", "cell_type_scarches_pred"))
        qry_final_key = str(qry_info.get("final_label_key", "cell_type_scarches_final"))
        qry_conf_key = str(qry_info.get("confidence_key", "scarches_confidence"))
        qry_latent_key = str(qry_info.get("latent_key", "X_scanvi"))

        if ref_latent_key not in ad_ref.obsm:
            raise KeyError(f"Reference latent key missing for {branch_name}: {ref_latent_key}")
        if qry_latent_key not in ad_qry.obsm:
            raise KeyError(f"Query latent key missing for {branch_name}: {qry_latent_key}")
        if ref_label_key not in ad_ref.obs.columns:
            raise KeyError(f"Reference label key missing for {branch_name}: {ref_label_key}")
        for key in [qry_pred_key, qry_final_key, qry_conf_key]:
            if key not in ad_qry.obs.columns:
                raise KeyError(f"Query key missing for {branch_name}: {key}")

        target_latent_key = "X_scanvi"
        if ref_latent_key != target_latent_key:
            ad_ref.obsm[target_latent_key] = np.asarray(ad_ref.obsm[ref_latent_key])
        if qry_latent_key != target_latent_key:
            ad_qry.obsm[target_latent_key] = np.asarray(ad_qry.obsm[qry_latent_key])

        ad_ref.obs_names = [f"ref::{x}" for x in ad_ref.obs_names]
        ad_qry.obs_names = [f"qry::{x}" for x in ad_qry.obs_names]

        ad_ref.obs["data_source"] = "reference"
        ad_qry.obs["data_source"] = "query"

        for key in ["sample", "tissue", ref_label_key, qry_pred_key, qry_final_key, qry_conf_key]:
            if key not in ad_ref.obs.columns:
                ad_ref.obs[key] = pd.NA
            if key not in ad_qry.obs.columns:
                ad_qry.obs[key] = pd.NA

        merged = sc.concat(
            [ad_ref, ad_qry],
            axis=0,
            join="inner",
            merge="same",
            uns_merge="first",
            index_unique=None,
        )
        merged.layers["counts"] = sparse.vstack(
            [ad_ref.layers["counts"], ad_qry.layers["counts"]],
            format="csr",
        )
        _build_log1p_layer(merged)

        sc.pp.neighbors(merged, use_rep=target_latent_key, n_neighbors=30, random_state=42)
        sc.tl.umap(merged, random_state=42)

        _sanitize_object_cols(merged.obs)
        _sanitize_object_cols(merged.var)

        branch_out_dir.mkdir(parents=True, exist_ok=True)
        merged_input = branch_out_dir / f"merged_input_{branch_name}_{PIPELINE_TAG}.h5ad"
        merged.write_h5ad(merged_input, compression="gzip", compression_opts=9)

        cfg_kwargs = {
            "lineage_name": f"Stromal {branch_name}",
            "lineage_slug": f"stromal_{branch_name}",
            "version": PIPELINE_VERSION,
            "input_merged_h5ad": str(merged_input),
            "output_dir": str(branch_out_dir),
            "reference_label_key": ref_label_key,
            "query_pred_key": qry_pred_key,
            "query_final_key": qry_final_key,
            "query_conf_key": qry_conf_key,
            "marker_genes": BRANCH_MARKERS.get(branch_name, []),
            "source_key": "data_source",
            "reference_source_value": "reference",
            "query_source_value": "query",
            "batch_key": "sample",
            "tissue_key": "tissue",
            "latent_key": target_latent_key,
            "umap_key": "X_umap",
            "marker_layer_preferred": "log1p",
            "unlabeled": "Unknown",
        }
        return merged_input, cfg_kwargs
    finally:
        if ad_ref is not None:
            del ad_ref
        if ad_qry is not None:
            del ad_qry
        gc.collect()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run stromal branchwise merged scHPL on branch reference/query manifests"
    )
    parser.add_argument("--reference-manifest", type=Path, default=DEFAULT_REFERENCE_MANIFEST)
    parser.add_argument("--query-manifest", type=Path, default=DEFAULT_QUERY_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--branches",
        nargs="+",
        default=["endothelial", "fibroblast", "smc"],
        choices=["endothelial", "fibroblast", "smc"],
    )
    parser.add_argument("--allow-overwrite", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    ref_manifest = _load_json(args.reference_manifest)
    qry_manifest = _load_json(args.query_manifest)
    ref_branches = ref_manifest.get("branches", {})
    qry_branches = qry_manifest.get("branches", {})

    summary_rows: list[dict[str, Any]] = []
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for branch_name in args.branches:
        if branch_name not in ref_branches:
            raise KeyError(f"Reference manifest missing branch: {branch_name}")
        if branch_name not in qry_branches:
            raise KeyError(f"Query manifest missing branch: {branch_name}")

        branch_out_dir = args.output_dir / branch_name
        existing_items = list(branch_out_dir.iterdir()) if branch_out_dir.exists() else []
        if existing_items and not args.allow_overwrite:
            raise FileExistsError(
                f"Branch output directory already exists and is not empty: {branch_out_dir}\n"
                "Pass --allow-overwrite if you intentionally want to reuse it."
            )

        merged_input, cfg_kwargs = _merge_branch_reference_query(
            branch_name=branch_name,
            ref_info=ref_branches[branch_name],
            qry_info=qry_branches[branch_name],
            branch_out_dir=branch_out_dir,
        )

        summary_row = {
            "branch": branch_name,
            "merged_input_h5ad": str(merged_input),
            "reference_h5ad": str(ref_branches[branch_name]["reference_h5ad"]),
            "query_h5ad": str(qry_branches[branch_name]["query_h5ad"]),
            "reference_label_key": cfg_kwargs["reference_label_key"],
            "query_pred_key": cfg_kwargs["query_pred_key"],
            "query_final_key": cfg_kwargs["query_final_key"],
            "query_conf_key": cfg_kwargs["query_conf_key"],
            "latent_key": cfg_kwargs["latent_key"],
            "umap_key": cfg_kwargs["umap_key"],
            "output_dir": str(branch_out_dir),
        }
        summary_rows.append(summary_row)

        print("=" * 80)
        print(f"Stromal branch preflight: {branch_name}")
        print("=" * 80)
        print(json.dumps(summary_row, indent=2))

        if args.preflight_only:
            continue

        cfg = LineageMergeSchPLConfig(**cfg_kwargs)
        run_pipeline(cfg)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.output_dir / "branchwise_schpl_inputs.csv", index=False)
    with open(args.output_dir / "branchwise_schpl_manifest.json", "w") as f:
        json.dump(summary_rows, f, indent=2)

    print("=" * 80)
    print("Stromal branchwise merge + scHPL complete")
    print("=" * 80)
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
