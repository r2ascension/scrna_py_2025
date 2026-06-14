#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Merge epithelial reference and no-COVID query outputs for downstream scHPL."""

from __future__ import annotations

import argparse
import gc
import json
import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lineage_merge_schpl_core_20260402 import sanitize_object_columns

DEFAULT_VERSION = "v1_0_nocovid_20260523"
DEFAULT_REF_H5AD = (
    PROJECT_ROOT
    / "data"
    / "py"
    / "0317"
    / "epithelial_v2_7_HOTFIX"
    / "REF"
    / "epithelial_scanvi_v2_7_HOTFIX_REF_final.h5ad"
)
DEFAULT_QUERY_H5AD = (
    PROJECT_ROOT
    / "data"
    / "py"
    / "0523"
    / "epithelial_scarches_query_v1_2_nocovid_20260523"
    / "epithelial_query_mapped_v1_2.h5ad"
)
DEFAULT_OUTPUT_PREFIX = "epithelial_merge_input"
DEFAULT_OUTPUT_FILENAME = "epithelial_reference_plus_query_merged_nocovid_20260523.h5ad"
REQUIRED_REF_OBS = ["cell_type_L3", "sample", "dataset", "tissue"]
REQUIRED_QRY_OBS = [
    "scanvi_fine_pred",
    "scanvi_fine_pred_filtered",
    "scanvi_fine_conf",
    "sample",
    "dataset",
    "tissue",
]
REQUIRED_SHARED_OBSM = ["X_scanvi_fine", "X_umap_fine"]
OPTIONAL_SHARED_OBSM = ["X_scanvi_major", "X_umap_major"]


def _has_nullable_strings(df: pd.DataFrame) -> bool:
    return any(str(dtype).startswith("string") for dtype in df.dtypes)


def derive_output_dir(version: str, output_prefix: str = DEFAULT_OUTPUT_PREFIX) -> Path:
    date_match = re.search(r"(20\d{6})$", version)
    day_dir = date_match.group(1)[4:] if date_match else "adhoc"
    return PROJECT_ROOT / "data" / "py" / day_dir / f"{output_prefix}_{version}"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Merge epithelial reference and query for scHPL input")
    parser.add_argument("--reference-h5ad", default=str(DEFAULT_REF_H5AD), help="Reference h5ad path")
    parser.add_argument("--query-h5ad", default=str(DEFAULT_QUERY_H5AD), help="Mapped query h5ad path")
    parser.add_argument("--version", default=DEFAULT_VERSION, help="Version tag used in output dir and manifest")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional explicit output directory. Defaults to data/py/<MMDD>/epithelial_merge_input_<version>",
    )
    parser.add_argument("--output-h5ad", default=DEFAULT_OUTPUT_FILENAME, help="Merged output filename inside output dir")
    parser.add_argument("--allow-overwrite", action="store_true", help="Allow writing into a non-empty output directory")
    parser.add_argument("--preflight-only", action="store_true", help="Validate inputs and report merge plan without writing")
    return parser


def _ensure_keys(adata: sc.AnnData, required_obs: list[str], required_obsm: list[str], *, adata_name: str) -> None:
    missing_obs = [key for key in required_obs if key not in adata.obs.columns]
    missing_obsm = [key for key in required_obsm if key not in adata.obsm.keys()]
    if missing_obs or missing_obsm:
        raise KeyError(f"[{adata_name}] missing keys: obs={missing_obs}, obsm={missing_obsm}")


def _to_csr(matrix):
    return csr_matrix(matrix) if issparse(matrix) else csr_matrix(matrix)


def _keep_counts_only(adata: sc.AnnData, *, adata_name: str) -> None:
    if "counts" in adata.layers:
        counts = adata.layers["counts"]
    else:
        counts = adata.X
        print(f"[WARN] {adata_name} is missing a 'counts' layer; falling back to X")
    counts = _to_csr(counts)
    adata.layers.clear()
    adata.layers["counts"] = counts
    adata.X = counts.copy()


def _prune_obsm(adata: sc.AnnData, keep_keys: list[str]) -> None:
    for key in list(adata.obsm.keys()):
        if key not in keep_keys:
            del adata.obsm[key]


def validate_output_dir(path: Path, allow_overwrite: bool) -> None:
    if not path.exists():
        return
    existing = [p.name for p in path.iterdir()]
    if existing and not allow_overwrite:
        raise FileExistsError(
            f"Output directory already exists and is not empty: {path}\n"
            "Pass --allow-overwrite only if you intentionally want to reuse it."
        )


def main() -> None:
    args = build_arg_parser().parse_args()
    ref_path = Path(args.reference_h5ad).expanduser().resolve()
    qry_path = Path(args.query_h5ad).expanduser().resolve()
    out_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else derive_output_dir(args.version).resolve()
    out_path = out_dir / args.output_h5ad

    if not ref_path.exists():
        raise FileNotFoundError(f"Reference h5ad not found: {ref_path}")
    if not qry_path.exists():
        raise FileNotFoundError(f"Query h5ad not found: {qry_path}")
    validate_output_dir(out_dir, allow_overwrite=args.allow_overwrite)

    print("=" * 80)
    print("Epithelial reference + no-COVID query merge")
    print("=" * 80)
    print(f"Reference: {ref_path}")
    print(f"Query    : {qry_path}")
    print(f"Output   : {out_path}")
    print(f"Version  : {args.version}")

    print("\n[1/5] Loading reference and query")
    adata_ref = sc.read_h5ad(ref_path)
    adata_qry = sc.read_h5ad(qry_path)
    adata_ref.var_names_make_unique()
    adata_qry.var_names_make_unique()
    adata_ref.obs_names_make_unique()
    adata_qry.obs_names_make_unique()

    _ensure_keys(adata_ref, REQUIRED_REF_OBS, REQUIRED_SHARED_OBSM, adata_name="reference")
    _ensure_keys(adata_qry, REQUIRED_QRY_OBS, REQUIRED_SHARED_OBSM, adata_name="query")

    optional_shared = [key for key in OPTIONAL_SHARED_OBSM if key in adata_ref.obsm and key in adata_qry.obsm]
    keep_obsm = REQUIRED_SHARED_OBSM + optional_shared
    common_genes = [gene for gene in adata_ref.var_names if gene in adata_qry.var_names]
    if not common_genes:
        raise ValueError("Reference and query do not share any genes")

    print(f"  Reference shape : {adata_ref.shape}")
    print(f"  Query shape     : {adata_qry.shape}")
    print(f"  Common genes    : {len(common_genes):,}/{adata_ref.n_vars:,} ({len(common_genes) / max(adata_ref.n_vars, 1):.1%})")
    print(f"  Shared obsm     : {keep_obsm}")

    if args.preflight_only:
        print("[OK] Preflight complete; --preflight-only set, skipping write.")
        return

    print("\n[2/5] Harmonizing genes and storage")
    adata_ref = adata_ref[:, common_genes].copy()
    adata_qry = adata_qry[:, common_genes].copy()

    if "cell_type_L3" not in adata_qry.obs.columns:
        adata_qry.obs["cell_type_L3"] = pd.NA
    if "scanvi_fine_pred" not in adata_ref.obs.columns:
        adata_ref.obs["scanvi_fine_pred"] = pd.NA
    if "scanvi_fine_pred_filtered" not in adata_ref.obs.columns:
        adata_ref.obs["scanvi_fine_pred_filtered"] = pd.NA
    if "scanvi_fine_conf" not in adata_ref.obs.columns:
        adata_ref.obs["scanvi_fine_conf"] = float("nan")

    _keep_counts_only(adata_ref, adata_name="reference")
    _keep_counts_only(adata_qry, adata_name="query")
    _prune_obsm(adata_ref, keep_obsm)
    _prune_obsm(adata_qry, keep_obsm)
    adata_ref.uns = {}
    adata_qry.uns = {}
    for key in list(adata_ref.obsp.keys()):
        del adata_ref.obsp[key]
    for key in list(adata_qry.obsp.keys()):
        del adata_qry.obsp[key]

    print("\n[3/5] Prefixing obs names to preserve provenance")
    adata_ref.obs_names = [f"ref::{x}" for x in adata_ref.obs_names]
    adata_qry.obs_names = [f"qry::{x}" for x in adata_qry.obs_names]

    print("\n[4/5] Concatenating")
    adata_merged = sc.concat(
        {"reference": adata_ref, "query": adata_qry},
        axis=0,
        join="inner",
        merge="unique",
        label="data_source",
    )
    sanitize_object_columns(adata_merged.obs, "obs")
    sanitize_object_columns(adata_merged.var, "var")

    expected_n = adata_ref.n_obs + adata_qry.n_obs
    if adata_merged.n_obs != expected_n:
        raise ValueError(f"Merged n_obs mismatch: {adata_merged.n_obs} != {expected_n}")
    for key in keep_obsm:
        if key not in adata_merged.obsm:
            raise KeyError(f"Merged output lost obsm['{key}']")
        if adata_merged.obsm[key].shape[0] != adata_merged.n_obs:
            raise ValueError(f"Merged obsm['{key}'] row mismatch: {adata_merged.obsm[key].shape[0]} != {adata_merged.n_obs}")

    manifest = {
        "version": args.version,
        "reference_h5ad": str(ref_path),
        "query_h5ad": str(qry_path),
        "output_h5ad": str(out_path),
        "common_genes": len(common_genes),
        "reference_cells": int(adata_ref.n_obs),
        "query_cells": int(adata_qry.n_obs),
        "latent_key": "X_scanvi_fine",
        "umap_key": "X_umap_fine",
        "kept_obsm": keep_obsm,
        "reference_label_key": "cell_type_L3",
        "query_pred_key": "scanvi_fine_pred",
        "query_final_key": "scanvi_fine_pred_filtered",
        "query_conf_key": "scanvi_fine_conf",
        "counts_layer": "counts",
    }
    adata_merged.uns["merge_manifest"] = manifest

    print("\n[5/5] Writing outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    if (
        hasattr(ad, "settings")
        and hasattr(ad.settings, "allow_write_nullable_strings")
        and (_has_nullable_strings(adata_merged.obs) or _has_nullable_strings(adata_merged.var))
        and not ad.settings.allow_write_nullable_strings
    ):
        print(
            "[INFO] Detected pandas StringDtype columns in merged AnnData; "
            "enabling anndata.settings.allow_write_nullable_strings=True before write"
        )
        ad.settings.allow_write_nullable_strings = True
    try:
        adata_merged.write_h5ad(out_path, compression="gzip")
    except RuntimeError as e:
        err_msg = str(e)
        if "allow_write_nullable_strings" not in err_msg:
            raise
        if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
            print(
                "[WARN] anndata refused nullable string write; enabling "
                "anndata.settings.allow_write_nullable_strings=True and retrying"
            )
            ad.settings.allow_write_nullable_strings = True
            adata_merged.write_h5ad(out_path, compression="gzip")
        else:
            raise
    manifest_path = out_dir / f"merge_manifest_{args.version}.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print(f"[OK] merged h5ad   : {out_path}")
    print(f"[OK] merge manifest: {manifest_path}")

    del adata_ref, adata_qry, adata_merged
    gc.collect()

    print("\n" + "=" * 80)
    print("Merge complete")
    print("=" * 80)


if __name__ == "__main__":
    main()
