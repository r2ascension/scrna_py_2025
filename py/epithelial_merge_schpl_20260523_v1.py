#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Thin, production-safe epithelial scHPL entrypoint for the no-COVID rerun."""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import scanpy as sc

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lineage_merge_schpl_core_20260402 import LineageMergeSchPLConfig, run_pipeline

DEFAULT_VERSION = "v1_0_nocovid_20260523"
DEFAULT_INPUT_MERGED_H5AD = (
    PROJECT_ROOT
    / "data"
    / "py"
    / "0523"
    / "epithelial_merge_input_v1_0_nocovid_20260523"
    / "epithelial_reference_plus_query_merged_nocovid_20260523.h5ad"
)
DEFAULT_OUTPUT_PREFIX = "epithelial_merge_schpl"
DEFAULT_MARKERS = [
    "EPCAM",
    "KRT5",
    "KRT8",
    "KRT17",
    "FOXJ1",
    "SCGB1A1",
    "MUC1",
    "MUC5AC",
]


def derive_output_dir(version: str, output_prefix: str = DEFAULT_OUTPUT_PREFIX) -> Path:
    date_match = re.search(r"(20\d{6})$", version)
    day_dir = date_match.group(1)[4:] if date_match else "adhoc"
    return PROJECT_ROOT / "data" / "py" / day_dir / f"{output_prefix}_{version}"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Hardened epithelial scHPL entrypoint")
    parser.add_argument(
        "--input-h5ad",
        default=str(DEFAULT_INPUT_MERGED_H5AD),
        help="Merged reference+query h5ad used as scHPL input.",
    )
    parser.add_argument(
        "--version",
        default=DEFAULT_VERSION,
        help="Run version tag; also used to derive the default output directory.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional explicit output directory. By default this is derived from --version.",
    )
    parser.add_argument(
        "--allow-overwrite",
        action="store_true",
        help="Allow running into a non-empty output directory. Off by default for safety.",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Validate config and input structure, then exit without running the full pipeline.",
    )
    return parser


def build_config(args: argparse.Namespace) -> LineageMergeSchPLConfig:
    input_h5ad = Path(args.input_h5ad).expanduser().resolve()
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else derive_output_dir(args.version).resolve()
    )
    return LineageMergeSchPLConfig(
        lineage_name="Epithelial",
        lineage_slug="epithelial",
        version=args.version,
        input_merged_h5ad=str(input_h5ad),
        output_dir=str(output_dir),
        reference_label_key="cell_type_L3",
        query_pred_key="scanvi_fine_pred",
        query_final_key="scanvi_fine_pred_filtered",
        query_conf_key="scanvi_fine_conf",
        marker_genes=list(DEFAULT_MARKERS),
        source_key="data_source",
        reference_source_value="reference",
        query_source_value="query",
        batch_key="sample",
        tissue_key="tissue",
        latent_key="X_scanvi_fine",
        umap_key="X_umap_fine",
        marker_layer_preferred="counts",
        unlabeled="Unknown",
        schpl_classifier="knn",
        schpl_dimred=True,
        schpl_use_re=True,
        schpl_n_neighbors=50,
        schpl_fn=0.5,
        schpl_rej_threshold=0.5,
        schpl_qc_n_neighbors=30,
        schpl_qc_leiden_resolution=0.5,
        schpl_novel_rej_threshold=0.50,
        schpl_novel_min_cells=50,
        random_seed=42,
        dpi=300,
        figure_format="pdf",
    )


def validate_config(cfg: LineageMergeSchPLConfig, allow_overwrite: bool = False) -> None:
    input_path = Path(cfg.input_merged_h5ad)
    output_dir = Path(cfg.output_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_path}")
    if input_path.suffix != ".h5ad":
        raise ValueError(f"Input must be an .h5ad file, got: {input_path}")

    if cfg.reference_source_value == cfg.query_source_value:
        raise ValueError("reference_source_value and query_source_value must differ")
    if cfg.query_pred_key == cfg.query_final_key:
        raise ValueError("query_pred_key and query_final_key must differ")

    if output_dir.exists():
        existing_items = [p.name for p in output_dir.iterdir()]
        if existing_items and not allow_overwrite:
            raise FileExistsError(
                f"Output directory already exists and is not empty: {output_dir}\n"
                "Pass --allow-overwrite only if you intentionally want to reuse it."
            )
        if existing_items and allow_overwrite:
            warnings.warn(
                f"Output directory is non-empty and will be reused because --allow-overwrite was set: {output_dir}",
                stacklevel=2,
            )

    adata = None
    try:
        adata = sc.read_h5ad(input_path, backed="r")
        required_obs = [
            cfg.source_key,
            cfg.batch_key,
            cfg.tissue_key,
            cfg.reference_label_key,
            cfg.query_pred_key,
            cfg.query_final_key,
            cfg.query_conf_key,
        ]
        missing_obs = [key for key in required_obs if key not in adata.obs.columns]
        if missing_obs:
            raise KeyError(f"Missing obs keys in input h5ad: {missing_obs}")

        required_obsm = [cfg.latent_key, cfg.umap_key]
        missing_obsm = [key for key in required_obsm if key not in adata.obsm.keys()]
        if missing_obsm:
            raise KeyError(f"Missing obsm keys in input h5ad: {missing_obsm}")

        source_values = {str(x) for x in adata.obs[cfg.source_key].astype(str).unique()}
        expected_sources = {cfg.reference_source_value, cfg.query_source_value}
        if source_values != expected_sources:
            raise ValueError(f"{cfg.source_key} must contain exactly {sorted(expected_sources)}, got {sorted(source_values)}")

        src = adata.obs[cfg.source_key].astype(str)
        ref_mask = src.eq(cfg.reference_source_value)
        qry_mask = src.eq(cfg.query_source_value)
        n_ref = int(ref_mask.sum())
        n_qry = int(qry_mask.sum())
        if n_ref <= 0 or n_qry <= 0:
            raise ValueError(f"Source split failed: reference={n_ref}, query={n_qry}")

        latent_shape = tuple(adata.obsm[cfg.latent_key].shape)
        umap_shape = tuple(adata.obsm[cfg.umap_key].shape)
        if latent_shape[0] != adata.n_obs:
            raise ValueError(f"obsm['{cfg.latent_key}'] row count mismatch: {latent_shape[0]} != {adata.n_obs}")
        if umap_shape[0] != adata.n_obs:
            raise ValueError(f"obsm['{cfg.umap_key}'] row count mismatch: {umap_shape[0]} != {adata.n_obs}")

        present_markers = [gene for gene in cfg.marker_genes if gene in adata.var_names]
        if not present_markers:
            raise ValueError("None of the configured marker genes were found in var_names")
        missing_markers = [gene for gene in cfg.marker_genes if gene not in adata.var_names]
        if missing_markers:
            warnings.warn(f"Missing marker genes ({len(missing_markers)}/{len(cfg.marker_genes)}): {missing_markers}", stacklevel=2)

        print("=" * 80)
        print("Epithelial scHPL preflight")
        print("=" * 80)
        print(f"Input              : {input_path}")
        print(f"Output             : {output_dir}")
        print(f"Version            : {cfg.version}")
        print(f"Reference / Query  : {n_ref:,} / {n_qry:,}")
        print(f"Latent key         : {cfg.latent_key} {latent_shape}")
        print(f"UMAP key           : {cfg.umap_key} {umap_shape}")
        print(f"Marker layer       : {cfg.marker_layer_preferred}")
        print(f"Markers present    : {len(present_markers)}/{len(cfg.marker_genes)}")
        if missing_markers:
            print(f"Markers missing    : {missing_markers}")
        print("=" * 80)
    finally:
        if adata is not None and getattr(adata, "file", None) is not None:
            adata.file.close()


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    cfg = build_config(args)
    validate_config(cfg, allow_overwrite=args.allow_overwrite)
    if args.preflight_only:
        print("[OK] Preflight completed; --preflight-only set, skipping run_pipeline().")
        return
    run_pipeline(cfg)


if __name__ == "__main__":
    main()
