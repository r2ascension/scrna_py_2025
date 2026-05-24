#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import scanpy as sc

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lineage_merge_schpl_core_20260402 import LineageMergeSchPLConfig, run_pipeline


def _split_marker_genes(values: list[str] | None) -> list[str]:
    if not values:
        return []
    genes: list[str] = []
    for value in values:
        if "," in value:
            genes.extend(item.strip() for item in value.split(",") if item.strip())
        else:
            item = value.strip()
            if item:
                genes.append(item)
    # keep order, remove duplicates
    return list(dict.fromkeys(genes))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generic lineage merge + scHPL CLI wrapper"
    )
    parser.add_argument("--lineage-name", required=True, help="Display name, e.g. 'T/NK'")
    parser.add_argument("--lineage-slug", required=True, help="Short slug, e.g. 'tnk'")
    parser.add_argument("--version", required=True, help="Run version tag")
    parser.add_argument("--input-h5ad", required=True, help="Merged reference+query h5ad")
    parser.add_argument("--output-dir", required=True, help="Output directory for scHPL run")
    parser.add_argument("--reference-label-key", required=True)
    parser.add_argument("--query-pred-key", required=True)
    parser.add_argument("--query-final-key", required=True)
    parser.add_argument("--query-conf-key", required=True)
    parser.add_argument("--latent-key", required=True)
    parser.add_argument("--umap-key", required=True)
    parser.add_argument("--marker-layer-preferred", default="log1p")
    parser.add_argument("--marker-genes", nargs="*", default=[])
    parser.add_argument("--source-key", default="data_source")
    parser.add_argument("--reference-source-value", default="reference")
    parser.add_argument("--query-source-value", default="query")
    parser.add_argument("--batch-key", default="sample")
    parser.add_argument("--tissue-key", default="tissue")
    parser.add_argument("--unlabeled", default="Unknown")
    parser.add_argument("--allow-overwrite", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    return parser


def build_config(args: argparse.Namespace) -> LineageMergeSchPLConfig:
    return LineageMergeSchPLConfig(
        lineage_name=args.lineage_name,
        lineage_slug=args.lineage_slug,
        version=args.version,
        input_merged_h5ad=str(Path(args.input_h5ad).expanduser().resolve()),
        output_dir=str(Path(args.output_dir).expanduser().resolve()),
        reference_label_key=args.reference_label_key,
        query_pred_key=args.query_pred_key,
        query_final_key=args.query_final_key,
        query_conf_key=args.query_conf_key,
        marker_genes=_split_marker_genes(args.marker_genes),
        source_key=args.source_key,
        reference_source_value=args.reference_source_value,
        query_source_value=args.query_source_value,
        batch_key=args.batch_key,
        tissue_key=args.tissue_key,
        latent_key=args.latent_key,
        umap_key=args.umap_key,
        marker_layer_preferred=args.marker_layer_preferred,
        unlabeled=args.unlabeled,
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

    if output_dir.exists():
        existing_items = [p.name for p in output_dir.iterdir()]
        if existing_items and not allow_overwrite:
            raise FileExistsError(
                f"Output directory already exists and is not empty: {output_dir}\n"
                "Pass --allow-overwrite only when you intentionally want to reuse it."
            )
        if existing_items and allow_overwrite:
            warnings.warn(
                f"Output directory is non-empty and will be reused: {output_dir}",
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
        missing_obs = [k for k in required_obs if k not in adata.obs.columns]
        if missing_obs:
            raise KeyError(f"Missing obs keys in input h5ad: {missing_obs}")

        required_obsm = [cfg.latent_key, cfg.umap_key]
        missing_obsm = [k for k in required_obsm if k not in adata.obsm.keys()]
        if missing_obsm:
            raise KeyError(f"Missing obsm keys in input h5ad: {missing_obsm}")

        n_obs = int(adata.n_obs)
        latent_shape = tuple(adata.obsm[cfg.latent_key].shape)
        umap_shape = tuple(adata.obsm[cfg.umap_key].shape)
        if latent_shape[0] != n_obs:
            raise ValueError(
                f"obsm['{cfg.latent_key}'] row count {latent_shape[0]} does not match n_obs={n_obs}"
            )
        if umap_shape[0] != n_obs:
            raise ValueError(
                f"obsm['{cfg.umap_key}'] row count {umap_shape[0]} does not match n_obs={n_obs}"
            )

        src = adata.obs[cfg.source_key].astype(str)
        ref_mask = src.eq(cfg.reference_source_value)
        qry_mask = src.eq(cfg.query_source_value)
        n_ref = int(ref_mask.sum())
        n_qry = int(qry_mask.sum())
        if n_ref <= 0 or n_qry <= 0:
            raise ValueError(
                f"Expected both reference and query cells in '{cfg.source_key}', got ref={n_ref}, qry={n_qry}"
            )

        ref_labels = adata.obs.loc[ref_mask, cfg.reference_label_key]
        if int(ref_labels.isna().sum()) == n_ref:
            raise ValueError(
                f"All reference labels are NA for '{cfg.reference_label_key}'"
            )

        for key in [cfg.query_pred_key, cfg.query_final_key, cfg.query_conf_key]:
            qry_vals = adata.obs.loc[qry_mask, key]
            if int(qry_vals.isna().sum()) == n_qry:
                raise ValueError(f"All query values are NA for key '{key}'")

        marker_hits = [g for g in cfg.marker_genes if g in adata.var_names]
        if cfg.marker_genes and not marker_hits:
            warnings.warn(
                "None of the configured marker genes were found in var_names; "
                "marker panel will be disabled in downstream output.",
                stacklevel=2,
            )

        print("=" * 80)
        print("Generic lineage schpl preflight")
        print("=" * 80)
        print(f"Input              : {input_path}")
        print(f"Output             : {output_dir}")
        print(f"Lineage            : {cfg.lineage_name} ({cfg.lineage_slug})")
        print(f"Version            : {cfg.version}")
        print(f"Reference / Query  : {n_ref:,} / {n_qry:,}")
        print(f"Reference label key: {cfg.reference_label_key}")
        print(f"Query pred/final   : {cfg.query_pred_key} / {cfg.query_final_key}")
        print(f"Confidence key     : {cfg.query_conf_key}")
        print(f"Latent / UMAP      : {cfg.latent_key} {latent_shape} / {cfg.umap_key} {umap_shape}")
        print(f"Marker genes found : {len(marker_hits)}/{len(cfg.marker_genes)}")
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
