#!/usr/bin/env python3
"""Export sample-level feature families for normal-airway site classification."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from normal_airway_ml_common_20260527 import (
    DEFAULT_OUTPUT_ROOT,
    aggregate_composition_features,
    aggregate_latent_summary_features,
    aggregate_pseudobulk_pca_features,
    build_sample_metadata,
    deep_get,
    ensure_counts_layer,
    ensure_dir,
    join_feature_blocks,
    load_config,
    make_toy_airway_adata,
    prepare_obs_contract,
    timestamp_slug,
    write_json,
)


def run_feature_export(
    adata: ad.AnnData,
    cfg: dict[str, Any],
    output_dir: Path | str,
    cell_type_key_override: str | None = None,
) -> dict[str, Any]:
    output_dir = ensure_dir(output_dir)
    feature_dir = ensure_dir(Path(output_dir) / "features")

    counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
    counts_info = ensure_counts_layer(adata, counts_layer)
    prepared_obs, contract = prepare_obs_contract(
        adata.obs.copy(),
        list(adata.obsm.keys()),
        cfg,
        cell_type_key_override=cell_type_key_override or deep_get(cfg, "feature_export", "cell_type_key", default=None),
    )
    adata.obs = prepared_obs

    sample_meta = build_sample_metadata(adata.obs, contract["analysis_mask"])
    epsilon = float(deep_get(cfg, "feature_export", "epsilon", default=1e-6))
    composition = aggregate_composition_features(adata.obs, contract["analysis_mask"], epsilon=epsilon)
    latent = aggregate_latent_summary_features(
        adata,
        adata.obs,
        contract["analysis_mask"],
        contract.get("latent_key"),
        max_dims=int(deep_get(cfg, "feature_export", "latent_max_dims", default=16)),
    )
    pseudobulk = aggregate_pseudobulk_pca_features(
        adata,
        adata.obs,
        contract["analysis_mask"],
        counts_layer=counts_layer,
        min_cells_per_group=int(deep_get(cfg, "feature_export", "min_cells_per_sample_celltype", default=10)),
        min_samples_per_celltype=int(deep_get(cfg, "feature_export", "min_samples_per_celltype", default=3)),
        top_genes=int(deep_get(cfg, "feature_export", "pseudobulk_top_genes", default=200)),
        n_pcs=int(deep_get(cfg, "feature_export", "pseudobulk_n_pcs", default=5)),
    )

    feature_table = join_feature_blocks(
        sample_meta,
        [composition["wide"], latent.get("wide", pd.DataFrame()), pseudobulk.get("wide", pd.DataFrame())],
    )
    feature_table.index.name = "sample"

    long_blocks = []
    if not composition["long"].empty:
        comp_long = composition["long"].copy()
        comp_long["feature_family"] = "composition"
        long_blocks.append(comp_long)
    sample_features_long = pd.concat(long_blocks, ignore_index=True) if long_blocks else pd.DataFrame()

    feature_table.to_csv(feature_dir / "sample_features_wide.tsv.gz", sep="\t", compression="gzip", index_label="sample_id")
    sample_meta.to_csv(feature_dir / "sample_metadata.tsv", sep="\t", index=False)
    composition["counts"].to_csv(feature_dir / "composition_counts.tsv", sep="\t")
    composition["fractions"].to_csv(feature_dir / "composition_fractions.tsv", sep="\t")
    if not sample_features_long.empty:
        sample_features_long.to_csv(feature_dir / "sample_features_long.tsv.gz", sep="\t", index=False, compression="gzip")

    manifest = {
        "output_dir": str(feature_dir),
        "n_samples": int(feature_table.shape[0]),
        "n_feature_columns": int(feature_table.shape[1] - sample_meta.shape[1]),
        "counts_info": counts_info,
        "resolved_columns": {
            k: contract.get(k)
            for k in ["sample_col", "sample_unit_col", "dataset_col", "study_col", "batch_col", "tissue_col", "condition_col", "cell_type_col", "latent_key"]
        },
        "analysis_n_cells": int(contract["n_analysis_cells"]),
        "feature_family_status": {
            "composition": "available",
            "latent_summary": latent.get("status", "skipped"),
            "pseudobulk_pca": pseudobulk.get("status", "skipped"),
            "core_discovery_input": "available",
            "cnmf_usage": "pending",
            "communication": "pending",
        },
        "site_counts": sample_meta["site_label"].value_counts().to_dict() if not sample_meta.empty else {},
    }
    write_json(manifest, feature_dir / "feature_manifest.json")
    return manifest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export normal-airway ML feature families")
    parser.add_argument("--input-h5ad", default=None, help="Input h5ad path")
    parser.add_argument("--config-yaml", default=None, help="Config YAML path")
    parser.add_argument("--output-dir", default=None, help="Explicit output dir")
    parser.add_argument("--toy", action="store_true", help="Use internally generated toy AnnData")
    parser.add_argument("--cell-type-key", default=None, help="Override cell type key")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    cfg = load_config(args.config_yaml)
    output_dir = Path(args.output_dir) if args.output_dir else ensure_dir(Path(deep_get(cfg, "run", "output_root", default=str(DEFAULT_OUTPUT_ROOT))) / f"normal_airway_ml_features_{timestamp_slug()}")
    if args.toy:
        adata = make_toy_airway_adata(cfg)
    else:
        input_h5ad = Path(args.input_h5ad or deep_get(cfg, "data_contract", "input_h5ad", default=""))
        if not input_h5ad.exists():
            raise FileNotFoundError(f"Input h5ad not found: {input_h5ad}")
        adata = ad.read_h5ad(input_h5ad)
    manifest = run_feature_export(adata=adata, cfg=cfg, output_dir=output_dir, cell_type_key_override=args.cell_type_key)
    print(pd.Series(manifest).to_json(force_ascii=False, indent=2))


if __name__ == "__main__":
    main()
