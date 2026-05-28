#!/usr/bin/env python3
"""Toy + real-subset smoke tests for the normal-airway ML MVP."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import anndata as ad

from normal_airway_ml_audit_20260527 import run_audit
from normal_airway_ml_common_20260527 import (
    DEFAULT_OUTPUT_ROOT,
    deep_get,
    ensure_dir,
    load_config,
    load_real_smoke_subset,
    make_toy_airway_adata,
    timestamp_slug,
    write_json,
)
from normal_airway_ml_feature_export_20260527 import run_feature_export
from normal_airway_ml_train_sample_level_20260527 import run_training


def build_smoke_input(smoke_mode: str, cfg: dict[str, Any]) -> tuple[ad.AnnData, dict[str, Any]]:
    seed = int(deep_get(cfg, "run", "random_seed", default=20260527))
    if smoke_mode == "toy":
        adata = make_toy_airway_adata(cfg, seed=seed)
        return adata, {"source_mode": "toy", "shape": [int(adata.n_obs), int(adata.n_vars)]}
    if smoke_mode == "real":
        input_h5ad = Path(deep_get(cfg, "smoke", "real_input_h5ad", default=deep_get(cfg, "data_contract", "input_h5ad", default="")))
        if not input_h5ad.exists():
            raise FileNotFoundError(f"Real smoke input h5ad not found: {input_h5ad}")
        adata, summary = load_real_smoke_subset(input_h5ad, cfg=cfg, seed=seed)
        return adata, {"source_mode": "real", **summary}
    raise ValueError(f"Unsupported smoke_mode: {smoke_mode}")


def run_smoke(smoke_mode: str, cfg: dict[str, Any], output_dir: Path | str) -> dict[str, Any]:
    output_dir = ensure_dir(output_dir)
    artifact_dir = ensure_dir(Path(output_dir) / "artifacts")
    audit_dir = ensure_dir(Path(output_dir) / "audit")
    feature_dir = ensure_dir(Path(output_dir) / "feature_export")
    train_dir = ensure_dir(Path(output_dir) / "train")

    try:
        adata, input_summary = build_smoke_input(smoke_mode, cfg)
    except ValueError as exc:
        summary = {
            "smoke_mode": smoke_mode,
            "status": "blocked",
            "reason": str(exc),
            "configured_real_input_h5ad": deep_get(cfg, "smoke", "real_input_h5ad", default=None),
        }
        write_json(summary, Path(output_dir) / "smoke_summary.json")
        return summary

    input_h5ad = artifact_dir / f"normal_airway_ml_smoke_input_{smoke_mode}.h5ad"
    adata.write_h5ad(input_h5ad, compression="gzip")

    audit_summary = run_audit(input_h5ad=input_h5ad, cfg=cfg, output_dir=audit_dir)
    feature_manifest = run_feature_export(adata=adata, cfg=cfg, output_dir=feature_dir)
    model_manifest = run_training(feature_dir=Path(feature_dir) / "features", cfg=cfg, output_dir=train_dir)

    summary = {
        "smoke_mode": smoke_mode,
        "status": "ok",
        "input_summary": input_summary,
        "input_h5ad": str(input_h5ad),
        "audit_summary_json": str(Path(audit_dir) / "audit_summary.json"),
        "feature_manifest_json": str(Path(feature_dir) / "features" / "feature_manifest.json"),
        "model_manifest_json": str(Path(train_dir) / "model" / "model_manifest.json"),
        "audit_analysis_cells": audit_summary.get("n_analysis_cells"),
        "feature_samples": feature_manifest.get("n_samples"),
        "trained_models": model_manifest.get("model_names"),
    }
    write_json(summary, Path(output_dir) / "smoke_summary.json")
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run normal-airway ML smoke tests")
    parser.add_argument("--smoke-mode", choices=["toy", "real"], default="toy")
    parser.add_argument("--config-yaml", default=None, help="Config YAML path")
    parser.add_argument("--output-dir", default=None, help="Explicit output dir")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    cfg = load_config(args.config_yaml)
    output_dir = Path(args.output_dir) if args.output_dir else ensure_dir(Path(deep_get(cfg, "run", "output_root", default=str(DEFAULT_OUTPUT_ROOT))) / f"normal_airway_ml_smoke_{args.smoke_mode}_{timestamp_slug()}")
    summary = run_smoke(smoke_mode=args.smoke_mode, cfg=cfg, output_dir=output_dir)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
