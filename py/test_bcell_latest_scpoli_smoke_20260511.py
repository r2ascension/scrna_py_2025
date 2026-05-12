#!/usr/bin/env python3
"""Smoke test scPoli on the latest B-cell h5ad.

This is intentionally small and conservative: it verifies that the newly copied
scHPL/pertpy/scPoli/scFoundation/scLong environment can train and classify with
scPoli on the latest B-cell reference object without modifying the source h5ad.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy import sparse

from scarches.models import scPoli


DEFAULT_H5AD = Path(
    "/home/h2048/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508/"
    "bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad"
)
DEFAULT_OUTDIR = Path("/home/h2048/temp/bcell_scpoli_smoke_20260511")


def set_cpu_safe(seed: int) -> None:
    """Force CPU path and deterministic-ish seeds for this smoke test."""
    os.environ.setdefault("PYTHONNOUSERSITE", "1")
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    # The host CUDA stack may be unavailable/stale; force scPoli trainer to CPU.
    torch.cuda.is_available = lambda: False  # type: ignore[assignment]
    torch.cuda.device_count = lambda: 0  # type: ignore[assignment]


def balanced_indices(obs: pd.DataFrame, label_key: str, max_per_label: int, seed: int) -> np.ndarray:
    labels = obs[label_key].astype(str)
    keep = ~(labels.str.lower().isin(["unknown", "none", "nan", "na", ""] ))
    labels = labels[keep]
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for label, idx in labels.groupby(labels).groups.items():
        idx_arr = np.array(list(idx))
        if len(idx_arr) > max_per_label:
            idx_arr = rng.choice(idx_arr, size=max_per_label, replace=False)
        selected.append(idx_arr)
    out = np.concatenate(selected)
    rng.shuffle(out)
    return out


def ensure_counts_x(adata: ad.AnnData, counts_layer: str) -> None:
    if counts_layer not in adata.layers:
        raise KeyError(f"counts layer not found: {counts_layer!r}; layers={list(adata.layers.keys())}")
    x = adata.layers[counts_layer].copy()
    if sparse.issparse(x):
        x = x.tocsr().astype(np.float32)
    else:
        x = np.asarray(x, dtype=np.float32)
    adata.X = x


def main() -> None:
    parser = argparse.ArgumentParser(description="Smoke test scPoli on latest B-cell h5ad")
    parser.add_argument("--input", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--label-key", default="cell_type_expert")
    parser.add_argument("--condition-key", default="study")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--max-per-label", type=int, default=120)
    parser.add_argument("--epochs", type=int, default=6)
    parser.add_argument("--pretraining-epochs", type=int, default=4)
    parser.add_argument("--latent-dim", type=int, default=8)
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260511)
    args = parser.parse_args()

    set_cpu_safe(args.seed)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "model").mkdir(exist_ok=True)

    t0 = time.time()
    adata_full = ad.read_h5ad(args.input)
    if args.label_key not in adata_full.obs:
        raise KeyError(f"label key missing: {args.label_key}")
    if args.condition_key not in adata_full.obs:
        raise KeyError(f"condition key missing: {args.condition_key}")

    selected = balanced_indices(adata_full.obs, args.label_key, args.max_per_label, args.seed)
    adata = adata_full[selected].copy()
    del adata_full

    ensure_counts_x(adata, args.counts_layer)
    adata.obs[args.label_key] = adata.obs[args.label_key].astype(str).astype("category")
    adata.obs[args.condition_key] = adata.obs[args.condition_key].astype(str).astype("category")

    # Keep only labels with at least two cells after sampling, because prototype
    # covariance is undefined for singleton classes.
    counts = adata.obs[args.label_key].astype(str).value_counts()
    valid_labels = counts[counts >= 2].index.tolist()
    adata = adata[adata.obs[args.label_key].astype(str).isin(valid_labels)].copy()
    adata.obs[args.label_key] = adata.obs[args.label_key].astype(str).astype("category")
    adata.obs[args.condition_key] = adata.obs[args.condition_key].astype(str).astype("category")

    labeled_indices = list(range(adata.n_obs))
    print("scPoli B-cell smoke test")
    print("input", args.input)
    print("outdir", args.outdir)
    print("shape", adata.shape)
    print("label counts")
    print(adata.obs[args.label_key].astype(str).value_counts().to_string())
    print("condition counts")
    print(adata.obs[args.condition_key].astype(str).value_counts().to_string())

    model = scPoli(
        adata=adata,
        condition_keys=args.condition_key,
        cell_type_keys=args.label_key,
        labeled_indices=labeled_indices,
        hidden_layer_sizes=[args.hidden],
        latent_dim=args.latent_dim,
        embedding_dims=8,
        recon_loss="nb",
        use_mmd=False,
        use_bn=False,
        use_ln=True,
        dr_rate=0.05,
    )
    model.train(
        n_epochs=args.epochs,
        pretraining_epochs=args.pretraining_epochs,
        batch_size=args.batch_size,
        train_frac=0.9,
        use_early_stopping=False,
        reload_best=False,
        unlabeled_prototype_training=False,
        monitor=True,
        monitor_only_val=False,
        seed=args.seed,
    )

    latent = model.get_latent(adata, mean=True)
    adata.obsm["X_scPoli_smoke"] = latent
    pred = model.classify(adata, get_prob=True, scale_uncertainties=True)[args.label_key]
    adata.obs["scpoli_pred"] = pd.Categorical(pred["preds"])
    adata.obs["scpoli_uncert"] = pred["uncert"].astype(float)
    adata.obs["scpoli_correct"] = (
        adata.obs["scpoli_pred"].astype(str).values == adata.obs[args.label_key].astype(str).values
    )

    confusion = pd.crosstab(
        adata.obs[args.label_key].astype(str),
        adata.obs["scpoli_pred"].astype(str),
        rownames=[args.label_key],
        colnames=["scpoli_pred"],
    )
    per_label = pd.DataFrame(
        {
            "n": adata.obs[args.label_key].astype(str).value_counts().sort_index(),
            "accuracy": adata.obs.groupby(args.label_key, observed=True)["scpoli_correct"].mean().sort_index(),
        }
    )

    adata_path = args.outdir / "bcell_latest_scpoli_smoke.h5ad"
    confusion_path = args.outdir / "bcell_latest_scpoli_confusion.tsv"
    per_label_path = args.outdir / "bcell_latest_scpoli_per_label.tsv"
    summary_path = args.outdir / "bcell_latest_scpoli_summary.json"

    adata.write_h5ad(adata_path)
    confusion.to_csv(confusion_path, sep="\t")
    per_label.to_csv(per_label_path, sep="\t")
    model.save(str(args.outdir / "model"), overwrite=True, save_anndata=False)

    logs = {k: [float(x) for x in v] for k, v in getattr(model.trainer, "logs", {}).items()}
    summary = {
        "input_h5ad": str(args.input),
        "output_dir": str(args.outdir),
        "adata_output": str(adata_path),
        "model_dir": str(args.outdir / "model"),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "label_key": args.label_key,
        "condition_key": args.condition_key,
        "counts_layer": args.counts_layer,
        "epochs": args.epochs,
        "pretraining_epochs": args.pretraining_epochs,
        "latent_dim": args.latent_dim,
        "hidden": args.hidden,
        "batch_size": args.batch_size,
        "seed": args.seed,
        "label_counts": adata.obs[args.label_key].astype(str).value_counts().to_dict(),
        "condition_counts": adata.obs[args.condition_key].astype(str).value_counts().to_dict(),
        "overall_self_classification_accuracy": float(adata.obs["scpoli_correct"].mean()),
        "uncertainty_min": float(adata.obs["scpoli_uncert"].min()),
        "uncertainty_median": float(adata.obs["scpoli_uncert"].median()),
        "uncertainty_max": float(adata.obs["scpoli_uncert"].max()),
        "latent_mean_abs": float(np.mean(np.abs(latent))),
        "runtime_sec": round(time.time() - t0, 3),
        "trainer_logs": logs,
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nDONE")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
