#!/usr/bin/env python3
"""Lightweight scPoli parameter optimization for latest B-cell h5ad.

The previous run was a smoke test. This script keeps the same source h5ad but
uses a larger balanced subset, a fixed stratified labeled/holdout split, and a
small parameter grid focused on better label-transfer accuracy plus less study
separation in the latent space.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import shutil
import time
import traceback
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy import sparse
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    silhouette_score,
)
from sklearn.model_selection import train_test_split

from scarches.models import scPoli


DEFAULT_H5AD = Path(
    "/home/h2048/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508/"
    "bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad"
)
DEFAULT_OUTDIR = Path("/home/h2048/temp/bcell_scpoli_paramopt_20260512")


@dataclass(frozen=True)
class Config:
    name: str
    condition_key: str = "study"
    hidden_layer_sizes: tuple[int, ...] = (128,)
    latent_dim: int = 16
    embedding_dims: int = 16
    dr_rate: float = 0.05
    use_mmd: bool = False
    mmd_on: str = "z"
    beta: float = 1.0
    use_bn: bool = False
    use_ln: bool = True
    epochs: int = 18
    pretraining_epochs: int = 9
    eta: float = 1.5
    lr: float = 1e-3
    alpha_epoch_anneal: int = 20
    weight_decay: float = 0.02
    train_frac: float = 0.9
    batch_size: int = 128
    use_stratified_sampling: bool = False


CONFIGS: list[Config] = [
    Config(
        name="longer_baseline_study_h128_lat16_eta1p5",
        hidden_layer_sizes=(128,),
        latent_dim=16,
        embedding_dims=16,
        epochs=18,
        pretraining_epochs=9,
        eta=1.5,
        weight_decay=0.02,
    ),
    Config(
        name="capacity_study_h256_128_lat24_eta2",
        hidden_layer_sizes=(256, 128),
        latent_dim=24,
        embedding_dims=16,
        dr_rate=0.08,
        epochs=20,
        pretraining_epochs=10,
        eta=2.0,
        weight_decay=0.02,
    ),
    Config(
        name="mmd_study_beta0p2_h128_lat16_eta1p5",
        hidden_layer_sizes=(128,),
        latent_dim=16,
        embedding_dims=16,
        use_mmd=True,
        beta=0.2,
        epochs=18,
        pretraining_epochs=9,
        eta=1.5,
        weight_decay=0.02,
    ),
    Config(
        name="mmd_study_beta0p5_h256_128_lat24_eta2",
        hidden_layer_sizes=(256, 128),
        latent_dim=24,
        embedding_dims=16,
        dr_rate=0.08,
        use_mmd=True,
        beta=0.5,
        epochs=20,
        pretraining_epochs=10,
        eta=2.0,
        weight_decay=0.02,
    ),
    Config(
        name="batch_condition_h128_lat16_eta2",
        condition_key="batch",
        hidden_layer_sizes=(128,),
        latent_dim=16,
        embedding_dims=12,
        epochs=18,
        pretraining_epochs=9,
        eta=2.0,
        weight_decay=0.02,
        use_stratified_sampling=True,
    ),
    Config(
        name="study_mmd_beta0p2_stratified_h128_lat16_eta2",
        condition_key="study",
        hidden_layer_sizes=(128,),
        latent_dim=16,
        embedding_dims=16,
        use_mmd=True,
        beta=0.2,
        epochs=18,
        pretraining_epochs=9,
        eta=2.0,
        weight_decay=0.02,
        use_stratified_sampling=True,
    ),
    Config(
        name="batch_mmd_beta0p2_h128_lat16_eta2",
        condition_key="batch",
        hidden_layer_sizes=(128,),
        latent_dim=16,
        embedding_dims=12,
        use_mmd=True,
        beta=0.2,
        epochs=18,
        pretraining_epochs=9,
        eta=2.0,
        weight_decay=0.02,
        use_stratified_sampling=True,
    ),
]


def set_cpu_safe(seed: int) -> None:
    os.environ.setdefault("PYTHONNOUSERSITE", "1")
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.is_available = lambda: False  # type: ignore[assignment]
    torch.cuda.device_count = lambda: 0  # type: ignore[assignment]
    torch.cuda.get_rng_state_all = lambda: []  # type: ignore[assignment]


def balanced_indices(obs: pd.DataFrame, label_key: str, max_per_label: int, seed: int) -> np.ndarray:
    labels = obs[label_key].astype(str)
    keep = ~(labels.str.lower().isin(["unknown", "none", "nan", "na", ""]))
    labels = labels[keep]
    rng = np.random.default_rng(seed)
    selected: list[np.ndarray] = []
    for _, idx in labels.groupby(labels).groups.items():
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


def stratified_labeled_holdout(obs: pd.DataFrame, label_key: str, test_size: float, seed: int) -> tuple[list[int], list[int]]:
    labels = obs[label_key].astype(str).values
    all_idx = np.arange(obs.shape[0])
    train_idx, test_idx = train_test_split(
        all_idx,
        test_size=test_size,
        random_state=seed,
        stratify=labels,
    )
    return sorted(map(int, train_idx)), sorted(map(int, test_idx))


def safe_silhouette(x: np.ndarray, labels: np.ndarray) -> float | None:
    labels = np.asarray(labels).astype(str)
    vc = pd.Series(labels).value_counts()
    keep_labels = vc[vc >= 2].index
    mask = np.isin(labels, keep_labels)
    labels2 = labels[mask]
    if mask.sum() < 4 or len(np.unique(labels2)) < 2 or len(np.unique(labels2)) >= mask.sum():
        return None
    try:
        return float(silhouette_score(x[mask], labels2, metric="euclidean"))
    except Exception:
        return None


def within_label_condition_asw(x: np.ndarray, label: np.ndarray, condition: np.ndarray) -> float | None:
    label = np.asarray(label).astype(str)
    condition = np.asarray(condition).astype(str)
    vals: list[float] = []
    for lab in sorted(np.unique(label)):
        mask = label == lab
        if mask.sum() < 8:
            continue
        cond_counts = pd.Series(condition[mask]).value_counts()
        keep_conditions = cond_counts[cond_counts >= 2].index
        mask2 = mask & np.isin(condition, keep_conditions)
        if mask2.sum() < 8 or len(np.unique(condition[mask2])) < 2:
            continue
        try:
            vals.append(abs(float(silhouette_score(x[mask2], condition[mask2], metric="euclidean"))))
        except Exception:
            continue
    if not vals:
        return None
    return float(np.mean(vals))


def composite_score(metrics: dict[str, Any]) -> float:
    macro_f1 = metrics.get("holdout_macro_f1") or 0.0
    bal_acc = metrics.get("holdout_balanced_accuracy") or 0.0
    label_asw = metrics.get("label_silhouette")
    label_score = ((label_asw + 1.0) / 2.0) if label_asw is not None else 0.5
    batch_asw = metrics.get("within_label_condition_asw_abs")
    batch_score = 1.0 - min(abs(batch_asw), 1.0) if batch_asw is not None else 0.5
    return float(0.45 * macro_f1 + 0.25 * bal_acc + 0.20 * label_score + 0.10 * batch_score)


def run_config(
    adata_base: ad.AnnData,
    cfg: Config,
    label_key: str,
    counts_layer: str,
    labeled_idx: list[int],
    holdout_idx: list[int],
    outdir: Path,
    seed: int,
    save_model: bool = False,
) -> dict[str, Any]:
    run_dir = outdir / cfg.name
    run_dir.mkdir(parents=True, exist_ok=True)
    adata = adata_base.copy()
    ensure_counts_x(adata, counts_layer)
    adata.obs[label_key] = adata.obs[label_key].astype(str).astype("category")
    adata.obs[cfg.condition_key] = adata.obs[cfg.condition_key].astype(str).astype("category")
    adata.obs["scpoli_split"] = "unlabeled_holdout"
    adata.obs.iloc[labeled_idx, adata.obs.columns.get_loc("scpoli_split")] = "labeled_train"
    adata.obs["scpoli_split"] = adata.obs["scpoli_split"].astype("category")

    started = time.time()
    print(f"\n## RUN {cfg.name}", flush=True)
    print(json.dumps(asdict(cfg), ensure_ascii=False), flush=True)

    model = scPoli(
        adata=adata,
        condition_keys=cfg.condition_key,
        cell_type_keys=label_key,
        labeled_indices=labeled_idx,
        hidden_layer_sizes=list(cfg.hidden_layer_sizes),
        latent_dim=cfg.latent_dim,
        embedding_dims=cfg.embedding_dims,
        recon_loss="nb",
        use_mmd=cfg.use_mmd,
        mmd_on=cfg.mmd_on,
        beta=cfg.beta,
        use_bn=cfg.use_bn,
        use_ln=cfg.use_ln,
        dr_rate=cfg.dr_rate,
    )
    model.train(
        n_epochs=cfg.epochs,
        pretraining_epochs=cfg.pretraining_epochs,
        batch_size=cfg.batch_size,
        train_frac=cfg.train_frac,
        use_early_stopping=False,
        reload_best=False,
        prototype_training=True,
        unlabeled_prototype_training=False,
        monitor=True,
        monitor_only_val=False,
        seed=seed,
        eta=cfg.eta,
        lr=cfg.lr,
        alpha_epoch_anneal=cfg.alpha_epoch_anneal,
        weight_decay=cfg.weight_decay,
        use_stratified_sampling=cfg.use_stratified_sampling,
    )

    latent = model.get_latent(adata, mean=True)
    pred = model.classify(adata, get_prob=True, scale_uncertainties=True)[label_key]
    y_true = adata.obs[label_key].astype(str).values
    y_pred = np.asarray(pred["preds"]).astype(str)
    holdout = np.asarray(holdout_idx, dtype=int)
    train = np.asarray(labeled_idx, dtype=int)

    metrics: dict[str, Any] = {
        "config": asdict(cfg),
        "name": cfg.name,
        "status": "ok",
        "runtime_sec": round(time.time() - started, 3),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "train_n": int(len(train)),
        "holdout_n": int(len(holdout)),
        "train_accuracy": float(accuracy_score(y_true[train], y_pred[train])),
        "holdout_accuracy": float(accuracy_score(y_true[holdout], y_pred[holdout])),
        "holdout_balanced_accuracy": float(balanced_accuracy_score(y_true[holdout], y_pred[holdout])),
        "holdout_macro_f1": float(f1_score(y_true[holdout], y_pred[holdout], average="macro", zero_division=0)),
        "holdout_weighted_f1": float(f1_score(y_true[holdout], y_pred[holdout], average="weighted", zero_division=0)),
        "label_silhouette": safe_silhouette(latent, y_true),
        "condition_silhouette": safe_silhouette(latent, adata.obs[cfg.condition_key].astype(str).values),
        "within_label_condition_asw_abs": within_label_condition_asw(
            latent,
            y_true,
            adata.obs[cfg.condition_key].astype(str).values,
        ),
        "condition_ari_vs_label": float(adjusted_rand_score(y_true, adata.obs[cfg.condition_key].astype(str).values)),
        "uncertainty_holdout_median": float(np.median(np.asarray(pred["uncert"])[holdout])),
        "uncertainty_all_median": float(np.median(np.asarray(pred["uncert"]))),
    }
    metrics["composite_score"] = composite_score(metrics)
    logs = getattr(model.trainer, "logs", {})
    metrics["final_epoch_cvae_loss"] = float(logs["epoch_cvae_loss"][-1]) if "epoch_cvae_loss" in logs else None
    metrics["final_val_cvae_loss"] = float(logs["val_cvae_loss"][-1]) if "val_cvae_loss" in logs else None
    metrics["trainer_logs"] = {k: [float(x) for x in v] for k, v in logs.items()}

    adata.obs[f"{cfg.name}_pred"] = pd.Categorical(y_pred)
    adata.obs[f"{cfg.name}_uncert"] = np.asarray(pred["uncert"], dtype=float)
    adata.obsm[f"X_scPoli_{cfg.name}"] = latent

    labels = sorted(pd.unique(pd.Series(y_true)))
    conf = pd.DataFrame(
        confusion_matrix(y_true[holdout], y_pred[holdout], labels=labels),
        index=labels,
        columns=labels,
    )
    conf.to_csv(run_dir / "holdout_confusion.tsv", sep="\t")
    (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8")

    # Keep a compact h5ad for every run with obs + latent for downstream visual checks.
    slim = ad.AnnData(
        X=sparse.csr_matrix((adata.n_obs, 0)),
        obs=adata.obs[[label_key, cfg.condition_key, "scpoli_split", f"{cfg.name}_pred", f"{cfg.name}_uncert"]].copy(),
        obsm={f"X_scPoli_{cfg.name}": latent},
    )
    slim.write_h5ad(run_dir / "latent_predictions_only.h5ad")
    if save_model:
        model.save(str(run_dir / "model"), overwrite=True, save_anndata=False)

    print(
        f"DONE {cfg.name}: holdout_acc={metrics['holdout_accuracy']:.3f}, "
        f"macroF1={metrics['holdout_macro_f1']:.3f}, composite={metrics['composite_score']:.3f}",
        flush=True,
    )
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description="Optimize scPoli parameters on latest B-cell h5ad")
    parser.add_argument("--input", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--label-key", default="cell_type_expert")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--max-per-label", type=int, default=220)
    parser.add_argument("--holdout-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=20260512)
    parser.add_argument("--only", nargs="*", default=None, help="Optional config names to run")
    args = parser.parse_args()

    set_cpu_safe(args.seed)
    if args.outdir.exists():
        shutil.rmtree(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    full = ad.read_h5ad(args.input)
    selected = balanced_indices(full.obs, args.label_key, args.max_per_label, args.seed)
    adata = full[selected].copy()
    del full

    # Drop singleton labels for prototype covariance safety.
    counts = adata.obs[args.label_key].astype(str).value_counts()
    valid_labels = counts[counts >= 5].index.tolist()
    adata = adata[adata.obs[args.label_key].astype(str).isin(valid_labels)].copy()
    for key in [args.label_key, "study", "batch", "dataset", "tissue"]:
        if key in adata.obs:
            adata.obs[key] = adata.obs[key].astype(str).astype("category")

    labeled_idx, holdout_idx = stratified_labeled_holdout(adata.obs, args.label_key, args.holdout_frac, args.seed)
    manifest = {
        "input_h5ad": str(args.input),
        "outdir": str(args.outdir),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "label_key": args.label_key,
        "counts_layer": args.counts_layer,
        "max_per_label": args.max_per_label,
        "holdout_frac": args.holdout_frac,
        "seed": args.seed,
        "label_counts": adata.obs[args.label_key].astype(str).value_counts().to_dict(),
        "study_counts": adata.obs["study"].astype(str).value_counts().to_dict() if "study" in adata.obs else {},
        "batch_counts": adata.obs["batch"].astype(str).value_counts().to_dict() if "batch" in adata.obs else {},
        "train_n": len(labeled_idx),
        "holdout_n": len(holdout_idx),
        "configs": [asdict(c) for c in CONFIGS],
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    pd.DataFrame({"cell": adata.obs_names, "split": "holdout"}).assign(
        split=lambda x: np.where(np.isin(np.arange(adata.n_obs), holdout_idx), "holdout", "labeled_train"),
        label=adata.obs[args.label_key].astype(str).values,
        study=adata.obs["study"].astype(str).values if "study" in adata.obs else "",
        batch=adata.obs["batch"].astype(str).values if "batch" in adata.obs else "",
    ).to_csv(args.outdir / "fixed_split.tsv", sep="\t", index=False)

    requested = set(args.only) if args.only else None
    results: list[dict[str, Any]] = []
    for cfg in CONFIGS:
        if requested is not None and cfg.name not in requested:
            continue
        try:
            results.append(
                run_config(
                    adata,
                    cfg,
                    args.label_key,
                    args.counts_layer,
                    labeled_idx,
                    holdout_idx,
                    args.outdir,
                    args.seed,
                    save_model=True,
                )
            )
        except Exception as exc:
            err = {
                "name": cfg.name,
                "config": asdict(cfg),
                "status": "error",
                "error": repr(exc),
                "traceback": traceback.format_exc(),
            }
            results.append(err)
            run_dir = args.outdir / cfg.name
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "metrics.json").write_text(json.dumps(err, indent=2, ensure_ascii=False), encoding="utf-8")
            print(f"ERROR {cfg.name}: {exc!r}", flush=True)

    result_df = pd.DataFrame([{k: v for k, v in r.items() if k not in {"trainer_logs", "config", "traceback"}} for r in results])
    if "composite_score" in result_df:
        result_df = result_df.sort_values(["status", "composite_score"], ascending=[True, False])
    result_df.to_csv(args.outdir / "scpoli_paramopt_results.tsv", sep="\t", index=False)

    ok = [r for r in results if r.get("status") == "ok"]
    if ok:
        best = max(ok, key=lambda r: r.get("composite_score", -np.inf))
        (args.outdir / "best_config.json").write_text(json.dumps(best, indent=2, ensure_ascii=False), encoding="utf-8")
        # Save the best full model by reusing the already trained model output is not possible here,
        # so the best run keeps latent/predictions; rerun with --only best_name and save_model=True if needed.
        print("\nBEST", best["name"], json.dumps({
            "holdout_accuracy": best["holdout_accuracy"],
            "holdout_macro_f1": best["holdout_macro_f1"],
            "holdout_balanced_accuracy": best["holdout_balanced_accuracy"],
            "label_silhouette": best["label_silhouette"],
            "within_label_condition_asw_abs": best["within_label_condition_asw_abs"],
            "composite_score": best["composite_score"],
        }, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
