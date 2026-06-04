#!/usr/bin/env python3
"""Contrast- and stability-oriented core-difference discovery for normal-airway site comparisons.

This module complements the original sample-level classification MVP by adding:
1. sample-unit-aware pseudobulk aggregation within each cell type,
2. pairwise and one-vs-rest tissue contrasts,
3. repeated resampling stability selection across multiple model families, and
4. shared / site-specific core-gene summaries.
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
from scipy import stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

PY_ROOT = Path(__file__).resolve().parents[1]
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from normal_airway_ml_common_20260527 import (  # noqa: E402
    DEFAULT_OUTPUT_ROOT,
    _matrix_take_rows,
    deep_get,
    ensure_counts_layer,
    ensure_dir,
    load_config,
    prepare_obs_contract,
    timestamp_slug,
    write_json,
)

try:  # pragma: no cover - optional
    from xgboost import XGBClassifier
except Exception:  # pragma: no cover
    XGBClassifier = None


CORE_METHOD_ORDER = [
    "lasso_logistic",
    "elastic_net_logistic",
    "ridge_logistic",
    "linear_svm",
    "random_forest",
    "lda",
    "naive_bayes",
    "xgboost",
]


def slugify_token(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(text).strip()).strip("_") or "unknown"


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.ones_like(p, dtype=float)
    finite = np.isfinite(p)
    if not np.any(finite):
        return out
    idx = np.flatnonzero(finite)
    ranked_order = idx[np.argsort(p[idx])]
    ranked_p = p[ranked_order]
    n = ranked_p.size
    adj = ranked_p * n / np.arange(1, n + 1, dtype=float)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out[ranked_order] = np.clip(adj, 0.0, 1.0)
    return out


def normalize_log_cpm(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    library = counts.sum(axis=1, keepdims=True)
    library[library <= 0] = 1.0
    return np.log2((counts / library) * 1e6 + 1.0)


def build_site_contrasts(sample_meta: pd.DataFrame, cfg: dict[str, Any]) -> list[dict[str, Any]]:
    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    site_order = [str(x) for x in deep_get(core_cfg, "site_order", default=["nasal", "sinus", "bronchus", "lung_parenchyma"])]
    available_sites = [site for site in site_order if site in set(sample_meta["site_label"].astype(str))]
    min_samples_per_group = int(deep_get(core_cfg, "min_samples_per_group", default=2))
    min_total_samples = int(deep_get(core_cfg, "min_total_samples", default=4))

    contrasts: list[dict[str, Any]] = []
    if bool(deep_get(core_cfg, "pairwise", default=True)):
        for i, site_a in enumerate(available_sites):
            for site_b in available_sites[i + 1 :]:
                mask = sample_meta["site_label"].astype(str).isin([site_a, site_b])
                sub = sample_meta.loc[mask]
                n_a = int((sub["site_label"].astype(str) == site_a).sum())
                n_b = int((sub["site_label"].astype(str) == site_b).sum())
                if min(n_a, n_b) < min_samples_per_group or int(sub.shape[0]) < min_total_samples:
                    continue
                contrasts.append(
                    {
                        "contrast_id": f"pairwise__{site_a}__vs__{site_b}",
                        "contrast_type": "pairwise",
                        "positive_label": site_b,
                        "negative_label": site_a,
                        "included_sites": [site_a, site_b],
                    }
                )
    if bool(deep_get(core_cfg, "one_vs_rest", default=True)):
        for site in available_sites:
            n_pos = int((sample_meta["site_label"].astype(str) == site).sum())
            n_neg = int(sample_meta.shape[0] - n_pos)
            if min(n_pos, n_neg) < min_samples_per_group or int(sample_meta.shape[0]) < min_total_samples:
                continue
            contrasts.append(
                {
                    "contrast_id": f"one_vs_rest__{site}",
                    "contrast_type": "one_vs_rest",
                    "positive_label": site,
                    "negative_label": "rest",
                    "included_sites": [site, "rest"],
                }
            )
    return contrasts


def aggregate_pseudobulk_by_celltype(adata: ad.AnnData, prepared_obs: pd.DataFrame, contract: dict[str, Any], cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    counts_layer = str(deep_get(core_cfg, "counts_layer", default=deep_get(cfg, "data_contract", "counts_layer", default="counts")))
    ensure_counts_layer(adata, counts_layer)
    min_cells = int(deep_get(core_cfg, "min_cells_per_sample_celltype", default=deep_get(cfg, "feature_export", "min_cells_per_sample_celltype", default=10)))
    sample_group_col = str(contract.get("sample_unit_col", "sample_unit_resolved"))

    counts_matrix = adata.layers[counts_layer]
    frame = prepared_obs.loc[
        contract["analysis_mask"],
        [
            sample_group_col,
            "sample_resolved",
            "dataset_resolved",
            "study_resolved",
            "batch_resolved",
            "site_label",
            "site_axis",
            "condition_resolved",
            "cell_type_resolved",
        ],
    ].copy()
    frame["row_idx"] = np.flatnonzero(np.asarray(contract["analysis_mask"], dtype=bool))

    pb_by_celltype: dict[str, dict[str, Any]] = {}
    group_sizes = (
        frame.groupby(["cell_type_resolved", sample_group_col], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    valid_groups = group_sizes[group_sizes["n_cells"] >= min_cells]
    for celltype, sub_groups in valid_groups.groupby("cell_type_resolved", observed=True):
        counts_rows: list[np.ndarray] = []
        meta_rows: list[dict[str, Any]] = []
        valid_sample_units = sub_groups[sample_group_col].astype(str).tolist()
        for sample_unit in valid_sample_units:
            cell_rows = frame.loc[
                (frame["cell_type_resolved"].astype(str) == str(celltype))
                & (frame[sample_group_col].astype(str) == str(sample_unit))
            ].copy()
            if int(cell_rows.shape[0]) < min_cells:
                continue
            idx = cell_rows["row_idx"].to_numpy(dtype=int)
            counts_rows.append(np.asarray(_matrix_take_rows(counts_matrix, idx).sum(axis=0), dtype=np.float32))
            meta_rows.append(
                {
                    "sample_unit": str(sample_unit),
                    "sample": str(cell_rows["sample_resolved"].astype(str).value_counts().index[0]),
                    "dataset": str(cell_rows["dataset_resolved"].astype(str).value_counts().index[0]),
                    "study": str(cell_rows["study_resolved"].astype(str).value_counts().index[0]),
                    "batch": str(cell_rows["batch_resolved"].astype(str).value_counts().index[0]),
                    "site_label": str(cell_rows["site_label"].astype(str).value_counts().index[0]),
                    "site_axis": str(cell_rows["site_axis"].astype(str).value_counts().index[0]),
                    "condition": str(cell_rows["condition_resolved"].astype(str).value_counts().index[0]),
                    "cell_type": str(celltype),
                    "n_cells": int(cell_rows.shape[0]),
                }
            )
        if not counts_rows:
            continue
        counts = np.vstack(counts_rows).astype(np.float32)
        meta = pd.DataFrame(meta_rows).set_index("sample_unit", drop=False)
        pb_by_celltype[str(celltype)] = {
            "counts": counts,
            "log_cpm": normalize_log_cpm(counts),
            "sample_meta": meta,
        }
    return pb_by_celltype


def build_binary_method(method_name: str, cfg: dict[str, Any], seed: int):
    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    max_iter = int(deep_get(core_cfg, "logistic_max_iter", default=3000))
    elastic_net_l1_ratio = float(deep_get(core_cfg, "elastic_net_l1_ratio", default=0.5))
    if method_name == "lasso_logistic":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        penalty="l1",
                        solver="saga",
                        class_weight="balanced",
                        max_iter=max_iter,
                        random_state=seed,
                    ),
                ),
            ]
        )
    if method_name == "elastic_net_logistic":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        penalty="elasticnet",
                        solver="saga",
                        l1_ratio=elastic_net_l1_ratio,
                        class_weight="balanced",
                        max_iter=max_iter,
                        random_state=seed,
                    ),
                ),
            ]
        )
    if method_name == "ridge_logistic":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        penalty="l2",
                        solver="lbfgs",
                        class_weight="balanced",
                        max_iter=max_iter,
                        random_state=seed,
                    ),
                ),
            ]
        )
    if method_name == "linear_svm":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                ("clf", LinearSVC(class_weight="balanced", random_state=seed)),
            ]
        )
    if method_name == "random_forest":
        return RandomForestClassifier(
            n_estimators=int(deep_get(core_cfg, "random_forest_estimators", default=300)),
            random_state=seed,
            n_jobs=-1,
            class_weight="balanced_subsample",
        )
    if method_name == "lda":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                ("clf", LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")),
            ]
        )
    if method_name == "naive_bayes":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                ("clf", GaussianNB()),
            ]
        )
    if method_name == "xgboost":
        if XGBClassifier is None:
            return None
        return XGBClassifier(
            n_estimators=int(deep_get(core_cfg, "xgboost_estimators", default=200)),
            max_depth=int(deep_get(core_cfg, "xgboost_max_depth", default=4)),
            learning_rate=float(deep_get(core_cfg, "xgboost_learning_rate", default=0.05)),
            subsample=0.9,
            colsample_bytree=0.9,
            objective="binary:logistic",
            eval_metric="logloss",
            random_state=seed,
            n_jobs=1,
        )
    raise KeyError(f"Unsupported core discovery method: {method_name}")


def extract_binary_feature_scores(model, method_name: str, feature_names: list[str]) -> pd.DataFrame:
    if hasattr(model, "named_steps"):
        clf = model.named_steps.get("clf", model)
    else:
        clf = model
    if method_name in {"lasso_logistic", "elastic_net_logistic", "ridge_logistic", "linear_svm"}:
        coef = np.asarray(getattr(clf, "coef_", np.zeros((1, len(feature_names)))), dtype=float)
        score = np.mean(np.abs(coef), axis=0)
    elif method_name == "lda":
        if hasattr(clf, "coef_"):
            coef = np.asarray(clf.coef_, dtype=float)
        else:
            coef = np.asarray(getattr(clf, "scalings_", np.zeros((len(feature_names), 1))), dtype=float).T
        score = np.mean(np.abs(coef), axis=0)
    elif method_name == "naive_bayes":
        theta = np.asarray(getattr(clf, "theta_", np.zeros((1, len(feature_names)))), dtype=float)
        var = np.asarray(getattr(clf, "var_", np.ones_like(theta)), dtype=float)
        centered = theta - theta.mean(axis=0, keepdims=True)
        score = np.mean(np.abs(centered) / (np.sqrt(var) + 1e-8), axis=0)
    else:
        score = np.asarray(getattr(clf, "feature_importances_", np.zeros(len(feature_names))), dtype=float)
    return pd.DataFrame({"gene": feature_names, "score": score}).sort_values("score", ascending=False).reset_index(drop=True)


def generate_group_balanced_splits(sample_meta: pd.DataFrame, y: pd.Series, cfg: dict[str, Any]) -> list[tuple[np.ndarray, np.ndarray, str]]:
    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    dataset_col = str(deep_get(core_cfg, "dataset_group_column", default="dataset"))
    if dataset_col not in sample_meta.columns:
        return []
    group_series = sample_meta[dataset_col].astype(str)
    classes = sorted(y.astype(int).unique().tolist())
    if classes != [0, 1]:
        return []

    groups_by_class = {cls: pd.Index(group_series.loc[y.astype(int) == cls].unique()) for cls in classes}
    if any(len(groups) < 2 for groups in groups_by_class.values()):
        return []

    repeats = int(deep_get(core_cfg, "stability_repeats", default=12))
    test_size = float(deep_get(core_cfg, "test_size", default=0.34))
    rng = np.random.default_rng(int(deep_get(core_cfg, "random_seed", default=20260527)))
    splits: list[tuple[np.ndarray, np.ndarray, str]] = []
    seen: set[tuple[str, ...]] = set()
    attempts = 0
    max_attempts = max(repeats * 20, 40)
    while len(splits) < repeats and attempts < max_attempts:
        attempts += 1
        test_groups_selected: list[str] = []
        for cls in classes:
            groups = groups_by_class[cls]
            holdout_n = max(1, int(round(len(groups) * test_size)))
            holdout_n = min(holdout_n, len(groups) - 1)
            chosen = rng.choice(groups.to_numpy(), size=holdout_n, replace=False)
            test_groups_selected.extend([str(x) for x in chosen.tolist()])
        test_groups = tuple(sorted(set(test_groups_selected)))
        if not test_groups or test_groups in seen:
            continue
        test_mask = group_series.isin(test_groups).to_numpy()
        train_idx = np.flatnonzero(~test_mask)
        test_idx = np.flatnonzero(test_mask)
        if train_idx.size == 0 or test_idx.size == 0:
            continue
        if y.iloc[train_idx].nunique() < 2 or y.iloc[test_idx].nunique() < 2:
            continue
        seen.add(test_groups)
        splits.append((train_idx, test_idx, "group_holdout"))
    return splits


def generate_binary_splits(sample_meta: pd.DataFrame, y: pd.Series, cfg: dict[str, Any]) -> list[tuple[np.ndarray, np.ndarray, str]]:
    group_splits = generate_group_balanced_splits(sample_meta, y, cfg)
    if group_splits:
        return group_splits

    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    repeats = int(deep_get(core_cfg, "stability_repeats", default=12))
    test_size = float(deep_get(core_cfg, "test_size", default=0.34))
    seed = int(deep_get(core_cfg, "random_seed", default=20260527))
    splitter = StratifiedShuffleSplit(n_splits=repeats, test_size=test_size, random_state=seed)
    splits: list[tuple[np.ndarray, np.ndarray, str]] = []
    try:
        dummy = np.zeros((len(y), 1), dtype=float)
        for train_idx, test_idx in splitter.split(dummy, y.astype(int).to_numpy()):
            if y.iloc[train_idx].nunique() < 2 or y.iloc[test_idx].nunique() < 2:
                continue
            splits.append((train_idx, test_idx, "stratified_shuffle"))
    except ValueError:
        pass
    return splits


def validation_row(y_true: np.ndarray, y_pred: np.ndarray, method_name: str, split_name: str) -> dict[str, Any]:
    return {
        "method": method_name,
        "split": split_name,
        "n_samples": int(len(y_true)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
    }


def run_stability_selection(
    X: pd.DataFrame,
    y: pd.Series,
    sample_meta: pd.DataFrame,
    cfg: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    requested_methods = deep_get(core_cfg, "methods", default=CORE_METHOD_ORDER)
    if not isinstance(requested_methods, list):
        requested_methods = CORE_METHOD_ORDER
    methods = [m for m in requested_methods if m in CORE_METHOD_ORDER]
    top_k = int(deep_get(core_cfg, "model_top_k_genes", default=40))
    seed = int(deep_get(core_cfg, "random_seed", default=20260527))
    feature_names = X.columns.astype(str).tolist()
    splits = generate_binary_splits(sample_meta, y, cfg)
    split_kind_counts = pd.Series([split_kind for _, _, split_kind in splits], dtype="string").value_counts().to_dict() if splits else {}
    split_metadata: dict[str, Any] = {
        "split_strategy": "group_holdout" if int(split_kind_counts.get("group_holdout", 0)) > 0 else (
            "stratified_shuffle" if int(split_kind_counts.get("stratified_shuffle", 0)) > 0 else "none"
        ),
        "n_splits_total": int(len(splits)),
        "n_splits_group_holdout": int(split_kind_counts.get("group_holdout", 0)),
        "n_splits_stratified_shuffle": int(split_kind_counts.get("stratified_shuffle", 0)),
    }

    selection_counts: dict[str, pd.Series] = {method: pd.Series(0, index=feature_names, dtype=float) for method in methods}
    fit_counts: dict[str, int] = {method: 0 for method in methods}
    full_score_frames: list[pd.DataFrame] = []
    validation_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []

    for split_idx, (train_idx, test_idx, split_kind) in enumerate(splits, start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx].astype(int)
        y_test = y.iloc[test_idx].astype(int)
        for method_offset, method in enumerate(methods, start=1):
            model = build_binary_method(method, cfg, seed + split_idx * 100 + method_offset)
            if model is None:
                status_rows.append({"method": method, "status": "skipped", "reason": "unavailable"})
                continue
            try:
                model.fit(X_train, y_train)
                y_pred = np.asarray(model.predict(X_test), dtype=int)
                validation_rows.append(validation_row(y_test.to_numpy(), y_pred, method, f"{split_kind}_{split_idx:02d}"))
                scores = extract_binary_feature_scores(model, method, feature_names)
                chosen = scores.head(min(top_k, scores.shape[0]))["gene"].astype(str)
                selection_counts[method].loc[chosen] += 1.0
                fit_counts[method] += 1
            except Exception as exc:
                status_rows.append({"method": method, "status": "failed", "reason": str(exc)})

    for method_offset, method in enumerate(methods, start=1):
        model = build_binary_method(method, cfg, seed + 5000 + method_offset)
        if model is None:
            continue
        try:
            model.fit(X, y.astype(int))
            score_df = extract_binary_feature_scores(model, method, feature_names).rename(columns={"score": f"full_score__{method}"})
            full_score_frames.append(score_df)
            status_rows.append({"method": method, "status": "ok", "fits": int(fit_counts.get(method, 0))})
        except Exception as exc:
            status_rows.append({"method": method, "status": "failed_full_fit", "reason": str(exc)})

    selection_df = pd.DataFrame(index=pd.Index(feature_names, name="gene"))
    for method in methods:
        denom = max(fit_counts.get(method, 0), 1)
        selection_df[f"selection_freq__{method}"] = selection_counts[method].reindex(feature_names).fillna(0.0) / float(denom)
        selection_df[f"fit_count__{method}"] = float(fit_counts.get(method, 0))
    if methods:
        selection_freq_cols = [f"selection_freq__{method}" for method in methods]
        selection_df["selection_freq_mean"] = selection_df[selection_freq_cols].mean(axis=1)
        selection_df["selected_methods_n"] = (selection_df[selection_freq_cols] > 0).sum(axis=1)
    else:
        selection_df["selection_freq_mean"] = 0.0
        selection_df["selected_methods_n"] = 0

    for score_df in full_score_frames:
        selection_df = selection_df.join(score_df.set_index("gene"), how="left")

    validation_df = pd.DataFrame(validation_rows)
    status_df = pd.DataFrame(status_rows).drop_duplicates()
    return selection_df.reset_index(), validation_df, status_df, split_metadata


def compute_direction_consistency(log_cpm: np.ndarray, sample_meta: pd.DataFrame, pos_mask: np.ndarray, neg_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dataset_series = sample_meta["dataset"].astype(str)
    global_effect = log_cpm[pos_mask].mean(axis=0) - log_cpm[neg_mask].mean(axis=0)
    support = np.zeros(log_cpm.shape[1], dtype=float)
    effective = np.zeros(log_cpm.shape[1], dtype=float)
    for dataset in dataset_series.unique().tolist():
        ds_mask = dataset_series.eq(str(dataset)).to_numpy()
        ds_pos = ds_mask & pos_mask
        ds_neg = ds_mask & neg_mask
        if int(ds_pos.sum()) < 1 or int(ds_neg.sum()) < 1:
            continue
        ds_effect = log_cpm[ds_pos].mean(axis=0) - log_cpm[ds_neg].mean(axis=0)
        nonzero = np.abs(ds_effect) > 1e-12
        effective += nonzero.astype(float)
        support += ((np.sign(ds_effect) == np.sign(global_effect)) & nonzero).astype(float)
    consistency = np.full(log_cpm.shape[1], np.nan, dtype=float)
    np.divide(support, effective, out=consistency, where=effective > 0)
    return consistency.astype(float), effective.astype(float)


def compute_gene_level_stats(
    counts: np.ndarray,
    sample_meta: pd.DataFrame,
    contrast_spec: dict[str, Any],
    cfg: dict[str, Any],
    gene_names: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    positive_label = str(contrast_spec["positive_label"])
    negative_label = str(contrast_spec["negative_label"])
    contrast_type = str(contrast_spec["contrast_type"])

    site_series = sample_meta["site_label"].astype(str)
    if contrast_type == "pairwise":
        included_mask = site_series.isin([positive_label, negative_label]).to_numpy()
        contrast_meta = sample_meta.loc[included_mask].copy()
        contrast_counts = counts[included_mask]
        y = contrast_meta["site_label"].astype(str).eq(positive_label).astype(int)
    else:
        contrast_meta = sample_meta.copy()
        contrast_counts = counts.copy()
        y = contrast_meta["site_label"].astype(str).eq(positive_label).astype(int)

    pos_mask = y.to_numpy(dtype=int) == 1
    neg_mask = ~pos_mask
    log_cpm = normalize_log_cpm(contrast_counts)

    mean_counts_pos = contrast_counts[pos_mask].mean(axis=0)
    mean_counts_neg = contrast_counts[neg_mask].mean(axis=0)
    mean_log_pos = log_cpm[pos_mask].mean(axis=0)
    mean_log_neg = log_cpm[neg_mask].mean(axis=0)
    effect = mean_log_pos - mean_log_neg

    var_pos = np.var(log_cpm[pos_mask], axis=0, ddof=1) if int(pos_mask.sum()) > 1 else np.zeros(log_cpm.shape[1], dtype=float)
    var_neg = np.var(log_cpm[neg_mask], axis=0, ddof=1) if int(neg_mask.sum()) > 1 else np.zeros(log_cpm.shape[1], dtype=float)
    pooled_var = ((max(int(pos_mask.sum()) - 1, 0) * var_pos) + (max(int(neg_mask.sum()) - 1, 0) * var_neg)) / max(int(pos_mask.sum() + neg_mask.sum()) - 2, 1)
    cohen_d = np.divide(effect, np.sqrt(np.maximum(pooled_var, 1e-12)), out=np.zeros_like(effect), where=np.sqrt(np.maximum(pooled_var, 1e-12)) > 0)

    test_res = stats.ttest_ind(log_cpm[pos_mask], log_cpm[neg_mask], axis=0, equal_var=False, nan_policy="omit")
    p_values = np.nan_to_num(np.asarray(test_res.pvalue, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
    adj_p = bh_adjust(p_values)

    direction_consistency, effective_datasets = compute_direction_consistency(log_cpm, contrast_meta, pos_mask, neg_mask)
    total_counts = contrast_counts.sum(axis=0)
    detected_samples = (contrast_counts > 0).sum(axis=0)
    detected_pos = (contrast_counts[pos_mask] > 0).mean(axis=0)
    detected_neg = (contrast_counts[neg_mask] > 0).mean(axis=0)

    stats_df = pd.DataFrame(
        {
            "gene": gene_names,
            "mean_counts_positive": mean_counts_pos,
            "mean_counts_negative": mean_counts_neg,
            "mean_log_cpm_positive": mean_log_pos,
            "mean_log_cpm_negative": mean_log_neg,
            "effect_log_cpm": effect,
            "log2_fc": np.log2(mean_counts_pos + 1.0) - np.log2(mean_counts_neg + 1.0),
            "cohen_d": cohen_d,
            "p_value": p_values,
            "adj_p": adj_p,
            "total_count": total_counts,
            "detected_samples": detected_samples,
            "detected_fraction_positive": detected_pos,
            "detected_fraction_negative": detected_neg,
            "direction_consistency": direction_consistency,
            "n_effective_datasets": effective_datasets,
        }
    )

    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    min_total_count = float(deep_get(core_cfg, "min_total_count_per_gene", default=10.0))
    min_detected_samples = int(deep_get(core_cfg, "min_detected_samples_per_gene", default=2))
    candidate_mask = (stats_df["total_count"] >= min_total_count) & (stats_df["detected_samples"] >= min_detected_samples)
    max_genes_for_models = deep_get(core_cfg, "max_genes_for_models", default=300)
    if max_genes_for_models in {None, "", "null"}:
        max_genes_for_models = None
    else:
        max_genes_for_models = int(max_genes_for_models)

    model_candidates = stats_df.loc[candidate_mask].copy()
    model_candidates["model_rank_abs_effect"] = model_candidates["cohen_d"].abs().rank(method="average", ascending=False)
    model_candidates["model_rank_significance"] = model_candidates["adj_p"].rank(method="average", ascending=True)
    model_candidates = model_candidates.sort_values(["model_rank_abs_effect", "model_rank_significance", "gene"], ascending=[True, True, True])
    if max_genes_for_models is not None and int(model_candidates.shape[0]) > max_genes_for_models:
        model_candidates = model_candidates.head(max_genes_for_models)

    model_gene_list = model_candidates["gene"].astype(str).tolist()
    if model_gene_list:
        X_model = pd.DataFrame(log_cpm, index=contrast_meta.index.astype(str), columns=gene_names).loc[:, model_gene_list]
        selection_df, validation_df, status_df, split_metadata = run_stability_selection(X_model, y, contrast_meta, cfg)
    else:
        selection_df = pd.DataFrame({"gene": gene_names, "selection_freq_mean": 0.0, "selected_methods_n": 0})
        validation_df = pd.DataFrame()
        status_df = pd.DataFrame()
        split_metadata = {
            "split_strategy": "none",
            "n_splits_total": 0,
            "n_splits_group_holdout": 0,
            "n_splits_stratified_shuffle": 0,
        }

    stats_df = stats_df.merge(selection_df, on="gene", how="left")
    selection_cols = [col for col in stats_df.columns if col.startswith("selection_freq__") or col.startswith("fit_count__") or col.startswith("full_score__")]
    fill_zero_cols = [col for col in selection_cols + ["selection_freq_mean", "selected_methods_n"] if col in stats_df.columns]
    if fill_zero_cols:
        stats_df[fill_zero_cols] = stats_df[fill_zero_cols].fillna(0.0)

    weights = deep_get(core_cfg, "consensus_weights", default={}) or {}
    weight_effect = float(weights.get("effect", 0.35))
    weight_significance = float(weights.get("significance", 0.25))
    weight_stability = float(weights.get("stability", 0.25))
    weight_dataset = float(weights.get("dataset_consistency", 0.15))

    stats_df["effect_rank_pct"] = stats_df["cohen_d"].abs().rank(pct=True, method="average").fillna(0.0)
    stats_df["significance_rank_pct"] = (-np.log10(stats_df["adj_p"].clip(lower=1e-300))).rank(pct=True, method="average").fillna(0.0)
    stats_df["dataset_consistency_score"] = np.where(stats_df["n_effective_datasets"] >= 2, stats_df["direction_consistency"].fillna(0.5), 0.5)
    stats_df["consensus_score"] = (
        weight_effect * stats_df["effect_rank_pct"]
        + weight_significance * stats_df["significance_rank_pct"]
        + weight_stability * stats_df["selection_freq_mean"].fillna(0.0)
        + weight_dataset * stats_df["dataset_consistency_score"].fillna(0.5)
    )

    stable_cfg = deep_get(core_cfg, "stable_gene", default={}) or {}
    stable_mask = candidate_mask.reindex(stats_df.index, fill_value=False)
    stable_mask &= stats_df["adj_p"] <= float(stable_cfg.get("max_adj_p", 0.1))
    stable_mask &= stats_df["log2_fc"].abs() >= float(stable_cfg.get("min_abs_log2_fc", 0.5))
    stable_mask &= stats_df["selection_freq_mean"].fillna(0.0) >= float(stable_cfg.get("min_selection_freq", 0.35))
    stable_mask &= stats_df["selected_methods_n"].fillna(0.0) >= int(stable_cfg.get("min_selected_methods", 1))
    if bool((stats_df["n_effective_datasets"] >= 2).any()):
        stable_mask &= np.where(
            stats_df["n_effective_datasets"] >= 2,
            stats_df["direction_consistency"].fillna(0.0) >= float(stable_cfg.get("min_direction_consistency", 0.6)),
            True,
        )
    stats_df["stable_core"] = stable_mask.astype(bool)
    stats_df["stable_single_method"] = (stats_df["stable_core"] & (stats_df["selected_methods_n"].fillna(0.0) == 1.0)).astype(bool)

    stats_df["contrast_id"] = str(contrast_spec["contrast_id"])
    stats_df["contrast_type"] = contrast_type
    stats_df["positive_label"] = positive_label
    stats_df["negative_label"] = negative_label
    stats_df["n_samples_positive"] = int(pos_mask.sum())
    stats_df["n_samples_negative"] = int(neg_mask.sum())
    stats_df["enriched_site"] = np.where(
        stats_df["log2_fc"] >= 0,
        positive_label,
        negative_label if contrast_type == "pairwise" else "rest",
    )
    stats_df = stats_df.sort_values(["stable_core", "consensus_score", "cohen_d"], ascending=[False, False, False]).reset_index(drop=True)
    return stats_df, validation_df, status_df, split_metadata


def summarize_core_tables(stable_core_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if stable_core_df.empty:
        site_level = pd.DataFrame(
            columns=[
                "cell_type",
                "gene",
                "enriched_site",
                "n_supporting_contrasts",
                "supporting_contrasts",
                "best_consensus_score",
                "mean_selection_freq",
                "min_adj_p",
                "mean_abs_log2_fc",
            ]
        )
        gene_support = pd.DataFrame(
            columns=[
                "cell_type",
                "gene",
                "n_target_sites",
                "total_supporting_contrasts",
                "supported_sites",
                "best_consensus_score",
                "mean_selection_freq",
            ]
        )
        shared_core = gene_support.copy()
        site_specific = pd.DataFrame(
            columns=[
                "cell_type",
                "gene",
                "enriched_site",
                "n_supporting_contrasts",
                "supporting_contrasts",
                "best_consensus_score",
                "mean_selection_freq",
                "min_adj_p",
                "mean_abs_log2_fc",
                "n_target_sites",
            ]
        )
        return site_level, gene_support, shared_core, site_specific

    site_level = (
        stable_core_df.loc[stable_core_df["enriched_site"].astype(str) != "rest"]
        .groupby(["cell_type", "gene", "enriched_site"], observed=True)
        .agg(
            n_supporting_contrasts=("contrast_id", "nunique"),
            supporting_contrasts=("contrast_id", lambda s: "|".join(sorted(set(map(str, s))))),
            best_consensus_score=("consensus_score", "max"),
            mean_selection_freq=("selection_freq_mean", "mean"),
            min_adj_p=("adj_p", "min"),
            mean_abs_log2_fc=("log2_fc", lambda s: float(np.mean(np.abs(np.asarray(s, dtype=float))))),
        )
        .reset_index()
        .sort_values(["n_supporting_contrasts", "best_consensus_score", "mean_abs_log2_fc"], ascending=[False, False, False])
    )
    gene_support = (
        site_level.groupby(["cell_type", "gene"], observed=True)
        .agg(
            n_target_sites=("enriched_site", "nunique"),
            total_supporting_contrasts=("n_supporting_contrasts", "sum"),
            supported_sites=("enriched_site", lambda s: "|".join(sorted(set(map(str, s))))),
            best_consensus_score=("best_consensus_score", "max"),
            mean_selection_freq=("mean_selection_freq", "mean"),
        )
        .reset_index()
        .sort_values(["n_target_sites", "total_supporting_contrasts", "best_consensus_score"], ascending=[False, False, False])
    )
    shared_core = gene_support.loc[gene_support["n_target_sites"] >= 2].copy()
    site_specific = site_level.merge(gene_support[["cell_type", "gene", "n_target_sites"]], on=["cell_type", "gene"], how="left")
    site_specific = site_specific.loc[site_specific["n_target_sites"] == 1].copy()
    return site_level, gene_support, shared_core, site_specific


def summarize_stable_method_support(stable_core_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    distribution_cols = ["selected_methods_n", "n_genes", "fraction_of_stable"]
    by_contrast_cols = ["cell_type", "contrast_id", "contrast_type", "n_stable_genes", "n_single_method_stable", "single_method_ratio"]
    if stable_core_df.empty:
        return pd.DataFrame(columns=distribution_cols), pd.DataFrame(columns=by_contrast_cols)

    stable = stable_core_df.copy()
    stable["selected_methods_n"] = stable["selected_methods_n"].fillna(0.0).astype(int)
    distribution_df = (
        stable.groupby("selected_methods_n", observed=True)
        .size()
        .rename("n_genes")
        .reset_index()
        .sort_values("selected_methods_n")
    )
    total = float(distribution_df["n_genes"].sum())
    distribution_df["fraction_of_stable"] = np.where(total > 0, distribution_df["n_genes"] / total, 0.0)

    by_contrast_df = (
        stable.groupby(["cell_type", "contrast_id", "contrast_type"], observed=True)
        .agg(
            n_stable_genes=("gene", "nunique"),
            n_single_method_stable=("selected_methods_n", lambda s: int((np.asarray(s, dtype=int) == 1).sum())),
        )
        .reset_index()
    )
    by_contrast_df["single_method_ratio"] = np.where(
        by_contrast_df["n_stable_genes"] > 0,
        by_contrast_df["n_single_method_stable"] / by_contrast_df["n_stable_genes"],
        np.nan,
    )
    return distribution_df, by_contrast_df


def run_core_discovery(
    input_h5ad: Path | str,
    cfg: dict[str, Any],
    output_dir: Path | str,
    cell_type_key_override: str | None = None,
) -> dict[str, Any]:
    output_dir = ensure_dir(output_dir)
    pseudobulk_dir = ensure_dir(Path(output_dir) / "pseudobulk")
    pairwise_dir = ensure_dir(Path(output_dir) / "pairwise")
    one_vs_rest_dir = ensure_dir(Path(output_dir) / "one_vs_rest")
    summary_dir = ensure_dir(Path(output_dir) / "summaries")

    input_h5ad = Path(input_h5ad)
    if not input_h5ad.exists():
        raise FileNotFoundError(f"Core discovery input h5ad not found: {input_h5ad}")

    adata = ad.read_h5ad(input_h5ad)
    counts_layer = str(deep_get(cfg, "core_discovery", "counts_layer", default=deep_get(cfg, "data_contract", "counts_layer", default="counts")))
    counts_info = ensure_counts_layer(adata, counts_layer)
    prepared_obs, contract = prepare_obs_contract(
        adata.obs.copy(),
        list(adata.obsm.keys()),
        cfg,
        cell_type_key_override=cell_type_key_override or deep_get(cfg, "core_discovery", "cell_type_key", default=deep_get(cfg, "feature_export", "cell_type_key", default=None)),
    )
    adata.obs = prepared_obs

    pseudobulk_by_celltype = aggregate_pseudobulk_by_celltype(adata, adata.obs, contract, cfg)
    gene_names = adata.var_names.astype(str).tolist()
    core_cfg = deep_get(cfg, "core_discovery", default={}) or {}
    export_pseudobulk_tables = bool(deep_get(core_cfg, "export_pseudobulk_tables", default=True))

    celltype_rows: list[dict[str, Any]] = []
    contrast_summary_rows: list[dict[str, Any]] = []
    stable_core_frames: list[pd.DataFrame] = []
    validation_frames: list[pd.DataFrame] = []
    method_status_frames: list[pd.DataFrame] = []
    stable_core_columns: list[str] | None = None

    for celltype, payload in pseudobulk_by_celltype.items():
        counts = payload["counts"]
        sample_meta = payload["sample_meta"].copy()
        celltype_token = slugify_token(celltype)
        celltype_rows.append(
            {
                "cell_type": celltype,
                "n_samples": int(sample_meta.shape[0]),
                "n_sites": int(sample_meta["site_label"].astype(str).nunique()),
                "n_datasets": int(sample_meta["dataset"].astype(str).nunique()),
                "min_cells_per_sample": int(sample_meta["n_cells"].min()),
                "max_cells_per_sample": int(sample_meta["n_cells"].max()),
            }
        )
        sample_meta.to_csv(pseudobulk_dir / f"{celltype_token}__sample_metadata.tsv", sep="\t", index=False)
        if export_pseudobulk_tables:
            pd.DataFrame(counts, index=sample_meta.index.astype(str), columns=gene_names).to_csv(
                pseudobulk_dir / f"{celltype_token}__counts.tsv.gz", sep="\t", compression="gzip", index_label="sample_unit"
            )

        contrasts = build_site_contrasts(sample_meta, cfg)
        for contrast_spec in contrasts:
            stats_df, validation_df, status_df, split_metadata = compute_gene_level_stats(counts, sample_meta, contrast_spec, cfg, gene_names)
            if stable_core_columns is None:
                stable_core_columns = ["cell_type", *stats_df.columns.astype(str).tolist()]
            contrast_token = slugify_token(contrast_spec["contrast_id"])
            target_dir = pairwise_dir if contrast_spec["contrast_type"] == "pairwise" else one_vs_rest_dir
            result_path = target_dir / f"{celltype_token}__{contrast_token}.tsv.gz"
            stats_df.to_csv(result_path, sep="\t", index=False, compression="gzip")
            if not validation_df.empty:
                validation_df = validation_df.copy()
                validation_df.insert(0, "cell_type", celltype)
                validation_df.insert(1, "contrast_id", contrast_spec["contrast_id"])
                validation_frames.append(validation_df)
            if not status_df.empty:
                status_df = status_df.copy()
                status_df.insert(0, "cell_type", celltype)
                status_df.insert(1, "contrast_id", contrast_spec["contrast_id"])
                method_status_frames.append(status_df)
            stable_rows = stats_df.loc[stats_df["stable_core"]].copy()
            stable_single_method_n = int(stable_rows["selected_methods_n"].fillna(0.0).eq(1.0).sum()) if not stable_rows.empty else 0
            stable_single_method_ratio = float(stable_single_method_n / stable_rows.shape[0]) if not stable_rows.empty else None
            if not stable_rows.empty:
                stable_rows.insert(0, "cell_type", celltype)
                stable_core_frames.append(stable_rows)

            contrast_summary_rows.append(
                {
                    "cell_type": celltype,
                    "contrast_id": contrast_spec["contrast_id"],
                    "contrast_type": contrast_spec["contrast_type"],
                    "positive_label": contrast_spec["positive_label"],
                    "negative_label": contrast_spec["negative_label"],
                    "split_strategy": split_metadata.get("split_strategy"),
                    "n_splits_total": int(split_metadata.get("n_splits_total", 0)),
                    "n_splits_group_holdout": int(split_metadata.get("n_splits_group_holdout", 0)),
                    "n_splits_stratified_shuffle": int(split_metadata.get("n_splits_stratified_shuffle", 0)),
                    "n_samples_total": int(stats_df[["n_samples_positive", "n_samples_negative"]].iloc[0].sum()) if not stats_df.empty else int(sample_meta.shape[0]),
                    "n_samples_positive": int(stats_df["n_samples_positive"].iloc[0]) if not stats_df.empty else None,
                    "n_samples_negative": int(stats_df["n_samples_negative"].iloc[0]) if not stats_df.empty else None,
                    "n_candidate_genes": int(stats_df.shape[0]),
                    "n_stable_core_genes": int(stats_df["stable_core"].sum()) if "stable_core" in stats_df.columns else 0,
                    "n_single_method_stable": stable_single_method_n,
                    "single_method_stable_ratio": stable_single_method_ratio,
                    "mean_balanced_accuracy": float(validation_df["balanced_accuracy"].mean()) if not validation_df.empty else None,
                    "mean_macro_f1": float(validation_df["macro_f1"].mean()) if not validation_df.empty else None,
                    "top_core_genes": "|".join(stable_rows.head(10)["gene"].astype(str).tolist()) if not stable_rows.empty else "",
                    "result_path": str(result_path),
                }
            )

    celltype_summary_columns = ["cell_type", "n_samples", "n_sites", "n_datasets", "min_cells_per_sample", "max_cells_per_sample"]
    contrast_summary_columns = [
        "cell_type",
        "contrast_id",
        "contrast_type",
        "positive_label",
        "negative_label",
        "split_strategy",
        "n_splits_total",
        "n_splits_group_holdout",
        "n_splits_stratified_shuffle",
        "n_samples_total",
        "n_samples_positive",
        "n_samples_negative",
        "n_candidate_genes",
        "n_stable_core_genes",
        "n_single_method_stable",
        "single_method_stable_ratio",
        "mean_balanced_accuracy",
        "mean_macro_f1",
        "top_core_genes",
        "result_path",
    ]
    celltype_summary_df = pd.DataFrame(celltype_rows, columns=celltype_summary_columns)
    if not celltype_summary_df.empty:
        celltype_summary_df = celltype_summary_df.sort_values(["n_samples", "cell_type"], ascending=[False, True])
    contrast_summary_df = pd.DataFrame(contrast_summary_rows, columns=contrast_summary_columns)
    if not contrast_summary_df.empty:
        contrast_summary_df = contrast_summary_df.sort_values(["n_stable_core_genes", "cell_type", "contrast_id"], ascending=[False, True, True])
    stable_core_df = pd.concat(stable_core_frames, ignore_index=True) if stable_core_frames else pd.DataFrame(columns=stable_core_columns or ["cell_type"])
    validation_df = pd.concat(validation_frames, ignore_index=True) if validation_frames else pd.DataFrame(columns=["cell_type", "contrast_id", "method", "split", "n_samples", "accuracy", "balanced_accuracy", "macro_f1"])
    method_status_df = pd.concat(method_status_frames, ignore_index=True) if method_status_frames else pd.DataFrame(columns=["cell_type", "contrast_id", "method", "status", "reason", "fits"])

    celltype_summary_df.to_csv(summary_dir / "celltype_pseudobulk_summary.tsv", sep="\t", index=False)
    contrast_summary_df.to_csv(summary_dir / "contrast_summary.tsv", sep="\t", index=False)
    stable_core_df.to_csv(summary_dir / "stable_core_genes.tsv", sep="\t", index=False)
    validation_df.to_csv(summary_dir / "stability_validation.tsv", sep="\t", index=False)
    method_status_df.to_csv(summary_dir / "stability_method_status.tsv", sep="\t", index=False)

    site_level_df, gene_support_df, shared_core_df, site_specific_df = summarize_core_tables(stable_core_df) if not stable_core_df.empty else (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    stable_method_distribution_df, stable_method_by_contrast_df = summarize_stable_method_support(stable_core_df)
    site_level_df.to_csv(summary_dir / "site_core_summary.tsv", sep="\t", index=False)
    gene_support_df.to_csv(summary_dir / "gene_support_summary.tsv", sep="\t", index=False)
    shared_core_df.to_csv(summary_dir / "shared_core_genes.tsv", sep="\t", index=False)
    site_specific_df.to_csv(summary_dir / "site_specific_core_genes.tsv", sep="\t", index=False)
    stable_method_distribution_df.to_csv(summary_dir / "stable_method_support_distribution.tsv", sep="\t", index=False)
    stable_method_by_contrast_df.to_csv(summary_dir / "stable_method_support_by_contrast.tsv", sep="\t", index=False)

    manifest = {
        "input_h5ad": str(input_h5ad),
        "output_dir": str(output_dir),
        "counts_info": counts_info,
        "resolved_columns": {
            k: contract.get(k)
            for k in ["sample_col", "sample_unit_col", "dataset_col", "study_col", "batch_col", "tissue_col", "condition_col", "cell_type_col", "latent_key"]
        },
        "n_analysis_cells": int(contract["n_analysis_cells"]),
        "n_celltypes_pseudobulked": int(len(pseudobulk_by_celltype)),
        "n_contrasts": int(contrast_summary_df.shape[0]) if not contrast_summary_df.empty else 0,
        "n_stable_core_rows": int(stable_core_df.shape[0]) if not stable_core_df.empty else 0,
        "n_shared_core_rows": int(shared_core_df.shape[0]) if not shared_core_df.empty else 0,
        "n_site_specific_rows": int(site_specific_df.shape[0]) if not site_specific_df.empty else 0,
        "methods": [m for m in deep_get(core_cfg, "methods", default=CORE_METHOD_ORDER) if m in CORE_METHOD_ORDER],
        "xgboost_available": bool(XGBClassifier is not None),
    }
    write_json(manifest, Path(output_dir) / "core_manifest.json")
    return manifest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run contrast/stability core discovery for normal-airway site comparisons")
    parser.add_argument("--input-h5ad", default=None, help="Input h5ad path")
    parser.add_argument("--config-yaml", default=None, help="Config YAML path")
    parser.add_argument("--output-dir", default=None, help="Explicit output directory")
    parser.add_argument("--cell-type-key", default=None, help="Override cell type key")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    cfg = load_config(args.config_yaml)
    input_h5ad = Path(args.input_h5ad or deep_get(cfg, "data_contract", "input_h5ad", default=""))
    output_dir = Path(args.output_dir) if args.output_dir else ensure_dir(Path(deep_get(cfg, "run", "output_root", default=str(DEFAULT_OUTPUT_ROOT))) / f"normal_airway_ml_core_discovery_{timestamp_slug()}")
    manifest = run_core_discovery(
        input_h5ad=input_h5ad,
        cfg=cfg,
        output_dir=output_dir,
        cell_type_key_override=args.cell_type_key,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
