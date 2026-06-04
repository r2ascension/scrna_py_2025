#!/usr/bin/env python3
"""Train sample-level normal-airway site classifiers with grouped validation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, f1_score
from sklearn.model_selection import GroupKFold, KFold, LeaveOneGroupOut, StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import LinearSVC

from normal_airway_ml_common_20260527 import DEFAULT_OUTPUT_ROOT, deep_get, ensure_dir, load_config, split_feature_metadata, timestamp_slug, write_json

try:  # pragma: no cover - optional
    from xgboost import XGBClassifier
except Exception:  # pragma: no cover
    XGBClassifier = None


MODEL_ORDER = [
    "logistic_regression",
    "ridge_logistic",
    "lasso_logistic",
    "elastic_net_logistic",
    "linear_svm",
    "lda",
    "naive_bayes",
    "random_forest",
    "xgboost",
]


def build_logistic_pipeline(model_name: str, cfg: dict[str, Any], seed: int) -> Pipeline:
    penalty = "l2"
    solver = "lbfgs"
    clf_kwargs: dict[str, Any] = {}
    if model_name == "lasso_logistic":
        penalty = "l1"
        solver = "saga"
    elif model_name == "elastic_net_logistic":
        penalty = "elasticnet"
        solver = "saga"
        clf_kwargs["l1_ratio"] = float(deep_get(cfg, "train", "elastic_net_l1_ratio", default=0.5))

    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    max_iter=int(deep_get(cfg, "train", "logistic_max_iter", default=2000)),
                    random_state=seed,
                    class_weight="balanced",
                    penalty=penalty,
                    solver=solver,
                    **clf_kwargs,
                ),
            ),
        ]
    )


def load_feature_bundle(feature_dir: Path | str) -> tuple[pd.DataFrame, pd.DataFrame]:
    feature_dir = Path(feature_dir)
    feature_table = pd.read_csv(feature_dir / "sample_features_wide.tsv.gz", sep="\t", index_col=0)
    sample_meta, feature_matrix = split_feature_metadata(feature_table)
    sample_meta.index = feature_table.index.astype(str)
    feature_matrix.index = feature_table.index.astype(str)
    return sample_meta, feature_matrix


def build_model(model_name: str, cfg: dict[str, Any], feature_names: list[str]):
    seed = int(deep_get(cfg, "train", "random_seed", default=20260527))
    if model_name in {"logistic_regression", "ridge_logistic", "lasso_logistic", "elastic_net_logistic"}:
        return build_logistic_pipeline(model_name, cfg, seed)
    if model_name == "linear_svm":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                ("clf", LinearSVC(class_weight="balanced", random_state=seed)),
            ]
        )
    if model_name == "lda":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                ("clf", LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")),
            ]
        )
    if model_name == "naive_bayes":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scaler", StandardScaler()),
                ("clf", GaussianNB()),
            ]
        )
    if model_name == "random_forest":
        return RandomForestClassifier(
            n_estimators=int(deep_get(cfg, "train", "random_forest_estimators", default=300)),
            random_state=seed,
            n_jobs=-1,
            class_weight="balanced_subsample",
        )
    if model_name == "xgboost":
        if XGBClassifier is None:
            return None
        return XGBClassifier(
            n_estimators=250,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.9,
            colsample_bytree=0.9,
            objective="multi:softprob",
            eval_metric="mlogloss",
            random_state=seed,
            n_jobs=1,
        )
    raise KeyError(f"Unsupported model_name: {model_name}")


def metric_row(y_true: np.ndarray, y_pred: np.ndarray, model_name: str, split_name: str) -> dict[str, Any]:
    labels = sorted({str(x) for x in y_true} | {str(x) for x in y_pred})
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    recalls = []
    for idx in range(cm.shape[0]):
        denom = float(cm[idx, :].sum())
        if denom > 0:
            recalls.append(float(cm[idx, idx]) / denom)
    balanced_accuracy = float(np.mean(recalls)) if recalls else 0.0
    return {
        "model": model_name,
        "split": split_name,
        "n_samples": int(len(y_true)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": balanced_accuracy,
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", labels=labels, zero_division=0)),
    }


def resolve_splitter(sample_meta: pd.DataFrame, y: pd.Series, cfg: dict[str, Any]):
    requested_splits = int(deep_get(cfg, "train", "n_splits", default=4))
    dataset_col = deep_get(cfg, "train", "dataset_group_column", default="dataset")
    if dataset_col in sample_meta.columns:
        groups = sample_meta[dataset_col].astype(str)
        n_groups = groups.nunique()
        if n_groups >= 2:
            n_splits = min(requested_splits, n_groups)
            if n_splits >= 2:
                candidate = GroupKFold(n_splits=n_splits)
                dummy = np.zeros((len(y), 1), dtype=float)
                feasible_splits = 0
                for train_idx, test_idx in candidate.split(dummy, y, groups.to_numpy()):
                    if y.iloc[train_idx].nunique() >= 2 and y.iloc[test_idx].nunique() >= 1:
                        feasible_splits += 1
                if feasible_splits >= 2:
                    return candidate, groups.to_numpy(), f"GroupKFold[{dataset_col}]"
    min_class = int(y.value_counts().min()) if not y.empty else 0
    if min_class >= 2 and len(y) >= 4:
        n_splits = min(requested_splits, min_class, len(y))
        if n_splits >= 2:
            return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=int(deep_get(cfg, "train", "random_seed", default=20260527))), None, "StratifiedKFold"
    n_splits = min(max(2, requested_splits), len(y))
    return KFold(n_splits=n_splits, shuffle=True, random_state=int(deep_get(cfg, "train", "random_seed", default=20260527))), None, "KFold"


def crossval_predictions(X: pd.DataFrame, y: pd.Series, sample_meta: pd.DataFrame, cfg: dict[str, Any], model_name: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    splitter, groups, splitter_name = resolve_splitter(sample_meta, y, cfg)
    pred_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    labels = sorted(y.unique().tolist())

    split_iter = splitter.split(X, y, groups) if groups is not None else splitter.split(X, y)
    for fold_id, (train_idx, test_idx) in enumerate(split_iter, start=1):
        if y.iloc[train_idx].nunique() < 2 or y.iloc[test_idx].nunique() < 1:
            continue
        model = build_model(model_name, cfg, list(X.columns))
        if model is None:
            break
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        model.fit(X_train, y_train)
        y_pred = pd.Series(model.predict(X_test), index=X_test.index, dtype="string")
        metric_rows.append(metric_row(y_test.to_numpy(), y_pred.to_numpy(), model_name, f"cv_fold_{fold_id}"))
        fold_df = pd.DataFrame(
            {
                "sample": X_test.index.astype(str),
                "truth": y_test.astype(str).values,
                "pred": y_pred.astype(str).values,
                "model": model_name,
                "split": f"cv_fold_{fold_id}",
                "dataset": sample_meta.loc[X_test.index, "dataset"].astype(str).values if "dataset" in sample_meta.columns else "dataset_missing",
            }
        )
        pred_rows.append(fold_df)

    predictions = pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame(columns=["sample", "truth", "pred", "model", "split", "dataset"])
    metrics = pd.DataFrame(metric_rows)
    summary = {
        "splitter": splitter_name,
        "n_splits": int(metrics["split"].nunique()) if not metrics.empty else 0,
        "labels": labels,
    }
    return predictions, metrics, summary


def leave_one_dataset_out_predictions(X: pd.DataFrame, y: pd.Series, sample_meta: pd.DataFrame, cfg: dict[str, Any], model_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "dataset" not in sample_meta.columns or sample_meta["dataset"].astype(str).nunique() < 2:
        return pd.DataFrame(), pd.DataFrame()
    logo = LeaveOneGroupOut()
    groups = sample_meta["dataset"].astype(str).to_numpy()
    pred_rows: list[pd.DataFrame] = []
    metric_rows: list[dict[str, Any]] = []
    for fold_id, (train_idx, test_idx) in enumerate(logo.split(X, y, groups), start=1):
        if y.iloc[train_idx].nunique() < 2 or y.iloc[test_idx].nunique() < 1:
            continue
        model = build_model(model_name, cfg, list(X.columns))
        if model is None:
            break
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        model.fit(X_train, y_train)
        y_pred = pd.Series(model.predict(X_test), index=X_test.index, dtype="string")
        dataset_name = safe_dataset_name(sample_meta.loc[X_test.index, "dataset"])
        metric_rows.append(metric_row(y_test.to_numpy(), y_pred.to_numpy(), model_name, f"leave_one_dataset_out::{dataset_name}"))
        pred_rows.append(
            pd.DataFrame(
                {
                    "sample": X_test.index.astype(str),
                    "truth": y_test.astype(str).values,
                    "pred": y_pred.astype(str).values,
                    "model": model_name,
                    "split": f"leave_one_dataset_out::{dataset_name}",
                    "dataset": sample_meta.loc[X_test.index, "dataset"].astype(str).values,
                }
            )
        )
    return (
        pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame(),
        pd.DataFrame(metric_rows),
    )


def safe_dataset_name(series: pd.Series) -> str:
    if series.empty:
        return "dataset_missing"
    return str(series.astype(str).value_counts().index[0])


def compute_feature_importance(model, feature_names: list[str], model_name: str) -> pd.DataFrame:
    if model_name in {"logistic_regression", "ridge_logistic", "lasso_logistic", "elastic_net_logistic", "linear_svm"}:
        clf = model.named_steps["clf"] if hasattr(model, "named_steps") else model
        coef = np.asarray(clf.coef_, dtype=float)
        importance = np.mean(np.abs(coef), axis=0)
    elif model_name == "lda":
        clf = model.named_steps["clf"] if hasattr(model, "named_steps") else model
        if hasattr(clf, "coef_"):
            coef = np.asarray(clf.coef_, dtype=float)
        else:
            coef = np.asarray(getattr(clf, "scalings_", np.zeros((len(feature_names), 1))), dtype=float).T
        importance = np.mean(np.abs(coef), axis=0)
    elif model_name == "naive_bayes":
        clf = model.named_steps["clf"] if hasattr(model, "named_steps") else model
        theta = np.asarray(getattr(clf, "theta_", np.zeros((1, len(feature_names)))), dtype=float)
        var = np.asarray(getattr(clf, "var_", np.ones_like(theta)), dtype=float)
        centered = theta - theta.mean(axis=0, keepdims=True)
        importance = np.mean(np.abs(centered) / (np.sqrt(var) + 1e-8), axis=0)
    else:
        importance = np.asarray(getattr(model, "feature_importances_", np.zeros(len(feature_names))), dtype=float)
    df = pd.DataFrame({"feature": feature_names, "importance": importance})
    return df.sort_values("importance", ascending=False).reset_index(drop=True)


def dataset_only_negative_control(sample_meta: pd.DataFrame, y: pd.Series, cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "dataset" not in sample_meta.columns:
        return pd.DataFrame(), pd.DataFrame()
    X = sample_meta.loc[:, ["dataset"]].copy()
    preprocessor = ColumnTransformer(
        transformers=[("dataset", OneHotEncoder(handle_unknown="ignore"), ["dataset"])],
        remainder="drop",
    )
    model = Pipeline(
        steps=[
            ("pre", preprocessor),
            ("clf", LogisticRegression(max_iter=int(deep_get(cfg, "train", "logistic_max_iter", default=2000)), class_weight="balanced")),
        ]
    )
    splitter, groups, splitter_name = resolve_splitter(sample_meta, y, cfg)
    pred_rows: list[pd.DataFrame] = []
    metric_rows: list[dict[str, Any]] = []
    split_iter = splitter.split(X, y, groups) if groups is not None else splitter.split(X, y)
    for fold_id, (train_idx, test_idx) in enumerate(split_iter, start=1):
        model.fit(X.iloc[train_idx], y.iloc[train_idx])
        y_pred = pd.Series(model.predict(X.iloc[test_idx]), index=X.iloc[test_idx].index)
        metric_rows.append(metric_row(y.iloc[test_idx].to_numpy(), y_pred.to_numpy(), "dataset_only_control", f"{splitter_name}_fold_{fold_id}"))
        pred_rows.append(
            pd.DataFrame(
                {
                    "sample": X.iloc[test_idx].index.astype(str),
                    "truth": y.iloc[test_idx].astype(str).values,
                    "pred": y_pred.astype(str).values,
                    "model": "dataset_only_control",
                    "split": f"{splitter_name}_fold_{fold_id}",
                    "dataset": sample_meta.loc[X.iloc[test_idx].index, "dataset"].astype(str).values,
                }
            )
        )
    return pd.concat(pred_rows, ignore_index=True), pd.DataFrame(metric_rows)


def permutation_control(X: pd.DataFrame, y: pd.Series, sample_meta: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    repeats = int(deep_get(cfg, "train", "permutation_repeats", default=10))
    rng = np.random.default_rng(int(deep_get(cfg, "train", "random_seed", default=20260527)))
    dataset_col = "dataset" if "dataset" in sample_meta.columns else None
    rows: list[dict[str, Any]] = []
    for rep in range(repeats):
        y_perm = y.copy()
        if dataset_col is not None:
            permuted = []
            for _, sub in sample_meta.groupby(dataset_col, observed=True):
                idx = sub.index
                vals = y.loc[idx].to_numpy().copy()
                rng.shuffle(vals)
                permuted.append(pd.Series(vals, index=idx))
            y_perm = pd.concat(permuted).reindex(y.index)
        else:
            vals = y_perm.to_numpy().copy()
            rng.shuffle(vals)
            y_perm = pd.Series(vals, index=y.index)
        preds, metrics, _ = crossval_predictions(X, y_perm.astype(str), sample_meta, cfg, "logistic_regression")
        if not metrics.empty:
            rows.append(
                {
                    "repeat": rep + 1,
                    "mean_accuracy": float(metrics["accuracy"].mean()),
                    "mean_balanced_accuracy": float(metrics["balanced_accuracy"].mean()),
                    "mean_macro_f1": float(metrics["macro_f1"].mean()),
                    "n_prediction_rows": int(preds.shape[0]),
                }
            )
    return pd.DataFrame(rows)


def run_training(feature_dir: Path | str, cfg: dict[str, Any], output_dir: Path | str) -> dict[str, Any]:
    output_dir = ensure_dir(output_dir)
    model_dir = ensure_dir(Path(output_dir) / "model")
    sample_meta, X = load_feature_bundle(feature_dir)
    label_column = str(deep_get(cfg, "train", "label_column", default="site_label"))
    if label_column not in sample_meta.columns:
        raise KeyError(f"Label column missing from sample metadata: {label_column}")
    y = sample_meta[label_column].astype(str)
    valid = y.notna() & (y != "")
    sample_meta = sample_meta.loc[valid].copy()
    X = X.loc[valid].copy()
    y = y.loc[valid].copy()
    X = X.select_dtypes(include=[np.number]).copy()

    feature_names = X.columns.astype(str).tolist()
    metrics_frames: list[pd.DataFrame] = []
    prediction_frames: list[pd.DataFrame] = []
    confusion_frames: list[pd.DataFrame] = []
    importance_frames: list[pd.DataFrame] = []
    model_summary_rows: list[dict[str, Any]] = []

    requested_models = deep_get(cfg, "train", "model_names", default=MODEL_ORDER)
    if not isinstance(requested_models, list):
        requested_models = MODEL_ORDER
    model_names = [m for m in requested_models if m in MODEL_ORDER]
    if not model_names:
        model_names = ["logistic_regression", "random_forest"]

    for model_name in model_names:
        model = build_model(model_name, cfg, feature_names)
        if model is None:
            model_summary_rows.append({"model": model_name, "status": "skipped", "reason": "xgboost_unavailable"})
            continue
        cv_preds, cv_metrics, split_info = crossval_predictions(X, y, sample_meta, cfg, model_name)
        lo_preds, lo_metrics = leave_one_dataset_out_predictions(X, y, sample_meta, cfg, model_name)
        combined_preds = pd.concat([cv_preds, lo_preds], ignore_index=True)
        combined_metrics = pd.concat([cv_metrics, lo_metrics], ignore_index=True)
        prediction_frames.append(combined_preds)
        metrics_frames.append(combined_metrics)

        model.fit(X, y)
        importance = compute_feature_importance(model, feature_names, model_name)
        importance["model"] = model_name
        importance_frames.append(importance)

        if not combined_preds.empty:
            cm = pd.crosstab(combined_preds["truth"], combined_preds["pred"]).reset_index().rename(columns={"truth": "truth_label"})
            cm.insert(0, "model", model_name)
            confusion_frames.append(cm)
        model_summary_rows.append(
            {
                "model": model_name,
                "status": "ok",
                "splitter": split_info["splitter"],
                "cv_accuracy_mean": float(cv_metrics["accuracy"].mean()) if not cv_metrics.empty else None,
                "cv_balanced_accuracy_mean": float(cv_metrics["balanced_accuracy"].mean()) if not cv_metrics.empty else None,
                "cv_macro_f1_mean": float(cv_metrics["macro_f1"].mean()) if not cv_metrics.empty else None,
                "lodo_accuracy_mean": float(lo_metrics["accuracy"].mean()) if not lo_metrics.empty else None,
            }
        )

    dataset_only_preds, dataset_only_metrics = dataset_only_negative_control(sample_meta, y, cfg)
    permutation_df = permutation_control(X, y, sample_meta, cfg)

    performance_summary = pd.DataFrame(model_summary_rows)
    validation_summary = pd.concat(metrics_frames + ([dataset_only_metrics] if not dataset_only_metrics.empty else []), ignore_index=True) if metrics_frames or not dataset_only_metrics.empty else pd.DataFrame()
    validation_predictions = pd.concat(prediction_frames + ([dataset_only_preds] if not dataset_only_preds.empty else []), ignore_index=True) if prediction_frames or not dataset_only_preds.empty else pd.DataFrame()
    feature_importance = pd.concat(importance_frames, ignore_index=True) if importance_frames else pd.DataFrame()
    confusion_summary = pd.concat(confusion_frames, ignore_index=True) if confusion_frames else pd.DataFrame()

    performance_summary.to_csv(model_dir / "performance_summary.tsv", sep="\t", index=False)
    validation_summary.to_csv(model_dir / "validation_summary.tsv", sep="\t", index=False)
    validation_predictions.to_csv(model_dir / "validation_predictions.tsv", sep="\t", index=False)
    feature_importance.to_csv(model_dir / "feature_importance.tsv", sep="\t", index=False)
    confusion_summary.to_csv(model_dir / "confusion_matrix.tsv", sep="\t", index=False)
    permutation_df.to_csv(model_dir / "permutation_summary.tsv", sep="\t", index=False)

    manifest = {
        "feature_dir": str(Path(feature_dir)),
        "output_dir": str(model_dir),
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "label_column": label_column,
        "labels": sorted(y.unique().tolist()),
        "dataset_counts": sample_meta["dataset"].astype(str).value_counts().to_dict() if "dataset" in sample_meta.columns else {},
        "model_names": model_names,
        "xgboost_available": bool(XGBClassifier is not None),
    }
    write_json(manifest, model_dir / "model_manifest.json")
    return manifest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train normal-airway site classifiers from exported feature tables")
    parser.add_argument("--feature-dir", required=True, help="Directory containing sample_features_wide.tsv.gz")
    parser.add_argument("--config-yaml", default=None, help="Config YAML path")
    parser.add_argument("--output-dir", default=None, help="Explicit output directory")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    cfg = load_config(args.config_yaml)
    output_dir = Path(args.output_dir) if args.output_dir else ensure_dir(Path(deep_get(cfg, "run", "output_root", default=str(DEFAULT_OUTPUT_ROOT))) / f"normal_airway_ml_train_{timestamp_slug()}")
    manifest = run_training(feature_dir=Path(args.feature_dir), cfg=cfg, output_dir=output_dir)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
