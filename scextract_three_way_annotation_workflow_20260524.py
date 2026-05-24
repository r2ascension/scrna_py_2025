#!/usr/bin/env python3
"""Three-way annotation workflow around scExtract, CellTypist, and ScType.

Notes
-----
The current upstream `yxwucq/scExtract` repository on GitHub main provides the
built-in literature-guided `auto_extract.auto_extract()` function, but does *not*
currently ship `add_celltypist_annotation` or `add_sctype_annotation` helpers on
that branch. Therefore this orchestrator integrates upstream `scExtract` at the
function level where available, and adds companion local function-based branches
for CellTypist and ScType-style annotation.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import re
import sys
import traceback
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import yaml
except ImportError:  # pragma: no cover - handled at runtime
    yaml = None

from scextract_annotation_bridge_20260524 import apply_bridge


WORKSPACE_ROOT = Path("/home/h2048")
DEFAULT_ENV_PYTHON = WORKSPACE_ROOT / "miniconda3" / "envs" / "scextract_env" / "bin" / "python"


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(payload: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def stringify(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def load_config(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if yaml is not None:
        cfg = yaml.safe_load(text)
        return cfg or {}
    try:
        return json.loads(text)
    except Exception as exc:  # pragma: no cover - runtime guard
        raise RuntimeError("PyYAML is not installed and config is not valid JSON.") from exc


def deep_get(obj: dict[str, Any], *keys: str, default: Any = None) -> Any:
    cur: Any = obj
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def now_stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def resolve_run_dir(cfg: dict[str, Any], cli_output_dir: str | None) -> Path:
    if cli_output_dir:
        return ensure_dir(Path(cli_output_dir))
    output_root = Path(deep_get(cfg, "run", "output_root", default=str(WORKSPACE_ROOT / "data" / "py" / datetime.now().strftime("%Y%m%d"))))
    prefix = deep_get(cfg, "run", "run_name_prefix", default="scextract_three_way")
    return ensure_dir(output_root / f"{prefix}_{now_stamp()}")


def resolve_scextract_python(cfg: dict[str, Any]) -> str:
    configured = deep_get(cfg, "scextract", "python_executable", default="auto")
    if configured and str(configured).lower() != "auto":
        return str(configured)
    if DEFAULT_ENV_PYTHON.exists():
        return str(DEFAULT_ENV_PYTHON)
    return sys.executable


def normalize_obs_strings(adata: ad.AnnData) -> None:
    for col in adata.obs.columns:
        dtype_name = getattr(adata.obs[col].dtype, "name", str(adata.obs[col].dtype))
        if dtype_name in {"object", "category", "string"}:
            adata.obs[col] = adata.obs[col].astype("string")


def is_count_like_matrix(matrix, max_sample: int = 2000) -> bool:
    if sp.issparse(matrix):
        data = np.asarray(matrix.data)
    else:
        data = np.asarray(matrix).ravel()
    if data.size == 0:
        return False
    if data.size > max_sample:
        rng = np.random.default_rng(42)
        idx = rng.choice(data.size, size=max_sample, replace=False)
        data = data[idx]
    data = data[np.isfinite(data)]
    if data.size == 0:
        return False
    return bool(np.all(data >= 0) and np.allclose(data, np.round(data), atol=1e-6))


def ensure_counts_layer(adata: ad.AnnData, counts_layer: str) -> dict[str, Any]:
    if counts_layer in adata.layers:
        return {"counts_layer": counts_layer, "created": False, "source": counts_layer}
    if is_count_like_matrix(adata.X):
        adata.layers[counts_layer] = adata.X.copy()
        return {"counts_layer": counts_layer, "created": True, "source": "X"}
    return {"counts_layer": counts_layer, "created": False, "source": None, "warning": "Missing counts layer and X is not count-like"}


def normalize_gene_names(adata: ad.AnnData) -> None:
    if adata.var_names.duplicated().sum() > 0:
        adata.var_names_make_unique()
    symbol_cols = ["gene_symbol", "gene_symbols", "symbol", "feature_name", "gene_name"]
    for col in symbol_cols:
        if col in adata.var.columns:
            vals = adata.var[col].astype(str).str.strip()
            non_empty = (vals != "").mean() if len(vals) else 0.0
            if non_empty > 0.5:
                adata.var["symbol_base"] = vals.values
                return
    adata.var["symbol_base"] = adata.var_names.astype(str).values


def preserve_full_raw_if_missing(adata: ad.AnnData, counts_layer: str) -> None:
    if adata.raw is not None and adata.raw.n_vars >= adata.n_vars:
        if "symbol_base" in adata.var.columns and "symbol_base" not in adata.raw.var.columns:
            common = adata.raw.var_names.intersection(adata.var_names)
            if len(common) > 0:
                adata.raw.var.loc[common, "symbol_base"] = adata.var.loc[common, "symbol_base"].values
        return
    if counts_layer not in adata.layers:
        raise KeyError(f"Required counts layer not found: {counts_layer}")
    adata.raw = ad.AnnData(X=adata.layers[counts_layer], obs=adata.obs.copy(), var=adata.var.copy())


def build_celltypist_input(adata: ad.AnnData, model_features: pd.Index | None, counts_layer: str) -> ad.AnnData:
    import scanpy as sc

    if adata.raw is not None:
        X = adata.raw.X
        var = adata.raw.var.copy()
    else:
        if counts_layer not in adata.layers:
            raise KeyError(f"Required counts layer not found: {counts_layer}")
        X = adata.layers[counts_layer]
        var = adata.var.copy()

    if "symbol_base" in var.columns:
        base_symbols = var["symbol_base"].astype(str).values
    elif "symbol_base" in adata.var.columns:
        mapping = adata.var["symbol_base"].astype(str).to_dict()
        base_symbols = np.array([mapping.get(g, str(g)) for g in var.index], dtype=object)
    else:
        base_symbols = var.index.astype(str).values

    base_index = pd.Index(base_symbols)
    if model_features is not None:
        keep = base_index.isin(model_features) & (~base_index.duplicated(keep="first"))
    else:
        keep = ~base_index.duplicated(keep="first")

    X_ct = X[:, keep].copy()
    var_ct = pd.DataFrame(index=base_index[keep])
    obs_ct = pd.DataFrame(index=adata.obs_names)
    ad_ct = ad.AnnData(X=X_ct, obs=obs_ct, var=var_ct)
    sc.pp.normalize_total(ad_ct, target_sum=1e4)
    sc.pp.log1p(ad_ct)
    return ad_ct


def extract_celltypist_confidence(pred) -> np.ndarray:
    if hasattr(pred, "probability_matrix") and pred.probability_matrix is not None:
        try:
            return np.asarray(pred.probability_matrix.max(axis=1)).ravel()
        except Exception:
            pass
    pred_df = getattr(pred, "predicted_labels", None)
    if isinstance(pred_df, pd.DataFrame):
        for col in ("conf_score", "confidence", "confidence_score", "prob"):
            if col in pred_df.columns:
                return pred_df[col].to_numpy()
    return np.ones(shape=(len(pred.predicted_labels),), dtype=float)


def series_numeric_like(series: pd.Series) -> bool:
    sample = series.astype(str).replace({"nan": None, "<NA>": None}).dropna().head(20)
    if sample.empty:
        return False
    return bool(sample.map(lambda x: re.fullmatch(r"-?\d+(\.\d+)?", x) is not None).all())


def infer_scextract_label_column(adata: ad.AnnData) -> str:
    candidates = []
    for col in ("leiden", "louvain"):
        if col in adata.obs.columns and not series_numeric_like(adata.obs[col]):
            candidates.append(col)
    if candidates:
        return candidates[0]
    for col in adata.obs.columns:
        if col.endswith("_rough") and not series_numeric_like(adata.obs[col]):
            return col
    raise RuntimeError("Unable to infer scExtract annotation label column from output AnnData")


def save_status(branch_dir: Path, branch_name: str, status: str, message: str = "", **extra) -> dict[str, Any]:
    payload = {"branch": branch_name, "status": status, "message": message, **extra}
    write_json(payload, branch_dir / f"{branch_name}_status.json")
    return payload


def ensure_batch_key(adata: ad.AnnData, batch_key: str, auto_create: bool) -> dict[str, Any]:
    if batch_key in adata.obs.columns:
        return {"batch_key": batch_key, "created": False}
    if auto_create:
        adata.obs[batch_key] = pd.Categorical(["batch_0"] * adata.n_obs)
        return {"batch_key": batch_key, "created": True}
    return {"batch_key": batch_key, "created": False, "warning": f"Missing batch key: {batch_key}"}


def ensure_cluster_key(adata: ad.AnnData, cluster_key: str, auto_compute: bool) -> dict[str, Any]:
    if cluster_key in adata.obs.columns:
        return {"cluster_key": cluster_key, "created": False}
    if not auto_compute:
        raise KeyError(f"Missing required cluster key: {cluster_key}")

    import scanpy as sc

    work_key = "X_pca"
    if work_key not in adata.obsm:
        n_comps = max(2, min(50, adata.n_obs - 1, adata.n_vars - 1))
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=42)
    n_neighbors = max(2, min(15, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca", random_state=42)
    if cluster_key == "louvain":
        sc.tl.louvain(adata, key_added=cluster_key, random_state=42)
    else:
        sc.tl.leiden(adata, key_added=cluster_key, random_state=42)
    return {"cluster_key": cluster_key, "created": True}


def prepare_preflight_adata(input_h5ad: Path, cfg: dict[str, Any], runtime_dir: Path) -> tuple[Path, dict[str, Any]]:
    adata = ad.read_h5ad(input_h5ad)
    normalize_obs_strings(adata)
    normalize_gene_names(adata)

    counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
    batch_key = deep_get(cfg, "data_contract", "batch_key", default="Batch")
    cluster_key = deep_get(cfg, "data_contract", "cluster_key", default="leiden")

    counts_info = ensure_counts_layer(adata, counts_layer)
    batch_info = ensure_batch_key(adata, batch_key, bool(deep_get(cfg, "data_contract", "auto_create_batch_key", default=True)))

    preserve_info: dict[str, Any]
    try:
        preserve_full_raw_if_missing(adata, counts_layer)
        preserve_info = {"raw_ready": True, "raw_n_vars": int(adata.raw.n_vars) if adata.raw is not None else 0}
    except Exception as exc:
        preserve_info = {"raw_ready": False, "warning": str(exc)}

    cluster_info = {"cluster_key": cluster_key, "created": False, "prepared": cluster_key in adata.obs.columns}
    preflight_path = runtime_dir / "preflight_input.h5ad"
    adata.write_h5ad(preflight_path, compression="gzip")
    return preflight_path, {
        "counts": counts_info,
        "batch": batch_info,
        "raw": preserve_info,
        "cluster": cluster_info,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "preflight_path": str(preflight_path),
    }


def ensure_scextract_repo_ready(repo_dir: Path) -> None:
    main_py = repo_dir / "main.py"
    config_py = repo_dir / "auto_extract" / "config.py"
    if not main_py.exists():
        raise FileNotFoundError(f"scExtract main.py not found: {main_py}")
    if not config_py.exists():
        raise FileNotFoundError(f"scExtract auto_extract/config.py not found: {config_py}")


def load_scextract_auto_extract(repo_dir: Path):
    repo_dir_str = str(repo_dir)
    if repo_dir_str not in sys.path:
        sys.path.insert(0, repo_dir_str)
    module = importlib.import_module("auto_extract.auto_extract")
    auto_extract_func = getattr(module, "auto_extract", None)
    if auto_extract_func is None:
        raise AttributeError("Could not resolve auto_extract.auto_extract from scExtract repository")
    return auto_extract_func


def run_scextract_branch(preflight_h5ad: Path, pdf_path: Path, cfg: dict[str, Any], branch_dir: Path) -> tuple[dict[str, Any], Path | None]:
    branch_name = "scExtract_builtin"
    output_name = deep_get(cfg, "scextract", "output_name", default="scextract_processed.h5ad")
    repo_dir = Path(deep_get(cfg, "scextract", "repo_dir", default=str(WORKSPACE_ROOT / "tools" / "scExtract")))
    ensure_scextract_repo_ready(repo_dir)

    if not os.getenv(deep_get(cfg, "deepseek", "env_key", default="DEEPSEEK_API_KEY")):
        return save_status(branch_dir, branch_name, "error", "Missing DEEPSEEK_API_KEY in environment"), None

    output_path = branch_dir / output_name
    function_log = branch_dir / "scExtract_builtin_function_call.json"
    write_json(
        {
            "function": "auto_extract.auto_extract",
            "repo_dir": str(repo_dir),
            "adata_path": str(preflight_h5ad),
            "pdf_path": str(pdf_path),
            "output_dir": str(branch_dir),
            "output_name": output_name,
        },
        function_log,
    )

    try:
        auto_extract_func = load_scextract_auto_extract(repo_dir)
        auto_extract_func(
            adata_path=str(preflight_h5ad),
            pdf_path=str(pdf_path),
            output_dir=str(branch_dir),
            output_name=output_name,
        )
    except Exception as exc:
        tb_path = branch_dir / "scExtract_builtin_traceback.txt"
        tb_path.write_text(traceback.format_exc(), encoding="utf-8")
        return save_status(
            branch_dir,
            branch_name,
            "error",
            f"scExtract function call failed: {exc}",
            traceback_file=str(tb_path),
            function="auto_extract.auto_extract",
            function_log=str(function_log),
        ), None

    log_path = branch_dir / "auto_extract.log"
    if not output_path.exists():
        return save_status(branch_dir, branch_name, "error", "Expected scExtract output h5ad not found", log_path=str(log_path), function_log=str(function_log)), None

    adata_out = ad.read_h5ad(output_path)
    label_col = infer_scextract_label_column(adata_out)
    return save_status(
        branch_dir,
        branch_name,
        "ok",
        f"Built-in scExtract function completed with label column {label_col}",
        output_h5ad=str(output_path),
        label_column=label_col,
        function="auto_extract.auto_extract",
        function_log=str(function_log),
        log_path=str(log_path),
        n_obs=int(adata_out.n_obs),
        n_vars=int(adata_out.n_vars),
    ), output_path


def run_celltypist_branch(preflight_h5ad: Path, cfg: dict[str, Any], branch_dir: Path) -> tuple[dict[str, Any], Path | None]:
    branch_name = "celltypist"
    if not bool(deep_get(cfg, "celltypist", "enabled", default=True)):
        return save_status(branch_dir, branch_name, "skipped", "CellTypist branch disabled in config"), None

    try:
        import celltypist
        from celltypist import models
    except Exception as exc:
        return save_status(branch_dir, branch_name, "error", f"CellTypist import failed: {exc}"), None

    adata = ad.read_h5ad(preflight_h5ad)
    counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
    counts_info = ensure_counts_layer(adata, counts_layer)
    normalize_gene_names(adata)
    preserve_full_raw_if_missing(adata, counts_layer)

    model_path = Path(deep_get(cfg, "celltypist", "model_path", default=""))
    if not model_path.exists():
        return save_status(branch_dir, branch_name, "error", f"CellTypist model not found: {model_path}"), None

    os.environ.setdefault("CELLTYPIST_FOLDER", str(model_path.parent))
    ct_model = models.Model.load(str(model_path))

    model_features = None
    for attr in ("features", "genes", "feature_names"):
        if hasattr(ct_model, attr):
            try:
                model_features = pd.Index([str(x) for x in getattr(ct_model, attr)])
                break
            except Exception:
                model_features = None

    ad_ct = build_celltypist_input(adata, model_features=model_features, counts_layer=counts_layer)
    majority_voting = bool(deep_get(cfg, "celltypist", "majority_voting", default=True))
    pred = celltypist.annotate(ad_ct, model=ct_model, majority_voting=majority_voting)
    pred_df = pred.predicted_labels.copy()
    pred_df.index = ad_ct.obs_names

    adata.obs["celltypist_annotation"] = pred_df["predicted_labels"].astype(str).reindex(adata.obs_names).values
    if majority_voting and "majority_voting" in pred_df.columns:
        adata.obs["celltypist_majority_voting"] = pred_df["majority_voting"].astype(str).reindex(adata.obs_names).values
    adata.obs["celltypist_confidence_scextract"] = extract_celltypist_confidence(pred)
    adata.uns["celltypist_model_path"] = str(model_path)
    adata.uns["celltypist_counts_layer_source"] = counts_info.get("source")

    output_path = branch_dir / "celltypist_annotation.h5ad"
    adata.write_h5ad(output_path, compression="gzip")
    return save_status(
        branch_dir,
        branch_name,
        "ok",
        f"CellTypist completed with model {model_path.name}",
        output_h5ad=str(output_path),
        model_path=str(model_path),
        n_unique_labels=int(pd.Series(adata.obs["celltypist_annotation"].astype(str)).nunique()),
        mean_confidence=float(np.nanmean(pd.to_numeric(adata.obs["celltypist_confidence_scextract"], errors="coerce"))),
    ), output_path


def parse_marker_string(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    tokens = re.split(r"[,;/\s]+", text)
    return [tok.strip().upper() for tok in tokens if tok and tok.strip()]


def build_gene_index(adata: ad.AnnData) -> tuple[dict[str, int], Any]:
    if adata.raw is not None:
        var = adata.raw.var
        matrix = adata.raw.X
    else:
        var = adata.var
        matrix = adata.X
    if "symbol_base" in var.columns:
        genes = var["symbol_base"].astype(str).values
    else:
        genes = var.index.astype(str).values
    mapping: dict[str, int] = {}
    for idx, gene in enumerate(genes):
        key = gene.strip().upper()
        if key and key not in mapping:
            mapping[key] = idx
    return mapping, matrix


def run_sctype_branch(preflight_h5ad: Path, cfg: dict[str, Any], branch_dir: Path) -> tuple[dict[str, Any], Path | None]:
    branch_name = "sctype"
    if not bool(deep_get(cfg, "sctype", "enabled", default=True)):
        return save_status(branch_dir, branch_name, "skipped", "ScType branch disabled in config"), None

    try:
        import scanpy as sc
    except Exception as exc:
        return save_status(branch_dir, branch_name, "error", f"scanpy import failed: {exc}"), None

    adata = ad.read_h5ad(preflight_h5ad)
    counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
    cluster_key = deep_get(cfg, "data_contract", "cluster_key", default="leiden")
    db_path = Path(deep_get(cfg, "sctype", "db_path", default=""))
    tissue = stringify(deep_get(cfg, "sctype", "tissue", default="Lung"))
    unknown_min_score = float(deep_get(cfg, "sctype", "unknown_min_score", default=0.0))

    if not db_path.exists():
        return save_status(branch_dir, branch_name, "error", f"ScType DB not found: {db_path}"), None

    normalize_gene_names(adata)
    counts_info = ensure_counts_layer(adata, counts_layer)
    preserve_full_raw_if_missing(adata, counts_layer)
    try:
        cluster_info = ensure_cluster_key(adata, cluster_key, bool(deep_get(cfg, "data_contract", "auto_compute_cluster_key", default=True)))
    except Exception as exc:
        return save_status(branch_dir, branch_name, "error", f"Could not prepare cluster key {cluster_key}: {exc}"), None

    db = pd.read_excel(db_path)
    required_cols = {"tissueType", "cellName", "geneSymbolmore1", "geneSymbolmore2"}
    if not required_cols.issubset(db.columns):
        return save_status(branch_dir, branch_name, "error", f"Unexpected ScType DB columns: {list(db.columns)}"), None

    if tissue and tissue.lower() not in {"all", "auto"}:
        filtered = db[db["tissueType"].astype(str).str.contains(tissue, case=False, na=False)].copy()
        if filtered.empty:
            filtered = db.copy()
    else:
        filtered = db.copy()

    gene_index, matrix = build_gene_index(adata)
    cluster_labels = adata.obs[cluster_key].astype(str).to_numpy()
    unique_clusters = pd.Index(pd.unique(cluster_labels))
    cluster_means: dict[str, np.ndarray] = {}

    if sp.issparse(matrix):
        matrix = matrix.tocsr()

    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        if sp.issparse(matrix):
            mean_vec = np.asarray(matrix[mask].mean(axis=0)).ravel()
        else:
            mean_vec = np.asarray(matrix[mask]).mean(axis=0)
        cluster_means[str(cluster)] = np.asarray(mean_vec, dtype=float)

    score_rows: list[dict[str, Any]] = []
    for row in filtered.itertuples(index=False):
        pos_markers = parse_marker_string(getattr(row, "geneSymbolmore1", None))
        neg_markers = parse_marker_string(getattr(row, "geneSymbolmore2", None))
        pos_idx = [gene_index[g] for g in pos_markers if g in gene_index]
        neg_idx = [gene_index[g] for g in neg_markers if g in gene_index]
        if not pos_idx:
            continue
        for cluster, mean_vec in cluster_means.items():
            pos_score = float(np.mean(mean_vec[pos_idx])) if pos_idx else 0.0
            neg_score = float(np.mean(mean_vec[neg_idx])) if neg_idx else 0.0
            score_rows.append(
                {
                    "cluster": cluster,
                    "cell_name": getattr(row, "cellName", None),
                    "short_name": getattr(row, "shortName", None),
                    "tissueType": getattr(row, "tissueType", None),
                    "score": pos_score - neg_score,
                    "positive_markers_present": len(pos_idx),
                    "negative_markers_present": len(neg_idx),
                }
            )

    if not score_rows:
        return save_status(branch_dir, branch_name, "error", "No ScType marker genes overlapped the dataset"), None

    score_df = pd.DataFrame(score_rows)
    best_df = score_df.sort_values(["cluster", "score"], ascending=[True, False]).groupby("cluster", as_index=False).first()
    best_df["assigned_label"] = best_df["cell_name"].astype(str)
    best_df.loc[best_df["score"] <= unknown_min_score, "assigned_label"] = deep_get(cfg, "data_contract", "unknown_label", default="Unknown")
    cluster_to_label = dict(zip(best_df["cluster"].astype(str), best_df["assigned_label"].astype(str)))
    cluster_to_score = dict(zip(best_df["cluster"].astype(str), best_df["score"].astype(float)))

    adata.obs["sctype_annotation"] = pd.Series(cluster_labels, index=adata.obs_names).map(cluster_to_label).fillna(deep_get(cfg, "data_contract", "unknown_label", default="Unknown")).astype(str).values
    adata.obs["sctype_score"] = pd.Series(cluster_labels, index=adata.obs_names).map(cluster_to_score).astype(float).values
    adata.uns["sctype_db_path"] = str(db_path)
    adata.uns["sctype_tissue_requested"] = tissue
    adata.uns["sctype_cluster_key"] = cluster_key
    adata.uns["sctype_counts_layer_source"] = counts_info.get("source")

    score_df.to_csv(branch_dir / "sctype_cluster_scores.tsv", sep="\t", index=False)
    best_df.to_csv(branch_dir / "sctype_cluster_best.tsv", sep="\t", index=False)
    output_path = branch_dir / "sctype_annotation.h5ad"
    adata.write_h5ad(output_path, compression="gzip")
    return save_status(
        branch_dir,
        branch_name,
        "ok",
        f"ScType completed using {cluster_key} and tissue={tissue}",
        output_h5ad=str(output_path),
        db_path=str(db_path),
        cluster_key=cluster_key,
        cluster_key_created=bool(cluster_info.get("created", False)),
        n_unique_labels=int(pd.Series(adata.obs["sctype_annotation"].astype(str)).nunique()),
    ), output_path


def merge_branch_outputs(preflight_h5ad: Path, branch_outputs: dict[str, Path | None], cfg: dict[str, Any], merged_dir: Path) -> tuple[Path, dict[str, Any]]:
    merged = ad.read_h5ad(preflight_h5ad)
    aliases = deep_get(cfg, "aliases", default={})
    merge_summary: dict[str, Any] = {}

    sc_path = branch_outputs.get("scExtract_builtin")
    if sc_path and Path(sc_path).exists():
        sc_adata = ad.read_h5ad(sc_path)
        label_col = infer_scextract_label_column(sc_adata)
        alias_col = aliases.get("scExtract", "cell_type_scextract")
        merged.obs[alias_col] = sc_adata.obs[label_col].astype(str).reindex(merged.obs_names).values
        for col in ("Certainty", "Source", "Tissue"):
            if col in sc_adata.obs.columns:
                merged.obs[f"scextract_{col.lower()}"] = sc_adata.obs[col].astype(str).reindex(merged.obs_names).values
        merge_summary["scExtract"] = {"label_column": label_col, "alias_column": alias_col}

    ct_path = branch_outputs.get("celltypist")
    if ct_path and Path(ct_path).exists():
        ct_adata = ad.read_h5ad(ct_path)
        alias_col = aliases.get("celltypist_annotation", "cell_type_celltypist_scextract")
        merged.obs[alias_col] = ct_adata.obs["celltypist_annotation"].astype(str).reindex(merged.obs_names).values
        if "celltypist_majority_voting" in ct_adata.obs.columns:
            merged.obs["celltypist_majority_voting_scextract"] = ct_adata.obs["celltypist_majority_voting"].astype(str).reindex(merged.obs_names).values
        if "celltypist_confidence_scextract" in ct_adata.obs.columns:
            merged.obs["celltypist_confidence_scextract"] = pd.to_numeric(ct_adata.obs["celltypist_confidence_scextract"], errors="coerce").reindex(merged.obs_names).values
        merge_summary["celltypist"] = {"alias_column": alias_col}

    sctype_path = branch_outputs.get("sctype")
    if sctype_path and Path(sctype_path).exists():
        sc_type_adata = ad.read_h5ad(sctype_path)
        alias_col = aliases.get("sctype_annotation", "cell_type_sctype_scextract")
        merged.obs[alias_col] = sc_type_adata.obs["sctype_annotation"].astype(str).reindex(merged.obs_names).values
        if "sctype_score" in sc_type_adata.obs.columns:
            merged.obs["sctype_score"] = pd.to_numeric(sc_type_adata.obs["sctype_score"], errors="coerce").reindex(merged.obs_names).values
        merge_summary["sctype"] = {"alias_column": alias_col}

    bridge_summary = apply_bridge(
        merged,
        selected_label_source=deep_get(cfg, "run", "selected_label_source", default="compare_only"),
        unknown_label=deep_get(cfg, "data_contract", "unknown_label", default="Unknown"),
        min_cells_per_type=int(deep_get(cfg, "data_contract", "min_cells_per_type", default=10)),
        confidence_column="celltypist_confidence_scextract",
        min_confidence=deep_get(cfg, "celltypist", "min_confidence", default=None),
    )
    merge_summary["bridge"] = bridge_summary

    output_path = merged_dir / "scextract_three_way_merged.h5ad"
    merged.write_h5ad(output_path, compression="gzip")
    return output_path, merge_summary


def compute_pairwise_agreement(adata: ad.AnnData, columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for left, right in combinations(columns, 2):
        s1 = adata.obs[left].astype("string") if left in adata.obs.columns else pd.Series(pd.NA, index=adata.obs_names)
        s2 = adata.obs[right].astype("string") if right in adata.obs.columns else pd.Series(pd.NA, index=adata.obs_names)
        valid = s1.notna() & s2.notna() & (s1.astype(str) != "") & (s2.astype(str) != "")
        if valid.any():
            agreement = float((s1[valid].astype(str) == s2[valid].astype(str)).mean())
            compared = int(valid.sum())
        else:
            agreement = np.nan
            compared = 0
        rows.append({"left": left, "right": right, "n_compared": compared, "exact_match_rate": agreement})
    return pd.DataFrame(rows)


def compute_truth_benchmark(adata: ad.AnnData, columns: list[str], true_group_key: str) -> pd.DataFrame:
    if true_group_key not in adata.obs.columns:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    truth = adata.obs[true_group_key].astype(str)
    try:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    except Exception:
        adjusted_rand_score = None
        normalized_mutual_info_score = None

    for col in columns:
        if col not in adata.obs.columns:
            continue
        pred = adata.obs[col].astype(str)
        valid = truth.notna() & pred.notna() & (truth != "") & (pred != "")
        if not valid.any():
            continue
        truth_valid = truth[valid].astype(str)
        pred_valid = pred[valid].astype(str)
        exact = float((truth_valid.str.strip().str.lower() == pred_valid.str.strip().str.lower()).mean())
        row: dict[str, Any] = {"prediction_column": col, "n_compared": int(valid.sum()), "exact_match_rate": exact}
        if adjusted_rand_score is not None and normalized_mutual_info_score is not None:
            row["ari"] = float(adjusted_rand_score(truth_valid, pred_valid))
            row["nmi"] = float(normalized_mutual_info_score(truth_valid, pred_valid))
        rows.append(row)
    return pd.DataFrame(rows)


def build_summary_outputs(merged_h5ad: Path, cfg: dict[str, Any], summary_dir: Path) -> dict[str, Any]:
    adata = ad.read_h5ad(merged_h5ad)
    columns = [
        deep_get(cfg, "aliases", "scExtract", default="cell_type_scextract"),
        deep_get(cfg, "aliases", "celltypist_annotation", default="cell_type_celltypist_scextract"),
        deep_get(cfg, "aliases", "sctype_annotation", default="cell_type_sctype_scextract"),
    ]
    columns = [col for col in columns if col in adata.obs.columns]

    branch_rows = []
    for col in columns:
        series = adata.obs[col].astype("string")
        missing = int(series.isna().sum() + (series.astype(str) == "").sum())
        branch_rows.append(
            {
                "column": col,
                "n_cells": int(adata.n_obs),
                "n_missing": missing,
                "n_unique_labels": int(series.dropna().astype(str).replace("", pd.NA).dropna().nunique()),
            }
        )
        counts = series.astype(str).value_counts().rename_axis("label").reset_index(name="count")
        counts.to_csv(summary_dir / f"{col}_counts.tsv", sep="\t", index=False)

    pairwise = compute_pairwise_agreement(adata, columns)
    pairwise.to_csv(summary_dir / "pairwise_agreement.tsv", sep="\t", index=False)

    true_group_key = deep_get(cfg, "data_contract", "true_group_key", default=None)
    truth_df = compute_truth_benchmark(adata, columns, true_group_key) if true_group_key else pd.DataFrame()
    if not truth_df.empty:
        truth_df.to_csv(summary_dir / "truth_benchmark.tsv", sep="\t", index=False)

    branch_df = pd.DataFrame(branch_rows)
    branch_df.to_csv(summary_dir / "branch_summary.tsv", sep="\t", index=False)
    summary = {
        "columns": columns,
        "branch_summary_tsv": str(summary_dir / "branch_summary.tsv"),
        "pairwise_agreement_tsv": str(summary_dir / "pairwise_agreement.tsv"),
        "truth_benchmark_tsv": str(summary_dir / "truth_benchmark.tsv") if not truth_df.empty else None,
    }
    write_json(summary, summary_dir / "summary_manifest.json")
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run three-way scExtract/CellTypist/ScType annotation workflow.")
    parser.add_argument("--adata", required=True, help="Input AnnData h5ad path")
    parser.add_argument("--pdf", required=True, help="Source article PDF/txt path")
    parser.add_argument("--config-yaml", required=True, help="Workflow YAML config")
    parser.add_argument("--output-dir", default=None, help="Explicit run directory (optional)")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    config_path = Path(args.config_yaml)
    cfg = load_config(config_path)
    run_dir = resolve_run_dir(cfg, args.output_dir)
    runtime_dir = ensure_dir(run_dir / "runtime")
    summary_dir = ensure_dir(run_dir / "summary")
    merged_dir = ensure_dir(run_dir / "merged")

    branch_dirs = {
        "scExtract_builtin": ensure_dir(run_dir / "branches" / "scExtract_builtin"),
        "celltypist": ensure_dir(run_dir / "branches" / "celltypist"),
        "sctype": ensure_dir(run_dir / "branches" / "sctype"),
    }

    manifest: dict[str, Any] = {
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "input_h5ad": str(Path(args.adata).resolve()),
        "input_pdf": str(Path(args.pdf).resolve()),
        "config_yaml": str(config_path.resolve()),
        "run_dir": str(run_dir.resolve()),
        "scextract_python": resolve_scextract_python(cfg),
        "deepseek_env_key": deep_get(cfg, "deepseek", "env_key", default="DEEPSEEK_API_KEY"),
        "deepseek_env_present": bool(os.getenv(deep_get(cfg, "deepseek", "env_key", default="DEEPSEEK_API_KEY"))),
    }
    write_json(manifest, run_dir / "run_manifest.json")

    preflight_h5ad, preflight_summary = prepare_preflight_adata(Path(args.adata), cfg, runtime_dir)
    manifest["preflight"] = preflight_summary
    write_json(manifest, run_dir / "run_manifest.json")

    statuses: list[dict[str, Any]] = []
    branch_outputs: dict[str, Path | None] = {}
    branch_runners = [
        ("scExtract_builtin", lambda: run_scextract_branch(preflight_h5ad, Path(args.pdf), cfg, branch_dirs["scExtract_builtin"])),
        ("celltypist", lambda: run_celltypist_branch(preflight_h5ad, cfg, branch_dirs["celltypist"])),
        ("sctype", lambda: run_sctype_branch(preflight_h5ad, cfg, branch_dirs["sctype"])),
    ]

    for branch_name, runner in branch_runners:
        try:
            status, output_path = runner()
        except Exception as exc:  # pragma: no cover - runtime guard
            tb_path = branch_dirs[branch_name] / f"{branch_name}_traceback.txt"
            tb_path.write_text(traceback.format_exc(), encoding="utf-8")
            status = save_status(branch_dirs[branch_name], branch_name, "error", str(exc), traceback_file=str(tb_path))
            output_path = None
        statuses.append(status)
        branch_outputs[branch_name] = output_path

    merged_h5ad, merge_summary = merge_branch_outputs(preflight_h5ad, branch_outputs, cfg, merged_dir)
    summary_manifest = build_summary_outputs(merged_h5ad, cfg, summary_dir)

    manifest.update(
        {
            "finished_at": datetime.now().isoformat(timespec="seconds"),
            "statuses": statuses,
            "branch_outputs": {k: (str(v) if v is not None else None) for k, v in branch_outputs.items()},
            "merged_h5ad": str(merged_h5ad),
            "merge_summary": merge_summary,
            "summary_manifest": summary_manifest,
        }
    )
    write_json(manifest, run_dir / "run_manifest.json")

    pd.DataFrame(statuses).to_csv(run_dir / "branch_status.tsv", sep="\t", index=False)
    print(json.dumps({
        "run_dir": str(run_dir),
        "merged_h5ad": str(merged_h5ad),
        "status_tsv": str(run_dir / "branch_status.tsv"),
        "summary_manifest": str(summary_dir / "summary_manifest.json"),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
