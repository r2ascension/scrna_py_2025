#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run cNMF independently for every L3 cell type.

This script reuses `/home/h2048/script/py/cnmf_helper_20260419_v1_1.py`, so the
per-celltype runs keep the same filtering contract and output schema as the
existing lineage-level `cnmf_full` runs.

Default output layout:
    <RUN_ROOT>/<lineage>/cnmf_by_celltype/<safe_L3>/cnmf_full/

Key env vars:
  CNMF_L3_LINEAGES=bcell,tnk
  CNMF_L3_CELLTYPES=Memory_B,Plasma_IgA
  CNMF_L3_LIMIT=1
  CNMF_L3_FORCE=1
  CNMF_L3_PREFLIGHT_ONLY=1
  CNMF_L3_N_ITER=100
  CNMF_L3_NUM_HVG=3000
"""

from __future__ import annotations

import gc
import importlib.util
import json
import os
import re
import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import anndata as ad
import numpy as np
import pandas as pd


RUN_ROOT = Path(os.environ.get("CNMF_L3_RUN_ROOT", "/home/h2048/output/program_full_parallel_methods_20260507"))
HELPER_PATH = Path(os.environ.get("CNMF_L3_HELPER_PY", "/home/h2048/script/py/cnmf_helper_20260419_v1_1.py"))
OUTPUT_SUBDIR = os.environ.get("CNMF_L3_OUTPUT_SUBDIR", "cnmf_by_celltype")
METHOD_DIRNAME = "cnmf_full"

LOW_INFO_PATTERN = re.compile(
    r"^(?:MT-|RPS[0-9A-Z]*$|RPL[0-9A-Z]*$|MRPS[0-9A-Z]*$|MRPL[0-9A-Z]*$|"
    r"(?:RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$|(?:AC|AL|AP|BX|Z)[0-9]+[.]|"
    r"RP[0-9]+-|CTD-|CTB-|CTC-|LOC[0-9]+)|-OT[0-9]+$",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True)
class LineageConfig:
    lineage: str
    h5ad_path: str
    celltype_col: str = "cell_type_L3"
    batch_col: str = "dataset"
    sample_col: str = "sample"
    tissue_col: str = "tissue"


LINEAGE_CONFIGS: List[LineageConfig] = [
    LineageConfig("bcell", "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/bcell_tissue_comparison_final.h5ad"),
    LineageConfig("stromal_smc", "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_1_neuronlike_20260414/stromal_smc_tissue_comparison_final.h5ad"),
    LineageConfig("stromal_fibroblast", "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_fibroblast_tissue_comparison_final.h5ad"),
    LineageConfig("stromal_endothelial", "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_1_rm_choir_20260414/stromal_endothelial_tissue_comparison_final.h5ad"),
    LineageConfig("tnk", "/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0/tnk_tissue_comparison_final.h5ad"),
    LineageConfig("myeloid", "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416/myeloid_tissue_comparison_final.h5ad"),
    LineageConfig("epithelial", "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/epithelial_tissue_comparison_final.h5ad"),
]


def split_env(name: str) -> List[str]:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return []
    return [x.strip() for x in raw.split(",") if x.strip()]


def truthy(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name, "1" if default else "0").strip().lower()
    return raw in {"1", "true", "yes", "y"}


def int_env(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


def safe_name(x: Any, max_len: int = 180) -> str:
    out = str(x).strip()
    for ch in ["/", "\\", " ", "|", ":", ";", ",", "\t", "'", '"']:
        out = out.replace(ch, "_")
    out = re.sub(r"_+", "_", out).strip("_")
    return (out or "NA")[:max_len]


def pick_primary_value(values: Iterable[Any]) -> Optional[str]:
    series = pd.Series(list(values), dtype="object").astype(str).str.strip()
    series = series[~series.isin(["", "NA", "nan", "None"])]
    if series.empty:
        return None
    return str(series.value_counts().index[0])


def load_helper(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"cNMF helper not found: {path}")
    spec = importlib.util.spec_from_file_location("pa_cnmf_helper_l3", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import cNMF helper: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_json(obj: Dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def selected_configs() -> List[LineageConfig]:
    lineages = set(split_env("CNMF_L3_LINEAGES"))
    if not lineages:
        return LINEAGE_CONFIGS
    return [cfg for cfg in LINEAGE_CONFIGS if cfg.lineage in lineages]


def build_cnmf_config(lineage: str, celltype: str) -> Dict[str, Any]:
    cfg: Dict[str, Any] = {
        "n_workers": int_env("CNMF_L3_N_WORKERS", 1),
        "gene_exclusion_config": {"lineage_context": f"{lineage} | {celltype}"},
    }
    if "CNMF_L3_N_ITER" in os.environ:
        cfg["n_iter"] = int_env("CNMF_L3_N_ITER", 100)
    if "CNMF_L3_NUM_HVG" in os.environ:
        cfg["num_hvg"] = int_env("CNMF_L3_NUM_HVG", 3000)
    if "CNMF_L3_DENSITY_THRESHOLD" in os.environ:
        try:
            cfg["density_threshold"] = float(os.environ["CNMF_L3_DENSITY_THRESHOLD"])
        except Exception:
            pass
    return cfg


def celltype_method_output_dir(lineage: str, safe_celltype: str, method_dirname: str = METHOD_DIRNAME) -> Path:
    return RUN_ROOT / lineage / OUTPUT_SUBDIR / safe_celltype / method_dirname


def legacy_cnmf_output_dir(lineage: str, safe_celltype: str) -> Path:
    return RUN_ROOT / lineage / OUTPUT_SUBDIR / f"cnmf_{safe_celltype}" / METHOD_DIRNAME


def is_legacy_cnmf_dir(celltype_root: Path, status: Dict[str, Any]) -> bool:
    if str(status.get("path_layout", "")).strip() == "celltype_method":
        return False
    return celltype_root.name.startswith("cnmf_")


def iter_cnmf_l3_status_records(root: Path = RUN_ROOT) -> List[tuple[Path, Dict[str, Any]]]:
    candidates = set(root.glob(f"*/{OUTPUT_SUBDIR}/*/{METHOD_DIRNAME}/cnmf_l3_status.json"))
    candidates.update(root.glob(f"*/{OUTPUT_SUBDIR}/cnmf_*/{METHOD_DIRNAME}/cnmf_l3_status.json"))
    ranked: List[tuple[tuple[str, str], bool, bool, Path, Dict[str, Any]]] = []
    for path in sorted(candidates):
        status = read_json(path)
        if not status:
            continue
        lineage = str(status.get("lineage") or path.parts[-5])
        celltype_root = path.parent.parent
        safe_celltype = str(status.get("safe_celltype") or "").strip()
        if not safe_celltype:
            safe_celltype = celltype_root.name[5:] if celltype_root.name.startswith("cnmf_") else celltype_root.name
        celltype_l3 = str(status.get("celltype_l3") or "").strip()
        key = (lineage, celltype_l3 or safe_celltype)
        ranked.append(
            (
                key,
                str(status.get("status") or "").strip() != "ok",
                is_legacy_cnmf_dir(celltype_root, status),
                path,
                status,
            )
        )
    ranked.sort(key=lambda x: (x[0][0], x[0][1], x[1], x[2], str(x[3])))
    selected: List[tuple[Path, Dict[str, Any]]] = []
    seen: set[tuple[str, str]] = set()
    for key, _, _, path, status in ranked:
        if key in seen:
            continue
        seen.add(key)
        selected.append((path, status))
    return selected


def write_cell_metadata(adata_sub: ad.AnnData, cfg: LineageConfig, out_dir: Path) -> None:
    keep_cols = [cfg.celltype_col, cfg.sample_col, cfg.tissue_col, cfg.batch_col, "condition", "cell_type_L2", "cell_type_L1"]
    keep_cols = [c for c in keep_cols if c in adata_sub.obs.columns]
    meta = adata_sub.obs[keep_cols].copy()
    meta.insert(0, "cell", adata_sub.obs_names.astype(str))
    meta.to_csv(out_dir / "cell_metadata_for_l3_cnmf.csv", index=False)


def filter_low_info_export_tables(out_dir: Path) -> Dict[str, Any]:
    """Apply the same extra low-information filter used by the unit LLM launcher."""
    gep_dir = out_dir / "gep_gene_tables"
    manifest: Dict[str, Any] = {"filtered_tables": [], "low_info_removed_rows": 0}
    if not gep_dir.exists():
        return manifest
    for score_path in sorted(gep_dir.glob("gep_gene_scores_k*.tsv")):
        try:
            df = pd.read_csv(score_path, sep="\t")
        except Exception:
            continue
        if "gene" not in df.columns or "gep" not in df.columns:
            continue
        before = len(df)
        keep = ~df["gene"].astype(str).str.contains(LOW_INFO_PATTERN, regex=True, na=False)
        df2 = df.loc[keep].copy()
        removed = before - len(df2)
        if removed > 0:
            if "rank" in df2.columns:
                df2["source_rank"] = df2["rank"]
                df2["rank"] = df2.groupby("gep").cumcount() + 1
            df2.to_csv(score_path, sep="\t", index=False)
            k_match = re.search(r"k(\d+)", score_path.name)
            if k_match:
                wide_path = gep_dir / f"gep_top_genes_k{k_match.group(1)}.csv"
                wide = {
                    gep: sub.sort_values("rank")["gene"].astype(str).tolist()
                    for gep, sub in df2.groupby("gep", sort=False)
                }
                if wide:
                    max_len = max(len(v) for v in wide.values())
                    wide_df = pd.DataFrame({k: v + [""] * (max_len - len(v)) for k, v in wide.items()})
                    wide_df.to_csv(wide_path, index=True)
            manifest["filtered_tables"].append(str(score_path))
            manifest["low_info_removed_rows"] += int(removed)
    (out_dir / "l3_cnmf_low_information_filter_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return manifest


def summarize_cnmf_output(out_dir: Path) -> Dict[str, Any]:
    gep_dir = out_dir / "gep_gene_tables"
    score_files = sorted(gep_dir.glob("gep_gene_scores_k*.tsv")) if gep_dir.exists() else []
    unit_count = 0
    low_info_rows = 0
    for path in score_files:
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if "gep" in df.columns:
            unit_count += int(df["gep"].astype(str).nunique())
        if "gene" in df.columns:
            low_info_rows += int(df["gene"].astype(str).str.contains(LOW_INFO_PATTERN, regex=True, na=False).sum())
    rec = read_json(out_dir / "k_selection_recommendation.json")
    return {
        "n_score_files": len(score_files),
        "n_k_gep_units": unit_count,
        "low_info_rows_in_exported_scores": low_info_rows,
        "recommended_k": rec.get("recommended_k"),
        "score_files": [str(p) for p in score_files],
    }


def collect_audit(root: Path = RUN_ROOT) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for path, status in iter_cnmf_l3_status_records(root):
        rows.append(
            {
                "lineage": status.get("lineage"),
                "celltype_l1": status.get("celltype_l1"),
                "celltype_l2": status.get("celltype_l2"),
                "celltype_l3": status.get("celltype_l3"),
                "safe_celltype": status.get("safe_celltype"),
                "status": status.get("status"),
                "n_cells": status.get("n_cells"),
                "n_samples": status.get("n_samples"),
                "n_tissues": status.get("n_tissues"),
                "recommended_k": status.get("recommended_k"),
                "n_score_files": status.get("n_score_files"),
                "n_k_gep_units": status.get("n_k_gep_units"),
                "low_info_rows_in_exported_scores": status.get("low_info_rows_in_exported_scores"),
                "elapsed_min": status.get("elapsed_min"),
                "output_dir": status.get("output_dir"),
                "error": status.get("error"),
                "status_json": str(path),
            }
        )
    return pd.DataFrame(rows)


def write_global_audit() -> Path:
    audit = collect_audit(RUN_ROOT)
    out = RUN_ROOT / "cnmf_l3_celltype_audit_20260519.tsv"
    if not audit.empty:
        audit = audit.sort_values(["lineage", "celltype_l3"])
    audit.to_csv(out, sep="\t", index=False)
    print(f"[AUDIT] {out} rows={len(audit)}", flush=True)
    return out


def run_one(helper, adata_all: ad.AnnData, cfg: LineageConfig, celltype: str, stats: Dict[str, Any]) -> Dict[str, Any]:
    safe_ct = safe_name(celltype)
    out_dir = celltype_method_output_dir(cfg.lineage, safe_ct)
    out_dir.mkdir(parents=True, exist_ok=True)
    status_path = out_dir / "cnmf_l3_status.json"
    prev = read_json(status_path)
    if not truthy("CNMF_L3_FORCE", False) and prev.get("status") == "ok" and (out_dir / "run_summary.json").exists():
        print(f"[SKIP] {cfg.lineage} / {celltype} existing ok", flush=True)
        return prev

    status: Dict[str, Any] = {
        "status": "planned" if truthy("CNMF_L3_PREFLIGHT_ONLY", False) else "started",
        "lineage": cfg.lineage,
        "celltype_l3": celltype,
        "safe_celltype": safe_ct,
        "path_layout": "celltype_method",
        "n_cells": int(stats["n_cells"]),
        "n_samples": int(stats["n_samples"]),
        "n_tissues": int(stats["n_tissues"]),
        "output_dir": str(out_dir),
        "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "error": None,
    }
    write_json(status, status_path)
    if truthy("CNMF_L3_PREFLIGHT_ONLY", False):
        return status

    t0 = time.time()
    try:
        mask = adata_all.obs[cfg.celltype_col].astype(str).to_numpy() == str(celltype)
        adata_sub = adata_all[mask, :].copy()
        if "cell_type_L1" in adata_sub.obs.columns:
            primary_l1 = pick_primary_value(adata_sub.obs["cell_type_L1"])
            if primary_l1:
                status["celltype_l1"] = primary_l1
        if "cell_type_L2" in adata_sub.obs.columns:
            primary_l2 = pick_primary_value(adata_sub.obs["cell_type_L2"])
            if primary_l2:
                status["celltype_l2"] = primary_l2
        write_json(status, status_path)
        write_cell_metadata(adata_sub, cfg, out_dir)
        run_name = safe_name(f"{cfg.lineage}_{safe_ct}_cnmf_l3_20260519")
        print(f"[RUN] {cfg.lineage} / {celltype} cells={adata_sub.n_obs} genes={adata_sub.n_vars} out={out_dir}", flush=True)
        result = helper.run_cnmf_full(
            adata=adata_sub,
            output_dir=out_dir,
            run_name=run_name,
            k_range=None,
            celltype_col=cfg.celltype_col,
            batch_col=cfg.batch_col if cfg.batch_col in adata_sub.obs.columns else None,
            use_batch_hvg=not truthy("CNMF_L3_DISABLE_BATCH_HVG", False),
            cnmf_config=build_cnmf_config(cfg.lineage, celltype),
            viz_config={},
        )
        extra_filter = filter_low_info_export_tables(out_dir)
        summary = summarize_cnmf_output(out_dir)
        status.update(
            {
                "status": "ok" if bool(result.get("success")) else "error",
                "elapsed_min": round((time.time() - t0) / 60.0, 3),
                "ended_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "run_name": run_name,
                "recommended_k": summary.get("recommended_k"),
                "n_score_files": summary.get("n_score_files"),
                "n_k_gep_units": summary.get("n_k_gep_units"),
                "low_info_rows_in_exported_scores": summary.get("low_info_rows_in_exported_scores"),
                "low_info_removed_rows_post_export": extra_filter.get("low_info_removed_rows", 0),
                "result_success": bool(result.get("success")),
                "error": None if bool(result.get("success")) else result.get("error", "run_cnmf_full returned success=False"),
            }
        )
        if status["status"] == "error":
            print(f"[ERROR] {cfg.lineage} / {celltype}: {status['error']}", flush=True)
        else:
            print(
                f"[OK] {cfg.lineage} / {celltype} units={status['n_k_gep_units']} "
                f"recommended_k={status['recommended_k']} low_info={status['low_info_rows_in_exported_scores']}",
                flush=True,
            )
        del adata_sub, result
        gc.collect()
    except Exception as exc:
        status.update(
            {
                "status": "error",
                "elapsed_min": round((time.time() - t0) / 60.0, 3),
                "ended_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "error": str(exc),
                "traceback": traceback.format_exc(),
            }
        )
        print(f"[ERROR] {cfg.lineage} / {celltype}: {exc}", flush=True)
    write_json(status, status_path)
    return status


def main() -> int:
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    helper = load_helper(HELPER_PATH)
    selected_celltypes = set(split_env("CNMF_L3_CELLTYPES"))
    min_cells = int_env("CNMF_L3_MIN_CELLS", 50)
    min_samples = int_env("CNMF_L3_MIN_SAMPLES", 1)
    limit = int_env("CNMF_L3_LIMIT", 0)
    selected_runs = 0
    inventory_rows: List[Dict[str, Any]] = []
    print(f"[BOOT] RUN_ROOT={RUN_ROOT}", flush=True)
    print(f"[BOOT] helper={HELPER_PATH}", flush=True)
    print(f"[BOOT] min_cells={min_cells} min_samples={min_samples} limit={limit} preflight={truthy('CNMF_L3_PREFLIGHT_ONLY', False)} force={truthy('CNMF_L3_FORCE', False)}", flush=True)

    for cfg in selected_configs():
        h5ad_path = Path(cfg.h5ad_path)
        if not h5ad_path.exists():
            raise FileNotFoundError(f"Missing h5ad for {cfg.lineage}: {h5ad_path}")
        print(f"\n[LOAD] {cfg.lineage} | {h5ad_path}", flush=True)
        adata_all = ad.read_h5ad(h5ad_path)
        for col in [cfg.celltype_col, cfg.sample_col]:
            if col not in adata_all.obs.columns:
                raise KeyError(f"{cfg.lineage} missing obs column {col}")
        obs = adata_all.obs.copy()
        obs["__celltype"] = obs[cfg.celltype_col].astype(str)
        obs["__sample"] = obs[cfg.sample_col].astype(str) if cfg.sample_col in obs.columns else "NA"
        obs["__tissue"] = obs[cfg.tissue_col].astype(str) if cfg.tissue_col in obs.columns else "NA"
        stats_df = (
            obs.loc[obs["__celltype"].notna() & (obs["__celltype"] != "")]
            .groupby("__celltype", observed=False)
            .agg(n_cells=("__celltype", "size"), n_samples=("__sample", pd.Series.nunique), n_tissues=("__tissue", pd.Series.nunique))
            .reset_index()
            .rename(columns={"__celltype": "celltype_l3"})
            .sort_values(["n_cells", "celltype_l3"], ascending=[False, True])
        )
        stats_df.insert(0, "lineage", cfg.lineage)
        stats_df["h5ad_path"] = str(h5ad_path)
        stats_df["eligible"] = (stats_df["n_cells"] >= min_cells) & (stats_df["n_samples"] >= min_samples)
        if selected_celltypes:
            stats_df["selected_by_env"] = stats_df["celltype_l3"].isin(selected_celltypes)
        else:
            stats_df["selected_by_env"] = True
        lineage_plan_dir = RUN_ROOT / cfg.lineage / OUTPUT_SUBDIR
        lineage_plan_dir.mkdir(parents=True, exist_ok=True)
        stats_df.to_csv(lineage_plan_dir / f"{cfg.lineage}_cnmf_l3_plan.tsv", sep="\t", index=False)
        inventory_rows.extend(stats_df.to_dict(orient="records"))
        run_df = stats_df.loc[stats_df["eligible"] & stats_df["selected_by_env"]].copy()
        print(f"[PLAN] {cfg.lineage} total={len(stats_df)} eligible={int(stats_df['eligible'].sum())} selected={len(run_df)}", flush=True)
        for rec in run_df.to_dict(orient="records"):
            if limit > 0 and selected_runs >= limit:
                break
            run_one(helper, adata_all, cfg, str(rec["celltype_l3"]), rec)
            selected_runs += 1
            write_global_audit()
        del adata_all
        gc.collect()
        if limit > 0 and selected_runs >= limit:
            break

    inventory = pd.DataFrame(inventory_rows)
    inventory_out = RUN_ROOT / "cnmf_l3_celltype_inventory_20260519.tsv"
    inventory.to_csv(inventory_out, sep="\t", index=False)
    audit_out = write_global_audit()
    print(f"[DONE] selected_runs={selected_runs} inventory={inventory_out} audit={audit_out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())