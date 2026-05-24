#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit L3 cell-type cNMF outputs and optional per-unit LLM outputs."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd


RUN_ROOT = Path(os.environ.get("CNMF_L3_RUN_ROOT", "/home/h2048/output/program_full_parallel_methods_20260507"))
OUTPUT_SUBDIR = os.environ.get("CNMF_L3_OUTPUT_SUBDIR", "cnmf_by_celltype")
METHOD_DIRNAME = "cnmf_full"
STAMP = os.environ.get("CNMF_L3_AUDIT_STAMP", "20260519")

LOW_INFO_PATTERN = re.compile(
    r"^(?:MT-|RPS[0-9A-Z]*$|RPL[0-9A-Z]*$|MRPS[0-9A-Z]*$|MRPL[0-9A-Z]*$|"
    r"(?:RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$|(?:AC|AL|AP|BX|Z)[0-9]+[.]|"
    r"RP[0-9]+-|CTD-|CTB-|CTC-|LOC[0-9]+)|-OT[0-9]+$",
    flags=re.IGNORECASE,
)


def read_json(path: Path) -> Dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def safe_get(d: Dict[str, Any], key: str, default: Any = None) -> Any:
    value = d.get(key, default)
    if isinstance(value, (list, tuple)) and len(value) == 1:
        return value[0]
    return value


def is_legacy_cnmf_dir(celltype_root: Path, status: Dict[str, Any]) -> bool:
    if str(status.get("path_layout", "")).strip() == "celltype_method":
        return False
    return celltype_root.name.startswith("cnmf_")


def iter_cnmf_l3_status_records() -> List[tuple[Path, Dict[str, Any]]]:
    candidates = set(RUN_ROOT.glob(f"*/{OUTPUT_SUBDIR}/*/{METHOD_DIRNAME}/cnmf_l3_status.json"))
    candidates.update(RUN_ROOT.glob(f"*/{OUTPUT_SUBDIR}/cnmf_*/{METHOD_DIRNAME}/cnmf_l3_status.json"))
    ranked: List[tuple[tuple[str, str], bool, bool, Path, Dict[str, Any]]] = []
    for path in sorted(candidates):
        status = read_json(path)
        if not status:
            continue
        lineage = str(safe_get(status, "lineage", path.parts[-5]))
        celltype_root = path.parent.parent
        safe_celltype = str(safe_get(status, "safe_celltype", "")).strip()
        if not safe_celltype:
            safe_celltype = celltype_root.name[5:] if celltype_root.name.startswith("cnmf_") else celltype_root.name
        celltype_l3 = str(safe_get(status, "celltype_l3", "")).strip()
        key = (lineage, celltype_l3 or safe_celltype)
        ranked.append(
            (
                key,
                str(safe_get(status, "status", "")).strip() != "ok",
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


def audit_score_tables(cnmf_dir: Path) -> Dict[str, Any]:
    gep_dir = cnmf_dir / "gep_gene_tables"
    score_files = sorted(gep_dir.glob("gep_gene_scores_k*.tsv")) if gep_dir.exists() else []
    unit_count = 0
    low_info_rows = 0
    top_genes_files = 0
    malformed_files: List[str] = []
    for score_path in score_files:
        try:
            df = pd.read_csv(score_path, sep="\t")
        except Exception:
            malformed_files.append(str(score_path))
            continue
        if "gep" in df.columns:
            unit_count += int(df["gep"].astype(str).nunique())
        if "gene" in df.columns:
            low_info_rows += int(df["gene"].astype(str).str.contains(LOW_INFO_PATTERN, regex=True, na=False).sum())
        else:
            malformed_files.append(str(score_path))
    if gep_dir.exists():
        top_genes_files = len(list(gep_dir.glob("gep_top_genes_k*.csv")))
    return {
        "n_score_files": len(score_files),
        "n_top_genes_files": top_genes_files,
        "n_k_gep_units": unit_count,
        "low_info_rows_in_scores": low_info_rows,
        "malformed_score_files": ";".join(malformed_files),
    }


def audit_llm(lineage_dir: Path, safe_celltype: str) -> Dict[str, Any]:
    llm_root = lineage_dir / OUTPUT_SUBDIR / "llm_parallel" / "cnmf" / "units"
    if not llm_root.exists():
        return {"llm_units": 0, "llm_ok": 0, "llm_error": 0, "llm_queued": 0, "llm_other": 0}
    status_files = sorted(llm_root.glob(f"cnmf_l3_{safe_celltype}_*/cnmf_l3_{safe_celltype}_*_LLM_status.json"))
    counts = {"ok": 0, "error": 0, "queued": 0, "other": 0}
    for path in status_files:
        x = read_json(path)
        status = str(x.get("status", "other"))
        if status == "ok":
            counts["ok"] += 1
        elif status == "error":
            counts["error"] += 1
        elif status.startswith("queued") or status.startswith("skipped"):
            counts["queued"] += 1
        else:
            counts["other"] += 1
    return {
        "llm_units": len(status_files),
        "llm_ok": counts["ok"],
        "llm_error": counts["error"],
        "llm_queued": counts["queued"],
        "llm_other": counts["other"],
    }


def collect_rows() -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for status_path, status in iter_cnmf_l3_status_records():
        cnmf_dir = status_path.parent
        lineage = str(safe_get(status, "lineage", status_path.parts[-5]))
        safe_celltype = str(safe_get(status, "safe_celltype", status_path.parent.parent.name.replace("cnmf_", "", 1)))
        row = {
            "lineage": lineage,
            "celltype_l3": safe_get(status, "celltype_l3"),
            "safe_celltype": safe_celltype,
            "status": safe_get(status, "status"),
            "n_cells": safe_get(status, "n_cells"),
            "n_samples": safe_get(status, "n_samples"),
            "n_tissues": safe_get(status, "n_tissues"),
            "recommended_k": safe_get(status, "recommended_k"),
            "elapsed_min": safe_get(status, "elapsed_min"),
            "low_info_removed_rows_post_export": safe_get(status, "low_info_removed_rows_post_export"),
            "output_dir": str(cnmf_dir),
            "error": safe_get(status, "error"),
            "status_json": str(status_path),
            "has_run_summary": (cnmf_dir / "run_summary.json").exists(),
            "has_k_metrics": (cnmf_dir / "k_stability_metrics.json").exists(),
            "has_recommendation": (cnmf_dir / "k_selection_recommendation.json").exists(),
            "n_visualizations": len(list((cnmf_dir / "visualizations").glob("*.png"))) if (cnmf_dir / "visualizations").exists() else 0,
        }
        row.update(audit_score_tables(cnmf_dir))
        row.update(audit_llm(RUN_ROOT / lineage, safe_celltype))
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    return df.sort_values(["lineage", "celltype_l3"], na_position="last")


def main() -> int:
    df = collect_rows()
    detail_path = RUN_ROOT / f"cnmf_l3_final_audit_{STAMP}.tsv"
    df.to_csv(detail_path, sep="\t", index=False)
    if df.empty:
        summary = {"n_rows": 0, "detail_path": str(detail_path)}
    else:
        by_status = df["status"].fillna("NA").value_counts().to_dict()
        summary = {
            "n_rows": int(len(df)),
            "status_counts": {str(k): int(v) for k, v in by_status.items()},
            "n_ok": int((df["status"] == "ok").sum()),
            "n_error": int((df["status"] == "error").sum()),
            "n_started": int((df["status"] == "started").sum()),
            "n_planned": int((df["status"] == "planned").sum()),
            "total_k_gep_units": int(pd.to_numeric(df.get("n_k_gep_units", 0), errors="coerce").fillna(0).sum()),
            "total_score_files": int(pd.to_numeric(df.get("n_score_files", 0), errors="coerce").fillna(0).sum()),
            "low_info_rows_in_scores": int(pd.to_numeric(df.get("low_info_rows_in_scores", 0), errors="coerce").fillna(0).sum()),
            "llm_units": int(pd.to_numeric(df.get("llm_units", 0), errors="coerce").fillna(0).sum()),
            "llm_ok": int(pd.to_numeric(df.get("llm_ok", 0), errors="coerce").fillna(0).sum()),
            "llm_error": int(pd.to_numeric(df.get("llm_error", 0), errors="coerce").fillna(0).sum()),
            "llm_queued": int(pd.to_numeric(df.get("llm_queued", 0), errors="coerce").fillna(0).sum()),
            "detail_path": str(detail_path),
        }
    summary_path = RUN_ROOT / f"cnmf_l3_final_audit_summary_{STAMP}.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[OK] detail={detail_path}")
    print(f"[OK] summary={summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())