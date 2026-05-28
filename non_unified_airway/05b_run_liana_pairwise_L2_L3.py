#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from non_unified_airway.communication import (
    DEFAULT_COMM_OUT_DIR,
    merge_method_side_results,
    module_available,
    standardize_interaction_columns,
    write_status,
)
from non_unified_airway.common import write_json, write_tsv


def write_table(path: Path, table: pd.DataFrame) -> None:
    kwargs = {"compression": "gzip"} if str(path).endswith(".gz") else {}
    write_tsv(path, table, **kwargs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LIANA+ pairwise communication for prepared L2/L3 inputs.")
    parser.add_argument("--prepared-manifest", type=Path, default=DEFAULT_COMM_OUT_DIR / "prepared_inputs_manifest.tsv")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_COMM_OUT_DIR / "liana")
    parser.add_argument("--resource-name", default="consensus")
    parser.add_argument("--groupby", default="comm_celltype")
    parser.add_argument("--expr-prop", type=float, default=0.1)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def _extract_liana_result(adata) -> pd.DataFrame:
    candidate_keys = [
        "liana_res",
        "liana_result",
        "rank_aggregate",
        "liana_rank_aggregate",
        "liana",
    ]
    for key in candidate_keys:
        if key not in adata.uns:
            continue
        obj = adata.uns[key]
        if isinstance(obj, pd.DataFrame):
            return obj.copy()
        try:
            return pd.DataFrame(obj)
        except Exception:
            continue
    raise KeyError(f"Could not find a LIANA result table in AnnData.uns; keys={list(adata.uns.keys())}")


def _run_liana_side(adata, side: str, *, groupby: str, resource_name: str, expr_prop: float) -> pd.DataFrame:
    import liana as li  # type: ignore
    import scanpy as sc

    side_data = adata[adata.obs["contrast_side"].astype(str) == side].copy()
    if side_data.n_obs == 0:
        return pd.DataFrame()
    if side_data.obs[groupby].astype(str).nunique() < 2:
        return pd.DataFrame()

    if "counts" in side_data.layers:
        side_data.X = side_data.layers["counts"].copy()
        sc.pp.normalize_total(side_data, target_sum=1e4)
        sc.pp.log1p(side_data)

    try:
        li.mt.rank_aggregate(
            side_data,
            groupby=groupby,
            resource_name=resource_name,
            expr_prop=expr_prop,
            use_raw=False,
            verbose=True,
        )
    except TypeError:
        li.mt.rank_aggregate(
            side_data,
            groupby=groupby,
            resource_name=resource_name,
            expr_prop=expr_prop,
            verbose=True,
        )

    result = _extract_liana_result(side_data)
    return standardize_interaction_columns(result)


def run_liana_pairwise(args: argparse.Namespace) -> dict[str, object]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.prepared_manifest, sep="\t")
    runnable = manifest.loc[manifest["status"].astype(str) == "ready"].copy()

    if not module_available("liana"):
        payload = {
            "status": "skipped_missing_dependency",
            "dependency": "liana",
            "n_inputs": int(runnable.shape[0]),
        }
        write_json(args.out_dir / "run_summary.json", payload)
        return payload

    import anndata as ad

    status_rows: list[dict[str, object]] = []
    for row in runnable.to_dict(orient="records"):
        contrast_id = str(row["contrast_id"])
        level = str(row["level"])
        run_dir = args.out_dir / contrast_id / level
        run_dir.mkdir(parents=True, exist_ok=True)
        adata_path = Path(str(row.get("h5ad_subset", "")))
        if not adata_path.exists():
            status = {
                "contrast_id": contrast_id,
                "level": level,
                "status": "blocked_missing_h5ad_subset",
                "h5ad_subset": str(adata_path),
            }
            write_status(run_dir / "liana_status.json", status)
            status_rows.append(status)
            continue

        result_path = run_dir / "liana_delta.tsv.gz"
        if result_path.exists() and not args.force:
            status = {
                "contrast_id": contrast_id,
                "level": level,
                "status": "skipped_exists",
                "result": str(result_path),
            }
            write_status(run_dir / "liana_status.json", status)
            status_rows.append(status)
            continue

        try:
            adata = ad.read_h5ad(adata_path)
            left_df = _run_liana_side(
                adata,
                "left",
                groupby=args.groupby,
                resource_name=args.resource_name,
                expr_prop=args.expr_prop,
            )
            right_df = _run_liana_side(
                adata,
                "right",
                groupby=args.groupby,
                resource_name=args.resource_name,
                expr_prop=args.expr_prop,
            )
            write_table(run_dir / "liana_left.tsv.gz", left_df)
            write_table(run_dir / "liana_right.tsv.gz", right_df)
            delta_df = merge_method_side_results(
                left_df,
                right_df,
                method="liana",
                left_label=str(row["left_label"]),
                right_label=str(row["right_label"]),
                contrast_id=contrast_id,
                level=level,
            )
            if not delta_df.empty:
                rank_cols = [col for col in ["aggregate_rank_left", "aggregate_rank_right", "magnitude_rank_left", "magnitude_rank_right"] if col in delta_df.columns]
                for column in rank_cols:
                    delta_df[column] = pd.to_numeric(delta_df[column], errors="coerce")
            write_table(result_path, delta_df)
            status = {
                "contrast_id": contrast_id,
                "level": level,
                "status": "ok",
                "n_left": int(left_df.shape[0]),
                "n_right": int(right_df.shape[0]),
                "n_delta": int(delta_df.shape[0]),
                "result": str(result_path),
            }
        except Exception as exc:  # pragma: no cover - runtime path
            status = {
                "contrast_id": contrast_id,
                "level": level,
                "status": "error",
                "error": str(exc),
            }
        write_status(run_dir / "liana_status.json", status)
        status_rows.append(status)

    status_df = pd.DataFrame(status_rows)
    write_table(args.out_dir / "liana_status.tsv", status_df)
    summary = {
        "status": "ok",
        "n_inputs": int(runnable.shape[0]),
        "n_ok": int(status_df["status"].eq("ok").sum()) if not status_df.empty else 0,
        "n_error": int(status_df["status"].eq("error").sum()) if not status_df.empty else 0,
        "n_blocked": int(status_df["status"].str.startswith("blocked").sum()) if not status_df.empty else 0,
    }
    write_json(args.out_dir / "run_summary.json", summary)
    return summary


def main() -> None:
    args = parse_args()
    summary = run_liana_pairwise(args)
    print(
        "[non_unified_airway][05b] "
        f"inputs={summary.get('n_inputs', 0)} ok={summary.get('n_ok', 0)} error={summary.get('n_error', 0)}"
    )


if __name__ == "__main__":
    main()
