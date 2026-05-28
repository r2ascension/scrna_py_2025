#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from non_unified_airway.communication import (
    DEFAULT_COMM_OUT_DIR,
    flatten_cellphonedb_table,
    write_status,
)
from non_unified_airway.common import write_json, write_tsv


def write_table(path: Path, table: pd.DataFrame) -> None:
    kwargs = {"compression": "gzip"} if str(path).endswith(".gz") else {}
    write_tsv(path, table, **kwargs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CellPhoneDB pairwise communication for prepared L2/L3 inputs.")
    parser.add_argument("--prepared-manifest", type=Path, default=DEFAULT_COMM_OUT_DIR / "prepared_inputs_manifest.tsv")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_COMM_OUT_DIR / "cellphonedb")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def _cellphonedb_available() -> tuple[bool, str]:
    executable = shutil.which("cellphonedb")
    if executable:
        return True, executable
    return False, "cellphonedb"


def _write_side_metadata(full_meta: pd.DataFrame, side: str, out_path: Path) -> pd.DataFrame:
    side_meta = full_meta.loc[full_meta["contrast_side"].astype(str) == side, ["cell_id", "comm_celltype"]].copy()
    side_meta = side_meta.rename(columns={"comm_celltype": "cell_type"})
    side_meta.to_csv(out_path, sep="\t", index=False)
    return side_meta


def _write_side_h5ad(source_h5ad: Path, side_meta: pd.DataFrame, side_h5ad: Path) -> None:
    import anndata as ad

    ad.settings.allow_write_nullable_strings = True

    adata = ad.read_h5ad(source_h5ad)
    side_ids = side_meta["cell_id"].astype(str).tolist()
    side_subset = adata[side_ids].copy()
    if getattr(side_subset, "X", None) is None and "counts" in side_subset.layers:
        side_subset.X = side_subset.layers["counts"].copy()
    elif "counts" in side_subset.layers:
        side_subset.X = side_subset.layers["counts"].copy()
    side_subset.obs["cell_type"] = side_meta.set_index("cell_id").loc[side_subset.obs_names, "cell_type"].astype(str).values
    side_subset.write_h5ad(side_h5ad, compression="gzip")


def _run_cpdb(meta_path: Path, h5ad_path: Path, out_dir: Path, *, threads: int, iterations: int) -> None:
    cmd = [
        "cellphonedb",
        "method",
        "statistical_analysis",
        str(meta_path),
        str(h5ad_path),
        "--output-path",
        str(out_dir),
        "--threads",
        str(threads),
        "--iterations",
        str(iterations),
    ]
    subprocess.run(cmd, check=True)


def _load_long_result(side_dir: Path, table_name: str, value_name: str) -> pd.DataFrame:
    table_path = side_dir / table_name
    if not table_path.exists():
        return pd.DataFrame()
    wide = pd.read_csv(table_path, sep="\t")
    long = flatten_cellphonedb_table(wide, value_name=value_name)
    long["interaction_key"] = (
        long["sender"].astype(str)
        + "|"
        + long["receiver"].astype(str)
        + "|"
        + long["ligand"].astype(str)
        + "|"
        + long["receptor"].astype(str)
    )
    return long


def _build_delta(left_long: pd.DataFrame, right_long: pd.DataFrame, *, contrast_id: str, level: str, left_label: str, right_label: str) -> pd.DataFrame:
    left = left_long[["interaction_key", "sender", "receiver", "ligand", "receptor", "mean"]].rename(columns={"mean": "left_mean"}) if not left_long.empty and "mean" in left_long.columns else pd.DataFrame(columns=["interaction_key", "sender", "receiver", "ligand", "receptor", "left_mean"])
    right = right_long[["interaction_key", "sender", "receiver", "ligand", "receptor", "mean"]].rename(columns={"mean": "right_mean"}) if not right_long.empty and "mean" in right_long.columns else pd.DataFrame(columns=["interaction_key", "sender", "receiver", "ligand", "receptor", "right_mean"])
    merged = left.merge(right, on=["interaction_key", "sender", "receiver", "ligand", "receptor"], how="outer")
    merged["left_mean"] = pd.to_numeric(merged.get("left_mean"), errors="coerce")
    merged["right_mean"] = pd.to_numeric(merged.get("right_mean"), errors="coerce")
    merged["delta_score_right_minus_left"] = merged["right_mean"].fillna(0) - merged["left_mean"].fillna(0)
    merged["method"] = "cellphonedb"
    merged["contrast_id"] = contrast_id
    merged["level"] = level
    merged["left_label"] = left_label
    merged["right_label"] = right_label
    return merged.sort_values("delta_score_right_minus_left", ascending=False).reset_index(drop=True)


def run_cellphonedb_pairwise(args: argparse.Namespace) -> dict[str, object]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.prepared_manifest, sep="\t")
    runnable = manifest.loc[manifest["status"].astype(str) == "ready"].copy()
    available, executable = _cellphonedb_available()
    if not available:
        payload = {
            "status": "skipped_missing_dependency",
            "dependency": executable,
            "n_inputs": int(runnable.shape[0]),
        }
        write_json(args.out_dir / "run_summary.json", payload)
        return payload

    status_rows: list[dict[str, object]] = []
    for row in runnable.to_dict(orient="records"):
        contrast_id = str(row["contrast_id"])
        level = str(row["level"])
        run_dir = args.out_dir / contrast_id / level
        run_dir.mkdir(parents=True, exist_ok=True)

        source_h5ad = Path(str(row.get("h5ad_subset", "")))
        meta_path = Path(str(row.get("cell_metadata_tsv", "")))
        if not source_h5ad.exists() or not meta_path.exists():
            status = {
                "contrast_id": contrast_id,
                "level": level,
                "status": "blocked_missing_inputs",
                "h5ad_subset": str(source_h5ad),
                "cell_metadata_tsv": str(meta_path),
            }
            write_status(run_dir / "cellphonedb_status.json", status)
            status_rows.append(status)
            continue

        result_path = run_dir / "cellphonedb_delta.tsv.gz"
        if result_path.exists() and not args.force:
            status = {
                "contrast_id": contrast_id,
                "level": level,
                "status": "skipped_exists",
                "result": str(result_path),
            }
            write_status(run_dir / "cellphonedb_status.json", status)
            status_rows.append(status)
            continue

        try:
            full_meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
            side_results: dict[str, pd.DataFrame] = {}
            for side in ["left", "right"]:
                side_dir = run_dir / side
                side_dir.mkdir(parents=True, exist_ok=True)
                side_meta_path = side_dir / "meta.tsv"
                side_meta = _write_side_metadata(full_meta, side, side_meta_path)
                side_h5ad = side_dir / "adata_side.h5ad"
                _write_side_h5ad(source_h5ad, side_meta, side_h5ad)
                _run_cpdb(side_meta_path, side_h5ad, side_dir, threads=args.threads, iterations=args.iterations)
                long_means = _load_long_result(side_dir, "means.txt", "mean")
                side_results[side] = long_means
                write_table(side_dir / "means_long.tsv.gz", long_means)

            delta_df = _build_delta(
                side_results.get("left", pd.DataFrame()),
                side_results.get("right", pd.DataFrame()),
                contrast_id=contrast_id,
                level=level,
                left_label=str(row["left_label"]),
                right_label=str(row["right_label"]),
            )
            write_table(result_path, delta_df)
            status = {
                "contrast_id": contrast_id,
                "level": level,
                "status": "ok",
                "n_left": int(side_results.get("left", pd.DataFrame()).shape[0]),
                "n_right": int(side_results.get("right", pd.DataFrame()).shape[0]),
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
        write_status(run_dir / "cellphonedb_status.json", status)
        status_rows.append(status)

    status_df = pd.DataFrame(status_rows)
    write_table(args.out_dir / "cellphonedb_status.tsv", status_df)
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
    summary = run_cellphonedb_pairwise(args)
    print(
        "[non_unified_airway][05c_v2] "
        f"inputs={summary.get('n_inputs', 0)} ok={summary.get('n_ok', 0)} error={summary.get('n_error', 0)}"
    )


if __name__ == "__main__":
    main()