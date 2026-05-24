#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Filter an AnnData H5AD by a removed-cell manifest.

Used for 2026-05-08 full reruns where outlier clusters were identified in the
R tissue-comparison output and the upstream scVI/scANVI input H5AD needs to be
recreated with those cells removed before fresh training.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd


def read_removed_cells(path: Path, column: str | None = None) -> pd.Index:
    df = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else None, engine="python")
    candidates = [column] if column else []
    candidates += ["cell", "cell_barcode", "obs_name", "barcode", "cell_id"]
    for col in candidates:
        if col and col in df.columns:
            vals = df[col].astype(str).str.strip()
            vals = vals[vals.ne("") & vals.ne("nan")]
            return pd.Index(vals.unique())
    raise KeyError(f"Could not detect removed-cell column in {path}; columns={list(df.columns)}")


def match_removed_to_obs(removed: pd.Index, obs_names: pd.Index) -> tuple[pd.Index, str]:
    direct = removed.intersection(obs_names)
    best = (direct, "exact")

    if len(direct) < len(removed):
        obs_stripped = pd.Index([x.split("::")[-1] for x in obs_names.astype(str)])
        rem_stripped = pd.Index([x.split("::")[-1] for x in removed.astype(str)])
        matched_stripped = rem_stripped.intersection(obs_stripped)
        if len(matched_stripped) > len(best[0]):
            # Map stripped IDs back to original obs names.
            stripped_to_obs = pd.Series(obs_names.astype(str), index=obs_stripped.astype(str))
            best = (pd.Index(stripped_to_obs.loc[matched_stripped].astype(str)), "strip_prefix")

    return best


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter an H5AD by removed cells")
    parser.add_argument("--input-h5ad", required=True, type=Path)
    parser.add_argument("--removed-cells", required=True, type=Path)
    parser.add_argument("--output-h5ad", required=True, type=Path)
    parser.add_argument("--removed-column", default=None)
    parser.add_argument("--min-match-rate", type=float, default=0.95)
    parser.add_argument("--summary-json", type=Path, default=None)
    parser.add_argument("--matched-cells-tsv", type=Path, default=None)
    parser.add_argument(
        "--compression",
        choices=["none", "gzip"],
        default="none",
        help="H5AD write compression. Default none is much faster for large rerun inputs.",
    )
    parser.add_argument("--compression-level", type=int, default=1)
    args = parser.parse_args()

    removed = read_removed_cells(args.removed_cells, column=args.removed_column)
    if len(removed) == 0:
        raise RuntimeError(f"No removed cells loaded from {args.removed_cells}")

    print(f"[filter] input h5ad: {args.input_h5ad}")
    print(f"[filter] removed manifest: {args.removed_cells} ({len(removed):,} cells)")
    adata = ad.read_h5ad(args.input_h5ad)
    obs_names = pd.Index(adata.obs_names.astype(str))
    matched, mode = match_removed_to_obs(removed, obs_names)
    match_rate = len(matched) / max(1, len(removed))
    print(f"[filter] match mode={mode} matched={len(matched):,}/{len(removed):,} rate={match_rate:.4f}")
    if match_rate < args.min_match_rate:
        missing = removed.difference(matched)
        raise RuntimeError(
            f"Removed-cell match rate too low ({match_rate:.4f} < {args.min_match_rate:.4f}). "
            f"First missing IDs: {missing[:10].tolist()}"
        )

    keep_mask = ~obs_names.isin(matched)
    filtered = adata[keep_mask].copy()
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    if args.compression == "gzip":
        filtered.write_h5ad(args.output_h5ad, compression="gzip", compression_opts=args.compression_level)
    else:
        filtered.write_h5ad(args.output_h5ad)
    print(f"[filter] wrote: {args.output_h5ad} shape={filtered.shape}")

    if args.matched_cells_tsv:
        args.matched_cells_tsv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"cell": matched.astype(str)}).to_csv(args.matched_cells_tsv, sep="\t", index=False)
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(
            json.dumps(
                {
                    "input_h5ad": str(args.input_h5ad),
                    "removed_cells_manifest": str(args.removed_cells),
                    "output_h5ad": str(args.output_h5ad),
                    "input_n_obs": int(adata.n_obs),
                    "input_n_vars": int(adata.n_vars),
                    "removed_requested": int(len(removed)),
                    "removed_matched": int(len(matched)),
                    "match_rate": float(match_rate),
                    "match_mode": mode,
                    "output_n_obs": int(filtered.n_obs),
                    "output_n_vars": int(filtered.n_vars),
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
