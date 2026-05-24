#!/usr/bin/env python3
from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Rebuild a final h5ad so that X/layers use full-gene matrices from a "
            "reference h5ad while preserving obs/obsm/uns from the final h5ad."
        )
    )
    p.add_argument("--final-h5ad", required=True, help="Existing final h5ad with final obs/obsm/uns")
    p.add_argument("--reference-h5ad", required=True, help="Reference h5ad containing full-gene X/layers")
    p.add_argument("--output-h5ad", required=True, help="Output full-gene h5ad path")
    p.add_argument("--final-cellid-col", default="", help="Optional obs column in final h5ad used as match key; default uses obs_names")
    p.add_argument("--reference-cellid-col", default="", help="Optional obs column in reference h5ad used as match key; default uses obs_names")
    p.add_argument("--expr-layer-candidate", action="append", dest="expr_layers", default=[], help="Preferred normalized-expression layer name in reference h5ad; can be repeated")
    p.add_argument("--counts-layer-candidate", action="append", dest="counts_layers", default=[], help="Preferred counts layer name in reference h5ad; can be repeated")
    p.add_argument("--expr-layer-output-name", default="", help="Optional layer name to use when storing expression matrix in output; default keeps source name except raw.X -> log1p")
    p.add_argument("--compression", default="gzip", help="h5ad compression (default: gzip)")
    return p.parse_args()


def get_match_ids(adata: ad.AnnData, col: str) -> pd.Index:
    if col:
        if col not in adata.obs.columns:
            raise KeyError(f"Match column not found in obs: {col}")
        return pd.Index(adata.obs[col].astype(str), name=col)
    return pd.Index(adata.obs_names.astype(str), name="obs_names")


def copy_matrix(mat):
    return mat.copy() if hasattr(mat, "copy") else np.array(mat, copy=True)


def expm1_matrix(mat):
    if sparse.issparse(mat):
        out = mat.copy()
        out.data = np.expm1(out.data)
        rounded = np.rint(out.data)
        if np.allclose(out.data, rounded, atol=1e-6, rtol=0):
            out.data = rounded.astype(np.int64, copy=False)
        return out

    out = np.expm1(np.asarray(mat))
    rounded = np.rint(out)
    if np.allclose(out, rounded, atol=1e-6, rtol=0):
        out = rounded.astype(np.int64, copy=False)
    return out


def resolve_matrix_source(adata: ad.AnnData, candidate: str):
    name = candidate.strip()

    if name in {"", "X"}:
        return copy_matrix(adata.X), "X", adata.var.copy()

    if name in {"raw", "raw.X", "raw_x"}:
        if adata.raw is None:
            raise KeyError("raw is not available")
        return copy_matrix(adata.raw.X), "raw.X", adata.raw.var.copy()

    if name in {"raw.X_expm1", "raw_expm1", "expm1(raw.X)", "expm1_raw.X", "expm1_raw_x"}:
        if adata.raw is None:
            raise KeyError("raw is not available")
        return expm1_matrix(adata.raw.X), "raw.X_expm1", adata.raw.var.copy()

    if name in adata.layers:
        return copy_matrix(adata.layers[name]), name, adata.var.copy()

    raise KeyError(f"candidate not found: {candidate}")


def pick_matrix(adata: ad.AnnData, candidates: list[str], target_n_vars: int | None = None):
    failures: list[str] = []
    for candidate in candidates:
        try:
            mat, source_name, var_df = resolve_matrix_source(adata, candidate)
        except Exception as exc:  # pragma: no cover - diagnostic path
            failures.append(f"{candidate}: {exc}")
            continue

        if target_n_vars is not None and mat.shape[1] != target_n_vars:
            failures.append(
                f"{candidate}: n_vars={mat.shape[1]} does not match target_n_vars={target_n_vars}"
            )
            continue

        return mat, source_name, var_df

    raise ValueError(
        "No usable matrix source found. Candidates tried: " + "; ".join(failures or candidates)
    )


def main() -> None:
    args = parse_args()
    expr_candidates = args.expr_layers or ["log1p", "data", "X"]
    counts_candidates = args.counts_layers or ["counts", "raw_counts", "raw.X_expm1"]

    ad.settings.allow_write_nullable_strings = True

    final_path = Path(args.final_h5ad)
    ref_path = Path(args.reference_h5ad)
    out_path = Path(args.output_h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    final = ad.read_h5ad(final_path, backed="r")
    ref = ad.read_h5ad(ref_path, backed="r")

    final_ids = get_match_ids(final, args.final_cellid_col)
    ref_ids = get_match_ids(ref, args.reference_cellid_col)

    if ref_ids.has_duplicates:
        dup = ref_ids[ref_ids.duplicated()].unique().tolist()[:10]
        raise ValueError(f"Reference match ids contain duplicates: {dup}")
    if final_ids.has_duplicates:
        dup = final_ids[final_ids.duplicated()].unique().tolist()[:10]
        raise ValueError(f"Final match ids contain duplicates: {dup}")

    indexer = ref_ids.get_indexer(final_ids)
    missing_mask = indexer < 0
    if missing_mask.any():
        missing = final_ids[missing_mask].tolist()[:10]
        raise ValueError(
            f"{missing_mask.sum()} final cells are missing in reference h5ad; examples: {missing}"
        )

    ref_sub = ref[indexer, :].to_memory()

    expr_mat, expr_source, expr_var = pick_matrix(ref_sub, expr_candidates)
    counts_mat, counts_source, counts_var = pick_matrix(
        ref_sub,
        counts_candidates,
        target_n_vars=expr_mat.shape[1],
    )

    if not expr_var.index.equals(counts_var.index):
        raise ValueError(
            "Expression and counts sources resolved to different gene sets; "
            f"expr_source={expr_source}, counts_source={counts_source}"
        )

    expr_output_name = args.expr_layer_output_name.strip() or (
        "log1p" if expr_source == "raw.X" else expr_source
    )

    out = ad.AnnData(X=expr_mat, obs=final.obs.copy(), var=expr_var.copy())
    out.obs_names = final.obs_names.copy()
    out.var_names = expr_var.index.copy()

    if "gene_symbol" not in out.var.columns:
        if "symbol_base" in out.var.columns:
            out.var["gene_symbol"] = out.var["symbol_base"].astype(str)
        else:
            out.var["gene_symbol"] = out.var_names.astype(str)
    out.var["hvg_in_original_final"] = out.var_names.isin(final.var_names)

    out.layers["counts"] = counts_mat
    if expr_output_name and expr_output_name != "X":
        out.layers[expr_output_name] = expr_mat.copy()

    for key in final.obsm.keys():
        out.obsm[key] = final.obsm[key].copy()
    for key in final.obsp.keys():
        out.obsp[key] = final.obsp[key].copy()

    out.uns = deepcopy(final.uns)
    out.uns["X_layer"] = expr_output_name or "X"
    out.uns["tc_full_gene_reference_h5ad"] = str(ref_path.resolve())
    out.uns["tc_full_gene_reference_match_key"] = args.reference_cellid_col or "__obs_names__"
    out.uns["tc_full_gene_final_match_key"] = args.final_cellid_col or "__obs_names__"
    out.uns["tc_full_gene_counts_source"] = counts_source
    out.uns["tc_full_gene_expr_source"] = expr_source
    out.uns["tc_full_gene_expr_layer_name"] = expr_output_name or "X"

    out.write_h5ad(out_path, compression=args.compression)

    if getattr(final, "isbacked", False):
        final.file.close()
    if getattr(ref, "isbacked", False):
        ref.file.close()

    check = ad.read_h5ad(out_path, backed="r")
    print("[OK] wrote", out_path)
    print("[OK] final cells matched to reference:", check.n_obs)
    print("[OK] full-gene shape:", check.shape)
    print("[OK] layers:", list(check.layers.keys()))
    print("[OK] first obs:", list(check.obs_names[:3]))
    print("[OK] first var:", list(check.var_names[:5]))
    check.file.close()


if __name__ == "__main__":
    main()
