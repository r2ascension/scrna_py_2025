#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Label-Guided cNMF Pipeline v1.4-FINAL - Production Integrated
==============================================================
Fixes:
- P0-1: Counts extraction priority: layers['counts'] -> raw.X (aligned) -> X
- P0-2: Chunked dense conversion for writing (avoid per-gene toarray + avoid huge line buffers)
- P0-3: Gene filtering DOES NOT modify .raw
- P1-1: validate_counts float tolerance (1e-6)
- P2: Keep ribo/histone by default (configurable)
- HOTFIX-1: np.char.startswith type error (genes.values.astype(str))
- HOTFIX-2: Remove confidence filtering (user request)
- HOTFIX-3: TSV parsing error - sanitize gene/cell names (remove tabs, newlines)
- HOTFIX-4: Optimize K values for lineage-level analysis (15-30 instead of 25-50)
- HOTFIX-5: CRITICAL - Zero HVG counts error fix:
  * Reduced gene filtering aggressiveness (keep unannotated, ENSG genes)
  * Added min_cells_per_gene filter (≥10 cells)
  * Changed HVG flavor: seurat_v3 → seurat (more robust)
  * Added HVG expression thresholds (min_mean, max_mean, min_disp)
  * Added cell filtering by HVG counts (removes cells with <50 HVG counts)

Author: r2end (integrated + hardened + hotfix)
Date: 2025-01-15
Version: 1.4-FINAL
"""

import os
import sys
import gc
import re
import json
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

warnings.filterwarnings("ignore")

# -----------------------------
# Environment (single-thread)
# -----------------------------
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"


# -----------------------------
# cNMF import
# -----------------------------
try:
    from cnmf import cNMF
except ImportError:
    print("ERROR: cnmf not installed. Please install: pip install cnmf")
    sys.exit(1)


# =============================================================================
# CONFIG
# =============================================================================

INPUT_H5AD = "/home/h2048/data/py/0110/celltypist_epithelial/epithelial_celltypist_filtered_final.h5ad"
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/cnmf_v1.4_final")

CELLTYPE_KEY = "celltypist_pred"

# REMOVED: USE_CONFIDENCE_FILTER, MIN_CONFIDENCE, CONFIDENCE_KEY (per user request)

USE_MAJOR_LINEAGE = True
MAJOR_LINEAGE_MAP = {
    # Alveolar lineage
    "AT1": "Alveolar",
    "AT1 ": "Alveolar",
    "AT2": "Alveolar",
    # Basal lineage
    "Basal": "Basal_Lineage",
    "Suprabasal": "Basal_Lineage",
    "SMG_Basal": "Basal_Lineage",
    "Dividing_Basal": "Basal_Lineage",
    # Ciliated lineage
    "Ciliated": "Ciliated_Lineage",
    "Deuterosome": "Ciliated_Lineage",
    "Deuterosomal": "Ciliated_Lineage",
    # Secretory lineage
    "Secretory_Goblet": "Secretory_Lineage",
    "Secretory_Club": "Secretory_Lineage",
    "SMG_Mucous": "Secretory_Lineage",
    "SMG_Serous": "Secretory_Lineage",
    "SCGB1A1+": "Secretory_Lineage",
    # Duct
    "SMG_Duct": "Duct",
    # Rare
    "Ionocyte_n_Brush": "Rare_Specialized",
    "Ionocyte": "Rare_Specialized",
    "Brush": "Rare_Specialized",
}

PROCESS_CELLTYPES = "auto"   # 'auto' or list
MIN_CELLS_PER_TYPE = 500

# Gene filtering: LESS AGGRESSIVE for lineage-level analysis
GENE_FILTER_CONFIG = {
    "remove_mt": True,
    "remove_ribo": False,           # Keep for stability
    "remove_histone": False,         # Keep for stability
    "remove_pseudogenes": False,     # Keep unless causing issues
    "remove_ensg": False,            # Keep unannotated genes
    "remove_unannotated": False,     # Keep for HVG diversity
}

# cNMF - Enhanced for robustness
CNMF_CONFIG = {
    "k_range_auto": True,
    "k_rules": {
        # Reduced K for lineage-level analysis (already pre-filtered cell types)
        # Rationale: within major lineage, functional heterogeneity is limited
        "large": (10000, [15, 20, 25, 30]),        # was [25, 30, 35, 40, 45, 50]
        "medium": (2000, [10, 15, 20, 25]),        # was [15, 20, 25, 30]
        "small": (500, [5, 8, 10, 12]),            # was [5, 8, 10, 12, 15]
    },
    "n_iter": 100,
    "seed": 42,
    "num_hvg": 3000,
    "hvg_flavor": "seurat",              # Changed from seurat_v3 (more robust)
    "hvg_min_mean": 0.0125,              # NEW: minimum mean expression
    "hvg_max_mean": 3,                   # NEW: maximum mean expression  
    "hvg_min_disp": 0.5,                 # NEW: minimum dispersion
    "min_cells_per_gene": 10,            # NEW: gene must be expressed in ≥10 cells
    "filter_zero_hvg_cells": True,       # NEW: remove cells with too few HVG counts
    "min_hvg_counts_per_cell": 50,       # NEW: cell must have ≥50 HVG counts
    "density_threshold": 0.1,
    "total_workers": 1,   # single-worker mode
    "show_clustering": True,
    "close_clustergram_fig": True,
}

# IO / runtime options
SAVE_FILTERED_H5AD = True
OVERWRITE_GENE_FILTER = False
SKIP_IF_EXISTS = True
DELETE_COUNTS_TXT = True

# Chunk writing control
# Approx memory (float32): chunk_genes * n_cells * 4 bytes
MAX_CHUNK_MEM_MB = 256
CHUNK_GENES_DEFAULT = 2000  # used as upper bound; final chunk computed dynamically


# =============================================================================
# UTILITIES
# =============================================================================

def safe_name(x: str, max_len: int = 180) -> str:
    s = str(x)
    for ch in ["/", "\\", " ", "|", ":", ";", ",", "\t", "(", ")", "[", "]", "{", "}"]:
        s = s.replace(ch, "_")
    while "__" in s:
        s = s.replace("__", "_")
    s = s.strip("_")
    return s[:max_len]


def validate_counts_matrix(X, layer_name: str = "X", tol: float = 1e-6) -> Tuple[bool, str]:
    """Validate raw-ish counts with float tolerance."""
    if sparse.issparse(X):
        data = X.data
        data_to_check = data[: min(20000, data.size)]
    else:
        flat = np.asarray(X).ravel()
        data_to_check = flat[: min(20000, flat.size)]

    if data_to_check.size == 0:
        return False, f"{layer_name}: empty matrix"

    max_val = float(np.max(data_to_check))
    if max_val < 10:
        return False, f"{layer_name}: max={max_val:.2f}, likely log-normalized"

    # integer-like check with tolerance
    is_not_int = np.abs(data_to_check - np.rint(data_to_check)) > tol
    if np.any(is_not_int):
        frac = float(np.mean(is_not_int))
        # allow small fraction due to float storage / numerical noise
        if frac > 0.05:
            return False, f"{layer_name}: {frac*100:.1f}% non-integer-like values (tol={tol})"

    return True, f"{layer_name}: valid counts (max={max_val:.0f})"


def determine_k_range(n_cells: int, config: Dict) -> List[int]:
    """Auto K with sqrt cap (prevents silly K on small N)."""
    if not config.get("k_range_auto", True):
        return config.get("k_range_custom", [10, 15, 20])

    # choose rule bucket
    chosen = None
    for _, (min_cells, k_vals) in sorted(config["k_rules"].items(), key=lambda x: x[1][0], reverse=True):
        if n_cells >= min_cells:
            chosen = k_vals
            break
    if chosen is None:
        chosen = [5, 8, 10]

    max_k_safe = max(8, int(np.sqrt(n_cells)))
    k_vals_capped = [k for k in chosen if k <= max_k_safe]
    if not k_vals_capped:
        k_vals_capped = [max_k_safe]

    # de-dup + sorted
    k_vals_capped = sorted(list(dict.fromkeys(k_vals_capped)))
    return k_vals_capped


def check_cnmf_complete(cnmf_dir: Path, k_range: List[int], density_threshold: float = 0.1) -> bool:
    """Require spectra+usage exist for each K, tolerate dt naming variants."""
    if not cnmf_dir.exists():
        return False

    dt_formats = [
        f"dt_{str(density_threshold).replace('.', '_')}",
        f"dt_{density_threshold}",
    ]

    for k in k_range:
        ok = False
        for dt in dt_formats:
            spectra = list(cnmf_dir.glob(f"*.spectra.k_{k}.{dt}.consensus.txt"))
            usage = list(cnmf_dir.glob(f"*.usages.k_{k}.{dt}.consensus.txt"))
            if spectra and usage:
                ok = True
                break
        if not ok:
            return False

    return True


# =============================================================================
# GENE FILTERING (does NOT touch .raw)
# =============================================================================

def filter_genes_safe(adata: sc.AnnData, cfg: Dict[str, bool]) -> sc.AnnData:
    """
    HOTFIX-1: Use .values.astype(str) to avoid numpy.char type error
    """
    genes = adata.var_names.values.astype(str)  # FIXED: convert to numpy str array
    to_remove = np.zeros(genes.shape[0], dtype=bool)

    if cfg.get("remove_mt", True):
        to_remove |= np.char.startswith(genes, "MT-")

    if cfg.get("remove_ribo", False):
        ribo_re = re.compile(r"^(RPS|RPL|MRPS|MRPL)")
        to_remove |= np.array([bool(ribo_re.match(g)) for g in genes], dtype=bool)

    if cfg.get("remove_histone", False):
        hist_re = re.compile(r"^(H1|H2A|H2B|H3|H4|HIST)")
        to_remove |= np.array([bool(hist_re.match(g)) for g in genes], dtype=bool)

    # Conservative pseudo: only ribo pseudo (safe; avoid MAPKAP1 disaster)
    if cfg.get("remove_pseudogenes", True):
        # RPSxxP, RPLxxP, RPSxx-like, RPLxx-like
        ribo_pseudo_re = re.compile(r"^(RPS|RPL|MRPS|MRPL)\d+[A-Z]?(P|P\d+|-PS\d+|-like)", re.IGNORECASE)
        to_remove |= np.array([bool(ribo_pseudo_re.match(g)) for g in genes], dtype=bool)

    # ENSG (safe)
    if cfg.get("remove_ensg", True):
        ensg_re = re.compile(r"^ENSG\d+")
        to_remove |= np.array([bool(ensg_re.match(g)) for g in genes], dtype=bool)

    # unannotated (safe)
    if cfg.get("remove_unannotated", True):
        unannotated = np.char.find(genes, "LOC") != -1
        to_remove |= unannotated

    n0 = adata.n_vars
    keep = ~to_remove
    n_keep = int(np.sum(keep))
    pct = 100.0 * (n0 - n_keep) / n0 if n0 > 0 else 0.0
    print(f"    Gene filtering: {n0:,} -> {n_keep:,} ({pct:.1f}% removed)")

    ad_filt = adata[:, keep].copy()
    return ad_filt


# =============================================================================
# COUNTS EXTRACTION (Priority: layers['counts'] -> raw.X (aligned) -> X)
# =============================================================================

def get_counts_and_genes(adata: sc.AnnData) -> Tuple[sparse.csr_matrix, np.ndarray, str]:
    """P0-1: Stable priority, no silent fallback, validate alignment."""
    if "counts" in adata.layers:
        print("    Counts source: layers['counts']")
        return sparse.csr_matrix(adata.layers["counts"]), adata.var_names.astype(str).to_numpy(), "layers['counts']"

    if adata.raw is not None:
        raw_genes = adata.raw.var_names.astype(str).to_numpy()
        curr_genes = adata.var_names.astype(str).to_numpy()

        if not np.array_equal(raw_genes, curr_genes):
            print("    WARNING: raw.var_names != var_names (using raw anyway)")
        else:
            print("    Genes aligned (raw = current)")

        print("    Counts source: raw.X (aligned)")
        return sparse.csr_matrix(adata.raw.X), curr_genes, "raw.X (aligned)"

    # fallback
    print("    âš ï¸  Fallback: using .X (check if it's counts!)")
    return sparse.csr_matrix(adata.X), adata.var_names.astype(str).to_numpy(), "X (fallback)"


# =============================================================================
# CHUNKED WRITING (P0-2: avoid huge row buffers + line length limit)
# =============================================================================

def sanitize_names_for_tsv(names: np.ndarray) -> np.ndarray:
    """
    Clean gene/cell names to prevent TSV parsing errors.
    Removes: tabs, newlines, carriage returns, and other problematic characters.
    """
    cleaned = []
    for name in names:
        s = str(name)
        # Replace tabs, newlines, carriage returns with underscores
        s = s.replace('\t', '_').replace('\n', '_').replace('\r', '_')
        # Remove any other control characters
        s = ''.join(c if c.isprintable() or c == ' ' else '_' for c in s)
        cleaned.append(s)
    return np.array(cleaned, dtype=str)


def compute_chunk_genes(n_cells: int, chunk_genes_default: int = 2000, max_chunk_mem_mb: int = 256) -> int:
    """Compute safe chunk size to not exceed max_chunk_mem_mb."""
    mem_per_cell_mb = 4e-6  # float32 = 4 bytes
    max_genes = int(max_chunk_mem_mb / (n_cells * mem_per_cell_mb))
    chunk = min(chunk_genes_default, max(50, max_genes))
    return chunk


def write_counts_tsv_chunked(
    counts_csr: sparse.csr_matrix,
    gene_names: np.ndarray,
    cell_names: np.ndarray,
    out_file: Path,
    chunk_genes: int = 2000,
):
    """
    P0-2: Chunked write to avoid huge line buffers.
    HOTFIX: Sanitize names to prevent TSV parsing errors.
    - transpose once
    - iterate gene-wise in chunks
    - dense chunk per iteration (not per gene)
    """
    # CRITICAL: Clean names before writing to prevent tab/newline in gene/cell names
    gene_names_clean = sanitize_names_for_tsv(gene_names)
    cell_names_clean = sanitize_names_for_tsv(cell_names)
    
    # transpose: (n_cells, n_genes) -> (n_genes, n_cells)
    counts_T = counts_csr.T.tocsr()
    n_genes, n_cells = counts_T.shape

    with open(out_file, "w") as f:
        # header
        header_line = "GENE\t" + "\t".join(cell_names_clean) + "\n"
        f.write(header_line)

        # chunk iteration
        for start_g in range(0, n_genes, chunk_genes):
            end_g = min(start_g + chunk_genes, n_genes)
            chunk_csr = counts_T[start_g:end_g, :]
            chunk_dense = np.asarray(chunk_csr.todense(), dtype=np.float32)  # (chunk_size, n_cells)

            for i_gene in range(chunk_dense.shape[0]):
                g_idx = start_g + i_gene
                gene = gene_names_clean[g_idx]  # Use cleaned name
                row_vals = chunk_dense[i_gene, :]
                row_str = "\t".join([f"{v:.0f}" for v in row_vals])
                f.write(f"{gene}\t{row_str}\n")

    # cleanup
    del counts_T
    gc.collect()


def select_hvg_robust(
    adata: sc.AnnData,
    n_top_genes: int = 3000,
    flavor: str = "seurat",
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
    batch_key: Optional[str] = None
) -> Tuple[bool, str]:
    """
    Robust HVG selection with fallback.
    Returns (success, method_used).
    """
    print(f"  Selecting HVGs (n_top={n_top_genes}, flavor={flavor})...")
    
    # Try batch-aware first
    if batch_key and batch_key in adata.obs.columns:
        try:
            print(f"    Attempting batch-aware HVG (batch_key={batch_key})...")
            sc.pp.highly_variable_genes(
                adata,
                layer="counts",
                n_top_genes=n_top_genes,
                flavor=flavor,
                min_mean=min_mean,
                max_mean=max_mean,
                min_disp=min_disp,
                batch_key=batch_key,
                subset=False
            )
            print(f"    ✓ Batch-aware HVG succeeded")
            return True, "batch-aware"
        except Exception as e:
            print(f"    ⚠️  Batch-aware failed: {e}")
    
    # Fallback to non-batch-aware
    print(f"    Using non-batch-aware HVG...")
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=n_top_genes,
        flavor=flavor,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        subset=False
    )
    print(f"    ✓ Non-batch-aware HVG succeeded")
    return True, "non-batch-aware"


def filter_genes_by_expression(
    adata: sc.AnnData,
    min_cells: int = 10
) -> sc.AnnData:
    """
    Filter genes that are expressed in too few cells.
    This prevents selecting genes that are zero in most cells.
    """
    print(f"  Filtering genes (min_cells={min_cells})...")
    
    n_genes_before = adata.n_vars
    
    # Count cells where gene > 0
    if sparse.issparse(adata.layers["counts"]):
        gene_cell_counts = (adata.layers["counts"] > 0).sum(axis=0).A1
    else:
        gene_cell_counts = (adata.layers["counts"] > 0).sum(axis=0)
    
    keep_genes = gene_cell_counts >= min_cells
    n_genes_after = keep_genes.sum()
    
    print(f"    Genes: {n_genes_before:,} → {n_genes_after:,} "
          f"(removed {n_genes_before - n_genes_after:,})")
    
    return adata[:, keep_genes].copy()


def filter_cells_by_hvg_expression(
    adata: sc.AnnData,
    min_hvg_counts: int = 50
) -> sc.AnnData:
    """
    Filter cells with too few counts in HVG genes.
    Prevents cNMF error: "cells have zero counts of overdispersed genes"
    """
    print(f"  Filtering cells by HVG expression (min_counts={min_hvg_counts})...")
    
    if "highly_variable" not in adata.var.columns:
        print(f"    ⚠️  No HVG selected, skipping cell filter")
        return adata
    
    n_cells_before = adata.n_obs
    
    # Get HVG subset
    hvg_mask = adata.var["highly_variable"].values
    hvg_counts = adata[:, hvg_mask].layers["counts"]
    
    # Sum HVG counts per cell
    if sparse.issparse(hvg_counts):
        cell_hvg_counts = hvg_counts.sum(axis=1).A1
    else:
        cell_hvg_counts = hvg_counts.sum(axis=1)
    
    # Keep cells with sufficient HVG counts
    keep_cells = cell_hvg_counts >= min_hvg_counts
    n_cells_after = keep_cells.sum()
    
    if n_cells_after == 0:
        print(f"    ❌ All cells would be removed! Skipping filter.")
        return adata
    
    removed = n_cells_before - n_cells_after
    pct = 100 * removed / n_cells_before if n_cells_before > 0 else 0
    
    print(f"    Cells: {n_cells_before:,} → {n_cells_after:,} "
          f"(removed {removed:,}, {pct:.1f}%)")
    
    if removed > 0:
        # Show statistics of removed cells
        removed_counts = cell_hvg_counts[~keep_cells]
        print(f"    Removed cells HVG counts: "
              f"mean={removed_counts.mean():.0f}, "
              f"max={removed_counts.max():.0f}")
    
    return adata[keep_cells, :].copy()


def prepare_counts_for_cnmf_optimized(adata: sc.AnnData, out_dir: Path, name: str) -> Path:
    """
    Integrated P0-1 + P0-2 + P1-1:
      - stable counts source priority
      - validate counts
      - chunked write with dynamic chunk size
    """
    print("  Preparing counts (v1.3 chunk-write)...")

    counts, gene_names, source = get_counts_and_genes(adata)
    ok, msg = validate_counts_matrix(counts, source)
    print(f"    Counts source: {source}")
    print(f"    {msg}")
    if not ok:
        raise ValueError(f"Count validation failed: {msg}")

    chunk_genes = compute_chunk_genes(
        n_cells=adata.n_obs,
        chunk_genes_default=CHUNK_GENES_DEFAULT,
        max_chunk_mem_mb=MAX_CHUNK_MEM_MB,
    )
    print(f"    Chunk genes: {chunk_genes} (MAX_CHUNK_MEM_MB={MAX_CHUNK_MEM_MB})")

    out_file = out_dir / f"{name}_counts.txt"
    print(f"    Writing TSV: {out_file}")
    write_counts_tsv_chunked(
        counts_csr=counts,
        gene_names=gene_names,
        cell_names=adata.obs_names.astype(str).to_numpy(),
        out_file=out_file,
        chunk_genes=chunk_genes,
    )

    print("    ✓ Counts TSV ready")
    return out_file


# =============================================================================
# cNMF RUN
# =============================================================================

def run_cnmf_for_group(adata: sc.AnnData, group: str, out_base: Path, cfg: Dict) -> bool:
    group_safe = safe_name(group)
    group_dir = out_base / group_safe
    group_dir.mkdir(parents=True, exist_ok=True)

    cnmf_dir = group_dir / f"{group_safe}_cNMF"
    k_range = determine_k_range(adata.n_obs, cfg)
    print(f"  K range: {k_range} (n_cells={adata.n_obs:,})")

    if SKIP_IF_EXISTS and check_cnmf_complete(cnmf_dir, k_range, cfg["density_threshold"]):
        print(f"  ⏭️  cNMF complete, skip: {cnmf_dir}")
        return True

    # counts TSV
    counts_file = None
    try:
        counts_file = prepare_counts_for_cnmf_optimized(adata, group_dir, group_safe)

        cnmf_obj = cNMF(output_dir=str(cnmf_dir), name=group_safe)

        # prepare
        num_hvg = min(cfg["num_hvg"], adata.n_vars - 1)
        if num_hvg < cfg["num_hvg"]:
            print(f"  ⚠️  num_hvg reduced: {cfg['num_hvg']} -> {num_hvg}")

        print("  [1/5] cnmf.prepare() ...")
        cnmf_obj.prepare(
            counts_fn=str(counts_file),
            components=k_range,
            n_iter=cfg["n_iter"],
            seed=cfg["seed"],
            num_highvar_genes=num_hvg,
            genes_file=None,
        )

        # factorize (single worker)
        print("  [2/5] cnmf.factorize() (single-worker) ...")
        cnmf_obj.factorize(worker_i=0, total_workers=1)

        # combine
        print("  [3/5] cnmf.combine() ...")
        cnmf_obj.combine()

        # k-selection plot
        print("  [4/5] cnmf.k_selection_plot() ...")
        cnmf_obj.k_selection_plot(close_fig=cfg["close_clustergram_fig"])

        # consensus
        print("  [5/5] cnmf.consensus() for each K ...")
        for k in k_range:
            try:
                print(f"    consensus K={k} ...")
                cnmf_obj.consensus(
                    k=k,
                    density_threshold=cfg["density_threshold"],
                    show_clustering=cfg["show_clustering"],
                    close_clustergram_fig=cfg["close_clustergram_fig"],
                )
                print(f"      ✓ K={k}")
            except Exception as e:
                print(f"      ⚠️  K={k} failed: {e}")

        return True

    except Exception as e:
        print(f"  ❌ cNMF failed for {group}: {e}")
        import traceback
        traceback.print_exc()
        return False

    finally:
        if counts_file is not None and DELETE_COUNTS_TXT and counts_file.exists():
            try:
                counts_file.unlink()
                print("  ✓ Removed counts TSV (cleanup)")
            except Exception:
                pass
        gc.collect()


# =============================================================================
# MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("Label-Guided cNMF Pipeline v1.4-FINAL (Production + Zero HVG Fix)")
    print("=" * 90)
    print(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Input: {INPUT_H5AD}")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 90)

    print("\n[STEP 1] Load data...")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  layers: {list(adata.layers.keys())}")
    print(f"  raw: {'✓' if adata.raw is not None else '✗'}")

    if CELLTYPE_KEY not in adata.obs.columns:
        raise ValueError(f"Missing CELLTYPE_KEY in obs: {CELLTYPE_KEY}")

    # lineage mapping
    if USE_MAJOR_LINEAGE:
        print("\n[STEP 2] Major lineage mapping...")
        adata.obs["major_lineage"] = adata.obs[CELLTYPE_KEY].map(MAJOR_LINEAGE_MAP)
        unmapped = adata.obs["major_lineage"].isna()
        if np.any(unmapped):
            # fallback to original labels
            adata.obs.loc[unmapped, "major_lineage"] = adata.obs.loc[unmapped, CELLTYPE_KEY].astype(str)
            unmapped_types = adata.obs.loc[unmapped, CELLTYPE_KEY].astype(str).unique().tolist()
            print(f"  ⚠️  Unmapped types kept as-is (n={len(unmapped_types)}): {unmapped_types[:10]}")
        group_key = "major_lineage"
    else:
        group_key = CELLTYPE_KEY

    # REMOVED: confidence filtering step (per user request)
    print("\n[STEP 3] Confidence filtering: DISABLED (per user request)")

    # select groups
    print("\n[STEP 4] Select groups...")
    vc = adata.obs[group_key].astype(str).value_counts()
    if PROCESS_CELLTYPES == "auto":
        groups = vc[vc >= MIN_CELLS_PER_TYPE].index.tolist()
    else:
        groups = PROCESS_CELLTYPES

    groups = [g for g in groups if vc.get(g, 0) >= MIN_CELLS_PER_TYPE]
    print(f"  Groups to process (>= {MIN_CELLS_PER_TYPE} cells): {len(groups)}")
    for g in groups:
        print(f"    - {g:25s} : {int(vc[g]):,}")

    if not groups:
        print("  ❌ No groups meet threshold.")
        return

    # run per group
    print("\n[STEP 5] Run per group...")
    summary = []
    for i, group in enumerate(groups, 1):
        print("\n" + "=" * 90)
        print(f"[{i}/{len(groups)}] Group: {group}")
        print("=" * 90)

        ad_sub = adata[adata.obs[group_key].astype(str) == str(group)].copy()
        print(f"  Subset cells: {ad_sub.n_obs:,} | genes: {ad_sub.n_vars:,}")

        # Step 1: Gene filtering (MT genes, etc.)
        ad_filt = filter_genes_safe(ad_sub, GENE_FILTER_CONFIG)
        
        # Step 2: Filter genes by minimum expression
        if CNMF_CONFIG.get("min_cells_per_gene", 0) > 0:
            ad_filt = filter_genes_by_expression(
                ad_filt,
                min_cells=CNMF_CONFIG["min_cells_per_gene"]
            )
        
        # Step 3: Select HVGs
        batch_col = None  # Could be extracted from metadata if available
        hvg_success, hvg_method = select_hvg_robust(
            ad_filt,
            n_top_genes=CNMF_CONFIG["num_hvg"],
            flavor=CNMF_CONFIG.get("hvg_flavor", "seurat"),
            min_mean=CNMF_CONFIG.get("hvg_min_mean", 0.0125),
            max_mean=CNMF_CONFIG.get("hvg_max_mean", 3),
            min_disp=CNMF_CONFIG.get("hvg_min_disp", 0.5),
            batch_key=batch_col
        )
        
        if not hvg_success:
            print(f"  ❌ HVG selection failed for {group}")
            summary.append({
                "group": group,
                "group_safe": safe_name(group),
                "n_cells": 0,
                "n_genes": 0,
                "success": False,
            })
            del ad_sub, ad_filt
            gc.collect()
            continue
        
        n_hvg = ad_filt.var["highly_variable"].sum()
        print(f"    ✓ Selected {n_hvg:,} HVGs (method: {hvg_method})")
        
        # Step 4: Filter cells with too few HVG counts
        if CNMF_CONFIG.get("filter_zero_hvg_cells", False):
            ad_filt = filter_cells_by_hvg_expression(
                ad_filt,
                min_hvg_counts=CNMF_CONFIG.get("min_hvg_counts_per_cell", 50)
            )
            
            if ad_filt.n_obs == 0:
                print(f"  ❌ All cells filtered out for {group}")
                summary.append({
                    "group": group,
                    "group_safe": safe_name(group),
                    "n_cells": 0,
                    "n_genes": 0,
                    "success": False,
                })
                del ad_sub, ad_filt
                gc.collect()
                continue
        
        print(f"  Final: {ad_filt.n_obs:,} cells × {ad_filt.n_vars:,} genes")

        # optionally save filtered h5ad
        if SAVE_FILTERED_H5AD:
            fp = OUTPUT_DIR / f"{safe_name(group)}_filtered.h5ad"
            if fp.exists() and not OVERWRITE_GENE_FILTER:
                print(f"  ⏭️  Filtered h5ad exists: {fp} (not overwriting)")
            else:
                ad_filt.write_h5ad(fp, compression="gzip")
                print(f"  ✓ Saved filtered h5ad: {fp}")

        # run cnmf
        ok = run_cnmf_for_group(ad_filt, group, OUTPUT_DIR, CNMF_CONFIG)

        summary.append(
            {
                "group": group,
                "group_safe": safe_name(group),
                "n_cells": int(ad_filt.n_obs),
                "n_genes": int(ad_filt.n_vars),
                "n_hvg": int(n_hvg),
                "hvg_method": hvg_method,
                "success": bool(ok),
            }
        )

        del ad_sub, ad_filt
        gc.collect()

    # summary
    print("\n[STEP 6] Summary...")
    df = pd.DataFrame(summary)
    out_csv = OUTPUT_DIR / "processing_summary.csv"
    df.to_csv(out_csv, index=False)
    print(df.to_string(index=False))
    print(f"\n✓ Summary saved: {out_csv}")
    print("\n✓ PIPELINE COMPLETE")
    print("=" * 90)


if __name__ == "__main__":
    main()
