#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cNMF Helper Module v1.1
=======================
Standalone helper module for cNMF (consensus Non-negative Matrix Factorization)
gene expression program analysis in scRNA-seq workflows.

Consolidates all utility, preparation, execution, post-processing, and
visualization functions from the batch production pipeline.

Usage:
    from cnmf_helper_20260419_v1_1 import *
    # or import specific functions

Author: r2end
Date: 2026-04-23
Version: v1.1

Changelog v1.1 (2026-04-23)
-----------------------------
P0-1 [BUGFIX] _worker_fn: total_workers was hardcoded to 1 in multiprocessing
     mode; each worker was processing the full task instead of its shard.
     Fixed: pass n_workers as total_workers argument to factorize().

P0-2 [BUGFIX] add_usage_to_adata: np.nanargmax raises ValueError on all-NaN
     rows. Fixed: pre-mask all-NaN rows before calling nanargmax.
     Also changed gep_dominant output from 0-based float to 1-based Int64
     (consistent with GEP_1...GEP_k naming convention elsewhere).

P1-3 [BUGFIX] load_gep_matrix: heuristic shape truncation silently corrupted
     data when matrix dimensions did not match K. Fixed: raise ValueError
     on ambiguous shape instead of slicing.

P1-4 [IMPROVE] prepare_cnmf_inputs: replaced adata.copy() (full object copy)
     with view-based sc.AnnData() constructor for both HVG counts and TP10K,
     reducing peak memory. Docstring now explicitly documents adata.raw mutation.

P1-5 [BUGFIX] plot_clustergram: was calling pdist(dist_mat) which computes
     distances between distance-profile vectors (double metric). Fixed:
     use squareform(dist_mat) -> linkage() for standard precomputed-distance
     hierarchical clustering.

P1-6 [IMPROVE] select_hvg_robust: expanded from 2-tier to 3-tier fallback:
     batch-aware seurat_v3 -> non-batch seurat_v3 -> non-batch cell_ranger.
     The cell_ranger flavor is most tolerant of sparse/edge-case data.

P2-7 [IMPROVE] setup_logger: refactored from root-logger manipulation
     (logging.basicConfig + clear root handlers) to named logger with
     propagate=False. Non-invasive when used as an imported library module.
"""

import os
import sys
import gc
import re
import json
import time
import warnings
import logging
import importlib.util
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from datetime import datetime
from multiprocessing import Process

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import entropy

warnings.filterwarnings('ignore')

# Force single-threaded NMF (must be set before cnmf import)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

try:
    from cnmf import cNMF
    CNMF_AVAILABLE = True
except ImportError:
    CNMF_AVAILABLE = False
    print("[WARN] cNMF not installed. Install: pip install cnmf")

_GENE_EXCLUSION_HELPER_PATH = Path(__file__).with_name('program_gene_exclusion_helper_20260505_v1.py')
_gene_exclusion_spec = importlib.util.spec_from_file_location(
    'program_gene_exclusion_helper_20260505_v1',
    _GENE_EXCLUSION_HELPER_PATH,
)
if _gene_exclusion_spec is None or _gene_exclusion_spec.loader is None:
    raise ImportError(f"Cannot load gene exclusion helper: {_GENE_EXCLUSION_HELPER_PATH}")
_gene_exclusion_module = importlib.util.module_from_spec(_gene_exclusion_spec)
sys.modules[_gene_exclusion_spec.name] = _gene_exclusion_module
_gene_exclusion_spec.loader.exec_module(_gene_exclusion_module)
apply_gene_exclusion_to_adata = _gene_exclusion_module.apply_gene_exclusion_to_adata

# ============================================================================
# DEFAULT CONFIGURATION
# ============================================================================

# K range presets (override per analysis)
# Typical GEP counts in scRNA-seq: 4-8 programs per cell type
K_RANGE_LARGE  = [4, 5, 6, 7, 8]   # >= 50k cells
K_RANGE_MEDIUM = [4, 5, 6, 7, 8]   # 10k - 50k cells
K_RANGE_SMALL  = [3, 4, 5, 6]      # < 10k cells

DEFAULT_CNMF_CONFIG = {
    'n_iter'             : 100,
    'seed'               : 42,
    'num_hvg'            : 3000,
    'hvg_flavor'         : 'seurat_v3',
    'density_threshold'  : 0.1,
    'show_clustering'    : True,
    'close_clustergram'  : True,
    'n_workers'          : 1,       # 1 = single-worker (safest); >1 = multiprocessing
    'exclude_technical'  : False,   # True = remove cell_cycle / stress IEGs from HVG
    'gene_exclusion_config': None,
}

DEFAULT_VIZ_CONFIG = {
    'dpi'                  : 300,
    'figure_format'        : 'png',
    'local_density_bins'   : 50,
    'local_density_k_nbrs' : 20,
    'cluster_method'       : 'average',
    'cluster_metric'       : 'euclidean',
}

DEFAULT_K_SELECT_CONFIG = {
    'elbow_sensitivity'    : 0.05,
    'optimal_sparsity_min' : 0.3,
    'optimal_sparsity_max' : 0.7,
}

# Technical gene lists
TECHNICAL_GENES = {
    'cell_cycle': [
        'MKI67', 'TOP2A', 'PCNA', 'CCNA2', 'CCNB1', 'CCNB2',
        'CDK1', 'AURKA', 'AURKB', 'CENPF', 'CENPE',
    ],
    'stress_ieg': [
        'FOS', 'JUN', 'JUNB', 'JUND', 'EGR1', 'EGR2', 'EGR3',
        'ATF3', 'DUSP1', 'HSPA1A', 'HSPA1B', 'HSP90AA1',
    ],
}

# Candidate column name lists (inferred in order)
CELLTYPE_COL_CANDIDATES = [
    'cell_type_final_l3', 'cell_type_final_l2', 'cell_type_final_l1',
    'scanvi_predictions', 'cell_type_scanvi_filt', 'scanvi_labels',
    'cell_type', 'celltype',
]
BATCH_COL_CANDIDATES = [
    'dataset', 'batch', 'sample_id', 'orig.ident', 'Sample',
]

# ============================================================================
# LOGGING
# ============================================================================

def setup_logger(name: str = 'cnmf_helper',
                 log_file: Optional[Path] = None,
                 level: int = logging.INFO) -> logging.Logger:
    """
    Create or retrieve a named logger writing to stdout (+ optional file).

    Uses a named logger (not root) so it does not interfere with the caller's
    logging configuration. Calling this multiple times with the same name is
    safe — handlers are not duplicated.

    Parameters
    ----------
    name     : Logger name (default 'cnmf_helper').
    log_file : If given, also write to this file.
    level    : Logging level (default INFO).

    Returns
    -------
    logging.Logger
    """
    log = logging.getLogger(name)
    log.setLevel(level)

    # Avoid duplicate handlers on repeated calls
    if log.handlers:
        return log

    fmt     = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    sh      = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    log.addHandler(sh)

    if log_file is not None:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        log.addHandler(fh)

    log.propagate = False   # don't bubble to root logger
    return log

logger = setup_logger()


# ============================================================================
# SECTION 1 — UTILITY FUNCTIONS
# ============================================================================

def safe_mkdir(directory: Path, description: str = '') -> bool:
    """
    Create directory (+ parents) and verify existence.

    Returns True on success, False on failure.
    """
    try:
        directory.mkdir(exist_ok=True, parents=True)
        if not directory.is_dir():
            logger.error(f"[WARN] Failed to create {description or directory}")
            return False
        return True
    except Exception as e:
        logger.error(f"[WARN] Cannot create {description or directory}: {e}")
        return False


def safe_name(x: str, max_len: int = 180) -> str:
    """
    Sanitize a string for use as directory / file name.
    Replaces path-unsafe characters with underscores.
    """
    s = str(x)
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t', "'", '"']:
        s = s.replace(ch, '_')
    return s[:max_len]


def first_glob_match(base: Path,
                     patterns: List[str]) -> Optional[Path]:
    """
    Return the first file in *base* matching any glob pattern in *patterns*.
    Returns None if no match found.
    """
    if base is None or not isinstance(base, Path):
        return None
    for pat in patterns:
        try:
            matches = sorted(base.glob(pat))
        except Exception:
            continue
        if matches:
            return matches[0]
    return None


def determine_k_range(n_cells: int) -> List[int]:
    """
    Select K range preset based on cell count.

    >= 50 000  → K_RANGE_LARGE
    10 000-50 000 → K_RANGE_MEDIUM
    <  10 000  → K_RANGE_SMALL
    """
    if n_cells >= 50_000:
        return list(K_RANGE_LARGE)
    elif n_cells >= 10_000:
        return list(K_RANGE_MEDIUM)
    else:
        return list(K_RANGE_SMALL)


def infer_column(adata: sc.AnnData,
                 candidates: List[str]) -> Optional[str]:
    """
    Return the first candidate column name present in adata.obs.
    Returns None if none found.
    """
    for col in candidates:
        if col in adata.obs.columns:
            return col
    return None


def detect_technical_genes(var_names: pd.Index) -> Dict[str, List[str]]:
    """
    Identify technical genes (cell_cycle, stress/IEG, ribosomal, mitochondrial)
    present in *var_names*.

    Returns dict: category -> list of matching gene names.
    """
    detected: Dict[str, List[str]] = {
        'cell_cycle'   : [],
        'stress_ieg'   : [],
        'ribosomal'    : [],
        'mitochondrial': [],
    }

    var_upper = [str(g).upper() for g in var_names]
    var_list  = list(var_names)

    for gene in TECHNICAL_GENES['cell_cycle']:
        if gene in var_upper:
            detected['cell_cycle'].append(var_list[var_upper.index(gene)])

    for gene in TECHNICAL_GENES['stress_ieg']:
        if gene in var_upper:
            detected['stress_ieg'].append(var_list[var_upper.index(gene)])

    for gene in var_list:
        gu = str(gene).upper()
        if gu.startswith('RPL') or gu.startswith('RPS'):
            detected['ribosomal'].append(gene)
        elif gu.startswith('MT-'):
            detected['mitochondrial'].append(gene)

    return detected


# ============================================================================
# SECTION 2 — HVG SELECTION
# ============================================================================

def select_hvg_robust(
    adata          : sc.AnnData,
    num_hvg        : int             = 3000,
    batch_key      : Optional[str]   = None,
    hvg_flavor     : str             = 'seurat_v3',
    exclude_tech   : bool            = False,
) -> Tuple[List[str], str, Dict[str, List[str]]]:
    """
    Robust HVG selection with batch-aware → non-batch-aware fallback.

    Always reads from layers['counts'].
    Caps num_hvg to adata.n_vars - 1 automatically.

    Parameters
    ----------
    adata       : AnnData with layers['counts']
    num_hvg     : number of HVGs to select
    batch_key   : obs column for batch-aware selection (None = non-batch)
    hvg_flavor  : 'seurat_v3' | 'cell_ranger' | 'seurat'
    exclude_tech: if True, remove cell_cycle / stress IEGs from final HVG list

    Returns
    -------
    hvg_genes   : list of selected HVG gene names
    method_used : 'batch-aware' | 'non-batch-aware'
    tech_in_hvg : dict of detected technical genes per category
    """
    if 'counts' not in adata.layers:
        raise ValueError("[WARN] layers['counts'] not found. Cannot select HVGs.")

    num_hvg = min(num_hvg, adata.n_vars - 1)
    method_used = 'unknown'

    logger.info(f"[INFO] Selecting {num_hvg} HVGs (flavor={hvg_flavor})...")

    # 3-tier fallback: batch-aware → non-batch-aware (same flavor) → cell_ranger
    # Tier 1: batch-aware seurat_v3
    try:
        if batch_key and batch_key in adata.obs.columns:
            logger.info(f"  Tier-1: batch-aware HVG (batch_key={batch_key}, flavor={hvg_flavor})...")
            sc.pp.highly_variable_genes(
                adata,
                layer       = 'counts',
                n_top_genes = num_hvg,
                flavor      = hvg_flavor,
                batch_key   = batch_key,
                subset      = False,
            )
            method_used = f'batch-aware-{hvg_flavor}'
            logger.info(f"  [OK] Tier-1 succeeded")
        else:
            raise ValueError("No batch_key — skip Tier-1")
    except Exception as e1:
        logger.warning(f"  Tier-1 failed ({e1}), trying Tier-2...")
        # Tier 2: non-batch-aware, same flavor
        try:
            sc.pp.highly_variable_genes(
                adata,
                layer       = 'counts',
                n_top_genes = num_hvg,
                flavor      = hvg_flavor,
                batch_key   = None,
                subset      = False,
            )
            method_used = f'non-batch-{hvg_flavor}'
            logger.info(f"  [OK] Tier-2 succeeded")
        except Exception as e2:
            logger.warning(f"  Tier-2 failed ({e2}), trying Tier-3 (cell_ranger)...")
            # Tier 3: cell_ranger flavor — most tolerant of sparse / edge-case data
            try:
                sc.pp.highly_variable_genes(
                    adata,
                    layer       = 'counts',
                    n_top_genes = num_hvg,
                    flavor      = 'cell_ranger',
                    batch_key   = None,
                    subset      = False,
                )
                method_used = 'non-batch-cell_ranger'
                logger.info(f"  [OK] Tier-3 (cell_ranger) succeeded")
            except Exception as e3:
                raise RuntimeError(
                    f"All 3 HVG tiers failed.\n"
                    f"  Tier-1: {e1}\n  Tier-2: {e2}\n  Tier-3: {e3}"
                ) from e3

    hvg_mask  = adata.var['highly_variable'].values
    hvg_genes = adata.var_names[hvg_mask].tolist()

    logger.info(f"  Selected {len(hvg_genes)} HVGs (method={method_used})")

    tech_in_hvg = detect_technical_genes(adata.var_names[hvg_mask])
    n_tech = sum(len(v) for v in tech_in_hvg.values())
    if n_tech > 0:
        for cat, genes in tech_in_hvg.items():
            if genes:
                logger.warning(f"  [WARN] Technical genes in HVG ({cat}): {len(genes)}")

    if exclude_tech:
        tech_flat = set(g for gl in tech_in_hvg.values() for g in gl)
        before = len(hvg_genes)
        hvg_genes = [g for g in hvg_genes if g not in tech_flat]
        logger.info(f"  Removed {before - len(hvg_genes)} technical genes; {len(hvg_genes)} remain")

    adata.uns['hvg_method'] = method_used

    return hvg_genes, method_used, tech_in_hvg


# ============================================================================
# SECTION 3 — INPUT PREPARATION
# ============================================================================

def prepare_cnmf_inputs(
    adata    : sc.AnnData,
    hvg_genes: List[str],
    output_dir: Path,
) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Prepare three cNMF input files:
      1. hvg_counts.h5ad  — HVG-subset raw counts (cNMF counts_fn)
      2. tp10k.h5ad       — full-gene TP10K matrix  (cNMF tpm_fn)
      3. hvg_genes.txt    — HVG gene list            (cNMF genes_file)

    SIDE EFFECT: adata.raw is overwritten with full-gene counts.
    Caller should be aware that adata is mutated in place.

    Memory strategy:
      - .raw shares adata.layers['counts'] (no copy)
      - HVG counts: view-based AnnData constructor (no full adata.copy())
      - TP10K: minimal AnnData from counts.copy(), written then deleted
      - Objects are deleted in order: HVG first, TP10K second

    Returns (hvg_counts_path, tp10k_path, hvg_txt_path)
    or      (None, None, None) on failure.
    """
    logger.info("[INFO] Preparing cNMF input files...")
    logger.info("  [NOTE] This call overwrites adata.raw.")

    try:
        # --- .raw: share full-gene counts (no copy) ---
        logger.info("  Setting .raw to full-gene counts...")
        adata.raw = sc.AnnData(
            X   = adata.layers['counts'],  # shared reference, no .copy()
            obs = adata.obs.copy(),
            var = adata.var.copy(),
        )
        logger.info(f"  [OK] .raw set ({adata.n_vars} genes)")

        hvg_mask     = adata.var_names.isin(hvg_genes)
        n_hvg_actual = int(hvg_mask.sum())

        # --- HVG gene list (cheap, write first) ---
        hvg_txt_path = output_dir / 'cnmf_input_hvg_genes.txt'
        with open(hvg_txt_path, 'w') as fh:
            for gene in adata.var_names[hvg_mask]:
                fh.write(f'{gene}\n')
        logger.info(f"  [OK] Saved HVG list: {hvg_txt_path.name} ({n_hvg_actual} genes)")

        # --- HVG counts: view-based constructor, avoids full adata.copy() ---
        logger.info(f"  Creating HVG counts subset ({n_hvg_actual} genes)...")
        adata_hvg = sc.AnnData(
            X   = adata.layers['counts'][:, hvg_mask],  # CSR view, minimal memory
            obs = adata.obs.copy(),
            var = adata.var[hvg_mask].copy(),
        )
        hvg_counts_path = output_dir / 'cnmf_input_hvg_counts.h5ad'
        adata_hvg.write_h5ad(hvg_counts_path, compression='gzip')
        logger.info(f"  [OK] Saved HVG counts: {hvg_counts_path.name}")
        del adata_hvg
        gc.collect()

        # --- TP10K (full-gene): minimal constructor, normalize, write, delete ---
        logger.info("  Creating TP10K matrix (full-gene)...")
        adata_tp10k = sc.AnnData(
            X   = adata.layers['counts'].copy(),
            obs = adata.obs.copy(),
            var = adata.var.copy(),
        )
        sc.pp.normalize_total(adata_tp10k, target_sum=1e4)
        tp10k_path = output_dir / 'cnmf_input_tp10k.h5ad'
        adata_tp10k.write_h5ad(tp10k_path, compression='gzip')
        logger.info(f"  [OK] Saved TP10K: {tp10k_path.name}")
        del adata_tp10k
        gc.collect()

        return hvg_counts_path, tp10k_path, hvg_txt_path

    except Exception as e:
        logger.error(f"  [WARN] Input preparation failed: {e}")
        return None, None, None


# ============================================================================
# SECTION 4 — CNMF EXECUTION
# ============================================================================

def _worker_fn(output_dir: str, name: str, worker_i: int, total_workers: int) -> None:
    """Target function for multiprocessing workers.
    
    CRITICAL: total_workers must match the actual pool size so each worker
    processes a distinct shard. Passing total_workers=1 here would cause all
    workers to re-run the full task and produce conflicting outputs.
    """
    try:
        obj = cNMF(output_dir=output_dir, name=name)
        obj.factorize(worker_i=worker_i, total_workers=total_workers)
    except Exception as e:
        print(f"[WARN] Worker {worker_i}/{total_workers} failed: {e}", flush=True)
        sys.exit(1)


def run_factorize(
    output_dir   : str,
    name         : str,
    n_workers    : int = 1,
) -> bool:
    """
    Run cNMF factorization.

    n_workers=1 (default): single-worker mode (safest, recommended).
    n_workers>1: multiprocessing with Process pool.

    Returns True on success.
    """
    if n_workers <= 1:
        logger.info("  Factorizing (single-worker mode)...")
        try:
            obj = cNMF(output_dir=output_dir, name=name)
            obj.factorize(worker_i=0, total_workers=1)
            logger.info("  [OK] Factorization complete")
            return True
        except Exception as e:
            logger.error(f"  [WARN] Factorization failed: {e}")
            return False
    else:
        logger.info(f"  Factorizing ({n_workers} workers)...")
        procs = []
        for wi in range(n_workers):
            p = Process(target=_worker_fn, args=(output_dir, name, wi, n_workers))
            p.start()
            procs.append(p)

        all_ok = True
        for i, p in enumerate(procs):
            p.join()
            if p.exitcode != 0:
                logger.error(f"  [WARN] Worker {i} failed (exitcode={p.exitcode})")
                all_ok = False
            else:
                logger.info(f"  [OK] Worker {i} done")

        return all_ok


def run_cnmf_pipeline(
    hvg_counts_path : Path,
    tp10k_path      : Path,
    hvg_txt_path    : Path,
    k_range         : List[int],
    output_dir      : Path,
    name            : str,
    cnmf_config     : Optional[Dict] = None,
) -> bool:
    """
    Execute the full 4-step cNMF pipeline:
      Step 1: prepare
      Step 2: factorize
      Step 3: combine
      Step 4: consensus (per K)

    Parameters
    ----------
    hvg_counts_path : Path to HVG raw counts h5ad (counts_fn)
    tp10k_path      : Path to full-gene TP10K h5ad (tpm_fn)
    hvg_txt_path    : Path to HVG gene list txt (genes_file)
    k_range         : list of K values to evaluate
    output_dir      : base output directory; cnmf outputs go to output_dir/cnmf_output/
    name            : cNMF run name (sanitized via safe_name before passing)
    cnmf_config     : override dict (falls back to DEFAULT_CNMF_CONFIG)

    Returns True if all 4 steps complete without fatal error.
    """
    if not CNMF_AVAILABLE:
        logger.error("[WARN] cNMF not installed. Cannot run pipeline.")
        return False

    cfg = {**DEFAULT_CNMF_CONFIG, **(cnmf_config or {})}

    cnmf_out = output_dir / 'cnmf_output'
    if not safe_mkdir(cnmf_out, 'cNMF output'):
        return False

    logger.info("=" * 70)
    logger.info(f"[INFO] cNMF Pipeline: {name}")
    logger.info(f"  K range: {k_range}")
    logger.info("=" * 70)

    try:
        cnmf_obj = cNMF(output_dir=str(cnmf_out), name=name)

        # --- Step 1: prepare ---
        logger.info("\nStep 1/4: prepare")
        cnmf_obj.prepare(
            counts_fn   = str(hvg_counts_path),
            tpm_fn      = str(tp10k_path),
            genes_file  = str(hvg_txt_path),
            components  = k_range,
            n_iter      = cfg['n_iter'],
            seed        = cfg['seed'],
        )
        logger.info("  [OK] prepare complete")

        # --- Step 2: factorize ---
        logger.info("\nStep 2/4: factorize")
        t0 = time.time()
        ok = run_factorize(
            output_dir = str(cnmf_out),
            name       = name,
            n_workers  = cfg['n_workers'],
        )
        if not ok:
            logger.error("  [WARN] Factorization failed — aborting pipeline")
            return False
        logger.info(f"  [OK] factorize complete ({(time.time()-t0)/60:.1f} min)")

        # --- Step 3: combine ---
        logger.info("\nStep 3/4: combine")
        cnmf_obj.combine()
        logger.info("  [OK] combine complete")

        # --- Step 4: consensus ---
        logger.info("\nStep 4/4: consensus")
        n_ok = 0
        for k in k_range:
            try:
                cnmf_obj.consensus(
                    k                    = k,
                    density_threshold    = cfg['density_threshold'],
                    show_clustering      = cfg['show_clustering'],
                    close_clustergram_fig= cfg['close_clustergram'],
                )
                logger.info(f"  [OK] K={k}")
                n_ok += 1
            except Exception as e:
                logger.warning(f"  [WARN] K={k} consensus failed: {e}")
            finally:
                plt.close('all')

        logger.info(f"\n[OK] cNMF pipeline done ({n_ok}/{len(k_range)} K values)")
        return n_ok > 0

    except Exception as e:
        logger.error(f"[WARN] cNMF pipeline fatal error: {e}")
        plt.close('all')
        return False


# ============================================================================
# SECTION 5 — GEP / USAGE MATRIX LOADING
# ============================================================================

def find_gep_spectra_file(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
) -> Optional[Path]:
    """
    Locate the spectra score file for a given K.

    Handles both dt_0_1 and dt_0.1 filename conventions.
    """
    return first_glob_match(
        cnmf_run_dir,
        [
            f'{name}.gene_spectra_score.k_{k}.dt_0_1.txt',
            f'{name}.gene_spectra_score.k_{k}.dt_0.1.txt',
            f'{name}.gene_spectra_score.k_{k}.dt_*.txt',
            f'{name}.gene_spectra_score.k_{k}.*.txt',
            f'{name}.spectra.k_{k}.dt_0_1.consensus.txt',
            f'{name}.spectra.k_{k}.dt_0.1.consensus.txt',
            f'{name}.spectra.k_{k}.dt_*.consensus.txt',
            f'{name}.spectra.k_{k}.*consensus.txt',
        ],
    )


def load_gep_score_dataframe(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
) -> Optional[pd.DataFrame]:
    """
    Load spectra scores as a standardized DataFrame with
    rows = GEP_1..GEP_k and columns = gene names.
    """
    spectra_file = find_gep_spectra_file(cnmf_run_dir, name, k)

    if spectra_file is None or not spectra_file.exists():
        logger.warning(f"  [WARN] Spectra file not found (K={k})")
        return None

    logger.info(f"  Loading spectra: {spectra_file.name}")

    try:
        df = pd.read_csv(spectra_file, sep='\t', index_col=0)

        if df.shape[0] == k:
            gep_df = df.copy()
        elif df.shape[1] == k:
            gep_df = df.T.copy()
        else:
            raise ValueError(
                f"Spectra matrix shape {df.shape} is inconsistent with K={k}. "
                f"Expected one dimension == {k}. "
                f"Check that the correct K and file are being used."
            )

        gep_df = gep_df.apply(pd.to_numeric, errors='coerce').fillna(0.0)

        if gep_df.shape[0] != k or gep_df.ndim != 2 or gep_df.shape[0] < 2:
            logger.warning(f"  [WARN] Invalid GEP score DataFrame shape: {gep_df.shape}")
            return None

        gep_df.index = [f'GEP_{i+1}' for i in range(k)]
        gep_df.columns = [str(col) for col in gep_df.columns]

        logger.info(f"  [OK] GEP score DataFrame loaded: {gep_df.shape}")
        return gep_df

    except Exception as e:
        logger.warning(f"  [WARN] Failed to load spectra: {e}")
        return None

def load_gep_matrix(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
) -> Optional[np.ndarray]:
    """
    Load gene spectra score matrix for a given K.

    Handles both dt_0_1 and dt_0.1 filename conventions.
    Returns ndarray of shape (k, n_genes) or None.
    """
    gep_df = load_gep_score_dataframe(cnmf_run_dir, name, k)
    if gep_df is None:
        return None

    mat = np.asarray(gep_df.values, dtype=float)
    if not np.isfinite(mat).all():
        mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
    return mat


def load_usage_matrix(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
) -> Optional[pd.DataFrame]:
    """
    Load per-cell GEP usage matrix for a given K.

    Returns DataFrame (cells × GEPs) with numeric values, or None.
    """
    usage_file = first_glob_match(
        cnmf_run_dir,
        [
            f'{name}.usages.k_{k}.dt_0_1.consensus.txt',
            f'{name}.usages.k_{k}.dt_0.1.consensus.txt',
            f'{name}.usages.k_{k}.dt_*.consensus.txt',
            f'{name}.usages.k_{k}.*consensus.txt',
        ],
    )

    if usage_file is None or not usage_file.exists():
        logger.warning(f"  [WARN] Usage file not found (K={k})")
        return None

    logger.info(f"  Loading usage: {usage_file.name}")

    try:
        df = pd.read_csv(usage_file, sep='\t', index_col=0)
        df = df.apply(pd.to_numeric, errors='coerce').fillna(0.0)
        logger.info(f"  [OK] Usage matrix loaded: {df.shape}")
        return df
    except Exception as e:
        logger.warning(f"  [WARN] Failed to load usage: {e}")
        return None


def load_top_genes(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
    n_top        : int = 50,
) -> Optional[pd.DataFrame]:
    """
    Load top-genes table (*.top_genes.k_{k}.dt_*.txt) for a given K.

    Returns DataFrame or None.
    """
    top_file = first_glob_match(
        cnmf_run_dir,
        [
            f'{name}.top_genes.k_{k}.dt_0_1.txt',
            f'{name}.top_genes.k_{k}.dt_0.1.txt',
            f'{name}.top_genes.k_{k}.dt_*.txt',
        ],
    )

    if top_file is None or not top_file.exists():
        logger.warning(f"  [WARN] Top-genes file not found (K={k})")
        return None

    try:
        df = pd.read_csv(top_file, sep='\t', index_col=0)
        logger.info(f"  [OK] Top-genes loaded: {df.shape}")
        return df.iloc[:n_top] if n_top < len(df) else df
    except Exception as e:
        logger.warning(f"  [WARN] Failed to load top_genes: {e}")
        return None


def add_usage_to_adata(
    adata        : sc.AnnData,
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
    obsm_key     : Optional[str] = None,
) -> bool:
    """
    Attach GEP usage scores to adata.obsm and adata.obs.

    Usage is aligned on cell barcodes (intersection).
    Adds:
      adata.obsm[obsm_key]    — usage matrix (n_cells × k), NaN for missing cells
      adata.obs['gep_dominant'] — dominant GEP index per cell

    Parameters
    ----------
    obsm_key : key name in obsm; defaults to f'X_cnmf_k{k}'
    """
    if obsm_key is None:
        obsm_key = f'X_cnmf_k{k}'

    usage_df = load_usage_matrix(cnmf_run_dir, name, k)
    if usage_df is None:
        return False

    common = adata.obs_names.intersection(usage_df.index)
    if len(common) == 0:
        logger.warning(f"  [WARN] No overlapping barcodes between adata and usage matrix")
        return False

    logger.info(f"  [OK] Aligning usage: {len(common)}/{adata.n_obs} cells")

    # Full matrix with NaN for missing cells
    usage_full = pd.DataFrame(
        np.full((adata.n_obs, usage_df.shape[1]), np.nan),
        index   = adata.obs_names,
        columns = usage_df.columns,
    )
    usage_full.loc[common, :] = usage_df.loc[common, :].values

    adata.obsm[obsm_key] = usage_full.values.astype(float)

    # Dominant GEP per cell (1-based, matching GEP_1...GEP_k naming convention).
    # CRITICAL: np.nanargmax raises on all-NaN rows, so mask those first.
    vals     = usage_full.values
    all_nan  = np.all(np.isnan(vals), axis=1)
    dom_arr  = np.full(vals.shape[0], pd.NA, dtype=object)
    valid    = ~all_nan
    if valid.any():
        dom_arr[valid] = np.nanargmax(vals[valid], axis=1) + 1  # 1-based
    adata.obs[f'gep_dominant_k{k}'] = pd.array(dom_arr, dtype='Int64')

    logger.info(f"  [OK] Usage added: obsm['{obsm_key}'], obs['gep_dominant_k{k}']")
    return True


# ============================================================================
# SECTION 6 — DISTANCE / SIMILARITY MATRICES
# ============================================================================

def compute_gep_distance_matrix(
    gep_matrix : np.ndarray,
    metric     : str = 'euclidean',
) -> Optional[np.ndarray]:
    """
    Compute pairwise GEP distance matrix from spectra score matrix.

    Parameters
    ----------
    gep_matrix : (k, n_genes) array
    metric     : distance metric passed to scipy.spatial.distance.pdist

    Returns
    -------
    (k, k) distance matrix or None.
    """
    if gep_matrix is None or gep_matrix.size == 0:
        return None

    try:
        dist_vec = pdist(gep_matrix, metric=metric)
        if dist_vec.size == 0:
            return None
        dist_mat = squareform(dist_vec)
        return dist_mat
    except Exception as e:
        logger.warning(f"  [WARN] GEP distance failed: {e}")
        return None


def get_gep_distance_matrix(
    cnmf_output_dir : Path,
    name            : str,
    k               : int,
) -> Optional[np.ndarray]:
    """
    Convenience wrapper: load spectra → compute distance matrix.
    """
    cnmf_run_dir = cnmf_output_dir / name
    gep = load_gep_matrix(cnmf_run_dir, name, k)
    if gep is None:
        return None
    return compute_gep_distance_matrix(
        gep,
        metric=DEFAULT_VIZ_CONFIG['cluster_metric'],
    )


# ============================================================================
# SECTION 7 — K-VALUE SELECTION METRICS
# ============================================================================

def compute_k_metrics(
    cnmf_output_dir : Path,
    name            : str,
    k               : int,
) -> Dict[str, Optional[float]]:
    """
    Compute four K-selection metrics for a single K value:

    1. mean_usage_entropy     — Shannon entropy per cell (higher = more diffuse)
    2. usage_sparsity_gini    — Gini coefficient of total GEP usage (higher = more specific)
    3. max_usage_concentration— mean of per-cell max usage value
    4. mean_pairwise_distance — mean inter-GEP Euclidean distance

    Returns dict with keys above; values are float or None if unavailable.
    """
    result: Dict[str, Optional[float]] = {
        'mean_usage_entropy'     : None,
        'usage_sparsity_gini'    : None,
        'max_usage_concentration': None,
        'mean_pairwise_distance' : None,
    }

    cnmf_run_dir = cnmf_output_dir / name

    # --- Usage matrix metrics ---
    usage_df = load_usage_matrix(cnmf_run_dir, name, k)
    if usage_df is not None:
        mat = usage_df.values

        # 1. Entropy
        ents = []
        for i in range(mat.shape[0]):
            u = mat[i, :]
            s = u.sum()
            if s > 0:
                ents.append(float(entropy(u / s)))
        if ents:
            result['mean_usage_entropy'] = float(np.mean(ents))

        # 2. Gini
        total = mat.sum(axis=0)
        if total.sum() > 0:
            sv = np.sort(total)
            n  = len(sv)
            idx = np.arange(1, n + 1)
            gini = (2 * np.sum(idx * sv)) / (n * sv.sum()) - (n + 1) / n
            result['usage_sparsity_gini'] = float(gini)

        # 3. Max concentration
        result['max_usage_concentration'] = float(np.mean(mat.max(axis=1)))

    # --- GEP distance ---
    dist_mat = get_gep_distance_matrix(cnmf_output_dir, name, k)
    if dist_mat is not None:
        tri = dist_mat[np.triu_indices_from(dist_mat, k=1)]
        tri = tri[np.isfinite(tri)]
        if tri.size > 0:
            result['mean_pairwise_distance'] = float(np.mean(tri))

    return result


def calculate_k_stability_metrics(
    cnmf_output_dir : Path,
    name            : str,
    k_range         : List[int],
    output_dir      : Path,
    viz_config      : Optional[Dict] = None,
) -> Dict[int, Dict[str, Optional[float]]]:
    """
    Compute stability metrics for all K values, save JSON, generate summary plot,
    and write K-selection recommendation.

    Returns
    -------
    metrics : {k: {metric_name: value, ...}, ...}
    """
    logger.info("\n" + "=" * 70)
    logger.info("[INFO] K-Value Stability Metrics")
    logger.info("=" * 70)

    metrics: Dict[int, Dict[str, Optional[float]]] = {}

    for k in k_range:
        logger.info(f"\n--- K={k} ---")
        m = compute_k_metrics(cnmf_output_dir, name, k)
        metrics[k] = m

        for key, val in m.items():
            if val is not None:
                logger.info(f"  {key}: {val:.4f}")
            else:
                logger.warning(f"  {key}: N/A")

    # Save JSON
    metrics_file = output_dir / 'k_stability_metrics.json'
    with open(metrics_file, 'w') as fh:
        json.dump(metrics, fh, indent=2)
    logger.info(f"\n[OK] Saved metrics: {metrics_file.name}")

    # K recommendation
    recommendation = generate_k_recommendation(metrics, k_range)
    rec_file = output_dir / 'k_selection_recommendation.json'
    with open(rec_file, 'w') as fh:
        json.dump(recommendation, fh, indent=2)
    logger.info(f"[OK] Saved recommendation: {rec_file.name}")
    logger.info(f"  Recommended K: {recommendation.get('recommended_k')} "
                f"(confidence: {recommendation.get('confidence')})")

    # Summary plot
    try:
        plot_k_metrics_summary(metrics, recommendation, output_dir, viz_config)
    except Exception as e:
        logger.warning(f"[WARN] K metrics plot failed: {e}")

    return metrics


def generate_k_recommendation(
    metrics : Dict[int, Dict[str, Optional[float]]],
    k_range : List[int],
) -> Dict[str, Any]:
    """
    Recommend a K value from stability metrics using three criteria:
      (a) entropy elbow (first point where delta-entropy < threshold)
      (b) distance stabilization
      (c) optimal sparsity Gini (closest to 0.5)

    Returns dict with recommended_k, confidence, reasoning, alternative_k.
    """
    cfg = DEFAULT_K_SELECT_CONFIG

    recommendation: Dict[str, Any] = {
        'recommended_k': None,
        'confidence'   : 'low',
        'reasoning'    : [],
        'alternative_k': [],
    }

    valid_k   = [k for k in k_range if metrics.get(k, {}).get('mean_usage_entropy') is not None]
    if len(valid_k) < 3:
        recommendation['reasoning'].append('Insufficient K values with valid metrics')
        return recommendation

    entropies  = [metrics[k]['mean_usage_entropy']      for k in valid_k]
    distances  = [metrics[k].get('mean_pairwise_distance', 0) or 0 for k in valid_k]
    sparsities = [metrics[k].get('usage_sparsity_gini', 0)    or 0 for k in valid_k]

    candidates = []

    # (a) Entropy elbow
    d_entropy = np.diff(entropies)
    elbow_idx = next((i for i, d in enumerate(d_entropy) if d < cfg['elbow_sensitivity']), None)
    if elbow_idx is not None and elbow_idx > 0:
        ek = valid_k[elbow_idx]
        candidates.append(ek)
        recommendation['reasoning'].append(f'Entropy plateaus at K={ek}')

    # (b) Distance stabilization
    d_dist = np.diff(distances)
    stable_idx = next((i for i, d in enumerate(d_dist) if abs(d) < cfg['elbow_sensitivity']), None)
    if stable_idx is not None and stable_idx > 0:
        dk = valid_k[stable_idx]
        candidates.append(dk)
        recommendation['reasoning'].append(f'GEP separation stabilizes at K={dk}')

    # (c) Gini closest to 0.5
    gini_scores = [1.0 - abs(s - 0.5) for s in sparsities]
    best_gini_idx = int(np.argmax(gini_scores))
    bk = valid_k[best_gini_idx]
    candidates.append(bk)
    recommendation['reasoning'].append(f'Optimal usage sparsity (Gini≈0.5) at K={bk}')

    if candidates:
        from collections import Counter
        counts = Counter(candidates)
        recommended_k = counts.most_common(1)[0][0]
        recommendation['recommended_k'] = recommended_k
        recommendation['confidence']    = 'high' if counts[recommended_k] >= 2 else 'medium'
        recommendation['alternative_k'] = sorted(
            set(candidates) - {recommended_k}
        )
    else:
        # Fallback: middle K
        mid_k = valid_k[len(valid_k) // 2]
        recommendation['recommended_k'] = mid_k
        recommendation['confidence']    = 'low'
        recommendation['reasoning'].append(f'Fallback: middle K={mid_k}')

    return recommendation


# ============================================================================
# SECTION 8 — VISUALIZATIONS
# ============================================================================

def plot_k_metrics_summary(
    metrics        : Dict[int, Dict[str, Optional[float]]],
    recommendation : Dict[str, Any],
    output_dir     : Path,
    viz_config     : Optional[Dict] = None,
) -> None:
    """
    2×2 summary figure:
      [0,0] Usage entropy vs K
      [0,1] Mean pairwise GEP distance vs K
      [1,0] Usage sparsity (Gini) vs K
      [1,1] Text recommendation panel
    """
    cfg = {**DEFAULT_VIZ_CONFIG, **(viz_config or {})}
    k_values = sorted(metrics.keys())

    entropies  = [metrics[k].get('mean_usage_entropy')      for k in k_values]
    distances  = [metrics[k].get('mean_pairwise_distance')  for k in k_values]
    sparsities = [metrics[k].get('usage_sparsity_gini')     for k in k_values]
    rec_k      = recommendation.get('recommended_k')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('K-Value Selection Metrics', fontsize=14, weight='bold')

    def _plot_metric(ax, k_vals, y_vals, color, ylabel, title):
        pairs = [(k, y) for k, y in zip(k_vals, y_vals) if y is not None]
        if not pairs:
            ax.set_visible(False)
            return
        ks, ys = zip(*pairs)
        ax.plot(ks, ys, 'o-', linewidth=2, markersize=8, color=color)
        if rec_k in ks:
            ri = list(ks).index(rec_k)
            ax.axvline(rec_k, color='red', linestyle='--', linewidth=2, alpha=0.7)
            ax.plot(rec_k, ys[ri], 'r*', markersize=18, label=f'Recommended K={rec_k}')
            ax.legend(fontsize=9)
        ax.set_xlabel('K', fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, weight='bold')
        ax.grid(alpha=0.3)

    _plot_metric(axes[0, 0], k_values, entropies,  'steelblue',
                 'Mean Usage Entropy', 'Usage Diversity')
    _plot_metric(axes[0, 1], k_values, distances,  'coral',
                 'Mean Pairwise Distance', 'GEP Separation')
    _plot_metric(axes[1, 0], k_values, sparsities, 'seagreen',
                 'Usage Sparsity (Gini)', 'Usage Specificity')

    # Text panel
    ax = axes[1, 1]
    ax.axis('off')
    lines = ['K-Selection Recommendation', '=' * 34, '']
    if rec_k:
        lines += [
            f'Recommended K : {rec_k}',
            f'Confidence    : {recommendation.get("confidence", "unknown").upper()}',
            '',
            'Reasoning:',
        ]
        for r in recommendation.get('reasoning', []):
            lines.append(f'  - {r}')
        alt = recommendation.get('alternative_k', [])
        if alt:
            lines.append(f'\nAlternative K: {alt}')
    else:
        lines.append('No recommendation available')

    ax.text(0.05, 0.95, '\n'.join(lines), transform=ax.transAxes,
            fontsize=9, va='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.35))

    plt.tight_layout()
    out = output_dir / f'k_stability_summary.{cfg["figure_format"]}'
    plt.savefig(out, dpi=cfg['dpi'], bbox_inches='tight')
    plt.close()
    logger.info(f"  [OK] Saved: {out.name}")


def plot_local_density_histogram(
    cnmf_output_dir : Path,
    name            : str,
    k               : int,
    viz_dir         : Path,
    viz_config      : Optional[Dict] = None,
) -> bool:
    """
    Plot histogram of mean k-nearest-neighbor distances across NMF spectra.
    Used to verify density threshold for consensus clustering.
    """
    cfg = {**DEFAULT_VIZ_CONFIG, **(viz_config or {})}
    logger.info(f"  Local density histogram (K={k})...")

    try:
        dist_mat = get_gep_distance_matrix(cnmf_output_dir, name, k)
        if dist_mat is None:
            logger.warning(f"  [WARN] No distance matrix for K={k}")
            return False

        n     = dist_mat.shape[0]
        nbrs  = min(cfg['local_density_k_nbrs'], n - 1)
        if nbrs < 1:
            return False

        mean_dists = []
        for i in range(n):
            d = np.asarray(dist_mat[i, :], dtype=float)
            d = d[np.isfinite(d)]
            if d.size == 0:
                continue
            knn = np.sort(d)[1: nbrs + 1]
            if knn.size:
                mean_dists.append(float(np.mean(knn)))

        if not mean_dists:
            return False

        arr    = np.array(mean_dists)
        med    = np.median(arr)
        mad    = np.median(np.abs(arr - med))
        thresh = med + 2 * mad
        n_above = int(np.sum(arr > thresh))

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(arr, bins=cfg['local_density_bins'],
                color='steelblue', edgecolor='black', alpha=0.7)
        ax.axvline(thresh, color='red', linestyle='--', linewidth=2,
                   label=f'Threshold = {thresh:.3f}')
        ax.set_xlabel('Mean distance to k nearest neighbors', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(f'Local Density Diagnostic (K={k})', fontsize=12, weight='bold')
        ax.text(0.98, 0.98,
                f'{n_above}/{len(arr)} ({100*n_above/len(arr):.0f}%) above threshold',
                transform=ax.transAxes, fontsize=9,
                va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.legend(fontsize=9)
        plt.tight_layout()

        out = viz_dir / f'local_density_k{k}.{cfg["figure_format"]}'
        plt.savefig(out, dpi=cfg['dpi'], bbox_inches='tight')
        plt.close()
        logger.info(f"    [OK] {out.name}")
        return True

    except Exception as e:
        logger.warning(f"  [WARN] local_density failed (K={k}): {e}")
        plt.close('all')
        return False


def plot_clustergram(
    cnmf_output_dir : Path,
    name            : str,
    k               : int,
    viz_dir         : Path,
    viz_config      : Optional[Dict] = None,
) -> bool:
    """
    Plot hierarchical clustergram of GEP pairwise distances.
    Reveals redundant or poorly separated GEPs.
    """
    cfg = {**DEFAULT_VIZ_CONFIG, **(viz_config or {})}
    logger.info(f"  Clustergram (K={k})...")

    try:
        dist_mat = get_gep_distance_matrix(cnmf_output_dir, name, k)
        if dist_mat is None:
            return False

        # Use squareform to convert the full distance matrix to condensed form,
        # then pass directly to linkage as precomputed distances.
        # Do NOT call pdist(dist_mat) here — that would compute distances between
        # distance-profile vectors (double metric), not between original GEPs.
        dv       = squareform(dist_mat, checks=False)
        link_mat = linkage(dv, method=cfg['cluster_method'])

        fig = plt.figure(figsize=(12, 10))

        ax_d = fig.add_axes([0.15, 0.75, 0.7, 0.20])
        dend = dendrogram(link_mat, ax=ax_d, color_threshold=0,
                          above_threshold_color='steelblue')
        ax_d.set_xticks([])
        ax_d.set_yticks([])
        for sp in ax_d.spines.values():
            sp.set_visible(False)

        order = dend['leaves']
        dm_ord = dist_mat[order, :][:, order]

        ax_h = fig.add_axes([0.15, 0.15, 0.70, 0.60])
        im   = ax_h.imshow(dm_ord, cmap='viridis', aspect='auto', interpolation='nearest')

        labels = [f'GEP_{i+1}' for i in order]
        ax_h.set_xticks(range(len(labels)))
        ax_h.set_yticks(range(len(labels)))
        ax_h.set_xticklabels(labels, rotation=90, fontsize=8)
        ax_h.set_yticklabels(labels, fontsize=8)

        cb_ax = fig.add_axes([0.87, 0.15, 0.02, 0.60])
        cb    = plt.colorbar(im, cax=cb_ax)
        cb.set_label('Distance', rotation=270, labelpad=20, fontsize=10)

        fig.suptitle(f'GEP Similarity Clustergram (K={k})',
                     fontsize=14, weight='bold', y=0.98)

        out = viz_dir / f'clustergram_k{k}.{cfg["figure_format"]}'
        plt.savefig(out, dpi=cfg['dpi'], bbox_inches='tight')
        plt.close()
        logger.info(f"    [OK] {out.name}")
        return True

    except Exception as e:
        logger.warning(f"  [WARN] clustergram failed (K={k}): {e}")
        plt.close('all')
        return False


def plot_gep_usage_heatmap(
    cnmf_output_dir : Path,
    name            : str,
    k               : int,
    adata           : sc.AnnData,
    celltype_col    : str,
    viz_dir         : Path,
    viz_config      : Optional[Dict] = None,
) -> bool:
    """
    Heatmap of mean GEP usage per cell type (rows = cell types, cols = GEPs).
    Both axes are hierarchically clustered.
    """
    cfg = {**DEFAULT_VIZ_CONFIG, **(viz_config or {})}
    logger.info(f"  GEP usage heatmap (K={k})...")

    if not celltype_col:
        logger.warning("  [WARN] celltype_col is None, skipping usage heatmap")
        return False

    try:
        cnmf_run_dir = cnmf_output_dir / name
        usage_df = load_usage_matrix(cnmf_run_dir, name, k)
        if usage_df is None:
            return False

        common = adata.obs_names.intersection(usage_df.index)
        if len(common) == 0:
            logger.warning("  [WARN] No overlapping barcodes")
            return False

        usage_aligned = usage_df.loc[common].copy()
        ct_labels     = adata.obs.loc[common, celltype_col].astype(str)
        usage_aligned['_celltype'] = ct_labels.values

        mean_mat = usage_aligned.groupby('_celltype').mean()
        mean_mat = mean_mat.fillna(0.0)

        n_rows, n_cols = mean_mat.shape
        if n_rows == 0 or n_cols == 0:
            return False

        # Hierarchical ordering
        def _hclust_order(mat: np.ndarray) -> List[int]:
            if mat.shape[0] < 2:
                return list(range(mat.shape[0]))
            try:
                return list(dendrogram(linkage(mat, method='average'), no_plot=True)['leaves'])
            except Exception:
                return list(range(mat.shape[0]))

        row_ord = _hclust_order(mean_mat.values)
        col_ord = _hclust_order(mean_mat.values.T)
        ordered = mean_mat.iloc[row_ord, col_ord]

        csv_out = viz_dir / f'gep_usage_mean_by_celltype_k{k}.csv'
        ordered.to_csv(csv_out)

        fig, ax = plt.subplots(figsize=(max(12, n_cols * 0.6), max(4, n_rows * 0.4)))
        im = ax.imshow(ordered.values, cmap='RdYlBu_r', aspect='auto', interpolation='nearest')

        ax.set_xticks(range(n_cols))
        ax.set_yticks(range(n_rows))
        ax.set_xticklabels(ordered.columns, rotation=90, fontsize=9)
        ax.set_yticklabels(ordered.index, fontsize=9)
        ax.set_xlabel('GEP', fontsize=11)
        ax.set_ylabel('Cell Type', fontsize=11)
        ax.set_title(f'Mean GEP Usage by Cell Type (K={k})', fontsize=12, weight='bold')

        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label('Mean Usage', rotation=270, labelpad=20, fontsize=10)

        plt.tight_layout()
        out = viz_dir / f'gep_usage_heatmap_k{k}.{cfg["figure_format"]}'
        plt.savefig(out, dpi=cfg['dpi'], bbox_inches='tight')
        plt.close()
        logger.info(f"    [OK] {out.name}; {csv_out.name}")
        return True

    except Exception as e:
        logger.warning(f"  [WARN] usage heatmap failed (K={k}): {e}")
        plt.close('all')
        return False


def plot_umap_gep_usage(
    adata           : sc.AnnData,
    k               : int,
    viz_dir         : Path,
    obsm_key        : Optional[str]  = None,
    n_geps_to_show  : int            = 8,
    umap_key        : str            = 'X_umap',
    viz_config      : Optional[Dict] = None,
) -> bool:
    """
    Plot per-GEP usage on UMAP (grid of sub-panels).

    Parameters
    ----------
    obsm_key       : key in adata.obsm containing usage matrix;
                     defaults to f'X_cnmf_k{k}'
    n_geps_to_show : how many GEPs to visualize (top columns by total usage)
    umap_key       : adata.obsm key with 2D coordinates
    """
    cfg = {**DEFAULT_VIZ_CONFIG, **(viz_config or {})}

    if obsm_key is None:
        obsm_key = f'X_cnmf_k{k}'

    if obsm_key not in adata.obsm:
        logger.warning(f"  [WARN] {obsm_key} not in adata.obsm. Run add_usage_to_adata() first.")
        return False

    if umap_key not in adata.obsm:
        logger.warning(f"  [WARN] {umap_key} not in adata.obsm.")
        return False

    usage_mat = adata.obsm[obsm_key]  # (n_cells, k)
    coords    = adata.obsm[umap_key]  # (n_cells, 2)

    # Select top GEPs by total usage
    col_totals = np.nansum(usage_mat, axis=0)
    top_idx    = np.argsort(col_totals)[::-1][:n_geps_to_show]

    n_show = len(top_idx)
    ncols  = min(4, n_show)
    nrows  = int(np.ceil(n_show / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes_flat = np.array(axes).flatten() if nrows * ncols > 1 else [axes]

    for i, gep_i in enumerate(top_idx):
        ax     = axes_flat[i]
        scores = usage_mat[:, gep_i]
        finite = np.isfinite(scores)

        sc_plot = ax.scatter(
            coords[finite, 0], coords[finite, 1],
            c     = scores[finite],
            cmap  = 'viridis',
            s     = 1,
            alpha = 0.6,
            rasterized=False,
        )
        plt.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(f'GEP {gep_i + 1}', fontsize=10, weight='bold')
        ax.set_xlabel('UMAP 1', fontsize=8)
        ax.set_ylabel('UMAP 2', fontsize=8)
        ax.tick_params(labelsize=7)

    # Hide unused panels
    for j in range(n_show, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle(f'GEP Usage on UMAP (K={k}, top {n_show} GEPs)',
                 fontsize=13, weight='bold')
    plt.tight_layout()

    out = viz_dir / f'umap_gep_usage_k{k}.{cfg["figure_format"]}'
    plt.savefig(out, dpi=cfg['dpi'], bbox_inches='tight')
    plt.close()
    logger.info(f"  [OK] {out.name}")
    return True


def generate_all_visualizations(
    cnmf_output_dir : Path,
    name            : str,
    k_range         : List[int],
    adata           : sc.AnnData,
    celltype_col    : Optional[str],
    output_dir      : Path,
    viz_config      : Optional[Dict] = None,
) -> Dict[int, Dict[str, bool]]:
    """
    Generate all three standard visualizations for each K:
      - local_density histogram
      - GEP clustergram
      - GEP usage heatmap by cell type

    Returns nested dict: {k: {plot_name: True/False}}
    """
    logger.info("\n" + "=" * 70)
    logger.info("[INFO] Generating visualizations")
    logger.info("=" * 70)

    viz_dir = output_dir / 'visualizations'
    if not safe_mkdir(viz_dir, 'visualization directory'):
        return {}

    results: Dict[int, Dict[str, bool]] = {}

    for k in k_range:
        logger.info(f"\n--- K={k} ---")
        kr: Dict[str, bool] = {}

        kr['local_density'] = plot_local_density_histogram(
            cnmf_output_dir, name, k, viz_dir, viz_config)

        kr['clustergram']   = plot_clustergram(
            cnmf_output_dir, name, k, viz_dir, viz_config)

        kr['usage_heatmap'] = plot_gep_usage_heatmap(
            cnmf_output_dir, name, k, adata, celltype_col, viz_dir, viz_config)

        n_ok = sum(kr.values())
        logger.info(f"  {n_ok}/3 plots succeeded")

        results[k] = kr
        plt.close('all')
        gc.collect()

    total_ok = sum(sum(v.values()) for v in results.values())
    total    = len(k_range) * 3
    logger.info(f"\n[OK] Visualization summary: {total_ok}/{total} plots")

    return results


# ============================================================================
# SECTION 9 — FULL PIPELINE WRAPPER
# ============================================================================

def run_cnmf_full(
    adata          : sc.AnnData,
    output_dir     : Path,
    run_name       : str,
    k_range        : Optional[List[int]] = None,
    celltype_col   : Optional[str]       = None,
    batch_col      : Optional[str]       = None,
    use_batch_hvg  : bool                = True,
    cnmf_config    : Optional[Dict]      = None,
    viz_config     : Optional[Dict]      = None,
) -> Dict[str, Any]:
    """
    Complete end-to-end cNMF workflow:
      1. HVG selection (batch-aware → fallback)
      2. Input preparation (hvg_counts, tp10k, hvg_genes.txt)
      3. cNMF pipeline (prepare → factorize → combine → consensus)
      4. K-value stability metrics + recommendation
            5. Export per-GEP top-gene tables (wide + long/score)
            6. All standard visualizations

    Parameters
    ----------
    adata         : AnnData with layers['counts']
    output_dir    : output directory (created if absent)
    run_name      : cNMF run identifier (will be sanitized)
    k_range       : list of K values; auto-selected from cell count if None
    celltype_col  : obs column for usage heatmap; auto-inferred if None
    batch_col     : obs column for batch-aware HVG; auto-inferred if None
    use_batch_hvg : attempt batch-aware HVG selection
    cnmf_config   : override dict for DEFAULT_CNMF_CONFIG
    viz_config    : override dict for DEFAULT_VIZ_CONFIG

    Returns
    -------
    result dict with keys: success, k_range, hvg_method, metrics, recommendation,
                           visualizations, paths
    """
    cfg    = {**DEFAULT_CNMF_CONFIG, **(cnmf_config or {})}
    result: Dict[str, Any] = {
        'success'      : False,
        'run_name'     : run_name,
        'k_range'      : [],
        'hvg_method'   : None,
        'tech_genes'   : {},
        'gene_exclusion': {},
        'metrics'      : {},
        'recommendation': {},
        'visualizations': {},
        'paths'        : {},
    }

    if not CNMF_AVAILABLE:
        result['error'] = 'cNMF not installed'
        return result

    if not safe_mkdir(output_dir, f'run {run_name}'):
        result['error'] = 'Cannot create output_dir'
        return result

    # Auto-infer columns
    if celltype_col is None:
        celltype_col = infer_column(adata, CELLTYPE_COL_CANDIDATES)
        if celltype_col:
            logger.info(f"[INFO] Auto-detected celltype_col: {celltype_col}")
        else:
            logger.warning("[WARN] No celltype column found; usage heatmap will be skipped")

    if batch_col is None:
        batch_col = infer_column(adata, BATCH_COL_CANDIDATES)
        if batch_col:
            logger.info(f"[INFO] Auto-detected batch_col: {batch_col}")

    # Gene exclusion (shared contract with R-side helpers)
    try:
        filtered_adata, exclusion_packet, exclusion_manifest = apply_gene_exclusion_to_adata(
            adata,
            output_dir=output_dir,
            gene_exclusion_config=cfg.get('gene_exclusion_config') or None,
            prefix='gene_exclusion',
            celltype_col=celltype_col,
        )
        adata = filtered_adata
        result['gene_exclusion'] = exclusion_manifest
        result['paths']['gene_exclusion_summary_csv'] = exclusion_manifest.get('summary_csv')
        result['paths']['gene_exclusion_audit_csv'] = exclusion_manifest.get('audit_csv')
        result['paths']['gene_exclusion_manifest_json'] = exclusion_manifest.get('manifest_json')
        logger.info(
            "[INFO] Gene exclusion kept %d/%d features (removed %d)",
            exclusion_manifest.get('n_kept_features', adata.n_vars),
            exclusion_manifest.get('n_total_features', adata.n_vars),
            exclusion_manifest.get('n_excluded_features', 0),
        )
    except Exception as e:
        result['error'] = f'Gene exclusion failed: {e}'
        logger.error(f"[WARN] {result['error']}")
        return result

    # K range
    if k_range is None:
        k_range = determine_k_range(adata.n_obs)
        logger.info(f"[INFO] Auto K range: {k_range}")
    result['k_range'] = k_range

    name_safe = safe_name(run_name)

    # 1. HVG selection
    try:
        hvg_genes, method, tech = select_hvg_robust(
            adata,
            num_hvg      = cfg['num_hvg'],
            batch_key    = batch_col if use_batch_hvg else None,
            hvg_flavor   = cfg['hvg_flavor'],
            exclude_tech = cfg['exclude_technical'],
        )
        result['hvg_method'] = method
        result['tech_genes'] = {k: len(v) for k, v in tech.items()}
    except Exception as e:
        result['error'] = f'HVG selection failed: {e}'
        logger.error(f"[WARN] {result['error']}")
        return result

    # 2. Prepare inputs
    hvg_counts, tp10k, hvg_txt = prepare_cnmf_inputs(adata, hvg_genes, output_dir)
    if hvg_counts is None:
        result['error'] = 'Input preparation failed'
        return result

    result['paths']['hvg_counts'] = str(hvg_counts)
    result['paths']['tp10k']      = str(tp10k)
    result['paths']['hvg_txt']    = str(hvg_txt)

    # 3. cNMF pipeline
    ok = run_cnmf_pipeline(
        hvg_counts_path = hvg_counts,
        tp10k_path      = tp10k,
        hvg_txt_path    = hvg_txt,
        k_range         = k_range,
        output_dir      = output_dir,
        name            = name_safe,
        cnmf_config     = cfg,
    )

    if not ok:
        result['error'] = 'cNMF pipeline failed'
        return result

    cnmf_output_dir = output_dir / 'cnmf_output'
    result['paths']['cnmf_output'] = str(cnmf_output_dir)

    # 4. K metrics
    metrics = calculate_k_stability_metrics(
        cnmf_output_dir = cnmf_output_dir,
        name            = name_safe,
        k_range         = k_range,
        output_dir      = output_dir,
        viz_config      = viz_config,
    )
    result['metrics']        = metrics
    result['recommendation'] = json.loads(
        (output_dir / 'k_selection_recommendation.json').read_text()
    )

    # 5. Export per-GEP gene tables
    gep_gene_tables = export_gep_gene_tables(
        cnmf_output_dir = cnmf_output_dir,
        name            = name_safe,
        k_range         = k_range,
        output_dir      = output_dir,
        n_top           = 50,
    )
    if gep_gene_tables:
        result['paths']['gep_gene_tables'] = gep_gene_tables

    # 6. Visualizations
    viz_results = generate_all_visualizations(
        cnmf_output_dir = cnmf_output_dir,
        name            = name_safe,
        k_range         = k_range,
        adata           = adata,
        celltype_col    = celltype_col,
        output_dir      = output_dir,
        viz_config      = viz_config,
    )
    result['visualizations'] = viz_results
    result['success']        = True

    logger.info("\n" + "=" * 70)
    logger.info(f"[OK] cNMF run complete: {run_name}")
    logger.info(f"  Recommended K: {result['recommendation'].get('recommended_k')}")
    logger.info("=" * 70)

    # Save summary
    summary_path = output_dir / 'run_summary.json'
    with open(summary_path, 'w') as fh:
        json.dump(result, fh, indent=2, default=str)

    return result


# ============================================================================
# SECTION 10 — CONVENIENCE / DOWNSTREAM UTILITIES
# ============================================================================

def extract_top_genes_per_gep(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
    n_top        : int = 30,
) -> Optional[pd.DataFrame]:
    """
    Return a DataFrame (n_top × k) where each column is a GEP and
    values are the top gene names ranked by spectra score.

    Useful for: LLM interpretation prompts, manual annotation.
    """
    gep_df = load_gep_score_dataframe(cnmf_run_dir, name, k)
    if gep_df is None:
        return None

    records = {}
    for gep_name, scores in gep_df.iterrows():
        top_series = scores.sort_values(ascending=False).head(n_top)
        records[gep_name] = top_series.index.tolist()

    return pd.DataFrame(records)


def save_gep_gene_table(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
    output_dir   : Path,
    n_top        : int = 50,
) -> Optional[Path]:
    """
    Save top genes per GEP to CSV: output_dir / gep_top_genes_k{k}.csv
    """
    df = extract_top_genes_per_gep(cnmf_run_dir, name, k, n_top=n_top)
    if df is None:
        logger.warning(f"  [WARN] Cannot extract top genes (K={k})")
        return None

    out = output_dir / f'gep_top_genes_k{k}.csv'
    df.to_csv(out, index=True)
    logger.info(f"  [OK] Saved top genes: {out.name}")
    return out


def extract_top_genes_with_scores_per_gep(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
    n_top        : int = 30,
) -> Optional[pd.DataFrame]:
    """
    Return a long-format DataFrame with ranked genes and spectra scores per GEP.

    Columns: k, gep, rank, gene, score
    """
    gep_df = load_gep_score_dataframe(cnmf_run_dir, name, k)
    if gep_df is None:
        return None

    rows: List[Dict[str, Any]] = []
    for gep_name, scores in gep_df.iterrows():
        top_series = scores.sort_values(ascending=False).head(n_top)
        for rank, (gene, score) in enumerate(top_series.items(), start=1):
            rows.append({
                'k'    : int(k),
                'gep'  : str(gep_name),
                'rank' : int(rank),
                'gene' : str(gene),
                'score': float(score),
            })

    return pd.DataFrame(rows, columns=['k', 'gep', 'rank', 'gene', 'score'])


def save_gep_gene_score_table(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
    output_dir   : Path,
    n_top        : int = 50,
) -> Optional[Path]:
    """
    Save long-format ranked genes with spectra scores to TSV:
    output_dir / gep_gene_scores_k{k}.tsv
    """
    df = extract_top_genes_with_scores_per_gep(cnmf_run_dir, name, k, n_top=n_top)
    if df is None:
        logger.warning(f"  [WARN] Cannot extract scored top genes (K={k})")
        return None

    out = output_dir / f'gep_gene_scores_k{k}.tsv'
    df.to_csv(out, sep='\t', index=False)
    logger.info(f"  [OK] Saved scored top genes: {out.name}")
    return out


def export_gep_gene_tables(
    cnmf_output_dir : Path,
    name            : str,
    k_range         : List[int],
    output_dir      : Path,
    n_top           : int = 50,
) -> Dict[str, Any]:
    """
    Export per-K GEP gene tables in both wide and long formats.

    Returns a JSON-friendly path manifest with keys:
      directory, wide_by_k, long_by_k
    """
    gep_tables_dir = output_dir / 'gep_gene_tables'
    if not safe_mkdir(gep_tables_dir, 'GEP gene table directory'):
        return {}

    cnmf_run_dir = cnmf_output_dir / name
    manifest: Dict[str, Any] = {
        'directory': str(gep_tables_dir),
        'wide_by_k': {},
        'long_by_k': {},
    }

    for k in k_range:
        wide_path = save_gep_gene_table(
            cnmf_run_dir=cnmf_run_dir,
            name=name,
            k=k,
            output_dir=gep_tables_dir,
            n_top=n_top,
        )
        if wide_path is not None:
            manifest['wide_by_k'][str(k)] = str(wide_path)

        long_path = save_gep_gene_score_table(
            cnmf_run_dir=cnmf_run_dir,
            name=name,
            k=k,
            output_dir=gep_tables_dir,
            n_top=n_top,
        )
        if long_path is not None:
            manifest['long_by_k'][str(k)] = str(long_path)

    return manifest


def build_llm_gep_prompt(
    cnmf_run_dir : Path,
    name         : str,
    k            : int,
    n_top        : int = 30,
    context      : str = '',
) -> str:
    """
    Build a structured LLM interpretation prompt for cNMF GEPs.

    Returns a text prompt string ready for submission to an LLM API.

    Parameters
    ----------
    context : biological context string (e.g. 'B cell scRNA-seq, CRSwNP')
    """
    top_df = extract_top_genes_per_gep(cnmf_run_dir, name, k, n_top=n_top)
    if top_df is None:
        return ''

    lines = [
        f'You are a computational biologist interpreting Gene Expression Programs (GEPs)',
        f'from consensus Non-negative Matrix Factorization (cNMF) of scRNA-seq data.',
        '',
    ]
    if context:
        lines += [f'Biological context: {context}', '']

    lines += [
        f'Below are the top {n_top} genes for each of the {k} GEPs.',
        'For each GEP, provide:',
        '  1. A concise biological name (< 6 words)',
        '  2. Key biological processes or cell states',
        '  3. Confidence (high / medium / low)',
        '',
        '--- GEP Gene Lists ---',
    ]

    for col in top_df.columns:
        genes = ', '.join(top_df[col].dropna().tolist())
        lines.append(f'{col}: {genes}')

    lines += [
        '',
        'Respond as JSON:',
        '{"GEP_1": {"name": "...", "processes": "...", "confidence": "..."}, ...}',
    ]

    return '\n'.join(lines)


def print_run_summary(result: Dict[str, Any]) -> None:
    """Print a human-readable summary of run_cnmf_full() output."""
    print("=" * 60)
    print(f"cNMF Run Summary: {result.get('run_name', 'unknown')}")
    print("=" * 60)
    print(f"  Status       : {'SUCCESS' if result.get('success') else 'FAILED'}")
    if not result.get('success'):
        print(f"  Error        : {result.get('error', 'unknown')}")
        return

    print(f"  K range      : {result.get('k_range')}")
    print(f"  HVG method   : {result.get('hvg_method')}")

    rec   = result.get('recommendation', {})
    rec_k = rec.get('recommended_k')
    conf  = rec.get('confidence', 'unknown')
    print(f"  Recommended K: {rec_k} (confidence: {conf})")

    if rec.get('alternative_k'):
        print(f"  Alternative K: {rec['alternative_k']}")

    print()
    print("  Reasoning:")
    for r in rec.get('reasoning', []):
        print(f"    - {r}")

    viz = result.get('visualizations', {})
    total_viz = sum(sum(v.values()) for v in viz.values())
    total_all = len(viz) * 3
    print(f"\n  Plots        : {total_viz}/{total_all} generated")
    print("=" * 60)


# ============================================================================
# MODULE SELF-CHECK
# ============================================================================

if __name__ == '__main__':
    print("[INFO] cnmf_helper_20260419_v1_1.py loaded successfully")
    print(f"  cNMF available : {CNMF_AVAILABLE}")
    print("  Exported functions:")
    exported = [
        'setup_logger', 'safe_mkdir', 'safe_name', 'first_glob_match',
        'determine_k_range', 'infer_column', 'detect_technical_genes',
        'select_hvg_robust',
        'prepare_cnmf_inputs',
        'run_factorize', 'run_cnmf_pipeline',
        'find_gep_spectra_file', 'load_gep_score_dataframe',
        'load_gep_matrix', 'load_usage_matrix', 'load_top_genes',
        'add_usage_to_adata',
        'compute_gep_distance_matrix', 'get_gep_distance_matrix',
        'compute_k_metrics', 'calculate_k_stability_metrics',
        'generate_k_recommendation',
        'plot_k_metrics_summary', 'plot_local_density_histogram',
        'plot_clustergram', 'plot_gep_usage_heatmap', 'plot_umap_gep_usage',
        'generate_all_visualizations',
        'run_cnmf_full',
        'extract_top_genes_per_gep', 'extract_top_genes_with_scores_per_gep',
        'save_gep_gene_table', 'save_gep_gene_score_table', 'export_gep_gene_tables',
        'build_llm_gep_prompt', 'print_run_summary',
    ]
    for fn in exported:
        print(f"    - {fn}()")
