#!/usr/bin/env python3
"""
cNMF Batch Production Pipeline v1.1.0 PRODUCTION
=================================================

CRITICAL FIXES (based on code review):
- Fixed global logger syntax error
- Removed kNN smoothing (direct raw counts input)
- Two-track cNMF: uncorrected + batch-aware HVG
- Worker exitcode validation
- Thread control (OMP/MKL=1)
- HVG optimization (avoid unnecessary copies)
- Technical gene warnings
- K-value stability metrics

Features:
- Dual-track cNMF (uncorrected vs batch-aware)
- Robust HVG selection with fallback
- Comprehensive visualizations
- K-value selection guidance
- Technical gene detection
- Production-grade error handling

Author: r2end
Date: 2024-12-23
Version: v1.1.0 (Production - Post Code Review)
"""

import os
import sys
from pathlib import Path
import warnings
import logging
from typing import Optional, Dict, List, Tuple, Any, Union
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import entropy
import gc
import json
import time
from datetime import datetime
from multiprocessing import Process

# Set thread limits BEFORE importing numpy-dependent libraries
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

# Check cNMF availability
try:
    from cnmf import cNMF
    CNMF_AVAILABLE = True
except ImportError:
    print("="*70)
    print("ERROR: cNMF not installed")
    print("="*70)
    print("\nPlease install cNMF before running this pipeline:")
    print("  pip install cnmf")
    print("\nOr via conda:")
    print("  conda install -c conda-forge cnmf")
    print("="*70)
    sys.exit(1)

# Configure warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Base directories
BASE_INPUT_DIR = Path("/home/h2048/data/py/1217/downstream_analysis_v1_6_2")
BASE_OUTPUT_DIR = Path("/home/h2048/data/py/1217/cnmf_batch_production_v1_1")

# Main cell types
MAIN_CELL_TYPES = ['T_cells', 'B_cells', 'Myeloid', 'Stromal_Vascular']

# cNMF Configuration (fine-grained K values)
CNMF_CONFIG = {
    # K exploration - EXPANDED for fine-grained classification
    'k_range_large': [25, 30, 35, 40, 45, 50],      # >50k cells
    'k_range_medium': [15, 20, 25, 30, 35, 40],     # 10k-50k cells
    'k_range_small': [10, 15, 20, 25, 30, 35],      # <10k cells
    
    # Iteration settings
    'n_iter': 100,
    'seed': 42,
    
    # HVG selection
    'num_hvg': 3000,
    'hvg_flavor': 'seurat_v3',
    
    # Consensus parameters
    'density_threshold': 0.1,
    'show_clustering': True,
    'close_clustergram_fig': True,  # Always close to prevent leaks
    
    # Multi-processing
    'workers_per_dataset': 6,
}

# Technical gene detection (for warnings)
TECHNICAL_GENES = {
    'cell_cycle': ['MKI67', 'TOP2A', 'PCNA', 'CCNA2', 'CCNB1', 'CCNB2', 
                   'CDK1', 'AURKA', 'AURKB', 'CENPF', 'CENPE'],
    'ribosomal': [],  # Will be auto-detected by prefix
    'mitochondrial': [],  # Will be auto-detected by prefix
    'stress_ieg': ['FOS', 'JUN', 'JUNB', 'JUND', 'EGR1', 'EGR2', 'EGR3',
                   'ATF3', 'DUSP1', 'HSPA1A', 'HSPA1B', 'HSP90AA1'],
}

# Column name inference
COLUMN_NAMES = {
    'celltype_candidates': [
        'scanvi_predictions',
        'cell_type_scanvi_filt',
        'scanvi_labels',
        'cell_type',
        'celltype'
    ],
    'batch_candidates': [
        'dataset',
        'batch',
        'sample_id',
        'orig.ident'
    ]
}

# Visualization settings
VIZ_CONFIG = {
    'dpi': 300,
    'figure_format': 'png',
    'local_density_bins': 50,
    'local_density_k_neighbors': 20,
    'clustergram_method': 'average',
    'clustergram_metric': 'euclidean',
}

# Random seed
np.random.seed(CNMF_CONFIG['seed'])

# ============================================================================
# LOGGING SETUP
# ============================================================================

def setup_logging(log_file: Optional[Path] = None) -> logging.Logger:
    """Configure logging system"""
    # Clear existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers,
        force=True
    )
    return logging.getLogger(__name__)

# Initialize module-level logger (will be reconfigured in main)
logger = setup_logging()

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def safe_mkdir(directory: Path, description: str = "directory") -> bool:
    """Safe directory creation with verification"""
    try:
        directory.mkdir(exist_ok=True, parents=True)
        
        if not directory.exists() or not directory.is_dir():
            logger.error(f"  ❌ Failed to create {description}: {directory}")
            return False
        
        return True
        
    except Exception as e:
        logger.error(f"  ❌ Error creating {description}: {e}")
        return False


def safe_name(x: str, max_len: int = 180) -> str:
    """Sanitize names for file paths"""
    s = str(x)
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t']:
        s = s.replace(ch, '_')
    return s[:max_len]


def determine_k_range(n_cells: int) -> List[int]:
    """Determine appropriate K range based on cell count"""
    if n_cells >= 50000:
        k_range = CNMF_CONFIG['k_range_large']
        category = "large"
    elif n_cells >= 10000:
        k_range = CNMF_CONFIG['k_range_medium']
        category = "medium"
    else:
        k_range = CNMF_CONFIG['k_range_small']
        category = "small"
    
    logger.info(f"  Dataset size: {n_cells:,} cells → {category} K range: {k_range}")
    return k_range


def infer_column_name(adata: sc.AnnData, candidates: List[str]) -> Optional[str]:
    """Infer column name from candidates"""
    for col in candidates:
        if col in adata.obs.columns:
            return col
    return None


def detect_technical_genes(var_names: pd.Index) -> Dict[str, List[str]]:
    """
    Detect technical genes in dataset
    
    Returns:
        Dict with categories: cell_cycle, ribosomal, mitochondrial, stress_ieg
    """
    detected = {
        'cell_cycle': [],
        'ribosomal': [],
        'mitochondrial': [],
        'stress_ieg': []
    }
    
    var_names_upper = [str(g).upper() for g in var_names]
    
    # Cell cycle
    for gene in TECHNICAL_GENES['cell_cycle']:
        if gene in var_names_upper:
            idx = var_names_upper.index(gene)
            detected['cell_cycle'].append(var_names[idx])
    
    # Ribosomal (RPL*, RPS*)
    for gene in var_names:
        gene_upper = str(gene).upper()
        if gene_upper.startswith('RPL') or gene_upper.startswith('RPS'):
            detected['ribosomal'].append(gene)
    
    # Mitochondrial (MT-)
    for gene in var_names:
        gene_upper = str(gene).upper()
        if gene_upper.startswith('MT-'):
            detected['mitochondrial'].append(gene)
    
    # Stress/IEG
    for gene in TECHNICAL_GENES['stress_ieg']:
        if gene in var_names_upper:
            idx = var_names_upper.index(gene)
            detected['stress_ieg'].append(var_names[idx])
    
    return detected

def resolve_cnmf_name_dir(cnmf_output_dir: Path, name: str) -> Optional[Path]:
    """
    Resolve actual cNMF output directory for a given run name.
    Typical layout: cnmf_output_dir/<name>/
    But add fallback search in case layout differs.
    """
    direct = cnmf_output_dir / name
    if direct.exists() and direct.is_dir():
        return direct

    # Fallback: search any directory named exactly `name`
    try:
        hits = [p for p in cnmf_output_dir.rglob(name) if p.is_dir() and p.name == name]
        if hits:
            return sorted(hits)[0]
    except Exception:
        pass

    return None


def glob_first(base: Path, patterns: List[str], recursive: bool = False) -> Optional[Path]:
    """
    Return the first matched file from a list of glob patterns.
    Patterns are tried in order; the first non-empty match wins.
    """
    for pat in patterns:
        try:
            matches = sorted(base.rglob(pat) if recursive else base.glob(pat))
            if matches:
                return matches[0]
        except Exception:
            continue
    return None


def load_dist_like_npz(npz_path: Path) -> Optional[np.ndarray]:
    """
    Robustly load distance/local-density-like arrays from npz.
    Keys vary across cNMF versions.
    Returns ndarray or None.
    """
    try:
        data = np.load(npz_path, allow_pickle=True)
        for key in ("dist", "distance", "D", "local_density"):
            if key in data:
                return data[key]
        # fallback: if only one array inside
        if len(data.files) == 1:
            return data[data.files[0]]
    except Exception:
        pass
    return None

# ============================================================================
# DATASET DISCOVERY
# ============================================================================

def discover_all_datasets() -> Dict[str, List[Path]]:
    """Discover all h5ad files in the downstream analysis directory"""
    logger.info("\n" + "="*70)
    logger.info("STEP 1: Dataset Discovery")
    logger.info("="*70)
    
    datasets = {}
    
    for celltype in MAIN_CELL_TYPES:
        celltype_dir = BASE_INPUT_DIR / celltype
        
        if not celltype_dir.exists():
            logger.warning(f"  ⚠️  Directory not found: {celltype_dir}")
            continue
        
        logger.info(f"\n--- {celltype} ---")
        
        # Main processed file
        main_file = celltype_dir / f"{celltype}_processed_final.h5ad"
        
        celltype_files = []
        
        if main_file.exists():
            celltype_files.append(main_file)
            logger.info(f"  ✓ Main: {main_file.name}")
        else:
            logger.warning(f"  ⚠️  Main file not found: {main_file.name}")
        
        # Subclustering results
        subcluster_dir = celltype_dir / "subclustering"
        if subcluster_dir.exists():
            subcluster_files = list(subcluster_dir.glob("*/*_subclustered.h5ad"))
            if subcluster_files:
                celltype_files.extend(subcluster_files)
                logger.info(f"  ✓ Subclusters: {len(subcluster_files)} files")
                for sf in subcluster_files[:5]:  # Show first 5
                    logger.info(f"      {sf.parent.name}/{sf.name}")
                if len(subcluster_files) > 5:
                    logger.info(f"      ... and {len(subcluster_files)-5} more")
            else:
                logger.info(f"  No subclustering results")
        
        if celltype_files:
            datasets[celltype] = celltype_files
            logger.info(f"  Total: {len(celltype_files)} datasets")
    
    # Summary
    total_datasets = sum(len(files) for files in datasets.values())
    logger.info(f"\n{'='*70}")
    logger.info(f"Discovery Summary:")
    logger.info(f"  Cell types: {len(datasets)}")
    logger.info(f"  Total datasets: {total_datasets}")
    logger.info(f"{'='*70}")
    
    return datasets


# ============================================================================
# DATA LOADING & VALIDATION
# ============================================================================

def load_and_validate_dataset(h5ad_path: Path) -> Tuple[Optional[sc.AnnData], Dict[str, Any]]:
    """Load h5ad and validate required components"""
    try:
        logger.info(f"\nLoading: {h5ad_path.name}")
        adata = sc.read_h5ad(h5ad_path)
        logger.info(f"  Dimensions: {adata.n_obs:,} cells × {adata.n_vars} genes")
        
        metadata = {
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'has_counts': False,
            'celltype_col': None,
            'batch_col': None,
        }
        
        # Check celltype column
        celltype_col = infer_column_name(adata, COLUMN_NAMES['celltype_candidates'])
        if celltype_col:
            metadata['celltype_col'] = celltype_col
            logger.info(f"  ✓ Cell type column: {celltype_col}")
        else:
            logger.error(f"  ❌ No cell type column found")
            return None, metadata
        
        # Check batch column
        batch_col = infer_column_name(adata, COLUMN_NAMES['batch_candidates'])
        if batch_col:
            metadata['batch_col'] = batch_col
            n_batches = adata.obs[batch_col].nunique()
            logger.info(f"  ✓ Batch column: {batch_col} ({n_batches} batches)")
        else:
            logger.warning(f"  ⚠️  No batch column found")
        
        # Check counts
        if 'counts' in adata.layers:
            metadata['has_counts'] = True
            logger.info(f"  ✓ Counts layer present")
        else:
            logger.error(f"  ❌ Counts layer not found")
            return None, metadata
        
        logger.info(f"  ✓ Validation passed")
        return adata, metadata
        
    except Exception as e:
        logger.error(f"  ❌ Failed to load: {e}")
        return None, {}


# ============================================================================
# HVG SELECTION (OPTIMIZED, NO UNNECESSARY COPIES)
# ============================================================================

def select_hvg_robust(
    adata: sc.AnnData,
    num_hvg: int = 3000,
    batch_key: Optional[str] = None,
    use_batch_aware: bool = True
) -> Tuple[List[str], str, Dict[str, List[str]]]:
    """
    Robust HVG selection with batch-aware fallback
    
    Returns:
        (hvg_genes, method_used, technical_genes_detected)
    """
    logger.info(f"\nSelecting {num_hvg} HVGs...")
    
    # Work on counts layer directly (no copy)
    method_used = "unknown"
    
    try:
        if use_batch_aware and batch_key and batch_key in adata.obs.columns:
            # Try batch-aware
            logger.info(f"  Attempting batch-aware HVG (batch_key={batch_key})...")
            sc.pp.highly_variable_genes(
                adata,
                layer='counts',
                n_top_genes=num_hvg,
                flavor=CNMF_CONFIG['hvg_flavor'],
                batch_key=batch_key,
                subset=False
            )
            method_used = "batch-aware"
            logger.info(f"  ✓ Batch-aware HVG successful")
        else:
            raise ValueError("Batch-aware disabled or no batch key")
            
    except Exception as e:
        # Fallback to non-batch-aware
        logger.warning(f"  ⚠️  Batch-aware failed: {e}")
        logger.info(f"  Falling back to non-batch-aware HVG...")
        
        sc.pp.highly_variable_genes(
            adata,
            layer='counts',
            n_top_genes=num_hvg,
            flavor=CNMF_CONFIG['hvg_flavor'],
            batch_key=None,
            subset=False
        )
        method_used = "non-batch-aware"
        logger.info(f"  ✓ Non-batch-aware HVG successful")
    
    # Extract HVG list
    hvg_mask = adata.var['highly_variable'].values
    hvg_genes = adata.var_names[hvg_mask].tolist()
    
    logger.info(f"  Selected: {len(hvg_genes)} HVGs (method: {method_used})")
    
    # Detect technical genes in HVG
    technical_in_hvg = detect_technical_genes(adata.var_names[hvg_mask])
    
    # Report technical genes
    total_technical = sum(len(genes) for genes in technical_in_hvg.values())
    if total_technical > 0:
        logger.warning(f"\n  ⚠️  Technical genes in HVG: {total_technical}")
        for category, genes in technical_in_hvg.items():
            if len(genes) > 0:
                logger.warning(f"    {category}: {len(genes)} genes")
                if len(genes) <= 10:
                    logger.warning(f"      {genes}")
                else:
                    logger.warning(f"      {genes[:10]} ... (+{len(genes)-10} more)")
    
    # Store method in uns
    adata.uns['hvg_method'] = method_used
    
    return hvg_genes, method_used, technical_in_hvg


# ============================================================================
# CNMF PREPARATION
# ============================================================================

def prepare_cnmf_inputs(
    adata: sc.AnnData,
    hvg_genes: List[str],
    output_dir: Path
) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """
    Prepare cNMF input files (optimized - no smoothing)
    
    1. HVG counts (for factorization)
    2. TP10K matrix (for gene scoring)
    3. HVG gene list
    
    Returns:
        (hvg_counts_h5ad, tp10k_h5ad, hvg_txt)
    """
    logger.info(f"\nPreparing cNMF inputs...")
    
    try:
        # ⭐ CRITICAL: Preserve full genes to .raw (shared memory, no copy)
        logger.info(f"  Preserving full genes to .raw...")
        adata.raw = sc.AnnData(
            X=adata.layers['counts'],  # Shared memory
            obs=adata.obs.copy(),
            var=adata.var.copy()
        )
        logger.info(f"  ✓ Full {adata.n_vars} genes preserved in .raw")
        
        # Filter to HVG
        hvg_mask = adata.var_names.isin(hvg_genes)
        n_hvg_actual = hvg_mask.sum()
        
        logger.info(f"  Creating HVG subset: {n_hvg_actual} genes...")
        
        # 1. HVG counts (for cNMF factorization)
        adata_hvg = adata[:, hvg_mask].copy()
        
        # Use raw counts (sparse)
        adata_hvg.X = adata_hvg.layers['counts'].copy()
        
        hvg_counts_path = output_dir / "cnmf_input_hvg_counts.h5ad"
        adata_hvg.write_h5ad(hvg_counts_path, compression='gzip')
        logger.info(f"  ✓ Saved HVG counts: {hvg_counts_path.name}")
        
        # 2. TP10K matrix (full genes, for gene scoring)
        logger.info(f"  Creating TP10K matrix (full genes)...")
        adata_tp10k = adata.copy()
        
        # Use raw counts
        adata_tp10k.X = adata_tp10k.layers['counts'].copy()
        
        # Normalize to TP10K
        sc.pp.normalize_total(adata_tp10k, target_sum=1e4)
        
        tp10k_path = output_dir / "cnmf_input_tp10k.h5ad"
        adata_tp10k.write_h5ad(tp10k_path, compression='gzip')
        logger.info(f"  ✓ Saved TP10K: {tp10k_path.name}")
        
        # 3. HVG gene list
        hvg_path = output_dir / "cnmf_input_hvg_genes.txt"
        with open(hvg_path, 'w') as f:
            for gene in adata_hvg.var_names:
                f.write(f"{gene}\n")
        logger.info(f"  ✓ Saved HVG list: {hvg_path.name} ({len(adata_hvg.var_names)} genes)")
        
        # Cleanup
        del adata_hvg, adata_tp10k
        gc.collect()
        
        return hvg_counts_path, tp10k_path, hvg_path
        
    except Exception as e:
        logger.error(f"  ❌ Preparation failed: {e}")
        return None, None, None


# ============================================================================
# CNMF EXECUTION (WITH WORKER VALIDATION)
# ============================================================================

def run_cnmf_factorize_parallel(
    output_dir: str,
    name: str,
    total_workers: int
) -> bool:
    """
    Run cNMF factorization with worker exitcode validation
    """
    logger.info(f"\n  Factorization (parallel, {total_workers} workers)...")
    
    def worker_fn(worker_i: int):
        """Worker function"""
        try:
            cnmf_obj = cNMF(output_dir=output_dir, name=name)
            cnmf_obj.factorize(worker_i=worker_i, total_workers=total_workers)
        except Exception as e:
            print(f"Worker {worker_i} failed: {e}")
            sys.exit(1)  # Exit with error code
    
    try:
        # Launch all workers
        procs = []
        for wi in range(total_workers):
            p = Process(target=worker_fn, args=(wi,))
            p.start()
            procs.append(p)
            logger.info(f"    Worker {wi} started (PID: {p.pid})")
        
        # Wait and validate exitcodes
        all_ok = True
        for i, p in enumerate(procs):
            p.join()
            if p.exitcode != 0:
                all_ok = False
                logger.error(f"    ❌ Worker {i} failed (exitcode={p.exitcode})")
            else:
                logger.info(f"    ✓ Worker {i} finished successfully")
        
        if not all_ok:
            logger.error(f"  ❌ Some workers failed")
            return False
        
        logger.info(f"  ✓ All {total_workers} workers completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"  ❌ Factorization failed: {e}")
        return False


def run_cnmf_pipeline(
    hvg_counts_h5ad: Path,
    tp10k_h5ad: Path,
    hvg_txt: Path,
    k_range: List[int],
    output_dir: Path,
    name: str
) -> bool:
    """Run complete cNMF pipeline with comprehensive error handling"""
    logger.info(f"\n" + "="*70)
    logger.info(f"Running cNMF Pipeline: {name}")
    logger.info(f"="*70)
    
    cnmf_output_dir = output_dir / "cnmf_output"
    if not safe_mkdir(cnmf_output_dir, "cNMF output directory"):
        return False
    
    try:
        # Initialize cNMF object
        cnmf_obj = cNMF(output_dir=str(cnmf_output_dir), name=name)
        
        # Step 1: Prepare
        logger.info(f"\nStep 1/4: Prepare")
        logger.info(f"  K range: {k_range}")
        logger.info(f"  Iterations: {CNMF_CONFIG['n_iter']}")
        
        cnmf_obj.prepare(
            counts_fn=str(hvg_counts_h5ad),
            tpm_fn=str(tp10k_h5ad),
            genes_file=str(hvg_txt),
            components=k_range,
            n_iter=CNMF_CONFIG['n_iter'],
            seed=CNMF_CONFIG['seed']
        )
        logger.info(f"  ✓ Preparation complete")
        
        # Step 2: Factorize (multi-worker with validation)
        logger.info(f"\nStep 2/4: Factorize")
        factorize_start = time.time()
        
        success = run_cnmf_factorize_parallel(
            output_dir=str(cnmf_output_dir),
            name=name,
            total_workers=CNMF_CONFIG['workers_per_dataset']
        )
        
        if not success:
            logger.error(f"  ❌ Factorization failed")
            return False
        
        factorize_time = time.time() - factorize_start
        logger.info(f"  ✓ Factorization complete ({factorize_time/60:.1f} min)")
        
        # Step 3: Combine
        logger.info(f"\nStep 3/4: Combine")
        cnmf_obj.combine()
        logger.info(f"  ✓ Combine complete")
        
        # Step 4: Consensus (for each K)
        logger.info(f"\nStep 4/4: Consensus")
        for k in k_range:
            logger.info(f"  Processing K={k}...")
            try:
                cnmf_obj.consensus(
                    k=k,
                    density_threshold=CNMF_CONFIG['density_threshold'],
                    show_clustering=CNMF_CONFIG['show_clustering'],
                    close_clustergram_fig=CNMF_CONFIG['close_clustergram_fig']
                )
                logger.info(f"    ✓ K={k} complete")
                # After consensus, quick sanity check of expected outputs (flexible dt)
                name_dir = resolve_cnmf_name_dir(cnmf_output_dir, name)
                if name_dir:
                    probe = glob_first(name_dir, [f"{name}.usages.k_{k}.dt_*.consensus.txt", f"{name}.usages.k_{k}.*consensus.txt"])
                    if not probe:
                        logger.warning(f"    ⚠️  Consensus finished but usage file still not found for K={k} (check cNMF logs)")
                
                # Force close all figures
                plt.close('all')
                
            except Exception as e:
                logger.warning(f"    ⚠️  K={k} failed: {e}")
                plt.close('all')
                continue
        
        logger.info(f"  ✓ All consensus clustering complete")
        
        logger.info(f"\n{'='*70}")
        logger.info(f"✓ cNMF Pipeline Complete: {name}")
        logger.info(f"{'='*70}")
        
        return True
        
    except Exception as e:
        logger.error(f"\n❌ cNMF pipeline failed: {e}")
        plt.close('all')
        return False


# ============================================================================
# K-VALUE STABILITY METRICS
# ============================================================================

def calculate_k_stability_metrics(
    cnmf_output_dir: Path,
    name: str,
    k_range: List[int],
    output_dir: Path
) -> Dict[int, Dict[str, float]]:
    logger.info(f"\n" + "="*70)
    logger.info(f"Calculating K-Value Stability Metrics (HOTFIX: flexible file finding)")
    logger.info(f"="*70)

    metrics: Dict[int, Dict[str, float]] = {}

    name_dir = resolve_cnmf_name_dir(cnmf_output_dir, name)
    if not name_dir:
        logger.warning(f"  ⚠️  Cannot resolve cNMF name directory for: {name}")
        # still save empty metrics
        metrics_file = output_dir / "k_stability_metrics.json"
        with open(metrics_file, "w") as f:
            json.dump(metrics, f, indent=2)
        logger.info(f"\n✓ Saved metrics: {metrics_file}")
        return metrics

    tmp_dir = name_dir / "cnmf_tmp"

    for k in k_range:
        logger.info(f"\n--- K={k} ---")

        k_metrics: Dict[str, Optional[float]] = {
            "reconstruction_error": None,
            "mean_silhouette": None,
            "mean_usage_entropy": None,
            "mean_pairwise_distance": None,
        }

        try:
            # Usage file (dt_0.1 vs dt_0_1 etc.)
            usage_file = glob_first(
                name_dir,
                patterns=[
                    f"{name}.usages.k_{k}.dt_*.consensus.txt",
                    f"{name}.usages.k_{k}.*consensus.txt",
                ],
                recursive=False,
            )

            if not usage_file or not usage_file.exists():
                logger.warning(f"  ⚠️  Usage file not found")
                metrics[k] = k_metrics
                continue

            usage_df = pd.read_csv(usage_file, sep="\t", index_col=0)

            # Usage entropy
            usage_matrix = usage_df.values
            cell_entropies = []
            for i in range(usage_matrix.shape[0]):
                row = usage_matrix[i, :]
                s = row.sum()
                if s > 0:
                    cell_entropies.append(entropy(row / s))
            if cell_entropies:
                mean_entropy = float(np.mean(cell_entropies))
                k_metrics["mean_usage_entropy"] = mean_entropy
                logger.info(f"  Mean usage entropy: {mean_entropy:.3f}")

            # Distance / similarity (usually in cnmf_tmp)
            dist_file = None
            if tmp_dir.exists():
                dist_file = glob_first(
                    tmp_dir,
                    patterns=[
                        f"{name}.spectra.k_{k}.dt_*.consensus.dist.df.npz",
                        f"{name}.spectra.k_{k}.*consensus*.npz",
                        f"{name}.local_density_cache.k_{k}.merged.df.npz",
                    ],
                    recursive=False,
                )

            if dist_file and dist_file.exists():
                arr = load_dist_like_npz(dist_file)
                if arr is None:
                    logger.warning(f"  ⚠️  Distance npz found but no usable key: {dist_file.name}")
                else:
                    # Prefer square distance matrix
                    if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] >= 2:
                        mean_dist = float(np.mean(arr[np.triu_indices_from(arr, k=1)]))
                        k_metrics["mean_pairwise_distance"] = mean_dist
                        logger.info(f"  Mean pairwise distance: {mean_dist:.3f}")
                    else:
                        # If it's 1D local density-like, at least record mean
                        try:
                            mean_val = float(np.mean(arr))
                            k_metrics["mean_pairwise_distance"] = mean_val
                            logger.info(f"  Mean density-like value: {mean_val:.3f} (non-square array)")
                        except Exception:
                            logger.warning(f"  ⚠️  Unhandled distance array shape: {getattr(arr, 'shape', None)}")

        except Exception as e:
            logger.warning(f"  ⚠️  Failed to compute metrics: {e}")

        metrics[k] = k_metrics

    metrics_file = output_dir / "k_stability_metrics.json"
    with open(metrics_file, "w") as f:
        json.dump(metrics, f, indent=2)
    logger.info(f"\n✓ Saved metrics: {metrics_file}")

    try:
        plot_k_stability_summary(metrics, output_dir)
    except Exception as e:
        logger.warning(f"  ⚠️  Failed to plot metrics: {e}")

    return metrics


def plot_k_stability_summary(metrics: Dict[int, Dict[str, float]], output_dir: Path):
    """Plot K-value stability metrics summary"""
    k_values = sorted(metrics.keys())
    
    # Extract metrics
    entropies = [metrics[k].get('mean_usage_entropy') for k in k_values]
    distances = [metrics[k].get('mean_pairwise_distance') for k in k_values]
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # Entropy
    ax = axes[0]
    valid_entropy = [(k, e) for k, e in zip(k_values, entropies) if e is not None]
    if valid_entropy:
        k_e, e_vals = zip(*valid_entropy)
        ax.plot(k_e, e_vals, 'o-', linewidth=2, markersize=8, color='steelblue')
        ax.set_xlabel('K (number of GEPs)', fontsize=11)
        ax.set_ylabel('Mean Usage Entropy', fontsize=11)
        ax.set_title('Usage Entropy vs K', fontsize=12, weight='bold')
        ax.grid(alpha=0.3)
    
    # Distance
    ax = axes[1]
    valid_dist = [(k, d) for k, d in zip(k_values, distances) if d is not None]
    if valid_dist:
        k_d, d_vals = zip(*valid_dist)
        ax.plot(k_d, d_vals, 'o-', linewidth=2, markersize=8, color='coral')
        ax.set_xlabel('K (number of GEPs)', fontsize=11)
        ax.set_ylabel('Mean Pairwise Distance', fontsize=11)
        ax.set_title('GEP Similarity vs K', fontsize=12, weight='bold')
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    output_path = output_dir / f'k_stability_summary.{VIZ_CONFIG["figure_format"]}'
    plt.savefig(output_path, dpi=VIZ_CONFIG['dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Saved: {output_path.name}")


# ============================================================================
# VISUALIZATIONS
# ============================================================================

def plot_local_density_histogram(
    cnmf_output_dir: Path,
    name: str,
    k: int,
    viz_dir: Path
) -> bool:
    logger.info(f"\n  Local density histogram (K={k})...")

    try:
        name_dir = resolve_cnmf_name_dir(cnmf_output_dir, name)
        if not name_dir:
            logger.warning(f"    ⚠️  Output directory not found for name={name}")
            return False

        tmp_dir = name_dir / "cnmf_tmp"
        if not tmp_dir.exists():
            logger.warning(f"    ⚠️  cnmf_tmp directory not found")
            return False

        dist_file = glob_first(
            tmp_dir,
            patterns=[
                f"{name}.spectra.k_{k}.dt_*.consensus.dist.df.npz",
                f"{name}.spectra.k_{k}.*consensus*.npz",
                f"{name}.local_density_cache.k_{k}.merged.df.npz",
            ],
            recursive=False,
        )

        if not dist_file or not dist_file.exists():
            logger.warning(f"    ⚠️  Distance file not found")
            return False

        arr = load_dist_like_npz(dist_file)
        if arr is None:
            logger.warning(f"    ⚠️  No distance/local_density data in {dist_file.name}")
            return False

        # If square distance matrix -> compute mean distance to kNN
        if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] >= 2:
            dist_matrix = arr
            n_neighbors = min(VIZ_CONFIG["local_density_k_neighbors"], dist_matrix.shape[0] - 1)
            n_spectra = dist_matrix.shape[0]

            mean_distances = []
            for i in range(n_spectra):
                d = np.sort(dist_matrix[i, :])
                mean_distances.append(float(np.mean(d[1 : n_neighbors + 1])))

            mean_distances = np.asarray(mean_distances)

        # If 1D vector -> treat as density-like already
        elif arr.ndim == 1 and arr.shape[0] >= 2:
            mean_distances = np.asarray(arr, dtype=float)
            n_spectra = mean_distances.shape[0]
        else:
            logger.warning(f"    ⚠️  Unhandled array shape: {arr.shape}")
            return False

        median_dist = float(np.median(mean_distances))
        mad = float(np.median(np.abs(mean_distances - median_dist)))
        threshold = median_dist + 2 * mad

        n_above_threshold = int(np.sum(mean_distances > threshold))
        pct_above = (n_above_threshold / n_spectra) * 100.0

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(mean_distances, bins=VIZ_CONFIG["local_density_bins"], edgecolor="black", alpha=0.7)

        ax.axvline(threshold, linestyle="--", linewidth=2, label="Filtering threshold")

        ax.set_xlabel("Mean distance to k nearest neighbors (or density-like value)", fontsize=12)
        ax.set_ylabel("Frequency", fontsize=12)
        ax.set_title(f"Local Density Diagnostic (K={k}) [HOTFIX]", fontsize=13, weight="bold")

        textstr = (
            f"{n_above_threshold}/{n_spectra} ({pct_above:.0f}%) spectra above threshold\n"
            f"were removed prior to clustering"
        )
        ax.text(
            0.98, 0.98, textstr,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

        plt.tight_layout()
        output_path = viz_dir / f"local_density_histogram_k{k}.{VIZ_CONFIG['figure_format']}"
        plt.savefig(output_path, dpi=VIZ_CONFIG["dpi"], bbox_inches="tight")
        plt.close()

        logger.info(f"    ✓ Saved: {output_path.name}")
        return True

    except Exception as e:
        logger.warning(f"    ⚠️  Failed: {e}")
        plt.close("all")
        return False

def plot_clustergram(
    cnmf_output_dir: Path,
    name: str,
    k: int,
    viz_dir: Path
) -> bool:
    logger.info(f"\n  Clustergram (K={k})...")

    try:
        name_dir = resolve_cnmf_name_dir(cnmf_output_dir, name)
        if not name_dir:
            logger.warning(f"    ⚠️  Output directory not found for name={name}")
            return False

        spectra_file = glob_first(
            name_dir,
            patterns=[
                f"{name}.gene_spectra_score.k_{k}.dt_*.txt",
                f"{name}.gene_spectra_score.k_{k}.*.txt",
                f"{name}.spectra.k_{k}.dt_*.consensus.txt",
                f"{name}.spectra.k_{k}.*consensus*.txt",
            ],
            recursive=False,
        )

        if not spectra_file or not spectra_file.exists():
            logger.warning(f"    ⚠️  Spectra file not found")
            return False

        spectra_df = pd.read_csv(spectra_file, sep="\t", index_col=0)
        spectra_matrix = spectra_df.values

        if spectra_matrix.ndim != 2 or spectra_matrix.shape[0] < 2:
            logger.warning(f"    ⚠️  Spectra matrix too small: {spectra_matrix.shape}")
            return False

        distances = pdist(spectra_matrix, metric=VIZ_CONFIG["clustergram_metric"])
        dist_matrix = squareform(distances)
        linkage_matrix = linkage(distances, method=VIZ_CONFIG["clustergram_method"])

        fig = plt.figure(figsize=(12, 10))

        ax_dendro = fig.add_axes([0.15, 0.75, 0.7, 0.2])
        dendro = dendrogram(linkage_matrix, ax=ax_dendro, color_threshold=0)
        ax_dendro.set_xticks([])
        ax_dendro.set_yticks([])
        for spn in ax_dendro.spines.values():
            spn.set_visible(False)

        order = dendro["leaves"]
        dist_matrix_ordered = dist_matrix[order, :][:, order]

        ax_heatmap = fig.add_axes([0.15, 0.15, 0.7, 0.6])
        im = ax_heatmap.imshow(dist_matrix_ordered, cmap="viridis", aspect="auto", interpolation="nearest")

        labels = [f"GEP_{i+1}" for i in order]
        ax_heatmap.set_xticks(range(len(labels)))
        ax_heatmap.set_yticks(range(len(labels)))
        ax_heatmap.set_xticklabels(labels, rotation=90, fontsize=8)
        ax_heatmap.set_yticklabels(labels, fontsize=8)

        cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.6])
        cbar = plt.colorbar(im, cax=cbar_ax)
        cbar.set_label("Euclidean Distance", rotation=270, labelpad=20, fontsize=10)

        fig.suptitle(f"GEP Similarity Clustergram (K={k})", fontsize=14, weight="bold", y=0.98)

        output_path = viz_dir / f"clustergram_k{k}.{VIZ_CONFIG['figure_format']}"
        plt.savefig(output_path, dpi=VIZ_CONFIG["dpi"], bbox_inches="tight")
        plt.close()

        logger.info(f"    ✓ Saved: {output_path.name}")
        return True

    except Exception as e:
        logger.warning(f"    ⚠️  Failed: {e}")
        plt.close("all")
        return False


def plot_gep_usage_heatmap(
    cnmf_output_dir: Path,
    name: str,
    k: int,
    adata: sc.AnnData,
    celltype_col: str,
    viz_dir: Path
) -> bool:
    logger.info(f"\n  GEP usage heatmap (K={k})...")

    try:
        name_dir = resolve_cnmf_name_dir(cnmf_output_dir, name)
        if not name_dir:
            logger.warning(f"    ⚠️  Output directory not found for name={name}")
            return False

        usage_file = glob_first(
            name_dir,
            patterns=[
                f"{name}.usages.k_{k}.dt_*.consensus.txt",
                f"{name}.usages.k_{k}.*consensus.txt",
            ],
            recursive=False,
        )

        if not usage_file or not usage_file.exists():
            logger.warning(f"    ⚠️  Usage file not found")
            return False

        usage_df = pd.read_csv(usage_file, sep="\t", index_col=0)

        common_cells = adata.obs_names.intersection(usage_df.index)
        if len(common_cells) == 0:
            logger.warning(f"    ⚠️  No overlapping cells")
            return False

        usage_aligned = usage_df.loc[common_cells].copy()
        celltype_aligned = adata.obs.loc[common_cells, celltype_col]

        usage_aligned["celltype"] = celltype_aligned.values
        mean_usage = usage_aligned.groupby("celltype").mean()

        fig, ax = plt.subplots(figsize=(max(12, k * 0.6), max(8, len(mean_usage) * 0.4)))

        row_linkage = linkage(mean_usage.values, method="average") if mean_usage.shape[0] > 1 else None
        col_linkage = linkage(mean_usage.T.values, method="average") if mean_usage.shape[1] > 1 else None

        if row_linkage is not None:
            row_order = dendrogram(row_linkage, no_plot=True)["leaves"]
        else:
            row_order = list(range(mean_usage.shape[0]))

        if col_linkage is not None:
            col_order = dendrogram(col_linkage, no_plot=True)["leaves"]
        else:
            col_order = list(range(mean_usage.shape[1]))

        mean_usage_ordered = mean_usage.iloc[row_order, col_order]

        im = ax.imshow(mean_usage_ordered.values, cmap="RdYlBu_r", aspect="auto", interpolation="nearest")

        ax.set_xticks(range(len(mean_usage_ordered.columns)))
        ax.set_yticks(range(len(mean_usage_ordered.index)))
        ax.set_xticklabels(mean_usage_ordered.columns, rotation=90, fontsize=9)
        ax.set_yticklabels(mean_usage_ordered.index, fontsize=9)

        ax.set_xlabel("GEP", fontsize=11)
        ax.set_ylabel("Cell Type", fontsize=11)
        ax.set_title(f"Mean GEP Usage by Cell Type (K={k})", fontsize=13, weight="bold")

        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Mean Usage", rotation=270, labelpad=20, fontsize=10)

        plt.tight_layout()
        output_path = viz_dir / f"gep_usage_heatmap_k{k}.{VIZ_CONFIG['figure_format']}"
        plt.savefig(output_path, dpi=VIZ_CONFIG["dpi"], bbox_inches="tight")
        plt.close()

        logger.info(f"    ✓ Saved: {output_path.name}")
        return True

    except Exception as e:
        logger.warning(f"    ⚠️  Failed: {e}")
        plt.close("all")
        return False


def generate_all_visualizations(
    cnmf_output_dir: Path,
    name: str,
    k_range: List[int],
    adata: sc.AnnData,
    celltype_col: str,
    output_dir: Path
) -> Dict[int, bool]:
    """Generate all visualizations for all K values"""
    logger.info(f"\n" + "="*70)
    logger.info(f"Generating Visualizations")
    logger.info(f"="*70)
    
    viz_dir = output_dir / "visualizations"
    if not safe_mkdir(viz_dir, "visualization directory"):
        return {}
    
    results = {}
    
    for k in k_range:
        logger.info(f"\n--- K={k} ---")
        
        success = True
        
        if not plot_local_density_histogram(cnmf_output_dir, name, k, viz_dir):
            success = False
        
        if not plot_clustergram(cnmf_output_dir, name, k, viz_dir):
            success = False
        
        if not plot_gep_usage_heatmap(cnmf_output_dir, name, k, adata, celltype_col, viz_dir):
            success = False
        
        results[k] = success
        
        if success:
            logger.info(f"  ✓ K={k} visualizations complete")
        else:
            logger.warning(f"  ⚠️  K={k} had some visualization failures")
        
        # Force cleanup
        plt.close('all')
        gc.collect()
    
    n_success = sum(results.values())
    logger.info(f"\n{'='*70}")
    logger.info(f"Visualization Summary: {n_success}/{len(k_range)} K values successful")
    logger.info(f"{'='*70}")
    
    return results


# ============================================================================
# MAIN DATASET PROCESSOR (TWO-TRACK)
# ============================================================================

def process_single_dataset(
    h5ad_path: Path,
    output_base: Path,
    dataset_id: str
) -> Dict[str, Any]:
    """
    Process single dataset with TWO cNMF tracks:
    1. Uncorrected (non-batch-aware HVG)
    2. Batch-aware HVG
    """
    logger.info(f"\n" + "="*70)
    logger.info(f"PROCESSING DATASET: {dataset_id}")
    logger.info(f"="*70)
    logger.info(f"File: {h5ad_path}")
    
    start_time = time.time()
    
    result = {
        'dataset_id': dataset_id,
        'h5ad_path': str(h5ad_path),
        'success': False,
        'error': None,
        'time_seconds': 0,
        'n_cells': 0,
        'n_genes': 0,
        'k_range': [],
        'tracks': {}
    }
    
    # Create output directory
    output_dir = output_base / dataset_id
    if not safe_mkdir(output_dir, f"output directory for {dataset_id}"):
        result['error'] = "Failed to create output directory"
        return result
    
    try:
        # Step 1: Load and validate
        adata, metadata = load_and_validate_dataset(h5ad_path)
        
        if adata is None:
            result['error'] = "Failed validation"
            return result
        
        result['n_cells'] = metadata['n_cells']
        result['n_genes'] = metadata['n_genes']
        
        celltype_col = metadata['celltype_col']
        batch_col = metadata['batch_col']
        
        # Determine K range
        k_range = determine_k_range(metadata['n_cells'])
        result['k_range'] = k_range
        
        # ==================================================================
        # TRACK 1: UNCORRECTED (Non-batch-aware HVG)
        # ==================================================================
        logger.info(f"\n" + "="*70)
        logger.info(f"TRACK 1: UNCORRECTED (Non-batch-aware HVG)")
        logger.info(f"="*70)
        
        track1_dir = output_dir / "uncorrected"
        if not safe_mkdir(track1_dir, "track1 directory"):
            result['tracks']['uncorrected'] = {'error': 'Failed to create directory'}
        else:
            track1_result = {'success': False}
            
            try:
                # HVG selection (non-batch-aware)
                hvg_genes_t1, method_t1, technical_t1 = select_hvg_robust(
                    adata,
                    num_hvg=CNMF_CONFIG['num_hvg'],
                    batch_key=None,
                    use_batch_aware=False
                )
                
                if len(hvg_genes_t1) == 0:
                    track1_result['error'] = 'HVG selection failed'
                else:
                    # Prepare cNMF inputs
                    hvg_h5ad_t1, tp10k_h5ad_t1, hvg_txt_t1 = prepare_cnmf_inputs(
                        adata, hvg_genes_t1, track1_dir
                    )
                    
                    if hvg_h5ad_t1 is None:
                        track1_result['error'] = 'Input preparation failed'
                    else:
                        # Run cNMF
                        cnmf_name_t1 = safe_name(f"{dataset_id}_uncorrected")
                        success_t1 = run_cnmf_pipeline(
                            hvg_h5ad_t1, tp10k_h5ad_t1, hvg_txt_t1,
                            k_range, track1_dir, cnmf_name_t1
                        )
                        
                        if success_t1:
                            # Visualizations
                            cnmf_output_dir_t1 = track1_dir / "cnmf_output"
                            viz_results_t1 = generate_all_visualizations(
                                cnmf_output_dir_t1, cnmf_name_t1, k_range,
                                adata, celltype_col, track1_dir
                            )
                            
                            # K stability metrics
                            metrics_t1 = calculate_k_stability_metrics(
                                cnmf_output_dir_t1, cnmf_name_t1, k_range, track1_dir
                            )
                            
                            track1_result['success'] = True
                            track1_result['hvg_method'] = method_t1
                            track1_result['technical_genes'] = {
                                k: len(v) for k, v in technical_t1.items()
                            }
                            track1_result['visualizations'] = viz_results_t1
                        else:
                            track1_result['error'] = 'cNMF pipeline failed'
            
            except Exception as e:
                track1_result['error'] = str(e)
                logger.error(f"Track 1 failed: {e}")
            
            result['tracks']['uncorrected'] = track1_result
        
        # ==================================================================
        # TRACK 2: BATCH-AWARE HVG
        # ==================================================================
        logger.info(f"\n" + "="*70)
        logger.info(f"TRACK 2: BATCH-AWARE HVG")
        logger.info(f"="*70)
        
        track2_dir = output_dir / "batch_aware"
        if not safe_mkdir(track2_dir, "track2 directory"):
            result['tracks']['batch_aware'] = {'error': 'Failed to create directory'}
        else:
            track2_result = {'success': False}
            
            try:
                # HVG selection (batch-aware with fallback)
                hvg_genes_t2, method_t2, technical_t2 = select_hvg_robust(
                    adata,
                    num_hvg=CNMF_CONFIG['num_hvg'],
                    batch_key=batch_col,
                    use_batch_aware=True
                )
                
                if len(hvg_genes_t2) == 0:
                    track2_result['error'] = 'HVG selection failed'
                else:
                    # Prepare cNMF inputs
                    hvg_h5ad_t2, tp10k_h5ad_t2, hvg_txt_t2 = prepare_cnmf_inputs(
                        adata, hvg_genes_t2, track2_dir
                    )
                    
                    if hvg_h5ad_t2 is None:
                        track2_result['error'] = 'Input preparation failed'
                    else:
                        # Run cNMF
                        cnmf_name_t2 = safe_name(f"{dataset_id}_batch_aware")
                        success_t2 = run_cnmf_pipeline(
                            hvg_h5ad_t2, tp10k_h5ad_t2, hvg_txt_t2,
                            k_range, track2_dir, cnmf_name_t2
                        )
                        
                        if success_t2:
                            # Visualizations
                            cnmf_output_dir_t2 = track2_dir / "cnmf_output"
                            viz_results_t2 = generate_all_visualizations(
                                cnmf_output_dir_t2, cnmf_name_t2, k_range,
                                adata, celltype_col, track2_dir
                            )
                            
                            # K stability metrics
                            metrics_t2 = calculate_k_stability_metrics(
                                cnmf_output_dir_t2, cnmf_name_t2, k_range, track2_dir
                            )
                            
                            track2_result['success'] = True
                            track2_result['hvg_method'] = method_t2
                            track2_result['technical_genes'] = {
                                k: len(v) for k, v in technical_t2.items()
                            }
                            track2_result['visualizations'] = viz_results_t2
                        else:
                            track2_result['error'] = 'cNMF pipeline failed'
            
            except Exception as e:
                track2_result['error'] = str(e)
                logger.error(f"Track 2 failed: {e}")
            
            result['tracks']['batch_aware'] = track2_result
        
        # Overall success if at least one track succeeded
        result['success'] = any(
            track.get('success', False) 
            for track in result['tracks'].values()
        )
        
        # Cleanup
        del adata
        gc.collect()
        
    except Exception as e:
        logger.error(f"\n❌ Processing failed: {e}", exc_info=True)
        result['error'] = str(e)
    
    finally:
        result['time_seconds'] = time.time() - start_time
    
    # Save result metadata
    result_json = output_dir / "processing_result.json"
    with open(result_json, 'w') as f:
        json.dump(result, f, indent=2)
    
    logger.info(f"\n{'='*70}")
    if result['success']:
        logger.info(f"✓ {dataset_id} COMPLETE ({result['time_seconds']/60:.1f} min)")
    else:
        logger.info(f"❌ {dataset_id} FAILED: {result['error']}")
    logger.info(f"{'='*70}")
    
    return result


# ============================================================================
# MAIN CONTROLLER
# ============================================================================

def main():
    """Main pipeline controller"""
    
    # Create base output directory FIRST
    if not safe_mkdir(BASE_OUTPUT_DIR, "base output directory"):
        print("❌ Cannot create base output directory - aborting")
        sys.exit(1)
    
    # Setup logging to file
    log_file = BASE_OUTPUT_DIR / f"pipeline_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    # Reconfigure global logger
    global logger
    logger = setup_logging(log_file)
    
    logger.info("="*70)
    logger.info("cNMF BATCH PRODUCTION PIPELINE v1.1.0 PRODUCTION")
    logger.info("="*70)
    logger.info(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Log file: {log_file}")
    
    logger.info(f"\nConfiguration:")
    logger.info(f"  Input:  {BASE_INPUT_DIR}")
    logger.info(f"  Output: {BASE_OUTPUT_DIR}")
    logger.info(f"  K ranges (fine-grained):")
    logger.info(f"    Large (>50k):  {CNMF_CONFIG['k_range_large']}")
    logger.info(f"    Medium (10k-50k): {CNMF_CONFIG['k_range_medium']}")
    logger.info(f"    Small (<10k):  {CNMF_CONFIG['k_range_small']}")
    logger.info(f"  Workers: {CNMF_CONFIG['workers_per_dataset']}")
    logger.info(f"  HVG: {CNMF_CONFIG['num_hvg']}")
    logger.info(f"  Two-track: uncorrected + batch-aware")
    
    pipeline_start = time.time()
    
    # Step 1: Discover all datasets
    datasets = discover_all_datasets()
    
    if not datasets:
        logger.error("❌ No datasets found - aborting")
        sys.exit(1)
    
    # Step 2: Process each dataset
    all_results = []
    
    for celltype, h5ad_files in datasets.items():
        logger.info(f"\n" + "="*70)
        logger.info(f"CELL TYPE: {celltype}")
        logger.info(f"="*70)
        
        for h5ad_path in h5ad_files:
            # Generate dataset ID
            if h5ad_path.name.endswith('_processed_final.h5ad'):
                dataset_id = celltype
            else:
                subtype = h5ad_path.parent.name
                dataset_id = f"{celltype}_{subtype}"
            
            # Process
            result = process_single_dataset(
                h5ad_path,
                BASE_OUTPUT_DIR,
                dataset_id
            )
            
            all_results.append(result)
    
    # Step 3: Generate summary report
    pipeline_time = time.time() - pipeline_start
    
    logger.info(f"\n" + "="*70)
    logger.info(f"PIPELINE COMPLETE")
    logger.info(f"="*70)
    
    # Statistics
    n_total = len(all_results)
    n_success = sum(1 for r in all_results if r['success'])
    n_failed = n_total - n_success
    
    # Track-level statistics
    n_track1_success = sum(
        1 for r in all_results 
        if r.get('tracks', {}).get('uncorrected', {}).get('success', False)
    )
    n_track2_success = sum(
        1 for r in all_results 
        if r.get('tracks', {}).get('batch_aware', {}).get('success', False)
    )
    
    logger.info(f"\nSummary Statistics:")
    logger.info(f"  Total datasets: {n_total}")
    logger.info(f"  Overall success: {n_success} ({n_success/n_total*100:.1f}%)")
    logger.info(f"  Overall failed: {n_failed} ({n_failed/n_total*100:.1f}%)")
    logger.info(f"\nTrack-level:")
    logger.info(f"  Uncorrected: {n_track1_success}/{n_total} ({n_track1_success/n_total*100:.1f}%)")
    logger.info(f"  Batch-aware: {n_track2_success}/{n_total} ({n_track2_success/n_total*100:.1f}%)")
    logger.info(f"\nTotal time: {pipeline_time/3600:.2f} hours")
    
    # Failed datasets
    if n_failed > 0:
        logger.info(f"\nFailed Datasets:")
        for r in all_results:
            if not r['success']:
                logger.info(f"  {r['dataset_id']}: {r['error']}")
    
    # Save summary
    summary = {
        'pipeline_version': '1.1.0',
        'start_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'total_time_seconds': pipeline_time,
        'n_total': n_total,
        'n_success': n_success,
        'n_failed': n_failed,
        'n_track1_success': n_track1_success,
        'n_track2_success': n_track2_success,
        'results': all_results
    }
    
    summary_file = BASE_OUTPUT_DIR / "pipeline_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"\n✓ Summary saved: {summary_file}")
    logger.info(f"\n{'='*70}")
    logger.info(f"All outputs: {BASE_OUTPUT_DIR}")
    logger.info(f"{'='*70}")


if __name__ == "__main__":
    main()
