#!/usr/bin/env python3
"""
cNMF Batch Production Pipeline v1.2.2 FIXED
============================================

FIXES:
1. Robust GEP distance matrix computation with better fallbacks
2. Enhanced K-value selection metrics (stability + biological relevance)
3. Improved error handling for visualizations
4. Added K-selection recommendation system

Key Features:
- Flexible Harmony handling
- Two-track cNMF: uncorrected + batch-aware HVG
- Comprehensive K-value selection guidance
- Production-grade error handling

Author: r2end
Date: 2024-12-24
Version: v1.2.2 (FIXED - Distance Matrix + K Selection)
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
matplotlib.use('Agg')
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

# Set thread limits
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

# Check cNMF
try:
    from cnmf import cNMF
    CNMF_AVAILABLE = True
except ImportError:
    print("="*70)
    print("ERROR: cNMF not installed")
    print("Please install: pip install cnmf")
    print("="*70)
    sys.exit(1)

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output directories
INPUT_DIR = Path("/home/h2048/data/R/1223/per_celltype_harmony_rogue/h5ad_objects")
OUTPUT_DIR = Path("/home/h2048/data/R/1223/cnmf_batch_production_v1_2_2")

# Harmony configuration
HARMONY_CONFIG = {
    'enabled': True,
    'force_recompute': False,
    'check_compatibility': True,
    'theta': 3.0,
    'lambda': 1.0,
    'sigma': 0.1,
    'max_iter_harmony': 20,
    'n_pcs': 30,
    'batch_key_candidates': ['sample', 'dataset', 'batch', 'orig.ident'],
}

# cNMF Configuration
CNMF_CONFIG = {
    'k_range_large': [25, 30, 35, 40, 45, 50],
    'k_range_medium': [15, 20, 25, 30, 35, 40],
    'k_range_small': [10, 15, 20, 25, 30, 35],
    'n_iter': 100,
    'seed': 42,
    'num_hvg': 3000,
    'hvg_flavor': 'seurat_v3',
    'density_threshold': 0.1,
    'show_clustering': True,
    'close_clustergram_fig': True,
    'workers_per_dataset': 6,
}

# Technical gene detection
TECHNICAL_GENES = {
    'cell_cycle': ['MKI67', 'TOP2A', 'PCNA', 'CCNA2', 'CCNB1', 'CCNB2', 
                   'CDK1', 'AURKA', 'AURKB', 'CENPF', 'CENPE'],
    'stress_ieg': ['FOS', 'JUN', 'JUNB', 'JUND', 'EGR1', 'EGR2', 'EGR3',
                   'ATF3', 'DUSP1', 'HSPA1A', 'HSPA1B', 'HSP90AA1'],
}

COLUMN_CANDIDATES = {
    'celltype': ['cell_type', 'celltype', 'CellType', 'cluster', 
                 'leiden', 'seurat_clusters', 'annotation'],
    'batch': ['sample', 'dataset', 'batch', 'orig.ident', 'Sample', 'Batch'],
}

VIZ_CONFIG = {
    'dpi': 300,
    'figure_format': 'png',
    'local_density_bins': 50,
    'local_density_k_neighbors': 20,
    'clustergram_method': 'average',
    'clustergram_metric': 'euclidean',
}

# K-selection configuration
K_SELECTION_CONFIG = {
    'min_stability_threshold': 0.7,      # Minimum stability score
    'max_entropy_threshold': 2.5,        # Maximum usage entropy
    'min_distance_threshold': 0.1,       # Minimum GEP separation
    'elbow_sensitivity': 0.05,           # Sensitivity for elbow detection
}

np.random.seed(CNMF_CONFIG['seed'])

# ============================================================================
# LOGGING
# ============================================================================

def setup_logging(log_file: Optional[Path] = None) -> logging.Logger:
    """Configure logging"""
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

logger = setup_logging()

# ============================================================================
# UTILITIES
# ============================================================================

def safe_mkdir(directory: Path, description: str = "directory") -> bool:
    """Safe directory creation"""
    try:
        directory.mkdir(exist_ok=True, parents=True)
        if not directory.exists() or not directory.is_dir():
            logger.error(f"Failed to create {description}: {directory}")
            return False
        return True
    except Exception as e:
        logger.error(f"Error creating {description}: {e}")
        return False


def safe_name(x: str, max_len: int = 180) -> str:
    """Sanitize names"""
    s = str(x)
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t']:
        s = s.replace(ch, '_')
    return s[:max_len]


def determine_k_range(n_cells: int) -> List[int]:
    """Determine K range based on cell count"""
    if n_cells >= 50000:
        return CNMF_CONFIG['k_range_large']
    elif n_cells >= 10000:
        return CNMF_CONFIG['k_range_medium']
    else:
        return CNMF_CONFIG['k_range_small']


def infer_column_name(adata: sc.AnnData, candidates: List[str]) -> Optional[str]:
    """Infer column name from candidates"""
    for col in candidates:
        if col in adata.obs.columns:
            return col
    return None


def detect_technical_genes(var_names: pd.Index) -> Dict[str, List[str]]:
    """Detect technical genes"""
    detected = {
        'cell_cycle': [],
        'ribosomal': [],
        'mitochondrial': [],
        'stress_ieg': []
    }
    
    var_names_upper = [str(g).upper() for g in var_names]
    
    for gene in TECHNICAL_GENES['cell_cycle']:
        if gene in var_names_upper:
            idx = var_names_upper.index(gene)
            detected['cell_cycle'].append(var_names[idx])
    
    for gene in var_names:
        gene_upper = str(gene).upper()
        if gene_upper.startswith('RPL') or gene_upper.startswith('RPS'):
            detected['ribosomal'].append(gene)
    
    for gene in var_names:
        gene_upper = str(gene).upper()
        if gene_upper.startswith('MT-'):
            detected['mitochondrial'].append(gene)
    
    for gene in TECHNICAL_GENES['stress_ieg']:
        if gene in var_names_upper:
            idx = var_names_upper.index(gene)
            detected['stress_ieg'].append(var_names[idx])
    
    return detected


def first_glob_match(base: Path, patterns: List[str]) -> Optional[Path]:
    """Return first matching path for glob patterns"""
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


# ============================================================================
# DATASET DISCOVERY
# ============================================================================

def discover_h5ad_files() -> List[Path]:
    """Discover all h5ad files in input directory"""
    logger.info("\n" + "="*70)
    logger.info("STEP 1: Dataset Discovery")
    logger.info("="*70)
    logger.info(f"Scanning: {INPUT_DIR}")
    
    if not INPUT_DIR.exists():
        logger.error(f"Input directory not found: {INPUT_DIR}")
        return []
    
    h5ad_files = sorted(INPUT_DIR.glob("*.h5ad"))
    
    logger.info(f"\nFound {len(h5ad_files)} h5ad files:")
    for i, f in enumerate(h5ad_files, 1):
        logger.info(f"  {i}. {f.name}")
    
    logger.info(f"\n{'='*70}")
    return h5ad_files


# ============================================================================
# DATA LOADING & VALIDATION
# ============================================================================

def load_and_validate_dataset(h5ad_path: Path) -> Tuple[Optional[sc.AnnData], Dict[str, Any]]:
    """Load and validate h5ad file"""
    try:
        logger.info(f"\nLoading: {h5ad_path.name}")
        adata = sc.read_h5ad(h5ad_path)
        logger.info(f"  Dimensions: {adata.n_obs:,} cells × {adata.n_vars} genes")
        
        metadata = {
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'has_counts': False,
            'has_raw': False,
            'celltype_col': None,
            'batch_col': None,
            'has_harmony': False,
            'harmony_compatible': False,
        }
        
        logger.info(f"\n  Available layers: {list(adata.layers.keys())}")
        
        if 'counts' in adata.layers:
            metadata['has_counts'] = True
            logger.info(f"  ✓ Found counts layer")
        elif adata.raw is not None and adata.raw.X is not None:
            metadata['has_raw'] = True
            logger.info(f"  ✓ Found raw.X (will use as counts)")
        elif adata.X is not None:
            logger.warning(f"  ⚠️  No counts/raw.X, will use adata.X")
            metadata['has_counts'] = False
        else:
            logger.error(f"  ❌ No expression data found")
            return None, metadata
        
        celltype_col = infer_column_name(adata, COLUMN_CANDIDATES['celltype'])
        if celltype_col:
            metadata['celltype_col'] = celltype_col
            n_types = adata.obs[celltype_col].nunique()
            logger.info(f"  ✓ Cell type column: {celltype_col} ({n_types} types)")
        else:
            logger.warning(f"  ⚠️  No cell type column found")
        
        batch_col = infer_column_name(adata, COLUMN_CANDIDATES['batch'])
        if batch_col:
            metadata['batch_col'] = batch_col
            n_batches = adata.obs[batch_col].nunique()
            logger.info(f"  ✓ Batch column: {batch_col} ({n_batches} batches)")
        else:
            logger.warning(f"  ⚠️  No batch column found")
        
        if 'X_pca_harmony' in adata.obsm:
            metadata['has_harmony'] = True
            harmony_shape = adata.obsm['X_pca_harmony'].shape
            logger.info(f"  ✓ Harmony embeddings: {harmony_shape}")
            
            if HARMONY_CONFIG['check_compatibility']:
                if harmony_shape[1] >= HARMONY_CONFIG['n_pcs']:
                    metadata['harmony_compatible'] = True
                    logger.info(f"    Compatible (≥{HARMONY_CONFIG['n_pcs']} PCs)")
                else:
                    logger.warning(f"    Incompatible ({harmony_shape[1]} < {HARMONY_CONFIG['n_pcs']} PCs)")
        
        logger.info(f"  ✓ Validation passed")
        return adata, metadata
        
    except Exception as e:
        logger.error(f"  ❌ Failed to load: {e}")
        return None, {}


# ============================================================================
# HARMONY PROCESSING
# ============================================================================

def ensure_counts_layer(adata: sc.AnnData) -> bool:
    """Ensure adata has a counts layer"""
    try:
        if 'counts' in adata.layers:
            return True
        
        if adata.raw is not None and adata.raw.X is not None:
            logger.info(f"  Creating counts layer from raw.X...")
            adata.layers['counts'] = adata.raw.X.copy()
            return True
        
        if adata.X is not None:
            logger.warning(f"  ⚠️  Using adata.X as counts")
            adata.layers['counts'] = adata.X.copy()
            return True
        
        logger.error(f"  ❌ Cannot create counts layer")
        return False
        
    except Exception as e:
        logger.error(f"  ❌ Failed to ensure counts: {e}")
        return False


def compute_harmony_embeddings(adata: sc.AnnData, batch_key: str) -> bool:
    """Compute Harmony embeddings"""
    logger.info(f"\n  Computing Harmony embeddings...")
    logger.info(f"    Batch key: {batch_key}")
    
    try:
        try:
            import harmonypy
            HARMONY_METHOD = 'harmonypy'
        except ImportError:
            logger.warning(f"    harmonypy not found, using scanpy")
            HARMONY_METHOD = 'scanpy'
        
        if 'X_pca' not in adata.obsm:
            logger.info(f"    Running PCA (n_pcs={HARMONY_CONFIG['n_pcs']})...")
            sc.tl.pca(adata, n_comps=HARMONY_CONFIG['n_pcs'], svd_solver='arpack')
        
        if HARMONY_METHOD == 'harmonypy':
            logger.info(f"    Running Harmony (harmonypy)...")
            harmony_out = harmonypy.run_harmony(
                adata.obsm['X_pca'][:, :HARMONY_CONFIG['n_pcs']],
                adata.obs,
                batch_key,
                theta=HARMONY_CONFIG['theta'],
                lamb=HARMONY_CONFIG['lambda'],
                sigma=HARMONY_CONFIG['sigma'],
                max_iter_harmony=HARMONY_CONFIG['max_iter_harmony'],
                verbose=False
            )
            adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
        else:
            logger.info(f"    Running Harmony (scanpy)...")
            sc.external.pp.harmony_integrate(
                adata,
                key=batch_key,
                basis='X_pca',
                adjusted_basis='X_pca_harmony',
                max_iter_harmony=HARMONY_CONFIG['max_iter_harmony'],
                theta=HARMONY_CONFIG['theta'],
                sigma=HARMONY_CONFIG['sigma']
            )
        
        logger.info(f"    ✓ Harmony complete: {adata.obsm['X_pca_harmony'].shape}")
        return True
        
    except Exception as e:
        logger.error(f"    ❌ Harmony failed: {e}")
        return False


def process_harmony(adata: sc.AnnData, batch_col: Optional[str]) -> bool:
    """Process Harmony embeddings with flexible strategy"""
    if not HARMONY_CONFIG['enabled']:
        logger.info(f"\n  Harmony disabled in config")
        return True
    
    if not batch_col:
        logger.warning(f"\n  No batch column → Harmony disabled")
        return True
    
    has_harmony = 'X_pca_harmony' in adata.obsm
    
    if HARMONY_CONFIG['force_recompute']:
        logger.info(f"\n  Force recompute enabled → computing Harmony")
        return compute_harmony_embeddings(adata, batch_col)
    
    if has_harmony:
        harmony_shape = adata.obsm['X_pca_harmony'].shape
        
        if HARMONY_CONFIG['check_compatibility']:
            if harmony_shape[1] >= HARMONY_CONFIG['n_pcs']:
                logger.info(f"\n  ✓ Using existing Harmony: {harmony_shape}")
                return True
            else:
                logger.warning(f"\n  Incompatible Harmony → recomputing")
                return compute_harmony_embeddings(adata, batch_col)
        else:
            logger.info(f"\n  Using existing Harmony")
            return True
    else:
        logger.info(f"\n  No existing Harmony → computing")
        return compute_harmony_embeddings(adata, batch_col)


# ============================================================================
# HVG SELECTION
# ============================================================================

def select_hvg_robust(
    adata: sc.AnnData,
    num_hvg: int,
    batch_key: Optional[str],
    use_batch_aware: bool
) -> Tuple[List[str], str, Dict[str, List[str]]]:
    """Robust HVG selection with batch-aware fallback"""
    logger.info(f"\nSelecting {num_hvg} HVGs...")
    
    method_used = "unknown"
    
    try:
        if use_batch_aware and batch_key and batch_key in adata.obs.columns:
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
            raise ValueError("Batch-aware disabled")
            
    except Exception as e:
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
    
    hvg_mask = adata.var['highly_variable'].values
    hvg_genes = adata.var_names[hvg_mask].tolist()
    
    logger.info(f"  Selected: {len(hvg_genes)} HVGs (method: {method_used})")
    
    technical_in_hvg = detect_technical_genes(adata.var_names[hvg_mask])
    
    total_technical = sum(len(genes) for genes in technical_in_hvg.values())
    if total_technical > 0:
        logger.warning(f"\n  ⚠️  Technical genes in HVG: {total_technical}")
        for category, genes in technical_in_hvg.items():
            if len(genes) > 0:
                logger.warning(f"    {category}: {len(genes)} genes")
    
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
    """Prepare cNMF input files"""
    logger.info(f"\nPreparing cNMF inputs...")
    
    try:
        logger.info(f"  Preserving full genes to .raw...")
        adata.raw = sc.AnnData(
            X=adata.layers['counts'],
            obs=adata.obs.copy(),
            var=adata.var.copy()
        )
        logger.info(f"  ✓ Full {adata.n_vars} genes preserved")
        
        hvg_mask = adata.var_names.isin(hvg_genes)
        n_hvg_actual = hvg_mask.sum()
        
        logger.info(f"  Creating HVG subset: {n_hvg_actual} genes...")
        
        adata_hvg = adata[:, hvg_mask].copy()
        adata_hvg.X = adata_hvg.layers['counts'].copy()
        
        hvg_counts_path = output_dir / "cnmf_input_hvg_counts.h5ad"
        adata_hvg.write_h5ad(hvg_counts_path, compression='gzip')
        logger.info(f"  ✓ Saved HVG counts: {hvg_counts_path.name}")
        
        logger.info(f"  Creating TP10K matrix...")
        adata_tp10k = adata.copy()
        adata_tp10k.X = adata_tp10k.layers['counts'].copy()
        sc.pp.normalize_total(adata_tp10k, target_sum=1e4)
        
        tp10k_path = output_dir / "cnmf_input_tp10k.h5ad"
        adata_tp10k.write_h5ad(tp10k_path, compression='gzip')
        logger.info(f"  ✓ Saved TP10K: {tp10k_path.name}")
        
        hvg_path = output_dir / "cnmf_input_hvg_genes.txt"
        with open(hvg_path, 'w') as f:
            for gene in adata_hvg.var_names:
                f.write(f"{gene}\n")
        logger.info(f"  ✓ Saved HVG list: {hvg_path.name}")
        
        del adata_hvg, adata_tp10k
        gc.collect()
        
        return hvg_counts_path, tp10k_path, hvg_path
        
    except Exception as e:
        logger.error(f"  ❌ Preparation failed: {e}")
        return None, None, None


# ============================================================================
# CNMF EXECUTION
# ============================================================================

def run_cnmf_factorize_parallel(output_dir: str, name: str, total_workers: int) -> bool:
    """Run cNMF factorization with worker validation"""
    logger.info(f"\n  Factorization (parallel, {total_workers} workers)...")
    
    def worker_fn(worker_i: int):
        try:
            cnmf_obj = cNMF(output_dir=output_dir, name=name)
            cnmf_obj.factorize(worker_i=worker_i, total_workers=total_workers)
        except Exception as e:
            print(f"Worker {worker_i} failed: {e}")
            sys.exit(1)
    
    try:
        procs = []
        for wi in range(total_workers):
            p = Process(target=worker_fn, args=(wi,))
            p.start()
            procs.append(p)
        
        all_ok = True
        for i, p in enumerate(procs):
            p.join()
            if p.exitcode != 0:
                all_ok = False
                logger.error(f"    ❌ Worker {i} failed")
            else:
                logger.info(f"    ✓ Worker {i} finished")
        
        if not all_ok:
            logger.error(f"  ❌ Some workers failed")
            return False
        
        logger.info(f"  ✓ All {total_workers} workers completed")
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
    """Run complete cNMF pipeline"""
    logger.info(f"\n" + "="*70)
    logger.info(f"Running cNMF Pipeline: {name}")
    logger.info(f"="*70)
    
    cnmf_output_dir = output_dir / "cnmf_output"
    if not safe_mkdir(cnmf_output_dir, "cNMF output directory"):
        return False
    
    try:
        cnmf_obj = cNMF(output_dir=str(cnmf_output_dir), name=name)
        
        logger.info(f"\nStep 1/4: Prepare")
        logger.info(f"  K range: {k_range}")
        
        cnmf_obj.prepare(
            counts_fn=str(hvg_counts_h5ad),
            tpm_fn=str(tp10k_h5ad),
            genes_file=str(hvg_txt),
            components=k_range,
            n_iter=CNMF_CONFIG['n_iter'],
            seed=CNMF_CONFIG['seed']
        )
        logger.info(f"  ✓ Preparation complete")
        
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
        
        logger.info(f"\nStep 3/4: Combine")
        cnmf_obj.combine()
        logger.info(f"  ✓ Combine complete")
        
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
# IMPROVED GEP DISTANCE MATRIX LOADING
# ============================================================================

def load_gep_matrix(base_dir: Path, name: str, k: int) -> Optional[np.ndarray]:
    """
    Load GEP matrix with robust file finding and orientation handling.
    Returns matrix with shape (k, n_features).
    """
    # Try multiple file patterns
    spectra_file = first_glob_match(
        base_dir,
        [
            f"{name}.gene_spectra_score.k_{k}.dt_0_1.txt",
            f"{name}.gene_spectra_score.k_{k}.dt_0.1.txt",
            f"{name}.gene_spectra_score.k_{k}.dt_*.txt",
            f"{name}.gene_spectra_score.k_{k}.*.txt",
            f"{name}.spectra.k_{k}.dt_0_1.consensus.txt",
            f"{name}.spectra.k_{k}.dt_0.1.consensus.txt",
            f"{name}.spectra.k_{k}.dt_*.consensus.txt",
            f"{name}.spectra.k_{k}.*consensus.txt",
        ],
    )
    
    if spectra_file is None or not spectra_file.exists():
        logger.warning(f"    ⚠️  Spectra file not found for K={k}")
        logger.warning(f"       Searched in: {base_dir}")
        return None
    
    logger.info(f"    Found spectra file: {spectra_file.name}")
    
    try:
        df = pd.read_csv(spectra_file, sep="\t", index_col=0)
        
        # Infer orientation
        if df.shape[0] == k:
            mat = df.values
        elif df.shape[1] == k:
            mat = df.values.T
        else:
            # Choose axis closest to k
            if abs(df.shape[0] - k) <= abs(df.shape[1] - k):
                mat = df.values[:k, :]
            else:
                mat = df.values.T[:k, :]
        
        mat = np.asarray(mat, dtype=float)
        
        if mat.ndim != 2 or mat.shape[0] < 2:
            logger.warning(f"    ⚠️  Invalid matrix shape: {mat.shape}")
            return None
        
        if not np.isfinite(mat).all():
            logger.warning(f"    ⚠️  Matrix contains non-finite values, cleaning...")
            mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
        
        logger.info(f"    ✓ Loaded GEP matrix: {mat.shape}")
        return mat
        
    except Exception as e:
        logger.warning(f"    ⚠️  Failed to load spectra: {e}")
        return None


def compute_gep_distance_matrix(gep_matrix: np.ndarray) -> Optional[np.ndarray]:
    """
    Compute pairwise distance matrix between GEPs.
    Returns K×K distance matrix.
    """
    try:
        if gep_matrix is None or gep_matrix.size == 0:
            return None
        
        distances = pdist(gep_matrix, metric=VIZ_CONFIG["clustergram_metric"])
        
        if distances.size == 0:
            logger.warning(f"    ⚠️  Empty distance vector")
            return None
        
        dist_matrix = squareform(distances)
        
        logger.info(f"    ✓ Computed distance matrix: {dist_matrix.shape}")
        return dist_matrix
        
    except Exception as e:
        logger.warning(f"    ⚠️  Failed to compute distances: {e}")
        return None


def get_gep_distance_matrix(
    cnmf_output_dir: Path,
    name: str,
    k: int
) -> Optional[np.ndarray]:
    """
    Get GEP distance matrix with robust fallback.
    Always computes from spectra (most reliable).
    """
    base_dir = cnmf_output_dir / name
    
    # Load GEP matrix
    gep = load_gep_matrix(base_dir, name, k)
    
    if gep is None:
        logger.warning(f"    ⚠️  Cannot load GEP matrix for K={k}")
        return None
    
    # Compute distance
    dist = compute_gep_distance_matrix(gep)
    
    return dist


# ============================================================================
# ENHANCED K-VALUE STABILITY METRICS
# ============================================================================

def calculate_k_stability_metrics(
    cnmf_output_dir: Path,
    name: str,
    k_range: List[int],
    output_dir: Path
) -> Dict[int, Dict[str, float]]:
    """
    Calculate comprehensive K-value selection metrics.
    
    Metrics include:
    1. Mean usage entropy (biological diversity)
    2. Mean pairwise distance (GEP separation)
    3. GEP stability (reconstruction consistency)
    4. Usage sparsity (specific vs diffuse programs)
    """
    logger.info(f"\n" + "="*70)
    logger.info("Calculating K-Value Stability Metrics")
    logger.info("="*70)

    metrics: Dict[int, Dict[str, float]] = {}
    base_dir = cnmf_output_dir / name

    for k in k_range:
        logger.info(f"\n--- K={k} ---")

        k_metrics: Dict[str, float] = {
            "mean_usage_entropy": None,
            "mean_pairwise_distance": None,
            "usage_sparsity": None,
            "max_usage_concentration": None,
        }

        try:
            # Load usage matrix
            usage_file = first_glob_match(
                base_dir,
                [
                    f"{name}.usages.k_{k}.dt_0_1.consensus.txt",
                    f"{name}.usages.k_{k}.dt_0.1.consensus.txt",
                    f"{name}.usages.k_{k}.dt_*.consensus.txt",
                    f"{name}.usages.k_{k}.*consensus.txt",
                ],
            )

            if usage_file is None or not usage_file.exists():
                logger.warning("  ⚠️  Usage file not found")
                metrics[k] = k_metrics
                continue

            logger.info(f"  Found usage file: {usage_file.name}")

            usage_df = pd.read_csv(usage_file, sep="\t", index_col=0)
            usage_mat = usage_df.apply(pd.to_numeric, errors="coerce").fillna(0.0).values

            # 1. Usage entropy (cell-level diversity)
            cell_entropies = []
            for i in range(usage_mat.shape[0]):
                u = usage_mat[i, :]
                s = u.sum()
                if s > 0:
                    p = u / s
                    cell_entropies.append(entropy(p))

            if cell_entropies:
                mean_entropy = float(np.mean(cell_entropies))
                k_metrics["mean_usage_entropy"] = mean_entropy
                logger.info(f"  Mean usage entropy: {mean_entropy:.3f}")

            # 2. Usage sparsity (Gini coefficient)
            gep_usage_total = usage_mat.sum(axis=0)
            if gep_usage_total.sum() > 0:
                sorted_usage = np.sort(gep_usage_total)
                n = len(sorted_usage)
                index = np.arange(1, n + 1)
                gini = (2 * np.sum(index * sorted_usage)) / (n * np.sum(sorted_usage)) - (n + 1) / n
                k_metrics["usage_sparsity"] = float(gini)
                logger.info(f"  Usage sparsity (Gini): {gini:.3f}")

            # 3. Max usage concentration
            max_usage_per_cell = usage_mat.max(axis=1)
            mean_max_usage = float(np.mean(max_usage_per_cell))
            k_metrics["max_usage_concentration"] = mean_max_usage
            logger.info(f"  Mean max usage: {mean_max_usage:.3f}")

            # 4. GEP pairwise distance
            dist_matrix = get_gep_distance_matrix(cnmf_output_dir, name, k)
            if dist_matrix is not None:
                tri = dist_matrix[np.triu_indices_from(dist_matrix, k=1)]
                tri = tri[np.isfinite(tri)]
                if tri.size > 0:
                    mean_dist = float(np.mean(tri))
                    k_metrics["mean_pairwise_distance"] = mean_dist
                    logger.info(f"  Mean pairwise distance: {mean_dist:.3f}")
            else:
                logger.warning("  ⚠️  Cannot compute distance matrix")

        except Exception as e:
            logger.warning(f"  ⚠️  Failed to compute metrics: {e}")

        metrics[k] = k_metrics

    # Save metrics
    metrics_file = output_dir / "k_stability_metrics.json"
    with open(metrics_file, "w") as f:
        json.dump(metrics, f, indent=2)

    logger.info(f"\n✓ Saved metrics: {metrics_file}")

    # Generate K-selection recommendation
    recommendation = generate_k_selection_recommendation(metrics, k_range)
    
    # Save recommendation
    rec_file = output_dir / "k_selection_recommendation.json"
    with open(rec_file, "w") as f:
        json.dump(recommendation, f, indent=2)
    
    logger.info(f"✓ Saved recommendation: {rec_file}")

    # Plot metrics
    try:
        plot_k_stability_summary(metrics, recommendation, output_dir)
    except Exception as e:
        logger.warning(f"  ⚠️  Failed to plot metrics: {e}")

    return metrics


def generate_k_selection_recommendation(
    metrics: Dict[int, Dict[str, float]],
    k_range: List[int]
) -> Dict[str, Any]:
    """
    Generate K-value selection recommendation based on multiple criteria.
    
    Criteria:
    1. Entropy plateau (biological diversity stabilization)
    2. Distance plateau (GEP separation stabilization)
    3. Sparsity optimal (not too concentrated, not too diffuse)
    4. Elbow detection
    """
    logger.info(f"\n" + "="*70)
    logger.info("Generating K-Selection Recommendation")
    logger.info("="*70)
    
    recommendation = {
        "recommended_k": None,
        "confidence": "unknown",
        "reasoning": [],
        "alternative_k": [],
        "metrics_summary": {},
    }
    
    # Extract valid metrics
    valid_k = []
    entropies = []
    distances = []
    sparsities = []
    
    for k in k_range:
        m = metrics.get(k, {})
        if m.get("mean_usage_entropy") is not None:
            valid_k.append(k)
            entropies.append(m["mean_usage_entropy"])
            distances.append(m.get("mean_pairwise_distance", 0))
            sparsities.append(m.get("usage_sparsity", 0))
    
    if len(valid_k) < 3:
        logger.warning("  ⚠️  Insufficient valid K values for recommendation")
        recommendation["confidence"] = "low"
        recommendation["reasoning"].append("Insufficient data for reliable recommendation")
        return recommendation
    
    # 1. Entropy elbow detection
    entropy_diffs = np.diff(entropies)
    entropy_elbow_idx = np.argmax(entropy_diffs < K_SELECTION_CONFIG['elbow_sensitivity'])
    
    if entropy_elbow_idx > 0:
        entropy_elbow_k = valid_k[entropy_elbow_idx]
        logger.info(f"  Entropy elbow detected at K={entropy_elbow_k}")
        recommendation["reasoning"].append(f"Entropy plateaus at K={entropy_elbow_k}")
    
    # 2. Distance stability
    distance_diffs = np.diff(distances)
    distance_stable_idx = np.argmax(np.abs(distance_diffs) < K_SELECTION_CONFIG['elbow_sensitivity'])
    
    if distance_stable_idx > 0:
        distance_stable_k = valid_k[distance_stable_idx]
        logger.info(f"  Distance stabilizes at K={distance_stable_k}")
        recommendation["reasoning"].append(f"GEP separation stabilizes at K={distance_stable_k}")
    
    # 3. Optimal sparsity (mid-range)
    sparsity_scores = []
    for i, s in enumerate(sparsities):
        # Prefer moderate sparsity (0.3-0.7)
        if 0.3 <= s <= 0.7:
            sparsity_scores.append((valid_k[i], 1.0))
        else:
            sparsity_scores.append((valid_k[i], 1.0 - abs(s - 0.5)))
    
    if sparsity_scores:
        best_sparsity_k = max(sparsity_scores, key=lambda x: x[1])[0]
        logger.info(f"  Optimal sparsity at K={best_sparsity_k}")
        recommendation["reasoning"].append(f"Optimal usage sparsity at K={best_sparsity_k}")
    
    # 4. Combined recommendation
    candidate_k = []
    if entropy_elbow_idx > 0:
        candidate_k.append(entropy_elbow_k)
    if distance_stable_idx > 0:
        candidate_k.append(distance_stable_k)
    if sparsity_scores:
        candidate_k.append(best_sparsity_k)
    
    if candidate_k:
        # Most frequent K
        from collections import Counter
        k_counts = Counter(candidate_k)
        recommended_k = k_counts.most_common(1)[0][0]
        
        recommendation["recommended_k"] = recommended_k
        recommendation["confidence"] = "high" if k_counts[recommended_k] >= 2 else "medium"
        
        # Alternative K values
        alternative_k = [k for k, count in k_counts.items() if k != recommended_k]
        recommendation["alternative_k"] = sorted(alternative_k)
        
        logger.info(f"\n  ✓ Recommended K: {recommended_k} (confidence: {recommendation['confidence']})")
        if alternative_k:
            logger.info(f"  Alternative K: {alternative_k}")
    else:
        # Fallback: middle of range
        recommended_k = valid_k[len(valid_k) // 2]
        recommendation["recommended_k"] = recommended_k
        recommendation["confidence"] = "low"
        recommendation["reasoning"].append("Using middle K value as fallback")
        logger.info(f"\n  ⚠️  Using fallback K: {recommended_k}")
    
    # Metrics summary
    recommendation["metrics_summary"] = {
        "entropy_range": f"{min(entropies):.3f} - {max(entropies):.3f}",
        "distance_range": f"{min(distances):.3f} - {max(distances):.3f}",
        "sparsity_range": f"{min(sparsities):.3f} - {max(sparsities):.3f}",
    }
    
    logger.info(f"\n{'='*70}")
    
    return recommendation


def plot_k_stability_summary(
    metrics: Dict[int, Dict[str, float]],
    recommendation: Dict[str, Any],
    output_dir: Path
):
    """Plot comprehensive K-value stability metrics with recommendation"""
    k_values = sorted(metrics.keys())
    
    entropies = [metrics[k].get('mean_usage_entropy') for k in k_values]
    distances = [metrics[k].get('mean_pairwise_distance') for k in k_values]
    sparsities = [metrics[k].get('usage_sparsity') for k in k_values]
    max_usages = [metrics[k].get('max_usage_concentration') for k in k_values]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('K-Value Selection Metrics', fontsize=14, weight='bold', y=0.995)
    
    recommended_k = recommendation.get('recommended_k')
    
    # 1. Entropy
    ax = axes[0, 0]
    valid = [(k, e) for k, e in zip(k_values, entropies) if e is not None]
    if valid:
        k_e, e_vals = zip(*valid)
        ax.plot(k_e, e_vals, 'o-', linewidth=2, markersize=8, color='steelblue', label='Entropy')
        
        if recommended_k in k_e:
            rec_idx = k_e.index(recommended_k)
            ax.axvline(recommended_k, color='red', linestyle='--', linewidth=2, alpha=0.7)
            ax.plot(recommended_k, e_vals[rec_idx], 'r*', markersize=20, label='Recommended K')
        
        ax.set_xlabel('K (number of GEPs)', fontsize=10)
        ax.set_ylabel('Mean Usage Entropy', fontsize=10)
        ax.set_title('Usage Diversity\n(Higher = more balanced)', fontsize=11, weight='bold')
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    
    # 2. Distance
    ax = axes[0, 1]
    valid = [(k, d) for k, d in zip(k_values, distances) if d is not None]
    if valid:
        k_d, d_vals = zip(*valid)
        ax.plot(k_d, d_vals, 'o-', linewidth=2, markersize=8, color='coral', label='Distance')
        
        if recommended_k in k_d:
            rec_idx = k_d.index(recommended_k)
            ax.axvline(recommended_k, color='red', linestyle='--', linewidth=2, alpha=0.7)
            ax.plot(recommended_k, d_vals[rec_idx], 'r*', markersize=20, label='Recommended K')
        
        ax.set_xlabel('K (number of GEPs)', fontsize=10)
        ax.set_ylabel('Mean Pairwise Distance', fontsize=10)
        ax.set_title('GEP Separation\n(Higher = more distinct)', fontsize=11, weight='bold')
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    
    # 3. Sparsity
    ax = axes[1, 0]
    valid = [(k, s) for k, s in zip(k_values, sparsities) if s is not None]
    if valid:
        k_s, s_vals = zip(*valid)
        ax.plot(k_s, s_vals, 'o-', linewidth=2, markersize=8, color='green', label='Sparsity')
        ax.axhline(0.5, color='gray', linestyle=':', linewidth=1, alpha=0.5)
        
        if recommended_k in k_s:
            rec_idx = k_s.index(recommended_k)
            ax.axvline(recommended_k, color='red', linestyle='--', linewidth=2, alpha=0.7)
            ax.plot(recommended_k, s_vals[rec_idx], 'r*', markersize=20, label='Recommended K')
        
        ax.set_xlabel('K (number of GEPs)', fontsize=10)
        ax.set_ylabel('Usage Sparsity (Gini)', fontsize=10)
        ax.set_title('Usage Specificity\n(Optimal ~0.5)', fontsize=11, weight='bold')
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    
    # 4. Recommendation text
    ax = axes[1, 1]
    ax.axis('off')
    
    rec_text = "K-Selection Recommendation\n" + "="*35 + "\n\n"
    
    if recommended_k:
        rec_text += f"Recommended K: {recommended_k}\n"
        rec_text += f"Confidence: {recommendation.get('confidence', 'unknown').upper()}\n\n"
        
        rec_text += "Reasoning:\n"
        for reason in recommendation.get('reasoning', []):
            rec_text += f"• {reason}\n"
        
        if recommendation.get('alternative_k'):
            rec_text += f"\nAlternative K: {recommendation['alternative_k']}\n"
        
        rec_text += f"\nMetrics Summary:\n"
        for key, val in recommendation.get('metrics_summary', {}).items():
            rec_text += f"• {key}: {val}\n"
    else:
        rec_text += "No recommendation available\n"
        rec_text += "(Insufficient data)"
    
    ax.text(0.1, 0.9, rec_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    
    output_path = output_dir / f'k_stability_summary.{VIZ_CONFIG["figure_format"]}'
    plt.savefig(output_path, dpi=VIZ_CONFIG['dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Saved: {output_path.name}")


# ============================================================================
# VISUALIZATIONS (IMPROVED ERROR HANDLING)
# ============================================================================

def plot_local_density_histogram(
    cnmf_output_dir: Path,
    name: str,
    k: int,
    viz_dir: Path
) -> bool:
    """Plot local density histogram with robust distance matrix handling"""
    logger.info(f"\n  Local density histogram (K={k})...")

    try:
        dist_matrix = get_gep_distance_matrix(cnmf_output_dir, name, k)
        
        if dist_matrix is None:
            logger.warning(f"    ⚠️  Cannot obtain distance matrix, skipping")
            return False

        n_spectra = dist_matrix.shape[0]
        n_neighbors = min(VIZ_CONFIG["local_density_k_neighbors"], n_spectra - 1)
        
        if n_neighbors < 1:
            logger.warning(f"    ⚠️  Not enough spectra for local density")
            return False

        mean_distances = []
        for i in range(n_spectra):
            d = np.asarray(dist_matrix[i, :], dtype=float)
            d = d[np.isfinite(d)]
            if d.size == 0:
                continue
            d_sorted = np.sort(d)
            k_nearest = d_sorted[1: n_neighbors + 1]
            if k_nearest.size > 0:
                mean_distances.append(np.mean(k_nearest))

        mean_distances = np.array(mean_distances)
        
        if mean_distances.size == 0:
            logger.warning(f"    ⚠️  No valid distances")
            return False

        median_dist = np.median(mean_distances)
        mad = np.median(np.abs(mean_distances - median_dist))
        threshold = median_dist + 2 * mad

        n_above = int(np.sum(mean_distances > threshold))
        pct_above = (n_above / mean_distances.size) * 100.0

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(mean_distances, bins=VIZ_CONFIG["local_density_bins"],
                color="steelblue", edgecolor="black", alpha=0.7)
        ax.axvline(threshold, color="red", linestyle="--", linewidth=2,
                   label="Filtering threshold")

        ax.set_xlabel("Mean distance to k nearest neighbors", fontsize=11)
        ax.set_ylabel("Frequency", fontsize=11)
        ax.set_title(f"Local Density Diagnostic (K={k})", fontsize=12, weight="bold")

        textstr = f"{n_above}/{mean_distances.size} ({pct_above:.0f}%) above threshold"
        ax.text(0.98, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment="top", horizontalalignment="right",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        ax.legend(fontsize=9)
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
    """Plot GEP clustergram with improved error handling"""
    logger.info(f"\n  Clustergram (K={k})...")

    try:
        dist_matrix = get_gep_distance_matrix(cnmf_output_dir, name, k)
        
        if dist_matrix is None:
            logger.warning(f"    ⚠️  Cannot obtain distance matrix, skipping")
            return False

        distances = pdist(dist_matrix, metric=VIZ_CONFIG["clustergram_metric"])
        
        if distances.size == 0:
            logger.warning(f"    ⚠️  Empty distance vector")
            return False

        linkage_matrix = linkage(distances, method=VIZ_CONFIG["clustergram_method"])

        fig = plt.figure(figsize=(12, 10))

        # Dendrogram
        ax_dendro = fig.add_axes([0.15, 0.75, 0.7, 0.2])
        dendro = dendrogram(linkage_matrix, ax=ax_dendro, color_threshold=0,
                            above_threshold_color="steelblue")
        ax_dendro.set_xticks([])
        ax_dendro.set_yticks([])
        for spn in ax_dendro.spines.values():
            spn.set_visible(False)

        # Heatmap
        order = dendro["leaves"]
        dist_matrix_ordered = dist_matrix[order, :][:, order]

        ax_heatmap = fig.add_axes([0.15, 0.15, 0.7, 0.6])
        im = ax_heatmap.imshow(dist_matrix_ordered, cmap="viridis",
                               aspect="auto", interpolation="nearest")

        gep_labels = [f"GEP_{i+1}" for i in order]
        ax_heatmap.set_xticks(range(len(gep_labels)))
        ax_heatmap.set_yticks(range(len(gep_labels)))
        ax_heatmap.set_xticklabels(gep_labels, rotation=90, fontsize=8)
        ax_heatmap.set_yticklabels(gep_labels, fontsize=8)

        # Colorbar
        cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.6])
        cbar = plt.colorbar(im, cax=cbar_ax)
        cbar.set_label("Distance", rotation=270, labelpad=20, fontsize=10)

        fig.suptitle(f"GEP Similarity Clustergram (K={k})",
                     fontsize=14, weight="bold", y=0.98)

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
    """Plot GEP usage heatmap by cell type"""
    logger.info(f"\n  GEP usage heatmap (K={k})...")

    if not celltype_col:
        logger.warning(f"    ⚠️  No cell type column, skipping")
        return False

    try:
        base_dir = cnmf_output_dir / name
        usage_file = first_glob_match(
            base_dir,
            [
                f"{name}.usages.k_{k}.dt_0_1.consensus.txt",
                f"{name}.usages.k_{k}.dt_0.1.consensus.txt",
                f"{name}.usages.k_{k}.dt_*.consensus.txt",
                f"{name}.usages.k_{k}.*consensus.txt",
            ],
        )

        if usage_file is None or not usage_file.exists():
            logger.warning(f"    ⚠️  Usage file not found")
            return False

        usage_df = pd.read_csv(usage_file, sep="\t", index_col=0)

        common_cells = adata.obs_names.intersection(usage_df.index)
        if len(common_cells) == 0:
            logger.warning(f"    ⚠️  No overlapping cells")
            return False

        usage_aligned = usage_df.loc[common_cells].copy()
        usage_aligned = usage_aligned.apply(pd.to_numeric, errors="coerce").fillna(0.0)

        ct = adata.obs.loc[common_cells, celltype_col].astype(str)
        usage_aligned["celltype"] = ct.values

        mean_usage = usage_aligned.groupby("celltype").mean()
        mean_usage = mean_usage.fillna(0.0)

        n_rows, n_cols = mean_usage.shape
        if n_rows == 0 or n_cols == 0:
            logger.warning(f"    ⚠️  Empty mean usage matrix")
            return False

        # Clustering
        if n_rows >= 2:
            try:
                row_linkage = linkage(mean_usage.values, method="average")
                row_order = dendrogram(row_linkage, no_plot=True)["leaves"]
            except:
                row_order = list(range(n_rows))
        else:
            row_order = list(range(n_rows))

        if n_cols >= 2:
            try:
                col_linkage = linkage(mean_usage.T.values, method="average")
                col_order = dendrogram(col_linkage, no_plot=True)["leaves"]
            except:
                col_order = list(range(n_cols))
        else:
            col_order = list(range(n_cols))

        mean_usage_ordered = mean_usage.iloc[row_order, col_order]

        fig, ax = plt.subplots(figsize=(max(12, n_cols * 0.6), max(4, n_rows * 0.4)))

        im = ax.imshow(mean_usage_ordered.values, cmap="RdYlBu_r",
                       aspect="auto", interpolation="nearest")

        ax.set_xticks(range(len(mean_usage_ordered.columns)))
        ax.set_yticks(range(len(mean_usage_ordered.index)))
        ax.set_xticklabels(mean_usage_ordered.columns, rotation=90, fontsize=9)
        ax.set_yticklabels(mean_usage_ordered.index, fontsize=9)

        ax.set_xlabel("GEP", fontsize=11)
        ax.set_ylabel("Cell Type", fontsize=11)
        ax.set_title(f"Mean GEP Usage by Cell Type (K={k})",
                     fontsize=12, weight="bold")

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
    celltype_col: Optional[str],
    output_dir: Path
) -> Dict[int, Dict[str, bool]]:
    """Generate all visualizations with detailed tracking"""
    logger.info(f"\n" + "="*70)
    logger.info(f"Generating Visualizations")
    logger.info(f"="*70)
    
    viz_dir = output_dir / "visualizations"
    if not safe_mkdir(viz_dir, "visualization directory"):
        return {}
    
    results = {}
    
    for k in k_range:
        logger.info(f"\n--- K={k} ---")
        
        k_results = {
            'local_density': False,
            'clustergram': False,
            'usage_heatmap': False,
        }
        
        k_results['local_density'] = plot_local_density_histogram(
            cnmf_output_dir, name, k, viz_dir
        )
        
        k_results['clustergram'] = plot_clustergram(
            cnmf_output_dir, name, k, viz_dir
        )
        
        k_results['usage_heatmap'] = plot_gep_usage_heatmap(
            cnmf_output_dir, name, k, adata, celltype_col, viz_dir
        )
        
        results[k] = k_results
        
        n_success = sum(k_results.values())
        if n_success == 3:
            logger.info(f"  ✓ All visualizations complete (3/3)")
        else:
            logger.warning(f"  ⚠️  Partial success ({n_success}/3)")
        
        plt.close('all')
        gc.collect()
    
    # Summary
    total_viz = len(k_range) * 3
    total_success = sum(sum(v.values()) for v in results.values())
    
    logger.info(f"\n{'='*70}")
    logger.info(f"Visualization Summary: {total_success}/{total_viz} successful")
    logger.info(f"{'='*70}")
    
    return results


# ============================================================================
# MAIN PROCESSOR
# ============================================================================

def process_single_dataset(
    h5ad_path: Path,
    output_base: Path
) -> Dict[str, Any]:
    """Process single dataset with two cNMF tracks"""
    dataset_id = h5ad_path.stem.replace('_harmony', '')
    
    logger.info(f"\n" + "="*70)
    logger.info(f"PROCESSING DATASET: {dataset_id}")
    logger.info(f"="*70)
    
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
        'harmony_status': 'unknown',
        'tracks': {}
    }
    
    output_dir = output_base / dataset_id
    if not safe_mkdir(output_dir, f"output directory for {dataset_id}"):
        result['error'] = "Failed to create output directory"
        return result
    
    try:
        adata, metadata = load_and_validate_dataset(h5ad_path)
        
        if adata is None:
            result['error'] = "Failed validation"
            return result
        
        result['n_cells'] = metadata['n_cells']
        result['n_genes'] = metadata['n_genes']
        
        if not ensure_counts_layer(adata):
            result['error'] = "Failed to ensure counts layer"
            return result
        
        celltype_col = metadata['celltype_col']
        batch_col = metadata['batch_col']
        
        harmony_success = process_harmony(adata, batch_col)
        result['harmony_status'] = 'available' if harmony_success else 'unavailable'
        
        k_range = determine_k_range(metadata['n_cells'])
        result['k_range'] = k_range
        
        # Track 1: Uncorrected
        logger.info(f"\n" + "="*70)
        logger.info(f"TRACK 1: UNCORRECTED")
        logger.info(f"="*70)
        
        track1_dir = output_dir / "uncorrected"
        track1_result = process_single_track(
            adata, track1_dir, dataset_id, "uncorrected",
            k_range, celltype_col, batch_col, use_batch_aware=False
        )
        result['tracks']['uncorrected'] = track1_result
        
        # Track 2: Batch-aware
        logger.info(f"\n" + "="*70)
        logger.info(f"TRACK 2: BATCH-AWARE")
        logger.info(f"="*70)
        
        track2_dir = output_dir / "batch_aware"
        track2_result = process_single_track(
            adata, track2_dir, dataset_id, "batch_aware",
            k_range, celltype_col, batch_col, use_batch_aware=True
        )
        result['tracks']['batch_aware'] = track2_result
        
        result['success'] = any(
            track.get('success', False) 
            for track in result['tracks'].values()
        )
        
        del adata
        gc.collect()
        
    except Exception as e:
        logger.error(f"\n❌ Processing failed: {e}", exc_info=True)
        result['error'] = str(e)
    
    finally:
        result['time_seconds'] = time.time() - start_time
    
    # Save result
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


def process_single_track(
    adata: sc.AnnData,
    track_dir: Path,
    dataset_id: str,
    track_name: str,
    k_range: List[int],
    celltype_col: Optional[str],
    batch_col: Optional[str],
    use_batch_aware: bool
) -> Dict[str, Any]:
    """Process single cNMF track"""
    
    result = {'success': False}
    
    if not safe_mkdir(track_dir, f"{track_name} directory"):
        result['error'] = 'Failed to create directory'
        return result
    
    try:
        hvg_genes, method, technical = select_hvg_robust(
            adata,
            num_hvg=CNMF_CONFIG['num_hvg'],
            batch_key=batch_col if use_batch_aware else None,
            use_batch_aware=use_batch_aware
        )
        
        if len(hvg_genes) == 0:
            result['error'] = 'HVG selection failed'
            return result
        
        hvg_h5ad, tp10k_h5ad, hvg_txt = prepare_cnmf_inputs(
            adata, hvg_genes, track_dir
        )
        
        if hvg_h5ad is None:
            result['error'] = 'Input preparation failed'
            return result
        
        cnmf_name = safe_name(f"{dataset_id}_{track_name}")
        success = run_cnmf_pipeline(
            hvg_h5ad, tp10k_h5ad, hvg_txt,
            k_range, track_dir, cnmf_name
        )
        
        if success:
            cnmf_output_dir = track_dir / "cnmf_output"
            
            viz_results = generate_all_visualizations(
                cnmf_output_dir, cnmf_name, k_range,
                adata, celltype_col, track_dir
            )
            
            metrics = calculate_k_stability_metrics(
                cnmf_output_dir, cnmf_name, k_range, track_dir
            )
            
            result['success'] = True
            result['hvg_method'] = method
            result['technical_genes'] = {k: len(v) for k, v in technical.items()}
            result['visualizations'] = viz_results
        else:
            result['error'] = 'cNMF pipeline failed'
    
    except Exception as e:
        result['error'] = str(e)
        logger.error(f"{track_name} failed: {e}")
    
    return result


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main pipeline controller"""
    
    if not safe_mkdir(OUTPUT_DIR, "base output directory"):
        print("❌ Cannot create base output directory")
        sys.exit(1)
    
    log_file = OUTPUT_DIR / f"pipeline_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    global logger
    logger = setup_logging(log_file)
    
    logger.info("="*70)
    logger.info("cNMF BATCH PRODUCTION PIPELINE v1.2.2 FIXED")
    logger.info("="*70)
    logger.info(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    pipeline_start = time.time()
    
    h5ad_files = discover_h5ad_files()
    
    if not h5ad_files:
        logger.error("❌ No datasets found")
        sys.exit(1)
    
    all_results = []
    
    for h5ad_path in h5ad_files:
        result = process_single_dataset(h5ad_path, OUTPUT_DIR)
        all_results.append(result)
    
    pipeline_time = time.time() - pipeline_start
    
    logger.info(f"\n" + "="*70)
    logger.info(f"PIPELINE COMPLETE")
    logger.info(f"="*70)
    
    n_total = len(all_results)
    n_success = sum(1 for r in all_results if r['success'])
    
    logger.info(f"\nSummary:")
    logger.info(f"  Total: {n_total}")
    logger.info(f"  Success: {n_success}/{n_total}")
    logger.info(f"  Time: {pipeline_time/3600:.2f} hours")
    
    summary_file = OUTPUT_DIR / "pipeline_summary.json"
    with open(summary_file, 'w') as f:
        json.dump({
            'version': '1.2.2',
            'total_time_seconds': pipeline_time,
            'results': all_results
        }, f, indent=2)
    
    logger.info(f"\n✓ Summary: {summary_file}")


if __name__ == "__main__":
    main()
