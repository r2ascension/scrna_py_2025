#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Label-Guided cNMF Pipeline v1.1 HOTFIX
=======================================

CRITICAL FIXES from code review:
- P0-1: Fixed factorize() parallel logic (single worker mode)
- P0-2: Complete safe_name() sanitization for all paths
- P0-3: Sparse matrix streaming write (no dense conversion)
- P1-4: Strict counts validation from layers['counts']
- P1-5: Robust "skip if exists" check for all K values
- P1-6: num_hvg upper bound protection
- P1-7: K range cap for small cell types
- P1-8: Robust dt_0.1/dt_0_1 file matching

Pipeline for running cNMF on cell type-specific subsets with gene filtering

Author: r2end
Date: 2025-01-14
Version: 1.1 - HOTFIX for production stability
"""

import os
import sys
import gc
import re
import json
import warnings
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import sparse

warnings.filterwarnings('ignore')

# Check cNMF availability
try:
    from cnmf import cNMF
    CNMF_AVAILABLE = True
except ImportError:
    print("="*80)
    print("ERROR: cNMF not installed")
    print("Please install: pip install cnmf")
    print("="*80)
    sys.exit(1)

# Set single-threaded mode for cNMF
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input file
INPUT_H5AD = "/home/h2048/data/py/0110/celltypist_epithelial/epithelial_scanvi_final.h5ad"

# Output directory
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/cnmf_by_celltype_v1_1")

# ===== Cell Type Selection =====
CELLTYPE_KEY = 'scanvi_predictions'  # Column name containing cell type labels

# Specify which cell types to process
PROCESS_CELLTYPES = 'auto'  # 'auto' or list like ['Basal', 'Goblet', 'Ciliated']
MIN_CELLS_PER_TYPE = 500     # Minimum cells required to run cNMF

# ⭐ NEW: Confidence gating (prevent label errors from contaminating GEPs)
USE_CONFIDENCE_FILTER = True
CONFIDENCE_KEY = 'scanvi_confidence'  # Optional: filter low-confidence cells
MIN_CONFIDENCE = 0.5  # Only use cells with confidence > 0.5

# ===== Gene Filtering Configuration =====
GENE_FILTER_CONFIG = {
    'remove_mt': True,           # Remove MT- genes
    'remove_ribo': True,         # Remove RPS/RPL/MRPS/MRPL
    'remove_histone': True,      # Remove H1/H2A/H2B/H3/H4/HIST
    'remove_pseudogenes': True,  # Remove *P* genes
    'remove_ensg': True,         # Remove ENSG* genes
    'remove_unannotated': True,  # Remove AC/AL/RP/CTD/LINC/LOC/etc
}

# ===== cNMF Configuration =====
CNMF_CONFIG = {
    # K range (number of GEPs to infer)
    'k_range_auto': True,  # Automatically determine based on cell count
    'k_range_custom': None,  # Or specify: [10, 15, 20, 25, 30]
    
    # K range rules (if auto) - ⭐ FIXED: Lower cap for small cell types
    'k_rules': {
        'large': (10000, [25, 30, 35, 40, 45, 50]),   # >10k cells
        'medium': (2000, [15, 20, 25, 30, 35]),       # 2k-10k cells
        'small': (500, [10, 15, 20, 25]),             # 500-2k cells (capped at 25)
    },
    
    # cNMF parameters
    'n_iter': 100,              # Number of NMF iterations
    'seed': 42,                 # Random seed
    'num_hvg': 3000,            # Number of highly variable genes
    'hvg_flavor': 'seurat_v3',  # HVG selection method
    'density_threshold': 0.1,   # Minimum usage density
    
    # ⭐ FIXED: Single worker mode (or set to 1 for production)
    'total_workers': 1,         # Use 1 for single-worker mode
    'use_multiprocessing': False,  # Set True if you want parallel workers
    
    # Visualization
    'show_clustering': True,
    'close_clustergram_fig': True,
}

# ===== Processing Options =====
SAVE_SUBSET_H5AD = True      # Save filtered subset h5ad for each cell type
SKIP_IF_EXISTS = True        # Skip if cNMF output already exists
OVERWRITE_GENE_FILTER = False  # Force re-filter even if filtered h5ad exists

# Random seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("Label-Guided cNMF Pipeline v1.1 HOTFIX")
print("="*80)
print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Input file: {INPUT_H5AD}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Cell type key: {CELLTYPE_KEY}")
print(f"⭐ Using single-worker mode for cNMF")
print("="*80)


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def safe_name(x: str, max_len: int = 180) -> str:
    """
    Sanitize names for file paths and cNMF identifiers.
    
    Parameters
    ----------
    x : str
        Original name
    max_len : int
        Maximum length
    
    Returns
    -------
    sanitized : str
        Safe name without special characters
    """
    s = str(x)
    # Remove all problematic characters
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t', '(', ')', '[', ']', '{', '}']:
        s = s.replace(ch, '_')
    # Remove multiple consecutive underscores
    while '__' in s:
        s = s.replace('__', '_')
    # Remove leading/trailing underscores
    s = s.strip('_')
    return s[:max_len]


def validate_counts_matrix(X, layer_name: str = "X") -> Tuple[bool, str]:
    """
    Validate that matrix contains raw counts.
    
    Parameters
    ----------
    X : array-like
        Matrix to validate
    layer_name : str
        Name of the layer (for error messages)
    
    Returns
    -------
    is_valid : bool
        Whether matrix is valid counts
    message : str
        Validation message
    """
    # Get max value
    if sparse.issparse(X):
        max_val = X.data.max() if X.data.size > 0 else 0
        # Sample for decimal check
        sample_data = X.data[:min(10000, X.data.size)]
    else:
        max_val = X.max()
        sample_data = X.flatten()[:min(10000, X.size)]
    
    # Check 1: Max value should be > 10 (typical for counts)
    if max_val < 10:
        return False, f"{layer_name} max={max_val:.2f}, likely log-normalized"
    
    # Check 2: Should not have many decimals (counts are integers)
    has_decimals = np.any(sample_data % 1 != 0)
    if has_decimals:
        decimal_frac = np.mean(sample_data % 1 != 0)
        if decimal_frac > 0.1:  # More than 10% have decimals
            return False, f"{layer_name} has {decimal_frac*100:.1f}% non-integer values"
    
    return True, f"{layer_name} appears to be valid counts (max={max_val:.0f})"


# ============================================================================
# GENE FILTERING FUNCTIONS
# ============================================================================

def filter_low_quality_genes(
    adata: sc.AnnData,
    config: Dict[str, bool],
    use_raw: bool = True,
    verbose: bool = True
) -> sc.AnnData:
    """
    Filter out low-quality genes from AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    config : dict
        Gene filtering configuration
    use_raw : bool
        Whether to filter adata.raw (if exists)
    verbose : bool
        Print detailed information
    
    Returns
    -------
    adata : AnnData
        Filtered AnnData object
    """
    if verbose:
        print("\n" + "="*80)
        print("Gene Filtering")
        print("="*80)
        print("⚠️  Note: This removes technical genes but may remove biologically")
        print("    relevant programs (proliferation, translation, stress)")
    
    # Get gene list
    if use_raw and adata.raw is not None:
        all_genes = adata.raw.var_names.tolist()
        filter_raw = True
        if verbose:
            print(f"\nFiltering adata.raw")
    else:
        all_genes = adata.var_names.tolist()
        filter_raw = False
        if verbose:
            print(f"\nFiltering adata.var")
    
    if verbose:
        print(f"Total genes before filtering: {len(all_genes):,}")
    
    genes_to_remove = []
    
    # 1. Mitochondrial genes
    if config['remove_mt']:
        mt_genes = [g for g in all_genes if g.startswith('MT-')]
        genes_to_remove.extend(mt_genes)
        if verbose:
            print(f"\n  1. Mitochondrial genes (MT-): {len(mt_genes)}")
            if len(mt_genes) > 0:
                print(f"     Examples: {', '.join(mt_genes[:5])}")
    
    # 2. Ribosomal genes
    if config['remove_ribo']:
        ribo_pattern = re.compile(r'^(RPS|RPL|MRPS|MRPL)')
        ribo_genes = [g for g in all_genes if ribo_pattern.match(g)]
        genes_to_remove.extend(ribo_genes)
        if verbose:
            print(f"  2. Ribosomal genes (RPS/RPL/MRPS/MRPL): {len(ribo_genes)}")
            if len(ribo_genes) > 0:
                print(f"     Examples: {', '.join(ribo_genes[:5])}")
    
    # 3. Histone genes
    if config['remove_histone']:
        histone_pattern = re.compile(r'^(H1|H2A|H2B|H3|H4|HIST)')
        histone_genes = [g for g in all_genes if histone_pattern.match(g)]
        genes_to_remove.extend(histone_genes)
        if verbose:
            print(f"  3. Histone genes (H1/H2A/H2B/H3/H4/HIST): {len(histone_genes)}")
            if len(histone_genes) > 0:
                print(f"     Examples: {', '.join(histone_genes[:5])}")
    
    # 4. Pseudogenes (ribosomal pseudogenes only for safety)
    if config['remove_pseudogenes']:
        pseudo_pattern = re.compile(r'^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$')
        pseudo_genes = [g for g in all_genes if pseudo_pattern.match(g)]
        genes_to_remove.extend(pseudo_genes)
        if verbose:
            print(f"  4. Ribosomal pseudogenes (RPS*P*/RPL*P*): {len(pseudo_genes)}")
            if len(pseudo_genes) > 0:
                print(f"     Examples: {', '.join(pseudo_genes[:5])}")
    
    # 5. ENSG genes
    if config['remove_ensg']:
        ensg_pattern = re.compile(r'^ENSG[0-9]+')
        ensg_genes = [g for g in all_genes if ensg_pattern.match(g)]
        genes_to_remove.extend(ensg_genes)
        if verbose:
            print(f"  5. ENSG unannotated genes: {len(ensg_genes)}")
            if len(ensg_genes) > 0:
                print(f"     Examples: {', '.join(ensg_genes[:5])}")
    
    # 6. Unannotated transcripts
    if config['remove_unannotated']:
        unannotated_pattern = re.compile(
            r'^(AC|AL|AP|BX|Z)[0-9]+\.|'
            r'^RP[0-9]+-|'
            r'^CTD-|^CTB-|^CTC-|'
            r'^LINC[0-9]+|'
            r'-AS[0-9]+$|'
            r'-OT[0-9]+$|'
            r'^LOC[0-9]+'
        )
        unannotated_genes = [g for g in all_genes if unannotated_pattern.search(g)]
        genes_to_remove.extend(unannotated_genes)
        if verbose:
            print(f"  6. Unannotated transcripts (AC/AL/RP/CTD/LINC/LOC/etc): {len(unannotated_genes)}")
            if len(unannotated_genes) > 0:
                print(f"     Examples: {', '.join(unannotated_genes[:5])}")
    
    # Remove duplicates
    genes_to_remove = list(set(genes_to_remove))
    genes_to_keep = [g for g in all_genes if g not in genes_to_remove]
    
    if verbose:
        print(f"\n{'─'*80}")
        print(f"Total unique genes to remove: {len(genes_to_remove):,}")
        print(f"Genes to keep: {len(genes_to_keep):,}")
        print(f"Percentage retained: {len(genes_to_keep)/len(all_genes)*100:.1f}%")
    
    # Apply filtering
    if filter_raw:
        if verbose:
            print(f"\nFiltering adata.raw...")
        raw_adata = adata.raw.to_adata()
        raw_adata_filtered = raw_adata[:, genes_to_keep].copy()
        adata.raw = raw_adata_filtered
        if verbose:
            print(f"  ✓ adata.raw filtered: {adata.raw.n_vars:,} genes remain")
        del raw_adata, raw_adata_filtered
    else:
        if verbose:
            print(f"\nFiltering adata.var...")
        adata = adata[:, genes_to_keep].copy()
        if verbose:
            print(f"  ✓ adata filtered: {adata.n_vars:,} genes remain")
    
    if verbose:
        print("\n✓ Gene filtering complete")
    
    return adata


# ============================================================================
# cNMF PREPARATION AND EXECUTION
# ============================================================================

def determine_k_range(n_cells: int, config: Dict) -> List[int]:
    """
    Determine appropriate K range based on cell count.
    ⭐ FIXED: Lower K cap for small cell types
    
    Parameters
    ----------
    n_cells : int
        Number of cells
    config : dict
        cNMF configuration
    
    Returns
    -------
    k_range : list
        List of K values to test
    """
    if not config['k_range_auto'] and config['k_range_custom'] is not None:
        return config['k_range_custom']
    
    # Auto-determine based on cell count
    for threshold, k_range in sorted(config['k_rules'].items(), 
                                     key=lambda x: x[1][0], reverse=True):
        min_cells, k_vals = config['k_rules'][threshold]
        if n_cells >= min_cells:
            # ⭐ Additional safety: cap K based on sqrt(n_cells)
            max_k_safe = max(10, int(np.sqrt(n_cells)))
            k_vals_capped = [k for k in k_vals if k <= max_k_safe]
            if len(k_vals_capped) == 0:
                k_vals_capped = [max_k_safe]
            return k_vals_capped
    
    # Fallback for very small datasets
    return [5, 10, 15]


def prepare_counts_for_cnmf_sparse(
    adata: sc.AnnData, 
    output_dir: Path, 
    name: str
) -> Path:
    """
    Prepare counts matrix for cNMF using sparse-aware writing.
    ⭐ FIXED: No dense conversion, stream write chunks
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    output_dir : Path
        Output directory
    name : str
        Dataset name (sanitized)
    
    Returns
    -------
    counts_file : Path
        Path to saved counts file
    """
    print(f"\n  Preparing counts matrix for cNMF (sparse mode)...")
    
    # ⭐ FIXED: Strict counts source priority
    if 'counts' in adata.layers:
        counts = adata.layers['counts']
        gene_names = adata.var_names
        print(f"    Using adata.layers['counts']: {counts.shape}")
    elif adata.raw is not None:
        counts = adata.raw.X
        gene_names = adata.raw.var_names
        print(f"    Using adata.raw.X: {counts.shape}")
    else:
        # Last resort: assume adata.X is counts
        counts = adata.X
        gene_names = adata.var_names
        print(f"    ⚠️  Using adata.X (assuming counts): {counts.shape}")
    
    # ⭐ CRITICAL: Validate counts
    is_valid, msg = validate_counts_matrix(counts, "Selected layer")
    print(f"    {msg}")
    if not is_valid:
        raise ValueError(
            f"Selected matrix does not appear to be raw counts.\n"
            f"  {msg}\n"
            f"  Please ensure adata.layers['counts'] contains raw integer counts."
        )
    
    # Convert to float32 sparse (cNMF compatible, memory efficient)
    if sparse.issparse(counts):
        if counts.dtype != np.float32:
            print(f"    Converting sparse matrix to float32...")
            counts = counts.astype(np.float32)
    else:
        print(f"    Converting dense to sparse float32...")
        counts = sparse.csr_matrix(counts, dtype=np.float32)
    
    counts_file = output_dir / f"{name}_counts.txt"
    print(f"    Saving to: {counts_file}")
    
    # ⭐ FIXED: Stream write without dense conversion
    # cNMF expects genes × cells format (transpose)
    print(f"    Transposing and writing (genes × cells)...")
    counts_T = counts.T.tocsr()  # Transpose to genes × cells
    
    # Write header
    with open(counts_file, 'w') as f:
        # Header: gene names in first column, cell names as columns
        f.write('\t' + '\t'.join(adata.obs_names) + '\n')
        
        # Write genes row by row (efficient for CSR)
        for gene_idx in range(counts_T.shape[0]):
            gene_row = counts_T.getrow(gene_idx).toarray().flatten()
            f.write(gene_names[gene_idx])
            for val in gene_row:
                f.write(f'\t{val:.6g}')  # 6 significant digits
            f.write('\n')
            
            # Progress
            if (gene_idx + 1) % 1000 == 0:
                print(f"      {gene_idx+1}/{counts_T.shape[0]} genes written...")
    
    print(f"    ✓ Counts saved: {counts_T.shape[0]} genes × {counts_T.shape[1]} cells")
    
    # Clean up
    del counts, counts_T
    gc.collect()
    
    return counts_file


def check_cnmf_complete(cnmf_dir: Path, k_range: List[int], density_threshold: float = 0.1) -> bool:
    """
    Check if cNMF has completed for all K values.
    ⭐ FIXED: Robust check for both dt formats
    
    Parameters
    ----------
    cnmf_dir : Path
        cNMF output directory
    k_range : list
        Expected K values
    density_threshold : float
        Density threshold used
    
    Returns
    -------
    is_complete : bool
        Whether all K values have completed
    """
    if not cnmf_dir.exists():
        return False
    
    # Check for both possible dt formats
    dt_formats = [
        f"dt_{str(density_threshold).replace('.', '_')}",
        f"dt_{density_threshold}"
    ]
    
    for k in k_range:
        found = False
        for dt_str in dt_formats:
            spectra_pattern = f"*.spectra.k_{k}.{dt_str}.consensus.txt"
            usage_pattern = f"*.usages.k_{k}.{dt_str}.consensus.txt"
            
            spectra_files = list(cnmf_dir.glob(spectra_pattern))
            usage_files = list(cnmf_dir.glob(usage_pattern))
            
            if len(spectra_files) > 0 and len(usage_files) > 0:
                found = True
                break
        
        if not found:
            return False
    
    return True


def run_cnmf_for_celltype(
    adata: sc.AnnData,
    cell_type: str,
    output_dir: Path,
    config: Dict
) -> bool:
    """
    Run complete cNMF workflow for a single cell type.
    ⭐ FIXED: Single-worker mode, safe naming, robust checks
    
    Parameters
    ----------
    adata : AnnData
        Subset AnnData for this cell type
    cell_type : str
        Cell type name (original)
    config : dict
        cNMF configuration
    output_dir : Path
        Output directory
    
    Returns
    -------
    success : bool
        Whether cNMF completed successfully
    """
    # ⭐ FIXED: Sanitize name consistently
    cell_type_safe = safe_name(cell_type)
    
    print(f"\n" + "="*80)
    print(f"Running cNMF for: {cell_type}")
    print(f"  Safe name: {cell_type_safe}")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print("="*80)
    
    # Create cell type-specific directory
    ct_output_dir = output_dir / cell_type_safe
    ct_output_dir.mkdir(parents=True, exist_ok=True)
    
    # ⭐ FIXED: Use safe name consistently
    cnmf_dir = ct_output_dir / f"{cell_type_safe}_cNMF"
    
    # Determine K range
    k_range = determine_k_range(adata.n_obs, config)
    print(f"\n  K range: {k_range}")
    print(f"  (Selected based on {adata.n_obs:,} cells)")
    
    # Check if already exists and complete
    if SKIP_IF_EXISTS and check_cnmf_complete(cnmf_dir, k_range, config['density_threshold']):
        print(f"\n⏭️  cNMF output complete for all K values, skipping: {cnmf_dir}")
        return True
    
    # Prepare counts file
    try:
        counts_file = prepare_counts_for_cnmf_sparse(adata, ct_output_dir, cell_type_safe)
    except Exception as e:
        print(f"\n❌ Failed to prepare counts: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Initialize cNMF object
    print(f"\n  Initializing cNMF...")
    cnmf_obj = cNMF(
        output_dir=str(cnmf_dir),
        name=cell_type_safe  # ⭐ FIXED: Use safe name
    )
    
    # Run cNMF pipeline
    try:
        # Step 1: Prepare
        print(f"\n  [Step 1/5] Preparing cNMF...")
        
        # ⭐ FIXED: Cap num_hvg to available genes
        num_hvg = min(config['num_hvg'], adata.n_vars - 1)
        if num_hvg < config['num_hvg']:
            print(f"    ⚠️  Reduced num_hvg from {config['num_hvg']} to {num_hvg} (limited by available genes)")
        
        cnmf_obj.prepare(
            counts_fn=str(counts_file),
            components=k_range,
            n_iter=config['n_iter'],
            seed=config['seed'],
            num_highvar_genes=num_hvg,
            genes_file=None,
        )
        print(f"    ✓ Preparation complete")
        
        # Step 2: Factorize
        print(f"\n  [Step 2/5] Running NMF factorization...")
        total_workers = config['total_workers']
        
        if total_workers == 1 or not config['use_multiprocessing']:
            # ⭐ FIXED: Single-worker mode (default, safest)
            print(f"    Using single-worker mode...")
            print(f"    This may take 15-45 minutes depending on dataset size...")
            cnmf_obj.factorize(worker_i=0, total_workers=1)
            print(f"    ✓ Factorization complete")
        else:
            # Multi-worker mode (advanced users only)
            print(f"    Using multi-worker mode ({total_workers} workers)...")
            print(f"    ⚠️  Note: This requires proper multiprocessing setup")
            from multiprocessing import Process
            
            def worker_fn(worker_i: int):
                try:
                    worker_obj = cNMF(output_dir=str(cnmf_dir), name=cell_type_safe)
                    worker_obj.factorize(worker_i=worker_i, total_workers=total_workers)
                except Exception as e:
                    print(f"      Worker {worker_i} failed: {e}")
                    sys.exit(1)
            
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
                    print(f"      ❌ Worker {i} failed")
                else:
                    print(f"      ✓ Worker {i} finished")
            
            if not all_ok:
                raise RuntimeError("Some workers failed during factorization")
        
        # Step 3: Combine
        print(f"\n  [Step 3/5] Combining results...")
        cnmf_obj.combine()
        print(f"    ✓ Results combined")
        
        # Step 4: K selection
        print(f"\n  [Step 4/5] Computing K selection metrics...")
        cnmf_obj.k_selection_plot(
            close_fig=config['close_clustergram_fig']
        )
        print(f"    ✓ K selection plot generated")
        
        # Step 5: Consensus for each K
        print(f"\n  [Step 5/5] Computing consensus for each K...")
        for k in k_range:
            print(f"    Computing consensus for K={k}...")
            try:
                cnmf_obj.consensus(
                    k=k,
                    density_threshold=config['density_threshold'],
                    show_clustering=config['show_clustering'],
                    close_clustergram_fig=config['close_clustergram_fig']
                )
                print(f"      ✓ K={k} complete")
            except Exception as e:
                print(f"      ⚠️  K={k} failed: {e}")
        
        print(f"\n✓ cNMF pipeline complete for {cell_type}")
        return True
        
    except Exception as e:
        print(f"\n❌ cNMF failed for {cell_type}: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Clean up counts file to save space
        if counts_file.exists():
            counts_file.unlink()
            print(f"  ✓ Cleaned up counts file")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Main pipeline execution."""
    
    # ========================================================================
    # STEP 1: Load annotated data
    # ========================================================================
    
    print("\n" + "="*80)
    print("STEP 1: Loading Annotated Data")
    print("="*80)
    
    print(f"\nLoading: {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    
    print(f"\n📊 Data Summary:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    
    # Check cell type column
    if CELLTYPE_KEY not in adata.obs.columns:
        raise ValueError(f"Cell type key '{CELLTYPE_KEY}' not found in adata.obs")
    
    # ⭐ NEW: Optional confidence filtering
    if USE_CONFIDENCE_FILTER and CONFIDENCE_KEY in adata.obs.columns:
        print(f"\n⭐ Applying confidence filter (>{MIN_CONFIDENCE})...")
        n_before = adata.n_obs
        adata = adata[adata.obs[CONFIDENCE_KEY] > MIN_CONFIDENCE].copy()
        n_after = adata.n_obs
        print(f"  Removed {n_before - n_after:,} low-confidence cells ({100*(n_before-n_after)/n_before:.1f}%)")
        print(f"  Remaining: {n_after:,} cells")
    
    print(f"\n📋 Cell Type Distribution:")
    celltype_counts = adata.obs[CELLTYPE_KEY].value_counts()
    for ct, count in celltype_counts.items():
        pct = 100 * count / adata.n_obs
        print(f"  {ct:30s}: {count:7,} ({pct:5.2f}%)")
    
    # ========================================================================
    # STEP 2: Determine which cell types to process
    # ========================================================================
    
    print("\n" + "="*80)
    print("STEP 2: Selecting Cell Types to Process")
    print("="*80)
    
    if PROCESS_CELLTYPES == 'auto':
        celltypes_to_process = celltype_counts[celltype_counts >= MIN_CELLS_PER_TYPE].index.tolist()
        print(f"\n⭐ Auto mode: Processing cell types with ≥{MIN_CELLS_PER_TYPE} cells")
    else:
        celltypes_to_process = PROCESS_CELLTYPES
        print(f"\n⭐ Manual mode: Processing specified cell types")
    
    print(f"\n📋 Cell types to process ({len(celltypes_to_process)}):")
    for ct in celltypes_to_process:
        count = celltype_counts.get(ct, 0)
        if count < MIN_CELLS_PER_TYPE:
            print(f"  ⚠️  {ct:30s}: {count:7,} (below threshold, will skip)")
        else:
            print(f"  ✓  {ct:30s}: {count:7,}")
    
    celltypes_to_process = [ct for ct in celltypes_to_process 
                           if celltype_counts.get(ct, 0) >= MIN_CELLS_PER_TYPE]
    
    if len(celltypes_to_process) == 0:
        print(f"\n❌ No cell types meet the minimum cell count threshold ({MIN_CELLS_PER_TYPE})")
        return
    
    print(f"\n✓ Will process {len(celltypes_to_process)} cell types")
    
    # ========================================================================
    # STEP 3: Process each cell type
    # ========================================================================
    
    print("\n" + "="*80)
    print("STEP 3: Processing Each Cell Type")
    print("="*80)
    
    results_summary = []
    
    for idx, cell_type in enumerate(celltypes_to_process, 1):
        print(f"\n{'='*80}")
        print(f"Processing {idx}/{len(celltypes_to_process)}: {cell_type}")
        print(f"{'='*80}")
        
        # ⭐ FIXED: Use safe name for file operations
        cell_type_safe = safe_name(cell_type)
        subset_h5ad_file = OUTPUT_DIR / f"{cell_type_safe}_filtered.h5ad"
        
        # Subset to this cell type
        print(f"\nSubsetting to {cell_type}...")
        adata_subset = adata[adata.obs[CELLTYPE_KEY] == cell_type].copy()
        print(f"  ✓ Subset: {adata_subset.n_obs:,} cells × {adata_subset.n_vars:,} genes")
        
        # Check if filtered version already exists
        if subset_h5ad_file.exists() and not OVERWRITE_GENE_FILTER:
            print(f"\n⏭️  Filtered h5ad exists, loading: {subset_h5ad_file}")
            adata_filtered = sc.read_h5ad(subset_h5ad_file)
        else:
            # Apply gene filtering
            adata_filtered = filter_low_quality_genes(
                adata_subset,
                GENE_FILTER_CONFIG,
                use_raw=True,
                verbose=True
            )
            
            # Save filtered subset
            if SAVE_SUBSET_H5AD:
                print(f"\nSaving filtered subset: {subset_h5ad_file}")
                adata_filtered.write_h5ad(subset_h5ad_file, compression='gzip')
                print(f"  ✓ Saved")
        
        print(f"\n📊 Filtered Data Summary:")
        print(f"  Cells: {adata_filtered.n_obs:,}")
        if adata_filtered.raw is not None:
            print(f"  Genes (raw): {adata_filtered.raw.n_vars:,}")
        print(f"  Genes (HVG): {adata_filtered.n_vars:,}")
        
        # Run cNMF
        success = run_cnmf_for_celltype(
            adata_filtered,
            cell_type,  # Original name for display
            OUTPUT_DIR,
            CNMF_CONFIG
        )
        
        results_summary.append({
            'cell_type': cell_type,
            'cell_type_safe': cell_type_safe,
            'n_cells': adata_filtered.n_obs,
            'n_genes': adata_filtered.raw.n_vars if adata_filtered.raw is not None else adata_filtered.n_vars,
            'success': success
        })
        
        # Clean up
        del adata_subset, adata_filtered
        gc.collect()
        
        print(f"\n✓ Completed processing: {cell_type}")
    
    # ========================================================================
    # STEP 4: Generate Summary Report
    # ========================================================================
    
    print("\n" + "="*80)
    print("STEP 4: Generating Summary Report")
    print("="*80)
    
    # Save summary
    summary_df = pd.DataFrame(results_summary)
    summary_file = OUTPUT_DIR / "processing_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"\n✓ Summary saved: {summary_file}")
    
    # Print summary
    print(f"\n📊 Processing Summary:")
    print(f"  Total cell types processed: {len(results_summary)}")
    print(f"  Successful: {sum(r['success'] for r in results_summary)}")
    print(f"  Failed: {sum(not r['success'] for r in results_summary)}")
    
    print(f"\n📋 Details:")
    for r in results_summary:
        status = "✓" if r['success'] else "❌"
        print(f"  {status} {r['cell_type']:30s}: {r['n_cells']:7,} cells, {r['n_genes']:7,} genes")
    
    # ========================================================================
    # COMPLETE
    # ========================================================================
    
    print("\n" + "="*80)
    print("✓ PIPELINE COMPLETE")
    print("="*80)
    print(f"Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"\nNext steps:")
    print(f"  1. Review K-selection plots in each cell type subdirectory")
    print(f"  2. Run cnmf_results_analysis_v1_1.py to compare K values")
    print(f"  3. Choose optimal K for each cell type")
    print(f"  4. Extract and analyze GEP signatures")
    print("="*80)


if __name__ == '__main__':
    main()
