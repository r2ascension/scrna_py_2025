#!/usr/bin/env python3
"""
scANVI BBKNN Analysis - Standalone Script v1.6.1
=================================================

Performs BBKNN clustering on already-processed scANVI results.
This script is separated from the main pipeline to avoid environment conflicts.

Usage:
------
# Single cell type
python scanvi_bbknn_analysis_v1_6_1.py --cell_type T_cells

# Multiple cell types
python scanvi_bbknn_analysis_v1_6_1.py --cell_types T_cells B_cells Myeloid

# Custom input directory
python scanvi_bbknn_analysis_v1_6_1.py \
    --input_dir /path/to/downstream_analysis \
    --cell_types T_cells

Requirements:
-------------
- bbknn
- scanpy
- numpy, pandas
- matplotlib

Author: r2end
Date: 2024-12-17
Version: v1.6.1
"""

import argparse
import sys
from pathlib import Path
import warnings
import logging
from typing import Optional, Dict, List, Tuple, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import gc
import json

# Check bbknn availability
try:
    import bbknn
    BBKNN_VERSION = bbknn.__version__ if hasattr(bbknn, '__version__') else 'unknown'
    BBKNN_AVAILABLE = True
except ImportError:
    print("❌ ERROR: bbknn not installed")
    print("Install: pip install bbknn")
    sys.exit(1)

# ============================================================================
# PACKAGE VERSION CHECKING
# ============================================================================

def get_package_version(package_name: str) -> Optional[str]:
    """Get package version string"""
    try:
        if package_name == 'scanpy':
            return sc.__version__
        elif package_name == 'pandas':
            return pd.__version__
        elif package_name == 'numpy':
            return np.__version__
        elif package_name == 'bbknn' and BBKNN_AVAILABLE:
            return BBKNN_VERSION
        else:
            import importlib
            mod = importlib.import_module(package_name)
            return getattr(mod, '__version__', 'unknown')
    except Exception:
        return None


def check_package_versions() -> Dict[str, Any]:
    """
    Check versions of critical packages
    
    Returns:
    --------
    dict with version info
    """
    versions = {}
    
    versions['scanpy'] = get_package_version('scanpy')
    versions['pandas'] = get_package_version('pandas')
    versions['numpy'] = get_package_version('numpy')
    versions['bbknn'] = get_package_version('bbknn')
    
    return {'versions': versions}


# Global package info
PACKAGE_INFO = check_package_versions()

# ⭐ Selective warning filtering (keep important warnings visible)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)
# Keep UserWarning and RuntimeWarning visible

# ============================================================================
# LOGGING
# ============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Log package versions at startup
logger.info("="*70)
logger.info("Package Versions")
logger.info("="*70)
for pkg, ver in PACKAGE_INFO['versions'].items():
    if ver:
        logger.info(f"  {pkg}: {ver}")
    else:
        logger.warning(f"  {pkg}: version unknown")
logger.info("="*70)

# ============================================================================
# CONFIGURATION
# ============================================================================

# ⭐ Cell type to filename mapping (based on actual data structure)
CELLTYPE_FILE_MAPPING = {
    'tcell': 'adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad',
    'bcell': 'adata_bcell_scvi_celltypist_scanvi_final_v2.h5ad',
    'myeloid': 'adata_myeloid_scvi_celltypist_scanvi_final_v2.h5ad',
    'stromal_vascular': 'adata_stromal_vascular_FINAL.h5ad',
    # Also support uppercase versions
    'T_cells': 'adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad',
    'B_cells': 'adata_bcell_scvi_celltypist_scanvi_final_v2.h5ad',
    'Myeloid': 'adata_myeloid_scvi_celltypist_scanvi_final_v2.h5ad',
    'Stromal_Vascular': 'adata_stromal_vascular_FINAL.h5ad',
}

# Default data directory
DEFAULT_DATA_DIR = Path('/home/h2048/data/core_data')
DEFAULT_OUTPUT_DIR = Path('/home/h2048/data/py/1217/downstream_analysis_v1_6_1')

BBKNN_PARAMS = {
    'batch_key': 'dataset',  # Will auto-detect if not found
    'n_pcs': 50,
    'neighbors_within_batch': 3,
    'metric': 'euclidean',
    'trim': None,
    'neighbors_key': 'neighbors_bbknn',
    'leiden_resolutions': [0.5, 1.0, 1.5],
    'default_leiden_res': 1.0,
    'umap_min_dist': 0.3,
    'umap_spread': 1.0,
    # Safety parameters
    'min_cells_per_batch': 10,
    'min_batches': 2,
    'auto_adjust_neighbors': True,
}

# ============================================================================
# COLUMN INFERENCE
# ============================================================================

def infer_columns(adata: sc.AnnData) -> Dict[str, Optional[str]]:
    """
    Infer column names from adata
    Replicates logic from main pipeline
    """
    cols = {
        'celltype_col': None,
        'batch_col': None,
        'disease_col': None,
    }
    
    # Celltype
    for col in ['scanvi_predictions', 'cell_type_scanvi_filt', 'cell_type', 'celltype']:
        if col in adata.obs.columns:
            cols['celltype_col'] = col
            logger.info(f"  ✓ Celltype column: {col}")
            break
    
    # Batch
    for col in ['dataset', 'batch', 'sample_id', 'orig.ident']:
        if col in adata.obs.columns:
            cols['batch_col'] = col
            logger.info(f"  ✓ Batch column: {col}")
            break
    
    # Disease (optional)
    for col in ['disease_status', 'Disease', 'condition', 'group']:
        if col in adata.obs.columns:
            cols['disease_col'] = col
            logger.info(f"  ✓ Disease column: {col}")
            break
    
    return cols


# ============================================================================
# BBKNN ANALYSIS
# ============================================================================

def run_bbknn_standalone(
    adata: sc.AnnData,
    cols: Dict[str, Optional[str]],
    output_dir: Path,
    params: Dict
) -> bool:
    """
    Standalone BBKNN analysis
    
    Parameters:
    -----------
    adata : AnnData
        Input data (from main pipeline final.h5ad)
    cols : dict
        Column mapping
    output_dir : Path
        Output directory for BBKNN results
    params : dict
        BBKNN parameters
    
    Returns:
    --------
    success : bool
    """
    logger.info("\n" + "="*70)
    logger.info("BBKNN CLUSTERING ANALYSIS")
    logger.info("="*70)
    
    bbknn_dir = output_dir / 'bbknn'
    bbknn_dir.mkdir(exist_ok=True, parents=True)
    
    # ===== STEP 1: Find batch key =====
    batch_key = None
    if params['batch_key'] in adata.obs.columns:
        batch_key = params['batch_key']
    else:
        batch_key = cols.get('batch_col')
    
    if batch_key is None or batch_key not in adata.obs.columns:
        logger.error("❌ No batch key found")
        return False
    
    logger.info(f"\nBatch key: {batch_key}")
    
    # ===== STEP 2: Batch-disease audit =====
    if cols['disease_col'] and cols['disease_col'] in adata.obs.columns:
        logger.warning("\n" + "="*70)
        logger.warning("⚠️  BBKNN BATCH-DISEASE AUDIT")
        logger.warning("="*70)
        
        batch_disease = pd.crosstab(
            adata.obs[batch_key],
            adata.obs[cols['disease_col']],
            margins=True
        )
        
        logger.warning("\nBatch × Disease contingency:")
        logger.warning(f"\n{batch_disease.to_string()}")
        
        # Check confounding
        confounded_batches = []
        for batch in batch_disease.index[:-1]:
            disease_counts = batch_disease.loc[batch, :]
            disease_counts = disease_counts[disease_counts.index != 'All']
            if (disease_counts > 0).sum() == 1:
                confounded_batches.append(batch)
        
        if confounded_batches:
            logger.warning(f"\n⚠️  WARNING: {len(confounded_batches)} batches are disease-specific!")
            logger.warning(f"   Examples: {confounded_batches[:5]}")
            logger.warning("   BBKNN may remove real disease biology!")
        
        logger.warning("="*70)
    
    # ===== STEP 3: Analyze batch composition =====
    batch_counts = adata.obs[batch_key].value_counts()
    n_batches = len(batch_counts)
    n_valid_batches = n_batches
    
    logger.info(f"\nBatch composition:")
    logger.info(f"  Total batches: {n_batches}")
    logger.info(f"  Cells per batch:")
    for batch, count in batch_counts.head(10).items():
        logger.info(f"    {batch}: {count:,} cells")
    if n_batches > 10:
        logger.info(f"    ... and {n_batches - 10} more batches")
    
    # Safety check
    if n_batches < params['min_batches']:
        logger.error(f"❌ Insufficient batches ({n_batches} < {params['min_batches']})")
        return False
    
    # ===== STEP 4: Filter small batches =====
    small_batches = batch_counts[batch_counts < params['min_cells_per_batch']]
    
    if len(small_batches) > 0:
        logger.warning(f"\n⚠️  Found {len(small_batches)} small batches (< {params['min_cells_per_batch']} cells):")
        for batch, count in small_batches.items():
            logger.warning(f"    {batch}: {count} cells")
        
        n_valid_batches = n_batches - len(small_batches)
        
        if n_valid_batches < params['min_batches']:
            logger.error(f"❌ Too few valid batches after filtering ({n_valid_batches} < {params['min_batches']})")
            return False
        
        logger.info(f"\n  Filtering to {n_valid_batches} valid batches...")
        valid_batches = batch_counts[batch_counts >= params['min_cells_per_batch']].index
        adata_filtered = adata[adata.obs[batch_key].isin(valid_batches)].copy()
        logger.info(f"  Retained: {adata_filtered.n_obs:,} cells ({adata_filtered.n_obs / adata.n_obs * 100:.1f}%)")
    else:
        adata_filtered = adata.copy()
        logger.info(f"  ✓ All batches have sufficient cells")
    
    # ===== STEP 5: Adjust neighbors =====
    min_batch_size = adata_filtered.obs[batch_key].value_counts().min()
    neighbors_within = params['neighbors_within_batch']
    
    if params['auto_adjust_neighbors'] and neighbors_within >= min_batch_size:
        original_neighbors = neighbors_within
        neighbors_within = max(1, min_batch_size - 1)
        logger.warning(f"\n⚠️  Auto-adjusting neighbors_within_batch:")
        logger.warning(f"    {original_neighbors} → {neighbors_within} (min batch size: {min_batch_size})")
    
    # ===== STEP 6: Backup UMAP =====
    # ⭐ FIX: Preserve original UMAP (from scANVI pipeline)
    if 'X_umap' in adata.obsm and 'X_umap_scanvi' not in adata.obsm:
        adata.obsm['X_umap_scanvi'] = np.asarray(adata.obsm['X_umap']).copy()
        logger.info("\n  ✓ Backed up X_umap → X_umap_scanvi")
    elif 'X_umap_scanvi' in adata.obsm:
        logger.info("\n  ✓ X_umap_scanvi already exists (preserved)")
    
    # ===== STEP 7: PCA if needed =====
    need_pca = ('X_pca' not in adata_filtered.obsm) or (adata_filtered.obsm['X_pca'].shape[1] < params['n_pcs'])
    if need_pca:
        logger.info("  Computing PCA...")
        use_hvg = ('highly_variable' in adata_filtered.var.columns) and np.any(adata_filtered.var['highly_variable'].values)
        try:
            # ⭐ FALLBACK: Handle different scanpy PCA API versions
            try:
                sc.pp.pca(adata_filtered, n_comps=params['n_pcs'], svd_solver='arpack', use_highly_variable=use_hvg)
            except TypeError:
                # Older scanpy versions may not support use_highly_variable
                if use_hvg:
                    adata_hvg = adata_filtered[:, adata_filtered.var['highly_variable']].copy()
                    sc.pp.pca(adata_hvg, n_comps=params['n_pcs'], svd_solver='arpack')
                    adata_filtered.obsm['X_pca'] = adata_hvg.obsm['X_pca']
                else:
                    sc.pp.pca(adata_filtered, n_comps=params['n_pcs'], svd_solver='arpack')
            logger.info(f"  ✓ PCA computed: {adata_filtered.obsm['X_pca'].shape}")
        except Exception as e:
            logger.error(f"  ❌ PCA failed: {e}")
            return False
    else:
        logger.info(f"  ✓ Using existing PCA: {adata_filtered.obsm['X_pca'].shape}")
    
    # ===== STEP 8: Run BBKNN =====
    logger.info(f"\n  Running BBKNN (neighbors_within_batch={neighbors_within})...")
    
    try:
        bbknn.bbknn(
            adata_filtered,
            batch_key=batch_key,
            neighbors_within_batch=neighbors_within,
            n_pcs=params['n_pcs'],
            metric=params['metric'],
            trim=params['trim'],
            key_added=params['neighbors_key'],
            copy=False,
        )
        logger.info("  ✓ BBKNN graph constructed")
        
    except ValueError as e:
        logger.error(f"  ❌ BBKNN failed: {e}")
        return False
    except Exception as e:
        logger.error(f"  ❌ Unexpected error: {e}")
        return False
    
    # ===== STEP 9: UMAP =====
    logger.info("  Computing UMAP...")
    try:
        sc.tl.umap(
            adata_filtered,
            neighbors_key=params['neighbors_key'],
            min_dist=params['umap_min_dist'],
            spread=params['umap_spread'],
            random_state=42,
        )
        logger.info("  ✓ UMAP computed")
    except Exception as e:
        logger.error(f"  ❌ UMAP failed: {e}")
        return False
    
    # ===== STEP 10: Transfer to original adata =====
    # ⭐ FIX: Save BBKNN UMAP to separate key, preserve original
    cell_mask = adata.obs_names.isin(adata_filtered.obs_names)
    n_transferred = cell_mask.sum()
    
    logger.info(f"\n  Transferring results to original adata:")
    logger.info(f"    Cells: {n_transferred:,} / {adata.n_obs:,} ({n_transferred/adata.n_obs*100:.1f}%)")
    
    # Create BBKNN UMAP array
    umap_bbknn = np.full((adata.n_obs, 2), np.nan)
    umap_bbknn[cell_mask] = adata_filtered.obsm['X_umap']
    
    # ⭐ CRITICAL: Save as X_umap_bbknn, keep X_umap intact
    adata.obsm['X_umap_bbknn'] = umap_bbknn
    
    # ===== STEP 11: Leiden clustering =====
    logger.info("\n  Leiden clustering...")
    for res in params['leiden_resolutions']:
        key = f'leiden_bbknn_res{res}'
        try:
            sc.tl.leiden(
                adata_filtered,
                neighbors_key=params['neighbors_key'],
                resolution=res,
                key_added=key,
                flavor='igraph',
                n_iterations=2,
                directed=False
            )
        except TypeError:
            sc.tl.leiden(
                adata_filtered,
                neighbors_key=params['neighbors_key'],
                resolution=res,
                key_added=key,
                n_iterations=2
            )
        
        # Transfer to full adata (with fallback for filtered cells)
        # ⭐ FIXED: Handle Categorical type mismatch
        try:
            # Get filtered values as strings first (avoids categorical mismatch)
            filtered_values = adata_filtered.obs[key].astype(str).values
            
            # Create full Series with 'filtered_out' as default
            leiden_full = pd.Series('filtered_out', index=adata.obs_names, dtype='object')
            
            # Assign filtered values using .loc
            leiden_full.loc[cell_mask] = filtered_values
            
            # Convert to categorical if original was categorical (preserves type)
            if isinstance(adata_filtered.obs[key].dtype, pd.CategoricalDtype):
                # Get all unique values (including 'filtered_out')
                all_categories = list(leiden_full.unique())
                # Ensure 'filtered_out' is first (for clarity)
                if 'filtered_out' in all_categories:
                    all_categories.remove('filtered_out')
                    all_categories = ['filtered_out'] + sorted(all_categories)
                else:
                    all_categories = sorted(all_categories)
                
                leiden_full = pd.Categorical(leiden_full, categories=all_categories)
            
            adata.obs[key] = leiden_full
            
        except Exception as e:
            logger.warning(f"    ⚠️  Failed to transfer {key}: {e}")
            logger.warning("    Using string type as fallback...")
            # Fallback: use string type
            leiden_full = pd.Series('filtered_out', index=adata.obs_names, dtype='object')
            leiden_full.loc[cell_mask] = adata_filtered.obs[key].astype(str).values
            adata.obs[key] = leiden_full
        
        logger.info(f"    Res {res}: {adata_filtered.obs[key].nunique()} clusters")
    
    default_key = f'leiden_bbknn_res{params["default_leiden_res"]}'
    if default_key in adata.obs.columns:
        adata.obs['leiden_bbknn'] = adata.obs[default_key]
    
    # ===== STEP 12: Save filtered data =====
    filtered_path = bbknn_dir / 'adata_bbknn_filtered.h5ad'
    adata_filtered.write_h5ad(filtered_path, compression='gzip')
    logger.info(f"\n  ✓ Saved filtered data: {filtered_path}")
    
    # ===== STEP 13: Comparison plots =====
    logger.info("\n  Creating comparison plots...")
    
    try:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        adata_plot = adata_filtered
        
        # Original UMAP
        if 'X_umap_scanvi' in adata.obsm and cols['celltype_col']:
            adata_orig_subset = adata[cell_mask].copy()
            
            sc.pl.umap(adata_orig_subset, color=cols['celltype_col'], ax=axes[0, 0], show=False, 
                       title='scANVI UMAP', frameon=False, size=50)
            sc.pl.umap(adata_orig_subset, color=batch_key, ax=axes[0, 1], show=False,
                       title='scANVI UMAP (Batch)', frameon=False, size=50)
            
            del adata_orig_subset
        else:
            axes[0, 0].axis('off')
            axes[0, 1].axis('off')
        
        # BBKNN results
        if cols['celltype_col']:
            sc.pl.umap(adata_plot, color=cols['celltype_col'], ax=axes[0, 2], show=False,
                       title='BBKNN UMAP', frameon=False, size=50)
        else:
            axes[0, 2].axis('off')
        
        sc.pl.umap(adata_plot, color=batch_key, ax=axes[1, 0], show=False,
                   title='BBKNN UMAP (Batch)', frameon=False, size=50)
        
        for i, res in enumerate([0.5, 1.0]):
            if i < 2:
                key = f'leiden_bbknn_res{res}'
                if key in adata_plot.obs.columns:
                    sc.pl.umap(adata_plot, color=key, ax=axes[1, i+1], show=False,
                              title=f'BBKNN Leiden (res={res})', frameon=False, size=50)
        
        plt.tight_layout()
        plt.savefig(bbknn_dir / 'bbknn_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("  ✓ Saved: bbknn_comparison.png")
        
    except Exception as e:
        logger.warning(f"  ⚠️  Plot generation failed: {e}")
    
    # ===== STEP 14: Save summary =====
    summary = {
        'total_cells': int(adata.n_obs),
        'cells_used': int(n_transferred),
        'cells_filtered': int(adata.n_obs - n_transferred),
        'total_batches': int(n_batches),
        'batches_used': int(n_valid_batches),
        'batches_filtered': int(len(small_batches)),
        'neighbors_within_batch_used': int(neighbors_within),
        'neighbors_within_batch_requested': int(params['neighbors_within_batch']),
        'min_batch_size': int(min_batch_size),
    }
    
    with open(bbknn_dir / 'bbknn_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info("\n" + "="*70)
    logger.info("✓ BBKNN COMPLETE")
    logger.info("="*70)
    logger.info(f"Cells used: {n_transferred:,} / {adata.n_obs:,} ({n_transferred/adata.n_obs*100:.1f}%)")
    logger.info(f"Batches used: {n_valid_batches} / {n_batches}")
    if len(small_batches) > 0:
        logger.info(f"Filtered batches: {list(small_batches.index)}")
    logger.info(f"\nResults:")
    logger.info(f"  Original UMAP: X_umap (preserved)")
    logger.info(f"  scANVI backup: X_umap_scanvi (if exists)")
    logger.info(f"  BBKNN UMAP: X_umap_bbknn (new)")
    logger.info(f"  Filtered data: {filtered_path}")
    logger.info("="*70)
    
    # Cleanup
    del adata_filtered
    gc.collect()
    
    return True


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def process_celltype(
    input_dir: Path,
    celltype: str,
    params: Dict,
    data_dir: Optional[Path] = None,
    output_dir: Optional[Path] = None
) -> bool:
    """
    Process one cell type
    
    Parameters:
    -----------
    input_dir : Path
        Base input directory (e.g., downstream_analysis_v1_6_1) - for processed files
    celltype : str
        Cell type name (e.g., tcell, T_cells)
    params : dict
        BBKNN parameters
    data_dir : Path, optional
        Directory containing raw data files (default: DEFAULT_DATA_DIR)
    output_dir : Path, optional
        Output directory (default: input_dir)
    
    Returns:
    --------
    success : bool
    """
    logger.info("\n" + "="*70)
    logger.info(f"PROCESSING: {celltype}")
    logger.info("="*70)
    
    # Determine data source
    data_dir = data_dir or DEFAULT_DATA_DIR
    output_dir = output_dir or input_dir
    
    # Try to find the input file
    final_h5ad = None
    
    # Strategy 1: Try processed file from downstream analysis
    celltype_dir = input_dir / celltype
    processed_file = celltype_dir / f'{celltype}_processed_final.h5ad'
    if processed_file.exists():
        final_h5ad = processed_file
        logger.info(f"  Using processed file from downstream analysis")
    
    # Strategy 2: Try raw file from core_data
    if final_h5ad is None:
        filename = CELLTYPE_FILE_MAPPING.get(celltype)
        if filename:
            raw_file = data_dir / filename
            if raw_file.exists():
                final_h5ad = raw_file
                logger.info(f"  Using raw file from core_data")
            else:
                logger.warning(f"  Expected file not found: {raw_file}")
        else:
            logger.warning(f"  No filename mapping for celltype: {celltype}")
    
    # Strategy 3: Try alternative naming (case-insensitive)
    if final_h5ad is None:
        # Try with different case variations (skip if already checked in Strategy 2)
        celltype_lower = celltype.lower()
        for alt_name, alt_file in CELLTYPE_FILE_MAPPING.items():
            if alt_name.lower() == celltype_lower and alt_name != celltype:
                alt_path = data_dir / alt_file
                if alt_path.exists():
                    final_h5ad = alt_path
                    logger.info(f"  Using file with alternative naming: {alt_file}")
                    break
    
    if final_h5ad is None or not final_h5ad.exists():
        logger.error(f"❌ File not found for celltype: {celltype}")
        logger.error(f"  Tried:")
        logger.error(f"    - Processed: {processed_file}")
        if celltype in CELLTYPE_FILE_MAPPING:
            logger.error(f"    - Raw: {data_dir / CELLTYPE_FILE_MAPPING[celltype]}")
        return False
    
    logger.info(f"\nLoading: {final_h5ad}")
    
    try:
        adata = sc.read_h5ad(final_h5ad)
        logger.info(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars} genes")
        
        # Infer columns
        logger.info("\nInferring column names...")
        cols = infer_columns(adata)
        
        if not cols['celltype_col']:
            logger.warning("⚠️  No celltype column found (optional for BBKNN)")
        
        if not cols['batch_col']:
            logger.error("❌ No batch column found (required for BBKNN)")
            return False
        
        # Run BBKNN (use output_dir for results)
        celltype_output_dir = output_dir / celltype
        celltype_output_dir.mkdir(exist_ok=True, parents=True)
        success = run_bbknn_standalone(adata, cols, celltype_output_dir, params)
        
        if success:
            # Save updated adata with BBKNN results
            output_path = celltype_output_dir / f'{celltype}_with_bbknn.h5ad'
            adata.write_h5ad(output_path, compression='gzip')
            logger.info(f"\n  ✓ Saved: {output_path}")
        
        return success
        
    except Exception as e:
        logger.error(f"\n❌ ERROR: {e}", exc_info=True)
        return False
    finally:
        if 'adata' in locals():
            del adata
        gc.collect()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='BBKNN analysis for scANVI downstream results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single cell type
  python scanvi_bbknn_analysis_v1_6_1.py --cell_type T_cells
  
  # Multiple cell types
  python scanvi_bbknn_analysis_v1_6_1.py --cell_types T_cells B_cells Myeloid
  
  # Custom directory and parameters
  python scanvi_bbknn_analysis_v1_6_1.py \\
      --input_dir /path/to/downstream_analysis \\
      --cell_types T_cells \\
      --neighbors 5 \\
      --min_cells 20
        """
    )
    
    parser.add_argument(
        '--input_dir',
        type=str,  # ⭐ FIX: Use str instead of Path for better compatibility
        default=str(DEFAULT_OUTPUT_DIR),
        help='Base input directory for processed files (default: downstream_analysis_v1_6_1)'
    )
    
    parser.add_argument(
        '--data_dir',
        type=str,  # ⭐ FIX: Use str instead of Path for better compatibility
        default=str(DEFAULT_DATA_DIR),
        help='Directory containing raw data files (default: /home/h2048/data/core_data)'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,  # ⭐ FIX: Use str instead of Path for better compatibility
        default=None,
        help='Output directory (default: same as input_dir)'
    )
    
    parser.add_argument(
        '--cell_type',
        type=str,
        help='Single cell type to process'
    )
    
    parser.add_argument(
        '--cell_types',
        nargs='+',
        help='Multiple cell types to process'
    )
    
    parser.add_argument(
        '--neighbors',
        type=int,
        default=3,
        help='neighbors_within_batch (default: 3)'
    )
    
    parser.add_argument(
        '--min_cells',
        type=int,
        default=10,
        help='min_cells_per_batch (default: 10)'
    )
    
    parser.add_argument(
        '--min_batches',
        type=int,
        default=2,
        help='min_batches (default: 2)'
    )
    
    args = parser.parse_args()
    
    # Determine cell types to process
    if args.cell_type:
        cell_types = [args.cell_type]
    elif args.cell_types:
        cell_types = args.cell_types
    else:
        # Default: all available
        cell_types = ['tcell', 'bcell', 'myeloid', 'stromal_vascular']
    
    # Update parameters
    params = BBKNN_PARAMS.copy()
    params['neighbors_within_batch'] = args.neighbors
    params['min_cells_per_batch'] = args.min_cells
    params['min_batches'] = args.min_batches
    
    # ⭐ FIX: Convert string arguments to Path objects
    args.input_dir = Path(args.input_dir)
    args.data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir) if args.output_dir else args.input_dir
    
    # Log configuration
    logger.info("="*70)
    logger.info("scANVI BBKNN ANALYSIS - v1.6.1")
    logger.info("="*70)
    logger.info(f"\nInput directory (processed): {args.input_dir}")
    logger.info(f"Data directory (raw): {args.data_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cell types: {cell_types}")
    logger.info(f"\nBBKNN parameters:")
    logger.info(f"  neighbors_within_batch: {params['neighbors_within_batch']}")
    logger.info(f"  min_cells_per_batch: {params['min_cells_per_batch']}")
    logger.info(f"  min_batches: {params['min_batches']}")
    logger.info(f"\nPackage versions:")
    for pkg, ver in PACKAGE_INFO['versions'].items():
        if ver:
            logger.info(f"  {pkg}: {ver}")
    
    # Process each cell type
    results = {}
    for ct in cell_types:
        success = process_celltype(
            args.input_dir, 
            ct, 
            params,
            data_dir=args.data_dir,
            output_dir=output_dir
        )
        results[ct] = success
    
    # Summary
    logger.info("\n" + "="*70)
    logger.info("BBKNN ANALYSIS COMPLETE")
    logger.info("="*70)
    
    for ct, success in results.items():
        status = "✓ SUCCESS" if success else "❌ FAILED"
        logger.info(f"  {ct}: {status}")
    
    successful = sum(results.values())
    logger.info(f"\nTotal: {successful}/{len(results)} successful")
    logger.info("="*70)


if __name__ == "__main__":
    main()
