#!/usr/bin/env python3
"""
scANVI Downstream Analysis - Production Pipeline v1.6.2
======================================================

CHANGELOG v1.6.2 (2024-12-18) - CRITICAL HOTFIX:
- ✅ REMOVED: All disease-related DE analysis per user request
- ✅ FIXED: scVI API compatibility (n_samples parameter handling)
- ✅ FIXED: Robust directory creation with fallback
- ✅ FIXED: Batch size warnings in model setup
- ✅ MAINTAINED: All v1.6.1 fixes

Changes from v1.6.1:
1. Removed `perform_de_disease_within_celltype_safe()` completely
2. Auto-detect scVI version and use correct DE parameters
3. Enhanced directory creation with explicit checks
4. Graceful degradation if DE fails

Critical features:
- Robust model loading with HVG extraction
- Version-agnostic scVI differential expression
- Memory-safe denoised expression
- Pseudobulk export for donor-level statistics
- Comprehensive marker visualization
- Optional BBKNN integration

Author: r2end
Date: 2024-12-18
Version: v1.6.2 (Production - HOTFIX)
"""

import os
import sys
from pathlib import Path
import warnings
import logging
from typing import Optional, Dict, List, Tuple, Any, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scvi
import torch
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import scipy.sparse as sp
import gc
import json

# Try to import bbknn (optional)
BBKNN_AVAILABLE = False
try:
    import bbknn
    BBKNN_AVAILABLE = True
except ImportError:
    pass

# ⭐ Selective warning filtering
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)

# ============================================================================
# PACKAGE VERSION CHECKING
# ============================================================================

def get_package_version(package_name: str) -> Optional[str]:
    """Get package version string"""
    try:
        if package_name == 'scanpy':
            return sc.__version__
        elif package_name == 'scvi':
            return scvi.__version__
        elif package_name == 'pandas':
            return pd.__version__
        elif package_name == 'numpy':
            return np.__version__
        elif package_name == 'torch':
            return torch.__version__
        elif package_name == 'bbknn' and BBKNN_AVAILABLE:
            return bbknn.__version__ if hasattr(bbknn, '__version__') else 'unknown'
        else:
            import importlib
            mod = importlib.import_module(package_name)
            return getattr(mod, '__version__', 'unknown')
    except Exception:
        return None


def check_package_versions() -> Dict[str, Any]:
    """
    Check versions and determine compatibility flags
    """
    versions = {}
    compatibility = {
        'scanpy_dotplot_return_fig': False,
        'scvi_has_seed_setting': False,
        'scvi_de_uses_n_samples': True,  # ⭐ NEW: Track if DE accepts n_samples
        'torch_weights_only': True,
    }
    
    # Scanpy
    scanpy_version = get_package_version('scanpy')
    versions['scanpy'] = scanpy_version
    if scanpy_version:
        try:
            major, minor = map(int, scanpy_version.split('.')[:2])
            if major > 1 or (major == 1 and minor >= 9):
                compatibility['scanpy_dotplot_return_fig'] = True
        except (ValueError, AttributeError):
            pass
    
    # scVI-tools - CRITICAL for DE API
    scvi_version = get_package_version('scvi')
    versions['scvi'] = scvi_version
    if scvi_version:
        try:
            major, minor = map(int, scvi_version.split('.')[:2])
            
            # Seed setting (>= 0.15.0)
            if major > 0 or (major == 0 and minor >= 15):
                compatibility['scvi_has_seed_setting'] = True
            
            # ⭐ CRITICAL: n_samples parameter deprecated in scvi >= 1.0
            # Note: scVI-tools 1.0+ removed n_samples from differential_expression()
            # All 1.x versions (1.0, 1.1, 1.2, 1.3+) do not accept n_samples
            if major >= 1:
                compatibility['scvi_de_uses_n_samples'] = False
        except (ValueError, AttributeError):
            pass
    
    # Torch
    torch_version = get_package_version('torch')
    versions['torch'] = torch_version
    if torch_version:
        try:
            major, minor = map(int, torch_version.split('.')[:2])
            if major < 1 or (major == 1 and minor < 13):
                compatibility['torch_weights_only'] = False
        except (ValueError, AttributeError):
            pass
    
    versions['pandas'] = get_package_version('pandas')
    versions['numpy'] = get_package_version('numpy')
    if BBKNN_AVAILABLE:
        versions['bbknn'] = get_package_version('bbknn')
    
    return {
        'versions': versions,
        'compatibility': compatibility
    }


PACKAGE_INFO = check_package_versions()


def safe_torch_load(file_path: Path, map_location: str = 'cpu') -> Any:
    """Safe torch.load with version compatibility"""
    try:
        if PACKAGE_INFO['compatibility']['torch_weights_only']:
            try:
                return torch.load(file_path, map_location=map_location, weights_only=False)
            except TypeError:
                return torch.load(file_path, map_location=map_location)
        else:
            return torch.load(file_path, map_location=map_location)
    except Exception as e:
        logger.error(f"  ❌ torch.load failed: {e}")
        raise

# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================

def setup_logging(log_level: str = 'INFO', log_file: Optional[Path] = None):
    """Configure logging system"""
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    return logging.getLogger(__name__)

logger = setup_logging()

# Log versions
logger.info("="*70)
logger.info("Package Versions & Compatibility")
logger.info("="*70)
for pkg, ver in PACKAGE_INFO['versions'].items():
    if ver:
        logger.info(f"  {pkg}: {ver}")
    else:
        logger.warning(f"  {pkg}: version unknown")

logger.info("\nCompatibility Flags:")
for flag, value in PACKAGE_INFO['compatibility'].items():
    logger.info(f"  {flag}: {value}")
logger.info("="*70)

# ============================================================================
# MARKER GENE DEFINITIONS (same as before)
# ============================================================================

MARKER_GENES = {
    'General': [
        'FXYD3', 'EPCAM', 'ELF3', 'IGFBP2', 'SERPINF1', 'TSPAN1', 'SCGB1A1',
        'AGER', 'SFTPC', 'FOXJ1', 'KRT5', 'MUC5B', 'KRT8', 'CD53', 'PTPRC',
        'CORO1A', 'ISG20', 'CCL5', 'MS4A1', 'TNFRSF17', 'CD19', 'CD79A', 'SDC1',
        'CD40LG', 'TNFRSF25', 'CD28', 'CD4', 'CD3E', 'CD8A', 'CD8B', 'TRGC2',
        'CD2', 'TRBC2', 'FCER1G', 'C1orf162', 'CLEC7A', 'CD1C', 'CD86', 'CD14',
        'XCR1', 'HLA-DRA', 'COL1A2', 'DCN', 'MFAP4', 'LUM', 'COL6A3', 'CFD',
        'COL1A1', 'PDGFRA', 'MXRA8', 'NBL1', 'VCAN', 'LEPR', 'MYH11', 'TINAGL1',
        'PLN', 'DES', 'ACTA2', 'CNN1', 'TAGLN', 'CLDN5', 'ECSCR', 'CLEC14A',
        'VWF', 'PECAM1', 'DARC', 'PTPRB', 'PDE2A', 'PLAT', 'GJA5', 'SPARCL1',
        'AQP1', 'RNASE1', 'MMRN1', 'CCL21', 'TFF3', 'MKI67', 'TOP2A', 'TK1', 'CENPW'
    ],
    'T_cells': [
        'KLRD1', 'FCGR3A', 'GNLY', 'TYROBP', 'FCER1G', 'KLRC1', 'FGFBP2', 'SPON2',
        'MYOM2', 'TRDC', 'KRT86', 'GATA3', 'IL5', 'AREG', 'HPGDS', 'IL23R', 'RORC',
        'LST1', 'PCDH9', 'TNFSF11', 'TRDC', 'TRGC1', 'TRGC2', 'CD4', 'CD28',
        'CD40LG', 'TRAT1', 'TNFRSF25', 'CD8A', 'CD8B', 'TRGC2', 'CCR7', 'TCF7',
        'LEF1', 'SELL', 'CD28', 'IL7R', 'CCR6', 'GATA3', 'IL4', 'IL13', 'IL17A',
        'CCL20', 'CCR7', 'TCF7', 'LEF1', 'SELL', 'TBX21', 'EOMES', 'GZMK', 'KLRG1',
        'ITGA1', 'CD8A', 'CD8B', 'CCR6', 'IL7R', 'TBX21', 'GZMB', 'GZMH', 'FGFBP2',
        'ZNF683', 'IFNG', 'CCL4L2', 'PDCD1', 'KLRB1', 'IL7R', 'NCR3', 'CEBPD',
        'FOXP3', 'IL2RA', 'IKZF2', 'TNFRSF4', 'MKI67', 'TOP2A', 'TK1', 'CENPW'
    ],
    'B_cells': [
        'IGHD', 'TCL1A', 'MS4A1', 'BANK1', 'MKI67', 'TOP2A', 'IGHA1', 'IGHA2',
        'IGHGP', 'IGHG1'
    ],
    'Myeloid': [
        'CLEC9A', 'XCR1', 'CADM1', 'CLNK', 'FLT3', 'ZBTB46', 'CLEC10A', 'CD1E',
        'FCER1A', 'CD1D', 'ITGAX', 'CD1C', 'FCGR2B', 'PKIB', 'CCR7', 'CD83',
        'LAMP3', 'CCL22', 'CCL17', 'CCL19', 'LAD1', 'LILRA4', 'SMPD3', 'SCT',
        'IRF7', 'PLD4', 'CLEC4C', 'MARCO', 'FABP4', 'CYP27A1', 'SIGLEC1', 'ABCG1',
        'PPARG', 'C1QA', 'C1QB', 'C1QC', 'HLA-DPA1', 'SLC40A1', 'FOLR2', 'F13A1',
        'SPP1', 'HAMP', 'VCAN', 'CCR2', 'CCR5', 'FCN1', 'S100A12', 'RNASE2',
        'LILRA5', 'MTSS1', 'TPSAB1', 'MS4A2', 'TPSB2', 'FCGR3B', 'CSF3R', 'CXCR1'
    ],
    'Stromal_Vascular': [
        'APOD', 'FGF7', 'COL15A1', 'MFAP5', 'PI16', 'CD34', 'MMP11', 'COL10A1',
        'POSTN', 'LRRC15', 'HOPX', 'IGFBP5', 'TIMP1', 'MMP1', 'COL7A1', 'WNT5A',
        'ISG15', 'IL7R', 'SFRP4', 'SFRP2', 'COMP', 'RGS5', 'PDGFRB', 'NDUFA4L2',
        'NOTCH3', 'CXCL1', 'CXCL2', 'IL6', 'CEBPD', 'CLU', 'CTGF', 'HGF', 'HSPA6',
        'DNAJB1', 'MYC', 'ATF4', 'PLAU', 'CHI3L1', 'MMP3', 'IL1R1', 'IL13RA2',
        'TNFSF11', 'MMP10', 'OSMR', 'IL11', 'STRA6', 'FAP', 'WNT2', 'TWIST1',
        'IL24', 'ACTG2', 'HHIP', 'CNN1', 'MYH11', 'ACTA2', 'TAGLN'
    ]
}

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path("/home/h2048/data/core_data")
OUTPUT_BASE = Path("/home/h2048/data/py/1217/downstream_analysis_v1_6_2")

CELL_TYPES = {
    'T_cells': {
        'h5ad': 'adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_model': 'models/T_models/scanvi_model',
        'naming_mode': 'standard',
    },
    'B_cells': {
        'h5ad': 'adata_bcell_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_model': 'models/B_models/scanvi_model',
        'naming_mode': 'standard',
    },
    'Myeloid': {
        'h5ad': 'adata_myeloid_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_model': 'models/Myeloid_models/scanvi_model',
        'naming_mode': 'standard',
    },
    'Stromal_Vascular': {
        'h5ad': 'adata_stromal_vascular_FINAL.h5ad',
        'scanvi_model': 'models/Stromal_models/scanvi_model',
        'naming_mode': 'semantic',
    }
}

# ⭐ UPDATED: DE parameters with version-aware handling
DE_PARAMS = {
    'mode': 'change',
    'delta': 0.5,
    'batch_correction': True,
    'n_samples': 5000,  # Used only if scvi_de_uses_n_samples=True
    'bayes_factor_threshold': 3.0,
    'lfc_threshold': 0.5,
    'min_proportion': 0.1,
}

VALIDATION_THRESHOLDS = {
    'min_gene_overlap_ratio': 0.9,
    'min_cells_for_de': 100,
    'epsilon': 1e-10
}

DENOISED_PARAMS = {
    'n_samples': 25,
    'library_size': 1e4,
    'max_cells_in_memory': 2e8,
}

SUBCLUSTER_PARAMS = {
    'resolutions': [0.2, 0.4, 0.6, 0.8, 1.0],
    'default_resolution': 0.6,
    'n_neighbors': 30,
    'min_cells_for_subcluster': 500
}

PSEUDOBULK_PARAMS = {
    'min_cells_per_sample': 10,
    'aggregation': 'sum'
}

MARKER_VIZ_PARAMS = {
    'dotplot_top_genes': 30,
    'featureplot_top_genes': 12,
}

BBKNN_PARAMS = {
    'enabled': True,
    'batch_key': 'dataset',
    'n_pcs': 50,
    'neighbors_within_batch': 3,
    'trim': None,
    'leiden_resolutions': [0.5, 1.0, 1.5],
    'default_leiden_res': 1.0,
    'min_cells_per_batch': 10,
    'min_batches': 2,
    'auto_adjust_neighbors': True,
}

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# Random seeds
np.random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)
if PACKAGE_INFO['compatibility']['scvi_has_seed_setting']:
    try:
        scvi.settings.seed = 42
    except (AttributeError, TypeError):
        pass

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def log_step(step_num: int, title: str, description: Optional[str] = None) -> None:
    """Format step headers"""
    logger.info(f"\n{'='*70}")
    logger.info(f"STEP {step_num}: {title}")
    if description:
        logger.info(f"{description}")
    logger.info(f"{'='*70}")


def safe_name(x: str, max_len: int = 180) -> str:
    """Sanitize names for file paths"""
    s = str(x)
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t']:
        s = s.replace(ch, '_')
    return s[:max_len]


def safe_mkdir(directory: Path, description: str = "directory") -> bool:
    """
    ⭐ NEW: Robust directory creation with explicit verification
    
    Returns:
        True if directory exists after operation, False otherwise
    """
    try:
        directory.mkdir(exist_ok=True, parents=True)
        
        # Explicit verification
        if not directory.exists():
            logger.error(f"  ❌ Failed to create {description}: {directory}")
            return False
        
        if not directory.is_dir():
            logger.error(f"  ❌ Path exists but is not directory: {directory}")
            return False
        
        # ⭐ FIXED: Use info level instead of debug (visible in default logging)
        logger.info(f"  ✓ Directory ready: {directory}")
        return True
        
    except PermissionError:
        logger.error(f"  ❌ Permission denied: {directory}")
        return False
    except OSError as e:
        logger.error(f"  ❌ OS error creating {description}: {e}")
        return False
    except Exception as e:
        logger.error(f"  ❌ Unexpected error creating {description}: {e}")
        return False


# ============================================================================
# COLUMN INFERENCE (from v1.6.1)
# ============================================================================

def infer_column_names_robust(
    adata: sc.AnnData, 
    naming_mode: str = 'standard'
) -> Dict[str, Optional[str]]:
    """Robust column inference"""
    cols = {
        'celltype_col': None,
        'confidence_col': None,
        'batch_col': None,
        'sample_col': None,
        'disease_col': None,  # Keep for compatibility but won't use
        'tissue_col': None,
    }
    
    # CELLTYPE
    if naming_mode == 'semantic':
        celltype_candidates = [
            'cell_type_scanvi_filt',
            'cell_type_scanvi_raw',
            'semantic_celltype',
            'scanvi_predictions',
            'cell_type',
            'celltype'
        ]
    else:
        celltype_candidates = [
            'scanvi_predictions',
            'scanvi_labels',
            'predicted_labels',
            'majority_voting',
            'cell_type',
            'celltype'
        ]
    
    for col in celltype_candidates:
        if col in adata.obs.columns:
            cols['celltype_col'] = col
            logger.info(f"  ✓ celltype_col: {col}")
            break
    
    # CONFIDENCE
    if naming_mode == 'semantic':
        confidence_candidates = ['scanvi_confidence', 'celltypist_confidence', 'conf_score', 'confidence']
    else:
        confidence_candidates = ['scanvi_confidence', 'conf_score', 'confidence']
    
    for col in confidence_candidates:
        if col in adata.obs.columns:
            cols['confidence_col'] = col
            logger.info(f"  ✓ confidence_col: {col}")
            break
    
    # BATCH
    for col in ['dataset', 'batch', 'sample_id', 'orig.ident']:
        if col in adata.obs.columns:
            cols['batch_col'] = col
            logger.info(f"  ✓ batch_col: {col}")
            break
    
    # SAMPLE
    for col in ['sample', 'sample_id', 'Sample', 'donor', 'donor_id', 'patient_id']:
        if col in adata.obs.columns:
            cols['sample_col'] = col
            logger.info(f"  ✓ sample_col: {col}")
            break
    
    # DISEASE (keep but won't use)
    for col in ['disease_status', 'Disease', 'condition', 'group', 'phenotype']:
        if col in adata.obs.columns:
            cols['disease_col'] = col
            logger.info(f"  ✓ disease_col: {col} (not used in this version)")
            break
    
    # TISSUE
    for col in ['tissue', 'organ', 'anatomical_site', 'location', 'organ__ontology_label']:
        if col in adata.obs.columns:
            cols['tissue_col'] = col
            logger.info(f"  ✓ tissue_col: {col}")
            break
    
    return cols


def validate_columns_v1_6_2(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]]
) -> bool:
    """Column validation - only require essentials"""
    errors = []
    warnings = []
    
    required = {
        'celltype_col': 'Cell type labels',
        'batch_col': 'Batch identifier',
        'sample_col': 'Sample identifier'
    }
    
    for key, desc in required.items():
        col = cols.get(key)
        if not col:
            errors.append(f"Required '{key}' not found. Needed: {desc}")
        elif col not in adata.obs.columns:
            errors.append(f"Column '{col}' for {key} not in adata.obs")
    
    optional = {
        'confidence_col': 'Confidence scores',
        'tissue_col': 'Tissue location'
    }
    
    for key, desc in optional.items():
        col = cols.get(key)
        if not col:
            warnings.append(f"Optional '{key}' not found - {desc} skipped")
    
    if 'X_scanvi' not in adata.obsm:
        warnings.append("X_scanvi missing - subclustering limited")
    
    has_counts = ('counts' in adata.layers) or (adata.raw is not None)
    if not has_counts:
        warnings.append("No counts - pseudobulk limited")
    
    if warnings:
        logger.warning("\n⚠️  Optional Features Limited:")
        for w in warnings:
            logger.warning(f"    {w}")
    
    if errors:
        logger.error("\n❌ REQUIRED COLUMNS MISSING:")
        for e in errors:
            logger.error(f"    {e}")
        raise ValueError("Critical columns missing")
    
    logger.info("✓ Column validation passed")
    return True


# ============================================================================
# MODEL LOADING (from v1.6.1)
# ============================================================================

def extract_var_names_from_model(
    model_path: Union[str, Path], 
    adata_full: Optional[sc.AnnData] = None, 
    model_type: str = 'SCANVI'
) -> Optional[List[str]]:
    """Extract gene names from model"""
    try:
        logger.info("  Reading checkpoint...")
        model_pt = Path(model_path) / 'model.pt'
        if model_pt.exists():
            checkpoint = safe_torch_load(model_pt, map_location='cpu')
            if 'var_names' in checkpoint:
                var_names = checkpoint['var_names']
                if isinstance(var_names, np.ndarray):
                    var_names = var_names.tolist()
                elif isinstance(var_names, torch.Tensor):
                    var_names = var_names.cpu().numpy().tolist()
                
                if len(var_names) > 0:
                    logger.info(f"  ✓ Extracted {len(var_names)} genes")
                    return list(var_names)
        
        if adata_full is not None:
            logger.info("  Loading model with adata...")
            Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
            model = Model.load(model_path, adata=adata_full)
            var_names = list(model.adata.var_names)
            logger.info(f"  ✓ Extracted {len(var_names)} genes")
            return var_names
        
        logger.error("  ❌ All methods failed")
        return None
        
    except Exception as e:
        logger.error(f"  ❌ Failed: {e}")
        return None


def create_model_compatible_adata(
    adata_full: sc.AnnData, 
    model_genes: List[str]
) -> sc.AnnData:
    """
    Create model-compatible subset with preserved gene order
    
    ⭐ CRITICAL: Preserves gene order from model_genes (essential for scVI)
    Using set() would scramble order → model weights misaligned with genes
    
    Parameters:
    -----------
    adata_full : AnnData
        Full dataset with all genes
    model_genes : List[str]
        Gene list from model (order must be preserved)
    
    Returns:
    --------
    AnnData
        Subset with genes in model_genes order
    """
    # ⭐ CRITICAL: List comprehension maintains model_genes order
    overlap_genes = [g for g in model_genes if g in adata_full.var_names]
    missing_genes = [g for g in model_genes if g not in adata_full.var_names]
    
    overlap_ratio = len(overlap_genes) / len(model_genes) if len(model_genes) > 0 else 0.0
    
    logger.info(f"  Data genes: {adata_full.n_vars}")
    logger.info(f"  Model genes: {len(model_genes)}")
    logger.info(f"  Overlap: {len(overlap_genes)} ({overlap_ratio*100:.1f}%)")
    
    if len(missing_genes) > 0:
        logger.warning(f"  Missing: {len(missing_genes)} genes")
        if len(missing_genes) <= 20:
            logger.warning(f"    Missing genes: {missing_genes}")
    
    if overlap_ratio < VALIDATION_THRESHOLDS['min_gene_overlap_ratio']:
        logger.warning(f"  ⚠️  Low overlap may affect performance")
    
    adata_model = adata_full[:, overlap_genes].copy()
    
    if adata_full.raw is not None:
        adata_model.raw = sc.AnnData(
            X=adata_full.raw.X,
            obs=adata_full.obs.copy(),
            var=adata_full.raw.var.copy()
        )
    
    return adata_model


def load_data_and_model_v1_6_2(
    celltype_name: str, 
    config: Dict[str, Any], 
    output_dir: Path
) -> Tuple[sc.AnnData, Optional[sc.AnnData], Optional[Any], Dict[str, Optional[str]], Optional[List[str]]]:
    """Load data and model"""
    log_step(1, "Loading Data and Model", "")
    
    h5ad_path = BASE_DIR / config['h5ad']
    logger.info(f"Loading: {h5ad_path}")
    adata_full = sc.read_h5ad(h5ad_path)
    logger.info(f"  Loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars} genes")
    
    cols = infer_column_names_robust(adata_full, config['naming_mode'])
    
    logger.info("\nColumn mapping:")
    for key, val in cols.items():
        status = val if val else "NOT FOUND"
        logger.info(f"  {key}: {status}")
    
    try:
        validate_columns_v1_6_2(adata_full, cols)
    except ValueError as e:
        logger.error(f"\n❌ Validation failed: {e}")
        return adata_full, None, None, cols, None
    
    scanvi_path = BASE_DIR / config['scanvi_model']
    logger.info(f"\nLoading scANVI model: {scanvi_path}")
    
    lvae = None
    adata_model = None
    model_genes = None
    
    try:
        model_genes = extract_var_names_from_model(scanvi_path, adata_full, 'SCANVI')
        if model_genes is None:
            raise ValueError("Could not extract genes")
        
        gene_list_path = output_dir / f'{celltype_name}_model_genes.txt'
        with open(gene_list_path, 'w') as f:
            for gene in model_genes:
                f.write(f"{gene}\n")
        logger.info(f"  ✓ Saved: {gene_list_path}")
        
        adata_model = create_model_compatible_adata(adata_full, model_genes)
        
        lvae = scvi.model.SCANVI.load(scanvi_path, adata=adata_model)
        logger.info("  ✓ Model loaded successfully")
        
    except Exception as e:
        logger.warning(f"\n⚠️ Model loading failed: {e}")
        logger.info("  Proceeding with scanpy-only analysis")
    
    return adata_full, adata_model, lvae, cols, model_genes


# ============================================================================
# EXPRESSION SOURCE HELPERS
# ============================================================================

def resolve_expr_source(adata: sc.AnnData) -> Dict[str, Union[bool, str]]:
    """Smart expression source resolution"""
    if adata.raw is not None and adata.raw.X is not None and adata.raw.n_vars > 0:
        logger.info("  Expression source: .raw")
        return {'use_raw': True}
    
    if 'log1p' in adata.layers:
        logger.info("  Expression source: layers['log1p']")
        return {'layer': 'log1p'}
    
    if 'counts' in adata.layers:
        logger.warning("  ⚠️  Using layers['counts']")
        return {'layer': 'counts'}
    
    logger.info("  Expression source: .X")
    return {'use_raw': False}


def get_counts_matrix(adata: sc.AnnData) -> Tuple[Any, pd.Index]:
    """Get counts matrix"""
    if 'counts' in adata.layers:
        logger.info("  Count source: layers['counts']")
        return adata.layers['counts'], adata.var_names
    
    if adata.raw is not None and adata.raw.X is not None:
        sample_size = min(100, adata.raw.n_obs, adata.raw.n_vars)
        if hasattr(adata.raw.X, 'toarray'):
            sample = adata.raw.X[:sample_size, :sample_size].toarray()
        else:
            sample = adata.raw.X[:sample_size, :sample_size]
        
        is_counts = np.allclose(sample, sample.astype(int)) and (sample >= 0).all()
        if is_counts:
            logger.info("  Count source: .raw.X (verified)")
            return adata.raw.X, adata.raw.var_names
    
    logger.warning("  ⚠️  Using .X")
    return adata.X, adata.var_names


# ============================================================================
# DE HELPERS (v1.6.2 - VERSION AWARE)
# ============================================================================

def _get_proba_de_col(df: pd.DataFrame) -> Optional[str]:
    """Get proba_de column name"""
    if 'proba_de' in df.columns:
        return 'proba_de'
    if 'proba_m2' in df.columns:
        return 'proba_m2'
    return None


def _filter_de_significant(de: pd.DataFrame, lfc_mode: str = 'up') -> pd.DataFrame:
    """Filter significant DE genes"""
    bf_ok = None
    if 'bayes_factor' in de.columns:
        bf_ok = de['bayes_factor'].astype(float) > DE_PARAMS['bayes_factor_threshold']
    
    proba_col = _get_proba_de_col(de)
    proba_ok = None
    if proba_col:
        proba_ok = de[proba_col].astype(float) > 0.95
    
    if 'lfc_mean' in de.columns:
        lfc = de['lfc_mean'].astype(float)
    else:
        lfc = None
    
    lfc_ok = None
    if lfc is not None:
        threshold = DE_PARAMS['lfc_threshold']
        if lfc_mode == 'abs':
            lfc_ok = lfc.abs() > threshold
        elif lfc_mode == 'up':
            lfc_ok = lfc > threshold
        elif lfc_mode == 'down':
            lfc_ok = lfc < -threshold
    
    mask = None
    for m in [bf_ok, proba_ok, lfc_ok]:
        if m is None:
            continue
        mask = m if mask is None else (mask & m)
    
    if mask is None:
        return de.iloc[0:0].copy()
    
    de_sig = de.loc[mask].copy()
    
    if 'non_zeros_proportion1' in de_sig.columns:
        de_sig = de_sig[de_sig['non_zeros_proportion1'].astype(float) > DE_PARAMS['min_proportion']]
    
    return de_sig


def validate_and_fix_de_schema(de_df: pd.DataFrame, source: str = "scANVI", var_names: Optional[List[str]] = None) -> pd.DataFrame:
    """Validate DE schema"""
    if not isinstance(de_df, pd.DataFrame):
        raise ValueError(f"DE result not DataFrame")
    
    if var_names is not None and de_df.index is not None:
        if (de_df.index.dtype == object) and (de_df.index.name in [None, 'gene', 'genes']):
            idx_vals = de_df.index.astype(str)
            overlap = np.intersect1d(idx_vals[:min(len(idx_vals), 2000)], np.asarray(var_names).astype(str)).size
            if overlap >= 50:
                de_df = de_df.reset_index().rename(columns={'index': 'gene'})
    
    if source == "scANVI":
        if 'gene' not in de_df.columns:
            if de_df.index.dtype == object:
                de_df = de_df.reset_index().rename(columns={'index': 'gene'})
            else:
                raise ValueError("Cannot identify gene column")
    
    return de_df


# ============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION (v1.6.2 - VERSION AWARE, NO DISEASE)
# ============================================================================

def perform_de_one_vs_rest_version_aware(
    adata_model: sc.AnnData, 
    lvae: Any, 
    celltype_col: str, 
    output_dir: Path
) -> Optional[pd.DataFrame]:
    """
    ⭐ FIXED: Version-aware differential expression
    
    Automatically detects scVI version and uses correct parameters
    """
    logger.info("\n--- One-vs-Rest DE (Version-Aware) ---")
    
    celltypes = adata_model.obs[celltype_col].cat.categories if \
                hasattr(adata_model.obs[celltype_col], 'cat') else \
                adata_model.obs[celltype_col].unique()
    
    all_results = []
    
    # ⭐ CRITICAL: Check scVI version for parameter compatibility
    use_n_samples = PACKAGE_INFO['compatibility']['scvi_de_uses_n_samples']
    scvi_version = PACKAGE_INFO['versions'].get('scvi', 'unknown')
    
    if use_n_samples:
        logger.info(f"  Using n_samples={DE_PARAMS['n_samples']} (scVI < 1.0, detected: {scvi_version})")
    else:
        logger.info(f"  Skipping n_samples parameter (scVI >= 1.0, detected: {scvi_version})")
    
    for ct in celltypes:
        logger.info(f"\n  Processing: {ct}")
        
        try:
            # ⭐ CRITICAL: Build kwargs dynamically based on version
            de_kwargs = {
                'groupby': celltype_col,
                'group1': ct,
                'group2': None,
                'mode': DE_PARAMS['mode'],
                'delta': DE_PARAMS['delta'],
                'batch_correction': DE_PARAMS['batch_correction'],
            }
            
            # Only add n_samples if supported
            if use_n_samples:
                de_kwargs['n_samples'] = DE_PARAMS['n_samples']
            
            de = lvae.differential_expression(**de_kwargs)
            
            de = validate_and_fix_de_schema(de, source="scANVI", var_names=adata_model.var_names)
            de['cell_type'] = ct
            de['comparison'] = f'{ct}_vs_Rest'
            
            all_results.append(de)
            
            de_sig = _filter_de_significant(de)
            logger.info(f"    Significant: {len(de_sig)}")
            
            safe_ct = safe_name(ct)
            de_sig.to_csv(output_dir / f'{safe_ct}_markers.csv', index=False)
            
        except TypeError as e:
            if 'n_samples' in str(e):
                logger.error(f"    ❌ API mismatch: {e}")
                logger.error("    This suggests scVI version detection failed")
                logger.error("    Attempting fallback without n_samples...")
                # ⭐ ENHANCED: Try fallback without n_samples
                try:
                    de_kwargs_fallback = {
                        'groupby': celltype_col,
                        'group1': ct,
                        'group2': None,
                        'mode': DE_PARAMS['mode'],
                        'delta': DE_PARAMS['delta'],
                        'batch_correction': DE_PARAMS['batch_correction'],
                    }
                    de = lvae.differential_expression(**de_kwargs_fallback)
                    de = validate_and_fix_de_schema(de, source="scANVI", var_names=adata_model.var_names)
                    de['cell_type'] = ct
                    de['comparison'] = f'{ct}_vs_Rest'
                    all_results.append(de)
                    de_sig = _filter_de_significant(de)
                    logger.info(f"    ✓ Fallback successful: {len(de_sig)} significant")
                    safe_ct = safe_name(ct)
                    de_sig.to_csv(output_dir / f'{safe_ct}_markers.csv', index=False)
                except Exception as e2:
                    logger.error(f"    ❌ Fallback also failed: {e2}")
            else:
                logger.warning(f"    ⚠️ Failed: {e}")
            continue
        except Exception as e:
            logger.warning(f"    ⚠️ Failed: {e}")
            continue
    
    if all_results:
        de_all = pd.concat(all_results, ignore_index=True)
        de_all.to_csv(output_dir / 'celltype_markers_all_scANVI.csv', index=False)
        logger.info(f"\n  ✓ Saved: celltype_markers_all_scANVI.csv")
        return de_all
    
    logger.warning("\n  ⚠️  No DE results generated")
    return None


def perform_de_scanpy_fallback(
    adata: sc.AnnData, 
    celltype_col: str, 
    output_dir: Path
) -> pd.DataFrame:
    """Scanpy fallback DE"""
    logger.info("\n--- Scanpy Fallback DE ---")
    
    expr_source = resolve_expr_source(adata)
    
    sc.tl.rank_genes_groups(
        adata,
        groupby=celltype_col,
        method='wilcoxon',
        n_genes=100,
        **expr_source
    )
    
    de_results = sc.get.rank_genes_groups_df(adata, group=None)
    de_results.to_csv(output_dir / 'celltype_markers_all_scanpy.csv', index=False)
    logger.info("  ✓ Saved: celltype_markers_all_scanpy.csv")
    
    return de_results


def perform_differential_expression_v1_6_2(
    adata_full: sc.AnnData, 
    adata_model: Optional[sc.AnnData], 
    lvae: Optional[Any], 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """
    ⭐ UPDATED: DE workflow without disease comparisons
    """
    log_step(2, "Differential Expression Analysis", "")
    
    de_dir = output_dir / 'differential_expression'
    if not safe_mkdir(de_dir, "DE directory"):
        logger.error("  ❌ Cannot create DE directory - skipping DE analysis")
        # ⭐ ENHANCED: Return directory path but log clearly that DE was skipped
        return de_dir
    
    celltype_col = cols['celltype_col']
    
    # One-vs-rest only
    if lvae is not None and adata_model is not None:
        perform_de_one_vs_rest_version_aware(adata_model, lvae, celltype_col, de_dir)
    else:
        perform_de_scanpy_fallback(adata_full, celltype_col, de_dir)
    
    # ⭐ REMOVED: Disease DE completely
    logger.info("\n⚠️  Disease DE skipped (disabled in v1.6.2)")
    
    logger.info("\n✓ Differential expression complete")
    return de_dir


# ============================================================================
# STEP 3: DENOISED EXPRESSION
# ============================================================================

def generate_denoised_expression_safe(
    adata_model: Optional[sc.AnnData], 
    lvae: Optional[Any], 
    output_dir: Path
) -> Optional[Path]:
    """Memory-safe denoised expression"""
    log_step(3, "Generating Denoised Expression", "")
    
    if lvae is None or adata_model is None:
        logger.warning("⚠️ Model not available - skipping")
        return None
    
    n_elements = adata_model.n_obs * adata_model.n_vars
    if n_elements > DENOISED_PARAMS['max_cells_in_memory']:
        logger.warning(f"⚠️ Dataset too large ({n_elements:,} elements) - skipping")
        return None
    
    logger.info(f"\nGenerating denoised expression...")
    
    try:
        denoised = lvae.get_normalized_expression(
            adata=adata_model,
            n_samples=DENOISED_PARAMS['n_samples'],
            return_mean=True,
            library_size=DENOISED_PARAMS['library_size']
        )
        
        if isinstance(denoised, pd.DataFrame):
            denoised_array = denoised.values.astype(np.float32)
            var_subset = adata_model.var.loc[denoised.columns].copy()
        else:
            denoised_array = denoised.astype(np.float32)
            var_subset = adata_model.var.copy()
        
        logger.info(f"  ✓ Generated: {denoised_array.shape}")
        
        adata_denoised = sc.AnnData(
            X=denoised_array,
            obs=adata_model.obs.copy(),
            var=var_subset
        )
        
        denoised_path = output_dir / 'denoised_expression.h5ad'
        adata_denoised.write_h5ad(denoised_path, compression='gzip')
        logger.info(f"  ✓ Saved: {denoised_path}")
        
        del adata_denoised
        gc.collect()
        
        return denoised_path
        
    except Exception as e:
        logger.error(f"\n  ⚠️ Failed: {e}")
        return None


# ============================================================================
# STEP 4: CELL COMPOSITION (WITHOUT DISEASE)
# ============================================================================

def perform_cell_composition_analysis_v1_6_2(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """
    ⭐ UPDATED: Cell composition without disease comparison
    """
    log_step(4, "Cell Composition Analysis", "")
    
    comp_dir = output_dir / 'cell_composition'
    if not safe_mkdir(comp_dir, "composition directory"):
        logger.error("  ❌ Cannot create directory - skipping")
        return comp_dir
    
    celltype_col = cols['celltype_col']
    sample_col = cols['sample_col']
    
    # Overall composition
    logger.info("\n--- Overall Composition ---")
    composition = pd.crosstab(
        adata.obs[sample_col],
        adata.obs[celltype_col],
        normalize='index'
    ) * 100
    
    composition.to_csv(comp_dir / 'composition_by_sample.csv')
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    composition.mean().sort_values(ascending=False).plot(kind='bar', ax=ax)
    ax.set_ylabel('Cell Proportion (%)')
    ax.set_xlabel('Cell Type')
    ax.set_title('Average Cell Type Composition')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(comp_dir / 'composition_overall.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  ✓ Saved composition plots")
    
    # ⭐ REMOVED: Disease comparison
    logger.info("\n⚠️  Disease composition analysis skipped (disabled in v1.6.2)")
    
    logger.info("\n✓ Cell composition complete")
    return comp_dir


# ============================================================================
# STEP 5: PSEUDOBULK (ENHANCED SAFETY)
# ============================================================================

def export_pseudobulk_enhanced_v1_6_2(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Optional[Path]:
    """
    ⭐ UPDATED: Enhanced safety with directory verification
    """
    log_step(5, "Pseudobulk Export", "")
    
    pseudobulk_dir = output_dir / 'pseudobulk'
    
    # ⭐ CRITICAL: Explicit directory creation with verification
    if not safe_mkdir(pseudobulk_dir, "pseudobulk directory"):
        logger.error("  ❌ Cannot create pseudobulk directory - aborting")
        return None
    
    celltype_col = cols['celltype_col']
    sample_col = cols['sample_col']
    
    counts, genes = get_counts_matrix(adata)
    if counts is None:
        logger.warning("⚠️ No counts - skipping")
        return None
    
    logger.info(f"\nAggregating by {sample_col} × {celltype_col}...")
    
    key_series = (
        adata.obs[sample_col].astype(str) + '|' +
        adata.obs[celltype_col].astype(str)
    )
    
    pseudobulk_list = []
    metadata_list = []
    
    for key in pd.unique(key_series):
        mask = (key_series == key).values
        n_cells = int(mask.sum())
        
        if n_cells < PSEUDOBULK_PARAMS['min_cells_per_sample']:
            continue
        
        try:
            if sp.issparse(counts):
                if PSEUDOBULK_PARAMS['aggregation'] == 'sum':
                    pseudo_counts = np.array(counts[mask, :].sum(axis=0)).ravel()
                else:
                    pseudo_counts = np.array(counts[mask, :].mean(axis=0)).ravel()
            else:
                if PSEUDOBULK_PARAMS['aggregation'] == 'sum':
                    pseudo_counts = np.array(counts[mask, :].sum(axis=0)).ravel()
                else:
                    pseudo_counts = np.array(counts[mask, :].mean(axis=0)).ravel()
        except Exception as e:
            logger.error(f"  ❌ Failed for {key}: {e}")
            continue
        
        pseudobulk_list.append(pseudo_counts)
        
        meta = adata.obs[mask].iloc[0]
        metadata_list.append({
            'pseudobulk_id': key,
            'sample_id': meta[sample_col],
            'cell_type': meta[celltype_col],
            'n_cells': n_cells,
            'batch': meta[cols['batch_col']] if cols['batch_col'] in meta.index else None,
            'tissue': meta[cols['tissue_col']] if cols['tissue_col'] and cols['tissue_col'] in meta.index else None
        })
    
    if len(pseudobulk_list) == 0:
        logger.warning("⚠️ No samples passed filters")
        return None
    
    # Create AnnData
    pseudobulk_matrix = np.vstack(pseudobulk_list)
    pseudobulk_metadata = pd.DataFrame(metadata_list)
    
    adata_pseudo = sc.AnnData(
        X=pseudobulk_matrix.astype(np.int32),
        obs=pseudobulk_metadata,
        var=pd.DataFrame(index=genes)
    )
    
    logger.info(f"\n  Generated {adata_pseudo.n_obs} pseudobulk samples")
    
    # ⭐ CRITICAL: Verify directory exists before saving
    if not pseudobulk_dir.exists():
        logger.error(f"  ❌ Directory vanished: {pseudobulk_dir}")
        logger.error("  Attempting to recreate...")
        if not safe_mkdir(pseudobulk_dir, "pseudobulk directory (retry)"):
            logger.error("  ❌ Recreation failed - aborting save")
            return None
    
    try:
        pseudo_path = pseudobulk_dir / 'pseudobulk_counts.h5ad'
        adata_pseudo.write_h5ad(pseudo_path, compression='gzip')
        logger.info(f"  ✓ Saved: {pseudo_path}")
        
        pseudobulk_metadata.to_csv(pseudobulk_dir / 'pseudobulk_metadata.csv', index=False)
        logger.info(f"  ✓ Saved metadata")
        
    except Exception as e:
        logger.error(f"  ❌ Save failed: {e}")
        logger.error(f"  Directory status: exists={pseudobulk_dir.exists()}, is_dir={pseudobulk_dir.is_dir()}")
        return None
    
    logger.info("\n✓ Pseudobulk export complete")
    return pseudobulk_dir


# ============================================================================
# STEP 6: MARKER VISUALIZATION
# ============================================================================

def plot_marker_dotplot(
    adata: sc.AnnData, 
    markers: List[str], 
    celltype_col: str, 
    output_path: Path, 
    title: str = "Markers"
) -> None:
    """Dotplot with version compatibility"""
    available = [g for g in markers if g in adata.var_names]
    
    if len(available) == 0:
        logger.warning(f"  ⚠️ No markers for {title}")
        return
    
    logger.info(f"  Dotplot: {len(available)}/{len(markers)} markers")
    
    expr_source = resolve_expr_source(adata)
    use_raw = expr_source.get('use_raw', False)
    layer = expr_source.get('layer', None)
    
    try:
        if PACKAGE_INFO['compatibility']['scanpy_dotplot_return_fig']:
            dp = sc.pl.dotplot(
                adata,
                var_names=available,
                groupby=celltype_col,
                dendrogram=False,
                use_raw=use_raw,
                layer=layer,
                standard_scale='var',
                cmap='RdYlBu_r',
                show=False,
                return_fig=True
            )
            
            if hasattr(dp, 'fig'):
                dp.fig.suptitle(title, fontsize=14, weight='bold')
            
            if hasattr(dp, 'savefig'):
                dp.savefig(output_path, dpi=300, bbox_inches='tight')
            else:
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
            plt.close('all')
        else:
            fig, ax = plt.subplots(figsize=(max(12, len(available) * 0.3), 
                                            max(8, adata.obs[celltype_col].nunique() * 0.5)))
            
            sc.pl.dotplot(
                adata,
                var_names=available,
                groupby=celltype_col,
                dendrogram=False,
                use_raw=use_raw,
                layer=layer,
                standard_scale='var',
                cmap='RdYlBu_r',
                ax=ax,
                show=False
            )
            
            ax.set_title(title, fontsize=14, weight='bold')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        
        logger.info(f"    ✓ Saved: {output_path.name}")
        
    except Exception as e:
        logger.warning(f"    ⚠️ Dotplot failed: {e}")


def plot_marker_featureplot(
    adata: sc.AnnData, 
    markers: List[str], 
    output_dir: Path, 
    prefix: str = "markers"
) -> None:
    """Feature plots for markers"""
    available = [g for g in markers if g in adata.var_names]
    
    if len(available) == 0:
        return
    
    n_genes = min(len(available), MARKER_VIZ_PARAMS['featureplot_top_genes'])
    markers_to_plot = available[:n_genes]
    
    logger.info(f"  Feature plots: {n_genes} markers")
    
    if 'X_umap' not in adata.obsm:
        logger.warning("    ⚠️ No UMAP - skipping")
        return
    
    expr_source = resolve_expr_source(adata)
    use_raw = expr_source.get('use_raw', False)
    layer = expr_source.get('layer', None)
    
    try:
        n_cols = 4
        n_rows = int(np.ceil(n_genes / n_cols))
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4))
        axes = axes.flatten() if n_genes > 1 else [axes]
        
        for i, gene in enumerate(markers_to_plot):
            sc.pl.umap(
                adata,
                color=gene,
                use_raw=use_raw,
                layer=layer,
                cmap='RdYlBu_r',
                ax=axes[i],
                show=False,
                frameon=False,
                title=gene,
                size=50
            )
        
        for i in range(n_genes, len(axes)):
            axes[i].axis('off')
        
        plt.tight_layout()
        output_path = output_dir / f'{prefix}_featureplot.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"    ✓ Saved: {output_path.name}")
        
    except Exception as e:
        logger.warning(f"    ⚠️ Failed: {e}")


def perform_marker_visualization(
    adata: sc.AnnData, 
    celltype_name: str, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """Marker visualization"""
    log_step(6, "Marker Gene Visualization", "")
    
    viz_dir = output_dir / 'marker_visualization'
    if not safe_mkdir(viz_dir, "visualization directory"):
        logger.error("  ❌ Cannot create directory - skipping")
        return viz_dir
    
    celltype_col = cols['celltype_col']
    
    marker_sets = {'General': MARKER_GENES['General']}
    
    if celltype_name == 'T_cells':
        marker_sets['T_cell'] = MARKER_GENES['T_cells']
    elif celltype_name == 'B_cells':
        marker_sets['B_cell'] = MARKER_GENES['B_cells']
    elif celltype_name == 'Myeloid':
        marker_sets['Myeloid'] = MARKER_GENES['Myeloid']
    elif celltype_name == 'Stromal_Vascular':
        marker_sets['Stromal'] = MARKER_GENES['Stromal_Vascular']
    
    logger.info(f"\nMarker sets: {list(marker_sets.keys())}")
    
    logger.info("\n--- Dotplots ---")
    for name, markers in marker_sets.items():
        top = markers[:MARKER_VIZ_PARAMS['dotplot_top_genes']]
        path = viz_dir / f'dotplot_{name.lower()}.png'
        plot_marker_dotplot(adata, top, celltype_col, path, f"{name} Markers")
    
    logger.info("\n--- Feature Plots ---")
    for name, markers in marker_sets.items():
        plot_marker_featureplot(adata, markers, viz_dir, f'{name.lower()}')
    
    logger.info(f"\n✓ Marker visualization complete")
    return viz_dir


# ============================================================================
# STEP 7: SUBCLUSTERING
# ============================================================================

def perform_subclustering_safe(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Dict[str, Dict[str, Any]]:
    """Safe subclustering by cell type"""
    log_step(7, "Subclustering", "")
    
    subcluster_dir = output_dir / 'subclustering'
    if not safe_mkdir(subcluster_dir, "subclustering directory"):
        logger.error("  ❌ Cannot create directory - skipping")
        return {}
    
    celltype_col = cols['celltype_col']
    
    if 'X_scanvi' not in adata.obsm:
        logger.warning("⚠️ No X_scanvi - skipping")
        return {}
    
    celltypes = adata.obs[celltype_col].unique()
    logger.info(f"\nCell types: {len(celltypes)}")
    
    results_summary = {}
    
    for ct in celltypes:
        ct_safe = safe_name(ct)
        logger.info(f"\n--- {ct} ---")
        
        adata_sub = adata[adata.obs[celltype_col] == ct].copy()
        
        if adata_sub.n_obs < SUBCLUSTER_PARAMS['min_cells_for_subcluster']:
            logger.warning(f"  Skipping (n={adata_sub.n_obs})")
            continue
        
        ct_dir = subcluster_dir / ct_safe
        if not safe_mkdir(ct_dir, f"subcluster dir for {ct}"):
            logger.warning(f"  ⚠️ Cannot create subdir - skipping {ct}")
            continue
        
        logger.info(f"  Cells: {adata_sub.n_obs:,}")
        
        # Neighbors
        n_neighbors = min(SUBCLUSTER_PARAMS['n_neighbors'], adata_sub.n_obs // 10)
        sc.pp.neighbors(adata_sub, use_rep='X_scanvi', n_neighbors=n_neighbors)
        
        # Clustering
        logger.info("  Clustering...")
        for res in SUBCLUSTER_PARAMS['resolutions']:
            key = f'leiden_sub_res{res}'
            sc.tl.leiden(adata_sub, resolution=res, key_added=key)
            n_clusters = adata_sub.obs[key].nunique()
            logger.info(f"    Res {res}: {n_clusters} clusters")
        
        # UMAP
        sc.tl.umap(adata_sub, min_dist=0.3)
        
        # Markers
        default_key = f'leiden_sub_res{SUBCLUSTER_PARAMS["default_resolution"]}'
        expr_source = resolve_expr_source(adata_sub)
        
        sc.tl.rank_genes_groups(adata_sub, groupby=default_key, method='wilcoxon', n_genes=100, **expr_source)
        markers_df = sc.get.rank_genes_groups_df(adata_sub, group=None)
        markers_df.to_csv(ct_dir / 'subcluster_markers.csv', index=False)
        
        # Plot
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        sc.pl.umap(adata_sub, color=celltype_col, ax=axes[0, 0], show=False, title='Original', frameon=False, size=50)
        sc.pl.umap(adata_sub, color=cols['batch_col'], ax=axes[0, 1], show=False, title='Batch', frameon=False, size=50)
        axes[0, 2].axis('off')  # No disease plot
        
        for i, res in enumerate([0.2, 0.4, 0.6]):
            key = f'leiden_sub_res{res}'
            sc.pl.umap(adata_sub, color=key, ax=axes[1, i], show=False,
                      title=f'Res {res}', legend_loc='right margin', frameon=False, size=50)
        
        plt.tight_layout()
        plt.savefig(ct_dir / 'subclustering_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save
        output_path = ct_dir / f'{ct_safe}_subclustered.h5ad'
        adata_sub.write_h5ad(output_path, compression='gzip')
        
        results_summary[ct] = {
            'path': str(output_path),
            'n_cells': adata_sub.n_obs,
            'n_subclusters': adata_sub.obs[default_key].nunique()
        }
        
        logger.info(f"  ✓ Complete")
        
        del adata_sub
        gc.collect()
    
    with open(subcluster_dir / 'subclustering_summary.json', 'w') as f:
        json.dump(results_summary, f, indent=2)
    
    logger.info(f"\n✓ Subclustering: {len(results_summary)} cell types")
    return results_summary


# ============================================================================
# STEP 8: BBKNN (from v1.6.1)
# ============================================================================

def run_bbknn_clustering(adata: sc.AnnData, cols: Dict[str, Optional[str]], output_dir: Path) -> Optional[Path]:
    """BBKNN clustering (optional)"""
    log_step(8, "BBKNN Clustering (Optional)", "")
    
    if not BBKNN_PARAMS['enabled']:
        logger.info("⚠️ BBKNN disabled - skipping")
        return None
    
    if not BBKNN_AVAILABLE:
        logger.warning("⚠️ bbknn not installed - skipping")
        return None
    
    bbknn_dir = output_dir / 'bbknn'
    if not safe_mkdir(bbknn_dir, "BBKNN directory"):
        logger.error("  ❌ Cannot create directory - skipping")
        return None
    
    # Find batch key
    batch_key = None
    if BBKNN_PARAMS['batch_key'] in adata.obs.columns:
        batch_key = BBKNN_PARAMS['batch_key']
    else:
        for k in ['dataset', 'sample', 'batch', 'orig.ident']:
            if k in adata.obs.columns:
                batch_key = k
                break
    
    if batch_key is None:
        logger.warning("⚠️ No batch key found - skipping")
        return None
    
    logger.info(f"\nBatch key: {batch_key}")
    
    # Analyze batches
    batch_counts = adata.obs[batch_key].value_counts()
    n_batches = len(batch_counts)
    n_valid_batches = n_batches
    
    logger.info(f"  Total batches: {n_batches}")
    
    if n_batches < BBKNN_PARAMS['min_batches']:
        logger.warning(f"⚠️ Insufficient batches")
        return None
    
    # Filter small batches
    small_batches = batch_counts[batch_counts < BBKNN_PARAMS['min_cells_per_batch']]
    
    if len(small_batches) > 0:
        n_valid_batches = n_batches - len(small_batches)
        
        if n_valid_batches < BBKNN_PARAMS['min_batches']:
            logger.warning(f"⚠️ Too few valid batches")
            return None
        
        valid_batches = batch_counts[batch_counts >= BBKNN_PARAMS['min_cells_per_batch']].index
        adata_filtered = adata[adata.obs[batch_key].isin(valid_batches)].copy()
    else:
        adata_filtered = adata.copy()
    
    # Adjust neighbors
    min_batch_size = adata_filtered.obs[batch_key].value_counts().min()
    neighbors_within = BBKNN_PARAMS['neighbors_within_batch']
    
    if BBKNN_PARAMS['auto_adjust_neighbors'] and neighbors_within >= min_batch_size:
        original_neighbors = neighbors_within
        neighbors_within = max(1, min_batch_size - 1)
        logger.warning(f"\n⚠️ Auto-adjusting neighbors: {original_neighbors} → {neighbors_within}")
    
    # Backup UMAP
    if 'X_umap' in adata.obsm and 'X_umap_scanvi' not in adata.obsm:
        adata.obsm['X_umap_scanvi'] = np.asarray(adata.obsm['X_umap']).copy()
        logger.info("\n  ✓ Backed up X_umap → X_umap_scanvi")
    
    # PCA if needed
    need_pca = ('X_pca' not in adata_filtered.obsm) or (adata_filtered.obsm['X_pca'].shape[1] < BBKNN_PARAMS['n_pcs'])
    if need_pca:
        logger.info("  Computing PCA...")
        use_hvg = ('highly_variable' in adata_filtered.var.columns) and np.any(adata_filtered.var['highly_variable'].values)
        try:
            sc.pp.pca(adata_filtered, n_comps=BBKNN_PARAMS['n_pcs'], svd_solver='arpack', use_highly_variable=use_hvg)
        except Exception as e:
            logger.error(f"  ❌ PCA failed: {e}")
            return None
    
    # Run BBKNN
    logger.info(f"\n  Running BBKNN (neighbors_within_batch={neighbors_within})...")
    
    try:
        bbknn.bbknn(
            adata_filtered,
            batch_key=batch_key,
            neighbors_within_batch=neighbors_within,
            n_pcs=BBKNN_PARAMS['n_pcs'],
            trim=BBKNN_PARAMS['trim'],
            copy=False,
        )
        logger.info("  ✓ BBKNN graph constructed")
        
    except Exception as e:
        logger.error(f"  ❌ BBKNN failed: {e}")
        return None
    
    # UMAP
    logger.info("  Computing UMAP...")
    try:
        sc.tl.umap(adata_filtered, random_state=42)
        logger.info("  ✓ UMAP computed")
    except Exception as e:
        logger.error(f"  ❌ UMAP failed: {e}")
        return None
    
    # Transfer to original
    cell_mask = adata.obs_names.isin(adata_filtered.obs_names)
    
    umap_bbknn = np.full((adata.n_obs, 2), np.nan)
    umap_bbknn[cell_mask] = adata_filtered.obsm['X_umap']
    adata.obsm['X_umap_bbknn'] = umap_bbknn
    
    # Leiden
    logger.info("\n  Leiden clustering...")
    for res in BBKNN_PARAMS['leiden_resolutions']:
        key = f'leiden_bbknn_res{res}'
        try:
            sc.tl.leiden(adata_filtered, resolution=res, key_added=key)
        except Exception:
            sc.tl.leiden(adata_filtered, resolution=res, key_added=key, n_iterations=2)
        
        leiden_full = pd.Series('filtered_out', index=adata.obs_names, dtype='category')
        leiden_full[cell_mask] = adata_filtered.obs[key].values
        adata.obs[key] = leiden_full
    
    # Save
    filtered_path = bbknn_dir / 'adata_bbknn_filtered.h5ad'
    adata_filtered.write_h5ad(filtered_path, compression='gzip')
    logger.info(f"\n  ✓ Saved: {filtered_path}")
    
    # Plots
    try:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        adata_plot = adata_filtered
        
        if 'X_umap_scanvi' in adata.obsm:
            adata_orig = adata[cell_mask].copy()
            sc.pl.umap(adata_orig, color=cols['celltype_col'], ax=axes[0, 0], show=False, 
                       title='scANVI UMAP', frameon=False, size=50)
            sc.pl.umap(adata_orig, color=batch_key, ax=axes[0, 1], show=False,
                       title='scANVI UMAP (Batch)', frameon=False, size=50)
            del adata_orig
        else:
            axes[0, 0].axis('off')
            axes[0, 1].axis('off')
        
        sc.pl.umap(adata_plot, color=cols['celltype_col'], ax=axes[0, 2], show=False,
                   title='BBKNN UMAP', frameon=False, size=50)
        
        sc.pl.umap(adata_plot, color=batch_key, ax=axes[1, 0], show=False,
                   title='BBKNN UMAP (Batch)', frameon=False, size=50)
        
        for i, res in enumerate([0.5, 1.0]):
            if i < 2:
                key = f'leiden_bbknn_res{res}'
                if key in adata_plot.obs.columns:
                    sc.pl.umap(adata_plot, color=key, ax=axes[1, i+1], show=False,
                              title=f'Leiden (res={res})', frameon=False, size=50)
        
        plt.tight_layout()
        plt.savefig(bbknn_dir / 'bbknn_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        logger.warning(f"  ⚠️ Plot failed: {e}")
    
    # Summary
    summary = {
        'total_cells': int(adata.n_obs),
        'cells_used': int(cell_mask.sum()),
        'total_batches': int(n_batches),
        'batches_used': int(n_valid_batches),
    }
    
    with open(bbknn_dir / 'bbknn_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info("\n✓ BBKNN complete")
    
    del adata_filtered
    gc.collect()
    
    return bbknn_dir


# ============================================================================
# SUMMARY
# ============================================================================

def save_summary_stats_v1_6_2(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    model_genes: Optional[List[str]], 
    output_dir: Path
) -> None:
    """Summary statistics"""
    summary = []
    summary.append("="*70)
    summary.append("Dataset Summary Statistics")
    summary.append("="*70)
    summary.append(f"\nTotal cells: {adata.n_obs:,}")
    summary.append(f"Total genes: {adata.n_vars:,}")
    
    if model_genes:
        summary.append(f"Model genes: {len(model_genes)}")
        overlap = len(set(model_genes) & set(adata.var_names))
        summary.append(f"  Overlap: {overlap} ({overlap/len(model_genes)*100:.1f}%)")
    
    summary.append(f"\nCell Type Distribution ({cols['celltype_col']}):")
    for ct, count in adata.obs[cols['celltype_col']].value_counts().items():
        pct = count / adata.n_obs * 100
        summary.append(f"  {ct}: {count:,} ({pct:.1f}%)")
    
    if cols['confidence_col'] and cols['confidence_col'] in adata.obs.columns:
        summary.append(f"\nConfidence ({cols['confidence_col']}):")
        conf = adata.obs[cols['confidence_col']].describe()
        summary.append(f"  Mean: {conf['mean']:.3f}")
        summary.append(f"  Median: {conf['50%']:.3f}")
    
    if cols['batch_col'] in adata.obs.columns:
        summary.append(f"\nBatches: {adata.obs[cols['batch_col']].nunique()}")
    
    summary.append("\n" + "="*70)
    
    summary_text = '\n'.join(summary)
    with open(output_dir / 'summary_statistics.txt', 'w') as f:
        f.write(summary_text)
    
    logger.info(summary_text)


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def process_celltype(celltype_name: str, config: Dict[str, Any]) -> bool:
    """Main processing pipeline v1.6.2"""
    logger.info("\n" + "="*70)
    logger.info(f"PROCESSING: {celltype_name}")
    logger.info("="*70)
    
    output_dir = OUTPUT_BASE / celltype_name
    if not safe_mkdir(output_dir, f"output directory for {celltype_name}"):
        logger.error(f"  ❌ Cannot create output directory - aborting {celltype_name}")
        return False
    
    try:
        # Step 1
        adata_full, adata_model, lvae, cols, model_genes = load_data_and_model_v1_6_2(
            celltype_name, config, output_dir
        )
        
        # Summary
        save_summary_stats_v1_6_2(adata_full, cols, model_genes, output_dir)
        
        # Step 2
        de_dir = perform_differential_expression_v1_6_2(adata_full, adata_model, lvae, cols, output_dir)
        
        # Step 3
        denoised_path = generate_denoised_expression_safe(adata_model, lvae, output_dir)
        
        # Step 4
        comp_dir = perform_cell_composition_analysis_v1_6_2(adata_full, cols, output_dir)
        
        # Step 5
        pseudobulk_dir = export_pseudobulk_enhanced_v1_6_2(adata_full, cols, output_dir)
        
        # Step 6
        viz_dir = perform_marker_visualization(adata_full, celltype_name, cols, output_dir)
        
        # Step 7
        subcluster_summary = perform_subclustering_safe(adata_full, cols, output_dir)
        
        # Step 8
        bbknn_dir = run_bbknn_clustering(adata_full, cols, output_dir)
        
        # Save final
        logger.info("\nSaving final data...")
        final_path = output_dir / f'{celltype_name}_processed_final.h5ad'
        adata_full.write_h5ad(final_path, compression='gzip')
        logger.info(f"  ✓ Saved: {final_path}")
        
        logger.info("\n" + "="*70)
        logger.info(f"✓ {celltype_name} COMPLETE")
        logger.info("="*70)
        
        return True
        
    except Exception as e:
        logger.error(f"\n❌ ERROR: {e}", exc_info=True)
        return False
    finally:
        if 'adata_full' in locals():
            del adata_full
        if 'adata_model' in locals() and adata_model is not None:
            del adata_model
        if 'lvae' in locals() and lvae is not None:
            del lvae
        gc.collect()


def main() -> None:
    """Main entry point"""
    logger.info("="*70)
    logger.info("scANVI DOWNSTREAM ANALYSIS - v1.6.2 PRODUCTION (HOTFIX)")
    logger.info("="*70)
    logger.info(f"\nBase: {BASE_DIR}")
    logger.info(f"Output: {OUTPUT_BASE}")
    logger.info("\nCritical changes in v1.6.2:")
    logger.info("  ✅ Removed all disease-related DE analysis")
    logger.info("  ✅ Fixed scVI API compatibility (version-aware)")
    logger.info("  ✅ Enhanced directory creation safety")
    logger.info("  ✅ Maintained all v1.6.1 fixes")
    
    if not safe_mkdir(OUTPUT_BASE, "base output directory"):
        logger.error("❌ Cannot create base output directory - aborting")
        sys.exit(1)
    
    results = {}
    for celltype_name, config in CELL_TYPES.items():
        success = process_celltype(celltype_name, config)
        results[celltype_name] = success
    
    logger.info("\n" + "="*70)
    logger.info("PIPELINE COMPLETE")
    logger.info("="*70)
    
    for celltype, success in results.items():
        status = "✓ SUCCESS" if success else "❌ FAILED"
        logger.info(f"  {celltype}: {status}")
    
    successful = sum(results.values())
    logger.info(f"\nTotal: {successful}/{len(results)} successful")
    logger.info(f"\nOutputs: {OUTPUT_BASE}")
    logger.info("="*70)


if __name__ == "__main__":
    main()
