#!/usr/bin/env python3
"""
scANVI Downstream Analysis - Production Pipeline v1.6
======================================================

CHANGELOG v1.6.1 (2024-12-17) - HOTFIX:
- ✅ P0-1: FIXED: Gene order preservation (critical for scVI model compatibility)
- ✅ P0-2: FIXED: BBKNN n_valid_batches initialization
- ✅ P0-3: FIXED: BBKNN saves to X_umap_bbknn (preserves X_umap)
- ✅ P0-4: FIXED: Dotplot using return_fig pattern
- ✅ P0-Bio-1: FIXED: Disease DE explicit 2-group comparison
- ✅ P0-Bio-2: ADDED: Batch-disease audit in BBKNN
- ✅ P1-1: FIXED: Marker gene corrections (MFSD2A, ATF4, CD1C)
- ✅ P1-2: FIXED: Bidirectional LFC filtering (abs/up/down modes)
- ✅ P1-3: FIXED: Complete random seed setup (numpy, torch, scvi)
- ✅ P1-4: FIXED: Selective warning filtering

CHANGELOG v1.6 (2024-12-17):
- ✅ FIXED: confidence_col now optional (warning instead of error)
- ✅ FIXED: Single-sample group safety in DE analysis
- ✅ FIXED: Robust column inference for both naming modes
- ✅ ADDED: Optional BBKNN clustering (preserves existing UMAP)
- ✅ IMPROVED: Better error messages and validation
- ✅ MAINTAINED: All v1.5 features (marker visualization)

Critical features:
1. Robust model loading with HVG extraction
2. Single-sample safety checks in DE
3. Smart column inference (standard vs semantic modes)
4. Memory-safe denoised expression
5. Pseudobulk export for donor-level statistics
6. Comprehensive marker visualization
7. Optional BBKNN integration

Author: r2end
Date: 2024-12-17
Version: v1.6.1 (Production - HOTFIX)
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
import torch  # Add this import
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

# ⭐ FIX P1-4: Selective warning filtering
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)
# Keep UserWarning and RuntimeWarning visible

# ============================================================================
# PACKAGE VERSION CHECKING AND FALLBACK
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
    Check versions of critical packages and determine fallback strategies
    
    Returns:
    --------
    dict with version info and compatibility flags
    """
    versions = {}
    compatibility = {
        'scanpy_dotplot_return_fig': False,
        'scvi_has_seed_setting': False,
        'torch_weights_only': True,  # Newer torch versions require weights_only
    }
    
    # Scanpy version check
    scanpy_version = get_package_version('scanpy')
    versions['scanpy'] = scanpy_version
    if scanpy_version:
        try:
            # Check if dotplot supports return_fig (scanpy >= 1.9.0)
            major, minor = map(int, scanpy_version.split('.')[:2])
            if major > 1 or (major == 1 and minor >= 9):
                compatibility['scanpy_dotplot_return_fig'] = True
        except (ValueError, AttributeError):
            pass
    
    # scvi-tools version check
    scvi_version = get_package_version('scvi')
    versions['scvi'] = scvi_version
    if scvi_version:
        try:
            # Check if scvi has settings.seed (scvi-tools >= 0.15.0)
            major, minor = map(int, scvi_version.split('.')[:2])
            if major > 0 or (major == 0 and minor >= 15):
                compatibility['scvi_has_seed_setting'] = True
        except (ValueError, AttributeError):
            pass
    
    # Torch version check
    torch_version = get_package_version('torch')
    versions['torch'] = torch_version
    if torch_version:
        try:
            # torch.load weights_only parameter (torch >= 1.13.0)
            major, minor = map(int, torch_version.split('.')[:2])
            if major < 1 or (major == 1 and minor < 13):
                compatibility['torch_weights_only'] = False
        except (ValueError, AttributeError):
            pass
    
    # Other packages
    versions['pandas'] = get_package_version('pandas')
    versions['numpy'] = get_package_version('numpy')
    if BBKNN_AVAILABLE:
        versions['bbknn'] = get_package_version('bbknn')
    
    return {
        'versions': versions,
        'compatibility': compatibility
    }


# Global compatibility flags (set during initialization)
PACKAGE_INFO = check_package_versions()


def safe_torch_load(file_path: Path, map_location: str = 'cpu') -> Any:
    """
    ⭐ FALLBACK: Safe torch.load with version compatibility
    
    Handles differences in torch.load API across versions
    Note: Must be defined after PACKAGE_INFO initialization
    """
    try:
        if PACKAGE_INFO['compatibility']['torch_weights_only']:
            try:
                return torch.load(file_path, map_location=map_location, weights_only=False)
            except TypeError:
                # Fallback for older torch versions
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

# Initialize logger
logger = setup_logging()

# Log package versions
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
# MARKER GENE DEFINITIONS
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
    'Epithelial': [
        'TP63', 'KRT5', 'SCGB1A1', 'SERPINB3', 'SCGB3A2', 'SCGB3A1', 'TCN1',
        'ASRGL1', 'FOXJ1', 'RSPH1', 'PIFO', 'BEST4', 'C20orf85', 'C9orf24',  # ← NOTE: C20orf85 may also be CIMIP1
        'MUC5AC', 'SPDEF', 'LYPD2', 'ITLN1', 'ASCL1', 'GRP', 'POU2F3', 'ASCL2',
        'CFTR', 'FOXI1', 'ASCL3', 'BSND', 'IGF1', 'CLCNKB', 'AGER', 'RTKN2',
        'CLIC5', 'SPOCK2', 'TIMP3', 'SFTPC', 'LAMP3', 'MFSD2A', 'C8orf4',  # ← FIXED: MFSD2A
        'C11orf96', 'VIM', 'SOX9', 'KRT14', 'MYH11', 'ACTA2', 'DMBT1', 'RNASE1',
        'MUC5B', 'SPDEF', 'LYZ', 'LTF', 'SFTPB', 'SCGB3A2', 'SFTA2'
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
    'Myeloid': [
        'CLEC9A', 'XCR1', 'CADM1', 'CLNK', 'FLT3', 'ZBTB46', 'CLEC10A', 'CD1E',
        'FCER1A', 'CD1D', 'ITGAX', 'CD1C', 'FCGR2B', 'PKIB', 'CCR7', 'CD83',  # ← FIXED: CD1C
        'LAMP3', 'CCL22', 'CCL17', 'CCL19', 'LAD1', 'LILRA4', 'SMPD3', 'SCT',
        'IRF7', 'PLD4', 'CLEC4C', 'MARCO', 'FABP4', 'CYP27A1', 'SIGLEC1', 'ABCG1',
        'PPARG', 'C1QA', 'C1QB', 'C1QC', 'HLA-DPA1', 'SLC40A1', 'FOLR2', 'F13A1',
        'SPP1', 'HAMP', 'VCAN', 'CCR2', 'CCR5', 'FCN1', 'S100A12', 'RNASE2',
        'LILRA5', 'MTSS1', 'TPSAB1', 'MS4A2', 'TPSB2', 'FCGR3B', 'CSF3R', 'CXCR1'
    ],
    'B_cells': [
        'IGHD', 'TCL1A', 'MS4A1', 'BANK1', 'MKI67', 'TOP2A', 'IGHA1', 'IGHA2',
        'IGHGP', 'IGHG1'
    ],
    'Endothelial': [
        'S100B', 'ALDH1A1', 'GJA4', 'HEY1', 'DKK2', 'IGFBP3', 'EFNB2', 'ACKR1',
        'VWF', 'RGCC', 'VWA1', 'IL7R', 'FCN3', 'MT1M', 'TFF3', 'LYVE1', 'MMRN1',
        'CCL21', 'MYL9', 'ACTA2', 'TINAGL1', 'NOTCH3', 'LAMC3'
    ],
    'Fibroblast': [
        'APOD', 'FGF7', 'COL15A1', 'MFAP5', 'PI16', 'CD34', 'MMP11', 'COL10A1',
        'POSTN', 'LRRC15', 'HOPX', 'IGFBP5', 'TIMP1', 'MMP1', 'COL7A1', 'WNT5A',
        'ISG15', 'IL7R', 'SFRP4', 'SFRP2', 'COMP', 'RGS5', 'PDGFRB', 'NDUFA4L2',
        'NOTCH3', 'CXCL1', 'CXCL2', 'IL6', 'CEBPD', 'CLU', 'CTGF', 'HGF', 'HSPA6',
        'DNAJB1', 'MYC', 'ATF4', 'PLAU', 'CHI3L1', 'MMP3', 'IL1R1', 'IL13RA2',  # ← FIXED: ATF4
        'TNFSF11', 'MMP10', 'OSMR', 'IL11', 'STRA6', 'FAP', 'WNT2', 'TWIST1',
        'IL24', 'ACTG2', 'HHIP', 'CNN1', 'MYH11', 'ACTA2', 'TAGLN'
    ],
    'SMC': [
        'RGS5', 'CD36', 'NOTCH3', 'SCIN', 'EPAS1', 'HLA-C', 'IGKC', 'PTP4A3',
        'FN1', 'COL18A1', 'WFDC1', 'IGHG4', 'IGHG1', 'RERGL', 'PLN', 'SORBS2',
        'DSTN', 'TSC22D1', 'BCAM', 'C11orf96', 'FILIP1', 'WTIP', 'NRGN',
        'TMEM176B', 'ANGPTL1', 'CFH', 'FHL1', 'VCAN', 'GGT5', 'C1S', 'COL6A3',
        'STEAP4', 'ADGRL3'
    ],
    'Stromal_Vascular': [
        'APOD', 'FGF7', 'COL15A1', 'MFAP5', 'PI16', 'CD34', 'MMP11', 'COL10A1',
        'POSTN', 'LRRC15', 'HOPX', 'IGFBP5', 'TIMP1', 'MMP1', 'COL7A1', 'WNT5A',
        'ISG15', 'IL7R', 'SFRP4', 'SFRP2', 'COMP', 'RGS5', 'PDGFRB', 'NDUFA4L2',
        'NOTCH3', 'CXCL1', 'CXCL2', 'IL6', 'CEBPD', 'CLU', 'CTGF', 'HGF', 'HSPA6',
        'DNAJB1', 'MYC', 'AFT4', 'PLAU', 'CHI3L1', 'MMP3', 'IL1R1', 'IL13RA2',
        'TNFSF11', 'MMP10', 'OSMR', 'IL11', 'STRA6', 'FAP', 'WNT2', 'TWIST1',
        'IL24', 'ACTG2', 'HHIP', 'CNN1', 'MYH11', 'ACTA2', 'TAGLN', 'KRT18',
        'SLPI', 'UPK3B', 'MSLN', 'CALB2', 'WT1', 'KLK11', 'ITLN1', 'WSB1',
        'DDX17', 'CTNNB1', 'RBP1', 'STAR', 'STMN1', 'CXCL12', 'CD74', 'HLA-DRB1',
        'HLA-DRA', 'ADAMDEC1', 'CCL8', 'APOE', 'APOC1', 'LIMCH1', 'A2M', 'ADH1B',
        'PRG4', 'CRTAC1', 'CXCL14', 'VSTM2A', 'SOX6', 'COL4A5', 'COL4A6', 'TSLP',
        'FRZB', 'BMP5', 'BMP2', 'CPM', 'F3', 'RGS5', 'CD36', 'NOTCH3', 'SCIN',
        'EPAS1', 'HLA-C', 'IGKC', 'PTP4A3', 'FN1', 'COL18A1', 'WFDC1', 'IGHG4',
        'IGHG1', 'RERGL', 'PLN', 'SORBS2', 'DSTN', 'TSC22D1', 'BCAM', 'C11orf96',
        'FILIP1', 'WTIP', 'NRGN', 'TMEM176B', 'ANGPTL1', 'CFH', 'FHL1', 'VCAN',
        'GGT5', 'C1S', 'COL6A3', 'STEAP4', 'ADGRL3'
    ]
}

# ⭐ NEW: Add synonym mapping for robust marker detection
MARKER_GENE_SYNONYMS = {
    'MF5D2A': 'MFSD2A',
    'AFT4': 'ATF4',
    'CDIC': 'CD1C',
    'C20orf85': ['C20orf85', 'CIMIP1', 'LLC1'],  # Multiple possible names
}

# ============================================================================
# MARKER GENE VALIDATION
# ============================================================================

def validate_marker_genes(adata: sc.AnnData, marker_dict: Dict[str, List[str]]) -> Dict[str, Any]:
    """
    Pre-flight sanity check for marker genes
    
    Returns report with:
    - missing_genes: genes not in adata
    - possible_typos: genes with similar names in adata
    - duplicates: genes appearing multiple times
    """
    all_markers = []
    for markers in marker_dict.values():
        all_markers.extend(markers)
    
    missing = [g for g in all_markers if g not in adata.var_names]
    found = [g for g in all_markers if g in adata.var_names]
    duplicates = [g for g in set(all_markers) if all_markers.count(g) > 1]
    
    # Try synonym mapping for missing genes
    corrected = {}
    for gene in missing:
        if gene in MARKER_GENE_SYNONYMS:
            synonym = MARKER_GENE_SYNONYMS[gene]
            if isinstance(synonym, str) and synonym in adata.var_names:
                corrected[gene] = synonym
            elif isinstance(synonym, list):
                for s in synonym:
                    if s in adata.var_names:
                        corrected[gene] = s
                        break
    
    report = {
        'total_markers': len(all_markers),
        'found': len(found),
        'missing': len(missing),
        'duplicates': len(duplicates),
        'corrected': len(corrected),
        'missing_genes': missing,
        'duplicate_genes': duplicates,
        'corrections': corrected
    }
    
    logger.info(f"\n--- Marker Gene Validation ---")
    logger.info(f"  Total markers: {report['total_markers']}")
    logger.info(f"  Found: {report['found']} ({report['found']/report['total_markers']*100:.1f}%)")
    logger.info(f"  Missing: {report['missing']}")
    
    if corrected:
        logger.info(f"  Auto-corrected: {len(corrected)}")
        for old, new in corrected.items():
            logger.info(f"    {old} → {new}")
    
    if duplicates:
        logger.warning(f"  ⚠️  Duplicates: {duplicates}")
    
    return report

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path("/home/h2048/data/core_data")
OUTPUT_BASE = Path("/home/h2048/data/py/1217/downstream_analysis_v1_6")

# Cell types configuration
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

# Analysis parameters
DE_PARAMS = {
    'mode': 'change',
    'delta': 0.5,
    'batch_correction': True,
    'n_samples': 5000,
    'bayes_factor_threshold': 3.0,
    'lfc_threshold': 0.5,
    'min_proportion': 0.1,
    'proba_de_threshold': 0.95,
    'min_samples_per_group': 2,  # ⭐ NEW: Minimum samples for DE
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
    'save_separate': True,
    'save_in_layers': False
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
    'expression_cutoff': 0.1,
}

# ⭐ NEW: BBKNN Parameters
BBKNN_PARAMS = {
    'enabled': True,  # Toggle BBKNN on/off
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
    'min_cells_per_batch': 10,  # Minimum cells required per batch
    'min_batches': 2,  # Minimum batches required to run BBKNN
    'auto_adjust_neighbors': True,  # Auto-reduce neighbors_within_batch if needed
}

# Visualization
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# ⭐ FIX P1-3: Complete random seed setup with version fallback
np.random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)
# ⭐ FALLBACK: Only set scvi seed if supported
if PACKAGE_INFO['compatibility']['scvi_has_seed_setting']:
    try:
        scvi.settings.seed = 42
    except (AttributeError, TypeError):
        pass  # Fallback for older scvi-tools versions

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


# ============================================================================
# ENHANCED COLUMN INFERENCE (v1.6)
# ============================================================================

def infer_column_names_robust(
    adata: sc.AnnData, 
    naming_mode: str = 'standard'
) -> Dict[str, Optional[str]]:
    """
    ⭐ ENHANCED: Robust column inference supporting both naming modes
    
    Parameters:
    -----------
    adata : AnnData
    naming_mode : str
        'standard' for T/B/Myeloid
        'semantic' for Epithelial/Stromal
    
    Returns:
    --------
    dict: Column name mapping with None for missing optional columns
    """
    cols = {
        'celltype_col': None,
        'confidence_col': None,
        'batch_col': None,
        'sample_col': None,
        'disease_col': None,
        'tissue_col': None
    }
    
    # ===== CELLTYPE COLUMN =====
    if naming_mode == 'semantic':
        celltype_candidates = [
            'cell_type_scanvi_filt',      # Priority 1
            'cell_type_scanvi_raw',       # Priority 2
            'semantic_celltype',          # Priority 3
            'scanvi_predictions',         # Fallback
            'cell_type',                  # Generic
            'celltype'
        ]
    else:  # standard
        celltype_candidates = [
            'scanvi_predictions',         # Priority 1
            'scanvi_labels',              # Priority 2
            'predicted_labels',           # CellTypist original
            'majority_voting',            # CellTypist
            'cell_type',                  # Generic
            'celltype'
        ]
    
    for col in celltype_candidates:
        if col in adata.obs.columns:
            cols['celltype_col'] = col
            logger.info(f"  ✓ Found celltype_col: {col}")
            break
    
    # ===== CONFIDENCE COLUMN =====
    if naming_mode == 'semantic':
        confidence_candidates = [
            'scanvi_confidence',          # Priority 1
            'celltypist_confidence',      # Priority 2
            'conf_score',                 # CellTypist original
            'confidence'                  # Generic
        ]
    else:  # standard
        confidence_candidates = [
            'scanvi_confidence',          # Priority 1
            'conf_score',                 # Priority 2
            'confidence'                  # Generic
        ]
    
    for col in confidence_candidates:
        if col in adata.obs.columns:
            cols['confidence_col'] = col
            logger.info(f"  ✓ Found confidence_col: {col}")
            break
    
    # ===== BATCH COLUMN =====
    batch_candidates = ['dataset', 'batch', 'sample_id', 'orig.ident']
    for col in batch_candidates:
        if col in adata.obs.columns:
            cols['batch_col'] = col
            logger.info(f"  ✓ Found batch_col: {col}")
            break
    
    # ===== SAMPLE COLUMN =====
    sample_candidates = ['sample', 'sample_id', 'Sample', 'donor', 'donor_id', 'patient_id']
    for col in sample_candidates:
        if col in adata.obs.columns:
            cols['sample_col'] = col
            logger.info(f"  ✓ Found sample_col: {col}")
            break
    
    # ===== DISEASE COLUMN (Optional) =====
    disease_candidates = ['disease_status', 'Disease', 'condition', 'group', 'phenotype']
    for col in disease_candidates:
        if col in adata.obs.columns:
            cols['disease_col'] = col
            logger.info(f"  ✓ Found disease_col: {col}")
            break
    
    # ===== TISSUE COLUMN (Optional) =====
    tissue_candidates = ['tissue', 'organ', 'anatomical_site', 'location', 'organ__ontology_label']
    for col in tissue_candidates:
        if col in adata.obs.columns:
            cols['tissue_col'] = col
            logger.info(f"  ✓ Found tissue_col: {col}")
            break
    
    return cols


def validate_columns_v1_6(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]]
) -> bool:
    """
    ⭐ FIXED: Strict validation with confidence_col as OPTIONAL
    
    Raises ValueError only for truly required columns:
    - celltype_col
    - batch_col  
    - sample_col
    
    Warnings only for:
    - confidence_col (optional, used for QC plots)
    - disease_col (optional, needed for disease DE)
    - tissue_col (optional, needed for tissue comparison)
    """
    errors = []
    warnings = []
    
    # Required columns
    required = {
        'celltype_col': 'Cell type labels (e.g., scanvi_predictions or cell_type_scanvi_filt)',
        'batch_col': 'Batch/dataset identifier (e.g., dataset, batch)',
        'sample_col': 'Sample/donor identifier (e.g., sample_id, donor)'
    }
    
    for key, desc in required.items():
        col = cols.get(key)
        if not col:
            errors.append(f"Required column '{key}' not found. Needed: {desc}")
        elif col not in adata.obs.columns:
            errors.append(f"Column '{col}' specified for {key} not found in adata.obs")
    
    # Optional columns (warning only)
    optional = {
        'confidence_col': 'Confidence scores (e.g., scanvi_confidence, conf_score)',
        'disease_col': 'Disease status (e.g., disease_status, condition)',
        'tissue_col': 'Tissue location (e.g., tissue, organ)'
    }
    
    for key, desc in optional.items():
        col = cols.get(key)
        if not col:
            warnings.append(f"Optional column '{key}' not found. {desc} will be skipped")
        elif col not in adata.obs.columns:
            warnings.append(f"Column '{col}' for {key} specified but not in adata.obs")
    
    # Check for X_scanvi
    if 'X_scanvi' not in adata.obsm:
        warnings.append("X_scanvi not in adata.obsm - subclustering will be limited")
    
    # Check for counts
    has_counts = ('counts' in adata.layers) or (adata.raw is not None)
    if not has_counts:
        warnings.append("No counts matrix - pseudobulk/DE may be limited")
    
    # Log warnings
    if warnings:
        logger.warning("\n⚠️  Optional Features Limited:")
        for w in warnings:
            logger.warning(f"    {w}")
    
    # Raise errors
    if errors:
        logger.error("\n❌ REQUIRED COLUMNS MISSING:")
        for e in errors:
            logger.error(f"    {e}")
        logger.error(f"\nAvailable columns: {list(adata.obs.columns)[:20]}...")
        raise ValueError("Critical columns missing")
    
    logger.info("✓ Column validation passed")
    return True


# ============================================================================
# MODEL LOADING (from v1.5)
# ============================================================================

def extract_var_names_from_model(
    model_path: Union[str, Path], 
    adata_full: Optional[sc.AnnData] = None, 
    model_type: str = 'SCANVI'
) -> Optional[List[str]]:
    """Extract gene names from model (multiple methods)"""
    import torch
    
    try:
        # Method A: Direct checkpoint
        logger.info("  Method A: Reading checkpoint...")
        model_pt = Path(model_path) / 'model.pt'
        if model_pt.exists():
            # ⭐ FALLBACK: Use safe_torch_load for version compatibility
            checkpoint = safe_torch_load(model_pt, map_location='cpu')
            if 'var_names' in checkpoint:
                var_names = checkpoint['var_names']
                if isinstance(var_names, np.ndarray):
                    var_names = var_names.tolist()
                elif isinstance(var_names, torch.Tensor):
                    var_names = var_names.cpu().numpy().tolist()
                
                if len(var_names) > 0:
                    logger.info(f"  ✓ Extracted {len(var_names)} genes from checkpoint")
                    return list(var_names)
        
        # Method B: Load with adata
        if adata_full is not None:
            logger.info("  Method B: Loading model with adata...")
            Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
            model = Model.load(model_path, adata=adata_full)
            var_names = list(model.adata.var_names)
            logger.info(f"  ✓ Extracted {len(var_names)} genes from model.adata")
            return var_names
        
        # Method C: Legacy
        logger.info("  Method C: Legacy method...")
        Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
        model = Model.load(model_path, adata=None)
        vn_dict = model.get_var_names()
        if isinstance(vn_dict, dict):
            var_names = list(list(vn_dict.values())[0])
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
    ⭐ FIXED: Preserve gene order from model_genes (critical for scVI)
    
    Using set() scrambles order → model weights misaligned with genes
    """
    # ⭐ CRITICAL FIX: List comprehension maintains model_genes order
    overlap_genes = [g for g in model_genes if g in adata_full.var_names]
    missing_genes = [g for g in model_genes if g not in adata_full.var_names]
    
    overlap_ratio = len(overlap_genes) / len(model_genes) if len(model_genes) > 0 else 0.0
    
    logger.info(f"  Data genes: {adata_full.n_vars}")
    logger.info(f"  Model genes: {len(model_genes)}")
    logger.info(f"  Overlap: {len(overlap_genes)} ({overlap_ratio*100:.1f}%)")
    
    if len(missing_genes) > 0:
        logger.warning(f"  Missing: {len(missing_genes)} genes")
        if len(missing_genes) <= 20:
            logger.warning(f"    {missing_genes}")
    
    if overlap_ratio < VALIDATION_THRESHOLDS['min_gene_overlap_ratio']:
        logger.warning(f"  ⚠️  Low overlap may affect performance")
    
    # Subset with preserved order
    adata_model = adata_full[:, overlap_genes].copy()
    
    # Preserve .raw if exists
    if adata_full.raw is not None:
        adata_model.raw = sc.AnnData(
            X=adata_full.raw.X,
            obs=adata_full.obs.copy(),
            var=adata_full.raw.var.copy()
        )
    
    return adata_model


def load_data_and_model_v1_6(
    celltype_name: str, 
    config: Dict[str, Any], 
    output_dir: Path
) -> Tuple[sc.AnnData, Optional[sc.AnnData], Optional[Any], Dict[str, Optional[str]], Optional[List[str]]]:
    """Load data and model with robust column inference"""
    log_step(1, "Loading Data and Model", "")
    
    # Load data
    h5ad_path = BASE_DIR / config['h5ad']
    logger.info(f"Loading: {h5ad_path}")
    adata_full = sc.read_h5ad(h5ad_path)
    logger.info(f"  Loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars} genes")
    
    # ⭐ ENHANCED: Robust column inference
    cols = infer_column_names_robust(adata_full, config['naming_mode'])
    
    logger.info("\nColumn mapping:")
    for key, val in cols.items():
        status = val if val else "NOT FOUND"
        logger.info(f"  {key}: {status}")
    
    # ⭐ FIXED: New validation
    try:
        validate_columns_v1_6(adata_full, cols)
    except ValueError as e:
        logger.error(f"\n❌ Validation failed: {e}")
        return adata_full, None, None, cols, None
    
    # Load model
    scanvi_path = BASE_DIR / config['scanvi_model']
    logger.info(f"\nLoading scANVI model: {scanvi_path}")
    
    lvae = None
    adata_model = None
    model_genes = None
    
    try:
        model_genes = extract_var_names_from_model(scanvi_path, adata_full, 'SCANVI')
        if model_genes is None:
            raise ValueError("Could not extract genes")
        
        # Save gene list
        gene_list_path = output_dir / f'{celltype_name}_model_genes.txt'
        with open(gene_list_path, 'w') as f:
            for gene in model_genes:
                f.write(f"{gene}\n")
        logger.info(f"  ✓ Saved: {gene_list_path}")
        
        # Create compatible subset
        adata_model = create_model_compatible_adata(adata_full, model_genes)
        
        # Load model
        lvae = scvi.model.SCANVI.load(scanvi_path, adata=adata_model)
        logger.info("  ✓ Model loaded successfully")
        
    except Exception as e:
        logger.warning(f"\n⚠️ Model loading failed: {e}")
        logger.info("  Proceeding with scanpy-only analysis")
    
    return adata_full, adata_model, lvae, cols, model_genes


# ============================================================================
# EXPRESSION SOURCE HELPERS (from v1.5)
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
        logger.warning("  ⚠️  Using layers['counts'] (may need normalization)")
        return {'layer': 'counts'}
    
    logger.info("  Expression source: .X")
    return {'use_raw': False}


def get_counts_matrix(adata: sc.AnnData) -> Tuple[Any, pd.Index]:
    """Get counts matrix for pseudobulk"""
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
    
    logger.warning("  ⚠️  Using .X (may not be counts)")
    return adata.X, adata.var_names


# ============================================================================
# DE HELPERS WITH SINGLE-SAMPLE SAFETY (v1.6)
# ============================================================================

def _get_proba_de_col(df: pd.DataFrame) -> Optional[str]:
    """Get proba_de column name (version-agnostic)"""
    if 'proba_de' in df.columns:
        return 'proba_de'
    if 'proba_m2' in df.columns:
        return 'proba_m2'
    return None


def _compute_log2fc_from_means(df: pd.DataFrame, eps: float = 1e-8) -> Optional[np.ndarray]:
    """Compute log2FC from mean columns"""
    pairs = [
        ('raw_normalized_mean1', 'raw_normalized_mean2'),
        ('raw_mean1', 'raw_mean2'),
        ('scale1', 'scale2'),
    ]
    for a, b in pairs:
        if a in df.columns and b in df.columns:
            x1 = df[a].astype(float).to_numpy()
            x2 = df[b].astype(float).to_numpy()
            return np.log2((x1 + eps) / (x2 + eps))
    return None


def _filter_de_significant(de: pd.DataFrame, lfc_mode: str = 'up') -> pd.DataFrame:
    """
    ⭐ FIXED: Support bidirectional LFC filtering
    
    Parameters:
    -----------
    lfc_mode : str
        'up': only up-regulated (lfc > threshold)
        'down': only down-regulated (lfc < -threshold)
        'abs': both directions (|lfc| > threshold)
    """
    # Bayes factor
    bf_ok = None
    if 'bayes_factor' in de.columns:
        bf_ok = de['bayes_factor'].astype(float) > DE_PARAMS['bayes_factor_threshold']
    
    # Proba DE
    proba_col = _get_proba_de_col(de)
    proba_ok = None
    if proba_col:
        proba_ok = de[proba_col].astype(float) > DE_PARAMS['proba_de_threshold']
    
    # LFC with mode support
    if 'lfc_mean' in de.columns:
        lfc = de['lfc_mean'].astype(float)
    else:
        lfc_array = _compute_log2fc_from_means(de)
        if lfc_array is not None:
            lfc = pd.Series(lfc_array, index=de.index, name='lfc_approx')
        else:
            lfc = None
    
    # ⭐ FIXED: Apply mode-specific LFC filter
    lfc_ok = None
    if lfc is not None:
        threshold = DE_PARAMS['lfc_threshold']
        if lfc_mode == 'abs':
            lfc_ok = lfc.abs() > threshold
        elif lfc_mode == 'up':
            lfc_ok = lfc > threshold
        elif lfc_mode == 'down':
            lfc_ok = lfc < -threshold
        else:
            raise ValueError(f"Invalid lfc_mode: {lfc_mode}")
    
    # Combine filters
    mask = None
    for m in [bf_ok, proba_ok, lfc_ok]:
        if m is None:
            continue
        mask = m if mask is None else (mask & m)
    
    if mask is None:
        return de.iloc[0:0].copy()
    
    de_sig = de.loc[mask].copy()
    
    # Expression proportion filter
    if 'non_zeros_proportion1' in de_sig.columns:
        de_sig = de_sig[de_sig['non_zeros_proportion1'].astype(float) > DE_PARAMS['min_proportion']]
    
    # Add computed LFC if needed
    if 'lfc_mean' not in de_sig.columns and lfc is not None:
        if isinstance(lfc, pd.Series):
            de_sig = de_sig.join(lfc.rename('lfc_approx_log2'), how='left')
    
    return de_sig


def validate_and_fix_de_schema(de_df: pd.DataFrame, source: str = "scANVI", var_names: Optional[List[str]] = None) -> pd.DataFrame:
    """Validate DE schema"""
    if not isinstance(de_df, pd.DataFrame):
        raise ValueError(f"DE result not DataFrame (got {type(de_df)})")
    
    # Reset index if gene-like
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
        
        if ('bayes_factor' not in de_df.columns) and (_get_proba_de_col(de_df) is None):
            raise ValueError(f"Missing bayes_factor or proba_de. Columns: {list(de_df.columns)[:30]}")
    
    return de_df


def check_group_sample_counts(adata: sc.AnnData, groupby: str, sample_col: str, min_samples: int = 2) -> Tuple[bool, Dict[str, int]]:
    """
    ⭐ NEW: Check if each group has enough samples for DE
    
    Returns:
    --------
    valid : bool
        True if all groups have >= min_samples
    sample_counts : dict
        Group -> sample count mapping
    """
    if groupby not in adata.obs.columns or sample_col not in adata.obs.columns:
        return False, {}
    
    sample_counts = {}
    for group in adata.obs[groupby].unique():
        mask = adata.obs[groupby] == group
        n_samples = adata.obs.loc[mask, sample_col].nunique()
        sample_counts[str(group)] = n_samples
    
    valid = all(n >= min_samples for n in sample_counts.values())
    return valid, sample_counts


# ============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION (v1.6 - SINGLE-SAMPLE SAFE)
# ============================================================================

def perform_de_one_vs_rest_explicit(
    adata_model: sc.AnnData, 
    lvae: Any, 
    celltype_col: str, 
    output_dir: Path
) -> Optional[pd.DataFrame]:
    """One-vs-rest DE with explicit loop"""
    logger.info("\n--- One-vs-Rest DE ---")
    
    celltypes = adata_model.obs[celltype_col].cat.categories if \
                hasattr(adata_model.obs[celltype_col], 'cat') else \
                adata_model.obs[celltype_col].unique()
    
    all_results = []
    
    for ct in celltypes:
        logger.info(f"\n  Processing: {ct}")
        
        try:
            de = lvae.differential_expression(
                groupby=celltype_col,
                group1=ct,
                group2=None,
                mode=DE_PARAMS['mode'],
                delta=DE_PARAMS['delta'],
                batch_correction=DE_PARAMS['batch_correction'],
                n_samples=DE_PARAMS['n_samples']
            )
            
            de = validate_and_fix_de_schema(de, source="scANVI", var_names=adata_model.var_names)
            de['cell_type'] = ct
            de['comparison'] = f'{ct}_vs_Rest'
            
            all_results.append(de)
            
            de_sig = _filter_de_significant(de)
            logger.info(f"    Significant: {len(de_sig)}")
            
            safe_ct = safe_name(ct)
            de_sig.to_csv(output_dir / f'{safe_ct}_markers.csv', index=False)
            
        except Exception as e:
            logger.warning(f"    ⚠️ Failed: {e}")
            continue
    
    if all_results:
        de_all = pd.concat(all_results, ignore_index=True)
        de_all.to_csv(output_dir / 'celltype_markers_all_scANVI.csv', index=False)
        logger.info(f"\n  ✓ Saved: celltype_markers_all_scANVI.csv")
        return de_all
    
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


def perform_de_disease_within_celltype_safe(
    adata_for_de: sc.AnnData, 
    lvae: Optional[Any], 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> None:
    """
    ⭐ FIXED: Explicit 2-group disease comparison with clear labeling
    """
    if not cols['disease_col'] or cols['disease_col'] not in adata_for_de.obs.columns:
        logger.warning("\n⚠️ No disease column - skipping")
        return
    
    if not cols['sample_col'] or cols['sample_col'] not in adata_for_de.obs.columns:
        logger.warning("\n⚠️ No sample column - cannot verify sample counts")
        return
    
    logger.info("\n--- Disease DE (Within Cell Type) ---")
    logger.warning("⚠️ NOTE: Cell-level DE - see pseudobulk for donor-level")
    
    celltype_col = cols['celltype_col']
    disease_col = cols['disease_col']
    sample_col = cols['sample_col']
    
    for ct in adata_for_de.obs[celltype_col].unique():
        logger.info(f"\n  {ct}...")
        adata_ct = adata_for_de[adata_for_de.obs[celltype_col] == ct].copy()
        
        if adata_ct.n_obs < VALIDATION_THRESHOLDS['min_cells_for_de']:
            logger.warning(f"    Skipping (n={adata_ct.n_obs})")
            continue
        
        # ⭐ FIXED: Enforce exactly 2 disease states
        disease_states = sorted(pd.unique(adata_ct.obs[disease_col].astype(str)))
        
        if len(disease_states) < 2:
            logger.warning(f"    Skipping (only 1 disease state: {disease_states})")
            continue
        elif len(disease_states) > 2:
            logger.warning(f"    Skipping (>2 disease states: {disease_states})")
            logger.warning(f"    For multi-group, please pre-merge disease_col or use pairwise")
            continue
        
        # Sample count check
        valid, sample_counts = check_group_sample_counts(
            adata_ct, 
            disease_col, 
            sample_col, 
            min_samples=DE_PARAMS['min_samples_per_group']
        )
        
        if not valid:
            logger.warning(f"    Skipping - insufficient samples:")
            for group, count in sample_counts.items():
                logger.warning(f"      {group}: {count} samples")
            continue
        
        # ⭐ FIXED: Explicit group1/group2
        state1, state2 = disease_states
        logger.info(f"    Comparing: {state1} vs {state2}")
        logger.info(f"    Sample counts: {sample_counts}")
        
        safe_ct = safe_name(ct)
        comparison_name = f'{safe_name(state1)}_vs_{safe_name(state2)}'
        
        if lvae is not None:
            try:
                # ⭐ CRITICAL: Explicit group1 and group2
                de = lvae.differential_expression(
                    adata=adata_ct,
                    groupby=disease_col,
                    group1=state1,  # ← Explicit
                    group2=state2,  # ← Explicit
                    mode=DE_PARAMS['mode'],
                    delta=0.25,
                    batch_correction=True,
                    n_samples=DE_PARAMS['n_samples']
                )
                
                de = validate_and_fix_de_schema(de, source="scANVI", var_names=adata_ct.var_names)
                
                # ⭐ FIXED: Clear filename with comparison
                de.to_csv(output_dir / f'{safe_ct}_{comparison_name}_scANVI.csv', index=False)
                
                # ⭐ FIXED: Use 'abs' mode for disease DE (up and down)
                de_sig = _filter_de_significant(de, lfc_mode='abs')
                logger.info(f"    ✓ Significant: {len(de_sig)} (bidirectional)")
                
            except Exception as e:
                logger.warning(f"    ⚠️ scANVI DE failed: {e}")
        else:
            # Scanpy fallback
            try:
                expr_source = resolve_expr_source(adata_ct)
                sc.tl.rank_genes_groups(adata_ct, groupby=disease_col, method='wilcoxon', **expr_source)
                de = sc.get.rank_genes_groups_df(adata_ct, group=None)
                
                de.to_csv(output_dir / f'{safe_ct}_{comparison_name}_scanpy.csv', index=False)
                logger.info(f"    ✓ Scanpy DE done")
                
            except Exception as e:
                logger.warning(f"    ⚠️ Scanpy DE failed: {e}")


def perform_differential_expression(
    adata_full: sc.AnnData, 
    adata_model: Optional[sc.AnnData], 
    lvae: Optional[Any], 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """Complete DE workflow"""
    log_step(2, "Differential Expression Analysis", "")
    
    de_dir = output_dir / 'differential_expression'
    de_dir.mkdir(exist_ok=True, parents=True)
    
    celltype_col = cols['celltype_col']
    adata_for_de = adata_model if adata_model is not None else adata_full
    
    # One-vs-rest
    if lvae is not None and adata_model is not None:
        perform_de_one_vs_rest_explicit(adata_model, lvae, celltype_col, de_dir)
    else:
        perform_de_scanpy_fallback(adata_full, celltype_col, de_dir)
    
    # ⭐ FIXED: Disease DE with safety checks
    perform_de_disease_within_celltype_safe(adata_for_de, lvae, cols, de_dir)
    
    logger.info("\n✓ Differential expression complete")
    return de_dir


# ============================================================================
# STEP 3: DENOISED EXPRESSION (from v1.5)
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
        logger.warning(f"⚠️ Dataset too large - skipping")
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
# STEP 4: CELL COMPOSITION (from v1.5)
# ============================================================================

def perform_cell_composition_analysis(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """Cell composition analysis"""
    log_step(4, "Cell Composition Analysis", "")
    
    comp_dir = output_dir / 'cell_composition'
    comp_dir.mkdir(exist_ok=True, parents=True)
    
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
    
    # Disease comparison
    if cols['disease_col'] and cols['disease_col'] in adata.obs.columns:
        logger.info("\n--- Disease Composition ---")
        
        sample_disease = adata.obs[[sample_col, cols['disease_col']]].drop_duplicates()
        sample_disease = sample_disease.set_index(sample_col)
        composition_with_disease = composition.join(sample_disease)
        
        comp_by_disease = composition_with_disease.groupby(cols['disease_col']).mean()
        comp_by_disease.to_csv(comp_dir / 'composition_by_disease.csv')
        
        # Stats
        comp_stats = []
        disease_states = adata.obs[cols['disease_col']].unique()
        
        if len(disease_states) == 2:
            state1, state2 = disease_states
            
            for ct in composition.columns:
                group1 = composition_with_disease[
                    composition_with_disease[cols['disease_col']] == state1
                ][ct]
                group2 = composition_with_disease[
                    composition_with_disease[cols['disease_col']] == state2
                ][ct]
                
                if len(group1) > 0 and len(group2) > 0:
                    u_stat, p_val = mannwhitneyu(group1, group2)
                    
                    comp_stats.append({
                        'cell_type': ct,
                        f'{state1}_mean': group1.mean(),
                        f'{state2}_mean': group2.mean(),
                        'fold_change': group1.mean() / (group2.mean() + 1e-10),
                        'p_value': p_val
                    })
            
            if comp_stats:
                comp_stats_df = pd.DataFrame(comp_stats)
                comp_stats_df['p_adj'] = multipletests(comp_stats_df['p_value'], method='fdr_bh')[1]
                comp_stats_df = comp_stats_df.sort_values('p_adj')
                comp_stats_df.to_csv(comp_dir / 'composition_disease_statistics.csv', index=False)
                
                sig = comp_stats_df[comp_stats_df['p_adj'] < 0.05]
                if len(sig) > 0:
                    logger.info(f"\n  Significant changes: {len(sig)}")
    
    logger.info("\n✓ Cell composition complete")
    return comp_dir


# ============================================================================
# STEP 5: PSEUDOBULK (from v1.5 - with bug fixes)
# ============================================================================

def export_pseudobulk_enhanced(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Optional[Path]:
    """Enhanced pseudobulk export"""
    log_step(5, "Pseudobulk Export", "")
    
    pseudobulk_dir = output_dir / 'pseudobulk'
    pseudobulk_dir.mkdir(exist_ok=True, parents=True)
    
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
        
        # Aggregate
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
            'disease_status': meta[cols['disease_col']] if cols['disease_col'] and cols['disease_col'] in meta.index else None,
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
    
    pseudo_path = pseudobulk_dir / 'pseudobulk_counts.h5ad'
    adata_pseudo.write_h5ad(pseudo_path, compression='gzip')
    pseudobulk_metadata.to_csv(pseudobulk_dir / 'pseudobulk_metadata.csv', index=False)
    
    logger.info(f"  ✓ Saved: {pseudo_path}")
    logger.info("\n⚠️ For tissue comparison, use DESeq2/edgeR on this data")
    
    return pseudobulk_dir


# ============================================================================
# STEP 6: MARKER VISUALIZATION (from v1.5)
# ============================================================================

def plot_marker_dotplot(
    adata: sc.AnnData, 
    markers: List[str], 
    celltype_col: str, 
    output_path: Path, 
    title: str = "Markers"
) -> None:
    """
    ⭐ FIXED: Use scanpy's return_fig pattern instead of ax
    """
    available = [g for g in markers if g in adata.var_names]
    
    if len(available) == 0:
        logger.warning(f"  ⚠️ No markers for {title}")
        return
    
    logger.info(f"  Dotplot: {len(available)}/{len(markers)} markers")
    
    expr_source = resolve_expr_source(adata)
    use_raw = expr_source.get('use_raw', False)
    layer = expr_source.get('layer', None)
    
    try:
        # ⭐ FALLBACK: Use return_fig if supported, otherwise use ax
        if PACKAGE_INFO['compatibility']['scanpy_dotplot_return_fig']:
            # Newer scanpy: use return_fig
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
            
            # Set title on the figure
            if hasattr(dp, 'fig'):
                dp.fig.suptitle(title, fontsize=14, weight='bold')
            
            # Save using dotplot object's method
            if hasattr(dp, 'savefig'):
                dp.savefig(output_path, dpi=300, bbox_inches='tight')
            else:
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
            
            plt.close('all')
        else:
            # Older scanpy: use ax parameter
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
        # Try fallback with matrixplot
        try:
            logger.info(f"    Trying matrixplot fallback...")
            sc.pl.matrixplot(
                adata,
                var_names=available,
                groupby=celltype_col,
                dendrogram=False,
                use_raw=use_raw,
                layer=layer,
                standard_scale='var',
                cmap='RdYlBu_r',
                show=False,
                save=str(output_path)
            )
            logger.info(f"    ✓ Matrixplot saved")
        except Exception as e2:
            logger.error(f"    ❌ Both dotplot and matrixplot failed: {e2}")


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
    viz_dir.mkdir(exist_ok=True, parents=True)
    
    celltype_col = cols['celltype_col']
    
    # Select markers
    marker_sets = {'General': MARKER_GENES['General']}
    
    if celltype_name == 'T_cells':
        marker_sets['T_cell'] = MARKER_GENES['T_cells']
    elif celltype_name == 'B_cells':
        marker_sets['B_cell'] = MARKER_GENES['B_cells']
    elif celltype_name == 'Myeloid':
        marker_sets['Myeloid'] = MARKER_GENES['Myeloid']
    elif celltype_name == 'Stromal_Vascular':
        marker_sets['Stromal'] = MARKER_GENES['Stromal_Vascular']
        marker_sets['Fibroblast'] = MARKER_GENES['Fibroblast']
        marker_sets['SMC'] = MARKER_GENES['SMC']
        marker_sets['Endothelial'] = MARKER_GENES['Endothelial']
    
    logger.info(f"\nMarker sets: {list(marker_sets.keys())}")
    
    # Dotplots
    logger.info("\n--- Dotplots ---")
    for name, markers in marker_sets.items():
        top = markers[:MARKER_VIZ_PARAMS['dotplot_top_genes']]
        path = viz_dir / f'dotplot_{name.lower()}.png'
        plot_marker_dotplot(adata, top, celltype_col, path, f"{name} Markers")
    
    # Feature plots
    logger.info("\n--- Feature Plots ---")
    for name, markers in marker_sets.items():
        plot_marker_featureplot(adata, markers, viz_dir, f'{name.lower()}')
    
    logger.info(f"\n✓ Marker visualization complete")
    return viz_dir


# ============================================================================
# STEP 7: SUBCLUSTERING (from v1.5)
# ============================================================================

def perform_subclustering_safe(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Dict[str, Dict[str, Any]]:
    """Safe subclustering by cell type"""
    log_step(7, "Subclustering", "")
    
    subcluster_dir = output_dir / 'subclustering'
    subcluster_dir.mkdir(exist_ok=True, parents=True)
    
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
        ct_dir.mkdir(exist_ok=True, parents=True)
        
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
        
        if cols['disease_col'] and cols['disease_col'] in adata_sub.obs.columns:
            sc.pl.umap(adata_sub, color=cols['disease_col'], ax=axes[0, 2], show=False, title='Disease', frameon=False, size=50)
        else:
            axes[0, 2].axis('off')
        
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
# STEP 8: BBKNN CLUSTERING (NEW in v1.6)
# ============================================================================

def run_bbknn_clustering(adata: sc.AnnData, cols: Dict[str, Optional[str]], output_dir: Path) -> Optional[Path]:
    """
    ⭐ FIXED v1.6.1:
    - P0-2: Initialize n_valid_batches properly
    - P0-3: Save to X_umap_bbknn (not X_umap)
    - P0-Bio-2: Add batch-disease audit
    """
    log_step(8, "BBKNN Clustering (Optional)", "")
    
    if not BBKNN_PARAMS['enabled']:
        logger.info("⚠️ BBKNN disabled - skipping")
        return None
    
    if not BBKNN_AVAILABLE:
        logger.warning("⚠️ bbknn not installed - skipping")
        return None
    
    bbknn_dir = output_dir / 'bbknn'
    bbknn_dir.mkdir(exist_ok=True, parents=True)
    
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
    
    # ⭐ NEW: Batch-disease audit (P0-Bio-2)
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
        logger.warning(f"\n{batch_disease}")
        
        # Check if batch is confounded with disease
        confounded_batches = []
        for batch in batch_disease.index[:-1]:  # Exclude 'All' row
            disease_counts = batch_disease.loc[batch, :]
            disease_counts = disease_counts[disease_counts.index != 'All']
            if (disease_counts > 0).sum() == 1:
                confounded_batches.append(batch)
        
        if confounded_batches:
            logger.warning(f"\n⚠️  WARNING: {len(confounded_batches)} batches are disease-specific!")
            logger.warning(f"   {confounded_batches[:5]}")
            logger.warning("   BBKNN may remove real disease biology!")
            logger.warning("   Consider: Use scVI latent space for disease analysis")
        
        logger.warning("="*70)
    
    # Analyze batch composition
    batch_counts = adata.obs[batch_key].value_counts()
    n_batches = len(batch_counts)
    
    # ⭐ FIX P0-2: Initialize n_valid_batches here
    n_valid_batches = n_batches
    
    logger.info(f"\nBatch composition:")
    logger.info(f"  Total batches: {n_batches}")
    
    # Safety checks...
    if n_batches < BBKNN_PARAMS['min_batches']:
        logger.warning(f"⚠️ Insufficient batches")
        return None
    
    # Filter small batches
    small_batches = batch_counts[batch_counts < BBKNN_PARAMS['min_cells_per_batch']]
    
    if len(small_batches) > 0:
        n_valid_batches = n_batches - len(small_batches)  # ← Update here
        
        if n_valid_batches < BBKNN_PARAMS['min_batches']:
            logger.warning(f"⚠️ Too few valid batches after filtering")
            return None
        
        valid_batches = batch_counts[batch_counts >= BBKNN_PARAMS['min_cells_per_batch']].index
        adata_filtered = adata[adata.obs[batch_key].isin(valid_batches)].copy()
    else:
        adata_filtered = adata.copy()
    
    # ===== STEP 5: Adjust neighbors_within_batch if needed =====
    min_batch_size = adata_filtered.obs[batch_key].value_counts().min()
    neighbors_within = BBKNN_PARAMS['neighbors_within_batch']
    
    if BBKNN_PARAMS['auto_adjust_neighbors'] and neighbors_within >= min_batch_size:
        original_neighbors = neighbors_within
        neighbors_within = max(1, min_batch_size - 1)  # At least 1, at most (min_batch_size - 1)
        logger.warning(f"\n⚠️ Auto-adjusting neighbors_within_batch:")
        logger.warning(f"    {original_neighbors} → {neighbors_within} (min batch size: {min_batch_size})")
    
    # ===== STEP 6: Backup existing UMAP =====
    if 'X_umap' in adata.obsm and 'X_umap_scanvi' not in adata.obsm:
        adata.obsm['X_umap_scanvi'] = np.asarray(adata.obsm['X_umap']).copy()
        logger.info("\n  ✓ Backed up X_umap → X_umap_scanvi")
    
    # ===== STEP 7: PCA if needed =====
    need_pca = ('X_pca' not in adata_filtered.obsm) or (adata_filtered.obsm['X_pca'].shape[1] < BBKNN_PARAMS['n_pcs'])
    if need_pca:
        logger.info("  Computing PCA...")
        use_hvg = ('highly_variable' in adata_filtered.var.columns) and np.any(adata_filtered.var['highly_variable'].values)
        try:
            sc.pp.pca(adata_filtered, n_comps=BBKNN_PARAMS['n_pcs'], svd_solver='arpack', use_highly_variable=use_hvg)
        except Exception as e:
            logger.error(f"  ❌ PCA failed: {e}")
            return None
    
    # ===== STEP 8: Run BBKNN (with error handling) =====
    logger.info(f"\n  Running BBKNN (neighbors_within_batch={neighbors_within})...")
    
    try:
        bbknn.bbknn(
            adata_filtered,
            batch_key=batch_key,
            neighbors_within_batch=neighbors_within,
            n_pcs=BBKNN_PARAMS['n_pcs'],
            metric=BBKNN_PARAMS['metric'],
            trim=BBKNN_PARAMS['trim'],
            key_added=BBKNN_PARAMS['neighbors_key'],
            copy=False,
        )
        logger.info("  ✓ BBKNN graph constructed")
        
    except ValueError as e:
        logger.error(f"  ❌ BBKNN failed: {e}")
        logger.error("  Possible causes:")
        logger.error("    - neighbors_within_batch too large for smallest batch")
        logger.error("    - Incompatible batch structure")
        logger.error("  Suggestion: Try reducing min_cells_per_batch or neighbors_within_batch")
        return None
    except Exception as e:
        logger.error(f"  ❌ Unexpected BBKNN error: {e}")
        return None
    
    # ===== STEP 9: UMAP =====
    logger.info("  Computing UMAP...")
    try:
        sc.tl.umap(
            adata_filtered,
            neighbors_key=BBKNN_PARAMS['neighbors_key'],
            min_dist=BBKNN_PARAMS['umap_min_dist'],
            spread=BBKNN_PARAMS['umap_spread'],
            random_state=42,
        )
        logger.info("  ✓ UMAP computed")
    except Exception as e:
        logger.error(f"  ❌ UMAP failed: {e}")
        return None
    
    # ===== STEP 10: Transfer results to original adata =====
    # Transfer UMAP and neighbors to original adata (for cells that were kept)
    cell_mask = adata.obs_names.isin(adata_filtered.obs_names)
    n_transferred = cell_mask.sum()
    
    # Create BBKNN UMAP array
    umap_bbknn = np.full((adata.n_obs, 2), np.nan)
    umap_bbknn[cell_mask] = adata_filtered.obsm['X_umap']
    
    # ⭐ CRITICAL FIX: Save as X_umap_bbknn, keep X_umap intact
    adata.obsm['X_umap_bbknn'] = umap_bbknn
    
    # Transfer neighbors (only for included cells - this is a limitation)
    # Note: Neighbors graph cannot be easily transferred with NaN cells
    # So we save filtered results separately
    
    # ===== STEP 11: Leiden clustering =====
    logger.info("\n  Leiden clustering...")
    for res in BBKNN_PARAMS['leiden_resolutions']:
        key = f'leiden_bbknn_res{res}'
        try:
            sc.tl.leiden(
                adata_filtered,
                neighbors_key=BBKNN_PARAMS['neighbors_key'],
                resolution=res,
                key_added=key,
                flavor='igraph',
                n_iterations=2,
                directed=False
            )
        except TypeError:
            sc.tl.leiden(
                adata_filtered,
                neighbors_key=BBKNN_PARAMS['neighbors_key'],
                resolution=res,
                key_added=key,
                n_iterations=2
            )
        
        # Transfer to full adata
        leiden_full = pd.Series('filtered_out', index=adata.obs_names, dtype='category')
        leiden_full[cell_mask] = adata_filtered.obs[key].values
        adata.obs[key] = leiden_full
    
    default_key = f'leiden_bbknn_res{BBKNN_PARAMS["default_leiden_res"]}'
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
        
        # Use filtered adata for plotting (to avoid NaN issues)
        adata_plot = adata_filtered
        
        # Row 1: Original vs BBKNN UMAP
        if 'X_umap_scanvi' in adata.obsm:
            # Plot original UMAP (from full adata, filtered cells)
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
        sc.pl.umap(adata_plot, color=cols['celltype_col'], ax=axes[0, 2], show=False,
                   title='BBKNN UMAP', frameon=False, size=50)
        
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
        logger.warning(f"  ⚠️ Plot generation failed: {e}")
    
    # ===== STEP 14: Save summary =====
    summary = {
        'total_cells': int(adata.n_obs),
        'cells_used': int(n_transferred),
        'cells_filtered': int(adata.n_obs - n_transferred),
        'total_batches': int(n_batches),
        'batches_used': int(n_valid_batches),
        'batches_filtered': int(len(small_batches)),
        'neighbors_within_batch_used': int(neighbors_within),
        'neighbors_within_batch_requested': int(BBKNN_PARAMS['neighbors_within_batch']),
        'min_batch_size': int(min_batch_size),
    }
    
    with open(bbknn_dir / 'bbknn_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info("\n" + "="*70)
    logger.info("✓ BBKNN COMPLETE")
    logger.info("="*70)
    logger.info(f"Results:")
    logger.info(f"  scANVI UMAP: X_umap (preserved)")
    logger.info(f"  BBKNN UMAP: X_umap_bbknn (new)")
    logger.info(f"  scANVI backup: X_umap_scanvi")
    logger.info("="*70)
    
    # Cleanup
    del adata_filtered
    gc.collect()
    
    return bbknn_dir



# ============================================================================
# SUMMARY
# ============================================================================

def save_summary_stats_v1_6(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    model_genes: Optional[List[str]], 
    output_dir: Path
) -> None:
    """
    ⭐ FIXED: Safe summary with None checks
    """
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
    
    # Cell types
    summary.append(f"\nCell Type Distribution ({cols['celltype_col']}):")
    for ct, count in adata.obs[cols['celltype_col']].value_counts().items():
        pct = count / adata.n_obs * 100
        summary.append(f"  {ct}: {count:,} ({pct:.1f}%)")
    
    # ⭐ FIXED: Check if confidence_col exists
    if cols['confidence_col'] and cols['confidence_col'] in adata.obs.columns:
        summary.append(f"\nConfidence ({cols['confidence_col']}):")
        conf = adata.obs[cols['confidence_col']].describe()
        summary.append(f"  Mean: {conf['mean']:.3f}")
        summary.append(f"  Median: {conf['50%']:.3f}")
        low_conf = (adata.obs[cols['confidence_col']] < 0.5).sum()
        summary.append(f"  Low (<0.5): {low_conf:,} ({low_conf/adata.n_obs*100:.1f}%)")
    else:
        summary.append(f"\n⚠️ Confidence column not available")
    
    # Batches
    if cols['batch_col'] in adata.obs.columns:
        summary.append(f"\nBatches: {adata.obs[cols['batch_col']].nunique()}")
    
    # Disease
    if cols['disease_col'] and cols['disease_col'] in adata.obs.columns:
        summary.append(f"\nDisease Status:")
        for status, count in adata.obs[cols['disease_col']].value_counts().items():
            pct = count / adata.n_obs * 100
            summary.append(f"  {status}: {count:,} ({pct:.1f}%)")
    
    summary.append("\n" + "="*70)
    
    summary_text = '\n'.join(summary)
    with open(output_dir / 'summary_statistics.txt', 'w') as f:
        f.write(summary_text)
    
    logger.info(summary_text)


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def process_celltype(celltype_name: str, config: Dict[str, Any]) -> bool:
    """Main processing pipeline v1.6"""
    logger.info("\n" + "="*70)
    logger.info(f"PROCESSING: {celltype_name}")
    logger.info("="*70)
    
    output_dir = OUTPUT_BASE / celltype_name
    output_dir.mkdir(exist_ok=True, parents=True)
    
    try:
        # Step 1: Load
        adata_full, adata_model, lvae, cols, model_genes = load_data_and_model_v1_6(
            celltype_name, config, output_dir
        )
        
        # Summary
        save_summary_stats_v1_6(adata_full, cols, model_genes, output_dir)
        
        # Step 2: DE
        de_dir = perform_differential_expression(adata_full, adata_model, lvae, cols, output_dir)
        
        # Step 3: Denoised
        denoised_path = generate_denoised_expression_safe(adata_model, lvae, output_dir)
        
        # Step 4: Composition
        comp_dir = perform_cell_composition_analysis(adata_full, cols, output_dir)
        
        # Step 5: Pseudobulk
        pseudobulk_dir = export_pseudobulk_enhanced(adata_full, cols, output_dir)
        
        # Step 6: Markers
        viz_dir = perform_marker_visualization(adata_full, celltype_name, cols, output_dir)
        
        # Step 7: Subclustering
        subcluster_summary = perform_subclustering_safe(adata_full, cols, output_dir)
        
        # ⭐ NEW Step 8: BBKNN (optional)
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
    logger.info("scANVI DOWNSTREAM ANALYSIS - v1.6.1 PRODUCTION (HOTFIX)")
    logger.info("="*70)
    logger.info(f"\nBase: {BASE_DIR}")
    logger.info(f"Output: {OUTPUT_BASE}")
    logger.info("\nCritical fixes in v1.6.1:")
    logger.info("  ✅ Gene order preservation (P0-1)")
    logger.info("  ✅ BBKNN fixes (P0-2, P0-3)")
    logger.info("  ✅ Disease DE explicit comparison (P0-Bio-1)")
    logger.info("  ✅ Marker gene corrections (P1-1)")
    logger.info("  ✅ Bidirectional LFC filtering (P1-2)")
    logger.info("  ✅ Complete random seed (P1-3)")
    logger.info("  ✅ Package version fallback support")
    
    # Log compatibility flags
    logger.info("\nCompatibility flags:")
    compat = PACKAGE_INFO['compatibility']
    logger.info(f"  scanpy dotplot return_fig: {compat['scanpy_dotplot_return_fig']}")
    logger.info(f"  scvi seed setting: {compat['scvi_has_seed_setting']}")
    logger.info(f"  torch weights_only: {compat['torch_weights_only']}")
    
    OUTPUT_BASE.mkdir(exist_ok=True, parents=True)
    
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
