#!/usr/bin/env python3
"""
scANVI Downstream Analysis - Production Pipeline v1.5 (Marker Visualization)
=============================================================================

New in v1.5:
- ✓ Removed tissue-based subsetting
- ✓ Added comprehensive marker gene visualization (dotplot + featureplot)
- ✓ Preserved all cell type subsetting functionality

Critical features maintained from v1.4:
1. ✓ Robust model loading with proper gene set matching
2. ✓ Explicit DE loop without comparison string dependency
3. ✓ Schema validation for DE results
4. ✓ Smart expression source resolution
5. ✓ Memory-safe denoised expression (no layers write)
6. ✓ Pseudobulk export for donor-level statistics
7. ✓ HVG extraction utilities
8. ✓ Subclustering by cell type

Author: r2end
Date: 2024-12-17
Version: v1.5
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
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import gc
import json

warnings.filterwarnings('ignore')

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
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    return logging.getLogger(__name__)

# Initialize logger
logger = setup_logging()

# ============================================================================
# MARKER GENE DEFINITIONS (from R script)
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
        'ASRGL1', 'FOXJ1', 'RSPH1', 'PIFO', 'BEST4', 'C20orf85', 'C9orf24',
        'MUC5AC', 'SPDEF', 'LYPD2', 'ITLN1', 'ASCL1', 'GRP', 'POU2F3', 'ASCL2',
        'CFTR', 'FOXI1', 'ASCL3', 'BSND', 'IGF1', 'CLCNKB', 'AGER', 'RTKN2',
        'CLIC5', 'SPOCK2', 'TIMP3', 'SFTPC', 'LAMP3', 'MF5D2A', 'C8orf4',
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
        'FCER1A', 'CD1D', 'ITGAX', 'CDIC', 'FCGR2B', 'PKIB', 'CCR7', 'CD83',
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
        'DNAJB1', 'MYC', 'AFT4', 'PLAU', 'CHI3L1', 'MMP3', 'IL1R1', 'IL13RA2',
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
    ],
    
    'Immune': [
        'CLEC9A', 'XCR1', 'CADM1', 'CLNK', 'FLT3', 'ZBTB46', 'CLEC10A', 'CD1E',
        'FCER1A', 'CD1D', 'ITGAX', 'CDIC', 'FCGR2B', 'PKIB', 'CCR7', 'CD83',
        'LAMP3', 'CCL22', 'CCL17', 'CCL19', 'LAD1', 'LILRA4', 'SMPD3', 'SCT',
        'IRF7', 'PLD4', 'CLEC4C', 'MARCO', 'FABP4', 'CYP27A1', 'SIGLEC1', 'ABCG1',
        'PPARG', 'C1QA', 'C1QB', 'C1QC', 'HLA-DPA1', 'SLC40A1', 'FOLR2', 'F13A1',
        'SPP1', 'HAMP', 'VCAN', 'CCR2', 'CCR5', 'FCN1', 'S100A12', 'RNASE2',
        'LILRA5', 'MTSS1', 'TPSAB1', 'MS4A2', 'TPSB2', 'FCGR3B', 'CSF3R', 'CXCR1',
        'KLRD1', 'FCGR3A', 'GNLY', 'TYROBP', 'FCER1G', 'KLRC1', 'FGFBP2', 'SPON2',
        'MYOM2', 'TRDC', 'KRT86', 'GATA3', 'IL5', 'AREG', 'HPGDS', 'IL23R', 'RORC',
        'LST1', 'PCDH9', 'TNFSF11', 'TRDC', 'TRGC1', 'TRGC2', 'CD4', 'CD28',
        'CD40LG', 'TRAT1', 'TNFRSF25', 'CD8A', 'CD8B', 'TRGC2', 'CCR7', 'TCF7',
        'LEF1', 'SELL', 'CD28', 'IL7R', 'CCR6', 'GATA3', 'IL4', 'IL13', 'IL17A',
        'CCL20', 'CCR7', 'TCF7', 'LEF1', 'SELL', 'TBX21', 'EOMES', 'GZMK', 'KLRG1',
        'ITGA1', 'CD8A', 'CD8B', 'CCR6', 'IL7R', 'TBX21', 'GZMB', 'GZMH', 'FGFBP2',
        'ZNF683', 'IFNG', 'CCL4L2', 'PDCD1', 'KLRB1', 'IL7R', 'NCR3', 'CEBPD',
        'FOXP3', 'IL2RA', 'IKZF2', 'TNFRSF4', 'MKI67', 'TOP2A', 'TK1', 'CENPW',
        'IGHD', 'TCL1A', 'MS4A1', 'BANK1', 'MKI67', 'TOP2A', 'IGHA1', 'IGHA2',
        'IGHGP', 'IGHG1'
    ]
}

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path("/home/h2048/data/core_data")
OUTPUT_BASE = Path("/home/h2048/data/py/1217/downstream_analysis_v1_5")

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
    # Filtering thresholds
    'bayes_factor_threshold': 3.0,
    'lfc_threshold': 0.5,
    'min_proportion': 0.1
}

# Validation thresholds
VALIDATION_THRESHOLDS = {
    'min_gene_overlap_ratio': 0.9,  # Minimum gene overlap ratio for model compatibility
    'min_cells_for_de': 100,  # Minimum cells for DE analysis
    'epsilon': 1e-10  # Small epsilon for numerical stability
}

DENOISED_PARAMS = {
    'n_samples': 25,
    'library_size': 1e4,
    'max_cells_in_memory': 2e8,
    'save_separate': True,
    'save_in_layers': False  # ⭐ Never write to layers by default
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

# Marker visualization parameters
MARKER_VIZ_PARAMS = {
    'dotplot_top_genes': 30,  # Top N markers per category
    'featureplot_top_genes': 12,  # Top N for feature plots
    'expression_cutoff': 0.1,  # Minimum expression fraction for display
}

# Visualization
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
np.random.seed(42)

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


def validate_columns_strict(adata: sc.AnnData, cols: Dict[str, Optional[str]]) -> bool:
    """
    Strict validation with clear error messages
    
    FIXED: Make confidence_col optional (warning instead of error)
    """
    errors = []
    warnings = []
    
    # Required columns (must exist)
    required_cols = ['celltype_col', 'batch_col', 'sample_col']
    for col_name in required_cols:
        col = cols.get(col_name)
        if not col:
            errors.append(f"Missing config: {col_name}")
        elif col not in adata.obs.columns:
            errors.append(f"Column '{col}' ({col_name}) not found in adata.obs")
    
    # ⭐ FIXED: Optional columns (warning only)
    optional_cols = ['confidence_col', 'disease_col', 'tissue_col']
    for col_name in optional_cols:
        col = cols.get(col_name)
        if not col:
            warnings.append(f"Optional column not found: {col_name}")
        elif col not in adata.obs.columns:
            warnings.append(f"Column '{col}' ({col_name}) specified but not in adata.obs")
    
    # Log warnings
    if warnings:
        logger.warning("\n⚠️  Optional columns missing:")
        for warning in warnings:
            logger.warning(f"    {warning}")
        logger.warning("  Analysis will continue but some features may be limited")
    
    # Raise errors
    if errors:
        raise ValueError("\n".join(errors))
    
    return True


def create_model_compatible_adata(
    adata_full: sc.AnnData, 
    model_genes: List[str],
    min_overlap_ratio: float = None
) -> sc.AnnData:
    """
    Create subset with genes matching model training set
    Uses shared memory via .raw for efficiency
    
    Parameters:
    -----------
    adata_full : AnnData
        Full dataset
    model_genes : List[str]
        List of genes used in model training
    min_overlap_ratio : float, optional
        Minimum overlap ratio (default from VALIDATION_THRESHOLDS)
    
    Returns:
    --------
    AnnData
        Model-compatible subset
    """
    if min_overlap_ratio is None:
        min_overlap_ratio = VALIDATION_THRESHOLDS['min_gene_overlap_ratio']
    
    overlap_genes = list(set(model_genes) & set(adata_full.var_names))
    
    logger.info(f"  Total genes in data: {adata_full.n_vars}")
    logger.info(f"  Model training genes: {len(model_genes)}")
    overlap_ratio = len(overlap_genes) / len(model_genes) if len(model_genes) > 0 else 0.0
    logger.info(f"  Overlap: {len(overlap_genes)} ({overlap_ratio*100:.1f}%)")
    
    if overlap_ratio < min_overlap_ratio:
        warning_msg = f"  ⚠️  Warning: <{min_overlap_ratio*100:.0f}% gene overlap may affect model performance"
        logger.warning(warning_msg)
    
    adata_model = adata_full[:, overlap_genes].copy()
    
    # Preserve .raw if it exists (for full gene access later)
    if adata_full.raw is not None:
        adata_model.raw = sc.AnnData(
            X=adata_full.raw.X,
            obs=adata_full.obs.copy(),
            var=adata_full.raw.var.copy()
        )
    
    return adata_model


def validate_and_fix_de_schema(de_df: pd.DataFrame, source: str = "scANVI") -> pd.DataFrame:
    """
    Validate and fix DE result schema
    Ensures consistent column names across methods
    
    Parameters:
    -----------
    de_df : DataFrame
        Differential expression results
    source : str
        Source method ("scANVI" or "scanpy")
    
    Returns:
    --------
    DataFrame
        Validated and fixed DataFrame
    """
    if 'gene' not in de_df.columns and 'proba_de' not in de_df.columns:
        if de_df.index.name is None:
            de_df['gene'] = de_df.index
        else:
            de_df = de_df.reset_index().rename(columns={'index': 'gene'})
    
    required_cols = {'gene'}
    if source == "scANVI":
        required_cols.update({'lfc_mean', 'bayes_factor', 'proba_de'})
    
    missing = required_cols - set(de_df.columns)
    if missing:
        raise ValueError(f"DE schema validation failed. Missing: {missing}")
    
    return de_df

# ============================================================================
# UTILITY: HVG EXTRACTION FROM MODEL
# ============================================================================

def extract_var_names_from_model(
    model_path: Union[str, Path], 
    adata_full: Optional[sc.AnnData] = None, 
    model_type: str = 'SCANVI'
) -> Optional[List[str]]:
    """
    Extract gene names (training var_names) from saved model
    
    IMPROVED VERSION (v1.3 - Based on diagnostic results):
    - Method A: Direct checkpoint read (fastest, no model loading)
    - Method B: Load with adata (most reliable)
    - Method C: Legacy registry extraction (fallback)
    
    Parameters:
    -----------
    model_path : Path or str
        Path to saved model directory
    adata_full : AnnData, optional
        Full dataset for model loading fallback (recommended)
    model_type : str
        'SCVI' or 'SCANVI'
    
    Returns:
    --------
    list of gene names, or None if extraction fails
    """
    import torch
    
    try:
        # Method A: Direct checkpoint read (FAST - from diagnostic Method 4)
        logger.info("  Method A: Reading var_names from checkpoint...")
        try:
            model_pt = Path(model_path) / 'model.pt'
            if model_pt.exists():
                checkpoint = torch.load(model_pt, map_location='cpu', weights_only=False)
                
                # Check if var_names exists in checkpoint
                if 'var_names' in checkpoint:
                    var_names = checkpoint['var_names']
                    if isinstance(var_names, np.ndarray):
                        var_names = var_names.tolist()
                    elif isinstance(var_names, torch.Tensor):
                        var_names = var_names.cpu().numpy().tolist()
                    
                    if len(var_names) > 0:
                        logger.info(f"  ✓ Extracted {len(var_names)} genes from checkpoint")
                        logger.debug(f"    First 5: {var_names[:5]}")
                        logger.debug(f"    Last 5: {var_names[-5:]}")
                        return list(var_names)
        except (FileNotFoundError, KeyError, RuntimeError) as e:
            logger.warning(f"  ⚠️ Checkpoint read failed: {e}")
        except Exception as e:
            logger.error(f"  ❌ Unexpected error in Method A: {e}", exc_info=True)
        
        # Method B: Load model with adata (RELIABLE - from diagnostic Method 3/5)
        if adata_full is not None:
            logger.info("  Method B: Loading model with adata...")
            try:
                Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
                model = Model.load(model_path, adata=adata_full)
                
                var_names = list(model.adata.var_names)
                logger.info(f"  ✓ Extracted {len(var_names)} genes from model.adata")
                logger.debug(f"    First 5: {var_names[:5]}")
                logger.debug(f"    Last 5: {var_names[-5:]}")
                
                return var_names
                
            except (FileNotFoundError, ValueError, RuntimeError) as e:
                logger.warning(f"  ⚠️ Model loading failed: {e}")
            except Exception as e:
                logger.error(f"  ❌ Unexpected error in Method B: {e}", exc_info=True)
        else:
            logger.warning("  ⚠️ Method B skipped (no adata_full provided)")
        
        # Method C: Legacy registry extraction (for compatibility)
        logger.info("  Method C: Legacy registry extraction...")
        Model = scvi.model.SCANVI if model_type == 'SCANVI' else scvi.model.SCVI
        
        try:
            model = Model.load(model_path, adata=None)
            vn_dict = model.get_var_names()
            
            if isinstance(vn_dict, dict):
                var_names = list(list(vn_dict.values())[0])
                logger.info(f"  ✓ Extracted {len(var_names)} genes via legacy method")
                return var_names
        except (FileNotFoundError, ValueError, RuntimeError) as e:
            logger.warning(f"  ⚠️ Legacy method failed: {e}")
        except Exception as e:
            logger.error(f"  ❌ Unexpected error in Method C: {e}", exc_info=True)
        
        logger.error("  ❌ All extraction methods failed")
        return None
        
    except Exception as e:
        logger.error(f"  ❌ Failed to extract var_names: {e}", exc_info=True)
        return None


def save_hvg_list(genes: List[str], output_path: Union[str, Path]) -> None:
    """Save gene list to text file for reproducibility"""
    with open(output_path, 'w') as f:
        for gene in genes:
            f.write(f"{gene}\n")
    logger.info(f"  ✓ Saved gene list: {output_path}")


# ============================================================================
# UTILITY: EXPRESSION SOURCE RESOLUTION
# ============================================================================

def resolve_expr_source(adata: sc.AnnData) -> Dict[str, Union[bool, str]]:
    """
    Smart resolution of expression matrix source
    Avoids blind use_raw=True that may fail or use wrong matrix
    
    Priority:
    1. adata.raw (if exists and has data)
    2. adata.layers['log1p']
    3. adata.layers['counts'] (with warning)
    4. adata.X (fallback)
    
    Returns:
    --------
    dict
        Dict for scanpy functions: {'use_raw': bool} or {'layer': str}
    """
    # Check .raw
    if adata.raw is not None:
        if adata.raw.X is not None and adata.raw.n_vars > 0:
            logger.info(f"  Expression source: .raw ({adata.raw.n_vars} genes)")
            return {'use_raw': True}
    
    # Check layers
    if 'log1p' in adata.layers:
        logger.info("  Expression source: layers['log1p']")
        return {'layer': 'log1p'}
    
    if 'counts' in adata.layers:
        logger.warning("  ⚠️  Expression source: layers['counts'] (raw counts - may need normalization)")
        return {'layer': 'counts'}
    
    # Fallback to .X
    logger.info("  Expression source: .X (assuming log-normalized)")
    return {'use_raw': False}


def get_counts_matrix(adata: sc.AnnData) -> Tuple[Any, pd.Index]:
    """
    Smart extraction of count matrix for pseudobulk
    
    Priority:
    1. adata.layers['counts']
    2. adata.raw.X (if looks like counts)
    3. adata.X (with warning)
    
    Returns:
    --------
    tuple
        (counts_matrix, gene_names)
    """
    if 'counts' in adata.layers:
        logger.info("  Count source: layers['counts']")
        counts = adata.layers['counts']
        genes = adata.var_names
        return counts, genes
    
    if adata.raw is not None and adata.raw.X is not None:
        # Check if raw.X looks like counts (integers, no negatives)
        # Use small sample to avoid memory issues
        sample_size = min(100, adata.raw.n_obs, adata.raw.n_vars)
        if hasattr(adata.raw.X, 'toarray'):
            sample = adata.raw.X[:sample_size, :sample_size].toarray()
        else:
            sample = adata.raw.X[:sample_size, :sample_size]
        
        is_counts = np.allclose(sample, sample.astype(int)) and (sample >= 0).all()
        
        if is_counts:
            logger.info("  Count source: .raw.X (verified as counts)")
            counts = adata.raw.X
            genes = adata.raw.var_names
            return counts, genes
    
    logger.warning("  ⚠️  Count source: .X (may not be raw counts)")
    counts = adata.X
    genes = adata.var_names
    return counts, genes


# ============================================================================
# STEP 1: DATA & MODEL LOADING (ROBUST)
# ============================================================================

def infer_column_names(
    adata: sc.AnnData, 
    naming_mode: str = 'standard'
) -> Dict[str, Optional[str]]:
    """
    Infer column names from adata.obs with mode-specific strategies
    
    naming_mode:
    - 'standard': celltype_scanvi, scanvi_confidence, etc.
    - 'semantic': semantic_celltype, semantic_confidence, etc.
    """
    cols = {
        'celltype_col': None,
        'confidence_col': None,
        'batch_col': None,
        'sample_col': None,
        'disease_col': None,
        'tissue_col': None
    }
    
    # Celltype and confidence
    if naming_mode == 'semantic':
        celltype_candidates = ['semantic_celltype', 'cell_type', 'celltype']
        confidence_candidates = ['semantic_confidence','scanvi_confidence','celltypist_confidence','conf_score','confidence']
    else:
        celltype_candidates = ['celltype_scanvi', 'cell_type', 'celltype']
        confidence_candidates = ['scanvi_confidence','confidence']
    
    for col in celltype_candidates:
        if col in adata.obs.columns:
            cols['celltype_col'] = col
            break
    
    for col in confidence_candidates:
        if col in adata.obs.columns:
            cols['confidence_col'] = col
            break
    
    # Batch
    batch_candidates = ['dataset', 'batch', 'sample_id']
    for col in batch_candidates:
        if col in adata.obs.columns:
            cols['batch_col'] = col
            break
    
    # Sample (for pseudobulk)
    sample_candidates = ['sample', 'sample_id', 'donor', 'batch']
    for col in sample_candidates:
        if col in adata.obs.columns:
            cols['sample_col'] = col
            break
    
    # Disease
    disease_candidates = ['disease_status', 'condition', 'group']
    for col in disease_candidates:
        if col in adata.obs.columns:
            cols['disease_col'] = col
            break
    
    # Tissue
    tissue_candidates = ['tissue', 'tissue_type', 'anatomical_location']
    for col in tissue_candidates:
        if col in adata.obs.columns:
            cols['tissue_col'] = col
            break
    
    return cols


def load_data_and_model_production(
    celltype_name: str, 
    config: Dict[str, Any], 
    output_dir: Path
) -> Tuple[sc.AnnData, Optional[sc.AnnData], Optional[Any], Dict[str, Optional[str]], Optional[List[str]]]:
    """
    Robust data and model loading with gene set matching
    
    Returns:
    --------
    adata_full : Full dataset (all genes)
    adata_model : Model-compatible subset (HVG only)
    lvae : scANVI model (or None if loading failed)
    cols : Column name mapping
    model_genes : List of genes used in model training
    """
    log_step(1, "Loading Data and Model", "")
    
    # Load data
    h5ad_path = BASE_DIR / config['h5ad']
    logger.info(f"Loading: {h5ad_path}")
    
    try:
        adata_full = sc.read_h5ad(h5ad_path)
        logger.info(f"  Loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars} genes")
    except (FileNotFoundError, OSError) as e:
        logger.error(f"  ❌ Failed to load data: {e}")
        raise
    except Exception as e:
        logger.error(f"  ❌ Unexpected error loading data: {e}", exc_info=True)
        raise
    
    # Infer columns
    cols = infer_column_names(adata_full, config['naming_mode'])
    
    logger.info("\nColumn mapping:")
    for key, val in cols.items():
        logger.info(f"  {key}: {val}")
    
    try:
        validate_columns_strict(adata_full, cols)
    except ValueError as e:
        logger.error(f"\n❌ Column validation failed: {e}")
        return adata_full, None, None, cols, None
    
    # Try to load model
    scanvi_path = BASE_DIR / config['scanvi_model']
    logger.info(f"\nAttempting scANVI model load: {scanvi_path}")
    
    lvae = None
    adata_model = None
    model_genes = None
    
    try:
        # Step 1: Extract genes from model (FIXED: now passes adata_full)
        logger.info("\n--- Extracting training genes from model ---")
        model_genes = extract_var_names_from_model(
            scanvi_path, 
            adata_full=adata_full,  # Pass full adata for Method B fallback
            model_type='SCANVI'
        )
        
        if model_genes is None:
            raise ValueError("Could not extract gene list from model")
        
        # Save gene list for reproducibility
        gene_list_path = output_dir / f'{celltype_name}_model_genes.txt'
        save_hvg_list(model_genes, gene_list_path)
        
        # Step 2: Create model-compatible adata
        logger.info("\n--- Creating model-compatible subset ---")
        adata_model = create_model_compatible_adata(adata_full, model_genes)
        
        # Step 3: Load model with matched adata
        logger.info("\n--- Loading scANVI model ---")
        lvae = scvi.model.SCANVI.load(scanvi_path, adata=adata_model)
        
        logger.info("  ✓ Model loaded successfully")
        logger.info(f"  Model adata: {adata_model.n_obs:,} cells × {adata_model.n_vars} genes")
        
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        logger.warning(f"\n⚠️ Model loading failed: {e}")
        logger.info("  Proceeding with scanpy-only analysis")
        lvae = None
        adata_model = None
    except Exception as e:
        logger.error(f"\n❌ Unexpected error in model loading: {e}", exc_info=True)
        logger.info("  Proceeding with scanpy-only analysis")
        lvae = None
        adata_model = None
    
    return adata_full, adata_model, lvae, cols, model_genes


# ============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION (ROBUST)
# ============================================================================

def perform_de_one_vs_rest_explicit(
    adata_model: sc.AnnData, 
    lvae: Any, 
    celltype_col: str, 
    output_dir: Path
) -> Optional[pd.DataFrame]:
    """
    Explicit loop DE without comparison string dependency
    Most robust approach per code review
    """
    logger.info("\n--- One-vs-Rest DE (Explicit Loop) ---")
    
    celltypes = adata_model.obs[celltype_col].cat.categories if \
                hasattr(adata_model.obs[celltype_col], 'cat') else \
                adata_model.obs[celltype_col].unique()
    
    all_results = []
    
    for ct in celltypes:
        logger.info(f"\n  Processing: {ct}")
        
        try:
            # Explicit group1 vs None (all others)
            de = lvae.differential_expression(
                groupby=celltype_col,
                group1=ct,
                group2=None,
                mode=DE_PARAMS['mode'],
                delta=DE_PARAMS['delta'],
                batch_correction=DE_PARAMS['batch_correction'],
                n_samples=DE_PARAMS['n_samples']
            )
            
            # Validate and fix schema
            de = validate_and_fix_de_schema(de, source="scANVI")
            
            # Add tracking columns
            de['cell_type'] = ct
            de['comparison'] = f'{ct}_vs_Rest'
            
            all_results.append(de)
            
            # Filter significant
            de_sig = de[
                (de['bayes_factor'] > DE_PARAMS['bayes_factor_threshold']) &
                (de['lfc_mean'] > DE_PARAMS['lfc_threshold'])
            ]
            
            if 'non_zeros_proportion1' in de.columns:
                de_sig = de_sig[de_sig['non_zeros_proportion1'] > DE_PARAMS['min_proportion']]
            
            logger.info(f"    Significant markers: {len(de_sig)}")
            
            # Save per-celltype
            safe_name = str(ct).replace('/', '_').replace(' ', '_')
            de_sig.to_csv(output_dir / f'{safe_name}_markers.csv', index=False)
            
        except (ValueError, RuntimeError) as e:
            logger.warning(f"    ⚠️ DE failed for {ct}: {e}")
            continue
        except Exception as e:
            logger.error(f"    ❌ Unexpected error in DE for {ct}: {e}", exc_info=True)
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
    """Scanpy fallback with proper expression source"""
    logger.info("\n--- Scanpy Fallback DE ---")
    
    expr_source = resolve_expr_source(adata)
    logger.info(f"  Expression source: {expr_source}")
    
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


def perform_de_disease_within_celltype(
    adata_for_de: sc.AnnData, 
    lvae: Optional[Any], 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> None:
    """
    Disease vs Healthy DE within each cell type
    ⚠️ Note: Cell-level DE inflates p-values
         Pseudobulk recommended for final inference
    """
    if not cols['disease_col'] or cols['disease_col'] not in adata_for_de.obs.columns:
        logger.warning("\n⚠️ No disease column - skipping")
        return
    
    logger.info("\n--- Disease vs Healthy DE (Within Cell Type) ---")
    logger.warning("⚠️ NOTE: Cell-level statistics - see pseudobulk for donor-level inference")
    
    celltype_col = cols['celltype_col']
    disease_col = cols['disease_col']
    
    for ct in adata_for_de.obs[celltype_col].unique():
        logger.info(f"\n  {ct}...")
        adata_ct = adata_for_de[adata_for_de.obs[celltype_col] == ct].copy()
        
        if adata_ct.n_obs < VALIDATION_THRESHOLDS['min_cells_for_de']:
            logger.warning(f"    Skipping (n={adata_ct.n_obs})")
            continue
        
        disease_states = adata_ct.obs[disease_col].unique()
        if len(disease_states) < 2:
            logger.warning(f"    Skipping (only {disease_states})")
            continue
        
        if lvae is not None:
            try:
                de = lvae.differential_expression(
                    adata=adata_ct,
                    groupby=disease_col,
                    mode=DE_PARAMS['mode'],
                    delta=0.25,
                    batch_correction=True,
                    n_samples=DE_PARAMS['n_samples']
                )
                
                de = validate_and_fix_de_schema(de, source="scANVI")
                
                safe_name = str(ct).replace('/', '_').replace(' ', '_')
                de.to_csv(output_dir / f'{safe_name}_disease_vs_healthy_scANVI.csv', index=False)
                
                n_sig = len(de[de['bayes_factor'] > 3])
                logger.info(f"    ✓ {n_sig} significant (BF>3)")
                
            except (ValueError, RuntimeError) as e:
                logger.warning(f"    ⚠️ Failed: {e}")
            except Exception as e:
                logger.error(f"    ❌ Unexpected error: {e}", exc_info=True)
        else:
            expr_source = resolve_expr_source(adata_ct)
            sc.tl.rank_genes_groups(adata_ct, groupby=disease_col, method='wilcoxon', **expr_source)
            de = sc.get.rank_genes_groups_df(adata_ct, group=None)
            safe_name = str(ct).replace('/', '_').replace(' ', '_')
            de.to_csv(output_dir / f'{safe_name}_disease_vs_healthy_scanpy.csv', index=False)
            logger.info(f"    ✓ Scanpy DE done")


def perform_differential_expression(
    adata_full: sc.AnnData, 
    adata_model: Optional[sc.AnnData], 
    lvae: Optional[Any], 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """
    Differential Expression Analysis - PSEUDOBULK ONLY
    
    Cell-level DE is DISABLED due to pseudo-replication issues.
    For tissue comparison, use donor-level pseudobulk (Step 5).
    """
    log_step(2, "Differential Expression Analysis", "")
    
    de_dir = output_dir / 'differential_expression'
    de_dir.mkdir(exist_ok=True, parents=True)
    
    logger.warning("\n" + "="*70)
    logger.warning("WARNING: CELL-LEVEL DE DISABLED (By Design)")
    logger.warning("="*70)
    logger.warning("\nReason:")
    logger.warning("  - Cell-level statistics have pseudo-replication issues")
    logger.warning("  - P-values can be inflated by orders of magnitude")
    logger.warning("  - Not suitable for tissue comparison")
    logger.warning("  - Single-sample groups cause crashes")
    
    logger.info("\nRecommended Workflow for Tissue Comparison:")
    logger.info("  1. This pipeline generates pseudobulk counts (Step 5)")
    logger.info("  2. Use DESeq2/edgeR in R for donor-level statistics")
    logger.info("  3. Properly controls for batch and donor effects")
    
    logger.info("\nDocumentation:")
    logger.info(f"  README: pseudobulk/README_TISSUE_DE.md")
    logger.info(f"  Template: pseudobulk/tissue_DE_DESeq2.R")
    
    # Create instruction file
    readme_content = """# Tissue Comparison via Pseudobulk DE

## Why Skip Cell-Level DE?

Problems with cell-level statistics:
- Cells from same donor are NOT independent
- P-values severely underestimated
- Cannot properly control for donor effects
- Single-sample groups cause crashes

Solution: Donor-level pseudobulk with DESeq2/edgeR

## Workflow

Step 1: Run this pipeline to generate pseudobulk_counts.h5ad

Step 2: Use tissue_DE_DESeq2.R template for tissue comparison

Step 3: Interpret with FDR < 0.05 and |log2FC| > 1

## Reference

Squair et al. (2021). Confronting false discoveries in single-cell DE.
Nat Commun 12, 5692.
"""
    
    readme_path = de_dir / 'README_TISSUE_DE.md'
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    
    logger.info(f"\nCreated: {readme_path}")
    logger.info("="*70)
    
    return de_dir


# ============================================================================
# STEP 3: DENOISED EXPRESSION (MEMORY-SAFE)
# ============================================================================

def generate_denoised_expression_safe(
    adata_model: Optional[sc.AnnData], 
    lvae: Optional[Any], 
    output_dir: Path
) -> Optional[Path]:
    """Memory-safe denoised expression - NEVER writes to layers"""
    log_step(3, "Generating Denoised Expression", "")
    
    if lvae is None or adata_model is None:
        logger.warning("⚠️ Model not available - skipping")
        return None
    
    n_elements = adata_model.n_obs * adata_model.n_vars
    if n_elements > DENOISED_PARAMS['max_cells_in_memory']:
        logger.warning(f"⚠️ Dataset too large ({n_elements:.2e} > {DENOISED_PARAMS['max_cells_in_memory']:.2e})")
        logger.warning("  Skipping to prevent memory issues")
        return None
    
    logger.info(f"\nGenerating denoised expression...")
    logger.info(f"  n_samples: {DENOISED_PARAMS['n_samples']}")
    logger.info(f"  library_size: {DENOISED_PARAMS['library_size']}")
    
    try:
        denoised = lvae.get_normalized_expression(
            adata=adata_model,
            n_samples=DENOISED_PARAMS['n_samples'],
            return_mean=True,
            library_size=DENOISED_PARAMS['library_size']
        )
        
        denoised = denoised.astype(np.float32)
        
        logger.info(f"  ✓ Generated: {denoised.shape}")
        logger.info(f"  Memory: {denoised.nbytes / 1e9:.2f} GB")
        
        # Save as separate file ONLY
        logger.info("\n  Saving as separate h5ad...")
        adata_denoised = sc.AnnData(
            X=denoised,
            obs=adata_model.obs.copy(),
            var=adata_model.var.copy()
        )
        
        denoised_path = output_dir / 'denoised_expression.h5ad'
        adata_denoised.write_h5ad(denoised_path, compression='gzip')
        logger.info(f"  ✓ Saved: {denoised_path}")
        
        del adata_denoised
        gc.collect()
        
        return denoised_path
        
    except (ValueError, RuntimeError, MemoryError) as e:
        logger.error(f"\n  ⚠️ Failed: {e}", exc_info=True)
        return None
    except Exception as e:
        logger.error(f"\n  ❌ Unexpected error: {e}", exc_info=True)
        return None


# ============================================================================
# STEP 4: CELL COMPOSITION ANALYSIS
# ============================================================================

def perform_cell_composition_analysis(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """Cell composition without tissue-specific subsetting"""
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
        
        # Statistical testing
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
                        'fold_change': group1.mean() / (group2.mean() + VALIDATION_THRESHOLDS['epsilon']),
                        'p_value': p_val
                    })
            
            if comp_stats:
                comp_stats_df = pd.DataFrame(comp_stats)
                comp_stats_df['p_adj'] = multipletests(comp_stats_df['p_value'], method='fdr_bh')[1]
                comp_stats_df = comp_stats_df.sort_values('p_adj')
                comp_stats_df.to_csv(comp_dir / 'composition_disease_statistics.csv', index=False)
                
                sig = comp_stats_df[comp_stats_df['p_adj'] < 0.05]
                if len(sig) > 0:
                    logger.info(f"\n  Significant changes (FDR < 0.05): {len(sig)}")
    
    logger.info("\n✓ Cell composition complete")
    return comp_dir


# ============================================================================
# STEP 5: PSEUDOBULK EXPORT (ENHANCED)
# ============================================================================

def export_pseudobulk_enhanced(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Optional[Path]:
    """
    Enhanced pseudobulk export with metadata preservation
    Critical for donor-level statistical inference
    """
    log_step(5, "Pseudobulk Export (Donor-Level)", "")
    
    pseudobulk_dir = output_dir / 'pseudobulk'
    pseudobulk_dir.mkdir(exist_ok=True, parents=True)
    
    celltype_col = cols['celltype_col']
    sample_col = cols['sample_col']
    
    # Get counts
    counts, genes = get_counts_matrix(adata)
    
    if counts is None:
        logger.warning("⚠️ No counts matrix - skipping")
        return None
    
    logger.info(f"\nAggregating by {sample_col} × {celltype_col}...")
    
    # Create aggregation key
    adata.obs['pseudobulk_key'] = (
        adata.obs[sample_col].astype(str) + '|' +
        adata.obs[celltype_col].astype(str)
    )
    
    pseudobulk_list = []
    metadata_list = []
    
    for key in adata.obs['pseudobulk_key'].unique():
        mask = adata.obs['pseudobulk_key'] == key
        n_cells = mask.sum()
        
        if n_cells < PSEUDOBULK_PARAMS['min_cells_per_sample']:
            continue
        
        # ⭐ FIX: Convert pandas Series to numpy array for sparse matrix indexing
        # Use np.asarray for compatibility with both pandas Series and numpy arrays
        mask_array = np.asarray(mask, dtype=bool)
        
        # Aggregate - memory-efficient for sparse matrices
        try:
            if PSEUDOBULK_PARAMS['aggregation'] == 'sum':
                aggregated = counts[mask_array, :].sum(axis=0)
            else:  # mean
                aggregated = counts[mask_array, :].mean(axis=0)
            
            # Convert to 1D numpy array efficiently
            if hasattr(aggregated, 'A1'):  # scipy sparse matrix result
                pseudo_counts = aggregated.A1
            elif hasattr(aggregated, 'A'):  # scipy sparse matrix (older)
                pseudo_counts = np.asarray(aggregated).flatten()
            elif hasattr(aggregated, 'flatten'):  # numpy array
                pseudo_counts = aggregated.flatten()
            else:  # fallback
                pseudo_counts = np.asarray(aggregated).flatten()
            
            # Validate shape consistency
            if len(pseudo_counts) != len(genes):
                logger.warning(f"  ⚠️ Shape mismatch for {key}: expected {len(genes)}, got {len(pseudo_counts)}")
                continue
                
        except (IndexError, ValueError, MemoryError) as e:
            logger.error(f"  ❌ Failed to aggregate for key {key}: {e}")
            continue
        except Exception as e:
            logger.error(f"  ❌ Unexpected error aggregating {key}: {e}", exc_info=True)
            continue
        
        pseudobulk_list.append(pseudo_counts)
        
        # Metadata
        meta = adata.obs[mask].iloc[0]
        metadata_list.append({
            'pseudobulk_id': key,
            'sample_id': meta[sample_col],
            'cell_type': meta[celltype_col],
            'n_cells': n_cells,
            'batch': meta[cols['batch_col']] if cols['batch_col'] in meta else None,
            'disease_status': meta[cols['disease_col']] if cols['disease_col'] and cols['disease_col'] in meta else None,
            'tissue': meta[cols['tissue_col']] if cols['tissue_col'] and cols['tissue_col'] in meta else None
        })
    
    if len(pseudobulk_list) == 0:
        logger.warning("⚠️ No pseudobulk samples passed filters")
        return None
    
    # Create pseudobulk AnnData
    # Validate all arrays have same shape before stacking
    if len(pseudobulk_list) > 0:
        expected_shape = len(pseudobulk_list[0])
        for i, arr in enumerate(pseudobulk_list):
            if len(arr) != expected_shape:
                logger.error(f"  ❌ Shape mismatch: sample {i} has {len(arr)} genes, expected {expected_shape}")
                raise ValueError(f"Inconsistent gene counts in pseudobulk samples")
    
    # Determine appropriate dtype before stacking (memory-efficient)
    dtype = np.int32
    if len(pseudobulk_list) > 0:
        # Check first few samples to determine if int64 is needed
        sample_size = min(10, len(pseudobulk_list))
        for arr in pseudobulk_list[:sample_size]:
            if arr.size > 0:
                arr_max = arr.max()
                arr_min = arr.min()
                if arr_max > np.iinfo(np.int32).max or arr_min < np.iinfo(np.int32).min:
                    logger.warning(f"  ⚠️ Values exceed int32 range, using int64")
                    dtype = np.int64
                    break
    
    pseudobulk_matrix = np.vstack(pseudobulk_list)
    pseudobulk_metadata = pd.DataFrame(metadata_list)
    
    adata_pseudo = sc.AnnData(
        X=pseudobulk_matrix.astype(dtype),
        obs=pseudobulk_metadata,
        var=pd.DataFrame(index=genes)
    )
    
    logger.info(f"\n  Generated {adata_pseudo.n_obs} pseudobulk samples")
    logger.info(f"  Genes: {adata_pseudo.n_vars}")
    
    # Save
    pseudo_path = pseudobulk_dir / 'pseudobulk_counts.h5ad'
    adata_pseudo.write_h5ad(pseudo_path, compression='gzip')
    logger.info(f"  ✓ Saved: {pseudo_path}")
    
    # Cell count summary
    ct_cell_counts = pseudobulk_metadata.groupby('cell_type')['n_cells'].agg(['sum', 'mean', 'median'])
    ct_cell_counts.to_csv(pseudobulk_dir / 'cell_counts_summary.csv')
    logger.info(f"  ✓ Saved: cell_counts_summary.csv")
    
    return pseudobulk_dir


# ============================================================================
# STEP 6: MARKER VISUALIZATION (NEW IN v1.5)
# ============================================================================

def plot_marker_dotplot(
    adata: sc.AnnData, 
    markers: List[str], 
    celltype_col: str, 
    output_path: Path, 
    title: str = "Marker Genes"
) -> None:
    """
    Create dotplot for marker genes
    
    Parameters:
    -----------
    adata : AnnData
        Dataset with expression data
    markers : list
        List of marker genes
    celltype_col : str
        Column name for cell types
    output_path : Path
        Output file path
    title : str
        Plot title
    """
    # Filter to genes present in dataset
    available_markers = [g for g in markers if g in adata.var_names]
    
    if len(available_markers) == 0:
        logger.warning(f"  ⚠️ No markers found in dataset for {title}")
        return
    
    logger.info(f"  Creating dotplot: {len(available_markers)}/{len(markers)} markers available")
    
    # Resolve expression source
    expr_source = resolve_expr_source(adata)
    use_raw = expr_source.get('use_raw', False)
    layer = expr_source.get('layer', None)
    
    try:
        # Create figure
        fig, ax = plt.subplots(figsize=(max(12, len(available_markers) * 0.3), 
                                        max(8, adata.obs[celltype_col].nunique() * 0.5)))
        
        sc.pl.dotplot(
            adata,
            var_names=available_markers,
            groupby=celltype_col,
            dendrogram=False,
            use_raw=use_raw,
            layer=layer,
            standard_scale='var',  # Z-score normalization
            cmap='RdYlBu_r',
            ax=ax,
            show=False
        )
        
        ax.set_title(title, fontsize=14, weight='bold')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"    ✓ Saved: {output_path.name}")
        
    except (ValueError, RuntimeError) as e:
        logger.warning(f"    ⚠️ Dotplot failed: {e}")
    except Exception as e:
        logger.error(f"    ❌ Unexpected error in dotplot: {e}", exc_info=True)


def plot_marker_featureplot(
    adata: sc.AnnData, 
    markers: List[str], 
    output_dir: Path, 
    prefix: str = "markers"
) -> None:
    """
    Create feature plots (UMAP) for top marker genes
    
    Parameters:
    -----------
    adata : AnnData
        Dataset with UMAP coordinates
    markers : list
        List of marker genes
    output_dir : Path
        Output directory
    prefix : str
        Filename prefix
    """
    # Filter to genes present in dataset
    available_markers = [g for g in markers if g in adata.var_names]
    
    if len(available_markers) == 0:
        logger.warning(f"  ⚠️ No markers found for feature plot")
        return
    
    # Limit to top N genes for visualization
    n_genes = min(len(available_markers), MARKER_VIZ_PARAMS['featureplot_top_genes'])
    markers_to_plot = available_markers[:n_genes]
    
    logger.info(f"  Creating feature plots: {n_genes} markers")
    
    # Check if UMAP exists
    if 'X_umap' not in adata.obsm:
        logger.warning("    ⚠️ No UMAP coordinates - skipping feature plots")
        return
    
    # Resolve expression source
    expr_source = resolve_expr_source(adata)
    use_raw = expr_source.get('use_raw', False)
    layer = expr_source.get('layer', None)
    
    try:
        # Create grid of plots
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
        
        # Hide unused subplots
        for i in range(n_genes, len(axes)):
            axes[i].axis('off')
        
        plt.tight_layout()
        output_path = output_dir / f'{prefix}_featureplot.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"    ✓ Saved: {output_path.name}")
        
    except (ValueError, RuntimeError) as e:
        logger.warning(f"    ⚠️ Feature plot failed: {e}")
    except Exception as e:
        logger.error(f"    ❌ Unexpected error in feature plot: {e}", exc_info=True)


def perform_marker_visualization(
    adata: sc.AnnData, 
    celltype_name: str, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Path:
    """
    Comprehensive marker visualization using predefined gene lists
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    celltype_name : str
        Cell type name (for selecting appropriate markers)
    cols : dict
        Column name mapping
    output_dir : Path
        Output directory
    """
    log_step(6, "Marker Gene Visualization", "")
    
    viz_dir = output_dir / 'marker_visualization'
    viz_dir.mkdir(exist_ok=True, parents=True)
    
    celltype_col = cols['celltype_col']
    
    # Select appropriate marker sets based on cell type
    marker_sets = {
        'General': MARKER_GENES['General']
    }
    
    # Add cell-type-specific markers
    if celltype_name == 'T_cells':
        marker_sets['T_cell_specific'] = MARKER_GENES['T_cells']
    elif celltype_name == 'B_cells':
        marker_sets['B_cell_specific'] = MARKER_GENES['B_cells']
    elif celltype_name == 'Myeloid':
        marker_sets['Myeloid_specific'] = MARKER_GENES['Myeloid']
    elif celltype_name == 'Stromal_Vascular':
        marker_sets['Stromal_specific'] = MARKER_GENES['Stromal_Vascular']
        marker_sets['Fibroblast'] = MARKER_GENES['Fibroblast']
        marker_sets['SMC'] = MARKER_GENES['SMC']
        marker_sets['Endothelial'] = MARKER_GENES['Endothelial']
    
    # Always include relevant general categories
    if any(x in celltype_name.lower() for x in ['epithelial', 'stromal']):
        marker_sets['Epithelial'] = MARKER_GENES['Epithelial']
    
    logger.info(f"\nGenerating visualizations for {len(marker_sets)} marker sets:")
    for name in marker_sets:
        logger.info(f"  - {name} ({len(marker_sets[name])} genes)")
    
    # Create dotplots
    logger.info("\n--- Creating Dotplots ---")
    for set_name, markers in marker_sets.items():
        # Limit markers for cleaner visualization
        top_markers = markers[:MARKER_VIZ_PARAMS['dotplot_top_genes']]
        
        output_path = viz_dir / f'dotplot_{set_name.lower()}.png'
        plot_marker_dotplot(
            adata, 
            top_markers, 
            celltype_col, 
            output_path, 
            title=f"{set_name} Markers"
        )
    
    # Create feature plots
    logger.info("\n--- Creating Feature Plots ---")
    for set_name, markers in marker_sets.items():
        plot_marker_featureplot(
            adata, 
            markers, 
            viz_dir, 
            prefix=f'{set_name.lower()}'
        )
    
    # Create comprehensive summary plot
    logger.info("\n--- Creating Summary Overview ---")
    try:
        # Select top markers from each set
        all_top_markers = []
        for markers in marker_sets.values():
            available = [g for g in markers if g in adata.var_names]
            all_top_markers.extend(available[:5])  # Top 5 from each set
        
        # Remove duplicates while preserving order
        seen = set()
        all_top_markers = [x for x in all_top_markers if not (x in seen or seen.add(x))]
        
        if len(all_top_markers) > 0:
            summary_path = viz_dir / 'summary_all_markers.png'
            plot_marker_dotplot(
                adata, 
                all_top_markers, 
                celltype_col, 
                summary_path, 
                title=f"{celltype_name} - Top Markers Summary"
            )
    except (ValueError, RuntimeError) as e:
        logger.warning(f"  ⚠️ Summary plot failed: {e}")
    except Exception as e:
        logger.error(f"  ❌ Unexpected error in summary plot: {e}", exc_info=True)
    
    logger.info(f"\n✓ Marker visualization complete")
    logger.info(f"  Output: {viz_dir}")
    
    return viz_dir


# ============================================================================
# STEP 7: SUBCLUSTERING (PRESERVED FROM v1.4)
# ============================================================================

def perform_subclustering_safe(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    output_dir: Path
) -> Dict[str, Dict[str, Any]]:
    """
    Safe subclustering by cell type (NO tissue subsetting)
    Preserved functionality from v1.4
    """
    log_step(7, "Subclustering by Cell Type", "")
    
    subcluster_dir = output_dir / 'subclustering'
    subcluster_dir.mkdir(exist_ok=True, parents=True)
    
    celltype_col = cols['celltype_col']
    
    celltypes = adata.obs[celltype_col].unique()
    logger.info(f"\nCell types to subcluster: {len(celltypes)}")
    
    results_summary = {}
    
    for ct in celltypes:
        ct_safe = str(ct).replace('/', '_').replace(' ', '_')
        logger.info(f"\n--- {ct} ---")
        
        # Filter
        adata_sub = adata[adata.obs[celltype_col] == ct].copy()
        
        if adata_sub.n_obs < SUBCLUSTER_PARAMS['min_cells_for_subcluster']:
            logger.warning(f"  Skipping (n={adata_sub.n_obs} < {SUBCLUSTER_PARAMS['min_cells_for_subcluster']})")
            continue
        
        # Create directory
        ct_dir = subcluster_dir / ct_safe
        ct_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info(f"  Cells: {adata_sub.n_obs:,}")
        
        # Neighbors on scANVI latent
        n_neighbors = min(SUBCLUSTER_PARAMS['n_neighbors'], adata_sub.n_obs // 10)
        sc.pp.neighbors(adata_sub, use_rep='X_scanvi', n_neighbors=n_neighbors)
        
        # Multi-resolution clustering
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
        
        # Figures
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
        
        # Store summary only
        results_summary[ct] = {
            'path': str(output_path),
            'n_cells': adata_sub.n_obs,
            'n_subclusters': adata_sub.obs[default_key].nunique()
        }
        
        logger.info(f"  ✓ Complete")
        
        del adata_sub
        gc.collect()
    
    # Save summary
    with open(subcluster_dir / 'subclustering_summary.json', 'w') as f:
        json.dump(results_summary, f, indent=2)
    
    logger.info(f"\n{'='*70}")
    logger.info(f"✓ Subclustering complete: {len(results_summary)} cell types")
    logger.info(f"{'='*70}")
    
    return results_summary


# ============================================================================
# SUMMARY GENERATION
# ============================================================================

def save_summary_stats(
    adata: sc.AnnData, 
    cols: Dict[str, Optional[str]], 
    model_genes: Optional[List[str]], 
    output_dir: Path
) -> None:
    """Enhanced summary with model info"""
    summary = []
    summary.append("="*70)
    summary.append("Dataset Summary Statistics")
    summary.append("="*70)
    summary.append(f"\nTotal cells: {adata.n_obs:,}")
    summary.append(f"Total genes: {adata.n_vars:,}")
    
    if model_genes:
        summary.append(f"Model training genes: {len(model_genes)}")
        overlap = len(set(model_genes) & set(adata.var_names))
        summary.append(f"  Overlap with data: {overlap} ({overlap/len(model_genes)*100:.1f}%)")
    
    # Cell types
    summary.append(f"\nCell Type Distribution ({cols['celltype_col']}):")
    for ct, count in adata.obs[cols['celltype_col']].value_counts().items():
        pct = count / adata.n_obs * 100
        summary.append(f"  {ct}: {count:,} ({pct:.1f}%)")
    
    # Confidence
    if cols['confidence_col'] and cols['confidence_col'] in adata.obs.columns:
        summary.append(f"\nscANVI Confidence:")
        conf = adata.obs[cols['confidence_col']].describe()
        summary.append(f"  Mean: {conf['mean']:.3f}")
        summary.append(f"  Median: {conf['50%']:.3f}")
        low_conf = (adata.obs[cols['confidence_col']] < 0.5).sum()
        summary.append(f"  Low (<0.5): {low_conf:,} ({low_conf/adata.n_obs*100:.1f}%)")
    else:
        summary.append(f"\n⚠️ Confidence column not found")
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
    """Main processing pipeline with marker visualization"""
    logger.info("\n" + "="*70)
    logger.info(f"PROCESSING: {celltype_name}")
    logger.info("="*70)
    
    output_dir = OUTPUT_BASE / celltype_name
    output_dir.mkdir(exist_ok=True, parents=True)
    
    adata_full = None
    adata_model = None
    lvae = None
    
    try:
        # Step 1: Load
        adata_full, adata_model, lvae, cols, model_genes = load_data_and_model_production(
            celltype_name, config, output_dir
        )
        
        # Summary
        save_summary_stats(adata_full, cols, model_genes, output_dir)
        
        # Step 2: DE
        de_dir = perform_differential_expression(adata_full, adata_model, lvae, cols, output_dir)
        
        # Step 3: Denoised
        denoised_path = generate_denoised_expression_safe(adata_model, lvae, output_dir)
        
        # Step 4: Composition
        comp_dir = perform_cell_composition_analysis(adata_full, cols, output_dir)
        
        # Step 5: Pseudobulk
        pseudobulk_dir = export_pseudobulk_enhanced(adata_full, cols, output_dir)
        
        # Step 6: Marker Visualization (NEW)
        viz_dir = perform_marker_visualization(adata_full, celltype_name, cols, output_dir)
        
        # Step 7: Subclustering
        subcluster_summary = perform_subclustering_safe(adata_full, cols, output_dir)
        
        # Save final
        logger.info("\nSaving final data...")
        final_path = output_dir / f'{celltype_name}_processed_final.h5ad'
        adata_full.write_h5ad(final_path, compression='gzip')
        logger.info(f"  ✓ Saved: {final_path}")
        
        logger.info("\n" + "="*70)
        logger.info(f"✓ {celltype_name} COMPLETE")
        logger.info("="*70)
        
        return True
        
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        logger.error(f"\n❌ ERROR: {e}", exc_info=True)
        return False
    except Exception as e:
        logger.error(f"\n❌ UNEXPECTED ERROR: {e}", exc_info=True)
        return False
    finally:
        # Ensure cleanup even if errors occur
        if adata_full is not None:
            del adata_full
        if adata_model is not None:
            del adata_model
        if lvae is not None:
            del lvae
        gc.collect()


def main() -> None:
    """Main entry point"""
    logger.info("="*70)
    logger.info("scANVI DOWNSTREAM ANALYSIS - v1.5 (Marker Visualization)")
    logger.info("="*70)
    logger.info(f"\nBase: {BASE_DIR}")
    logger.info(f"Output: {OUTPUT_BASE}")
    logger.info("\nNew in v1.5:")
    logger.info("  ✓ Comprehensive marker gene visualization")
    logger.info("  ✓ No tissue-based subsetting")
    logger.info("  ✓ Preserved cell type subclustering")
    
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
