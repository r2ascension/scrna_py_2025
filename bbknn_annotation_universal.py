#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bbknn_annotation_universal.py - FIXED VERSION

Universal Cell Type Annotation Pipeline for BBKNN-integrated Data

✅ FIXES APPLIED (Priority 1 + Priority 2):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Priority 1 (Critical - FIXED):
  ✓ Cell type parameter validation added
  ✓ Standardized cluster ID handling across all operations
  ✓ UTF-8 encoding for all CSV read/write operations
  ✓ Clear error messages for missing input files

Priority 2 (Important - FIXED):
  ✓ DE analysis reports excluded small clusters
  ✓ Memory-efficient balanced subsampling using boolean masks
  ✓ Empty heatmap batch detection and skipping
  ✓ Improved color palette for >20 cell types
  ✓ Progress bars with tqdm for long operations
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Features:
1) Support for multiple cell types with specific marker sets
2) Flexible resolution selection from clustering results
3) Three-stage workflow: pre_annotation → annotation → post_annotation
4) Automatic marker checking and visualization
5) Standardized output structure
6) Robust fallbacks for DE matrix, cluster keys, and Scanpy version compatibility
7) Balanced subsampling for very large datasets (plotting only)
8) CLI overrides for all key parameters

Author: Clinical-Bioinformatics Team
Version: v1.2 (Fixed)
Date: 2025-11-11
"""

import sys
import os
from pathlib import Path
import warnings
import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import time
from tqdm import tqdm

warnings.filterwarnings('ignore')

# Try to import packaging for version checks
try:
    from packaging import version as pkg_version
    HAS_PACKAGING = True
except ImportError:
    HAS_PACKAGING = False
    warnings.warn("⚠️  'packaging' module not found. Version checks disabled.\n"
                  "   Install with: pip install packaging")

# ==================== USER CONFIGURATIONS (can be overridden by CLI) ====================

# Select cell type to annotate
CELL_TYPE = "All"  # Options: "All", "Epithelial", "T", "Myeloid", "B", "Fibroblast", "Endothelial", "SMC"

# Select clustering resolution to use for annotation
# "default" will use "leiden_bbknn"; numbers will use "leiden_bbknn_res{RESOLUTION}"
RESOLUTION = "leiden_bbknn_res1.8"

# Run mode selection:
# - "pre_annotation":  Generate marker analysis and annotation template
# - "post_annotation": Apply annotation and generate cell type visualizations
# - "full":            Complete workflow (both stages; requires Annotation.csv if proceeding to post-)
RUN_MODE = "pre_annotation"

# Input / Output defaults (can be overridden by CLI)
INPUT_H5AD_PATH = '/home/h2048/data/py/1125/bbknn_output_optimized/adata_bbknn_integrated_with_subclusters.h5ad'  # will be resolved by _resolve_input_path if None
OUTPUT_DIR = '/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/Annotation'       # will be derived from CELL_TYPE + RESOLUTION if None

# Annotation file (created in pre_annotation, used in post_annotation)
ANNOTATION_FILE = "/home/h2048/data/py/1029/bbknn_annotation/Epithelial/res_leiden_bbknn_res2.4/Annotation.csv"  # Format: Cluster,CellType (no header)

# Marker/DE settings
RUN_FIND_MARKERS = True
MARKER_MIN_PCT = 0.25
MARKER_LOGFC_THRESHOLD = 0.25
TOP_N_MARKERS = 10        # top markers per cluster in heatmap
MIN_CELLS_PER_CLUSTER = 0  # skip DE for tiny clusters
MAX_CELLS_PER_GROUP_PLOTTING = 0  # balanced subsampling per cluster for heavy plotting

# Visualization settings
GENERATE_DOTPLOT = True
GENERATE_HEATMAP = True
HEATMAP_BATCH_SIZE = 25  # genes per heatmap batch
DPI = 300
FIGURE_FORMAT = "pdf"
VERBOSE = True
USE_PROGRESS_BARS = True  # ✅ NEW: Enable progress bars
np.random.seed(0)

# Cell type colors (will auto-generate for unknown labels)
CELL_TYPE_COLORS = {
    "Epithelial": "#E41A1C",
    "Fibroblast": "#377EB8",
    "T": "#4DAF4A",
    "Myeloid": "#984EA3",
    "Endothelial": "#FF7F00",
    "SMC": "#FFFF33",
    "B": "#A65628",
    "Plasma": "#F781BF",
    "NK": "#66C2A5",
    "Proliferation": "#FC8D62",
    "Unknown": "#CCCCCC"
}

# ==================== MARKER GENE CONFIGURATIONS ====================

MARKER_CONFIGS = {
    "All": {
        "description": "Major cell type annotation (Epithelial, Immune, Stromal)",
        "markers": [
            # Epithelial
            'EPCAM', 'KRT8', 'KRT18', 'KRT19',
            'SCGB1A1', 'MUC5B', 'FOXJ1', 'KRT5',
            # Pan-immune
            'PTPRC', 'CD53', 'CORO1A',
            # T
            'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A',
            # B
            'CD19', 'MS4A1', 'CD79A', 'CD79B',
            # Plasma
            'MZB1', 'SDC1', 'JCHAIN',
            # Myeloid
            'CD68', 'CD14', 'CD163', 'CD1C', 'CLEC9A',
            # NK
            'NCAM1', 'FCGR3A', 'GNLY',
            # Fibroblast
            'COL1A1', 'COL1A2', 'DCN', 'LUM',
            # SMC
            'ACTA2', 'TAGLN', 'MYH11',
            # Endothelial
            'PECAM1', 'VWF', 'CDH5',
            # Proliferation
            'MKI67', 'TOP2A'
        ]
    },
    "Epithelial": {
        "description": "Epithelial cell subtype annotation",
        "markers": [
            'TP63','KRT5',
            'SCGB1A1','SERPINB3','SCGB3A2','SCGB3A1','TCN1','ASRGL1','BPIFB1',
            'FOXJ1','RSPH1','PIFO','BEST4','C20orf85','C9orf24','CYP2F1',
            'MUC5AC','SPDEF','LYPD2','ITLN1',
            'ASCL1','GRP',
            'POU2F3','ASCL2',
            'CFTR','FOXI1','ASCL3','BSND','IGF1','CLCNKB','PDE1C',
            'AGER','RTKN2','CLIC5','SPOCK2','TIMP3',
            'SFTPC','LAMP3','MF5D2A','C8orf4','C11orf96',
            'VIM','SOX9',
            'KRT14','MYH11','ACTA2','MYLK',
            'DMBT1','RNASE1',
            'MUC5B','SPDEF',
            'LYZ','LTF','PIP','CCL28',
            'SFTPB','SCGB3A2','SFTA2'
        ]
    },
    "T": {
        "description": "T cell subtype annotation",
        "markers": [
            # Pan T
            'CD3D', 'CD3E', 'CD3G',
            # CD4
            'CD4', 'IL7R',
            # CD8
            'CD8A', 'CD8B',
            # Naive
            'CCR7', 'SELL', 'LEF1', 'TCF7',
            # Effector/Memory
            'GZMA', 'GZMB', 'GZMK', 'PRF1',
            # Treg
            'FOXP3', 'IL2RA', 'CTLA4',
            # Activation
            'IFNG', 'TNF', 'IL2',
            # Exhaustion
            'PDCD1', 'LAG3', 'TIGIT', 'HAVCR2',
            # Resident memory
            'CD69', 'ITGAE', 'CXCR6',
            # Proliferation
            'MKI67', 'TOP2A'
        ]
    },
    "Myeloid": {
        "description": "Myeloid cell subtype annotation",
        "markers": [
            'CLEC9A','XCR1','CADM1','CLNK','FLT3','ZBTB46',
            'CLEC10A','CD1E','FCER1A','CD1D','ITGAX','CDIC','FCGR2B','PKIB',
            'CCR7','CD83','LAMP3','CCL22','CCL17','CCL19','LAD1',
            'LILRA4','SMPD3','SCT','IRF7','PLD4','CLEC4C',
            'MARCO','FABP4','CYP27A1','SIGLEC1','ABCG1','PPARG',
            'C1QA','C1QB','C1QC','HLA-DPA1','SLC40A1',
            'FOLR2','F13A1',
            'SPP1','HAMP','VCAN','CCR2','CCR5',
            'FCN1',
            'S100A12','RNASE2',
            'LILRA5','MTSS1',
            'TPSAB1','MS4A2','TPSB2',
            'FCGR3B','CSF3R','CXCR1'
        ]
    },
    "B": {
        "description": "B cell subtype annotation",
        "markers": [
            # Pan B
            'CD19', 'MS4A1', 'CD79A', 'CD79B',
            # Naive B
            'IGHD', 'TCL1A', 'FCER2', 'IL4R',
            # Memory B
            'CD27', 'TNFRSF13B', 'AIM2',
            # GC B
            'BCL6', 'AICDA', 'MEF2B', 'RGS13',
            # Plasma
            'MZB1', 'SDC1', 'PRDM1', 'XBP1', 'JCHAIN', 'IGHA1', 'IGHG1',
            # Activation
            'CD86', 'CD80', 'CD40',
            # Proliferation
            'MKI67', 'TOP2A'
        ]
    },
    "Fibroblast": {
        "description": "Fibroblast subtype annotation",
        "markers": [
            'APOD','FGF7','COL15A1','MFAP5','PI16','CD34',
            'MMP11','COL10A1','POSTN','LRRC15','HOPX','IGFBP5','TIMP1','MMP1','COL7A1','WNT5A','ISG15','IL7R','SFRP4','SFRP2','COMP','RGS5','PDGFRB','NDUFA4L2','NOTCH3',
            'CXCL1','CXCL2','IL6','CEBPD','CLU','CTGF','HGF','HSPA6','DNAJB1','MYC','AFT4','PLAU','CHI3L1','MMP3','IL1R1','IL13RA2','TNFSF11','MMP10','OSMR','IL11','STRA6','FAP','WNT2','TWIST1','IL24',
            'ACTG2','HHIP','CNN1',
            'MYH11','ACTA2','TAGLN',
            'KRT18','SLPI','UPK3B','MSLN','CALB2','WT1','KLK11','ITLN1',
            'WSB1','DDX17','CTNNB1',
            'RBP1','STAR','STMN1',
            'CXCL12','CD74','HLA-DRB1','HLA-DRA',
            'ADAMDEC1','CCL8','APOE','APOC1',
            'LIMCH1','A2M','ADH1B',
            'PRG4','CRTAC1',
            'CXCL14','VSTM2A','SOX6','COL4A5','COL4A6','TSLP','FRZB','BMP5','BMP2','CPM','F3'
        ]
    },
    "Endothelial": {
        "description": "Endothelial cell subtype annotation",
        "markers": [
            'S100B','ALDH1A1','CDH19','GFRA3','NRXN1','SCN7A',
            'GJA4','HEY1','DKK2','IGFBP3','EFNB2',
            'ACKR1','VWF',
            'RGCC','VWA1','IL7R','FCN3','MT1M',
            'TFF3','LYVE1','MMRN1','CCL21',
            'MYL9','ACTA2','TINAGL1','NOTCH3','LAMC3','PDGFRB','RGS5','LPL','HIGD1B','STEAP4','NDUFA4L2'
        ]
    },
    "SMC": {
        "description": "Smooth muscle cell subtype annotation",
        "markers": [
            # Smooth muscle
            'ACTA2', 'TAGLN', 'MYH11', 'CNN1', 'MYLK', 'TPM1', 'TPM2', 'CARMN',
            # Pericyte
            'RGS5', 'PDGFRB', 'NOTCH3', 'MCAM', 'KCNJ8', 'ABCC9',
            # Vascular SMC
            'MYH11', 'ACTA2', 'TAGLN',
            # Proliferation
            'MKI67', 'TOP2A'
        ]
    }
}

# ==================== LOGGING UTILS ====================

def log_msg(msg: str):
    """Print log message if VERBOSE is True."""
    if VERBOSE:
        print(msg)

def log_step(step_num, step_name):
    """Print nicely formatted step header."""
    log_msg("\n" + "=" * 70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("=" * 70)

# ==================== VALIDATION FUNCTIONS (✅ FIX 1) ====================

def validate_cell_type(cell_type: str, marker_configs: dict) -> bool:
    """
    Validate cell type parameter at startup
    
    Args:
        cell_type: User-specified cell type
        marker_configs: Dictionary of available configurations
        
    Returns:
        True if valid
        
    Raises:
        ValueError: If cell type is invalid
    """
    if cell_type not in marker_configs:
        available = ', '.join(sorted(marker_configs.keys()))
        raise ValueError(
            f"\n❌ Invalid CELL_TYPE: '{cell_type}'\n"
            f"Available options: {available}\n"
            f"Please update CELL_TYPE in the script or use --cell-type flag"
        )
    return True

# ==================== CLUSTER ID STANDARDIZATION (✅ FIX 2) ====================

def standardize_cluster_id(cluster_value):
    """
    Standardize cluster ID format for consistent mapping across all operations.
    This is critical for reliable annotation mapping.
    
    Handles various input types:
    - int: 0 → "0"
    - float: 1.0 → "1", 1.5 → "1.5"
    - str: "  2  " → "2", "cluster_3" → "cluster_3"
    - categorical: converts to string
    
    Args:
        cluster_value: Original cluster identifier (any type)
        
    Returns:
        Standardized string representation
    """
    # Convert to string and remove leading/trailing whitespace
    s = str(cluster_value).strip()
    
    # If it looks like a float integer (e.g., "1.0"), convert to int string
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except (ValueError, AttributeError):
        pass
    
    return s

# ==================== PATH RESOLUTION (✅ FIX 4) ====================

def _resolve_input_path(input_path, cell_type):
    """
    Resolve input file path with clear error messages.
    
    Args:
        input_path: User-provided path (may be None)
        cell_type: Cell type for default path construction
        
    Returns:
        Validated file path
        
    Raises:
        FileNotFoundError: If no valid path found with helpful suggestions
    """
    # If explicit path provided, validate and return
    if input_path is not None:
        p = Path(input_path)
        if not p.exists():
            raise FileNotFoundError(
                f"\n❌ Input file not found: {input_path}\n"
                f"Please check the path and try again"
            )
        return str(p)
    
    # Try default locations
    base = Path(f"/home/h2048/data/py/1029/bbknn_celltype_analysis/{cell_type}")
    candidates = [
        base / "output_universal" / f"{cell_type.lower()}_bbknn_integrated.h5ad",
        base / "output_optimized" / f"{cell_type.lower()}_bbknn_integrated.h5ad",
        base / "output" / f"{cell_type.lower()}_bbknn_integrated.h5ad",
        base / f"{cell_type.lower()}_bbknn_integrated.h5ad",
    ]
    
    # Try each candidate
    for p in candidates:
        if p.exists():
            return str(p)
    
    # None found - provide helpful error message
    raise FileNotFoundError(
        f"\n❌ No BBKNN-integrated data found for cell type: {cell_type}\n"
        f"\n📂 Searched locations:\n" +
        '\n'.join(f"  {i+1}. {p}" for i, p in enumerate(candidates)) +
        f"\n\n💡 Please either:\n"
        f"  1. Run BBKNN integration first to create the input file, or\n"
        f"  2. Specify input path explicitly with --input flag\n"
        f"\nExample: --input /path/to/your/data.h5ad"
    )

def _default_output_dir(cell_type: str, resolution: str) -> str:
    """Generate default output directory path"""
    base = f"/home/h2048/data/py/1029/bbknn_annotation/{cell_type}"
    return f"{base}/res_{resolution}"

# ==================== CLUSTER KEY RESOLUTION ====================

def resolve_cluster_key(adata, preferred_key: str) -> str:
    """
    Resolve a best-matching cluster key from adata.obs.
    Tries exact match, common variants, and then a fuzzy 'contains' search for 'leiden'.
    """
    obs_cols = list(adata.obs.columns)
    if preferred_key in obs_cols:
        return preferred_key

    candidates = []
    if preferred_key.startswith("leiden_bbknn_res"):
        res = preferred_key.replace("leiden_bbknn_res", "")
        # allow both "resX" and "res_X" conventions
        candidates += [
            f"leiden_bbknn_res{res}",
            f"leiden_bbknn_res_{res}",
            f"leiden_res{res}",
            f"leiden_res_{res}",
            f"leiden_{res}",
            "leiden_bbknn",
            "leiden"
        ]
    elif preferred_key == "leiden_bbknn":
        candidates += ["leiden_bbknn", "leiden", "leiden_res_1.0", "leiden_res1.0"]

    for c in candidates:
        if c in obs_cols:
            return c
    
    # last resort: any column containing 'leiden'
    leiden_cols = [c for c in obs_cols if 'leiden' in c.lower()]
    if leiden_cols:
        return leiden_cols[0]
    
    raise KeyError(
        f"❌ Cluster key '{preferred_key}' not found in data.\n"
        f"Available clustering columns: {leiden_cols if leiden_cols else 'None found'}\n"
        f"Please run clustering first or specify correct cluster key"
    )

# ==================== MARKER AVAILABILITY (CASE-INSENSITIVE) ====================

def _build_upper_to_original(var_names) -> dict:
    """
    Build a mapping from UPPERCASE gene symbol to the first original var_name
    to support case-insensitive marker presence checks without mutating adata.var_names.
    """
    upper = [str(g).upper() for g in var_names]
    first_idx = {}
    for i, g in enumerate(upper):
        if g not in first_idx:
            first_idx[g] = i
    return {g: var_names[idx] for g, idx in first_idx.items()}

def normalize_gene_symbols(var_names):
    """
    Return a deduplicated uppercase Index (for internal checks only).
    Does not modify adata.var_names.
    """
    up = pd.Index([str(g).upper() for g in var_names])
    _, first_idx = np.unique(up, return_index=True)
    return up[first_idx]

# ==================== DE MATRIX PICKER ====================

def _pick_anndata_for_de(adata):
    """
    Decide which matrix to use for differential expression (DE).
    Preference: layers['counts'] → raw → X.
    Returns: (adata_for_de, use_raw_flag, de_source_label)
    """
    if 'counts' in adata.layers:
        adata_tmp = adata.copy()
        adata_tmp.raw = None  # avoid Scanpy auto-using raw
        adata_tmp.X = adata.layers['counts']  # use counts explicitly
        return adata_tmp, False, "counts"
    if adata.raw is not None:
        return adata, True, "raw"
    return adata, False, "X"

# ==================== TINY GROUP FILTER & SUBSAMPLING (✅ FIX 6) ====================

def _filter_tiny_groups(adata, groupby: str, min_cells=MIN_CELLS_PER_CLUSTER):
    """
    Remove cells from groups (clusters) that are too small for stable DE.
    Returns filtered AnnData and the kept group names.
    """
    vc = adata.obs[groupby].astype(str).value_counts()
    keep_groups = set(vc[vc >= min_cells].index.astype(str))
    mask = adata.obs[groupby].astype(str).isin(keep_groups)
    return adata[mask].copy(), sorted(list(keep_groups))

def balanced_subsample_by_group(adata, groupby: str, per_group=MAX_CELLS_PER_GROUP_PLOTTING, 
                               random_state=0):
    """
    Memory-efficient balanced subsampling per group for heavy plotting tasks.
    
    ✅ FIXED: Uses boolean mask instead of growing list for better memory efficiency.
    
    Args:
        adata: AnnData object
        groupby: Column name to group by
        per_group: Maximum cells per group (None to skip subsampling)
        random_state: Random seed for reproducibility
        
    Returns:
        Subsampled AnnData object
    """
    if per_group is None or per_group <= 0:
        return adata
    
    np.random.seed(random_state)
    
    # ✅ Use boolean mask instead of growing list (memory efficient)
    mask = np.zeros(adata.n_obs, dtype=bool)
    
    for g, group_idx in adata.obs.groupby(groupby, observed=True).indices.items():
        n_take = min(per_group, len(group_idx))
        if n_take > 0:
            # Randomly select indices
            selected = np.random.choice(group_idx, size=n_take, replace=False)
            mask[selected] = True
    
    n_selected = mask.sum()
    pct = n_selected / adata.n_obs * 100
    log_msg(f"   Balanced subsample: {n_selected:,} / {adata.n_obs:,} cells ({pct:.1f}%)")
    
    return adata[mask].copy()

# ==================== SCANPY RGG TO DATAFRAME (COMPATIBLE) ====================

def _rgg_to_dataframe(adata):
    """
    Robustly convert rank_genes_groups to a tidy DataFrame compatible across Scanpy versions.
    Required columns (normalized): cluster, gene, avg_log2FC, pvals, pvals_adj, pct.1, pct.2
    """
    try:
        import scanpy as sc
        df = sc.get.rank_genes_groups_df(adata, group=None)
        # unify field names
        rename_map = {
            'group': 'cluster',
            'names': 'gene',
            'logfoldchanges': 'avg_log2FC',
            'pvals': 'pvals',
            'pvals_adj': 'pvals_adj',
            'pts': 'pct.1',
            'pts_rest': 'pct.2'
        }
        for k, v in rename_map.items():
            if k in df.columns and v not in df.columns:
                df = df.rename(columns={k: v})
        # Ensure required columns
        for col in ['cluster', 'gene', 'avg_log2FC', 'pct.1', 'pct.2', 'pvals', 'pvals_adj']:
            if col not in df.columns:
                df[col] = np.nan
        return df
    except Exception:
        # manual extraction for older Scanpy versions
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        rows = []
        for g in groups:
            names = result['names'][g]
            scores = result.get('scores', {}).get(g, [np.nan]*len(names))
            lfc = result.get('logfoldchanges', {})
            lfc = lfc.get(g, [np.nan]*len(names)) if isinstance(lfc, dict) else [np.nan]*len(names)
            pvals = result.get('pvals', {}).get(g, [np.nan]*len(names))
            pvals_adj = result.get('pvals_adj', {}).get(g, [np.nan]*len(names))
            pts = result.get('pts', {}).get(g, [np.nan]*len(names))
            pts_rest = result.get('pts_rest', {}).get(g, [np.nan]*len(names))
            for i, gene in enumerate(names):
                rows.append({
                    'cluster': g,
                    'gene': gene,
                    'avg_log2FC': lfc[i],
                    'pvals': pvals[i],
                    'pvals_adj': pvals_adj[i],
                    'pct.1': pts[i],
                    'pct.2': pts_rest[i]
                })
        return pd.DataFrame(rows)

# ==================== COLOR PALETTE GENERATOR (✅ FIX 8) ====================

def generate_color_palette(n_colors: int, base_palette: dict = None) -> list:
    """
    Generate a color palette for n_colors, extending beyond predefined colors if needed.
    
    ✅ FIXED: Better handling for >20 cell types using colormap interpolation.
    
    Args:
        n_colors: Number of colors needed
        base_palette: Optional base color dictionary
        
    Returns:
        List of color hex codes
    """
    if base_palette is None:
        base_palette = CELL_TYPE_COLORS
    
    colors = []
    base_colors = list(base_palette.values())
    
    if n_colors <= len(base_colors):
        return base_colors[:n_colors]
    
    # Use base colors first
    colors.extend(base_colors)
    
    # Generate additional colors using tab20 + tab20b + tab20c for variety
    remaining = n_colors - len(colors)
    
    if remaining > 0:
        # Combine multiple colormaps for better variety
        cmaps = [plt.cm.tab20, plt.cm.tab20b, plt.cm.tab20c, plt.cm.Set3]
        cmap_colors = []
        
        for cmap in cmaps:
            n_cmap = cmap.N if hasattr(cmap, 'N') else 20
            for i in range(n_cmap):
                cmap_colors.append(mcolors.to_hex(cmap(i)))
        
        # Add as many as needed
        colors.extend(cmap_colors[:remaining])
    
    return colors[:n_colors]

# ==================== CORE FUNCTIONS ====================

def load_integrated_data(h5ad_path, cluster_key_spec):
    """
    Load BBKNN-integrated AnnData object and resolve cluster key.
    
    ✅ Includes Scanpy version checking
    """
    log_step(1, "Loading Integrated Data")

    try:
        import scanpy as sc
        sc_version = sc.__version__
        log_msg(f"   Scanpy version: {sc_version}")
        
        # ✅ Version check (P2 Fix #9)
        if HAS_PACKAGING:
            if pkg_version.parse(sc_version) < pkg_version.parse("1.9.0"):
                warnings.warn(
                    "\n⚠️  Scanpy < 1.9.0 detected. Some features may have compatibility issues.\n"
                    "   Recommended: pip install --upgrade scanpy\n"
                    f"   Current version: {sc_version}"
                )
    except ImportError:
        raise ImportError(
            "❌ scanpy is required but not found.\n"
            "Install with: pip install scanpy"
        )

    if not os.path.exists(h5ad_path):
        raise FileNotFoundError(f"❌ Input file not found: {h5ad_path}")

    log_msg(f"\n   Loading data from:\n   {h5ad_path}")
    
    # Use progress bar if enabled
    if USE_PROGRESS_BARS:
        with tqdm(total=1, desc="Loading h5ad", unit="file") as pbar:
            adata = sc.read_h5ad(h5ad_path)
            pbar.update(1)
    else:
        adata = sc.read_h5ad(h5ad_path)

    log_msg(f"\n   ✓ Data loaded successfully")
    log_msg(f"   Cells: {adata.n_obs:,}")
    log_msg(f"   Genes: {adata.n_vars:,}")

    # Determine cluster key based on resolution selection spec
    if cluster_key_spec == "default":
        cluster_key_candidate = "leiden_bbknn"
    else:
        cluster_key_candidate = f"leiden_bbknn_res{cluster_key_spec}"

    # Resolve best cluster key in the object
    resolved_key = resolve_cluster_key(adata, cluster_key_candidate)
    log_msg(f"\n   Using clustering: {resolved_key}")

    cluster_counts = adata.obs[resolved_key].value_counts().sort_index()
    log_msg(f"   Number of clusters: {len(cluster_counts)}")
    log_msg(f"   Cluster size range: {cluster_counts.min():,} - {cluster_counts.max():,} cells")

    return adata, resolved_key


def check_marker_availability(adata, markers, cell_type):
    """
    Check availability of marker genes (case-insensitive) in dataset.
    """
    log_step(2, f"Checking Marker Genes for {cell_type} Annotation")

    description = MARKER_CONFIGS[cell_type]["description"]
    log_msg(f"   Annotation purpose: {description}")
    log_msg(f"   Total markers: {len(markers)}")

    upper_to_original = _build_upper_to_original(adata.var_names)
    marker_up = [g.upper() for g in markers]

    available_markers_orig = []
    missing_markers = []
    for g in marker_up:
        if g in upper_to_original:
            available_markers_orig.append(upper_to_original[g])
        else:
            missing_markers.append(g)

    log_msg(f"   Available: {len(available_markers_orig)} "
            f"({len(available_markers_orig)/len(markers)*100:.1f}%)")

    if missing_markers:
        log_msg(f"\n   Missing markers ({len(missing_markers)}):")
        for gene in missing_markers[:15]:
            log_msg(f"     - {gene}")
        if len(missing_markers) > 15:
            log_msg(f"     ... and {len(missing_markers)-15} more")

    if len(available_markers_orig) < len(markers) * 0.5:
        log_msg("\n   ⚠️  Warning: Less than 50% of markers available!")
        log_msg("   Consider:")
        log_msg("     1. Check gene naming (SYMBOL vs ENSEMBL)")
        log_msg("     2. Verify species match (human vs mouse)")
        log_msg("     3. Update marker list for your dataset")

    return available_markers_orig


def generate_dotplot(adata, available_markers, cluster_key, output_dir):
    """
    Generate marker gene DotPlot with explicit use_raw selection.
    """
    log_step(3, "Generating Marker Gene DotPlot")

    if not available_markers:
        log_msg("   No markers available for dotplot")
        return

    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    try:
        import scanpy as sc
        
        # Balanced subsample for plotting if dataset is huge
        adata_plot = balanced_subsample_by_group(
            adata, groupby=cluster_key, per_group=MAX_CELLS_PER_GROUP_PLOTTING, random_state=0
        )

        log_msg(f"\n   Creating dotplot with {len(available_markers)} markers...")

        n_clusters = adata_plot.obs[cluster_key].nunique()
        fig_width = max(12, len(available_markers) * 0.35)
        fig_height = max(8, n_clusters * 0.5)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sc.pl.dotplot(
            adata_plot,
            var_names=available_markers,
            groupby=cluster_key,
            ax=ax,
            show=False,
            standard_scale='var',
            cmap='Reds',
            use_raw=(adata_plot.raw is not None)
        )
        ax.set_title(f"Marker Gene Expression • {CELL_TYPE} • {cluster_key}\n"
                    f"n={adata_plot.n_obs:,} cells (balanced subsample)")

        output_file = fig_dir / f"dotplot_markers_{CELL_TYPE}.{FIGURE_FORMAT}"
        plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
        plt.close()

        log_msg(f"   ✓ Dotplot saved to: {output_file}")

    except Exception as e:
        log_msg(f"   ✗ Failed to generate dotplot: {e}")
        import traceback
        if VERBOSE:
            log_msg(traceback.format_exc())


def find_all_markers(adata, cluster_key, output_dir):
    """
    Find marker genes for all clusters using scanpy.tl.rank_genes_groups.
    
    ✅ FIXED (P2 Fix #5): Reports excluded clusters explicitly
    """
    log_step(4, "Finding Cluster Marker Genes")

    import scanpy as sc

    log_msg(f"   Cluster key: {cluster_key}")
    log_msg(f"   Parameters:")
    log_msg(f"     - min_pct: {MARKER_MIN_PCT}")
    log_msg(f"     - logfc_threshold: {MARKER_LOGFC_THRESHOLD}")
    log_msg(f"     - min_cells_per_cluster: {MIN_CELLS_PER_CLUSTER}")

    try:
        # Pick DE matrix source
        adata_de, use_raw_flag, de_src = _pick_anndata_for_de(adata)
        log_msg(f"\n   DE matrix source: {de_src} (use_raw={use_raw_flag})")

        # ✅ Check cluster sizes BEFORE filtering
        original_cluster_counts = adata_de.obs[cluster_key].value_counts()
        log_msg(f"\n   Original cluster distribution:")
        log_msg(f"     Total clusters: {len(original_cluster_counts)}")
        log_msg(f"     Size range: {original_cluster_counts.min():,} - "
                f"{original_cluster_counts.max():,} cells")

        # Filter tiny groups
        adata_de2, kept_groups = _filter_tiny_groups(adata_de, cluster_key, MIN_CELLS_PER_CLUSTER)

        # ✅ Report excluded clusters (P2 Fix #5)
        all_clusters = set(original_cluster_counts.index.astype(str))
        kept_clusters_set = set(kept_groups)
        excluded_clusters = sorted(all_clusters - kept_clusters_set, 
                                  key=lambda x: (len(x), x))

        if excluded_clusters:
            log_msg(f"\n   ⚠️  Excluded {len(excluded_clusters)} small clusters "
                    f"(< {MIN_CELLS_PER_CLUSTER} cells):")
            for cl in excluded_clusters:
                n_cells = original_cluster_counts[cl]
                log_msg(f"     Cluster {cl}: {n_cells} cells")
            log_msg(f"\n   ⚠️  These clusters will have NO marker genes in output")
            log_msg(f"   Consider lowering --min-cells-per-cluster if needed")

        if len(kept_groups) == 0:
            log_msg(f"\n   ❌ No clusters pass the min-cells threshold; skipping DE.")
            log_msg(f"   Recommendation: Lower --min-cells-per-cluster (current: {MIN_CELLS_PER_CLUSTER})")
            return None

        log_msg(f"\n   Proceeding with {len(kept_groups)} clusters for DE analysis...")

        # Run DE with progress bar
        log_msg("\n   Running differential expression analysis...")
        if USE_PROGRESS_BARS:
            with tqdm(total=1, desc="Rank genes", unit="analysis") as pbar:
                sc.tl.rank_genes_groups(
                    adata_de2,
                    groupby=cluster_key,
                    method='wilcoxon',
                    use_raw=use_raw_flag,
                    pts=True
                )
                pbar.update(1)
        else:
            sc.tl.rank_genes_groups(
                adata_de2,
                groupby=cluster_key,
                method='wilcoxon',
                use_raw=use_raw_flag,
                pts=True
            )
        log_msg("     ✓ Analysis complete")

        # Extract results to DataFrame
        log_msg("\n   Extracting marker genes...")
        df = _rgg_to_dataframe(adata_de2)
        
        # Ensure required columns
        need_cols = ['cluster', 'gene', 'avg_log2FC', 'pct.1', 'pct.2', 'pvals', 'pvals_adj']
        for c in need_cols:
            if c not in df.columns:
                df[c] = np.nan

        # Thresholding
        df_thr = df[(df['pct.1'] >= MARKER_MIN_PCT) &
                    (df['avg_log2FC'].fillna(0).abs() >= MARKER_LOGFC_THRESHOLD)].copy()

        if len(df_thr) > 0:
            # ✅ UTF-8 encoding (P1 Fix #3)
            output_file = Path(output_dir) / 'cluster_markers.csv'
            df_thr.to_csv(output_file, index=False, encoding='utf-8')
            log_msg(f"   ✓ Marker genes saved to: {output_file}")
            log_msg(f"   Total markers (post-threshold): {len(df_thr)}")

            # Display top markers per cluster
            log_msg("\n   Top markers per cluster:")
            for cluster in sorted(df_thr['cluster'].astype(str).unique(), key=lambda x: (len(x), x)):
                cluster_markers = (
                    df_thr[df_thr['cluster'].astype(str) == str(cluster)]
                    .sort_values('avg_log2FC', ascending=False)
                    .head(5)
                )
                top_genes = ', '.join(cluster_markers['gene'].astype(str).values)
                log_msg(f"     Cluster {cluster}: {top_genes}")
        else:
            log_msg("\n   ⚠️  No significant markers found after thresholding")
            log_msg(f"   Consider adjusting thresholds:")
            log_msg(f"     - Lower min_pct (current: {MARKER_MIN_PCT})")
            log_msg(f"     - Lower logfc_threshold (current: {MARKER_LOGFC_THRESHOLD})")
            df_thr = None

        # record DE matrix source on adata for the final report
        adata.uns['__de_source__'] = de_src

        return df_thr

    except Exception as e:
        log_msg(f"\n   ❌ Failed to find markers: {e}")
        import traceback
        if VERBOSE:
            log_msg(traceback.format_exc())
        return None


def generate_marker_heatmaps(adata, cluster_key, marker_df, output_dir):
    """
    Generate heatmaps of top marker genes per cluster.
    
    ✅ FIXED (P2 Fix #7): Skips empty batches
    """
    log_step(5, "Generating Marker Gene Heatmaps")

    import scanpy as sc

    if marker_df is None or len(marker_df) == 0:
        log_msg("   No markers available for heatmap")
        return

    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    try:
        log_msg(f"   Selecting top {TOP_N_MARKERS} markers per cluster...")

        top_markers_per_cluster = (
            marker_df.sort_values('avg_log2FC', ascending=False)
            .groupby('cluster', as_index=False)
            .head(TOP_N_MARKERS)
        )

        top_marker_genes = top_markers_per_cluster['gene'].astype(str).unique().tolist()
        log_msg(f"   Total unique markers: {len(top_marker_genes)}")

        if len(top_marker_genes) == 0:
            log_msg("   ⚠️  No genes to plot")
            return

        # Balanced subsample for plotting only
        adata_plot = balanced_subsample_by_group(
            adata, groupby=cluster_key, per_group=MAX_CELLS_PER_GROUP_PLOTTING, random_state=0
        )

        n_batches = int(np.ceil(len(top_marker_genes) / HEATMAP_BATCH_SIZE))
        log_msg(f"\n   Generating {n_batches} heatmap(s)...")

        # ✅ Use progress bar (P2 Fix #9)
        batch_range = range(n_batches)
        if USE_PROGRESS_BARS:
            batch_range = tqdm(batch_range, desc="Heatmap batches", unit="batch")

        for batch_idx in batch_range:
            start_idx = batch_idx * HEATMAP_BATCH_SIZE
            end_idx = min((batch_idx + 1) * HEATMAP_BATCH_SIZE, len(top_marker_genes))
            batch_genes = top_marker_genes[start_idx:end_idx]

            # ✅ Skip empty batches (P2 Fix #7)
            if len(batch_genes) == 0:
                log_msg(f"     Batch {batch_idx + 1}: Empty, skipping")
                continue

            if not USE_PROGRESS_BARS:
                log_msg(f"     Batch {batch_idx + 1}: {len(batch_genes)} genes")

            fig, ax = plt.subplots(figsize=(12, max(8, len(batch_genes) * 0.3)))
            sc.pl.heatmap(
                adata_plot,
                var_names=batch_genes,
                groupby=cluster_key,
                ax=ax,
                show=False,
                cmap='RdBu_r',
                standard_scale='var',
                dendrogram=True,
                use_raw=(adata_plot.raw is not None)
            )
            ax.set_title(f"Top Marker Genes • {CELL_TYPE} • {cluster_key}\n"
                        f"Batch {batch_idx + 1}/{n_batches} • n={adata_plot.n_obs:,} cells")

            output_file = fig_dir / f"heatmap_top_markers_batch{batch_idx + 1}.{FIGURE_FORMAT}"
            plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
            plt.close()

        log_msg(f"\n   ✓ Heatmaps saved to: {fig_dir}")

    except Exception as e:
        log_msg(f"\n   ✗ Failed to generate heatmaps: {e}")
        import traceback
        if VERBOSE:
            log_msg(traceback.format_exc())


def create_annotation_template(adata, cluster_key, output_dir):
    """
    Create annotation template CSV file (no header).
    
    ✅ FIXED (P1 Fix #2): Standardizes cluster IDs and creates reference file
    ✅ FIXED (P1 Fix #3): UTF-8 encoding
    """
    log_step(6, "Creating Annotation Template")

    # ✅ Get unique clusters and standardize (P1 Fix #2)
    raw_clusters = adata.obs[cluster_key].unique()
    clusters = sorted([standardize_cluster_id(c) for c in raw_clusters], 
                     key=lambda x: (len(x), x))

    # ✅ Check for duplicate clusters after standardization
    if len(clusters) != len(set(clusters)):
        warnings.warn(
            "\n⚠️  Warning: Some clusters have identical standardized IDs.\n"
            "   This may cause annotation conflicts. Check cluster_id_reference.csv"
        )

    template_df = pd.DataFrame({
        'Cluster': clusters,
        'CellType': ['Unknown'] * len(clusters)
    })

    # ✅ Save with UTF-8 encoding (P1 Fix #3)
    template_file = Path(output_dir) / "Annotation_template.csv"
    template_df.to_csv(template_file, index=False, header=False, encoding='utf-8')

    log_msg(f"   ✓ Annotation template created: {template_file}")
    log_msg(f"   Number of clusters: {len(clusters)}")
    log_msg(f"   Cluster IDs: {', '.join(str(c) for c in clusters[:10])}")
    if len(clusters) > 10:
        log_msg(f"                ... and {len(clusters)-10} more")

    # ✅ Create a reference file showing original vs standardized IDs (P1 Fix #2)
    reference_df = pd.DataFrame({
        'Original': raw_clusters,
        'Standardized': [standardize_cluster_id(c) for c in raw_clusters]
    })
    reference_file = Path(output_dir) / "cluster_id_reference.csv"
    reference_df.to_csv(reference_file, index=False, encoding='utf-8')
    log_msg(f"   ✓ Cluster ID reference saved: {reference_file}")
    log_msg(f"      (Use this to verify cluster ID mapping)")

    log_msg("\n   📋 Next steps:")
    log_msg("   1. Review dotplot and marker heatmaps in figures/")
    log_msg("   2. Copy Annotation_template.csv to Annotation.csv:")
    log_msg(f"      cp {template_file} {Path(output_dir) / 'Annotation.csv'}")
    log_msg("   3. Edit Annotation.csv to assign cell types")
    log_msg("   4. Run script again with RUN_MODE='post_annotation'")
    log_msg("\n   💡 Template format (NO header, format: Cluster,CellType):")
    log_msg("   0,Basal")
    log_msg("   1,Secretory")
    log_msg("   2,Ciliated")
    log_msg("   ...")



def apply_annotation(adata, cluster_key, annotation_file, output_dir):
    """
    Apply cell type annotation from CSV file (no header).
    
    ✅ FIXED (P1 Fix #2): Robust cluster ID standardization and mapping
    ✅ FIXED (P1 Fix #3): UTF-8 encoding
    """
    log_step(7, "Applying Cell Type Annotation")

    annotation_path = Path(output_dir) / annotation_file
    if not annotation_path.exists():
        raise FileNotFoundError(
            f"\n❌ Annotation file not found: {annotation_path}\n"
            f"\nPlease create this file by:\n"
            f"  1. Copy Annotation_template.csv to Annotation.csv:\n"
            f"     cp {Path(output_dir) / 'Annotation_template.csv'} {annotation_path}\n"
            f"  2. Edit Annotation.csv with cell type assignments\n"
            f"  3. Ensure format: Cluster,CellType (NO header)\n"
            f"\nExample content:\n"
            f"  0,Basal\n"
            f"  1,Secretory\n"
            f"  2,Ciliated"
        )

    log_msg(f"   Reading annotation from:\n   {annotation_path}")

    # ✅ Read annotation with UTF-8 encoding and force string type (P1 Fix #3)
    try:
        annotation_df = pd.read_csv(
            annotation_path, 
            header=None, 
            names=['Cluster', 'CellType'],
            encoding='utf-8',
            dtype=str  # Force string type to avoid type mismatches
        )
    except Exception as e:
        raise ValueError(
            f"\n❌ Failed to read annotation file: {e}\n"
            f"Please check file format and encoding"
        )

    # ✅ Standardize cluster IDs in annotation (P1 Fix #2)
    annotation_df['Cluster'] = annotation_df['Cluster'].apply(standardize_cluster_id)
    annotation_df['CellType'] = annotation_df['CellType'].str.strip()

    # ✅ Check for empty cell types
    empty_mask = annotation_df['CellType'].str.len() == 0
    if empty_mask.any():
        empty_clusters = annotation_df.loc[empty_mask, 'Cluster'].tolist()
        raise ValueError(
            f"\n❌ Empty cell type assignments found for clusters: {empty_clusters}\n"
            f"Please assign cell types to all clusters in {annotation_file}"
        )

    # ✅ Check for duplicate cluster entries
    duplicates = annotation_df['Cluster'].duplicated()
    if duplicates.any():
        dup_clusters = annotation_df.loc[duplicates, 'Cluster'].tolist()
        warnings.warn(
            f"\n⚠️  Duplicate cluster entries found: {dup_clusters}\n"
            f"   Keeping first occurrence for each cluster"
        )
        annotation_df = annotation_df[~duplicates].copy()

    log_msg(f"   Annotations loaded: {len(annotation_df)} clusters")

    # Create mapping dictionary
    cluster_to_celltype = dict(zip(annotation_df['Cluster'], annotation_df['CellType']))

    # ✅ Standardize adata cluster IDs (P1 Fix #2)
    adata.obs['_cluster_std'] = adata.obs[cluster_key].apply(standardize_cluster_id)

    # Apply mapping
    adata.obs['cell_type'] = adata.obs['_cluster_std'].map(cluster_to_celltype)

    # ✅ Check unmapped clusters
    unmapped_mask = adata.obs['cell_type'].isna()
    n_unmapped = unmapped_mask.sum()

    if n_unmapped > 0:
        unmapped_clusters = adata.obs.loc[unmapped_mask, '_cluster_std'].unique()
        warnings.warn(
            f"\n⚠️  Warning: {n_unmapped:,} cells in {len(unmapped_clusters)} "
            f"clusters could not be mapped:\n"
            f"   Unmapped clusters: {sorted(unmapped_clusters)}\n"
            f"   Check if these clusters exist in {annotation_file}\n"
            f"   These cells will be labeled as 'Unknown'"
        )
        # Fill unmapped as 'Unknown'
        adata.obs['cell_type'] = adata.obs['cell_type'].fillna('Unknown')

    # Clean up temporary column
    adata.obs.drop('_cluster_std', axis=1, inplace=True)

    # Summary
    log_msg("\n   Cell type distribution:")
    celltype_counts = adata.obs['cell_type'].value_counts().sort_values(ascending=False)
    for celltype, count in celltype_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"     {celltype}: {count:,} cells ({pct:.1f}%)")

    # ✅ Validation: Check if high proportion of 'Unknown'
    if 'Unknown' in celltype_counts.index:
        unknown_count = celltype_counts['Unknown']
        if unknown_count > adata.n_obs * 0.05:  # More than 5% unknown
            warnings.warn(
                f"\n⚠️  High proportion of 'Unknown' cells: {unknown_count:,} "
                f"({unknown_count/adata.n_obs*100:.1f}%)\n"
                f"   Consider reviewing annotation completeness"
            )

    # Save a copy of the final mapping with UTF-8
    mapping_out = Path(output_dir) / "Annotation_applied.csv"
    annotation_df.to_csv(mapping_out, index=False, encoding='utf-8')
    log_msg(f"\n   ✓ Annotation mapping saved: {mapping_out}")

    return adata


def generate_celltype_visualizations(adata, cluster_key, output_dir):
    """
    Generate cell type UMAP and additional comparison/proportion plots.
    
    ✅ FIXED (P2 Fix #8): Better color palette for >20 cell types
    """
    log_step(8, "Generating Cell Type Visualizations")

    import scanpy as sc

    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ✅ Generate color palette with improved handling for many cell types (P2 Fix #8)
    unique_celltypes = sorted(adata.obs['cell_type'].astype(str).unique(), 
                             key=lambda x: (len(x), x))
    log_msg(f"   Cell types: {len(unique_celltypes)}")

    # Generate colors using improved palette generator
    all_colors = generate_color_palette(len(unique_celltypes), CELL_TYPE_COLORS)
    
    celltype_colors = {}
    for i, celltype in enumerate(unique_celltypes):
        if celltype in CELL_TYPE_COLORS:
            celltype_colors[celltype] = CELL_TYPE_COLORS[celltype]
        else:
            celltype_colors[celltype] = all_colors[i]

    # 1) Cell type UMAP
    try:
        log_msg("\n   Generating cell type UMAP...")
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(
            adata,
            color='cell_type',
            ax=ax,
            show=False,
            title=f'{CELL_TYPE} - Cell Type Annotation\n'
                  f'Resolution: {RESOLUTION} • n={adata.n_obs:,} cells',
            palette=celltype_colors,
            legend_loc='right margin',
            frameon=False
        )

        output_file = fig_dir / f"umap_celltype_{CELL_TYPE}_res{RESOLUTION}.{FIGURE_FORMAT}"
        plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
        plt.close()
        log_msg(f"   ✓ Cell type UMAP saved: {output_file}")

    except Exception as e:
        log_msg(f"   ✗ Failed to generate cell type UMAP: {e}")

    # 2) Side-by-side: clusters and cell types
    try:
        log_msg("\n   Generating side-by-side comparison...")
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))

        sc.pl.umap(adata, color=cluster_key, ax=axes[0], show=False,
                   title=f'Clusters ({cluster_key})', frameon=False)
        sc.pl.umap(adata, color='cell_type', ax=axes[1], show=False,
                   title='Cell Type Annotation', palette=celltype_colors, frameon=False)

        plt.tight_layout()
        output_file = fig_dir / f"umap_comparison_{CELL_TYPE}_res{RESOLUTION}.{FIGURE_FORMAT}"
        plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
        plt.close()
        log_msg(f"   ✓ Comparison plot saved: {output_file}")

    except Exception as e:
        log_msg(f"   ✗ Failed to generate comparison plot: {e}")

    # 3) Cell type proportions bar plot
    try:
        log_msg("\n   Generating cell type proportions plot...")
        celltype_counts = adata.obs['cell_type'].value_counts()

        fig, ax = plt.subplots(figsize=(max(10, len(celltype_counts) * 0.6), 6))
        bars = ax.bar(range(len(celltype_counts)), celltype_counts.values)

        # Color bars
        for bar, celltype in zip(bars, celltype_counts.index):
            bar.set_color(celltype_colors.get(celltype, '#CCCCCC'))

        ax.set_xticks(range(len(celltype_counts)))
        ax.set_xticklabels(celltype_counts.index, rotation=45, ha='right')
        ax.set_ylabel('Number of Cells', fontsize=12)
        ax.set_title(f'Cell Type Distribution - {CELL_TYPE}\n'
                    f'Resolution: {RESOLUTION} • Total: {adata.n_obs:,} cells', 
                    fontsize=14)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add count labels on bars
        for i, (celltype, count) in enumerate(celltype_counts.items()):
            pct = count / adata.n_obs * 100
            ax.text(i, count, f'{count:,}\n({pct:.1f}%)', 
                   ha='center', va='bottom', fontsize=8)

        plt.tight_layout()
        output_file = fig_dir / f"barplot_celltype_{CELL_TYPE}_res{RESOLUTION}.{FIGURE_FORMAT}"
        plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
        plt.close()
        log_msg(f"   ✓ Proportion plot saved: {output_file}")

    except Exception as e:
        log_msg(f"   ✗ Failed to generate proportion plot: {e}")

    # 4) Cell type by batch (if batch info available)
    batch_keys = ['dataset', 'batch', 'sample']
    batch_key = None
    for key in batch_keys:
        if key in adata.obs.columns:
            batch_key = key
            break

    if batch_key:
        try:
            log_msg(f"\n   Generating cell type by batch plot (batch_key: {batch_key})...")

            crosstab = pd.crosstab(
                adata.obs[batch_key],
                adata.obs['cell_type'],
                normalize='index'
            ) * 100

            fig, ax = plt.subplots(figsize=(max(12, len(crosstab.columns) * 0.8),
                                            max(8, len(crosstab.index) * 0.6)))

            crosstab.plot(kind='bar', stacked=True, ax=ax,
                          color=[celltype_colors.get(ct, '#CCCCCC') for ct in crosstab.columns])

            ax.set_ylabel('Percentage (%)', fontsize=12)
            ax.set_xlabel(batch_key.capitalize(), fontsize=12)
            ax.set_title(f'Cell Type Distribution by {batch_key.capitalize()}\n'
                        f'{CELL_TYPE} • Resolution: {RESOLUTION}', fontsize=14)
            ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()

            output_file = fig_dir / f"barplot_celltype_by_batch_{CELL_TYPE}_res{RESOLUTION}.{FIGURE_FORMAT}"
            plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
            plt.close()
            log_msg(f"   ✓ Batch distribution plot saved: {output_file}")

        except Exception as e:
            log_msg(f"   ✗ Failed to generate batch plot: {e}")

    log_msg(f"\n   ✓ All visualizations saved to: {fig_dir}")


def save_annotated_data(adata, output_dir):
    """
    Save annotated AnnData object to h5ad with gzip compression.
    """
    log_step(9, "Saving Annotated Data")

    output_file = Path(output_dir) / f"{CELL_TYPE.lower()}_annotated_res{RESOLUTION}.h5ad"

    log_msg(f"   Saving to:\n   {output_file}")
    
    if USE_PROGRESS_BARS:
        with tqdm(total=1, desc="Saving h5ad", unit="file") as pbar:
            adata.write_h5ad(output_file, compression='gzip')
            pbar.update(1)
    else:
        adata.write_h5ad(output_file, compression='gzip')
    
    log_msg("   ✓ Data saved successfully")

    # Report file size
    file_size_mb = output_file.stat().st_size / (1024 * 1024)
    log_msg(f"   File size: {file_size_mb:.1f} MB")

    return output_file


def generate_summary_report(adata, cluster_key, markers, output_dir, run_mode, processing_time=None):
    """
    Generate analysis summary report with configuration, data overview, and outputs.
    """
    from datetime import datetime
    try:
        import scanpy as sc
        sc_ver = sc.__version__
    except Exception:
        sc_ver = "unknown"

    log_step(10, "Generating Summary Report")

    report = []
    report.append("=" * 70)
    report.append(f"{CELL_TYPE} Cell Type Annotation Summary")
    report.append("=" * 70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if processing_time:
        report.append(f"Processing time: {processing_time:.1f} seconds ({processing_time/60:.1f} minutes)")
    report.append("")
    
    report.append("[ Configuration ]")
    report.append(f"  Cell type: {CELL_TYPE}")
    report.append(f"  Annotation purpose: {MARKER_CONFIGS[CELL_TYPE]['description']}")
    report.append(f"  Resolution: {RESOLUTION}")
    report.append(f"  Cluster key: {cluster_key}")
    report.append(f"  Run mode: {run_mode}")
    report.append(f"  Scanpy version: {sc_ver}")
    de_src = adata.uns.get('__de_source__', 'unknown')
    report.append(f"  DE matrix: {de_src}")
    report.append("")

    report.append("[ Data Overview ]")
    report.append(f"  Cells: {adata.n_obs:,}")
    report.append(f"  Genes: {adata.n_vars:,}")
    report.append("")

    if cluster_key in adata.obs.columns:
        report.append(f"[ Clustering Results ({cluster_key}) ]")
        cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
        report.append(f"  Number of clusters: {len(cluster_counts)}")
        report.append(f"  Cluster size range: {cluster_counts.min():,} - {cluster_counts.max():,} cells")
        report.append("")

    if 'cell_type' in adata.obs.columns:
        report.append("[ Cell Type Annotation ]")
        celltype_counts = adata.obs['cell_type'].value_counts().sort_values(ascending=False)
        report.append(f"  Number of cell types: {len(celltype_counts)}")
        report.append("")
        report.append("  Cell type distribution:")
        for celltype, count in celltype_counts.items():
            pct = count / adata.n_obs * 100
            report.append(f"    {celltype}: {count:,} cells ({pct:.1f}%)")
        report.append("")

    report.append("[ Marker Genes ]")
    report.append(f"  Total markers queried: {len(markers)}")
    # availability re-check
    upper_to_original = _build_upper_to_original(adata.var_names)
    marker_up = [g.upper() for g in markers]
    available = [g for g in marker_up if g in upper_to_original]
    report.append(f"  Available in dataset: {len(available)} ({len(available)/len(markers)*100:.1f}%)")
    report.append("")

    report.append("[ Output Files ]")
    fig_dir = Path(output_dir) / 'figures'
    if RUN_MODE in ["pre_annotation", "full"]:
        report.append(f"  - Dotplot: {fig_dir / f'dotplot_markers_{CELL_TYPE}.{FIGURE_FORMAT}'}")
        if RUN_FIND_MARKERS:
            report.append(f"  - Markers: {Path(output_dir) / 'cluster_markers.csv'}")
            report.append(f"  - Heatmaps: {fig_dir / f'heatmap_top_markers_*.{FIGURE_FORMAT}'}")
        report.append(f"  - Template: {Path(output_dir) / 'Annotation_template.csv'}")
        report.append(f"  - ID Reference: {Path(output_dir) / 'cluster_id_reference.csv'}")
    
    if RUN_MODE in ["post_annotation", "full"]:
        anno_h5ad = Path(output_dir) / f'{CELL_TYPE.lower()}_annotated_res{RESOLUTION}.h5ad'
        if anno_h5ad.exists():
            report.append(f"  - Annotated data: {anno_h5ad}")
            report.append(f"  - Cell type UMAPs: {fig_dir / f'umap_celltype_*.{FIGURE_FORMAT}'}")
            report.append(f"  - Comparison plot: {fig_dir / f'umap_comparison_*.{FIGURE_FORMAT}'}")
            report.append(f"  - Proportion plots: {fig_dir / f'barplot_celltype_*.{FIGURE_FORMAT}'}")

    report.append("")
    report.append("[ Fixes Applied ]")
    report.append("  ✓ Priority 1 (Critical):")
    report.append("    • Cell type parameter validation")
    report.append("    • Standardized cluster ID handling")
    report.append("    • UTF-8 encoding for all CSV operations")
    report.append("    • Improved file path error messages")
    report.append("")
    report.append("  ✓ Priority 2 (Important):")
    report.append("    • DE analysis exclusion reporting")
    report.append("    • Memory-efficient balanced subsampling")
    report.append("    • Empty heatmap batch handling")
    report.append("    • Better color palette generation")
    report.append("    • Progress bars for long operations")
    report.append("")
    report.append("=" * 70)
    report.append("Analysis completed successfully")
    report.append("=" * 70)

    report_text = '\n'.join(report)

    # Save report with UTF-8 encoding
    report_path = Path(output_dir) / "annotation_summary.txt"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report_text)

    log_msg(f"\n   ✓ Report saved to: {report_path}")
    log_msg("\n" + report_text)



# ==================== CLI ARG PARSER ====================

def _parse_args():
    """Parse command line arguments"""
    p = argparse.ArgumentParser(
        description="Universal BBKNN Cell Type Annotation Pipeline (FIXED)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Pre-annotation only (generate markers and template)
  python %(prog)s --cell-type All --run-mode pre_annotation
  
  # Post-annotation only (after manually creating Annotation.csv)
  python %(prog)s --cell-type All --run-mode post_annotation
  
  # Full workflow (requires Annotation.csv for post steps)
  python %(prog)s --cell-type All --run-mode full
  
  # Custom resolution and paths
  python %(prog)s --cell-type T --resolution 1.0 --input custom.h5ad --output results/
  
  # Disable DE and progress bars
  python %(prog)s --no-de --no-progress
        """
    )
    
    # Required/Main arguments
    p.add_argument("--cell-type", default=CELL_TYPE,
                   choices=list(MARKER_CONFIGS.keys()),
                   help="Cell type scope (default: %(default)s)")
    
    p.add_argument("--resolution", default=RESOLUTION,
                   help="Clustering resolution (default: '%(default)s')")
    
    p.add_argument("--run-mode", default=RUN_MODE,
                   choices=["pre_annotation", "post_annotation", "full"],
                   help="Run mode (default: %(default)s)")
    
    # I/O arguments
    p.add_argument("--input", default=None,
                   help="Override input .h5ad path")
    
    p.add_argument("--output", default=None,
                   help="Override output directory")
    
    # DE parameters
    p.add_argument("--no-de", action="store_true",
                   help="Disable FindAllMarkers (skip DE analysis)")
    
    p.add_argument("--min-cells-per-cluster", type=int, default=MIN_CELLS_PER_CLUSTER,
                   help=f"Minimum cells per cluster for DE (default: {MIN_CELLS_PER_CLUSTER})")
    
    p.add_argument("--plot-cells-per-cluster", type=int, default=MAX_CELLS_PER_GROUP_PLOTTING,
                   help=f"Balanced subsample size per cluster for plotting (default: {MAX_CELLS_PER_GROUP_PLOTTING})")
    
    # Other options
    p.add_argument("--no-progress", action="store_true",
                   help="Disable progress bars")
    
    p.add_argument("--quiet", action="store_true",
                   help="Suppress verbose logging")
    
    p.add_argument("--version", action="version",
                   version="%(prog)s v1.2 (FIXED)")
    
    return p.parse_args()


# ==================== MAIN ====================

def main():
    """Main pipeline execution"""
    print("\n" + "="*70)
    print(f"Universal Cell Type Annotation Pipeline - FIXED VERSION")
    print("="*70)
    print("✅ All Priority 1 + Priority 2 fixes applied")
    print("="*70)

    # Parse CLI args and override globals when provided
    args = _parse_args()
    
    global CELL_TYPE, RESOLUTION, RUN_MODE, INPUT_H5AD_PATH, OUTPUT_DIR
    global RUN_FIND_MARKERS, MIN_CELLS_PER_CLUSTER, MAX_CELLS_PER_GROUP_PLOTTING
    global USE_PROGRESS_BARS, VERBOSE

    CELL_TYPE = args.cell_type
    RESOLUTION = args.resolution
    RUN_MODE = args.run_mode
    
    if args.no_de:
        RUN_FIND_MARKERS = False
    
    MIN_CELLS_PER_CLUSTER = args.min_cells_per_cluster
    MAX_CELLS_PER_GROUP_PLOTTING = args.plot_cells_per_cluster
    USE_PROGRESS_BARS = not args.no_progress
    VERBOSE = not args.quiet

    # ✅ Validate cell type (P1 Fix #1)
    try:
        validate_cell_type(CELL_TYPE, MARKER_CONFIGS)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)

    # ✅ Resolve input path with clear errors (P1 Fix #4)
    if args.input:
        INPUT_H5AD_PATH = args.input
    else:
        try:
            INPUT_H5AD_PATH = _resolve_input_path(None, CELL_TYPE)
        except FileNotFoundError as e:
            print(str(e), file=sys.stderr)
            sys.exit(1)

    # Resolve output directory
    if args.output:
        OUTPUT_DIR = args.output
    else:
        OUTPUT_DIR = _default_output_dir(CELL_TYPE, RESOLUTION)

    # Compute cluster key spec from RESOLUTION
    cluster_key_spec = "default" if RESOLUTION == "default" else RESOLUTION

    # Display run mode
    print(f"\n🔄 Run Configuration")
    print("="*70)
    print(f"Cell Type:    {CELL_TYPE}")
    print(f"Resolution:   {RESOLUTION}")
    print(f"Run Mode:     {RUN_MODE.upper()}")
    print(f"Input:        {INPUT_H5AD_PATH}")
    print(f"Output:       {OUTPUT_DIR}")
    print(f"Progress Bar: {'Enabled' if USE_PROGRESS_BARS else 'Disabled'}")
    print("="*70)

    if RUN_MODE == "pre_annotation":
        print("\n📋 Mode: PRE-ANNOTATION")
        print("   Will perform:")
        print("   ✓ Load data")
        print("   ✓ Check markers")
        print("   ✓ Generate DotPlot")
        if RUN_FIND_MARKERS:
            print("   ✓ Find cluster markers")
            print("   ✓ Generate marker heatmaps")
        print("   ✓ Create Annotation_template.csv")
        print("   ✗ Skip annotation application")
        print("   ✗ Skip cell type UMAP generation")
        print("\n   📝 Next step: Create Annotation.csv from template")

    elif RUN_MODE == "post_annotation":
        print("\n🏷️  Mode: POST-ANNOTATION")
        print("   Will perform:")
        print("   ✓ Load data")
        print("   ✓ Apply annotation from Annotation.csv")
        print("   ✓ Generate cell type UMAPs")
        print("   ✓ Save annotated data")
        print("   ✗ Skip marker analysis")
        print("\n   📋 Requirement: Annotation.csv must exist")

    else:  # full
        print("\n🔄 Mode: FULL")
        print("   Will perform all steps:")
        print("   ✓ Complete pre-annotation workflow")
        print("   ✓ Complete post-annotation workflow")
        print("\n   ⚠️  Note: Requires Annotation.csv for post-annotation steps")

    print("="*70)

    start_time = time.time()

    # Environment check
    log_msg("\n📦 Checking Python environment...")
    try:
        import scanpy as sc
        sc_version = sc.__version__
        log_msg(f"   ✓ scanpy: {sc_version}")
        
        # Version check
        if HAS_PACKAGING and pkg_version.parse(sc_version) < pkg_version.parse("1.9.0"):
            log_msg(f"   ⚠️  Scanpy version < 1.9.0 may have compatibility issues")
    except ImportError as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        print("\n Installation: pip install scanpy", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_msg(f"\n📂 Output directory: {output_dir}")

    # Resolve markers for current CELL_TYPE
    MARKERS = MARKER_CONFIGS[CELL_TYPE]["markers"]

    # Load data and resolve cluster key
    adata = None
    try:
        adata, resolved_cluster_key = load_integrated_data(INPUT_H5AD_PATH, cluster_key_spec)
    except Exception as e:
        print(f"\n❌ Failed to load/resolve data: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # ==================== PRE-ANNOTATION WORKFLOW ====================
    if RUN_MODE in ["pre_annotation", "full"]:
        try:
            # Check markers
            available_markers = check_marker_availability(adata, MARKERS, CELL_TYPE)

            # DotPlot
            if GENERATE_DOTPLOT:
                generate_dotplot(adata, available_markers, resolved_cluster_key, output_dir)

            # Find all markers (DE)
            marker_df = None
            if RUN_FIND_MARKERS:
                marker_df = find_all_markers(adata, resolved_cluster_key, output_dir)

            # Heatmaps
            if GENERATE_HEATMAP and marker_df is not None:
                generate_marker_heatmaps(adata, resolved_cluster_key, marker_df, output_dir)

            # Create template
            create_annotation_template(adata, resolved_cluster_key, output_dir)

            if RUN_MODE == "pre_annotation":
                processing_time = time.time() - start_time
                generate_summary_report(adata, resolved_cluster_key, MARKERS, output_dir, 
                                       RUN_MODE, processing_time)

                print("\n" + "="*70)
                print("✅ PRE-ANNOTATION COMPLETED")
                print("="*70)
                print("\n📝 Next steps:")
                print(f"1. Review results in: {output_dir}/figures/")
                print(f"2. Copy Annotation_template.csv to Annotation.csv:")
                print(f"   cp {output_dir}/Annotation_template.csv {output_dir}/Annotation.csv")
                print(f"3. Edit Annotation.csv with cell type assignments")
                print(f"4. Run again with --run-mode post_annotation")
                print("="*70)
                return
        
        except Exception as e:
            print(f"\n❌ Pre-annotation workflow failed: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            sys.exit(1)

    # ==================== POST-ANNOTATION WORKFLOW ====================
    if RUN_MODE in ["post_annotation", "full"]:
        try:
            # In full mode, require Annotation.csv to proceed to post steps
            if RUN_MODE == "full":
                annotation_path = Path(output_dir) / ANNOTATION_FILE
                if not annotation_path.exists():
                    print("\n" + "="*70)
                    print("⚠️  ANNOTATION FILE NOT FOUND")
                    print("="*70)
                    print(f"\nCannot proceed with post-annotation workflow.")
                    print(f"Annotation file not found: {annotation_path}")
                    print("\nPre-annotation completed. Please:")
                    print(f"1. Copy template: cp {output_dir}/Annotation_template.csv {annotation_path}")
                    print("2. Edit Annotation.csv with cell type assignments")
                    print("3. Run script again")
                    print("="*70)

                    # Generate pre-annotation report
                    processing_time = time.time() - start_time
                    generate_summary_report(adata, resolved_cluster_key, MARKERS, output_dir, 
                                           "pre_annotation", processing_time)
                    return

            # Apply annotation
            adata = apply_annotation(adata, resolved_cluster_key, ANNOTATION_FILE, output_dir)

            # Visualizations
            generate_celltype_visualizations(adata, resolved_cluster_key, output_dir)

            # Save annotated data
            save_annotated_data(adata, output_dir)
        
        except Exception as e:
            print(f"\n❌ Post-annotation workflow failed: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            sys.exit(1)

    # ==================== FINAL REPORT ====================
    processing_time = time.time() - start_time
    generate_summary_report(adata, resolved_cluster_key, MARKERS, output_dir, RUN_MODE, processing_time)

    print("\n" + "="*70)
    print(f"✅ Pipeline completed in {processing_time:.1f} seconds ({processing_time/60:.1f} minutes)")
    print("="*70)
    print(f"\n📂 All results saved to: {output_dir}")
    print("="*70)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️  Pipeline interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n\n❌ Fatal error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)