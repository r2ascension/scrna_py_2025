# %% [markdown]
# # T/NK Cell Annotation Validation: Dotplot & Feature Plot
#
# **Purpose**: Validate L2/L3/L4 T cell annotations using marker gene dotplots and UMAP feature plots
# **Data**: `adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad`
# **Date**: 2026-02-19
#
# **Annotation Hierarchy**:
# - L2: NK cells / CD4 T cells / CD8 T cells
# - L3: Functional subtypes (TRM, TEMRA, Naive, Memory, Th1, Th17, Tfh, Treg, etc.)
# - L4: Descriptive subcluster labels (marker-based)
# - L4_pipeline: Original pipeline cluster names (key column, read-only)
#

# %% [markdown]
# ## 0. Configuration

# %%
# ============================================================================
# CONFIGURATION
# ============================================================================

INPUT_H5AD  = "/home/h2048/data/py/0129/tnk_analysis_unified/results/subcluster_unified_v2_20260129/adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"
OUTPUT_DIR  = "/home/h2048/data/py/0212/tcell_annotation_viz"

# Column in adata.obs that holds the fine-grained pipeline cluster names
# (e.g. 'tcm_naive_helper_t_cells_c0') -- read-only source key
CLUSTER_KEY = 'cell_type_L3'

# Batch key (used for neighbours fallback if UMAP is missing)
BATCH_KEY   = 'dataset'

DPI = 300

# USE_RAW is the intent; USE_RAW_EFFECTIVE is set after load (True only when .raw actually exists)
USE_RAW = True

print("Configuration loaded")
print(f"  Input : {INPUT_H5AD}")
print(f"  Output: {OUTPUT_DIR}")

# %% [markdown]
# ## 1. Imports & Setup

# %%
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import gc

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=DPI, facecolor='white', frameon=False)

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
fig_dir = output_dir / 'figures'
fig_dir.mkdir(exist_ok=True)

# Set figdir early so stacked_violin save= writes to the right place
sc.settings.figdir = str(fig_dir)

print(f"scanpy : {sc.__version__}")
print(f"Output : {output_dir}")
print(f"figdir : {sc.settings.figdir}")

# %% [markdown]
# ## 2. Marker Gene Definitions
#
# **Critical rule**: every key in `NK_L3_MARKERS` / `CD4_L3_MARKERS` / `CD8_L3_MARKERS`
# must be **identical** to the corresponding value in `L4_TO_L3` (and appear in `*_L3_ORDER`).
# Any mismatch causes `make_dotplot()` to silently return empty gene lists.

# %%
# ============================================================================
# L3-LEVEL MARKERS
# Keys must exactly match L4_TO_L3 values.
# ============================================================================

# --- NK Cells ---
# Notes:
#   'Inflammatory NK cells': S100A8/S100A9 removed (myeloid contamination signal);
#                            replaced with ISG/IFN-response markers
#   'Resting NK cells':      TRGC2 removed (TCR gamma-constant, belongs to T/gammadelta);
#                            replaced with KLRD1/LTB/MAL
NK_L3_MARKERS = {
    'Exhausted NK cells':    ['TOX', 'HAVCR2', 'TIGIT', 'LAG3', 'PDCD1', 'NCAM1', 'NCR1', 'CXCR6', 'GZMK', 'KLRC1'],
    'Cytotoxic NK cells':    ['FCGR3A', 'NKG7', 'GNLY', 'PRF1', 'FGFBP2', 'GZMB', 'GZMH', 'FCER1G', 'TYROBP', 'KLRD1'],
    'Resting NK cells':      ['KLRC1', 'KLRD1', 'CXCR4', 'LTB', 'MAL'],
    'Inflammatory NK cells': ['IFITM3', 'ISG15', 'IFI6', 'IFIT1', 'MX1', 'STAT1', 'HLA-DRA', 'HLA-DRB1'],
    'ILC3 NK cells':         ['KIT', 'IL7R', 'IL1R1', 'IL18R1', 'RORA', 'AHR', 'CCR6', 'IL23R', 'ZBTB16'],
}

# --- CD4 T Cells ---
CD4_L3_MARKERS = {
    'Naive CD4 T cells':      ['CCR7', 'SELL', 'TCF7', 'LEF1', 'LTB', 'MAL', 'IL7R', 'FCMR', 'CD27', 'BACH2'],
    'Memory CD4 T cells':     ['IL7R', 'LTB', 'CCR6', 'S100A4', 'IL32', 'RGS1', 'CXCR4', 'LGALS3', 'ICOS'],
    'Th1 CD4 T cells':        ['IFNG', 'TBX21', 'CXCR3', 'GZMK', 'TNF', 'NKG7', 'PHLDA1', 'GADD45G'],
    'Th17 CD4 T cells':       ['CCR6', 'IL23R', 'RORC', 'RORA', 'KLRB1', 'IL4I1', 'RUNX2'],
    'Tfh-like CD4 T cells':   ['CXCR5', 'PDCD1', 'ICOS', 'SH2D1A', 'CXCL13', 'TOX2', 'MAF', 'TIGIT'],
    'Regulatory CD4 T cells': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TIGIT', 'TNFRSF18', 'BATF'],
}

# --- CD8 T Cells ---
# Memory CD8 now includes CCR7/SELL/LEF1 as 'naive/TCM catch':
# if a cluster named 'TEM' but showing high CCR7/SELL -> reclassify L3 to Naive/TCM
CD8_L3_MARKERS = {
    'Memory CD8 T cells':        ['IL7R', 'TCF7', 'LEF1', 'CCR7', 'SELL', 'LTB', 'CD27', 'MAL', 'TOB1'],
    'TEMRA CD8 T cells':         ['KLRG1', 'CX3CR1', 'FGFBP2', 'NKG7', 'GNLY', 'PRF1', 'GZMB', 'GZMH', 'S1PR5', 'ZEB2', 'ADGRG1'],
    'TRM CD8 T cells':           ['ITGA1', 'ITGAE', 'CXCR6', 'ZNF683', 'RUNX3', 'CD69', 'RGS1', 'GZMK', 'KLRC1'],
    'Gamma-delta CD8 T cells':   ['TRDC', 'TRGC1', 'TRGC2', 'TRDV1', 'TRGV9', 'NKG7', 'ZNF683', 'ITGAE'],
    'MAIT CD8 T cells':          ['KLRB1', 'ZBTB16', 'SLC4A10', 'IL18R1', 'IL18RAP', 'CCR6', 'NKG7', 'TRAC'],
}

# ============================================================================
# L4-LEVEL MARKERS (per pipeline cluster)
# Contamination clusters are excluded here; handled in Section 8.
# ============================================================================

L4_MARKERS = {
    # NK -------------------------------------------------------------------
    'cd16-_nk_cells_c0':           ['TOX', 'NCAM1', 'NCR1', 'HAVCR2', 'TIGIT', 'GZMK', 'KLRC1'],
    'cd16plus_nk_cells_c0':        ['FCGR3A', 'NKG7', 'GNLY', 'PRF1', 'FGFBP2', 'GZMB', 'S1PR5'],
    'cd16plus_nk_cells_c1':        ['CX3CR1', 'FCGR3A', 'FGFBP2', 'PRF1', 'GZMB', 'TBX21', 'KLRD1'],
    'cd16plus_nk_cells_c2':        ['AREG', 'XCL1', 'XCL2', 'CCL3', 'CCL4', 'GZMH', 'NKG7'],
    'nk_cells_c0':                 ['IFITM3', 'ISG15', 'HLA-DRA', 'HLA-DRB1', 'STAT1', 'MX1', 'IFI6'],
    'nk_cells_c1':                 ['CXCR4', 'KLRC1', 'KLRD1', 'LTB', 'MAL'],
    'ilc3_c0':                     ['KIT', 'IL7R', 'IL1R1', 'IL18R1', 'RORA', 'AHR'],
    # CD4 ------------------------------------------------------------------
    'tcm_naive_helper_t_cells_c0':            ['CCR7', 'SELL', 'TCF7', 'LEF1', 'LTB', 'IL7R', 'CD27'],
    'cd4plus_tem_effector_helper_t_cells_c0': ['CCR7', 'SELL', 'MAL', 'LEF1', 'BACH2', 'TCF7', 'IL7R'],
    'tem_effector_helper_t_cells_c0':         ['IL7R', 'S100A4', 'IL32', 'LGALS3', 'RGS1', 'TNFAIP3'],
    'tem_effector_helper_t_cells_c1':         ['ICOS', 'CD40LG', 'TNFSF8', 'CCR6', 'IL7R', 'RGS1'],
    'type_1_helper_t_cells_c0':               ['IFNG', 'TBX21', 'CXCR3', 'TNF', 'GZMK', 'PHLDA1', 'GADD45G'],
    'type_17_helper_t_cells_c0':              ['CCR6', 'IL23R', 'RORC', 'RORA', 'KLRB1', 'IL4I1'],
    'follicular_helper_t_cells_c0':           ['CXCR5', 'TCF7', 'IL6R', 'SH2D1A', 'ICOS', 'CD27'],
    'follicular_helper_t_cells_c1':           ['CXCR5', 'PDCD1', 'TOX2', 'ICOS', 'TIGIT', 'CXCL13'],
    'regulatory_t_cells_c0':                  ['CXCL13', 'CXCR5', 'PDCD1', 'TOX2', 'TIGIT', 'ICOS'],
    'cd4plus_regulatory_t_cells_c0':          ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TIGIT', 'TNFRSF18'],
    'regulatory_t_cells_c1':                  ['FOXP3', 'IL2RA', 'CTLA4', 'BATF', 'TIGIT', 'TNFRSF18'],
    'regulatory_t_cells_c2':                  ['FOXP3', 'CTLA4', 'IL2RA', 'PDCD1', 'TOX2', 'CXCR5'],
    # CD8 ------------------------------------------------------------------
    'cd8plus_tem_effector_helper_t_cells_c0': ['IL7R', 'TCF7', 'LEF1', 'CCR7', 'SELL', 'TOB1', 'LTB'],
    'cd8plus_tem_temra_cytotoxic_t_cells_c0': ['KLRG1', 'FGFBP2', 'GZMH', 'PRF1', 'NKG7', 'CX3CR1', 'ZEB2'],
    'cd8plus_tem_temra_cytotoxic_t_cells_c1': ['CX3CR1', 'KLRG1', 'FGFBP2', 'ADGRG1', 'PRF1', 'GNLY', 'S1PR5'],
    'tem_temra_cytotoxic_t_cells_c0':         ['ZEB2', 'KLRG1', 'FGFBP2', 'GZMH', 'PRF1', 'NKG7', 'S1PR5'],
    'tem_temra_cytotoxic_t_cells_c1':         ['CX3CR1', 'S1PR5', 'KLRG1', 'FGFBP2', 'GNLY', 'PRF1', 'NKG7'],
    'cd8plus_tem_trm_cytotoxic_t_cells_c0':   ['ITGAE', 'GZMK', 'IFNG', 'CXCR6', 'CD69', 'ZNF683', 'RUNX3'],
    'cd8plus_tem_trm_cytotoxic_t_cells_c1':   ['HLA-DPA1', 'HLA-DPB1', 'HLA-DRB1', 'CD74', 'IFITM1', 'ISG15'],
    'tem_trm_cytotoxic_t_cells_c0':           ['CD69', 'RGS1', 'ITGA1', 'CXCR6', 'ZNF683', 'RUNX3', 'GZMK'],
    'tem_trm_cytotoxic_t_cells_c1':           ['KLRG1', 'CD27', 'LTB', 'TRAT1', 'SH2D1A', 'IL7R'],
    'trm_cytotoxic_t_cells_c0':               ['ITGA1', 'CXCR6', 'CD69', 'RGS1', 'ZNF683', 'RUNX3', 'GZMK'],
    'trm_cytotoxic_t_cells_c2':               ['MT1E', 'MT1X', 'SPRY1', 'CXCR6', 'ITGA1', 'ZNF683'],
    'cd8plus_trm_cytotoxic_t_cells_c0':       ['ITGA1', 'CXCR6', 'RUNX3', 'JAML', 'KLRC1', 'ZNF683'],
    'gamma-delta_t_cells_c0':                 ['TRDC', 'TRGC1', 'TRGV9', 'TRDV1', 'NKG7', 'ZNF683'],
    'cd8plus_gamma-delta_t_cells_c0':         ['TRDC', 'TRGV9', 'TRDV1', 'ITGAE', 'ZNF683', 'CXCL13'],
    'mait_cells_c0':                          ['KLRB1', 'ZBTB16', 'IL18R1', 'SLC4A10', 'CCR6', 'NKG7'],
}

# Clusters with strong epithelial contamination -- handled separately in Section 8
CONTAMINATION_CLUSTERS = [
    'trm_cytotoxic_t_cells_c1',
    'cd8plus_trm_cytotoxic_t_cells_c1',
    'cd8plus_trm_cytotoxic_t_cells_c2',
    'mait_cells_c1',
]
CONTAMINATION_T_MARKERS   = ['TRAC', 'CD3D', 'CD8A', 'LTB', 'ITGAE', 'KLRB1', 'ZBTB16']
CONTAMINATION_EPI_MARKERS = ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'SCGB1A1',
                              'WFDC2', 'SLPI', 'BPIFB1', 'BPIFA1', 'STATH',
                              'LTF', 'DMBT1', 'MSMB', 'TG']

# ============================================================================
# L4 -> L3 -> L2 Mappings
# ============================================================================

L4_TO_L3 = {
    'cd16-_nk_cells_c0':                       'Exhausted NK cells',
    'cd16plus_nk_cells_c0':                    'Cytotoxic NK cells',
    'cd16plus_nk_cells_c1':                    'Cytotoxic NK cells',
    'cd16plus_nk_cells_c2':                    'Cytotoxic NK cells',
    'nk_cells_c0':                             'Inflammatory NK cells',
    'nk_cells_c1':                             'Resting NK cells',
    'ilc3_c0':                                 'ILC3 NK cells',
    'tcm_naive_helper_t_cells_c0':             'Naive CD4 T cells',
    'cd4plus_tem_effector_helper_t_cells_c0':  'Naive CD4 T cells',
    'tem_effector_helper_t_cells_c0':          'Memory CD4 T cells',
    'tem_effector_helper_t_cells_c1':          'Memory CD4 T cells',
    'type_1_helper_t_cells_c0':                'Th1 CD4 T cells',
    'type_17_helper_t_cells_c0':               'Th17 CD4 T cells',
    'follicular_helper_t_cells_c0':            'Tfh-like CD4 T cells',
    'follicular_helper_t_cells_c1':            'Tfh-like CD4 T cells',
    'regulatory_t_cells_c0':                   'Tfh-like CD4 T cells',
    'cd4plus_regulatory_t_cells_c0':           'Regulatory CD4 T cells',
    'regulatory_t_cells_c1':                   'Regulatory CD4 T cells',
    'regulatory_t_cells_c2':                   'Regulatory CD4 T cells',
    'cd8plus_tem_effector_helper_t_cells_c0':  'Memory CD8 T cells',
    'cd8plus_tem_temra_cytotoxic_t_cells_c0':  'TEMRA CD8 T cells',
    'cd8plus_tem_temra_cytotoxic_t_cells_c1':  'TEMRA CD8 T cells',
    'tem_temra_cytotoxic_t_cells_c0':          'TEMRA CD8 T cells',
    'tem_temra_cytotoxic_t_cells_c1':          'TEMRA CD8 T cells',
    'cd8plus_tem_trm_cytotoxic_t_cells_c0':    'TRM CD8 T cells',
    'cd8plus_tem_trm_cytotoxic_t_cells_c1':    'TRM CD8 T cells',
    'tem_trm_cytotoxic_t_cells_c0':            'TRM CD8 T cells',
    'tem_trm_cytotoxic_t_cells_c1':            'TRM CD8 T cells',
    'trm_cytotoxic_t_cells_c0':                'TRM CD8 T cells',
    'trm_cytotoxic_t_cells_c1':                'TRM CD8 T cells',
    'trm_cytotoxic_t_cells_c2':                'TRM CD8 T cells',
    'cd8plus_trm_cytotoxic_t_cells_c0':        'TRM CD8 T cells',
    'cd8plus_trm_cytotoxic_t_cells_c1':        'TRM CD8 T cells',
    'cd8plus_trm_cytotoxic_t_cells_c2':        'TRM CD8 T cells',
    'gamma-delta_t_cells_c0':                  'Gamma-delta CD8 T cells',
    'cd8plus_gamma-delta_t_cells_c0':          'Gamma-delta CD8 T cells',
    'mait_cells_c0':                           'MAIT CD8 T cells',
    'mait_cells_c1':                           'MAIT CD8 T cells',
}

L3_TO_L2 = {
    'Exhausted NK cells':     'NK cells',
    'Cytotoxic NK cells':     'NK cells',
    'Inflammatory NK cells':  'NK cells',
    'Resting NK cells':       'NK cells',
    'ILC3 NK cells':          'NK cells',
    'Naive CD4 T cells':      'CD4 T cells',
    'Memory CD4 T cells':     'CD4 T cells',
    'Th1 CD4 T cells':        'CD4 T cells',
    'Th17 CD4 T cells':       'CD4 T cells',
    'Tfh-like CD4 T cells':   'CD4 T cells',
    'Regulatory CD4 T cells': 'CD4 T cells',
    'Memory CD8 T cells':     'CD8 T cells',
    'TEMRA CD8 T cells':      'CD8 T cells',
    'TRM CD8 T cells':        'CD8 T cells',
    'Gamma-delta CD8 T cells':'CD8 T cells',
    'MAIT CD8 T cells':       'CD8 T cells',
}

# ============================================================================
# L4 Descriptive Label Mapping (pipeline cluster name -> publication label)
#
# Note on nk_cells_c0: S100A8/A9 removed from label (myeloid contamination
# signal); ISG/IFN-response signature used instead for consistency with
# NK_L3_MARKERS definition above.
#
# Note on contamination clusters (trm_c1, cd8plus_trm_c1/c2, mait_c1):
# retained in mapping with 'Epithelial-associated/like' labels; downstream
# filtering/removal decision left to Section 8 QC.
# ============================================================================

L4_TO_L4_LABEL = {
    # NK
    'cd16-_nk_cells_c0':                       'GZMK+ Exhausted NK CD16- NK cells',
    'cd16plus_nk_cells_c0':                    'HAVCR2+ Cytotoxic NK CD16+ NK cells',
    'cd16plus_nk_cells_c1':                    'CX3CR1+ Cytotoxic NK CD16+ NK cells',
    'cd16plus_nk_cells_c2':                    'AREG+ Cytokine-producing NK CD16+ NK cells',
    'nk_cells_c0':                             'ISG/IFN-response HLA-II+ Inflammatory NK cells',
    'nk_cells_c1':                             'CXCR4+ Resting NK NK cells',
    'ilc3_c0':                                 'KIT+ IL7R+ ILC3 NK ILC cells',
    # CD4
    'tcm_naive_helper_t_cells_c0':             'CCR7+ SELL+ Naive-like CD4 TCM/Naive T cells',
    'cd4plus_tem_effector_helper_t_cells_c0':  'MAL+ LEF1+ Naive-like CD4 TCM T cells',
    'tem_effector_helper_t_cells_c0':          'SLC2A3+ Inflammatory CD4 TEM T cells',
    'tem_effector_helper_t_cells_c1':          'ICOS+ CD40LG+ Activated CD4 TEM T cells',
    'type_1_helper_t_cells_c0':                'IFNG+ TNF+ Activated CD4 Th1 T cells',
    'type_17_helper_t_cells_c0':               'IL23R+ CCR6+ Activated CD4 Th17 T cells',
    'follicular_helper_t_cells_c0':            'TCF7+ Resting CD4 Tfh T cells',
    'follicular_helper_t_cells_c1':            'TIGIT+ PD-1+ CD4 Tfh T cells',
    'regulatory_t_cells_c0':                   'PD-1+ CXCL13+ CD4 Follicular regulatory T cells',
    'cd4plus_regulatory_t_cells_c0':           'CXCR5+ CTLA4+ PD-1+ CD4 Follicular regulatory T cells',
    'regulatory_t_cells_c1':                   'IL2RA+ ICOS+ Activated CD4 Regulatory T cells',
    'regulatory_t_cells_c2':                   'CXCR5+ CTLA4+ PD-1+ CD4 Follicular regulatory T cells',
    # CD8
    'cd8plus_tem_effector_helper_t_cells_c0':  'IL7R+ TCF7+ Naive-like CD8 TCM T cells',
    'cd8plus_tem_temra_cytotoxic_t_cells_c0':  'MHC-II+ Cytotoxic CD8 TEMRA T cells',
    'cd8plus_tem_temra_cytotoxic_t_cells_c1':  'CX3CR1+ Terminal-effector Cytotoxic CD8 TEMRA T cells',
    'tem_temra_cytotoxic_t_cells_c0':          'ZEB2+ Terminal-effector Cytotoxic CD8 TEMRA T cells',
    'tem_temra_cytotoxic_t_cells_c1':          'CX3CR1+ Terminal-effector Cytotoxic CD8 TEMRA T cells',
    'cd8plus_tem_trm_cytotoxic_t_cells_c0':    'Exhausted Cytotoxic CD8 TRM T cells',
    'cd8plus_tem_trm_cytotoxic_t_cells_c1':    'MHC-II+ Activated Cytotoxic CD8 TRM T cells',
    'tem_trm_cytotoxic_t_cells_c0':            'MHC-II+ Activated Cytotoxic CD8 TRM T cells',
    'tem_trm_cytotoxic_t_cells_c1':            'IFNG+ Metabolic-activated Cytotoxic CD8 TRM T cells',
    'cd8plus_trm_cytotoxic_t_cells_c0':        'Resting Cytotoxic CD8 TRM T cells',
    'cd8plus_trm_cytotoxic_t_cells_c1':        'Epithelial-associated Cytotoxic CD8 TRM T cells',
    'cd8plus_trm_cytotoxic_t_cells_c2':        'Epithelial-like Cytotoxic CD8 TRM T cells',
    'trm_cytotoxic_t_cells_c0':                'Resting Cytotoxic CD8 TRM T cells',
    'trm_cytotoxic_t_cells_c1':                'CXCR4+ Activated Cytotoxic CD8 TRM T cells',
    'trm_cytotoxic_t_cells_c2':                'Metallothionein-high Cytotoxic CD8 TRM T cells',
    'cd8plus_gamma-delta_t_cells_c0':          'CXCL13+ ENTPD1+ TRM-like CD8 gammadelta T cells',
    'gamma-delta_t_cells_c0':                  'KLRC1+ TRDV1+ TRM-like CD8 gammadelta T cells',
    'mait_cells_c0':                           'IL23R+ CCR6+ Type-17-like CD8 MAIT T cells',
    'mait_cells_c1':                           'Epithelial-associated Type-17-like CD8 MAIT T cells',
}

NK_L3_ORDER  = ['Exhausted NK cells', 'Resting NK cells', 'ILC3 NK cells',
                'Inflammatory NK cells', 'Cytotoxic NK cells']
CD4_L3_ORDER = ['Naive CD4 T cells', 'Memory CD4 T cells', 'Th1 CD4 T cells',
                'Th17 CD4 T cells', 'Tfh-like CD4 T cells', 'Regulatory CD4 T cells']
CD8_L3_ORDER = ['Memory CD8 T cells', 'TEMRA CD8 T cells', 'TRM CD8 T cells',
                'Gamma-delta CD8 T cells', 'MAIT CD8 T cells']

print("Marker dictionaries defined")
print(f"  NK L3  : {list(NK_L3_MARKERS.keys())}")
print(f"  CD4 L3 : {list(CD4_L3_MARKERS.keys())}")
print(f"  CD8 L3 : {list(CD8_L3_MARKERS.keys())}")
print(f"  L4     : {len(L4_MARKERS)} clusters")
print(f"  L4 labels: {len(L4_TO_L4_LABEL)} entries")

# %% [markdown]
# ## 3. Load Data

# %%
print("Loading h5ad...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"  Shape : {adata.shape}")
print(f"  obs   : {list(adata.obs.columns)}")
print(f"  obsm  : {list(adata.obsm.keys())}")
print(f"  layers: {list(adata.layers.keys()) if adata.layers else 'None'}")
print(f"  .raw  : {adata.raw is not None}")

# -----------------------------------------------------------------------
# USE_RAW_EFFECTIVE: set True only when .raw actually exists
# All plotting calls must use USE_RAW_EFFECTIVE, not the raw USE_RAW flag
# -----------------------------------------------------------------------
USE_RAW_EFFECTIVE = USE_RAW and (adata.raw is not None)
print(f"\nUSE_RAW_EFFECTIVE = {USE_RAW_EFFECTIVE}")
if USE_RAW and not adata.raw:
    print("  WARNING: .raw is None -- using adata.X for all visualizations")

# -----------------------------------------------------------------------
# UMAP fallback: check neighbours graph before calling sc.tl.umap()
# -----------------------------------------------------------------------
if 'X_umap' not in adata.obsm:
    print("\nX_umap not found, computing...")
    has_neighbors = ('neighbors' in adata.uns) or ('connectivities' in adata.obsp)
    if not has_neighbors:
        print("  No neighbours graph -- building now")
        if 'X_scanvi' in adata.obsm:
            sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=30)
        elif 'X_scvi' in adata.obsm:
            sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=30)
        else:
            sc.pp.pca(adata)
            sc.pp.neighbors(adata, n_pcs=50, n_neighbors=30)
    sc.tl.umap(adata)
    print("  UMAP computed")

# %%
# ============================================================================
# Apply L4 / L3 / L2 annotations
# New columns: cell_type_L4_new, cell_type_L3_new, cell_type_L2_new
# Original CLUSTER_KEY column ('cell_type_L3') is never modified.
# ============================================================================

print(f"All unique values in '{CLUSTER_KEY}':")
current_clusters = sorted(adata.obs[CLUSTER_KEY].dropna().unique().tolist())
print(f"  {len(current_clusters)} clusters")
for c in current_clusters:
    print(f"    '{c}'")

# Derive all three levels from the read-only pipeline cluster key
adata.obs['cell_type_L4_new'] = adata.obs[CLUSTER_KEY].map(L4_TO_L4_LABEL)
adata.obs['cell_type_L3_new'] = adata.obs[CLUSTER_KEY].map(L4_TO_L3)
adata.obs['cell_type_L2_new'] = adata.obs['cell_type_L3_new'].map(L3_TO_L2)

n_mapped_l4 = adata.obs['cell_type_L4_new'].notna().sum()
n_mapped_l3 = adata.obs['cell_type_L3_new'].notna().sum()
n_mapped_l2 = adata.obs['cell_type_L2_new'].notna().sum()
print(f"\nMapped L4: {n_mapped_l4:,}/{adata.n_obs:,} ({n_mapped_l4/adata.n_obs*100:.1f}%)")
print(f"Mapped L3: {n_mapped_l3:,}/{adata.n_obs:,} ({n_mapped_l3/adata.n_obs*100:.1f}%)")
print(f"Mapped L2: {n_mapped_l2:,}/{adata.n_obs:,} ({n_mapped_l2/adata.n_obs*100:.1f}%)")

# Report any pipeline clusters missing from L4_TO_L4_LABEL or L4_TO_L3
unmapped_l4 = adata.obs.loc[adata.obs['cell_type_L4_new'].isna(), CLUSTER_KEY].value_counts()
unmapped_l3 = adata.obs.loc[adata.obs['cell_type_L3_new'].isna(), CLUSTER_KEY].value_counts()
if len(unmapped_l4):
    print("\nUnmapped L4 clusters (add to L4_TO_L4_LABEL):")
    for k, v in unmapped_l4.items():
        print(f"  '{k}': {v} cells")
if len(unmapped_l3):
    print("\nUnmapped L3 clusters (add to L4_TO_L3):")
    for k, v in unmapped_l3.items():
        print(f"  '{k}': {v} cells")

print("\nL2 cell counts:"); print(adata.obs['cell_type_L2_new'].value_counts())
print("\nL3 cell counts:"); print(adata.obs['cell_type_L3_new'].value_counts())
print("\nL4 cell counts:"); print(adata.obs['cell_type_L4_new'].value_counts())

# %%
# ============================================================================
# Gene availability check
# ============================================================================

def _avail(adata, use_raw_eff):
    if use_raw_eff and adata.raw is not None:
        return set(adata.raw.var_names), '.raw'
    return set(adata.var_names), 'adata.var'

def check_genes(gene_dict, adata, use_raw_eff):
    avail, src = _avail(adata, use_raw_eff)
    missing_all, valid = [], {}
    for grp, genes in gene_dict.items():
        valid[grp] = [g for g in genes if g in avail]
        missing_all += [g for g in genes if g not in avail]
    if missing_all:
        print(f"  Missing from {src}: {sorted(set(missing_all))}")
    else:
        print(f"  All found in {src}")
    return valid

def check_flat(genes, adata, use_raw_eff):
    avail, src = _avail(adata, use_raw_eff)
    valid   = [g for g in genes if g in avail]
    missing = [g for g in genes if g not in avail]
    if missing:
        print(f"  Missing from {src}: {missing}")
    return valid

print("NK L3...");  NK_L3_VALID  = check_genes(NK_L3_MARKERS,  adata, USE_RAW_EFFECTIVE)
print("CD4 L3..."); CD4_L3_VALID = check_genes(CD4_L3_MARKERS, adata, USE_RAW_EFFECTIVE)
print("CD8 L3..."); CD8_L3_VALID = check_genes(CD8_L3_MARKERS, adata, USE_RAW_EFFECTIVE)
print("L4...");     L4_VALID     = check_genes(L4_MARKERS,     adata, USE_RAW_EFFECTIVE)

# %% [markdown]
# ## 4. UMAP Overview

# %%
fig, axes = plt.subplots(1, 4, figsize=(32, 7))

sc.pl.umap(adata, color='cell_type_L2_new', title='L2: Major Lineage',
           ax=axes[0], show=False, legend_loc='on data', legend_fontsize=7)
sc.pl.umap(adata, color='cell_type_L3_new', title='L3: Functional Subtype',
           ax=axes[1], show=False, legend_loc='right margin', legend_fontsize=6)
sc.pl.umap(adata, color='cell_type_L4_new', title='L4: Descriptive Label',
           ax=axes[2], show=False, legend_loc='right margin', legend_fontsize=4)
sc.pl.umap(adata, color=CLUSTER_KEY, title='L4 Pipeline: Cluster Key',
           ax=axes[3], show=False, legend_loc='right margin', legend_fontsize=5)

plt.tight_layout()
out = fig_dir / 'umap_L2_L3_L4.pdf'
fig.savefig(out, bbox_inches='tight', dpi=DPI)
plt.show()
print(f"Saved: {out.name}")

# %% [markdown]
# ## 5. Dotplot Helpers

# %%
def make_dotplot(adata, l2_filter, obs_key, marker_dict_valid, order,
                 title, filename, use_raw_eff):
    """
    Dotplot with y-axis = L3 labels and x-axis = markers bracketed by L3 group.
    marker_dict_valid keys must exactly match values in order.
    """
    if l2_filter is not None:
        mask = adata.obs['cell_type_L2_new'] == l2_filter
        adata_sub = adata[mask].copy()
    else:
        adata_sub = adata.copy()
    print(f"  [{l2_filter}] {adata_sub.n_obs:,} cells")

    available = set(adata_sub.obs[obs_key].dropna().unique())
    order_valid = [o for o in order if o in available]
    missing = [o for o in order if o not in available]
    if missing:
        print(f"  L3 labels missing from data: {missing}")

    genes_flat, var_labels, var_positions, pos = [], [], [], 0
    for grp in order_valid:
        genes = marker_dict_valid.get(grp, [])
        if not genes:
            print(f"  No valid genes for '{grp}', skipping bracket")
            continue
        var_labels.append(grp)
        var_positions.append([pos, pos + len(genes) - 1])
        genes_flat.extend(genes)
        pos += len(genes)

    if not genes_flat:
        print("  No genes to plot")
        return

    adata_sub.obs[obs_key] = pd.Categorical(adata_sub.obs[obs_key], categories=order_valid)
    figsize = (max(8, len(genes_flat) * 0.35 + 3),
               max(4, len(order_valid) * 0.55 + 2))

    dp = sc.pl.dotplot(
        adata_sub, var_names=genes_flat, groupby=obs_key,
        categories_order=order_valid, use_raw=use_raw_eff,
        standard_scale='var', color_map='RdBu_r',
        var_group_labels=var_labels, var_group_positions=var_positions,
        title=title, figsize=figsize, return_fig=True, show=False,
    )
    dp.savefig(fig_dir / filename, bbox_inches='tight', dpi=DPI)
    dp.show()
    print(f"  Saved: {filename}")


def make_dotplot_L4(adata, l2_filter, l4_order, l4_dict_valid,
                    title, filename, use_raw_eff):
    """Dotplot with y-axis = raw pipeline cluster names."""
    if l2_filter is not None:
        mask = adata.obs['cell_type_L2_new'] == l2_filter
        adata_sub = adata[mask].copy()
    else:
        adata_sub = adata.copy()
    print(f"  [{l2_filter}] {adata_sub.n_obs:,} cells")

    available = set(adata_sub.obs[CLUSTER_KEY].dropna().unique())
    order_valid = [c for c in l4_order if c in available]
    missing = [c for c in l4_order if c not in available]
    if missing:
        print(f"  Clusters not found: {missing}")

    genes_flat, var_labels, var_positions, pos = [], [], [], 0
    for cluster in order_valid:
        genes = l4_dict_valid.get(cluster, [])
        if not genes:
            continue
        short = (cluster
                 .replace('_cytotoxic_t_cells_', '_')
                 .replace('_helper_t_cells_', '_')
                 .replace('_t_cells_', '_'))
        var_labels.append(short)
        var_positions.append([pos, pos + len(genes) - 1])
        genes_flat.extend(genes)
        pos += len(genes)

    if not genes_flat:
        print("  No genes to plot")
        return

    adata_sub.obs[CLUSTER_KEY] = pd.Categorical(adata_sub.obs[CLUSTER_KEY],
                                                 categories=order_valid)
    figsize = (max(10, len(genes_flat) * 0.3 + 3),
               max(5, len(order_valid) * 0.5 + 2))

    dp = sc.pl.dotplot(
        adata_sub, var_names=genes_flat, groupby=CLUSTER_KEY,
        categories_order=order_valid, use_raw=use_raw_eff,
        standard_scale='var', color_map='RdBu_r',
        var_group_labels=var_labels, var_group_positions=var_positions,
        title=title, figsize=figsize, return_fig=True, show=False,
    )
    dp.savefig(fig_dir / filename, bbox_inches='tight', dpi=DPI)
    dp.show()
    print(f"  Saved: {filename}")


print("Dotplot helpers defined")

# %% [markdown]
# ## 6. Dotplots: L3 Functional Subtypes

# %%
make_dotplot(adata, 'NK cells', 'cell_type_L3_new',
             NK_L3_VALID, NK_L3_ORDER,
             'NK Cells: L3 Marker Genes', 'dotplot_NK_L3.pdf', USE_RAW_EFFECTIVE)

# %%
make_dotplot(adata, 'CD4 T cells', 'cell_type_L3_new',
             CD4_L3_VALID, CD4_L3_ORDER,
             'CD4 T Cells: L3 Marker Genes', 'dotplot_CD4_L3.pdf', USE_RAW_EFFECTIVE)

# %%
make_dotplot(adata, 'CD8 T cells', 'cell_type_L3_new',
             CD8_L3_VALID, CD8_L3_ORDER,
             'CD8 T Cells: L3 Marker Genes', 'dotplot_CD8_L3.pdf', USE_RAW_EFFECTIVE)

# %% [markdown]
# ## 7. Dotplots: L4 Pipeline Clusters

# %%
NK_L4_ORDER = ['cd16-_nk_cells_c0',
               'cd16plus_nk_cells_c0', 'cd16plus_nk_cells_c1', 'cd16plus_nk_cells_c2',
               'nk_cells_c0', 'nk_cells_c1', 'ilc3_c0']
make_dotplot_L4(adata, 'NK cells', NK_L4_ORDER,
                {k: L4_VALID.get(k, []) for k in NK_L4_ORDER},
                'NK Cells: L4 Cluster Markers', 'dotplot_NK_L4.pdf', USE_RAW_EFFECTIVE)

# %%
CD4_L4_ORDER = [
    'tcm_naive_helper_t_cells_c0', 'cd4plus_tem_effector_helper_t_cells_c0',
    'tem_effector_helper_t_cells_c0', 'tem_effector_helper_t_cells_c1',
    'type_1_helper_t_cells_c0', 'type_17_helper_t_cells_c0',
    'follicular_helper_t_cells_c0', 'follicular_helper_t_cells_c1',
    'regulatory_t_cells_c0',
    'cd4plus_regulatory_t_cells_c0', 'regulatory_t_cells_c1', 'regulatory_t_cells_c2',
]
make_dotplot_L4(adata, 'CD4 T cells', CD4_L4_ORDER,
                {k: L4_VALID.get(k, []) for k in CD4_L4_ORDER},
                'CD4 T Cells: L4 Cluster Markers', 'dotplot_CD4_L4.pdf', USE_RAW_EFFECTIVE)

# %%
# Contamination clusters excluded here; see Section 8
CD8_L4_ORDER_CLEAN = [
    'cd8plus_tem_effector_helper_t_cells_c0',
    'cd8plus_tem_temra_cytotoxic_t_cells_c0', 'cd8plus_tem_temra_cytotoxic_t_cells_c1',
    'tem_temra_cytotoxic_t_cells_c0', 'tem_temra_cytotoxic_t_cells_c1',
    'cd8plus_tem_trm_cytotoxic_t_cells_c0', 'cd8plus_tem_trm_cytotoxic_t_cells_c1',
    'tem_trm_cytotoxic_t_cells_c0', 'tem_trm_cytotoxic_t_cells_c1',
    'trm_cytotoxic_t_cells_c0', 'trm_cytotoxic_t_cells_c2',
    'cd8plus_trm_cytotoxic_t_cells_c0',
    'gamma-delta_t_cells_c0', 'cd8plus_gamma-delta_t_cells_c0',
    'mait_cells_c0',
]
make_dotplot_L4(adata, 'CD8 T cells', CD8_L4_ORDER_CLEAN,
                {k: L4_VALID.get(k, []) for k in CD8_L4_ORDER_CLEAN},
                'CD8 T Cells: L4 Cluster Markers (Clean)', 'dotplot_CD8_L4.pdf', USE_RAW_EFFECTIVE)

# %% [markdown]
# ## 8. Contamination / Doublet Check: T + Epithelial Dual Panel
#
# These clusters express both T cell and epithelial markers simultaneously,
# suggesting doublets or ambient RNA.
# - Top rows: clean TRM reference
# - Bottom rows: suspect clusters
# - Left bracket: T cell markers | Right bracket: epithelial markers
#
# Interpretation guide:
# - T high + Epi low  -> normal TRM/MAIT (keep)
# - T high + Epi high -> likely doublet (flag)
# - T low  + Epi high -> likely ambient RNA or misclassified epithelial (remove)

# %%
t_panel   = check_flat(CONTAMINATION_T_MARKERS,   adata, USE_RAW_EFFECTIVE)
epi_panel = check_flat(CONTAMINATION_EPI_MARKERS, adata, USE_RAW_EFFECTIVE)

# Clean reference clusters shown above contamination clusters for visual comparison
reference_clean = ['trm_cytotoxic_t_cells_c0', 'cd8plus_trm_cytotoxic_t_cells_c0', 'mait_cells_c0']
contam_order = [c for c in reference_clean + CONTAMINATION_CLUSTERS
                if c in set(adata.obs[CLUSTER_KEY].dropna().unique())]

mask = adata.obs[CLUSTER_KEY].isin(contam_order)
adata_contam = adata[mask].copy()
adata_contam.obs[CLUSTER_KEY] = pd.Categorical(adata_contam.obs[CLUSTER_KEY],
                                                categories=contam_order)
print(f"{adata_contam.n_obs:,} cells, {len(contam_order)} clusters")

genes_dual = t_panel + epi_panel
t_end  = len(t_panel) - 1
ep_end = len(t_panel) + len(epi_panel) - 1

dp = sc.pl.dotplot(
    adata_contam,
    var_names=genes_dual, groupby=CLUSTER_KEY,
    categories_order=contam_order,
    use_raw=USE_RAW_EFFECTIVE,
    standard_scale='var', color_map='RdBu_r',
    var_group_labels=['T cell markers', 'Epithelial markers'],
    var_group_positions=[[0, t_end], [t_end + 1, ep_end]],
    title='Contamination Check: T vs Epithelial Markers\n(Clean TRM top, suspect clusters bottom)',
    figsize=(max(12, len(genes_dual) * 0.4 + 3),
             max(5, len(contam_order) * 0.6 + 2)),
    return_fig=True, show=False,
)
out = fig_dir / 'dotplot_contamination_check.pdf'
dp.savefig(out, bbox_inches='tight', dpi=DPI)
dp.show()
print(f"Saved: {out.name}")

# %% [markdown]
# ## 9. Cluster Sanity Check Table: Naive / TEMRA / TRM / Epithelial Scores
#
# Computes mean expression score per cluster per panel.
# Sort by `Naive/TCM_mean` to identify clusters with high CCR7/SELL co-expression.
# Any cluster currently labelled TEM/effector but scoring high on Naive/TCM panel
# -> candidate for L3 reclassification.

# %%
SANITY_PANELS = {
    'Naive/TCM':   ['CCR7', 'SELL', 'TCF7', 'LEF1', 'LTB', 'IL7R', 'CD27', 'MAL'],
    'TEMRA':       ['KLRG1', 'CX3CR1', 'FGFBP2', 'PRF1', 'GNLY', 'NKG7', 'S1PR5', 'ZEB2'],
    'TRM':         ['ITGA1', 'ITGAE', 'CXCR6', 'CD69', 'RGS1', 'ZNF683', 'RUNX3', 'GZMK'],
    'Epithelial':  ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'SCGB1A1', 'WFDC2', 'SLPI', 'BPIFA1'],
}

avail_genes, _ = _avail(adata, USE_RAW_EFFECTIVE)

rows = []
for cluster in sorted(adata.obs[CLUSTER_KEY].dropna().unique()):
    mask = adata.obs[CLUSTER_KEY] == cluster
    adata_c = adata[mask]
    row = {'cluster': cluster,
           'n_cells': int(mask.sum()),
           'L4_label':  L4_TO_L4_LABEL.get(cluster, 'UNMAPPED'),
           'L3_current': L4_TO_L3.get(cluster, 'UNMAPPED')}
    for pname, pgenes in SANITY_PANELS.items():
        valid = [g for g in pgenes if g in avail_genes]
        if not valid:
            row[f'{pname}_mean'] = np.nan
            continue
        if USE_RAW_EFFECTIVE and adata.raw is not None:
            raw_idx = [list(adata.raw.var_names).index(g) for g in valid]
            X = adata_c.raw.X[:, raw_idx]
        else:
            gene_idx = [list(adata.var_names).index(g) for g in valid]
            X = adata_c.X[:, gene_idx]
        if sp.issparse(X):
            X = X.toarray()
        row[f'{pname}_mean'] = float(np.mean(X))
    rows.append(row)

sanity_df = pd.DataFrame(rows).sort_values('Naive/TCM_mean', ascending=False)
cols = ['cluster', 'n_cells', 'L4_label', 'L3_current',
        'Naive/TCM_mean', 'TEMRA_mean', 'TRM_mean', 'Epithelial_mean']
print(sanity_df[cols].to_string(index=False))

sanity_df.to_csv(output_dir / 'cluster_sanity_scores.csv', index=False)
print("\nSaved: cluster_sanity_scores.csv")

# Clusters to consider reclassifying to Naive/TCM
NAIVE_THRESHOLD = 0.3   # review actual distribution first
candidates = sanity_df[sanity_df['Naive/TCM_mean'] > NAIVE_THRESHOLD]
if len(candidates):
    print(f"\nCandidates for L3 reclassification to Naive/TCM (score > {NAIVE_THRESHOLD}):")
    print(candidates[['cluster', 'L4_label', 'L3_current', 'Naive/TCM_mean']].to_string(index=False))

# %% [markdown]
# ## 10. Feature Plots

# %%
def plot_feature_grid(adata, genes, title, filename, use_raw_eff, ncols=5):
    avail, _ = _avail(adata, use_raw_eff)
    genes = [g for g in genes if g in avail]
    if not genes:
        print(f"No valid genes for {title}")
        return
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
    axes = axes.flatten()
    for i, gene in enumerate(genes):
        sc.pl.umap(adata, color=gene, ax=axes[i], show=False,
                   use_raw=use_raw_eff, frameon=False,
                   colorbar_loc=None, color_map='viridis', vmin=0)
        axes[i].set_title(gene, fontsize=10, fontweight='bold')
    for i in range(len(genes), len(axes)):
        axes[i].set_visible(False)
    plt.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()
    out = fig_dir / filename
    fig.savefig(out, bbox_inches='tight', dpi=DPI)
    plt.show()
    print(f"Saved: {filename}")

print("plot_feature_grid defined")

# %%
LINEAGE_GENES = [
    'CD3D', 'CD4', 'CD8A', 'CD8B', 'GNLY', 'NKG7', 'FCGR3A', 'NCAM1',
    'CCR7', 'SELL', 'TCF7', 'LEF1',
    'GZMB', 'GZMK', 'GZMH', 'PRF1',
    'ITGA1', 'CXCR6', 'CD69', 'RGS1',
    'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT',
    'FOXP3', 'IL2RA', 'CTLA4',
    'ZBTB16', 'TRDV1', 'IL23R', 'CX3CR1', 'KLRG1',
]
plot_feature_grid(adata, LINEAGE_GENES, 'T/NK Lineage Markers',
                  'featureplot_lineage.pdf', USE_RAW_EFFECTIVE)

# %%
NK_GENES = ['FCGR3A', 'GNLY', 'NKG7', 'PRF1', 'GZMB', 'FGFBP2',
            'TOX', 'NCAM1', 'NCR1', 'HAVCR2',
            'IFITM3', 'ISG15', 'HLA-DRA', 'STAT1',
            'KLRC1', 'KLRD1', 'CXCR4', 'LTB',
            'KIT', 'IL7R', 'IL1R1', 'RORA', 'CX3CR1', 'AREG']
plot_feature_grid(adata, NK_GENES, 'NK Cell Subtype Markers',
                  'featureplot_NK.pdf', USE_RAW_EFFECTIVE)

# %%
CD4_GENES = [
    'CCR7', 'SELL', 'TCF7', 'MAL', 'LEF1', 'BACH2',
    'RGS1', 'LGALS3', 'ICOS', 'CD40LG',
    'IFNG', 'TBX21', 'CXCR3', 'TNF',
    'IL23R', 'CCR6', 'RORC', 'KLRB1',
    'CXCR5', 'TOX2', 'PDCD1', 'TIGIT', 'CXCL13',
    'FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'BATF',
]
plot_feature_grid(adata, CD4_GENES, 'CD4 T Cell Subtype Markers',
                  'featureplot_CD4.pdf', USE_RAW_EFFECTIVE)

# %%
CD8_GENES = [
    'IL7R', 'TCF7', 'CD27', 'CCR7', 'SELL', 'TOB1',
    'KLRG1', 'CX3CR1', 'FGFBP2', 'GZMH', 'S1PR5', 'ZEB2', 'ADGRG1',
    'ITGA1', 'ITGAE', 'CXCR6', 'CD69', 'RGS1', 'ZNF683', 'RUNX3', 'GZMK',
    'TRDC', 'TRDV1', 'TRGV9',
    'KLRB1', 'ZBTB16', 'IL18R1', 'SLC4A10',
]
plot_feature_grid(adata, CD8_GENES, 'CD8 T Cell Subtype Markers',
                  'featureplot_CD8.pdf', USE_RAW_EFFECTIVE)

# %% [markdown]
# ## 11. Stacked Violin

# %%
VIOLIN_GENES = [
    'CD3D', 'CD4', 'CD8A', 'GNLY', 'NKG7', 'FCGR3A',
    'CCR7', 'SELL', 'TCF7', 'LEF1',
    'GZMB', 'GZMK', 'GZMH', 'PRF1',
    'ITGA1', 'CXCR6', 'CD69',
    'PDCD1', 'HAVCR2', 'TIGIT', 'LAG3',
    'FOXP3', 'IL2RA', 'CTLA4',
    'KLRG1', 'CX3CR1',
    'CXCR5', 'IFNG', 'IL23R', 'ZBTB16',
]

avail_v, _ = _avail(adata, USE_RAW_EFFECTIVE)
VIOLIN_GENES = [g for g in VIOLIN_GENES if g in avail_v]

ALL_L3_ORDER = NK_L3_ORDER + CD4_L3_ORDER + CD8_L3_ORDER
avail_L3 = set(adata.obs['cell_type_L3_new'].dropna().unique())
ALL_L3_ORDER = [l for l in ALL_L3_ORDER if l in avail_L3]

adata.obs['cell_type_L3_ordered'] = pd.Categorical(
    adata.obs['cell_type_L3_new'], categories=ALL_L3_ORDER
)

# sc.settings.figdir is already set to fig_dir (Section 1)
# scanpy prefixes 'stacked_violins_' to the save= argument
# -> final file: {fig_dir}/stacked_violins_tcell_all_L3.pdf
sc.pl.stacked_violin(
    adata,
    var_names=VIOLIN_GENES,
    groupby='cell_type_L3_ordered',
    use_raw=USE_RAW_EFFECTIVE,
    figsize=(16, max(10, len(ALL_L3_ORDER) * 0.6)),
    dendrogram=False,
    show=False,
    save='_tcell_all_L3.pdf',
)

# Confirm file was written
expected = fig_dir / 'stacked_violins_tcell_all_L3.pdf'
if expected.exists():
    print(f"Saved: {expected.name}")
else:
    written = sorted(fig_dir.glob('stacked_violins*.pdf'))
    print(f"Actual filename(s): {[p.name for p in written]}")

# %% [markdown]
# ## 12. Summary & Export

# %%
# -----------------------------------------------------------------------
# Annotation summary table: L2 / L3 / L4 / pipeline cluster / cell count
# -----------------------------------------------------------------------
summary = (
    adata.obs
    .groupby(['cell_type_L2_new', 'cell_type_L3_new', 'cell_type_L4_new', CLUSTER_KEY])
    .size()
    .reset_index(name='n_cells')
)
summary['pct_total'] = (summary['n_cells'] / adata.n_obs * 100).round(2)
summary.columns = ['L2', 'L3', 'L4', 'L4_pipeline', 'n_cells', 'pct_total']
print(summary.to_string(index=False))
summary.to_csv(output_dir / 'annotation_summary.csv', index=False)

# -----------------------------------------------------------------------
# Per-cell annotation export (all four levels)
# -----------------------------------------------------------------------
export = adata.obs[['cell_type_L2_new', 'cell_type_L3_new',
                     'cell_type_L4_new', CLUSTER_KEY]].copy()
export.columns = ['L2_updated', 'L3_updated', 'L4_updated', 'L4_pipeline']
export.to_csv(output_dir / 'tcell_annotations_updated.csv')

print("\n" + "=" * 60)
print("DONE")
print(f"Figures ({len(list(fig_dir.glob('*.pdf')))}) in: {fig_dir}")
for f in sorted(fig_dir.glob('*.pdf')):
    print(f"  {f.name}")