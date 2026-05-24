# %% [markdown]
# # T/NK Cell scVI + scANVI Reference Model (scArches-ready)
# Version: v1.2 (2026-03-15)
#
# **Purpose**:
#   1. Apply refined L4/L3/L2 annotation labels
#   2. Remove epithelial-contaminated doublet clusters (4 clusters)
#   3. Train scVI with categorical covariates (encode_covariates=True for scArches)
#   4. Train scANVI using L3 as supervision labels
#   5. Save reference model + var_names.csv for future scArches query mapping
#
# **v1.2 fixes vs v1.1**:
#   [BUG]  Sec 5: sample-only counts check missed global Inf/NaN values; HVG then
#          crashed with "input data contains infinity" inside pd.cut(means, bins=...)
#   [FIX1] Global finite check via X.data (not sample) for all candidate layers
#   [FIX2] Prefer 'raw_counts' over 'counts' when both present and valid
#          ('counts' may have been silently overwritten mid-pipeline)
#   [FIX3] sanitize_sparse_counts() replaces non-finite entries, eliminates zeros
#   [FIX4] filter_genes(min_cells=3) + filter_cells(min_genes=200) BEFORE HVG
#          stabilizes seurat_v3 loess (near-singularity on zero/near-constant genes)
#   [FIX5] 4-tier HVG fallback with explicit flavor at every level:
#          seurat_v3(batch,counts) -> seurat_v3(non-batch,counts) ->
#          seurat(batch,log1p) -> seurat(non-batch,log1p)
#   [FIX6] Final assert: counts must be fully finite and non-negative before HVG
#
# **v1.1 fixes vs v1.0**:
#   [P0-1] PIPELINE_KEY separated from CLUSTER_KEY -- prevents doublet filter failure
#          caused by overwriting cell_type_L3 column before mask_clean was built
#   [P0-2] sc.pl.embedding basis='umap' (was 'X_umap', caused KeyError on obsm lookup)
#   [P0-3] matplotlib.use('Agg') moved before pyplot import
#   [P1-1] .raw.X now stores log1p (was counts; dotplot/featureplot dynamic range fix)
#   [P1-2] counts reconstruction checks get_indexer for -1 and validates distribution
#   [P1-4] n_pcs removed from sc.pp.neighbors (use_rep makes it redundant / version trap)
#   [P1-5] n_samples_per_label adapted to min class size (avoids ILC3/gdT oversampling)
#   [P2-2] Removed sc.read(model.pt) (nonsense); training loss curve saved BEFORE del
#   [P2-3] uns records raw_x_semantics, counts_layer_semantics, scanvi_label_categories
#   [self] leiden wrapped in try-except; scanvi_label_categories saved to uns for query
#   [self] n_unlabeled check now uses PIPELINE_KEY (valid after P0-1 fix)
#
# **Input** : adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad
# **Output** :
#   - adata_tnk_scanvi_ref_20260315_v1_2.h5ad
#   - tnk_scvi_ref_model/ (scVI model dir, includes var_names.csv)
#   - tnk_scanvi_ref_model/ (scANVI model dir)
#
# **Expected data structure (output)**:
#   .X               : log1p normalized (HVG, float32 CSR)
#   .layers['counts']: raw counts (HVG, float32 CSR)
#   .layers['log1p'] : log1p normalized (HVG, float32 CSR)
#   .raw.X           : log1p normalized (FULL GENE SPACE, for feature/dotplots)
#   .raw.var          : full gene var table (minimal fields)
#   .obsm['X_scvi']  : scVI latent (n_cells x N_LATENT)
#   .obsm['X_scanvi']: scANVI latent (n_cells x N_LATENT)
#   .obsm['X_umap']  : UMAP 2D
#   .obs key columns : cell_type_L2/L3/L4/L4_pipeline, scanvi_label,
#                      scanvi_pred_L3, scanvi_pred_prob, leiden_r0.5, leiden_r1.0
#
# **Memory** : ~20-28 GB peak (HVG-only training, full-gene .raw shared memory)
# **Runtime**: ~1.5-2 h (GPU V100)
#
# **Key QRM rules applied**:
#   - QRM 15.1 : reference uses from_scvi_model (query must use load_query_data)
#   - QRM 15.2 : SCANVI_LABELS_KEY constant = 'scanvi_label'; query must match
#   - QRM 15.4 : covariates computed BEFORE gene alignment
#   - QRM 15.5 : full-gene .raw saved BEFORE HVG subset
#   - QRM 13   : .raw before scoring; uns_keys in order; category dtype before write

# %% [markdown]
# ## 0. Configuration

# %%
# ============================================================================
# CONFIGURATION  (edit here only)
# ============================================================================

INPUT_H5AD  = "/home/h2048/data/py/0129/tnk_analysis_unified/results/subcluster_unified_v2_20260129/adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"
OUTPUT_DIR  = "/home/h2048/data/py/0315/tnk_scarches_ref"

# Source pipeline cluster column in input h5ad (values like 'trm_cytotoxic_t_cells_c0')
# This column is READ-ONLY; do NOT overwrite it.
CLUSTER_KEY = "cell_type_L3"

# New frozen copy of pipeline clusters (written to obs, then used for all lookups)
# Separating this from CLUSTER_KEY prevents the P0-1 bug where L3 overwrite
# breaks doublet filtering that depends on original pipeline cluster names.
PIPELINE_KEY = "cell_type_L4_pipeline"

# Batch key for scVI/scANVI training and HVG selection
# Use 'sample' (individual samples) for scVI -- finer-grained correction.
# 'dataset' is reserved for BBKNN only (graph-based integration needs coarser batches).
BATCH_KEY   = "sample"

# scANVI label key -- MUST be kept identical in all future query scripts (QRM 15.2)
SCANVI_LABELS_KEY   = "scanvi_label"
UNLABELED_CATEGORY  = "Unknown"

# scVI / scANVI training parameters
N_HVG         = 6000
N_LATENT      = 75
N_HIDDEN      = 128
N_LAYERS      = 2
DROPOUT       = 0.1
SCVI_EPOCHS   = 400
SCANVI_EPOCHS = 200
BATCH_SIZE    = 256    # reduce to 128 if GPU OOM

# Categorical covariates to encode (checked for existence in obs before use)
CANDIDATE_CAT_COVARIATES  = ['tissue', 'disease_status', 'sex']
CANDIDATE_CONT_COVARIATES = []   # e.g. ['log10_n_counts'] if available

# UMAP / plotting
DPI         = 300
N_NEIGHBORS = 50

# Random seed
SEED = 42

print("Configuration loaded")
print(f"  Input       : {INPUT_H5AD}")
print(f"  Output      : {OUTPUT_DIR}")
print(f"  CLUSTER_KEY : {CLUSTER_KEY}  (read-only source)")
print(f"  PIPELINE_KEY: {PIPELINE_KEY} (frozen copy for all lookups)")

# %% [markdown]
# ## 1. Imports & Setup
#
# P0-3 fix: matplotlib.use('Agg') MUST come before pyplot import to take effect.

# %%
import os, gc, warnings
warnings.filterwarnings('ignore')

# Set torch/BLAS threads BEFORE importing torch (prevents oversubscription)
os.environ['OMP_NUM_THREADS']      = '8'
os.environ['MKL_NUM_THREADS']      = '8'
os.environ['OPENBLAS_NUM_THREADS'] = '8'

import matplotlib
matplotlib.use('Agg')           # P0-3 fix: before pyplot import
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scipy.sparse as sparse
import scanpy as sc
import scvi
from pathlib import Path

sc.settings.verbosity  = 2
sc.settings.n_jobs     = 16    # P2-1: reduced from 48 to avoid CPU oversubscription
sc.settings.set_figure_params(dpi=DPI, facecolor='white', frameon=False)

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
fig_dir    = output_dir / 'figures'
fig_dir.mkdir(exist_ok=True)
sc.settings.figdir = str(fig_dir)

scvi.settings.seed = SEED
np.random.seed(SEED)

print(f"scanpy  : {sc.__version__}")
print(f"scvi    : {scvi.__version__}")
print(f"Output  : {output_dir}")

# %% [markdown]
# ## 2. Annotation Mappings
#
# Keys are PIPELINE_KEY values (original pipeline cluster names).
# L4_TO_L4_LABEL : refined publication labels (markers in parentheses)
# L4_TO_L3       : pipeline cluster -> L3 functional subtype (scANVI supervision)
# L3_TO_L2       : L3 -> L2 major lineage

# %%
# ============================================================================
# CONTAMINATION CLUSTERS (epithelial doublets) -- removed before all training
# Identified by co-expression of T/NK markers + epithelial markers
# ============================================================================
CONTAMINATION_CLUSTERS = [
    'trm_cytotoxic_t_cells_c1',          # WFDC2/SLPI high
    'cd8plus_trm_cytotoxic_t_cells_c1',  # Epi-associated
    'cd8plus_trm_cytotoxic_t_cells_c2',  # BPIFA1/STATH high
    'mait_cells_c1',                     # BPIFA1/STATH high
]

# ============================================================================
# L4 clean labels  (cell_type_L4)
# No gene markers -- suitable for figure legends, dotplot groupby, publication panels.
# Two NK Cytotoxic CD16+ clusters distinguished by dominant marker protein name.
# ============================================================================
L4_TO_L4_LABEL = {
    # NK ---------------------------------------------------------------
    'cd16-_nk_cells_c0':                       'NK Exhausted CD16-',
    'cd16plus_nk_cells_c0':                    'NK Cytotoxic CD16+ HAVCR2hi',
    'cd16plus_nk_cells_c1':                    'NK Cytotoxic CD16+ CX3CR1hi',
    'cd16plus_nk_cells_c2':                    'NK AREG-high CD16+',
    'nk_cells_c0':                             'NK ISG/IFN-high',
    'nk_cells_c1':                             'NK Resting',
    'ilc3_c0':                                 'ILC3',
    # CD4 --------------------------------------------------------------
    'tcm_naive_helper_t_cells_c0':             'CD4 Naive/TCM',
    'cd4plus_tem_effector_helper_t_cells_c0':  'CD4 TCM',
    'tem_effector_helper_t_cells_c0':          'CD4 TEM Inflammatory',
    'tem_effector_helper_t_cells_c1':          'CD4 TEM Activated',
    'type_1_helper_t_cells_c0':                'CD4 Th1',
    'type_17_helper_t_cells_c0':               'CD4 Th17',
    'follicular_helper_t_cells_c0':            'CD4 Tfh Resting',
    'follicular_helper_t_cells_c1':            'CD4 Tfh PD-1-high',
    'regulatory_t_cells_c0':                   'CD4 Tfh/Tfr-like',
    'cd4plus_regulatory_t_cells_c0':           'CD4 Tfr-like',
    'regulatory_t_cells_c1':                   'CD4 Treg Activated',
    'regulatory_t_cells_c2':                   'CD4 Tfr-like',
    # CD8 --------------------------------------------------------------
    'cd8plus_tem_effector_helper_t_cells_c0':  'CD8 Naive/TCM',
    'cd8plus_tem_temra_cytotoxic_t_cells_c0':  'CD8 TEMRA',
    'cd8plus_tem_temra_cytotoxic_t_cells_c1':  'CD8 TEMRA Terminal CX3CR1hi',
    'tem_temra_cytotoxic_t_cells_c0':          'CD8 TEMRA Terminal ZEB2hi',
    'tem_temra_cytotoxic_t_cells_c1':          'CD8 TEMRA Terminal CX3CR1hi',
    'cd8plus_tem_trm_cytotoxic_t_cells_c0':    'CD8 TRM Activated',
    'cd8plus_tem_trm_cytotoxic_t_cells_c1':    'CD8 TRM HLA-II/ISG-high',
    'tem_trm_cytotoxic_t_cells_c0':            'CD8 TRM Activated',
    'tem_trm_cytotoxic_t_cells_c1':            'CD8 TRM-like',
    'cd8plus_trm_cytotoxic_t_cells_c0':        'CD8 TRM Resting',
    'trm_cytotoxic_t_cells_c0':                'CD8 TRM Resting',
    'trm_cytotoxic_t_cells_c2':                'CD8 TRM MT-high',
    'cd8plus_gamma-delta_t_cells_c0':          'gammadelta TRM-like CXCL13hi',
    'gamma-delta_t_cells_c0':                  'gammadelta TRM-like TRDV1hi',
    'mait_cells_c0':                           'MAIT Type17-like',
    # Contamination (audit only; cells removed before training)
    'trm_cytotoxic_t_cells_c1':               'CD8 TRM Epi-associated [doublet]',
    'cd8plus_trm_cytotoxic_t_cells_c1':       'CD8 TRM Epi-associated [doublet]',
    'cd8plus_trm_cytotoxic_t_cells_c2':       'CD8 TRM Epi-associated [doublet]',
    'mait_cells_c1':                          'MAIT Epi-associated [doublet]',
}

# ============================================================================
# L4 labels WITH marker genes  (cell_type_L4_markers)
# For internal QC, dotplot titles, supplementary tables.
# Pattern: "{clean_label} ({key markers})"
# ============================================================================
L4_TO_L4_LABEL_MARKERS = {
    # NK ---------------------------------------------------------------
    'cd16-_nk_cells_c0':                       'NK Exhausted CD16- (TOX/GZMK/NCR1)',
    'cd16plus_nk_cells_c0':                    'NK Cytotoxic CD16+ HAVCR2hi (HAVCR2/FGFBP2/GZMB)',
    'cd16plus_nk_cells_c1':                    'NK Cytotoxic CD16+ CX3CR1hi (CX3CR1/FGFBP2)',
    'cd16plus_nk_cells_c2':                    'NK AREG-high CD16+ (AREG/XCL1/XCL2)',
    'nk_cells_c0':                             'NK ISG/IFN-high (ISG15/IFIT1/MX1; HLA-DRA)',
    'nk_cells_c1':                             'NK Resting (CXCR4/KLRC1/KLRD1)',
    'ilc3_c0':                                 'ILC3 (KIT/IL7R/IL1R1/IL18R1)',
    # CD4 --------------------------------------------------------------
    'tcm_naive_helper_t_cells_c0':             'CD4 Naive/TCM (CCR7/SELL/TCF7)',
    'cd4plus_tem_effector_helper_t_cells_c0':  'CD4 TCM (MAL/LEF1/TCF7)',
    'tem_effector_helper_t_cells_c0':          'CD4 TEM Inflammatory (SLC2A3/RGS1)',
    'tem_effector_helper_t_cells_c1':          'CD4 TEM Activated (ICOS/CD40LG)',
    'type_1_helper_t_cells_c0':                'CD4 Th1 (IFNG/TNF/TBX21)',
    'type_17_helper_t_cells_c0':               'CD4 Th17 (IL23R/CCR6/RORC)',
    'follicular_helper_t_cells_c0':            'CD4 Tfh Resting (CXCR5/TCF7/IL6R)',
    'follicular_helper_t_cells_c1':            'CD4 Tfh PD-1-high (PDCD1/TIGIT/CXCL13)',
    'regulatory_t_cells_c0':                   'CD4 Tfh/Tfr-like (PDCD1/CXCL13/TOX2)',
    'cd4plus_regulatory_t_cells_c0':           'CD4 Tfr-like (CXCR5/CTLA4/PDCD1)',
    'regulatory_t_cells_c1':                   'CD4 Treg Activated (FOXP3/IL2RA/ICOS)',
    'regulatory_t_cells_c2':                   'CD4 Tfr-like (FOXP3/CTLA4/CXCR5/PDCD1)',
    # CD8 --------------------------------------------------------------
    'cd8plus_tem_effector_helper_t_cells_c0':  'CD8 Naive/TCM (IL7R/TCF7/CCR7/SELL)',
    'cd8plus_tem_temra_cytotoxic_t_cells_c0':  'CD8 TEMRA (KLRG1/FGFBP2/GZMH)',
    'cd8plus_tem_temra_cytotoxic_t_cells_c1':  'CD8 TEMRA Terminal CX3CR1hi (CX3CR1/KLRG1/FGFBP2)',
    'tem_temra_cytotoxic_t_cells_c0':          'CD8 TEMRA Terminal ZEB2hi (ZEB2/KLRG1/FGFBP2)',
    'tem_temra_cytotoxic_t_cells_c1':          'CD8 TEMRA Terminal CX3CR1hi (CX3CR1/KLRG1/S1PR5)',
    'cd8plus_tem_trm_cytotoxic_t_cells_c0':    'CD8 TRM Activated (ITGAE/GZMK/IFNG)',
    'cd8plus_tem_trm_cytotoxic_t_cells_c1':    'CD8 TRM HLA-II/ISG-high (CD74/HLA-DPA1/ISG15)',
    'tem_trm_cytotoxic_t_cells_c0':            'CD8 TRM Activated (CD69/RGS1/CXCR6)',
    'tem_trm_cytotoxic_t_cells_c1':            'CD8 TRM-like (KLRG1/IL7R/TRAT1)',
    'cd8plus_trm_cytotoxic_t_cells_c0':        'CD8 TRM Resting (ITGA1/CXCR6/ZNF683)',
    'trm_cytotoxic_t_cells_c0':                'CD8 TRM Resting (ITGA1/CXCR6/ZNF683)',
    'trm_cytotoxic_t_cells_c2':                'CD8 TRM MT-high (MT1E/MT1X)',
    'cd8plus_gamma-delta_t_cells_c0':          'gammadelta TRM-like CXCL13hi (CXCL13/ENTPD1/TRDC)',
    'gamma-delta_t_cells_c0':                  'gammadelta TRM-like TRDV1hi (TRDV1/KLRC1/TRGV9)',
    'mait_cells_c0':                           'MAIT Type17-like (KLRB1/ZBTB16/IL18R1/CCR6)',
    # Contamination (audit only; cells removed before training)
    'trm_cytotoxic_t_cells_c1':               'CD8 TRM Epi-associated (WFDC2/SLPI/SCGB1A1) [doublet]',
    'cd8plus_trm_cytotoxic_t_cells_c1':       'CD8 TRM Epi-associated (WFDC2/SLPI) [doublet]',
    'cd8plus_trm_cytotoxic_t_cells_c2':       'CD8 TRM Epi-associated (BPIFA1/STATH/WFDC2) [doublet]',
    'mait_cells_c1':                          'MAIT Epi-associated (BPIFA1/STATH) [doublet]',
}

# ============================================================================
# L4 -> L3 mapping (scANVI supervision level)
# Contamination clusters NOT included -- they are removed before training
# ============================================================================
L4_TO_L3 = {
    'cd16-_nk_cells_c0':                       'NK Exhausted',
    'cd16plus_nk_cells_c0':                    'NK Cytotoxic',
    'cd16plus_nk_cells_c1':                    'NK Cytotoxic',
    'cd16plus_nk_cells_c2':                    'NK Cytotoxic',
    'nk_cells_c0':                             'NK Inflammatory',
    'nk_cells_c1':                             'NK Resting',
    'ilc3_c0':                                 'ILC3',
    'tcm_naive_helper_t_cells_c0':             'CD4 Naive/TCM',
    'cd4plus_tem_effector_helper_t_cells_c0':  'CD4 Naive/TCM',
    'tem_effector_helper_t_cells_c0':          'CD4 TEM',
    'tem_effector_helper_t_cells_c1':          'CD4 TEM',
    'type_1_helper_t_cells_c0':                'CD4 Th1',
    'type_17_helper_t_cells_c0':               'CD4 Th17',
    'follicular_helper_t_cells_c0':            'CD4 Tfh',
    'follicular_helper_t_cells_c1':            'CD4 Tfh',
    'regulatory_t_cells_c0':                   'CD4 Tfh',    # Tfr-like under Tfh umbrella
    'cd4plus_regulatory_t_cells_c0':           'CD4 Treg',
    'regulatory_t_cells_c1':                   'CD4 Treg',
    'regulatory_t_cells_c2':                   'CD4 Treg',
    'cd8plus_tem_effector_helper_t_cells_c0':  'CD8 Naive/TCM',
    'cd8plus_tem_temra_cytotoxic_t_cells_c0':  'CD8 TEMRA',
    'cd8plus_tem_temra_cytotoxic_t_cells_c1':  'CD8 TEMRA',
    'tem_temra_cytotoxic_t_cells_c0':          'CD8 TEMRA',
    'tem_temra_cytotoxic_t_cells_c1':          'CD8 TEMRA',
    'cd8plus_tem_trm_cytotoxic_t_cells_c0':    'CD8 TRM',
    'cd8plus_tem_trm_cytotoxic_t_cells_c1':    'CD8 TRM',
    'tem_trm_cytotoxic_t_cells_c0':            'CD8 TRM',
    'tem_trm_cytotoxic_t_cells_c1':            'CD8 TRM',
    'cd8plus_trm_cytotoxic_t_cells_c0':        'CD8 TRM',
    'trm_cytotoxic_t_cells_c0':                'CD8 TRM',
    'trm_cytotoxic_t_cells_c2':                'CD8 TRM',
    'cd8plus_gamma-delta_t_cells_c0':          'CD8 Gamma-delta',
    'gamma-delta_t_cells_c0':                  'CD8 Gamma-delta',
    'mait_cells_c0':                           'CD8 MAIT',
}

L3_TO_L2 = {
    'NK Exhausted':     'NK cells',
    'NK Cytotoxic':     'NK cells',
    'NK Inflammatory':  'NK cells',
    'NK Resting':       'NK cells',
    'ILC3':             'NK cells',
    'CD4 Naive/TCM':    'CD4 T cells',
    'CD4 TEM':          'CD4 T cells',
    'CD4 Th1':          'CD4 T cells',
    'CD4 Th17':         'CD4 T cells',
    'CD4 Tfh':          'CD4 T cells',
    'CD4 Treg':         'CD4 T cells',
    'CD8 Naive/TCM':    'CD8 T cells',
    'CD8 TEMRA':        'CD8 T cells',
    'CD8 TRM':          'CD8 T cells',
    'CD8 Gamma-delta':  'CD8 T cells',
    'CD8 MAIT':         'CD8 T cells',
}

print("Annotation mappings defined")
print(f"  L4 labels  : {len(L4_TO_L4_LABEL)} entries (incl. 4 doublet audit entries)")
print(f"  L3 labels  : {len(L4_TO_L3)} clean clusters -> {len(set(L4_TO_L3.values()))} L3 types")
print(f"  Doublet clusters to remove: {CONTAMINATION_CLUSTERS}")

# %% [markdown]
# ## 3. Load Data + Apply Labels + Remove Doublets
#
# P0-1 fix:
#   1. Freeze pipeline cluster column FIRST as PIPELINE_KEY
#   2. Derive all new annotation columns from PIPELINE_KEY (not CLUSTER_KEY)
#   3. Doublet removal uses PIPELINE_KEY -- guaranteed to see original cluster names

# %%
print("Loading h5ad...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"  Shape     : {adata.shape}")
print(f"  obs cols  : {list(adata.obs.columns)}")
print(f"  obsm      : {list(adata.obsm.keys())}")
print(f"  layers    : {list(adata.layers.keys()) if adata.layers else 'None'}")
print(f"  .raw      : {adata.raw is not None}")

# %%
# Step 1: Freeze pipeline cluster column (CLUSTER_KEY -> PIPELINE_KEY)
# PIPELINE_KEY is the only column we use for lookups from this point on.
adata.obs[PIPELINE_KEY] = adata.obs[CLUSTER_KEY].astype(str).copy()

# Step 2: Derive all annotation columns from PIPELINE_KEY
adata.obs['cell_type_L4']         = adata.obs[PIPELINE_KEY].map(L4_TO_L4_LABEL)
adata.obs['cell_type_L4_markers'] = adata.obs[PIPELINE_KEY].map(L4_TO_L4_LABEL_MARKERS)
adata.obs['cell_type_L3']         = adata.obs[PIPELINE_KEY].map(L4_TO_L3)
adata.obs['cell_type_L2']         = adata.obs['cell_type_L3'].map(L3_TO_L2)

# Audit: contamination cell counts BEFORE removal
n_before = adata.n_obs
print(f"\nBefore doublet removal: {n_before:,} cells")
for c in CONTAMINATION_CLUSTERS:
    n = (adata.obs[PIPELINE_KEY] == c).sum()   # P0-1 fix: use PIPELINE_KEY
    print(f"  {c}: {n:,} cells (flagged for removal)")

# Step 3: Remove doublet clusters (use PIPELINE_KEY not CLUSTER_KEY)
mask_clean = ~adata.obs[PIPELINE_KEY].isin(CONTAMINATION_CLUSTERS)
adata = adata[mask_clean].copy()
n_after = adata.n_obs
print(f"\nAfter doublet removal : {n_after:,} cells "
      f"({n_before - n_after:,} removed, {(n_before-n_after)/n_before*100:.1f}%)")

# Sanity check: all PIPELINE_KEY values should now be clean (non-contamination)
remaining_contam = adata.obs[PIPELINE_KEY].isin(CONTAMINATION_CLUSTERS).sum()
if remaining_contam > 0:
    raise ValueError(f"Doublet removal failed: {remaining_contam} contamination cells remain")
print("  Doublet removal verified: OK")

# Sanity check: all cells should have valid L3 labels (no NaN)
n_unmapped_l3 = adata.obs['cell_type_L3'].isna().sum()
n_unmapped_l2 = adata.obs['cell_type_L2'].isna().sum()
if n_unmapped_l3 > 0:
    missing_keys = adata.obs.loc[adata.obs['cell_type_L3'].isna(), PIPELINE_KEY].value_counts()
    print(f"\nWARNING: {n_unmapped_l3} cells with unmapped L3 (add to L4_TO_L3):")
    for k, v in missing_keys.items():
        print(f"  '{k}': {v} cells")
else:
    print("  All clusters mapped to L3: OK")

if n_unmapped_l2 > 0:
    missing_l2 = adata.obs.loc[adata.obs['cell_type_L2'].isna(), 'cell_type_L3'].value_counts()
    print(f"\nWARNING: {n_unmapped_l2} cells with unmapped L2 (add to L3_TO_L2):")
    for k, v in missing_l2.items():
        print(f"  '{k}': {v} cells")

print("\nL2 composition:"); print(adata.obs['cell_type_L2'].value_counts())
print("\nL3 composition:"); print(adata.obs['cell_type_L3'].value_counts())

# %% [markdown]
# ## 4. Covariate Check  (QRM 15.4 -- MUST be before gene alignment)

# %%
# ============================================================================
# All covariate computation must happen on FULL GENE SPACE (this section).
# If you need MT%, cell cycle scores, or stress scores, compute them here.
# ============================================================================

cat_covariates  = [c for c in CANDIDATE_CAT_COVARIATES  if c in adata.obs.columns]
cont_covariates = [c for c in CANDIDATE_CONT_COVARIATES if c in adata.obs.columns]

print(f"Available categorical covariates : {cat_covariates}")
print(f"Available continuous  covariates : {cont_covariates}")

assert BATCH_KEY in adata.obs.columns, f"BATCH_KEY '{BATCH_KEY}' not found in obs"
print(f"Batch key '{BATCH_KEY}': {adata.obs[BATCH_KEY].nunique()} batches")

# Remove small batches (<3 cells) to prevent scVI/BBKNN instability
batch_counts  = adata.obs[BATCH_KEY].value_counts()
small_batches = batch_counts[batch_counts < 3].index.tolist()
if small_batches:
    print(f"\nRemoving {len(small_batches)} small batches: {small_batches}")
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    print(f"  After small-batch filter: {adata.n_obs:,} cells")
else:
    print("  No small batches found")

# Sanitize categorical covariate dtypes (prevent literal 'nan' strings)
for col in cat_covariates + [BATCH_KEY]:
    adata.obs[col] = adata.obs[col].astype('string').fillna('Unknown').astype('category')

# Covariate summary
print("\nCovariate summary:")
for col in cat_covariates:
    print(f"  {col}: {adata.obs[col].nunique()} categories")

# %% [markdown]
# ## 5. Data Integrity: Layer Validation, Sanitization, and Pre-Filtering
#
# v1.2 rewrite -- fixes the "input data contains infinity" crash in HVG:
#   - Global finite check via X.data (not a sample)
#   - Prefer 'raw_counts' over 'counts' when both present
#   - sanitize_sparse_counts() clears Inf/NaN in-place
#   - filter_genes + filter_cells BEFORE HVG to stabilize seurat_v3 loess

# %%
# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

def _layer_flat_data(X):
    """Return 1-D array of all STORED (non-structural-zero) values."""
    if sparse.issparse(X):
        return X.data          # only non-zero stored values -- O(nnz), not O(n*p)
    return np.asarray(X, dtype=np.float32).ravel()


def summarize_layer(X, name):
    """Print global stats for a layer and return (finite, is_int, min, max)."""
    arr = _layer_flat_data(X)
    if arr.size == 0:
        print(f"  {name}: EMPTY (all structural zeros or zero matrix)")
        return True, True, 0.0, 0.0
    finite   = bool(np.isfinite(arr).all())
    n_inf    = int(np.isinf(arr).sum())
    n_nan    = int(np.isnan(arr).sum())
    is_int   = bool(np.all(arr == np.floor(arr)))
    min_val  = float(arr.min())
    max_val  = float(arr.max())
    print(f"  {name}: finite={finite}, inf={n_inf}, nan={n_nan}, "
          f"integer={is_int}, min={min_val:.2f}, max={max_val:.2f}")
    return finite, is_int, min_val, max_val


def sanitize_sparse_counts(X):
    """
    Replace non-finite entries with 0 in-place (COO-safe CSR conversion).
    Returns (cleaned_csr, n_replaced).
    """
    X = sparse.csr_matrix(X, dtype=np.float32, copy=True)
    bad   = ~np.isfinite(X.data)
    n_bad = int(bad.sum())
    if n_bad > 0:
        X.data[bad] = 0.0
        X.eliminate_zeros()
    return X, n_bad


# ----------------------------------------------------------------------------
# Step 5a: Identify and select best raw-count layer (global check)
# ----------------------------------------------------------------------------
layers_present = list(adata.layers.keys()) if adata.layers else []
print(f"Layers present: {layers_present}")
print("\nGlobal layer diagnostics (FULL matrix, not sample):")

candidates = []   # (layer_name, finite, is_int, min_val, max_val)
for lname in ['raw_counts', 'counts']:   # raw_counts preferred (may be more pristine)
    if lname in layers_present:
        finite, is_int, min_val, max_val = summarize_layer(adata.layers[lname], lname)
        candidates.append((lname, finite, is_int, min_val, max_val))

if not candidates:
    # Last resort: reconstruct from .raw
    if adata.raw is not None:
        print("  No raw-count layer found; reconstructing from .raw...")
        indexer  = adata.raw.var_names.get_indexer(adata.var_names)
        n_miss   = int((indexer == -1).sum())
        if n_miss > 0:
            raise ValueError(
                f"Cannot reconstruct counts from .raw: {n_miss} genes not found "
                f"in .raw.var -- .raw may be log-transformed or misaligned"
            )
        adata.layers['counts'] = sparse.csr_matrix(
            adata.raw.X[:, indexer], dtype=np.float32
        )
        summarize_layer(adata.layers['counts'], 'counts (from .raw)')
        candidates = [('counts', True, True, 0.0, 1.0)]  # assume ok, sanitize below
    else:
        raise ValueError(
            "No raw-count layer ('raw_counts', 'counts') and no .raw. "
            "Input h5ad must contain at least one raw count layer."
        )

# Choose: prefer raw_counts > counts; require finite + integer + non-negative
chosen = None
for lname, finite, is_int, min_val, max_val in candidates:
    if finite and is_int and min_val >= 0:
        chosen = lname
        break

if chosen is None:
    # No fully valid layer -- pick best candidate and sanitize
    chosen = candidates[0][0]
    print(f"\nWARNING: no fully valid raw-count layer found; sanitizing '{chosen}'")
else:
    print(f"\nSelected raw-count layer: '{chosen}'")

# Materialize as 'counts' (float32 CSR), sanitizing non-finite values
adata.layers['counts'], n_bad_counts = sanitize_sparse_counts(adata.layers[chosen])
if n_bad_counts > 0:
    print(f"  Sanitized {n_bad_counts} non-finite entries in '{chosen}' -> set to 0")

# Final assertion: counts must be fully finite and non-negative
counts_flat = _layer_flat_data(adata.layers['counts'])
if not np.isfinite(counts_flat).all():
    raise ValueError("counts layer still contains non-finite values after sanitization")
if np.any(counts_flat < 0):
    raise ValueError("counts layer contains negative values -- not valid raw counts")
print("  counts layer: finite=True, non-negative=True (assertion passed)")

# ----------------------------------------------------------------------------
# Step 5b: Ensure log1p layer and set .X
# ----------------------------------------------------------------------------
if 'log1p' in layers_present:
    adata.layers['log1p'], n_bad_log1p = sanitize_sparse_counts(adata.layers['log1p'])
    if n_bad_log1p > 0:
        print(f"  Sanitized {n_bad_log1p} non-finite entries in 'log1p'")
    adata.X = sparse.csr_matrix(adata.layers['log1p'], dtype=np.float32)
    print("Using existing 'log1p' layer for .X")
else:
    print("Computing log1p from counts layer...")
    adata.X = sparse.csr_matrix(adata.layers['counts'], dtype=np.float32)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = sparse.csr_matrix(adata.X, dtype=np.float32)
    print("  log1p computed from counts and stored")

# Also print log1p diagnostics to confirm it looks reasonable
print("\nlog1p layer check:")
summarize_layer(adata.layers['log1p'], 'log1p')

# ----------------------------------------------------------------------------
# Step 5c: Pre-filter genes and cells BEFORE HVG
# (stabilizes seurat_v3 loess; removes near-constant / zero genes)
# This filter is on full gene space -- .raw (Sec 7) will capture post-filter
# full gene space, which is the correct reference space for this dataset.
# ----------------------------------------------------------------------------
print("\nPre-filtering (stabilize HVG loess)...")
n_vars_before = adata.n_vars
n_obs_before  = adata.n_obs
sc.pp.filter_genes(adata, min_cells=3)    # remove near-absent genes
sc.pp.filter_cells(adata, min_genes=200)  # remove degenerate cells
print(f"  Genes : {n_vars_before:,} -> {adata.n_vars:,} "
      f"(removed {n_vars_before - adata.n_vars:,})")
print(f"  Cells : {n_obs_before:,} -> {adata.n_obs:,} "
      f"(removed {n_obs_before - adata.n_obs:,})")

# %% [markdown]
# ## 6. HVG Selection (4-tier explicit fallback)
#
# v1.2: explicit flavor at every tier; last two tiers run on 'log1p' not 'counts'
# Tier 1: seurat_v3 batch-aware  on counts  (best)
# Tier 2: seurat_v3 non-batch    on counts  (if batch causes instability)
# Tier 3: seurat    batch-aware  on log1p   (if seurat_v3 still fails)
# Tier 4: seurat    non-batch    on log1p   (guaranteed to work fallback)

# %%
print("Selecting highly variable genes...")
hvg_method = None

try:
    sc.pp.highly_variable_genes(
        adata, layer='counts', n_top_genes=N_HVG,
        batch_key=BATCH_KEY, flavor='seurat_v3', subset=False,
    )
    hvg_method = 'seurat_v3 batch-aware (counts)'

except Exception as e1:
    print(f"  Tier 1 failed (seurat_v3 batch-aware counts): {e1}")
    try:
        sc.pp.highly_variable_genes(
            adata, layer='counts', n_top_genes=N_HVG,
            flavor='seurat_v3', subset=False,
        )
        hvg_method = 'seurat_v3 non-batch (counts)'

    except Exception as e2:
        print(f"  Tier 2 failed (seurat_v3 non-batch counts): {e2}")
        try:
            sc.pp.highly_variable_genes(
                adata, layer='log1p', n_top_genes=N_HVG,
                batch_key=BATCH_KEY, flavor='seurat', subset=False,
            )
            hvg_method = 'seurat batch-aware (log1p)'

        except Exception as e3:
            print(f"  Tier 3 failed (seurat batch-aware log1p): {e3}")
            sc.pp.highly_variable_genes(
                adata, layer='log1p', n_top_genes=N_HVG,
                flavor='seurat', subset=False,
            )
            hvg_method = 'seurat non-batch (log1p)'

n_hvg = int(adata.var['highly_variable'].sum())
print(f"  {n_hvg} HVGs selected  [{hvg_method}]")
if n_hvg == 0:
    raise ValueError("No HVGs selected -- check input data")

# %% [markdown]
# ## 7. Save Full-Gene .raw  (QRM 15.5 -- BEFORE HVG subset)
#
# P1-1 fix: .raw.X stores log1p so use_raw=True in dotplot/featureplot has
# proper dynamic range (not raw count magnitudes).
# "Full gene" here means post-filter (Sec 5c) but pre-HVG subset -- this is
# the correct reference space: clean, normalized, covering all expressed genes.

# %%
print("Saving full-gene .raw with log1p values (shared memory, no .copy())...")
# Use adata.layers['log1p'] which is the same object as adata.X at this point.
# Shared memory: no extra copy of the full-gene matrix.
adata.raw = sc.AnnData(
    X   = adata.layers['log1p'],   # P1-1 fix: log1p (not counts)
    obs = adata.obs[[]].copy(),    # minimal obs (no metadata columns)
    var = adata.var[['gene_ids', 'feature_types']].copy()
          if all(c in adata.var.columns for c in ['gene_ids', 'feature_types'])
          else adata.var[[]].copy(),  # minimal var (strip HVG columns)
)
print(f"  .raw saved: {adata.raw.n_obs:,} cells x {adata.raw.n_vars:,} genes (log1p)")

# %% [markdown]
# ## 8. HVG Subset

# %%
adata_hvg = adata[:, adata.var['highly_variable']].copy()
print(f"HVG subset shape: {adata_hvg.shape}")

# adata[:, hvg_mask].copy() does NOT carry .raw in all AnnData versions.
# adata_hvg.raw = adata.raw FAILS: the setter only accepts AnnData, not Raw.
# (QRM 15.5: AnnData setter不接受Raw对象)
# Correct fix: reconstruct a fresh AnnData from adata.raw's components.
adata_hvg.raw = sc.AnnData(
    X   = adata.raw.X,          # shared reference to the log1p full-gene matrix
    obs = adata_hvg.obs[[]].copy(),   # minimal obs matching HVG subset cell order
    var = adata.raw.var.copy(),       # full-gene var (post-filter, pre-HVG-subset)
)
print(f"  .raw: {adata_hvg.raw.n_vars:,} genes (full space, log1p) -- reconstructed from Raw")

# Enforce float32 CSR for memory efficiency
adata_hvg.layers['counts'] = sparse.csr_matrix(adata_hvg.layers['counts'], dtype=np.float32)
adata_hvg.X                = sparse.csr_matrix(adata_hvg.X,                dtype=np.float32)
if 'log1p' in adata_hvg.layers:
    adata_hvg.layers['log1p'] = sparse.csr_matrix(adata_hvg.layers['log1p'], dtype=np.float32)

# Free full-gene adata (raw.X is still referenced by adata_hvg.raw, so not freed)
del adata
gc.collect()
print("  Full adata freed (full-gene counts shared in .raw)")

# %% [markdown]
# ## 9. scVI Setup and Training (encode_covariates=True for scArches)

# %%
# Prepare SCANVI label column (all cells are clean after doublet removal)
adata_hvg.obs[SCANVI_LABELS_KEY] = (
    adata_hvg.obs['cell_type_L3'].astype(str).astype('category')
)

# Print per-class cell counts (needed for P1-5 adaptive n_samples_per_label)
print(f"scANVI label distribution ('{SCANVI_LABELS_KEY}'):")
label_counts = adata_hvg.obs[SCANVI_LABELS_KEY].value_counts()
print(label_counts)
min_class_size = int(label_counts.min())
print(f"\nSmallest class: '{label_counts.idxmin()}' ({min_class_size} cells)")

# %%
# Convert all label / annotation columns to category dtype before setup_anndata
label_cols = [
    'cell_type_L2', 'cell_type_L3', 'cell_type_L4', PIPELINE_KEY,
    SCANVI_LABELS_KEY, BATCH_KEY,
] + cat_covariates
for col in label_cols:
    if col in adata_hvg.obs.columns:
        adata_hvg.obs[col] = adata_hvg.obs[col].astype('category')

# %%
# scVI setup_anndata
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer                        = 'counts',
    batch_key                    = BATCH_KEY,
    categorical_covariate_keys   = cat_covariates if cat_covariates else None,
    continuous_covariate_keys    = cont_covariates if cont_covariates else None,
)
print("scVI setup_anndata complete")
print(f"  Batch key       : {BATCH_KEY} ({adata_hvg.obs[BATCH_KEY].nunique()} batches)")
print(f"  Cat covariates  : {cat_covariates}")
print(f"  Cont covariates : {cont_covariates}")
print(f"  encode_covariates=True, deeply_inject_covariates=False (scArches standard)")

# %%
# Build scVI model
# encode_covariates=True + deeply_inject_covariates=False: required for scArches surgery
vae = scvi.model.SCVI(
    adata_hvg,
    n_latent                 = N_LATENT,
    n_hidden                 = N_HIDDEN,
    n_layers                 = N_LAYERS,
    dropout_rate             = DROPOUT,
    gene_likelihood          = 'nb',
    encode_covariates        = True,   # scArches requirement
    deeply_inject_covariates = False,  # scArches standard
    use_layer_norm           = 'both',
    use_batch_norm           = 'none',
)

# %%
print(f"\nTraining scVI ({SCVI_EPOCHS} max epochs, batch_size={BATCH_SIZE})...")
vae.train(
    max_epochs               = SCVI_EPOCHS,
    batch_size               = BATCH_SIZE,
    train_size               = 0.9,
    early_stopping           = True,
    early_stopping_patience  = 30,
    plan_kwargs              = {'lr': 1e-3},
)
print("scVI training complete")

# %%
# Extract latent + save model BEFORE training history is lost (P2-2 fix)
adata_hvg.obsm['X_scvi'] = vae.get_latent_representation()
print(f"X_scvi shape: {adata_hvg.obsm['X_scvi'].shape}")

# Plot scVI training loss curve (must happen BEFORE del vae)
try:
    fig_loss, ax_loss = plt.subplots(figsize=(8, 4))
    train_hist = vae.history.get('train_loss_epoch', None)
    val_hist   = vae.history.get('elbo_validation', None)
    if train_hist is not None:
        ax_loss.plot(train_hist.values.flatten(), label='train ELBO')
    if val_hist is not None:
        ax_loss.plot(val_hist.values.flatten(), label='val ELBO')
    ax_loss.set_xlabel('Epoch'); ax_loss.set_ylabel('ELBO'); ax_loss.legend()
    ax_loss.set_title('scVI Training Loss')
    fig_loss.savefig(fig_dir / 'scvi_training_loss.pdf', bbox_inches='tight', dpi=DPI)
    plt.close(fig_loss)
    print("Saved: scvi_training_loss.pdf")
except Exception as e:
    print(f"  scVI loss plot skipped: {e}")

# Save scVI model + var_names.csv (critical for scArches query gene alignment)
scvi_model_dir = str(output_dir / 'tnk_scvi_ref_model')
vae.save(scvi_model_dir, overwrite=True)

var_names_path = Path(scvi_model_dir) / 'var_names.csv'
pd.Series(adata_hvg.var_names.tolist()).to_csv(
    var_names_path, index=False, header=False
)
print(f"scVI model saved : {scvi_model_dir}/")
print(f"var_names.csv    : {var_names_path} ({len(adata_hvg.var_names)} genes)")

# %% [markdown]
# ## 10. scANVI Training (L3 supervision)

# %%
# ============================================================================
# Build scANVI from scVI using from_scvi_model.
# This is the REFERENCE training path -- correct here.
# QUERY scripts MUST use SCANVI.load_query_data(adata, ref_scanvi) (QRM 15.1)
# ============================================================================
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category = UNLABELED_CATEGORY,
    labels_key         = SCANVI_LABELS_KEY,
)

n_labeled   = (adata_hvg.obs[SCANVI_LABELS_KEY] != UNLABELED_CATEGORY).sum()
n_unlabeled = (adata_hvg.obs[SCANVI_LABELS_KEY] == UNLABELED_CATEGORY).sum()
print(f"scANVI setup:")
print(f"  labels_key    : {SCANVI_LABELS_KEY}")
print(f"  labeled cells : {n_labeled:,}")
print(f"  unlabeled     : {n_unlabeled:,}  (should be 0 after doublet removal)")
print(f"  classes       : {sorted(adata_hvg.obs[SCANVI_LABELS_KEY].cat.categories.tolist())}")

# %%
# P1-5 fix: adaptive n_samples_per_label to avoid oversampling rare classes
# Use 80% of smallest class size, capped at 100
n_samples_per_label = max(10, min(100, int(min_class_size * 0.8)))
print(f"\nAdaptive n_samples_per_label: {n_samples_per_label} "
      f"(min class={min_class_size}, 80% cap at 100)")

print(f"\nTraining scANVI ({SCANVI_EPOCHS} max epochs)...")
lvae.train(
    max_epochs              = SCANVI_EPOCHS,
    batch_size              = BATCH_SIZE,
    train_size              = 0.9,
    early_stopping          = True,
    early_stopping_patience = 20,
    n_samples_per_label     = n_samples_per_label,
)
print("scANVI training complete")

# %%
# Extract scANVI outputs
adata_hvg.obsm['X_scanvi']       = lvae.get_latent_representation()
adata_hvg.obs['scanvi_pred_L3']  = lvae.predict()
scanvi_proba                      = lvae.predict(soft=True)
adata_hvg.obs['scanvi_pred_prob'] = scanvi_proba.max(axis=1).values.astype(np.float32)

print(f"X_scanvi shape: {adata_hvg.obsm['X_scanvi'].shape}")
concordance = (
    adata_hvg.obs['scanvi_pred_L3'].astype(str) ==
    adata_hvg.obs[SCANVI_LABELS_KEY].astype(str)
).mean()
print(f"\nOverall concordance (pred vs label): {concordance:.3f}")
print(f"\nPrediction confidence:\n{adata_hvg.obs['scanvi_pred_prob'].describe()}")

low_conf = (adata_hvg.obs['scanvi_pred_prob'] < 0.5).sum()
print(f"\nLow-confidence cells (<0.5): {low_conf:,} ({low_conf/adata_hvg.n_obs*100:.1f}%)")

# Concordance breakdown by L3 class
print("\nPer-class concordance:")
for l3_class in sorted(adata_hvg.obs[SCANVI_LABELS_KEY].cat.categories):
    mask = adata_hvg.obs[SCANVI_LABELS_KEY].astype(str) == l3_class
    if mask.sum() == 0:
        continue
    acc = (adata_hvg.obs.loc[mask, 'scanvi_pred_L3'].astype(str) == l3_class).mean()
    print(f"  {l3_class}: {acc:.3f}  (n={mask.sum():,})")

# Plot scANVI training loss (must happen BEFORE del lvae)
try:
    fig_loss2, ax_loss2 = plt.subplots(figsize=(8, 4))
    train_hist2 = lvae.history.get('train_loss_epoch', None)
    if train_hist2 is not None:
        ax_loss2.plot(train_hist2.values.flatten(), label='train ELBO')
    ax_loss2.set_xlabel('Epoch'); ax_loss2.set_ylabel('ELBO'); ax_loss2.legend()
    ax_loss2.set_title('scANVI Training Loss')
    fig_loss2.savefig(fig_dir / 'scanvi_training_loss.pdf', bbox_inches='tight', dpi=DPI)
    plt.close(fig_loss2)
    print("Saved: scanvi_training_loss.pdf")
except Exception as e:
    print(f"  scANVI loss plot skipped: {e}")

# Save scANVI model
scanvi_model_dir = str(output_dir / 'tnk_scanvi_ref_model')
lvae.save(scanvi_model_dir, overwrite=True)
print(f"\nscANVI model saved: {scanvi_model_dir}/")

# Free GPU memory (del AFTER history/model already saved)
del vae, lvae
gc.collect()
try:
    import torch
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        print("GPU memory released")
except Exception:
    pass

# %% [markdown]
# ## 11. Neighbors + UMAP (on X_scanvi)

# %%
print("Computing neighbors on X_scanvi...")
sc.pp.neighbors(
    adata_hvg,
    use_rep     = 'X_scanvi',
    n_neighbors = N_NEIGHBORS,
    # n_pcs NOT set -- P1-4 fix: use_rep makes n_pcs redundant/version-dependent
)

sc.tl.umap(adata_hvg, min_dist=0.3, spread=1.0)
print(f"UMAP done: {adata_hvg.obsm['X_umap'].shape}")

# %%
# Leiden clustering (exploratory; for visualization/QC only)
try:
    for res in [0.5, 1.0]:
        sc.tl.leiden(adata_hvg, resolution=res, key_added=f'leiden_r{res}')
        print(f"  leiden_r{res}: {adata_hvg.obs[f'leiden_r{res}'].nunique()} clusters")
except Exception as e:
    print(f"  Leiden skipped (leidenalg not available): {e}")

# %% [markdown]
# ## 12. Visualization
#
# P0-2 fix: sc.pl.embedding uses basis='umap' (NOT 'X_umap').
# Passing 'X_umap' causes scanpy to look for obsm['X_X_umap'] -> KeyError.

# %%
# Overview UMAPs
fig, axes = plt.subplots(2, 3, figsize=(24, 14))

sc.pl.embedding(adata_hvg, basis='umap', color='cell_type_L2',   # P0-2 fix
                title='L2: Major Lineage', ax=axes[0, 0], show=False,
                legend_loc='on data', legend_fontsize=8)
sc.pl.embedding(adata_hvg, basis='umap', color='cell_type_L3',
                title='L3: Functional Subtype (ground truth)', ax=axes[0, 1], show=False,
                legend_loc='right margin', legend_fontsize=6)
sc.pl.embedding(adata_hvg, basis='umap', color='scanvi_pred_L3',
                title='L3: scANVI Prediction', ax=axes[0, 2], show=False,
                legend_loc='right margin', legend_fontsize=6)
sc.pl.embedding(adata_hvg, basis='umap', color=BATCH_KEY,
                title=f'Batch ({BATCH_KEY})', ax=axes[1, 0], show=False,
                legend_loc='right margin', legend_fontsize=5)
sc.pl.embedding(adata_hvg, basis='umap', color='scanvi_pred_prob',
                title='scANVI Prediction Confidence', ax=axes[1, 1], show=False,
                color_map='RdYlGn', vmin=0, vmax=1)
if cat_covariates:
    sc.pl.embedding(adata_hvg, basis='umap', color=cat_covariates[0],
                    title=f'Covariate: {cat_covariates[0]}', ax=axes[1, 2], show=False,
                    legend_loc='right margin', legend_fontsize=6)
elif 'leiden_r1.0' in adata_hvg.obs.columns:
    sc.pl.embedding(adata_hvg, basis='umap', color='leiden_r1.0',
                    title='Leiden r=1.0', ax=axes[1, 2], show=False,
                    legend_loc='on data', legend_fontsize=6)
else:
    axes[1, 2].set_visible(False)

plt.tight_layout()
# Rasterize scatter points post-hoc (rasterized= kwarg conflicts in scanpy 1.11+)
for ax in fig.axes:
    for coll in ax.collections:
        coll.set_rasterized(True)
fig.savefig(fig_dir / 'umap_overview.pdf', bbox_inches='tight', dpi=DPI)
plt.close(fig)
print("Saved: umap_overview.pdf")

# %%
# Key marker feature plots (use .raw = full-gene log1p space)
MARKER_GENES = [
    'CD3D', 'CD4', 'CD8A', 'GNLY', 'NKG7', 'FCGR3A',  # lineage
    'CCR7', 'SELL', 'TCF7',                              # naive/memory
    'GZMB', 'PRF1', 'GZMK',                             # cytotoxic
    'ITGA1', 'CXCR6', 'CD69',                           # TRM
    'PDCD1', 'HAVCR2', 'FOXP3',                         # exhaustion/Treg
    'KLRG1', 'CX3CR1',                                  # TEMRA
    'KIT', 'ZBTB16', 'TRDC',                            # ILC3/MAIT/gdT
]
use_raw_eff = adata_hvg.raw is not None
avail_raw   = set(adata_hvg.raw.var_names) if use_raw_eff else set(adata_hvg.var_names)
marker_valid = [g for g in MARKER_GENES if g in avail_raw]
missing_mks  = [g for g in MARKER_GENES if g not in avail_raw]
if missing_mks:
    print(f"Marker genes not in .raw: {missing_mks}")

ncols = 5
nrows = int(np.ceil(len(marker_valid) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
axes = axes.flatten()
for i, gene in enumerate(marker_valid):
    sc.pl.embedding(adata_hvg, basis='umap', color=gene,   # P0-2 fix
                    ax=axes[i], show=False, use_raw=use_raw_eff,
                    color_map='viridis', vmin=0)
    axes[i].set_title(gene, fontsize=10, fontweight='bold')
for i in range(len(marker_valid), len(axes)):
    axes[i].set_visible(False)
plt.suptitle('T/NK Key Marker Genes (scANVI UMAP, use_raw=True)', fontsize=13, y=1.01)
plt.tight_layout()
for ax in fig.axes:
    for coll in ax.collections:
        coll.set_rasterized(True)
fig.savefig(fig_dir / 'featureplot_markers.pdf', bbox_inches='tight', dpi=DPI)
plt.close(fig)
print("Saved: featureplot_markers.pdf")

# %% [markdown]
# ## 13. uns Metadata + Category dtype Conversion (QRM 13)

# %%
# Write uns in deterministic order (QRM 13 -- uns_keys order trap)
# P2-3 fix: add raw_x_semantics, counts_layer_semantics, scanvi_label_categories

scanvi_label_categories = sorted(
    adata_hvg.obs[SCANVI_LABELS_KEY].cat.categories.tolist()
)

# chosen_count_layer is set in Sec 5 (v1.2+); guard for older script versions
if 'chosen_count_layer' not in dir():
    chosen_count_layer = 'counts'   # fallback: assume counts was used

adata_hvg.uns['pipeline_params'] = {
    'script_version'           : 'v1.2',
    'n_hvg'                    : int(N_HVG),
    'hvg_method'               : hvg_method,
    'batch_key'                : BATCH_KEY,
    'cat_covariates'           : cat_covariates,
    'cont_covariates'          : cont_covariates,
    'n_latent'                 : int(N_LATENT),
    'n_hidden'                 : int(N_HIDDEN),
    'n_layers'                 : int(N_LAYERS),
    'dropout'                  : float(DROPOUT),
    'scvi_epochs'              : int(SCVI_EPOCHS),
    'scanvi_epochs'            : int(SCANVI_EPOCHS),
    'batch_size'               : int(BATCH_SIZE),
    'n_samples_per_label'      : int(n_samples_per_label),
    'seed'                     : int(SEED),
    'scanvi_labels_key'        : SCANVI_LABELS_KEY,
    'unlabeled_category'       : UNLABELED_CATEGORY,
    'doublets_removed'         : CONTAMINATION_CLUSTERS,
    'scvi_model_dir'           : scvi_model_dir,
    'scanvi_model_dir'         : scanvi_model_dir,
    # P2-3 fix: data semantics recorded for downstream use
    'raw_x_semantics'          : 'log1p_normalized_full_gene_space',
    'counts_layer_semantics'   : 'raw_integer_counts_HVG_subset',
    'x_semantics'              : 'log1p_normalized_HVG_subset',
    'counts_source_layer'      : chosen_count_layer,   # v1.2: which layer was used as counts
    # Self fix: label categories for query script validation (QRM 15.2)
    'scanvi_label_categories'  : scanvi_label_categories,
}

adata_hvg.uns['annotation_levels'] = {
    'L2'             : 'cell_type_L2  (NK cells / CD4 T cells / CD8 T cells)',
    'L3'             : 'cell_type_L3  (functional subtype, scANVI supervision label)',
    'L4'             : 'cell_type_L4  (clean label, no gene markers; use for legends/figures)',
    'L4_markers'     : 'cell_type_L4_markers  (label + key marker genes; for QC/supplements)',
    'L4_pipeline'    : f'{PIPELINE_KEY}  (original pipeline cluster key, read-only)',
}

# Convert all label/category columns to category dtype before write (QRM 13)
cat_cols_to_convert = [
    'cell_type_L2', 'cell_type_L3', 'cell_type_L4', 'cell_type_L4_markers',
    PIPELINE_KEY, SCANVI_LABELS_KEY, 'scanvi_pred_L3', BATCH_KEY,
] + cat_covariates
for col in cat_cols_to_convert:
    if col in adata_hvg.obs.columns:
        adata_hvg.obs[col] = adata_hvg.obs[col].astype('category')

print("uns metadata written")
print(f"  scanvi_label_categories ({len(scanvi_label_categories)}): {scanvi_label_categories}")
print("Category dtype conversion: OK")

# %% [markdown]
# ## 14. Save Output h5ad

# %%
out_h5ad = output_dir / 'adata_tnk_scanvi_ref_20260315_v1_2.h5ad'
print(f"Writing {out_h5ad.name}...")
print(f"  .X shape   : {adata_hvg.shape}  (HVG, log1p)")
print(f"  .raw shape : {adata_hvg.raw.n_obs:,} x {adata_hvg.raw.n_vars:,}  (full gene, log1p)")
print(f"  .layers    : {list(adata_hvg.layers.keys())}")
print(f"  .obsm      : {list(adata_hvg.obsm.keys())}")
adata_hvg.write_h5ad(out_h5ad, compression='gzip', compression_opts=9)
print(f"Saved: {out_h5ad}")

# %% [markdown]
# ## 15. Final Summary

# %%
print("=" * 70)
print("PIPELINE COMPLETE  (v1.2)")
print("=" * 70)

print(f"\nOutput h5ad  : {out_h5ad}")
print(f"scVI model   : {scvi_model_dir}/")
print(f"scANVI model : {scanvi_model_dir}/")
print(f"var_names    : {Path(scvi_model_dir) / 'var_names.csv'}")
print(f"Figures      : {fig_dir}/  ({len(list(fig_dir.glob('*.pdf')))} PDFs)")

print(f"\nAnnData structure:")
print(f"  .X               : {adata_hvg.shape} HVG log1p (float32 CSR)")
print(f"  .layers['counts']: HVG raw counts (float32 CSR)")
print(f"  .layers['log1p'] : HVG log1p (float32 CSR)")
print(f"  .raw.X           : {adata_hvg.raw.n_vars:,} genes full-space log1p")

print(f"\n.obs annotation columns:")
for col in ['cell_type_L2', 'cell_type_L3', 'cell_type_L4', 'cell_type_L4_markers',
            PIPELINE_KEY, SCANVI_LABELS_KEY,
            'scanvi_pred_L3', 'scanvi_pred_prob']:
    if col in adata_hvg.obs.columns:
        n_cat = adata_hvg.obs[col].nunique()
        print(f"  {col:<35}: {n_cat} categories")

print(f"\n.obsm keys:")
for key in adata_hvg.obsm.keys():
    print(f"  {key}: {adata_hvg.obsm[key].shape}")

print(f"\nscArches query checklist (for future query scripts):")
print(f"  1. Load var_names.csv  -> {Path(scvi_model_dir) / 'var_names.csv'}")
print(f"  2. Compute covariates BEFORE gene alignment (QRM 15.4)")
print(f"  3. Save full-gene .raw (log1p) BEFORE HVG subset (QRM 15.5)")
print(f"  4. Use COO scatter for gene alignment (QRM 15.6)")
print(f"  5. SCVI.load(scvi_dir, extend_categories=True)")
print(f"     -> SCVI.prepare_query_anndata() -> SCVI.load_query_data()")
print(f"  6. SCANVI.load(scanvi_dir, extend_categories=True)")
print(f"     -> SCANVI.prepare_query_anndata() -> SCANVI.load_query_data()")
print(f"     (NOT from_scvi_model -- QRM 15.1)")
print(f"  7. Add obs column '{SCANVI_LABELS_KEY}' = '{UNLABELED_CATEGORY}' before SCANVI.load (QRM 15.2)")
print(f"  8. Expected label space: {scanvi_label_categories}")
print(f"  9. Multi-subset query: merge ALL subsets into ONE object before mapping (QRM 15.3)")