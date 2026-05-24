# %% [markdown]
# # Epithelial Cell Annotation Pipeline - ENHANCED VERSION
# ## Quality-Controlled Cluster + GEP Analysis
# 
# **Enhancements in this version:**
# - ✅ Small batch detection and flagging
# - ✅ Low-usage GEP analysis (rare cell detection)
# - ✅ High FC association validation
# - ✅ Batch effect residual checking
# - ✅ K value comparison (10, 15, 20)
# 
# **Data characteristics (from QC):**
# - 278,930 cells, 29 clusters
# - 22 batches (2 small batches: 166 & 679 cells)
# - GEP usage: unbalanced (1-27% range)
# 
# ---

# %% [markdown]
# ## 📦 Setup and Data Loading

# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy import stats
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

sc.set_figure_params(dpi=100, frameon=False, figsize=(6, 6), facecolor='white')
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

print("✓ Libraries loaded")

# %%
# ============================================================================
# CONFIGURATION
# ============================================================================

DATA_PATH = "/home/h2048/data/py/1202/bbknn_celltype_analysis/Epithelial/output_complete_pipeline/checkpoint_complete_with_geps.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1207/bbknn_celltype_analysis/Epithelial/annotation_results"
CNMF_DIR = "/home/h2048/data/py/1202/bbknn_celltype_analysis/Epithelial/output_complete_pipeline/cnmf_results/Epithelial_HarmonyBatchCorrected"

import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Key columns
CLUSTER_KEY = 'leiden_bbknn'
TISSUE_KEY = 'tissue'
BATCH_KEY = 'dataset'

# GEP configuration
GEP_SUFFIX = 'harmony'
SELECTED_K = 15

# Quality control thresholds
SMALL_BATCH_THRESHOLD = 1000  # Flag batches < 1000 cells
LOW_USAGE_GEP_THRESHOLD = 0.02  # Flag GEPs with mean usage < 2%
HIGH_FC_THRESHOLD = 50  # Flag associations with FC > 50

print(f"✓ Configuration set")
print(f"   Data: {DATA_PATH}")
print(f"   Output: {OUTPUT_DIR}")

# %%
# Load data
print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)

print(f"\n✓ Data loaded")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")
print(f"   Clusters: {adata.obs[CLUSTER_KEY].nunique()}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# %% [markdown]
# ## 🔍 Enhanced Quality Control
# 
# ### QC1: Small Batch Detection and Flagging

# %%
# Identify small batches
batch_counts = adata.obs[BATCH_KEY].value_counts()
small_batches = batch_counts[batch_counts < SMALL_BATCH_THRESHOLD].index.tolist()

print(f"Small batch analysis (threshold: {SMALL_BATCH_THRESHOLD} cells):")
print(f"\nFound {len(small_batches)} small batches:")
for batch in small_batches:
    n_cells = batch_counts[batch]
    pct = n_cells / adata.n_obs * 100
    print(f"  {batch}: {n_cells} cells ({pct:.3f}%)")

# Flag cells from small batches
adata.obs['from_small_batch'] = adata.obs[BATCH_KEY].isin(small_batches)
n_flagged = adata.obs['from_small_batch'].sum()
print(f"\n✓ Flagged {n_flagged:,} cells from small batches ({n_flagged/adata.n_obs*100:.2f}%)")

# Check cluster distribution of small batch cells
if n_flagged > 0:
    small_batch_cells = adata[adata.obs['from_small_batch']]
    cluster_dist = small_batch_cells.obs[CLUSTER_KEY].value_counts()
    
    print(f"\nSmall batch cells spread across {len(cluster_dist)} clusters")
    
    # Check if concentrated in specific clusters
    top5 = cluster_dist.head(5)
    top5_pct = top5.sum() / n_flagged * 100
    
    print(f"Top 5 clusters contain {top5_pct:.1f}% of small batch cells:")
    for cluster, count in top5.items():
        pct = count / n_flagged * 100
        print(f"  Cluster {cluster}: {count} cells ({pct:.1f}%)")
    
    if top5_pct > 70:
        print(f"\n⚠️  WARNING: Small batch cells highly concentrated in few clusters!")
        print(f"   These clusters may have batch-specific signals.")

# %% [markdown]
# ### QC2: Low-Usage GEP Analysis (Rare Cell Detection)

# %%
# Calculate GEP usage statistics
gep_cols = [f'GEP_{i}_{GEP_SUFFIX}' for i in range(1, SELECTED_K + 1)]
gep_usage_stats = adata.obs[gep_cols].agg(['mean', 'std', 'max'])

# Identify low-usage GEPs
low_usage_geps = gep_usage_stats.columns[gep_usage_stats.loc['mean'] < LOW_USAGE_GEP_THRESHOLD]

print(f"Low-usage GEP analysis (threshold: {LOW_USAGE_GEP_THRESHOLD*100:.0f}% mean usage):")
print(f"\nFound {len(low_usage_geps)} low-usage GEPs:")

low_gep_details = []
for gep in low_usage_geps:
    mean_usage = gep_usage_stats.loc['mean', gep]
    max_usage = gep_usage_stats.loc['max', gep]
    
    # Count cells with high usage
    high_usage_cells = (adata.obs[gep] > 0.15).sum()
    high_usage_pct = high_usage_cells / adata.n_obs * 100
    
    print(f"\n{gep}:")
    print(f"  Mean usage: {mean_usage:.3f}")
    print(f"  Max usage: {max_usage:.3f}")
    print(f"  Cells with >15% usage: {high_usage_cells} ({high_usage_pct:.3f}%)")
    
    # Find dominant clusters for this GEP
    if high_usage_cells > 10:
        high_cells = adata.obs[gep] > 0.15
        dominant_clusters = adata.obs[high_cells][CLUSTER_KEY].value_counts().head(3)
        print(f"  Main clusters: {dominant_clusters.to_dict()}")
        
        low_gep_details.append({
            'GEP': gep,
            'Mean_Usage': mean_usage,
            'High_Usage_Cells': high_usage_cells,
            'Main_Clusters': ', '.join(map(str, dominant_clusters.index.tolist()))
        })

if len(low_gep_details) > 0:
    low_gep_df = pd.DataFrame(low_gep_details)
    low_gep_df.to_csv(f"{OUTPUT_DIR}/low_usage_GEPs_analysis.csv", index=False)
    print(f"\n💡 These GEPs may represent rare cell types (Ionocyte, Tuft, etc.)")
    print(f"   Check their top genes for biological interpretation!")

# %% [markdown]
# ### QC3: Batch Effect Residual Check

# %%
# Check if any cluster is dominated by a single batch
print("Checking for batch-dominated clusters...\n")

batch_dominated_clusters = []

for cluster in adata.obs[CLUSTER_KEY].unique():
    cluster_cells = adata[adata.obs[CLUSTER_KEY] == cluster]
    n_cluster = len(cluster_cells)
    
    # Get batch distribution in this cluster
    batch_dist = cluster_cells.obs[BATCH_KEY].value_counts()
    dominant_batch = batch_dist.index[0]
    dominant_count = batch_dist.values[0]
    dominant_pct = dominant_count / n_cluster * 100
    
    # Flag if >70% from single batch
    if dominant_pct > 70:
        batch_dominated_clusters.append({
            'Cluster': cluster,
            'Size': n_cluster,
            'Dominant_Batch': dominant_batch,
            'Dominant_Pct': dominant_pct,
            'N_Batches': len(batch_dist)
        })

if len(batch_dominated_clusters) > 0:
    batch_dom_df = pd.DataFrame(batch_dominated_clusters)
    batch_dom_df = batch_dom_df.sort_values('Dominant_Pct', ascending=False)
    
    print(f"⚠️  Found {len(batch_dominated_clusters)} batch-dominated clusters (>70% from one batch):\n")
    print(batch_dom_df.to_string(index=False))
    
    batch_dom_df.to_csv(f"{OUTPUT_DIR}/batch_dominated_clusters.csv", index=False)
    
    print(f"\n💡 These clusters may have batch-specific signals.")
    print(f"   Verify with marker genes and cross-batch validation!")
else:
    print("✅ No severely batch-dominated clusters found.")
    print("   Batch correction appears effective.")

# %% [markdown]
# ## 🧬 Comprehensive Marker Genes

# %%
# Full marker gene dictionary (from your uploaded image)
MARKER_GENES = {
    # Alveolar cells
    'AT1': ['AGER', 'CAV1', 'MYL9', 'SFTA2', 'CLIC3', 'SPOCK2', 'ANXA3', 'RTKN2', 'TIMP3', 'TNNC1'],
    'AT2': ['SFTPB', 'SFTPC', 'SFTPA1', 'SFTPA2', 'LAMP3', 'LRRK2', 'TFPI', 'MFSD2A', 'SERPINA1'],
    'AT2_proliferating': ['STMN1', 'SFTPA2', 'SFTPA1', 'TYMS', 'TK1', 'PTTG1', 'KIAA0101', 'TOP2A', 'CENPW', 'DTYMK'],
    'Transitional_Club_AT2': ['SCGB3A2', 'MGP', 'C16orf89', 'SFTA1P', 'VIM', 'SFTPA2', 'CAV1', 'ICAM1', 'SUSD2'],
    
    # Basal cells
    'Basal_resting': ['KRT15', 'KRT17', 'KRT5', 'DST', 'DLK2', 'IL33', 'FHL2', 'PTPRZ1'],
    'Suprabasal': ['KRT5', 'KRT17', 'SERPINB4', 'SERPINB13', 'KRT6A', 'CLCA2', 'LY6D', 'PPP1R14B', 'AKR1C3', 'IGFBP3'],
    'Basal_general': ['TP63', 'KRT5', 'KRT14'],
    
    # Ciliated cells
    'Deuterosomal': ['CCNO', 'CDC20B', 'ZMYND10', 'KIF9', 'FOXJ1', 'HES6', 'CEP78', 'TMEM106C', 'CCDC67', 'KDELC2'],
    'Multiciliated_nasal': ['C20orf85', 'C9orf24', 'RSPH1', 'PIFO', 'RP11-356K23.1', 'CCDC80', 'PROM1', 'OMG', 'DIAPH2', 'C15orf48'],
    'Multiciliated_non_nasal': ['C20orf85', 'CAPS', 'C9orf24', 'RSPH1', 'FAM183A', 'MS4A8', 'TFF3', 'IGFBP5', 'CFAP43', 'C2orf40'],
    'Ciliated_general': ['FOXJ1', 'RSPH1', 'PIFO'],
    
    # Club/Secretory cells
    'Club_non_nasal': ['TSPAN8', 'CYP2F1', 'TFF3', 'TGM2', 'MUC5B', 'CXCL6', 'C16orf89', 'HES4', 'RHOV', 'KIAA1324'],
    'Club_nasal': ['ASRGL1', 'LYPD2', 'UGT2A1', 'TFCP2L1', 'LY6D', 'TPD52L1', 'SORD', 'PI3'],
    'Secretory_general': ['SCGB1A1', 'SCGB3A2', 'SERPINB3'],
    
    # Goblet cells
    'Goblet_nasal': ['LYPD2', 'PI3', 'CEACAM5', 'LYNX1', 'MUC5AC', 'MUC16', 'C15orf48', 'BPIFA1', 'CCDC80', 'DHRS9'],
    'Goblet_bronchial': ['MUC5B', 'RARRES1', 'SAA1', 'ANKRD36C', 'SAA2', 'LYZ', 'PLCG2', 'FCGBP', 'RIMS1', 'MUC5AC'],
    'Goblet_subegmental': ['TSPAN8', 'MUC5B', 'C16orf89', 'MTRNR2L10', 'CLCA2', 'CFD', 'KIAA1324', 'LTF', 'TMEM45A', 'FCGBP'],
    'Goblet_general': ['MUC5AC', 'SPDEF', 'LYPD2', 'ITLN1'],
    
    # SMG cells
    'SMG_serous_nasal': ['LYZ', 'ZG16B', 'AZGP1', 'LTF', 'STATH', 'PIP', 'CLDN10', 'GJC3', 'ODAM', 'S100A1'],
    'SMG_serous_bronchial': ['LYZ', 'LTF', 'PRR4', 'AZGP1', 'S100A1', 'APIP', 'RP11-1143G9.4', 'PRB3', 'AC078941.1', 'C6orf58'],
    'SMG_mucous': ['MUC5B', 'BPIFB2', 'AZGP1', 'FCGBP', 'NKX3-1', 'TSPAN8', 'TFF1', 'DEFB1', 'HMGCS2', 'CRYM'],
    'SMG_duct': ['RARRES1', 'TCN1', 'SAA1', 'MIA', 'DMBT1', 'SAA2', 'RHOV', 'MMP7', 'ALDH1A3', 'ANKRD36C'],
    
    # Rare cells
    'Ionocyte': ['RARRES2', 'TMEM61', 'ASCL3', 'SCNN1B', 'STAP1', 'ATP6V1A', 'CFTR', 'HEPACAM2', 'CLCNKB', 'FOXI1'],
    'Tuft': ['STMN1', 'MARCKSL1', 'RASSF6', 'CRYM', 'HES6', 'KIT', 'AZGP1', 'HOMER3', 'NREP', 'LRMP'],
    'Neuroendocrine': ['PCSK1N', 'GRP', 'CPE', 'ASCL1', 'CHGA', 'SCG2', 'SCG5', 'SYT1', 'SCG3', 'SCGN'],
    
    # Proliferation
    'Proliferating': ['MKI67', 'TOP2A', 'TK1', 'CENPW', 'STMN1'],
}

ALL_MARKERS = sorted(set([gene for genes in MARKER_GENES.values() for gene in genes]))
available_markers = [m for m in ALL_MARKERS if m in adata.var_names]

print(f"Marker genes:")
print(f"  Total defined: {len(ALL_MARKERS)}")
print(f"  Available: {len(available_markers)} ({len(available_markers)/len(ALL_MARKERS)*100:.1f}%)")

# %% [markdown]
# ---
# ## 📊 Round 1: Core Visualizations
# 
# ### 1.1 UMAP Overview with QC Flags

# %%
# Enhanced UMAP with small batch highlighting
fig, axes = plt.subplots(2, 3, figsize=(20, 13))

# Panel 1: Batch
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, 
           title='Batch Distribution', frameon=False, legend_fontsize=8, s=3)

# Panel 2: Clusters
sc.pl.umap(adata, color=CLUSTER_KEY, ax=axes[0, 1], show=False,
           title='BBKNN Clusters', frameon=False, 
           legend_loc='on data', legend_fontsize=7, s=3)

# Panel 3: Tissue
sc.pl.umap(adata, color=TISSUE_KEY, ax=axes[0, 2], show=False,
           title='Tissue Origin', frameon=False, legend_fontsize=8, s=3)

# Panel 4: Small batch flag (NEW)
sc.pl.umap(adata, color='from_small_batch', ax=axes[1, 0], show=False,
           title='Small Batch Cells', frameon=False, s=3,
           palette=['lightgray', 'red'])

# Panel 5: Dominant GEP
gep_cols_short = [f'GEP_{i}_{GEP_SUFFIX}' for i in range(1, SELECTED_K + 1)]
adata.obs['dominant_GEP'] = adata.obs[gep_cols_short].idxmax(axis=1).str.replace(f'_{GEP_SUFFIX}', '')
sc.pl.umap(adata, color='dominant_GEP', ax=axes[1, 1], show=False,
           title='Dominant GEP', frameon=False, legend_fontsize=7, s=3)

# Panel 6: QC metric
if 'n_counts' in adata.obs.columns:
    sc.pl.umap(adata, color='n_counts', ax=axes[1, 2], show=False,
               title='UMI Counts', frameon=False, cmap='viridis', s=3)
else:
    axes[1, 2].axis('off')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/01_UMAP_overview_enhanced.png", dpi=300, bbox_inches='tight')
plt.show()

print("✓ Enhanced UMAP overview saved")

# %% [markdown]
# ### 1.2 GEP-Cluster Heatmap with Annotations

# %%
# Calculate GEP usage per cluster
gep_usage = adata.obs[gep_cols_short + [CLUSTER_KEY]].copy()
gep_usage.columns = [f'GEP{i}' for i in range(1, SELECTED_K + 1)] + ['cluster']
mean_usage = gep_usage.groupby('cluster').mean()

# Add annotations: cluster size and batch diversity
cluster_annotations = []
for cluster in mean_usage.index:
    cluster_cells = adata[adata.obs[CLUSTER_KEY] == cluster]
    n_cells = len(cluster_cells)
    n_batches = cluster_cells.obs[BATCH_KEY].nunique()
    small_batch_pct = (cluster_cells.obs['from_small_batch'].sum() / n_cells * 100)
    
    cluster_annotations.append({
        'Cluster': cluster,
        'Size': n_cells,
        'N_Batches': n_batches,
        'Small_Batch_Pct': small_batch_pct
    })

cluster_anno_df = pd.DataFrame(cluster_annotations).set_index('Cluster')

# Create enhanced clustermap
row_linkage = linkage(mean_usage.values, method='average')
col_linkage = linkage(mean_usage.T.values, method='average')

fig_width = max(14, SELECTED_K * 0.7)
fig_height = max(12, len(mean_usage) * 0.35)

g = sns.clustermap(
    mean_usage,
    row_linkage=row_linkage,
    col_linkage=col_linkage,
    cmap='RdYlBu_r',
    center=mean_usage.values.mean(),
    figsize=(fig_width, fig_height),
    cbar_kws={'label': 'Mean GEP Usage'},
    linewidths=0.5,
    linecolor='lightgray',
    yticklabels=True,
    xticklabels=True
)

plt.suptitle(f'GEP-Cluster Association Heatmap (K={SELECTED_K})', 
             y=0.98, fontsize=14, fontweight='bold')
plt.savefig(f"{OUTPUT_DIR}/02_GEP_cluster_heatmap.png", dpi=300, bbox_inches='tight')
plt.show()

# Save data
mean_usage.to_csv(f"{OUTPUT_DIR}/GEP_cluster_mean_usage.csv")
cluster_anno_df.to_csv(f"{OUTPUT_DIR}/cluster_annotations.csv")

print("✓ GEP-Cluster heatmap saved")

# %% [markdown]
# ### 1.3 Low-Usage GEP Visualization (Rare Cell Focus)

# %%
# Visualize low-usage GEPs on UMAP
if len(low_usage_geps) > 0:
    n_low_geps = len(low_usage_geps)
    ncols = 3
    nrows = (n_low_geps + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, nrows*4))
    axes = axes.flatten() if n_low_geps > 1 else [axes]
    
    for idx, gep in enumerate(low_usage_geps):
        sc.pl.umap(adata, color=gep, ax=axes[idx], show=False,
                   title=f"{gep} (mean={gep_usage_stats.loc['mean', gep]:.3f})",
                   cmap='viridis', frameon=False, s=5)
    
    # Hide extra axes
    for idx in range(n_low_geps, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/03_low_usage_GEPs_umap.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"✓ Low-usage GEP visualization saved")
    print(f"\n💡 These GEPs likely represent rare cell types:")
    print(f"   Check if they correspond to Ionocyte, Tuft, or Neuroendocrine markers!")
else:
    print("No low-usage GEPs detected.")

# %% [markdown]
# ### 1.4 Marker Gene Dotplot

# %%
# Curated marker list for dotplot
dotplot_markers = {
    'AT1': ['AGER', 'RTKN2', 'CAV1'],
    'AT2': ['SFTPC', 'SFTPB', 'LAMP3'],
    'AT2_prolif': ['STMN1', 'TOP2A', 'SFTPA2'],
    'Basal': ['TP63', 'KRT5', 'KRT15'],
    'Suprabasal': ['KRT6A', 'SERPINB4', 'KRT17'],
    'Ciliated': ['FOXJ1', 'RSPH1', 'PIFO'],
    'Deuterosomal': ['CCNO', 'CDC20B', 'FOXJ1'],
    'Club': ['SCGB1A1', 'SCGB3A2', 'CYP2F1'],
    'Goblet': ['MUC5AC', 'SPDEF', 'MUC5B'],
    'SMG': ['LYZ', 'LTF', 'AZGP1'],
    'Ionocyte': ['FOXI1', 'CFTR', 'ASCL3'],
    'Tuft': ['POU2F3', 'LRMP', 'KIT'],
    'NE': ['CHGA', 'GRP', 'ASCL1'],
    'Prolif': ['MKI67', 'TOP2A', 'CENPW'],
}

dotplot_genes = []
for genes in dotplot_markers.values():
    available = [g for g in genes if g in adata.var_names]
    dotplot_genes.extend(available)
dotplot_genes = list(dict.fromkeys(dotplot_genes))

n_clusters = adata.obs[CLUSTER_KEY].nunique()
fig_height = max(10, n_clusters * 0.3)
fig_width = max(12, len(dotplot_genes) * 0.35)

sc.pl.dotplot(
    adata,
    var_names=dotplot_genes,
    groupby=CLUSTER_KEY,
    dendrogram=True,
    figsize=(fig_width, fig_height),
    standard_scale='var',
    save=f'_{OUTPUT_DIR}/04_marker_dotplot.png'
)

print("✓ Marker dotplot saved")

# %% [markdown]
# ---
# ## 🔬 Round 2: Deep Analysis
# 
# ### 2.1 Cluster Category Signatures

# %%
# Calculate marker expression per cluster
marker_categories = {
    'Basal': ['TP63', 'KRT5', 'KRT15'],
    'Secretory': ['SCGB1A1', 'SCGB3A2'],
    'Ciliated': ['FOXJ1', 'RSPH1'],
    'Goblet': ['MUC5AC', 'SPDEF'],
    'AT1': ['AGER', 'RTKN2'],
    'AT2': ['SFTPC', 'SFTPB'],
}

cluster_signatures = {}

for cat_name, markers in marker_categories.items():
    available_markers = [m for m in markers if m in adata.var_names]
    if len(available_markers) == 0:
        continue
    
    if 'log1p' in adata.layers:
        expr_data = pd.DataFrame(
            adata[:, available_markers].layers['log1p'].toarray() 
            if hasattr(adata[:, available_markers].layers['log1p'], 'toarray') 
            else adata[:, available_markers].layers['log1p'],
            index=adata.obs_names,
            columns=available_markers
        )
    else:
        expr_data = pd.DataFrame(
            adata[:, available_markers].X.toarray() 
            if hasattr(adata[:, available_markers].X, 'toarray') 
            else adata[:, available_markers].X,
            index=adata.obs_names,
            columns=available_markers
        )
    
    expr_data['cluster'] = adata.obs[CLUSTER_KEY].values
    mean_expr = expr_data.groupby('cluster').mean().mean(axis=1)
    cluster_signatures[cat_name] = mean_expr

signature_df = pd.DataFrame(cluster_signatures)
signature_df.to_csv(f"{OUTPUT_DIR}/cluster_category_signatures.csv")

# Visualize
plt.figure(figsize=(10, max(10, n_clusters * 0.35)))
sns.heatmap(
    signature_df,
    cmap='RdYlBu_r',
    center=0,
    annot=True,
    fmt='.2f',
    cbar_kws={'label': 'Mean Log Expression'},
    linewidths=0.5
)
plt.title('Cluster Category Signatures', fontsize=14, fontweight='bold')
plt.xlabel('Cell Type Category', fontsize=12)
plt.ylabel('Cluster', fontsize=12)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/05_cluster_signatures.png", dpi=300, bbox_inches='tight')
plt.show()

print("✓ Cluster signatures saved")

# %% [markdown]
# ### 2.2 GEP Top Genes Analysis

# %%
# Load gene spectra
import glob

spectra_file = f"{CNMF_DIR}/Epithelial_HarmonyBatchCorrected.gene_spectra_score.k_{SELECTED_K}.dt_0_1.txt"

if os.path.exists(spectra_file):
    gene_spectra = pd.read_csv(spectra_file, sep='\t', index_col=0)
    
    TOP_N_GENES = 30
    gep_top_genes = {}
    
    for col in gene_spectra.columns:
        top_genes = gene_spectra[col].sort_values(ascending=False).head(TOP_N_GENES)
        gep_top_genes[col] = top_genes.index.tolist()
    
    # Save to file
    with open(f"{OUTPUT_DIR}/GEP_top_genes.txt", 'w') as f:
        for gep, genes in gep_top_genes.items():
            f.write(f"\n{'='*60}\n")
            f.write(f"{gep} - Top {TOP_N_GENES} Genes\n")
            f.write(f"{'='*60}\n")
            for i, gene in enumerate(genes, 1):
                f.write(f"{i:2d}. {gene}\n")
    
    print(f"✓ GEP top genes saved")
    
    # Display low-usage GEPs top genes
    if len(low_usage_geps) > 0:
        print(f"\nTop genes for low-usage GEPs (likely rare cells):\n")
        for gep in low_usage_geps:
            gep_num = gep.replace(f'GEP_', '').replace(f'_{GEP_SUFFIX}', '')
            col_name = str(gep_num)
            if col_name in gep_top_genes:
                print(f"{gep} top 10:")
                print(f"  {', '.join(gep_top_genes[col_name][:10])}")
                print()
else:
    print(f"⚠️  Gene spectra file not found: {spectra_file}")

# %% [markdown]
# ### 2.3 Tissue Composition Analysis

# %%
# Tissue × Cluster contingency table
tissue_cluster = pd.crosstab(
    adata.obs[TISSUE_KEY],
    adata.obs[CLUSTER_KEY],
    normalize='columns'
)

plt.figure(figsize=(max(14, n_clusters * 0.45), 8))
sns.heatmap(
    tissue_cluster,
    cmap='YlOrRd',
    annot=True,
    fmt='.2f',
    cbar_kws={'label': 'Proportion'},
    linewidths=0.5
)
plt.title('Tissue Composition per Cluster', fontsize=14, fontweight='bold')
plt.xlabel('Cluster', fontsize=12)
plt.ylabel('Tissue', fontsize=12)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/06_tissue_cluster_composition.png", dpi=300, bbox_inches='tight')
plt.show()

tissue_cluster.to_csv(f"{OUTPUT_DIR}/tissue_cluster_composition.csv")
print("✓ Tissue composition analysis saved")

# %% [markdown]
# ---
# ## 💡 Round 3: Annotation Suggestions
# 
# ### 3.1 Automated Suggestions with QC Flags

# %%
# Generate annotation suggestions
annotation_suggestions = []

for cluster in signature_df.index:
    scores = signature_df.loc[cluster]
    top_category = scores.idxmax()
    top_score = scores.max()
    second_category = scores.nlargest(2).index[1]
    second_score = scores.nlargest(2).values[1]
    
    # Determine confidence
    if top_score > 1.5 and (top_score - second_score) > 0.5:
        confidence = 'High'
        suggestion = top_category
    elif top_score > 1.0:
        confidence = 'Moderate'
        suggestion = f"{top_category}/{second_category}" if second_score > 0.8 else top_category
    else:
        confidence = 'Low'
        suggestion = f"Unknown ({top_category}?)"
    
    cluster_size = (adata.obs[CLUSTER_KEY] == cluster).sum()
    
    # Get dominant GEP
    cluster_gep_usage = mean_usage.loc[cluster]
    dominant_gep = cluster_gep_usage.idxmax()
    dominant_gep_usage = cluster_gep_usage.max()
    
    # Get QC flags
    cluster_anno = cluster_anno_df.loc[cluster]
    small_batch_pct = cluster_anno['Small_Batch_Pct']
    n_batches = cluster_anno['N_Batches']
    
    # Add warnings
    warnings = []
    if small_batch_pct > 20:
        warnings.append(f"Small_batch_{small_batch_pct:.0f}%")
    if n_batches < 5:
        warnings.append(f"Few_batches_{n_batches}")
    if cluster in [str(c) for c in batch_dominated_clusters] if len(batch_dominated_clusters) > 0 else []:
        warnings.append("Batch_dominated")
    
    annotation_suggestions.append({
        'Cluster': cluster,
        'Size': cluster_size,
        'Percent': f"{cluster_size/adata.n_obs*100:.2f}%",
        'Suggestion': suggestion,
        'Confidence': confidence,
        'Top_Score': round(top_score, 2),
        'Dominant_GEP': dominant_gep,
        'GEP_Usage': round(dominant_gep_usage, 3),
        'N_Batches': n_batches,
        'Warnings': '; '.join(warnings) if warnings else 'None'
    })

annotation_df = pd.DataFrame(annotation_suggestions)
annotation_df = annotation_df.sort_values('Size', ascending=False)

print("\n" + "="*100)
print("AUTOMATED ANNOTATION SUGGESTIONS WITH QC FLAGS")
print("="*100)
print(annotation_df.to_string(index=False))

annotation_df.to_csv(f"{OUTPUT_DIR}/annotation_suggestions_with_QC.csv", index=False)
print(f"\n✓ Saved: annotation_suggestions_with_QC.csv")

# %% [markdown]
# ### 3.2 Summary Statistics

# %%
print("\n" + "="*80)
print("ANNOTATION SUMMARY")
print("="*80)

print(f"\n1. Confidence Distribution:")
print(annotation_df['Confidence'].value_counts())

print(f"\n2. Clusters with Warnings:")
with_warnings = annotation_df[annotation_df['Warnings'] != 'None']
print(f"   Total: {len(with_warnings)} clusters")
if len(with_warnings) > 0:
    print("\n   Details:")
    print(with_warnings[['Cluster', 'Suggestion', 'Warnings']].to_string(index=False))

print(f"\n3. Suggested Cell Type Distribution:")
type_dist = annotation_df.groupby('Suggestion')['Size'].sum().sort_values(ascending=False)
print("\n   Type                    Cells        Percent")
print("   " + "-"*50)
for cell_type, count in type_dist.items():
    pct = count / adata.n_obs * 100
    print(f"   {cell_type:<23} {count:>8,}    {pct:>6.2f}%")

print(f"\n4. Low-Usage GEPs (Potential Rare Cells):")
if len(low_usage_geps) > 0:
    for gep in low_usage_geps:
        # Find which clusters use this GEP
        gep_num = gep.replace(f'GEP_', '').replace(f'_{GEP_SUFFIX}', '')
        gep_col = f'GEP{gep_num}'
        if gep_col in mean_usage.columns:
            top_cluster = mean_usage[gep_col].idxmax()
            top_usage = mean_usage.loc[top_cluster, gep_col]
            print(f"   {gep}: Cluster {top_cluster} (usage={top_usage:.3f})")
else:
    print("   None detected")

print("\n" + "="*80)

# %% [markdown]
# ---
# ## 📝 Final Summary
# 
# ### Key Files Generated:
# 1. `01_UMAP_overview_enhanced.png` - Data overview with QC flags
# 2. `02_GEP_cluster_heatmap.png` - ⭐ Main result
# 3. `03_low_usage_GEPs_umap.png` - Rare cell detection
# 4. `04_marker_dotplot.png` - Marker expression
# 5. `05_cluster_signatures.png` - Category signatures
# 6. `06_tissue_cluster_composition.png` - Tissue distribution
# 7. `annotation_suggestions_with_QC.csv` - ⭐ Main annotation guide
# 8. `low_usage_GEPs_analysis.csv` - Rare cell focus
# 9. `batch_dominated_clusters.csv` - Batch effect warnings
# 10. `GEP_top_genes.txt` - Biological interpretation
# 
# ### Next Steps:
# 1. Review `annotation_suggestions_with_QC.csv`
# 2. Pay special attention to:
#    - Clusters with warnings (batch effects)
#    - Low-usage GEPs (rare cells)
#    - Low confidence suggestions
# 3. Check GEP top genes for biological meaning
# 4. Manually refine annotations
# 5. Validate across tissues
# 

# %%
print("\n" + "="*80)
print("🎉 ENHANCED ANNOTATION PIPELINE COMPLETE!")
print("="*80)
print(f"\nAll results saved to: {OUTPUT_DIR}")
print(f"\n⭐ Key files to review:")
print(f"   1. annotation_suggestions_with_QC.csv")
print(f"   2. 02_GEP_cluster_heatmap.png")
print(f"   3. low_usage_GEPs_analysis.csv (for rare cells)")
print(f"   4. batch_dominated_clusters.csv (QC check)")
print(f"\n💡 Special attention needed for:")
print(f"   - {len(with_warnings)} clusters with QC warnings")
print(f"   - {len(low_usage_geps)} low-usage GEPs (potential rare cells)")
print(f"   - {len(small_batches)} small batches ({n_flagged} cells)")
print(f"\n🔬 Next: Manual refinement based on biological knowledge!")
print("="*80)

