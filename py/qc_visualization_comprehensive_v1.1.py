# %% [markdown]
# # Comprehensive Quality Control Visualization
# ## Post-integration QC metrics visualization by cluster and sample
# 
# **Quality Control Criteria:**
# - nFeature_RNA: 200-6000
# - MT%: <20
# - RP%: <40
# - DoubletFinder: 0.06 threshold
# - DecontX: <0.25

# %%
# Import required libraries
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

# Set plotting parameters
sc.set_figure_params(dpi=100, dpi_save=300, frameon=False, figsize=(8, 6))
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

print("Scanpy version:", sc.__version__)
print("Setup complete!")

# %% [markdown]
# ## 1. Load Data and Configuration

# %%
# File path
data_path = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
output_dir = "/home/h2048/data/py/1203/qc_visualization_output"


# Create output directory
import os
os.makedirs(output_dir, exist_ok=True)

# Load data
print(f"Loading data from: {data_path}")
adata = sc.read_h5ad(data_path)

print(f"\nData loaded successfully!")
print(f"Total cells: {adata.n_obs:,}")
print(f"Total genes: {adata.n_vars:,}")
print(f"\nAvailable metadata columns:")
print(adata.obs.columns.tolist())

# %%
# Check QC metrics availability
qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rp']
alternative_names = {
    'n_genes_by_counts': ['nFeature_RNA', 'n_genes'],
    'total_counts': ['nCount_RNA', 'n_counts'],
    'pct_counts_mt': ['percent.mt', 'percent_mt', 'pct_mt'],
    'pct_counts_rp': ['percent.rp', 'percent_rp', 'pct_rp']
}

# Map metric names
metric_mapping = {}
for standard_name, alternatives in alternative_names.items():
    if standard_name in adata.obs.columns:
        metric_mapping[standard_name] = standard_name
    else:
        for alt in alternatives:
            if alt in adata.obs.columns:
                metric_mapping[standard_name] = alt
                break

print("\nQC metrics mapping:")
for std, actual in metric_mapping.items():
    print(f"  {std} -> {actual}")

# Check for doublet and contamination scores
doublet_cols = [col for col in adata.obs.columns if 'doublet' in col.lower()]
decontx_cols = [col for col in adata.obs.columns if 'decontx' in col.lower() or 'contamination' in col.lower()]

print(f"\nDoublet-related columns: {doublet_cols}")
print(f"DecontX-related columns: {decontx_cols}")

# %%
# Standardize column names for easier access
adata.obs['n_genes'] = adata.obs[metric_mapping.get('n_genes_by_counts', 'n_genes_by_counts')]
adata.obs['n_counts'] = adata.obs[metric_mapping.get('total_counts', 'total_counts')]
adata.obs['mt_pct'] = adata.obs[metric_mapping.get('pct_counts_mt', 'pct_counts_mt')]
adata.obs['rp_pct'] = adata.obs[metric_mapping.get('pct_counts_rp', 'pct_counts_rp')]

# Add log-transformed counts for visualization
adata.obs['log_n_genes'] = np.log10(adata.obs['n_genes'] + 1)
adata.obs['log_n_counts'] = np.log10(adata.obs['n_counts'] + 1)

print("Standardized column names created!")

# Display basic statistics
print("\n=== QC Metrics Summary ===")
summary_stats = adata.obs[['n_genes', 'n_counts', 'mt_pct', 'rp_pct']].describe()
print(summary_stats)

# %% [markdown]
# ## 2. Overall QC Distribution

# %%
# Create comprehensive QC overview
fig = plt.figure(figsize=(16, 12))
gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

# 1. nFeature distribution
ax1 = fig.add_subplot(gs[0, 0])
ax1.hist(adata.obs['n_genes'], bins=100, color='steelblue', alpha=0.7, edgecolor='black')
ax1.axvline(200, color='red', linestyle='--', linewidth=2, label='Min: 200')
ax1.axvline(6000, color='red', linestyle='--', linewidth=2, label='Max: 6000')
ax1.set_xlabel('Number of Genes', fontsize=12)
ax1.set_ylabel('Cell Count', fontsize=12)
ax1.set_title('Gene Count Distribution', fontsize=13, fontweight='bold')
ax1.legend(fontsize=10)

# 2. nCount distribution
ax2 = fig.add_subplot(gs[0, 1])
ax2.hist(adata.obs['n_counts'], bins=100, color='forestgreen', alpha=0.7, edgecolor='black')
ax2.set_xlabel('Total Counts (UMI)', fontsize=12)
ax2.set_ylabel('Cell Count', fontsize=12)
ax2.set_title('UMI Count Distribution', fontsize=13, fontweight='bold')

# 3. MT% distribution
ax3 = fig.add_subplot(gs[0, 2])
ax3.hist(adata.obs['mt_pct'], bins=100, color='salmon', alpha=0.7, edgecolor='black')
ax3.axvline(20, color='red', linestyle='--', linewidth=2, label='Threshold: 20%')
ax3.set_xlabel('Mitochondrial %', fontsize=12)
ax3.set_ylabel('Cell Count', fontsize=12)
ax3.set_title('MT% Distribution', fontsize=13, fontweight='bold')
ax3.legend(fontsize=10)

# 4. RP% distribution
ax4 = fig.add_subplot(gs[1, 0])
ax4.hist(adata.obs['rp_pct'], bins=100, color='mediumpurple', alpha=0.7, edgecolor='black')
ax4.axvline(40, color='red', linestyle='--', linewidth=2, label='Threshold: 40%')
ax4.set_xlabel('Ribosomal %', fontsize=12)
ax4.set_ylabel('Cell Count', fontsize=12)
ax4.set_title('RP% Distribution', fontsize=13, fontweight='bold')
ax4.legend(fontsize=10)

# 5. nGenes vs nCounts
ax5 = fig.add_subplot(gs[1, 1])
scatter = ax5.scatter(adata.obs['n_counts'], adata.obs['n_genes'], 
                     c=adata.obs['mt_pct'], cmap='RdYlBu_r', 
                     alpha=0.5, s=5, rasterized=True)
ax5.set_xlabel('Total Counts (UMI)', fontsize=12)
ax5.set_ylabel('Number of Genes', fontsize=12)
ax5.set_title('Counts vs Genes (colored by MT%)', fontsize=13, fontweight='bold')
cbar = plt.colorbar(scatter, ax=ax5)
cbar.set_label('MT%', fontsize=10)

# 6. nGenes vs MT%
ax6 = fig.add_subplot(gs[1, 2])
ax6.scatter(adata.obs['n_genes'], adata.obs['mt_pct'], 
           alpha=0.3, s=5, color='coral', rasterized=True)
ax6.axhline(20, color='red', linestyle='--', linewidth=2)
ax6.axvline(200, color='red', linestyle='--', linewidth=2)
ax6.axvline(6000, color='red', linestyle='--', linewidth=2)
ax6.set_xlabel('Number of Genes', fontsize=12)
ax6.set_ylabel('Mitochondrial %', fontsize=12)
ax6.set_title('Genes vs MT%', fontsize=13, fontweight='bold')

# 7. nGenes vs RP%
ax7 = fig.add_subplot(gs[2, 0])
ax7.scatter(adata.obs['n_genes'], adata.obs['rp_pct'], 
           alpha=0.3, s=5, color='mediumpurple', rasterized=True)
ax7.axhline(40, color='red', linestyle='--', linewidth=2)
ax7.axvline(200, color='red', linestyle='--', linewidth=2)
ax7.axvline(6000, color='red', linestyle='--', linewidth=2)
ax7.set_xlabel('Number of Genes', fontsize=12)
ax7.set_ylabel('Ribosomal %', fontsize=12)
ax7.set_title('Genes vs RP%', fontsize=13, fontweight='bold')

# 8. MT% vs RP%
ax8 = fig.add_subplot(gs[2, 1])
ax8.scatter(adata.obs['mt_pct'], adata.obs['rp_pct'], 
           alpha=0.3, s=5, color='teal', rasterized=True)
ax8.axhline(40, color='red', linestyle='--', linewidth=2)
ax8.axvline(20, color='red', linestyle='--', linewidth=2)
ax8.set_xlabel('Mitochondrial %', fontsize=12)
ax8.set_ylabel('Ribosomal %', fontsize=12)
ax8.set_title('MT% vs RP%', fontsize=13, fontweight='bold')

# 9. Summary statistics table
ax9 = fig.add_subplot(gs[2, 2])
ax9.axis('off')
summary_text = f"""
QC Summary Statistics
{'='*30}
Total Cells: {adata.n_obs:,}
Total Genes: {adata.n_vars:,}

Genes per Cell:
  Mean: {adata.obs['n_genes'].mean():.0f}
  Median: {adata.obs['n_genes'].median():.0f}

UMI per Cell:
  Mean: {adata.obs['n_counts'].mean():.0f}
  Median: {adata.obs['n_counts'].median():.0f}

MT% per Cell:
  Mean: {adata.obs['mt_pct'].mean():.2f}%
  Median: {adata.obs['mt_pct'].median():.2f}%

RP% per Cell:
  Mean: {adata.obs['rp_pct'].mean():.2f}%
  Median: {adata.obs['rp_pct'].median():.2f}%
"""
ax9.text(0.1, 0.5, summary_text, fontsize=11, fontfamily='monospace',
        verticalalignment='center')

plt.suptitle('Overall QC Metrics Distribution', fontsize=16, fontweight='bold', y=0.98)
plt.savefig(f"{output_dir}/00_overall_qc_distribution.pdf", bbox_inches='tight', dpi=300)
plt.savefig(f"{output_dir}/00_overall_qc_distribution.png", bbox_inches='tight', dpi=300)
plt.show()

print("Overall QC distribution plot saved!")

# %% [markdown]
# ## 3. QC Visualization by Cluster

# %%
# Check cluster column
cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower()]
print(f"Available cluster columns: {cluster_cols}")

# Use the most appropriate cluster column
if 'cluster' in adata.obs.columns:
    cluster_col = 'cluster'
elif 'clusters' in adata.obs.columns:
    cluster_col = 'clusters'
elif 'leiden' in adata.obs.columns:
    cluster_col = 'leiden'
elif len(cluster_cols) > 0:
    cluster_col = cluster_cols[0]
else:
    cluster_col = None
    print("Warning: No cluster column found!")

if cluster_col:
    print(f"\nUsing cluster column: {cluster_col}")
    print(f"Number of clusters: {adata.obs[cluster_col].nunique()}")
    print(f"Cluster distribution:")
    print(adata.obs[cluster_col].value_counts().sort_index())

# %%
# Violin plots by cluster
if cluster_col:
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # QC thresholds for reference lines
    thresholds = {
        'n_genes': [200, 6000],
        'mt_pct': [20],
        'rp_pct': [40]
    }
    
    # 1. Number of genes
    sc.pl.violin(adata, keys='n_genes', groupby=cluster_col, rotation=45,
                ax=axes[0, 0], show=False)
    axes[0, 0].axhline(200, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[0, 0].axhline(6000, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[0, 0].set_title('Number of Genes by Cluster', fontsize=13, fontweight='bold')
    axes[0, 0].set_ylabel('Number of Genes', fontsize=11)
    
    # 2. Total counts
    sc.pl.violin(adata, keys='n_counts', groupby=cluster_col, rotation=45,
                ax=axes[0, 1], show=False)
    axes[0, 1].set_title('Total UMI Counts by Cluster', fontsize=13, fontweight='bold')
    axes[0, 1].set_ylabel('Total Counts', fontsize=11)
    
    # 3. MT%
    sc.pl.violin(adata, keys='mt_pct', groupby=cluster_col, rotation=45,
                ax=axes[1, 0], show=False)
    axes[1, 0].axhline(20, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[1, 0].set_title('Mitochondrial % by Cluster', fontsize=13, fontweight='bold')
    axes[1, 0].set_ylabel('MT %', fontsize=11)
    
    # 4. RP%
    sc.pl.violin(adata, keys='rp_pct', groupby=cluster_col, rotation=45,
                ax=axes[1, 1], show=False)
    axes[1, 1].axhline(40, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[1, 1].set_title('Ribosomal % by Cluster', fontsize=13, fontweight='bold')
    axes[1, 1].set_ylabel('RP %', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/01_qc_violin_by_cluster.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/01_qc_violin_by_cluster.png", bbox_inches='tight', dpi=300)
    plt.show()
    
    print("QC violin plots by cluster saved!")

# %%
# Box plots by cluster with statistical summary
if cluster_col:
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    metrics = ['n_genes', 'n_counts', 'mt_pct', 'rp_pct']
    titles = ['Number of Genes', 'Total UMI Counts', 'Mitochondrial %', 'Ribosomal %']
    colors = ['steelblue', 'forestgreen', 'salmon', 'mediumpurple']
    
    for idx, (metric, title, color) in enumerate(zip(metrics, titles, colors)):
        ax = axes[idx // 2, idx % 2]
        
        # Create box plot
        adata.obs.boxplot(column=metric, by=cluster_col, ax=ax, 
                         patch_artist=True, showfliers=False)
        
        # Customize appearance
        ax.set_title(f'{title} by Cluster', fontsize=13, fontweight='bold')
        ax.set_xlabel('Cluster', fontsize=11)
        ax.set_ylabel(title, fontsize=11)
        ax.tick_params(axis='x', rotation=45)
        
        # Add threshold lines
        if metric == 'n_genes':
            ax.axhline(200, color='red', linestyle='--', linewidth=1, alpha=0.7)
            ax.axhline(6000, color='red', linestyle='--', linewidth=1, alpha=0.7)
        elif metric == 'mt_pct':
            ax.axhline(20, color='red', linestyle='--', linewidth=1, alpha=0.7)
        elif metric == 'rp_pct':
            ax.axhline(40, color='red', linestyle='--', linewidth=1, alpha=0.7)
    
    plt.suptitle('')  # Remove default title
    plt.tight_layout()
    plt.savefig(f"{output_dir}/02_qc_boxplot_by_cluster.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/02_qc_boxplot_by_cluster.png", bbox_inches='tight', dpi=300)
    plt.show()
    
    print("QC box plots by cluster saved!")

# %%
# Ridge plots by cluster (distribution density)
if cluster_col:
    import joypy
    from matplotlib import cm
    
    metrics = ['n_genes', 'mt_pct', 'rp_pct']
    titles = ['Gene Count Distribution', 'MT% Distribution', 'RP% Distribution']
    
    for metric, title in zip(metrics, titles):
        fig, axes = joypy.joyplot(
            adata.obs, 
            column=metric, 
            by=cluster_col,
            figsize=(10, 8),
            title=f"{title} by Cluster",
            colormap=cm.viridis,
            alpha=0.7,
            linewidth=1.5
        )
        
        plt.xlabel(title, fontsize=12)
        plt.savefig(f"{output_dir}/03_qc_ridge_{metric}_by_cluster.pdf", 
                   bbox_inches='tight', dpi=300)
        plt.show()
    
    print("QC ridge plots by cluster saved!")

# %% [markdown]
# ## 4. QC Visualization by Sample

# %%
# Check sample column
sample_cols = [col for col in adata.obs.columns if 'sample' in col.lower()]
print(f"Available sample columns: {sample_cols}")

# Use the most appropriate sample column
if 'sample' in adata.obs.columns:
    sample_col = 'sample'
elif 'sample_id' in adata.obs.columns:
    sample_col = 'sample_id'
elif 'orig.ident' in adata.obs.columns:
    sample_col = 'orig.ident'
elif len(sample_cols) > 0:
    sample_col = sample_cols[0]
else:
    sample_col = None
    print("Warning: No sample column found!")

if sample_col:
    print(f"\nUsing sample column: {sample_col}")
    print(f"Number of samples: {adata.obs[sample_col].nunique()}")
    print(f"\nSample distribution:")
    print(adata.obs[sample_col].value_counts())

# %%
# Violin plots by sample
if sample_col:
    n_samples = adata.obs[sample_col].nunique()
    
    # Adjust figure size based on number of samples
    figsize = (max(16, n_samples * 0.8), 12)
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    # 1. Number of genes
    sc.pl.violin(adata, keys='n_genes', groupby=sample_col, rotation=90,
                ax=axes[0, 0], show=False)
    axes[0, 0].axhline(200, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[0, 0].axhline(6000, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[0, 0].set_title('Number of Genes by Sample', fontsize=13, fontweight='bold')
    axes[0, 0].set_ylabel('Number of Genes', fontsize=11)
    
    # 2. Total counts
    sc.pl.violin(adata, keys='n_counts', groupby=sample_col, rotation=90,
                ax=axes[0, 1], show=False)
    axes[0, 1].set_title('Total UMI Counts by Sample', fontsize=13, fontweight='bold')
    axes[0, 1].set_ylabel('Total Counts', fontsize=11)
    
    # 3. MT%
    sc.pl.violin(adata, keys='mt_pct', groupby=sample_col, rotation=90,
                ax=axes[1, 0], show=False)
    axes[1, 0].axhline(20, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[1, 0].set_title('Mitochondrial % by Sample', fontsize=13, fontweight='bold')
    axes[1, 0].set_ylabel('MT %', fontsize=11)
    
    # 4. RP%
    sc.pl.violin(adata, keys='rp_pct', groupby=sample_col, rotation=90,
                ax=axes[1, 1], show=False)
    axes[1, 1].axhline(40, color='red', linestyle='--', linewidth=1, alpha=0.7)
    axes[1, 1].set_title('Ribosomal % by Sample', fontsize=13, fontweight='bold')
    axes[1, 1].set_ylabel('RP %', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/04_qc_violin_by_sample.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/04_qc_violin_by_sample.png", bbox_inches='tight', dpi=300)
    plt.show()
    
    print("QC violin plots by sample saved!")

# %%
# Heatmap of QC metrics by sample
if sample_col:
    # Calculate mean values for each sample
    qc_summary = adata.obs.groupby(sample_col).agg({
        'n_genes': 'mean',
        'n_counts': 'mean',
        'mt_pct': 'mean',
        'rp_pct': 'mean'
    }).round(2)
    
    # Normalize for better visualization
    qc_summary_norm = (qc_summary - qc_summary.min()) / (qc_summary.max() - qc_summary.min())
    
    # Create heatmap
    fig, axes = plt.subplots(1, 2, figsize=(16, max(8, len(qc_summary) * 0.4)))
    
    # Raw values heatmap
    sns.heatmap(qc_summary, annot=True, fmt='.1f', cmap='YlOrRd', 
               cbar_kws={'label': 'Mean Value'}, ax=axes[0], linewidths=0.5)
    axes[0].set_title('Mean QC Metrics by Sample (Raw Values)', 
                     fontsize=13, fontweight='bold')
    axes[0].set_xlabel('QC Metric', fontsize=11)
    axes[0].set_ylabel('Sample', fontsize=11)
    
    # Normalized heatmap
    sns.heatmap(qc_summary_norm, annot=True, fmt='.2f', cmap='viridis', 
               cbar_kws={'label': 'Normalized Value'}, ax=axes[1], linewidths=0.5)
    axes[1].set_title('Mean QC Metrics by Sample (Normalized)', 
                     fontsize=13, fontweight='bold')
    axes[1].set_xlabel('QC Metric', fontsize=11)
    axes[1].set_ylabel('Sample', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/05_qc_heatmap_by_sample.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/05_qc_heatmap_by_sample.png", bbox_inches='tight', dpi=300)
    plt.show()
    
    print("QC heatmap by sample saved!")
    
    # Save summary table
    qc_summary.to_csv(f"{output_dir}/qc_summary_by_sample.csv")
    print(f"QC summary table saved to: {output_dir}/qc_summary_by_sample.csv")

# %% [markdown]
# ## 5. Cell Composition Analysis

# %%
# Cell count bar plots
if cluster_col and sample_col:
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 1. Cells per cluster
    cluster_counts = adata.obs[cluster_col].value_counts().sort_index()
    axes[0].bar(range(len(cluster_counts)), cluster_counts.values, 
               color='steelblue', alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('Cluster', fontsize=12)
    axes[0].set_ylabel('Cell Count', fontsize=12)
    axes[0].set_title('Cell Distribution by Cluster', fontsize=13, fontweight='bold')
    axes[0].set_xticks(range(len(cluster_counts)))
    axes[0].set_xticklabels(cluster_counts.index, rotation=45)
    axes[0].grid(axis='y', alpha=0.3)
    
    # Add count labels on bars
    for i, v in enumerate(cluster_counts.values):
        axes[0].text(i, v + max(cluster_counts.values) * 0.01, 
                    f'{v:,}', ha='center', fontsize=9)
    
    # 2. Cells per sample
    sample_counts = adata.obs[sample_col].value_counts()
    axes[1].barh(range(len(sample_counts)), sample_counts.values, 
                color='forestgreen', alpha=0.7, edgecolor='black')
    axes[1].set_ylabel('Sample', fontsize=12)
    axes[1].set_xlabel('Cell Count', fontsize=12)
    axes[1].set_title('Cell Distribution by Sample', fontsize=13, fontweight='bold')
    axes[1].set_yticks(range(len(sample_counts)))
    axes[1].set_yticklabels(sample_counts.index)
    axes[1].grid(axis='x', alpha=0.3)
    
    # Add count labels on bars
    for i, v in enumerate(sample_counts.values):
        axes[1].text(v + max(sample_counts.values) * 0.01, i, 
                    f'{v:,}', va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/06_cell_counts.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/06_cell_counts.png", bbox_inches='tight', dpi=300)
    plt.show()
    
    print("Cell count plots saved!")

# %%
# Stacked bar plot: cluster composition per sample
if cluster_col and sample_col:
    # Create crosstab
    composition = pd.crosstab(adata.obs[sample_col], adata.obs[cluster_col], 
                             normalize='index') * 100
    
    # Plot
    fig, ax = plt.subplots(figsize=(14, max(8, len(composition) * 0.4)))
    composition.plot(kind='barh', stacked=True, ax=ax, 
                    colormap='tab20', edgecolor='white', linewidth=0.5)
    
    ax.set_xlabel('Percentage (%)', fontsize=12)
    ax.set_ylabel('Sample', fontsize=12)
    ax.set_title('Cluster Composition by Sample', fontsize=13, fontweight='bold')
    ax.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left', 
             ncol=1, fontsize=9)
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/07_cluster_composition_by_sample.pdf", 
               bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/07_cluster_composition_by_sample.png", 
               bbox_inches='tight', dpi=300)
    plt.show()
    
    # Save composition table
    composition.to_csv(f"{output_dir}/cluster_composition_by_sample.csv")
    print("Cluster composition plot and table saved!")

# %% [markdown]
# ## 6. Doublet and Contamination Analysis (if available)

# %%
# Check for doublet information
if len(doublet_cols) > 0:
    doublet_col = doublet_cols[0]
    print(f"Using doublet column: {doublet_col}")
    
    # If it's a score column
    if adata.obs[doublet_col].dtype in ['float64', 'float32']:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Distribution
        axes[0].hist(adata.obs[doublet_col], bins=100, color='coral', 
                    alpha=0.7, edgecolor='black')
        axes[0].axvline(0.06, color='red', linestyle='--', linewidth=2, 
                       label='Threshold: 0.06')
        axes[0].set_xlabel('Doublet Score', fontsize=12)
        axes[0].set_ylabel('Cell Count', fontsize=12)
        axes[0].set_title('Doublet Score Distribution', fontsize=13, fontweight='bold')
        axes[0].legend(fontsize=10)
        
        # By cluster
        if cluster_col:
            sc.pl.violin(adata, keys=doublet_col, groupby=cluster_col, 
                        rotation=45, ax=axes[1], show=False)
            axes[1].axhline(0.06, color='red', linestyle='--', linewidth=1, alpha=0.7)
            axes[1].set_title('Doublet Score by Cluster', fontsize=13, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/08_doublet_analysis.pdf", bbox_inches='tight', dpi=300)
        plt.savefig(f"{output_dir}/08_doublet_analysis.png", bbox_inches='tight', dpi=300)
        plt.show()
        
        print("Doublet analysis plots saved!")
else:
    print("No doublet information found in dataset.")

# %%
# Check for contamination (DecontX) information
if len(decontx_cols) > 0:
    decontx_col = decontx_cols[0]
    print(f"Using contamination column: {decontx_col}")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Distribution
    axes[0].hist(adata.obs[decontx_col], bins=100, color='mediumpurple', 
                alpha=0.7, edgecolor='black')
    axes[0].axvline(0.25, color='red', linestyle='--', linewidth=2, 
                   label='Threshold: 0.25')
    axes[0].set_xlabel('Contamination Score', fontsize=12)
    axes[0].set_ylabel('Cell Count', fontsize=12)
    axes[0].set_title('DecontX Contamination Distribution', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=10)
    
    # By cluster
    if cluster_col:
        sc.pl.violin(adata, keys=decontx_col, groupby=cluster_col, 
                    rotation=45, ax=axes[1], show=False)
        axes[1].axhline(0.25, color='red', linestyle='--', linewidth=1, alpha=0.7)
        axes[1].set_title('Contamination Score by Cluster', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/09_contamination_analysis.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/09_contamination_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()
    
    print("Contamination analysis plots saved!")
else:
    print("No contamination information found in dataset.")

# %% [markdown]
# ## 7. Generate Summary Report

# %%
# Generate comprehensive summary statistics
summary_report = []
summary_report.append("="*80)
summary_report.append("QUALITY CONTROL SUMMARY REPORT")
summary_report.append("="*80)
summary_report.append("")

# Dataset overview
summary_report.append("[1] Dataset Overview")
summary_report.append("-" * 80)
summary_report.append(f"Total Cells: {adata.n_obs:,}")
summary_report.append(f"Total Genes: {adata.n_vars:,}")
if cluster_col:
    summary_report.append(f"Number of Clusters: {adata.obs[cluster_col].nunique()}")
if sample_col:
    summary_report.append(f"Number of Samples: {adata.obs[sample_col].nunique()}")
summary_report.append("")

# QC metrics statistics
summary_report.append("[2] QC Metrics Statistics")
summary_report.append("-" * 80)

metrics_info = [
    ('n_genes', 'Genes per Cell', [200, 6000]),
    ('n_counts', 'UMI Counts per Cell', None),
    ('mt_pct', 'Mitochondrial %', [20]),
    ('rp_pct', 'Ribosomal %', [40])
]

for metric, label, thresholds in metrics_info:
    values = adata.obs[metric]
    summary_report.append(f"\n{label}:")
    summary_report.append(f"  Mean: {values.mean():.2f}")
    summary_report.append(f"  Median: {values.median():.2f}")
    summary_report.append(f"  Std: {values.std():.2f}")
    summary_report.append(f"  Min: {values.min():.2f}")
    summary_report.append(f"  Max: {values.max():.2f}")
    
    if thresholds:
        if len(thresholds) == 2:
            passed = ((values >= thresholds[0]) & (values <= thresholds[1])).sum()
            summary_report.append(f"  Cells passing QC ({thresholds[0]}-{thresholds[1]}): "
                                f"{passed:,} ({passed/len(values)*100:.1f}%)")
        else:
            passed = (values < thresholds[0]).sum()
            summary_report.append(f"  Cells passing QC (<{thresholds[0]}): "
                                f"{passed:,} ({passed/len(values)*100:.1f}%)")

summary_report.append("")

# Cluster statistics
if cluster_col:
    summary_report.append("[3] Cluster Statistics")
    summary_report.append("-" * 80)
    cluster_stats = adata.obs.groupby(cluster_col).agg({
        'n_genes': ['count', 'mean', 'std'],
        'mt_pct': 'mean',
        'rp_pct': 'mean'
    }).round(2)
    summary_report.append(cluster_stats.to_string())
    summary_report.append("")

# Sample statistics
if sample_col:
    summary_report.append("[4] Sample Statistics")
    summary_report.append("-" * 80)
    sample_stats = adata.obs.groupby(sample_col).agg({
        'n_genes': ['count', 'mean', 'std'],
        'mt_pct': 'mean',
        'rp_pct': 'mean'
    }).round(2)
    summary_report.append(sample_stats.to_string())
    summary_report.append("")

summary_report.append("="*80)
summary_report.append(f"Report generated: {pd.Timestamp.now()}")
summary_report.append("="*80)

# Print to console
report_text = "\n".join(summary_report)
print(report_text)

# Save to file
with open(f"{output_dir}/QC_SUMMARY_REPORT.txt", 'w') as f:
    f.write(report_text)

print(f"\nSummary report saved to: {output_dir}/QC_SUMMARY_REPORT.txt")

# %% [markdown]
# ## 8. Export Summary Tables

# %%
# Export detailed statistics tables
export_dir = f"{output_dir}/tables"
os.makedirs(export_dir, exist_ok=True)

# 1. Overall statistics
overall_stats = adata.obs[['n_genes', 'n_counts', 'mt_pct', 'rp_pct']].describe()
overall_stats.to_csv(f"{export_dir}/01_overall_qc_statistics.csv")

# 2. Statistics by cluster
if cluster_col:
    cluster_stats = adata.obs.groupby(cluster_col)[[
        'n_genes', 'n_counts', 'mt_pct', 'rp_pct'
    ]].describe().T
    cluster_stats.to_csv(f"{export_dir}/02_qc_statistics_by_cluster.csv")

# 3. Statistics by sample
if sample_col:
    sample_stats = adata.obs.groupby(sample_col)[[
        'n_genes', 'n_counts', 'mt_pct', 'rp_pct'
    ]].describe().T
    sample_stats.to_csv(f"{export_dir}/03_qc_statistics_by_sample.csv")

# 4. Cell counts
if cluster_col:
    cluster_counts = adata.obs[cluster_col].value_counts().to_frame('cell_count')
    cluster_counts.to_csv(f"{export_dir}/04_cell_counts_by_cluster.csv")

if sample_col:
    sample_counts = adata.obs[sample_col].value_counts().to_frame('cell_count')
    sample_counts.to_csv(f"{export_dir}/05_cell_counts_by_sample.csv")

# 5. Cluster composition by sample
if cluster_col and sample_col:
    composition_abs = pd.crosstab(adata.obs[sample_col], adata.obs[cluster_col])
    composition_pct = pd.crosstab(adata.obs[sample_col], adata.obs[cluster_col], 
                                 normalize='index') * 100
    
    composition_abs.to_csv(f"{export_dir}/06_cluster_composition_absolute.csv")
    composition_pct.to_csv(f"{export_dir}/07_cluster_composition_percentage.csv")

print(f"\nAll summary tables exported to: {export_dir}/")
print("\nExported files:")
for file in sorted(os.listdir(export_dir)):
    print(f"  - {file}")

# %% [markdown]
# ## 9. Summary
# 
# All QC visualizations and summary tables have been generated and saved to the output directory.
# 
# **Output Structure:**
# ```
# qc_visualization_output/
# ├── 00_overall_qc_distribution.pdf/png
# ├── 01_qc_violin_by_cluster.pdf/png
# ├── 02_qc_boxplot_by_cluster.pdf/png
# ├── 03_qc_ridge_*_by_cluster.pdf
# ├── 04_qc_violin_by_sample.pdf/png
# ├── 05_qc_heatmap_by_sample.pdf/png
# ├── 06_cell_counts.pdf/png
# ├── 07_cluster_composition_by_sample.pdf/png
# ├── 08_doublet_analysis.pdf/png (if available)
# ├── 09_contamination_analysis.pdf/png (if available)
# ├── QC_SUMMARY_REPORT.txt
# └── tables/
#     ├── 01_overall_qc_statistics.csv
#     ├── 02_qc_statistics_by_cluster.csv
#     ├── 03_qc_statistics_by_sample.csv
#     ├── 04_cell_counts_by_cluster.csv
#     ├── 05_cell_counts_by_sample.csv
#     ├── 06_cluster_composition_absolute.csv
#     └── 07_cluster_composition_percentage.csv
# ```

# %%
print("\n" + "="*80)
print("QC VISUALIZATION PIPELINE COMPLETED!")
print("="*80)
print(f"\nAll outputs saved to: {output_dir}/")
print(f"\nGenerated {len([f for f in os.listdir(output_dir) if f.endswith(('.pdf', '.png'))])} plots")
print(f"Generated {len(os.listdir(export_dir))} summary tables")
print("\nPipeline finished successfully!")


