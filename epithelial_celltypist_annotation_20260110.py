# %% [markdown]
# # CellTypist Automated Cell Type Annotation for Epithelial Cells
# 
# **Purpose:**
# - Load epithelial cell dataset with manual annotations
# - Run CellTypist automated annotation with customizable model
# - Compare predictions with manual annotations
# - Generate comprehensive visualizations
# 
# **Author:** r2end  
# **Date:** 2025-01-10  
# **Version:** 1.0
# 
# **Key Features:**
# - Memory-optimized (use adata.raw for full gene access)
# - Gene matching with base symbols (no -1/-2 suffixes)
# - Majority voting for robust predictions
# - Index-aligned result writing
# - Comprehensive quality control

# %% [markdown]
# ## 📋 Configuration Parameters

# %%
import os
import sys
import gc
import warnings
from pathlib import Path
from datetime import datetime
import time

import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ============================================================================
# INPUT/OUTPUT CONFIGURATION
# ============================================================================

# --- Input Data ---
INPUT_H5AD = "/home/h2048/data/core2/adata_epithelial_with_manual_annotations_raw.h5ad"

# --- Output Directory ---
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/celltypist_epithelial")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# CELLTYPIST MODEL CONFIGURATION (⬅️ MODIFY HERE)
# ============================================================================

# --- Available Models ---
# To see all available models, run: models.models_description()
# Common options:
#   - 'Immune_All_Low.pkl'           # Pan-immune, lower resolution
#   - 'Immune_All_High.pkl'          # Pan-immune, higher resolution
#   - 'Developing_human_lung.pkl'    # Lung-specific
#   - 'Healthy_adult_lung.pkl'       # Adult lung (if available)
#   - 'Human_Lung_Atlas.pkl'         # Human lung atlas (if available)

CELLTYPIST_MODEL = "/home/h2048/data/source/reference/celltypist_models/Cells_Lung_Airway.pkl"  # ⬅️ CHANGE THIS

# --- CellTypist Parameters ---
CELLTYPIST_MAJORITY_VOTING = True        # Use majority voting (recommended: True)
CELLTYPIST_OVER_CLUSTERING = None        # Clustering resolution (None = auto)

# ============================================================================
# ANNOTATION COMPARISON
# ============================================================================

# --- Manual Annotation Key ---
MANUAL_ANNOTATION_KEY = "Manual_Annotation"      # Column name with manual annotations

# ============================================================================
# VISUALIZATION SETTINGS
# ============================================================================

FIGURE_DPI = 300
FIGURE_FORMAT = "pdf"                    # Options: "pdf", "png"

# ============================================================================
# PERFORMANCE SETTINGS
# ============================================================================

N_JOBS = 48                              # Parallel processing threads

# ============================================================================
# APPLY SETTINGS
# ============================================================================

sc.settings.n_jobs = N_JOBS
sc.settings.figdir = OUTPUT_DIR / "figures"
sc.settings.figdir.mkdir(exist_ok=True)

print("="*80)
print("CellTypist Automated Annotation Pipeline")
print("="*80)
print(f"Input:  {INPUT_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print(f"Model:  {CELLTYPIST_MODEL}")
print(f"Date:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# %% [markdown]
# ## 📂 Step 1: Load Data

# %%
print("\n" + "="*80)
print("STEP 1: Loading Data")
print("="*80)

start_time = time.time()

print(f"\nLoading: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

elapsed = time.time() - start_time
print(f"✓ Loaded in {elapsed:.1f}s")

print(f"\n📊 Dataset Summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Batches: {adata.obs['batch'].nunique() if 'batch' in adata.obs else 'N/A'}")

# %%
# Check data structure
print(f"\n🔍 Data Structure Check:")
print(f"  .X shape: {adata.X.shape}")
print(f"  .X type: {type(adata.X)}")
print(f"  .layers: {list(adata.layers.keys())}")
print(f"  .raw: {'✓ Present' if adata.raw is not None else '✗ Missing'}")

if adata.raw is not None:
    print(f"  .raw.X shape: {adata.raw.X.shape}")
    print(f"  .raw genes: {adata.raw.n_vars:,}")

# %%
# Check manual annotations
print(f"\n🏷️  Manual Annotations Check:")

if MANUAL_ANNOTATION_KEY in adata.obs.columns:
    manual_types = adata.obs[MANUAL_ANNOTATION_KEY].value_counts()
    print(f"✓ Found manual annotations: '{MANUAL_ANNOTATION_KEY}'")
    print(f"  Cell types: {len(manual_types)}")
    print(f"\nTop cell types:")
    print(manual_types.head(10))
else:
    print(f"⚠️  Manual annotation key '{MANUAL_ANNOTATION_KEY}' not found")
    print(f"   Available columns: {list(adata.obs.columns[:10])}...")

# %% [markdown]
# ## 🧬 Step 2: Prepare Gene Names

# %%
print("\n" + "="*80)
print("STEP 2: Gene Name Preparation")
print("="*80)

# Check if symbol_base already exists
if 'symbol_base' not in adata.var.columns:
    print("\n⚠️  'symbol_base' not found, creating from var_names...")
    
    # Remove -1, -2, -3 suffixes to get base symbols
    adata.var['symbol_base'] = adata.var_names.str.replace(r'-\d+$', '', regex=True)
    print("✓ Created symbol_base column")
else:
    print("✓ Found existing symbol_base column")

# %%
# Patch symbol_base to .raw if needed
if adata.raw is not None and 'symbol_base' not in adata.raw.var.columns:
    print("\n🔧 Patching symbol_base to adata.raw...")
    
    # Match by var_names
    common = adata.var_names.intersection(adata.raw.var_names)
    
    if len(common) > 0:
        adata.raw.var.loc[common, 'symbol_base'] = adata.var.loc[common, 'symbol_base'].values
        print(f"✓ Patched {len(common):,} genes")
    else:
        print("⚠️  No common genes found, creating from raw var_names...")
        adata.raw.var['symbol_base'] = adata.raw.var_names.str.replace(r'-\d+$', '', regex=True)

# %%
# Gene name statistics
print(f"\n📊 Gene Name Statistics:")
print(f"  Total genes (main): {adata.n_vars:,}")

if 'symbol_base' in adata.var.columns:
    n_unique = adata.var['symbol_base'].nunique()
    n_dup = adata.var['symbol_base'].duplicated().sum()
    print(f"  Unique base symbols: {n_unique:,}")
    print(f"  Duplicated symbols: {n_dup:,}")
    
    if n_dup > 0:
        print(f"\n  Top duplicated symbols:")
        dup_symbols = adata.var.loc[adata.var['symbol_base'].duplicated(keep=False), 'symbol_base']
        print(dup_symbols.value_counts().head())

if adata.raw is not None and 'symbol_base' in adata.raw.var.columns:
    n_unique_raw = adata.raw.var['symbol_base'].nunique()
    print(f"  Total genes (raw): {adata.raw.n_vars:,}")
    print(f"  Unique symbols (raw): {n_unique_raw:,}")

# %% [markdown]
# ## 📦 Step 3: Load CellTypist Model

# %%
print("\n" + "="*80)
print("STEP 3: Loading CellTypist Model")
print("="*80)

print(f"\n📥 Loading model: {CELLTYPIST_MODEL}")

try:
    # Try to load from local cache first
    model = models.Model.load(model=CELLTYPIST_MODEL)
    print(f"✓ Model loaded from cache")
except:
    # Download if not available
    print(f"⬇️  Model not in cache, downloading...")
    models.download_models(model=CELLTYPIST_MODEL, force_update=False)
    model = models.Model.load(model=CELLTYPIST_MODEL)
    print(f"✓ Model downloaded and loaded")

# %%
# Model information
print(f"\n📋 Model Information:")
print(f"  Model name: {model.model_name}")
print(f"  Cell types: {len(model.cell_types)}")
print(f"  Features: {len(model.features)}")

print(f"\n  Model cell types:")
for ct in sorted(model.cell_types)[:20]:  # Show first 20
    print(f"    - {ct}")
if len(model.cell_types) > 20:
    print(f"    ... and {len(model.cell_types) - 20} more")

# %% [markdown]
# ## 🔍 Step 4: Gene Matching & Feature Preparation

# %%
print("\n" + "="*80)
print("STEP 4: Gene Matching & Feature Preparation")
print("="*80)

# Get model features
model_features = set(model.features)
print(f"\n📊 Feature Overlap Analysis:")
print(f"  Model features: {len(model_features):,}")

# %%
# ⭐ CRITICAL: Use symbol_base from .raw for matching (no -1/-2 suffixes)
print(f"\n🔍 Building CellTypist input with base symbol matching...")

# Use .raw if available (full gene set), otherwise use main adata
source_adata = adata.raw if adata.raw is not None else adata

# Get base symbols without duplicates
sym = pd.Index(source_adata.var['symbol_base']).astype(str)

# Find overlap: genes in model AND not duplicated in our data
keep_mask = sym.isin(model_features) & ~sym.duplicated(keep='first')
n_overlap = keep_mask.sum()

print(f"  Our data: {len(sym):,} genes")
print(f"  Overlap: {n_overlap:,} genes")
print(f"  Coverage: {n_overlap/len(model_features)*100:.1f}%")

if n_overlap < 1000:
    print(f"\n  ⚠️  WARNING: Low overlap ({n_overlap} genes)")
    print(f"  Model may not perform well. Consider:")
    print(f"    - Using a different model")
    print(f"    - Checking gene name format")

# %%
# Build CellTypist input (only overlapping genes, deduplicated)
print(f"\n🔧 Building CellTypist AnnData...")

X_ct = source_adata.X[:, keep_mask].copy()
var_ct = pd.DataFrame(index=sym[keep_mask])

adata_celltypist = sc.AnnData(
    X=X_ct,
    obs=pd.DataFrame(index=adata.obs_names),
    var=var_ct
)

print(f"✓ CellTypist input: {adata_celltypist.n_obs:,} cells × {adata_celltypist.n_vars:,} genes")

# Memory cleanup
del X_ct, var_ct
gc.collect()

# %% [markdown]
# ## 🔬 Step 5: Normalize Data for CellTypist

# %%
print("\n" + "="*80)
print("STEP 5: Data Normalization")
print("="*80)

print(f"\n📊 Normalizing to 10,000 counts per cell...")
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
print(f"✓ Normalized")

print(f"\n📊 Log1p transformation...")
sc.pp.log1p(adata_celltypist)
print(f"✓ Log-transformed")

print(f"\n✓ Data ready for CellTypist prediction")

# %% [markdown]
# ## 🚀 Step 6: Run CellTypist Prediction

# %%
print("\n" + "="*80)
print("STEP 6: CellTypist Prediction")
print("="*80)

print(f"\n🚀 Running CellTypist...")
print(f"  Model: {CELLTYPIST_MODEL}")
print(f"  Majority voting: {CELLTYPIST_MAJORITY_VOTING}")
print(f"  Over-clustering: {CELLTYPIST_OVER_CLUSTERING}")

start_time = time.time()

predictions = celltypist.annotate(
    adata_celltypist,
    model=model,
    majority_voting=CELLTYPIST_MAJORITY_VOTING,
    over_clustering=CELLTYPIST_OVER_CLUSTERING
)

elapsed = time.time() - start_time
print(f"\n✓ Prediction complete: {int(elapsed//60)}m {int(elapsed%60)}s")

# Memory cleanup
del adata_celltypist
gc.collect()

# %% [markdown]
# ## 📝 Step 7: Extract and Align Results

# %%
print("\n" + "="*80)
print("STEP 7: Extracting Results")
print("="*80)

# ⭐ CRITICAL: Use majority voting results if enabled, otherwise use direct predictions
print(f"\n🔍 Extracting predictions...")

pred_df = predictions.predicted_labels

# Ensure index alignment (CRITICAL for correct cell mapping)
pred_df = pred_df.reindex(adata.obs_names)

print(f"✓ Predictions extracted and aligned")
print(f"  Total cells: {len(pred_df):,}")
print(f"  Non-null: {pred_df.notna().sum():,}")

# %%
# Add predictions to main adata
if CELLTYPIST_MAJORITY_VOTING:
    # Majority voting results (more robust)
    adata.obs['celltypist_pred'] = pred_df['majority_voting'].values
    adata.obs['celltypist_conf'] = pred_df['conf_score'].values
    
    # Also keep original predictions
    adata.obs['celltypist_pred_original'] = pred_df['predicted_labels'].values
    
    print(f"\n✓ Added columns:")
    print(f"  - celltypist_pred (majority voting)")
    print(f"  - celltypist_conf (confidence score)")
    print(f"  - celltypist_pred_original (direct prediction)")
else:
    # Direct predictions only
    adata.obs['celltypist_pred'] = pred_df['predicted_labels'].values
    adata.obs['celltypist_conf'] = pred_df['conf_score'].values
    
    print(f"\n✓ Added columns:")
    print(f"  - celltypist_pred (direct prediction)")
    print(f"  - celltypist_conf (confidence score)")

# %%
# Prediction summary
print(f"\n📊 Prediction Summary:")

pred_counts = adata.obs['celltypist_pred'].value_counts()
print(f"  Predicted cell types: {len(pred_counts)}")
print(f"\nTop predictions:")
print(pred_counts.head(15))

# %%
# Confidence score statistics
print(f"\n📊 Confidence Score Statistics:")

conf_stats = adata.obs['celltypist_conf'].describe()
print(conf_stats)

# Flag low-confidence predictions
LOW_CONF_THRESHOLD = 0.5
low_conf = (adata.obs['celltypist_conf'] < LOW_CONF_THRESHOLD).sum()
print(f"\n⚠️  Low confidence predictions (<{LOW_CONF_THRESHOLD}): {low_conf:,} ({low_conf/adata.n_obs*100:.1f}%)")

# %% [markdown]
# ## 📊 Step 8: Compare with Manual Annotations

# %%
print("\n" + "="*80)
print("STEP 8: Comparing with Manual Annotations")
print("="*80)

if MANUAL_ANNOTATION_KEY in adata.obs.columns:
    print(f"\n🔍 Comparing predictions with manual annotations...")
    
    # Create comparison dataframe
    comparison = pd.DataFrame({
        'manual': adata.obs[MANUAL_ANNOTATION_KEY],
        'predicted': adata.obs['celltypist_pred'],
        'confidence': adata.obs['celltypist_conf']
    })
    
    # Remove NaN values
    comparison = comparison.dropna(subset=['manual', 'predicted'])
    
    print(f"✓ Comparison ready")
    print(f"  Valid pairs: {len(comparison):,}")
    
    # %%
    # Confusion matrix
    from sklearn.metrics import confusion_matrix, classification_report
    
    print(f"\n📊 Classification Metrics:")
    
    # Overall agreement
    agreement = (comparison['manual'] == comparison['predicted']).sum()
    total = len(comparison)
    accuracy = agreement / total * 100
    
    print(f"  Overall agreement: {agreement:,} / {total:,} ({accuracy:.1f}%)")
    
    # %%
    # Detailed classification report
    print(f"\n📋 Classification Report:")
    
    try:
        report = classification_report(
            comparison['manual'],
            comparison['predicted'],
            output_dict=False,
            zero_division=0
        )
        print(report)
    except Exception as e:
        print(f"  ⚠️  Could not generate full report: {e}")
        print(f"  (Too many classes or class mismatch)")
    
    # %%
    # Crosstab: Manual vs Predicted
    print(f"\n📊 Crosstab (Manual vs Predicted):")
    
    crosstab = pd.crosstab(
        comparison['manual'],
        comparison['predicted'],
        margins=True,
        margins_name='Total'
    )
    
    print(crosstab)
    
    # Save crosstab
    crosstab_file = OUTPUT_DIR / "comparison_crosstab.csv"
    crosstab.to_csv(crosstab_file)
    print(f"\n✓ Crosstab saved: {crosstab_file}")

else:
    print(f"\n⚠️  Skipping comparison (manual annotations not found)")

# %% [markdown]
# ## 📈 Step 9: Visualization

# %%
print("\n" + "="*80)
print("STEP 9: Visualization")
print("="*80)

# Check if UMAP coordinates exist
if 'X_umap' in adata.obsm:
    print(f"\n✓ UMAP coordinates found")
    has_umap = True
elif 'X_umap_scvi' in adata.obsm:
    print(f"\n✓ scVI UMAP coordinates found, using as default")
    adata.obsm['X_umap'] = adata.obsm['X_umap_scvi']
    has_umap = True
else:
    print(f"\n⚠️  UMAP coordinates not found")
    print(f"   Available embeddings: {list(adata.obsm.keys())}")
    has_umap = False

# %%
# Visualization 1: Predicted cell types on UMAP
if has_umap:
    print(f"\n📊 Plotting predicted cell types...")
    
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(
        adata,
        color='celltypist_pred',
        title=f'CellTypist Predictions ({CELLTYPIST_MODEL})',
        legend_loc='right margin',
        frameon=False,
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(
        sc.settings.figdir / f"umap_celltypist_predictions.{FIGURE_FORMAT}",
        dpi=FIGURE_DPI,
        bbox_inches='tight'
    )
    plt.show()
    print(f"✓ Saved: umap_celltypist_predictions.{FIGURE_FORMAT}")

# %%
# Visualization 2: Confidence scores on UMAP
if has_umap:
    print(f"\n📊 Plotting confidence scores...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color='celltypist_conf',
        title='CellTypist Confidence Scores',
        cmap='RdYlGn',
        vmin=0,
        vmax=1,
        frameon=False,
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(
        sc.settings.figdir / f"umap_confidence_scores.{FIGURE_FORMAT}",
        dpi=FIGURE_DPI,
        bbox_inches='tight'
    )
    plt.show()
    print(f"✓ Saved: umap_confidence_scores.{FIGURE_FORMAT}")

# %%
# Visualization 3: Manual vs Predicted (if available)
if has_umap and MANUAL_ANNOTATION_KEY in adata.obs.columns:
    print(f"\n📊 Plotting manual vs predicted comparison...")
    
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    # Manual annotations
    sc.pl.umap(
        adata,
        color=MANUAL_ANNOTATION_KEY,
        title='Manual Annotations',
        legend_loc='right margin',
        frameon=False,
        ax=axes[0],
        show=False
    )
    
    # Predicted annotations
    sc.pl.umap(
        adata,
        color='celltypist_pred',
        title='CellTypist Predictions',
        legend_loc='right margin',
        frameon=False,
        ax=axes[1],
        show=False
    )
    
    plt.tight_layout()
    plt.savefig(
        sc.settings.figdir / f"umap_manual_vs_predicted.{FIGURE_FORMAT}",
        dpi=FIGURE_DPI,
        bbox_inches='tight'
    )
    plt.show()
    print(f"✓ Saved: umap_manual_vs_predicted.{FIGURE_FORMAT}")

# %%
# Visualization 4: Cell type composition barplot
print(f"\n📊 Plotting cell type composition...")

fig, ax = plt.subplots(figsize=(12, 6))

pred_counts = adata.obs['celltypist_pred'].value_counts()
pred_counts.plot(kind='bar', ax=ax, color='steelblue')

ax.set_xlabel('Predicted Cell Type', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_title(f'Cell Type Composition (CellTypist: {CELLTYPIST_MODEL})', fontsize=14)
ax.tick_params(axis='x', rotation=45, labelsize=10)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(
    sc.settings.figdir / f"barplot_cell_composition.{FIGURE_FORMAT}",
    dpi=FIGURE_DPI,
    bbox_inches='tight'
)
plt.show()
print(f"✓ Saved: barplot_cell_composition.{FIGURE_FORMAT}")

# %%
# Visualization 5: Confidence distribution by cell type
print(f"\n📊 Plotting confidence distribution by cell type...")

fig, ax = plt.subplots(figsize=(14, 6))

# Prepare data
conf_by_type = adata.obs[['celltypist_pred', 'celltypist_conf']].copy()
conf_by_type = conf_by_type.sort_values('celltypist_pred')

# Violin plot
import seaborn as sns
sns.violinplot(
    data=conf_by_type,
    x='celltypist_pred',
    y='celltypist_conf',
    ax=ax,
    inner='box',
    palette='Set2'
)

ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Low confidence threshold')
ax.set_xlabel('Predicted Cell Type', fontsize=12)
ax.set_ylabel('Confidence Score', fontsize=12)
ax.set_title('Confidence Score Distribution by Cell Type', fontsize=14)
ax.tick_params(axis='x', rotation=45, labelsize=9)
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(
    sc.settings.figdir / f"violin_confidence_by_celltype.{FIGURE_FORMAT}",
    dpi=FIGURE_DPI,
    bbox_inches='tight'
)
plt.show()
print(f"✓ Saved: violin_confidence_by_celltype.{FIGURE_FORMAT}")

# %% [markdown]
# ## 💾 Step 10: Save Results

# %%
print("\n" + "="*80)
print("STEP 10: Saving Results")
print("="*80)

# Save annotated adata
output_h5ad = OUTPUT_DIR / "adata_epithelial_celltypist_annotated.h5ad"

print(f"\n💾 Saving annotated data...")
adata.write_h5ad(output_h5ad, compression='gzip', compression_opts=9)

file_size = output_h5ad.stat().st_size / 1e9
print(f"✓ Saved: {output_h5ad}")
print(f"  Size: {file_size:.2f} GB")

# %%
# Export predictions to CSV
print(f"\n📄 Exporting predictions to CSV...")

predictions_df = adata.obs[[
    'celltypist_pred',
    'celltypist_conf'
]].copy()

if CELLTYPIST_MAJORITY_VOTING:
    predictions_df['celltypist_pred_original'] = adata.obs['celltypist_pred_original']

if MANUAL_ANNOTATION_KEY in adata.obs.columns:
    predictions_df[MANUAL_ANNOTATION_KEY] = adata.obs[MANUAL_ANNOTATION_KEY]

if 'batch' in adata.obs.columns:
    predictions_df['batch'] = adata.obs['batch']

predictions_csv = OUTPUT_DIR / "celltypist_predictions.csv"
predictions_df.to_csv(predictions_csv)

print(f"✓ Saved: {predictions_csv}")

# %%
# Save summary statistics
print(f"\n📊 Saving summary statistics...")

summary_file = OUTPUT_DIR / "annotation_summary.txt"

with open(summary_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("CellTypist Annotation Summary\n")
    f.write("="*80 + "\n\n")
    
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input: {INPUT_H5AD}\n")
    f.write(f"Model: {CELLTYPIST_MODEL}\n")
    f.write(f"Majority Voting: {CELLTYPIST_MAJORITY_VOTING}\n\n")
    
    f.write(f"Dataset:\n")
    f.write(f"  Total cells: {adata.n_obs:,}\n")
    f.write(f"  Total genes: {adata.n_vars:,}\n\n")
    
    f.write(f"Predictions:\n")
    f.write(f"  Cell types predicted: {adata.obs['celltypist_pred'].nunique()}\n")
    f.write(f"  Mean confidence: {adata.obs['celltypist_conf'].mean():.3f}\n")
    f.write(f"  Median confidence: {adata.obs['celltypist_conf'].median():.3f}\n\n")
    
    f.write(f"Cell Type Composition:\n")
    for ct, count in adata.obs['celltypist_pred'].value_counts().items():
        pct = count / adata.n_obs * 100
        f.write(f"  {ct}: {count:,} ({pct:.1f}%)\n")

print(f"✓ Saved: {summary_file}")

# %% [markdown]
# ## ✅ Pipeline Complete

# %%
print("\n" + "="*80)
print("PIPELINE COMPLETE")
print("="*80)

print(f"\n📂 Output Directory: {OUTPUT_DIR}")
print(f"\n📄 Generated Files:")
print(f"  1. adata_epithelial_celltypist_annotated.h5ad  (Annotated data)")
print(f"  2. celltypist_predictions.csv                  (Prediction table)")
print(f"  3. annotation_summary.txt                      (Summary statistics)")
if MANUAL_ANNOTATION_KEY in adata.obs.columns:
    print(f"  4. comparison_crosstab.csv                     (Manual vs Predicted)")

print(f"\n📊 Figures (in figures/ subdirectory):")
if has_umap:
    print(f"  - umap_celltypist_predictions.{FIGURE_FORMAT}")
    print(f"  - umap_confidence_scores.{FIGURE_FORMAT}")
    if MANUAL_ANNOTATION_KEY in adata.obs.columns:
        print(f"  - umap_manual_vs_predicted.{FIGURE_FORMAT}")
print(f"  - barplot_cell_composition.{FIGURE_FORMAT}")
print(f"  - violin_confidence_by_celltype.{FIGURE_FORMAT}")

print(f"\n✅ All tasks completed successfully!")

# %%
# Final data check
print(f"\n🔍 Final Data Check:")
print(f"  adata.obs columns added:")
print(f"    - celltypist_pred")
print(f"    - celltypist_conf")
if CELLTYPIST_MAJORITY_VOTING:
    print(f"    - celltypist_pred_original")

print(f"\n  Example predictions:")
print(adata.obs[['celltypist_pred', 'celltypist_conf']].head(10))

# %%
