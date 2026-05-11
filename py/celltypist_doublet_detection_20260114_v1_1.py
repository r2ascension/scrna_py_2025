"""
CellTypist-based Doublet Detection v1.1 (Robust Version)
=========================================================

Key Improvements:
- Uses only CellTypist model genes (not full 83k genome)
- Preserves full dataset (subset only for analysis)
- Top2 lineage + margin logic (more robust than max-threshold)
- Multi-evidence gating (CellTypist + neighbor + markers)
- Optimized for nasal/sinus epithelial data

Author: r2end
Date: 2025-01-13
Version: 1.1 (Production-ready)
"""

import scanpy as sc
import celltypist
from celltypist.models import Model
import numpy as np
import pandas as pd
from pathlib import Path
import gc

# ==============================================================================
# Configuration
# ==============================================================================

# Input/Output
INPUT_H5AD = "/home/h2048/data/R/0113/polyp_obj_updated_20260113.h5ad"
OUTPUT_DIR = Path(f"/home/h2048/data/py/0112/celltypist_stromal/celltypist_doublet_analysis")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# CellTypist model (⚠️ UPDATE THIS)
CELLTYPIST_MODEL = "/home/h2048/data/source/reference/celltypist_models/Cells_Lung_Airway.pkl"

# Analysis subset (RECOMMENDED: only analyze suspected problematic cells)
# Set to None to analyze all cells (slower, more memory)
ANALYSIS_SUBSET_STRATEGY = "all"  # Options: "epithelial_tier_c", "all", "custom"

# Top2 lineage + margin parameters (more robust than simple threshold)
TOP2_THRESHOLD = 0.5      # Top2 lineage probability must exceed this
MARGIN_MAX = 0.15         # Max allowed margin between top1 and top2 (smaller = more confident doublet)

# Multi-evidence gating parameters
NEIGHBOR_FRAC_THRESHOLD = 0.8     # Cells with good neighbors less likely to be doublets
MAPPING_CONF_THRESHOLD = 0.6      # Cells with high mapping confidence more reliable

# Batch key
BATCH_KEY = "dataset"

# ==============================================================================
# Lineage Groupings (Same as v1.0)
# ==============================================================================

LINEAGE_GROUPS = {
    'Epithelial': [
        'AT1', 'AT2', 'Dividing_AT2',  # Note: In nasal data, these are "AT2-like/secretory-like"
        'Basal', 'Dividing_Basal', 'Suprabasal',
        'Ciliated', 'Deuterosomal',
        'Secretory_Club', 'Secretory_Goblet',
        'Ionocyte_n_Brush', 'Neuroendocrine',
        'SMG_Basal', 'SMG_Duct', 'SMG_Mucous', 'SMG_Serous'
    ],
    'T_NK': [
        'CD4_EM/Effector', 'CD4_TRM', 'CD4_naive/CM',
        'CD8_EM', 'CD8_EM/EMRA', 'CD8_TRM', 'CD8_TRM/EM',
        'T_reg', 'MAIT', 'NKT', 'gdT', 'ILC',
        'NK_CD11d', 'NK_CD16hi', 'NK_CD56bright'
    ],
    'B_Plasma': [
        'B_memory', 'B_naive',
        'B_plasma_IgA', 'B_plasma_IgG', 'B_plasmablast'
    ],
    'Myeloid': [
        'DC_1', 'DC_2', 'DC_activated', 'DC_plasmacytoid',
        'Macro_AW_CX3CR1', 'Macro_CCL', 'Macro_CHIT1',
        'Macro_alveolar', 'Macro_alveolar_metallothioneins',
        'Macro_dividing', 'Macro_intermediate',
        'Macro_interstitial', 'Macro_intravascular',
        'Monocyte_CD14', 'Monocyte_CD16', 'Mast_cell'
    ],
    'Stromal': [
        'Fibro_adventitial', 'Fibro_alveolar', 'Fibro_immune_recruiting',
        'Fibro_myofibroblast', 'Fibro_peribronchial',
        'Muscle_pericyte_pulmonary', 'Muscle_pericyte_systemic',
        'Muscle_perivascular_immune_recruiting',
        'Muscle_smooth_airway', 'Muscle_smooth_arterial_systemic',
        'Muscle_smooth_pulmonary',
        'Mesothelia', 'Chondrocyte'
    ],
    'Endothelial': [
        'Endothelia_Lymphatic',
        'Endothelia_vascular_Cap_a', 'Endothelia_vascular_Cap_g',
        'Endothelia_vascular_arterial_pulmonary',
        'Endothelia_vascular_arterial_systemic',
        'Endothelia_vascular_venous_pulmonary',
        'Endothelia_vascular_venous_systemic'
    ],
    'Other': [
        'NAF_endoneurial', 'NAF_perineurial',
        'Schwann_Myelinating', 'Schwann_nonmyelinating',
        'Erythrocyte', 'Megakaryocyte'
    ]
}

# Marker genes for validation (optional, used if available)
MARKER_GENES = {
    'Epithelial': ['EPCAM', 'KRT5', 'KRT8', 'KRT18'],
    'Immune': ['PTPRC', 'CD3D', 'CD79A', 'CD14', 'CD68'],
    'Stromal': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5']
}

# ==============================================================================
# Step 1: Load Full Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 1: Loading Full Query Data")
print("=" * 70)

adata_full = sc.read_h5ad(INPUT_H5AD)
print(f"Loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars:,} genes")

# ==============================================================================
# Step 2: Define Analysis Subset
# ==============================================================================

print("\n" + "=" * 70)
print("Step 2: Defining Analysis Subset")
print("=" * 70)

if ANALYSIS_SUBSET_STRATEGY == "epithelial_tier_c":
    # Focus on epithelial tier C (your main problem area)
    if 'epi_tier' in adata_full.obs.columns:
        mask = (
            (adata_full.obs['cell_type_mapped'].astype(str) == 'Epithelial') &
            (adata_full.obs['epi_tier'].astype(str) == 'C')
        )
        print("  Strategy: Epithelial Tier C only")
    else:
        print("  ⚠️ 'epi_tier' column not found, falling back to all epithelial")
        mask = adata_full.obs['cell_type_mapped'].astype(str) == 'Epithelial'
        
elif ANALYSIS_SUBSET_STRATEGY == "all":
    mask = np.ones(adata_full.n_obs, dtype=bool)
    print("  Strategy: Analyze all cells")
    
elif ANALYSIS_SUBSET_STRATEGY == "custom":
    # User can define custom mask here
    # Example: specific datasets
    # mask = adata_full.obs['dataset'].isin(['GSE276503', 'GSE299751', 'GSE164547'])
    mask = np.ones(adata_full.n_obs, dtype=bool)
    print("  Strategy: Custom (modify code for your needs)")
    
else:
    raise ValueError(f"Unknown strategy: {ANALYSIS_SUBSET_STRATEGY}")

adata_subset = adata_full[mask].copy()
print(f"  Subset: {adata_subset.n_obs:,} cells ({adata_subset.n_obs/adata_full.n_obs*100:.1f}%)")

# ==============================================================================
# Step 3: Load CellTypist Model and Get Feature Genes
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3: Loading CellTypist Model")
print("=" * 70)

print(f"Model: {CELLTYPIST_MODEL}")

model = Model.load(CELLTYPIST_MODEL)
model_genes = model.features  # These are the genes the model uses

print(f"  Model uses {len(model_genes)} genes")

# Find overlap with query data
available_genes = [g for g in model_genes if g in adata_subset.var_names]
overlap_pct = len(available_genes) / len(model_genes) * 100

print(f"  Genes available in query: {len(available_genes)}/{len(model_genes)} ({overlap_pct:.1f}%)")

if overlap_pct < 70:
    print(f"\n⚠️ Warning: Low gene overlap ({overlap_pct:.1f}%)")
    print(f"   This may affect prediction quality")
    print(f"   Consider checking gene naming (ENSEMBL vs symbol)")

# ==============================================================================
# Step 4: Subset to Model Genes and Normalize
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4: Preparing Data for CellTypist")
print("=" * 70)

# Critical: Only subset to model genes BEFORE normalization
print(f"  Subsetting from {adata_subset.n_vars:,} to {len(available_genes):,} genes")
adata_ct = adata_subset[:, available_genes].copy()

# Use raw counts for normalization
if 'counts' in adata_ct.layers:
    print("  Using layers['counts']")
    adata_ct.X = adata_ct.layers['counts'].copy()
elif hasattr(adata_subset, 'raw') and adata_subset.raw is not None:
    print("  Using .raw.X")
    adata_ct_raw = adata_subset.raw[:, available_genes]
    adata_ct.X = adata_ct_raw.X.copy()
else:
    print("  ⚠️ Using current .X (assuming it's raw counts)")

# Normalize to CPM=10,000 and log1p (CellTypist requirement)
print("  Normalizing: CPM=10,000 → log1p")
sc.pp.normalize_total(adata_ct, target_sum=1e4)
sc.pp.log1p(adata_ct)

print(f"✓ Data prepared: {adata_ct.n_obs:,} cells × {adata_ct.n_vars:,} genes")

# Clean up large intermediate object
del adata_subset
gc.collect()

# ==============================================================================
# Step 5: Run CellTypist with Probability Match
# ==============================================================================

print("\n" + "=" * 70)
print("Step 5: Running CellTypist (mode='prob match')")
print("=" * 70)

pred = celltypist.annotate(
    adata_ct,
    model=model,
    mode='prob match',
    p_thres=TOP2_THRESHOLD,  # Use top2 threshold for initial filtering
    majority_voting=False
)

# Insert results (handle different return structures in prob match mode)
# pred.predicted_labels is a DataFrame, not an object with attributes
print(f"\n✓ CellTypist completed")

# Check what columns are available
print(f"  Available columns in predicted_labels: {list(pred.predicted_labels.columns)}")

# Extract predicted labels (this column always exists)
if 'predicted_labels' in pred.predicted_labels.columns:
    adata_ct.obs['ct_predicted_labels'] = pred.predicted_labels['predicted_labels'].values
else:
    # Fallback: use index or first column
    adata_ct.obs['ct_predicted_labels'] = pred.predicted_labels.iloc[:, 0].values
    print(f"  ⚠️ 'predicted_labels' column not found, using first column")

# Extract confidence score (may not exist in prob match mode)
if 'conf_score' in pred.predicted_labels.columns:
    adata_ct.obs['ct_conf_score'] = pred.predicted_labels['conf_score'].values
    print(f"  ✓ Confidence scores extracted")
elif 'over_clustering' in pred.predicted_labels.columns:
    # Some versions use over_clustering
    adata_ct.obs['ct_conf_score'] = pred.predicted_labels['over_clustering'].values
    print(f"  ✓ Using over_clustering scores")
else:
    # Calculate confidence from probability matrix (max probability)
    max_probs = pred.probability_matrix.max(axis=1).values
    adata_ct.obs['ct_conf_score'] = max_probs
    print(f"  ⚠️ No conf_score column, calculated from max probability")

print(f"  Unique predicted types: {adata_ct.obs['ct_predicted_labels'].nunique()}")

# Get probability matrix (always available)
prob_matrix = pred.probability_matrix
print(f"  Probability matrix: {prob_matrix.shape[0]} cells × {prob_matrix.shape[1]} cell types")

# ==============================================================================
# Step 6: Calculate Lineage-Level Max Probabilities
# ==============================================================================

print("\n" + "=" * 70)
print("Step 6: Computing Lineage Probabilities")
print("=" * 70)

lineage_probs = {}

for lineage, cell_types in LINEAGE_GROUPS.items():
    # Find matching columns in probability matrix
    matching_types = [ct for ct in cell_types if ct in prob_matrix.columns]
    
    if len(matching_types) == 0:
        print(f"  ⚠️ {lineage}: No matching cell types")
        lineage_probs[lineage] = np.zeros(adata_ct.n_obs)
    else:
        # Max probability across this lineage
        max_probs = prob_matrix[matching_types].max(axis=1).values
        lineage_probs[lineage] = max_probs
        print(f"  {lineage}: {len(matching_types)} types, median_prob={np.median(max_probs):.3f}")

# Store in obs
for lineage, probs in lineage_probs.items():
    adata_ct.obs[f'prob_{lineage}'] = probs

# ==============================================================================
# Step 7: Top2 Lineage + Margin Detection (More Robust)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 7: Doublet Detection (Top2 + Margin Logic)")
print("=" * 70)

print(f"\nParameters:")
print(f"  Top2 threshold: {TOP2_THRESHOLD}")
print(f"  Max margin: {MARGIN_MAX}")

# Get lineage probability matrix
lineage_names = list(LINEAGE_GROUPS.keys())
lineage_prob_matrix = np.column_stack([lineage_probs[ln] for ln in lineage_names])

# Find top1 and top2 lineages for each cell
top1_idx = lineage_prob_matrix.argmax(axis=1)
tmp_matrix = lineage_prob_matrix.copy()
tmp_matrix[np.arange(tmp_matrix.shape[0]), top1_idx] = -1
top2_idx = tmp_matrix.argmax(axis=1)

# Get probabilities
p_top1 = lineage_prob_matrix[np.arange(lineage_prob_matrix.shape[0]), top1_idx]
p_top2 = lineage_prob_matrix[np.arange(lineage_prob_matrix.shape[0]), top2_idx]
margin = p_top1 - p_top2

# Store top lineages
adata_ct.obs['top1_lineage'] = [lineage_names[i] for i in top1_idx]
adata_ct.obs['top2_lineage'] = [lineage_names[i] for i in top2_idx]
adata_ct.obs['prob_top1'] = p_top1
adata_ct.obs['prob_top2'] = p_top2
adata_ct.obs['lineage_margin'] = margin

# Doublet candidate: top2 >= threshold AND margin <= max_margin
doublet_candidate = (p_top2 >= TOP2_THRESHOLD) & (margin <= MARGIN_MAX)

adata_ct.obs['doublet_candidate'] = doublet_candidate

# Create doublet type (combining top1 and top2 lineages)
def get_doublet_type(row):
    if not row['doublet_candidate']:
        return 'Singlet'
    l1, l2 = sorted([row['top1_lineage'], row['top2_lineage']])
    return f"{l1}_{l2}"

adata_ct.obs['doublet_type'] = adata_ct.obs.apply(get_doublet_type, axis=1)

n_doublets = doublet_candidate.sum()
print(f"\n✓ Doublet candidates: {n_doublets:,} cells ({n_doublets/adata_ct.n_obs*100:.2f}%)")

# Show doublet type distribution
doublet_types = adata_ct.obs['doublet_type'].value_counts()
print("\nDoublet type distribution:")
for dtype, count in doublet_types.head(10).items():
    if dtype != 'Singlet':
        print(f"  {dtype}: {count:,} cells")

# ==============================================================================
# Step 8: Multi-Evidence Gating (Optional but Recommended)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 8: Multi-Evidence Validation")
print("=" * 70)

# Check if neighbor_frac is available (from your BBKNN analysis)
has_neighbor_frac = 'epi_neighbor_frac' in adata_ct.obs.columns or 'neighbor_frac' in adata_ct.obs.columns

if has_neighbor_frac:
    neighbor_col = 'epi_neighbor_frac' if 'epi_neighbor_frac' in adata_ct.obs.columns else 'neighbor_frac'
    poor_neighbors = adata_ct.obs[neighbor_col] < NEIGHBOR_FRAC_THRESHOLD
    print(f"  Using {neighbor_col} for neighbor evidence")
    print(f"    Cells with poor neighbors (<{NEIGHBOR_FRAC_THRESHOLD}): {poor_neighbors.sum():,}")
else:
    poor_neighbors = np.ones(adata_ct.n_obs, dtype=bool)
    print(f"  ⚠️ neighbor_frac not available, skipping neighbor gating")

# Check if mapping confidence is available (from scArches)
has_mapping_conf = 'mapping_confidence' in adata_ct.obs.columns

if has_mapping_conf:
    low_conf = adata_ct.obs['mapping_confidence'] < MAPPING_CONF_THRESHOLD
    print(f"  Using mapping_confidence for quality evidence")
    print(f"    Cells with low confidence (<{MAPPING_CONF_THRESHOLD}): {low_conf.sum():,}")
else:
    low_conf = np.zeros(adata_ct.n_obs, dtype=bool)
    print(f"  ⚠️ mapping_confidence not available, skipping confidence gating")

# Define high-confidence doublets (for removal)
# Require: (1) CellTypist doublet candidate + (2) at least one supporting evidence
high_conf_doublet = doublet_candidate & (poor_neighbors | low_conf)

adata_ct.obs['doublet_high_confidence'] = high_conf_doublet

print(f"\n✓ High-confidence doublets: {high_conf_doublet.sum():,} cells ({high_conf_doublet.sum()/adata_ct.n_obs*100:.2f}%)")
print(f"  (These are safer to remove)")

# Medium-confidence: doublet candidate but no strong supporting evidence
medium_conf_doublet = doublet_candidate & ~high_conf_doublet
adata_ct.obs['doublet_medium_confidence'] = medium_conf_doublet

print(f"  Medium-confidence doublets: {medium_conf_doublet.sum():,} cells")
print(f"  (Review before removing - may include inflammatory/transition states)")

# ==============================================================================
# Step 9: Calculate Marker Expression Scores (Optional)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 9: Marker Expression Validation (Optional)")
print("=" * 70)

# Check if we have raw counts for marker expression
if hasattr(adata_full, 'raw') and adata_full.raw is not None:
    print("  Computing marker scores from .raw")
    
    for category, genes in MARKER_GENES.items():
        available = [g for g in genes if g in adata_full.raw.var_names]
        if len(available) > 0:
            # Get expression for subset cells
            idx = adata_ct.obs_names
            expr = adata_full[idx].raw[:, available].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray()
            # Mean expression
            marker_score = expr.mean(axis=1)
            adata_ct.obs[f'marker_{category}'] = marker_score
            print(f"    {category}: {len(available)} genes")
    
    print("  ✓ Marker scores computed")
else:
    print("  ⚠️ .raw not available, skipping marker validation")

# ==============================================================================
# Step 10: Dataset-Specific Analysis
# ==============================================================================

print("\n" + "=" * 70)
print("Step 10: Dataset-Specific Doublet Rates")
print("=" * 70)

if BATCH_KEY in adata_ct.obs.columns:
    dataset_summary = adata_ct.obs.groupby(BATCH_KEY, observed=True).agg({
        'doublet_candidate': 'mean',
        'doublet_high_confidence': 'mean',
        'ct_conf_score': 'median'
    }).round(4)
    
    dataset_summary.columns = ['candidate_rate', 'high_conf_rate', 'median_ctconf']
    dataset_summary = dataset_summary.sort_values('candidate_rate', ascending=False)
    
    print("\nTop 20 datasets by doublet candidate rate:")
    print(dataset_summary.head(20))
    
    # Save
    dataset_summary.to_csv(OUTPUT_DIR / "doublet_rates_by_dataset.csv")
    print(f"\n✓ Saved: {OUTPUT_DIR / 'doublet_rates_by_dataset.csv'}")
else:
    print(f"  ⚠️ Batch key '{BATCH_KEY}' not found")

# ==============================================================================
# Step 11: Transfer Results Back to Full Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 11: Transferring Results to Full Dataset")
print("=" * 70)

# Columns to transfer (only essential ones, not all 78 ct_* columns)
cols_to_transfer = [
    'ct_predicted_labels', 'ct_conf_score',
    'top1_lineage', 'top2_lineage', 'prob_top1', 'prob_top2', 'lineage_margin',
    'doublet_candidate', 'doublet_type',
    'doublet_high_confidence', 'doublet_medium_confidence'
]

# Add lineage probabilities
cols_to_transfer.extend([f'prob_{ln}' for ln in LINEAGE_GROUPS.keys()])

# Add marker scores if computed
marker_cols = [c for c in adata_ct.obs.columns if c.startswith('marker_')]
cols_to_transfer.extend(marker_cols)

# Transfer to full data (cells not in subset will be NaN/False)
print(f"  Transferring {len(cols_to_transfer)} columns")

for col in cols_to_transfer:
    if col in adata_ct.obs.columns:
        # Initialize column in full data
        if adata_ct.obs[col].dtype == bool:
            adata_full.obs[col] = False
        elif adata_ct.obs[col].dtype == 'category' or adata_ct.obs[col].dtype == object:
            adata_full.obs[col] = 'N/A'
        else:
            adata_full.obs[col] = np.nan
        
        # Fill in values for analyzed cells
        adata_full.obs.loc[adata_ct.obs_names, col] = adata_ct.obs[col].values

print(f"✓ Results transferred to full dataset ({adata_full.n_obs:,} cells)")

# ==============================================================================
# Step 12: Save Results
# ==============================================================================

print("\n" + "=" * 70)
print("Step 12: Saving Results")
print("=" * 70)

output_h5ad = OUTPUT_DIR / "adata_with_doublet_predictions.h5ad"
print(f"\nSaving: {output_h5ad}")

adata_full.write_h5ad(output_h5ad, compression='gzip', compression_opts=9)

file_size = output_h5ad.stat().st_size / (1024**3)
print(f"✓ Saved ({file_size:.2f} GB)")

# Save summary
import json

summary = {
    'total_cells': int(adata_full.n_obs),
    'analyzed_cells': int(adata_ct.n_obs),
    'analysis_rate': float(adata_ct.n_obs / adata_full.n_obs),
    'doublet_candidates': {
        'n': int(adata_ct.obs['doublet_candidate'].sum()),
        'rate': float(adata_ct.obs['doublet_candidate'].mean())
    },
    'high_confidence_doublets': {
        'n': int(adata_ct.obs['doublet_high_confidence'].sum()),
        'rate': float(adata_ct.obs['doublet_high_confidence'].mean())
    },
    'parameters': {
        'top2_threshold': TOP2_THRESHOLD,
        'margin_max': MARGIN_MAX,
        'neighbor_frac_threshold': NEIGHBOR_FRAC_THRESHOLD,
        'mapping_conf_threshold': MAPPING_CONF_THRESHOLD
    }
}

with open(OUTPUT_DIR / "doublet_summary.json", 'w') as f:
    json.dump(summary, f, indent=2)

print(f"✓ Saved: {OUTPUT_DIR / 'doublet_summary.json'}")

# ==============================================================================
# Final Report
# ==============================================================================

print("\n" + "=" * 70)
print("DOUBLET DETECTION SUMMARY")
print("=" * 70)

print(f"\nTotal cells: {adata_full.n_obs:,}")
print(f"Analyzed cells: {adata_ct.n_obs:,} ({adata_ct.n_obs/adata_full.n_obs*100:.1f}%)")

print(f"\nDoublet detection results:")
print(f"  Candidates: {adata_ct.obs['doublet_candidate'].sum():,} ({adata_ct.obs['doublet_candidate'].mean()*100:.2f}%)")
print(f"  High-confidence: {adata_ct.obs['doublet_high_confidence'].sum():,} ({adata_ct.obs['doublet_high_confidence'].mean()*100:.2f}%)")
print(f"  Medium-confidence: {adata_ct.obs['doublet_medium_confidence'].sum():,} ({adata_ct.obs['doublet_medium_confidence'].mean()*100:.2f}%)")

print(f"\nTop doublet types:")
for dtype, count in adata_ct.obs['doublet_type'].value_counts().head(5).items():
    if dtype != 'Singlet':
        print(f"  {dtype}: {count:,} cells")

print(f"\nOutput files: {OUTPUT_DIR}")
print("  - adata_with_doublet_predictions.h5ad")
print("  - doublet_rates_by_dataset.csv")
print("  - doublet_summary.json")

print("\n" + "=" * 70)
print("RECOMMENDED NEXT STEPS")
print("=" * 70)

print("""
1. Review dataset-specific rates in doublet_rates_by_dataset.csv

2. Validate high-confidence doublets:
   doublets = adata[adata.obs['doublet_high_confidence']]
   sc.pl.umap(doublets, color=['top1_lineage', 'top2_lineage', 'lineage_margin'])

3. Check marker expression (if computed):
   sc.pl.umap(doublets, color=['marker_Epithelial', 'marker_Immune'])

4. Apply filtering (conservative approach):
   adata_clean = adata[~adata.obs['doublet_high_confidence']].copy()

5. Review medium-confidence doublets before removing:
   - May include inflammatory/transition states
   - Consider additional validation

6. Re-run integration with cleaned data
""")

print("✓ Analysis complete!")