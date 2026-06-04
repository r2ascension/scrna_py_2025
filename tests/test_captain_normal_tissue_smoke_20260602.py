#!/usr/bin/env python3
"""Smoke test CAPTAIN gene vocabulary coverage with normal tissue airway genes.

Uses bare dict lookups (avoids slow GeneVocab initialization with 60K entries).
For the full GeneVocab load test, see earlier terminal output in WORKLOG.
"""
import json
import pickle
import time

# 1. Load CAPTAIN vocab as dict (fast!)
t0 = time.time()
with open('/home/h2048/CAPTAIN/token_dict/vocab.json') as f:
    token2idx = json.load(f)
t1 = time.time()
print(f'CAPTAIN vocab size: {len(token2idx)} (loaded in {t1-t0:.3f}s)')

# Special tokens
for tok in ['<pad>', '<cls>', '<eoc>']:
    print(f'  Special: {tok}={token2idx.get(tok, "N/A")}')

# 2. Test airway/epithelial marker genes
test_genes = [
    'TP53', 'EGFR', 'KRT5', 'MUC5AC', 'FOXJ1', 'SFTPC', 'SCGB1A1',
    'CDH1', 'VIM', 'ACTA2', 'EPCAM', 'KRT8', 'KRT18', 'SOX2', 'SOX9',
    'DNAH5', 'FOXI1', 'ASCL1', 'CHGA', 'PDPN', 'COL1A1', 'PECAM1',
    'PTPRC', 'CD3E', 'CD4', 'CD8A', 'CD19', 'MS4A1', 'NCAM1', 'MKI67',
]
print('\n--- Airway/Epithelial marker genes ---')
found = sum(1 for g in test_genes if g in token2idx)
missing = [g for g in test_genes if g not in token2idx]
if missing:
    print(f'  MISSING: {missing}')
print(f'  Found: {found}/{len(test_genes)}')

# 3. Cell-type specific marker genes
tissue_markers = {
    'Basal': ['KRT5', 'KRT14', 'TP63'],
    'Club': ['SCGB1A1', 'SCGB3A1'],
    'Goblet': ['MUC5AC', 'MUC5B', 'SPDEF'],
    'Ciliated': ['FOXJ1', 'DNAH5', 'TPPP3'],
    'Ionocyte': ['FOXI1', 'CFTR', 'ASCL3'],
    'AT1': ['AGER', 'PDPN', 'HOPX'],
    'AT2': ['SFTPC', 'SFTPB', 'LAMP3'],
    'Neuroendocrine': ['ASCL1', 'CHGA', 'SYP'],
    'Basal_Inflammatory': ['KRT5', 'CXCL8', 'IL1B'],
}
print('\n--- Cell-type marker coverage ---')
for ct, genes in tissue_markers.items():
    ok = sum(1 for g in genes if g in token2idx)
    bad = [g for g in genes if g not in token2idx]
    status = 'OK' if not bad else f'MISSING: {bad}'
    print(f'  {ct}: {ok}/{len(genes)} {status}')

# 4. Protein token dict (CAPTAIN's unique CSP vocabulary)
with open('/home/h2048/CAPTAIN/token_dict/csp_token_dict.pickle', 'rb') as f:
    csp_tokens = pickle.load(f)
print(f'\nCAPTAIN CSP proteins: {len(csp_tokens)}')
print(f'  Sample: {sorted(csp_tokens.keys())[:10]}')

# 5. Check normal tissue gene overlap (fast set intersection)
gene_df_col0 = set()
with open('/home/h2048/data/coredata0524/all_genes_union.csv') as f:
    next(f)  # skip header
    for line in f:
        gene_df_col0.add(line.strip().split(',')[0].upper().strip())
captain_set = set(k.upper().strip() for k in token2idx)
overlap = gene_df_col0 & captain_set
print(f'\nGene overlap: {len(overlap)}/{len(gene_df_col0)} normal genes in CAPTAIN ({100*len(overlap)/max(len(gene_df_col0),1):.1f}%)')

# 6. Environment
import torch
print(f'\ntorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'GPU: {torch.cuda.get_device_name(0)}')

print(f'\n=== SMOKE TEST PASSED ===')
