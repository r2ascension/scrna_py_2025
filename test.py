import anndata as ad
import h5py

# 1. 检查文件完整性
h5ad_path = '/home/h2048/data/py/1029/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad'

# 方法A: 尝试直接读取
try:
    adata = ad.read_h5ad(h5ad_path)
    print("✓ Python anndata 可以正常读取")
    print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    print(f"Layers: {list(adata.layers.keys())}")
except Exception as e:
    print(f"✗ Python 读取失败: {e}")

# 方法B: 检查HDF5文件结构0
with h5py.File(h5ad_path, 'r') as f:
    print("\n=== HDF5 Structure ===")
    print(f"Keys: {list(f.keys())}")
    if 'layers' in f:
        print(f"Layers: {list(f['layers'].keys())}")
        if 'counts' in f['layers']:
            counts_dataset = f['layers']['counts']
            print(f"Counts compression: {counts_dataset.compression}")
            print(f"Counts shape: {counts_dataset.shape}")