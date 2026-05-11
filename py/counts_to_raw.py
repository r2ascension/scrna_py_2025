#!/usr/bin/env python3
import sys
import anndata as ad

in_path  = '/home/h2048/data/py/1204/Tcell_scANVI/adata_tcell_scANVI_annotated.h5ad'
out_path = '/home/h2048/data/py/1204/Tcell_scANVI/adata_tcell_scANVI_annotated_raw.h5ad'

print(f"Reading: {in_path}")
adata = ad.read_h5ad(in_path)

# 检查 layers 里有没有 counts/raw_counts
layer_keys = list(adata.layers.keys())
print("Available layers:", layer_keys)

raw_mat = None
if "counts" in layer_keys:
    print("Using adata.layers['counts'] as raw counts")
    raw_mat = adata.layers["counts"]
elif "raw_counts" in layer_keys:
    print("Using adata.layers['raw_counts'] as raw counts")
    raw_mat = adata.layers["raw_counts"]
else:
    raise ValueError("No 'counts' or 'raw_counts' layer found in adata.layers")

# 构建一个只包含 raw counts 的 AnnData，用来填充 .raw
# obs 轴默认和 adata 一致，var 用原来的 var
raw_adata = ad.AnnData(X=raw_mat, var=adata.var)

# 覆盖 .raw
adata.raw = raw_adata
print("Set adata.raw from selected layer.")

# 写出
print(f"Writing: {out_path}")
adata.write_h5ad(out_path)
print("Done.")
