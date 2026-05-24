#!/usr/bin/env python3
import anndata as ad
import numpy as np
import os

# ======== 路径配置（按你给的写死）========
FULL_PATH = "/home/h2048/data/py/1204/Tcell_pure_scvi/starcat_analysis/adata_tcell_scvi_starcat_annotated.h5ad"
SCANVI_PATH = "/home/h2048/data/py/1204/Tcell_scANVI/adata_tcell_scANVI_annotated.h5ad"

# 写回同一个文件的话，设为 SCANVI_PATH；也可以改成新文件名
OUTPUT_PATH = "/home/h2048/data/py/1204/Tcell_scANVI/adata_tcell_scANVI_annotated_fullRaw.h5ad"


def main():
    print("Reading full (starcat) file:")
    print("  ", FULL_PATH)
    adata_full = ad.read_h5ad(FULL_PATH)

    print("Reading scANVI file:")
    print("  ", SCANVI_PATH)
    adata_scanvi = ad.read_h5ad(SCANVI_PATH)

    # ============ 1. 对齐细胞顺序（你说细胞完全一致，但仍然稳妥检查一下） ============
    if not np.array_equal(adata_full.obs_names, adata_scanvi.obs_names):
        print("Cell order differs, reindexing full to match scANVI ...")
        # 按 scANVI 的 obs 顺序重排 full
        missing = np.setdiff1d(adata_scanvi.obs_names, adata_full.obs_names)
        if len(missing) > 0:
            raise ValueError(f"Cells in scANVI not found in full object: {len(missing)} cells (e.g. {missing[:5]})")
        adata_full = adata_full[adata_scanvi.obs_names, :].copy()
    else:
        print("Cell order is identical between full and scANVI.")

    # ============ 2. 从 full 取完整的 counts 矩阵 ============
    layer_keys = list(adata_full.layers.keys())
    print("Full object layers:", layer_keys)

    if "counts" in layer_keys:
        print("Using adata_full.layers['counts'] as full counts.")
        full_counts = adata_full.layers["counts"]
    elif "raw_counts" in layer_keys:
        print("Using adata_full.layers['raw_counts'] as full counts.")
        full_counts = adata_full.layers["raw_counts"]
    elif adata_full.raw is not None:
        print("No 'counts'/'raw_counts' layer; using adata_full.raw.X as full counts.")
        full_counts = adata_full.raw.X
    else:
        print("No 'counts'/'raw_counts'/raw found; using adata_full.X as full counts (⚠️ 可能是归一化矩阵).")
        full_counts = adata_full.X

    full_var = adata_full.var.copy()
    full_var_names = full_var.index.to_numpy()

    print(f"Full counts shape: {full_counts.shape[0]} cells x {full_counts.shape[1]} genes")

    # ============ 3. 写入 scANVI 的 raw（全基因） ============
    print("Setting adata_scanvi.raw with full counts and full var ...")
    adata_scanvi.raw = ad.AnnData(X=full_counts, var=full_var)

    # ============ 4. 写入 scANVI 的 layers['counts']（按 scANVI HVG 子集） ============
    scanvi_var_names = adata_scanvi.var_names.to_numpy()
    print(f"scANVI HVG genes: {len(scanvi_var_names)}")

    # 确保 HVG 都在 full 里
    missing_genes = np.setdiff1d(scanvi_var_names, full_var_names)
    if len(missing_genes) > 0:
        raise ValueError(
            f"{len(missing_genes)} scANVI genes not found in full object (e.g. {missing_genes[:5]}). "
            "Check that the two files come from the same reference."
        )

    # 按 scANVI 基因顺序在 full_counts 里取列
    # full_var.index 是 pandas Index，get_indexer 保证顺序一致
    indexer = full_var.index.get_indexer(scanvi_var_names)
    if (indexer < 0).any():
        raise ValueError("Some scANVI genes could not be indexed in full_var.")

    print("Subsetting full counts to scANVI HVG gene set for layers['counts'] ...")
    # full_counts 是 (n_cells x n_full_genes)，取 subset 列
    hvg_counts = full_counts[:, indexer]

    print(f"hvg_counts shape (for layers['counts']): {hvg_counts.shape[0]} cells x {hvg_counts.shape[1]} genes")

    adata_scanvi.layers["counts"] = hvg_counts

    # ============ 5. 保存 ============
    print("Writing updated scANVI AnnData to:")
    print("  ", OUTPUT_PATH)
    # 如果原文件存在并且 OUTPUT_PATH == SCANVI_PATH，就覆盖
    if os.path.exists(OUTPUT_PATH) and OUTPUT_PATH != SCANVI_PATH:
        print("⚠️  WARNING: output file exists and will be overwritten.")

    adata_scanvi.write_h5ad(OUTPUT_PATH)
    print("Done.")


if __name__ == "__main__":
    main()
