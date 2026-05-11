library(Seurat)
library(patchwork)
library(ggplot2)

obj <- seurat_obj1
DefaultAssay(obj) <- "RNA"


seurat_obj1 <- FindNeighbors(
  seurat_obj1,
  reduction = "scanvi", # ⭐ 使用scANVI降维
  dims = 1:75, # 使用全部75个维度
  k.param = 30, # 邻居数量
  # prune.SNN = 1/15,          # SNN修剪阈值
  verbose = TRUE
)

seurat_obj1 <- FindClusters(
  seurat_obj1,
  resolution = 1,
  algorithm = 4, # 4 = Leiden算法（推荐）
  random.seed = 42,
  verbose = FALSE
)


# Reductions(obj)
# # 建议确认你要用的分组列是否存在
# c("dataset","sample","batch","tissue","tissue_sampling_method","cell_type","scanvi_label","leiden_3") %in% colnames(obj@meta.data)

# (p1 <- DimPlot(obj, reduction="umap_scvi",   group.by="dataset", raster=TRUE) + ggtitle("scVI UMAP: dataset")) |
# (p2 <- DimPlot(obj, reduction="umap_scanvi", group.by="dataset", raster=TRUE) + ggtitle("scANVI UMAP: dataset"))

# (p3 <- DimPlot(obj, reduction="umap_scvi",   group.by="sample", raster=TRUE) + ggtitle("scVI UMAP: sample")) |
# (p4 <- DimPlot(obj, reduction="umap_scanvi", group.by="sample", raster=TRUE) + ggtitle("scANVI UMAP: sample"))

group_key <- "seurat_clusters" # 或 "scanvi_label" / "ann_level_2" / "leiden_3"

p1 <- DimPlot(
  obj,
  reduction = "umap_scanvi",
  group.by = group_key,
  label = FALSE,
  repel = TRUE
) +
  ggtitle("umap")

p1
summary(seurat_obj1$scanvi_confidence)
getwd()
Marker_Fibroblast <- c(
  'APOD',
  'FGF7',
  'COL15A1',
  'MFAP5',
  'PI16',
  'CD34',
  'MMP11',
  'COL10A1',
  'POSTN',
  'LRRC15',
  'HOPX',
  'IGFBP5',
  'TIMP1',
  'MMP1',
  'COL7A1',
  'WNT5A',
  'ISG15',
  'IL7R',
  'SFRP4',
  'SFRP2',
  'COMP',
  'RGS5',
  'PDGFRB',
  'NDUFA4L2',
  'NOTCH3',
  'CXCL1',
  'CXCL2',
  'IL6',
  'CEBPD',
  'CLU',
  'CTGF',
  'HGF',
  'HSPA6',
  'DNAJB1',
  'MYC',
  'AFT4',
  'PLAU',
  'CHI3L1',
  'MMP3',
  'IL1R1',
  'IL13RA2',
  'TNFSF11',
  'MMP10',
  'OSMR',
  'IL11',
  'STRA6',
  'FAP',
  'WNT2',
  'TWIST1',
  'IL24',
  'ACTG2',
  'HHIP',
  'CNN1',
  'MYH11',
  'ACTA2',
  'TAGLN',
  'KRT18',
  'SLPI',
  'UPK3B',
  'MSLN',
  'CALB2',
  'WT1',
  'KLK11',
  'ITLN1',
  'WSB1',
  'DDX17',
  'CTNNB1',
  'RBP1',
  'STAR',
  'STMN1',
  'CXCL12',
  'CD74',
  'HLA.DRB1',
  'HLA.DRA',
  'ADAMDEC1',
  'CCL8',
  'APOE',
  'APOC1',
  'LIMCH1',
  'A2M',
  'ADH1B',
  'PRG4',
  'CRTAC1',
  'CXCL14',
  'VSTM2A',
  'SOX6',
  'COL4A5',
  'COL4A6',
  'TSLP',
  'FRZB',
  'BMP5',
  'BMP2',
  'CPM',
  'F3'
)


seurat_obj1 <- NormalizeData(seurat_obj1) #归一化
seurat_obj1 <- FindVariableFeatures(
  seurat_obj1,
  selection.method = "vst",
  nfeatures = 2000
) #寻找变异基因

seurat_obj1 <- ScaleData(seurat_obj1) #标准化


print("Generating feature plots...")
pdf("umap_FeaturePlot.pdf", width = 8, height = 8)
for (marker in Marker_Fibroblast) {
  if (marker %in% rownames(seurat_obj1)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(seurat_obj1, features = marker, raster = TRUE))
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()
