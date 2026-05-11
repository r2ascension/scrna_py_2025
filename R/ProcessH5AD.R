# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
setwd("/home/h2048/data/R/1221")
library(reticulate)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
# Source the SCNT module functions (GetSeurat function)
# Make sure the SCNT.R file path is correct
# source("path/to/SCNT.R")

# ============================================================
# Configuration: Modify these parameters for your data
# ============================================================

# File paths for the two h5ad files
h5ad_file1 <- "/home/h2048/data/py/The integrated Human Lung Cell Atlas.h5ad"
h5ad_file2 <- "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"

# Metadata for each dataset
# Dataset 1 metadata
dataset1_name <- "Kerstin_B_Meyer_2021"
tissue_method1 <- "brush" # e.g., "FACS", "Dissociation", "Biopsy"

# Dataset 2 metadata
dataset2_name <- "Waradon_Sungnak_2020"
tissue_method2 <- "brush"

# ============================================================
# Step 1: Read first h5ad file
# ============================================================

# 检查counts slot (应该是raw counts)
counts <- GetAssayData(seurat_obj1, slot = "counts")
max_val <- max(counts)
mean_val <- mean(counts[counts > 0])

cat(sprintf("Max: %.2f\n", max_val))
cat(sprintf("Mean (non-zero): %.2f\n", mean_val))

if (max_val > 100) {
  cat("✓ HAS RAW COUNTS!\n")
} else if (max_val < 20) {
  cat("✗ LOG-TRANSFORMED\n")
} else {
  cat("⚠️  UNCLEAR\n")
}

# ============================================================
# Step 2: Read second h5ad file
# ============================================================

cat("\nReading second h5ad file...\n")
seurat_obj2 <- GetSeurat(h5ad_path = h5ad_file2, debug = TRUE)

# Add metadata columns for dataset 2
seurat_obj2$dataset <- dataset2_name
seurat_obj2$tissue_sampling_method <- tissue_method2

# Check the metadata
cat("\nSecond dataset metadata preview:\n")
print(head(seurat_obj2@meta.data))
cat(sprintf("\nSecond dataset: %d cells\n", ncol(seurat_obj2)))

# ============================================================
# Step 3: (Optional) Merge the two Seurat objects
# ============================================================

cat("\nMerging datasets...\n")
seurat_merged <- merge(
  x = seurat_obj1,
  y = seurat_obj2,
  add.cell.ids = c(dataset1_name, dataset2_name),
  project = "Merged_Analysis"
)

# Verify the merged object
cat("\nMerged dataset summary:\n")
print(table(seurat_merged$dataset))
print(table(seurat_merged$tissue_sampling_method))
cat(sprintf("\nTotal cells after merging: %d\n", ncol(seurat_merged)))

# Check metadata columns
cat("\nMetadata columns in merged object:\n")
print(colnames(seurat_merged@meta.data))

# ============================================================
# Step 4: Save the results
# ============================================================
cat("Reading first h5ad file...\n")
seurat_obj1 <- GetSeurat(h5ad_path = h5ad_file1, debug = TRUE)

seurat_obj1 <- subset(seurat_obj1, subset = tissue_sampling_method == 'brush')
table(seurat_obj1$tissue)
# Add metadata columns for dataset 1
# Without aggregation (keep duplicates with suffix)
DefaultAssay(seurat_obj1) <- 'RNA.symbol'
str(seurat_obj1@assays)
head(rownames(seurat_obj1), 20)
seurat_obj1@assays[[RNA]] <- NULL
seurat_obj1@assays$RNA <- seurat_obj1@assays$RNA.symbol
getwd()
run_tissue_comparison_analysis(
  seurat_obj1,
  cell_anno_col = "ann_level_1",
  tissue_col = "tissue",
  sample_col = "sample",
  min_cell_per_sample = 3,
  min_sample_per_tissue = 3,
  run_gsva = TRUE,
  run_go = TRUE,
  output_dir = "./analysis_results",
  species = "Homo sapiens"
)


seurat_obj1 <- NormalizeData(seurat_obj1) #归一化
seurat_obj1 <- FindVariableFeatures(
  seurat_obj1,
  selection.method = "vst",
  nfeatures = 4000
) #寻找变异基因
seurat_obj1 <- ScaleData(seurat_obj1) #标准化

seurat_obj1 <- RunPCA(seurat_obj1, npcs = 30)

getwd()
seurat_obj1 <- RunHarmony(
  obj = seurat_obj1,
  group.by.vars = 'sample',
  reduction.use = 'pca',
  dims.use = 1:30,
  # reduction.save = "harmony",
  verbose = TRUE
)
seurat_obj1 <- RunUMAP(
  object = seurat_obj1,
  reduction = "harmony",
  dims = 1:30,
  # reduction.name = 'umap_harmony',
  # reduction.key  = "umaph_",
  seed.use = 42,
  min.dist = 0.5,
  verbose = FALSE
)

cat("Neighbors + clustering...\n")
seurat_obj1 <- FindNeighbors(
  seurat_obj1,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)
seurat_obj1 <- FindClusters(seurat_obj1, res = 1, algorithm = 4, seed = 42)

pdf("qc_umap_harmony.pdf", width = 12, height = 9)
print(
  DimPlot(
    seurat_obj1,
    group.by = 'RNA.symbol_snn_res.1',
    label = TRUE,
    raster = TRUE
  ) +
    ggtitle(paste0("UMAP (Harmony) - ", 'RNA.symbol_snn_res.1'))
)
print(
  DimPlot(seurat_obj1, group.by = 'dataset', raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'dataset'))
)
dev.off()

# 0) 选择聚类/分组作为身份（按需改成你的meta列名，如 "cell_type"）
Idents(seurat_obj1) <- "RNA.symbol_snn_res.1"
marker_gene_1222 <- c(
  "LDB3",
  "SORBS1",
  "B7Z9B7",
  "MYL4",
  "SLMAP",
  "MCAM",
  "CAVIN3",
  "COL21A1",
  "TGFB1I1",
  "TNS1",
  "COL3A1",
  "ITGA1",
  "LPP",
  "SORBS3",
  "FERMT2",
  "PDLIM7",
  "RSU1",
  "SUSD2",
  "HEL-117",
  "PKP2",
  "EHD2",
  "SEPTIN10",
  "ATP6V0A1",
  "FLNC",
  "HEL32",
  "PPP1R12B",
  "PDLIM3",
  "PURB",
  "HEL114",
  "PALLD",
  "HNRNPD",
  "SSR1",
  "MAP4",
  "TLN1",
  "MYO1C"
)

# 仅保留对象中存在的基因
marker_gene_1222 <- intersect(rownames(seurat_obj1), marker_gene_1222)
marker_gene_1222
pdf("feature_plots_key_genes.pdf", width = 10, height = 8)

for (gene in marker_gene_1222) {
  p <- FeaturePlot(
    object = seurat_obj1,
    features = gene,
    reduction = "umap",
    raster = FALSE,
    pt.size = 0.5
  ) +
    ggtitle(paste("Expression of", gene)) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

  print(p)
}

dev.off()


pdf("qc_umap_harmony.pdf", width = 12, height = 9)
print(
  DimPlot(
    seurat_obj1,
    group.by = 'RNA.symbol_snn_res.1',
    label = TRUE,
    raster = TRUE
  ) +
    ggtitle(paste0("UMAP (Harmony) - ", 'RNA.symbol_snn_res.1'))
)
print(
  DimPlot(seurat_obj1, group.by = 'dataset', raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'dataset'))
)
print(
  DimPlot(seurat_obj1, group.by = 'tissue', raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'tissue'))
)
print(
  DimPlot(seurat_obj1, group.by = 'ann_finest_level', raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'ann_finest_level'))
)
print(
  DimPlot(seurat_obj1, group.by = 'ann_level_2', raster = TRUE) +
    ggtitle(paste0("UMAP (Harmony) - ", 'ann_level_2'))
)
dev.off()

saveRDS(seurat_obj1, 'Brush.rds')

# ggsave("marker_heatmap.pdf", p, width = 10, height = 12)
pdf("epi_marker.pdf", width = 12, height = 8, onefile = TRUE)
print(p)
dev.off()

# 确保有常见 QC 指标（如已存在会跳过）
if (!"percent.mt" %in% colnames(seurat_obj1@meta.data)) {
  seurat_obj1[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj1,
    pattern = "^MT-"
  )
}
if (!"percent.rb" %in% colnames(seurat_obj1@meta.data)) {
  seurat_obj1[["percent.rb"]] <- PercentageFeatureSet(
    seurat_obj1,
    pattern = "^(RPL|RPS)"
  )
}

cluster_col <- "RNA.symbol_snn_res.1"

pdf("QC_by_cluster.pdf", width = 12, height = 8, onefile = TRUE)

# 1) Violin：按 cluster 看 nCount / nFeature / mt / rb
p1 <- VlnPlot(
  seurat_obj1,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"),
  group.by = cluster_col,
  pt.size = 0.05,
  ncol = 2
)
print(p1)

# 2) Scatter：nCount vs nFeature（按 cluster 上色）
print(FeatureScatter(
  seurat_obj1,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = cluster_col
))

# 3) Scatter：percent.mt vs nCount / nFeature（帮助定位低质/高线粒体群）
print(FeatureScatter(
  seurat_obj1,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = cluster_col
))
print(FeatureScatter(
  seurat_obj1,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt",
  group.by = cluster_col
))

dev.off()

p0 <- VlnPlot(
  seurat_obj1,
  features = c('PTPRC', 'LST1', 'FCGR3B'),
  group.by = cluster_col,
  pt.size = 0.05,
  ncol = 2
)
print(p0)

pdf(
  "QC_by_cluster_oneplot_per_page.pdf",
  width = 12,
  height = 6,
  onefile = TRUE
)

qc_feats <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
for (f in qc_feats) {
  p <- VlnPlot(
    seurat_obj1,
    features = f,
    group.by = cluster_col,
    pt.size = 0
  )
  print(p) # 每次 print 自动新开一页
}

dev.off()

seurat_obj1$dataset <- dataset1_name
seurat_obj1$tissue_sampling_method <- tissue_method1

# Check the metadata
cat("\nFirst dataset metadata preview:\n")
print(head(seurat_obj1@meta.data))
cat(sprintf("\nFirst dataset: %d cells\n", ncol(seurat_obj1)))
seurat_obj1$tissue[seurat_obj1$tissue == 'nasal cavity'] <- 'nose'
# Save individual Seurat objects
seurat_obj1$sample <- seurat_obj1$sample_id
table(seurat_obj1$COVID_status, seurat_obj1$tissue)
table(seurat_obj1$disease)
seurat_obj1 <- subset(
  seurat_obj1,
  subset = tissue == "nasal cavity" &
    COVID_status %in% c("COVID+", "Post-COVID")
)

# seurat_obj1$tissue <- 'nose'
str(seurat_obj1)
saveRDS(seurat_obj1, file = "stromal_vascular_FINAL.rds")
seurat_obj2$tissue <- 'nose'
table(seurat_obj2$Location)
saveRDS(seurat_obj2, file = "Waradon_Sungnak_2020.rds")

# Save merged object
saveRDS(seurat_merged, file = "seurat_merged.rds")

# Optional: Save back to h5ad format using GetH5ad function
# GetH5ad(seurat_merged, "merged_output.h5ad", mode = "sc", assay = "RNA")

cat("\n=== Processing completed successfully ===\n")

# ============================================================
# Alternative: Read and process multiple datasets in a loop
# ============================================================

# If you have more than 2 datasets, you can use this approach:
process_multiple_h5ad <- function(file_list, metadata_list) {
  # file_list: character vector of h5ad file paths
  # metadata_list: list of lists containing metadata for each file
  # Example: list(
  #   list(dataset = "DS1", tissue_sampling_method = "FACS"),
  #   list(dataset = "DS2", tissue_sampling_method = "Biopsy")
  # )

  seurat_list <- list()

  for (i in seq_along(file_list)) {
    cat(sprintf(
      "\nProcessing file %d/%d: %s\n",
      i,
      length(file_list),
      basename(file_list[i])
    ))

    # Read h5ad file
    seurat_obj <- GetSeurat(h5ad_path = file_list[i], debug = FALSE)

    # Add metadata
    for (meta_key in names(metadata_list[[i]])) {
      seurat_obj[[meta_key]] <- metadata_list[[i]][[meta_key]]
    }

    seurat_list[[i]] <- seurat_obj
  }

  # Merge all objects
  if (length(seurat_list) > 1) {
    merged_obj <- Reduce(function(x, y) merge(x, y), seurat_list)
  } else {
    merged_obj <- seurat_list[[1]]
  }

  return(list(individual = seurat_list, merged = merged_obj))
}

# Example usage of the batch processing function:
# file_paths <- c("file1.h5ad", "file2.h5ad", "file3.h5ad")
# metadata_info <- list(
#   list(dataset = "CRSwNP_Site1", tissue_sampling_method = "FACS", disease = "CRSwNP"),
#   list(dataset = "CRSwNP_Site2", tissue_sampling_method = "Biopsy", disease = "CRSwNP"),
#   list(dataset = "Control_Site1", tissue_sampling_method = "FACS", disease = "Control")
# )
# results <- process_multiple_h5ad(file_paths, metadata_info)
# seurat_merged <- results$merged
