library(BiocManager)
library(devtools)
# Sys.chmod("C:/Program Files/R/R-4.3.3/library",'777')
BiocManager::install('clusterProfiler', force = TRUE)
BiocManager::install("KEGGREST")
BiocManager::install('CellChat')
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DropletUtils")
BiocManager::install("Startrac")
BiocManager::install("SeuratDisk")
devtools::install_github('ropensci/magick')
devtools::install_github("YuLab-SMU/clusterProfiler")
devtools::install_github("aertslab/SCopeLoomR")
install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
# install.packages("Sesuit")
packageurl <-
  install.packages('ggalluvial')
remotes::install_github("PaulingLiu/ROGUE")
devtools::install_github("sajuukLyu/ggunchull", type = "source")
BiocManager::install("speckle")

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("mojaveazure/seurat-disk")

library(pypr)
installed.packages()[, c("Package", "LibPath")]
installed.packages()[, c("Package", "LibPath")]
remotes::install_github("PaulingLiu/ROGUE")
library(kBET)
usethis::edit_r_environ()

install.packages('magick')
devtools::install_github('vertesy/Seurat.utils')
remotes::install_github("LTLA/scuttle")
BiocManager::install("scuttle")
install.packages("hdf5r")
install.packages("remotes")
remotes::install_github("MarioniLab/DropletUtils")
BiocManager::install("Sesuit")
install.packages("ggrastr")

library(loupeR)
loupeR::setup()
library(Seurat.utils)
install.packages("Seurat.utils")
install.packages("reticulate")
install.packages("cowplot")
install.packages("DoubletFinder")

devtools::install_github("satijalab/seurat-data")
install.packages("tictoc")
SeuratData::InstallData("pbmc3k")
install.packages(
  "E:/idm/PaulingLiu-ROGUE-6e1c8f9.tar.gz",
  repos = NULL,
  type = "source"
)
SeuratData::InstallData("ifnb")
install.packages('NMF')
library("devtools")
install_github("Danko-Lab/BayesPrism/BayesPrism")
BiocManager::install("apeglm")
devtools::install_github("cellgeni/sceasy")
install.packages('apeglm')
install.packages('BiocManager')
BiocManager::install('SCENIC')

scenic_infercnv_packages <- c(
  "SCENIC",
  "SCopeLoomR",
  "AUCell",
  "RcisTarget",
  "doParallel",
  "infercnv",
  "Matrix",
  "ComplexHeatmap",
  "circlize",
  "HiddenMarkov",
  "foreach",
  "doSNOW",
  "parallel",
  "future",
  "BiocParallel",
  "arrow"
)
library(devtools)
install_github("campbio/celda")
for (packages in scenic_infercnv_packages) {
  BiocManager::install(packages)
}
install_github("Danko-Lab/BayesPrism/BayesPrism")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DropletUtils")
library(DropletUtils)

devtools::install_github("pengminshi/mrtree")
library(devtools)
devtools::install_local("/home/h2048/temp/ecotyper.zip")
devtools::install_local("/home/h2048/temp/BayesPrism.zip")
getOption('timeout')
# [1] 60
options(timeout = 6000)


# Install devtools (if not already installed)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install SCNT from GitHub
devtools::install_github("746443qjb/SCNT")
install.packages("bbknnR")
options(download.file.method = 'libcurl')
options(url.method = 'libcurl')
BiocManager::install(c(
  'BiocGenerics',
  'DelayedArray',
  'DelayedMatrixStats',
  'limma',
  'lme4',
  'S4Vectors',
  'SingleCellExperiment',
  'SummarizedExperiment',
  'batchelor',
  'HDF5Array',
  'ggrastr'
))
BiocManager::install("bnprks/BPCells/r")
install.packages('jsonlite')
BiocManager::install('Penghuihuang2000/BLEND')
devtools::install_github('Penghuihuang2000/BLEND')
install.packages("curl")
devtools::install_local('/home/h2048/temp/ggplot2_4.0.1.tar.gz')
devtools::install_local('/home/h2048/temp/ggstats_0.12.0.tar.gz')
devtools::install_local('/home/h2048/temp/GGally_2.4.0.tar.gz')
devtools::install_local('/home/h2048/temp/corpcor_1.6.10.tar.gz')


devtools::install_local('/home/h2048/temp/yulab.utils.zip')
devtools::install_local('/home/h2048/temp/GOSemSim.zip')
devtools::install_github("YuLab-SMU/enrichit")

devtools::install_local('/home/h2048/temp/fanyi_0.1.0.tar.gz')

devtools::install_local('/home/h2048/temp/enrichit-1.zip')
devtools::install_local('/home/h2048/temp/clusterProfiler-1.zip')

## SCENIC需要一些依赖包，先安装好
# install.package("BiocManager")
BiocManager::install(c(
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "zoo",
  "mixtools",
  "rbokeh",
  "DT",
  "NMF",
  "pheatmap",
  "R2HTML",
  "Rtsne",
  "doMC",
  "doRNG",
  "scRNAseq"
))
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")
BiocManager::install("Seurat")
#check
library(SCENIC)
packageVersion("SCENIC")


install.packages("reticulate")

BiocManager::install("TOAST")
devtools::install_github("YuLab-SMU/tigeR")
install.packages(c(
  "RColorBrewer",
  "cluster",
  "circlize",
  "cowplot",
  "data.table",
  "doParallel",
  "ggplot2",
  "grid",
  "reshape2",
  "viridis",
  "config",
  "argparse",
  "colorspace",
  "plyr"
))
BiocManager::install("ComplexHeatmap")
BiocManager::install("Biobase")
BiocManager::install("NMF")
