# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(rsvd)
library(harmony)
library(presto)
Sys.setenv(RETICULATE_PYTHON = "C:/Users/崔填祎/AppData/Local/Programs/Python/Python311/python.exe")
library(reticulate)
library(DoubletFinder)
library(DirichletReg)
library(speckle)
library(tidyverse)
library(MAST)
require(DirichletReg)
require(speckle)
require(ggplot2)
require(cowplot)
require(edgeR)
require(reshape2)
require(pheatmap)
require(tidyr)
library(patchwork)
Markers <- c('FXYD3', 'EPCAM', 'ELF3', 'SERPINF1', 'TSPAN1',
             'SCGB1A1', 'AGER', 'SFTPC', 'FOXJ1', 'KRT5', 'MUC5B', 'KRT8',#Epithelial
             'CD53', 'PTPRC', 'CORO1A', 'CCL5',
             'MS4A1', 'TNFRSF17', 'CD19', 'CD79A', #B
             'CD40LG', 'TNFRSF25', 'CD28', 'CD4', 'CD3E', 'CD8A', 'CD8B', 
             'TRGC2', 'CD2', 'TRBC2', #T
             'FCER1G', 'C1orf162', 'CLEC7A', 'CD1C', 'CD86', 'CD14', 'XCR1', 'HLA-DRA',#Myeloid
             'COL1A2', 'DCN', 'MFAP4', 'LUM', 'COL6A3', 'CFD', 'COL1A1', 'PDGFRA', 
             'MXRA8', 'LEPR',
             'MYH11', 'TINAGL1', 'PLN', 'DES', 'ACTA2', 'CNN1', 'TAGLN', #Fibroblast & SMC
             'CLDN5', 'ECSCR', 'CLEC14A', 'VWF', 'PECAM1', 'DARC', 'PTPRB', 'PDE2A', 
             'PLAT', 'GJA5', 'SPARCL1', 'AQP1', 'MMRN1', 'CCL21', #Endothelial
             'MKI67', 'TOP2A', 'TK1', 'CENPW')#Proliferation

Markers <- c('FXYD3', 'EPCAM', 'ELF3', 'IGFBP2', 'SERPINF1', 'TSPAN1',
             'SCGB1A1', 'AGER', 'SFTPC', 'FOXJ1', 'KRT5', 'MUC5B', 'KRT8',
             'CD53', 'PTPRC', 'CORO1A', 'ISG20', 'CCL5',
             'MS4A1', 'TNFRSF17', 'CD19', 'CD79A', 'SDC1',
             'CD40LG', 'TNFRSF25', 'CD28', 'CD4', 'CD3E', 'CD8A', 'CD8B', 
             'TRGC2', 'CD2', 'TRBC2',
             'FCER1G', 'C1orf162', 'CLEC7A', 'CD1C', 'CD86', 'CD14', 'XCR1', 'HLA-DRA',
             'COL1A2', 'DCN', 'MFAP4', 'LUM', 'COL6A3', 'CFD', 'COL1A1', 'PDGFRA', 
             'MXRA8', 'NBL1', 'VCAN', 'LEPR',
             'MYH11', 'TINAGL1', 'PLN', 'DES', 'ACTA2', 'CNN1', 'TAGLN',
             'CLDN5', 'ECSCR', 'CLEC14A', 'VWF', 'PECAM1', 'DARC', 'PTPRB', 'PDE2A', 
             'PLAT', 'GJA5', 'SPARCL1', 'AQP1', 'RNASE1', 'MMRN1', 'CCL21', 'TFF3',
             'MKI67', 'TOP2A', 'TK1', 'CENPW')
Marker_Myeloid <- c('CLEC9A','XCR1','CADM1','CLNK','FLT3','ZBTB46',
                    'CLEC10A','CD1E','FCER1A','CD1D','ITGAX','CDIC','FCGR2B','PKIB',
                    'CCR7','CD83','LAMP3','CCL22','CCL17','CCL19','LAD1',
                    'LILRA4','SMPD3','SCT','IRF7','PLD4','CLEC4C',
                    'MARCO','FABP4','CYP27A1','SIGLEC1','ABCG1','PPARG',
                    'C1QA','C1QB','C1QC','HLA-DPA1','SLC40A1',
                    'FOLR2','F13A1',
                    'SPP1','HAMP','VCAN','CCR2','CCR5',
                    'FCN1',
                    'S100A12','RNASE2',
                    'LILRA5','MTSS1',
                    'TPSAB1','MS4A2','TPSB2',
                    'FCGR3B','CSF3R','CXCR1')
Marker_Endothelial <- c('S100B','ALDH1A1',
                        'GJA4','HEY1','DKK2','IGFBP3','EFNB2',
                        'ACKR1','VWF',
                        'RGCC','VWA1','IL7R','FCN3','MT1M',
                        'TFF3','LYVE1','MMRN1','CCL21',
                        'MYL9','ACTA2','TINAGL1','NOTCH3','LAMC3')
Marker_Fibroblast <- c('APOD','FGF7','COL15A1','MFAP5','PI16','CD34',
                       'MMP11','COL10A1','POSTN','LRRC15','HOPX','IGFBP5','TIMP1','MMP1','COL7A1','WNT5A','ISG15','IL7R','SFRP4','SFRP2','COMP','RGS5','PDGFRB','NDUFA4L2','NOTCH3',
                       'CXCL1','CXCL2','IL6','CEBPD','CLU','CTGF','HGF','HSPA6','DNAJB1','MYC','AFT4','PLAU','CHI3L1','MMP3','IL1R1','IL13RA2','TNFSF11','MMP10','OSMR','IL11','STRA6','FAP','WNT2','TWIST1','IL24',
                       'ACTG2','HHIP','CNN1',
                       'MYH11','ACTA2','TAGLN',
                       'KRT18','SLPI','UPK3B','MSLN','CALB2','WT1','KLK11','ITLN1',
                       'WSB1','DDX17','CTNNB1',
                       'RBP1','STAR','STMN1',
                       'CXCL12','CD74','HLA.DRB1','HLA.DRA',
                       'ADAMDEC1','CCL8','APOE','APOC1',
                       'LIMCH1','A2M','ADH1B',
                       'PRG4','CRTAC1',
                       'CXCL14','VSTM2A','SOX6','COL4A5','COL4A6','TSLP','FRZB','BMP5','BMP2','CPM','F3')
Marker_Epithelial <- c('TP63','KRT5',
                       'SCGB1A1','SERPINB3','SCGB3A2','SCGB3A1','TCN1','ASRGL1',
                       'FOXJ1','RSPH1','PIFO','BEST4','C20orf85','C9orf24',
                       'MUC5AC','SPDEF','LYPD2','ITLN1',
                       'ASCL1','GRP',
                       'POU2F3','ASCL2','CFTR','FOXI1','ASCL3','BSND','IGF1','CLCNKB',
                       'AGER','RTKN2','CLIC5','SPOCK2','TIMP3',
                       'SFTPC','LAMP3','MF5D2A','C8orf4','C11orf96',
                       'VIM','SOX9',
                       'KRT14','MYH11','ACTA2',
                       'DMBT1','RNASE1',
                       'MUC5B','SPDEF',
                       'LYZ','LTF',
                       'SFTPB','SCGB3A2','SFTA2')
Marker_T <- c('KLRD1','FCGR3A','GNLY','TYROBP','FCER1G','KLRC1','FGFBP2','SPON2','MYOM2',
               'TRDC','KRT86',
               'GATA3','IL5','AREG','HPGDS',
               'IL23R','RORC','LST1','PCDH9','TNFSF11',
               'TRDC','TRGC1','TRGC2',
               'CD4','CD28','CD40LG','TRAT1','TNFRSF25',
               'CD8A','CD8B','TRGC2',
               'CCR7','TCF7','LEF1','SELL',
               'CD28','IL7R','CCR6',
               'GATA3','IL4','IL13',
               'IL17A','CCL20',
               'CCR7','TCF7','LEF1','SELL',
               'TBX21','EOMES','GZMK','KLRG1',
               'ITGA1','CD8A','CD8B','CCR6','IL7R',
               'TBX21','GZMB','GZMH','FGFBP2',
               'ZNF683','IFNG','CCL4L2','PDCD1',
               'KLRB1','IL7R','NCR3','CEBPD',
               'FOXP3','IL2RA','IKZF2','TNFRSF4',
               'MKI67','TOP2A','TK1','CENPW')
Marker_SMC <- c('RGS5','CD36','NOTCH3',
                 'SCIN','EPAS1','HLA.C','IGKC','PTP4A3','FN1','COL18A1','WFDC1','IGHG4','IGHG1',
                 'RERGL','PLN','SORBS2','DSTN','TSC22D1','BCAM','C11orf96','FILIP1','WTIP','NRGN',
                 'TMEM176B','ANGPTL1','CFH','FHL1','VCAN','GGT5','C1S','COL6A3','STEAP4','ADGRL3')
Marker_B <- c('IGHD','TCL1A',
               'MS4A1','BANK1',
               'MKI67','TOP2A',
               'IGHA1','IGHA2',
               'IGHGP','IGHG1')
Marker_Stromal <- c('APOD','FGF7','COL15A1','MFAP5','PI16','CD34',
                       'MMP11','COL10A1','POSTN','LRRC15','HOPX','IGFBP5','TIMP1','MMP1','COL7A1','WNT5A','ISG15','IL7R','SFRP4','SFRP2','COMP','RGS5','PDGFRB','NDUFA4L2','NOTCH3',
                       'CXCL1','CXCL2','IL6','CEBPD','CLU','CTGF','HGF','HSPA6','DNAJB1','MYC','AFT4','PLAU','CHI3L1','MMP3','IL1R1','IL13RA2','TNFSF11','MMP10','OSMR','IL11','STRA6','FAP','WNT2','TWIST1','IL24',
                       'ACTG2','HHIP','CNN1',
                       'MYH11','ACTA2','TAGLN',
                       'KRT18','SLPI','UPK3B','MSLN','CALB2','WT1','KLK11','ITLN1',
                       'WSB1','DDX17','CTNNB1',
                       'RBP1','STAR','STMN1',
                       'CXCL12','CD74','HLA.DRB1','HLA.DRA',
                       'ADAMDEC1','CCL8','APOE','APOC1',
                       'LIMCH1','A2M','ADH1B',
                       'PRG4','CRTAC1',
                       'CXCL14','VSTM2A','SOX6','COL4A5','COL4A6','TSLP','FRZB','BMP5','BMP2','CPM','F3',
                      'RGS5','CD36','NOTCH3',
                      'SCIN','EPAS1','HLA.C','IGKC','PTP4A3','FN1','COL18A1','WFDC1','IGHG4','IGHG1',
                      'RERGL','PLN','SORBS2','DSTN','TSC22D1','BCAM','C11orf96','FILIP1','WTIP','NRGN',
                      'TMEM176B','ANGPTL1','CFH','FHL1','VCAN','GGT5','C1S','COL6A3','STEAP4','ADGRL3')
Marker_Immune <- c('CLEC9A','XCR1','CADM1','CLNK','FLT3','ZBTB46',
                    'CLEC10A','CD1E','FCER1A','CD1D','ITGAX','CDIC','FCGR2B','PKIB',
                    'CCR7','CD83','LAMP3','CCL22','CCL17','CCL19','LAD1',
                    'LILRA4','SMPD3','SCT','IRF7','PLD4','CLEC4C',
                    'MARCO','FABP4','CYP27A1','SIGLEC1','ABCG1','PPARG',
                    'C1QA','C1QB','C1QC','HLA-DPA1','SLC40A1',
                    'FOLR2','F13A1',
                    'SPP1','HAMP','VCAN','CCR2','CCR5',
                    'FCN1',
                    'S100A12','RNASE2',
                    'LILRA5','MTSS1',
                    'TPSAB1','MS4A2','TPSB2',
                    'FCGR3B','CSF3R','CXCR1',
                   'KLRD1','FCGR3A','GNLY','TYROBP','FCER1G','KLRC1','FGFBP2','SPON2','MYOM2',
                   'TRDC','KRT86',
                   'GATA3','IL5','AREG','HPGDS',
                   'IL23R','RORC','LST1','PCDH9','TNFSF11',
                   'TRDC','TRGC1','TRGC2',
                   'CD4','CD28','CD40LG','TRAT1','TNFRSF25',
                   'CD8A','CD8B','TRGC2',
                   'CCR7','TCF7','LEF1','SELL',
                   'CD28','IL7R','CCR6',
                   'GATA3','IL4','IL13',
                   'IL17A','CCL20',
                   'CCR7','TCF7','LEF1','SELL',
                   'TBX21','EOMES','GZMK','KLRG1',
                   'ITGA1','CD8A','CD8B','CCR6','IL7R',
                   'TBX21','GZMB','GZMH','FGFBP2',
                   'ZNF683','IFNG','CCL4L2','PDCD1',
                   'KLRB1','IL7R','NCR3','CEBPD',
                   'FOXP3','IL2RA','IKZF2','TNFRSF4',
                   'MKI67','TOP2A','TK1','CENPW',
                   'IGHD','TCL1A',
                   'MS4A1','BANK1',
                   'MKI67','TOP2A',
                   'IGHA1','IGHA2',
                   'IGHGP','IGHG1')
setwd("E:/R/0301/")
#' 优化的降维和聚类函数，简化为只使用Harmony进行批次校正
#' @param seurat_obj Seurat对象
#' @param cell_type 要分析的细胞类型
#' @param markers 预设的marker基因列表，默认NULL
#' @param output_dir 输出目录，默认"."
#' @param harmony 是否使用harmony校正批次效应，默认TRUE
#' @param resolution 聚类分辨率，默认2
#' @param cell_type_col 细胞类型所在列名，默认"Annotation"
#' @param harmony_vars 用于harmony批次校正的变量，默认c("sample", "tissue")
#' @param max_cells_for_heatmap 热图最大细胞数，默认3000
#' @return 处理后的Seurat对象
run_subset_analysis <- function(seurat_obj, 
                                cell_type, 
                                markers = NULL, 
                                output_dir = ".", 
                                harmony = TRUE,
                                resolution = 2,
                                cell_type_col = "Annotation",
                                harmony_vars = c("sample", "tissue"),
                                max_cells_for_heatmap = 3000) {
  
  # 设置输出目录
  if(!is.null(markers)) markers <- unique(markers)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_prefix <- file.path(output_dir, paste0(cell_type, "_"))
  
  # 检查是否已存在分析文件
  processed_rds <- paste0(output_prefix, "processed.rds")
  analyzed_rds <- paste0(output_prefix, "analyzed.rds")
  
  message(paste0("Starting analysis for cell type: ", cell_type))
  
  # 如果存在已处理的RDS，直接加载
  if(file.exists(processed_rds)) {
    message(paste0("Loading existing processed data from: ", processed_rds))
    subset_obj <- readRDS(processed_rds)
  } else {
    # 打印调试信息
    message(paste0("Looking for cells where ", cell_type_col, " = '", cell_type, "'"))
    message(paste0("Total cells before subsetting: ", ncol(seurat_obj)))
    cell_count_pre <- sum(seurat_obj@meta.data[[cell_type_col]] == cell_type, na.rm = TRUE)
    message(paste0("Cells matching this condition: ", cell_count_pre))
    
    # 直接通过元数据选择细胞
    cells_to_keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cell_type_col]] == cell_type]
    
    if(length(cells_to_keep) == 0) {
      stop(paste0("No cells found for cell type: ", cell_type))
    }
    
    # 使用cells参数进行子集选择
    subset_obj <- subset(seurat_obj, cells = cells_to_keep)
    
    # 检查子集中的细胞数量
    cell_count <- ncol(subset_obj)
    message(paste0("Cells after subsetting: ", cell_count))
    
    if(cell_count < 10) {
      warning(paste0("Cell type '", cell_type, "' has only ", cell_count, " cells. Analysis may not be reliable."))
      if(cell_count == 0) {
        stop(paste0("No cells found for cell type: ", cell_type))
      }
    }
    
    message(paste0("Extracted ", cell_count, " cells for analysis."))
    
    # 运行降维和聚类 - 使用全部细胞
    subset_obj <- run_dim_reduction(subset_obj, 
                                    harmony = harmony, 
                                    vars_use = harmony_vars,
                                    resolution = resolution)
    
    # 保存降维和聚类结果
    message(paste0("Saving processed data to: ", processed_rds))
    saveRDS(subset_obj, processed_rds)
  }
  
  # 如果存在已分析的RDS，直接加载
  if(file.exists(analyzed_rds)) {
    message(paste0("Loading existing analysis from: ", analyzed_rds))
    return(readRDS(analyzed_rds))
  }
  
  # 运行分析 - 使用优化版本的run_analysis函数
  subset_obj <- run_optimized_analysis(subset_obj, 
                                       markers = markers, 
                                       output_prefix = output_prefix,
                                       max_cells_for_heatmap = max_cells_for_heatmap)
  
  # 保存分析结果
  message(paste0("Saving analyzed results to: ", analyzed_rds))
  saveRDS(subset_obj, analyzed_rds)
  
  message(paste0("Analysis completed for cell type: ", cell_type))
  return(subset_obj)
}

#' 优化的数据分析函数 - 内存友好版本，仅在可视化时抽样
#'
#' @param seurat_obj Seurat对象
#' @param markers 预设的marker基因列表，默认NULL
#' @param output_prefix 输出文件前缀
#' @param dot_plot_cols 点图的颜色设置，默认c("lightgrey", "red")
#' @param find_markers_params 寻找marker基因的参数列表
#' @param max_cells_for_heatmap 热图最大细胞数，默认3000
#' @return 处理后的Seurat对象
run_optimized_analysis <- function(seurat_obj, 
                                   markers = NULL, 
                                   output_prefix = "",
                                   dot_plot_cols = c("lightgrey", "red"),
                                   find_markers_params = list(only.pos = TRUE, 
                                                              min.pct = 0.25, 
                                                              logfc.threshold = 0.25),
                                   max_cells_for_heatmap = 3000) {
  
  # 确保输出目录存在
  output_dir <- dirname(output_prefix)
  if(output_dir != "" && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 首先输出UMAP聚类图
  message("Generating UMAP cluster plots...")
  tryCatch({
    pdf(paste0(output_prefix, "umap_clusters.pdf"), width = 10, height = 8)
    print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
            ggtitle("UMAP Clustering"))
    
    # 如果数据中有组织或条件信息，也绘制按组织/条件着色的UMAP
    if("Tissue" %in% colnames(seurat_obj@meta.data)) {
      print(DimPlot(seurat_obj, reduction = "umap", group.by = "Tissue") + 
              ggtitle("UMAP by Tissue"))
    }
    
    if("condition" %in% colnames(seurat_obj@meta.data)) {
      print(DimPlot(seurat_obj, reduction = "umap", group.by = "condition") + 
              ggtitle("UMAP by Condition"))
    }
    
    if("tissue" %in% colnames(seurat_obj@meta.data)) {
      print(DimPlot(seurat_obj, reduction = "umap", group.by = "tissue") + 
              ggtitle("UMAP by Tissue"))
    }
    dev.off()
  }, error = function(e) {
    message(paste0("Error generating UMAP plots: ", e$message))
  })
  
  # Marker基因分析 - 使用全部细胞
  if(!is.null(markers)) {
    # 过滤掉不在数据中的marker基因
    valid_markers <- markers[markers %in% rownames(seurat_obj)]
    if(length(valid_markers) > 0) {
      message(paste("Analyzing", length(valid_markers), "marker genes..."))
      
      # 生成Feature plots
      tryCatch({
        feature_plot_width <- min(10, 5 * ceiling(sqrt(length(valid_markers))))
        feature_plot_height <- min(8, 4 * ceiling(length(valid_markers) / ceiling(sqrt(length(valid_markers)))))
        
        pdf(paste0(output_prefix, "umap_FeaturePlot.pdf"), width = feature_plot_width, height = feature_plot_height)
        # 每页最多显示6个基因
        for(i in seq(1, length(valid_markers), by = 6)) {
          end_idx <- min(i + 5, length(valid_markers))
          if(i <= end_idx) {  # 防止索引错误
            current_markers <- valid_markers[i:end_idx]
            p <- FeaturePlot(seurat_obj, 
                             features = current_markers,
                             raster = TRUE,
                             ncol = min(3, length(current_markers)))
            print(p)
          }
        }
        dev.off()
      }, error = function(e) {
        message(paste0("Error generating feature plots: ", e$message))
      })
      
      # DotPlot - 使用全部细胞
      tryCatch({
        dot_plot_width <- max(10, length(valid_markers) * 0.3)
        pdf(paste0(output_prefix, "markers_dotplot.pdf"), width = dot_plot_width, height = 10)
        print(DotPlot(seurat_obj, 
                      features = valid_markers, 
                      group.by = "seurat_clusters",
                      cols = dot_plot_cols,
                      dot.scale = 8) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        dev.off()
      }, error = function(e) {
        message(paste0("Error generating dot plot: ", e$message))
      })
      
      # 计算每个cluster中每个marker基因的表达统计 - 使用全部细胞进行计算
      message("Calculating marker gene expression statistics...")
      tryCatch({
        results <- data.frame()
        
        for(cluster in unique(Idents(seurat_obj))) {
          cells_in_cluster <- WhichCells(seurat_obj, idents = cluster)
          
          for(gene in valid_markers) {
            expr_data <- GetAssayData(seurat_obj, layer = "data")[gene, cells_in_cluster]
            mean_expr <- mean(expr_data)
            mean_expr_nonzero <- if(sum(expr_data > 0) > 0) mean(expr_data[expr_data > 0]) else 0
            
            pct_expr <- sum(expr_data > 0) / length(expr_data) * 100
            
            # 使用rbind而不是列表提高效率
            results <- rbind(results, data.frame(
              Cluster = cluster,
              Gene = gene,
              Mean_Expression = round(mean_expr, 3),
              Mean_Expression_NonZero = round(mean_expr_nonzero, 3),
              Percent_Expressing = round(pct_expr, 2)
            ))
          }
        }
        
        write.csv(results, paste0(output_prefix, "cluster_marker_stats.csv"), row.names = FALSE)
      }, error = function(e) {
        message(paste0("Error calculating marker statistics: ", e$message))
      })
    } else {
      warning("None of the provided markers were found in the dataset.")
    }
  }
  
  # 找marker基因 - 使用全部细胞进行计算
  cluster_markers_file <- paste0(output_prefix, "cluster_markers.csv")
  if(!file.exists(cluster_markers_file)) {
    message("Finding cluster-specific marker genes...")
    tryCatch({
      # 使用do.call允许传递自定义参数列表
      markers_args <- c(list(object = seurat_obj), find_markers_params)
      markers <- do.call(FindAllMarkers, markers_args)
      
      write.csv(markers, cluster_markers_file, row.names = FALSE)
    }, error = function(e) {
      message(paste0("Error finding markers: ", e$message))
      # 创建空的结果文件避免重复尝试
      write.csv(data.frame(), cluster_markers_file, row.names = FALSE)
      return(seurat_obj)
    })
  } else {
    message(paste0("Loading existing marker genes from: ", cluster_markers_file))
    markers <- read.csv(cluster_markers_file)
  }
  
  # 生成热图 - 仅在热图可视化时进行抽样
  if(file.exists(cluster_markers_file) && file.size(cluster_markers_file) > 0) {
    if(!exists("markers")) markers <- read.csv(cluster_markers_file)
    
    if(nrow(markers) > 0) {
      message("Generating marker gene heatmaps...")
      
      # 对每个群集取前5个marker基因(减少基因数量)
      tryCatch({
        top_markers <- markers %>% 
          group_by(cluster) %>% 
          top_n(5, wt = avg_log2FC)  # 减少为每个cluster的前5个
        
        marker_genes <- unique(top_markers$gene)
        
        if(length(marker_genes) > 0) {
          # 仅为热图创建抽样对象
          cell_count <- ncol(seurat_obj)
          
          if(cell_count > max_cells_for_heatmap) {
            message(paste0("Large dataset detected (", cell_count, " cells). Creating sampled object with ", 
                           max_cells_for_heatmap, " cells for heatmap visualization only."))
            
            # 按聚类进行分层抽样
            set.seed(42) # 设置随机种子确保可重复性
            
            # 获取每个聚类的细胞
            cells_by_cluster <- split(colnames(seurat_obj), seurat_obj$seurat_clusters)
            
            # 计算每个聚类应抽取的细胞数量（按比例）
            cells_per_cluster <- lapply(cells_by_cluster, function(cells) {
              n_sample <- ceiling(length(cells) / cell_count * max_cells_for_heatmap)
              if(length(cells) <= n_sample) return(cells)
              return(sample(cells, n_sample))
            })
            
            # 合并所有抽样细胞
            sampled_cells <- unlist(cells_per_cluster)
            
            # 创建子集用于热图
            heatmap_obj <- subset(seurat_obj, cells = sampled_cells)
            message(paste0("Created temporary object with ", ncol(heatmap_obj), " cells for heatmap visualization"))
          } else {
            heatmap_obj <- seurat_obj
          }
          
          # 仅缩放热图所需的基因
          heatmap_obj <- ScaleData(heatmap_obj, features = marker_genes)
          
          # 减小批次大小，每批次最多15个基因
          batch_size <- 15  # 热图中一次显示的基因数
          marker_batches <- split(marker_genes, ceiling(seq_along(marker_genes)/batch_size))
          
          for(i in seq_along(marker_batches)) {
            batch_genes <- marker_batches[[i]]
            if(length(batch_genes) > 0) {
              # 设置适当的高度
              pdf_height <- max(8, length(batch_genes) * 0.2)
              
              # 使用tryCatch捕获错误
              tryCatch({
                message(paste0("Creating heatmap batch ", i, " with ", length(batch_genes), " genes"))
                
                heatmap_file <- paste0(output_prefix, "marker_heatmap_batch_", i, ".pdf")
                pdf(heatmap_file, width = 10, height = pdf_height)
                
                # 使用raster=TRUE提高性能，减小点大小
                p <- DoHeatmap(heatmap_obj, 
                               features = batch_genes, 
                               group.by = "seurat_clusters",
                               size = 3,       
                               raster = TRUE)  # 使用光栅化提高性能
                
                print(p)
                dev.off()
                
                # 清理内存
                rm(p)
                gc()
                
              }, error = function(e) {
                message(paste0("Error generating heatmap batch ", i, ": ", e$message))
                # 如果文件已经创建但不完整，删除它
                if(file.exists(paste0(output_prefix, "marker_heatmap_batch_", i, ".pdf"))) {
                  file.remove(paste0(output_prefix, "marker_heatmap_batch_", i, ".pdf"))
                }
              })
            }
          }
          
          # 清理临时对象释放内存
          if(exists("heatmap_obj") && !identical(heatmap_obj, seurat_obj)) {
            rm(heatmap_obj)
            gc()
          }
        }
      }, error = function(e) {
        message(paste0("Error in heatmap preparation: ", e$message))
      })
    }
  }
  
  return(seurat_obj)
}
# 设置工作目录

# 1. 读取所有RDS文件
batch_files <- list.files("E:\\R\\0216\\rds", pattern = "\\.rds$", full.names = TRUE)
if(length(batch_files) == 0) {
  stop("No RDS files found!")
}

# 2. 读取数据集
merged_list <- lapply(batch_files, function(file) {
  cat(sprintf("Processing: %s\n", basename(file)))
  current_batch <- readRDS(file)
  DietSeurat(current_batch, dimreducs = NULL, graphs = NULL)
})

# 3. 获取所有基因
all_genes <- unique(unlist(lapply(merged_list, rownames)))
cat("Total unique genes:", length(all_genes), "\n")

# 4. 对每个样本进行基因填充
filled_list <- list()
for(i in seq_along(merged_list)) {
  obj <- merged_list[[i]]
  missing_genes <- setdiff(all_genes, rownames(obj))
  
  if(length(missing_genes) > 0) {
    cat(sprintf("Batch %d: Adding %d missing genes\n", i, length(missing_genes)))
    
    # 创建填充矩阵
    fill_mat <- Matrix::Matrix(0, 
                               nrow = length(missing_genes), 
                               ncol = ncol(obj), 
                               sparse = TRUE)
    rownames(fill_mat) <- missing_genes
    colnames(fill_mat) <- colnames(obj)
    
    # 添加到原矩阵
    new_counts <- rbind(GetAssayData(obj, layer = "counts"), fill_mat)
    new_counts <- new_counts[all_genes,]  # 确保基因顺序一致
    obj[["RNA"]] <- CreateAssayObject(counts = new_counts)
  }
  filled_list[[i]] <- obj
}

# 5. 处理每个批次
processed_list <- list()
for(i in seq_along(filled_list)) {
  cat(sprintf("\nProcessing batch %d/%d\n", i, length(filled_list)))
  
  obj <- filled_list[[i]]
  cat(sprintf("Cells: %d, Features: %d\n", ncol(obj), nrow(obj)))
  
  # 标准化
  cat("Performing normalization...\n")
  obj <- NormalizeData(obj)
  
  # 找变异基因
  cat("Finding variable features...\n")
  obj <- FindVariableFeatures(obj, 
                              selection.method = "vst",
                              nfeatures = 3000)
  
  processed_list[[i]] <- obj
  cat("Batch processing completed\n")
}

# 6. 合并所有批次
cat("\nMerging all batches...\n")
seurat_obj_main <- merge(x = processed_list[[1]], 
                         y = processed_list[-1],
                         add.cell.ids = paste0("Batch", seq_along(processed_list)))
seurat_obj_main <- subset(seurat_obj_main, subset = tissue_sampling_method !='brush'
                           & study!='Jain_Misharin_2021')
# 8. 清理中间对象节省内存
rm(merged_list, filled_list, processed_list)
gc()

# 获取所有unique的sample
samples <- unique(seurat_obj_main$sample)
samples <- samples[!is.na(samples)]  # 移除NA值

# 为每个sample创建一个列表来存储结果
sample_results <- list()

# 循环处理每个sample
for(sample_name in samples) {
  cat(sprintf("\nProcessing sample: %s\n", sample_name))
  
  # 提取当前sample的数据
  current_obj <- subset(seurat_obj_main, subset = sample == sample_name)
  
  # 获取细胞数量
  n_cells <- ncol(current_obj)
  if(n_cells < 50) {
    cat("Too few cells, skipping...\n")
    next
  }
  
  # 计算预期双细胞率
  doublet_rate <- min(0.075, n_cells / 1000 * 0.008)
  nExp_poi <- round(doublet_rate * n_cells)
  cat(sprintf("Expected doublets: %d (rate: %.3f)\n", nExp_poi, doublet_rate))
  
  # 数据准备
  current_obj <- NormalizeData(current_obj, verbose = TRUE)
  current_obj <- FindVariableFeatures(current_obj, verbose = TRUE)
  current_obj <- ScaleData(current_obj, verbose = TRUE)
  current_obj <- RunPCA(current_obj, npcs = 30, verbose = TRUE)
  
  # 参数优化和doublet检测
  tryCatch({
    sweep.res.list <- paramSweep(current_obj, PCs = 1:30, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    # 确保pK是数值
    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    if(is.na(pK)) pK <- 0.09
    cat(sprintf("Using pK: %.3f\n", pK))
    
    # 运行DoubletFinder
    current_obj <- doubletFinder(
      current_obj,
      PCs = 1:30,
      pN = 0.25,
      pK = pK,
      nExp = nExp_poi,
      sct = FALSE
    )
    
    # 获取classification列名
    DF.name <- colnames(current_obj@meta.data)[grep("DF.classifications", 
                                                    colnames(current_obj@meta.data))]
    
    # 保存doublet结果
    if(length(DF.name) > 0) {
      current_obj$doublet_status <- current_obj@meta.data[[DF.name[1]]]
      # 统计doublet数量
      doublet_count <- table(current_obj$doublet_status)
      cat("Doublet detection results:\n")
      print(doublet_count)
    } else {
      current_obj$doublet_status <- "Unknown"
      cat("No doublet classification found\n")
    }
    
  }, error = function(e) {
    warning(sprintf("DoubletFinder failed for sample %s: %s", sample_name, e$message))
    current_obj$doublet_status <- "Unknown"
  })
  
  # 保存结果
  sample_results[[sample_name]] <- current_obj$doublet_status
}

# 将所有结果整合回主对象
seurat_obj_main$doublet_status <- "Unknown"  # 初始化列
for(sample_name in names(sample_results)) {
  cells_in_sample <- WhichCells(seurat_obj_main, expression = sample == sample_name)
  seurat_obj_main$doublet_status[cells_in_sample] <- sample_results[[sample_name]]
}

# 打印总体统计
cat("\nOverall doublet detection results:\n")
print(table(seurat_obj_main$doublet_status))

cells_before <- ncol(seurat_obj_main)
seurat_obj_main <- subset(seurat_obj_main, subset = doublet_status != "Doublet")
cells_after <- ncol(seurat_obj_main)
cat(sprintf("Removed %d doublets (%.1f%%)\n", 
            cells_before - cells_after, 
            (cells_before - cells_after)/cells_before * 100))

# 定义可能的doublet组合
doublet_patterns <- list(
  "EPCAM-PECAM1-CD3D" = c("EPCAM", "PECAM1", "CD3D"),
  "PECAM1-CD79A" = c("PECAM1", "CD79A"),
  "CD79A-CD3D" = c("CD79A", "CD3D"),
  "EPCAM-KRT8-CD3D" = c("EPCAM", "KRT8", "CD3D"),
  "FOXJ1-COL1A2" = c("FOXJ1", "COL1A2"),
  "EPCAM-PECAM1" = c("EPCAM", "PECAM1"),
  "EPCAM-CD79A" = c("EPCAM", "CD79A")
)

#' Doublet识别和分析函数
#' @param seurat_obj Seurat对象
#' @param gene_combinations doublet基因组合列表
#' @param min_cells 最小细胞数阈值
#' @return 包含doublet分析结果的列表
analyze_doublets <- function(seurat_obj, gene_combinations, min_cells = 50) {
  # 初始化结果存储
  doublet_cells <- vector("list", length(gene_combinations) + 1)
  names(doublet_cells) <- c(names(gene_combinations), "CD14/CD68-EPCAM/CD79A")
  
  # 初始化doublet状态
  seurat_obj$doublet_status <- "Singlet" 
  
  # 按样本处理
  for(sample_name in unique(seurat_obj$sample)) {
    current_obj <- subset(seurat_obj, subset = sample == sample_name)
    n_cells <- ncol(current_obj)
    
    if(n_cells < min_cells) {
      cat(sprintf("Skipping sample %s: too few cells (%d)\n", sample_name, n_cells))
      next
    }
    
    # 处理常规doublet模式
    for(pattern_name in names(gene_combinations)) {
      genes <- gene_combinations[[pattern_name]]
      expr_mat <- GetAssayData(current_obj, layer = "data")[genes, ]
      threshold <- 0.4
      
      current_doublets <- colnames(expr_mat)[
        apply(expr_mat, 2, function(x) all(x > threshold))
      ]
      doublet_cells[[pattern_name]] <- c(
        doublet_cells[[pattern_name]], 
        current_doublets
      )
    }
    
    # 处理特殊组合
    group1_genes <- c("CD14", "CD68")
    group2_genes <- c("EPCAM", "CD79A")
    
    expr_mat_group1 <- GetAssayData(current_obj, layer = "data")[group1_genes, ]
    expr_mat_group2 <- GetAssayData(current_obj, layer = "data")[group2_genes, ]
    
    group1_expr <- apply(expr_mat_group1, 2, function(x) any(x > threshold))
    group2_expr <- apply(expr_mat_group2, 2, function(x) any(x > threshold))
    
    current_doublets <- colnames(current_obj)[group1_expr & group2_expr]
    doublet_cells[["CD14/CD68-EPCAM/CD79A"]] <- c(
      doublet_cells[["CD14/CD68-EPCAM/CD79A"]], 
      current_doublets
    )
  }
  
  # 标记doublets
  all_doublets <- unique(unlist(doublet_cells))
  seurat_obj$doublet_status[colnames(seurat_obj) %in% all_doublets] <- "Doublet"
  
  # 计算统计信息
  stats_df <- data.frame(
    Sample = unique(seurat_obj$sample),
    Total_Cells = tapply(seurat_obj$doublet_status, seurat_obj$sample, length),
    Doublets = tapply(seurat_obj$doublet_status == "Doublet", seurat_obj$sample, sum),
    stringsAsFactors = FALSE
  )
  
  stats_df$Singlets <- stats_df$Total_Cells - stats_df$Doublets
  stats_df$Doublet_Rate <- stats_df$Doublets / stats_df$Total_Cells * 100
  
  # 添加总计行
  total_row <- data.frame(
    Sample = "Total",
    Total_Cells = sum(stats_df$Total_Cells),
    Doublets = sum(stats_df$Doublets),
    Singlets = sum(stats_df$Singlets),
    Doublet_Rate = sum(stats_df$Doublets) / sum(stats_df$Total_Cells) * 100
  )
  stats_df <- rbind(stats_df, total_row)
  
  return(list(
    seurat_obj = seurat_obj,
    doublet_cells = doublet_cells,
    pattern_counts = sapply(doublet_cells, length),
    stats = stats_df
  ))
}

# 使用示例
results <- analyze_doublets(seurat_obj_main, doublet_patterns)
print(results$stats)

# 移除doublets并显示统计信息
cells_before <- ncol(results$seurat_obj)
seurat_obj_main <- subset(results$seurat_obj, subset = doublet_status != "Doublet")
cells_after <- ncol(seurat_obj_main)

# 打印移除结果
cat(sprintf("Removed %d doublets (%.1f%%)\n", 
            cells_before - cells_after, 
            (cells_before - cells_after)/cells_before * 100))

# 7. 输出最终结果信息
cat(sprintf("\nFinal merged object:\n"))
cat(sprintf("Total cells: %d\n", ncol(seurat_obj_main)))
cat(sprintf("Total features: %d\n", nrow(seurat_obj_main)))

all_marker_genes <- c(
  # T细胞标记基因
  "CD3E",     # CD3ε
  "CD3D",     # CD3δ
  "CD3G",     # CD3γ
  "CD5",      # CD5
  "TRAC",     # T细胞受体α链恒定区
  "TRBC1",    # T细胞受体β链恒定区1
  "TRBC2",    # T细胞受体β链恒定区2
  "TRGC1",    # T细胞受体γ链恒定区1
  "TRGC2",    # T细胞受体γ链恒定区2
  "TRDC",     # T细胞受体δ链恒定区
  "ITGAE",    # CD103，组织驻留T细胞
  "CD69",     # 组织驻留T细胞标记
  
  # 髓系细胞标记基因
  "ITGAX",    # CD11c，树突状细胞
  "ITGAM",    # CD11b，多种髓系细胞
  "CD14",     # 单核/巨噬细胞
  "CD68",     # 巨噬细胞
  "SIGLEC1",  # CD169，某些巨噬细胞亚群
  "MRC1",     # CD206，M2型巨噬细胞
  "MARCO",    # 肺泡巨噬细胞特异性标记
  "FUT4",     # CD15，中性粒细胞
  "CEACAM8",  # CD66b，中性粒细胞
  "MPO",      # 髓过氧化物酶，中性粒细胞
  
  # 功能相关基因
  "IFNG",     # 干扰素γ，T细胞激活
  "IL10",     # 白细胞介素10
  "LYZ",      # 溶菌酶，髓系细胞
  "CST3",     # 胱抑素C，髓系细胞
  "S100A8",   # 钙结合蛋白
  "S100A9",   # 钙结合蛋白
  "LAMP1",    # 溶酶体相关膜蛋白1
  
  # 组织特异性和其他相关基因
  "SPP1",     # 骨桥蛋白
  "EPCAM",    # 上皮细胞标记
  "PECAM1",   # 内皮细胞标记
  "PTPRC",    # CD45，所有白细胞
  "SELL"      # CD62L，L-选择素
)

# 检查基因是否存在于数据集中
genes_present <- all_marker_genes[all_marker_genes %in% rownames(seurat_obj_main)]
genes_missing <- all_marker_genes[!all_marker_genes %in% rownames(seurat_obj_main)]

# 创建输出文件夹（如果不存在）
dir.create("marker_gene_plots", showWarnings = FALSE)

# 打印存在和缺失的基因
cat("存在的基因数量:", length(genes_present), "\n")
cat("缺失的基因数量:", length(genes_missing), "\n")

# 保存基因存在信息到文本文件
sink("marker_gene_plots/gene_status.txt")
cat("存在的基因数量:", length(genes_present), "\n")
cat("存在的基因:", paste(genes_present, collapse=", "), "\n\n")
cat("缺失的基因数量:", length(genes_missing), "\n")
if(length(genes_missing) > 0) {
  cat("缺失的基因:", paste(genes_missing, collapse=", "), "\n")
}
sink()

# 确保基因在scale.data中可用
# 由于很多基因不在scale.data中，我们需要先对这些基因进行标准化
genes_to_scale <- genes_present[!genes_present %in% rownames(seurat_obj_main[["RNA"]]@scale.data)]
if(length(genes_to_scale) > 0) {
  seurat_obj_main <- ScaleData(seurat_obj_main, features = genes_to_scale)
  cat("已对缺失的", length(genes_to_scale), "个基因进行标准化\n")
}

# 创建多页PDF文件
if(length(genes_present) > 0) {
  # 使用英文标题避免中文编码问题
  # 1. 热图 - 所有基因在聚类中的表达情况
  pdf("marker_gene_plots/heatmap_all_markers.pdf", width = 12, height = 20)
  print(
    DoHeatmap(seurat_obj_main, features = genes_present, group.by = "seurat_clusters") +
      ggtitle("Expression of Myeloid and T-cell Markers in Clusters")
  )
  dev.off()
  
  # 2. 关键基因组合的特征图 (一页一个基因)
  key_genes <- c("CD3E", "CD3D", "CD5", "ITGAE", "TRAC", # T细胞标记
                 "ITGAM", "ITGAX", "CD14", "CD68", "MPO") # 髓系标记
  key_genes_present <- key_genes[key_genes %in% genes_present]
  
  pdf("marker_gene_plots/feature_plots_key_genes.pdf", width = 10, height = 8)
  for(gene in key_genes_present) {
    # 使用英文标题，减少样本量
    p <- FeaturePlot(seurat_obj_main, 
                     features = gene, 
                     reduction = "umap",
                     raster = FALSE,
                     pt.size = 0.5,
                     cells = sample(colnames(seurat_obj_main), min(20000, ncol(seurat_obj_main)))) + 
      ggtitle(paste("Expression of", gene)) +
      theme(plot.title = element_text(size = 16, face = "bold"))
    print(p)
  }
  dev.off()
  
  # 3. 使用单独的FeaturePlot
  t_cell_markers <- c("CD3E", "CD3D", "CD5", "ITGAE", "TRAC")
  myeloid_markers <- c("ITGAM", "ITGAX", "CD14", "CD68", "MPO")
  
  t_present <- t_cell_markers[t_cell_markers %in% genes_present]
  myeloid_present <- myeloid_markers[myeloid_markers %in% genes_present]
  
  if(length(t_present) > 0 && length(myeloid_present) > 0) {
    # 尝试使用基本的R绘图功能而不是复杂的patchwork布局
    pdf("marker_gene_plots/paired_feature_plots.pdf", width = 12, height = 8)
    
    for(t_gene in t_present) {
      for(m_gene in myeloid_present) {
        # 绘制T细胞标记
        par(mfrow=c(1,2)) # 设置为1行2列的布局
        
        # 获取表达数据
        t_expr <- GetAssayData(seurat_obj_main, layer = "data")[t_gene, ]
        m_expr <- GetAssayData(seurat_obj_main, layer = "data")[m_gene, ]
        umap_coords <- Embeddings(seurat_obj_main, reduction = "umap")
        
        # 随机抽样20000个细胞
        set.seed(42)
        sample_idx <- sample(1:ncol(seurat_obj_main), min(20000, ncol(seurat_obj_main)))
        
        # T细胞标记图
        plot(umap_coords[sample_idx,1], umap_coords[sample_idx,2], 
             pch=16, cex=0.5, col=colorRampPalette(c("grey", "blue"))(10)[cut(t_expr[sample_idx], breaks=10)],
             main=paste("Expression of", t_gene), xlab="UMAP1", ylab="UMAP2")
        
        # 髓系标记图
        plot(umap_coords[sample_idx,1], umap_coords[sample_idx,2], 
             pch=16, cex=0.5, col=colorRampPalette(c("grey", "red"))(10)[cut(m_expr[sample_idx], breaks=10)],
             main=paste("Expression of", m_gene), xlab="UMAP1", ylab="UMAP2")
        
        # 重置绘图参数
        par(mfrow=c(1,1))
        
        # 添加散点图比较两个基因的表达
        plot(t_expr[sample_idx], m_expr[sample_idx], 
             pch=16, cex=0.5, col=adjustcolor("black", alpha.f=0.5),
             main=paste(t_gene, "vs", m_gene, "Expression"),
             xlab=t_gene, ylab=m_gene)
        abline(h=0, v=0, lty=2, col="grey")
      }
    }
    dev.off()
  }
  
  # 4. 小提琴图展示不同聚类中的标记基因表达
  pdf("marker_gene_plots/violin_plots.pdf", width = 12, height = 8)
  
  # 为避免一次绘制太多小提琴图，分批绘制
  t_batches <- split(t_present, ceiling(seq_along(t_present)/4))
  m_batches <- split(myeloid_present, ceiling(seq_along(myeloid_present)/4))
  
  for(t_batch in t_batches) {
    if(length(t_batch) > 0) {
      p <- VlnPlot(seurat_obj_main, features = t_batch, group.by = "seurat_clusters", pt.size = 0, ncol = 2) +
        plot_annotation(title = "T-cell Markers Expression")
      print(p)
    }
  }
  
  for(m_batch in m_batches) {
    if(length(m_batch) > 0) {
      p <- VlnPlot(seurat_obj_main, features = m_batch, group.by = "seurat_clusters", pt.size = 0, ncol = 2) +
        plot_annotation(title = "Myeloid Markers Expression")
      print(p)
    }
  }
  dev.off()
  
  # 5. 散点图比较潜在双表型标记对的表达 - 使用基本R绘图
  pdf("marker_gene_plots/scatter_plots.pdf", width = 10, height = 10)
  # 选择主要的几对T细胞-髓系标记组合
  key_pairs <- list(
    c("CD3E", "ITGAM"),  # CD3-CD11b
    c("CD3E", "ITGAX"),  # CD3-CD11c
    c("CD5", "CD14"),    # CD5-CD14
    c("CD3E", "CD68"),   # CD3-CD68
    c("ITGAE", "CD68")   # CD103-CD68
  )
  
  for(pair in key_pairs) {
    if(all(pair %in% genes_present)) {
      # 提取基因表达数据
      gene1 <- pair[1]
      gene2 <- pair[2]
      
      # 获取表达数据和聚类信息
      expr_data <- GetAssayData(seurat_obj_main, layer = "data")[c(gene1, gene2), ]
      clusters <- seurat_obj_main$seurat_clusters
      
      # 随机抽样
      set.seed(42)
      sample_idx <- sample(1:ncol(seurat_obj_main), min(10000, ncol(seurat_obj_main)))
      
      # 使用不同颜色表示不同聚类
      cluster_colors <- rainbow(length(unique(clusters)))
      names(cluster_colors) <- unique(clusters)
      
      # 绘制散点图
      plot(expr_data[1, sample_idx], expr_data[2, sample_idx], 
           pch=16, cex=0.5, 
           col=cluster_colors[clusters[sample_idx]],
           main=paste(gene1, "vs", gene2, "Expression"),
           xlab=gene1, ylab=gene2)
      
      # 添加图例
      legend("topright", legend=names(cluster_colors), 
             col=cluster_colors, pch=16, cex=0.8, title="Clusters")
      
      # 添加参考线
      abline(h=0, v=0, lty=2, col="grey")
    }
  }
  dev.off()
}

cat("所有图表已保存到 'marker_gene_plots' 文件夹\n")

# 设置T细胞标记和髓系标记
t_cell_markers <- c("CD3E", "CD5", "CD3D")
myeloid_markers <- c("ITGAM", "ITGAX", "CD14", "CD68")

# 创建输出PDF
pdf("marker_gene_plots/t_myeloid_blend_plots.pdf", width = 24, height = 10)

# 为每个T细胞标记与每个髓系标记创建blend FeaturePlot
for (t_marker in t_cell_markers) {
  if (t_marker %in% rownames(seurat_obj_main)) {
    for (m_marker in myeloid_markers) {
      if (m_marker %in% rownames(seurat_obj_main)) {
        # 检查两个标记是否都存在
        cat(paste("绘制", t_marker, "vs", m_marker, "blend plot\n"))
        
        # 使用blend参数创建共表达图
        p <- FeaturePlot(seurat_obj_main, 
                         features = c(t_marker, m_marker),
                         blend = TRUE)
        
        print(p + plot_annotation(title = paste(t_marker, "vs", m_marker, "Co-expression")))
      }
    }
  }
}

dev.off()

cat("所有blend共表达图表已保存到marker_gene_plots/t_myeloid_blend_plots.pdf文件\n")

# 设置T细胞标记和髓系标记
t_cell_markers <- c("CD3E", "CD5", "CD3D")
myeloid_markers <- c("ITGAM", "ITGAX", "CD14", "CD68")

# 定义函数：识别并提取双表型细胞
identify_dual_phenotype_cells <- function(seurat_obj, 
                                          t_marker, 
                                          m_marker, 
                                          threshold = 0.3,
                                          visualize = TRUE) {
  # 检查标记是否存在
  if(!(t_marker %in% rownames(seurat_obj) && m_marker %in% rownames(seurat_obj))) {
    cat(paste("标记", t_marker, "或", m_marker, "不在数据集中\n"))
    return(NULL)
  }
  
  # 获取表达数据
  expr_data <- GetAssayData(seurat_obj, layer = "data")[c(t_marker, m_marker), ]
  
  # 识别双表型细胞
  dual_cells <- colnames(expr_data)[expr_data[t_marker, ] > threshold & 
                                      expr_data[m_marker, ] > threshold]
  
  # 打印结果
  cat(paste("找到", length(dual_cells), "个同时表达", t_marker, 
            "和", m_marker, "的细胞\n"))
  
  # 如果需要可视化
  if(visualize && length(dual_cells) > 0) {
    # 在UMAP上展示这些细胞
    p1 <- DimPlot(seurat_obj, 
                  cells.highlight = dual_cells,
                  cols.highlight = "purple",
                  sizes.highlight = 1,
                  reduction = "umap") +
      ggtitle(paste(t_marker, "+", m_marker, "+细胞"))
    
    # 在表达散点图上展示
    all_cells <- sample(colnames(seurat_obj), min(10000, ncol(seurat_obj)))
    expr_data_sample <- FetchData(seurat_obj, 
                                  vars = c(t_marker, m_marker, "seurat_clusters"), 
                                  cells = all_cells)
    
    dual_cells_in_sample <- intersect(dual_cells, all_cells)
    expr_data_sample$cell_type <- "Other"
    expr_data_sample$cell_type[rownames(expr_data_sample) %in% dual_cells_in_sample] <- "Dual Phenotype"
    
    p2 <- ggplot(expr_data_sample, aes_string(x = t_marker, y = m_marker)) +
      geom_point(aes(color = cell_type), alpha = 0.7, size = 1) +
      scale_color_manual(values = c("Dual Phenotype" = "purple", "Other" = "grey")) +
      theme_classic() +
      labs(title = paste(t_marker, "vs", m_marker, "表达")) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      geom_vline(xintercept = threshold, linetype = "dashed", color = "red")
    
    # 打印图表
    print(p1)
    print(p2)
  }
  
  # 返回双表型细胞的列表
  return(dual_cells)
}

# 处理所有标记组合并保存结果
pdf("marker_gene_plots/dual_phenotype_cells.pdf", width = 16, height = 8)

# 存储所有结果
all_dual_cells <- list()

# 处理每一对组合
for (t_marker in t_cell_markers) {
  for (m_marker in myeloid_markers) {
    # 定义组合名称
    combo_name <- paste(t_marker, m_marker, sep = "_")
    
    # 识别双表型细胞
    dual_cells <- identify_dual_phenotype_cells(seurat_obj_main, 
                                                t_marker, 
                                                m_marker, 
                                                threshold = 0.3,  # 可以调整阈值
                                                visualize = TRUE)
    
    # 保存结果
    all_dual_cells[[combo_name]] <- dual_cells
  }
}
dev.off()

# 保存所有双表型细胞的统计信息
sink("marker_gene_plots/dual_phenotype_stats.txt")
cat("双表型细胞统计信息:\n\n")

for (combo_name in names(all_dual_cells)) {
  cells <- all_dual_cells[[combo_name]]
  if (!is.null(cells)) {
    cat(paste(combo_name, "组合:", length(cells), "个细胞\n"))
    
    # 如果存在足够的双表型细胞，分析它们的聚类分布
    if (length(cells) >= 10) {
      cluster_counts <- table(seurat_obj_main$seurat_clusters[cells])
      cat("  聚类分布:\n")
      for (cluster in names(cluster_counts)) {
        cat(paste("    聚类", cluster, ":", cluster_counts[cluster], "个细胞\n"))
      }
    }
    cat("\n")
  }
}
sink()

# 创建一个新的元数据列标记双表型细胞
seurat_obj_main$dual_phenotype <- "None"

# 对每种组合的双表型细胞进行标记
for (combo_name in names(all_dual_cells)) {
  cells <- all_dual_cells[[combo_name]]
  if (!is.null(cells) && length(cells) > 0) {
    seurat_obj_main$dual_phenotype[colnames(seurat_obj_main) %in% cells] <- combo_name
  }
}

# 保存标记了双表型细胞的Seurat对象(可选)
# saveRDS(seurat_obj_main, "seurat_obj_with_dual_phenotype.rds")

# 创建一个合并的UMAP图，显示所有双表型细胞
pdf("marker_gene_plots/all_dual_phenotype_cells.pdf", width = 12, height = 10)

# 所有双表型细胞的合并列表
all_dual_cells_combined <- unique(unlist(all_dual_cells))

# 如果有双表型细胞，创建UMAP图
if (length(all_dual_cells_combined) > 0) {
  p <- DimPlot(seurat_obj_main,
               group.by = "dual_phenotype",
               cols = c("None" = "grey", 
                        setNames(rainbow(length(names(all_dual_cells))), 
                                 names(all_dual_cells)))) +
    ggtitle("所有双表型细胞在UMAP上的分布")
  print(p)
}

dev.off()
cat("原始细胞数:", ncol(seurat_obj_main), "\n")
seurat_obj_main <- subset(seurat_obj_main, subset = dual_phenotype == "None")
cat("剩余细胞数:", ncol(seurat_obj_main), "\n")
cat("双表型细胞分析完成，结果已保存到marker_gene_plots文件夹\n")

saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")




regression_outliers_by_sample <- function(seurat_obj, outliers_threshold = 0.999) {
  # 创建输出目录 [^1]
  dir.create("qaqc_plots", recursive = TRUE, showWarnings = FALSE)
  
  # 存储所有要保留的细胞
  all_cells_keep <- vector()
  
  # 按样本处理
  for(sample_name in unique(seurat_obj$sample)) {
    # 获取当前样本数据
    current_obj <- subset(seurat_obj, subset = sample == sample_name)
    
    # 获取metrics数据
    nCount <- current_obj$nCount_RNA
    nFeature <- current_obj$nFeature_RNA
    
    # 对数转换
    log_counts <- log10(nCount + 1)
    log_features <- log10(nFeature + 1)
    
    # 线性回归
    fit <- lm(log_features ~ log_counts)
    
    # 计算预测区间
    pred <- predict(fit, 
                    data.frame(log_counts = log_counts), 
                    interval = "prediction", 
                    level = outliers_threshold)
    
    # 标记异常值
    outliers <- log_features < pred[,"lwr"] | log_features > pred[,"upr"]
    
    # 可视化
    pdf(file.path("qaqc_plots", 
                  paste0("regression_outliers_", sample_name, ".pdf")),
        width = 10, height = 8)
    
    # 绘制散点图
    plot(log_counts, log_features,
         pch = 16, cex = 0.6,
         col = ifelse(outliers, "red", "grey50"),
         xlab = "log10(nCount_RNA)",
         ylab = "log10(nFeature_RNA)",
         main = paste("Regression Outliers Detection -", sample_name))
    
    # 添加拟合线和预测区间
    lines(sort(log_counts), pred[order(log_counts),"lwr"], col = "blue", lty = 2)
    lines(sort(log_counts), pred[order(log_counts),"upr"], col = "blue", lty = 2)
    lines(sort(log_counts), pred[order(log_counts),"fit"], col = "blue", lwd = 2)
    
    legend("topleft",
           c("Cells", "Outliers", "Fit", "Prediction Interval"),
           col = c("grey50", "red", "blue", "blue"),
           pch = c(16, 16, NA, NA),
           lty = c(NA, NA, 1, 2),
           lwd = c(NA, NA, 2, 1))
    
    dev.off()
    
    # 输出每个样本的结果
    cat(sprintf("Sample %s: Detected %d outliers (%.1f%%)\n", 
                sample_name,
                sum(outliers), 
                mean(outliers) * 100))
    
    # 收集要保留的细胞
    cells_keep <- colnames(current_obj)[!outliers]
    all_cells_keep <- c(all_cells_keep, cells_keep)
  }
  
  # 返回所有要保留的细胞
  return(all_cells_keep)
}

# 使用函数
cells_to_keep <- regression_outliers_by_sample(seurat_obj_main)
seurat_obj_main <- subset(seurat_obj_main, cells = cells_to_keep)

seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")

# 2. 先进行标准化和特征选择
print("Normalizing data and finding variable features...")
seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
# table(seurat_obj_main@meta.data$sample)
# seurat_obj_main <- subset(seurat_obj_main, subset =  sample!='SC144')
# clusters_to_remove <- c("93", "82")  # 根据需要添加
# seurat_obj_main <- subset(seurat_obj_main, 
#                               subset = RNA_snn_res.4 != "93" & RNA_snn_res.4 != "82")
seurat_obj_main <- NormalizeData(seurat_obj_main)
seurat_obj_main <- FindVariableFeatures(seurat_obj_main, selection.method = "vst", nfeatures = 4000)
seurat_obj_main <- ScaleData(seurat_obj_main)
seurat_obj_main <- RunPCA(seurat_obj_main,npcs = 40)
# # 
# seurat_obj_main <- RenameGenesSeurat(
#   obj = seurat_obj_main,
#   newnames_file = "E:/R/sources/10x/featrues/3/merged_all_features.csv",
#   remove_duplicates = TRUE
# )
# saveRDS(seurat_obj_main,"final_analyzed_seurat_obj_main.rds")
seurat_obj_main <- RunHarmony(
  object = seurat_obj_main,           
  group.by.vars = c("sample",'tissue_sampling_method'),   
  theta = c(1,1),      # Higher theta for more diverse clustering               
  lambda = c(7,7),   # Higher lambda to reduce overcorrection                
  sigma = 0.01,           # Lower sigma for tighter clusters
  nclust = 60,            # Increased number of clusters
  reduction.use = "pca",
  max_iter = 9, 
  early_stop = FALSE,
  dims = 1:40           # More iterations for better convergence
)

print("Calculating integration quality metrics...")

# 获取降维结果和元数据
embeddings <- Embeddings(seurat_obj_main, "harmony")
metadata <- seurat_obj_main@meta.data

# 过滤掉ann_level_2为NA的细胞
valid_cells <- !is.na(metadata$ann_level_2)
embeddings <- embeddings[valid_cells, ]
metadata <- metadata[valid_cells, ]

# 检查并处理无效值
embeddings_valid <- is.finite(rowSums(embeddings))
if(sum(!embeddings_valid) > 0) {
  cat(sprintf("Removed %d cells with NA/NaN/Inf values\n", sum(!embeddings_valid)))
  embeddings <- embeddings[embeddings_valid, ]
  metadata <- metadata[embeddings_valid, ]
}

# 设置采样参数
n_cells_sample <- min(10000, nrow(embeddings))
n_iterations <- 10
k <- 30

# 多次抽样计算LISI scores
ilisi_scores <- c()
clisi_scores <- c()

for(i in 1:n_iterations) {
  set.seed(i)
  cells_idx <- sample(nrow(embeddings), n_cells_sample)
  
  # 抽取对应的数据
  sampled_embeddings <- embeddings[cells_idx, ]
  sampled_metadata <- metadata[cells_idx, ]
  
  # 尝试计算LISI scores
  tryCatch({
    lisi_res <- lisi::compute_lisi(sampled_embeddings, 
                                   sampled_metadata, 
                                   c("sample", "ann_level_2"), 
                                   k)
    
    # 存储结果
    ilisi_scores[i] <- mean(lisi_res[, "sample"])
    clisi_scores[i] <- mean(lisi_res[, "ann_level_2"])
    
    cat(sprintf("Iteration %d completed: iLISI = %.3f, cLISI = %.3f\n", 
                i, ilisi_scores[i], clisi_scores[i]))
  }, error = function(e) {
    cat(sprintf("Error in iteration %d: %s\n", i, e$message))
  })
}

# 计算并输出结果
if(length(ilisi_scores) > 0) {
  cat("\nIntegration Quality Metrics:")
  cat("\n--------------------------")
  cat(sprintf("\nNumber of cells used: %d (after removing NA annotations)", nrow(embeddings)))
  cat(sprintf("\niLISI Score: %.3f ± %.3f", mean(ilisi_scores), sd(ilisi_scores)))
  cat(sprintf("\ncLISI Score: %.3f ± %.3f", mean(clisi_scores), sd(clisi_scores)))
  cat("\n\nInterpretation:")
  cat("\n- iLISI: Higher values (closer to N_batches) indicate better batch mixing")
  cat("\n- cLISI: Lower values (closer to 1) indicate better cell type separation")
  cat("\n--------------------------\n")
} else {
  cat("\nWarning: Could not compute LISI scores due to errors in all iterations\n")
}

# # 创建输出目录（如果不存在）
# if(!dir.exists("qc_plots")) {
#   dir.create("qc_plots")
# }
# 
# # 保存所有knee plots到一个PDF文件
# pdf("qc_plots/knee_plots_analysis.pdf", width = 12, height = 10)
# 
# # 1. UMI knee plot
# counts <- GetAssayData(seurat_obj_main, layer = "counts")
# total_umi <- Matrix::colSums(counts)
# umi_rank <- rank(-total_umi)
# 
# print(ggplot(data.frame(rank = umi_rank, UMI = total_umi), 
#              aes(x = rank, y = UMI)) +
#         geom_line(color = "blue") +
#         geom_point(size = 0.5, alpha = 0.5, color = "blue") +
#         scale_x_log10() +
#         scale_y_log10() +
#         theme_bw() +
#         labs(title = "Post-integration UMI Knee Plot",
#              x = "Cell rank",
#              y = "Total UMIs"))
# 
# # 2. 基因数量的knee plot
# genes_per_cell <- Matrix::colSums(counts > 0)
# gene_rank <- rank(-genes_per_cell)
# 
# print(ggplot(data.frame(rank = gene_rank, Genes = genes_per_cell), 
#              aes(x = rank, y = Genes)) +
#         geom_line(color = "blue") +
#         geom_point(size = 0.5, alpha = 0.5, color = "blue") +
#         scale_x_log10() +
#         scale_y_log10() +
#         theme_bw() +
#         labs(title = "Post-integration Gene Detection Knee Plot",
#              x = "Cell rank",
#              y = "Number of genes detected"))
# 
# # 3. PC肘部图
# pct_var <- seurat_obj_main[["pca"]]@stdev / sum(seurat_obj_main[["pca"]]@stdev) * 100
# cumsum_var <- cumsum(pct_var)
# 
# elbow_data <- data.frame(
#   PC = 1:length(pct_var),
#   variance = pct_var,
#   cumulative = cumsum_var
# )
# 
# print(ggplot(elbow_data, aes(x = PC)) +
#         geom_line(aes(y = variance), color = "blue") +
#         geom_point(aes(y = variance), color = "blue") +
#         geom_line(aes(y = cumulative), color = "red") +
#         geom_point(aes(y = cumulative), color = "red") +
#         theme_bw() +
#         labs(title = "PCA Elbow Plot",
#              x = "Principal Component",
#              y = "Percentage of Variance Explained") +
#         scale_y_continuous(
#           name = "Individual Variance Explained (%)",
#           sec.axis = sec_axis(~., name = "Cumulative Variance Explained (%)")
#         ))
# 
# # 4. Harmony scores knee plot - 使用采样
# set.seed(42)
# n_sample <- min(10000, ncol(seurat_obj_main))
# sample_idx <- sample(1:ncol(seurat_obj_main), n_sample)
# 
# harmony_scores <- seurat_obj_main[["harmony"]]@cell.embeddings[sample_idx, 1:20]
# harmony_dist <- dist(harmony_scores)
# harmony_scores_ordered <- sort(colMeans(as.matrix(harmony_dist)), decreasing = TRUE)
# 
# print(ggplot(data.frame(rank = 1:length(harmony_scores_ordered),
#                         score = harmony_scores_ordered), 
#              aes(x = rank, y = score)) +
#         geom_line(color = "blue") +
#         geom_point(size = 0.5, alpha = 0.5, color = "blue") +
#         theme_bw() +
#         labs(title = "Harmony Scores Distribution (10k cells sampled)",
#              x = "Cell rank",
#              y = "Average Harmony distance"))
# 
# dev.off()
# 
# # 保存统计信息到CSV
# stats_df <- data.frame(
#   Metric = c("Total_cells", "Median_UMI", "Median_genes", 
#              "PC_90percent_var", "Harmony_mean_dist"),
#   Value = c(
#     length(total_umi),
#     median(total_umi),
#     median(genes_per_cell),
#     which(cumsum_var > 90)[1],
#     mean(harmony_scores_ordered)
#   )
# )
# 
# write.csv(stats_df, "qc_plots/knee_plot_statistics.csv", row.names = FALSE)

# 直接运行UMAP - 优化参数设置和注释
print("Running UMAP and clustering...")
tryCatch({
  # 运行UMAP降维
  seurat_obj_main <- RunUMAP(seurat_obj_main, 
                           reduction = "harmony", 
                           dims = 1:30, 
                           n.neighbors = 20,
                           min.dist = 0.4,
                           n.epochs = 1600,
                           repulsion.strength = 1.1,
                           metric = "correlation",
                           verbose = FALSE)
  
  # 寻找邻居关系用于聚类
  seurat_obj_main <- FindNeighbors(seurat_obj_main, 
                                 reduction = "harmony", 
                                 dims = 1:30,
                                 k.param = 20,
                                 verbose = FALSE) 

  # 在多个分辨率下进行聚类
  seurat_obj_main <- FindClusters(seurat_obj_main,
                                resolution = 3,
                                random.seed = 42,
                                verbose = TRUE)
  
  print(paste("Completed clustering with", length(unique(seurat_obj_main$seurat_clusters)), "clusters"))
}, error = function(e) {
  stop(paste("Error in UMAP/clustering:", e$message))
})

# 创建可视化输出目录
vis_dir <- "visualization_results"
if(!dir.exists(vis_dir)) dir.create(vis_dir)

# 创建统一的可视化函数，提高代码复用性
generate_dimplot <- function(object, group.by, filename, title = NULL, 
                           split.by = NULL, cells.highlight = NULL, 
                           pt.size = 0.5, label = TRUE, label.size = 4,
                           width = 10, height = 8, raster = TRUE,
                           cols = NULL, return_plot = FALSE) {
  
  if(is.null(title)) title <- paste("UMAP by", group.by)
  
  p <- DimPlot(object, 
              reduction = "umap",
              group.by = group.by,
              split.by = split.by,
              cells.highlight = cells.highlight,
              cols = cols,
              pt.size = pt.size,
              label = label,
              label.size = label.size,
              raster = raster) + 
        ggtitle(title) + 
        theme(legend.text = element_text(size = 10))
  
  # 保存图形
  if(!is.null(filename)) {
    pdf(file.path(vis_dir, filename), width = width, height = height)
    print(p)
    dev.off()
  }
  
  if(return_plot) return(p)
}

# 绘制主UMAP聚类图
print("Generating visualization plots...")
Idents(seurat_obj_main) <- seurat_obj_main$seurat_clusters
generate_dimplot(seurat_obj_main, "seurat_clusters", "umap_clusters.pdf", 
                "UMAP Clustering")

# 主要元数据可视化 - 使用优化的函数
print("Generating metadata visualizations...")
pdf(file.path(vis_dir, "umap_anno_integration.pdf"), width = 15, height = 10)
plots <- lapply(c("study", "sample", "ann_level_2", "tissue"), function(meta_col) {
  if(meta_col %in% colnames(seurat_obj_main@meta.data)) {
    return(generate_dimplot(seurat_obj_main, meta_col, NULL, 
                         paste("UMAP by", meta_col), return_plot = TRUE))
  }
  return(NULL)
})
plots <- plots[!sapply(plots, is.null)]
for(p in plots) print(p)
dev.off()

# 创建分离的UMAP图 - 组织、样本和研究的高级可视化
metadata_categories <- list(
  tissue = list(col = "tissue", title = "Tissue"),
  sample = list(col = "sample", title = "Sample"),
  study = list(col = "study", title = "Study")
)

for(category_name in names(metadata_categories)) {
  category <- metadata_categories[[category_name]]
  col_name <- category$col
  
  if(col_name %in% colnames(seurat_obj_main@meta.data)) {
    unique_values <- unique(seurat_obj_main@meta.data[[col_name]])
    unique_values <- unique_values[!is.na(unique_values)]
    
    if(length(unique_values) > 0) {
      print(paste("Generating", category$title, "split visualizations..."))
      
      # 生成分离图形
      pdf(file.path(vis_dir, paste0(category_name, "_umap_split.pdf")), 
          width = 12, height = 12)
      Idents(seurat_obj_main) <- col_name
      
      # 如果值太多，添加分页
      values_per_page <- 12
      value_pages <- split(unique_values, ceiling(seq_along(unique_values)/values_per_page))
      
      for(page_values in value_pages) {
        for(value in page_values) {
          tryCatch({
            p <- DimPlot(seurat_obj_main, 
                        reduction = "umap",
                        cells.highlight = CellsByIdentities(object = seurat_obj_main, 
                                                          idents = value),
                        cols = "grey",
                        pt.size = 0.5,
                        label = FALSE) +
                  ggtitle(value) +
                  theme(legend.text = element_text(size = 12))
            print(p)
          }, error = function(e) {
            message(paste("Error plotting", category$title, value, ":", e$message))
          })
        }
      }
      dev.off()
    }
  }
}

# 处理ann_level_2的特殊情况
print("Generating cell type annotations visualizations...")
if("ann_level_2" %in% colnames(seurat_obj_main@meta.data)) {
  # 创建不含NA值的版本
  seurat_obj_main$ann_level_2_no_na <- seurat_obj_main$ann_level_2
  seurat_obj_main$ann_level_2_no_na[is.na(seurat_obj_main$ann_level_2_no_na)] <- "Unknown"
  Idents(seurat_obj_main) <- "ann_level_2_no_na"
  
  # 获取不含NA的唯一细胞类型
  ann_level_2 <- unique(seurat_obj_main$ann_level_2)
  ann_level_2 <- ann_level_2[!is.na(ann_level_2)]
  
  # 生成可视化
  if(length(ann_level_2) > 0) {
    pdf(file.path(vis_dir, "ann_level_2_umap_split.pdf"), width = 12, height = 12)
    
    # 批量处理，每页最多9个类型
    cell_type_pages <- split(ann_level_2, ceiling(seq_along(ann_level_2)/9))
    
    for(page_types in cell_type_pages) {
      for(anno in page_types) {
        tryCatch({
          p <- DimPlot(seurat_obj_main, 
                      reduction = "umap",
                      cells.highlight = CellsByIdentities(seurat_obj_main, idents = anno),
                      cols = "grey",
                      pt.size = 0.5,
                      raster = TRUE,
                      label = FALSE) +
                ggtitle(anno) +
                theme(legend.text = element_text(size = 12))
          print(p)
        }, error = function(e) {
          message(paste("Error plotting annotation", anno, ":", e$message))
        })
      }
    }
    dev.off()
  }
}

# 标记基因可视化 - 优化特征图生成
print("Generating marker gene feature plots...")

# 首先检查并移除标记基因中的重复项
Markers <- unique(Markers)
message(paste("Processing", length(Markers), "unique marker genes"))

# 函数：分批次创建特征图以优化性能
generate_feature_plots <- function(seurat_obj, markers, batch_size = 6, 
                                 output_file = "umap_FeaturePlot.pdf",
                                 width = 8, height = 8) {
  
  # 只保留对象中存在的基因
  valid_markers <- markers[markers %in% rownames(seurat_obj)]
  if(length(valid_markers) == 0) {
    message("No valid markers found in the dataset.")
    return(NULL)
  }
  
  message(paste("Creating feature plots for", length(valid_markers), "genes"))
  
  # 分批次
  marker_batches <- split(valid_markers, ceiling(seq_along(valid_markers)/batch_size))
  
  # 创建PDF
  pdf(file.path(vis_dir, output_file), width = width, height = height)
  
  # 进度跟踪
  total_batches <- length(marker_batches)
  pb <- txtProgressBar(min = 0, max = total_batches, style = 3)
  
  # 处理每个批次
  for(i in seq_along(marker_batches)) {
    batch_markers <- marker_batches[[i]]
    tryCatch({
      # 为每个基因单独创建特征图
      for(marker in batch_markers) {
        p <- FeaturePlot(seurat_obj, 
                        features = marker,
                        raster = TRUE,
                        pt.size = 0.5) +
             ggtitle(marker)
        print(p)
      }
    }, error = function(e) {
      message(paste("Error plotting batch", i, ":", e$message))
    })
    setTxtProgressBar(pb, i)
  }
  close(pb)
  dev.off()
  
  return(valid_markers)
}

# 生成标记基因特征图
plotted_markers <- generate_feature_plots(seurat_obj_main, Markers)

# 优化的点图生成
print("Generating DotPlot for marker genes...")
tryCatch({
  # 优化宽度计算
  dot_plot_width <- max(14, length(plotted_markers) * 0.25)
  pdf(file.path(vis_dir, "markers_dotplot.pdf"), width = dot_plot_width, height = 18)
  
  dot_plot <- DotPlot(seurat_obj_main, 
                     features = plotted_markers, 
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "red"),
                     dot.scale = 8) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  # 保存数据
  source_data <- dot_plot$data
  write.csv(source_data, file.path(vis_dir, "dotplot_data.csv"), row.names = FALSE)
  print(dot_plot)
  dev.off()
}, error = function(e) {
  message(paste("Error generating dot plot:", e$message))
})

# 识别并保存cluster特异性marker基因
print("Finding and analyzing cluster-specific marker genes...")
markers_file <- file.path(vis_dir, "cluster_markers.csv")

# 只有当标记文件不存在时才寻找标记
if(!file.exists(markers_file)) {
  tryCatch({
    Idents(seurat_obj_main) <- seurat_obj_main$seurat_clusters
    markers <- FindAllMarkers(seurat_obj_main, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25,
                             test.use = "wilcox",
                             verbose = TRUE)
    
    write.csv(markers, markers_file)
    message(paste("Saved marker genes to", markers_file))
  }, error = function(e) {
    message(paste("Error finding cluster markers:", e$message))
  })
} else {
  message(paste("Using existing marker gene file:", markers_file))
  markers <- read.csv(markers_file)
}

# 优化的热图生成
print("Generating marker heatmap...")
if(exists("markers") && nrow(markers) > 0) {
  tryCatch({
    # 为每个聚类选择top markers
    top10_markers <- markers %>% 
      group_by(cluster) %>% 
      top_n(10, wt = avg_log2FC) %>%
      ungroup()
    
    # 提取唯一基因并缩放
    marker_genes <- unique(top10_markers$gene)
    if(length(marker_genes) > 0) {
      message(paste("Scaling data for", length(marker_genes), "top marker genes"))
      seurat_obj_main <- ScaleData(seurat_obj_main, features = marker_genes, verbose = FALSE)
      
      # 将marker基因分成较小的批次以优化内存使用
      batch_size <- 25
      marker_batches <- split(marker_genes, ceiling(seq_along(marker_genes)/batch_size))
      
      # 创建热图
      pdf(file.path(vis_dir, "marker_heatmap_batches.pdf"), width = 32, height = 18)
      
      for(i in seq_along(marker_batches)) {
        message(paste("Processing heatmap batch", i, "of", length(marker_batches)))
        batch_genes <- marker_batches[[i]]
        
        tryCatch({
          p <- DoHeatmap(seurat_obj_main, 
                        features = batch_genes,
                        size = 24,
                        raster = TRUE) + 
              NoLegend() +
              ggtitle(paste("Marker Genes Batch", i)) +
              theme(axis.text.y = element_text(size = 12))
          print(p)
        }, error = function(e) {
          message(paste("Error generating heatmap batch", i, ":", e$message))
        })
        
        # 主动清理内存
        gc(verbose = FALSE)
      }
      dev.off()
    }
  }, error = function(e) {
    message(paste("Error in heatmap generation:", e$message))
  })
}

print("Data visualization complete!")

# 14. 可视化
pdf("umap_anno_integration.pdf", width = 15, height = 10)
p1 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "study", raster = TRUE, 
              pt.size = 0.5) + ggtitle("Batches")
p2 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "sample", raster = TRUE, 
              pt.size = 0.5) + ggtitle("sample")
p3 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "ann_level_2", raster = TRUE, 
              pt.size = 0.5) + ggtitle("ann_level_2")
p4 <- DimPlot(seurat_obj_main, reduction = "umap", group.by = "tissue", raster = TRUE, 
              pt.size = 0.5) + ggtitle("tissue")
p1
p2
p3
p4
dev.off()

tissues <- unique(seurat_obj_main$tissue)
pdf("tissue_umap_split.pdf", width = 12, height = 12)
Idents(seurat_obj_main) <- "tissue"  # 先设置identity
for(tissue_name in tissues) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = seurat_obj_main, 
                                                    idents = tissue_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(tissue_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

samples <- unique(seurat_obj_main$sample)
pdf("sample_umap_split.pdf", width = 12, height = 12)
Idents(seurat_obj_main) <- "sample"  # 先设置identity
for(sample_name in samples) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = seurat_obj_main, 
                                                    idents = sample_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(sample_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

studys <- unique(seurat_obj_main$study)
pdf("study_umap_split.pdf", width = 12, height = 12)
Idents(seurat_obj_main) <- "study"  # 先设置identity
for(study_name in studys) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(object = seurat_obj_main, 
                                                    idents = study_name),
                cols = "grey",
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(study_name) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()

# 先设置identities
# 先正确设置identities，排除NA值
seurat_obj_main$ann_level_2_no_na <- seurat_obj_main$ann_level_2
seurat_obj_main$ann_level_2_no_na[is.na(seurat_obj_main$ann_level_2_no_na)] <- "Unknown"
Idents(seurat_obj_main) <- "ann_level_2_no_na"

pdf("ann_level_2_umap_split.pdf", width = 12, height = 12)
# 获取唯一的细胞类型（不包括NA和Unknown）
ann_level_2 <- unique(seurat_obj_main$ann_level_2)
ann_level_2 <- ann_level_2[!is.na(ann_level_2)]  # 移除NA值

# 为每个细胞类型绘图
for(anno in ann_level_2) {
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = CellsByIdentities(seurat_obj_main, 
                                                    idents = anno),
                cols = "grey",
                pt.size = 0.5,
                raster = TRUE,
                label = FALSE) +
          ggtitle(anno) +
          theme(legend.text = element_text(size = 12)))
}
dev.off()



# 1. 首先检查 Markers 中是否有重复
print(length(Markers))
print(length(unique(Markers)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Markers[duplicated(Markers)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Markers <- unique(Markers)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot.pdf", width = 32, height = 18)
dot_plot <- DotPlot(seurat_obj_main, 
        features = Markers, 
        group.by = "seurat_clusters",
        split.by = NULL,
        cols = c("lightgrey", "red"),
        dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- dot_plot$data
# 查看或保存数据
View(source_data)
write.csv(source_data, "dotplot_data.csv")
dot_plot
dev.off()
getwd()
# saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")
# 2. 计算每个cluster中每个基因的表达情况
print("Calculating expression statistics...")

# 10. 识别cluster特异性marker基因
print("Finding cluster markers...")
Idents(seurat_obj_main) <- seurat_obj_main$seurat_clusters
markers <- FindAllMarkers(seurat_obj_main, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "cluster_markers.csv")

# table(subset(seurat_obj_main, subset = seurat_clusters =='85')@meta.data[["cell_type"]])

# 11. 生成热图
print("Generating heatmap...")
top10_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)


# 1. 先缩放marker基因
marker_genes <- unique(top10_markers$gene)
seurat_obj_main <- ScaleData(seurat_obj_main, features = marker_genes)

# 1. 将marker基因分成较小的批次
batch_size <- 25  # 每批50个基因
marker_batches <- split(marker_genes, ceiling(seq_along(marker_genes)/batch_size))

# 生成文件名
filename <- "marker_heatmap_batch.pdf"
pdf(filename, width = 32, height = 18)
# 2. 为每个批次生成热图
for(i in seq_along(marker_batches)) {
  print(paste("Processing batch", i, "of", length(marker_batches)))
  
  # 绘制热图
  print(DoHeatmap(seurat_obj_main, 
                  features = marker_batches[[i]],
                  size = 24) + 
          NoLegend() +
          ggtitle(paste("Marker Genes Batch", i))+
          theme(axis.text.y = element_text(size = 36))
  )
  # dev.off()
}
dev.off()
# 

# seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
# 1. 读取注释文件（无表头）
annotations <- read.csv("annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- seurat_obj_main$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
seurat_obj_main$Annotation <- new_idents
seurat_obj_main <- SetIdent(seurat_obj_main, value = new_idents)

# 4. 设置颜色
cell_colors <- c(
  "Epithelial" = "#E41A1C",   
  "Fibroblast" = "#377EB8",    
  "T" = "#4DAF4A",             
  "Myeloid" = "#984EA3",       
  "Endothelial" = "#FF7F00",   
  "SMC" = "#FFFF33",           
  "B" = "#A65628",             
  "Proliferation" = "#F781BF"  
)

# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(seurat_obj_main, 
        reduction = "umap",
        raster = TRUE,
        label = TRUE,
        pt.size = 0.5,
        label.size = 4,
        # cols = cell_colors
        ) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 6. 绘制分割版本
pdf("cell_types_umap_split.pdf", width = 12, height = 12)
# 获取所有细胞类型
cell_types <- unique(seurat_obj_main$Annotation)
cell_types <- cell_types[!is.na(cell_types)]  # 移除NA值

# 为每个细胞类型单独绘图
for(cell_type in cell_types) {
  # 创建文件名
  
  
  # 绘制该细胞类型的UMAP图
  
  print(DimPlot(seurat_obj_main, 
                reduction = "umap",
                cells.highlight = WhichCells(seurat_obj_main, idents = cell_type),
                cols.highlight = cell_colors[cell_type],
                cols = "grey",  # 其他细胞显示为灰色
                pt.size = 0.5,
                label = FALSE) +
          ggtitle(cell_type) +
          theme(legend.text = element_text(size = 12)))
  
}
dev.off()
# 方法1
table(seurat_obj_main@meta.data$tissue[seurat_obj_main@meta.data$seurat_clusters == "43"])
table(seurat_obj_main@meta.data$sample[seurat_obj_main@meta.data$percent.mt >= 15])
# 5. 保存修改后的 Seurat 对象
# saveRDS(seurat_obj_main, "final_analyzed_seurat_obj_main.rds")

print("Analysis completed successfully!")

results <- analyze_cell_proportions_complete(
  seurat_obj = seurat_obj_main,
  exclude_filter = "tissue_sampling_method != 'scraping' & study != 'Jain_Misharin_2021'",
  output_dir = "proportion_results"
)

table(seurat_obj_main@meta.data$seurat_clusters[seurat_obj_main@meta.data$Annotation=='Epithelial' & seurat_obj_main@meta.data$ann_level_2 =='Myeloid'])

# 假设我们选择 cell_type 为 "T_cells" 的细胞
# 1. 提取Undefined细胞子集
# 从主Seurat对象中筛选出Annotation标注为"Undefined"的细胞
Undefined_obj <- subset(seurat_obj_main, subset = Annotation == "Undefined")

# 2. 标准数据预处理
# 数据标准化，默认使用"LogNormalize"方法
Undefined_obj <- NormalizeData(Undefined_obj)

# 3. 寻找高变异基因
# 使用方差稳定化转换(VST)方法选择3000个高变基因，用于后续降维分析
Undefined_obj <- FindVariableFeatures(Undefined_obj, selection.method = "vst", nfeatures = 3000)

# 4. 数据缩放
# 对所有高变基因进行归一化处理，使其均值为0，方差为1
Undefined_obj <- ScaleData(Undefined_obj, features = VariableFeatures(Undefined_obj))

# 5. 主成分分析
# 使用高变基因进行主成分分析，降低数据维度
Undefined_obj <- RunPCA(Undefined_obj)

# 7. Harmony批次效应校正
# 使用Harmony对PCA结果进行批次校正，减少样本间和组织间的批次效应
Undefined_obj <- RunHarmony(
  object = Undefined_obj,           
  group.by.vars = c("sample",'tissue_sampling_method'),   
  theta = c(1,1),      # Higher theta for more diverse clustering               
  lambda = c(7,7),   # Higher lambda to reduce overcorrection                
  sigma = 0.1,           # Lower sigma for tighter clusters
  nclust = 60,            # Increased number of clusters
  reduction.use = "pca",
  max_iter = 20, 
  early_stop = TRUE,
  dims = 1:40           # More iterations for better convergence
)
Undefined_obj <- RunUMAP(Undefined_obj, 
                           reduction = "harmony", 
                           dims = 1:30, 
                           n.neighbors = 30,
                           min.dist = 0.3,
                           metric = "correlation")
# Find neighbors
Undefined_obj <- FindNeighbors(Undefined_obj, 
                                 reduction = "harmony", 
                                 dims = 1:30,
                                 k.param = 30) 

Undefined_obj <- FindClusters(Undefined_obj,
                                resolution = 2, # 多个分辨率
                                verbose = TRUE)

# 11. 可视化聚类结果
print("生成聚类可视化...")
# 设置默认使用分辨率为2的聚类结果
Idents(Undefined_obj) <- Undefined_obj$seurat_clusters

# 保存UMAP聚类图
pdf("Undefined_umap_clusters.pdf", width = 10, height = 8)
DimPlot(Undefined_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Undefined细胞聚类") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
dev.off()

print("Generating feature plots...")
pdf("umap_FeaturePlot_Undefined.pdf", width = 8, height = 8)
for(marker in Markers) {
  if(marker %in% rownames(Undefined_obj)) {
    print(paste("Processing marker:", marker))
    print(FeaturePlot(Undefined_obj, features = marker, raster = TRUE))
    
  } else {
    print(paste("Marker not found:", marker))
  }
}
dev.off()

# 12. 保存分析结果
# saveRDS(Undefined_obj, "Undefined_analyzed.rds")
print("Undefined细胞分析完成!")

# 15. 批次整合质量评估
print("评估批次整合质量...")
# 获取harmony降维结果和元数据
embeddings <- Embeddings(Undefined_obj, "harmony")
metadata <- Undefined_obj@meta.data

# 抽样计算LISI分数以评估批次整合质量
set.seed(42)
n_cells_sample <- min(5000, nrow(embeddings))
sample_idx <- sample(nrow(embeddings), n_cells_sample)

# 使用lisi包计算整合得分
if(requireNamespace("lisi", quietly = TRUE)) {
  lisi_res <- lisi::compute_lisi(embeddings[sample_idx, 1:30], 
                                 metadata[sample_idx, ], 
                                 c("sample", "tissue"), 
                                 30)
  
  # 保存整合质量指标
  cat(sprintf("\n批次整合质量评估:\n"))
  cat(sprintf("- 样本整合得分(iLISI): %.3f\n", mean(lisi_res[, "sample"])))
  cat(sprintf("- 组织保真度(cLISI): %.3f\n", mean(lisi_res[, "tissue"])))
  
  # 保存到文件
  integration_metrics <- data.frame(
    Metric = c("iLISI_sample", "cLISI_tissue"),
    Value = c(mean(lisi_res[, "sample"]), mean(lisi_res[, "tissue"]))
  )
  write.csv(integration_metrics, "Undefined_integration_metrics.csv", row.names = FALSE)
} else {
  warning("lisi包未安装，跳过批次整合质量评估")
}

# 16. 生成UMAP可视化图：按样本、组织等元数据
# 按样本(sample)分组绘制UMAP图
pdf("Undefined_umap_by_metadata.pdf", width =
      12, height = 10)
p1 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "sample", pt.size = 0.5) + 
  ggtitle("按样本分组")
p2 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "tissue", pt.size = 0.5) + 
  ggtitle("按组织分组")
p3 <- DimPlot(Undefined_obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5) + 
  ggtitle("聚类结果")

print(p1)
print(p2)
print(p3)
dev.off()
print("Analysis completed!")

# 1. 首先检查 Markers 中是否有重复
print(length(Markers))
print(length(unique(Markers)))
# 方法1：直接显示所有重复的元素
duplicated_markers <- Markers[duplicated(Markers)]
print("重复的 markers 是：")
print(duplicated_markers)

# 2. 如果发现有重复，可以移除重复项
Markers <- unique(Markers)

# 1. 生成点状图
print("Generating DotPlot...")
pdf("markers_dotplot_Undefined.pdf", width = 32, height = 18)
dot_plot <- DotPlot(Undefined_obj, 
                    features = Markers, 
                    group.by = "seurat_clusters",
                    split.by = NULL,
                    cols = c("lightgrey", "red"),
                    dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)
  )

source_data <- dot_plot$data
# 查看或保存数据
View(source_data)
write.csv(source_data, "dotplot_data_Undefined.csv")
dot_plot
dev.off()
getwd()
# 9. 保存分析结果
print("Saving analysis results...")
# saveRDS(Epithelial_obj, "final_analyzed_Epithelial_obj.rds")

# 1. 读取注释文件（无表头）
annotations <- read.csv("Undefined_annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- Undefined_obj$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)

# 3. 添加到meta.data并设置标识
Undefined_obj$Annotation <- new_idents
Undefined_obj <- SetIdent(Undefined_obj, value = new_idents)

# 4. 设置颜色
cell_colors <- c(
  "Epithelial" = "#E41A1C",   
  "Fibroblast" = "#377EB8",    
  "T" = "#4DAF4A",             
  "Myeloid" = "#984EA3",       
  "Endothelial" = "#FF7F00",   
  "SMC" = "#FFFF33",           
  "B" = "#A65628",             
  "Proliferation" = "#F781BF"  
)

# 5. 绘制UMAP图
pdf("Undefined_cell_types_umap.pdf", width = 8, height = 8)
DimPlot(Undefined_obj, 
        reduction = "umap",
        raster = TRUE,
        label = TRUE,
        pt.size = 0.5,
        label.size = 4,
        # cols = cell_colors
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()

# 获取Undefined_obj中细胞的索引
cells_to_update <- rownames(Undefined_obj@meta.data)

# 检查这些细胞是否都存在于主对象中
cells_in_main <- cells_to_update %in% rownames(seurat_obj_main@meta.data)
valid_cells_to_update <- cells_to_update[cells_in_main]

# 只更新主对象中存在的细胞的Annotation
seurat_obj_main@meta.data[valid_cells_to_update, "Annotation"] <- 
  Undefined_obj@meta.data[valid_cells_to_update, "Annotation"]

# 验证更新
table(seurat_obj_main@meta.data$Annotation)

Idents(seurat_obj_main) <- "Annotation"
# 5. 绘制UMAP图
pdf("cell_types_umap.pdf", width = 8, height = 8)
DimPlot(seurat_obj_main, 
        reduction = "umap",
        raster = TRUE,
        label = TRUE,
        pt.size = 0.5,
        label.size = 4,
        # cols = cell_colors
) +
  ggtitle("Cell Types") +
  theme(legend.text = element_text(size = 12))
dev.off()
rm(Undefined_obj)
gc()
# seurat_obj_main <- readRDS("final_analyzed_seurat_obj_main.rds")
annotations <- read.csv("annotation.csv", header = FALSE)
colnames(annotations) <- c("Cluster", "Annotation")

# 2. 创建新的标识（使用seurat_clusters匹配）
current_clusters <- seurat_obj_main$seurat_clusters  # 获取当前的cluster标识
new_idents <- annotations$Annotation[match(current_clusters, annotations$Cluster)]
names(new_idents) <- names(current_clusters)
cell_types <- unique(seurat_obj_main@meta.data$Annotation)
cell_types <- cell_types[!is.na(cell_types)]
cell_types <- cell_types[! cell_types %in% c('Proliferation','Epithelial','Endothelial','Stromal')]
cell_types
# rm(Endothelial_obj,Epithelial_obj,Fibroblast_obj,Myeloid_obj,SMC_obj,T_obj,B_obj)
for (cell_type in cell_types) {
  new_var_name <- paste(cell_type, 'obj', sep = "_")
  Marker_current <- get(paste0('Marker_',cell_type))
  assign(new_var_name, run_subset_analysis(
    seurat_obj_main,
    cell_type = cell_type,
    markers = Marker_current,
    output_dir = paste(cell_type, "analysis", sep = "_")
  ))
}
getwd()
#' 基于MASC的单细胞组织间两两比较分析函数（修复版）
#' 
#' 使用混合效应模型对所有组织对进行两两比较分析
#' 
#' @param seurat_obj Seurat对象
#' @param cell_type_col 细胞类型列名，默认"Annotation"
#' @param sample_col 样本/供体列名，默认"sample"，作为随机效应
#' @param contrast_col 组织/条件列名，默认"tissue" 
#' @param fixed_effects_cols 固定效应列名向量，默认NULL
#' @param exclude_filter 排除条件的表达式字符串
#' @param output_dir 输出目录
#' @param p_threshold p值显著性阈值，默认0.05
#' @param fdr_threshold FDR显著性阈值，默认0.1
#' @param min_cells 每种细胞类型的最小细胞数，默认10
#' @param min_prop 考虑分析的最小细胞比例，默认0.001
#' @param color_palette 用于绘图的颜色列表
#' @param save_models 是否保存MASC模型，默认FALSE
#' @param min_samples 每个组织的最小样本数，默认2
#' @param specific_pairs 指定要比较的组织对的列表，默认NULL表示比较所有可能的组织对
#' @param seed 随机数种子
#' @return 分析结果列表
scPairwiseMASCAnalysis <- function(seurat_obj, 
                                   cell_type_col = "Annotation", 
                                   sample_col = "sample", 
                                   contrast_col = "tissue", 
                                   fixed_effects_cols = NULL,
                                   exclude_filter = NULL, 
                                   output_dir = "pairwise_MASC_analysis",
                                   p_threshold = 0.05,
                                   fdr_threshold = 0.1,
                                   min_cells = 10,
                                   min_prop = 0.001,
                                   color_palette = NULL,
                                   save_models = FALSE,
                                   min_samples = 2,
                                   specific_pairs = NULL,
                                   seed = 42) {
  
  # 设置随机数种子
  set.seed(seed)
  
  # 创建输出目录
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 记录函数开始时间
  start_time <- Sys.time()
  message(paste("Starting pairwise MASC cell association analysis:", format(start_time)))
  
  #--------------------------------------------------
  # 1. 数据准备
  #--------------------------------------------------
  message("Preparing data...")
  
  # 确保必要的包已加载
  required_packages <- c("Seurat", "dplyr", "ggplot2", "reshape2", 
                         "tidyr", "lme4", "cowplot", "utils")
  
  # 检查并加载包
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed. Please install it."))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # 检查pheatmap包
  run_heatmap <- requireNamespace("pheatmap", quietly = TRUE)
  
  # 获取元数据
  meta_data <- seurat_obj@meta.data
  
  # 检查必要的列
  required_cols <- c(cell_type_col, sample_col, contrast_col)
  if (!is.null(fixed_effects_cols)) {
    required_cols <- c(required_cols, fixed_effects_cols)
  }
  
  if (!all(required_cols %in% colnames(meta_data))) {
    stop(paste("Missing required columns:", 
               paste(setdiff(required_cols, colnames(meta_data)), collapse = ", ")))
  }
  
  # 应用筛选条件
  if (!is.null(exclude_filter) && exclude_filter != "") {
    tryCatch({
      filter_expr <- parse(text = paste("subset(meta_data,", exclude_filter, ")"))
      filtered_meta <- eval(filter_expr)
      message(paste("Cells before filtering:", nrow(meta_data)))
      message(paste("Cells after filtering:", nrow(filtered_meta)))
    }, error = function(e) {
      message(paste("Error in filter expression:", e$message))
      message("Using all cells without filtering")
      filtered_meta <- meta_data
    })
  } else {
    filtered_meta <- meta_data
  }
  
  # 移除NA值
  na_cols <- c(cell_type_col, sample_col, contrast_col)
  if (!is.null(fixed_effects_cols)) {
    na_cols <- c(na_cols, fixed_effects_cols)
  }
  
  # 创建逻辑向量指示哪些行有NA值
  has_na <- apply(filtered_meta[, na_cols, drop = FALSE], 1, function(x) any(is.na(x)))
  if (any(has_na)) {
    message(paste("Removing", sum(has_na), "cells with NA values in key columns"))
    filtered_meta <- filtered_meta[!has_na, ]
  }
  
  # 转换因子列
  filtered_meta[[cell_type_col]] <- as.factor(filtered_meta[[cell_type_col]])
  filtered_meta[[sample_col]] <- as.factor(filtered_meta[[sample_col]])
  filtered_meta[[contrast_col]] <- as.factor(filtered_meta[[contrast_col]])
  
  if (!is.null(fixed_effects_cols)) {
    for (col in fixed_effects_cols) {
      filtered_meta[[col]] <- as.factor(filtered_meta[[col]])
    }
  }
  
  # 获取所有组织水平
  all_tissues <- levels(filtered_meta[[contrast_col]])
  n_tissues <- length(all_tissues)
  
  if (n_tissues < 2) {
    stop("Need at least 2 tissues for comparison")
  }
  
  message(paste("Found", n_tissues, "tissues/conditions:", paste(all_tissues, collapse=", ")))
  
  # 统计每个组织的样本数量（修复BUG：使用唯一样本计数）
  samples_by_tissue <- tapply(filtered_meta[[sample_col]], 
                              filtered_meta[[contrast_col]], 
                              function(x) length(unique(x)))
  
  message("Sample count per tissue:")
  for (tissue in names(samples_by_tissue)) {
    message(paste(" -", tissue, ":", samples_by_tissue[tissue], "samples"))
  }
  
  # 检查样本总数
  n_samples <- length(unique(filtered_meta[[sample_col]]))
  n_cell_types <- length(unique(filtered_meta[[cell_type_col]]))
  
  message(paste("Total unique samples/donors:", n_samples))
  message(paste("Number of cell types:", n_cell_types))
  
  # 检查样本数量是否足够进行统计分析
  if (n_samples < 4) {
    warning("Very few samples (<4). Statistical results may not be reliable.")
  }
  
  # 检查每个组织是否有足够样本
  tissues_with_enough_samples <- names(samples_by_tissue)[samples_by_tissue >= min_samples]
  if (length(tissues_with_enough_samples) < 2) {
    stop(paste("Not enough tissues with", min_samples, "or more samples for comparison"))
  }
  
  if (length(tissues_with_enough_samples) < n_tissues) {
    warning(paste("The following tissues have fewer than", min_samples, "samples and will be excluded:"))
    warning(paste(" ", setdiff(names(samples_by_tissue), tissues_with_enough_samples), collapse="\n  "))
    
    # 过滤掉样本不足的组织
    filtered_meta <- filtered_meta[filtered_meta[[contrast_col]] %in% tissues_with_enough_samples, ]
    filtered_meta[[contrast_col]] <- factor(filtered_meta[[contrast_col]])  # 重新设置因子水平
  }
  
  # 更新组织列表
  all_tissues <- levels(filtered_meta[[contrast_col]])
  n_tissues <- length(all_tissues)
  
  #--------------------------------------------------
  # 2. 计算细胞计数和比例（用于可视化和过滤）
  #--------------------------------------------------
  message("Calculating cell counts and proportions...")
  
  # 创建计数表
  counts_table <- table(filtered_meta[[cell_type_col]], filtered_meta[[sample_col]])
  counts_matrix <- as.matrix(counts_table)
  
  # 添加样本的对比信息
  sample_contrast_df <- unique(filtered_meta[, c(sample_col, contrast_col)])
  sample_to_contrast <- setNames(
    as.character(sample_contrast_df[[contrast_col]]), 
    as.character(sample_contrast_df[[sample_col]])
  )
  
  # 检查每种细胞类型的总细胞数
  cell_type_totals <- rowSums(counts_matrix)
  low_count_cell_types <- names(cell_type_totals)[cell_type_totals < min_cells]
  
  if (length(low_count_cell_types) > 0) {
    message(paste("Warning: The following cell types have fewer than", min_cells, "cells:"))
    message(paste(low_count_cell_types, collapse = ", "))
    message("These cell types will be excluded from analysis.")
    
    # 过滤低计数的细胞类型
    filtered_meta <- filtered_meta[!filtered_meta[[cell_type_col]] %in% low_count_cell_types, ]
    counts_matrix <- counts_matrix[cell_type_totals >= min_cells, ]
  }
  
  # 计算每个样本的细胞比例
  props_by_sample <- prop.table(counts_matrix, margin = 2)
  
  # 过滤低比例的细胞类型
  cell_type_mean_props <- rowMeans(props_by_sample)
  low_prop_cell_types <- names(cell_type_mean_props)[cell_type_mean_props < min_prop]
  
  if (length(low_prop_cell_types) > 0) {
    message(paste("Warning: The following cell types have mean proportion <", min_prop, ":"))
    message(paste(low_prop_cell_types, collapse = ", "))
    message("These cell types will be excluded from analysis.")
    
    # 过滤低比例的细胞类型
    filtered_meta <- filtered_meta[!filtered_meta[[cell_type_col]] %in% low_prop_cell_types, ]
    counts_matrix <- counts_matrix[cell_type_mean_props >= min_prop, ]
    props_by_sample <- prop.table(counts_matrix, margin = 2)
  }
  
  # 转换为长格式用于ggplot绘图
  props_df <- reshape2::melt(props_by_sample)
  colnames(props_df) <- c("CellType", "Sample", "Proportion")
  
  # 添加对比信息到props_df
  props_df$Contrast <- sample_to_contrast[as.character(props_df$Sample)]
  
  # 检查是否有NA值
  if (any(is.na(props_df$Contrast))) {
    warning("Some samples couldn't be matched to contrasts. Check sample_to_contrast mapping.")
    # 移除NA值
    props_df <- props_df[!is.na(props_df$Contrast), ]
  }
  
  # 保存比例数据
  write.csv(props_df, file.path(output_dir, "cell_proportions.csv"), row.names = FALSE)
  
  #--------------------------------------------------
  # 3. 可视化细胞组成
  #--------------------------------------------------
  message("Generating cell composition visualizations...")
  
  # 定义一致的颜色调色板
  if(is.null(color_palette)) {
    n_cell_types <- length(unique(props_df$CellType))
    
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector <- unlist(mapply(RColorBrewer::brewer.pal, 
                                  qual_col_pals$maxcolors, 
                                  rownames(qual_col_pals)))
      color_palette <- sample(col_vector, n_cell_types)
    } else {
      # 基础R颜色
      color_palette <- sample(colors(), n_cell_types)
    }
    
    names(color_palette) <- unique(props_df$CellType)
  }
  
  # 样本级别的堆叠柱状图
  p1 <- ggplot(props_df, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right") +
    labs(title = "Cell Type Composition by Sample",
         y = "Proportion", 
         x = "Sample", 
         fill = "Cell Type") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_palette)
  
  # 添加样本分组
  p1 <- p1 + facet_grid(. ~ Contrast, scales = "free_x", space = "free")
  
  # 保存图像
  ggsave(file.path(output_dir, "cell_composition_by_sample.pdf"), 
         p1, width = 14, height = 8)
  
  # 对比级别的堆叠柱状图
  props_by_contrast <- props_df %>%
    group_by(CellType, Contrast) %>%
    summarise(Proportion = mean(Proportion), .groups = "drop")
  
  p2 <- ggplot(props_by_contrast, aes(x = Contrast, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right") +
    labs(title = "Average Cell Type Composition by Contrast",
         y = "Proportion", 
         x = "Contrast Group", 
         fill = "Cell Type") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_palette)
  
  # 保存图像
  ggsave(file.path(output_dir, "cell_composition_by_contrast.pdf"), 
         p2, width = 10, height = 8)
  
  # 为每个细胞类型创建箱线图
  props_df %>%
    group_by(CellType) %>%
    group_walk(function(cell_data, key) {
      cell_type <- key$CellType
      p <- ggplot(cell_data, aes(x = Contrast, y = Proportion, fill = Contrast)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
        theme_bw() +
        labs(title = paste0(cell_type, " - Proportion by Condition"),
             y = "Proportion", 
             x = "Contrast Group") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 12, face = "bold"),
              legend.position = "none") +
        scale_y_continuous(labels = scales::percent)
      
      ggsave(file.path(output_dir, paste0(cell_type, "_boxplot.pdf")), 
             p, width = 8, height = 6)
    })
  
  #--------------------------------------------------
  # 4. 定义MASC函数
  #--------------------------------------------------
  
  # 检查是否在环境中找到MASC
  if (!exists("MASC", envir = .GlobalEnv)) {
    MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                     verbose = FALSE, save_models = FALSE, save_model_dir = NULL) {
      # Check inputs
      if (is.factor(dataset[[contrast]]) == FALSE) {
        stop("Specified contrast term is not coded as a factor in dataset")
      }
      
      # 将cluster转换为character
      cluster <- as.character(cluster)
      # 创建设计矩阵
      designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
      dataset <- cbind(designmat, dataset)
      
      # 创建输出列表
      res <- vector(mode = "list", length = length(unique(cluster)))
      names(res) <- attributes(designmat)$dimnames[[2]]
      
      # 创建模型公式
      if (!is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                              paste0("(1|", random_effects, ")", collapse = " + ")),
                            collapse = " + ")
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
        }
      } else if (!is.null(fixed_effects) && is.null(random_effects)) {
        model_rhs <- paste0(fixed_effects, collapse = " + ")
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
          stop("No random effects specified")
        }
      } else if (is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
        }
      } else {
        model_rhs <- "1" # 仅包含截距
        if (verbose == TRUE) {
          message(paste("Using null model:", "cluster ~", model_rhs))
          stop("No random or fixed effects specified")
        }
      }
      
      # 初始化列表存储每个簇的模型对象
      cluster_models <- vector(mode = "list",
                               length = length(attributes(designmat)$dimnames[[2]]))
      names(cluster_models) <- attributes(designmat)$dimnames[[2]]
      
      # 为每个簇运行嵌套混合效应模型
      for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
        test_cluster <- attributes(designmat)$dimnames[[2]][i]
        if (verbose == TRUE) {
          message(paste("Creating logistic mixed models for", test_cluster))
        }
        null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                       model_rhs), collapse = ""))
        full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                       model_rhs), collapse = ""))
        # 运行空模型和完整混合效应模型
        tryCatch({
          null_model <- lme4::glmer(formula = null_fm, data = dataset,
                                    family = binomial, nAGQ = 1, verbose = 0,
                                    control = glmerControl(optimizer = "bobyqa"))
          full_model <- lme4::glmer(formula = full_fm, data = dataset,
                                    family = binomial, nAGQ = 1, verbose = 0,
                                    control = glmerControl(optimizer = "bobyqa"))
          model_lrt <- anova(null_model, full_model)
          
          # 计算对比项beta的置信区间
          contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
          contrast_ci <- confint.merMod(full_model, method = "Wald",
                                        parm = contrast_lvl2)
          
          # 保存模型对象到列表
          cluster_models[[i]]$null_model <- null_model
          cluster_models[[i]]$full_model <- full_model
          cluster_models[[i]]$model_lrt <- model_lrt
          cluster_models[[i]]$confint <- contrast_ci
        }, error = function(e) {
          warning(paste("Model for cluster", test_cluster, "failed:", e$message))
          cluster_models[[i]]$error <- e$message
        })
      }
      
      # 组织结果为输出数据框
      output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                           size = colSums(designmat))
      
      # 检查哪些模型成功运行
      successful_models <- sapply(cluster_models, function(x) !is.null(x$model_lrt))
      
      if (any(successful_models)) {
        output$model.pvalue <- NA
        output$model.pvalue[successful_models] <- sapply(
          cluster_models[successful_models], 
          function(x) x$model_lrt[["Pr(>Chisq)"]][2]
        )
        
        # 添加FDR校正
        output$FDR <- p.adjust(output$model.pvalue, method = "BH")
        
        # 添加优势比和置信区间
        contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
        output[[paste(contrast_lvl2, "OR", sep = ".")]] <- NA
        output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- NA
        output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- NA
        
        for (i in which(successful_models)) {
          output[[paste(contrast_lvl2, "OR", sep = ".")]][i] <- 
            exp(fixef(cluster_models[[i]]$full_model)[[contrast_lvl2]])
          output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]][i] <- 
            exp(cluster_models[[i]]$confint[contrast_lvl2, "2.5 %"])
          output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]][i] <- 
            exp(cluster_models[[i]]$confint[contrast_lvl2, "97.5 %"])
        }
      }
      
      # 保存模型对象并返回结果
      if (save_models == TRUE) {
        if (is.null(save_model_dir)) {
          save_dir <- getwd()
        } else {
          save_dir <- save_model_dir
          dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
        }
        saveRDS(cluster_models, file = file.path(save_dir, "masc_models.rds"))
        message(paste("Models saved to", file.path(save_dir, "masc_models.rds")))
        return(list(results = output, models = cluster_models))
      } else {
        return(list(results = output, models = cluster_models))
      }
    }
  }
  
  #--------------------------------------------------
  # 5. 构建所有可能的组织对
  #--------------------------------------------------
  message("Building tissue/condition pairs for comparison...")
  
  # 生成所有可能的组织对组合
  if (is.null(specific_pairs)) {
    tissue_pairs <- combn(all_tissues, 2, simplify = FALSE)
    message(paste("Generated", length(tissue_pairs), "tissue pairs for comparison"))
  } else {
    # 验证指定的组织对
    tissue_pairs <- specific_pairs
    valid_pairs <- sapply(tissue_pairs, function(pair) {
      all(pair %in% all_tissues) && length(pair) == 2
    })
    
    if (!all(valid_pairs)) {
      invalid_indices <- which(!valid_pairs)
      warning(paste("Invalid tissue pairs specified at indices:", 
                    paste(invalid_indices, collapse=", ")))
      tissue_pairs <- tissue_pairs[valid_pairs]
    }
    
    message(paste("Using", length(tissue_pairs), "specified tissue pairs for comparison"))
  }
  
  if (length(tissue_pairs) == 0) {
    stop("No valid tissue pairs for comparison")
  }
  
  #--------------------------------------------------
  # 6. 对每个组织对运行MASC分析
  #--------------------------------------------------
  
  # 创建结果存储列表
  pairwise_results <- list()
  
  # 创建进度条
  pb <- utils::txtProgressBar(min = 0, max = length(tissue_pairs), style = 3)
  
  # 对每个组织对进行分析
  # 对每个组织对进行分析
  for (i in 1:length(tissue_pairs)) {
    pair <- tissue_pairs[[i]]
    tissue1 <- pair[1]
    tissue2 <- pair[2]
    
    comparison_name <- paste(tissue2, "vs", tissue1)
    message(paste("\nComparing", comparison_name, "(",i, "of", length(tissue_pairs),")"))
    
    # 创建只包含这两个组织的数据子集
    pair_meta <- filtered_meta[filtered_meta[[contrast_col]] %in% pair, ]
    
    # 修复: 计算每个组织的唯一样本数 - 使用简单的计数方法避免NA值
    sample_tissue_pairs <- unique(pair_meta[, c(sample_col, contrast_col)])
    samples_per_tissue <- table(sample_tissue_pairs[[contrast_col]])
    
    # 确保只有当前组织对
    samples_per_tissue <- samples_per_tissue[names(samples_per_tissue) %in% pair]
    
    message(paste(" Sample counts:", 
                  paste(names(samples_per_tissue), "=", samples_per_tissue, collapse=", ")))
    
    # 确保有足够样本
    if (length(samples_per_tissue) != 2 || any(samples_per_tissue < min_samples)) {
      skip_reason <- if(length(samples_per_tissue) != 2) 
        "missing data for one tissue" 
      else 
        "insufficient samples"
      warning(paste("Skipping", comparison_name, "-", skip_reason))
      next
    }
    
    # 重新编码对比因子以确保tissue1是参照水平
    pair_meta[[contrast_col]] <- factor(pair_meta[[contrast_col]], levels = pair)
    
    # 将数据转换为MASC所需格式
    masc_dataset <- pair_meta[, c(cell_type_col, sample_col, contrast_col)]
    
    # 添加固定效应列
    if (!is.null(fixed_effects_cols)) {
      masc_dataset <- cbind(masc_dataset, pair_meta[, fixed_effects_cols, drop = FALSE])
    }
    
    # 运行MASC分析
    pair_results <- tryCatch({
      MASC(dataset = masc_dataset,
           cluster = masc_dataset[[cell_type_col]],
           contrast = contrast_col,
           random_effects = sample_col,
           fixed_effects = fixed_effects_cols,
           verbose = TRUE,
           save_models = save_models,
           save_model_dir = ifelse(save_models, 
                                   file.path(output_dir, paste0("models_", gsub(" ", "_", comparison_name))), 
                                   NULL))
    }, error = function(e) {
      warning(paste("MASC analysis failed for", comparison_name, ":", e$message))
      return(NULL)
    })
    
    # 如果分析成功
    if (!is.null(pair_results)) {
      # 提取结果并添加组织信息
      results_df <- pair_results$results
      results_df$Tissue1 <- tissue1  # 参照水平
      results_df$Tissue2 <- tissue2  # 对比水平
      results_df$Comparison <- comparison_name
      
      # 将结果添加到列表
      pairwise_results[[comparison_name]] <- list(
        results = results_df,
        models = if(save_models) pair_results$models else NULL
      )
    }
    
    # 更新进度条
    utils::setTxtProgressBar(pb, i)
  }
  
  # 关闭进度条
  close(pb)
  
  # 检查是否有成功的分析
  if (length(pairwise_results) == 0) {
    stop("No successful comparisons were made")
  }
  
  #--------------------------------------------------
  # 7. 整合所有比较结果
  #--------------------------------------------------
  message("\nIntegrating all pairwise comparison results...")
  
  # 提取和组合所有结果
  all_results <- do.call(rbind, lapply(names(pairwise_results), function(comparison) {
    results_df <- pairwise_results[[comparison]]$results
    
    # 找到并统一重命名列
    or_cols <- grep("\\.OR$", colnames(results_df), value = TRUE)
    ci_lower_cols <- grep("\\.OR\\.95pct\\.ci\\.lower$", colnames(results_df), value = TRUE)
    ci_upper_cols <- grep("\\.OR\\.95pct\\.ci\\.upper$", colnames(results_df), value = TRUE)
    
    if (length(or_cols) > 0) {
      colnames(results_df)[colnames(results_df) == or_cols[1]] <- "OR"
      colnames(results_df)[colnames(results_df) == ci_lower_cols[1]] <- "CI_Lower"
      colnames(results_df)[colnames(results_df) == ci_upper_cols[1]] <- "CI_Upper"
    }
    
    # 添加比较名称作为行标识
    results_df$Comparison <- comparison
    return(results_df)
  }))
  
  # 重命名列以匹配细胞类型
  cluster_to_celltype <- function(cluster_name) {
    # 从"clusterCellType"中提取"CellType"
    gsub("^cluster", "", cluster_name)
  }
  
  all_results$CellType <- cluster_to_celltype(all_results$cluster)
  
  # 添加显著性标记
  all_results$Significant_FDR <- all_results$FDR < fdr_threshold
  all_results$Significant_P <- all_results$model.pvalue < p_threshold
  
  # 添加显著性星号
  all_results$Significance <- ifelse(all_results$FDR < 0.001, "***",
                                     ifelse(all_results$FDR < 0.01, "**",
                                            ifelse(all_results$FDR < 0.05, "*", 
                                                   ifelse(all_results$FDR < 0.1, ".", "ns"))))
  
  # 保存整合结果
  write.csv(all_results, file.path(output_dir, "all_pairwise_masc_results.csv"), row.names = FALSE)
  
  # 拆分结果为每个细胞类型的数据框
  cell_type_results <- split(all_results, all_results$CellType)
  
  # 为每个细胞类型保存结果
  for (cell_type in names(cell_type_results)) {
    write.csv(cell_type_results[[cell_type]], 
              file.path(output_dir, paste0(cell_type, "_pairwise_results.csv")), 
              row.names = FALSE)
  }
  
  #--------------------------------------------------
  # 8. 显示显著结果摘要
  #--------------------------------------------------
  
  # 找出显著的比较
  sig_results <- all_results[all_results$Significant_FDR, ]
  
  if (nrow(sig_results) > 0) {
    message(paste("Found", nrow(sig_results), "significant cell type differences (FDR <", fdr_threshold, ")"))
    
    # 按细胞类型分组摘要
    sig_by_celltype <- split(sig_results, sig_results$CellType)
    
    for (cell_type in names(sig_by_celltype)) {
      message(paste("\nCell type:", cell_type))
      ct_results <- sig_by_celltype[[cell_type]]
      
      for (i in 1:nrow(ct_results)) {
        row <- ct_results[i, ]
        direction <- ifelse(row$OR > 1, "enriched", "depleted")
        message(paste0("  - ", row$Comparison, ": ", 
                       direction, " (OR = ", round(row$OR, 2), 
                       ", FDR = ", format(row$FDR, digits = 3), ")"))
      }
    }
  } else {
    message("No significant cell type differences found")
  }
  
  #--------------------------------------------------
  # 9. 可视化结果
  #--------------------------------------------------
  message("\nGenerating visualizations for pairwise comparisons...")
  
  # 为每个细胞类型创建森林图
  for (cell_type in unique(all_results$CellType)) {
    # 提取该细胞类型的结果
    ct_results <- all_results[all_results$CellType == cell_type, ]
    
    # 跳过没有OR的结果
    ct_results <- ct_results[!is.na(ct_results$OR), ]
    
    if (nrow(ct_results) > 0) {
      # 计算对数优势比和置信区间
      ct_results$log_OR <- log(ct_results$OR)
      ct_results$log_CI_Lower <- log(ct_results$CI_Lower)
      ct_results$log_CI_Upper <- log(ct_results$CI_Upper)
      
      # 按对数优势比排序
      ct_results <- ct_results[order(ct_results$log_OR), ]
      
      # 创建比较因子，保持排序
      ct_results$Comparison <- factor(ct_results$Comparison, levels = ct_results$Comparison)
      
      # 森林图
      p_forest <- ggplot(ct_results, aes(x = log_OR, y = Comparison, color = Significant_FDR)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = log_CI_Lower, xmax = log_CI_Upper), height = 0.2) +
        scale_color_manual(values = c("FALSE" = "darkgray", "TRUE" = "red")) +
        labs(
          title = paste("MASC Odds Ratios for", cell_type),
          subtitle = paste("Red: FDR <", fdr_threshold),
          x = "Log Odds Ratio (95% CI)",
          y = "Comparison"
        ) +
        theme_bw() +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(face = ifelse(ct_results$Significant_FDR, "bold", "plain"))
        )
      
      ggsave(file.path(output_dir, paste0(cell_type, "_forest_plot.pdf")), p_forest, width = 10, height = 8)
    }
  }
  
  # 创建热图显示所有细胞类型和比较 - 使用修复的代码
  if (run_heatmap && nrow(all_results) > 0) {
    message("Generating heatmap of log odds ratios across comparisons...")
    
    tryCatch({
      # 准备热图数据
      # 为缺失的OR值添加占位符
      all_results$log_OR <- log(all_results$OR)
      all_results$log_OR[!is.finite(all_results$log_OR)] <- 0
      
      # 宽格式转换：行为细胞类型，列为比较
      heatmap_data <- all_results %>%
        select(CellType, Comparison, log_OR, Significant_FDR)
      
      # 检查是否至少有一行数据
      if(nrow(heatmap_data) > 0) {
        # 重标或截断极端值以便更好地可视化
        max_abs_log_or <- max(abs(heatmap_data$log_OR), na.rm = TRUE)
        cap_value <- min(max_abs_log_or, 3)  # 截断极值为±3
        heatmap_data$log_OR_capped <- pmin(pmax(heatmap_data$log_OR, -cap_value), cap_value)
        
        # 转换为宽格式矩阵
        heat_matrix <- tidyr::pivot_wider(
          heatmap_data, 
          id_cols = CellType,
          names_from = Comparison,
          values_from = log_OR_capped,
          values_fill = 0
        )
        
        # 提取矩阵
        rownames(heat_matrix) <- heat_matrix$CellType
        heat_matrix <- as.matrix(heat_matrix[, -1, drop = FALSE])
        
        # 确保至少有一列
        if(ncol(heat_matrix) > 0) {
          # 提取显著性矩阵
          sig_df <- tidyr::pivot_wider(
            heatmap_data, 
            id_cols = CellType,
            names_from = Comparison,
            values_from = Significant_FDR,
            values_fill = FALSE
          )
          
          rownames(sig_df) <- sig_df$CellType
          sig_df <- sig_df[, -1, drop = FALSE]  # 删除CellType列保留其他列
          
          # 确保行名匹配
          common_rows <- intersect(rownames(heat_matrix), rownames(sig_df))
          common_cols <- intersect(colnames(heat_matrix), colnames(sig_df))
          
          if(length(common_rows) > 0 && length(common_cols) > 0) {
            heat_matrix_subset <- heat_matrix[common_rows, common_cols, drop = FALSE]
            sig_df_subset <- as.matrix(sig_df[common_rows, common_cols, drop = FALSE])
            
            # 创建显著性标记矩阵
            sig_symbols <- matrix("", 
                                  nrow = nrow(heat_matrix_subset), 
                                  ncol = ncol(heat_matrix_subset),
                                  dimnames = dimnames(heat_matrix_subset))
            
            # 安全地填充显著性标记
            for(r in 1:nrow(sig_symbols)) {
              for(c in 1:ncol(sig_symbols)) {
                row_name <- rownames(sig_symbols)[r]
                col_name <- colnames(sig_symbols)[c]
                
                # 使用行名和列名而不是索引来确保安全
                if(!is.na(sig_df_subset[row_name, col_name]) && 
                   sig_df_subset[row_name, col_name]) {
                  sig_symbols[r, c] <- "*"
                }
              }
            }
            
            # 设置列注释 - 组织对
            comparison_split <- strsplit(colnames(heat_matrix_subset), " vs ")
            tissue1 <- sapply(comparison_split, function(x) if(length(x) > 1) x[2] else NA)
            tissue2 <- sapply(comparison_split, function(x) if(length(x) > 0) x[1] else NA)
            
            col_anno <- data.frame(
              Tissue1 = tissue1,
              Tissue2 = tissue2,
              row.names = colnames(heat_matrix_subset)
            )
            
            # 定义颜色
            tissue_colors <- setNames(
              rainbow(length(all_tissues)),
              all_tissues
            )
            
            anno_colors <- list(
              Tissue1 = tissue_colors,
              Tissue2 = tissue_colors
            )
            
            # 创建热图
            pdf(file.path(output_dir, "pairwise_log_OR_heatmap.pdf"), width = 12, height = 10)
            pheatmap::pheatmap(
              heat_matrix_subset,
              main = "Log Odds Ratios Across All Comparisons",
              annotation_col = col_anno,
              annotation_colors = anno_colors,
              display_numbers = sig_symbols,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              color = colorRampPalette(c("blue", "white", "red"))(100),
              breaks = seq(-cap_value, cap_value, length.out = 101)
            )
            dev.off()
            
            message("Heatmap created successfully")
          } else {
            message("No common cell types or comparisons found for heatmap after matching")
          }
        } else {
          message("No comparison columns available for heatmap")
        }
      } else {
        message("No data available for heatmap")
      }
    }, error = function(e) {
      warning(paste("Error generating heatmap:", e$message))
      message("Continuing with other visualizations...")
    })
  }
  # 为细胞类型在各组织中的丰度创建堆叠柱状图
  props_by_contrast %>%
    ggplot(aes(x = Contrast, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Cell Type Composition Across Tissues",
      y = "Proportion",
      x = "Tissue"
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_palette)
  
  ggsave(file.path(output_dir, "cell_composition_horizontal.pdf"), width = 10, height = 8)
  
  # 创建多面板图，显示显著差异的细胞类型比例
  if (nrow(sig_results) > 0) {
    # 获取显著差异的细胞类型
    sig_celltypes <- unique(sig_results$CellType)
    
    # 为每个显著的细胞类型绘制箱线图矩阵
    for (cell_type in sig_celltypes) {
      # 获取该细胞类型的比例数据
      ct_data <- props_df[props_df$CellType == cell_type, ]
      
      # 获取该细胞类型的显著对比
      sig_comparisons <- sig_results[sig_results$CellType == cell_type, ]
      
      # 创建箱线图
      p <- ggplot(ct_data, aes(x = Contrast, y = Proportion, fill = Contrast)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
        theme_bw() +
        labs(
          title = paste0(cell_type, " - Significant Differences"),
          subtitle = paste("Number of significant comparisons:", nrow(sig_comparisons)),
          y = "Proportion",
          x = "Tissue"
        ) +
        scale_y_continuous(labels = scales::percent) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold")
        )
      
      # 保存图像
      ggsave(file.path(output_dir, paste0(cell_type, "_significant_boxplot.pdf")), 
             p, width = 10, height = 6)
    }
  }
  
  # 记录分析结束时间
  end_time <- Sys.time()
  message(paste("\nPairwise MASC analysis completed, duration:", format(end_time - start_time)))
  
  # 返回结果
  return(list(
    counts = counts_matrix,
    proportions = props_df,
    all_results = all_results,
    pairwise_results = pairwise_results,
    significant_results = sig_results
  ))
}

results <- scPairwiseMASCAnalysis(
  seurat_obj = seurat_obj_main,
  cell_type_col = "Annotation",
  sample_col = "sample",
  contrast_col = "tissue",
  exclude_filter = "tissue_sampling_method != 'scraping'",
  fixed_effects_cols = NULL,  # 添加任何固定效应协变量
  output_dir = "pairwise_MASC_results",
  # 如果您需要少于2个样本也能进行比较，可以降低此阈值
  min_samples = 2
)

