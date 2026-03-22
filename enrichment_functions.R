# =================== 10. 添加富集分析函数 ===================
# 安装和加载所需的包
install_packages <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # List of packages to install
  bioc_packages <- c("clusterProfiler", "org.Hs.eg.db", "pathview", "DOSE")
  cran_packages <- c("ggplot2", "dplyr", "tidyr", "plotly", "htmlwidgets")
  
  # Install Bioconductor packages
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Install CRAN packages
  for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  cat("All required packages have been installed and loaded.\n")
}

# Function to create output directories
create_output_dirs <- function(prefix) {
  # Create main output directory
  main_dir <- file.path(getwd(), prefix)
  dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create subdirectories for GO and KEGG
  go_dir <- file.path(main_dir, "go")
  kegg_dir <- file.path(main_dir, "kegg")
  dir.create(go_dir, showWarnings = FALSE)
  dir.create(kegg_dir, showWarnings = FALSE)
  
  return(list(
    main_dir = main_dir,
    go_dir = go_dir,
    kegg_dir = kegg_dir
  ))
}

# Function to preprocess data from CSV or data frame
preprocess_data <- function(data_input) {
  # Check if input is a file path or data frame
  if (is.character(data_input)) {
    # Read the CSV file with proper column types
    data <- read.csv(data_input, stringsAsFactors = FALSE)
  } else if (is.data.frame(data_input)) {
    # Input is already a data frame
    data <- data_input
  } else {
    stop("Input must be either a CSV file path or a data frame")
  }
  
  # Verify expected columns
  expected_cols <- c("p_val", "logFC", "AveExpr", "p_val_adj", "gene_names")
  missing_cols <- setdiff(expected_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    # Try to map common column names
    col_mapping <- list(
      p_val = c("p_val", "p.value", "pvalue", "PValue", "p_value"),
      logFC = c("logFC", "avg_log2FC", "log2FC", "log2FoldChange", "logfoldchanges"),
      AveExpr = c("AveExpr", "baseMean", "AveExpr", "mean.expr"),
      p_val_adj = c("p_val_adj", "padj", "p.adjust", "p_adj", "adjusted.p.val", "FDR"),
      gene_names = c("gene_names", "gene", "Gene", "gene_id", "ID", "gene_name", "Symbol")
    )
    
    for (expected_col in missing_cols) {
      possible_names <- col_mapping[[expected_col]]
      for (alt_name in possible_names) {
        if (alt_name %in% colnames(data)) {
          data[[expected_col]] <- data[[alt_name]]
          cat(paste("Mapped column", alt_name, "to", expected_col, "\n"))
          break
        }
      }
    }
    
    # Check again after mapping
    missing_cols <- setdiff(expected_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
    }
  }
  
  # If gene_names column is missing but rownames contain gene names
  if (!"gene_names" %in% colnames(data) && is.character(rownames(data))) {
    data$gene_names <- rownames(data)
  }
  
  # Print data summary
  cat("\nData Summary:\n")
  cat("Total genes:", nrow(data), "\n")
  cat("Non-NA genes:", sum(complete.cases(data)), "\n")
  
  # Remove any NA values
  data <- na.omit(data)
  
  return(data)
}

# Function to prepare gene lists
prepare_gene_lists <- function(data, log2fc_cutoff, pval_cutoff) {
  # Filter DEGs
  deg_list <- data %>%
    dplyr::filter(abs(logFC) > log2fc_cutoff) %>%
    dplyr::filter(p_val_adj < pval_cutoff)
  
  # Split into up and down regulated genes
  up_genes <- deg_list %>% 
    dplyr::filter(logFC > log2fc_cutoff) %>%
    pull(gene_names)
  
  down_genes <- deg_list %>% 
    dplyr::filter(logFC < -log2fc_cutoff) %>%
    pull(gene_names)
  
  cat("\nDifferential Expression Analysis:\n")
  cat("Upregulated genes:", length(up_genes), "\n")
  cat("Downregulated genes:", length(down_genes), "\n")
  
  # Convert to ENTREZ IDs
  up_entrez <- mapIds(org.Hs.eg.db,
                      keys = up_genes,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")
  
  down_entrez <- mapIds(org.Hs.eg.db,
                        keys = down_genes,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
  
  # Prepare gene list for GSEA
  all_entrez <- mapIds(org.Hs.eg.db,
                       keys = data$gene_names,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  gene_list <- data$logFC
  names(gene_list) <- all_entrez
  gene_list <- gene_list[!is.na(names(gene_list))]
  
  cat("Mapped to ENTREZ IDs:\n")
  cat("Upregulated:", length(na.omit(up_entrez)), "\n")
  cat("Downregulated:", length(na.omit(down_entrez)), "\n")
  cat("Total genes for GSEA:", length(gene_list), "\n")
  
  return(list(
    up_entrez = na.omit(up_entrez),
    down_entrez = na.omit(down_entrez),
    gene_list = gene_list
  ))
}

# Function to process GO enrichment result
process_go_result <- function(ego_obj, ont) {
  if(is.null(ego_obj) || nrow(ego_obj@result) == 0) {
    return(NULL)
  }
  
  # Convert to data frame and add ONTOLOGY column first
  result_df <- data.frame(
    ONTOLOGY = rep(ont, nrow(ego_obj@result)),
    ID = ego_obj@result$ID,
    Description = ego_obj@result$Description,
    GeneRatio = ego_obj@result$GeneRatio,
    BgRatio = ego_obj@result$BgRatio,
    pvalue = ego_obj@result$pvalue,
    p.adjust = ego_obj@result$p.adjust,
    qvalue = ego_obj@result$qvalue,
    geneID = ego_obj@result$geneID,
    Count = ego_obj@result$Count,
    stringsAsFactors = FALSE
  )
  
  # Remove any NA rows
  result_df <- result_df[complete.cases(result_df), ]
  
  # Sort by p.adjust
  result_df <- result_df[order(result_df$p.adjust), ]
  
  return(result_df)
}

# Function to classify KEGG pathway
get_kegg_category <- function(description) {
  # Define standard KEGG categories and their keywords
  categories <- list(
    "Human Diseases" = c("cancer", "disease", "infection", "disorder", "diabetes", "obesity"),
    "Organismal Systems" = c("system", "immune", "endocrine", "circulatory", "digestive", "nervous", "aging"),
    "Environmental Information Processing" = c("signal", "signaling", "interaction", "transduction"),
    "Cellular Processes" = c("transport", "catabolism", "autophagy", "cycle", "death", "motility"),
    "Metabolism" = c("metabolism", "metabolic", "biosynthesis", "degradation"),
    "Genetic Information Processing" = c("transcription", "translation", "replication", "repair", "ribosome")
  )
  
  # Convert description to lower case for case-insensitive matching
  desc_lower <- tolower(description)
  
  # Check each category
  for(cat in names(categories)) {
    if(any(sapply(categories[[cat]], function(x) grepl(x, desc_lower)))) {
      return(cat)
    }
  }
  
  # If no category matches, return NA
  return("NA")
}

# Function to process enrichment results
process_enrichment_result <- function(enrichment_obj, type = "GO") {
  if(is.null(enrichment_obj) || nrow(enrichment_obj@result) == 0) {
    return(NULL)
  }
  
  # Convert to data frame
  result_df <- data.frame(
    ID = enrichment_obj@result$ID,
    Description = enrichment_obj@result$Description,
    GeneRatio = enrichment_obj@result$GeneRatio,
    BgRatio = enrichment_obj@result$BgRatio,
    pvalue = enrichment_obj@result$pvalue,
    p.adjust = enrichment_obj@result$p.adjust,
    qvalue = enrichment_obj@result$qvalue,
    geneID = enrichment_obj@result$geneID,
    Count = enrichment_obj@result$Count,
    stringsAsFactors = FALSE
  )
  
  # Add category information
  if(type == "GO") {
    result_df$Ontology <- sapply(strsplit(result_df$ID, ":"), `[`, 1)
  } else {
    result_df$MainClass <- sapply(result_df$Description, get_kegg_category)
  }
  
  return(result_df)
}

# Function to save GO visualizations
save_go_plots <- function(enrichment_obj, prefix, output_dir) {
  if(is.null(enrichment_obj) || nrow(enrichment_obj) == 0) {
    warning("No enrichment results for ", prefix)
    return()
  }
  
  # Filter for significant results first
  enrichment_obj <- enrichment_obj[enrichment_obj$p.adjust <= 0.05, ]
  
  if(nrow(enrichment_obj) == 0) {
    warning("No significant results (p.adjust <= 0.05) for ", prefix)
    return()
  }
  
  # Add numeric gene ratio for sorting
  enrichment_obj$GeneRatio_num <- sapply(strsplit(as.character(enrichment_obj$GeneRatio), "/"), 
                                         function(x) as.numeric(x[1])/as.numeric(x[2]))
  
  # Get top 10 terms for each ontology
  go_data <- do.call(rbind, lapply(unique(enrichment_obj$ONTOLOGY), function(ont) {
    ont_data <- subset(enrichment_obj, ONTOLOGY == ont)
    ont_data <- ont_data[order(-ont_data$GeneRatio_num), ]
    head(ont_data, 10)
  }))
  
  if(!is.null(go_data) && nrow(go_data) > 0) {
    # Create proper factor levels for ordering
    go_data$Description <- factor(go_data$Description,
                                  levels = rev(go_data$Description))
    go_data$ONTOLOGY <- factor(go_data$ONTOLOGY,
                               levels = c("BP", "CC", "MF"),
                               labels = c("Biological Process",
                                          "Cellular Component",
                                          "Molecular Function"))
    
    # Create combined plot
    p <- ggplot(go_data, 
                aes(x = GeneRatio_num, y = Description)) +
      geom_point(aes(size = Count, color = p.adjust)) +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
      scale_color_gradientn(
        colors = c("red", "purple", "blue"),
        limits = c(0, 0.05),
        breaks = seq(0, 0.05, by = 0.01),
        guide = guide_colorbar(nbin = 50)
      ) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8),
            strip.text.y = element_text(angle = 0)) +
      labs(title = paste(prefix, "GO Analysis"),
           x = "Gene Ratio",
           y = "GO Term",
           color = "Adjusted p-value",
           size = "Gene Count")
    
    # Save combined plot
    pdf(file.path(output_dir, paste0(prefix, "_GO_combined.pdf")),
        width = 12, height = max(15, nrow(go_data) * 0.3))
    print(p)
    dev.off()
    
    # Create individual ontology plots
    for(ont in unique(enrichment_obj$ONTOLOGY)) {
      ont_data <- subset(enrichment_obj, ONTOLOGY == ont)
      if(nrow(ont_data) > 0) {
        # Sort by Gene Ratio and take top 20
        ont_data <- ont_data[order(-ont_data$GeneRatio_num), ]
        ont_data <- head(ont_data, 20)
        
        # Dotplot
        p1 <- ggplot(ont_data, 
                     aes(x = GeneRatio_num, 
                         y = reorder(Description, GeneRatio_num))) +
          geom_point(aes(size = Count, color = p.adjust)) +
          scale_color_gradientn(
            colors = c("red", "purple", "blue"),
            limits = c(0, 0.05),
            breaks = seq(0, 0.05, by = 0.01),
            guide = guide_colorbar(nbin = 50)
          ) +
          theme_bw() +
          theme(axis.text.y = element_text(size = 8)) +
          labs(title = paste(prefix, ont, "GO Terms"),
               x = "Gene Ratio",
               y = "GO Term",
               color = "Adjusted p-value",
               size = "Gene Count")
        
        pdf(file.path(output_dir, paste0(prefix, "_", ont, "_dotplot.pdf")),
            width = 12, height = 8)
        print(p1)
        dev.off()
      }
    }
  }
  
  # Save results by ontology
  for(ont in unique(enrichment_obj$ONTOLOGY)) {
    ont_data <- subset(enrichment_obj, ONTOLOGY == ont)
    if(nrow(ont_data) > 0) {
      # Sort by p.adjust first, then by GeneRatio
      ont_data <- ont_data[order(ont_data$p.adjust, -ont_data$GeneRatio_num), ]
      
      # Save full results
      write.csv(ont_data,
                file.path(output_dir, paste0(prefix, "_", ont, ".csv")),
                row.names = FALSE)
      
      # Save top results without quotes
      write.table(head(ont_data, 50),
                  file.path(output_dir, paste0(prefix, "_", ont, "_top50.txt")),
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
}

# Function to save KEGG visualizations
save_kegg_plots <- function(kegg_data, prefix, output_dir) {
  if(is.null(kegg_data) || nrow(kegg_data) == 0) {
    warning("No enrichment results for ", prefix)
    return()
  }
  
  # Filter for significant results first
  kegg_data <- kegg_data[kegg_data$p.adjust <= 0.05, ]
  
  if(nrow(kegg_data) == 0) {
    warning("No significant results (p.adjust <= 0.05) for ", prefix)
    return()
  }
  
  # Check if this is GSEA result
  is_gsea <- "NES" %in% colnames(kegg_data)
  
  # Add necessary columns for visualization
  if(is_gsea) {
    kegg_data$GeneRatio_num <- abs(kegg_data$NES)
    kegg_data$size_metric <- kegg_data$setSize
  } else {
    kegg_data$GeneRatio_num <- sapply(strsplit(as.character(kegg_data$GeneRatio), "/"), 
                                      function(x) as.numeric(x[1])/as.numeric(x[2]))
    kegg_data$size_metric <- kegg_data$Count
  }
  
  # Add category information if not already present
  if(!"MainClass" %in% colnames(kegg_data)) {
    kegg_data$MainClass <- sapply(kegg_data$Description, get_kegg_category)
  }
  
  # Order categories
  category_order <- c("Human Diseases", "Organismal Systems", 
                      "Environmental Information Processing", "Cellular Processes",
                      "Metabolism", "Genetic Information Processing", "NA")
  
  # Get top pathways for each category
  kegg_data_top <- do.call(rbind, lapply(category_order, function(cat) {
    cat_data <- subset(kegg_data, MainClass == cat)
    if(nrow(cat_data) > 0) {
      if(is_gsea) {
        cat_data <- cat_data[order(-abs(cat_data$NES)), ]
      } else {
        cat_data <- cat_data[order(-cat_data$GeneRatio_num), ]
      }
      head(cat_data, 10)
    }
  }))
  
  if(!is.null(kegg_data_top) && nrow(kegg_data_top) > 0) {
    # Create proper factor levels for ordering
    kegg_data_top$Description <- factor(kegg_data_top$Description,
                                        levels = rev(kegg_data_top$Description))
    kegg_data_top$MainClass <- factor(kegg_data_top$MainClass,
                                      levels = category_order)
    
    # Create combined plot
    p <- ggplot(kegg_data_top, 
                aes(x = GeneRatio_num, y = Description)) +
      geom_point(aes(size = size_metric, color = p.adjust)) +
      facet_grid(MainClass ~ ., scales = "free_y", space = "free_y") +
      scale_color_gradientn(
        colors = c("red", "purple", "blue"),
        limits = c(0, 0.05),
        breaks = seq(0, 0.05, by = 0.01),
        guide = guide_colorbar(nbin = 50)
      ) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8),
            strip.text.y = element_text(angle = 0)) +
      labs(title = paste(prefix, "KEGG Pathways"),
           x = if(is_gsea) "Normalized Enrichment Score (absolute)" else "Gene Ratio",
           y = "Pathway",
           color = "Adjusted p-value",
           size = if(is_gsea) "Gene Set Size" else "Gene Count")
    
    # Save combined plot
    pdf(file.path(output_dir, paste0(prefix, "_KEGG_combined.pdf")),
        width = 12, height = max(15, nrow(kegg_data_top) * 0.3))
    print(p)
    dev.off()
    
    # Create individual category plots
    for(cat in category_order) {
      cat_data <- subset(kegg_data, MainClass == cat)
      if(nrow(cat_data) > 0) {
        # Sort and take top 20
        if(is_gsea) {
          cat_data <- cat_data[order(-abs(cat_data$NES)), ]
        } else {
          cat_data <- cat_data[order(-cat_data$GeneRatio_num), ]
        }
        cat_data <- head(cat_data, 20)
        
        # Dotplot
        p1 <- ggplot(cat_data, 
                     aes(x = GeneRatio_num, 
                         y = reorder(Description, GeneRatio_num))) +
          geom_point(aes(size = size_metric, color = p.adjust)) +
          scale_color_gradientn(
            colors = c("red", "purple", "blue"),
            limits = c(0, 0.05),
            breaks = seq(0, 0.05, by = 0.01),
            guide = guide_colorbar(nbin = 50)
          ) +
          theme_bw() +
          theme(axis.text.y = element_text(size = 8)) +
          labs(title = paste(prefix, cat, "Pathways"),
               x = if(is_gsea) "Normalized Enrichment Score (absolute)" else "Gene Ratio",
               y = "Pathway",
               color = "Adjusted p-value",
               size = if(is_gsea) "Gene Set Size" else "Gene Count")
        
        pdf(file.path(output_dir, paste0(prefix, "_", make.names(cat), "_dotplot.pdf")),
            width = 12, height = 8)
        print(p1)
        dev.off()
        
        # Save results
        if(is_gsea) {
          cat_data <- cat_data[order(cat_data$p.adjust, -abs(cat_data$NES)), ]
        } else {
          cat_data <- cat_data[order(cat_data$p.adjust, -cat_data$GeneRatio_num), ]
        }
        write.csv(cat_data,
                  file.path(output_dir, paste0(prefix, "_", make.names(cat), ".csv")),
                  row.names = FALSE)
      }
    }
  }
  
  # Save complete results
  kegg_data <- kegg_data[order(kegg_data$MainClass, kegg_data$p.adjust), ]
  write.csv(kegg_data, 
            file.path(output_dir, paste0(prefix, "_KEGG_all.csv")),
            row.names = FALSE)
}


# Main GO analysis function (works with both CSV file or data frame)
perform_go_analysis <- function(data_input, log2fc_cutoff = 1, pval_cutoff = 0.05, output_prefix = NULL) {
  # Determine output prefix
  if(is.null(output_prefix)) {
    if(is.character(data_input)) {
      output_prefix <- tools::file_path_sans_ext(basename(data_input))
    } else {
      output_prefix <- paste0("go_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    }
  }
  
  # Create output directories
  dirs <- create_output_dirs(output_prefix)
  output_dir <- dirs$go_dir
  
  # Process data
  data <- preprocess_data(data_input)
  gene_lists <- prepare_gene_lists(data, log2fc_cutoff, pval_cutoff)
  
  # Perform GO enrichment separately for each ontology
  cat("\nPerforming GO enrichment analysis...\n")
  
  # Define ontologies
  ontologies <- c("BP", "CC", "MF")
  ont_names <- c(BP = "Biological Process",
                 CC = "Cellular Component",
                 MF = "Molecular Function")
  
  # Initialize results lists
  go_up_list <- list()
  go_down_list <- list()
  
  # Perform analysis for each ontology
  for(ont in ontologies) {
    cat("\nAnalyzing", ont_names[ont], "terms...\n")
    
    # Up-regulated genes
    go_up_list[[ont]] <- tryCatch({
      ego_up <- enrichGO(gene = gene_lists$up_entrez,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = ont,
                         pAdjustMethod = "BH",
                         pvalueCutoff = pval_cutoff,
                         qvalueCutoff = 0.2,
                         readable = TRUE)
      
      cat("Up-regulated ", ont, " terms: ", nrow(ego_up@result), "\n")
      process_go_result(ego_up, ont)
      
    }, error = function(e) {
      warning("Error in up-regulated ", ont, " analysis: ", e$message)
      NULL
    })
    
    # Down-regulated genes
    go_down_list[[ont]] <- tryCatch({
      ego_down <- enrichGO(gene = gene_lists$down_entrez,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = ont,
                           pAdjustMethod = "BH",
                           pvalueCutoff = pval_cutoff,
                           qvalueCutoff = 0.2,
                           readable = TRUE)
      
      cat("Down-regulated ", ont, " terms: ", nrow(ego_down@result), "\n")
      process_go_result(ego_down, ont)
      
    }, error = function(e) {
      warning("Error in down-regulated ", ont, " analysis: ", e$message)
      NULL
    })
  }
  
  # Combine results
  go_up <- do.call(rbind, go_up_list)
  go_down <- do.call(rbind, go_down_list)
  
  # Save results and create plots
  if(!is.null(go_up) && nrow(go_up) > 0) {
    write.csv(go_up, file.path(output_dir, "upregulated_GO_full.csv"), 
              row.names = FALSE)
    save_go_plots(go_up, "GO_up", output_dir)
  }
  
  if(!is.null(go_down) && nrow(go_down) > 0) {
    write.csv(go_down, file.path(output_dir, "downregulated_GO_full.csv"), 
              row.names = FALSE)
    save_go_plots(go_down, "GO_down", output_dir)
  }
  
  # Create summary file
  summary_data <- data.frame(
    Category = c("Total Genes", "Upregulated", "Downregulated",
                 paste0("GO Terms ", ontologies, " (Up)"),
                 paste0("GO Terms ", ontologies, " (Down)")),
    Count = c(nrow(data),
              length(gene_lists$up_entrez),
              length(gene_lists$down_entrez),
              sapply(ontologies, function(ont) {
                if(is.null(go_up_list[[ont]])) 0 else nrow(go_up_list[[ont]])
              }),
              sapply(ontologies, function(ont) {
                if(is.null(go_down_list[[ont]])) 0 else nrow(go_down_list[[ont]])
              }))
  )
  write.csv(summary_data, file.path(output_dir, "analysis_summary.csv"), 
            row.names = FALSE)
  
  return(list(
    upregulated = go_up,
    downregulated = go_down,
    upregulated_by_ont = go_up_list,
    downregulated_by_ont = go_down_list
  ))
}

# Main KEGG analysis function (works with both CSV file or data frame)
perform_kegg_analysis <- function(data_input, log2fc_cutoff = 1, pval_cutoff = 0.05, output_prefix = NULL) {
  # Determine output prefix
  if(is.null(output_prefix)) {
    if(is.character(data_input)) {
      output_prefix <- tools::file_path_sans_ext(basename(data_input))
    } else {
      output_prefix <- paste0("kegg_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    }
  }
  
  # Create output directories
  dirs <- create_output_dirs(output_prefix)
  output_dir <- dirs$kegg_dir
  
  # Process data
  data <- preprocess_data(data_input)
  gene_lists <- prepare_gene_lists(data, log2fc_cutoff, pval_cutoff)
  
  # Perform KEGG enrichment
  cat("\nPerforming KEGG pathway analysis...\n")
  
  # Up-regulated genes analysis
  cat("\nAnalyzing up-regulated KEGG pathways...\n")
  kegg_up <- tryCatch({
    ekegg_up <- enrichKEGG(
      gene = gene_lists$up_entrez,
      organism = 'hsa',
      pvalueCutoff = pval_cutoff,
      pAdjustMethod = "BH",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500
    )
    result <- process_enrichment_result(ekegg_up, type = "KEGG")
    if(!is.null(result)) {
      # Filter by adjusted p-value
      result <- result[result$p.adjust < pval_cutoff, ]
      cat("Found", nrow(result), "significantly enriched pathways\n")
    }
    result
  }, error = function(e) {
    warning("Error in up-regulated KEGG analysis: ", e$message)
    NULL
  })
  
  # Down-regulated genes analysis
  cat("\nAnalyzing down-regulated KEGG pathways...\n")
  kegg_down <- tryCatch({
    ekegg_down <- enrichKEGG(
      gene = gene_lists$down_entrez,
      organism = 'hsa',
      pvalueCutoff = pval_cutoff,
      pAdjustMethod = "BH",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500
    )
    result <- process_enrichment_result(ekegg_down, type = "KEGG")
    if(!is.null(result)) {
      # Filter by adjusted p-value
      result <- result[result$p.adjust < pval_cutoff, ]
      cat("Found", nrow(result), "significantly enriched pathways\n")
    }
    result
  }, error = function(e) {
    warning("Error in down-regulated KEGG analysis: ", e$message)
    NULL
  })
  
  # GSEA analysis
  cat("\nPerforming KEGG GSEA analysis...\n")
  kegg_gsea <- tryCatch({
    ekegg_gsea <- gseKEGG(
      geneList = sort(gene_lists$gene_list, decreasing = TRUE),
      organism = 'hsa',
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = pval_cutoff,
      pAdjustMethod = "BH"
    )
    if(!is.null(ekegg_gsea) && nrow(ekegg_gsea@result) > 0) {
      result_df <- as.data.frame(ekegg_gsea)
      # Filter by adjusted p-value
      result_df <- result_df[result_df$p.adjust < pval_cutoff, ]
      if(nrow(result_df) > 0) {
        result_df$MainClass <- sapply(result_df$Description, get_kegg_category)
        cat("Found", nrow(result_df), "significant GSEA pathways\n")
        result_df
      } else {
        cat("No significant GSEA pathways found after multiple testing correction\n")
        NULL
      }
    } else {
      cat("No significant GSEA pathways found\n")
      NULL
    }
  }, error = function(e) {
    warning("Error in KEGG GSEA analysis: ", e$message)
    NULL
  })
  
  # Save results and plots
  cat("\nSaving results and generating plots...\n")
  
  if(!is.null(kegg_up) && nrow(kegg_up) > 0) {
    write.csv(kegg_up, file.path(output_dir, "upregulated_KEGG.csv"), 
              row.names = FALSE)
    save_kegg_plots(kegg_up, "KEGG_up", output_dir)
  } else {
    cat("No significant up-regulated KEGG pathways to save\n")
  }
  
  if(!is.null(kegg_down) && nrow(kegg_down) > 0) {
    write.csv(kegg_down, file.path(output_dir, "downregulated_KEGG.csv"), 
              row.names = FALSE)
    save_kegg_plots(kegg_down, "KEGG_down", output_dir)
  } else {
    cat("No significant down-regulated KEGG pathways to save\n")
  }
  
  if(!is.null(kegg_gsea) && nrow(kegg_gsea) > 0) {
    write.csv(kegg_gsea, file.path(output_dir, "GSEA_KEGG.csv"), 
              row.names = FALSE)
    save_kegg_plots(kegg_gsea, "GSEA", output_dir)
  } else {
    cat("No significant GSEA results to save\n")
  }
  
  # Create summary file
  summary_data <- data.frame(
    Category = c("Total Genes", "Upregulated", "Downregulated",
                 "KEGG Pathways (Up)", "KEGG Pathways (Down)", "GSEA Pathways"),
    Count = c(nrow(data),
              length(gene_lists$up_entrez),
              length(gene_lists$down_entrez),
              ifelse(!is.null(kegg_up), nrow(kegg_up), 0),
              ifelse(!is.null(kegg_down), nrow(kegg_down), 0),
              ifelse(!is.null(kegg_gsea), nrow(kegg_gsea), 0))
  )
  write.csv(summary_data, file.path(output_dir, "analysis_summary.csv"), 
            row.names = FALSE)
  
  # Print final summary
  cat("\nKEGG Analysis Summary:\n")
  cat("Up-regulated pathways:", ifelse(!is.null(kegg_up), nrow(kegg_up), 0), "\n")
  cat("Down-regulated pathways:", ifelse(!is.null(kegg_down), nrow(kegg_down), 0), "\n")
  cat("GSEA pathways:", ifelse(!is.null(kegg_gsea), nrow(kegg_gsea), 0), "\n")
  
  return(list(
    upregulated = kegg_up,
    downregulated = kegg_down,
    gsea = kegg_gsea
  ))
}

# =================== 11. 细胞类型间差异分析与富集分析集成 ===================
# 函数：进行差异分析并直接执行富集分析
run_enrichment_analysis <- function(seurat_obj, ident1, ident2, group.by = "Annotation", 
                                    output_prefix = NULL, log2fc_cutoff = 1, pval_cutoff = 0.05,
                                    test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25) {
  # 保存原始标识
  original_ident <- Idents(seurat_obj)
  
  # 确保组织列是因子类型
  if(!is.factor(seurat_obj@meta.data[[group.by]])) {
    seurat_obj@meta.data[[group.by]] <- as.factor(seurat_obj@meta.data[[group.by]])
  }
  
  # 设置标识
  Idents(seurat_obj) <- seurat_obj@meta.data[[group.by]]
  # 检查身份是否存在
  all_idents <- levels(Idents(seurat_obj))
  if(!(ident1 %in% all_idents) && !(ident1 %in% unique(Idents(seurat_obj)))) {
    missing <- ident1
    stop("The identity '", missing, "' is not present in the Seurat object")
  }
  if(!(ident2 %in% all_idents) && !(ident2 %in% unique(Idents(seurat_obj)))) {
    missing <- ident2
    stop("The identity '", missing, "' is not present in the Seurat object")
  }
  
  # 确定输出前缀
  if(is.null(output_prefix)) {
    output_prefix <- paste0(make.names(ident1), "_vs_", make.names(ident2), "_", 
                           format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
  }
  
  # 运行差异分析但不设置Idents
  print(paste("Running DE analysis between", ident1, "and", ident2))
  de_results <- FindMarkers(
    seurat_obj,
    group.by = group.by,
    ident.1 = ident1,
    ident.2 = ident2,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = test.use
  )
  
  # 准备结果为富集分析所需格式
  de_results$gene_names <- rownames(de_results)
  colnames(de_results)[c(2, 5)] <- c("logFC", "p_val_adj")
  de_results$AveExpr <- (de_results$pct.1 + de_results$pct.2) / 2
  
  # 保存DE结果
  csv_file <- paste0(output_prefix, ".csv")
  write.csv(de_results, csv_file, row.names = FALSE)
  
  # 运行GO分析
  print("Running GO enrichment analysis...")
  go_results <- perform_go_analysis(
    de_results, 
    log2fc_cutoff = log2fc_cutoff, 
    pval_cutoff = pval_cutoff,
    output_prefix = output_prefix
  )
  
  # 运行KEGG分析
  print("Running KEGG pathway analysis...")
  kegg_results <- perform_kegg_analysis(
    de_results, 
    log2fc_cutoff = log2fc_cutoff, 
    pval_cutoff = pval_cutoff,
    output_prefix = output_prefix
  )
  
  # 恢复原始身份
  Idents(seurat_obj) <- original_ident
  
  return(list(
    de_results = de_results,
    go_results = go_results,
    kegg_results = kegg_results
  ))
}

# 函数：运行多组比较
run_multiple_comparisons <- function(seurat_obj, group_list, reference_group = NULL, 
                                    group.by = "Annotation", 
                                    output_prefix = "comparison", 
                                    log2fc_cutoff = 1, pval_cutoff = 0.05) {
  # 检查输入
  if(!is.null(reference_group) && !(reference_group %in% group_list)) {
    stop("Reference group must be in the group list")
  }
  
  results_list <- list()
  
  if(is.null(reference_group)) {
    # 所有组之间两两比较
    comparisons <- combn(group_list, 2, simplify = FALSE)
    
    for(i in seq_along(comparisons)) {
      group_pair <- comparisons[[i]]
      group1 <- group_pair[1]
      group2 <- group_pair[2]
      
      # 生成输出前缀
      pair_prefix <- paste0(output_prefix, "_", make.names(group1), "_vs_", make.names(group2))
      
      # 运行比较
      cat(sprintf("\nRunning comparison %d/%d: %s vs %s\n", 
                i, length(comparisons), group1, group2))
      
      results <- run_enrichment_analysis(
        seurat_obj, 
        ident1 = group1, 
        ident2 = group2, 
        group.by = group.by,
        output_prefix = pair_prefix,
        log2fc_cutoff = log2fc_cutoff, 
        pval_cutoff = pval_cutoff
      )
      
      results_list[[paste(group1, "vs", group2)]] <- results
    }
  } else {
    # 所有组与参考组比较
    other_groups <- setdiff(group_list, reference_group)
    
    for(i in seq_along(other_groups)) {
      group <- other_groups[i]
      
      # 生成输出前缀
      pair_prefix <- paste0(output_prefix, "_", make.names(group), "_vs_", make.names(reference_group))
      
      # 运行比较
      cat(sprintf("\nRunning comparison %d/%d: %s vs %s (reference)\n", 
                i, length(other_groups), group, reference_group))
      
      results <- run_enrichment_analysis(
        seurat_obj, 
        ident1 = group, 
        ident2 = reference_group, 
        group.by = group.by,
        output_prefix = pair_prefix,
        log2fc_cutoff = log2fc_cutoff, 
        pval_cutoff = pval_cutoff
      )
      
      results_list[[paste(group, "vs", reference_group)]] <- results
    }
  }
  
  return(results_list)
}

# =================== 11. 细胞类型间差异分析与富集分析集成 ===================
# 函数：进行差异分析并直接执行富集分析
run_enrichment_analysis <- function(seurat_obj, ident1, ident2, group.by = "Annotation", 
                                    output_prefix = NULL, log2fc_cutoff = 1, pval_cutoff = 0.05,
                                    test.use = "MAST", min.pct = 0.1, logfc.threshold = 0.25) {
  # 记录原始标识
  original_ident <- Idents(seurat_obj)
  
  # 确定输出前缀
  if(is.null(output_prefix)) {
    output_prefix <- paste0(make.names(ident1), "_vs_", make.names(ident2), "_", 
                            format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
  }
  
  # 直接运行差异分析，使用group.by参数
  print(paste("Running DE analysis between", ident1, "and", ident2))
  de_results <- FindMarkers(
    seurat_obj,
    group.by = group.by,
    ident.1 = ident1,
    ident.2 = ident2,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = test.use
  )
  
  # 准备结果为富集分析所需格式
  de_results$gene_names <- rownames(de_results)
  colnames(de_results)[c(2, 5)] <- c("logFC", "p_val_adj")
  de_results$AveExpr <- (de_results$pct.1 + de_results$pct.2) / 2
  
  # 保存DE结果
  csv_file <- paste0(output_prefix, ".csv")
  write.csv(de_results, csv_file, row.names = FALSE)
  
  # 运行GO分析
  print("Running GO enrichment analysis...")
  go_results <- perform_go_analysis(
    de_results, 
    log2fc_cutoff = log2fc_cutoff, 
    pval_cutoff = pval_cutoff,
    output_prefix = output_prefix
  )
  
  # 运行KEGG分析
  print("Running KEGG pathway analysis...")
  kegg_results <- perform_kegg_analysis(
    de_results, 
    log2fc_cutoff = log2fc_cutoff, 
    pval_cutoff = pval_cutoff,
    output_prefix = output_prefix
  )
  
  # 恢复原始标识
  Idents(seurat_obj) <- original_ident
  
  return(list(
    de_results = de_results,
    go_results = go_results,
    kegg_results = kegg_results
  ))
}
# 函数：运行多组比较
run_multiple_comparisons <- function(seurat_obj, group_list, reference_group = NULL, 
                                    group.by = "Annotation", 
                                    output_prefix = "comparison", 
                                    log2fc_cutoff = 1, pval_cutoff = 0.05) {
  # 检查输入
  if(!is.null(reference_group) && !(reference_group %in% group_list)) {
    stop("Reference group must be in the group list")
  }
  
  results_list <- list()
  
  if(is.null(reference_group)) {
    # 所有组之间两两比较
    comparisons <- combn(group_list, 2, simplify = FALSE)
    
    for(i in seq_along(comparisons)) {
      group_pair <- comparisons[[i]]
      group1 <- group_pair[1]
      group2 <- group_pair[2]
      
      # 生成输出前缀
      pair_prefix <- paste0(output_prefix, "_", make.names(group1), "_vs_", make.names(group2))
      
      # 运行比较
      cat(sprintf("\nRunning comparison %d/%d: %s vs %s\n", 
                i, length(comparisons), group1, group2))
      
      results <- run_enrichment_analysis(
        seurat_obj, 
        ident1 = group1, 
        ident2 = group2, 
        group.by = group.by,
        output_prefix = pair_prefix,
        log2fc_cutoff = log2fc_cutoff, 
        pval_cutoff = pval_cutoff
      )
      
      results_list[[paste(group1, "vs", group2)]] <- results
    }
  } else {
    # 所有组与参考组比较
    other_groups <- setdiff(group_list, reference_group)
    
    for(i in seq_along(other_groups)) {
      group <- other_groups[i]
      
      # 生成输出前缀
      pair_prefix <- paste0(output_prefix, "_", make.names(group), "_vs_", make.names(reference_group))
      
      # 运行比较
      cat(sprintf("\nRunning comparison %d/%d: %s vs %s (reference)\n", 
                i, length(other_groups), group, reference_group))
      
      results <- run_enrichment_analysis(
        seurat_obj, 
        ident1 = group, 
        ident2 = reference_group, 
        group.by = group.by,
        output_prefix = pair_prefix,
        log2fc_cutoff = log2fc_cutoff, 
        pval_cutoff = pval_cutoff
      )
      
      results_list[[paste(group, "vs", reference_group)]] <- results
    }
  }
  
  return(results_list)
}

# 函数：对所有组织两两比较特定细胞类型
run_all_tissue_comparisons <- function(seurat_obj, cell_type, 
                                      cell_type_col = "Annotation",
                                      tissue_col = "tissue",
                                      output_dir = NULL,
                                      log2fc_cutoff = 1, 
                                      pval_cutoff = 0.05,
                                      min_cells = 50) {
  # 创建输出目录
  if(is.null(output_dir)) {
    output_dir <- file.path(getwd(), paste0(cell_type, "_tissue_comparisons"))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 提取指定细胞类型的细胞
  cat(paste0("Extracting ", cell_type, " cells...\n"))
  # 使用seurat_obj的元数据直接比较
  cells <- WhichCells(seurat_obj, cells = rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cell_type_col]] == cell_type])
  
  if(length(cells) == 0) {
    stop(paste0("No cells found with ", cell_type_col, " = '", cell_type, "'"))
  }
  
  # 创建子集
  subset_obj <- subset(seurat_obj, cells = cells)
  cat(paste0("Created subset with ", ncol(subset_obj), " cells\n"))
  
  # 获取所有组织类型
  tissues <- unique(subset_obj@meta.data[[tissue_col]])
  tissues <- tissues[!is.na(tissues)]  # 移除NA
  cat(paste0("Found ", length(tissues), " tissues: ", paste(tissues, collapse = ", "), "\n"))
  
  # 检查每个组织的细胞数量
  tissue_counts <- table(subset_obj@meta.data[[tissue_col]])
  valid_tissues <- names(tissue_counts)[tissue_counts >= min_cells]
  
  if(length(tissues) != length(valid_tissues)) {
    excluded <- setdiff(tissues, valid_tissues)
    cat(paste0("Excluding tissues with fewer than ", min_cells, " cells: ", 
               paste(excluded, collapse = ", "), "\n"))
    tissues <- valid_tissues
  }
  
  if(length(tissues) < 2) {
    stop("Need at least 2 tissues for comparison")
  }
  
  # 创建所有可能的组织对
  tissue_pairs <- combn(tissues, 2, simplify = FALSE)
  cat(paste0("Will perform ", length(tissue_pairs), " pairwise comparisons\n"))
  
  # 初始化结果列表
  all_results <- list()
  
  # 对每对组织进行比较
  for(i in seq_along(tissue_pairs)) {
    pair <- tissue_pairs[[i]]
    tissue1 <- pair[1]
    tissue2 <- pair[2]
    
    comparison_name <- paste0(tissue1, "_vs_", tissue2)
    cat(paste0("\nPerforming comparison ", i, "/", length(tissue_pairs), ": ", comparison_name, "\n"))
    
    # 输出前缀
    output_prefix <- file.path(output_dir, paste0(cell_type, "_", comparison_name))
    
    # 检查每个组织的细胞数量
    tissue1_cells <- sum(subset_obj@meta.data[[tissue_col]] == tissue1)
    tissue2_cells <- sum(subset_obj@meta.data[[tissue_col]] == tissue2)
    cat(paste0("  ", tissue1, ": ", tissue1_cells, " cells\n"))
    cat(paste0("  ", tissue2, ": ", tissue2_cells, " cells\n"))
    
    # 执行富集分析
    tryCatch({
      results <- run_enrichment_analysis(
        subset_obj, 
        ident1 = tissue1, 
        ident2 = tissue2, 
        group.by = tissue_col,
        output_prefix = output_prefix,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff
      )
      
      all_results[[comparison_name]] <- results
      cat(paste0("  Completed ", comparison_name, " analysis\n"))
      
    }, error = function(e) {
      cat(paste0("  ERROR in ", comparison_name, ": ", e$message, "\n"))
    })
  }
  
  # 创建汇总报告
  summary_file <- file.path(output_dir, "analysis_summary.txt")
  cat(paste0("Creating summary report at ", summary_file, "\n"))
  
  sink(summary_file)
  cat(paste0("# ", cell_type, " Tissue Comparison Analysis\n\n"))
  cat(paste0("Analysis date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
  cat(paste0("Cell type: ", cell_type, "\n"))
  cat(paste0("Total cells: ", ncol(subset_obj), "\n\n"))
  
  cat("## Tissue Counts\n\n")
  for(tissue in tissues) {
    count <- sum(subset_obj@meta.data[[tissue_col]] == tissue)
    cat(paste0("- ", tissue, ": ", count, " cells\n"))
  }
  
  cat("\n## Comparison Results\n\n")
  for(comparison in names(all_results)) {
    result <- all_results[[comparison]]
    de_count <- nrow(result$de_results)
    up_count <- sum(result$de_results$logFC > 0)
    down_count <- sum(result$de_results$logFC < 0)
    
    cat(paste0("### ", comparison, "\n"))
    cat(paste0("- Total DEGs: ", de_count, "\n"))
    cat(paste0("- Up-regulated: ", up_count, "\n"))
    cat(paste0("- Down-regulated: ", down_count, "\n"))
    
    # GO结果汇总
    if(!is.null(result$go_results$upregulated) && nrow(result$go_results$upregulated) > 0) {
      cat(paste0("- GO terms (up): ", nrow(result$go_results$upregulated), "\n"))
    } else {
      cat("- GO terms (up): 0\n")
    }
    
    if(!is.null(result$go_results$downregulated) && nrow(result$go_results$downregulated) > 0) {
      cat(paste0("- GO terms (down): ", nrow(result$go_results$downregulated), "\n"))
    } else {
      cat("- GO terms (down): 0\n")
    }
    
    # KEGG结果汇总
    if(!is.null(result$kegg_results$upregulated) && nrow(result$kegg_results$upregulated) > 0) {
      cat(paste0("- KEGG pathways (up): ", nrow(result$kegg_results$upregulated), "\n"))
    } else {
      cat("- KEGG pathways (up): 0\n")
    }
    
    if(!is.null(result$kegg_results$downregulated) && nrow(result$kegg_results$downregulated) > 0) {
      cat(paste0("- KEGG pathways (down): ", nrow(result$kegg_results$downregulated), "\n"))
    } else {
      cat("- KEGG pathways (down): 0\n")
    }
    
    if(!is.null(result$kegg_results$gsea) && nrow(result$kegg_results$gsea) > 0) {
      cat(paste0("- GSEA pathways: ", nrow(result$kegg_results$gsea), "\n"))
    } else {
      cat("- GSEA pathways: 0\n")
    }
    
    cat("\n")
  }
  
  sink()
  
  cat("\nAnalysis complete! Results saved to ", output_dir, "\n")
  return(all_results)
}

# 批量分析多个细胞类型的组织差异
run_all_celltypes_tissue_comparisons <- function(seurat_obj, 
                                               cell_types, 
                                               cell_type_col = "Annotation",
                                               tissue_col = "tissue",
                                               base_output_dir = "tissue_comparisons",
                                               log2fc_cutoff = 1, 
                                               pval_cutoff = 0.05,
                                               min_cells = 50) {
  # 创建主输出目录
  dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 初始化结果列表
  all_celltype_results <- list()
  
  # 遍历所有细胞类型
  for(cell_type in cell_types) {
    cat(paste0("\n=====================================\n"))
    cat(paste0("Analyzing ", cell_type, "\n"))
    cat(paste0("=====================================\n"))
    
    output_dir <- file.path(base_output_dir, cell_type)
    
    # 运行单个细胞类型的所有组织比较
    tryCatch({
      results <- run_all_tissue_comparisons(
        seurat_obj = seurat_obj,
        cell_type = cell_type,
        cell_type_col = cell_type_col,
        tissue_col = tissue_col,
        output_dir = output_dir,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff,
        min_cells = min_cells
      )
      
      all_celltype_results[[cell_type]] <- results
      cat(paste0("\nCompleted analysis for ", cell_type, "\n"))
      
    }, error = function(e) {
      cat(paste0("\nERROR processing ", cell_type, ": ", e$message, "\n"))
    })
  }
  
  # 创建汇总报告
  summary_file <- file.path(base_output_dir, "master_summary.txt")
  cat(paste0("\nCreating master summary report at ", summary_file, "\n"))
  
  sink(summary_file)
  cat("# Tissue Comparison Analysis - Master Summary\n\n")
  cat(paste0("Analysis date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
  cat(paste0("Cell types analyzed: ", paste(cell_types, collapse = ", "), "\n\n"))
  
  # 汇总每个细胞类型的结果
  for(cell_type in names(all_celltype_results)) {
    cat(paste0("## ", cell_type, "\n\n"))
    comparisons <- names(all_celltype_results[[cell_type]])
    cat(paste0("Completed ", length(comparisons), " tissue comparisons\n\n"))
    
    for(comparison in comparisons) {
      result <- all_celltype_results[[cell_type]][[comparison]]
      de_count <- nrow(result$de_results)
      
      cat(paste0("### ", comparison, "\n"))
      cat(paste0("- DEGs: ", de_count, "\n"))
      
      # Top GO和KEGG结果
      if(!is.null(result$go_results$upregulated) && nrow(result$go_results$upregulated) > 0) {
        top_go <- head(result$go_results$upregulated[order(result$go_results$upregulated$p.adjust),], 3)
        cat("- Top GO terms (up):\n")
        for(i in 1:nrow(top_go)) {
          cat(paste0("  * ", top_go$Description[i], " (p=", signif(top_go$p.adjust[i], 3), ")\n"))
        }
      }
      
      if(!is.null(result$kegg_results$upregulated) && nrow(result$kegg_results$upregulated) > 0) {
        top_kegg <- head(result$kegg_results$upregulated[order(result$kegg_results$upregulated$p.adjust),], 3)
        cat("- Top KEGG pathways (up):\n")
        for(i in 1:nrow(top_kegg)) {
          cat(paste0("  * ", top_kegg$Description[i], " (p=", signif(top_kegg$p.adjust[i], 3), ")\n"))
        }
      }
      
      cat("\n")
    }
    cat("\n")
  }
  
  sink()
  
  cat("\nAll analyses complete! Results saved to ", base_output_dir, "\n")
  return(all_celltype_results)
}

# 既有的差异基因分析CSV直接用于富集分析
run_enrichment_from_csv <- function(csv_file, log2fc_cutoff = 1, pval_cutoff = 0.05) {
  # 运行GO分析
  print("Running GO enrichment analysis...")
  go_results <- perform_go_analysis(csv_file, log2fc_cutoff, pval_cutoff)
  
  # 运行KEGG分析
  print("Running KEGG pathway analysis...")
  kegg_results <- perform_kegg_analysis(csv_file, log2fc_cutoff, pval_cutoff)
  
  return(list(
    go_results = go_results,
    kegg_results = kegg_results
  ))
}

run_specific_tissue_comparison <- function(seurat_obj, 
                                           cell_type, 
                                           tissue1, 
                                           tissue2, 
                                           cell_type_col = "Annotation", 
                                           tissue_col = "tissue", 
                                           output_dir = NULL, 
                                           log2fc_cutoff = 1, 
                                           pval_cutoff = 0.05,
                                           test.use = "MAST",
                                           min.pct = 0.1,
                                           logfc.threshold = 0.25) {
  
  # 创建输出目录
  if(is.null(output_dir)) {
    output_dir <- file.path(getwd(), paste0(cell_type, "_", tissue1, "_vs_", tissue2))
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 提取特定细胞类型的细胞 - 修改这部分以修复错误
  cat(paste0("Extracting ", cell_type, " cells...\n"))
  
  # 使用更安全的方法获取细胞类型列
  cell_type_data <- FetchData(seurat_obj, vars = cell_type_col)
  cells <- rownames(cell_type_data)[cell_type_data[,1] == cell_type]
  
  if(length(cells) == 0) {
    stop(paste0("No cells found with ", cell_type_col, " = '", cell_type, "'"))
  }
  
  # 创建细胞类型子集
  subset_obj <- subset(seurat_obj, cells = cells)
  cat(paste0("Created subset with ", ncol(subset_obj), " cells\n"))
  
  # 其余代码保持不变...
  # 检查指定的组织是否存在
  all_tissues <- unique(subset_obj@meta.data[[tissue_col]])
  all_tissues <- all_tissues[!is.na(all_tissues)]
  
  if(!(tissue1 %in% all_tissues)) {
    stop(paste0("Tissue '", tissue1, "' not found in the dataset. Available tissues: ", 
                paste(all_tissues, collapse = ", ")))
  }
  if(!(tissue2 %in% all_tissues)) {
    stop(paste0("Tissue '", tissue2, "' not found in the dataset. Available tissues: ", 
                paste(all_tissues, collapse = ", ")))
  }
  
  # 检查每个组织的细胞数量
  tissue1_cells <- sum(subset_obj@meta.data[[tissue_col]] == tissue1, na.rm = TRUE)
  tissue2_cells <- sum(subset_obj@meta.data[[tissue_col]] == tissue2, na.rm = TRUE)
  
  cat(paste0("Cells by tissue:\n",
             "  ", tissue1, ": ", tissue1_cells, " cells\n",
             "  ", tissue2, ": ", tissue2_cells, " cells\n"))
  
  if(tissue1_cells == 0) {
    stop(paste0("No ", cell_type, " cells found in tissue '", tissue1, "'"))
  }
  if(tissue2_cells == 0) {
    stop(paste0("No ", cell_type, " cells found in tissue '", tissue2, "'"))
  }
  
  # 确定输出前缀
  output_prefix <- file.path(output_dir, paste0(cell_type, "_", tissue1, "_vs_", tissue2))
  
  # 执行差异分析和富集分析
  cat(paste0("\nPerforming differential expression and enrichment analysis between ",
             tissue1, " and ", tissue2, " for ", cell_type, " cells...\n"))
  
  # 运行富集分析
  results <- run_enrichment_analysis(
    subset_obj, 
    ident1 = tissue1, 
    ident2 = tissue2, 
    group.by = tissue_col,
    output_prefix = output_prefix,
    log2fc_cutoff = log2fc_cutoff,
    pval_cutoff = pval_cutoff,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold
  )
  
  # 创建简单的分析摘要报告
  summary_file <- file.path(output_dir, "analysis_summary.txt")
  cat(paste0("Creating summary report at ", summary_file, "\n"))
  
  sink(summary_file)
  cat(paste0("# ", cell_type, " Cell Comparison: ", tissue1, " vs ", tissue2, "\n\n"))
  cat(paste0("Analysis date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
  cat(paste0("Parameters:\n",
             "- log2fc_cutoff: ", log2fc_cutoff, "\n",
             "- pval_cutoff: ", pval_cutoff, "\n",
             "- test.use: ", test.use, "\n",
             "- min.pct: ", min.pct, "\n",
             "- logfc.threshold: ", logfc.threshold, "\n\n"))
  
  cat(paste0("Cell counts:\n",
             "- ", cell_type, " cells in ", tissue1, ": ", tissue1_cells, "\n",
             "- ", cell_type, " cells in ", tissue2, ": ", tissue2_cells, "\n\n"))
  
  # 差异表达基因统计
  de_count <- nrow(results$de_results)
  up_count <- sum(results$de_results$logFC > 0)
  down_count <- sum(results$de_results$logFC < 0)
  
  cat(paste0("Differential expression results:\n",
             "- Total DEGs: ", de_count, "\n",
             "- Upregulated in ", tissue1, ": ", up_count, "\n",
             "- Downregulated in ", tissue1, " (upregulated in ", tissue2, "): ", down_count, "\n\n"))
  
  # GO富集结果
  if(!is.null(results$go_results$upregulated) && nrow(results$go_results$upregulated) > 0) {
    cat(paste0("GO enrichment for genes upregulated in ", tissue1, ":\n"))
    cat(paste0("- Total GO terms: ", nrow(results$go_results$upregulated), "\n"))
    
    # 显示前5个GO条目
    if(nrow(results$go_results$upregulated) > 0) {
      top_go <- head(results$go_results$upregulated[order(results$go_results$upregulated$p.adjust),], 5)
      cat("- Top GO terms:\n")
      for(i in 1:nrow(top_go)) {
        cat(paste0("  * ", top_go$Description[i], " (", top_go$ONTOLOGY[i], ", p.adj=", 
                   signif(top_go$p.adjust[i], 3), ", ", top_go$Count[i], " genes)\n"))
      }
    }
    cat("\n")
  } else {
    cat(paste0("No significant GO terms found for genes upregulated in ", tissue1, "\n\n"))
  }
  
  if(!is.null(results$go_results$downregulated) && nrow(results$go_results$downregulated) > 0) {
    cat(paste0("GO enrichment for genes upregulated in ", tissue2, ":\n"))
    cat(paste0("- Total GO terms: ", nrow(results$go_results$downregulated), "\n"))
    
    # 显示前5个GO条目
    if(nrow(results$go_results$downregulated) > 0) {
      top_go <- head(results$go_results$downregulated[order(results$go_results$downregulated$p.adjust),], 5)
      cat("- Top GO terms:\n")
      for(i in 1:nrow(top_go)) {
        cat(paste0("  * ", top_go$Description[i], " (", top_go$ONTOLOGY[i], ", p.adj=", 
                   signif(top_go$p.adjust[i], 3), ", ", top_go$Count[i], " genes)\n"))
      }
    }
    cat("\n")
  } else {
    cat(paste0("No significant GO terms found for genes upregulated in ", tissue2, "\n\n"))
  }
  
  # KEGG富集结果
  if(!is.null(results$kegg_results$upregulated) && nrow(results$kegg_results$upregulated) > 0) {
    cat(paste0("KEGG pathways for genes upregulated in ", tissue1, ":\n"))
    cat(paste0("- Total pathways: ", nrow(results$kegg_results$upregulated), "\n"))
    
    # 显示前5个KEGG通路
    if(nrow(results$kegg_results$upregulated) > 0) {
      top_kegg <- head(results$kegg_results$upregulated[order(results$kegg_results$upregulated$p.adjust),], 5)
      cat("- Top pathways:\n")
      for(i in 1:nrow(top_kegg)) {
        cat(paste0("  * ", top_kegg$Description[i], " (p.adj=", 
                   signif(top_kegg$p.adjust[i], 3), ", ", top_kegg$Count[i], " genes)\n"))
      }
    }
    cat("\n")
  } else {
    cat(paste0("No significant KEGG pathways found for genes upregulated in ", tissue1, "\n\n"))
  }
  
  if(!is.null(results$kegg_results$downregulated) && nrow(results$kegg_results$downregulated) > 0) {
    cat(paste0("KEGG pathways for genes upregulated in ", tissue2, ":\n"))
    cat(paste0("- Total pathways: ", nrow(results$kegg_results$downregulated), "\n"))
    
    # 显示前5个KEGG通路
    if(nrow(results$kegg_results$downregulated) > 0) {
      top_kegg <- head(results$kegg_results$downregulated[order(results$kegg_results$downregulated$p.adjust),], 5)
      cat("- Top pathways:\n")
      for(i in 1:nrow(top_kegg)) {
        cat(paste0("  * ", top_kegg$Description[i], " (p.adj=", 
                   signif(top_kegg$p.adjust[i], 3), ", ", top_kegg$Count[i], " genes)\n"))
      }
    }
    cat("\n")
  } else {
    cat(paste0("No significant KEGG pathways found for genes upregulated in ", tissue2, "\n\n"))
  }
  
  # GSEA结果
  if(!is.null(results$kegg_results$gsea) && nrow(results$kegg_results$gsea) > 0) {
    cat("KEGG GSEA results:\n")
    cat(paste0("- Total enriched pathways: ", nrow(results$kegg_results$gsea), "\n"))
    
    # 区分正负富集
    pos_gsea <- results$kegg_results$gsea[results$kegg_results$gsea$NES > 0, ]
    neg_gsea <- results$kegg_results$gsea[results$kegg_results$gsea$NES < 0, ]
    
    cat(paste0("- Pathways enriched in ", tissue1, ": ", nrow(pos_gsea), "\n"))
    cat(paste0("- Pathways enriched in ", tissue2, ": ", nrow(neg_gsea), "\n"))
    
    # 显示前几个GSEA结果
    if(nrow(pos_gsea) > 0) {
      top_pos <- head(pos_gsea[order(pos_gsea$p.adjust),], 5)
      cat(paste0("\nTop pathways enriched in ", tissue1, ":\n"))
      for(i in 1:nrow(top_pos)) {
        cat(paste0("  * ", top_pos$Description[i], " (NES=", 
                   round(top_pos$NES[i], 2), ", p.adj=", 
                   signif(top_pos$p.adjust[i], 3), ")\n"))
      }
    }
    
    if(nrow(neg_gsea) > 0) {
      top_neg <- head(neg_gsea[order(neg_gsea$p.adjust),], 5)
      cat(paste0("\nTop pathways enriched in ", tissue2, ":\n"))
      for(i in 1:nrow(top_neg)) {
        cat(paste0("  * ", top_neg$Description[i], " (NES=", 
                   round(top_neg$NES[i], 2), ", p.adj=", 
                   signif(top_neg$p.adjust[i], 3), ")\n"))
      }
    }
  } else {
    cat("No significant GSEA results found\n")
  }
  
  sink()
  
  # 完成
  cat(paste0("\nAnalysis complete! Results saved to ", output_dir, "\n"))
  return(results)
}

# =================== 12. 使用示例 ===================
# 对seurat_obj_main中"Epithelial"细胞类型的所有组织进行两两比较
if(FALSE) {  # 设为TRUE来执行
  # 加载必要的包
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  
  # 单个细胞类型的所有组织两两比较
  epithelial_results <- run_all_tissue_comparisons(
    seurat_obj = seurat_obj_main,
    cell_type = "Epithelial",  # 细胞类型
    cell_type_col = "Annotation",  # 细胞类型所在的列名
    tissue_col = "tissue",  # 组织信息所在的列名
    output_dir = "Epithelial_tissue_comparisons",
    log2fc_cutoff = 1.5,
    pval_cutoff = 0.05,
    min_cells = 50  # 每个组织至少需要的细胞数
  )
}

# 对多个细胞类型的所有组织进行两两比较
if(FALSE) {  # 设为TRUE来执行
  # 加载必要的包
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  
  # 批量分析多个细胞类型
  all_results <- run_all_celltypes_tissue_comparisons(
    seurat_obj = seurat_obj_main,
    cell_types = c("Epithelial", "Fibroblast", "T", "B", "Myeloid"),
    cell_type_col = "Annotation",
    tissue_col = "tissue",
    base_output_dir = "all_tissue_comparisons",
    log2fc_cutoff = 1,
    pval_cutoff = 0.05,
    min_cells = 50
  )
}

# 从CSV文件直接运行富集分析
if(FALSE) {  # 设为TRUE来执行
  # 加载必要的包
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
  
  csv_file <- "epi_Tissue-Lung_vs_Tissue-Sinus_2025-01-19-02-39-56.csv"
  enrichment_from_csv <- run_enrichment_from_csv(
    csv_file, 
    log2fc_cutoff = 2, 
    pval_cutoff = 0.05
  )
}