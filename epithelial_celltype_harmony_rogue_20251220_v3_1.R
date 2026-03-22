# ===== Per-Cell-Type Harmony Integration and ROGUE Analysis =====
# Version: v1.2.1 (Critical Fixes from Second Code Review)
# Author: r2end
# Date: 2024-12-20
#
# Purpose: 
# - Subset each cell type from Manual_Annotation
# - Perform Harmony batch correction (or PCA if single batch)
# - Cluster at appropriate resolution
# - Calculate cluster-level ROGUE purity scores (overall + per-sample)
#
# Workflow:
# 1. Load full epithelial Seurat object
# 2. For each cell type in Manual_Annotation:
#    - Subset cells (P0-1 FIX: explicit cell selection)
#    - Set DefaultAssay to "RNA" (P1-1 FIX)
#    - Normalize, find variable features, scale (P1-1 FIX: only HVG)
#      * IMPORTANT: This does NOT delete genes/cells, only feature selection
#    - PCA -> Harmony (if multi-batch) -> UMAP -> Clustering (P1-2 FIX: res=1.5)
#      * PCA uses HVG (feature selection), does NOT modify object
#    - Calculate ROGUE for each cluster:
#      * Overall cluster ROGUE (all samples combined)
#      * Per-sample ROGUE (each sample within cluster separately)
#      * ROGUE calculation uses matr.filter (removes low-quality cells/genes)
#        but this happens ONLY during computation, NOT to Seurat object
#
# New in v1.2.1 (Second Review Fixes):
# - P0-1: Encapsulated ROGUE calculation in calc_rogue_safe() with sparse->dense fallback
# - P0-2: Fixed per-sample boxplot statistical validity (points + cluster boxplot)
# - P0-3: Added BATCH_COL validation and NA sample filtering
# - P1-1: Explicitly set DefaultAssay to "RNA" after subset
# - P1-2: Track filtering ratios to assess ROGUE quality impact
# - P1-3: Unified ROGUE thresholds (overall & per-sample use same values)
# - Added extensive comments clarifying feature selection vs deletion
#
# New in v1.2:
# - Per-sample ROGUE calculation for batch-specific quality assessment
# - Identify problematic sample×cluster combinations
# - Compare overall vs per-sample ROGUE to detect batch heterogeneity
# - Additional visualizations: per-sample distributions, sample variance
#
# Key Fixes from First Review (v1.1):
# - P0-1: Fixed subset() NSE issue - use explicit cell selection
# - P0-2: Fixed UMAP plotting - use gridExtra instead of patchwork
# - P0-3: Fixed ROGUE OOM - extract counts once, downsample, try sparse first
# - P0-4: Fixed heatmap cluster collision - use cluster_rank instead of raw IDs
# - P1-1: ScaleData only on HVG (saves time/memory)
# - P1-2: Reduced resolution from 3 to 1.5 (appropriate for already-subset data)
# - P1-3: Added Harmony parameter compatibility warnings
# - P1-4: Added biological interpretation warnings
#
# Critical Understanding:
# - NO genes/cells are removed from Seurat object before or during Harmony
# - HVG selection is for PCA/Harmony computation, not deletion
# - ROGUE filtering (matr.filter) is ONLY for ROGUE calculation
# - All count data remains intact in the object throughout
#
# Dependencies: Seurat, harmony, ROGUE, tidyverse, gridExtra, pheatmap

# ===== Load Libraries =====
cat("Loading required libraries...\n")
required_packages <- c("Seurat", "harmony", "ROGUE", "tidyverse", "gridExtra", "pheatmap")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed. Install with: install.packages('%s')", pkg, pkg))
  }
}

suppressMessages({
  library(Seurat)
  library(harmony)
  library(ROGUE)
  library(tidyverse)
  library(gridExtra)
  library(pheatmap)
})

cat("✓ All libraries loaded successfully\n\n")

# ===== Configuration Parameters =====
INPUT_RDS <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "per_celltype_harmony_rogue"
ANNOTATION_COL <- "Manual_Annotation"  # Column to subset by
BATCH_COL <- "dataset"                   # Batch correction variable

# Harmony parameters
HARMONY_THETA <- 1
HARMONY_LAMBDA <- 7
HARMONY_SIGMA <- 0.1
HARMONY_NCLUST <- 30
HARMONY_MAX_ITER <- 20
HARMONY_DIMS <- 1:50

# Clustering parameters
N_PCS <- 50
UMAP_DIMS <- 1:30
UMAP_NEIGHBORS <- 30
NEIGHBOR_K <- 45
CLUSTER_RESOLUTION <- 1.5  # Reduced from 3 - within single cell type, lower res avoids over-fragmentation

# ROGUE parameters
MIN_CELLS_PER_CELLTYPE <- 100  # Minimum cells to process a cell type
MIN_CELLS_PER_CLUSTER <- 10    # Minimum cells for ROGUE calculation
MAX_CELLS_ROGUE <- 2000         # Downsample large clusters to prevent OOM (P0-3 fix)

# ROGUE calculation parameters (P1-2 & P1-3 fix: unified thresholds)
ROGUE_MIN_CELLS <- 10           # Minimum cells for matr.filter (same for overall & per-sample)
ROGUE_MIN_GENES <- 200          # Minimum genes for matr.filter (same for overall & per-sample)
ROGUE_FILTER_WARNING_THRESHOLD <- 0.7  # Warn if filtered_fraction < this value

# ===== Setup Output Directory =====
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "seurat_objects"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "rogue_results"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE)

# Helper function: Safe filename generation
safe_filename <- function(name) {
  # Replace spaces, slashes, and other problematic characters
  name <- gsub("[\\s/\\\\:*?\"<>|]", "_", name)
  # Remove leading/trailing underscores
  name <- gsub("^_+|_+$", "", name)
  # Collapse multiple underscores
  name <- gsub("_+", "_", name)
  return(name)
}

# Helper function: Safe ROGUE calculation with sparse->dense fallback
# Handles the entire chain: matr.filter -> SE_fun -> CalculateRogue
calc_rogue_safe <- function(mat, min.cells = 10, min.genes = 200, platform = "UMI") {
  
  # Step 1: Try sparse path (entire chain)
  sparse_result <- try({
    # Filter
    mat_filt <- matr.filter(mat, min.cells = min.cells, min.genes = min.genes)
    
    # Check basic validity
    if (ncol(mat_filt) < min.cells || nrow(mat_filt) < min.genes/2) {
      stop("Insufficient data after filtering")
    }
    
    # Calculate entropy
    ent <- SE_fun(mat_filt)
    
    # Validate entropy
    if (any(is.na(ent$entropy)) || any(is.infinite(ent$entropy))) {
      stop("Invalid entropy values")
    }
    
    # Calculate ROGUE
    rogue_val <- CalculateRogue(ent, platform = platform)
    
    # Validate ROGUE
    if (is.na(rogue_val) || is.infinite(rogue_val) || 
        rogue_val < 0 || rogue_val > 1) {
      stop(sprintf("Invalid ROGUE value: %.4f", rogue_val))
    }
    
    # Return results
    list(
      rogue_value = rogue_val,
      n_genes = nrow(mat_filt),
      n_cells = ncol(mat_filt),
      method = "sparse"
    )
  }, silent = TRUE)
  
  # If sparse succeeded, return
  if (!inherits(sparse_result, "try-error")) {
    return(sparse_result)
  }
  
  # Step 2: Dense fallback (entire chain)
  dense_result <- try({
    # Densify
    mat_dense <- as.matrix(mat)
    
    # Filter
    mat_filt <- matr.filter(mat_dense, min.cells = min.cells, min.genes = min.genes)
    
    # Check basic validity
    if (ncol(mat_filt) < min.cells || nrow(mat_filt) < min.genes/2) {
      stop("Insufficient data after filtering")
    }
    
    # Calculate entropy
    ent <- SE_fun(mat_filt)
    
    # Validate entropy
    if (any(is.na(ent$entropy)) || any(is.infinite(ent$entropy))) {
      stop("Invalid entropy values")
    }
    
    # Calculate ROGUE
    rogue_val <- CalculateRogue(ent, platform = platform)
    
    # Validate ROGUE
    if (is.na(rogue_val) || is.infinite(rogue_val) || 
        rogue_val < 0 || rogue_val > 1) {
      stop(sprintf("Invalid ROGUE value: %.4f", rogue_val))
    }
    
    # Return results
    list(
      rogue_value = rogue_val,
      n_genes = nrow(mat_filt),
      n_cells = ncol(mat_filt),
      method = "dense"
    )
  }, silent = TRUE)
  
  # If dense also failed, return NULL
  if (inherits(dense_result, "try-error")) {
    return(NULL)
  }
  
  return(dense_result)
}

cat("=== Per-Cell-Type Harmony Integration and ROGUE Analysis ===\n")
cat(sprintf("Input: %s\n", INPUT_RDS))
cat(sprintf("Output: %s\n\n", OUTPUT_DIR))

# ===== Load Data =====
cat("Loading Seurat object...\n")

# Check if input file exists
if (!file.exists(INPUT_RDS)) {
  stop(sprintf("Input file not found: %s", INPUT_RDS))
}

# Load and validate
seurat_obj <- readRDS(INPUT_RDS)

# Validate object structure
if (!inherits(seurat_obj, "Seurat")) {
  stop("Input file is not a Seurat object")
}

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("Seurat object does not contain 'RNA' assay")
}

if (!"counts" %in% names(seurat_obj@assays$RNA@layers)) {
  stop("RNA assay does not contain 'counts' layer")
}

cat("✓ Seurat object loaded successfully\n")
cat(sprintf("  Total cells: %s\n", format(ncol(seurat_obj), big.mark = ",")))
cat(sprintf("  Total genes: %s\n", format(nrow(seurat_obj), big.mark = ",")))
cat("\n")

# Check annotation column
if (!ANNOTATION_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Column '%s' not found in metadata", ANNOTATION_COL))
}

# Check batch column (P0-3 fix)
if (!BATCH_COL %in% colnames(seurat_obj@meta.data)) {
  stop(sprintf("Batch column '%s' not found in metadata", BATCH_COL))
}

if (all(is.na(seurat_obj@meta.data[[BATCH_COL]]))) {
  stop(sprintf("Batch column '%s' contains only NA values", BATCH_COL))
}

n_na_batch <- sum(is.na(seurat_obj@meta.data[[BATCH_COL]]))
if (n_na_batch > 0) {
  cat(sprintf("⚠ WARNING: %d cells have NA in batch column '%s' (will be excluded from per-sample analysis)\n", 
              n_na_batch, BATCH_COL))
}

# Get cell types
cell_types <- unique(seurat_obj@meta.data[[ANNOTATION_COL]])
cell_types <- cell_types[!is.na(cell_types)]
cell_type_counts <- table(seurat_obj@meta.data[[ANNOTATION_COL]])

cat(sprintf("Found %d cell types:\n", length(cell_types)))
print(cell_type_counts)
cat("\n")

# Filter cell types by minimum cell count
cell_types_filtered <- names(cell_type_counts[cell_type_counts >= MIN_CELLS_PER_CELLTYPE])
cat(sprintf("Processing %d cell types with >= %d cells\n\n", 
            length(cell_types_filtered), MIN_CELLS_PER_CELLTYPE))

# ===== Initialize Results Storage =====
all_rogue_results <- list()
all_rogue_per_sample_results <- list()
processing_summary <- list()

# ===== Process Each Cell Type =====
for (i in seq_along(cell_types_filtered)) {
  ct <- cell_types_filtered[i]
  start_time <- Sys.time()
  
  cat(sprintf("\n===== [%d/%d] Processing: %s =====\n", i, length(cell_types_filtered), ct))
  n_cells <- cell_type_counts[ct]
  cat(sprintf("Total cells: %s\n", format(n_cells, big.mark = ",")))
  
  tryCatch({
    
    # ===== 1. Subset Data =====
    cat("Step 1/7: Subsetting data...\n")
    
    # Use explicit cell selection to avoid NSE issues with Seurat::subset()
    cells_ct <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[ANNOTATION_COL]] == ct]
    if (length(cells_ct) == 0) {
      stop(sprintf("No cells found for cell type: %s", ct))
    }
    seurat_sub <- subset(seurat_obj, cells = cells_ct)
    
    # P1-1 FIX: Explicitly set default assay to RNA
    DefaultAssay(seurat_sub) <- "RNA"
    
    # Check batch distribution
    batch_counts <- table(seurat_sub@meta.data[[BATCH_COL]])
    n_batches <- length(batch_counts)
    cat(sprintf("  Batches: %d (range: %d-%d cells)\n", 
                n_batches, min(batch_counts), max(batch_counts)))
    
    use_harmony <- n_batches >= 2
    if (!use_harmony) {
      cat("  INFO: Only 1 batch, will skip Harmony and use PCA directly\n")
    }
    
    # ===== 2. Normalization and Feature Selection =====
    # IMPORTANT: This step does NOT remove genes/cells from the object
    # - FindVariableFeatures only MARKS HVG, does not delete non-HVG
    # - ScaleData only creates scaled values for HVG, does not delete genes
    # - The full gene count matrix remains intact for downstream use
    cat("Step 2/7: Normalization and feature selection...\n")
    seurat_sub <- NormalizeData(seurat_sub, verbose = FALSE)
    seurat_sub <- FindVariableFeatures(
      seurat_sub,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    )
    # Scale only HVG (saves time and memory)
    seurat_sub <- ScaleData(
      seurat_sub, 
      features = VariableFeatures(seurat_sub), 
      verbose = FALSE
    )
    
    # ===== 3. PCA =====
    # PCA uses HVG by default (feature selection, not deletion)
    cat("Step 3/7: PCA...\n")
    seurat_sub <- RunPCA(
      seurat_sub,
      npcs = N_PCS,
      verbose = FALSE
    )
    
    # ===== 4. Harmony Batch Correction (if multiple batches) =====
    reduction_use <- "pca"
    if (use_harmony) {
      cat("Step 4/7: Harmony integration...\n")
      
      # WARNING (P1-3): Harmony parameter names may vary across versions
      # If RunHarmony fails, check for: reduction.use vs reduction, dims.use vs dims
      # Current params tested on harmony v1.0+
      
      # WARNING (P1-4): Correcting by 'sample' may remove real biological differences
      # if disease/condition is highly correlated with sample
      # Consider comparing pre/post Harmony UMAP by disease status
      
      seurat_sub <- RunHarmony(
        object = seurat_sub,
        group.by.vars = BATCH_COL,
        theta = HARMONY_THETA,
        lambda = HARMONY_LAMBDA,
        sigma = HARMONY_SIGMA,
        nclust = HARMONY_NCLUST,
        reduction.use = "pca",
        max_iter = HARMONY_MAX_ITER,
        early_stop = TRUE,
        dims = HARMONY_DIMS,
        verbose = FALSE
      )
      
      # Verify Harmony reduction was created
      if ("harmony" %in% names(seurat_sub@reductions)) {
        reduction_use <- "harmony"
        cat("  ✓ Harmony reduction created successfully\n")
      } else {
        cat("  ⚠ Harmony failed, falling back to PCA\n")
        reduction_use <- "pca"
      }
    } else {
      cat("Step 4/7: Skipping Harmony (single batch)...\n")
    }
    
    # Verify the reduction exists
    if (!reduction_use %in% names(seurat_sub@reductions)) {
      stop(sprintf("Reduction '%s' not found in Seurat object", reduction_use))
    }
    
    # ===== 5. UMAP =====
    cat("Step 5/7: UMAP...\n")
    seurat_sub <- RunUMAP(
      seurat_sub,
      reduction = reduction_use,
      dims = UMAP_DIMS,
      n.neighbors = UMAP_NEIGHBORS,
      n.trees = 500,
      min.dist = 0.4,
      metric = "correlation",
      verbose = FALSE
    )
    
    # ===== 6. Clustering =====
    cat("Step 6/7: Finding neighbors and clustering...\n")
    seurat_sub <- FindNeighbors(
      seurat_sub,
      reduction = reduction_use,
      dims = UMAP_DIMS,
      k.param = NEIGHBOR_K,
      verbose = FALSE
    )
    
    seurat_sub <- FindClusters(
      seurat_sub,
      algorithm = 4,
      group.singletons = FALSE,
      resolution = CLUSTER_RESOLUTION,
      verbose = FALSE
    )
    
    # Get clustering column name
    cluster_cols <- grep("^RNA_snn_res", colnames(seurat_sub@meta.data), value = TRUE)
    if (length(cluster_cols) == 0) {
      stop("No clustering results found. Check FindClusters output.")
    }
    cluster_col <- cluster_cols[length(cluster_cols)]  # Use highest resolution
    n_clusters <- length(unique(seurat_sub@meta.data[[cluster_col]]))
    cat(sprintf("  Found %d clusters at resolution %.1f\n", n_clusters, CLUSTER_RESOLUTION))
    
    # ===== 7. Calculate ROGUE for Each Cluster =====
    # IMPORTANT: ROGUE calculation uses matr.filter which removes low-quality cells/genes
    # This happens ONLY during ROGUE computation and does NOT modify the Seurat object
    # The filtering can systematically raise ROGUE values by removing poor-quality cells
    cat("Step 7/7: Calculating ROGUE for each cluster (overall + per-sample)...\n")
    
    # Extract counts ONCE for all clusters (avoid repeated large matrix extraction)
    expr_all <- GetAssayData(seurat_sub, layer = "counts", assay = "RNA")
    
    cluster_rogue <- list()
    cluster_rogue_per_sample <- list()
    clusters <- sort(unique(seurat_sub@meta.data[[cluster_col]]))
    
    for (clust in clusters) {
      
      # Get cells for this cluster
      cells_keep <- rownames(seurat_sub@meta.data)[seurat_sub@meta.data[[cluster_col]] == clust]
      n_cells_clust <- length(cells_keep)
      
      cat(sprintf("  Cluster %s: %d cells... ", clust, n_cells_clust))
      
      if (n_cells_clust < MIN_CELLS_PER_CLUSTER) {
        cat("SKIPPED (too few cells)\n")
        next
      }
      
      # ===== Calculate Overall Cluster ROGUE =====
      overall_rogue <- tryCatch({
        # Downsample if too many cells (prevent memory explosion)
        cells_sampled <- cells_keep
        if (n_cells_clust > MAX_CELLS_ROGUE) {
          set.seed(42)
          cells_sampled <- sample(cells_keep, MAX_CELLS_ROGUE)
          cat(sprintf("(downsampled to %d) ", MAX_CELLS_ROGUE))
        }
        
        # Extract expression matrix for this cluster
        expr_subset <- expr_all[, cells_sampled, drop = FALSE]
        
        # Check if matrix is valid
        if (ncol(expr_subset) == 0 || nrow(expr_subset) == 0) {
          cat("SKIPPED (empty matrix)\n")
          NULL  # Changed from return(NULL)
        } else {
          n_cells_before <- ncol(expr_subset)
          n_genes_before <- nrow(expr_subset)
          
          # Use calc_rogue_safe for robust calculation (P0-1 fix)
          rogue_result <- calc_rogue_safe(
            expr_subset, 
            min.cells = ROGUE_MIN_CELLS, 
            min.genes = ROGUE_MIN_GENES, 
            platform = "UMI"
          )
          
          if (is.null(rogue_result)) {
            cat("SKIPPED (calc_rogue_safe failed)\n")
            NULL  # Changed from return(NULL)
          } else {
            # Calculate filtering ratio (P1-2 fix: track quality impact)
            n_cells_after <- rogue_result$n_cells
            n_genes_after <- rogue_result$n_genes
            filter_fraction <- n_cells_after / n_cells_before
            
            # Warn if too many cells filtered
            warning_flag <- ""
            if (filter_fraction < ROGUE_FILTER_WARNING_THRESHOLD) {
              warning_flag <- sprintf(" ⚠LOW_RETENTION(%.1f%%)", filter_fraction * 100)
            }
            
            cat(sprintf("Overall ROGUE = %.4f (cells:%d->%d, genes:%d->%d, method:%s)%s ", 
                       rogue_result$rogue_value, n_cells_before, n_cells_after, 
                       n_genes_before, n_genes_after, rogue_result$method, warning_flag))
            
            # Store overall results
            result <- data.frame(
              cell_type = ct,
              cluster = clust,
              sample = "ALL",  # Indicator for overall
              n_cells = n_cells_clust,
              n_cells_used = length(cells_sampled),
              n_cells_filtered = n_cells_after,
              n_genes_used = n_genes_after,
              filter_fraction = filter_fraction,
              rogue_value = rogue_result$rogue_value,
              calc_method = rogue_result$method,
              stringsAsFactors = FALSE
            )
            
            result  # Changed from return(result) - tryCatch returns last expression
          }
        }
      }, error = function(e) {
        cat(sprintf("ERROR: %s\n", conditionMessage(e)))
        NULL  # Changed from return(NULL)
      })
      
      if (!is.null(overall_rogue)) {
        cluster_rogue[[paste0(clust, "_ALL")]] <- overall_rogue
      }
      
      # ===== Calculate Per-Sample ROGUE =====
      # Get sample distribution for this cluster
      cluster_meta <- seurat_sub@meta.data[cells_keep, , drop = FALSE]
      samples_in_cluster <- unique(cluster_meta[[BATCH_COL]])
      
      # P0-3 fix: Remove NA samples
      samples_in_cluster <- samples_in_cluster[!is.na(samples_in_cluster)]
      n_samples <- length(samples_in_cluster)
      
      cat(sprintf("| Per-sample (%d samples): ", n_samples))
      
      sample_rogue_list <- list()
      for (samp in samples_in_cluster) {
        
        # Get cells for this cluster × sample
        cells_sample <- rownames(cluster_meta)[cluster_meta[[BATCH_COL]] == samp & 
                                                 !is.na(cluster_meta[[BATCH_COL]])]
        n_cells_sample <- length(cells_sample)
        
        # Skip if too few cells
        if (n_cells_sample < MIN_CELLS_PER_CLUSTER) {
          next
        }
        
        # FIX: Don't use return() in for loop - just assign NULL or result
        sample_rogue <- tryCatch({
          # Downsample if needed
          cells_sample_use <- cells_sample
          if (n_cells_sample > MAX_CELLS_ROGUE) {
            set.seed(42)
            cells_sample_use <- sample(cells_sample, MAX_CELLS_ROGUE)
          }
          
          # Extract expression matrix
          expr_sample <- expr_all[, cells_sample_use, drop = FALSE]
          
          if (ncol(expr_sample) == 0 || nrow(expr_sample) == 0) {
            NULL  # Changed from return(NULL)
          } else {
            n_cells_before <- ncol(expr_sample)
            n_genes_before <- nrow(expr_sample)
            
            # P1-3 fix: Use UNIFIED thresholds (same as overall)
            rogue_result <- calc_rogue_safe(
              expr_sample, 
              min.cells = ROGUE_MIN_CELLS, 
              min.genes = ROGUE_MIN_GENES, 
              platform = "UMI"
            )
            
            if (is.null(rogue_result)) {
              NULL  # Changed from return(NULL)
            } else {
              # Calculate filtering ratio
              n_cells_after <- rogue_result$n_cells
              n_genes_after <- rogue_result$n_genes
              filter_fraction <- n_cells_after / n_cells_before
              
              result_sample <- data.frame(
                cell_type = ct,
                cluster = clust,
                sample = samp,
                n_cells = n_cells_sample,
                n_cells_used = length(cells_sample_use),
                n_cells_filtered = n_cells_after,
                n_genes_used = n_genes_after,
                filter_fraction = filter_fraction,
                rogue_value = rogue_result$rogue_value,
                calc_method = rogue_result$method,
                stringsAsFactors = FALSE
              )
              
              result_sample  # Changed from return(result_sample)
            }
          }
        }, error = function(e) {
          NULL  # Changed from return(NULL)
        })
        
        if (!is.null(sample_rogue)) {
          sample_rogue_list[[paste0(clust, "_", samp)]] <- sample_rogue
        }
      }
      
      # Store per-sample results
      if (length(sample_rogue_list) > 0) {
        sample_rogue_df <- bind_rows(sample_rogue_list)
        cluster_rogue_per_sample[[as.character(clust)]] <- sample_rogue_df
        cat(sprintf("%d/%d samples OK\n", nrow(sample_rogue_df), n_samples))
      } else {
        cat("0 samples OK\n")
      }
    }
    
    # Clean up global counts matrix
    rm(expr_all)
    gc(verbose = FALSE)
    
    # Compile cluster-level results
    if (length(cluster_rogue) > 0) {
      ct_rogue_df <- bind_rows(cluster_rogue)
      all_rogue_results[[ct]] <- ct_rogue_df
      
      # Save per-cell-type overall ROGUE results
      safe_ct_name <- safe_filename(ct)
      write.csv(
        ct_rogue_df,
        file.path(OUTPUT_DIR, "rogue_results", sprintf("%s_rogue_overall.csv", safe_ct_name)),
        row.names = FALSE
      )
      
      cat(sprintf("\n  Successfully calculated overall ROGUE for %d/%d clusters\n", 
                  nrow(ct_rogue_df), n_clusters))
    } else {
      cat("\n  WARNING: No valid overall ROGUE values calculated\n")
    }
    
    # Compile and save per-sample results
    if (length(cluster_rogue_per_sample) > 0) {
      ct_rogue_per_sample_df <- bind_rows(cluster_rogue_per_sample)
      all_rogue_per_sample_results[[ct]] <- ct_rogue_per_sample_df
      
      # Save per-cell-type per-sample ROGUE results
      safe_ct_name <- safe_filename(ct)
      write.csv(
        ct_rogue_per_sample_df,
        file.path(OUTPUT_DIR, "rogue_results", sprintf("%s_rogue_per_sample.csv", safe_ct_name)),
        row.names = FALSE
      )
      
      cat(sprintf("  Successfully calculated per-sample ROGUE: %d cluster×sample combinations\n", 
                  nrow(ct_rogue_per_sample_df)))
    } else {
      cat("  WARNING: No valid per-sample ROGUE values calculated\n")
    }
    
    # ===== Save Processed Seurat Object =====
    safe_ct_name <- safe_filename(ct)
    saveRDS(
      seurat_sub,
      file.path(OUTPUT_DIR, "seurat_objects", sprintf("%s_harmony.rds", safe_ct_name))
    )
    
    # ===== Generate Per-Cell-Type Plots =====
    # UMAP plots
    pdf(file.path(OUTPUT_DIR, "plots", sprintf("%s_umap.pdf", safe_ct_name)), 
        width = 14, height = 5)
    p1 <- DimPlot(seurat_sub, reduction = "umap", group.by = BATCH_COL) +
      ggtitle(sprintf("%s - Batch Distribution", ct)) +
      theme(plot.title = element_text(face = "bold"))
    p2 <- DimPlot(seurat_sub, reduction = "umap", group.by = cluster_col, label = TRUE) +
      ggtitle(sprintf("%s - Clusters (res=%.1f)", ct, CLUSTER_RESOLUTION)) +
      theme(plot.title = element_text(face = "bold"))
    # Use gridExtra instead of patchwork
    gridExtra::grid.arrange(p1, p2, nrow = 1)
    dev.off()
    
    # ROGUE plots if available
    if (length(cluster_rogue) > 0) {
      # Overall ROGUE barplot
      pdf(file.path(OUTPUT_DIR, "plots", sprintf("%s_rogue_barplot_overall.pdf", safe_ct_name)), 
          width = 8, height = 6)
      
      ct_rogue_df <- ct_rogue_df %>%
        mutate(
          purity_class = case_when(
            rogue_value >= 0.9 ~ "High (≥0.9)",
            rogue_value >= 0.7 ~ "Moderate (0.7-0.9)",
            TRUE ~ "Low (<0.7)"
          )
        )
      
      p <- ggplot(ct_rogue_df, aes(x = reorder(cluster, rogue_value), y = rogue_value, 
                                   fill = purity_class)) +
        geom_col() +
        geom_text(aes(label = sprintf("%.3f", rogue_value)), 
                  hjust = -0.1, size = 3) +
        geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red", alpha = 0.5) +
        coord_flip() +
        scale_fill_manual(values = c("High (≥0.9)" = "forestgreen",
                                     "Moderate (0.7-0.9)" = "orange",
                                     "Low (<0.7)" = "firebrick")) +
        labs(title = sprintf("Overall ROGUE Cluster Purity - %s", ct),
             x = "Cluster",
             y = "ROGUE Value",
             fill = "Purity Class") +
        theme_bw() +
        theme(plot.title = element_text(face = "bold")) +
        ylim(0, 1.05)
      print(p)
      dev.off()
      
      # Per-sample ROGUE plot if available
      if (length(cluster_rogue_per_sample) > 0) {
        ct_rogue_per_sample_df <- bind_rows(cluster_rogue_per_sample)
        
        # P0-2 fix: Per-sample visualization using points + cluster-level boxplot
        # (Each cluster×sample has only 1 ROGUE value, so boxplot across samples makes sense)
        pdf(file.path(OUTPUT_DIR, "plots", sprintf("%s_rogue_per_sample.pdf", safe_ct_name)), 
            width = 12, height = 8)
        
        p_sample <- ggplot(ct_rogue_per_sample_df, 
                          aes(x = factor(cluster), y = rogue_value)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "grey80") +
          geom_point(aes(color = sample), alpha = 0.7, size = 3, 
                    position = position_jitter(width = 0.15, seed = 42)) +
          geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", 
                    color = "red", alpha = 0.5) +
          labs(title = sprintf("Per-Sample ROGUE Distribution - %s", ct),
               subtitle = "Boxplot shows variation across samples; each point = one sample in cluster",
               x = "Cluster",
               y = "ROGUE Value",
               color = "Sample") +
          theme_bw() +
          theme(plot.title = element_text(face = "bold"),
                legend.position = "right",
                axis.text.x = element_text(angle = 0, hjust = 0.5)) +
          ylim(0, 1.05)
        print(p_sample)
        dev.off()
      }
    }
    
    # Record processing success
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    processing_summary[[ct]] <- list(
      status = "SUCCESS",
      n_cells = n_cells,
      n_batches = n_batches,
      n_clusters = n_clusters,
      n_rogue_calculated = ifelse(length(cluster_rogue) > 0, nrow(bind_rows(cluster_rogue)), 0),
      reduction_used = reduction_use,
      time_seconds = elapsed_time
    )
    
    cat(sprintf("✓ %s processing complete (reduction: %s, time: %.1fs)\n", ct, reduction_use, elapsed_time))
    
    # Clean up
    rm(seurat_sub)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    end_time <- Sys.time()
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    cat(sprintf("\n✗ ERROR processing %s: %s\n", ct, conditionMessage(e)))
    processing_summary[[ct]] <- list(
      status = "FAILED",
      error = conditionMessage(e),
      time_seconds = elapsed_time
    )
  })
  
  cat("\n")
}

# ===== Compile Final Results =====
cat("\n===== Compiling Final Results =====\n")

# ===== Overall ROGUE Results =====
if (length(all_rogue_results) > 0) {
  
  # Combine all overall ROGUE results
  final_rogue_df <- bind_rows(all_rogue_results)
  
  # Add purity classification
  final_rogue_df <- final_rogue_df %>%
    mutate(
      purity_class = case_when(
        rogue_value >= 0.9 ~ "High (≥0.9)",
        rogue_value >= 0.7 ~ "Moderate (0.7-0.9)",
        TRUE ~ "Low (<0.7)"
      )
    ) %>%
    arrange(cell_type, desc(rogue_value))
  
  # Save combined overall results
  write.csv(
    final_rogue_df,
    file.path(OUTPUT_DIR, "all_celltype_cluster_rogue_overall.csv"),
    row.names = FALSE
  )
  
  cat(sprintf("Overall ROGUE values calculated: %d clusters\n", nrow(final_rogue_df)))
  cat(sprintf("Cell types with results: %d\n", length(unique(final_rogue_df$cell_type))))
  
  # ===== Summary Statistics (Overall) =====
  summary_stats <- final_rogue_df %>%
    group_by(cell_type) %>%
    summarise(
      n_clusters = n(),
      mean_rogue = mean(rogue_value),
      median_rogue = median(rogue_value),
      min_rogue = min(rogue_value),
      max_rogue = max(rogue_value),
      high_purity_pct = sum(rogue_value >= 0.9) / n() * 100,
      mean_filter_fraction = mean(filter_fraction, na.rm = TRUE),  # P1-2 addition
      low_retention_clusters = sum(filter_fraction < ROGUE_FILTER_WARNING_THRESHOLD, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_rogue))
  
  write.csv(
    summary_stats,
    file.path(OUTPUT_DIR, "rogue_summary_by_celltype_overall.csv"),
    row.names = FALSE
  )
  
  cat("\n=== Summary Statistics by Cell Type (Overall) ===\n")
  print(summary_stats, n = Inf)
  
  # Warn about clusters with low retention (P1-2 quality check)
  low_retention_total <- sum(final_rogue_df$filter_fraction < ROGUE_FILTER_WARNING_THRESHOLD, 
                             na.rm = TRUE)
  if (low_retention_total > 0) {
    cat(sprintf("\n⚠ WARNING: %d clusters have low cell retention (<%.0f%%) after ROGUE filtering\n",
                low_retention_total, ROGUE_FILTER_WARNING_THRESHOLD * 100))
    cat("  These ROGUE values may be inflated due to aggressive cell filtering.\n")
    cat("  Check 'filter_fraction' column in output files.\n")
  }
  
  # ===== Per-Sample ROGUE Results =====
  if (length(all_rogue_per_sample_results) > 0) {
    cat("\n===== Processing Per-Sample ROGUE Results =====\n")
    
    # Combine all per-sample results
    final_rogue_per_sample_df <- bind_rows(all_rogue_per_sample_results)
    
    # Add purity classification
    final_rogue_per_sample_df <- final_rogue_per_sample_df %>%
      mutate(
        purity_class = case_when(
          rogue_value >= 0.9 ~ "High (≥0.9)",
          rogue_value >= 0.7 ~ "Moderate (0.7-0.9)",
          TRUE ~ "Low (<0.7)"
        )
      ) %>%
      arrange(cell_type, cluster, sample)
    
    # Save combined per-sample results
    write.csv(
      final_rogue_per_sample_df,
      file.path(OUTPUT_DIR, "all_celltype_cluster_rogue_per_sample.csv"),
      row.names = FALSE
    )
    
    cat(sprintf("Per-sample ROGUE values: %d cluster×sample combinations\n", 
                nrow(final_rogue_per_sample_df)))
    
    # ===== Per-Sample Summary Statistics =====
    # By cell type and sample
    sample_summary <- final_rogue_per_sample_df %>%
      group_by(cell_type, sample) %>%
      summarise(
        n_clusters = n(),
        mean_rogue = mean(rogue_value),
        median_rogue = median(rogue_value),
        min_rogue = min(rogue_value),
        max_rogue = max(rogue_value),
        .groups = 'drop'
      ) %>%
      arrange(cell_type, desc(mean_rogue))
    
    write.csv(
      sample_summary,
      file.path(OUTPUT_DIR, "rogue_summary_by_celltype_and_sample.csv"),
      row.names = FALSE
    )
    
    # Identify problematic sample×cluster combinations (low ROGUE)
    problematic_combinations <- final_rogue_per_sample_df %>%
      filter(rogue_value < 0.7) %>%
      arrange(rogue_value) %>%
      select(cell_type, cluster, sample, n_cells, rogue_value)
    
    if (nrow(problematic_combinations) > 0) {
      write.csv(
        problematic_combinations,
        file.path(OUTPUT_DIR, "problematic_sample_cluster_combinations.csv"),
        row.names = FALSE
      )
      cat(sprintf("⚠ Identified %d problematic cluster×sample combinations (ROGUE < 0.7)\n", 
                  nrow(problematic_combinations)))
    }
    
    # ===== Compare Overall vs Per-Sample ROGUE =====
    # Merge overall and per-sample for comparison
    rogue_comparison <- final_rogue_df %>%
      select(cell_type, cluster, overall_rogue = rogue_value) %>%
      left_join(
        final_rogue_per_sample_df %>%
          group_by(cell_type, cluster) %>%
          summarise(
            mean_sample_rogue = mean(rogue_value),
            sd_sample_rogue = sd(rogue_value),
            min_sample_rogue = min(rogue_value),
            max_sample_rogue = max(rogue_value),
            n_samples = n(),
            .groups = 'drop'
          ),
        by = c("cell_type", "cluster")
      )
    
    write.csv(
      rogue_comparison,
      file.path(OUTPUT_DIR, "rogue_overall_vs_per_sample_comparison.csv"),
      row.names = FALSE
    )
    
    cat("\n=== Per-Sample Statistics Summary ===\n")
    cat(sprintf("Total samples analyzed: %d\n", 
                length(unique(final_rogue_per_sample_df$sample))))
    cat(sprintf("Average clusters per sample: %.1f\n", 
                nrow(final_rogue_per_sample_df) / length(unique(final_rogue_per_sample_df$sample))))
  } else {
    cat("\n⚠ No per-sample ROGUE results available\n")
  }
  
  # ===== Generate Combined Visualizations =====
  
  # 1. Box plot across cell types
  tryCatch({
    pdf(file.path(OUTPUT_DIR, "combined_rogue_boxplot.pdf"), width = 12, height = 6)
    p1 <- ggplot(final_rogue_df, aes(x = reorder(cell_type, rogue_value, FUN = median), 
                                     y = rogue_value, fill = cell_type)) +
      geom_boxplot(outlier.shape = 16, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
      geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red", alpha = 0.5) +
      coord_flip() +
      labs(title = "ROGUE Distribution Across Cell Types",
           subtitle = sprintf("Total: %d clusters from %d cell types", 
                             nrow(final_rogue_df), 
                             length(unique(final_rogue_df$cell_type))),
           x = "Cell Type",
           y = "ROGUE Value") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.y = element_text(size = 10),
            plot.title = element_text(face = "bold"))
    print(p1)
    dev.off()
    cat("✓ Generated combined_rogue_boxplot.pdf\n")
  }, error = function(e) {
    cat(sprintf("✗ Failed to generate boxplot: %s\n", conditionMessage(e)))
  })
  
  # 2. Heatmap: cell type x cluster rank (avoid cluster ID collision)
  tryCatch({
    pdf(file.path(OUTPUT_DIR, "combined_rogue_heatmap.pdf"), width = 14, height = 10)
    
    # Create cluster rank within each cell type (1 = highest ROGUE)
    rogue_rank_df <- final_rogue_df %>%
      group_by(cell_type) %>%
      arrange(desc(rogue_value), .by_group = TRUE) %>%
      mutate(cluster_rank = row_number()) %>%
      ungroup()
    
    # Prepare matrix with cluster_rank instead of cluster ID
    rogue_matrix <- rogue_rank_df %>%
      mutate(cluster_rank = as.character(cluster_rank)) %>%
      select(cell_type, cluster_rank, rogue_value) %>%
      pivot_wider(names_from = cluster_rank, values_from = rogue_value) %>%
      column_to_rownames("cell_type") %>%
      as.matrix()
    
    # Handle missing values
    rogue_matrix[is.na(rogue_matrix)] <- 0
    
    # Only create heatmap if matrix is not empty
    if (nrow(rogue_matrix) > 0 && ncol(rogue_matrix) > 0) {
      pheatmap::pheatmap(
        rogue_matrix,
        cluster_rows = TRUE,
        cluster_cols = FALSE,  # Keep cluster rank order
        color = colorRampPalette(c("firebrick", "white", "forestgreen"))(100),
        breaks = seq(0, 1, length.out = 101),
        border_color = "grey80",
        main = "ROGUE Values by Cell Type (Columns = Rank within Cell Type)",
        fontsize = 9,
        cellwidth = 15,
        cellheight = 15,
        display_numbers = TRUE,
        number_format = "%.2f",
        na_col = "grey90"
      )
    }
    dev.off()
    cat("✓ Generated combined_rogue_heatmap.pdf (cluster_rank based)\n")
  }, error = function(e) {
    cat(sprintf("✗ Failed to generate heatmap: %s\n", conditionMessage(e)))
  })
  
  # 3. Summary bar plot
  tryCatch({
    pdf(file.path(OUTPUT_DIR, "celltype_mean_rogue.pdf"), width = 10, height = 6)
    p3 <- ggplot(summary_stats, aes(x = reorder(cell_type, mean_rogue), y = mean_rogue)) +
      geom_col(aes(fill = mean_rogue), alpha = 0.8) +
      geom_errorbar(aes(ymin = min_rogue, ymax = max_rogue), width = 0.3, alpha = 0.5) +
      geom_text(aes(label = sprintf("%.3f", mean_rogue)), hjust = -0.1, size = 3.5) +
      geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red", alpha = 0.5) +
      coord_flip() +
      scale_fill_gradient2(low = "firebrick", mid = "orange", high = "forestgreen",
                           midpoint = 0.8) +
      labs(title = "Mean ROGUE by Cell Type",
           subtitle = "Error bars show min-max range across clusters",
           x = "Cell Type",
           y = "Mean ROGUE Value",
           fill = "Mean ROGUE") +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold")) +
      ylim(0, 1.05)
    print(p3)
    dev.off()
    cat("✓ Generated celltype_mean_rogue.pdf\n")
  }, error = function(e) {
    cat(sprintf("✗ Failed to generate mean ROGUE plot: %s\n", conditionMessage(e)))
  })
  
  # ===== Per-Sample ROGUE Visualizations =====
  if (length(all_rogue_per_sample_results) > 0 && exists("final_rogue_per_sample_df")) {
    cat("\n=== Generating Per-Sample Visualizations ===\n")
    
    # 4. Per-sample distribution by cell type (P0-2 fix: proper boxplot usage)
    tryCatch({
      pdf(file.path(OUTPUT_DIR, "per_sample_rogue_by_celltype.pdf"), width = 14, height = 8)
      
      # Each cell type gets a boxplot showing distribution across samples
      p4 <- ggplot(final_rogue_per_sample_df, 
                   aes(x = reorder(cell_type, rogue_value, FUN = median), 
                       y = rogue_value)) +
        geom_boxplot(aes(fill = cell_type), alpha = 0.7, outlier.shape = NA) +
        geom_point(aes(color = sample), alpha = 0.5, size = 2,
                  position = position_jitter(width = 0.2, seed = 42)) +
        geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", 
                  color = "red", alpha = 0.5) +
        coord_flip() +
        labs(title = "Per-Sample ROGUE Distribution by Cell Type",
             subtitle = "Boxplot = variation across all cluster×sample combinations",
             x = "Cell Type",
             y = "ROGUE Value",
             fill = "Cell Type",
             color = "Sample") +
        theme_bw() +
        theme(legend.position = "right",
              axis.text.y = element_text(size = 9),
              plot.title = element_text(face = "bold"))
      print(p4)
      dev.off()
      cat("✓ Generated per_sample_rogue_by_celltype.pdf\n")
    }, error = function(e) {
      cat(sprintf("✗ Failed to generate per-sample boxplot: %s\n", conditionMessage(e)))
    })
    
    # 5. Heatmap: Sample × Cell Type (better than the previous faceted plot)
    tryCatch({
      pdf(file.path(OUTPUT_DIR, "per_sample_rogue_heatmap.pdf"), width = 12, height = 10)
      
      # Aggregate per-sample data (mean ROGUE across clusters for each sample×celltype)
      sample_celltype_matrix <- final_rogue_per_sample_df %>%
        group_by(sample, cell_type) %>%
        summarise(mean_rogue = mean(rogue_value), .groups = 'drop') %>%
        pivot_wider(names_from = cell_type, values_from = mean_rogue) %>%
        column_to_rownames("sample") %>%
        as.matrix()
      
      sample_celltype_matrix[is.na(sample_celltype_matrix)] <- 0
      
      if (nrow(sample_celltype_matrix) > 0 && ncol(sample_celltype_matrix) > 0) {
        pheatmap::pheatmap(
          sample_celltype_matrix,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          color = colorRampPalette(c("firebrick", "white", "forestgreen"))(100),
          breaks = seq(0, 1, length.out = 101),
          border_color = "grey80",
          main = "Mean Per-Sample ROGUE: Sample × Cell Type",
          fontsize = 9,
          cellwidth = 20,
          cellheight = 15,
          display_numbers = TRUE,
          number_format = "%.2f",
          na_col = "grey90"
        )
      }
      dev.off()
      cat("✓ Generated per_sample_rogue_heatmap.pdf\n")
    }, error = function(e) {
      cat(sprintf("✗ Failed to generate sample heatmap: %s\n", conditionMessage(e)))
    })
    
    # 5. Sample variance plot (identify variable samples)
    tryCatch({
      sample_variance <- final_rogue_per_sample_df %>%
        group_by(sample) %>%
        summarise(
          n_clusters = n(),
          mean_rogue = mean(rogue_value),
          sd_rogue = sd(rogue_value),
          .groups = 'drop'
        ) %>%
        arrange(desc(sd_rogue))
      
      pdf(file.path(OUTPUT_DIR, "sample_rogue_variance.pdf"), width = 10, height = 6)
      p5 <- ggplot(sample_variance, 
                   aes(x = reorder(sample, sd_rogue), y = sd_rogue)) +
        geom_col(aes(fill = mean_rogue), alpha = 0.8) +
        geom_text(aes(label = sprintf("%.3f", sd_rogue)), hjust = -0.1, size = 3) +
        coord_flip() +
        scale_fill_gradient2(low = "firebrick", mid = "orange", high = "forestgreen",
                            midpoint = 0.8) +
        labs(title = "ROGUE Variance by Sample",
             subtitle = "High SD indicates inconsistent cluster quality within sample",
             x = "Sample",
             y = "Standard Deviation of ROGUE",
             fill = "Mean ROGUE") +
        theme_bw() +
        theme(plot.title = element_text(face = "bold"))
      print(p5)
      dev.off()
      cat("✓ Generated sample_rogue_variance.pdf\n")
    }, error = function(e) {
      cat(sprintf("✗ Failed to generate sample variance plot: %s\n", conditionMessage(e)))
    })
    
    # 6. Overall vs mean per-sample ROGUE comparison
    if (exists("rogue_comparison")) {
      tryCatch({
        pdf(file.path(OUTPUT_DIR, "overall_vs_sample_rogue_comparison.pdf"), 
            width = 10, height = 8)
        
        comparison_plot <- rogue_comparison %>%
          filter(!is.na(mean_sample_rogue))
        
        p6 <- ggplot(comparison_plot, 
                     aes(x = overall_rogue, y = mean_sample_rogue)) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
          geom_errorbar(aes(ymin = min_sample_rogue, ymax = max_sample_rogue), 
                       alpha = 0.3, width = 0.01) +
          geom_point(aes(color = cell_type, size = n_samples), alpha = 0.7) +
          geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
          labs(title = "Overall vs Per-Sample ROGUE Comparison",
               subtitle = "Points on diagonal = consistent across samples; error bars = min-max range",
               x = "Overall Cluster ROGUE",
               y = "Mean Per-Sample ROGUE",
               color = "Cell Type",
               size = "# Samples") +
          theme_bw() +
          theme(plot.title = element_text(face = "bold"),
                legend.position = "right") +
          coord_fixed()
        print(p6)
        dev.off()
        cat("✓ Generated overall_vs_sample_rogue_comparison.pdf\n")
      }, error = function(e) {
        cat(sprintf("✗ Failed to generate comparison plot: %s\n", conditionMessage(e)))
      })
    }
  }
  
} else {
  cat("WARNING: No ROGUE values were successfully calculated\n")
}

# ===== Processing Summary =====
cat("\n=== Processing Summary ===\n")
summary_df <- bind_rows(lapply(names(processing_summary), function(ct) {
  s <- processing_summary[[ct]]
  data.frame(
    cell_type = ct,
    status = s$status,
    n_cells = ifelse(is.null(s$n_cells), NA, s$n_cells),
    n_batches = ifelse(is.null(s$n_batches), NA, s$n_batches),
    n_clusters = ifelse(is.null(s$n_clusters), NA, s$n_clusters),
    n_rogue = ifelse(is.null(s$n_rogue_calculated), 0, s$n_rogue_calculated),
    reduction = ifelse(is.null(s$reduction_used), NA, s$reduction_used),
    time_sec = ifelse(is.null(s$time_seconds), NA, round(s$time_seconds, 1)),
    error_msg = ifelse(is.null(s$error), "", as.character(s$error)),
    stringsAsFactors = FALSE
  )
}))

# Sort by status and time
summary_df <- summary_df %>%
  arrange(status, desc(time_sec))

write.csv(summary_df, file.path(OUTPUT_DIR, "processing_summary.csv"), row.names = FALSE)
print(summary_df, row.names = FALSE)

# Calculate timing statistics
if (nrow(summary_df) > 0 && any(!is.na(summary_df$time_sec))) {
  total_time <- sum(summary_df$time_sec, na.rm = TRUE)
  mean_time <- mean(summary_df$time_sec, na.rm = TRUE)
  cat(sprintf("\nTiming Statistics:\n"))
  cat(sprintf("  Total processing time: %.1f seconds (%.1f minutes)\n", total_time, total_time/60))
  cat(sprintf("  Average time per cell type: %.1f seconds\n", mean_time))
}

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
cat(sprintf("Successfully processed: %d/%d cell types\n", 
            sum(summary_df$status == "SUCCESS"), nrow(summary_df)))

if (any(summary_df$status == "FAILED")) {
  failed_cts <- summary_df %>% filter(status == "FAILED") %>% pull(cell_type)
  cat(sprintf("\n⚠ Failed cell types (%d): %s\n", 
              length(failed_cts), paste(failed_cts, collapse = ", ")))
  cat("Check error messages in processing_summary.csv\n")
}