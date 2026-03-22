# Install required packages if not already installed
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  install.packages("hdf5r")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}

library(hdf5r)
library(Matrix)

# Function to convert a single h5 file
h5_to_10x <- function(h5_path, output_dir) {
  tryCatch({
    # Convert file paths to use forward slashes
    h5_path <- gsub("\\\\", "/", h5_path)
    output_dir <- gsub("\\\\", "/", output_dir)
    
    # Create output directory if it doesn't exist
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Open the H5 file
    h5_data <- H5File$new(h5_path, mode = "r")
    
    # Read the matrix data
    counts <- h5_data[["matrix/data"]]$read()
    indices <- h5_data[["matrix/indices"]]$read()
    indptr <- h5_data[["matrix/indptr"]]$read()
    
    # Read gene names and barcodes
    genes <- h5_data[["matrix/features/name"]]$read()
    gene_ids <- h5_data[["matrix/features/id"]]$read()
    barcodes <- h5_data[["matrix/barcodes"]]$read()
    
    # Create sparse matrix
    sparse_matrix <- sparseMatrix(
      i = indices + 1,
      p = indptr,
      x = counts,
      dims = c(length(genes), length(barcodes))
    )
    
    # Write files in 10x format
    # 1. features.tsv.gz (genes)
    genes_df <- data.frame(
      gene_ids = gene_ids,
      gene_names = genes,
      feature_types = rep("Gene Expression", length(genes))
    )
    
    features_path <- file.path(output_dir, "features.tsv.gz")
    features_path <- gsub("\\\\", "/", features_path)
    write.table(
      genes_df,
      file = gzfile(features_path),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    
    # 2. barcodes.tsv.gz
    barcodes_path <- file.path(output_dir, "barcodes.tsv.gz")
    barcodes_path <- gsub("\\\\", "/", barcodes_path)
    write.table(
      barcodes,
      file = gzfile(barcodes_path),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    
    # 3. matrix.mtx.gz
    matrix_path <- file.path(output_dir, "matrix.mtx")
    matrix_path <- gsub("\\\\", "/", matrix_path)
    Matrix::writeMM(sparse_matrix, file = matrix_path)
    
    # Compress matrix file
    R.utils::gzip(matrix_path, overwrite = TRUE)
    
    # Close the H5 file
    h5_data$close()
    
    message("Successfully converted: ", basename(h5_path))
    return(TRUE)
  }, error = function(e) {
    message("Error processing ", basename(h5_path), ": ", e$message)
    return(FALSE)
  })
}

# Function to process all h5 files in a directory
process_h5_directory <- function(input_dir, output_base_dir) {
  # Normalize paths
  input_dir <- normalizePath(input_dir, winslash = "/")
  output_base_dir <- normalizePath(output_base_dir, winslash = "/", mustWork = FALSE)
  
  # Get list of all h5 files
  h5_files <- list.files(
    path = input_dir,
    pattern = "\\.h5$",
    full.names = TRUE
  )
  
  if (length(h5_files) == 0) {
    stop("No .h5 files found in the input directory")
  }
  
  # Process each file
  results <- data.frame(
    file = basename(h5_files),
    status = character(length(h5_files)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(h5_files)) {
    h5_file <- h5_files[i]
    # Create output directory with same name as h5 file (without extension)
    file_name <- tools::file_path_sans_ext(basename(h5_file))
    output_dir <- file.path(output_base_dir, file_name)
    
    message("\nProcessing file ", i, " of ", length(h5_files), ": ", basename(h5_file))
    success <- h5_to_10x(h5_file, output_dir)
    
    results$status[i] <- if(success) "Success" else "Failed"
  }
  
  # Write processing report
  report_path <- file.path(output_base_dir, "conversion_report.csv")
  write.csv(results, report_path, row.names = FALSE)
  
  message("\nProcessing complete!")
  message("Total files processed: ", nrow(results))
  message("Successful conversions: ", sum(results$status == "Success"))
  message("Failed conversions: ", sum(results$status == "Failed"))
  message("Report saved to: ", report_path)
}

# Example usage:
process_h5_directory("E:/R/Source/20250116", "E:/R/Output/20250116")