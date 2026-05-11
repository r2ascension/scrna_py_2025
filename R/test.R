suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})


# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
# output_dir <- "/home/h2048/data/R/0103/ciliated"
# dir.create(output_dir,recursive=TRUE)
# setwd()
library(reticulate)
library(data.table)
library(harmony)
library(ggplot2)
library(patchwork)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()
Markers <- c(
  'FXYD3',
  'EPCAM',
  'ELF3',
  'SERPINF1',
  'TSPAN1',
  'SCGB1A1',
  'AGER',
  'SFTPC',
  'FOXJ1',
  'KRT5',
  'MUC5B',
  'KRT8',
  'CD53',
  'PTPRC',
  'CORO1A',
  'CCL5',
  'MS4A1',
  'TNFRSF17',
  'CD19',
  'CD79A',
  'CD40LG',
  'TNFRSF25',
  'CD28',
  'CD4',
  'CD3E',
  'CD8A',
  'CD8B',
  'TRGC2',
  'CD2',
  'TRBC2',
  'FCER1G',
  'C1orf162',
  'CLEC7A',
  'CD1C',
  'CD86',
  'CD14',
  'XCR1',
  'HLA-DRA',
  'COL1A2',
  'DCN',
  'MFAP4',
  'LUM',
  'COL6A3',
  'CFD',
  'COL1A1',
  'PDGFRA',
  'MXRA8',
  'LEPR',
  'MYH11',
  'TINAGL1',
  'PLN',
  'DES',
  'ACTA2',
  'CNN1',
  'TAGLN',
  'CLDN5',
  'ECSCR',
  'CLEC14A',
  'VWF',
  'PECAM1',
  'DARC',
  'PTPRB',
  'PDE2A',
  'PLAT',
  'GJA5',
  'SPARCL1',
  'AQP1',
  'MMRN1',
  'CCL21',
  'MKI67',
  'TOP2A',
  'TK1',
  'CENPW'
) #Proliferation
markers <- Markers # 你上面那串向量

output_dir <- '/home/h2048/data/R/0118/nose_n_airway/'
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

seurat_obj <- readRDS(
  '/home/h2048/data/source/1124/processed/seurat_without_partial_dual_contamination_doublets_2_1028.rds'
)
table(seurat_obj$tissue, seurat_obj$tissue_sampling_method)

table(
  seurat_obj@meta.data$tissue_sampling_method[
    seurat_obj@meta.data$tissue_level_2 == "inferior turbinate"
  ],
  useNA = "ifany"
)

seurat_obj <- subset(
  seurat_obj,
  subset = tissue %in% c('nose', 'respiratory airway')
)
saveRDS(seurat_obj, 'nose_n_airway.rds')


table(seurat_obj$dataset, seurat_obj$batch)
c(unique(seurat_obj$sample))
# 1) 基因是否存在（大小写不敏感匹配到对象真实gene名）
assay_use <- DefaultAssay(seurat_obj)
feat <- rownames(seurat_obj[[assay_use]])
idx <- match(toupper(markers), toupper(feat))

markers_present <- feat[idx[!is.na(idx)]]
markers_missing <- markers[is.na(idx)]

cat(sprintf(
  "Assay: %s\nMarkers in object: %d/%d\n",
  assay_use,
  length(markers_present),
  length(markers)
))
if (length(markers_missing)) {
  cat("Missing markers:\n")
  print(markers_missing)
}

# 2) 快速看表达（DotPlot：默认按当前Idents分组）
#    若还没NormalizeData，DotPlot会很不稳；这里做个轻量兜底
if (!"data" %in% SeuratObject::Layers(seurat_obj[[assay_use]])) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
}

p_dot <- DotPlot(seurat_obj, features = markers_present) +
  RotatedAxis()

print(p_dot)

# 3) 可选：挑前12个做FeaturePlot快速扫一眼
for (g in head(markers_present, 12)) {
  print(FeaturePlot(seurat_obj, features = g, raster = TRUE) + ggtitle(g))
}


METADATA_CSV <- "/home/h2048/data/source/reference/metadata_full_atlas_20251210.csv"
SAMPLE_COLUMN <- "sample" # Column name in both Seurat object and CSV

# Step 1: Read reference metadata CSV
cat("Reading reference metadata from CSV...\n")
ref_metadata <- read.csv(METADATA_CSV, stringsAsFactors = FALSE)
cat(sprintf(
  "Reference metadata dimensions: %d rows, %d columns\n",
  nrow(ref_metadata),
  ncol(ref_metadata)
))
cat(
  "Reference columns:",
  paste(colnames(ref_metadata), collapse = ", "),
  "\n\n"
)

# Step 2: Get current sample names from Seurat object
cat("Extracting sample information from Seurat object...\n")
current_samples <- unique(seurat_obj@meta.data[[SAMPLE_COLUMN]])
cat(sprintf(
  "Found %d unique samples in Seurat object\n",
  length(current_samples)
))
cat(
  "Sample preview:",
  paste(head(current_samples, 3), collapse = ", "),
  "...\n\n"
)

# Step 3: Check if samples need _seurat suffix matching
cat("Checking sample name matching...\n")
ref_samples <- unique(ref_metadata[[SAMPLE_COLUMN]])

# Try direct matching first
direct_match <- sum(current_samples %in% ref_samples)
cat(sprintf("Direct matches: %d/%d\n", direct_match, length(current_samples)))

# Try with _seurat suffix
suffix_match <- sum(paste0(current_samples, "_seurat") %in% ref_samples)
cat(sprintf(
  "Matches with '_seurat' suffix: %d/%d\n",
  suffix_match,
  length(current_samples)
))

# Determine matching strategy
if (direct_match >= suffix_match) {
  cat("\nUsing direct matching strategy\n")
  use_suffix <- FALSE
  seurat_obj@meta.data$sample_matched <- seurat_obj@meta.data[[SAMPLE_COLUMN]]
} else {
  cat("\nUsing '_seurat' suffix matching strategy\n")
  use_suffix <- TRUE
  seurat_obj@meta.data$sample_matched <- paste0(
    seurat_obj@meta.data[[SAMPLE_COLUMN]],
    "_seurat"
  )
}

# Step 4: Merge metadata
cat("\nMerging metadata...\n")

# Get columns to update (exclude the sample column itself)
update_columns <- setdiff(colnames(ref_metadata), SAMPLE_COLUMN)
cat("Columns to update:", paste(update_columns, collapse = ", "), "\n")

# Create a temporary data frame with cell barcodes and matched sample names
temp_df <- data.frame(
  cell_barcode = rownames(seurat_obj@meta.data),
  sample_matched = seurat_obj@meta.data$sample_matched,
  stringsAsFactors = FALSE
)

# Merge with reference metadata using base R merge
merged_metadata <- merge(
  temp_df,
  ref_metadata,
  by.x = "sample_matched",
  by.y = SAMPLE_COLUMN,
  all.x = TRUE,
  sort = FALSE
)

# Restore original row order
merged_metadata <- merged_metadata[
  match(temp_df$cell_barcode, merged_metadata$cell_barcode),
]

# Step 5: Update Seurat object metadata
cat("\nUpdating Seurat object metadata...\n")
for (col in update_columns) {
  if (col %in% colnames(merged_metadata)) {
    seurat_obj@meta.data[[col]] <- merged_metadata[[col]]

    # Count non-NA values
    non_na_count <- sum(!is.na(merged_metadata[[col]]))
    cat(sprintf(
      "  - %s: %d/%d cells updated\n",
      col,
      non_na_count,
      nrow(seurat_obj@meta.data)
    ))
  }
}

# Step 6: Clean up temporary column
seurat_obj@meta.data$sample_matched <- NULL

# Step 7: Summary
cat("\n===== Update Summary =====\n")
cat(sprintf("Total cells: %d\n", nrow(seurat_obj@meta.data)))
cat(sprintf("Metadata columns after update: %d\n", ncol(seurat_obj@meta.data)))
cat("\nNew/updated columns from CSV:\n")
for (col in update_columns) {
  if (col %in% colnames(seurat_obj@meta.data)) {
    na_count <- sum(is.na(seurat_obj@meta.data[[col]]))
    cat(sprintf(
      "  - %s: %d cells with data, %d NA\n",
      col,
      nrow(seurat_obj@meta.data) - na_count,
      na_count
    ))
  }
}

# Step 8: Check for unmatched samples
unmatched_samples <- setdiff(
  if (use_suffix) paste0(current_samples, "_seurat") else current_samples,
  ref_samples
)
if (length(unmatched_samples) > 0) {
  cat("\n⚠️ Warning: Following samples in Seurat object not found in CSV:\n")
  cat(paste(unmatched_samples, collapse = "\n"), "\n")
}

cat("\n✓ Metadata update completed!\n")

# Optional: View updated metadata
cat("\nPreview of updated metadata:\n")
print(head(seurat_obj@meta.data[, c(SAMPLE_COLUMN, update_columns)], 10))


# ===== Diagnostic: Check CSV Sample Names =====

# Read CSV and check sample format
ref_metadata <- read.csv(
  "/home/h2048/data/source/reference/metadata_full_atlas_20251210.csv",
  stringsAsFactors = FALSE
)

cat("=== CSV Sample Format Diagnosis ===\n")
cat("First 10 samples in CSV:\n")
print(head(ref_metadata$sample, 10))

cat(
  "\nSamples with '_seurat' suffix in CSV:",
  sum(grepl("_seurat$", ref_metadata$sample)),
  "/",
  nrow(ref_metadata),
  "\n"
)

cat(
  "\nSamples without '_seurat' suffix in CSV:",
  sum(!grepl("_seurat$", ref_metadata$sample)),
  "/",
  nrow(ref_metadata),
  "\n"
)

# Check Seurat object samples
seurat_samples <- unique(seurat_obj@meta.data$sample)
cat("\n=== Seurat Object Sample Format ===\n")
cat("First 10 samples in Seurat object:\n")
print(head(seurat_samples, 10))

cat(
  "\nSamples with '_seurat' suffix in Seurat:",
  sum(grepl("_seurat$", seurat_samples)),
  "/",
  length(seurat_samples),
  "\n"
)

# Try matching strategies
cat("\n=== Matching Test ===\n")

# Strategy 1: Direct match
direct_matches <- sum(seurat_samples %in% ref_metadata$sample)
cat(
  "Strategy 1 - Direct match: ",
  direct_matches,
  "/",
  length(seurat_samples),
  "\n"
)

# Strategy 2: Remove _seurat from Seurat samples
seurat_samples_clean <- gsub("_seurat$", "", seurat_samples)
match_remove_suffix <- sum(seurat_samples_clean %in% ref_metadata$sample)
cat(
  "Strategy 2 - Remove '_seurat' from Seurat samples: ",
  match_remove_suffix,
  "/",
  length(seurat_samples),
  "\n"
)

# Strategy 3: Add _seurat to CSV samples
csv_samples_with_suffix <- paste0(ref_metadata$sample, "_seurat")
match_add_suffix <- sum(seurat_samples %in% csv_samples_with_suffix)
cat(
  "Strategy 3 - Add '_seurat' to CSV samples: ",
  match_add_suffix,
  "/",
  length(seurat_samples),
  "\n"
)

cat("\n✓ Diagnosis complete. Choose the best matching strategy.\n")


setwd('/home/h2048/data/source/polyp')
METADATA_CSV <- "/home/h2048/data/source/reference/metadata_full_atlas_20251210.csv"
SAMPLE_COLUMN <- "sample"

# Step 1: Read reference metadata CSV
cat("Reading reference metadata from CSV...\n")
ref_metadata <- read.csv(METADATA_CSV, stringsAsFactors = FALSE)
cat(sprintf(
  "Reference metadata: %d rows, %d columns\n",
  nrow(ref_metadata),
  ncol(ref_metadata)
))

# Step 2: Create matching key by removing _seurat suffix
cat("\nPreparing sample matching...\n")
seurat_obj@meta.data$sample_for_matching <- gsub(
  "_seurat$",
  "",
  seurat_obj@meta.data[[SAMPLE_COLUMN]]
)

# Check matching success
unique_seurat_samples <- unique(seurat_obj@meta.data$sample_for_matching)
unique_csv_samples <- unique(ref_metadata[[SAMPLE_COLUMN]])
matched_count <- sum(unique_seurat_samples %in% unique_csv_samples)
cat(sprintf(
  "Matching result: %d/%d samples matched\n",
  matched_count,
  length(unique_seurat_samples)
))

# Step 3: Merge metadata
cat("\nMerging metadata...\n")
update_columns <- setdiff(colnames(ref_metadata), SAMPLE_COLUMN)
cat("Columns to update:", paste(update_columns, collapse = ", "), "\n\n")

# Create temporary dataframe
temp_df <- data.frame(
  cell_barcode = rownames(seurat_obj@meta.data),
  sample_for_matching = seurat_obj@meta.data$sample_for_matching,
  stringsAsFactors = FALSE
)

# Merge using base R
merged_metadata <- merge(
  temp_df,
  ref_metadata,
  by.x = "sample_for_matching",
  by.y = SAMPLE_COLUMN,
  all.x = TRUE,
  sort = FALSE
)

# Restore original row order
merged_metadata <- merged_metadata[
  match(temp_df$cell_barcode, merged_metadata$cell_barcode),
]

# Step 4: Update Seurat object metadata
cat("Updating Seurat object metadata...\n")
for (col in update_columns) {
  if (col %in% colnames(merged_metadata)) {
    seurat_obj@meta.data[[col]] <- merged_metadata[[col]]

    non_na_count <- sum(!is.na(merged_metadata[[col]]))
    cat(sprintf(
      "  - %-30s: %d/%d cells updated\n",
      col,
      non_na_count,
      nrow(seurat_obj@meta.data)
    ))
  }
}

# Step 5: Clean up temporary column
seurat_obj@meta.data$sample_for_matching <- NULL

# Step 6: Summary
cat("\n===== Update Summary =====\n")
cat(sprintf("Total cells: %d\n", nrow(seurat_obj@meta.data)))
cat(sprintf("Total metadata columns: %d\n", ncol(seurat_obj@meta.data)))

cat("\nUpdate statistics by column:\n")
for (col in update_columns) {
  if (col %in% colnames(seurat_obj@meta.data)) {
    na_count <- sum(is.na(seurat_obj@meta.data[[col]]))
    non_na_count <- nrow(seurat_obj@meta.data) - na_count
    percentage <- round(100 * non_na_count / nrow(seurat_obj@meta.data), 1)
    cat(sprintf(
      "  %-30s: %d cells (%s%%) with data\n",
      col,
      non_na_count,
      percentage
    ))
  }
}

# Step 7: Check for unmatched samples
unmatched_seurat <- setdiff(unique_seurat_samples, unique_csv_samples)
if (length(unmatched_seurat) > 0) {
  cat(sprintf(
    "\n⚠️  Warning: %d samples in Seurat not found in CSV:\n",
    length(unmatched_seurat)
  ))
  cat(paste(head(unmatched_seurat, 10), collapse = ", "), "\n")
  if (length(unmatched_seurat) > 10) {
    cat(sprintf("... and %d more\n", length(unmatched_seurat) - 10))
  }
}

unmatched_csv <- setdiff(unique_csv_samples, unique_seurat_samples)
if (length(unmatched_csv) > 0) {
  cat(sprintf(
    "\nℹ️  Info: %d samples in CSV not found in Seurat:\n",
    length(unmatched_csv)
  ))
  cat(paste(head(unmatched_csv, 10), collapse = ", "), "\n")
  if (length(unmatched_csv) > 10) {
    cat(sprintf("... and %d more\n", length(unmatched_csv) - 10))
  }
}

cat("\n✓ Metadata update completed!\n")

# Step 8: Preview updated metadata
cat("\nPreview of updated metadata (first 5 cells):\n")
preview_cols <- c(SAMPLE_COLUMN, head(update_columns, 5))
print(head(seurat_obj@meta.data[, preview_cols], 5))

# Optional: Save updated object
saveRDS(seurat_obj, "polyp_obj_updated_20260113.rds")
GetH5ad(seurat_obj, 'polyp_obj_updated_20260113.h5ad')
seurat_obj_polyp <- subset(seurat_obj, subset = condition %in% c('nasal polyp'))
saveRDS(seurat_obj_polyp, "polyp_obj_updated_purepolyp_20260113.rds")
table(seurat_obj$dataset, seurat_obj$condition)
seurat_obj$condition[
  seurat_obj$condition %in% c('nasal polyps')
] <- 'nasal polyp'
seurat_obj$tissue_level_2[
  seurat_obj$tissue_level_2 %in% c('nasal polyps')
] <- 'nasal polyp'
# Load required libraries
# library(CHOIR)
library(Seurat)
library(reticulate)
library(dplyr)
library(SCNT)
library(data.table)
output_dir <- '/home/h2048/data/R/0113'
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
library(reticulate)
library(harmony)
library(ggplot2)
# Specify conda environment by name
use_condaenv("bbknn_env", required = TRUE)
# Verify the environment
py_config()


# ===== Fix Metadata Types Before Writing H5AD =====

cat("Checking and fixing metadata types...\n\n")

# Step 1: Identify problematic columns
cat("=== Metadata Type Diagnosis ===\n")
meta_types <- sapply(seurat_obj@meta.data, class)
factor_cols <- names(meta_types[meta_types == "factor"])
numeric_cols <- names(meta_types[meta_types %in% c("numeric", "integer")])
character_cols <- names(meta_types[meta_types == "character"])

cat("Factor columns (", length(factor_cols), "):\n")
if (length(factor_cols) > 0) {
  cat(paste(head(factor_cols, 10), collapse = ", "), "\n")
  if (length(factor_cols) > 10) {
    cat("... and", length(factor_cols) - 10, "more\n")
  }
}

cat("\nNumeric columns (", length(numeric_cols), "):\n")
if (length(numeric_cols) > 0) {
  cat(paste(head(numeric_cols, 10), collapse = ", "), "\n")
}

cat("\nCharacter columns (", length(character_cols), "):\n")
if (length(character_cols) > 0) {
  cat(paste(head(character_cols, 10), collapse = ", "), "\n")
}

# Step 2: Check cellLabel specifically
if ("cellLabel" %in% colnames(seurat_obj@meta.data)) {
  cat("\n=== cellLabel Column Details ===\n")
  cat("Type:", class(seurat_obj@meta.data$cellLabel), "\n")
  cat("Sample values:\n")
  print(head(table(seurat_obj@meta.data$cellLabel), 10))
}

# Step 3: Fix all factor columns to character
cat("\n=== Converting Factor Columns to Character ===\n")
if (length(factor_cols) > 0) {
  for (col in factor_cols) {
    seurat_obj@meta.data[[col]] <- as.character(seurat_obj@meta.data[[col]])
    cat(sprintf("  ✓ Converted %s to character\n", col))
  }
} else {
  cat("No factor columns found.\n")
}

# Step 4: Handle NA values in character columns
cat("\n=== Handling NA Values ===\n")
for (col in colnames(seurat_obj@meta.data)) {
  if (is.character(seurat_obj@meta.data[[col]])) {
    na_count <- sum(is.na(seurat_obj@meta.data[[col]]))
    if (na_count > 0) {
      # Convert NA to "Unknown" for h5ad compatibility
      seurat_obj@meta.data[[col]][is.na(seurat_obj@meta.data[[
        col
      ]])] <- "Unknown"
      cat(sprintf(
        "  ✓ Replaced %d NA values in %s with 'Unknown'\n",
        na_count,
        col
      ))
    }
  }
}

# Step 5: Verify cellLabel is now character
if ("cellLabel" %in% colnames(seurat_obj@meta.data)) {
  cat("\n=== Final cellLabel Verification ===\n")
  cat("Type:", class(seurat_obj@meta.data$cellLabel), "\n")
  cat("Contains NA:", any(is.na(seurat_obj@meta.data$cellLabel)), "\n")
  cat("Unique values:", length(unique(seurat_obj@meta.data$cellLabel)), "\n")
}

cat("\n✓ Metadata type conversion completed!\n")

# Step 6: Try writing h5ad again
cat("\n=== Attempting to Write H5AD ===\n")
tryCatch(
  {
    GetH5ad(seurat_obj, 'polypobj_updated.h5ad')
    cat("✓ Successfully wrote polypobj_updated.h5ad\n")
  },
  error = function(e) {
    cat("❌ Error writing h5ad:\n")
    cat(conditionMessage(e), "\n")
    cat("\nTrying alternative approach...\n")

    # Alternative: Remove problematic columns temporarily
    cat("\nIdentifying problematic columns...\n")

    # Check for any remaining non-standard types
    problematic_cols <- c()
    for (col in colnames(seurat_obj@meta.data)) {
      col_class <- class(seurat_obj@meta.data[[col]])
      if (!col_class %in% c("character", "numeric", "integer", "logical")) {
        problematic_cols <- c(problematic_cols, col)
        cat(sprintf("  - %s: %s\n", col, col_class))
      }
    }

    if (length(problematic_cols) > 0) {
      cat("\nRemoving problematic columns and retrying...\n")
      seurat_obj_clean <- seurat_obj
      for (col in problematic_cols) {
        seurat_obj_clean@meta.data[[col]] <- NULL
        cat(sprintf("  ✓ Removed %s\n", col))
      }

      GetH5ad(seurat_obj_clean, 'polypobj_updated_clean.h5ad')
      cat(
        "✓ Successfully wrote polypobj_updated_clean.h5ad (with some columns removed)\n"
      )
      cat("Removed columns:", paste(problematic_cols, collapse = ", "), "\n")
    }
  }
)

# ===== Check Non-Unknown Metadata Values =====

# List of newly added metadata columns
new_metadata_cols <- c(
  "GEO",
  "dataset",
  "tissue",
  "tissue_level_2",
  "disease_level_1",
  "disease_level_2",
  "condition",
  "tissue_sampling_method",
  "filtered_or_raw",
  "target_cell",
  "frozen_or_fresh",
  "sex",
  "age",
  "platform",
  "reference_genome",
  "assay",
  "tissue_dissociation_protocol"
)

cat("===== Non-Unknown Metadata Analysis =====\n\n")

for (col in new_metadata_cols) {
  if (col %in% colnames(seurat_obj@meta.data)) {
    cat("========================================\n")
    cat("Column:", col, "\n")
    cat("----------------------------------------\n")

    # Get non-Unknown values
    values <- seurat_obj@meta.data[[col]]
    non_unknown <- values[values != "Unknown"]

    # Statistics
    total_cells <- length(values)
    unknown_cells <- sum(values == "Unknown")
    non_unknown_cells <- total_cells - unknown_cells
    percentage <- round(100 * non_unknown_cells / total_cells, 2)

    cat(sprintf("Total cells: %d\n", total_cells))
    cat(sprintf("Unknown: %d (%.2f%%)\n", unknown_cells, 100 - percentage))
    cat(sprintf("Non-Unknown: %d (%.2f%%)\n", non_unknown_cells, percentage))

    if (non_unknown_cells > 0) {
      cat("\nNon-Unknown value distribution:\n")
      value_table <- sort(table(non_unknown), decreasing = TRUE)

      # Show top 20 or all if fewer
      n_show <- min(20, length(value_table))
      for (i in 1:n_show) {
        val_name <- names(value_table)[i]
        val_count <- value_table[i]
        val_pct <- round(100 * val_count / total_cells, 2)
        cat(sprintf(
          "  - %-40s: %7d cells (%.2f%%)\n",
          val_name,
          val_count,
          val_pct
        ))
      }

      if (length(value_table) > 20) {
        cat(sprintf(
          "  ... and %d more unique values\n",
          length(value_table) - 20
        ))
      }

      # Show which samples have this data
      cat("\nSamples with non-Unknown values:\n")
      samples_with_data <- unique(seurat_obj@meta.data$sample[
        values != "Unknown"
      ])
      cat(sprintf("  Total: %d samples\n", length(samples_with_data)))
      cat("  Sample list:\n")
      cat(paste("   ", head(samples_with_data, 10), collapse = "\n"), "\n")
      if (length(samples_with_data) > 10) {
        cat(sprintf(
          "  ... and %d more samples\n",
          length(samples_with_data) - 10
        ))
      }
    } else {
      cat("\n⚠️  All values are 'Unknown'\n")
    }

    cat("\n")
  } else {
    cat("Column", col, "not found in metadata\n\n")
  }
}

# Summary: Show samples without metadata
cat("========================================\n")
cat("SUMMARY: Samples Missing Metadata\n")
cat("========================================\n\n")

# Get samples with all Unknown values
all_samples <- unique(seurat_obj@meta.data$sample)
samples_with_metadata <- unique(seurat_obj@meta.data$sample[
  seurat_obj@meta.data$GEO != "Unknown"
])
samples_without_metadata <- setdiff(all_samples, samples_with_metadata)

cat(sprintf("Total samples: %d\n", length(all_samples)))
cat(sprintf("Samples with metadata: %d\n", length(samples_with_metadata)))
cat(sprintf("Samples without metadata: %d\n", length(samples_without_metadata)))

if (length(samples_without_metadata) > 0) {
  cat("\nSamples without metadata (need to add to CSV):\n")
  cat(paste(samples_without_metadata, collapse = "\n"), "\n")
}

cat("\n✓ Analysis completed!\n")
table(seurat_obj$orig.ident[
  seurat_obj$sample %in% c('Jose_Ordovas_Montanes_2018')
])
seurat_obj$sample[
  seurat_obj$sample %in%
    c('Jose_Ordovas_Montanes_2018', 'Jose_Ordovas_Montanes_2018_1')
] <- seurat_obj$orig.ident[
  seurat_obj$sample %in%
    c('Jose_Ordovas_Montanes_2018', 'Jose_Ordovas_Montanes_2018_1')
]
table(seurat_obj$sample)


# ===== Update Sample Names from sample_id =====

# Configuration
OLD_SAMPLE <- 'Kerstin_B_Meyer_2021_covid'

# Step 1: Check if sample_id column exists
if (!"sample_id" %in% colnames(seurat_obj@meta.data)) {
  stop("❌ Error: 'sample_id' column not found in metadata")
}

# Step 2: Check current status
cat("=== Before Update ===\n")
mask <- seurat_obj$sample == OLD_SAMPLE
n_cells <- sum(mask)
cat(sprintf("Cells with sample='%s': %d\n", OLD_SAMPLE, n_cells))

if (n_cells > 0) {
  # Show unique sample_id values for these cells
  unique_sample_ids <- unique(seurat_obj$sample_id[mask])
  cat(sprintf("Unique sample_id values: %d\n", length(unique_sample_ids)))
  cat("Sample_id values:\n")
  print(table(seurat_obj$sample_id[mask]))

  # Step 3: Perform update
  cat("\n=== Performing Update ===\n")
  seurat_obj$sample[mask] <- seurat_obj$sample_id[mask]

  # Step 4: Verify update
  cat("\n=== After Update ===\n")
  cat(sprintf(
    "Cells with sample='%s': %d (should be 0)\n",
    OLD_SAMPLE,
    sum(seurat_obj$sample == OLD_SAMPLE)
  ))

  # Show new sample distribution for updated cells
  cat("\nNew sample values for updated cells:\n")
  print(table(seurat_obj$sample[mask]))

  cat("\n✓ Update completed successfully!\n")
} else {
  cat(sprintf("\n⚠️  Warning: No cells found with sample='%s'\n", OLD_SAMPLE))
}

# Optional: Check if there are other samples that might need similar updates
cat("\n=== Other Samples That Might Need Updates ===\n")
all_samples <- unique(seurat_obj$sample)
suspicious_samples <- all_samples[!grepl("^GSM|^HRR", all_samples)]
if (length(suspicious_samples) > 0) {
  cat("Samples without GSM/HRR prefix:\n")
  for (s in suspicious_samples) {
    n <- sum(seurat_obj$sample == s)
    cat(sprintf("  - %-50s: %d cells\n", s, n))
  }
} else {
  cat("All samples have standard GSM/HRR prefixes.\n")
}

table(seurat_obj$dataset[
  seurat_obj$sample %in%
    c(
      'Polyp15AEPI',
      'Polyp16BEPI',
      'Polyp17AEPI',
      'Polyp17BEPI',
      'Polyp19AEPI',
      'Polyp19BEPI',
      'Polyp11TOT',
      'Polyp12TOT',
      'Polyp1TOT',
      'Polyp2TOT',
      'Polyp3TOT',
      'Polyp4TOT',
      'Polyp5TOT',
      'Polyp6ATOT',
      'Polyp6BTOT',
      'Polyp7TOT',
      'Polyp8TOT',
      'Polyp9TOT'
    )
])
seurat_obj$dataset[
  seurat_obj$sample %in%
    c(
      'Polyp15AEPI',
      'Polyp16BEPI',
      'Polyp17AEPI',
      'Polyp17BEPI',
      'Polyp19AEPI',
      'Polyp19BEPI',
      'Polyp11TOT',
      'Polyp12TOT',
      'Polyp1TOT',
      'Polyp2TOT',
      'Polyp3TOT',
      'Polyp4TOT',
      'Polyp5TOT',
      'Polyp6ATOT',
      'Polyp6BTOT',
      'Polyp7TOT',
      'Polyp8TOT',
      'Polyp9TOT'
    )
] <- 'Jose_Ordovas_Montanes_2018'
covid <- c(
  'AP1-NB',
  'AP10-NB',
  'AP11-NB',
  'AP14-NB_1',
  'AP14-NB_2',
  'AP4-NB',
  'AP5-NB',
  'AP7-NB',
  'AP8-NB',
  'AP9-NB',
  'PC11-NB',
  'PC12-NB',
  'PC2-NB',
  'PC6-NB',
  'PC9-NB',
  'PP1-NB',
  'PP10-NB',
  'PP11-NB v2.0',
  'PP12-NB',
  'PP13-NB',
  'PP15-NB',
  'PP16-NB',
  'PP17-NB',
  'PP18-NB',
  'PP19-NB',
  'PP2-NB',
  'PP3-NB',
  'PP4-NB',
  'PP5-NB_1',
  'PP5-NB_2',
  'PP6-NB_1',
  'PP6-NB_2',
  'PP7-NB_v1.1',
  'PP8-NB',
  'PP9-NB'
)
seurat_obj$dataset[seurat_obj$sample %in% covid] <- 'Kerstin_B_Meyer_2021_covid'
table(seurat_obj$disease__ontology_label, seurat_obj$polyp)
table(seurat_obj$disease_level_1, seurat_obj$dataset)
seurat_obj$disease_level_1[
  seurat_obj$dataset %in% c('Kerstin_B_Meyer_2021_covid')
] <- seurat_obj$disease[seurat_obj$dataset %in% c('Kerstin_B_Meyer_2021_covid')]
seurat_obj$disease_level_1[seurat_obj$polyp %in% 'YES'] <- 'CRSwNP'
seurat_obj$disease_level_1[seurat_obj$polyp %in% 'NO'] <- 'CRSsNP'
seurat_obj_polyp <- subset(
  seurat_obj,
  subset = GEO %in% c('GSE202100', 'GSE235711')
)

library(MSigDB)
