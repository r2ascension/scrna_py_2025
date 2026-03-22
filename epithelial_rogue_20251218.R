# ===== ROGUE Analysis by Dataset and Cell Type =====
# Calculate cluster purity for each dataset-celltype combination

# Load required packages
suppressMessages(library(ROGUE))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

# ===== Configuration =====
INPUT_RDS <- "/home/h2048/data/R/1217/epithelial_bbknn_raw_20251217.rds"
OUTPUT_DIR <- "/home/h2048/data/R/1218/rogue_results/rogue_results"

DATASET_COL <- "dataset" # Dataset metadata column
CELLTYPE_COL <- "Manual_Annotation" # Cell type metadata column
MIN_CELLS_PER_GROUP <- 30 # Minimum cells for stable calculation

# ===== Load Data =====
seurat_obj <- readRDS(INPUT_RDS)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTPUT_DIR)
# ===== Get Valid Dataset-CellType Combinations =====
meta_data <- seurat_obj@meta.data

combo_counts <- meta_data %>%
  group_by(across(all_of(c(DATASET_COL, CELLTYPE_COL)))) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  filter(n_cells >= MIN_CELLS_PER_GROUP)

cat(sprintf("Total valid combinations: %d\n", nrow(combo_counts)))
cat(sprintf(
  "Datasets: %d, Cell types: %d\n",
  n_distinct(combo_counts[[DATASET_COL]]),
  n_distinct(combo_counts[[CELLTYPE_COL]])
))

# ===== Calculate ROGUE for Each Combination =====
cat("\n=== Calculating ROGUE values ===\n")

rogue_results_list <- list()
pb <- txtProgressBar(max = nrow(combo_counts), style = 3)

for (i in 1:nrow(combo_counts)) {
  dataset_id <- combo_counts[[DATASET_COL]][i]
  cell_type <- combo_counts[[CELLTYPE_COL]][i]
  n_cells <- combo_counts$n_cells[i]

  tryCatch(
    {
      # Subset cells for this dataset-celltype combination
      cells_keep <- meta_data[[DATASET_COL]] == dataset_id &
        meta_data[[CELLTYPE_COL]] == cell_type

      # Extract expression matrix
      expr_subset <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
      expr_subset <- expr_subset[, cells_keep]
      expr_subset <- as.matrix(expr_subset)

      # Filter low-abundance genes and cells
      expr_subset <- matr.filter(expr_subset, min.cells = 10, min.genes = 200)

      # Check matrix validity
      if (ncol(expr_subset) < MIN_CELLS_PER_GROUP || nrow(expr_subset) < 100) {
        rogue_val <- NA
      } else {
        # Calculate entropy
        ent_res <- SE_fun(expr_subset)

        # Check for invalid values
        if (any(is.na(ent_res$entropy)) || any(is.infinite(ent_res$entropy))) {
          rogue_val <- NA
        } else {
          # Calculate ROGUE
          rogue_val <- CalculateRogue(ent_res, platform = "UMI")
        }
      }

      rogue_results_list[[i]] <- data.frame(
        dataset = dataset_id,
        cell_type = cell_type,
        n_cells = n_cells,
        n_genes_used = ifelse(!is.na(rogue_val), nrow(expr_subset), NA),
        rogue_value = rogue_val
      )

      # Clean up
      rm(expr_subset)
      if (exists("ent_res")) {
        rm(ent_res)
      }
      gc(verbose = FALSE)
    },
    error = function(e) {
      rogue_results_list[[i]] <- data.frame(
        dataset = dataset_id,
        cell_type = cell_type,
        n_cells = n_cells,
        n_genes_used = NA,
        rogue_value = NA
      )
    }
  )

  setTxtProgressBar(pb, i)
}
close(pb)

# ===== Compile Results =====
rogue_df <- bind_rows(rogue_results_list)

# Report failed calculations
n_failed <- sum(is.na(rogue_df$rogue_value))
if (n_failed > 0) {
  cat(sprintf(
    "\n⚠️  %d/%d combinations failed (%.1f%%)\n",
    n_failed,
    nrow(rogue_df),
    100 * n_failed / nrow(rogue_df)
  ))
}

# Save long format
write.csv(
  rogue_df,
  file.path(OUTPUT_DIR, "rogue_values_long_format.csv"),
  row.names = FALSE
)

# ===== Convert to Wide Format (Matrix) =====
rogue_wide <- rogue_df %>%
  select(dataset, cell_type, rogue_value) %>%
  pivot_wider(names_from = cell_type, values_from = rogue_value)

write.csv(
  rogue_wide,
  file.path(OUTPUT_DIR, "rogue_values_wide_format.csv"),
  row.names = FALSE
)

# ===== Summary Statistics by Cell Type =====
celltype_summary <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  group_by(cell_type) %>%
  summarise(
    n_datasets = n(),
    total_cells = sum(n_cells),
    mean_rogue = mean(rogue_value),
    sd_rogue = sd(rogue_value),
    median_rogue = median(rogue_value),
    min_rogue = min(rogue_value),
    max_rogue = max(rogue_value),
    .groups = 'drop'
  ) %>%
  arrange(desc(median_rogue))

write.csv(
  celltype_summary,
  file.path(OUTPUT_DIR, "rogue_summary_by_celltype.csv"),
  row.names = FALSE
)

# ===== Summary Statistics by Dataset =====
dataset_summary <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  group_by(dataset) %>%
  summarise(
    n_celltypes = n(),
    total_cells = sum(n_cells),
    mean_rogue = mean(rogue_value),
    median_rogue = median(rogue_value),
    .groups = 'drop'
  ) %>%
  arrange(desc(median_rogue))

write.csv(
  dataset_summary,
  file.path(OUTPUT_DIR, "rogue_summary_by_dataset.csv"),
  row.names = FALSE
)

# ===== Visualizations =====

# 1. Boxplot by Cell Type (across datasets)
pdf(
  file.path(OUTPUT_DIR, "rogue_boxplot_by_celltype.pdf"),
  width = 12,
  height = 6
)
p1 <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  ggplot(aes(x = reorder(cell_type, rogue_value, median), y = rogue_value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red") +
  labs(
    title = "ROGUE Distribution by Cell Type (Across Datasets)",
    x = "Cell Type",
    y = "ROGUE Value",
    subtitle = sprintf("n = %d datasets", n_distinct(rogue_df$dataset))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  coord_cartesian(ylim = c(0.5, 1.0))
print(p1)
dev.off()

# 2. Boxplot by Dataset (across cell types)
pdf(
  file.path(OUTPUT_DIR, "rogue_boxplot_by_dataset.pdf"),
  width = 12,
  height = 6
)
p2 <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  ggplot(aes(x = reorder(dataset, rogue_value, median), y = rogue_value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightcoral", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red") +
  labs(
    title = "ROGUE Distribution by Dataset (Across Cell Types)",
    x = "Dataset",
    y = "ROGUE Value",
    subtitle = sprintf("n = %d cell types", n_distinct(rogue_df$cell_type))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  coord_cartesian(ylim = c(0.5, 1.0))
print(p2)
dev.off()

# 3. Heatmap
pdf(file.path(OUTPUT_DIR, "rogue_heatmap.pdf"), width = 14, height = 10)
rogue_matrix <- rogue_df %>%
  select(dataset, cell_type, rogue_value) %>%
  pivot_wider(names_from = cell_type, values_from = rogue_value) %>%
  column_to_rownames("dataset") %>%
  as.matrix()

if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(
    rogue_matrix,
    color = colorRampPalette(c("red", "yellow", "green"))(100),
    breaks = seq(0.5, 1.0, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    na_col = "grey90",
    main = "ROGUE Values: Dataset × Cell Type",
    fontsize = 8,
    cellwidth = 15,
    cellheight = 12
  )
} else {
  cat("pheatmap package not available, skipping heatmap\n")
}
dev.off()

# 4. Faceted plot: each cell type across datasets
pdf(
  file.path(OUTPUT_DIR, "rogue_faceted_by_celltype.pdf"),
  width = 14,
  height = 10
)
p3 <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  ggplot(aes(x = dataset, y = rogue_value, fill = dataset)) +
  geom_col() +
  geom_hline(
    yintercept = c(0.7, 0.9),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
  labs(
    title = "ROGUE Values by Cell Type and Dataset",
    x = "Dataset",
    y = "ROGUE Value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold")
  ) +
  coord_cartesian(ylim = c(0.5, 1.0))
print(p3)
dev.off()

# 5. Dot plot: cell type × dataset
pdf(file.path(OUTPUT_DIR, "rogue_dotplot.pdf"), width = 14, height = 8)
p4 <- rogue_df %>%
  filter(!is.na(rogue_value)) %>%
  ggplot(aes(x = cell_type, y = dataset, size = n_cells, color = rogue_value)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(
    low = "red",
    mid = "yellow",
    high = "green",
    midpoint = 0.8,
    limits = c(0.5, 1.0)
  ) +
  scale_size_continuous(range = c(2, 10), labels = scales::comma) +
  labs(
    title = "ROGUE Values: Cell Type × Dataset",
    x = "Cell Type",
    y = "Dataset",
    size = "Number of Cells",
    color = "ROGUE Value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )
print(p4)
dev.off()

# ===== Print Summaries =====
cat("\n=== Summary by Cell Type ===\n")
print(celltype_summary, n = Inf)

cat("\n=== Summary by Dataset ===\n")
print(dataset_summary, n = Inf)

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
cat(sprintf(
  "Successfully calculated: %d/%d combinations\n",
  sum(!is.na(rogue_df$rogue_value)),
  nrow(rogue_df)
))
