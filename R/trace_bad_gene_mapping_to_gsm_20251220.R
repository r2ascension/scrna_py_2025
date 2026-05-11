#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# =======================
# CONFIG (按你实际路径改)
# =======================
MAPPING_CSV <- "/home/h2048/data/R/1215/merge/cleaned_samples_standardized_output/global_gene_mapping/global_gene_mapping_complete.csv"
CLEANED_RDS_DIR <- "/home/h2048/data/R/1215/merge/cleaned_samples"  # 用于追踪“在哪个GSM出现”
ASSAY <- "RNA"

OUT_DIR <- "/home/h2048/data/R/1215/merge/cleaned_samples_standardized_output/global_gene_mapping/_trace_back"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 可选：把你那批“怀疑没正确对应”的ENSG列表放到这里（支持：
#  - 每行一个基因ID：ENSG00000...
#  - 或者你现在粘贴的那种三列CSV行：ENSG...,ENSG...,SYMBOL）
GENE_LIST_FILE <- file.path(dirname(MAPPING_CSV), "genes_to_check.txt")

# =======================
# Helpers
# =======================
read_mapping_csv <- function(path) {
  dt <- fread(path, header = TRUE, sep = ",", showProgress = FALSE)
  # 兼容无表头的情况（你贴出来看像无表头的3列）
  if (ncol(dt) == 3 && !all(c("original_name","official_symbol","mapping_source") %in% names(dt))) {
    setnames(dt, c("original_name","official_symbol","mapping_source"))
  }
  if (!all(c("original_name","official_symbol","mapping_source") %in% names(dt))) {
    stop("Mapping CSV must contain columns: original_name, official_symbol, mapping_source")
  }
  dt[]
}

parse_gene_list_file <- function(path) {
  if (!file.exists(path)) return(character())
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  # 如果是三列CSV行，就取第一列
  x <- sub(",.*$", "", x)
  x <- unique(x)
  x
}

get_features_strict <- function(obj, assay = "RNA") {
  tryCatch(
    Features(obj, assay = assay),
    error = function(e) rownames(obj[[assay]])
  )
}

map_ensembl_to_symbol <- function(ensembl_ids) {
  ensembl_ids <- unique(ensembl_ids)
  ensembl_ids <- ensembl_ids[grepl("^ENSG[0-9]+$", ensembl_ids)]
  if (length(ensembl_ids) == 0) return(setNames(character(), character()))
  m <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  as.character(m)
}

# =======================
# 1) Load mapping + pick targets
# =======================
cat("Loading mapping CSV...\n")
map_dt <- read_mapping_csv(MAPPING_CSV)

deprecated_genes <- unique(map_dt[grepl("^DEPRECATED-ENSG[0-9]+$", original_name), original_name])

user_genes <- parse_gene_list_file(GENE_LIST_FILE)

# 如果你不放 genes_to_check.txt，则只追踪 DEPRECATED-ENSG 那批
target_genes <- unique(c(deprecated_genes, user_genes))

cat(sprintf("Targets:\n  DEPRECATED-ENSG*: %d\n  user genes file: %d\n  total unique targets: %d\n",
            length(deprecated_genes), length(user_genes), length(target_genes)))

if (length(target_genes) == 0) {
  stop("No target genes found. Provide DEPRECATED-ENSG entries in mapping OR create genes_to_check.txt")
}

# =======================
# 2) Trace back to GSM (scan per-sample RDS features)
# =======================
cat("Scanning sample RDS files to trace targets back to GSM...\n")
rds_files <- list.files(CLEANED_RDS_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS found in CLEANED_RDS_DIR")

sample_names <- sub("\\.rds$", "", basename(rds_files))

hits_list <- vector("list", length(rds_files))
names(hits_list) <- sample_names

for (i in seq_along(rds_files)) {
  sn <- sample_names[i]
  cat(sprintf("  [%d/%d] %s\n", i, length(rds_files), sn))
  obj <- readRDS(rds_files[i])
  feats <- get_features_strict(obj, assay = ASSAY)
  hits <- intersect(feats, target_genes)
  if (length(hits) > 0) {
    hits_list[[i]] <- data.table(sample = sn, original_name = hits)
  } else {
    hits_list[[i]] <- NULL
  }
  rm(obj); gc()
}

hits_dt <- rbindlist(hits_list, use.names = TRUE, fill = TRUE)
if (nrow(hits_dt) == 0) {
  warning("No targets found in any sample RDS features. Check CLEANED_RDS_DIR / ASSAY / whether features were already standardized.")
}

# 输出：每个 gene 在哪些 sample 出现
gene_to_samples <- hits_dt[, .(
  n_samples = uniqueN(sample),
  samples = paste(sort(unique(sample)), collapse = ";")
), by = original_name][order(-n_samples, original_name)]

# 输出：每个 sample 命中多少 target gene
sample_summary <- hits_dt[, .(
  n_targets_hit = uniqueN(original_name),
  targets = paste(sort(unique(original_name)), collapse = ";")
), by = sample][order(-n_targets_hit, sample)]

fwrite(hits_dt, file.path(OUT_DIR, "target_gene__sample_hits.tsv"), sep = "\t")
fwrite(gene_to_samples, file.path(OUT_DIR, "target_gene__to_samples.tsv"), sep = "\t")
fwrite(sample_summary, file.path(OUT_DIR, "sample__target_gene_summary.tsv"), sep = "\t")

cat(sprintf("Trace-back outputs written to: %s\n", OUT_DIR))

# =======================
# 3) Validate whether mapping is “stuck” at SYMBOL but ENSEMBL can map to real symbol
# =======================
cat("Validating ENSEMBL->SYMBOL mapping for targets...\n")

sub_dt <- map_dt[original_name %in% target_genes]
sub_dt[, base_ensembl := sub("^DEPRECATED-", "", original_name)]  # DEPRECATED-ENSG... -> ENSG...

ensembl_symbols <- map_ensembl_to_symbol(sub_dt$base_ensembl)
sub_dt[, ensembl_symbol := unname(ensembl_symbols[base_ensembl])]

# 关键标记：csv里看起来“没变”（official_symbol == original_name）
# 但 ENSEMBL 实际能映射到一个不同的 SYMBOL
sub_dt[, likely_wrong := FALSE]
sub_dt[
  !is.na(ensembl_symbol) &
    nzchar(ensembl_symbol) &
    (official_symbol == original_name | is.na(official_symbol)) &
    (ensembl_symbol != base_ensembl),
  likely_wrong := TRUE
]

# 也把“DEPRECATED-ENSG”但 ENSEMBL 可映射的单独标注
sub_dt[, deprecated_but_mappable := grepl("^DEPRECATED-ENSG", original_name) & !is.na(ensembl_symbol) & nzchar(ensembl_symbol)]

fwrite(sub_dt, file.path(OUT_DIR, "target_genes__mapping_validation.tsv"), sep = "\t")

cat("Done.\n")
cat("Key files:\n")
cat(sprintf("  - %s\n", file.path(OUT_DIR, "target_gene__to_samples.tsv")))
cat(sprintf("  - %s\n", file.path(OUT_DIR, "target_genes__mapping_validation.tsv")))
