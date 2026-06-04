#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

suppressPackageStartupMessages({
  library(Seurat)
})

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

make_tiny_seurat <- function() {
  set.seed(42)
  genes <- c('MS4A1','CD79A','CD79B','BANK1','HLA-DRA','CD74','MKI67','TOP2A')
  cells <- paste0('cell_', seq_len(12))
  counts <- matrix(
    rpois(length(genes) * length(cells), lambda = 5),
    nrow = length(genes),
    dimnames = list(genes, cells)
  )
  seu <- CreateSeuratObject(counts = counts)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu$cell_type_final_l3 <- rep(c('Naive_B', 'Memory_B'), each = 6)
  seu$sample <- rep(c('S1', 'S2', 'S3'), length.out = 12)
  seu
}

seu <- make_tiny_seurat()
CELLTYPE_COL <- 'cell_type_final_l3'
SAMPLE_COL <- 'sample'

.captured_metacells <- NULL
.captured_network <- NULL
.captured_normalize_metacells <- NULL

old_opts <- options(hdwgcna.allow_global_metacells_override = TRUE)

assign(
  'MetacellsByGroups',
  function(seurat_obj,
           group.by,
           ident.group,
           k,
           reduction = NULL,
           dims = NULL,
           assay = NULL,
           slot = NULL,
           layer = NULL,
           mode = NULL,
           cells.use = NULL,
           min_cells = NULL,
           max_shared = NULL,
           target_metacells = NULL,
           max_iter = NULL,
           verbose = NULL,
           wgcna_name = NULL) {
    .captured_metacells <<- list(
      group.by = group.by,
      ident.group = ident.group,
      k = k,
      reduction = reduction,
      dims = dims,
      min_cells = min_cells,
      max_shared = max_shared,
      target_metacells = target_metacells,
      wgcna_name = wgcna_name
    )
    seurat_obj
  },
  envir = .GlobalEnv
)

assign(
  'NormalizeMetacells',
  function(seurat_obj, scale.factor = NULL, wgcna_name = NULL, ...) {
    .captured_normalize_metacells <<- list(
      scale.factor = scale.factor,
      wgcna_name = wgcna_name
    )
    seurat_obj
  },
  envir = .GlobalEnv
)

assign(
  'ConstructNetwork',
  function(seurat_obj,
           soft_power,
           min_power = NULL,
           tom_outdir = NULL,
           tom_name = NULL,
           consensus = NULL,
           overwrite_tom = NULL,
           wgcna_name = NULL,
           blocks = NULL,
           maxBlockSize = NULL,
           randomSeed = NULL,
           corType = NULL,
           consensusQuantile = NULL,
           networkType = NULL,
           TOMType = NULL,
           TOMDenom = NULL,
           scaleTOMs = NULL,
           calibrationQuantile = NULL,
           sampleForCalibration = NULL,
           sampleForCalibrationFactor = NULL,
           useDiskCache = NULL,
           chunkSize = NULL,
           deepSplit = NULL,
           pamStage = NULL,
           detectCutHeight = NULL,
           minModuleSize = NULL,
           mergeCutHeight = NULL,
           saveConsensusTOMs = NULL,
           ...) {
    .captured_network <<- list(
      soft_power = soft_power,
      corType = corType,
      networkType = networkType,
      TOMType = TOMType,
      deepSplit = deepSplit,
      detectCutHeight = detectCutHeight,
      minModuleSize = minModuleSize,
      mergeCutHeight = mergeCutHeight,
      wgcna_name = wgcna_name
    )
    seurat_obj
  },
  envir = .GlobalEnv
)

assign('GetModules', function(seurat_obj, wgcna_name) data.frame(module = c('turquoise', 'grey'), gene_name = c('MS4A1', 'CD79A')), envir = .GlobalEnv)

on.exit({
  options(old_opts)
  rm('MetacellsByGroups', envir = .GlobalEnv)
  rm('NormalizeMetacells', envir = .GlobalEnv)
  rm('ConstructNetwork', envir = .GlobalEnv)
  rm('GetModules', envir = .GlobalEnv)
}, add = TRUE)

res_meta <- hdwgcna_metacells(
  seurat_obj = seu,
  group_by = c('cell_type_final_l3', 'sample'),
  metacell_k = 25,
  metacell_target = 50000,
  max_shared = 10,
  target_metacells = 250,
  metacell_reduction = 'harmony',
  metacell_dims = 1:5,
  metacell_min_cells = 20,
  wgcna_name = 'hdWGCNA'
)

assert_true(inherits(res_meta, 'Seurat'), 'hdwgcna_metacells should return a Seurat object')
assert_true(identical(.captured_metacells$group.by, c('cell_type_final_l3', 'sample')), 'MetacellsByGroups should receive group.by')
assert_true(identical(.captured_metacells$ident.group, 'cell_type_final_l3'), 'MetacellsByGroups should derive ident.group from group_by[1]')
assert_true(identical(.captured_metacells$max_shared, 10L), 'MetacellsByGroups should receive max_shared')
assert_true(identical(.captured_metacells$target_metacells, 250L), 'MetacellsByGroups should receive target_metacells')
assert_true(identical(.captured_metacells$reduction, 'harmony'), 'MetacellsByGroups should receive reduction')
assert_true(identical(.captured_metacells$dims, 1:5), 'MetacellsByGroups should receive dims')
assert_true(identical(.captured_metacells$min_cells, 20L), 'MetacellsByGroups should receive min_cells')
assert_true(identical(.captured_normalize_metacells$scale.factor, 50000), 'NormalizeMetacells should receive scale.factor')
assert_true(identical(.captured_normalize_metacells$wgcna_name, 'hdWGCNA'), 'NormalizeMetacells should receive wgcna_name')

hdwgcna_enable_construct_metacells_singleton_compat()
construct_metacells_body <- paste(
  deparse(body(get('ConstructMetacells', envir = asNamespace('hdWGCNA')))),
  collapse = '\n'
)
assert_true(grepl('nn_map\\[new_chosen, , drop = FALSE\\]', construct_metacells_body), 'ConstructMetacells singleton shim should preserve matrix shape for new_chosen subsets')
assert_true(grepl('nn_map\\[chosen, , drop = FALSE\\]', construct_metacells_body), 'ConstructMetacells singleton shim should preserve matrix shape for chosen subsets')

hdwgcna_enable_metacells_by_groups_null_filter_compat()
metacells_by_groups_body <- paste(
  deparse(body(get('MetacellsByGroups', envir = asNamespace('hdWGCNA')))),
  collapse = '\n'
)
assert_true(grepl('if \\(length\\(remove\\) >= 1\\)', metacells_by_groups_body), 'MetacellsByGroups NULL-filter shim should drop singleton NULL entries from metacell_list')

res_net <- hdwgcna_construct_network(
  seurat_obj = seu,
  soft_power = 8,
  net_type = 'signed hybrid',
  tom_type = 'signed',
  cor_type = 'bicor',
  deep_split = 4,
  detect_cut_height = 0.995,
  min_module_size = 20,
  merge_cut_height = 0.25,
  wgcna_name = 'hdWGCNA'
)

assert_true(inherits(res_net, 'Seurat'), 'hdwgcna_construct_network should return a Seurat object')
assert_true(identical(.captured_network$soft_power, 8), 'ConstructNetwork should receive soft power')
assert_true(identical(.captured_network$corType, 'bicor'), 'ConstructNetwork should receive corType')
assert_true(identical(.captured_network$networkType, 'signed hybrid'), 'ConstructNetwork should receive networkType')
assert_true(identical(.captured_network$TOMType, 'signed'), 'ConstructNetwork should receive TOMType')
assert_true(identical(.captured_network$deepSplit, 4L), 'ConstructNetwork should receive deepSplit')
assert_true(identical(.captured_network$detectCutHeight, 0.995), 'ConstructNetwork should receive detectCutHeight')
assert_true(identical(.captured_network$minModuleSize, 20L), 'ConstructNetwork should receive minModuleSize')
assert_true(identical(.captured_network$mergeCutHeight, 0.25), 'ConstructNetwork should receive mergeCutHeight')

cat('hdwgcna parameter passthrough regression test passed.\n')
