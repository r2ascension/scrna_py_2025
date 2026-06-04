#!/usr/bin/env Rscript

source('/home/h2048/script/R/hdwgcna_covarnet_helpers_v1_1.R')

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

ref_modules <- data.frame(
  module = c('blue', 'blue', 'turquoise', 'turquoise'),
  gene_name = c('G1', 'G2', 'G3', 'G4'),
  stringsAsFactors = FALSE
)
query_modules <- data.frame(
  module = c('blue', 'brown', 'brown', 'turquoise'),
  gene_name = c('G1', 'G2', 'G5', 'G4'),
  stringsAsFactors = FALSE
)

ov <- hdwgcna_pairwise_module_overlap(
  reference_modules = ref_modules,
  query_modules = query_modules,
  reference_label = 'nose',
  query_label = 'bronchial'
)
assert_true(is.data.frame(ov), 'module overlap should return a data.frame')
assert_true(any(ov$reference_module == 'blue' & ov$query_module == 'blue' & ov$n_overlap == 1L), 'module overlap should detect overlapping blue modules')
assert_true(any(ov$reference_module == 'turquoise' & ov$query_module == 'turquoise' & ov$n_overlap == 1L), 'module overlap should detect overlapping turquoise modules')

ref_hubs <- data.frame(
  module = c('blue', 'blue', 'turquoise'),
  gene_name = c('G1', 'G2', 'G4'),
  hub_score = c(3, 2, 5),
  stringsAsFactors = FALSE
)
query_hubs <- data.frame(
  module = c('blue', 'brown', 'turquoise'),
  gene_name = c('G1', 'G2', 'G4'),
  hub_score = c(4, 1, 6),
  stringsAsFactors = FALSE
)

hub_ov <- hdwgcna_pairwise_hub_overlap(
  reference_hubs = ref_hubs,
  query_hubs = query_hubs,
  reference_label = 'nose',
  query_label = 'bronchial',
  top_n = 2
)
assert_true(is.data.frame(hub_ov), 'hub overlap should return a data.frame')
assert_true(any(hub_ov$reference_module == 'blue' & hub_ov$query_module == 'blue' & hub_ov$n_overlap_hubs == 1L), 'hub overlap should detect overlapping blue hubs')

mock_preservation <- list(
  preservation = list(
    observed = list(
      'nose.nose' = list(
        'inColumnsAlsoPresentIn.bronchial' = data.frame(
          moduleSize = c(12, 20),
          medianRank.pres = c(1.5, 2.0),
          row.names = c('blue', 'turquoise'),
          check.names = FALSE
        )
      )
    ),
    Z = list(
      'nose.nose' = list(
        'inColumnsAlsoPresentIn.bronchial' = data.frame(
          moduleSize = c(12, 20),
          Zsummary.pres = c(6.2, 1.9),
          row.names = c('blue', 'turquoise'),
          check.names = FALSE
        )
      )
    )
  )
)

pres_df <- hdwgcna_tidy_module_preservation(
  module_preservation_result = mock_preservation,
  reference_label = 'nose',
  query_label = 'bronchial'
)
assert_true(is.data.frame(pres_df), 'tidy preservation should return a data.frame')
assert_true(identical(sort(pres_df$module), c('blue', 'turquoise')), 'tidy preservation should keep module names from rownames')
assert_true(any(pres_df$module == 'blue' & abs(pres_df$Zsummary.pres - 6.2) < 1e-8), 'tidy preservation should keep Zsummary.pres values')
assert_true(all(pres_df$reference_network == 'nose'), 'tidy preservation should annotate reference network label')
assert_true(all(pres_df$query_network == 'bronchial'), 'tidy preservation should annotate query network label')

cat('hdwgcna tissue-preservation helper regression test passed.\n')
