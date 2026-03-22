#!/usr/bin/env Rscript
# ==============================================================================
# Test: Multiple Database Input
# ==============================================================================

library(clusterProfiler)
library(fanyi)

cat("\n=== Loading Multiple Enrichment Results ===\n\n")

# Paths
BASE_DIR <- "/home/h2048/data/R/0127/bcell_interpret_v2_4_1/reports"

# Load available databases
go_bp_enrich <- readRDS(file.path(BASE_DIR, "go_bp_enrich.rds"))
cat("✓ GO BP loaded\n")

hallmark_enrich <- tryCatch(
  {
    h <- readRDS(file.path(BASE_DIR, "hallmark_enrich.rds"))
    cat("✓ Hallmark loaded\n")
    h
  },
  error = function(e) {
    cat("✗ Hallmark not found\n")
    return(NULL)
  }
)

bcell_markers_enrich <- tryCatch(
  {
    b <- readRDS(file.path(BASE_DIR, "bcell_markers_enrich.rds"))
    cat("✓ B Cell Markers loaded\n")
    b
  },
  error = function(e) {
    cat("✗ B Cell Markers not found\n")
    return(NULL)
  }
)

# Setup
ctx <- paste(
  "Human normal tissue B cells (non-diseased baseline).",
  "Goal: annotate B-cell subclusters with high precision.",
  "Integrate evidence from GO pathways, hallmark signatures, and B-cell markers."
)

Sys.setenv(DEEPSEEK_API_KEY = "sk-ed1879cf6fa14b04aac9cb6c078a3d05")
fanyi::set_translate_option(
  key = Sys.getenv("DEEPSEEK_API_KEY"),
  source = "deepseek"
)

MODEL <- "deepseek-chat"

# ==============================================================================
# Test 1: Single Database (baseline)
# ==============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Test 1: Single Database - list(go_bp_enrich)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

y_single <- tryCatch(
  {
    clusterProfiler::interpret(
      list(go_bp_enrich),
      task = "annotation",
      context = ctx,
      model = MODEL,
      n_pathways = 10
    )
  },
  error = function(e) {
    cat("[ERROR]:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (!is.null(y_single)) {
  cat(sprintf("[OK] Single DB: %d clusters\n", length(y_single)))
  cat(sprintf("Fields: %s\n", paste(names(y_single[[1]]), collapse = ", ")))

  # Check cell_type content
  if (!is.null(y_single[[1]]$cell_type)) {
    cat(sprintf("\nExample cell_type: %s\n", y_single[[1]]$cell_type))
  }
  if (!is.null(y_single[[1]]$markers)) {
    markers <- y_single[[1]]$markers
    if (length(markers) > 0) {
      cat(sprintf(
        "Example markers: %s\n",
        paste(head(markers, 3), collapse = ", ")
      ))
    }
  }
}

# ==============================================================================
# Test 2: Two Databases
# ==============================================================================

if (!is.null(hallmark_enrich)) {
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("Test 2: Two Databases - list(go_bp_enrich, hallmark_enrich)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  y_double <- tryCatch(
    {
      clusterProfiler::interpret(
        list(go_bp_enrich, hallmark_enrich),
        task = "annotation",
        context = ctx,
        model = MODEL,
        n_pathways = 10
      )
    },
    error = function(e) {
      cat("[ERROR]:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(y_double)) {
    cat(sprintf("[OK] Two DBs: %d clusters\n", length(y_double)))
    cat(sprintf("Fields: %s\n", paste(names(y_double[[1]]), collapse = ", ")))

    # Compare with single DB
    if (!is.null(y_single)) {
      cat("\n--- Comparison with Single DB ---\n")

      # Compare cell_type
      single_type <- y_single[[1]]$cell_type
      double_type <- y_double[[1]]$cell_type

      cat(sprintf("Single DB cell_type: %s\n", single_type))
      cat(sprintf("Two DBs cell_type:   %s\n", double_type))

      if (single_type != double_type) {
        cat("⚠ Different annotations!\n")
      } else {
        cat("✓ Same annotation\n")
      }

      # Compare reasoning length
      single_len <- nchar(y_single[[1]]$reasoning)
      double_len <- nchar(y_double[[1]]$reasoning)

      cat(sprintf("\nReasoning length:\n"))
      cat(sprintf("  Single DB: %d chars\n", single_len))
      cat(sprintf("  Two DBs:   %d chars\n", double_len))

      if (double_len > single_len) {
        cat(sprintf(
          "  → %d%% more detailed with 2 DBs\n",
          round((double_len - single_len) / single_len * 100)
        ))
      }
    }

    saveRDS(y_double, "test_double_db.rds")
    cat("\nSaved to: test_double_db.rds\n")
  }
}

# ==============================================================================
# Test 3: Three Databases
# ==============================================================================

if (!is.null(hallmark_enrich) && !is.null(bcell_markers_enrich)) {
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("Test 3: Three Databases - GO BP + Hallmark + B Cell Markers\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  y_triple <- tryCatch(
    {
      clusterProfiler::interpret(
        list(go_bp_enrich, hallmark_enrich, bcell_markers_enrich),
        task = "annotation",
        context = ctx,
        model = MODEL,
        n_pathways = 10
      )
    },
    error = function(e) {
      cat("[ERROR]:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (!is.null(y_triple)) {
    cat(sprintf("[OK] Three DBs: %d clusters\n", length(y_triple)))

    # Compare
    if (!is.null(y_single) && !is.null(y_double)) {
      cat("\n--- Comparison ---\n")

      triple_type <- y_triple[[1]]$cell_type
      triple_len <- nchar(y_triple[[1]]$reasoning)

      cat("Annotation:\n")
      cat(sprintf("  1 DB:  %s\n", y_single[[1]]$cell_type))
      cat(sprintf("  2 DBs: %s\n", y_double[[1]]$cell_type))
      cat(sprintf("  3 DBs: %s\n", triple_type))

      cat("\nReasoning length:\n")
      cat(sprintf("  1 DB:  %d chars\n", nchar(y_single[[1]]$reasoning)))
      cat(sprintf("  2 DBs: %d chars\n", nchar(y_double[[1]]$reasoning)))
      cat(sprintf("  3 DBs: %d chars\n", triple_len))
    }

    saveRDS(y_triple, "test_triple_db.rds")
    cat("\nSaved to: test_triple_db.rds\n")
  }
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Input Format Tests:\n")
cat(sprintf(
  "  list(go_bp):                    %s\n",
  ifelse(!is.null(y_single), "✓", "✗")
))

if (!is.null(hallmark_enrich)) {
  cat(sprintf(
    "  list(go_bp, hallmark):          %s\n",
    ifelse(exists("y_double") && !is.null(y_double), "✓", "✗")
  ))
}

if (!is.null(bcell_markers_enrich)) {
  cat(sprintf(
    "  list(go_bp, hallmark, bcell):   %s\n",
    ifelse(exists("y_triple") && !is.null(y_triple), "✓", "✗")
  ))
}

cat("\nRecommendation:\n")

if (exists("y_triple") && !is.null(y_triple)) {
  cat("  → ✓ Multi-database input works!\n")
  cat("  → Consider using 2-3 complementary databases:\n")
  cat("      - GO BP for biological processes\n")
  cat("      - Hallmark for canonical pathways\n")
  cat("      - B Cell Markers for cell type specificity\n")
} else if (exists("y_double") && !is.null(y_double)) {
  cat("  → ✓ Two-database input works!\n")
  cat("  → Recommend: GO BP + Hallmark or GO BP + B Cell Markers\n")
} else {
  cat("  → Use single database (GO BP)\n")
}

cat("\n")

#!/usr/bin/env Rscript
# ==============================================================================
# Test: 5-Database interpret() Call
# ==============================================================================
# Purpose: Test interpret() with GO BP, GO MF, GO CC, Hallmark, KEGG
# Date: 2026-01-27
# Status: Testing for production readiness
# ==============================================================================

library(clusterProfiler)
library(fanyi)

cat("\n")
cat(
  "================================================================================\n"
)
cat("5-Database interpret() Test\n")
cat(
  "================================================================================\n\n"
)

# ==============================================================================
# Load Enrichment Results
# ==============================================================================

BASE_DIR <- "/home/h2048/data/R/0127/bcell_interpret_v2_4_1/reports"

cat("Loading enrichment results...\n\n")

# Load all available databases
enrichment_list <- list()
db_names <- c()

# GO BP
if (file.exists(file.path(BASE_DIR, "go_bp_enrich.rds"))) {
  go_bp <- readRDS(file.path(BASE_DIR, "go_bp_enrich.rds"))
  enrichment_list[[length(enrichment_list) + 1]] <- go_bp
  db_names <- c(db_names, "GO BP")
  cat(sprintf(
    "✓ GO BP loaded: %d clusters, %d terms\n",
    length(unique(go_bp@compareClusterResult$Cluster)),
    nrow(go_bp@compareClusterResult)
  ))
}

# GO MF
if (file.exists(file.path(BASE_DIR, "go_mf_enrich.rds"))) {
  go_mf <- readRDS(file.path(BASE_DIR, "go_mf_enrich.rds"))
  enrichment_list[[length(enrichment_list) + 1]] <- go_mf
  db_names <- c(db_names, "GO MF")
  cat(sprintf(
    "✓ GO MF loaded: %d clusters, %d terms\n",
    length(unique(go_mf@compareClusterResult$Cluster)),
    nrow(go_mf@compareClusterResult)
  ))
}

# GO CC
if (file.exists(file.path(BASE_DIR, "go_cc_enrich.rds"))) {
  go_cc <- readRDS(file.path(BASE_DIR, "go_cc_enrich.rds"))
  enrichment_list[[length(enrichment_list) + 1]] <- go_cc
  db_names <- c(db_names, "GO CC")
  cat(sprintf(
    "✓ GO CC loaded: %d clusters, %d terms\n",
    length(unique(go_cc@compareClusterResult$Cluster)),
    nrow(go_cc@compareClusterResult)
  ))
}

# Hallmark
if (file.exists(file.path(BASE_DIR, "hallmark_enrich.rds"))) {
  hallmark <- readRDS(file.path(BASE_DIR, "hallmark_enrich.rds"))
  enrichment_list[[length(enrichment_list) + 1]] <- hallmark
  db_names <- c(db_names, "Hallmark")
  cat(sprintf(
    "✓ Hallmark loaded: %d clusters, %d terms\n",
    length(unique(hallmark@compareClusterResult$Cluster)),
    nrow(hallmark@compareClusterResult)
  ))
}

# KEGG
if (file.exists(file.path(BASE_DIR, "msigdb_kegg_enrich.rds"))) {
  kegg <- readRDS(file.path(BASE_DIR, "msigdb_kegg_enrich.rds"))
  enrichment_list[[length(enrichment_list) + 1]] <- kegg
  db_names <- c(db_names, "KEGG")
  cat(sprintf(
    "✓ KEGG loaded: %d clusters, %d terms\n",
    length(unique(kegg@compareClusterResult$Cluster)),
    nrow(kegg@compareClusterResult)
  ))
}

cat(sprintf("\nTotal databases loaded: %d\n", length(enrichment_list)))
cat(sprintf("Databases: %s\n\n", paste(db_names, collapse = ", ")))

if (length(enrichment_list) == 0) {
  stop("No enrichment objects found! Check BASE_DIR path.")
}

# ==============================================================================
# Setup API
# ==============================================================================

cat("Setting up DeepSeek API...\n")

DEEPSEEK_API_KEY <- "sk-ed1879cf6fa14b04aac9cb6c078a3d05"
MODEL <- "deepseek-chat"

Sys.setenv(DEEPSEEK_API_KEY = DEEPSEEK_API_KEY)
fanyi::set_translate_option(
  key = DEEPSEEK_API_KEY,
  source = "deepseek"
)

# Test connection
test_response <- tryCatch(
  {
    fanyi::chat_request("test", model = MODEL)
  },
  error = function(e) {
    cat("[ERROR] API connection failed:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (is.null(test_response)) {
  stop("API connection failed")
}

cat("[OK] API connection successful\n\n")

# ==============================================================================
# Context
# ==============================================================================

ctx <- paste(
  "Human normal respiratory tract B cells from nasal cavity, sinus, bronchi, and lung.",
  "Goal: Comprehensive cell type annotation using multi-database evidence.",
  "Databases: GO Biological Process (pathways), GO Molecular Function (activities),",
  "GO Cellular Component (localization), Hallmark gene sets (signatures), KEGG pathways (networks).",
  "B cell populations: Memory, Naive, Germinal Center, Plasma, Age-associated.",
  "Focus: Activation states, class-switching, proliferation, differentiation stages."
)

# ==============================================================================
# Test: Multi-Database Annotation
# ==============================================================================

SEP <- paste(rep("=", 80), collapse = "")

cat(SEP, "\n")
cat(sprintf("Running interpret() with %d databases\n", length(enrichment_list)))
cat(SEP, "\n\n")

cat("Databases:\n")
for (i in seq_along(db_names)) {
  cat(sprintf("  %d. %s\n", i, db_names[i]))
}
cat("\n")

start_time <- Sys.time()

annotation_results <- tryCatch(
  {
    interpret(
      x = enrichment_list,
      task = "annotation",
      context = ctx,
      model = MODEL,
      n_pathways = 15
    )
  },
  error = function(e) {
    cat("\n[ERROR] interpret() failed:\n")
    cat("  Message:", conditionMessage(e), "\n")
    cat("  Details:", e$message, "\n")
    return(NULL)
  }
)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

# ==============================================================================
# Check Results
# ==============================================================================

cat("\n")
cat(SEP, "\n")
cat("RESULTS\n")
cat(SEP, "\n\n")

if (is.null(annotation_results)) {
  cat("[FAILED] No results returned\n")
  cat("Possible reasons:\n")
  cat("  - API error\n")
  cat("  - Input format issue\n")
  cat("  - Model error\n")
  quit(status = 1)
}

cat(sprintf("[SUCCESS] interpret() completed in %.1f seconds\n\n", elapsed))

# Check structure
cat("Result structure:\n")
cat(sprintf("  Type: %s\n", class(annotation_results)))
cat(sprintf("  Length: %d clusters\n", length(annotation_results)))
cat(sprintf(
  "  Names: %s\n",
  paste(head(names(annotation_results), 3), collapse = ", ")
))

if (length(annotation_results) > 0) {
  cat("\nFirst cluster structure:\n")
  first_cluster <- annotation_results[[1]]
  cat(sprintf("  Cluster name: %s\n", names(annotation_results)[1]))
  cat(sprintf("  Fields: %s\n", paste(names(first_cluster), collapse = ", ")))

  # Check for expected fields
  cat("\nField validation:\n")

  expected_fields <- c(
    "cell_type",
    "confidence",
    "reasoning",
    "regulatory_drivers",
    "markers",
    "refined_network",
    "network_evidence",
    "cluster",
    "network"
  )

  for (field in expected_fields) {
    status <- if (!is.null(first_cluster[[field]])) "✓" else "✗"
    cat(sprintf("  %s %s\n", status, field))
  }

  # Sample content
  if (!is.null(first_cluster$cell_type)) {
    cat("\nSample content:\n")
    cat(sprintf("  Cell Type: %s\n", first_cluster$cell_type))
    cat(sprintf("  Confidence: %s\n", first_cluster$confidence))

    if (!is.null(first_cluster$markers) && length(first_cluster$markers) > 0) {
      cat(sprintf(
        "  Markers (%d total): %s\n",
        length(first_cluster$markers),
        paste(head(first_cluster$markers, 3), collapse = ", ")
      ))
    }

    if (
      !is.null(first_cluster$regulatory_drivers) &&
        length(first_cluster$regulatory_drivers) > 0
    ) {
      cat(sprintf(
        "  Regulatory Drivers: %s\n",
        paste(first_cluster$regulatory_drivers, collapse = ", ")
      ))
    }

    if (!is.null(first_cluster$reasoning)) {
      cat(sprintf(
        "  Reasoning (first 100 chars): %s...\n",
        substr(first_cluster$reasoning, 1, 100)
      ))
    }
  }

  # Save results
  output_file <- "test_5db_annotation_results.rds"
  saveRDS(annotation_results, output_file)
  cat(sprintf("\n✓ Results saved to: %s\n", output_file))

  # Generate summary CSV
  summary_df <- data.frame(
    Cluster = character(),
    Cell_Type = character(),
    Confidence = character(),
    Num_Markers = integer(),
    Num_Regulators = integer(),
    stringsAsFactors = FALSE
  )

  for (cluster_name in names(annotation_results)) {
    cluster_result <- annotation_results[[cluster_name]]

    # Handle special case (no enrichment)
    if (
      !is.null(cluster_result$overview) &&
        is.null(cluster_result$cell_type) &&
        cluster_result$confidence == "None"
    ) {
      summary_df <- rbind(
        summary_df,
        data.frame(
          Cluster = cluster_name,
          Cell_Type = "No enrichment",
          Confidence = "None",
          Num_Markers = 0,
          Num_Regulators = 0,
          stringsAsFactors = FALSE
        )
      )
      next
    }

    summary_df <- rbind(
      summary_df,
      data.frame(
        Cluster = cluster_name,
        Cell_Type = cluster_result$cell_type %||% "Unknown",
        Confidence = cluster_result$confidence %||% "NA",
        Num_Markers = length(cluster_result$markers %||% c()),
        Num_Regulators = length(cluster_result$regulatory_drivers %||% c()),
        stringsAsFactors = FALSE
      )
    )
  }

  write.csv(summary_df, "test_5db_summary.csv", row.names = FALSE)
  cat("✓ Summary saved to: test_5db_summary.csv\n")

  # Print summary table
  cat("\nAnnotation Summary:\n")
  print(summary_df)
}

# ==============================================================================
# Compare with 3-Database Results (if available)
# ==============================================================================

cat("\n")
cat(SEP, "\n")
cat("COMPARISON (if 3-DB results available)\n")
cat(SEP, "\n\n")

if (file.exists("test_triple_db.rds")) {
  cat("Loading 3-database results for comparison...\n")
  triple_results <- readRDS("test_triple_db.rds")

  if (!is.null(triple_results) && length(triple_results) > 0) {
    first_triple <- triple_results[[1]]
    first_five <- annotation_results[[1]]

    cat("\nFirst cluster comparison:\n")
    cat(sprintf("  3-DB cell_type: %s\n", first_triple$cell_type %||% "N/A"))
    cat(sprintf("  5-DB cell_type: %s\n", first_five$cell_type %||% "N/A"))

    reasoning_3db <- nchar(first_triple$reasoning %||% "")
    reasoning_5db <- nchar(first_five$reasoning %||% "")

    cat(sprintf("\nReasoning detail:\n"))
    cat(sprintf("  3-DB: %d characters\n", reasoning_3db))
    cat(sprintf("  5-DB: %d characters\n", reasoning_5db))

    if (reasoning_5db > reasoning_3db) {
      cat(sprintf(
        "  → 5-DB is %d%% more detailed\n",
        round((reasoning_5db - reasoning_3db) / reasoning_3db * 100)
      ))
    }
  }
} else {
  cat(
    "No 3-database results found for comparison (test_triple_db.rds not found)\n"
  )
}

# ==============================================================================
# Final Summary
# ==============================================================================

cat("\n")
cat(SEP, "\n")
cat("FINAL SUMMARY\n")
cat(SEP, "\n\n")

cat(sprintf("✓ Test completed successfully\n"))
cat(sprintf(
  "✓ Tested with %d databases: %s\n",
  length(enrichment_list),
  paste(db_names, collapse = ", ")
))
cat(sprintf("✓ Annotated %d clusters\n", length(annotation_results)))
cat(sprintf("✓ Runtime: %.1f seconds\n", elapsed))

cat("\nRecommendation:\n")
if (length(enrichment_list) >= 5) {
  cat("  → 5-database input works successfully!\n")
  cat("  → Consider using all 5 for production analysis\n")
  cat("  → Benefits: More comprehensive evidence, richer annotations\n")
} else {
  cat(sprintf("  → Only %d database(s) available\n", length(enrichment_list)))
  cat("  → Consider adding more enrichment databases\n")
}

cat("\nOutput files:\n")
cat("  - test_5db_annotation_results.rds (full results)\n")
cat("  - test_5db_summary.csv (summary table)\n")

cat("\n")
cat(SEP, "\n")
cat("TEST COMPLETE\n")
cat(SEP, "\n")
