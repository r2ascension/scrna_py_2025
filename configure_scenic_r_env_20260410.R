#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Configure / validate SCENIC R environment
# ==============================================================================
# Usage:
#   Rscript configure_scenic_r_env_20260410.R
#   SCENIC_CHECK_ONLY=true Rscript configure_scenic_r_env_20260410.R
#
# Notes:
#   - Prefers system R libraries already present on this workstation
#   - Installs only missing packages
#   - Verifies cisTarget feather databases required by scenic_core_20260410.R
# ==============================================================================

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = max(600L, getOption("timeout", 60L)),
  scipen = 999
)

CHECK_ONLY <- tolower(Sys.getenv("SCENIC_CHECK_ONLY", "false")) %in% c("1", "true", "yes")
DATABASE_DIR <- if (dir.exists("/home/h2048/data/index_genome/cisTarget_databases_rscenic")) {
  "/home/h2048/data/index_genome/cisTarget_databases_rscenic"
} else {
  "/home/h2048/data/index_genome/cisTarget_databases"
}

cran_packages <- c(
  "Seurat",
  "data.table",
  "dplyr",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "visNetwork",
  "htmlwidgets"
)

bioc_packages <- c(
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "BiocParallel",
  "GSEABase"
)

github_packages <- c(
  "aertslab/SCopeLoomR" = "SCopeLoomR",
  "aertslab/SCENIC" = "SCENIC"
)

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

ensure_cran_package <- function(pkg) {
  if (is_installed(pkg)) return(invisible(TRUE))
  if (CHECK_ONLY) return(invisible(FALSE))
  install.packages(pkg)
  invisible(is_installed(pkg))
}

ensure_biocmanager <- function() {
  if (!is_installed("BiocManager")) {
    if (CHECK_ONLY) return(invisible(FALSE))
    install.packages("BiocManager")
  }
  invisible(is_installed("BiocManager"))
}

ensure_bioc_package <- function(pkg) {
  if (is_installed(pkg)) return(invisible(TRUE))
  if (!ensure_biocmanager()) return(invisible(FALSE))
  if (CHECK_ONLY) return(invisible(FALSE))
  BiocManager::install(pkg, ask = FALSE, update = FALSE)
  invisible(is_installed(pkg))
}

ensure_github_installer <- function() {
  if (!is_installed("remotes")) {
    if (CHECK_ONLY) return(invisible(FALSE))
    install.packages("remotes")
  }
  invisible(is_installed("remotes"))
}

ensure_github_package <- function(repo, pkg) {
  if (is_installed(pkg)) return(invisible(TRUE))
  if (!ensure_github_installer()) return(invisible(FALSE))
  if (CHECK_ONLY) return(invisible(FALSE))
  remotes::install_github(repo, upgrade = "never", dependencies = FALSE)
  invisible(is_installed(pkg))
}

cat("============================================================\n")
cat("SCENIC R environment preflight\n")
cat("============================================================\n")
cat(sprintf("R.version : %s\n", R.version.string))
cat(sprintf("R.home    : %s\n", R.home()))
cat(sprintf("Check only: %s\n", CHECK_ONLY))
cat(sprintf("DB dir    : %s\n", DATABASE_DIR))

cran_status <- vapply(cran_packages, ensure_cran_package, logical(1))
bioc_status <- vapply(bioc_packages, ensure_bioc_package, logical(1))
github_status <- vapply(seq_along(github_packages), function(i) {
  ensure_github_package(names(github_packages)[i], unname(github_packages[i]))
}, logical(1))
names(github_status) <- unname(github_packages)

cat("\n[CRAN]\n")
for (pkg in names(cran_status)) {
  cat(sprintf("- %s : %s\n", pkg, cran_status[[pkg]]))
}

cat("\n[Bioconductor]\n")
for (pkg in names(bioc_status)) {
  cat(sprintf("- %s : %s\n", pkg, bioc_status[[pkg]]))
}

cat("\n[GitHub]\n")
for (pkg in names(github_status)) {
  cat(sprintf("- %s : %s\n", pkg, github_status[[pkg]]))
}

required_db_patterns <- c("10kb", "500bp")
db_files <- list.files(DATABASE_DIR, pattern = "\\.feather$", full.names = TRUE)
db_ok <- length(db_files) >= 2 && all(vapply(required_db_patterns, function(pat) {
  any(grepl(pat, basename(db_files), ignore.case = TRUE))
}, logical(1)))

cat("\n[Databases]\n")
cat(sprintf("- feather files found : %d\n", length(db_files)))
cat(sprintf("- 10kb+500bp ready    : %s\n", db_ok))
if (length(db_files) > 0) {
  cat(sprintf("- examples            : %s\n", paste(utils::head(basename(db_files), 2), collapse = ", ")))
}

all_pkg_ok <- all(cran_status) && all(bioc_status) && all(github_status)
all_ok <- all_pkg_ok && db_ok

cat("\n============================================================\n")
cat(sprintf("SCENIC R environment ready: %s\n", all_ok))
cat("============================================================\n")

if (!all_ok) {
  quit(save = "no", status = 1)
}
