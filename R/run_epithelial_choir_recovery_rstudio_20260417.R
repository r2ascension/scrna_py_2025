#!/usr/bin/env Rscript
# ==============================================================================
# RStudio launcher for epithelial CHOIR recovery
# Open this file in RStudio and click Source (or run source() in Console).
# ==============================================================================

ROOT_DIR <- "/home/h2048"
RECOVERY_SCRIPT <- "/home/h2048/script/R/epithelial_tissue_comparison_v1_3_2_choir_recovery_20260416.R"
FORCE_START <- identical(Sys.getenv("FORCE_START", unset = "0"), "1")

message(sprintf("[%s] RStudio launcher started", format(Sys.time(), "%F %T")))
message(sprintf("[%s] ROOT_DIR=%s", format(Sys.time(), "%F %T"), ROOT_DIR))
message(sprintf("[%s] RECOVERY_SCRIPT=%s", format(Sys.time(), "%F %T"), RECOVERY_SCRIPT))

if (!file.exists(RECOVERY_SCRIPT)) {
  stop(sprintf("Recovery script not found: %s", RECOVERY_SCRIPT), call. = FALSE)
}

find_existing_pids <- function(pattern) {
  out <- tryCatch(
    system2("pgrep", c("-f", pattern), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  out <- trimws(out)
  out[nzchar(out)]
}

existing_pids <- find_existing_pids(RECOVERY_SCRIPT)
if (length(existing_pids) > 0 && !FORCE_START) {
  stop(
    paste0(
      "Detected existing epithelial recovery process(es): ",
      paste(existing_pids, collapse = ", "),
      "\nRefusing to start a duplicate run.\n",
      "If you really want to bypass this guard, set Sys.setenv(FORCE_START='1') and rerun."
    ),
    call. = FALSE
  )
}

setwd(ROOT_DIR)
message(sprintf("[%s] Working directory set to %s", format(Sys.time(), "%F %T"), getwd()))
message(sprintf("[%s] Sourcing recovery script...", format(Sys.time(), "%F %T")))

source(RECOVERY_SCRIPT, echo = TRUE, chdir = FALSE)

message(sprintf("[%s] Launcher finished", format(Sys.time(), "%F %T")))
