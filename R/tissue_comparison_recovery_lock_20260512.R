# Lightweight file-lock helpers for 20260511/20260512 tissue-comparison recovery wrappers.
# Purpose: allow parallel front-loading while the original sequential launcher is still active.

.tc_recovery_file_nonempty <- function(path) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path)) return(FALSE)
  if (!file.exists(path)) return(FALSE)
  size <- suppressWarnings(file.info(path)$size)
  isTRUE(!is.na(size) && size > 0)
}

.tc_recovery_all_nonempty <- function(paths) {
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  length(paths) > 0 && all(vapply(paths, .tc_recovery_file_nonempty, logical(1)))
}

.tc_recovery_lock_pid <- function(lock_path) {
  lines <- tryCatch(readLines(lock_path, warn = FALSE), error = function(e) character())
  if (length(lines) == 0) return(NA_integer_)
  pid_line <- grep("^pid=", lines, value = TRUE)
  pid_raw <- if (length(pid_line) > 0) sub("^pid=", "", pid_line[[1]]) else lines[[1]]
  suppressWarnings(as.integer(trimws(pid_raw)))
}

.tc_recovery_pid_alive <- function(pid) {
  if (is.na(pid) || pid <= 1L) return(FALSE)
  status <- suppressWarnings(system2("/bin/kill", c("-0", as.character(pid)), stdout = FALSE, stderr = FALSE))
  identical(status, 0L)
}

.tc_recovery_lock_is_ours <- function(lock_path, pid = Sys.getpid()) {
  identical(.tc_recovery_lock_pid(lock_path), as.integer(pid))
}

tc_recovery_register_lock <- function(output_dir,
                                      lock_name,
                                      completion_paths = character(),
                                      label = "tissue-comparison recovery",
                                      wait_interval_sec = 60,
                                      stale_after_sec = 48 * 3600,
                                      wait_timeout_sec = Inf) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock_path <- file.path(output_dir, lock_name)
  completion_paths <- as.character(completion_paths)

  if (.tc_recovery_all_nonempty(completion_paths)) {
    cat(sprintf("[INFO] %s already complete; exiting wrapper early.\n", label))
    cat(sprintf("[INFO] Completion paths: %s\n", paste(completion_paths, collapse = " | ")))
    quit(save = "no", status = 0, runLast = FALSE)
  }

  start_time <- Sys.time()
  repeat {
    if (!file.exists(lock_path)) break

    owner_pid <- .tc_recovery_lock_pid(lock_path)
    owner_alive <- .tc_recovery_pid_alive(owner_pid)
    lock_age <- suppressWarnings(as.numeric(difftime(Sys.time(), file.info(lock_path)$mtime, units = "secs")))
    lock_age <- ifelse(is.na(lock_age), Inf, lock_age)

    if (.tc_recovery_all_nonempty(completion_paths)) {
      cat(sprintf("[INFO] %s completed while waiting for lock; exiting wrapper early.\n", label))
      quit(save = "no", status = 0, runLast = FALSE)
    }

    if (!owner_alive || lock_age > stale_after_sec) {
      cat(sprintf(
        "[WARN] Removing stale %s lock: %s (owner_pid=%s, alive=%s, age_sec=%.0f)\n",
        label, lock_path, as.character(owner_pid), as.character(owner_alive), lock_age
      ))
      unlink(lock_path, force = TRUE)
      break
    }

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (is.finite(wait_timeout_sec) && elapsed > wait_timeout_sec) {
      stop(sprintf("Timed out waiting for %s lock: %s (owner_pid=%s)", label, lock_path, owner_pid))
    }

    cat(sprintf(
      "[INFO] %s lock is held by pid %s; waiting %s sec before recheck.\n",
      label, owner_pid, wait_interval_sec
    ))
    Sys.sleep(wait_interval_sec)
  }

  writeLines(
    c(
      sprintf("pid=%s", Sys.getpid()),
      sprintf("started=%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf("label=%s", label),
      sprintf("host=%s", Sys.info()[["nodename"]])
    ),
    lock_path
  )
  cat(sprintf("[OK] Acquired %s lock: %s\n", label, lock_path))

  force(lock_path)
  force(label)
  function() {
    if (file.exists(lock_path) && .tc_recovery_lock_is_ours(lock_path)) {
      unlink(lock_path, force = TRUE)
      cat(sprintf("[OK] Released %s lock: %s\n", label, lock_path))
    }
    invisible(TRUE)
  }
}
