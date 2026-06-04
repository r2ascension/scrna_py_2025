#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(data.table)
	library(igraph)
})

DEFAULTS <- list(
	deseq_root = "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/pseudobulk_de_L3",
	source_dir = "/home/h2048/data/source/reference/ppi",
	output_dir = "/home/h2048/data/R/20260531/epithelial_airway_deg_ppi_workflow_20260531",
	top_n = 150L,
	padj_threshold = 0.05,
	abs_lfc_threshold = 1,
	string_score_threshold = 700L,
	min_consensus_support = 2L,
	target_pairs = NULL
)

print_usage <- function(defaults) {
	cat(
		paste(
			"Usage:",
			"  Rscript epithelial_airway_deg_ppi_workflow_20260531.R [--deseq-root PATH] [--source-dir PATH] [--output-dir PATH]",
			"         [--top-n INT] [--padj-threshold FLOAT] [--abs-lfc-threshold FLOAT]",
			"         [--string-score-threshold INT] [--min-consensus-support INT] [--target-pairs PAIR1,PAIR2]",
			"",
			"Defaults:",
			paste0("  --deseq-root ", defaults$deseq_root),
			paste0("  --source-dir ", defaults$source_dir),
			paste0("  --output-dir ", defaults$output_dir),
			paste0("  --top-n ", defaults$top_n),
			paste0("  --padj-threshold ", defaults$padj_threshold),
			paste0("  --abs-lfc-threshold ", defaults$abs_lfc_threshold),
			paste0("  --string-score-threshold ", defaults$string_score_threshold),
			paste0("  --min-consensus-support ", defaults$min_consensus_support),
			sep = "\n"
		),
		"\n"
	)
}

parse_args <- function(defaults) {
	args <- commandArgs(trailingOnly = TRUE)
	opts <- defaults
	if (!length(args)) {
		return(opts)
	}
	if (any(args %in% c("-h", "--help"))) {
		print_usage(defaults)
		quit(save = "no", status = 0)
	}

	i <- 1L
	while (i <= length(args)) {
		key <- args[[i]]
		if (!startsWith(key, "--")) {
			stop("Unexpected argument: ", key)
		}
		if (i == length(args)) {
			stop("Missing value for argument: ", key)
		}
		value <- args[[i + 1L]]
		key <- sub("^--", "", key)
		switch(
			key,
			"deseq-root" = opts$deseq_root <- value,
			"source-dir" = opts$source_dir <- value,
			"output-dir" = opts$output_dir <- value,
			"top-n" = opts$top_n <- as.integer(value),
			"padj-threshold" = opts$padj_threshold <- as.numeric(value),
			"abs-lfc-threshold" = opts$abs_lfc_threshold <- as.numeric(value),
			"string-score-threshold" = opts$string_score_threshold <- as.integer(value),
			"min-consensus-support" = opts$min_consensus_support <- as.integer(value),
			"target-pairs" = {
				value <- trimws(value)
				opts$target_pairs <- if (!nzchar(value)) NULL else trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
			},
			stop("Unknown argument: --", key)
		)
		i <- i + 2L
	}
	opts
}

validate_config <- function(config) {
	if (is.na(config$top_n) || config$top_n < 1L) {
		stop("--top-n must be a positive integer")
	}
	if (is.na(config$padj_threshold) || config$padj_threshold <= 0 || config$padj_threshold > 1) {
		stop("--padj-threshold must be in (0, 1]")
	}
	if (is.na(config$abs_lfc_threshold) || config$abs_lfc_threshold < 0) {
		stop("--abs-lfc-threshold must be >= 0")
	}
	if (is.na(config$string_score_threshold) || config$string_score_threshold < 0) {
		stop("--string-score-threshold must be >= 0")
	}
	if (is.na(config$min_consensus_support) || config$min_consensus_support < 1L) {
		stop("--min-consensus-support must be a positive integer")
	}
	config
}

dir_create <- function(path) {
	dir.create(path, recursive = TRUE, showWarnings = FALSE)
	invisible(path)
}

normalize_token <- function(x) {
	x <- trimws(tolower(as.character(x)))
	x <- gsub("[^a-z0-9]+", "_", x, perl = TRUE)
	x <- gsub("^_+|_+$", "", x, perl = TRUE)
	x
}

site_group_from_token <- function(x) {
	token <- normalize_token(x)
	if (!nzchar(token)) {
		return(token)
	}
	if (grepl("sinus|nose|nasal", token, ignore.case = TRUE, perl = TRUE)) {
		return("nasal")
	}
	if (grepl("bronch|trachea|airway", token, ignore.case = TRUE, perl = TRUE)) {
		return("bronchial")
	}
	if (grepl("lung|alveol|parench", token, ignore.case = TRUE, perl = TRUE)) {
		return("lung")
	}
	token
}

canonical_pair <- function(left, right) {
	values <- c(normalize_token(left), normalize_token(right))
	pair_rank <- c(nasal = 1L, bronchial = 2L, lung = 3L)
	ranks <- pair_rank[values]
	ranks[is.na(ranks)] <- 99L
	ord <- order(ranks, values)
	paste(values[ord], collapse = "_vs_")
}

safe_numeric <- function(x) {
	suppressWarnings(as.numeric(x))
}

pick_col <- function(column_names, patterns, required = TRUE) {
	for (pattern in patterns) {
		hit <- column_names[grepl(pattern, column_names, ignore.case = TRUE, perl = TRUE)]
		if (length(hit)) {
			return(hit[[1L]])
		}
	}
	if (required) {
		stop("Could not find column matching patterns: ", paste(patterns, collapse = " | "))
	}
	NA_character_
}

empty_edge_table <- function() {
	data.table(
		gene_a = character(),
		gene_b = character(),
		weight = numeric(),
		evidence_count = integer(),
		source_db = character()
	)
}

canonicalize_edges <- function(dt, gene_a_col, gene_b_col, db_name, weight_col = NULL, evidence_col = NULL) {
	if (!nrow(dt)) {
		return(empty_edge_table())
	}
	out <- copy(dt)
	out[, gene_a := toupper(trimws(as.character(get(gene_a_col))))]
	out[, gene_b := toupper(trimws(as.character(get(gene_b_col))))]
	out <- out[nzchar(gene_a) & nzchar(gene_b) & gene_a != gene_b]
	if (!nrow(out)) {
		return(empty_edge_table())
	}
	out[, node1 := pmin(gene_a, gene_b)]
	out[, node2 := pmax(gene_a, gene_b)]
	if (!is.null(weight_col)) {
		out[, edge_weight := safe_numeric(get(weight_col))]
		out[is.na(edge_weight), edge_weight := 1]
	} else {
		out[, edge_weight := 1]
	}

	if (!is.null(evidence_col)) {
		collapsed <- out[, .(
			weight = max(edge_weight, na.rm = TRUE),
			evidence_count = uniqueN(get(evidence_col))
		), by = .(gene_a = node1, gene_b = node2)]
	} else {
		collapsed <- out[, .(
			weight = max(edge_weight, na.rm = TRUE),
			evidence_count = .N
		), by = .(gene_a = node1, gene_b = node2)]
	}
	collapsed[, source_db := db_name]
	setorder(collapsed, gene_a, gene_b)
	collapsed[]
}

read_deseq_results <- function(deseq_root, padj_threshold, abs_lfc_threshold) {
	files <- list.files(
		deseq_root,
		pattern = "DESeq2_results\\.csv$",
		recursive = TRUE,
		full.names = TRUE
	)
	if (!length(files)) {
		stop("No DESeq2_results.csv files found under: ", deseq_root)
	}

	root_norm <- normalizePath(deseq_root, winslash = "/", mustWork = TRUE)
	results <- vector("list", length(files))

	for (i in seq_along(files)) {
		file_path <- files[[i]]
		rel_path <- sub(
			paste0("^", root_norm, "/?"),
			"",
			normalizePath(file_path, winslash = "/", mustWork = TRUE)
		)
		parts <- strsplit(rel_path, "/", fixed = TRUE)[[1L]]
		if (length(parts) < 3L) {
			next
		}
		cell_type <- parts[[1L]]
		contrast <- parts[[2L]]
		contrast_parts <- strsplit(contrast, "_vs_", fixed = TRUE)[[1L]]
		if (length(contrast_parts) != 2L) {
			next
		}

		left_detail <- normalize_token(contrast_parts[[1L]])
		right_detail <- normalize_token(contrast_parts[[2L]])
		left_group <- site_group_from_token(left_detail)
		right_group <- site_group_from_token(right_detail)
		pair_label <- canonical_pair(left_group, right_group)

		dt <- fread(file_path, showProgress = FALSE)
		needed <- c("gene", "baseMean", "log2FoldChange", "padj")
		if (!all(needed %in% names(dt))) {
			warning("Skipping malformed DESeq2 file: ", file_path)
			next
		}

		dt[, gene_symbol := toupper(trimws(as.character(gene)))]
		dt[, baseMean := safe_numeric(baseMean)]
		dt[, log2FoldChange := safe_numeric(log2FoldChange)]
		dt[, padj_numeric := safe_numeric(padj)]
		dt <- dt[
			nzchar(gene_symbol) &
				!is.na(log2FoldChange) &
				!is.na(padj_numeric) &
				padj_numeric <= padj_threshold &
				abs(log2FoldChange) >= abs_lfc_threshold
		]
		if (!nrow(dt)) {
			next
		}

		dt[, abs_log2FoldChange := abs(log2FoldChange)]
		dt[, neg_log10_padj := -log10(padj_numeric)]
		dt[, direction_resolved := fifelse(log2FoldChange >= 0, "up", "down")]
		dt[, ppi_priority_score := abs_log2FoldChange * pmax(neg_log10_padj, 1)]
		dt[, `:=`(
			pair_label = pair_label,
			cell_type = cell_type,
			contrast = contrast,
			contrast_left = left_detail,
			contrast_right = right_detail,
			left_site_group = left_group,
			right_site_group = right_group,
			result_path = file_path
		)]

		results[[i]] <- dt[, .(
			pair_label,
			cell_type,
			contrast,
			contrast_left,
			contrast_right,
			left_site_group,
			right_site_group,
			gene = gene_symbol,
			gene_symbol,
			baseMean,
			log2FoldChange,
			abs_log2FoldChange,
			padj = padj_numeric,
			neg_log10_padj,
			direction_resolved,
			ppi_priority_score,
			result_path
		)]
	}

	out <- rbindlist(results, use.names = TRUE, fill = TRUE)
	if (!nrow(out)) {
		stop("All DESeq2 files were filtered out; no significant DEG rows survived.")
	}
	setorder(out, pair_label, cell_type, -ppi_priority_score, gene_symbol)
	out[]
}

collapse_seed_candidates <- function(deg_dt) {
	seed_dt <- deg_dt[, .(
		n_cell_types = uniqueN(cell_type),
		n_occurrences = .N,
		cell_types = paste(sort(unique(cell_type)), collapse = ";"),
		contrasts = paste(sort(unique(contrast)), collapse = ";"),
		mean_log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
		mean_abs_log2FoldChange = mean(abs_log2FoldChange, na.rm = TRUE),
		max_abs_log2FoldChange = max(abs_log2FoldChange, na.rm = TRUE),
		min_padj = min(padj, na.rm = TRUE),
		neg_log10_min_padj = max(neg_log10_padj, na.rm = TRUE),
		up_hits = sum(direction_resolved == "up", na.rm = TRUE),
		down_hits = sum(direction_resolved == "down", na.rm = TRUE),
		raw_priority_sum = sum(ppi_priority_score, na.rm = TRUE)
	), by = .(pair_label, gene_symbol)]

	seed_dt[, dominant_direction := fifelse(
		up_hits > down_hits,
		"up",
		fifelse(down_hits > up_hits, "down", "mixed")
	)]
	seed_dt[, direction_consistency := pmax(up_hits, down_hits) / pmax(n_occurrences, 1)]
	seed_dt[, string_input_gene := gene_symbol]
	seed_dt[, ppi_priority_score := raw_priority_sum * pmax(direction_consistency, 0.1)]
	seed_dt[, raw_priority_sum := NULL]

	setorder(seed_dt, pair_label, -ppi_priority_score, -n_cell_types, -max_abs_log2FoldChange, gene_symbol)
	seed_dt[, seed_rank := seq_len(.N), by = pair_label]
	seed_dt[]
}

probe_remote_file_size <- function(url, wget_bin) {
	if (!nzchar(wget_bin)) {
		return(NA_real_)
	}
	probe <- try(
		suppressWarnings(system2(
			wget_bin,
			args = c("--server-response", "--spider", "--header=Range: bytes=0-0", "-T", "60", url),
			stdout = TRUE,
			stderr = TRUE
		)),
		silent = TRUE
	)
	if (inherits(probe, "try-error")) {
		return(NA_real_)
	}
	content_range <- probe[grepl("content-range:", probe, ignore.case = TRUE, perl = TRUE)]
	if (length(content_range)) {
		total <- suppressWarnings(as.numeric(sub(".*?/([0-9]+)\\s*$", "\\1", content_range[[1L]], perl = TRUE)))
		if (!is.na(total) && total > 0) {
			return(total)
		}
	}
	length_line <- probe[grepl("^length:\\s*[0-9]+", trimws(probe), ignore.case = TRUE, perl = TRUE)]
	if (length(length_line)) {
		total <- suppressWarnings(as.numeric(sub("^.*?([0-9]+).*$", "\\1", trimws(length_line[[1L]]), perl = TRUE)))
		if (!is.na(total) && total > 0) {
			return(total)
		}
	}
	NA_real_
}

combine_binary_parts <- function(part_files, destfile) {
	out_con <- file(destfile, open = "wb")
	on.exit(close(out_con), add = TRUE)
	for (part_file in part_files) {
		in_con <- file(part_file, open = "rb")
		repeat {
			chunk <- readBin(in_con, what = "raw", n = 1024L * 1024L)
			if (!length(chunk)) {
				break
			}
			writeBin(chunk, out_con, useBytes = TRUE)
		}
		close(in_con)
	}
	invisible(destfile)
}

download_file_parallel_wget <- function(url, destfile, wget_bin, max_parts = 4L) {
	total_bytes <- probe_remote_file_size(url, wget_bin)
	if (is.na(total_bytes) || total_bytes <= 0) {
		return(list(ok = FALSE, error = "could not determine remote file size"))
	}
	if (total_bytes < (20 * 1024^2)) {
		return(list(ok = FALSE, error = "remote file below parallel-download threshold"))
	}

	part_count <- min(max_parts, max(1L, as.integer(floor(total_bytes / (5 * 1024^2)))))
	part_count <- max(part_count, 2L)
	part_count <- min(part_count, as.integer(total_bytes))
	if (part_count <= 1L) {
		return(list(ok = FALSE, error = "remote file too small for parallel download"))
	}

	part_dir <- paste0(destfile, ".parts")
	if (file.exists(destfile)) {
		unlink(destfile)
	}
	if (dir.exists(part_dir)) {
		unlink(part_dir, recursive = TRUE, force = TRUE)
	}
	dir_create(part_dir)

	base_size <- as.integer(total_bytes %/% part_count)
	remainder <- as.integer(total_bytes %% part_count)
	starts <- integer(part_count)
	ends <- integer(part_count)
	cursor <- 0L
	for (idx in seq_len(part_count)) {
		extra <- if (idx <= remainder) 1L else 0L
		part_size <- base_size + extra
		starts[[idx]] <- cursor
		ends[[idx]] <- cursor + part_size - 1L
		cursor <- cursor + part_size
	}

	part_specs <- data.table(
		part_id = seq_len(part_count),
		start_byte = starts,
		end_byte = ends,
		expected_size = ends - starts + 1L,
		part_file = file.path(part_dir, sprintf("part_%02d.bin", seq_len(part_count)))
	)

	results <- parallel::mclapply(seq_len(nrow(part_specs)), function(i) {
		row <- part_specs[i]
		args <- c(
			"-O", row$part_file,
			paste0("--header=Range:bytes=", row$start_byte, "-", row$end_byte),
			"--tries=3",
			"--timeout=300",
			"--waitretry=5",
			url
		)
		status <- suppressWarnings(system2(wget_bin, args = args, stdout = "", stderr = ""))
		size <- if (file.exists(row$part_file)) file.info(row$part_file)$size else NA_real_
		list(
			status = status,
			part_file = row$part_file,
			size = size,
			expected_size = row$expected_size
		)
	}, mc.cores = min(part_count, max_parts))

	status_ok <- vapply(results, function(x) identical(x$status, 0L), logical(1L))
	size_ok <- vapply(results, function(x) !is.na(x$size) && identical(as.numeric(x$size), as.numeric(x$expected_size)), logical(1L))
	if (!all(status_ok) || !all(size_ok)) {
		unlink(part_dir, recursive = TRUE, force = TRUE)
		if (file.exists(destfile)) {
			unlink(destfile)
		}
		return(list(
			ok = FALSE,
			error = sprintf(
				"parallel wget failed: status_ok=%s size_ok=%s",
				paste(status_ok, collapse = ","),
				paste(size_ok, collapse = ",")
			)
		))
	}

	combine_binary_parts(part_specs$part_file, destfile)
	unlink(part_dir, recursive = TRUE, force = TRUE)
	final_size <- if (file.exists(destfile)) file.info(destfile)$size else NA_real_
	if (is.na(final_size) || final_size != total_bytes) {
		if (file.exists(destfile)) {
			unlink(destfile)
		}
		return(list(ok = FALSE, error = sprintf("combined file size mismatch: got %s expected %s", final_size, total_bytes)))
	}
	list(ok = TRUE, bytes = final_size)
}

download_file_cached <- function(url, destfile, max_wget_tries = 3L, wget_timeout = 300L) {
	dir_create(dirname(destfile))
	ok_file <- paste0(destfile, ".ok")
	if (file.exists(destfile) && file.exists(ok_file) && isTRUE(file.info(destfile)$size > 0)) {
		return(list(path = destfile, downloaded = FALSE, bytes = file.info(destfile)$size, url = url))
	}

	last_error <- NULL
	old_timeout <- getOption("timeout")
	on.exit(options(timeout = old_timeout), add = TRUE)
	options(timeout = max(3600, old_timeout))
	wget_bin <- Sys.which("wget")
	if (nzchar(wget_bin)) {
		parallel_res <- download_file_parallel_wget(url, destfile, wget_bin)
		if (isTRUE(parallel_res$ok) && file.exists(destfile) && isTRUE(file.info(destfile)$size > 0)) {
			writeLines(c(
				paste0("url=", url),
				paste0("downloaded_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
				paste0("bytes=", file.info(destfile)$size),
				"backend=wget-parallel"
			), con = ok_file)
			return(list(path = destfile, downloaded = TRUE, bytes = file.info(destfile)$size, url = url))
		}
		last_error <- if (!is.null(parallel_res$error)) parallel_res$error else last_error

		wget_args <- c(
			"--continue",
			paste0("--tries=", as.integer(max_wget_tries)),
			paste0("--timeout=", as.integer(wget_timeout)),
			"--waitretry=5",
			"-O", destfile,
			url
		)
		wget_status <- suppressWarnings(system2(wget_bin, args = wget_args, stdout = "", stderr = ""))
		if (identical(wget_status, 0L) && file.exists(destfile) && isTRUE(file.info(destfile)$size > 0)) {
			writeLines(c(
				paste0("url=", url),
				paste0("downloaded_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
				paste0("bytes=", file.info(destfile)$size),
				"backend=wget"
			), con = ok_file)
			return(list(path = destfile, downloaded = TRUE, bytes = file.info(destfile)$size, url = url))
		}
		last_error <- paste0("wget exit status ", wget_status)
	}

	if (file.exists(destfile) && !file.exists(ok_file)) {
		unlink(destfile)
	}
	for (method in c("libcurl", "auto")) {
		res <- try(
			download.file(url, destfile = destfile, method = method, mode = "wb", quiet = FALSE),
			silent = TRUE
		)
		if (!inherits(res, "try-error") && file.exists(destfile) && isTRUE(file.info(destfile)$size > 0)) {
			writeLines(c(
				paste0("url=", url),
				paste0("downloaded_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
				paste0("bytes=", file.info(destfile)$size),
				paste0("backend=download.file:", method)
			), con = ok_file)
			return(list(path = destfile, downloaded = TRUE, bytes = file.info(destfile)$size, url = url))
		}
		last_error <- as.character(res)
		if (file.exists(destfile) && isTRUE(file.info(destfile)$size == 0)) {
			unlink(destfile)
		}
	}
	stop("Failed to download ", url, " -> ", destfile, " | ", last_error)
}

download_ppi_databases <- function(source_dir) {
	hint_spec <- list(
		database = "HINT_human_binary_hq",
		url = "https://hint.yulab.org/download-raw/2024-06/HomoSapiens_binary_hq.txt",
		local_path = file.path(source_dir, "hint", "2024_06", "HomoSapiens_binary_hq.txt")
	)
	hint_ok_file <- paste0(hint_spec$local_path, ".ok")
	partial_hint_size <- if (file.exists(hint_spec$local_path)) file.info(hint_spec$local_path)$size else NA_real_
	use_partial_hint <- file.exists(hint_spec$local_path) && !file.exists(hint_ok_file) && !is.na(partial_hint_size) && partial_hint_size >= (2 * 1024^2)
	if (use_partial_hint) {
		hint_download <- list(path = hint_spec$local_path, downloaded = FALSE, bytes = partial_hint_size, url = hint_spec$url)
		hint_ok <- TRUE
		hint_bytes <- partial_hint_size
		hint_mode <- "partial-local-after-timeout"
		hint_note <- "using existing partial HINT text file (>=2MB) after repeated remote timeouts"
	} else {
		hint_download <- try(
			download_file_cached(hint_spec$url, hint_spec$local_path, max_wget_tries = 1L, wget_timeout = 300L),
			silent = TRUE
		)
		hint_ok <- !inherits(hint_download, "try-error") && file.exists(hint_spec$local_path)
		hint_bytes <- if (hint_ok) file.info(hint_spec$local_path)$size else NA_real_
		hint_mode <- if (hint_ok) "downloaded-local" else "download-failed-continued"
		hint_note <- if (hint_ok) NA_character_ else as.character(hint_download)
	}
	manifest_dt <- rbindlist(list(
		data.table(
			database = "STRING_api",
			url = "https://string-db.org/api/tsv/network",
			local_path = NA_character_,
			bytes = NA_real_,
			downloaded = FALSE,
			access_mode = "api",
			status_note = NA_character_
		),
		data.table(
			database = "BioGRID_psicquic",
			url = "http://tyersrest.tyerslab.com:8805/psicquic/webservices/current/search/interactor/",
			local_path = NA_character_,
			bytes = NA_real_,
			downloaded = FALSE,
			access_mode = "psicquic-rest",
			status_note = NA_character_
		),
		data.table(
			database = hint_spec$database,
			url = hint_spec$url,
			local_path = hint_spec$local_path,
			bytes = hint_bytes,
			downloaded = hint_ok,
			access_mode = hint_mode,
			status_note = hint_note
		)
	), use.names = TRUE, fill = TRUE)

	list(
		manifest = manifest_dt,
		hint_file = if (hint_ok) hint_spec$local_path else NA_character_,
		string_api_url = "https://string-db.org/api/tsv/network",
		biogrid_psicquic_url = "http://tyersrest.tyerslab.com:8805/psicquic/webservices/current/search/interactor/"
	)
}

http_get_text <- function(url, timeout_seconds = 180) {
	resp <- httr2::request(url) |>
		httr2::req_timeout(timeout_seconds) |>
		httr2::req_user_agent("epithelial_airway_deg_ppi_workflow_20260531/1.0") |>
		httr2::req_perform()
	httr2::resp_body_string(resp)
}

build_string_master_edges <- function(seed_genes, min_score) {
	seed_genes <- unique(toupper(seed_genes))
	if (!length(seed_genes)) {
		return(empty_edge_table())
	}
	query_url <- paste0(
		"https://string-db.org/api/tsv/network?identifiers=",
		utils::URLencode(paste(seed_genes, collapse = "\n"), reserved = TRUE),
		"&species=9606",
		"&required_score=", as.integer(min_score),
		"&network_type=functional"
	)
	txt <- try(http_get_text(query_url, timeout_seconds = 180), silent = TRUE)
	if (inherits(txt, "try-error") || !nzchar(trimws(txt))) {
		return(empty_edge_table())
	}
	dt <- fread(text = txt, sep = "\t", quote = "", showProgress = FALSE)
	if (!nrow(dt) || !all(c("preferredName_A", "preferredName_B", "score") %in% names(dt))) {
		return(empty_edge_table())
	}
	dt[, gene_a := toupper(trimws(preferredName_A))]
	dt[, gene_b := toupper(trimws(preferredName_B))]
	dt <- dt[gene_a %in% seed_genes & gene_b %in% seed_genes]
	if (!nrow(dt)) {
		return(empty_edge_table())
	}
	dt[, string_score_api := pmax(safe_numeric(score), 0) * 1000]
	canonicalize_edges(dt, "gene_a", "gene_b", db_name = "STRING", weight_col = "string_score_api")
}

extract_psicquic_candidates <- function(field_value) {
	if (is.na(field_value) || !nzchar(trimws(as.character(field_value)))) {
		return(character())
	}
	pieces <- trimws(unlist(strsplit(as.character(field_value), "\\|", perl = TRUE), use.names = FALSE))
	values <- trimws(sub("^[^:]+:", "", pieces, perl = TRUE))
	values <- sub("\\(.*$", "", values, perl = TRUE)
	values <- toupper(gsub('"', "", values, fixed = TRUE))
	values <- values[
		nzchar(values) &
			grepl("[A-Z]", values, perl = TRUE) &
			!grepl("^[0-9]+$", values, perl = TRUE) &
			!grepl("^MI:", values, ignore.case = TRUE, perl = TRUE) &
			!grepl("^PUBMED", values, ignore.case = TRUE, perl = TRUE) &
			!grepl("^TAXID", values, ignore.case = TRUE, perl = TRUE)
	]
	unique(values)
}

pick_psicquic_symbol <- function(primary_field, alias_field, seed_genes) {
	candidates <- unique(c(
		extract_psicquic_candidates(primary_field),
		extract_psicquic_candidates(alias_field)
	))
	hit <- candidates[candidates %in% seed_genes]
	if (length(hit)) {
		return(hit[[1L]])
	}
	if (length(candidates)) {
		return(candidates[[1L]])
	}
	NA_character_
}

build_biogrid_master_edges <- function(seed_genes, chunk_size = 20L) {
	seed_genes <- unique(toupper(seed_genes))
	if (!length(seed_genes)) {
		return(empty_edge_table())
	}
	endpoint <- "http://tyersrest.tyerslab.com:8805/psicquic/webservices/current/search/interactor/"
	chunks <- split(seed_genes, ceiling(seq_along(seed_genes) / chunk_size))
	chunk_rows <- lapply(chunks, function(chunk_genes) {
		query <- paste(chunk_genes, collapse = " OR ")
		query_url <- paste0(endpoint, utils::URLencode(query, reserved = TRUE), "?format=tab25")
		txt <- try(http_get_text(query_url, timeout_seconds = 180), silent = TRUE)
		if (inherits(txt, "try-error") || !nzchar(trimws(txt))) {
			return(NULL)
		}
		dt <- fread(text = txt, sep = "\t", header = FALSE, fill = TRUE, quote = "", showProgress = FALSE)
		if (!nrow(dt)) {
			return(NULL)
		}
		for (col in paste0("V", 1:11)) {
			if (!col %in% names(dt)) {
				dt[, (col) := NA_character_]
			}
		}
		dt[, gene_a := mapply(pick_psicquic_symbol, V3, V5, MoreArgs = list(seed_genes = seed_genes), USE.NAMES = FALSE)]
		dt[, gene_b := mapply(pick_psicquic_symbol, V4, V6, MoreArgs = list(seed_genes = seed_genes), USE.NAMES = FALSE)]
		dt <- dt[
			grepl("taxid:9606", V10, fixed = TRUE) &
				grepl("taxid:9606", V11, fixed = TRUE) &
				gene_a %in% seed_genes &
				gene_b %in% seed_genes
		]
		if (!nrow(dt)) {
			return(NULL)
		}
		dt[, publication_key := V9]
		dt
	})
	out <- rbindlist(chunk_rows, use.names = TRUE, fill = TRUE)
	if (!nrow(out)) {
		return(empty_edge_table())
	}
	canonicalize_edges(out, "gene_a", "gene_b", db_name = "BioGRID", evidence_col = "publication_key")
}

build_hint_master_edges <- function(hint_file, seed_genes) {
	if (is.na(hint_file) || !file.exists(hint_file)) {
		return(empty_edge_table())
	}
	seed_genes <- unique(toupper(seed_genes))
	dt <- fread(hint_file, sep = "\t", fill = TRUE, showProgress = FALSE)
	cols <- names(dt)
	col_a <- pick_col(cols, c("^Gene_A$", "^Gene A$"))
	col_b <- pick_col(cols, c("^Gene_B$", "^Gene B$"))
	taxid_col <- pick_col(cols, c("^taxid$"), required = FALSE)
	quality_col <- pick_col(cols, c("^high_quality$"), required = FALSE)

	if (!is.na(taxid_col)) {
		dt <- dt[safe_numeric(get(taxid_col)) == 9606]
	}
	if (!is.na(quality_col)) {
		dt <- dt[tolower(as.character(get(quality_col))) %in% c("true", "t", "1")]
	}

	dt[, gene_a := toupper(trimws(as.character(get(col_a))))]
	dt[, gene_b := toupper(trimws(as.character(get(col_b))))]
	dt <- dt[gene_a %in% seed_genes & gene_b %in% seed_genes]
	if (!nrow(dt)) {
		return(empty_edge_table())
	}

	canonicalize_edges(dt, "gene_a", "gene_b", db_name = "HINT")
}

write_seed_inputs <- function(seed_dt, top_n, output_dir) {
	dir_create(output_dir)
	input_dir <- file.path(output_dir, "ppi_input", "pairs")
	dir_create(input_dir)

	all_seed_path <- file.path(output_dir, "seed_summary_by_pair.tsv")
	top_seed_path <- file.path(output_dir, sprintf("seed_summary_by_pair_top%d.tsv", top_n))
	fwrite(seed_dt, all_seed_path, sep = "\t")
	top_seed_dt <- seed_dt[seed_rank <= top_n]
	fwrite(top_seed_dt, top_seed_path, sep = "\t")

	manifest_rows <- list()
	pair_levels <- unique(seed_dt$pair_label)
	for (pair_label in pair_levels) {
		pair_value <- pair_label
		safe_pair <- gsub("[^A-Za-z0-9_]+", "_", pair_label, perl = TRUE)
		pair_all <- seed_dt[pair_label == pair_value, gene_symbol]
		pair_top <- top_seed_dt[pair_label == pair_value, gene_symbol]
		all_path <- file.path(input_dir, sprintf("%s__all_sig.txt", safe_pair))
		top_path <- file.path(input_dir, sprintf("%s__top%d.txt", safe_pair, top_n))
		writeLines(pair_all, con = all_path)
		writeLines(pair_top, con = top_path)
		manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
			pair_label = pair_label,
			file_kind = c("all_sig", sprintf("top%d", top_n)),
			n_genes = c(length(pair_all), length(pair_top)),
			file_path = c(all_path, top_path)
		)
	}
	manifest_dt <- rbindlist(manifest_rows, use.names = TRUE)
	fwrite(manifest_dt, file.path(output_dir, "seed_input_manifest.tsv"), sep = "\t")
	list(top_seed_dt = top_seed_dt, manifest_dt = manifest_dt)
}

plot_graph_bundle <- function(graph, nodes_dt, pair_label, db_name, out_dir) {
	dir_create(out_dir)
	rescale_numeric <- function(x, to = c(0, 1), fallback = mean(to)) {
		x <- as.numeric(x)
		if (!length(x)) {
			return(numeric())
		}
		keep <- is.finite(x)
		if (!any(keep)) {
			return(rep(fallback, length(x)))
		}
		x_min <- min(x[keep])
		x_max <- max(x[keep])
		if (isTRUE(all.equal(x_min, x_max))) {
			return(rep(mean(to), length(x)))
		}
		to[[1L]] + ((x - x_min) / (x_max - x_min)) * (to[[2L]] - to[[1L]])
	}
	prepare_plot_nodes <- function(graph_obj, all_nodes_dt) {
		plot_dt <- merge(
			data.table(name = V(graph_obj)$name, vertex_order = seq_along(V(graph_obj)$name)),
			all_nodes_dt,
			by = "name",
			all.x = TRUE
		)
		setorder(plot_dt, vertex_order)
		plot_dt[, vertex_order := NULL]
		plot_dt[]
	}
	draw_summary_panel <- function(direction_colors, summary_lines, top_hub_text) {
		plot.new()
		title(main = "network summary", cex.main = 1.1, font.main = 2)
		legend(
			"topleft",
			inset = 0.02,
			legend = c("up", "down", "mixed"),
			pch = 16,
			col = unname(direction_colors[c("up", "down", "mixed")]),
			bty = "n",
			cex = 0.95
		)
		text(
			x = 0.02,
			y = 0.78,
			labels = paste(summary_lines, collapse = "\n"),
			adj = c(0, 1),
			cex = 1.02
		)
		text(
			x = 0.02,
			y = 0.36,
			labels = paste0("Top hubs\n", top_hub_text),
			adj = c(0, 1),
			cex = 0.98
		)
	}

	direction_colors <- c(up = "#C44E52", down = "#4C72B0", mixed = "#8C8C8C")
	plot_paths <- c(
		file.path(out_dir, "network_plot.png"),
		file.path(out_dir, "network_plot.pdf")
	)
	total_nodes <- vcount(graph)
	total_edges <- ecount(graph)
	connected_seed_count <- sum(nodes_dt$degree > 0, na.rm = TRUE)
	isolated_count <- sum(nodes_dt$degree == 0, na.rm = TRUE)
	connected_seed_names <- nodes_dt[degree > 0, name]

	if (length(connected_seed_names)) {
		core_graph <- induced_subgraph(graph, vids = connected_seed_names)
	} else {
		core_graph <- make_empty_graph(n = 0L, directed = FALSE)
	}
	if (vcount(core_graph) > 60L) {
		core_ranked <- nodes_dt[name %in% V(core_graph)$name][order(-degree, seed_rank, -ppi_priority_score)]
		keep_core_names <- core_ranked$name[seq_len(min(60L, nrow(core_ranked)))]
		core_graph <- induced_subgraph(core_graph, vids = keep_core_names)
	}
	omitted_connected_count <- max(0L, connected_seed_count - vcount(core_graph))
	core_nodes <- if (vcount(core_graph)) prepare_plot_nodes(core_graph, nodes_dt) else nodes_dt[0]
	core_components <- if (vcount(core_graph)) components(core_graph) else NULL
	core_component_count <- if (is.null(core_components)) 0L else core_components$no
	core_largest_component <- if (is.null(core_components)) 0L else max(core_components$csize)
	top_hubs <- nodes_dt[degree > 0][order(-degree, seed_rank, -ppi_priority_score)][seq_len(min(5L, .N))]
	top_hub_text <- if (nrow(top_hubs)) {
		paste(sprintf("%s (%d)", top_hubs$gene_symbol, as.integer(top_hubs$degree)), collapse = "\n")
	} else {
		"none"
	}
	summary_lines <- c(
		sprintf("seed genes: %d", total_nodes),
		sprintf("edges: %d", total_edges),
		sprintf("connected seeds total: %d", connected_seed_count),
		sprintf("connected seeds shown: %d", vcount(core_graph)),
		sprintf("connected omitted by cap: %d", omitted_connected_count),
		sprintf("isolated seeds hidden: %d", isolated_count),
		sprintf("connected components: %d", core_component_count),
		sprintf("largest component: %d nodes", core_largest_component)
	)
	footer_line <- if (total_edges > 0L) {
		if (omitted_connected_count > 0L) {
			sprintf(
				"主图显示 %d/%d 个有边 seed；另有 %d 个有边 seed 因展示上限省略，%d 个孤立 seed 已移到统计说明。",
				vcount(core_graph),
				connected_seed_count,
				omitted_connected_count,
				isolated_count
			)
		} else {
			sprintf("主图仅显示有边的核心网络；其余 %d 个孤立 seed 已移到统计说明。", isolated_count)
		}
	} else {
		"当前阈值下未观察到 seed 内部 PPI 边；保留统计说明便于回查。"
	}

	for (plot_path in plot_paths) {
		if (grepl("\\.png$", plot_path, ignore.case = TRUE)) {
			png(plot_path, width = 2200, height = 1500, res = 220)
		} else {
			pdf(plot_path, width = 14, height = 9)
		}
		old_par <- par(no.readonly = TRUE)
		layout(matrix(c(1, 2), nrow = 1L), widths = c(4.1, 1.6))

		par(mar = c(1.8, 1.8, 4.2, 0.8))
		if (total_edges == 0L || vcount(core_graph) == 0L) {
			plot.new()
			title(main = sprintf("%s | %s", pair_label, db_name), cex.main = 1.7, font.main = 2)
			text(0.5, 0.62, labels = "No connected seed network", cex = 1.6, font = 2)
			text(
				0.5,
				0.42,
				labels = sprintf(
					"%d 个 seed 通过筛选，但当前 %s 条件下没有保留下来的内部互作边。",
					total_nodes,
					db_name
				),
				cex = 1.05
			)
			text(0.5, 0.25, labels = footer_line, cex = 0.95, col = "grey30")
		} else {
			core_priority <- log1p(pmax(core_nodes$ppi_priority_score, 0))
			core_degree <- pmax(core_nodes$degree, 0)
			vertex_size <- rescale_numeric(core_priority, to = c(9, 21), fallback = 12) +
				rescale_numeric(core_degree, to = c(0, 6), fallback = 2)
			vertex_colors <- direction_colors[ifelse(core_nodes$dominant_direction %in% names(direction_colors), core_nodes$dominant_direction, "mixed")]
			vertex_frame <- ifelse(core_degree >= 2, "#1F1F1F", "white")
			label_limit <- if (vcount(core_graph) <= 18L) {
				vcount(core_graph)
			} else if (vcount(core_graph) <= 35L) {
				18L
			} else {
				14L
			}
			labels <- ifelse(
				rank(-core_nodes$degree, ties.method = "first") <= label_limit | core_nodes$degree >= 3,
				core_nodes$gene_symbol,
				NA_character_
			)

			if (vcount(core_graph) == 1L) {
				coords <- matrix(c(0, 0), ncol = 2L)
			} else if (core_component_count > 1L || vcount(core_graph) <= 25L) {
				coords <- layout_with_kk(core_graph)
			} else {
				coords <- layout_with_fr(core_graph, niter = 2000)
			}
			coords <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

			edge_weight_log <- log1p(pmax(as.numeric(E(core_graph)$weight), 0))
			edge_width <- rescale_numeric(
				edge_weight_log,
				to = if (identical(db_name, "consensus")) c(4, 9.5) else c(2.2, 7.2),
				fallback = if (identical(db_name, "consensus")) 5.8 else 3.5
			)
			edge_alpha <- rescale_numeric(
				edge_weight_log,
				to = if (identical(db_name, "consensus")) c(0.68, 0.95) else c(0.38, 0.82),
				fallback = if (identical(db_name, "consensus")) 0.82 else 0.55
			)
			edge_base_color <- if (identical(db_name, "consensus")) "#667A95" else "#7A7A7A"
			edge_colors <- vapply(edge_alpha, function(alpha) {
				grDevices::adjustcolor(edge_base_color, alpha.f = alpha)
			}, FUN.VALUE = character(1L))

			plot(
				core_graph,
				layout = coords,
				vertex.label = labels,
				vertex.label.cex = if (vcount(core_graph) <= 20L) 0.95 else 0.8,
				vertex.label.family = "sans",
				vertex.label.color = "black",
				vertex.size = vertex_size,
				vertex.color = vertex_colors,
				vertex.frame.color = vertex_frame,
				vertex.frame.width = 1.2,
				edge.width = edge_width,
				edge.color = edge_colors,
				main = sprintf("%s | %s", pair_label, db_name),
				margin = 0.08
			)
			mtext(footer_line, side = 1, line = -1.2, cex = 0.82, col = "grey30")
		}

		par(mar = c(2.5, 0.5, 4.2, 0.8))
		draw_summary_panel(direction_colors, summary_lines, top_hub_text)
		par(old_par)
		dev.off()
	}
}

save_graph_bundle <- function(pair_label, db_name, edge_dt, pair_seed_dt, out_dir) {
	dir_create(out_dir)
	vertices <- unique(pair_seed_dt[, .(
		name = gene_symbol,
		gene_symbol,
		seed_rank,
		n_cell_types,
		n_occurrences,
		dominant_direction,
		direction_consistency,
		ppi_priority_score,
		mean_abs_log2FoldChange,
		max_abs_log2FoldChange,
		min_padj,
		neg_log10_min_padj
	)])
	setorder(vertices, seed_rank)

	if (!nrow(edge_dt)) {
		graph <- make_empty_graph(n = nrow(vertices), directed = FALSE)
		if (nrow(vertices)) {
			V(graph)$name <- vertices$name
		}
	} else {
		graph <- graph_from_data_frame(
			d = as.data.frame(edge_dt),
			directed = FALSE,
			vertices = as.data.frame(vertices)
		)
	}

	if (vcount(graph)) {
		deg <- degree(graph, mode = "all")
		btw <- betweenness(graph, directed = FALSE, normalized = FALSE)
		comp <- components(graph)
		metrics_dt <- data.table(
			name = names(deg),
			degree = as.numeric(deg),
			betweenness = as.numeric(btw),
			component_id = as.integer(comp$membership[names(deg)]),
			component_size = as.integer(comp$csize[comp$membership[names(deg)]])
		)
	} else {
		metrics_dt <- data.table(
			name = character(),
			degree = numeric(),
			betweenness = numeric(),
			component_id = integer(),
			component_size = integer()
		)
	}

	nodes_dt <- merge(vertices, metrics_dt, by = "name", all.x = TRUE)
	nodes_dt[is.na(degree), degree := 0]
	nodes_dt[is.na(betweenness), betweenness := 0]
	nodes_dt[is.na(component_id), component_id := 0L]
	nodes_dt[is.na(component_size), component_size := 1L]
	setorder(nodes_dt, -degree, seed_rank, -ppi_priority_score, gene_symbol)

	hub_dt <- copy(nodes_dt[degree > 0])[seq_len(min(50L, sum(nodes_dt$degree > 0)))]
	summary_dt <- data.table(
		pair_label = pair_label,
		database = db_name,
		n_seed_genes = nrow(vertices),
		n_nodes = vcount(graph),
		n_connected_nodes = sum(nodes_dt$degree > 0),
		n_edges = ecount(graph),
		n_isolated_seeds = sum(nodes_dt$degree == 0),
		density = if (vcount(graph) > 1L) edge_density(graph, loops = FALSE) else NA_real_,
		n_connected_components = if (vcount(graph)) components(graph)$no else 0L,
		largest_component_nodes = if (vcount(graph)) max(components(graph)$csize) else 0L,
		top_hub_gene = if (nrow(hub_dt)) hub_dt$gene_symbol[[1L]] else NA_character_,
		top_hub_degree = if (nrow(hub_dt)) hub_dt$degree[[1L]] else NA_real_
	)

	fwrite(edge_dt, file.path(out_dir, "edges.tsv"), sep = "\t")
	fwrite(nodes_dt, file.path(out_dir, "nodes.tsv"), sep = "\t")
	fwrite(hub_dt, file.path(out_dir, "hub_genes.tsv"), sep = "\t")
	fwrite(summary_dt, file.path(out_dir, "network_stats.tsv"), sep = "\t")
	plot_graph_bundle(graph, nodes_dt, pair_label, db_name, out_dir)
	summary_dt
}

build_consensus_edges <- function(edge_tables, min_support) {
	merged <- rbindlist(edge_tables, use.names = TRUE, fill = TRUE)
	if (!nrow(merged)) {
		return(empty_edge_table())
	}
	consensus <- merged[, .(
		weight = max(weight, na.rm = TRUE),
		evidence_count = sum(evidence_count, na.rm = TRUE),
		db_support = uniqueN(source_db),
		source_db = paste(sort(unique(source_db)), collapse = "|")
	), by = .(gene_a, gene_b)]
	consensus <- consensus[db_support >= min_support]
	setorder(consensus, -db_support, -weight, gene_a, gene_b)
	consensus[]
}

write_summary_report <- function(output_dir, config, db_manifest_dt, top_seed_dt, network_summary_dt) {
	pair_summary_dt <- top_seed_dt[, .(
		n_top_seed_genes = .N,
		top_seed = gene_symbol[[1L]],
		top_seed_score = ppi_priority_score[[1L]]
	), by = pair_label]
	pair_lines <- unlist(lapply(seq_len(nrow(pair_summary_dt)), function(i) {
		row <- pair_summary_dt[i]
		sprintf(
			"- `%s`: top%d seed genes = %d; top seed = %s (score %.2f)",
			row$pair_label,
			config$top_n,
			row$n_top_seed_genes,
			row$top_seed,
			row$top_seed_score
		)
	}))

	network_lines <- unlist(lapply(seq_len(nrow(network_summary_dt)), function(i) {
		row <- network_summary_dt[i]
		sprintf(
			"- `%s` / `%s`: nodes=%d, edges=%d, isolated=%d, top hub=%s (degree %.0f)",
			row$pair_label,
			row$database,
			row$n_nodes,
			row$n_edges,
			row$n_isolated_seeds,
			ifelse(is.na(row$top_hub_gene), "NA", row$top_hub_gene),
			ifelse(is.na(row$top_hub_degree), 0, row$top_hub_degree)
		)
	}))

	report_lines <- c(
		"# Epithelial airway DEG → PPI workflow (R-only)",
		"",
		sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
		sprintf("- DE root: `%s`", config$deseq_root),
		sprintf("- Source dir: `%s`", config$source_dir),
		sprintf("- Output dir: `%s`", config$output_dir),
		sprintf("- padj threshold: %.4f", config$padj_threshold),
		sprintf("- |log2FC| threshold: %.2f", config$abs_lfc_threshold),
		sprintf("- STRING combined score threshold: %d", config$string_score_threshold),
		sprintf("- topN per pair: %d", config$top_n),
		"",
		"## Guardrail",
		"",
		"- DEG 只是 PPI 数据库映射的候选输入，不是从 RNA 表达直接推断蛋白互作。",
		"- 这里输出的是数据库支持下的 seed-network / hub ranking，后续仍应回到表达层、文献和功能实验做验证。",
		"",
		"## Database access summary",
		""
	)
	report_lines <- c(
		report_lines,
		vapply(seq_len(nrow(db_manifest_dt)), function(i) {
			row <- db_manifest_dt[i]
			if (!is.na(row$local_path) && nzchar(row$local_path)) {
				line <- sprintf(
					"- `%s`: `%s` [%s] (%s bytes)",
					row$database,
					row$local_path,
					row$access_mode,
					format(as.numeric(row$bytes), scientific = FALSE)
				)
			} else {
				line <- sprintf(
					"- `%s`: `%s` [%s]",
					row$database,
					row$url,
					row$access_mode
				)
			}
			if (!is.na(row$status_note) && nzchar(row$status_note)) {
				line <- paste0(line, " | note: ", row$status_note)
			}
			line
		}, FUN.VALUE = character(1L)),
		"",
		"## Pair-level seed overview",
		"",
		pair_lines,
		"",
		"## Network summary",
		"",
		network_lines
	)
	writeLines(report_lines, con = file.path(output_dir, "PPI_WORKFLOW_SUMMARY.md"))
}

main <- function() {
	config <- validate_config(parse_args(DEFAULTS))
	dir_create(config$output_dir)

	cat("[1/7] Reading DESeq2 outputs...\n")
	deg_dt <- read_deseq_results(
		deseq_root = config$deseq_root,
		padj_threshold = config$padj_threshold,
		abs_lfc_threshold = config$abs_lfc_threshold
	)

	if (!is.null(config$target_pairs) && length(config$target_pairs)) {
		keep_pairs <- unique(normalize_token(config$target_pairs))
		deg_dt <- deg_dt[pair_label %in% keep_pairs]
		if (!nrow(deg_dt)) {
			stop("No DEG rows remain after filtering target pairs: ", paste(keep_pairs, collapse = ", "))
		}
	}
	fwrite(deg_dt, file.path(config$output_dir, "deg_sig_long_table.tsv"), sep = "\t")

	cat("[2/7] Collapsing seed candidates by pair...\n")
	seed_dt <- collapse_seed_candidates(deg_dt)
	fwrite(seed_dt, file.path(config$output_dir, "seed_summary_by_pair_all.tsv"), sep = "\t")

	cat("[3/7] Preparing STRING / BioGRID / HINT access...\n")
	db_paths <- download_ppi_databases(config$source_dir)
	fwrite(db_paths$manifest, file.path(config$output_dir, "db_manifest.tsv"), sep = "\t")

	cat("[4/7] Writing pair seed inputs...\n")
	seed_outputs <- write_seed_inputs(seed_dt, config$top_n, config$output_dir)
	top_seed_dt <- seed_outputs$top_seed_dt
	pair_labels <- unique(top_seed_dt$pair_label)
	union_seed_genes <- unique(top_seed_dt$gene_symbol)

	cat("[5/7] Building master interaction tables...\n")
	string_edges_all <- build_string_master_edges(
		union_seed_genes,
		min_score = config$string_score_threshold
	)
	biogrid_edges_all <- build_biogrid_master_edges(
		union_seed_genes
	)
	hint_edges_all <- build_hint_master_edges(
		db_paths$hint_file,
		union_seed_genes
	)

	master_manifest_dt <- data.table(
		database = c("STRING", "BioGRID", "HINT"),
		n_edges = c(nrow(string_edges_all), nrow(biogrid_edges_all), nrow(hint_edges_all)),
		n_unique_genes = c(
			uniqueN(c(string_edges_all$gene_a, string_edges_all$gene_b)),
			uniqueN(c(biogrid_edges_all$gene_a, biogrid_edges_all$gene_b)),
			uniqueN(c(hint_edges_all$gene_a, hint_edges_all$gene_b))
		)
	)
	fwrite(master_manifest_dt, file.path(config$output_dir, "master_edge_manifest.tsv"), sep = "\t")

	cat("[6/7] Building per-pair networks...\n")
	network_summary_rows <- list()
	edge_tables <- list(STRING = string_edges_all, BioGRID = biogrid_edges_all, HINT = hint_edges_all)

	for (pair_label in pair_labels) {
		cat("  - Pair:", pair_label, "\n")
		pair_value <- pair_label
		pair_seed_dt <- top_seed_dt[pair_label == pair_value]
		pair_genes <- unique(pair_seed_dt$gene_symbol)
		pair_dir <- file.path(config$output_dir, "networks", pair_label)
		dir_create(pair_dir)

		pair_edge_tables <- list()
		for (db_name in names(edge_tables)) {
			db_edges <- edge_tables[[db_name]][gene_a %in% pair_genes & gene_b %in% pair_genes]
			pair_edge_tables[[db_name]] <- db_edges
			db_dir <- file.path(pair_dir, db_name)
			summary_dt <- save_graph_bundle(pair_label, db_name, db_edges, pair_seed_dt, db_dir)
			network_summary_rows[[length(network_summary_rows) + 1L]] <- summary_dt
		}

		consensus_edges <- build_consensus_edges(pair_edge_tables, min_support = config$min_consensus_support)
		consensus_dir <- file.path(pair_dir, "consensus")
		summary_dt <- save_graph_bundle(pair_label, "consensus", consensus_edges, pair_seed_dt, consensus_dir)
		network_summary_rows[[length(network_summary_rows) + 1L]] <- summary_dt
	}

	network_summary_dt <- rbindlist(network_summary_rows, use.names = TRUE, fill = TRUE)
	fwrite(network_summary_dt, file.path(config$output_dir, "network_summary.tsv"), sep = "\t")

	cat("[7/7] Writing markdown summary...\n")
	write_summary_report(config$output_dir, config, db_paths$manifest, top_seed_dt, network_summary_dt)

	cat("Done. Key outputs:\n")
	cat("  -", file.path(config$output_dir, "db_manifest.tsv"), "\n")
	cat("  -", file.path(config$output_dir, sprintf("seed_summary_by_pair_top%d.tsv", config$top_n)), "\n")
	cat("  -", file.path(config$output_dir, "network_summary.tsv"), "\n")
	cat("  -", file.path(config$output_dir, "PPI_WORKFLOW_SUMMARY.md"), "\n")
}

main()
