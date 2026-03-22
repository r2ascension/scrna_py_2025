# interpret_agent_hotfix.R
# Usage:
#   source("interpret.txt")            # or put interpret.txt path below
#   source("interpret_agent_hotfix.R") # override functions in current session

# 1) 先加载原始源码（改成你的 interpret.txt 路径）
#    如果你主脚本里已经 source 过 interpret.txt，这行可以注释掉


# ---- helper: ensure list or NULL ----
.ensure_list <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) return(x)
  NULL
}

# 2) 覆盖 interpret_agent（最小侵入：只加 is.list 防护 + 合并时也防护）
interpret_agent <- function(x, context = NULL, n_pathways = 50,
                            model = "deepseek-chat", api_key = NULL,
                            add_ppi = FALSE, gene_fold_change = NULL) {

  if (missing(x)) stop("enrichment result 'x' is required.")

  res_list <- process_enrichment_input(x, n_pathways)
  if (length(res_list) == 0) return("No significant pathways found to interpret.")

  results <- lapply(names(res_list), function(name) {

    item <- res_list[[name]]
    df <- item$df
    original_genes <- item$genes

    fallback_mode <- FALSE
    pathway_text <- ""

    if (nrow(df) == 0) {
      if (!is.null(original_genes) && length(original_genes) > 0) {
        fallback_mode <- TRUE
        warning(sprintf("Cluster '%s': No enriched pathways. Falling back to gene-based interpretation. Confidence may be lower.", name))
        pathway_text <- paste("No significant pathways enriched.",
                              "Top Genes:", paste(head(original_genes, 50), collapse = ", "))
      } else {
        return(NULL)
      }
    } else {
      message(sprintf("Processing cluster '%s' with Agent 1: The Cleaner...", name))
      cols_to_keep <- intersect(c("ID","Description","GeneRatio","NES","p.adjust","pvalue","geneID"), names(df))
      pathway_text <- paste(
        apply(df[, cols_to_keep, drop = FALSE], 1, function(row) {
          paste(names(row), row, sep=": ", collapse=", ")
        }),
        collapse = "\n"
      )
    }

    # --- Step 1: Agent Cleaner ---
    cleaned_pathways <- pathway_text
    if (!fallback_mode) {
      clean_res <- .ensure_list(run_agent_cleaner(pathway_text, context, model, api_key))

      # ⭐ 关键修复：必须是 list 且包含 kept_pathways
      if (is.null(clean_res) || is.null(clean_res$kept_pathways)) {
        warning("Agent Cleaner failed or returned non-JSON/empty results. Falling back to using top pathways.")
        cleaned_pathways <- pathway_text
      } else {
        cleaned_pathways <- paste(
          "Selected Relevant Pathways (filtered by Agent Cleaner):",
          paste(clean_res$kept_pathways, collapse = ", "),
          "\nReasoning:", clean_res$reasoning,
          sep = "\n"
        )
      }
    }

    # --- Step 2: Agent Detective ---
    message(sprintf("Processing cluster '%s' with Agent 2: The Detective...", name))

    ppi_network_text <- NULL
    if (add_ppi) {
      all_genes <- if (fallback_mode) original_genes else unique(unlist(strsplit(df$geneID, "/")))
      if (length(all_genes) > 0) ppi_network_text <- .get_ppi_context_text(all_genes, x)
    }

    fc_text <- NULL
    if (!is.null(gene_fold_change)) {
      all_genes <- if (fallback_mode) original_genes else unique(unlist(strsplit(df$geneID, "/")))
      common_genes <- intersect(all_genes, names(gene_fold_change))
      if (length(common_genes) > 0) {
        fc_subset <- gene_fold_change[common_genes]
        fc_subset <- fc_subset[order(abs(fc_subset), decreasing = TRUE)]
        top_fc <- head(fc_subset, 20)
        fc_text <- paste(names(top_fc), round(top_fc, 2), sep=":", collapse=", ")
      }
    }

    detective_res <- .ensure_list(run_agent_detective(
      cleaned_pathways, ppi_network_text, fc_text, context, model, api_key, fallback_mode
    ))

    # --- Step 3: Agent Synthesizer ---
    message(sprintf("Processing cluster '%s' with Agent 3: The Storyteller...", name))

    final_res <- run_agent_synthesizer(cleaned_pathways, detective_res, context, model, api_key, fallback_mode)

    # ⭐ 关键修复：Synthesizer 若返回 character，则包装成 list
    if (!is.list(final_res)) {
      final_res <- list(
        cluster = name,
        overview = as.character(final_res),
        confidence = "Low",
        reasoning = "Failed to parse structured response from LLM (non-JSON)."
      )
      class(final_res) <- c("interpretation", "list")
      return(final_res)
    }

    # post-processing
    final_res$cluster <- name
    if (fallback_mode) final_res$data_source <- "gene_list_only"

    # ⭐ 关键修复：只有 detective_res 是 list 才合并字段
    if (!is.null(detective_res) && is.list(detective_res)) {
      final_res$regulatory_drivers <- detective_res$key_drivers
      final_res$refined_network <- detective_res$refined_network
      final_res$network_evidence <- detective_res$network_evidence
    }

    # refined_network -> igraph（沿用原逻辑）
    if (!is.null(final_res$refined_network)) {
      rn_df <- tryCatch({
        if (is.data.frame(final_res$refined_network)) {
          final_res$refined_network
        } else {
          do.call(rbind, lapply(final_res$refined_network, as.data.frame))
        }
      }, error = function(e) NULL)

      if (!is.null(rn_df) && nrow(rn_df) > 0) {
        colnames(rn_df)[colnames(rn_df) == "source"] <- "from"
        colnames(rn_df)[colnames(rn_df) == "target"] <- "to"
        if ("from" %in% names(rn_df) && "to" %in% names(rn_df)) {
          final_res$network <- igraph::graph_from_data_frame(rn_df, directed = FALSE)
        }
      }
    }

    final_res
  })

  names(results) <- names(res_list)
  if (length(results) == 1 && names(results)[1] == "Default") {
    results[[1]]
  } else {
    class(results) <- c("interpretation_list", "list")
    results
  }
}
