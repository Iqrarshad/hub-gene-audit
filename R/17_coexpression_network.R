# Network built from the analysed cohorts rather than from STRING.
#
# Edges are computed on residuals after removing neuronal, glial and immune
# content. Raw co-expression would rebuild the composition artefact:
# neuronal genes correlate with each other because they track neuron
# content, not because they interact.
#
# Both networks are thresholded to the same edge count so density cannot
# confound the comparison.
#
# Co-expression has its own biases toward highly expressed and high
# variance genes; these are tested alongside so that trading one bias for
# another is visible.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC); library(Matrix)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")

mcc_score <- function(g) {
  cl <- max_cliques(g, min = 1)
  s <- setNames(numeric(vcount(g)), V(g)$name)
  for (c in cl) { k <- length(c)
    if (k > 1) for (v in names(c)) s[v] <- s[v] + factorial(k - 1) }
  s[s == 0] <- 1; s
}

hub_select <- function(g, top_k = 10, min_methods = 4) {
  if (ecount(g) == 0) return(list(hubs = character(0), tab = NULL))
  vs <- V(g)$name
  tab <- tibble(gene = vs, degree = degree(g),
                betweenness = betweenness(g, normalized = TRUE),
                closeness = closeness(g, normalized = TRUE),
                eigenvector = eigen_centrality(g)$vector,
                pagerank = page_rank(g)$vector,
                mcc = mcc_score(g)[vs])
  mets <- setdiff(names(tab), "gene")
  top <- vapply(mets, function(m)
    as.integer(tab$gene %in% tab$gene[head(order(tab[[m]],
                                                 decreasing = TRUE), top_k)]),
    integer(nrow(tab)))
  tab$n_methods <- rowSums(top)
  list(hubs = tab$gene[tab$n_methods >= min_methods], tab = tab)
}

predictability_auc <- function(nominated, all_genes, predictor) {
  d <- tibble(gene = all_genes,
              y = as.integer(all_genes %in% nominated),
              x = predictor[match(all_genes, names(predictor))]) %>%
    filter(!is.na(x))
  if (sum(d$y) < 3 || sum(d$y) == nrow(d)) return(NA_real_)
  fit <- glm(y ~ x, data = d, family = binomial())
  as.numeric(pROC::auc(pROC::roc(d$y, predict(fit, type = "response"),
                                 quiet = TRUE)))
}

# --- Residualised co-expression -----------------------------------------
build_coexpression <- function(mat, genes, n_edges, label) {
  g <- intersect(genes, rownames(mat))
  log_msg(label, ": ", length(g), " of ", length(genes), " genes present")
  if (length(g) < 30) stop("Too few genes present in ", label)

  lm_ <- log2(mat + 1)

  # Composition covariates, from the same panels used for differential
  # expression, so the two halves of the study are consistent
  neuro  <- composition_score(mat, NEURONAL_PANEL, exclude = g,
                              label = paste(label, "neuronal"))
  glial  <- composition_score(mat, GLIAL_PANEL,  exclude = g,
                              label = paste(label, "glial"))
  immune <- composition_score(mat, IMMUNE_PANEL, exclude = g,
                              label = paste(label, "immune"))
  Z <- cbind(1, neuro, glial, immune)
  keep_z <- apply(Z, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  Z <- Z[, keep_z | seq_len(ncol(Z)) == 1, drop = FALSE]

  X <- t(lm_[g, , drop = FALSE])          # samples x genes
  ok <- complete.cases(X) & complete.cases(Z)
  X <- X[ok, , drop = FALSE]; Zs <- Z[ok, , drop = FALSE]

  # Residualise every gene on composition
  qrz <- qr(Zs)
  R <- X - Zs %*% qr.coef(qrz, X)
  log_msg(label, ": residualised on ", ncol(Zs) - 1, " composition axes, ",
          nrow(R), " samples")

  # Spearman on residuals
  C <- suppressWarnings(cor(R, method = "spearman"))
  diag(C) <- 0
  C[!is.finite(C)] <- 0

  # Threshold to the requested number of edges, matching STRING density
  ut <- which(upper.tri(C), arr.ind = TRUE)
  vals <- abs(C[ut])
  if (n_edges >= length(vals)) n_edges <- floor(length(vals) * 0.05)
  cut <- sort(vals, decreasing = TRUE)[n_edges]
  sel <- ut[vals >= cut, , drop = FALSE]

  edges <- tibble(from = colnames(C)[sel[, 1]],
                  to   = colnames(C)[sel[, 2]],
                  rho  = C[sel])
  log_msg(label, ": ", nrow(edges), " edges at |rho| >= ", round(cut, 3))

  gr <- graph_from_data_frame(edges %>% select(from, to),
                              directed = FALSE, vertices = g)
  list(graph = gr, edges = edges, expr = lm_[g, , drop = FALSE])
}

main_17 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene
  net   <- readRDS(P("rds", "network.rds"))
  score <- as.integer(net$chosen)

  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  global_deg <- setNames(norm$global_degree, norm$gene)

  # ---- STRING reference -------------------------------------------------
  sf <- file.path(CACHE_DIR, "string_density", paste0("t", score, ".rds"))
  se <- if (file.exists(sf)) readRDS(sf) else
        readRDS(P("rds", paste0("string_", score, ".rds")))
  if (is.null(se)) stop("No cached STRING edges.")
  se <- se %>% filter(from %in% genes, to %in% genes)
  sg <- graph_from_data_frame(se, directed = FALSE, vertices = genes)
  s_hub <- hub_select(sg)
  n_edges <- ecount(sg)
  log_msg("STRING network: ", n_edges, " edges, ", length(s_hub$hubs),
          " hubs: ", paste(s_hub$hubs, collapse = ", "))

  # ---- Co-expression from each cohort ------------------------------------
  cohorts <- get_validation_cohorts()
  results <- list()

  for (nm in names(cohorts)) {
    mat <- cohorts[[nm]]$matrix
    if (ncol(mat) < 100) { log_msg("Skipping ", nm, ", n = ", ncol(mat)); next }

    co <- tryCatch(build_coexpression(mat, genes, n_edges, nm),
                   error = function(e) { log_msg(nm, ": ",
                                                 conditionMessage(e)); NULL })
    if (is.null(co)) next

    h <- hub_select(co$graph)
    log_msg(nm, " co-expression hubs: ",
            paste(head(h$hubs, 12), collapse = ", "))

    # Alternative biases: does expression level or variance predict hubs?
    mean_expr <- setNames(rowMeans(co$expr), rownames(co$expr))
    var_expr  <- setNames(apply(co$expr, 1, sd), rownames(co$expr))

    results[[nm]] <- tibble(
      network = paste0("co-expression (", nm, ")"),
      n_samples = ncol(mat),
      n_edges = ecount(co$graph),
      n_hubs = length(h$hubs),
      auc_global_string_degree = predictability_auc(h$hubs,
                                    V(co$graph)$name, global_deg),
      auc_mean_expression = predictability_auc(h$hubs,
                                    V(co$graph)$name, mean_expr),
      auc_expression_variance = predictability_auc(h$hubs,
                                    V(co$graph)$name, var_expr),
      overlap_with_string_hubs = length(intersect(h$hubs, s_hub$hubs)),
      hubs = paste(head(h$hubs, 10), collapse = "; "))

    write_csv(co$edges, P("tables",
              paste0("coexpression_edges_", make.names(nm), ".csv")))
  }

  # STRING row, scored the same way
  s_mean <- NULL
  if (length(cohorts)) {
    m1 <- cohorts[[1]]$matrix
    gg <- intersect(genes, rownames(m1))
    s_mean <- setNames(rowMeans(log2(m1[gg, , drop = FALSE] + 1)), gg)
  }
  string_row <- tibble(
    network = "STRING (score 400)",
    n_samples = NA_integer_,
    n_edges = n_edges,
    n_hubs = length(s_hub$hubs),
    auc_global_string_degree = predictability_auc(s_hub$hubs, genes,
                                                  global_deg),
    auc_mean_expression = if (!is.null(s_mean))
      predictability_auc(s_hub$hubs, names(s_mean), s_mean) else NA_real_,
    auc_expression_variance = NA_real_,
    overlap_with_string_hubs = length(s_hub$hubs),
    hubs = paste(s_hub$hubs, collapse = "; "))

  out <- bind_rows(string_row, bind_rows(results))
  write_csv(out, P("tables", "coexpression_vs_string.csv"))

  log_msg("=================================================")
  log_msg("DATA-DERIVED NETWORK VERSUS STRING")
  log_msg("  auc_global_string_degree: can hub membership be predicted from")
  log_msg("  how well studied a gene is? 0.5 means no relationship.")
  log_msg("  The other two AUCs test whether a DIFFERENT bias replaces it.")
  log_msg("=================================================")
  print(as.data.frame(out %>% select(-hubs) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))))

  for (i in seq_len(nrow(out)))
    log_msg(out$network[i], ": ", out$hubs[i])

  co_rows <- out %>% filter(grepl("co-expression", network))
  if (nrow(co_rows)) {
    m <- mean(co_rows$auc_global_string_degree, na.rm = TRUE)
    log_msg("Mean AUC for co-expression networks: ", round(m, 3),
            " versus ", round(string_row$auc_global_string_degree, 3),
            " for STRING")
    if (!is.na(m) && m < 0.7) {
      log_msg("CONCLUSION: hubs from data-derived networks are not ",
              "predictable from annotation density. Building the network ",
              "from the cohort removes the bias that no metric could.")
      me <- mean(co_rows$auc_mean_expression, na.rm = TRUE)
      mv <- mean(co_rows$auc_expression_variance, na.rm = TRUE)
      log_msg("Check for a replacement bias: mean expression AUC ",
              round(me, 3), ", expression variance AUC ", round(mv, 3))
      if ((!is.na(me) && me > 0.8) || (!is.na(mv) && mv > 0.8)) {
        log_msg("WARNING: one bias has been replaced by another. Report ",
                "both AUCs rather than only the favourable one.")
      }
    } else {
      log_msg("CONCLUSION: co-expression hubs remain predictable from ",
              "STRING degree (AUC ", round(m, 3), "). Building the network ",
              "from data does NOT solve the problem, and that must be ",
              "reported.")
    }
  }

  saveRDS(out, P("rds", "coexpression_vs_string.rds"))
  log_msg("17_coexpression_network complete.")
  invisible(out)
}
