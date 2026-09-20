# End-to-end comparison of two complete pipelines.
#
#   conventional  DEG -> STRING -> centrality consensus -> hub list
#   proposed      composition-adjusted DEG -> cohort co-expression ->
#                 cross-cohort agreement weights -> Louvain modules
#
# Judged on annotation bias, leave-one-cohort-out reproducibility, and
# enrichment. The third criterion guards against a pipeline that is
# reproducibly finding nothing.
#
# The conventional pipeline uses STRING, which does not change when a
# cohort is dropped, so its leave-one-out agreement is 1.0 by construction
# and is not comparable.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC)
  library(clusterProfiler); library(org.Hs.eg.db)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")

EDGE_BUDGET <- 139
MIN_MODULE  <- 5

jaccard <- function(a, b) {
  u <- length(union(a, b)); if (!u) return(NA_real_)
  length(intersect(a, b)) / u
}

residual_cor <- function(mat, genes, label = "") {
  g <- intersect(genes, rownames(mat))
  lm_ <- log2(mat + 1)
  Z <- cbind(1,
             composition_score(mat, NEURONAL_PANEL, exclude = g, label = ""),
             composition_score(mat, GLIAL_PANEL,  exclude = g, label = ""),
             composition_score(mat, IMMUNE_PANEL, exclude = g, label = ""))
  Z <- Z[, c(TRUE, apply(Z[, -1, drop = FALSE], 2,
                         function(v) all(is.finite(v)) && sd(v) > 0)),
         drop = FALSE]
  X <- t(lm_[g, , drop = FALSE])
  ok <- complete.cases(X) & complete.cases(Z)
  R <- X[ok, , drop = FALSE] - Z[ok, , drop = FALSE] %*%
       qr.coef(qr(Z[ok, , drop = FALSE]), X[ok, , drop = FALSE])
  C <- suppressWarnings(cor(R, method = "spearman"))
  C[!is.finite(C)] <- 0; diag(C) <- 0
  C
}

# Edge weight = minimum absolute correlation across cohorts. An edge counts
# only if it holds in every cohort, which the selector comparison found
# removes annotation bias relative to mean correlation.
agreement_weights <- function(mats, genes) {
  A <- simplify2array(lapply(mats, function(C) C[genes, genes, drop = FALSE]))
  W <- apply(A, c(1, 2), function(v) min(abs(v), na.rm = TRUE))
  dimnames(W) <- list(genes, genes); diag(W) <- 0
  W[!is.finite(W)] <- 0
  W
}

threshold_to_budget <- function(M, budget) {
  ut <- which(upper.tri(M), arr.ind = TRUE)
  v <- M[ut]
  if (budget >= length(v)) budget <- floor(length(v) * 0.05)
  cut <- sort(v, decreasing = TRUE)[budget]
  M[M < cut] <- 0
  M
}

modules_of <- function(W, min_size = MIN_MODULE) {
  g <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                   diag = FALSE)
  g <- delete_vertices(g, which(degree(g) == 0))
  if (vcount(g) < min_size) return(list())
  cl <- cluster_louvain(g, weights = E(g)$weight)
  mods <- split(V(g)$name, membership(cl))
  mods[lengths(mods) >= min_size]
}

bias_auc <- function(nominated, all_genes, global_deg) {
  d <- tibble(gene = all_genes,
              y = as.integer(all_genes %in% nominated),
              x = global_deg[match(all_genes, names(global_deg))]) %>%
    filter(!is.na(x))
  if (sum(d$y) < 3 || sum(d$y) == nrow(d)) return(NA_real_)
  fit <- glm(y ~ log10(x + 1), data = d, family = binomial())
  as.numeric(pROC::auc(pROC::roc(d$y, predict(fit, type = "response"),
                                 quiet = TRUE)))
}

enrichment_of <- function(genes, universe) {
  ez <- suppressWarnings(suppressMessages(
    bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)))$ENTREZID
  uz <- suppressWarnings(suppressMessages(
    bitr(universe, "SYMBOL", "ENTREZID", org.Hs.eg.db)))$ENTREZID
  if (length(ez) < 5) return(list(n = 0, top = NA_character_))
  k <- tryCatch(enrichGO(ez, universe = uz, OrgDb = org.Hs.eg.db,
                         ont = "BP", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable = TRUE),
                error = function(e) NULL)
  if (is.null(k)) return(list(n = 0, top = NA_character_))
  d <- as.data.frame(k)
  list(n = nrow(d), top = if (nrow(d)) d$Description[1] else NA_character_)
}

main_22 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes_adj  <- inter$shared          # composition-robust
  genes_conv <- inter$conventional    # unadjusted intersection
  if (is.data.frame(genes_adj))  genes_adj  <- genes_adj$gene
  if (is.data.frame(genes_conv)) genes_conv <- genes_conv$gene
  log_msg("Input sets: adjusted ", length(genes_adj),
          ", conventional ", length(genes_conv))

  net  <- readRDS(P("rds", "network.rds"))
  conv_hubs <- net$consensus$gene
  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  global_deg <- setNames(norm$global_degree, norm$gene)

  micro <- readRDS(P("rds", "deg_microarray.rds"))
  rna   <- readRDS(P("rds", "deg_rnaseq.rds"))
  universe <- unique(c(unlist(lapply(micro, function(x)
                        if (!is.null(x)) x$gene)), rna$LGG_vs_Normal$gene))
  universe <- universe[!is.na(universe) & universe != ""]

  cohorts <- get_validation_cohorts()
  usable <- names(cohorts)[map_int(cohorts, ~ ncol(.x$matrix)) >= 150]
  if (length(usable) < 3) stop("Need three cohorts.")
  log_msg("Cohorts: ", paste(usable, collapse = ", "))

  # ---- PROPOSED pipeline, full and leave-one-out ------------------------
  run_proposed <- function(cohort_subset) {
    mats <- map(cohort_subset, ~ residual_cor(cohorts[[.x]]$matrix, genes_adj))
    common <- sort(Reduce(intersect, map(mats, rownames)))
    W <- threshold_to_budget(agreement_weights(mats, common), EDGE_BUDGET)
    mods <- modules_of(W)
    list(modules = mods, genes = unlist(mods, use.names = FALSE),
         W = W, common = common)
  }

  prop <- run_proposed(usable)
  log_msg("Proposed pipeline: ", length(prop$modules), " modules, ",
          length(prop$genes), " genes in modules")
  for (i in seq_along(prop$modules))
    log_msg("  module ", i, " (n=", length(prop$modules[[i]]), "): ",
            paste(head(prop$modules[[i]], 12), collapse = ", "))

  # ---- Reproducibility, leave one cohort out ----------------------------
  loo_prop <- map_dbl(usable, function(drop) {
    p2 <- run_proposed(setdiff(usable, drop))
    # module-level agreement: best-matching module Jaccard, averaged
    if (!length(p2$modules) || !length(prop$modules)) return(NA_real_)
    mean(map_dbl(prop$modules, function(m)
      max(map_dbl(p2$modules, ~ jaccard(m, .x)), na.rm = TRUE)),
      na.rm = TRUE)
  })
  log_msg("Proposed, leave-one-cohort-out module agreement: ",
          paste(round(loo_prop, 3), collapse = ", "),
          " (mean ", round(mean(loo_prop, na.rm = TRUE), 3), ")")

  # ---- CONVENTIONAL pipeline, leave one cohort out ----------------------
  # The conventional pipeline uses STRING, which does not depend on the
  # cohorts, so its output cannot change when a cohort is dropped. Its
  # apparent stability is therefore 1.0 by construction and is reported as
  # such rather than as a merit.
  log_msg("Conventional pipeline output is cohort independent by ",
          "construction: STRING topology does not change when a cohort is ",
          "removed. Its leave-one-out agreement is 1.0 trivially and is not ",
          "comparable to the proposed pipeline's.")

  # ---- Evaluation -------------------------------------------------------
  all_genes <- intersect(union(genes_adj, genes_conv), names(global_deg))

  e_conv <- enrichment_of(conv_hubs, universe)
  e_prop <- enrichment_of(prop$genes, universe)

  # largest module alone, as a like-for-like short list
  biggest <- if (length(prop$modules))
    prop$modules[[which.max(lengths(prop$modules))]] else character(0)
  e_big <- enrichment_of(biggest, universe)

  ev <- tibble(
    pipeline = c("conventional (STRING + cytoHubba hubs)",
                 "proposed (all module genes)",
                 "proposed (largest module)"),
    n_output = c(length(conv_hubs), length(prop$genes), length(biggest)),
    auc_annotation_bias = c(
      bias_auc(conv_hubs, all_genes, global_deg),
      bias_auc(prop$genes, all_genes, global_deg),
      bias_auc(biggest, all_genes, global_deg)),
    median_global_degree = c(
      median(global_deg[match(conv_hubs, names(global_deg))], na.rm = TRUE),
      median(global_deg[match(prop$genes, names(global_deg))], na.rm = TRUE),
      median(global_deg[match(biggest, names(global_deg))], na.rm = TRUE)),
    loo_reproducibility = c(NA_real_, mean(loo_prop, na.rm = TRUE),
                            mean(loo_prop, na.rm = TRUE)),
    n_enriched_GO_BP = c(e_conv$n, e_prop$n, e_big$n),
    top_term = c(e_conv$top, e_prop$top, e_big$top))

  write_csv(ev, P("tables", "pipeline_benchmark.csv"))
  write_csv(tibble(module = rep(seq_along(prop$modules),
                                lengths(prop$modules)),
                   gene = unlist(prop$modules, use.names = FALSE)),
            P("tables", "proposed_pipeline_modules.csv"))

  log_msg("=================================================")
  log_msg("END-TO-END PIPELINE COMPARISON")
  log_msg("  annotation bias: LOWER better, 0.5 ideal")
  log_msg("  reproducibility: HIGHER better, oracle ceiling near 0.51")
  log_msg("  enrichment: guards against reproducibly finding nothing")
  log_msg("=================================================")
  print(as.data.frame(ev %>% select(-top_term) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))))
  for (i in seq_len(nrow(ev)))
    log_msg("  ", ev$pipeline[i], " top term: ", ev$top_term[i])

  # ---- Verdict ----------------------------------------------------------
  b_conv <- ev$auc_annotation_bias[1]
  b_prop <- ev$auc_annotation_bias[2]
  rep_prop <- mean(loo_prop, na.rm = TRUE)
  enr_prop <- ev$n_enriched_GO_BP[2]

  log_msg("---")
  if (!is.na(b_prop) && !is.na(b_conv) && b_prop < b_conv - 0.15 &&
      !is.na(rep_prop) && rep_prop > 0.3 && enr_prop > 0) {
    log_msg("The proposed pipeline reduces annotation bias from ",
            round(b_conv, 3), " to ", round(b_prop, 3),
            ", reproduces across cohorts at ", round(rep_prop, 3),
            " against a ceiling near 0.51, and retains biological signal (",
            enr_prop, " enriched GO BP terms). All three criteria met.")
  } else if (enr_prop == 0) {
    log_msg("The proposed pipeline produces no enriched terms. Lower bias ",
            "and higher stability without biological signal means it is ",
            "reproducibly finding nothing. Report this and do not ",
            "recommend the pipeline.")
  } else {
    log_msg("Mixed result. bias ", round(b_prop, 3), " vs ",
            round(b_conv, 3), "; reproducibility ", round(rep_prop, 3),
            "; enriched terms ", enr_prop,
            ". Report all three figures and temper the recommendation.")
  }

  saveRDS(list(evaluation = ev, modules = prop$modules, loo = loo_prop),
          P("rds", "pipeline_benchmark.rds"))
  log_msg("22_pipeline_benchmark complete.")
  invisible(ev)
}
