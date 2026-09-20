# Eigenvector centrality on a network whose edge weights are cross-cohort
# reproducibility rather than correlation strength or database confidence.
#
# Eigenvector centrality on a database network is recursive degree and
# inherits any degree bias directly. Weighting edges by consistency across
# cohorts makes a gene central if it connects through reproducible edges to
# other genes that themselves connect through reproducible edges.
#
# Three weightings are compared against a plain mean-correlation control.
# Eigenvector centrality concentrates on dense subgraphs, so the
# correlation with plain degree and the component membership of the
# nominations are reported.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")

N_NOMINATE <- 8
EDGE_BUDGET <- 139     # matched to the STRING network for comparability

# --- Residualised correlation matrix per cohort -------------------------
residual_cor <- function(mat, genes, label) {
  g <- intersect(genes, rownames(mat))
  lm_ <- log2(mat + 1)
  neuro  <- composition_score(mat, NEURONAL_PANEL, exclude = g, label = "")
  glial  <- composition_score(mat, GLIAL_PANEL,  exclude = g, label = "")
  immune <- composition_score(mat, IMMUNE_PANEL, exclude = g, label = "")
  Z <- cbind(1, neuro, glial, immune)
  Z <- Z[, c(TRUE, apply(Z[, -1, drop = FALSE], 2,
                         function(v) all(is.finite(v)) && sd(v) > 0)),
         drop = FALSE]
  X <- t(lm_[g, , drop = FALSE])
  ok <- complete.cases(X) & complete.cases(Z)
  X <- X[ok, , drop = FALSE]; Z <- Z[ok, , drop = FALSE]
  R <- X - Z %*% qr.coef(qr(Z), X)
  C <- suppressWarnings(cor(R, method = "spearman"))
  C[!is.finite(C)] <- 0; diag(C) <- 0
  log_msg(label, ": ", length(g), " genes, ", nrow(R), " samples")
  C
}

# --- Reproducibility weights --------------------------------------------
build_weights <- function(mats, genes) {
  A <- simplify2array(lapply(mats, function(C) C[genes, genes, drop = FALSE]))
  n <- length(genes)

  rho_mean <- apply(A, c(1, 2), function(v) mean(v, na.rm = TRUE))
  rho_min  <- apply(A, c(1, 2), function(v) min(abs(v), na.rm = TRUE))
  z <- atanh(pmin(pmax(A, -0.999), 0.999))
  z_mean <- apply(z, c(1, 2), mean, na.rm = TRUE)
  z_sd   <- apply(z, c(1, 2), sd,   na.rm = TRUE)
  sign_agree <- apply(A, c(1, 2), function(v) {
    s <- sign(v); s <- s[s != 0]
    if (!length(s)) return(0)
    max(mean(s > 0), mean(s < 0))
  })

  W <- list(
    min     = rho_min,
    invvar  = abs(z_mean) / (1 + ifelse(is.na(z_sd), 0, z_sd)),
    concord = abs(rho_mean) * sign_agree,
    # reference: plain mean correlation, ignoring reproducibility
    meanrho = abs(rho_mean))

  W <- lapply(W, function(M) {
    dimnames(M) <- list(genes, genes); diag(M) <- 0
    M[!is.finite(M)] <- 0; M
  })
  W
}

threshold_to_budget <- function(M, budget) {
  ut <- which(upper.tri(M), arr.ind = TRUE)
  v <- M[ut]
  if (budget >= length(v)) budget <- floor(length(v) * 0.05)
  cut <- sort(v, decreasing = TRUE)[budget]
  Mt <- M; Mt[Mt < cut] <- 0
  Mt
}

eigen_scores <- function(M) {
  g <- graph_from_adjacency_matrix(M, mode = "undirected", weighted = TRUE,
                                   diag = FALSE)
  ev <- tryCatch(eigen_centrality(g, weights = E(g)$weight)$vector,
                 error = function(e) rep(NA_real_, vcount(g)))
  comp <- components(g)
  tibble(gene = V(g)$name, eigen = ev, degree = degree(g),
         wdegree = strength(g),
         component_size = comp$csize[comp$membership])
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

jaccard <- function(a, b) {
  u <- length(union(a, b)); if (!u) return(NA_real_)
  length(intersect(a, b)) / u
}

main_21 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene
  net   <- readRDS(P("rds", "network.rds"))
  conv  <- net$consensus$gene

  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  global_deg <- setNames(norm$global_degree, norm$gene)

  cohorts <- get_validation_cohorts()
  usable <- names(cohorts)[map_int(cohorts, ~ ncol(.x$matrix)) >= 150]
  log_msg("Cohorts: ", paste(usable, collapse = ", "))
  if (length(usable) < 3) stop("Need three cohorts for a reproducibility weight.")

  mats <- map(usable, ~ residual_cor(cohorts[[.x]]$matrix, genes, .x))
  names(mats) <- usable
  common <- sort(Reduce(intersect, map(mats, rownames)))
  log_msg("Common genes: ", length(common))

  W <- build_weights(mats, common)

  res <- map_dfr(names(W), function(nm) {
    Mt <- threshold_to_budget(W[[nm]], EDGE_BUDGET)
    sc <- eigen_scores(Mt)
    nom <- sc %>% arrange(desc(eigen)) %>% head(N_NOMINATE) %>% pull(gene)
    tibble(
      weight = nm,
      n_edges = sum(Mt > 0) / 2,
      auc_annotation_bias = bias_auc(nom, common, global_deg),
      median_global_degree = median(global_deg[match(nom, names(global_deg))],
                                    na.rm = TRUE),
      cor_eigen_degree = suppressWarnings(cor(sc$eigen, sc$degree,
                                              method = "spearman",
                                              use = "complete.obs")),
      median_component = median(sc$component_size[match(nom, sc$gene)]),
      overlap_cytoHubba = length(intersect(nom, conv)),
      genes = paste(nom, collapse = "; "))
  })

  write_csv(res, P("tables", "reproducibility_eigenvector.csv"))

  log_msg("=================================================")
  log_msg("EIGENVECTOR ON REPRODUCIBILITY-WEIGHTED NETWORKS")
  log_msg("  auc_annotation_bias LOWER is better; 0.5 ideal")
  log_msg("  benchmarks: cytoHubba 0.994, binomial 0.825,")
  log_msg("              co-expression hubs 0.576")
  log_msg("  cor_eigen_degree: if high, eigenvector is just degree again")
  log_msg("=================================================")
  print(as.data.frame(res %>% select(-genes) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))))
  for (i in seq_len(nrow(res)))
    log_msg("  ", res$weight[i], ": ", res$genes[i])

  # --- Leave-one-cohort-out stability of the nominations -----------------
  if (length(usable) >= 3) {
    stab <- map_dfr(names(W), function(nm) {
      full <- eigen_scores(threshold_to_budget(W[[nm]], EDGE_BUDGET)) %>%
        arrange(desc(eigen)) %>% head(N_NOMINATE) %>% pull(gene)
      js <- map_dbl(usable, function(drop) {
        sub <- mats[setdiff(usable, drop)]
        Ws <- build_weights(sub, common)[[nm]]
        n2 <- eigen_scores(threshold_to_budget(Ws, EDGE_BUDGET)) %>%
          arrange(desc(eigen)) %>% head(N_NOMINATE) %>% pull(gene)
        jaccard(n2, full)
      })
      tibble(weight = nm, loo_stability = mean(js, na.rm = TRUE))
    })
    write_csv(stab, P("tables", "reproducibility_eigenvector_stability.csv"))
    log_msg("Leave-one-cohort-out stability of nominations:")
    print(as.data.frame(stab %>% mutate(loo_stability =
                                          round(loo_stability, 3))))
  }

  # --- Verdict -----------------------------------------------------------
  best <- res[which.min(res$auc_annotation_bias), ]
  ref_coexp <- 0.576
  log_msg("---")
  log_msg("Best weighting: ", best$weight, " at AUC ",
          round(best$auc_annotation_bias, 3))

  if (!is.na(best$cor_eigen_degree) && best$cor_eigen_degree > 0.85) {
    log_msg("WARNING: eigenvector correlates with plain degree at ",
            round(best$cor_eigen_degree, 3),
            ". It is recovering degree, not propagating reproducibility.")
  }
  if (!is.na(best$auc_annotation_bias) &&
      best$auc_annotation_bias < ref_coexp - 0.05) {
    log_msg("This beats the co-expression hub benchmark of ", ref_coexp,
            ". Worth developing.")
  } else if (!is.na(best$auc_annotation_bias) &&
             best$auc_annotation_bias > ref_coexp + 0.05) {
    log_msg("This does NOT beat the co-expression hub benchmark of ",
            ref_coexp, ". Four scoring approaches have now landed between ",
            "0.73 and 0.83 while the substrate change reached 0.58. The ",
            "evidence says the scoring function is not where the problem ",
            "lives. Report that and stop.")
  } else {
    log_msg("Comparable to the co-expression benchmark. No clear gain.")
  }

  # Does the reproducibility weighting beat the plain-correlation control?
  mr <- res$auc_annotation_bias[res$weight == "meanrho"]
  bt <- best$auc_annotation_bias
  if (!is.na(mr) && !is.na(bt)) {
    log_msg("Reproducibility weighting versus plain mean correlation: ",
            round(bt, 3), " versus ", round(mr, 3),
            ifelse(bt < mr - 0.05,
                   ". The reproducibility term is doing real work.",
                   ". The reproducibility term adds little over plain correlation."))
  }

  saveRDS(list(results = res, weights = names(W)),
          P("rds", "reproducibility_eigenvector.rds"))
  log_msg("21_reproducibility_weighted_eigenvector complete.")
  invisible(res)
}
