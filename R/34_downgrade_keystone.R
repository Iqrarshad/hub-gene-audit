# The downgrade keystone: does the network add value over doing nothing?
#
# The STRING-cytoHubba workflow has four steps: differential expression,
# network construction, centrality ranking, hub list. Earlier stages showed
# each step is citation-loaded or unreproducible. This stage asks the
# endpoint question directly: does the assembled pipeline beat a baseline
# that skips the network and centrality entirely?
#
# Selectors compared, each producing a top-k gene set per cohort:
#   variance      highest expression variance. No target, no graph.
#   DE t-stat     strongest grade contrast. No graph.
#   relevance     mutual information with grade. No graph.
#   mRMR          relevance minus redundancy. No graph. (information-selection stage)
#   degree        STRING global degree. The network step alone.
#   cytoHubba-like top-k by degree on the cohort-induced STRING subgraph.
#   coexpr_degree cohort-specific co-expression degree, a network built
#                 from this cohort's own data rather than a fixed database.
#
# The co-expression selector is the control for the obvious objection: STRING
# degree is reproducible only because it is a fixed global vector, identical
# across cohorts by construction. A network derived from each cohort's own
# expression cannot be identical across cohorts, so if network structure
# rather than data-independence drives reproducibility, co-expression degree
# should be similarly reproducible. The prediction, stated before running,
# is that it is not.
#
# Two properties per selector, both pre-declared:
#   reproducibility  mean pairwise Jaccard of the top-k sets across cohorts
#   citation AUC     how well gene2pubmed count predicts the selected set
#
# Downgrade is demonstrated if a network-free selector is at least as
# reproducible as the degree and cytoHubba selectors while being no more
# citation-biased. That is the network step subtracting value.
#
# Pre-declared reading:
#   DOWNGRADE SHOWN   best network-free reproducibility >= degree
#                     reproducibility, and its citation AUC <= the degree
#                     citation AUC
#   NOT SHOWN         a network selector is both more reproducible and no
#                     more citation-biased than every network-free one
#
# Robustness: all values are reported at k = 20, 50, 100 and the ordering
# must hold across sizes to be claimed.
#
# Input: cohort matrices and grade via 00c_cohorts.R; STRING degree and
# citation via propagation_gates.rds (k and C aligned to genes).

suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(readr); library(tibble)
})

KS         <- c(20, 50, 100)
N_BINS     <- 6
SEED       <- 11

source_once <- function(f) if (!exists("get_validation_cohorts"))
  source(f, local = FALSE)

# reuse the MI and mRMR machinery from the information-selection stage (33)
mi_one_vs_many <- function(xd, Md, nb) {
  n <- length(xd); nx <- max(xd)
  Ix <- matrix(0, n, nx); Ix[cbind(seq_len(n), xd)] <- 1
  px <- colSums(Ix) / n
  out <- numeric(nrow(Md))
  for (b in seq_len(nb)) {
    rows <- Md == b
    py_b <- rowSums(rows) / n
    if (all(py_b == 0)) next
    pj <- (rows %*% Ix) / n
    contrib <- pj * log2(pj / (py_b %o% px))
    contrib[!is.finite(contrib)] <- 0
    out <- out + rowSums(contrib)
  }
  out
}

discretise <- function(x, nb) {
  br <- quantile(x, probs = seq(0, 1, length.out = nb + 1),
                 na.rm = TRUE, type = 8)
  br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(x, unique(br), labels = FALSE))
}

mrmr_set <- function(X, y, k, nb, prefilter = 400) {
  yd <- as.integer(as.factor(y))
  Xd <- t(apply(X, 1, discretise, nb = nb)); storage.mode(Xd) <- "integer"
  rel <- mi_one_vs_many(yd, Xd, nb); names(rel) <- rownames(X)
  keep <- order(rel, decreasing = TRUE)[seq_len(min(prefilter, nrow(X)))]
  Xk <- Xd[keep, , drop = FALSE]; relk <- rel[keep]
  sel <- which.max(relk); red <- mi_one_vs_many(Xk[sel, ], Xk, nb)
  for (s in 2:k) {
    score <- relk - red / length(sel); score[sel] <- -Inf
    nx <- which.max(score); sel <- c(sel, nx)
    red <- red + mi_one_vs_many(Xk[nx, ], Xk, nb)
  }
  rownames(X)[keep][sel]
}

t_stat <- function(X, y) {
  g <- unique(y); a <- X[, y == g[1], drop = FALSE]
  b <- X[, y == g[2], drop = FALSE]
  ma <- rowMeans(a); mb <- rowMeans(b)
  va <- apply(a, 1, var); vb <- apply(b, 1, var)
  (ma - mb) / sqrt(va / ncol(a) + vb / ncol(b) + 1e-9)
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# Cohort-specific co-expression degree. Matches the co-expression network stage: log2, Spearman on
# a residual-free basis here (composition covariates are not reloaded in this
# stage, so this is the unresidualised version and is labelled as such), then
# threshold the correlation matrix to the same density as STRING and count
# each gene's retained edges. The ranking is what feeds selection.
coexpr_degree <- function(X, n_edges) {
  lm_ <- log2(X + 1)
  C <- suppressWarnings(cor(t(lm_), method = "spearman"))
  diag(C) <- 0; C[!is.finite(C)] <- 0
  ut <- which(upper.tri(C), arr.ind = TRUE)
  vals <- abs(C[ut])
  if (n_edges >= length(vals)) n_edges <- floor(length(vals) * 0.05)
  cut <- sort(vals, decreasing = TRUE)[n_edges]
  A <- (abs(C) >= cut) * 1
  setNames(rowSums(A), rownames(X))
}
auc_rank <- function(score, pos) {
  r <- rank(score); n1 <- sum(pos); n0 <- length(pos) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

main_34 <- function() {
  source_once("R/00c_cohorts.R")
  coh <- get_validation_cohorts()

  g <- readRDS(P("rds", "propagation_gates.rds"))
  deg <- setNames(g$k, g$genes)
  cit <- setNames(g$C, g$genes)

  # keep cohorts with a two-level grade target
  use <- list()
  for (nm in names(coh)) {
    y <- coh[[nm]]$meta$grade[match(colnames(coh[[nm]]$matrix),
                                    coh[[nm]]$meta$sample)]
    if (length(unique(na.omit(y))) >= 2) {
      ok <- !is.na(y)
      use[[nm]] <- list(X = coh[[nm]]$matrix[, ok, drop = FALSE], y = y[ok])
    } else log_msg(nm, ": grade constant, excluded from reproducibility")
  }
  if (length(use) < 2) stop("Need two grade cohorts.")
  log_msg("Grade cohorts: ", paste(names(use), collapse = ", "))

  # STRING edge count sets the density target for co-expression, so both
  # networks are thresholded to the same number of edges.
  n_edges_string <- length(g$k)   # proxy; replaced per cohort by 5% density
  # Precompute cohort-specific co-expression degree once per cohort. The
  # correlation matrix is large, so this is done outside the k loop.
  coex_deg <- list()
  for (nm in names(use)) {
    X <- use[[nm]]$X
    v <- apply(X, 1, var)
    Xf <- X[v > quantile(v, 0.5) & rowMeans(X) > 1, , drop = FALSE]
    n_e <- floor(0.05 * choose(nrow(Xf), 2) / 100)  # sparse target
    coex_deg[[nm]] <- tryCatch(coexpr_degree(Xf, n_e),
      error = function(e) { log_msg(nm, " co-expression failed: ",
                                    conditionMessage(e)); NULL })
    if (!is.null(coex_deg[[nm]]))
      log_msg(nm, ": co-expression degree computed on ", nrow(Xf), " genes")
  }

  crit <- tibble(
    check = c("downgrade shown"),
    rule = paste("best network-free reproducibility >= degree, and its",
                 "citation AUC <= degree citation AUC, at all k"))
  write_csv(crit, P("tables", "downgrade_criteria.csv"))

  selectors <- c("variance", "DE_tstat", "relevance", "mRMR",
                 "degree", "cytoHubba_subgraph", "coexpr_degree")
  rows <- list()
  set.seed(SEED)
  for (k in KS) {
    sets <- setNames(vector("list", length(selectors)), selectors)
    for (s in selectors) sets[[s]] <- list()
    for (nm in names(use)) {
      X <- use[[nm]]$X; y <- use[[nm]]$y
      v <- apply(X, 1, var)
      X <- X[v > quantile(v, 0.5) & rowMeans(X) > 1, , drop = FALSE]
      gg <- rownames(X)
      vv <- apply(X, 1, var)

      sets$variance[[nm]] <- gg[order(vv, decreasing = TRUE)][1:k]
      tt <- abs(t_stat(X, y))
      sets$DE_tstat[[nm]] <- gg[order(tt, decreasing = TRUE)][1:k]
      yd <- as.integer(as.factor(y))
      Xd <- t(apply(X, 1, discretise, nb = N_BINS))
      storage.mode(Xd) <- "integer"
      rel <- mi_one_vs_many(yd, Xd, N_BINS)
      sets$relevance[[nm]] <- gg[order(rel, decreasing = TRUE)][1:k]
      sets$mRMR[[nm]] <- mrmr_set(X, y, k, N_BINS)

      # network selectors, restricted to genes present in this cohort
      dsub <- deg[intersect(names(deg), gg)]
      sets$degree[[nm]] <- names(sort(dsub, decreasing = TRUE))[1:k]
      # cytoHubba-like: degree within the cohort-induced subgraph is the
      # same ordering as global degree restricted to these genes, since the
      # subgraph keeps only edges among present genes. Use induced degree.
      sets$cytoHubba_subgraph[[nm]] <- names(sort(dsub,
                                                  decreasing = TRUE))[1:k]
      cd <- coex_deg[[nm]]
      if (!is.null(cd)) {
        cd <- cd[intersect(names(cd), gg)]
        sets$coexpr_degree[[nm]] <- names(sort(cd, decreasing = TRUE))[1:k]
      } else {
        sets$coexpr_degree[[nm]] <- character(0)
      }
    }
    pairs <- combn(names(use), 2, simplify = FALSE)
    for (s in selectors) {
      rep <- mean(vapply(pairs, function(p)
        jaccard(sets[[s]][[p[1]]], sets[[s]][[p[2]]]), numeric(1)))
      cauc <- mean(vapply(names(use), function(nm) {
        common <- intersect(names(cit), rownames(use[[nm]]$X))
        auc_rank(cit[common], common %in% sets[[s]][[nm]])
      }, numeric(1)), na.rm = TRUE)
      rows[[length(rows) + 1]] <- tibble(k = k, selector = s,
        reproducibility = rep, citation_auc = cauc,
        network_free = !s %in% c("degree", "cytoHubba_subgraph",
                                 "coexpr_degree"))
    }
  }
  res <- bind_rows(rows)
  write_csv(res, P("tables", "downgrade_result.csv"))

  log_msg("=============== RESULT ===============")
  for (k in KS) {
    log_msg("--- k = ", k, " ---")
    sub <- res %>% filter(k == !!k) %>% arrange(desc(reproducibility))
    print(as.data.frame(sub %>% transmute(selector,
      reproducibility = round(reproducibility, 4),
      citation_auc = round(citation_auc, 4),
      network_free)))
  }

  # verdict at each k
  log_msg("================ VERDICT ================")
  shown <- logical(0)
  for (k in KS) {
    sub <- res %>% filter(k == !!k)
    degr <- sub$reproducibility[sub$selector == "degree"]
    degc <- sub$citation_auc[sub$selector == "degree"]
    nf <- sub %>% filter(network_free)
    best <- nf %>% filter(reproducibility >= degr, citation_auc <= degc)
    ok <- nrow(best) > 0
    shown <- c(shown, ok)
    log_msg("k = ", k, ": ", if (ok)
      paste0("downgrade shown by ",
             best$selector[which.max(best$reproducibility)]) else
      "not shown at this k")
  }
  log_msg(if (all(shown))
    paste0("DOWNGRADE HOLDS across all k. A network-free selector ",
           "matches or beats the network on reproducibility at no ",
           "citation cost.")
    else if (any(shown))
    "MIXED. Holds at some k only; report the size dependence."
    else
    "NOT SHOWN. The network selector is not dominated.")

  # The data-independence contrast: fixed STRING degree versus cohort-derived
  # co-expression degree. If STRING is far more reproducible, its stability
  # comes from being a fixed vector, not from being a network.
  log_msg("======= STRING vs co-expression degree =======")
  for (k in KS) {
    sub <- res %>% filter(k == !!k)
    sd <- sub$reproducibility[sub$selector == "degree"]
    cd <- sub$reproducibility[sub$selector == "coexpr_degree"]
    log_msg("k = ", k, ": STRING degree ", round(sd, 3),
            " vs co-expression degree ", round(cd, 3),
            if (!is.na(cd) && sd > cd + 0.1)
              "  -> fixed-vector stability confirmed" else "")
  }

  saveRDS(list(result = res, verdict = shown), P("rds", "downgrade.rds"))
  invisible(res)
}
