# Information-theoretic gene selection, as a prior-check.
#
# Question: can a method that never touches STRING or gene2pubmed select
# disease-informative genes that are (a) decoupled from citation count and
# (b) reproducible across independent cohorts? If so, it is a candidate
# replacement for the network-centrality workflow, immune to citation bias
# by construction rather than by correction.
#
# Method: minimum-redundancy-maximum-relevance (mRMR). Rank genes by mutual
# information with the disease target, penalised by redundancy with genes
# already chosen. Relevance rewards information about the phenotype;
# redundancy penalises genes that merely restate an already-selected one.
# This is the opposite of centrality, which rewards connectedness. No graph
# is built and the citation vector never enters the selection.
#
# Targets per cohort:
#   CGGA-325, CGGA-693   grade (LGG vs HGG), which varies
#   TCGA-LGG             IDH status, since grade does not vary in this
#                        pipeline's TCGA-LGG (all lower grade)
#
# Criteria, fixed before running:
#   KILL   citation AUC of the selected set > 0.65. The method found
#          high-variance well-studied genes anyway.
#   KILL   cross-cohort set agreement no better than the cytoHubba hub-list
#          agreement of 0.061 (F4). Immune to bias but unstable, so there
#          is no method.
#   PASS   citation AUC <= 0.55 AND cross-cohort agreement clearly above
#          0.061.
#   Between the two: inconclusive, report as such.
#
# Comparators: the same selection by relevance alone (mutual information,
# no redundancy penalty) and by variance alone. If mRMR does not beat
# these, the redundancy penalty is not what is doing the work.
#
# Input: cohort matrices and grade via 00c_cohorts.R; citation degrees via
# propagation_gates.rds (C aligned to genes) for the audit comparison only,
# never used in selection.

suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(readr); library(tibble)
})

N_SELECT   <- 50
N_BINS     <- 6      # discretisation bins for mutual information
SEED       <- 11
KILL_CIT   <- 0.65
PASS_CIT   <- 0.55
CYTO_AGREE <- 0.061  # F4 hub-list cross-cohort agreement

# ---- mutual information on discretised expression ----------------------
discretise <- function(x, nb) {
  br <- quantile(x, probs = seq(0, 1, length.out = nb + 1),
                 na.rm = TRUE, type = 8)
  br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(x, unique(br), labels = FALSE))
}

# ---- mutual information -------------------------------------------------
# MI between one discretised vector and every row of a discretised matrix,
mi_one_vs_many <- function(xd, Md, nb) {
  n <- length(xd)
  # indicator of x-bin. Sized by x's own levels, since the target may have
  # a different number of levels than the gene bins.
  nx <- max(xd)
  Ix <- matrix(0, n, nx); Ix[cbind(seq_len(n), xd)] <- 1
  px <- colSums(Ix) / n
  out <- numeric(nrow(Md))
  for (b in seq_len(nb)) {
    rows <- Md == b                      # nrow x n logical
    py_b <- rowSums(rows) / n            # P(row = b)
    if (all(py_b == 0)) next
    pj <- (rows %*% Ix) / n              # P(row = b, x = .)
    contrib <- pj * log2(pj / (py_b %o% px))
    contrib[!is.finite(contrib)] <- 0
    out <- out + rowSums(contrib)
  }
  out
}

# ---- mRMR selection -----------------------------------------------------
mrmr <- function(X, y, k, nb, prefilter = 400) {
  # X: genes in rows, samples in columns. y: target over samples.
  # Classic Peng 2005 criterion: relevance minus mean redundancy against
  # the already-selected set. Redundancy is a mean, on the same scale as
  # relevance. Candidates are prefiltered to the top `prefilter` by
  # relevance, since a redundancy penalty over thousands of irrelevant
  # genes is wasted work and is standard practice to drop.
  yd <- if (is.numeric(y) && length(unique(y)) > nb)
    discretise(y, nb) else as.integer(as.factor(y))
  Xd <- t(apply(X, 1, discretise, nb = nb))
  storage.mode(Xd) <- "integer"

  rel <- mi_one_vs_many(yd, Xd, nb)
  names(rel) <- rownames(X)

  keep <- order(rel, decreasing = TRUE)[seq_len(min(prefilter, nrow(X)))]
  Xk <- Xd[keep, , drop = FALSE]
  relk <- rel[keep]
  rel_only <- rownames(X)[order(rel, decreasing = TRUE)][1:k]

  sel <- which.max(relk)
  red_sum <- mi_one_vs_many(Xk[sel, ], Xk, nb)   # redundancy to first pick
  for (step in 2:k) {
    score <- relk - red_sum / (length(sel))
    score[sel] <- -Inf
    nxt <- which.max(score)
    sel <- c(sel, nxt)
    red_sum <- red_sum + mi_one_vs_many(Xk[nxt, ], Xk, nb)
  }
  list(genes = rownames(X)[keep][sel], relevance = rel,
       relevance_only = rel_only)
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

auc_rank <- function(score, pos) {
  r <- rank(score); n1 <- sum(pos); n0 <- length(pos) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

main_33 <- function() {
  source("R/00c_cohorts.R", local = FALSE)
  coh <- get_validation_cohorts()
  if (length(coh) < 2)
    stop("Need at least two cohorts with grade for the reproducibility ",
         "test. Found ", length(coh), ".")

  # Citation vector, for the audit comparison only. Aligned by gene symbol.
  g <- readRDS(P("rds", "propagation_gates.rds"))
  cit <- setNames(g$C, g$genes)

  crit <- tibble(
    criterion = c("citation AUC", "cross-cohort agreement"),
    pass_if = c(paste("<=", PASS_CIT),
                paste("clearly >", CYTO_AGREE)),
    kill_if = c(paste(">", KILL_CIT), paste("<=", CYTO_AGREE)))
  write_csv(crit, P("tables", "infoselect_criteria.csv"))
  log_msg("Criteria written before selection")

  set.seed(SEED)
  sets_mrmr <- list(); sets_rel <- list(); sets_var <- list()
  for (nm in names(coh)) {
    X <- coh[[nm]]$matrix
    meta <- coh[[nm]]$meta
    y <- meta$grade[match(colnames(X), meta$sample)]
    if (length(unique(na.omit(y))) < 2) {
      idh <- meta$IDH[match(colnames(X), meta$sample)]
      if (length(unique(na.omit(idh))) >= 2) {
        y <- idh
        log_msg(nm, ": grade constant, using IDH as target")
      } else {
        log_msg(nm, ": no varying target, skipped")
        next
      }
    }
    ok <- !is.na(y)
    X <- X[, ok, drop = FALSE]; y <- y[ok]
    # Keep the genes expressed and variable enough to carry information.
    v <- apply(X, 1, var)
    X <- X[v > quantile(v, 0.5) & rowMeans(X) > 1, , drop = FALSE]
    log_msg(nm, ": ", ncol(X), " samples, ", nrow(X),
            " genes after variance filter")

    res <- mrmr(X, y, N_SELECT, N_BINS)
    sets_mrmr[[nm]] <- res$genes
    sets_rel[[nm]]  <- res$relevance_only
    sets_var[[nm]]  <- rownames(X)[order(v[rownames(X)],
                                         decreasing = TRUE)][1:N_SELECT]
  }

  if (length(sets_mrmr) < 2)
    stop("Fewer than two cohorts produced a selection.")

  # ---- citation AUC of each cohort's mRMR set -------------------------
  log_msg("=============== citation AUC ===============")
  cit_auc <- c()
  for (nm in names(sets_mrmr)) {
    common <- intersect(names(cit), rownames(coh[[nm]]$matrix))
    pos <- common %in% sets_mrmr[[nm]]
    a <- auc_rank(cit[common], pos)
    cit_auc[nm] <- a
    log_msg(nm, ": AUC of citation predicting the mRMR set = ",
            round(a, 4))
  }

  # ---- cross-cohort agreement -----------------------------------------
  log_msg("========== cross-cohort agreement ==========")
  pairs <- combn(names(sets_mrmr), 2, simplify = FALSE)
  agree <- function(sets) mean(vapply(pairs,
    function(p) jaccard(sets[[p[1]]], sets[[p[2]]]), numeric(1)))
  a_mrmr <- agree(sets_mrmr); a_rel <- agree(sets_rel)
  a_var <- agree(sets_var)
  log_msg("mRMR            ", round(a_mrmr, 4))
  log_msg("relevance only  ", round(a_rel, 4))
  log_msg("variance only   ", round(a_var, 4))
  log_msg("cytoHubba (F4)  ", CYTO_AGREE, "  reference")

  # ---- verdict ---------------------------------------------------------
  mean_cit <- mean(cit_auc, na.rm = TRUE)
  out <- tibble(
    measure = c("mean citation AUC", "mRMR agreement",
                "relevance-only agreement", "variance-only agreement",
                "cytoHubba agreement"),
    value = c(mean_cit, a_mrmr, a_rel, a_var, CYTO_AGREE))
  write_csv(out, P("tables", "infoselect_result.csv"))

  cit_pass <- mean_cit <= PASS_CIT
  cit_kill <- mean_cit > KILL_CIT
  rep_kill <- a_mrmr <= CYTO_AGREE
  rep_pass <- a_mrmr > CYTO_AGREE * 2

  log_msg("================== SUMMARY =================")
  print(as.data.frame(out %>% mutate(value = round(value, 4))))
  log_msg("citation: ", if (cit_kill) "KILL" else if (cit_pass) "PASS"
          else "inconclusive",
          " (", round(mean_cit, 3), ")")
  log_msg("reproducibility: ", if (rep_kill) "KILL" else if (rep_pass)
          "PASS" else "inconclusive", " (", round(a_mrmr, 3),
          " vs ", CYTO_AGREE, ")")

  if (cit_pass && rep_pass)
    log_msg("Both criteria pass. Information selection is a candidate ",
            "replacement worth building properly.")
  else if (cit_kill || rep_kill)
    log_msg("A kill criterion fired. Report as a clean negative in the ",
            "audit paper and do not pursue the method.")
  else
    log_msg("Inconclusive. Neither pass nor kill; needs the full method ",
            "or a larger selection before deciding.")

  saveRDS(list(sets_mrmr = sets_mrmr, sets_rel = sets_rel,
               sets_var = sets_var, citation_auc = cit_auc,
               agreement = c(mrmr = a_mrmr, rel = a_rel, var = a_var),
               result = out), P("rds", "infoselect.rds"))
  invisible(out)
}
