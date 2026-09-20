# Why does propagation recover degree rather than signal?
#
# The propagation gates found held-out recovery of 0.066 against a degree
# recovery of 0.61 to 0.88. This stage separates the candidate causes before
# that is read as a property of the data.
#
# D1  Graph autocorrelation of E. Correlation between a gene's effect and
#     the mean effect of its neighbours, against a degree-stratified
#     permutation null. This bounds what any propagation method can reach.
# D2  The same restricted to genes carrying real effect, and as
#     discrimination of the top decile, since a Spearman over the noise
#     bulk is driven toward zero mechanically.
# D3  Signal choice. Absolute logFC, signed logFC and significance are not
#     equally smooth on a graph.
# D4  Normalisation. Symmetric D^-1/2 A D^-1/2 carries a sqrt(degree)
#     factor that row normalisation D^-1 A does not. The propagation gates
#     used row normalisation; the earlier symmetric run used D^-1/2 A D^-1/2.
# D5  Recovery after degree is partialled out.
# D6  Calibration. Signals with graded graph smoothness are generated,
#     masked and recovered. Recovery is read against autocorrelation, and
#     the observed autocorrelation from D1 is located on that curve. A
#     single pass or fail control cannot say whether observed recovery is
#     low in absolute terms or low for the smoothness available.
#
# Reading, fixed before running:
#   D6 fails                        the test is invalid, stop and fix it
#   D1 within the null              signal has no graph structure, so
#                                   propagation cannot work on this
#                                   substrate, and the operator is not at
#                                   fault
#   D1 above null but D4 differs    normalisation was the problem
#   D1 above null and D2 strong     the metric was the problem
#
# Input: propagation_gates.rds from the propagation-gates stage (27), deg_rnaseq.rds for D3.

suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(readr); library(tibble)
})

N_PERM    <- 200
N_BINS    <- 20
TOP_FRAC  <- 0.10
MASK_FRAC <- 0.20
SEED      <- 11
ALPHA     <- 1

# ---- helpers ------------------------------------------------------------
auc_rank <- function(score, pos) {
  r <- rank(score)
  n1 <- sum(pos); n0 <- length(pos) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

neighbour_mean <- function(A, v) {
  s <- Matrix::rowSums(A)
  as.vector(A %*% v) / pmax(s, 1e-12)
}

# Permute within degree bins so the null keeps the degree distribution and
# destroys only the alignment between signal and graph position.
degree_binned_perm <- function(v, k, n_bins) {
  b <- cut(rank(k, ties.method = "first"), breaks = n_bins, labels = FALSE)
  out <- v
  for (i in unique(b)) {
    idx <- which(b == i)
    out[idx] <- v[sample(idx)]
  }
  out
}

main_30 <- function() {
  g <- readRDS(P("rds", "propagation_gates.rds"))
  A <- g$A; E <- g$E; k <- g$k; genes <- g$genes
  keep <- k > 0
  log_msg("Nodes: ", length(E), "; connected: ", sum(keep))

  # ---- D1: graph autocorrelation ---------------------------------------
  log_msg("=============== D1 autocorrelation ===============")
  nb <- neighbour_mean(A, E)
  obs <- cor(E[keep], nb[keep], method = "spearman")
  set.seed(SEED)
  nullv <- vapply(seq_len(N_PERM), function(i) {
    Ep <- degree_binned_perm(E, k, N_BINS)
    cor(Ep[keep], neighbour_mean(A, Ep)[keep], method = "spearman")
  }, numeric(1))
  z <- (obs - mean(nullv)) / sd(nullv)
  p_emp <- (1 + sum(nullv >= obs)) / (N_PERM + 1)
  log_msg("observed rho(E, neighbour mean) = ", round(obs, 4))
  log_msg("null mean = ", round(mean(nullv), 4),
          "; sd = ", round(sd(nullv), 4))
  log_msg("z = ", round(z, 2), "; empirical p = ", signif(p_emp, 3))

  # ---- D2: signal genes only, and top-decile discrimination ------------
  log_msg("============ D2 signal genes and AUC =============")
  thr <- quantile(E[keep], 1 - TOP_FRAC)
  top <- keep & E >= thr
  sig <- which(top)
  rho_sig <- cor(E[sig], nb[sig], method = "spearman")
  auc_top <- auc_rank(nb[keep], top[keep])
  log_msg("top decile threshold |effect| = ", signif(thr, 4),
          " on ", length(sig), " genes")
  log_msg("rho within the top decile = ", round(rho_sig, 4))
  log_msg("AUC of neighbour mean for top-decile membership = ",
          round(auc_top, 4))
  auc_deg <- auc_rank(k[keep], top[keep])
  log_msg("AUC of degree alone for the same = ", round(auc_deg, 4),
          "   <- the baseline to beat")

  # ---- D3: alternative signal definitions ------------------------------
  log_msg("=============== D3 signal choice =================")
  d3 <- tibble(signal = character(), rho_nb = numeric())
  de_f <- P("rds", "deg_rnaseq.rds")
  if (file.exists(de_f)) {
    rna <- readRDS(de_f)
    d <- rna$LGG_vs_Normal
    cand <- list()
    fc <- intersect(c("logFC_adj", "logFC_raw", "logFC"), names(d))[1]
    if (!is.na(fc)) {
      cand[["abs logFC"]] <- abs(d[[fc]])
      cand[["signed logFC"]] <- d[[fc]]
    }
    pc <- intersect(c("adj.P.Val", "padj", "P.Value", "pvalue"), names(d))[1]
    if (!is.na(pc)) cand[["minus log10 p"]] <- -log10(pmax(d[[pc]], 1e-300))
    for (nm in names(cand)) {
      v <- cand[[nm]][match(genes, d$gene)]
      v[!is.finite(v)] <- 0
      r <- cor(v[keep], neighbour_mean(A, v)[keep], method = "spearman")
      d3 <- bind_rows(d3, tibble(signal = nm, rho_nb = r))
      log_msg(nm, ": rho = ", round(r, 4))
    }
  } else {
    log_msg("deg_rnaseq.rds absent; D3 skipped")
  }

  # ---- D4 and D5: normalisation and degree-adjusted recovery -----------
  log_msg("========== D4 normalisation, D5 partialled ==========")
  set.seed(SEED)
  m <- sample(which(keep), round(MASK_FRAC * sum(keep)))
  Em <- E; Em[m] <- 0
  s <- Matrix::rowSums(A); s[s == 0] <- 1

  diffuse <- function(mode) {
    if (mode == "symmetric") {
      di <- 1 / sqrt(pmax(k, 1e-12))
      D <- Matrix::Diagonal(x = di)
      Aa <- Matrix::Diagonal(n = nrow(A)) +
        ALPHA * (Matrix::Diagonal(n = nrow(A)) - D %*% A %*% D)
      as.vector(Matrix::solve(as(Matrix::forceSymmetric(Aa),
                                 "symmetricMatrix"), Em))
    } else {
      # Row normalised random walk with restart, as in the propagation-gates stage.
      P_ <- Matrix::Diagonal(x = 1 / s) %*% A
      r <- 1 / (1 + ALPHA)
      x <- Em
      for (i in 1:500) {
        xn <- (1 - r) * as.vector(Matrix::crossprod(P_, x)) + r * Em
        if (sum(abs(xn - x)) < 1e-10) break
        x <- xn
      }
      x
    }
  }

  d4 <- tibble()
  for (mode in c("symmetric", "row")) {
    x <- diffuse(mode)
    rt <- cor(x[m], E[m], method = "spearman")
    rd <- cor(x[m], k[m], method = "spearman")
    # D5: rank residuals of both after removing degree.
    e1 <- residuals(lm(rank(x[m]) ~ rank(k[m])))
    e2 <- residuals(lm(rank(E[m]) ~ rank(k[m])))
    rp <- cor(e1, e2, method = "spearman")
    d4 <- bind_rows(d4, tibble(normalisation = mode, rec_true = rt,
                               rec_degree = rd, rec_partial = rp))
    log_msg(mode, ": true ", round(rt, 4), "; degree ", round(rd, 4),
            "; partialled ", round(rp, 4))
  }

  # ---- D6: control and calibration curve -------------------------------
  # A single pass/fail control is uninformative. Instead a family of signals
  # is built with graded graph smoothness, and recovery is read against
  # their autocorrelation. The observed autocorrelation from D1 is then
  # located on that curve, which says what recovery the real signal should
  # produce if propagation is working as intended.
  log_msg("======== D6 control and calibration curve ========")
  di <- 1 / sqrt(pmax(k, 1e-12))
  Dm <- Matrix::Diagonal(x = di)
  Aa <- Matrix::Diagonal(n = nrow(A)) +
    ALPHA * (Matrix::Diagonal(n = nrow(A)) - Dm %*% A %*% Dm)
  Aa <- as(Matrix::forceSymmetric(Aa), "symmetricMatrix")
  ch <- Matrix::Cholesky(Aa, LDL = FALSE, perm = TRUE)

  set.seed(SEED + 3)
  z0 <- rnorm(length(E)); z0[!keep] <- 0
  base_smooth <- as.vector(scale(neighbour_mean(A, z0)))
  noise <- rnorm(length(E))
  cal <- tibble()
  for (w in c(1, 0.8, 0.6, 0.4, 0.25, 0.15, 0.08)) {
    v <- w * base_smooth + (1 - w) * noise
    v[!keep] <- 0
    ac <- cor(v[keep], neighbour_mean(A, v)[keep], method = "spearman")
    vm <- v; vm[m] <- 0
    x <- as.vector(Matrix::solve(ch, vm, system = "A"))
    rc <- cor(x[m], v[m], method = "spearman")
    rd <- cor(x[m], k[m], method = "spearman")
    e1 <- residuals(lm(rank(x[m]) ~ rank(k[m])))
    e2 <- residuals(lm(rank(v[m]) ~ rank(k[m])))
    rp <- cor(e1, e2, method = "spearman")
    cal <- bind_rows(cal, tibble(w = w, autocorr = ac, recovery = rc,
                                 rec_degree = rd, rec_partial = rp))
    log_msg("w = ", w, ": autocorr ", round(ac, 4),
            "; recovery ", round(rc, 4), "; partialled ", round(rp, 4))
  }
  write_csv(cal, P("tables", "propagation_audit_calibration.csv"))

  ctrl <- max(cal$recovery)
  log_msg("Maximum control recovery = ", round(ctrl, 4))
  log_msg(ifelse(ctrl > 0.5, "Test is sensitive.",
                 "TEST IS BROKEN. D1 to D5 cannot be interpreted."))

  # Locate the observed autocorrelation on the curve.
  ok <- order(cal$autocorr)
  exp_rec <- approx(cal$autocorr[ok], cal$recovery[ok], xout = obs,
                    rule = 2)$y
  exp_par <- approx(cal$autocorr[ok], cal$rec_partial[ok], xout = obs,
                    rule = 2)$y
  log_msg("At the observed autocorrelation of ", round(obs, 4), ":")
  log_msg("  expected recovery ", round(exp_rec, 4),
          " versus observed ", round(d4$rec_true[1], 4))
  log_msg("  expected partialled ", round(exp_par, 4),
          " versus observed ", round(d4$rec_partial[1], 4))
  log_msg(ifelse(d4$rec_partial[1] >= 0.5 * exp_par,
                 "Real signal behaves like a graph-smooth signal of the ",
                 "Real signal underperforms a matched smooth signal. "),
          "same autocorrelation.")

  out <- tibble(
    diagnostic = c("D1 autocorrelation", "D1 null mean", "D1 z",
                   "D2 rho top decile", "D2 AUC neighbour mean",
                   "D2 AUC degree", "D6 max control recovery",
                   "D6 expected recovery at observed autocorr",
                   "D6 expected partialled at observed autocorr"),
    value = c(obs, mean(nullv), z, rho_sig, auc_top, auc_deg, ctrl,
              exp_rec, exp_par))
  write_csv(out, P("tables", "propagation_audit.csv"))
  write_csv(d4, P("tables", "propagation_audit_norm.csv"))
  if (nrow(d3)) write_csv(d3, P("tables", "propagation_audit_signal.csv"))
  log_msg("================== SUMMARY ==================")
  print(as.data.frame(out %>% mutate(value = signif(value, 4))))
  print(as.data.frame(d4 %>% mutate(across(where(is.numeric),
                                           ~ signif(.x, 4)))))
  saveRDS(list(summary = out, norm = d4, signal = d3,
               calibration = cal, null = nullv),
          P("rds", "propagation_audit.rds"))
  invisible(out)
}
