# Constrained diffusion: is orthogonalising the operator better than
# adjusting the output?
#
# An earlier run on this network with alpha = 1 established:
#   graph autocorrelation of the effect vector   0.2166
#   expected recovery at that autocorrelation    0.2587  (calibration curve)
#   raw held-out recovery                        0.0651
#   recovery after degree is partialled out      0.1953
#   AUC of degree for strong differential expression  0.4167, below chance
#
# Differential expression is anti-correlated with degree, symmetric
# diffusion injects a positive degree component, and the two oppose. This
# stage constrains the operator instead of correcting its output:
#
#   x* = A^-1 E - A^-1 C (Lam^-1 + C' A^-1 C)^-1 C' A^-1 E,  A = I + alpha L
#
# with Lam = lambda I. As lambda tends to infinity this is the exact
# projection onto {x : C'x = 0}. Woodbury gives every lambda from one
# factorisation. Constraint sets: degree, citation, and both.
#
# Criteria, fixed before running and taken from that calibration
# curve rather than chosen by hand:
#   S1  raw rec_true >= 0.20        success, 77 percent of the 0.2587
#                                   expected at this autocorrelation
#   S2  raw rec_true < 0.12         failure
#       between the two             inconclusive
#   S3  the operator must match or beat post-hoc partialling of degree
#       from the unconstrained output (0.1953). If it does not, adjusting
#       the output is simpler and the operator is unnecessary.
#
# S3 is the criterion that decides whether there is a method here.
#
# Input: propagation_gates.rds from the propagation-gates stage (27).

suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(readr); library(tibble)
})

ALPHAS    <- c(0.5, 1, 2)
LAMBDAS   <- c(10 ^ seq(-2, 6, length.out = 25), Inf)
MASK_FRAC <- 0.20
SEED      <- 11
SUCCESS   <- 0.20
FAILURE   <- 0.12

build_operator <- function(A, alpha) {
  k <- Matrix::rowSums(A)
  di <- 1 / sqrt(pmax(k, 1e-12))
  D <- Matrix::Diagonal(x = di)
  I <- Matrix::Diagonal(n = nrow(A))
  as(Matrix::forceSymmetric(I + alpha * (I - D %*% A %*% D)),
     "symmetricMatrix")
}

unit_centre <- function(v) {
  v <- v - mean(v)
  v / sqrt(sum(v ^ 2))
}

# Woodbury solution for a rank-k constraint at one lambda.
constrained <- function(x0, W, CtX0, CtW, lambda) {
  kk <- ncol(W)
  Linv <- if (is.infinite(lambda)) matrix(0, kk, kk)
          else diag(1 / lambda, kk, kk)
  x0 - as.vector(W %*% solve(Linv + CtW, CtX0))
}

main_31 <- function() {
  g <- readRDS(P("rds", "propagation_gates.rds"))
  A <- g$A; E <- g$E; C <- g$C; k <- g$k
  keep <- k > 0
  set.seed(SEED)
  m <- sample(which(keep), round(MASK_FRAC * sum(keep)))
  Em <- E; Em[m] <- 0
  log_msg("Masked ", length(m), " of ", sum(keep), " connected nodes")

  c_deg <- unit_centre(log1p(k))
  c_cit <- unit_centre(C)
  log_msg("rho(log degree, log citations) = ",
          round(cor(log1p(k), C, method = "spearman"), 4))

  sets <- list(degree = cbind(c_deg),
               citation = cbind(c_cit),
               both = cbind(c_deg, c_cit))

  crit <- tibble(
    criterion = c("S1 success", "S2 failure", "S3 beats post-hoc"),
    threshold = c(paste(">=", SUCCESS), paste("<", FAILURE),
                  ">= post-hoc partialled value"))
  write_csv(crit, P("tables", "constrained_criteria.csv"))
  log_msg("Criteria written before the sweep")

  rows <- list(); base <- list()
  for (alpha in ALPHAS) {
    log_msg("=== alpha = ", alpha, " ===")
    Aa <- build_operator(A, alpha)
    ch <- Matrix::Cholesky(Aa, LDL = FALSE, perm = TRUE)
    x0 <- as.vector(Matrix::solve(ch, Em, system = "A"))

    # Unconstrained reference and the post-hoc alternative.
    r_raw <- cor(x0[m], E[m], method = "spearman")
    e1 <- residuals(lm(rank(x0[m]) ~ rank(k[m])))
    e2 <- residuals(lm(rank(E[m]) ~ rank(k[m])))
    r_post <- cor(e1, e2, method = "spearman")
    log_msg("unconstrained raw ", round(r_raw, 4),
            "; post-hoc partialled ", round(r_post, 4))
    base[[length(base) + 1]] <- tibble(alpha = alpha, raw = r_raw,
                                       posthoc = r_post)

    for (nm in names(sets)) {
      Cm <- sets[[nm]]
      W <- as.matrix(Matrix::solve(ch, Cm, system = "A"))
      CtW <- crossprod(Cm, W)
      CtX0 <- crossprod(Cm, x0)
      for (lam in LAMBDAS) {
        x <- constrained(x0, W, CtX0, CtW, lam)
        rows[[length(rows) + 1]] <- tibble(
          alpha = alpha, constraint = nm, lambda = lam,
          rec_true = cor(x[m], E[m], method = "spearman"),
          rec_degree = cor(x[m], k[m], method = "spearman"),
          rec_citation = cor(x[m], C[m], method = "spearman"),
          cx_degree = sum(c_deg * x), cx_citation = sum(c_cit * x))
      }
    }
  }
  res <- bind_rows(rows); bl <- bind_rows(base)
  write_csv(res, P("tables", "constrained_recovery.csv"))
  write_csv(bl, P("tables", "constrained_baselines.csv"))

  log_msg("================= BASELINES =================")
  print(as.data.frame(bl %>% mutate(across(where(is.numeric),
                                           ~ signif(.x, 4)))))

  log_msg("=========== BEST PER CONSTRAINT ============")
  best <- res %>% group_by(alpha, constraint) %>%
    slice_max(rec_true, n = 1, with_ties = FALSE) %>% ungroup() %>%
    arrange(desc(rec_true))
  print(as.data.frame(best %>%
    transmute(alpha, constraint, lambda = signif(lambda, 3),
              rec_true = round(rec_true, 4),
              rec_degree = round(rec_degree, 4),
              rec_citation = round(rec_citation, 4))))

  top <- best$rec_true[1]
  post_at <- bl$posthoc[bl$alpha == best$alpha[1]]
  edge <- is.infinite(best$lambda[1]) ||
    best$lambda[1] >= max(LAMBDAS[is.finite(LAMBDAS)])
  if (edge)
    log_msg("NOTE: best lambda is at or beyond the sweep edge, so the ",
            "projection limit is optimal and lambda is not a real knob.")

  out <- crit %>% mutate(
    value = c(signif(top, 4), signif(top, 4), signif(post_at, 4)),
    pass = c(top >= SUCCESS, !(top < FAILURE), top >= post_at))
  write_csv(out, P("tables", "constrained_criteria.csv"))
  log_msg("================== SUMMARY =================")
  print(as.data.frame(out))

  if (top >= SUCCESS && top >= post_at) {
    log_msg("Constrained operator reaches the target and is not beaten ",
            "by post-hoc adjustment.")
  } else if (top < post_at) {
    log_msg("Post-hoc adjustment of the unconstrained output matches or ",
            "beats the operator. Report the simpler procedure.")
  } else if (top < FAILURE) {
    log_msg("Below the failure threshold. Constraining the operator does ",
            "not recover the signal.")
  } else {
    log_msg("Inconclusive. Move to per-degree-bin constraints.")
  }

  saveRDS(list(recovery = res, baselines = bl, criteria = out),
          P("rds", "constrained_recovery.rds"))
  invisible(res)
}
