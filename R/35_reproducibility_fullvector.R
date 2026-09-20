# Full-vector reproducibility of co-expression degree across cohorts.
#
# The co-expression comparison showed that the top-k hub lists share no genes
# between CGGA-325 and CGGA-693, while STRING degree top-k lists share
# nearly all. A top-k Jaccard of zero can arise either because the whole
# degree ranking is uncorrelated, or because the ranking is weakly
# correlated but the sharp top-k cut amplifies small differences.
#
# This distinguishes the two. It computes the Spearman correlation of the
# entire co-expression degree vector between the two cohorts, on the shared
# gene set, at several network densities. It also reports the same for
# STRING degree as the reference, which is identical across cohorts by
# construction and therefore correlates at 1.
#
# Reading:
#   full-vector rho near 0    the ranking itself does not replicate;
#                             the zero top-k overlap is real, not a cut
#                             artefact
#   full-vector rho moderate  the ranking partly replicates but hub
#                             extraction destroys it; state that the
#                             instability is in the thresholding step
#
# Standalone. Reads cohort matrices via 00c_cohorts.R and STRING degree
# via propagation_gates.rds. Writes coexpr_fullvector.csv.

suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

source_once <- function(f) if (!exists("get_validation_cohorts"))
  source(f, local = FALSE)

coexpr_degree_vec <- function(X, density) {
  lm_ <- log2(X + 1)
  C <- suppressWarnings(cor(t(lm_), method = "spearman"))
  diag(C) <- 0; C[!is.finite(C)] <- 0
  ut <- which(upper.tri(C), arr.ind = TRUE)
  vals <- abs(C[ut])
  n_edges <- floor(density * length(vals))
  cut <- sort(vals, decreasing = TRUE)[n_edges]
  A <- (abs(C) >= cut) * 1
  setNames(rowSums(A), rownames(X))
}

main_35 <- function() {
  source_once("R/00c_cohorts.R")
  coh <- get_validation_cohorts()

  use <- list()
  for (nm in names(coh)) {
    y <- coh[[nm]]$meta$grade[match(colnames(coh[[nm]]$matrix),
                                    coh[[nm]]$meta$sample)]
    if (length(unique(na.omit(y))) >= 2) {
      ok <- !is.na(y)
      X <- coh[[nm]]$matrix[, ok, drop = FALSE]
      v <- apply(X, 1, var)
      use[[nm]] <- X[v > quantile(v, 0.5) & rowMeans(X) > 1, , drop = FALSE]
    }
  }
  nms <- names(use)
  if (length(use) < 2) stop("Need two grade cohorts.")
  log_msg("Cohorts: ", paste(nms, collapse = ", "))

  ov <- intersect(rownames(use[[1]]), rownames(use[[2]]))
  log_msg("Shared genes: ", length(ov))

  g <- readRDS(P("rds", "propagation_gates.rds"))
  string_deg <- setNames(g$k, g$genes)
  sd_shared <- string_deg[intersect(names(string_deg), ov)]

  rows <- list()
  paired <- NULL
  for (density in c(0.0005, 0.002, 0.01, 0.05)) {
    d1 <- coexpr_degree_vec(use[[1]], density)[ov]
    d2 <- coexpr_degree_vec(use[[2]], density)[ov]
    rho_coex <- cor(d1, d2, method = "spearman")
    # Keep the paired vectors at the mid density for the scatter figure.
    if (density == 0.01) {
      sd1 <- string_deg[ov]; sd2 <- string_deg[ov]
      paired <- tibble(gene = ov,
                       coexpr_325 = as.numeric(d1),
                       coexpr_693 = as.numeric(d2),
                       string = as.numeric(sd1))
    }
    log_msg("density ", density,
            ": full-vector Spearman of co-expression degree = ",
            round(rho_coex, 4),
            " (n = ", length(ov), ")")
    rows[[length(rows) + 1]] <- tibble(
      density = density, measure = "coexpression_degree",
      spearman = rho_coex, n = length(ov))
  }

  # STRING degree is the same vector for both cohorts, so its cross-cohort
  # Spearman is 1 by construction. Reported to make the contrast explicit.
  rows[[length(rows) + 1]] <- tibble(
    density = NA_real_, measure = "string_degree_reference",
    spearman = 1.0, n = length(sd_shared))
  log_msg("STRING degree cross-cohort Spearman = 1 by construction ",
          "(the same global vector is used for every cohort)")

  out <- bind_rows(rows)
  write_csv(out, P("tables", "coexpr_fullvector.csv"))
  if (!is.null(paired))
    write_csv(paired, P("tables", "coexpr_paired_degree.csv"))
  log_msg("=================== SUMMARY ===================")
  print(as.data.frame(out %>% mutate(spearman = round(spearman, 4))))

  m <- mean(out$spearman[out$measure == "coexpression_degree"])
  log_msg(if (abs(m) < 0.2)
    paste0("Full-vector correlation is near zero (mean ", round(m, 3),
           "). The co-expression degree ranking does not replicate ",
           "across cohorts. The zero top-k overlap is not a cut artefact.")
    else
    paste0("Full-vector correlation is moderate (mean ", round(m, 3),
           "). The ranking partly replicates; the instability is in the ",
           "top-k hub extraction step. Report it that way."))

  saveRDS(out, P("rds", "coexpr_fullvector.rds"))
  invisible(out)
}
