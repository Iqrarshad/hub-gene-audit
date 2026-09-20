# Reproducibility at the edge level rather than the hub level.
#
# Hub lists contain six to eight genes and are a lossy summary of a network
# with over ten thousand possible gene pairs. Jaccard on eight-item sets is
# unforgiving: a single swap moves it by 0.14.
#
# Correlation matrices are compared directly across cohorts on
# composition-residualised expression, with overlap reported across a range
# of edge budgets.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")

# Residualised correlation matrix, same construction as the co-expression network stage (17)
residual_cor <- function(mat, genes, label) {
  g <- intersect(genes, rownames(mat))
  lm_ <- log2(mat + 1)
  neuro  <- composition_score(mat, NEURONAL_PANEL, exclude = g, label = "")
  glial  <- composition_score(mat, GLIAL_PANEL,  exclude = g, label = "")
  immune <- composition_score(mat, IMMUNE_PANEL, exclude = g, label = "")
  Z <- cbind(1, neuro, glial, immune)
  keep <- c(TRUE, apply(Z[, -1, drop = FALSE], 2,
                        function(v) all(is.finite(v)) && sd(v) > 0))
  Z <- Z[, keep, drop = FALSE]

  X <- t(lm_[g, , drop = FALSE])
  ok <- complete.cases(X) & complete.cases(Z)
  X <- X[ok, , drop = FALSE]; Z <- Z[ok, , drop = FALSE]
  R <- X - Z %*% qr.coef(qr(Z), X)

  C <- suppressWarnings(cor(R, method = "spearman"))
  C[!is.finite(C)] <- 0
  diag(C) <- 0
  log_msg(label, ": ", length(g), " genes, ", nrow(R), " samples")
  C
}

edge_vector <- function(C, genes) {
  C <- C[genes, genes, drop = FALSE]
  ut <- upper.tri(C)
  tibble(pair = outer(genes, genes, paste, sep = "|")[ut],
         rho = C[ut])
}

main_19 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene

  cohorts <- get_validation_cohorts()
  usable <- names(cohorts)[map_int(cohorts, ~ ncol(.x$matrix)) >= 150]
  log_msg("Cohorts: ", paste(usable, collapse = ", "))
  if (length(usable) < 2) stop("Need at least two cohorts.")

  # Common gene set across all cohorts, so vectors are comparable
  common <- reduce(map(cohorts[usable], ~ intersect(genes,
                                                    rownames(.x$matrix))),
                   intersect)
  common <- sort(common)
  n_pairs <- choose(length(common), 2)
  log_msg("Common genes: ", length(common), "; gene pairs: ", n_pairs)

  mats <- map(usable, function(nm)
    residual_cor(cohorts[[nm]]$matrix, common, nm))
  names(mats) <- usable

  vecs <- map(mats, ~ edge_vector(.x, common))

  # --- 1. Full correlation vector agreement ------------------------------
  pairs <- combn(usable, 2, simplify = FALSE)
  res <- map_dfr(pairs, function(p) {
    a <- vecs[[p[1]]]; b <- vecs[[p[2]]]
    stopifnot(identical(a$pair, b$pair))
    rho_all <- cor(a$rho, b$rho, method = "spearman")
    rho_pear <- cor(a$rho, b$rho)

    # top-k overlap at the STRING edge budget and beyond
    topk <- function(k) {
      ta <- a$pair[order(-abs(a$rho))][seq_len(k)]
      tb <- b$pair[order(-abs(b$rho))][seq_len(k)]
      length(intersect(ta, tb)) / k
    }
    tibble(cohort_a = p[1], cohort_b = p[2],
           spearman_all_pairs = rho_all,
           pearson_all_pairs = rho_pear,
           top139_overlap = topk(139),
           top500_overlap = topk(500),
           top1000_overlap = topk(1000),
           expected_top139 = round(139 / n_pairs, 4))
  })
  write_csv(res, P("tables", "edge_reproducibility.csv"))

  log_msg("=================================================")
  log_msg("EDGE-LEVEL REPRODUCIBILITY")
  log_msg("  Spearman across all ", n_pairs, " gene pairs, and overlap of")
  log_msg("  the strongest edges. Expected top-139 overlap by chance is ",
          round(139 / n_pairs, 4))
  log_msg("=================================================")
  print(as.data.frame(res %>% mutate(across(where(is.numeric),
                                            ~ round(.x, 3)))))

  # --- 2. Overlap as a function of threshold -----------------------------
  ks <- c(50, 139, 300, 500, 1000, 2000)
  curve <- map_dfr(pairs, function(p) {
    a <- vecs[[p[1]]]; b <- vecs[[p[2]]]
    map_dfr(ks, function(k) {
      ta <- a$pair[order(-abs(a$rho))][seq_len(k)]
      tb <- b$pair[order(-abs(b$rho))][seq_len(k)]
      tibble(cohort_a = p[1], cohort_b = p[2], k = k,
             overlap = length(intersect(ta, tb)) / k,
             fold_over_chance = (length(intersect(ta, tb)) / k) /
                                (k / n_pairs))
    })
  })
  write_csv(curve, P("tables", "edge_overlap_curve.csv"))
  log_msg("Overlap by edge budget:")
  print(as.data.frame(curve %>% mutate(across(where(is.numeric),
                                              ~ round(.x, 3)))))

  # --- 3. Edge agreement versus hub agreement ----------------------------
  hub_f <- P("tables", "grade_matched_hub_agreement.csv")
  hub_between <- NA_real_
  if (file.exists(hub_f)) {
    ha <- read_csv(hub_f, show_col_types = FALSE)
    bet <- ha %>% filter(grepl("\\[all\\]", set_a), grepl("\\[all\\]", set_b))
    hub_between <- mean(bet$hub_jaccard, na.rm = TRUE)
  }

  cmp <- tibble(
    level = c("hub list (6-8 genes)", "top 139 edges",
              "all gene pairs (rank correlation)"),
    agreement = c(hub_between,
                  mean(res$top139_overlap, na.rm = TRUE),
                  mean(res$spearman_all_pairs, na.rm = TRUE)))
  write_csv(cmp, P("tables", "edge_vs_hub_agreement.csv"))
  log_msg("Agreement by level of description:")
  print(as.data.frame(cmp %>% mutate(agreement = round(agreement, 3))))

  # --- Verdict -----------------------------------------------------------
  edge_rho <- mean(res$spearman_all_pairs, na.rm = TRUE)
  top139 <- mean(res$top139_overlap, na.rm = TRUE)

  log_msg("---")
  if (edge_rho >= 0.5 && !is.na(hub_between) && hub_between < 0.3) {
    log_msg("VERDICT: the network reproduces across cohorts at the edge ",
            "level (rank correlation ", round(edge_rho, 3),
            ") while hub lists do not (", round(hub_between, 3), "). ",
            "Thresholding a continuous correlation structure into a short ",
            "hub list destroys reproducible information. Recommend ",
            "reporting edge-level or module-level structure.")
  } else if (edge_rho < 0.3) {
    log_msg("VERDICT: edge-level agreement is also low (", round(edge_rho, 3),
            "). Co-expression structure does not transfer between these ",
            "cohorts. No model can recover what is not stable, and this ",
            "should be reported as a limit of the data rather than of the ",
            "method.")
  } else {
    log_msg("VERDICT: edge agreement is moderate (", round(edge_rho, 3),
            "). Partial transfer. Report the edge and hub figures side by ",
            "side and draw no strong conclusion in either direction.")
  }

  if (!is.na(top139)) {
    log_msg("Top-139 edge overlap ", round(top139, 3), " is ",
            round(top139 / (139 / n_pairs), 1),
            " times chance, but note that a high fold-change over a tiny ",
            "baseline can still be a small absolute overlap.")
  }

  saveRDS(list(res = res, curve = curve, comparison = cmp, mats = mats),
          P("rds", "edge_reproducibility.rds"))
  log_msg("19_edge_reproducibility complete.")
  invisible(res)
}
