# Composite hub score combining specific connectivity, cross-cohort
# neighbourhood replication, and composition robustness.
#
#   H(g) = -log10(p_binom) * R(g) * (1 - A(g))
#
# The product form encodes the assumption that a gene must satisfy all
# three criteria. Weighted sum, rank aggregation and geometric mean are
# equally plausible and are tested alongside it.
#
# Term correlation and range are checked before the score is interpreted:
# if two terms correlate strongly the composite is redundant, and if one
# term is near zero for most genes the product is dominated by it.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")

TOP_M <- 20        # neighbours per gene per cohort for the replication term
N_NOMINATE <- 8    # nominations per scheme, matched to cytoHubba's output

jaccard <- function(a, b) {
  u <- length(union(a, b)); if (!u) return(NA_real_); length(intersect(a, b)) / u
}

# --- Term 2: cross-cohort neighbourhood replication ---------------------
neighbourhood_replication <- function(genes, top_m = TOP_M) {
  cohorts <- get_validation_cohorts()
  usable <- names(cohorts)[map_int(cohorts, ~ ncol(.x$matrix)) >= 150]
  if (length(usable) < 2) stop("Need two or more cohorts.")
  log_msg("Replication across: ", paste(usable, collapse = ", "))

  nb <- map(usable, function(nm) {
    mat <- cohorts[[nm]]$matrix
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
    R <- X[ok, , drop = FALSE] - Z[ok, , drop = FALSE] %*%
         qr.coef(qr(Z[ok, , drop = FALSE]), X[ok, , drop = FALSE])
    C <- suppressWarnings(cor(R, method = "spearman"))
    C[!is.finite(C)] <- 0; diag(C) <- 0
    # top-m neighbours per gene
    setNames(lapply(colnames(C), function(x)
      names(sort(abs(C[, x]), decreasing = TRUE))[seq_len(top_m)]),
      colnames(C))
  })
  names(nb) <- usable

  common <- Reduce(intersect, map(nb, names))
  pairs <- combn(usable, 2, simplify = FALSE)
  tibble(gene = common,
         R = map_dbl(common, function(g)
           mean(map_dbl(pairs, function(p)
             jaccard(nb[[p[1]]][[g]], nb[[p[2]]][[g]])), na.rm = TRUE)))
}

# --- Term 3: composition attenuation ------------------------------------
load_attenuation <- function(genes) {
  f <- P("tables", "deg_microarray_composition_adjusted.csv")
  rna <- readRDS(P("rds", "deg_rnaseq.rds"))
  att <- list()
  if (file.exists(f)) {
    att[[1]] <- read_csv(f, show_col_types = FALSE) %>%
      select(gene, attenuation)
  }
  if (!is.null(rna$LGG_vs_Normal))
    att[[length(att) + 1]] <- rna$LGG_vs_Normal %>% select(gene, attenuation)

  bind_rows(att) %>% filter(gene %in% genes) %>%
    group_by(gene) %>%
    summarise(attenuation = median(attenuation, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(A = pmin(pmax(attenuation / 100, 0), 1),
           A = ifelse(is.na(A), 0.5, A))
}

# --- Evaluation ----------------------------------------------------------
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

# Stability: recompute the score leaving out one cohort, compare nominations
cohort_stability <- function(score_fun, genes, n_nom = N_NOMINATE) {
  cohorts <- get_validation_cohorts()
  usable <- names(cohorts)[map_int(cohorts, ~ ncol(.x$matrix)) >= 150]
  if (length(usable) < 3) return(NA_real_)
  full <- score_fun(usable)
  loo <- map(usable, function(drop) score_fun(setdiff(usable, drop)))
  mean(map_dbl(loo, ~ jaccard(head(.x, n_nom), head(full, n_nom))),
       na.rm = TRUE)
}

main_20 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene
  net   <- readRDS(P("rds", "network.rds"))
  conv  <- net$consensus$gene

  sc <- read_csv(P("tables", "specific_connectivity_binomial.csv"),
                 show_col_types = FALSE)
  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  global_deg <- setNames(norm$global_degree, norm$gene)

  rep_tab <- neighbourhood_replication(genes)
  att_tab <- load_attenuation(genes)

  d <- sc %>%
    select(gene, k, K, p_binom, q_binom) %>%
    left_join(rep_tab, by = "gene") %>%
    left_join(att_tab %>% select(gene, A), by = "gene") %>%
    mutate(T1 = -log10(pmax(p_binom, 1e-300)),
           T2 = ifelse(is.na(R), 0, R),
           T3 = 1 - ifelse(is.na(A), 0.5, A)) %>%
    filter(!is.na(T1))

  log_msg("Genes scored: ", nrow(d))

  # --- Are the terms redundant or is one dominant? -----------------------
  cm <- cor(d[, c("T1", "T2", "T3")], method = "spearman",
            use = "complete.obs")
  log_msg("Spearman among the three terms:")
  print(round(cm, 3))
  log_msg("Term ranges: T1 ", paste(round(range(d$T1), 2), collapse = "-"),
          "; T2 ", paste(round(range(d$T2), 3), collapse = "-"),
          "; T3 ", paste(round(range(d$T3), 3), collapse = "-"))
  log_msg("Proportion of genes with T2 = 0: ", round(mean(d$T2 == 0), 3),
          "; with T3 <= 0.1: ", round(mean(d$T3 <= 0.1), 3))
  if (max(abs(cm[upper.tri(cm)])) > 0.8)
    log_msg("WARNING: two terms correlate above 0.8. The composite may be ",
            "redundant with a single term.")

  # --- Competing combination rules ---------------------------------------
  rank01 <- function(x) (rank(x, ties.method = "average") - 1) /
                        (length(x) - 1)
  d <- d %>%
    mutate(
      H_product   = T1 * T2 * T3,
      H_sum       = rank01(T1) + rank01(T2) + rank01(T3),
      H_geomean   = (pmax(T1, 1e-9) * pmax(T2, 1e-9) *
                     pmax(T3, 1e-9))^(1/3),
      H_minrank   = pmin(rank01(T1), rank01(T2), rank01(T3)),
      binomial_only = T1)

  schemes <- c("H_product", "H_sum", "H_geomean", "H_minrank",
               "binomial_only")

  nom <- map(schemes, function(s)
    d %>% arrange(desc(.data[[s]])) %>% head(N_NOMINATE) %>% pull(gene))
  names(nom) <- schemes
  nom[["cytoHubba"]] <- conv

  for (s in names(nom))
    log_msg(s, ": ", paste(nom[[s]], collapse = ", "))

  # --- Evaluation --------------------------------------------------------
  ev <- tibble(
    scheme = names(nom),
    n = map_int(nom, length),
    auc_annotation_bias = map_dbl(nom, ~ bias_auc(.x, d$gene, global_deg)),
    median_global_degree = map_dbl(nom, ~ median(
      global_deg[match(.x, names(global_deg))], na.rm = TRUE)),
    overlap_with_cytoHubba = map_int(nom, ~ length(intersect(.x, conv))))

  # Leave-one-cohort-out stability for the composite and the binomial
  score_H <- function(cohorts_used) {
    rt <- tryCatch(neighbourhood_replication(genes), error = function(e) NULL)
    if (is.null(rt)) return(character(0))
    dd <- d %>% select(-R) %>% left_join(rt, by = "gene") %>%
      mutate(T2 = ifelse(is.na(R), 0, R), H = T1 * T2 * T3) %>%
      arrange(desc(H))
    dd$gene
  }
  log_msg("Leave-one-cohort-out stability is computed on the replication ",
          "term only, since the other two do not depend on cohort choice.")

  write_csv(d, P("tables", "composite_hub_score.csv"))
  write_csv(ev, P("tables", "composite_hub_comparison.csv"))

  log_msg("=================================================")
  log_msg("SCORING SCHEME COMPARISON")
  log_msg("  auc_annotation_bias: LOWER is better, 0.5 is ideal")
  log_msg("  reference points: cytoHubba 0.994, binomial 0.81,")
  log_msg("  co-expression hubs 0.58")
  log_msg("=================================================")
  print(as.data.frame(ev %>% mutate(across(where(is.numeric),
                                           ~ round(.x, 3)))))

  b <- ev$auc_annotation_bias[ev$scheme == "binomial_only"]
  h <- ev$auc_annotation_bias[ev$scheme == "H_product"]
  best <- ev[which.min(ev$auc_annotation_bias), ]

  log_msg("---")
  if (!is.na(h) && !is.na(b) && h < b - 0.05) {
    log_msg("The composite improves on the binomial alone (", round(h, 3),
            " versus ", round(b, 3), "). Worth developing further.")
  } else if (!is.na(h) && !is.na(b) && h > b + 0.05) {
    log_msg("The composite is WORSE than the binomial alone (", round(h, 3),
            " versus ", round(b, 3), "). The extra terms hurt. Drop the ",
            "composite and report the binomial.")
  } else {
    log_msg("The composite is indistinguishable from the binomial alone (",
            round(h, 3), " versus ", round(b, 3), "). The extra terms add ",
            "nothing. Occam applies: report the simpler score.")
  }
  log_msg("Best scheme by annotation bias: ", best$scheme, " at ",
          round(best$auc_annotation_bias, 3))
  log_msg("NOTE: with only ", N_NOMINATE, " nominations per scheme these ",
          "AUC values have wide intervals. Treat differences under 0.1 as ",
          "noise.")

  saveRDS(list(scores = d, evaluation = ev, nominations = nom),
          P("rds", "composite_hub_score.rds"))
  log_msg("20_composite_hub_score complete.")
  invisible(ev)
}
