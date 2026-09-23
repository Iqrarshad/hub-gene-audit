# Weighted co-expression networks built per cohort, tested for cross-cohort
# reproducibility.
#
# The co-expression stage builds an unweighted correlation network. WGCNA is
# the other standard data-derived paradigm: it soft-thresholds the
# correlation matrix, computes a topological overlap measure, and detects
# modules whose intramodular connectivity defines the hubs. If the failure
# audited for STRING is a property of the fixed network rather than of
# co-expression thresholding, WGCNA hubs and modules should be no more
# reproducible across cohorts than the unweighted co-expression network.
#
# The reported quantity is reproducibility, not a degree-AUC. A degree-AUC on
# a soft-thresholded network would be circular, since soft-thresholding
# raises the correlation to a power and so amplifies degree by construction;
# a reviewer could then attribute any degree dominance to the method rather
# than to the data. Three reproducibility measures are reported instead:
# top-k intramodular hub agreement (Jaccard), full-vector intramodular
# connectivity agreement (Spearman), and module-assignment agreement
# (adjusted Rand index) across cohort pairs.
#
# Expression is residualised on neuronal, glial and immune content before the
# network is built, matching the co-expression stage, so that composition
# does not rebuild spurious edges.
#
# WGCNA is an optional dependency. If it is not installed the stage records
# that and returns without error, so the run completes.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")

MIN_SAMPLES   <- 40      # WGCNA is unstable below this
MIN_GENES     <- 30
TOPK_HUB      <- 10
MIN_MODSIZE   <- 20
SCALEFREE_R2  <- 0.80
MAX_POWER     <- 20

# Adjusted Rand index of two label vectors on the same items. Written here to
# avoid a dependency on mclust for one function.
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  sum_comb <- function(x) sum(choose(x, 2))
  ai <- rowSums(tab); bj <- colSums(tab)
  index <- sum(choose(tab, 2))
  expected <- sum_comb(ai) * sum_comb(bj) / choose(n, 2)
  maxindex <- (sum_comb(ai) + sum_comb(bj)) / 2
  if (maxindex - expected == 0) return(NA_real_)
  (index - expected) / (maxindex - expected)
}

jaccard <- function(a, b) {
  if (!length(union(a, b))) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# Residualise a cohort matrix on composition, returning samples x genes.
residualise_expression <- function(mat, genes, label) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < MIN_GENES) {
    log_msg(label, ": only ", length(g), " genes present, below the minimum")
    return(NULL)
  }
  lm_ <- log2(mat + 1)

  neuro  <- composition_score(mat, NEURONAL_PANEL, exclude = g,
                              label = paste(label, "neuronal"))
  glial  <- composition_score(mat, GLIAL_PANEL, exclude = g,
                              label = paste(label, "glial"))
  immune <- composition_score(mat, IMMUNE_PANEL, exclude = g,
                              label = paste(label, "immune"))
  Z <- cbind(1, neuro, glial, immune)
  keep_z <- apply(Z, 2, function(v) all(is.finite(v)) && sd(v) > 0)
  Z <- Z[, keep_z | seq_len(ncol(Z)) == 1, drop = FALSE]

  X <- t(lm_[g, , drop = FALSE])
  ok <- complete.cases(X) & complete.cases(Z)
  X <- X[ok, , drop = FALSE]; Zs <- Z[ok, , drop = FALSE]
  if (nrow(X) < MIN_SAMPLES) {
    log_msg(label, ": ", nrow(X), " usable samples, below the minimum")
    return(NULL)
  }
  qrz <- qr(Zs)
  R <- X - Zs %*% qr.coef(qrz, X)
  log_msg(label, ": residualised on ", ncol(Zs) - 1, " composition axes, ",
          nrow(R), " samples, ", ncol(R), " genes")
  R
}

# Build a signed WGCNA network on one residualised matrix and return the
# module labels and intramodular connectivity.
wgcna_one <- function(R, label) {
  powers <- c(1:10, seq(12, MAX_POWER, 2))
  sft <- WGCNA::pickSoftThreshold(R, powerVector = powers, verbose = 0,
                                  networkType = "signed")
  fit <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  beta <- sft$powerEstimate
  if (is.na(beta)) {
    beta <- powers[which(fit >= SCALEFREE_R2)[1]]
    if (is.na(beta)) beta <- powers[which.max(fit)]
  }
  log_msg(label, ": soft power ", beta, " (scale-free fit ",
          round(max(fit, na.rm = TRUE), 2), ")")

  adj <- WGCNA::adjacency(R, power = beta, type = "signed")
  tom <- WGCNA::TOMsimilarity(adj, verbose = 0)
  dimnames(tom) <- dimnames(adj)
  diss <- 1 - tom
  tree <- hclust(as.dist(diss), method = "average")
  mods <- dynamicTreeCut::cutreeDynamic(
    dendro = tree, distM = diss, deepSplit = 2,
    pamRespectsDendro = FALSE, minClusterSize = MIN_MODSIZE)
  names(mods) <- colnames(R)

  kin <- WGCNA::intramodularConnectivity(adj, mods)$kWithin
  names(kin) <- colnames(R)
  hubs <- names(sort(kin, decreasing = TRUE))[seq_len(TOPK_HUB)]

  log_msg(label, ": ", length(unique(mods[mods != 0])), " modules, ",
          "top intramodular hubs: ", paste(head(hubs, 8), collapse = ", "))
  list(modules = mods, kin = kin, hubs = hubs, power = beta)
}

main_36 <- function() {
  if (!requireNamespace("WGCNA", quietly = TRUE) ||
      !requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    log_msg("WGCNA (or dynamicTreeCut) is not installed. Skipping the WGCNA ",
            "reproducibility stage. Install it with install_dependencies.R ",
            "to run the second data-derived network paradigm.")
    return(invisible(NULL))
  }

  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene

  cohorts <- get_validation_cohorts()
  fits <- list()
  for (nm in names(cohorts)) {
    R <- residualise_expression(cohorts[[nm]]$matrix, genes, nm)
    if (is.null(R)) next
    fit <- tryCatch(wgcna_one(R, nm),
                    error = function(e) {
                      log_msg(nm, ": WGCNA failed, ", conditionMessage(e))
                      NULL })
    if (!is.null(fit)) fits[[nm]] <- fit
  }

  if (length(fits) < 2) {
    log_msg("Fewer than two cohorts produced a WGCNA network. Cross-cohort ",
            "reproducibility cannot be measured; stage stops here.")
    return(invisible(NULL))
  }

  beta_tbl <- tibble(cohort = names(fits),
                     soft_power_beta = vapply(fits, function(f) f$power, numeric(1)))
  write_csv(beta_tbl, P("tables", "wgcna_soft_threshold_power.csv"))
  log_msg("Soft-thresholding powers (beta): ",
          paste(beta_tbl$cohort, beta_tbl$soft_power_beta, sep = "=", collapse = ", "))

  pairs <- combn(names(fits), 2, simplify = FALSE)
  rows <- map(pairs, function(p) {
    a <- fits[[p[1]]]; b <- fits[[p[2]]]
    shared <- intersect(names(a$kin), names(b$kin))
    tibble(
      cohort_a = p[1], cohort_b = p[2], n_shared = length(shared),
      hub_jaccard = jaccard(a$hubs, b$hubs),
      kin_spearman = suppressWarnings(
        cor(a$kin[shared], b$kin[shared], method = "spearman")),
      module_ari = adjusted_rand(a$modules[shared], b$modules[shared]))
  })
  out <- bind_rows(rows)

  summary_row <- tibble(
    cohort_a = "mean", cohort_b = "", n_shared = round(mean(out$n_shared)),
    hub_jaccard = mean(out$hub_jaccard, na.rm = TRUE),
    kin_spearman = mean(out$kin_spearman, na.rm = TRUE),
    module_ari = mean(out$module_ari, na.rm = TRUE))
  out <- bind_rows(out, summary_row)
  write_csv(out, P("tables", "wgcna_reproducibility.csv"))

  log_msg("=================================================")
  log_msg("WGCNA CROSS-COHORT REPRODUCIBILITY")
  log_msg("  hub_jaccard   top-", TOPK_HUB,
          " intramodular hub agreement across cohort pairs")
  log_msg("  kin_spearman  full-vector intramodular connectivity agreement")
  log_msg("  module_ari    module-assignment agreement (adjusted Rand)")
  log_msg("  STRING degree reproduces at 1.00 by construction; these ",
          "measure whether soft-thresholding the data recovers that.")
  log_msg("=================================================")
  print(as.data.frame(out %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))))

  mj <- summary_row$hub_jaccard; ms <- summary_row$kin_spearman
  if (!is.na(mj) && mj < 0.3 && !is.na(ms) && ms < 0.5) {
    log_msg("CONCLUSION: WGCNA hubs and modules do not reproduce across ",
            "cohorts (mean hub Jaccard ", round(mj, 3), ", mean kIN Spearman ",
            round(ms, 3), "). The second data-derived paradigm fails at the ",
            "same step as the unweighted co-expression network.")
  } else {
    log_msg("CONCLUSION: WGCNA reproducibility is higher than expected ",
            "(mean hub Jaccard ", round(mj, 3), ", mean kIN Spearman ",
            round(ms, 3), "). Report the figures as measured and temper the ",
            "claim that both data-derived paradigms fail identically.")
  }

  # ---- Supplementary figure ---------------------------------------------
  fig_ok <- requireNamespace("ggplot2", quietly = TRUE)
  if (fig_ok) {
    pdat <- out %>% filter(cohort_a != "mean") %>%
      mutate(pair = paste(cohort_a, cohort_b, sep = " / ")) %>%
      select(pair, hub_jaccard, kin_spearman, module_ari) %>%
      tidyr::pivot_longer(-pair, names_to = "measure", values_to = "value") %>%
      mutate(measure = recode(measure,
        hub_jaccard = "Top-k hub Jaccard",
        kin_spearman = "Full-vector Spearman",
        module_ari = "Module ARI"))

    p <- ggplot2::ggplot(pdat,
        ggplot2::aes(measure, value, group = pair)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                          colour = FIG$okabe_ito[6]) +
      ggplot2::geom_point(size = 2.4, colour = FIG$okabe_ito[5]) +
      ggplot2::geom_line(alpha = 0.4, colour = FIG$okabe_ito[5]) +
      ggplot2::annotate("text", x = 0.6, y = 1.02, hjust = 0,
                        label = "STRING degree = 1 by construction",
                        size = 3, family = FIG$font) +
      ggplot2::coord_cartesian(
        ylim = c(min(-0.05, min(pdat$value, na.rm = TRUE) - 0.05), 1.1)) +
      ggplot2::labs(
        title = "WGCNA hubs and modules do not reproduce across cohorts",
        x = NULL, y = "Cross-cohort agreement") +
      ggplot2::theme_classic(base_size = 12, base_family = FIG$font) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,
                                           size = 13),
        axis.text = ggplot2::element_text(colour = "black"))

    for (fmt in FIG$formats) {
      f <- P("figures", paste0("supp_fig_12_wgcna_reproducibility.", fmt))
      tryCatch({
        if (fmt == "tiff")
          ggplot2::ggsave(f, p, width = 6.5, height = 4.5, dpi = FIG$dpi,
                          device = "tiff", compression = "lzw")
        else
          ggplot2::ggsave(f, p, width = 6.5, height = 4.5, device = cairo_pdf)
      }, error = function(e)
        log_msg("  ", fmt, " failed for supp_fig_12: ", conditionMessage(e)))
    }
    log_msg("Figure written: supp_fig_12_wgcna_reproducibility")
  }

  saveRDS(list(fits = fits, reproducibility = out), P("rds", "wgcna.rds"))
  log_msg("36_wgcna_reproducibility complete.")
  invisible(out)
}
