# Composition-adjusted differential expression.
#
# Tumour-versus-normal comparisons in brain are confounded by cell type
# composition: normal tissue is neuron-rich, tumour tissue is not, so
# neuronal genes appear downregulated for compositional reasons.
#
# For each gene the conventional effect, the effect with a composition
# covariate, and the attenuation between them are reported. A gene is
# composition-robust only if significant in both models with the same sign.
#
# The score excludes the genes under test, and drops any marker correlating
# above 0.85 with a tested gene, so no gene is adjusted against itself.

suppressPackageStartupMessages({
  library(limma); library(dplyr); library(readr); library(tibble)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

NEURONAL_PANEL <- c("RBFOX3", "SYN1", "SYT1", "NEFL", "NEFM", "MAP2",
                    "GRIN1", "GAD1", "GAD2", "SLC17A6", "STMN2", "TUBB3",
                    "ENO2", "CAMK2A", "NRGN", "SYP", "DLG4", "SNCB")

# Also track glial and immune content: a tumour-vs-normal contrast is
# confounded by these too, though less severely than by neurons.
GLIAL_PANEL  <- c("GFAP", "AQP4", "S100B", "SLC1A2", "SLC1A3", "ALDH1L1")
IMMUNE_PANEL <- c("PTPRC", "CD68", "AIF1", "CSF1R", "ITGAM", "TMEM119")

# --- Score construction --------------------------------------------------
composition_score <- function(mat, panel, exclude = character(0),
                              label = "score") {
  present <- setdiff(intersect(panel, rownames(mat)), exclude)
  if (length(present) < 4) {
    warning(label, ": only ", length(present), " markers present.")
    if (length(present) == 0) return(rep(NA_real_, ncol(mat)))
  }
  lm_ <- if (max(mat, na.rm = TRUE) > 50) log2(mat[present, , drop = FALSE] + 1)
         else mat[present, , drop = FALSE]
  z <- t(scale(t(lm_)))
  colMeans(z, na.rm = TRUE)
}

# Drop markers too collinear with the gene being tested
safe_panel_for_gene <- function(mat, panel, gene, cutoff = 0.85) {
  present <- intersect(panel, rownames(mat))
  present <- setdiff(present, gene)
  if (!gene %in% rownames(mat) || length(present) == 0) return(present)
  g <- as.numeric(mat[gene, ])
  keep <- vapply(present, function(m) {
    r <- suppressWarnings(cor(g, as.numeric(mat[m, ]), method = "spearman",
                              use = "complete.obs"))
    is.na(r) || abs(r) < cutoff
  }, logical(1))
  present[keep]
}

# --- The core routine ----------------------------------------------------
# mat   : numeric matrix, symbols as rownames
# group : factor-able vector, must contain the two contrast levels
# adjust_for: which composition axes enter the model.
#   "neuronal" and "glial" are the defaults for tumour-vs-normal.
#   "immune" MUST be added for any tumour-vs-tumour contrast, because grade
#   comparisons differ in leukocyte infiltration rather than in neuron
#   content. Omitting it was an error in the first version of this script:
#   the immune score was computed and then left out of the design matrix,
#   so grade contrasts were adjusted for neurons and glia only.
compositional_deg <- function(mat, group, label,
                              contrast = c("Tumour", "Normal"),
                              extra_covariates = NULL,
                              adjust_for = c("neuronal", "glial")) {

  group <- factor(group, levels = contrast)
  keep_s <- !is.na(group)
  mat <- mat[, keep_s, drop = FALSE]; group <- droplevels(group[keep_s])

  tb <- table(group)
  log_msg(label, ": ", paste(names(tb), tb, sep = "=", collapse = ", "))
  if (length(tb) < 2 || any(tb < 3)) {
    warning(label, ": fewer than 3 samples in a group. Skipping.")
    return(NULL)
  }

  lmat <- if (max(mat, na.rm = TRUE) > 50) log2(mat + 1) else mat

  neuro  <- composition_score(mat, NEURONAL_PANEL, label = paste(label, "neuronal"))
  glial  <- composition_score(mat, GLIAL_PANEL,    label = paste(label, "glial"))
  immune <- composition_score(mat, IMMUNE_PANEL,   label = paste(label, "immune"))

  n_res <- length(levels(group)) + 1 + !is.null(extra_covariates)
  if (ncol(lmat) - n_res < 3) {
    warning(label, ": too few residual degrees of freedom for adjustment. ",
            "Unadjusted model only.")
    neuro <- NULL
  }

  # --- Model 1: conventional ---------------------------------------------
  d0 <- model.matrix(~ 0 + group); colnames(d0) <- levels(group)
  cm <- makeContrasts(contrasts = paste0(contrast[1], "-", contrast[2]),
                      levels = d0)
  f0 <- eBayes(contrasts.fit(lmFit(lmat, d0), cm), trend = TRUE)
  t0 <- topTable(f0, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("gene") %>%
    transmute(gene, logFC_raw = logFC, p_raw = P.Value, q_raw = adj.P.Val)

  if (is.null(neuro)) {
    out <- t0 %>% mutate(dataset = label, logFC_adj = NA_real_,
                         q_adj = NA_real_, attenuation = NA_real_,
                         composition_robust = NA)
    return(out)
  }

  # --- Model 2: composition adjusted -------------------------------------
  cov_all <- data.frame(neuronal = neuro, glial = glial, immune = immune)
  cov_df <- cov_all[, intersect(adjust_for, names(cov_all)), drop = FALSE]
  log_msg(label, ": adjusting for ", paste(names(cov_df), collapse = " + "))
  cov_df <- cov_df[, apply(cov_df, 2, function(v)
    !all(is.na(v)) && sd(v, na.rm = TRUE) > 0), drop = FALSE]
  if (!is.null(extra_covariates)) cov_df <- cbind(cov_df, extra_covariates)

  d1 <- model.matrix(~ 0 + group + ., data = cov_df)
  colnames(d1)[seq_along(levels(group))] <- levels(group)
  cm1 <- makeContrasts(contrasts = paste0(contrast[1], "-", contrast[2]),
                       levels = d1)
  f1 <- eBayes(contrasts.fit(lmFit(lmat, d1), cm1), trend = TRUE)
  t1 <- topTable(f1, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("gene") %>%
    transmute(gene, logFC_adj = logFC, p_adj = P.Value, q_adj = adj.P.Val)

  res <- full_join(t0, t1, by = "gene") %>%
    mutate(
      dataset = label,
      covariates = paste(names(cov_df), collapse = "+"),
      attenuation = 100 * (1 - abs(logFC_adj) / abs(logFC_raw)),
      sig_raw = abs(logFC_raw) >= THRESH$deg_logfc & q_raw < THRESH$deg_adjp,
      sig_adj = abs(logFC_adj) >= THRESH$deg_logfc & q_adj < THRESH$deg_adjp,
      same_sign = sign(logFC_raw) == sign(logFC_adj),
      composition_robust = sig_raw & sig_adj & same_sign,
      composition_driven = sig_raw & !sig_adj) %>%
    arrange(desc(composition_robust), q_adj)

  log_msg(label, ": ", sum(res$sig_raw, na.rm = TRUE), " DEGs unadjusted, ",
          sum(res$composition_robust, na.rm = TRUE), " composition-robust, ",
          sum(res$composition_driven, na.rm = TRUE), " composition-driven")

  res
}
