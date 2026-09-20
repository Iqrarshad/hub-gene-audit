# Main manuscript figures for the hub-gene audit (Briefings in Bioinformatics).
#
# Six figures, each geometry matched to its result, each with a bold centred
# title and bold panel labels. Colourblind-safe Okabe-Ito without yellow.
# Output: 1200 DPI TIFF (LZW) and PDF, Times New Roman, max width 174 mm.
#
# F1  hub membership is degree            AUC of degree predicting hubs
# F2  degree tracks citation              publication count vs degree
# F3  reproducibility is fixedness        paired scatter, the keystone
# F4  the recovery ceiling                calibration curve with observed
# F5  degree is anti-informative          AUC of predictors of DE
# F6  the selection ladder                reproducibility by selector
#
# Reads saved tables only; a figure whose input is absent is skipped.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tibble)
  library(tidyr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# Named colours, so removing a palette slot never shifts a figure.
CLR <- list(
  ink     = "#000000",
  blue    = "#0072B2",   # primary
  verm    = "#D55E00",   # emphasis / wrong-direction
  green   = "#009E73",   # positive
  sky     = "#56B4E9",
  orange  = "#E69F00",
  purple  = "#CC79A7",
  grey    = "#7F7F7F")

theme_ms <- function(base = 12) {
  theme_classic(base_size = base, base_family = FIG$font) +
    theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.6, colour = "black"),
      axis.text = element_text(colour = "black", size = base - 1),
      axis.title = element_text(size = base),
      plot.title = element_text(size = base + 2, face = "bold",
                                hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = base - 1, colour = "grey25",
                                   hjust = 0.5, margin = margin(b = 6)),
      plot.tag = element_text(size = base + 2, face = "bold"),
      legend.text = element_text(size = base - 1),
      legend.title = element_blank(),
      plot.margin = margin(12, 14, 10, 10))
}

save_ms <- function(p, name, width, height) {
  for (fmt in FIG$formats) {
    f <- P("figures", paste0(name, ".", fmt))
    tryCatch({
      if (fmt == "tiff")
        ggsave(f, p, width = width, height = height, dpi = FIG$dpi,
               device = "tiff", compression = "lzw")
      else if (fmt == "pdf")
        ggsave(f, p, width = width, height = height, device = cairo_pdf)
      else
        ggsave(f, p, width = width, height = height, dpi = FIG$dpi,
               device = fmt)
    }, error = function(e)
      log_msg("  ", fmt, " failed for ", name, ": ",
              conditionMessage(e)))
  }
  log_msg("Figure written: ", name)
}

have <- function(f) {
  ok <- file.exists(P("tables", f))
  if (!ok) log_msg("Missing ", f, "; figure skipped")
  ok
}
rd <- function(f) read_csv(P("tables", f), show_col_types = FALSE)

# ---- F3: reproducibility is network fixedness (keystone) ----------------
fig_reproducibility <- function() {
  if (!have("coexpr_paired_degree.csv")) return(invisible())
  d <- rd("coexpr_paired_degree.csv")
  da <- d[is.finite(d$coexpr_325) & is.finite(d$coexpr_693), ]
  db <- d[is.finite(d$string), ]
  rho <- suppressWarnings(cor(da$coexpr_325, da$coexpr_693,
                              method = "spearman"))

  pa <- ggplot(da, aes(coexpr_325, coexpr_693)) +
    geom_point(size = 0.9, alpha = 0.35, colour = CLR$verm) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("(Spearman %.2f, n = %d)", rho, nrow(da)),
             hjust = 1.05, vjust = 1.6, fontface = "italic", size = 3.8,
             family = FIG$font) +
    labs(x = "Co-expression degree, CGGA-325",
         y = "Co-expression degree, CGGA-693",
         title = "Data-derived network", tag = "A") +
    theme_ms()

  pb <- ggplot(db, aes(string, string)) +
    geom_point(size = 0.9, alpha = 0.35, colour = CLR$blue) +
    annotate("text", x = Inf, y = -Inf,
             label = sprintf("(Spearman 1.00, n = %d)", nrow(db)),
             hjust = 1.05, vjust = -1.0, fontface = "italic", size = 3.8,
             family = FIG$font) +
    labs(x = "STRING degree, any cohort",
         y = "STRING degree, any cohort",
         title = "Fixed database network", tag = "B") +
    theme_ms()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- patchwork::wrap_plots(pa, pb, nrow = 1) +
      patchwork::plot_annotation(
        title = "Hub reproducibility reflects network fixedness",
        theme = theme(plot.title = element_text(
          size = 15, face = "bold", hjust = 0.5, family = FIG$font)))
    save_ms(p, "main_fig_03_reproducibility_keystone", 9.2, 4.8)
  } else {
    save_ms(pa, "F3a_coexpression", 4.8, 4.8)
    save_ms(pb, "F3b_string", 4.8, 4.8)
    log_msg("patchwork absent; F3 saved as two panels")
  }
}

# ---- F1: hub membership is degree --------------------------------------
fig_hub_is_degree <- function() {
  if (!have("hubbias_predictability.csv")) return(invisible())
  h <- rd("hubbias_predictability.csv")
  row <- h[h$target == "conventional hub", , drop = FALSE]
  val <- row$auc[1]; lo <- row$ci_lo[1]; hi <- row$ci_hi[1]
  d <- tibble(x = "STRING degree", auc = val)
  p <- ggplot(d, aes(x, auc)) +
    geom_col(width = 0.42, fill = CLR$verm, colour = "black",
             linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = CLR$grey,
               linewidth = 0.7) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12,
                  linewidth = 0.6) +
    annotate("text", x = 1, y = hi,
             label = sprintf("AUC %.3f", val),
             vjust = -0.8, size = 4.2, family = FIG$font) +
    annotate("text", x = 1.42, y = 0.5, label = "chance", hjust = 0,
             vjust = -0.4, size = 3.4, family = FIG$font,
             colour = CLR$grey) +
    coord_cartesian(ylim = c(0, 1.08), clip = "off") +
    labs(x = NULL, y = "AUC",
         title = "Hub status is predictable from degree") +
    theme_ms() + theme(plot.margin = margin(12, 46, 10, 10))
  save_ms(p, "main_fig_01_hub_is_degree", 5.0, 4.8)
}

# ---- F2: degree tracks citation ----------------------------------------
fig_degree_is_citation <- function() {
  if (!have("degree_vs_papers.csv")) return(invisible())
  d <- rd("degree_vs_papers.csv")
  if (!all(c("ld", "lp") %in% names(d))) {
    log_msg("F2 columns ld/lp not found"); return(invisible())
  }
  r2 <- summary(lm(ld ~ lp, data = d))$r.squared
  p <- ggplot(d, aes(lp, ld)) +
    geom_point(alpha = 0.16, size = 0.7, colour = CLR$blue) +
    geom_smooth(method = "lm", se = FALSE, colour = CLR$verm,
                linewidth = 1.1) +
    annotate("text", x = Inf, y = -Inf,
             label = sprintf("Variance explained %.1f%%", 100 * r2),
             hjust = 1.05, vjust = -1.0, size = 4, family = FIG$font) +
    labs(x = "Publication count (log10)",
         y = "STRING degree (log10)",
         title = "Interactome degree tracks study intensity") +
    theme_ms()
  save_ms(p, "main_fig_02_degree_tracks_citation", 6.0, 5.0)
}

# ---- S11: robustness, degree tracks citation with text-mining removed --
# Companion to F2. The all-channels score in F2 includes the text-mining
# channel, which is itself literature-derived. This panel repeats the
# regression on the text-mining-excluded score to show the relationship is
# not merely the text-mining channel.
fig_degree_is_citation_notm <- function() {
  if (!have("degree_vs_papers_notextmining.csv")) return(invisible())
  d <- rd("degree_vs_papers_notextmining.csv")
  if (!all(c("ld", "lp") %in% names(d))) {
    log_msg("S11 columns ld/lp not found"); return(invisible())
  }
  r2 <- summary(lm(ld ~ lp, data = d))$r.squared
  p <- ggplot(d, aes(lp, ld)) +
    geom_point(alpha = 0.16, size = 0.7, colour = CLR$blue) +
    geom_smooth(method = "lm", se = FALSE, colour = CLR$verm,
                linewidth = 1.1) +
    annotate("text", x = Inf, y = -Inf,
             label = sprintf("Variance explained %.1f%%", 100 * r2),
             hjust = 1.05, vjust = -1.0, size = 4, family = FIG$font) +
    labs(x = "Publication count (log10)",
         y = "STRING degree (log10)",
         title = "Degree tracks study intensity without text mining") +
    theme_ms()
  save_ms(p, "supp_fig_11_degree_tracks_citation_no_textmining", 6.0, 5.0)
}

# ---- F4: the recovery ceiling ------------------------------------------
fig_ceiling <- function() {
  if (!have("propagation_audit.csv") ||
      !have("propagation_audit_calibration.csv")) return(invisible())
  aud <- rd("propagation_audit.csv")
  cal <- rd("propagation_audit_calibration.csv")
  obs <- aud$value[aud$diagnostic == "D1 autocorrelation"]
  exp_r <- aud$value[aud$diagnostic ==
                     "D6 expected recovery at observed autocorr"]
  p <- ggplot(cal, aes(autocorr, recovery)) +
    geom_line(colour = CLR$blue, linewidth = 1.1) +
    geom_point(size = 2.6, colour = CLR$blue) +
    geom_segment(x = obs, xend = obs, y = 0, yend = exp_r,
                 linetype = 2, colour = CLR$verm, linewidth = 0.7) +
    geom_segment(x = 0, xend = obs, y = exp_r, yend = exp_r,
                 linetype = 2, colour = CLR$verm, linewidth = 0.7) +
    annotate("point", x = obs, y = exp_r, size = 3, colour = CLR$verm) +
    annotate("text", x = obs, y = exp_r,
             label = paste0("  ceiling ", sprintf("%.2f", exp_r)),
             hjust = 0, vjust = 1.7, size = 3.9, family = FIG$font,
             colour = CLR$verm) +
    labs(x = "Graph autocorrelation of the signal",
         y = "Achievable held-out recovery",
         title = "Recovery is bounded by the substrate") +
    theme_ms()
  save_ms(p, "main_fig_04_recovery_ceiling", 6.2, 5.0)
}

# ---- F5: degree is anti-informative ------------------------------------
fig_anti_informative <- function() {
  if (!have("propagation_audit.csv")) return(invisible())
  aud <- rd("propagation_audit.csv")
  nb <- aud$value[aud$diagnostic == "D2 AUC neighbour mean"]
  dg <- aud$value[aud$diagnostic == "D2 AUC degree"]
  d <- tibble(
    predictor = factor(c("Neighbourhood mean", "Degree"),
                       levels = c("Neighbourhood mean", "Degree")),
    auc = c(nb, dg),
    fill = c(CLR$green, CLR$verm))
  p <- ggplot(d, aes(predictor, auc, fill = fill)) +
    geom_col(width = 0.54, colour = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = CLR$grey,
               linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.7,
              size = 4.2, family = FIG$font) +
    annotate("text", x = 2.42, y = 0.5, label = "chance", hjust = 0,
             vjust = -0.4, size = 3.4, family = FIG$font,
             colour = CLR$grey) +
    scale_fill_identity() +
    coord_cartesian(ylim = c(0, 0.75), clip = "off") +
    labs(x = NULL,
         y = "AUC for strong differential expression",
         title = "Degree points the wrong way") +
    theme_ms() + theme(plot.margin = margin(12, 46, 10, 10))
  save_ms(p, "main_fig_05_degree_anti_informative", 5.6, 4.8)
}

# ---- F6: the selection ladder ------------------------------------------
fig_ladder <- function() {
  if (!have("downgrade_result.csv")) return(invisible())
  d <- rd("downgrade_result.csv")
  lab <- c(variance = "Variance", DE_tstat = "DE t-statistic",
           relevance = "Relevance (MI)", mRMR = "mRMR",
           degree = "STRING degree",
           cytoHubba_subgraph = "cytoHubba",
           coexpr_degree = "Co-expression degree")
  s <- d %>% group_by(selector) %>%
    summarise(reproducibility = mean(reproducibility),
              network_free = network_free[1], .groups = "drop") %>%
    filter(selector %in% names(lab)) %>%
    mutate(name = lab[selector]) %>%
    arrange(reproducibility) %>%
    mutate(name = factor(name, levels = name))
  p <- ggplot(s, aes(reproducibility, name, fill = network_free)) +
    geom_col(width = 0.66, colour = "black", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", reproducibility)),
              hjust = -0.25, size = 3.8, family = FIG$font) +
    scale_fill_manual(values = c(`TRUE` = CLR$orange, `FALSE` = CLR$blue),
                      labels = c(`TRUE` = "Network-free",
                                 `FALSE` = "Network-based")) +
    scale_x_continuous(limits = c(0, 1.12),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Mean cross-cohort reproducibility", y = NULL,
         title = "Reproducibility rises as data is ignored") +
    theme_ms() + theme(legend.position = c(0.75, 0.25))
  save_ms(p, "main_fig_06_selection_ladder", 7.0, 5.0)
}

main_38 <- function() {
  dir.create(P("figures"), showWarnings = FALSE, recursive = TRUE)
  fig_hub_is_degree()
  fig_degree_is_citation()
  fig_degree_is_citation_notm()
  fig_reproducibility()
  fig_ceiling()
  fig_anti_informative()
  fig_ladder()
  log_msg("Main figures written to ", P("figures"), " at ", FIG$dpi, " dpi")
  invisible(NULL)
}
