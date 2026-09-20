# Figures for the propagation audit.
#
# Reads saved tables only and computes nothing. Any figure whose input
# table is absent is skipped with a message.
#
# Output: TIFF (LZW), PDF and PNG at the resolution set in FIG$dpi.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tibble)
  library(tidyr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

theme_pub52 <- function(base = 12) {
  theme_classic(base_size = base, base_family = FIG$font) +
    theme(
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.6, colour = "black"),
      axis.text = element_text(colour = "black", size = base - 1),
      legend.key = element_blank(),
      plot.title = element_text(size = base + 2, face = "bold",
                                hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = base - 1, colour = "grey25",
                                   hjust = 0.5),
      plot.tag = element_text(size = base + 2, face = "bold"))
}

save_fig52 <- function(p, name, width = 7, height = 5, title = NULL) {
  if (!is.null(title)) p <- p + labs(title = title)
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
    }, error = function(e) log_msg("  ", fmt, " failed: ",
                                   conditionMessage(e)))
  }
  log_msg("Figure saved: ", name)
}

have <- function(f) {
  ok <- file.exists(P("tables", f))
  if (!ok) log_msg("Missing ", f, "; figure skipped")
  ok
}

# ---- A: calibration curve ----------------------------------------------
fig_calibration <- function(aud, cal) {
  obs <- aud$value[aud$diagnostic == "D1 autocorrelation"]
  exp_r <- aud$value[aud$diagnostic ==
                     "D6 expected recovery at observed autocorr"]
  p <- ggplot(cal, aes(autocorr, recovery)) +
    geom_line(colour = FIG$okabe_ito[6], linewidth = 0.7) +
    geom_point(size = 2.4, colour = FIG$okabe_ito[6]) +
    geom_segment(aes(x = obs, xend = obs, y = -Inf, yend = exp_r),
                 linetype = 2, colour = FIG$okabe_ito[7],
                 linewidth = 0.5) +
    geom_segment(aes(x = -Inf, xend = obs, y = exp_r, yend = exp_r),
                 linetype = 2, colour = FIG$okabe_ito[7],
                 linewidth = 0.5) +
    annotate("text", x = obs, y = exp_r,
             label = paste0("observed autocorrelation ", round(obs, 3),
                            "\nexpected recovery ", round(exp_r, 3)),
             hjust = 1.05, vjust = -0.4, size = 3.1, family = FIG$font,
             colour = FIG$okabe_ito[7]) +
    labs(x = "Graph autocorrelation of the signal",
         y = "Held-out recovery (Spearman)") +
    theme_pub52()
  save_fig52(p, "supp_fig_07", 6.5, 4.6, title = "Recovery scales with autocorrelation")
}

# ---- B: recovery against the ceiling -----------------------------------
fig_recovery <- function(aud, bl, best) {
  exp_r <- aud$value[aud$diagnostic ==
                     "D6 expected recovery at observed autocorr"]
  a <- bl$alpha[which.max(bl$posthoc)]
  d <- tibble(
    method = factor(
      c("Unconstrained\npropagation", "Constrained\noperator",
        "Post-hoc degree\nadjustment", "Expected at this\nautocorrelation"),
      levels = c("Unconstrained\npropagation", "Constrained\noperator",
                 "Post-hoc degree\nadjustment",
                 "Expected at this\nautocorrelation")),
    value = c(bl$raw[bl$alpha == a],
              max(best$rec_true[best$alpha == a]),
              bl$posthoc[bl$alpha == a], exp_r),
    kind = c("observed", "observed", "observed", "ceiling"))
  p <- ggplot(d, aes(method, value, fill = kind)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.5,
              size = 3.2, family = FIG$font) +
    scale_fill_manual(values = c(observed = FIG$okabe_ito[3],
                                 ceiling = FIG$okabe_ito[5]),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(x = NULL, y = "Held-out recovery (Spearman)") +
    theme_pub52()
  save_fig52(p, "supp_fig_08", 6.8, 4.6, title = "Constraint does not beat adjustment")
}

# ---- C: what predicts strong differential expression -------------------
fig_auc <- function(aud) {
  d <- tibble(
    predictor = factor(c("Neighbourhood mean", "Degree"),
                       levels = c("Neighbourhood mean", "Degree")),
    auc = c(aud$value[aud$diagnostic == "D2 AUC neighbour mean"],
            aud$value[aud$diagnostic == "D2 AUC degree"]))
  p <- ggplot(d, aes(predictor, auc)) +
    geom_col(width = 0.5, fill = FIG$okabe_ito[3]) +
    geom_hline(yintercept = 0.5, linetype = 2,
               colour = FIG$okabe_ito[7], linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.5,
              size = 3.4, family = FIG$font) +
    annotate("text", x = 2.35, y = 0.5, label = "chance", hjust = 0,
             vjust = -0.4, size = 3, family = FIG$font,
             colour = FIG$okabe_ito[7]) +
    coord_cartesian(ylim = c(0, 0.75), clip = "off") +
    labs(x = NULL,
         y = "AUC for membership of the top decile by effect size") +
    theme_pub52() +
    theme(plot.margin = margin(6, 34, 6, 6))
  save_fig52(p, "supp_fig_09", 5.6, 4.4, title = "Degree fails to predict expression")
}

# ---- D: degree dependence by normalisation -----------------------------
fig_normalisation <- function(nrm) {
  d <- nrm %>%
    select(normalisation, rec_true, rec_degree, rec_partial) %>%
    pivot_longer(-normalisation, names_to = "measure",
                 values_to = "value") %>%
    mutate(measure = recode(measure,
             rec_true = "vs true values",
             rec_degree = "vs degree",
             rec_partial = "vs true, degree partialled"),
           normalisation = recode(normalisation,
             symmetric = "Symmetric", row = "Row"))
  p <- ggplot(d, aes(measure, value, fill = normalisation)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    scale_fill_manual(values = FIG$okabe_ito[c(3, 6)], name = NULL) +
    labs(x = NULL, y = "Spearman correlation of the diffused value") +
    theme_pub52() + theme(legend.position = "top")
  save_fig52(p, "supp_fig_10", 6.6, 4.6, title = "Normalisation shifts degree dependence")
}

main_32 <- function() {
  dir.create(P("figures"), showWarnings = FALSE, recursive = TRUE)

  if (have("propagation_audit.csv") &&
      have("propagation_audit_calibration.csv")) {
    aud <- read_csv(P("tables", "propagation_audit.csv"),
                    show_col_types = FALSE)
    cal <- read_csv(P("tables", "propagation_audit_calibration.csv"),
                    show_col_types = FALSE)
    fig_calibration(aud, cal)
    fig_auc(aud)

    if (have("constrained_baselines.csv") &&
        have("constrained_recovery.csv")) {
      bl <- read_csv(P("tables", "constrained_baselines.csv"),
                     show_col_types = FALSE)
      rec <- read_csv(P("tables", "constrained_recovery.csv"),
                      show_col_types = FALSE)
      best <- rec %>% group_by(alpha, constraint) %>%
        slice_max(rec_true, n = 1, with_ties = FALSE) %>% ungroup()
      fig_recovery(aud, bl, best)
    }
  }

  if (have("propagation_audit_norm.csv")) {
    nrm <- read_csv(P("tables", "propagation_audit_norm.csv"),
                    show_col_types = FALSE)
    fig_normalisation(nrm)
  }

  log_msg("Propagation figures written to ", P("figures"),
          " at ", FIG$dpi, " dpi")
  invisible(NULL)
}
