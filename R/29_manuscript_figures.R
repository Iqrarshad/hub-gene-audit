# Publication figures for the principal findings.
#
# The discovery analysis lives in the core stages. This stage covers the bias
# and literature results, which are the manuscript's main figures.
#
# Reads saved tables only and computes nothing, so figures can be restyled
# without rerunning any analysis. Any figure whose input table is absent is
# skipped with a message.
#
# Output: TIFF (LZW), PDF and PNG at the resolution set in FIG$dpi.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tibble)
  library(tidyr); library(scales)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

theme_pub <- function(base = 12) {
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

save_fig <- function(p, name, width = 7, height = 5, title = NULL) {
  if (!is.null(title))
    p <- p + labs(title = title)
  for (fmt in FIG$formats) {
    f <- P("figures", paste0(name, ".", fmt))
    ok <- tryCatch({
      if (fmt == "tiff")
        ggsave(f, p, width = width, height = height, dpi = FIG$dpi,
               device = "tiff", compression = "lzw")
      else if (fmt == "pdf")
        ggsave(f, p, width = width, height = height, device = cairo_pdf)
      else
        ggsave(f, p, width = width, height = height, dpi = FIG$dpi,
               device = fmt)
      TRUE
    }, error = function(e) { log_msg("  ", fmt, " failed: ",
                                     conditionMessage(e)); FALSE })
  }
  log_msg("Figure saved: ", name)
}

have <- function(f) {
  ok <- file.exists(P("tables", f))
  if (!ok) log_msg("  skipped, missing: ", f)
  ok
}

# --- Fig: published hub genes sit at the top of the degree distribution --
fig_literature_degree <- function() {
  if (!have("lit30_final_by_direction.csv")) return(invisible(NULL))
  d <- read_csv(P("tables", "lit30_final_by_direction.csv"),
                show_col_types = FALSE) %>%
    mutate(direction = factor(direction, levels = c("down", "mixed", "up"),
                              labels = c("Downregulated", "Mixed",
                                         "Upregulated")))
  p <- ggplot(d, aes(x = direction, y = median_pct, fill = direction)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 50, linetype = 2, colour = "grey40") +
    geom_text(aes(label = sprintf("%.1f", median_pct)), vjust = -0.5,
              size = 3.2, family = FIG$font) +
    scale_fill_manual(values = FIG$okabe_ito[c(3, 5, 7)], guide = "none") +
    scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
    labs(x = NULL,
         y = "Median interactome degree percentile") +
    theme_pub()
  save_fig(p, "supp_fig_01", 5.0, 4.0, title = "Published hubs sit at high degree")
}

# --- Fig: reproducibility falls as results are summarised ---------------
fig_reproducibility_gradient <- function() {
  f <- P("tables", "edge_overlap_curve.csv")
  cmp <- P("tables", "edge_vs_hub_agreement.csv")
  if (!file.exists(cmp)) { log_msg("  skipped, missing edge_vs_hub_agreement"); return(invisible(NULL)) }
  d <- read_csv(cmp, show_col_types = FALSE) %>%
    mutate(level = factor(level, levels = level))
  p <- ggplot(d, aes(x = level, y = agreement, group = 1)) +
    geom_line(colour = FIG$okabe_ito[6], linewidth = 0.8) +
    geom_point(size = 3, colour = FIG$okabe_ito[6]) +
    geom_text(aes(label = sprintf("%.3f", agreement)), vjust = -1,
              size = 3.2, family = FIG$font) +
    scale_y_continuous(limits = c(0, 0.8)) +
    labs(x = NULL, y = "Cross-cohort agreement") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_fig(p, "supp_fig_02", 6.5, 4.5, title = "Reproducibility falls with summarisation")
}

# --- Fig: composition attenuation, contrast as internal control ---------
fig_composition <- function() {
  f1 <- P("tables", "composition_adjustment_summary.csv")
  f2 <- P("tables", "rnaseq_contrast_attenuation.csv")
  parts <- list()
  if (file.exists(f2)) {
    parts[[1]] <- read_csv(f2, show_col_types = FALSE) %>%
      select(label = dataset, attenuation = median_attenuation)
  }
  if (file.exists(f1)) {
    parts[[2]] <- read_csv(f1, show_col_types = FALSE) %>%
      select(label = dataset, attenuation = median_attenuation)
  }
  if (!length(parts)) { log_msg("  skipped, no attenuation tables"); return(invisible(NULL)) }
  d <- bind_rows(parts) %>% filter(!is.na(attenuation)) %>%
    mutate(control = grepl("HGG_vs_LGG|grade", label, ignore.case = TRUE),
           label = reorder(label, attenuation))
  p <- ggplot(d, aes(x = attenuation, y = label, fill = control)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(`FALSE` = FIG$okabe_ito[7],
                                 `TRUE`  = FIG$okabe_ito[4]),
                      labels = c(`FALSE` = "Tumour vs normal",
                                 `TRUE`  = "Tumour vs tumour (control)"),
                      name = NULL) +
    labs(x = "Median attenuation after composition adjustment (%)", y = NULL) +
    theme_pub()
  save_fig(p, "supp_fig_03", 6.8, 4.6, title = "Composition adjustment attenuates signal")
}

# --- Fig: attempted corrections -----------------------------------------
fig_corrections <- function() {
  rows <- list()
  add <- function(f, mcol, vcol, tag) {
    if (!file.exists(P("tables", f))) return(NULL)
    d <- read_csv(P("tables", f), show_col_types = FALSE)
    if (!all(c(mcol, vcol) %in% names(d))) return(NULL)
    tibble(method = d[[mcol]], auc = d[[vcol]], source = tag)
  }
  rows[[1]] <- add("citation_penalised_comparison.csv", "method",
                   "auc_citation_bias", "citation")
  rows[[2]] <- add("coessentiality_comparison.csv", "method",
                   "auc_citation_bias", "coessentiality")
  rows[[3]] <- add("reproducibility_eigenvector.csv", "weight",
                   "auc_annotation_bias", "eigenvector")
  d <- bind_rows(rows) %>% filter(!is.na(auc))
  if (!nrow(d)) { log_msg("  skipped, no correction tables"); return(invisible(NULL)) }
  d <- d %>% distinct(method, .keep_all = TRUE) %>%
    mutate(method = reorder(method, -auc))
  p <- ggplot(d, aes(x = auc, y = method)) +
    geom_col(fill = FIG$okabe_ito[3], width = 0.7) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = FIG$okabe_ito[7]) +
    scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    labs(x = "AUC of prior attention predicting the nominated genes",
         y = NULL) +
    theme_pub()
  save_fig(p, "supp_fig_04", 6.8, 4.6, title = "Corrections retain citation bias")
}

# --- Fig: stability ceiling ---------------------------------------------
fig_ceiling <- function() {
  if (!have("stability_ceiling_grid.csv")) return(invisible(NULL))
  d <- read_csv(P("tables", "stability_ceiling_grid.csv"),
                show_col_types = FALSE) %>%
    group_by(retained_frac) %>%
    summarise(oracle = mean(oracle), random = mean(random), .groups = "drop") %>%
    pivot_longer(c(oracle, random), names_to = "method", values_to = "v") %>%
    mutate(method = factor(method, levels = c("oracle", "random"),
                           labels = c("Perfect method (ceiling)",
                                      "Random (floor)")))
  p <- ggplot(d, aes(x = retained_frac, y = v, colour = method,
                     group = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
    annotate("text", x = 0.55, y = 1.03, label = "Value assumed by convention",
             size = 3, family = FIG$font, colour = "grey40", hjust = 0) +
    scale_colour_manual(values = FIG$okabe_ito[c(6, 1)], name = NULL) +
    scale_y_continuous(limits = c(0, 1.1)) +
    labs(x = "Fraction of genes retained",
         y = "Split-half Jaccard") +
    theme_pub()
  save_fig(p, "supp_fig_05", 6.5, 4.5, title = "Stability sits below its ceiling")
}

# --- Fig: recurrence against the degree-matched null --------------------
fig_recurrence <- function() {
  if (!have("lit30_final_recurrence_null.csv")) return(invisible(NULL))
  d <- read_csv(P("tables", "lit30_final_recurrence_null.csv"),
                show_col_types = FALSE)
  obs <- d$observed_recur[1]
  p <- ggplot(d, aes(x = null, y = mean_recur)) +
    geom_col(fill = FIG$okabe_ito[3], width = 0.55) +
    geom_errorbar(aes(ymin = mean_recur - sd, ymax = mean_recur + sd),
                  width = 0.15) +
    geom_hline(yintercept = obs, colour = FIG$okabe_ito[7], linewidth = 0.8) +
    annotate("text", x = 0.6, y = obs + 2,
             label = paste("observed =", obs), hjust = 0, size = 3.2,
             family = FIG$font, colour = FIG$okabe_ito[7]) +
    labs(x = NULL, y = "Genes recurring in two or more studies") +
    theme_pub()
  save_fig(p, "supp_fig_06", 6, 4.5, title = "Recurrence matches the degree null")
}

main_29 <- function() {
  dir.create(P("figures"), showWarnings = FALSE, recursive = TRUE)
  log_msg("Formats: ", paste(FIG$formats, collapse = ", "),
          " at ", FIG$dpi, " dpi")
  fig_literature_degree()
  fig_recurrence()
  fig_reproducibility_gradient()
  fig_composition()
  fig_corrections()
  fig_ceiling()
  n <- length(list.files(P("figures")))
  log_msg("29_manuscript_figures complete. ", n, " files in figures/")
}
