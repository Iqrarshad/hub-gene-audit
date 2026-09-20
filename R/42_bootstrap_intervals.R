# BCa bootstrap intervals for the headline AUCs.
#
# Several AUCs are reported in the manuscript as point estimates only. This
# stage attaches a bias-corrected and accelerated (BCa) interval with 5000
# resamples to each, from the same cached tables the earlier stages wrote,
# and collects them in one file. Nothing upstream is recomputed and no
# earlier output is changed; this stage only reads and adds.
#
# Covered claims:
#   Sec 3.2  publication count predicts hub membership          (papers)
#   Sec 3.2  degree residual predicts hub membership (0.739)    (residual)
#   Sec 3.2  degree predicts hub membership over the interactome (degree)
#   Sec 3.3  degree and publication AUCs for the 131 published hubs
#   Sec 3.5  attention predicts which pathways are enriched (0.708)
#   Discussion  DepMap glioma-selective dependency of the hub set
#
# The candidate-network degree->hub AUC of 0.994 already carries a DeLong
# interval from stage 12 and is not recomputed here.
#
# collect_auc_inputs() is shared with the multi-seed stage (43) so both use
# the identical class labels and predictors.
#
# Output: tables/bootstrap_intervals.csv

source("R/00h_bootstrap.R")

N_BOOT_CI <- 5000    # resamples for every interval in this stage

# Assemble every (analysis, quantity, y, x) the interval stages evaluate.
# Returns a list of records; missing inputs are dropped with a message so
# the stage degrades rather than fails.
collect_auc_inputs <- function() {
  recs <- list()
  keep <- function(analysis, hub_set, quantity, universe, claim, y, x) {
    y <- as.integer(y)
    if (sum(y == 1L, na.rm = TRUE) < 3L ||
        length(unique(stats::na.omit(x))) < 2L) return(invisible())
    recs[[length(recs) + 1]] <<- list(analysis = analysis, hub_set = hub_set,
      quantity = quantity, universe = universe, claim = claim, y = y, x = x)
  }

  dvp <- read_table_safe("degree_vs_papers.csv")
  own <- consensus_hub_genes()
  pub <- published_hub_genes()

  if (!is.null(dvp) && all(c("gene", "ld", "lp", "degree_resid") %in% names(dvp))) {
    g <- toupper(dvp$gene)
    if (length(own) >= 3L) {
      y <- as.integer(g %in% toupper(own))
      keep("hub components (our hubs)", "our cytoHubba hubs",
           "degree predicts hub", "interactome (all channels)",
           "Sec 3.2 degree component", y, dvp$ld)
      keep("hub components (our hubs)", "our cytoHubba hubs",
           "publications predict hub", "interactome (all channels)",
           "Sec 3.2 publication AUC 0.948", y, dvp$lp)
      keep("hub components (our hubs)", "our cytoHubba hubs",
           "degree residual predicts hub", "interactome (all channels)",
           "Sec 3.2 residual AUC 0.739", y, dvp$degree_resid)
      if ("degree_from_papers" %in% names(dvp))
        keep("hub components (our hubs)", "our cytoHubba hubs",
             "citation-explained degree predicts hub",
             "interactome (all channels)",
             "Sec 3.2 citation-explained degree", y, dvp$degree_from_papers)
    }
    if (length(pub) >= 3L) {
      yp <- as.integer(g %in% toupper(pub))
      keep("hub components (published hubs)", "published hub genes",
           "degree predicts hub", "interactome (all channels)",
           "Sec 3.3 degree AUC 0.876", yp, dvp$ld)
      keep("hub components (published hubs)", "published hub genes",
           "publications predict hub", "interactome (all channels)",
           "Sec 3.3 publication AUC 0.821", yp, dvp$lp)
      keep("hub components (published hubs)", "published hub genes",
           "degree residual predicts hub", "interactome (all channels)",
           "Sec 3.3 residual", yp, dvp$degree_resid)
    }
  } else log_msg("  degree_vs_papers.csv absent; hub-component inputs skipped")

  enr <- read_table_safe("annotation_bias_enrichment.csv")
  if (!is.null(enr) && all(c("attention", "enriched") %in% names(enr)))
    keep("enrichment analogue", "consensus hub pathways",
         "attention predicts enriched pathway", "tested KEGG pathways",
         "Sec 3.5 enrichment AUC 0.708",
         as.integer(enr$enriched), as.numeric(enr$attention))
  else log_msg("  annotation_bias_enrichment.csv absent; enrichment skipped")

  dep <- read_table_safe("depmap_glioma_selectivity.csv")
  if (!is.null(dep) && all(c("gene", "selectivity") %in% names(dep))) {
    gd <- toupper(dep$gene)
    if (length(own) >= 3L)
      keep("depmap selectivity", "our cytoHubba hubs",
           "hub membership from glioma-selective dependency",
           "DepMap genes", "Discussion DepMap selectivity",
           as.integer(gd %in% toupper(own)), as.numeric(dep$selectivity))
    if (length(pub) >= 3L)
      keep("depmap selectivity", "published hub genes",
           "hub membership from glioma-selective dependency",
           "DepMap genes", "Discussion DepMap selectivity (published)",
           as.integer(gd %in% toupper(pub)), as.numeric(dep$selectivity))
  } else log_msg("  depmap_glioma_selectivity.csv absent; DepMap skipped")

  recs
}

main_42 <- function() {
  set.seed(SEED)
  log_msg("=================================================")
  log_msg("BCa BOOTSTRAP INTERVALS (", N_BOOT_CI, " resamples)")
  log_msg("=================================================")

  recs <- collect_auc_inputs()
  if (!length(recs)) {
    log_msg("No AUC inputs found. Run the upstream stages first ",
            "(24 for hub components, 37 for enrichment, 23 for DepMap).")
    return(invisible(NULL))
  }

  rows <- list()
  for (r in recs) {
    ci <- bca_auc_ci(r$y, r$x, R = N_BOOT_CI, seed = SEED)
    rows[[length(rows) + 1]] <- data.frame(
      analysis = r$analysis, hub_set = r$hub_set, quantity = r$quantity,
      universe = r$universe, n_positives = sum(r$y == 1L, na.rm = TRUE),
      estimate = round(unname(ci["estimate"]), 4),
      ci_lo = round(unname(ci["lo"]), 4),
      ci_hi = round(unname(ci["hi"]), 4),
      method = attr(ci, "method"), n_boot = unname(ci["R"]),
      seed = SEED, manuscript_claim = r$claim, stringsAsFactors = FALSE)
    log_msg("  ", r$analysis, " / ", r$quantity, ": ",
            round(unname(ci["estimate"]), 3), " [",
            round(unname(ci["lo"]), 3), ", ", round(unname(ci["hi"]), 3),
            "] (", attr(ci, "method"), ")")
  }

  out <- do.call(rbind, rows)
  utils::write.csv(out, P("tables", "bootstrap_intervals.csv"),
                   row.names = FALSE)
  log_msg("Wrote bootstrap_intervals.csv (", nrow(out), " rows)")
  print(out[, c("analysis", "quantity", "estimate", "ci_lo", "ci_hi",
                "method")])
  log_msg("Note: candidate-network degree->hub AUC 0.994 keeps its stage 12 ",
          "DeLong interval and is not re-bootstrapped here.")
  log_msg("42_bootstrap_intervals complete.")
  invisible(out)
}
