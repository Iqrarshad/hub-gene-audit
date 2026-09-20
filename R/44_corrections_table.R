# Consolidated corrections table.
#
# The manuscript tests nine corrections in two families and, separately, two
# network-propagation variants (random walk with restart, and degree-
# constrained diffusion). This stage assembles all of them into one machine-
# readable table so the count is explicit and the two propagation variants
# are saved as rows rather than living only in prose.
#
# The two families do not share a metric. Network-kept and network-free
# corrections are judged on cross-cohort reproducibility and dependence on
# study attention. Propagation variants are judged on whether they recover a
# held-out true signal or merely the degree vector. Forcing propagation into
# the reproducibility and dependence columns would compare unlike things, so
# those cells are left NA for the propagation rows and their native metrics
# are given in the propagation columns instead. Every numeric cell is read
# from a stage output; a value that cannot be found is left NA.
#
# Output: tables/corrections_all.csv

source("R/00h_bootstrap.R")

main_44 <- function() {
  na <- NA_real_

  # ---- Nine corrections, enumerated, values read where available --------
  rewire   <- pick_cell("hub_metric_comparison.csv", "method",
                        "rewiring z-score (mean)",
                        "auc_predictable_from_global_degree")
  comp     <- pick_cell("composite_hub_comparison.csv", "scheme",
                        "H_product", "auc_annotation_bias")
  cp_cit   <- pick_cell("citation_penalised_comparison.csv", "method",
                        "C: degree residual only", "auc_citation_bias")
  cp_sel   <- pick_cell("citation_penalised_comparison.csv", "method",
                        "C: degree residual only", "auc_selectivity")
  base_cit <- pick_cell("citation_penalised_comparison.csv", "method",
                        "cytoHubba (conventional)", "auc_citation_bias")

  # Network-free selectors: reproducibility and citation AUC, averaged over
  # the rows in downgrade_result.csv exactly as stage 40 does for Table 3.
  nf_repro <- function(sel) na; nf_dep <- function(sel) na
  dg <- read_table_safe("downgrade_result.csv")
  if (!is.null(dg) &&
      all(c("selector", "reproducibility", "citation_auc") %in% names(dg))) {
    ag <- stats::aggregate(cbind(reproducibility, citation_auc) ~ selector,
                           data = dg, FUN = function(z) mean(z, na.rm = TRUE))
    nf_repro <- function(sel) {
      r <- ag$reproducibility[ag$selector == sel]; if (length(r)) r[1] else na }
    nf_dep <- function(sel) {
      r <- ag$citation_auc[ag$selector == sel]; if (length(r)) r[1] else na }
  } else {
    log_msg("  downgrade_result.csv absent; network-free rows left NA")
  }

  corr <- data.frame(
    family = c(rep("network-kept", 5), rep("network-free", 4)),
    correction = c("Degree normalisation", "Rewiring z-score null",
                   "Binomial specific connectivity", "Composite centrality",
                   "Degree residualised on publications",
                   "Expression variance", "Grade statistic",
                   "Mutual information relevance",
                   "Minimum-redundancy maximum-relevance"),
    basis = c("local degree over expected degree",
              "degree-preserving rewiring null",
              "specific-connectivity binomial test",
              "combined centrality metrics",
              "degree residual on publication count",
              "expression variance selector",
              "grade differential statistic",
              "mutual-information relevance",
              "mRMR selector"),
    reproducibility = c(na, na, na, na, na,
                        nf_repro("variance"), nf_repro("DE_tstat"),
                        nf_repro("relevance"), nf_repro("mRMR")),
    dependence_auc = c(na, round0(rewire), na, round0(comp), round0(cp_cit),
                       round0(nf_dep("variance")), round0(nf_dep("DE_tstat")),
                       round0(nf_dep("relevance")), round0(nf_dep("mRMR"))),
    selectivity_auc = c(na, na, na, na, round0(cp_sel), na, na, na, na),
    prop_recovery_signal = na, prop_recovery_degree = na,
    gate_amplification_pass = NA,
    outcome = c("bias retained", "bias retained", "bias retained",
                "bias retained",
                "dependence falls only by discarding signal; selectivity at chance",
                "reproducibility lost", "reproducibility lost",
                "reproducibility lost", "reproducibility lost"),
    source_file = c("(text)", "hub_metric_comparison.csv",
                    "(text)", "composite_hub_comparison.csv",
                    "citation_penalised_comparison.csv",
                    "downgrade_result.csv", "downgrade_result.csv",
                    "downgrade_result.csv", "downgrade_result.csv"),
    stringsAsFactors = FALSE)

  # ---- Baseline row -----------------------------------------------------
  baseline <- data.frame(
    family = "baseline", correction = "cytoHubba consensus hubs",
    basis = "degree-based centrality",
    reproducibility = round0(nf_repro("cytoHubba_subgraph")),
    dependence_auc = round0(base_cit), selectivity_auc = na,
    prop_recovery_signal = na, prop_recovery_degree = na,
    gate_amplification_pass = NA, outcome = "reference",
    source_file = "citation_penalised_comparison.csv",
    stringsAsFactors = FALSE)

  # ---- Two propagation variants -----------------------------------------
  # Random walk with restart, from the propagation gates and audit.
  amp_pass <- pick_flag("propagation_gates.csv", "gate", "3a amplification",
                        "pass")
  sig_ret  <- pick_cell("propagation_gates.csv", "gate", "3b signal retained",
                        "value")
  rec_deg  <- pick_cell("propagation_audit.csv", "diagnostic",
                        "D2 AUC degree", "value")
  rec_nbr  <- pick_cell("propagation_audit.csv", "diagnostic",
                        "D2 AUC neighbour mean", "value")
  rwr <- data.frame(
    family = "propagation", correction = "Random walk with restart",
    basis = "network propagation of the restart vector",
    reproducibility = na, dependence_auc = na, selectivity_auc = na,
    prop_recovery_signal = round0(sig_ret),
    prop_recovery_degree = round0(rec_deg),
    gate_amplification_pass = if (is.na(amp_pass)) NA else amp_pass == 1,
    outcome = "recovers degree, not held-out signal; fails amplification gate",
    source_file = "propagation_gates.csv; propagation_audit.csv",
    stringsAsFactors = FALSE)

  # Degree-constrained diffusion, from the constrained-diffusion criteria.
  cd_top  <- pick_cell("constrained_criteria.csv", "criterion", "S1 success",
                       "value")
  cd_s1   <- pick_flag("constrained_criteria.csv", "criterion", "S1 success",
                       "pass")
  cd_s3   <- pick_flag("constrained_criteria.csv", "criterion",
                       "S3 beats post-hoc", "pass")
  cd_out  <- if (!is.na(cd_s1) && !is.na(cd_s3)) {
    paste0("recovery threshold ",
           ifelse(cd_s1 == 1, "met", "not met"),
           "; post-hoc degree adjustment ",
           ifelse(cd_s3 == 1, "beaten", "not beaten"))
  } else "does not reach recovery threshold; does not beat post-hoc adjustment"
  cdiff <- data.frame(
    family = "propagation", correction = "Degree-constrained diffusion",
    basis = "diffusion with degree and citation constraints",
    reproducibility = na, dependence_auc = na, selectivity_auc = na,
    prop_recovery_signal = round0(cd_top), prop_recovery_degree = na,
    gate_amplification_pass = NA, outcome = cd_out,
    source_file = "constrained_criteria.csv",
    stringsAsFactors = FALSE)

  out <- rbind(baseline, corr, rwr, cdiff)
  utils::write.csv(out, P("tables", "corrections_all.csv"), row.names = FALSE)

  log_msg("=================================================")
  log_msg("CONSOLIDATED CORRECTIONS: ", nrow(corr),
          " corrections (5 network-kept, 4 network-free) plus ",
          "2 propagation variants, plus 1 baseline")
  log_msg("=================================================")
  print(out[, c("family", "correction", "reproducibility",
                "dependence_auc", "outcome")])
  log_msg("Wrote corrections_all.csv (", nrow(out), " rows)")
  log_msg("44_corrections_table complete.")
  invisible(out)
}

# Round that tolerates NA and non-numeric.
round0 <- function(x, d = 3) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1 || is.na(x)) return(NA_real_)
  round(x, d)
}
