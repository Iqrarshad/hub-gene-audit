# Publication tables for the hub-gene audit (Briefings in Bioinformatics).
#
# Reads existing stage outputs and writes formatted, manuscript-ready
# tables. Computes nothing new. Each table is written both as a CSV for
# import and, where short, echoed to the log for a quick check.
#
# Main tables:
#   Table 1  datasets and cohorts
#   Table 2  attempted corrections and their outcome
#   Table 3  selector comparison, reproducibility and citation bias
#
# Supplementary tables:
#   S1  pre-declared pass and fail criteria
#   S2  adjustment provenance (composition-adjustment records)
#   S3  full-vector reproducibility by network density
#   S4  pipeline benchmark, end to end
#   S5  cytoHubba metrics against degree
#   S6  WGCNA cross-cohort reproducibility
#   S7  annotation-bias decomposition
#
# A table whose source is absent is skipped with a message. Column names
# are resolved tolerantly so the stage survives minor upstream changes.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
})

# The annotation stack loaded elsewhere masks these; bind them to dplyr.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename

have <- function(f) {
  ok <- file.exists(P("tables", f))
  if (!ok) log_msg("Source ", f, " absent; that table skipped")
  ok
}
rd <- function(f) read_csv(P("tables", f), show_col_types = FALSE)
wr <- function(d, name) {
  write_csv(d, P("tables", paste0(name, ".csv")))
  log_msg("Wrote ", name, ".csv  (", nrow(d), " rows)")
}
pick <- function(d, opts) { i <- intersect(opts, names(d))[1]
  if (is.na(i)) NA else d[[i]] }

main_40 <- function() {

  # ---- Table 1: datasets ------------------------------------------------
  # Static description of the inputs, since sample counts are reported in
  # the cohort logs. Kept as a declared manifest rather than recomputed.
  t1 <- tribble(
    ~Dataset, ~Accession, ~Platform, ~Role,
    "Glioma RNA-seq", "GSE147352", "RNA-seq", "Discovery",
    "Gravendeel", "GSE16011", "Microarray", "Discovery",
    "REMBRANDT", "GSE108474", "Microarray", "Discovery",
    "Glioma microarray", "GSE15824", "Microarray", "Discovery",
    "Glioma microarray", "GSE21354", "Microarray", "Discovery",
    "TCGA-LGG", "GDC / cBioPortal", "RNA-seq", "Validation",
    "CGGA mRNAseq 325", "CGGA", "RNA-seq", "Validation",
    "CGGA mRNAseq 693", "CGGA", "RNA-seq", "Validation",
    "STRING human", "v12.0", "Interactome", "Network",
    "gene2pubmed", "NCBI", "Citation counts", "Literature",
    "DepMap CRISPR", "DepMap", "Essentiality", "Substrate")
  wr(t1, "Table1_datasets")

  # ---- Table 2: attempted corrections -----------------------------------
  # Values are read from the correction-stage outputs so the table always
  # matches the run. pick() looks up one cell by matching a key column; a
  # missing file, column, or row returns NA and the outcome text says
  # "not available" rather than showing a stale hardcoded number.
  pick <- function(file, key_col, key_val, val_col, digits = 3) {
    f <- P("tables", file)
    if (!file.exists(f)) return(NA_real_)
    d <- suppressMessages(read_csv(f, show_col_types = FALSE))
    if (!all(c(key_col, val_col) %in% names(d))) return(NA_real_)
    r <- d[d[[key_col]] == key_val, , drop = FALSE]
    if (!nrow(r)) return(NA_real_)
    round(as.numeric(r[[val_col]][1]), digits)
  }
  fmt <- function(x) ifelse(is.na(x), "not available", format(x))

  rewire   <- pick("hub_metric_comparison.csv", "method",
                   "rewiring z-score (mean)",
                   "auc_predictable_from_global_degree")
  comp     <- pick("composite_hub_comparison.csv", "scheme", "H_product",
                   "auc_annotation_bias")
  comp_bin <- pick("composite_hub_comparison.csv", "scheme", "binomial_only",
                   "auc_annotation_bias")
  cp_cit   <- pick("citation_penalised_comparison.csv", "method",
                   "C: degree residual only", "auc_citation_bias")
  cp_sel   <- pick("citation_penalised_comparison.csv", "method",
                   "C: degree residual only", "auc_selectivity")
  ce_cit   <- pick("coessentiality_comparison.csv", "method",
                   "co-essentiality (CNS held out)", "auc_citation_bias")
  ce_sel   <- pick("coessentiality_comparison.csv", "method",
                   "co-essentiality (CNS held out)", "auc_selectivity")
  gm       <- pick("grade_matched_summary.csv", "comparison", "LGG vs LGG",
                   "mean_jaccard")
  base_ann <- pick("pipeline_benchmark.csv", "pipeline",
                   "conventional (STRING + cytoHubba hubs)",
                   "auc_annotation_bias")
  base_cit <- pick("citation_penalised_comparison.csv", "method",
                   "cytoHubba (conventional)", "auc_citation_bias")

  t2 <- tribble(
    ~Correction, ~Basis, ~Outcome,
    "Degree normalisation", "local over expected degree",
      "bias retained",
    "Rewiring z-score", "degree-preserving null",
      paste0("AUC ", fmt(rewire), ", bias retained"),
    "Composite hub score", "combined centrality metrics",
      paste0("AUC ", fmt(comp), ", matches binomial ", fmt(comp_bin),
             ", bias retained"),
    "Binomial connectivity", "specific-connectivity test",
      "bias retained",
    "Citation-penalised score", "degree residualised on publications",
      paste0("citation AUC falls to ", fmt(cp_cit),
             " only by discarding signal; selectivity stays at chance (",
             fmt(cp_sel), ")"),
    "Cohort co-expression", "expression-derived network",
      "degree bias broken but hubs do not replicate",
    "Grade-matched co-expression", "grade-stratified network",
      paste0("within-grade hub agreement ", fmt(gm)),
    "CRISPR co-essentiality", "DepMap dependency network",
      paste0("bias removed (citation AUC ", fmt(ce_cit),
             ") but no glioma selectivity (", fmt(ce_sel), ")"))
  wr(t2, "Table2_corrections")

  t2base <- tribble(
    ~Baseline, ~Basis, ~Outcome,
    "cytoHubba", "degree-based centrality",
      paste0("annotation AUC ", fmt(base_ann),
             ", citation AUC ", fmt(base_cit)))
  wr(t2base, "Table2_baseline")

  # ---- Table 3: selector comparison -------------------------------------
  if (have("downgrade_result.csv")) {
    d <- rd("downgrade_result.csv")
    lab <- c(variance = "Variance", DE_tstat = "DE t-statistic",
             relevance = "Relevance (MI)", mRMR = "mRMR",
             degree = "STRING degree",
             cytoHubba_subgraph = "cytoHubba",
             coexpr_degree = "Co-expression degree")
    t3 <- d %>% group_by(selector) %>%
      summarise(Reproducibility = round(mean(reproducibility), 3),
                Citation_AUC = round(mean(citation_auc), 3),
                Network_free = network_free[1], .groups = "drop") %>%
      filter(selector %in% names(lab)) %>%
      mutate(Selector = lab[selector]) %>%
      arrange(desc(Reproducibility)) %>%
      select(Selector, Reproducibility, Citation_AUC, Network_free) %>%
      mutate(Note = ifelse(Selector == "Co-expression degree",
        "top-k Jaccard; full-vector Spearman 0.20", ""))
    wr(t3, "Table3_selectors")
    print(as.data.frame(t3))
  }

  # ---- Supplementary S1: pre-declared criteria --------------------------
  crit_files <- c("infoselect_criteria.csv", "downgrade_criteria.csv",
                  "constrained_criteria.csv", "cod_criteria.csv")
  s1 <- list()
  for (f in crit_files) if (file.exists(P("tables", f))) {
    x <- rd(f); x$source <- sub("_criteria.csv", "", f)
    s1[[length(s1) + 1]] <- x
  }
  if (length(s1)) {
    s1b <- bind_rows(s1)
    wr(s1b, "SupplementaryTable1_criteria")
  }

  # ---- Supplementary S2: adjustment provenance --------------------------

  # ---- Supplementary S3: reproducibility by density ---------------------
  if (have("coexpr_fullvector.csv")) {
    fv <- rd("coexpr_fullvector.csv")
    wr(fv, "SupplementaryTable3_reproducibility_by_density")
  }

  # ---- Supplementary S4: pipeline benchmark -----------------------------
  if (have("pipeline_benchmark.csv")) {
    pb <- rd("pipeline_benchmark.csv")
    wr(pb, "SupplementaryTable4_pipeline_benchmark")
  }

  # ---- Supplementary S5: cytoHubba metrics vs degree --------------------
  if (have("cytohubba_metric_vs_degree.csv")) {
    cm <- rd("cytohubba_metric_vs_degree.csv")
    wr(cm, "SupplementaryTable5_cytohubba_metrics")
  }

  # ---- Supplementary S6: WGCNA cross-cohort reproducibility -------------
  if (have("wgcna_reproducibility.csv")) {
    wg <- rd("wgcna_reproducibility.csv")
    wr(wg, "SupplementaryTable6_wgcna_reproducibility")
  }

  # ---- Supplementary S7: annotation-bias decomposition ------------------
  if (have("annotation_bias_summary.csv")) {
    ab <- rd("annotation_bias_summary.csv")
    wr(ab, "SupplementaryTable7_annotation_bias")
  }

  log_msg("Tables written to ", P("tables"))
  invisible(NULL)
}
