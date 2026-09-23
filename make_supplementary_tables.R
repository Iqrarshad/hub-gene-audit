#!/usr/bin/env Rscript
# make_supplementary_tables.R
# -----------------------------------------------------------------------------
# Consolidate the pipeline's analysis outputs into the manuscript's
# Supplementary Table numbering (S1-S10). Run AFTER run_all.R has completed a
# full run, e.g.:
#
#   Rscript make_supplementary_tables.R --results-dir "D:\Method Paper\...\results_final"
#
# This ONLY copies already-validated CSVs into correctly-named files under
# <results-dir>/supplementary_tables/. It recomputes nothing, so it is safe to
# rerun. If the final manuscript renumbers a supplementary table, edit MAP below
# and rerun; the mapping is the single source of truth for the S1-S10 labels.
#
# NOTE: confirm this mapping against the final manuscript before Zenodo upload.
# The S4-S9 rows are unambiguous; S1/S2/S3/S10 depend on the final supplementary
# structure and are the ones to double-check.
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
rd <- Sys.getenv("GLIOMA_RESULTS_DIR", "")
i  <- which(args == "--results-dir")
if (length(i) && length(args) >= i + 1) rd <- args[i + 1]
if (rd == "") stop("Set --results-dir <path> or the GLIOMA_RESULTS_DIR env var.")

src       <- file.path(rd, "tables")
repo_root <- getwd()                       # published_hub_lists.csv lives here
out       <- file.path(rd, "supplementary_tables")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# manuscript label  <-  canonical source CSV
MAP <- list(
  "S1_cohort_characteristics"        = "Table1_datasets.csv",
  "S2_cross_source_replicated_genes" = "rgselect_genes.csv",
  "S3_published_hub_gene_studies"    = "published_hub_lists.csv",
  "S4_cytohubba_metrics_vs_degree"   = "cytohubba_metric_vs_degree.csv",
  "S5_edge_reproducibility"          = "edge_reproducibility.csv",
  "S6_coexpression_vs_string"        = "coexpression_vs_string.csv",
  "S7_wgcna_reproducibility"         = "wgcna_reproducibility.csv",
  "S8_annotation_enrichment_bias"    = "annotation_bias_enrichment.csv",
  "S9_string_threshold_sensitivity"  = "density_threshold_sweep.csv",
  "S10_leave_one_cohort_out"         = "rgselect_loco.csv")

done <- character(0); miss <- character(0)
for (nm in names(MAP)) {
  f    <- MAP[[nm]]
  cand <- c(file.path(src, f), file.path(repo_root, f))
  hit  <- cand[file.exists(cand)]
  dst  <- file.path(out, paste0("SupplementaryTable_", nm, ".csv"))
  if (length(hit)) {
    file.copy(hit[1], dst, overwrite = TRUE)
    done <- c(done, sprintf("%-34s <- %s", nm, f))
  } else {
    miss <- c(miss, sprintf("%-34s <- %s (NOT FOUND)", nm, f))
  }
}

cat(sprintf("Wrote %d supplementary tables to %s\n", length(done), out))
for (d in done) cat("  ", d, "\n")
if (length(miss)) {
  cat("\nMISSING - run the full pipeline first, or fix MAP:\n")
  for (m in miss) cat("  ", m, "\n")
}
