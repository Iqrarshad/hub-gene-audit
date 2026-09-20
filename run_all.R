# run_all.R ===============================================================
# Master runner.
#
#   Rscript run_all.R              every stage in order
#   Rscript run_all.R 07 08 09     selected stages only
#   Rscript run_all.R --list       print the stage table and exit
#   Rscript run_all.R --group core run a named group
#   Rscript run_all.R --clean      discard cached results, then run all
#
# Paths are not hardcoded. Set them with:
#   Rscript run_all.R --data-dir /path/to/data --results-dir /path/to/out
# or with GLIOMA_DATA_DIR and GLIOMA_RESULTS_DIR in the environment, or in
# config.local.R. See R/00_config.R for the resolution order.
#
# Stages are independent given their inputs, so a failed stage can be fixed
# and rerun without redoing the ones before it. Intermediate objects are
# cached as .rds under RESULTS_DIR/rds.
#
# STAGE GROUPS
#   core       01-06        data, DEG, enrichment, STRING network, hubs
#   validate   07-10        hub null models and centrality diagnostics
#   bias       11-19,36,39   hub-bias quantification and network diagnostics
#   methods    20-23        attempted corrections, all reported as failures
#   lit        24-28,37     literature-scale analysis and the propagation gates
#   propagate  30-31,33-35   propagation audit, information selection, baselines
#   figures    29,32,38,40   manuscript figures and tables
#   recommend  41           reproducibility-gated selection and audit card
#   intervals  42-44        BCa intervals, multi-seed robustness, corrections table
#
# Run 'Rscript run_all.R --list' for the authoritative stage table.
#
# =========================================================================

args <- commandArgs(trailingOnly = TRUE)

if (!file.exists("R/00_config.R"))
  stop("Run this from the project root, the folder containing R/.")

# Path flags are consumed here and passed to the config as environment
# variables, so the remaining arguments are stage selectors only.
take_flag <- function(a, name) {
  i <- which(a == name)
  if (!length(i)) return(list(args = a, value = NA_character_))
  if (i[1] == length(a)) stop(name, " requires a path")
  list(args = a[-(i[1]:(i[1] + 1))], value = a[i[1] + 1])
}
for (fl in c("--data-dir", "--results-dir")) {
  got <- take_flag(args, fl)
  args <- got$args
  if (!is.na(got$value)) {
    if (fl == "--data-dir") Sys.setenv(GLIOMA_DATA_DIR = got$value)
    else Sys.setenv(GLIOMA_RESULTS_DIR = got$value)
  }
}

STAGES <- list(
  "01" = list(file = "R/01_download.R", fn = "main_01", grp = "core",
              desc = "Load data, build cohort manifest"),
  "02" = list(file = "R/02_deg_microarray.R", fn = "main_02", grp = "core",
              desc = "Microarray DEG, composition adjusted"),
  "03" = list(file = "R/03_deg_rnaseq.R", fn = "main_03", grp = "core",
              desc = "RNA-seq DEG, three contrasts"),
  "04" = list(file = "R/04_intersect.R", fn = "main_04", grp = "core",
              desc = "DEG intersection and Jaccard"),
  "05" = list(file = "R/05_enrichment.R", fn = "main_05", grp = "core",
              desc = "KEGG and GO enrichment, BH corrected"),
  "06" = list(file = "R/06_network.R", fn = "main_06", grp = "core",
              desc = "STRING network and consensus hubs"),

  "07" = list(file = "R/07_hub_null.R", fn = "main_07", grp = "validate",
              desc = "Hub spike-in null"),
  "08" = list(file = "R/08_centrality_diagnostics.R", fn = "main_08", grp = "validate",
              desc = "Degree normalisation, physical network, redundancy"),
  "09" = list(file = "R/09_connectivity_null_and_hla.R", fn = "main_09", grp = "validate",
              desc = "Degree-enrichment null and HLA direction"),
  "10" = list(file = "R/10_specific_connectivity.R", fn = "main_10", grp = "validate",
              desc = "Binomial test for specific connectivity"),


  "11" = list(file = "R/11_surfaceome.R", fn = "main_11", grp = "bias",
              desc = "Surfaceome composition of hub nominations"),
  "12" = list(file = "R/12_hub_bias_quantified.R", fn = "main_12", grp = "bias",
              desc = "Hub bias quantified, five tests"),
  "13" = list(file = "R/13_rewiring_zscore.R", fn = "main_13", grp = "bias",
              desc = "Rewiring z-score, with degree floor"),
  "14" = list(file = "R/14_density_threshold.R", fn = "main_14", grp = "bias",
              desc = "Network density and hub stability"),
  "15" = list(file = "R/15_inference_level_stability.R", fn = "main_15", grp = "bias",
              desc = "Stability by level of inference"),
  "16" = list(file = "R/16_stability_ceiling.R", fn = "main_16", grp = "bias",
              desc = "Oracle ceiling for stability metrics"),
  "17" = list(file = "R/17_coexpression_network.R", fn = "main_17", grp = "bias",
              desc = "Cohort co-expression network versus STRING"),
  "18" = list(file = "R/18_grade_matched_coexpression.R", fn = "main_18", grp = "bias",
              desc = "Grade-matched co-expression comparison"),
  "19" = list(file = "R/19_edge_reproducibility.R", fn = "main_19", grp = "bias",
              desc = "Edge-level reproducibility across cohorts"),

  "20" = list(file = "R/20_composite_hub_score.R", fn = "main_20", grp = "methods",
              desc = "Composite hub score, four combination rules"),
  "21" = list(file = "R/21_reproducibility_weighted_eigenvector.R", fn = "main_21",
              grp = "methods",
              desc = "Eigenvector on reproducibility-weighted network"),
  "22" = list(file = "R/22_pipeline_benchmark.R", fn = "main_22", grp = "methods",
              desc = "End-to-end pipeline comparison"),
  "23" = list(file = "R/23_depmap_validation.R", fn = "main_23", grp = "methods",
              desc = "DOS viability and DepMap CRISPR validation"),

  "27" = list(file = "R/27_propagation_gates.R", fn = "main_27", grp = "lit",
              desc = "Propagation gates, transcriptome scale"),
  "24" = list(file = "R/24_literature_bias_gene2pubmed.R", fn = "main_24", grp = "lit",
              desc = "Literature bias measured via gene2pubmed"),
  "25" = list(file = "R/25_citation_penalised_score.R", fn = "main_25", grp = "lit",
              desc = "Citation-penalised hub score"),
  "26" = list(file = "R/26_coessentiality_network.R", fn = "main_26", grp = "lit",
              desc = "CRISPR co-essentiality network substrate"),
  "28" = list(file = "R/28_literature_scale_analysis.R", fn = "main_28", grp = "lit",
              desc = "Literature analysis, 30 studies"),
  "29" = list(file = "R/29_manuscript_figures.R", fn = "main_29", grp = "figures",
              desc = "Publication figures for the principal findings"),

  "30" = list(file = "R/30_propagation_audit.R", fn = "main_30", grp = "propagate",
              desc = "Why propagation recovers degree, with calibration"),
  "31" = list(file = "R/31_constrained_diffusion.R", fn = "main_31", grp = "propagate",
              desc = "Constrained diffusion against post-hoc adjustment"),
  "32" = list(file = "R/32_propagation_figures.R", fn = "main_32", grp = "figures",
              desc = "Figures for the propagation audit"),

  "33" = list(file = "R/33_information_selection.R", fn = "main_33", grp = "propagate",
              desc = "Information-theoretic gene selection, prior-check"),

  "34" = list(file = "R/34_downgrade_keystone.R", fn = "main_34", grp = "propagate",
              desc = "Does the network beat a network-free baseline"),
  "35" = list(file = "R/35_reproducibility_fullvector.R", fn = "main_35", grp = "propagate",
              desc = "Full-vector co-expression degree reproducibility"),
  "36" = list(file = "R/36_wgcna_reproducibility.R", fn = "main_36", grp = "bias",
              desc = "WGCNA per-cohort networks, cross-cohort reproducibility"),
  "37" = list(file = "R/37_annotation_bias.R", fn = "main_37", grp = "lit",
              desc = "Annotation bias, the enrichment analogue of degree bias"),

  "38" = list(file = "R/38_manuscript_figures_main.R", fn = "main_38", grp = "figures",
              desc = "Main manuscript figures for the audit"),
  "39" = list(file = "R/39_cytohubba_completeness.R", fn = "main_39", grp = "bias",
              desc = "All 11 cytoHubba metrics versus degree"),
  "40" = list(file = "R/40_tables.R", fn = "main_40", grp = "figures",
              desc = "Publication tables, main and supplementary"),
  "41" = list(file = "R/41_reproducibility_gated_selection.R", fn = "main_41",
              grp = "recommend",
              desc = "Reproducibility-gated selection and the audit card"),

  "42" = list(file = "R/42_bootstrap_intervals.R", fn = "main_42",
              grp = "intervals",
              desc = "BCa 5000-resample intervals for the headline AUCs"),
  "43" = list(file = "R/43_multiseed_robustness.R", fn = "main_43",
              grp = "intervals",
              desc = "Multi-seed robustness of the bootstrap intervals"),
  "44" = list(file = "R/44_corrections_table.R", fn = "main_44",
              grp = "intervals",
              desc = "Consolidated corrections table incl. two propagation variants")
)

if (length(args) && args[1] == "--list") {
  cat(sprintf("%-4s %-9s %s\n", "id", "group", "description"))
  for (s in names(STAGES))
    cat(sprintf("%-4s %-9s %s\n", s, STAGES[[s]]$grp, STAGES[[s]]$desc))
  quit(save = "no")
}

source("R/00_config.R")

# --clean removes derived results only. Source data under DATA_DIR is left
# alone so the STRING and GEO downloads are not repeated.
if (length(args) && args[1] == "--clean") {
  for (d in c("rds", "tables", "figures")) {
    p <- P(d)
    if (dir.exists(p)) {
      unlink(list.files(p, full.names = TRUE, recursive = TRUE),
             recursive = TRUE)
      log_msg("Cleared ", p)
    }
    dir.create(p, showWarnings = FALSE, recursive = TRUE)
  }
  args <- args[-1]
}

if (length(args) >= 2 && args[1] == "--group") {
  to_run <- names(STAGES)[vapply(STAGES, function(x) x$grp, "") == args[2]]
  if (!length(to_run)) stop("Unknown group: ", args[2])
} else if (length(args) == 0) {
  to_run <- names(STAGES)
} else {
  to_run <- args
  unknown <- setdiff(to_run, names(STAGES))
  if (length(unknown))
    stop("Unknown stage(s): ", paste(unknown, collapse = ", "),
         "\nRun 'Rscript run_all.R --list' to see available stages.")
}

# Helper functions defined in these files are needed by several stages
ALWAYS_SOURCE <- c("R/01_download.R")
for (f in ALWAYS_SOURCE) source(f)

log_msg("=== Run started. Stages: ", paste(to_run, collapse = ", "), " ===")
t0 <- Sys.time(); failed <- character(0)

for (s in to_run) {
  st <- STAGES[[s]]
  log_msg("--- Stage ", s, " [", st$grp, "]: ", st$desc, " ---")
  ts <- Sys.time()
  ok <- tryCatch({ source(st$file); do.call(st$fn, list()); TRUE },
                 error = function(e) {
                   log_msg("STAGE ", s, " FAILED: ", conditionMessage(e))
                   failed <<- c(failed, s); FALSE })
  log_msg("Stage ", s, if (ok) " ok" else " failed", " in ",
          round(difftime(Sys.time(), ts, units = "mins"), 2), " min")
}

log_msg("=== Run finished in ",
        round(difftime(Sys.time(), t0, units = "mins"), 2), " min ===")
if (length(failed)) {
  log_msg("Failed stages: ", paste(failed, collapse = ", "))
} else {
  log_msg("All stages completed.")
}

writeLines(capture.output(sessionInfo()), P("logs", "sessionInfo.txt"))
