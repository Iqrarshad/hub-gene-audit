# Multi-seed robustness of the bootstrap intervals.
#
# The AUC point estimates are deterministic given the data, so they do not
# move with the random seed. What can move is the resampled interval. This
# stage recomputes every interval from stage 42 under a set of seeds and
# reports how stable the interval endpoints are, so the reported CIs can be
# shown to be a property of the data rather than of one seed.
#
# It reuses collect_auc_inputs() from stage 42, so the labels and predictors
# are identical to the single-seed run. It is fully offline and reads only
# the cached tables.
#
# Scope note. This covers seed stability of the interval estimates. Seed
# stability of the STRING spike-in null (stage 07) is a separate question,
# because re-drawing its random backgrounds calls the STRING API; run stage
# 07 under different values of SEED for that, with an API budget in mind. It
# is deliberately not triggered here.
#
# Outputs:
#   tables/multiseed_robustness.csv        endpoint stability per quantity
#   tables/multiseed_robustness_raw.csv    every seed's interval

source("R/00h_bootstrap.R")
source("R/42_bootstrap_intervals.R")   # for collect_auc_inputs()

MS_SEEDS <- 20L      # number of seeds
MS_R     <- 2000L    # resamples per seed (lower than stage 42 for runtime)

main_43 <- function() {
  recs <- collect_auc_inputs()
  if (!length(recs)) {
    log_msg("No AUC inputs found; multi-seed robustness skipped.")
    return(invisible(NULL))
  }

  seeds <- SEED + seq_len(MS_SEEDS) - 1L
  log_msg("=================================================")
  log_msg("MULTI-SEED ROBUSTNESS: ", MS_SEEDS, " seeds x ", MS_R,
          " resamples")
  log_msg("=================================================")

  raw <- list(); summ <- list()
  for (r in recs) {
    est  <- fast_auc(r$y, r$x)
    los  <- numeric(0); his <- numeric(0); meths <- character(0)
    for (s in seeds) {
      ci <- bca_auc_ci(r$y, r$x, R = MS_R, seed = s)
      los  <- c(los, unname(ci["lo"]))
      his  <- c(his, unname(ci["hi"]))
      meths <- c(meths, attr(ci, "method"))
      raw[[length(raw) + 1]] <- data.frame(
        analysis = r$analysis, quantity = r$quantity, hub_set = r$hub_set,
        seed = s, estimate = round(est, 4),
        ci_lo = round(unname(ci["lo"]), 4),
        ci_hi = round(unname(ci["hi"]), 4),
        method = attr(ci, "method"), stringsAsFactors = FALSE)
    }
    summ[[length(summ) + 1]] <- data.frame(
      analysis = r$analysis, quantity = r$quantity, hub_set = r$hub_set,
      estimate = round(est, 4),
      lo_mean = round(mean(los, na.rm = TRUE), 4),
      lo_sd   = round(stats::sd(los, na.rm = TRUE), 4),
      lo_min  = round(min(los, na.rm = TRUE), 4),
      lo_max  = round(max(los, na.rm = TRUE), 4),
      hi_mean = round(mean(his, na.rm = TRUE), 4),
      hi_sd   = round(stats::sd(his, na.rm = TRUE), 4),
      hi_min  = round(min(his, na.rm = TRUE), 4),
      hi_max  = round(max(his, na.rm = TRUE), 4),
      n_seeds = length(seeds), resamples_per_seed = MS_R,
      bca_fraction = round(mean(meths == "bca"), 3),
      stringsAsFactors = FALSE)
    log_msg("  ", r$quantity, " (", r$hub_set, "): est ", round(est, 3),
            "; lo ", round(mean(los, na.rm = TRUE), 3), " +/- ",
            round(stats::sd(los, na.rm = TRUE), 3), "; hi ",
            round(mean(his, na.rm = TRUE), 3), " +/- ",
            round(stats::sd(his, na.rm = TRUE), 3))
  }

  s_out <- do.call(rbind, summ)
  r_out <- do.call(rbind, raw)
  utils::write.csv(s_out, P("tables", "multiseed_robustness.csv"),
                   row.names = FALSE)
  utils::write.csv(r_out, P("tables", "multiseed_robustness_raw.csv"),
                   row.names = FALSE)
  log_msg("Wrote multiseed_robustness.csv (", nrow(s_out),
          " quantities) and multiseed_robustness_raw.csv (", nrow(r_out),
          " rows)")
  print(s_out[, c("quantity", "estimate", "lo_mean", "lo_sd", "hi_mean",
                  "hi_sd")])
  log_msg("43_multiseed_robustness complete.")
  invisible(s_out)
}
