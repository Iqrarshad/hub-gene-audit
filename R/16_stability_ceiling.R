# Empirical ceiling for split-half stability metrics.
#
# Under 50% gene removal, roughly half the true positives are absent by
# construction, so the maximum achievable Jaccard is well below 1.0:
#
#   full set yields k hubs; a random half retains about k/2 of them; the
#   half-network still nominates about k of its own
#
# An oracle method that always returns the correct answer is simulated
# under the identical procedure, giving the ceiling. A random method gives
# the floor. Observed values are then rescaled onto that range.
#
# Interpreting a raw Jaccard against 1.0 rather than against the ceiling is
# a design error.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr); library(tidyr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

N_SPLITS <- 500

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (!u) return(NA_real_)
  length(intersect(a, b)) / u
}

# --- Oracle -------------------------------------------------------------
# Knows the truth. On any gene set it returns the true hubs present in that
# set. Cannot do better; this is the ceiling.
oracle_stability <- function(n_genes, k_hubs, frac = 0.5,
                             n_splits = N_SPLITS) {
  genes <- paste0("G", seq_len(n_genes))
  truth <- genes[seq_len(k_hubs)]
  full <- truth                       # oracle on the full set
  js <- map_dbl(seq_len(n_splits), function(i) {
    set.seed(SEED + i)
    sub <- sample(genes, floor(n_genes * frac))
    jaccard(intersect(truth, sub), full)
  })
  mean(js, na.rm = TRUE)
}

# --- Random floor -------------------------------------------------------
# Nominates k genes at random from whatever is present.
random_stability <- function(n_genes, k_hubs, frac = 0.5,
                             n_splits = N_SPLITS) {
  genes <- paste0("G", seq_len(n_genes))
  js <- map_dbl(seq_len(n_splits), function(i) {
    set.seed(SEED + 10000 + i)
    full <- sample(genes, k_hubs)
    sub  <- sample(genes, floor(n_genes * frac))
    half <- sample(sub, min(k_hubs, length(sub)))
    jaccard(half, full)
  })
  mean(js, na.rm = TRUE)
}

# --- A realistic intermediate -------------------------------------------
# A method that is right most of the time: recovers the true hubs present,
# but misses a proportion and adds a proportion of false ones.
noisy_stability <- function(n_genes, k_hubs, miss = 0.2, false_rate = 0.2,
                            frac = 0.5, n_splits = N_SPLITS) {
  genes <- paste0("G", seq_len(n_genes))
  truth <- genes[seq_len(k_hubs)]
  call_on <- function(sub, seed) {
    set.seed(seed)
    present <- intersect(truth, sub)
    keep <- present[runif(length(present)) > miss]
    n_false <- rbinom(1, length(present), false_rate)
    extras <- if (n_false > 0)
      sample(setdiff(sub, truth), min(n_false, length(setdiff(sub, truth))))
      else character(0)
    c(keep, extras)
  }
  full <- call_on(genes, SEED)
  js <- map_dbl(seq_len(n_splits), function(i) {
    set.seed(SEED + 20000 + i)
    sub <- sample(genes, floor(n_genes * frac))
    jaccard(call_on(sub, SEED + 30000 + i), full)
  })
  mean(js, na.rm = TRUE)
}

main_16 <- function() {
  # Match the observed design
  n_genes <- 148
  k_hubs  <- 8

  log_msg("Design: ", n_genes, " genes, ", k_hubs,
          " hubs, 50 percent removal, ", N_SPLITS, " splits")

  # --- Ceiling across a range of designs ---------------------------------
  grid <- expand_grid(n_genes = c(148, 300, 600),
                      k_hubs = c(5, 8, 10, 20),
                      frac = c(0.5, 0.7, 0.9))

  ceil <- pmap_dfr(grid, function(n_genes, k_hubs, frac) {
    tibble(n_genes = n_genes, k_hubs = k_hubs, retained_frac = frac,
           oracle = oracle_stability(n_genes, k_hubs, frac, 200),
           random = random_stability(n_genes, k_hubs, frac, 200))
  })
  write_csv(ceil, P("tables", "stability_ceiling_grid.csv"))

  log_msg("=== CEILING BY DESIGN ===")
  print(as.data.frame(ceil %>% mutate(across(where(is.numeric),
                                             ~ round(.x, 3)))))

  # --- The specific design used -----------------------------------------
  orc <- oracle_stability(n_genes, k_hubs, 0.5, N_SPLITS)
  rnd <- random_stability(n_genes, k_hubs, 0.5, N_SPLITS)
  noisy <- noisy_stability(n_genes, k_hubs, 0.2, 0.2, 0.5, N_SPLITS)

  log_msg("=================================================")
  log_msg("EMPIRICAL CEILING FOR THE OBSERVED DESIGN")
  log_msg("=================================================")
  log_msg("  oracle  (always correct)      : ", round(orc, 3))
  log_msg("  noisy   (20% miss, 20% false) : ", round(noisy, 3))
  log_msg("  random  (chance)              : ", round(rnd, 3))

  # --- Rescale the observed values --------------------------------------
  f <- P("tables", "inference_level_stability.csv")
  obs <- if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL

  if (!is.null(obs)) {
    obs <- obs %>%
      mutate(oracle_ceiling = orc, random_floor = rnd,
             normalised = (stability - rnd) / (orc - rnd),
             pct_of_ceiling = round(100 * stability / orc, 1))
    write_csv(obs, P("tables", "stability_normalised.csv"))

    log_msg("=== OBSERVED VALUES RESCALED AGAINST THE CEILING ===")
    print(as.data.frame(obs %>%
      select(level, stability, oracle_ceiling, random_floor,
             normalised, pct_of_ceiling) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))))
  }

  # --- Verdict -----------------------------------------------------------
  hub_obs <- if (!is.null(obs))
    obs$stability[grepl("gene-level", obs$level)][1] else NA_real_

  if (!is.na(hub_obs)) {
    log_msg("Hub selection observed ", round(hub_obs, 3),
            " against a ceiling of ", round(orc, 3), " and a floor of ",
            round(rnd, 3))
    if (hub_obs > 0.9 * orc) {
      log_msg("VERDICT: hub selection performs at essentially the maximum ",
              "this design permits. The instability conclusion from stages ",
              "26 and 27 was an artefact of comparing against 1.0 rather ",
              "than against the achievable ceiling. It must be withdrawn.")
      log_msg("The BIAS findings are unaffected: they rest on AUC 0.994 for ",
              "global degree predicting hub status, and on degree-preserving ",
              "rewiring, neither of which uses split-half stability.")
    } else if (hub_obs > 0.6 * orc) {
      log_msg("VERDICT: hub selection reaches ", round(100 * hub_obs / orc),
              " percent of the achievable ceiling. Performance is moderate ",
              "rather than poor. State the ceiling alongside the observed ",
              "value.")
    } else {
      log_msg("VERDICT: hub selection reaches only ",
              round(100 * hub_obs / orc), " percent of the ceiling. The ",
              "instability finding stands even after correcting for the ",
              "design limit.")
    }
  }

  saveRDS(list(ceiling_grid = ceil, oracle = orc, random = rnd,
               noisy = noisy, observed = obs),
          P("rds", "stability_ceiling.rds"))
  log_msg("16_stability_ceiling complete.")
  invisible(ceil)
}
