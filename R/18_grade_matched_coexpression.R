# Is cross-cohort disagreement biological or is it instability?
#
# CGGA is restricted to LGG only and rebuilt. If the LGG-only networks
# converge toward the pure-LGG TCGA cohort, the disagreement was grade
# composition. A within-grade HGG comparison runs alongside as a control.
#
# If the CGGA LGG-only networks agree with each other but not with TCGA,
# that indicates a platform or batch effect rather than instability, and
# both comparisons are reported so the distinction is visible.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")
source("R/17_coexpression_network.R")

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (!u) return(NA_real_)
  length(intersect(a, b)) / u
}

main_18 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene
  net   <- readRDS(P("rds", "network.rds"))
  score <- as.integer(net$chosen)

  sf <- file.path(CACHE_DIR, "string_density", paste0("t", score, ".rds"))
  se <- if (file.exists(sf)) readRDS(sf) else
        readRDS(P("rds", paste0("string_", score, ".rds")))
  n_edges <- se %>% filter(from %in% genes, to %in% genes) %>% nrow()
  log_msg("Edge budget matched to STRING: ", n_edges)

  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  global_deg <- setNames(norm$global_degree, norm$gene)

  cohorts <- get_validation_cohorts()

  # --- Build one network per cohort x grade stratum -----------------------
  strata <- list()
  for (nm in names(cohorts)) {
    co <- cohorts[[nm]]
    mat <- co$matrix; meta <- co$meta

    combos <- list(all = meta$sample)
    for (gr in intersect(c("LGG", "HGG"), unique(meta$grade))) {
      ids <- meta$sample[meta$grade == gr]
      if (length(ids) >= 100) combos[[gr]] <- ids
    }

    for (lab in names(combos)) {
      ids <- intersect(combos[[lab]], colnames(mat))
      if (length(ids) < 100) {
        log_msg("Skipping ", nm, " ", lab, ": n = ", length(ids))
        next
      }
      key <- paste0(nm, " [", lab, "]")
      cx <- tryCatch(build_coexpression(mat[, ids, drop = FALSE], genes,
                                        n_edges, key),
                     error = function(e) { log_msg(key, ": ",
                                                   conditionMessage(e)); NULL })
      if (is.null(cx)) next
      h <- hub_select(cx$graph)
      strata[[key]] <- list(cohort = nm, stratum = lab, n = length(ids),
                            hubs = h$hubs, graph = cx$graph,
                            edges = cx$edges)
      log_msg(key, " (n = ", length(ids), "): ",
              paste(head(h$hubs, 10), collapse = ", "))
    }
  }

  if (length(strata) < 3) stop("Too few strata built.")

  # --- Pairwise hub agreement --------------------------------------------
  keys <- names(strata)
  pairs <- combn(keys, 2, simplify = FALSE)
  agr <- map_dfr(pairs, function(p) {
    a <- strata[[p[1]]]; b <- strata[[p[2]]]
    tibble(set_a = p[1], set_b = p[2],
           stratum_a = a$stratum, stratum_b = b$stratum,
           same_stratum = a$stratum == b$stratum,
           same_cohort_family = grepl("CGGA", a$cohort) &&
                                grepl("CGGA", b$cohort),
           n_a = a$n, n_b = b$n,
           hub_jaccard = jaccard(a$hubs, b$hubs),
           shared = paste(intersect(a$hubs, b$hubs), collapse = "; "))
  })
  write_csv(agr, P("tables", "grade_matched_hub_agreement.csv"))

  log_msg("=================================================")
  log_msg("HUB AGREEMENT BETWEEN COHORT-STRATUM NETWORKS")
  log_msg("=================================================")
  print(as.data.frame(agr %>%
    select(set_a, set_b, same_stratum, hub_jaccard, shared) %>%
    mutate(hub_jaccard = round(hub_jaccard, 3))))

  # --- The question ------------------------------------------------------
  lgg_pairs <- agr %>% filter(stratum_a == "LGG", stratum_b == "LGG")
  mix_pairs <- agr %>% filter(stratum_a == "all", stratum_b == "all")
  hgg_pairs <- agr %>% filter(stratum_a == "HGG", stratum_b == "HGG")
  cross     <- agr %>% filter(!same_stratum)

  summ <- tibble(
    comparison = c("LGG vs LGG", "HGG vs HGG", "mixed vs mixed",
                   "different strata"),
    n_pairs = c(nrow(lgg_pairs), nrow(hgg_pairs), nrow(mix_pairs),
                nrow(cross)),
    mean_jaccard = c(mean(lgg_pairs$hub_jaccard, na.rm = TRUE),
                     mean(hgg_pairs$hub_jaccard, na.rm = TRUE),
                     mean(mix_pairs$hub_jaccard, na.rm = TRUE),
                     mean(cross$hub_jaccard, na.rm = TRUE)))
  write_csv(summ, P("tables", "grade_matched_summary.csv"))
  log_msg("Agreement by stratum matching:")
  print(as.data.frame(summ %>% mutate(mean_jaccard = round(mean_jaccard, 3))))

  # Specifically: does LGG-only CGGA move toward TCGA?
  tcga_keys <- keys[grepl("TCGA", keys)]
  cgga_lgg  <- keys[grepl("CGGA", keys) & grepl("\\[LGG\\]", keys)]
  cgga_all  <- keys[grepl("CGGA", keys) & grepl("\\[all\\]", keys)]

  if (length(tcga_keys) && length(cgga_lgg) && length(cgga_all)) {
    tk <- tcga_keys[1]
    j_lgg <- map_dbl(cgga_lgg, ~ jaccard(strata[[.x]]$hubs,
                                         strata[[tk]]$hubs))
    j_all <- map_dbl(cgga_all, ~ jaccard(strata[[.x]]$hubs,
                                         strata[[tk]]$hubs))
    log_msg("---")
    log_msg("Agreement of CGGA with ", tk, ":")
    log_msg("  CGGA mixed grade  : ", paste(round(j_all, 3), collapse = ", "),
            "  (mean ", round(mean(j_all, na.rm = TRUE), 3), ")")
    log_msg("  CGGA LGG only     : ", paste(round(j_lgg, 3), collapse = ", "),
            "  (mean ", round(mean(j_lgg, na.rm = TRUE), 3), ")")

    if (mean(j_lgg, na.rm = TRUE) > mean(j_all, na.rm = TRUE) + 0.1) {
      log_msg("CONCLUSION: restricting CGGA to LGG moves it TOWARD TCGA. ",
              "The disagreement in the co-expression comparison was grade composition, not ",
              "instability. Co-expression structure is grade specific.")
    } else if (mean(j_lgg, na.rm = TRUE) > 0.3) {
      log_msg("CONCLUSION: LGG-only networks agree moderately regardless of ",
              "cohort. Grade explains part of the disagreement.")
    } else {
      log_msg("CONCLUSION: restricting to LGG does NOT improve agreement ",
              "with TCGA. The disagreement is not grade composition. Either ",
              "it is platform or batch, or the method does not reproduce ",
              "across datasets. Report it as a limitation.")
      if (nrow(lgg_pairs) && mean(lgg_pairs$hub_jaccard, na.rm = TRUE) > 0.3) {
        log_msg("  Note: the two CGGA LGG networks DO agree with each other ",
                "(", round(mean(lgg_pairs$hub_jaccard, na.rm = TRUE), 3),
                "), which points to a platform or batch effect between ",
                "CGGA and TCGA rather than to instability.")
      }
    }
  }

  # --- Recurrent genes across all strata ---------------------------------
  allh <- table(unlist(map(strata, "hubs")))
  rec <- tibble(gene = names(allh), n_networks = as.integer(allh),
                pct = round(100 * as.integer(allh) / length(strata), 1)) %>%
    arrange(desc(n_networks))
  write_csv(rec, P("tables", "grade_matched_recurrent_hubs.csv"))
  log_msg("Genes recurring as hubs across the ", length(strata),
          " cohort-stratum networks:")
  print(as.data.frame(head(rec, 15)))

  # --- Annotation bias check, per stratum --------------------------------
  bias <- map_dfr(keys, function(k) {
    s <- strata[[k]]
    tibble(network = k, n = s$n, n_hubs = length(s$hubs),
           auc_global_degree = predictability_auc(s$hubs,
                                V(s$graph)$name, global_deg))
  })
  write_csv(bias, P("tables", "grade_matched_bias.csv"))
  log_msg("Annotation bias per stratum network (0.5 = none):")
  print(as.data.frame(bias %>% mutate(auc_global_degree =
                                        round(auc_global_degree, 3))))

  saveRDS(list(strata = map(strata, ~ .x[c("cohort", "stratum", "n", "hubs")]),
               agreement = agr, summary = summ, recurrent = rec, bias = bias),
          P("rds", "grade_matched.rds"))
  log_msg("18_grade_matched_coexpression complete.")
  invisible(agr)
}
