# Literature-scale analysis: 30 published hub gene lists across 21 cancer
# descriptions.
#
# Each observed hub gene is replaced by a background gene of similar
# degree, matched within bins, so each study keeps its own degree profile
# and no shared pool is imposed. Recurrence then arises only if
# degree-matched substitutes coincide across studies.
#
# Drawing each study's list from a shared pool does not work: that null
# depends only on pool size, not on which genes are in the pool.
#
# Gates, declared before the result:
#   L1  resolved background >= 1500
#   L2  every degree bin >= 100 genes
#   L3  matched substitutes within 15% of the observed median degree
#   L4  random null mean below the observed count, else the pool is
#       exhausted and that null is uninformative
#
# Both p_more and p_fewer are reported, since degree may over-predict as
# well as under-predict the observed overlap.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(httr); library(pROC)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; distinct <- dplyr::distinct
summarise <- dplyr::summarise; arrange <- dplyr::arrange

BACKGROUND_N  <- 3000     # requested; expect ~63% to resolve in STRING
MIN_BG        <- 1500     # on RESOLVED genes
MIN_BIN       <- 100
MAX_MATCH_ERR <- 0.15
N_NULL        <- 5000
N_BINS        <- 10
API_DELAY     <- 0.30

safe_key <- function(x) substr(gsub("[^A-Za-z0-9._-]", "_", x), 1, 100)

global_degree <- function(gene, score = 400) {
  cd <- file.path(CACHE_DIR, "string_global")
  dir.create(cd, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cd, paste0(safe_key(gene), "_", score, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- tryCatch(httr::GET("https://string-db.org/api/tsv/interaction_partners",
        query = list(identifiers = gene, species = 9606,
                     required_score = score, limit = 10000,
                     caller_identity = "glioma_lit30v3")),
        error = function(e) NULL)
  n <- NA_integer_
  if (!is.null(r) && httr::status_code(r) == 200) {
    d <- tryCatch(readr::read_tsv(I(httr::content(r, as = "text",
                                                  encoding = "UTF-8")),
                                  show_col_types = FALSE),
                  error = function(e) NULL)
    if (!is.null(d)) n <- nrow(d)
  }
  saveRDS(n, f); Sys.sleep(API_DELAY); n
}

fetch_degrees <- function(genes, label) {
  log_msg("Querying degree for ", length(genes), " ", label,
          " genes (cached entries are instant)")
  t0 <- Sys.time()
  map_dfr(seq_along(genes), function(i) {
    if (i %% 200 == 0)
      log_msg("  ", i, "/", length(genes), "  (",
              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
              " min)")
    v <- tryCatch(global_degree(genes[i]), error = function(e) NA_integer_)
    tibble(gene = genes[i], degree = v)
  }) %>% filter(!is.na(degree), degree > 0)
}

build_background <- function(n_target = BACKGROUND_N) {
  cache <- P("rds", paste0("background_degrees_", n_target, ".rds"))
  if (file.exists(cache)) {
    d <- readRDS(cache)
    log_msg("Background from cache: ", nrow(d), " resolved genes")
    return(d)
  }
  micro <- readRDS(P("rds", "deg_microarray.rds"))
  rna   <- readRDS(P("rds", "deg_rnaseq.rds"))
  uni <- unique(c(unlist(lapply(micro, function(x)
                    if (!is.null(x)) x$gene)), rna$LGG_vs_Normal$gene))
  uni <- uni[!is.na(uni) & uni != ""]
  uni <- uni[!grepl("///|\\s", uni)]
  uni <- uni[grepl("^[A-Za-z][A-Za-z0-9._-]{0,20}$", uni)]
  uni <- uni[!grepl("^(AK|AL|AF|BC|NM_|NR_|XM_|XR_|ENSG)[0-9]", uni)]
  log_msg("Transcriptome universe after cleaning: ", length(uni))

  set.seed(SEED)
  samp <- sample(uni, min(n_target, length(uni)))
  d <- fetch_degrees(samp, "background")
  log_msg("Requested ", length(samp), "; resolved in STRING ", nrow(d),
          " (", round(100 * nrow(d) / length(samp), 1), "%)")
  saveRDS(d, cache)
  d
}

main_28 <- function() {
  # The published hub lists are a supporting file shipped with the repo,
  # not a generated table. Look in the project root first, then tables.
  src <- if (file.exists("published_hub_lists.csv"))
    "published_hub_lists.csv" else P("tables", "published_hub_lists.csv")
  if (!file.exists(src))
    stop("published_hub_lists.csv not found in the project root or tables. ",
         "It ships with the repository.")
  lit <- read_csv(src, show_col_types = FALSE)
  long <- lit %>%
    mutate(gene = strsplit(hub_genes, ";")) %>% unnest(gene) %>%
    mutate(gene = toupper(trimws(gene))) %>%
    filter(gene != "", !is.na(gene), !grepl("///|\\s", gene),
           grepl("^[A-Z][A-Z0-9._-]{0,20}$", gene))

  log_msg("Studies: ", n_distinct(long$study_id),
          "; gene-study pairs: ", nrow(long),
          "; unique genes: ", n_distinct(long$gene),
          "; cancer descriptions: ", n_distinct(long$cancer))

  bg  <- build_background()
  hub <- fetch_degrees(unique(long$gene), "hub")
  log_msg("Background resolved: ", nrow(bg),
          "; hub genes resolved: ", nrow(hub), " of ",
          n_distinct(long$gene))
  log_msg("Median degree -- background ", median(bg$degree),
          ", hub genes ", median(hub$degree))

  L1 <- nrow(bg) >= MIN_BG
  log_msg("GATE L1 resolved background >= ", MIN_BG, ": ",
          ifelse(L1, "PASS", "FAIL"))

  brks <- unique(quantile(bg$degree, seq(0, 1, length.out = N_BINS + 1),
                          na.rm = TRUE))
  bin_of <- function(d) as.integer(cut(d, brks, include.lowest = TRUE))
  bg$bin <- bin_of(bg$degree)
  pools <- split(bg$gene, bg$bin)
  log_msg("Genes per degree bin: ",
          paste(names(lengths(pools)), lengths(pools), sep = "=",
                collapse = ", "))
  L2 <- min(lengths(pools)) >= MIN_BIN
  log_msg("GATE L2 min bin >= ", MIN_BIN, ": ", ifelse(L2, "PASS", "FAIL"))

  hb <- setNames(bin_of(hub$degree), hub$gene)
  hb[is.na(hb)] <- max(bg$bin, na.rm = TRUE)

  # ---- descriptive ------------------------------------------------------
  e_bg <- ecdf(bg$degree)
  pos <- long %>% inner_join(hub, by = "gene") %>%
    mutate(pct = 100 * e_bg(degree))

  by_dir <- pos %>% group_by(direction) %>%
    summarise(n_studies = n_distinct(study_id), n_genes = n(),
              median_degree = median(degree), median_pct = median(pct),
              pct_above_90 = 100 * mean(pct >= 90), .groups = "drop")
  write_csv(by_dir, P("tables", "lit30_final_by_direction.csv"))
  log_msg("=========== DEGREE PERCENTILE BY DIRECTION ===========")
  print(as.data.frame(by_dir %>% mutate(across(where(is.numeric),
                                               ~ round(.x, 1)))))

  up <- pos$pct[pos$direction == "up"]; dn <- pos$pct[pos$direction == "down"]
  if (length(up) > 5 && length(dn) > 5) {
    wt <- wilcox.test(up, dn)
    log_msg("Up vs down, Wilcoxon p = ", signif(wt$p.value, 3),
            "; medians ", round(median(up), 1), " vs ", round(median(dn), 1))
  }
  log_msg("Median percentile, all published hub genes: ",
          round(median(pos$pct), 1))

  by_study <- pos %>% group_by(study_id, cancer, direction) %>%
    summarise(n = n(), median_pct = median(pct), .groups = "drop") %>%
    arrange(desc(median_pct))
  write_csv(by_study, P("tables", "lit30_final_by_study.csv"))

  nonhub <- bg %>% filter(!gene %in% hub$gene)
  roc1 <- pROC::roc(c(rep(1, nrow(hub)), rep(0, nrow(nonhub))),
                    c(hub$degree, nonhub$degree), quiet = TRUE)
  auc <- as.numeric(pROC::auc(roc1)); ci <- as.numeric(pROC::ci.auc(roc1))
  log_msg("AUC of degree separating published hubs from background: ",
          round(auc, 3), " (95% CI ", round(ci[1], 3), " to ",
          round(ci[3], 3), ")")

  # ---- recurrence nulls -------------------------------------------------
  per <- split(long$gene, long$study_id)
  per <- lapply(per, function(g) unique(g[g %in% hub$gene]))
  per <- per[lengths(per) > 0]
  obs_tab <- table(unlist(per))
  obs_recur <- sum(obs_tab >= 2); obs_max <- max(obs_tab)
  log_msg("Observed: ", obs_recur, " genes in 2+ studies; maximum ",
          obs_max, " studies for a single gene")

  draw_matched <- function(gs) map_chr(gs, function(g) {
    b <- as.character(hb[[g]])
    if (is.na(b) || is.null(pools[[b]])) sample(bg$gene, 1)
    else sample(pools[[b]], 1)
  })
  draw_random <- function(gs) sample(bg$gene, length(gs))

  run_null <- function(drawer, label) {
    set.seed(SEED)
    st <- map_dfr(seq_len(N_NULL), function(i) {
      d <- lapply(per, function(gs) unique(drawer(gs)))
      t <- table(unlist(d))
      tibble(recur = sum(t >= 2), mx = if (length(t)) max(t) else 0L)
    })
    tibble(null = label, mean_recur = mean(st$recur), sd = sd(st$recur),
           q05 = quantile(st$recur, 0.05), q95 = quantile(st$recur, 0.95),
           mean_max = mean(st$mx),
           p_more = (sum(st$recur >= obs_recur) + 1) / (N_NULL + 1),
           p_fewer = (sum(st$recur <= obs_recur) + 1) / (N_NULL + 1),
           p_max = (sum(st$mx >= obs_max) + 1) / (N_NULL + 1))
  }

  log_msg("Running ", N_NULL, " replicates per null ...")
  res <- bind_rows(run_null(draw_matched, "degree-matched"),
                   run_null(draw_random,  "random")) %>%
    mutate(observed_recur = obs_recur, observed_max = obs_max)
  write_csv(res, P("tables", "lit30_final_recurrence_null.csv"))
  log_msg("================ RECURRENCE NULLS ================")
  log_msg("  p_more  = P(null >= observed), the one-sided test for excess")
  log_msg("  p_fewer = P(null <= observed), tests whether degree")
  log_msg("            OVER-predicts the observed overlap")
  print(as.data.frame(res %>% mutate(across(where(is.numeric),
                                            ~ round(.x, 3)))))

  set.seed(SEED)
  sub <- unlist(lapply(per, draw_matched))
  obs_med <- median(hub$degree[hub$gene %in% unlist(per)])
  sub_med <- median(bg$degree[match(sub, bg$gene)], na.rm = TRUE)
  err <- abs(obs_med - sub_med) / obs_med
  log_msg("Match quality: observed median degree ", round(obs_med),
          ", substitutes ", round(sub_med), " (error ",
          round(100 * err, 1), "%)")
  L3 <- err <= MAX_MATCH_ERR
  log_msg("GATE L3 match error <= ", 100 * MAX_MATCH_ERR, "%: ",
          ifelse(L3, "PASS", "FAIL"))

  rd <- res %>% filter(null == "random")
  L4 <- rd$mean_recur < obs_recur
  log_msg("GATE L4 random null not saturated: ",
          ifelse(L4, "PASS", "FAIL, pool exhausted"))

  rec <- tibble(gene = names(obs_tab), n_studies = as.integer(obs_tab)) %>%
    left_join(hub, by = "gene") %>%
    mutate(pct = round(100 * e_bg(degree), 1)) %>%
    arrange(desc(n_studies), desc(degree))
  write_csv(rec, P("tables", "lit30_final_recurrent_genes.csv"))
  log_msg("Most recurrent genes:")
  print(as.data.frame(head(rec, 15)))

  gates <- tibble(
    gate = c("L1 resolved background", "L2 min bin size",
             "L3 match quality", "L4 random null usable"),
    value = c(nrow(bg), min(lengths(pools)), round(err, 4),
              round(rd$mean_recur, 1)),
    threshold = c(MIN_BG, MIN_BIN, MAX_MATCH_ERR, obs_recur),
    pass = c(L1, L2, L3, L4))
  write_csv(gates, P("tables", "lit30_final_gates.csv"))
  log_msg("==================== GATES ====================")
  print(as.data.frame(gates))

  dm <- res %>% filter(null == "degree-matched")
  log_msg("---")
  if (!(L1 && L2 && L3)) {
    log_msg("Gates L1 to L3 not all cleared; the recurrence null is not ",
            "reportable. The direction analysis and degree percentiles do ",
            "not depend on it and stand.")
  } else if (dm$p_more > 0.05 && dm$p_fewer > 0.05) {
    log_msg("Observed recurrence is consistent with degree-matched ",
            "substitution (p_more = ", signif(dm$p_more, 3),
            ", p_fewer = ", signif(dm$p_fewer, 3),
            "). Cross-study overlap is accounted for by interactome ",
            "degree.")
  } else if (dm$p_fewer <= 0.05) {
    log_msg("Observed recurrence is BELOW the degree-matched null (",
            obs_recur, " against ", round(dm$mean_recur, 1),
            ", p_fewer = ", signif(dm$p_fewer, 3),
            "). Degree over-predicts the overlap. The overlap is fully ",
            "accounted for by degree, and if anything these studies share ",
            "fewer genes than their connectivity profile alone would ",
            "produce.")
  } else {
    log_msg("Observed recurrence EXCEEDS the degree-matched null (p_more = ",
            signif(dm$p_more, 3), "). Degree accounts for ",
            round(dm$mean_recur, 1), " of ", obs_recur,
            " but not all. State as 'largely explained by degree'.")
  }
  log_msg("Single-gene concentration: observed maximum ", obs_max,
          " studies, null mean ", round(dm$mean_max, 2),
          ", p_max = ", signif(dm$p_max, 3))

  saveRDS(list(long = long, hub = hub, bg = bg, null = res,
               by_direction = by_dir, recurrent = rec, gates = gates,
               auc = auc, auc_ci = ci),
          P("rds", "lit30_final.rds"))
  log_msg("28_literature_scale_analysis complete.")
  invisible(res)
}
