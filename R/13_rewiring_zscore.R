# Centrality scored against a degree-preserving null.
#
#   z = (observed - mean over rewirings) / sd over rewirings
#
# A gene central only because it has many edges sits near zero, since the
# rewired networks reproduce that centrality. This is standard network
# science; what is absent from this literature is its use in the DEG-to-hub
# workflow, where raw centrality remains the norm.
#
# Genes are ranked by mean z across measures rather than by maximum.
# Ranking by maximum cherry-picks each gene's most extreme metric.
#
# A degree floor sensitivity analysis follows, since a metric can move from
# favouring high-degree genes to favouring single-edge genes.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC); library(httr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
slice <- dplyr::slice; arrange <- dplyr::arrange

N_NULL <- 500          # rewirings for the z-score
Z_THRESHOLD <- 2       # nomination cutoff

mcc_score <- function(g) {
  cl <- max_cliques(g, min = 1)
  sc <- setNames(numeric(vcount(g)), V(g)$name)
  for (c in cl) {
    k <- length(c)
    if (k > 1) for (v in names(c)) sc[v] <- sc[v] + factorial(k - 1)
  }
  sc[sc == 0] <- 1
  sc
}

centrality_matrix <- function(g) {
  vs <- V(g)$name
  cbind(
    degree      = degree(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness   = closeness(g, normalized = TRUE),
    eigenvector = eigen_centrality(g)$vector,
    pagerank    = page_rank(g)$vector,
    mcc         = mcc_score(g)[vs])
}

string_edges <- function(genes, score, tag) {
  cd <- file.path(CACHE_DIR, "string_bias")
  dir.create(cd, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cd, paste0(tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- tryCatch(httr::POST("https://string-db.org/api/tsv/network",
        body = list(identifiers = paste(genes, collapse = "%0d"),
                    species = 9606, required_score = score,
                    caller_identity = "glioma_rewire"), encode = "form"),
        error = function(e) NULL)
  e <- NULL
  if (!is.null(r) && httr::status_code(r) == 200) {
    d <- tryCatch(readr::read_tsv(I(httr::content(r, as = "text",
                                                  encoding = "UTF-8")),
                                  show_col_types = FALSE),
                  error = function(e) NULL)
    if (!is.null(d) && nrow(d))
      e <- d %>% select(from = preferredName_A, to = preferredName_B) %>%
        filter(from != to) %>% distinct()
  }
  saveRDS(e, f); e
}

# --- The metric ----------------------------------------------------------
rewiring_zscores <- function(g, n_null = N_NULL) {
  obs <- centrality_matrix(g)
  vs  <- rownames(obs)
  mets <- colnames(obs)

  # Accumulate mean and sd across rewirings without storing every matrix
  sums <- array(0, dim = dim(obs), dimnames = dimnames(obs))
  sqs  <- sums

  log_msg("Running ", n_null, " degree-preserving rewirings ...")
  for (i in seq_len(n_null)) {
    gr <- rewire(g, keeping_degseq(niter = ecount(g) * 10))
    m <- centrality_matrix(gr)
    m <- m[vs, , drop = FALSE]
    sums <- sums + m
    sqs  <- sqs + m^2
    if (i %% 100 == 0) log_msg("  ", i, "/", n_null)
  }

  mu <- sums / n_null
  sdv <- sqrt(pmax(sqs / n_null - mu^2, 0))
  z <- (obs - mu) / ifelse(sdv < 1e-10, NA, sdv)

  out <- as_tibble(z) %>% mutate(gene = vs) %>%
    relocate(gene)
  names(out)[-1] <- paste0("z_", mets)

  # Composite: how many measures exceed the threshold, and the max z
  zm <- as.matrix(out[, -1, drop = FALSE])
  out$n_metrics_z <- rowSums(zm > Z_THRESHOLD, na.rm = TRUE)
  out$max_z <- apply(zm, 1, function(v) suppressWarnings(max(v, na.rm = TRUE)))
  out$max_z[!is.finite(out$max_z)] <- NA_real_
  out$mean_z <- rowMeans(zm, na.rm = TRUE)

  # Raw values for reference
  raw <- as_tibble(obs) %>% mutate(gene = vs) %>% relocate(gene)
  names(raw)[-1] <- paste0("raw_", mets)

  # Ranked by mean z across measures rather than by maximum; the maximum
  # selects each gene's most extreme metric.
  left_join(out, raw, by = "gene") %>% arrange(desc(mean_z))
}

# --- Evaluation ----------------------------------------------------------
predictability_auc <- function(nominated, all_genes, global_degree) {
  d <- tibble(gene = all_genes,
              is_nom = as.integer(all_genes %in% nominated),
              gd = global_degree[match(all_genes, names(global_degree))]) %>%
    filter(!is.na(gd))
  if (sum(d$is_nom) < 3 || sum(d$is_nom) == nrow(d)) return(NA_real_)
  fit <- glm(is_nom ~ log10(gd + 1), data = d, family = binomial())
  as.numeric(pROC::auc(pROC::roc(d$is_nom, predict(fit, type = "response"),
                                 quiet = TRUE)))
}

main_13 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene
  net   <- readRDS(P("rds", "network.rds"))
  score <- as.integer(net$chosen)

  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  sc <- read_csv(P("tables", "specific_connectivity_binomial.csv"),
                 show_col_types = FALSE)
  gd <- setNames(norm$global_degree, norm$gene)

  e <- string_edges(genes, score, "observed")
  if (is.null(e)) stop("STRING network unavailable.")
  g <- graph_from_data_frame(e, directed = FALSE)
  miss <- setdiff(genes, V(g)$name)
  if (length(miss)) g <- add_vertices(g, length(miss), name = miss)
  log_msg("Network: ", vcount(g), " nodes, ", ecount(g), " edges")

  z <- rewiring_zscores(g)
  write_csv(z, P("tables", "rewiring_zscores.csv"))

  log_msg("Top genes by MEAN rewiring z-score:")
  print(as.data.frame(z %>%
    select(gene, mean_z, max_z, n_metrics_z, raw_degree) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>% head(20)))
  log_msg("Where the conventional hubs land on mean z:")
  print(as.data.frame(z %>% filter(gene %in% net$consensus$gene) %>%
    select(gene, mean_z, max_z, raw_degree) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))))

  # --- Nomination sets ---------------------------------------------------
  conv <- net$consensus$gene
  binom_hits <- sc$gene[sc$significant]
  z_hits <- z$gene[!is.na(z$mean_z) & z$mean_z > Z_THRESHOLD]
  z_hits_max <- z$gene[!is.na(z$max_z) & z$max_z > Z_THRESHOLD]
  combined <- intersect(z_hits, binom_hits)

  log_msg("Nominations:")
  log_msg("  conventional centrality : ", length(conv), " -> ",
          paste(conv, collapse = ", "))
  log_msg("  global degree correction: ", length(binom_hits), " -> ",
          paste(binom_hits, collapse = ", "))
  log_msg("  rewiring z-score        : ", length(z_hits), " -> ",
          paste(head(z_hits, 15), collapse = ", "))
  log_msg("  both corrections        : ", length(combined), " -> ",
          paste(combined, collapse = ", "))

  # --- The evaluation ----------------------------------------------------
  ev <- tibble(
    method = c("raw centrality (conventional)",
               "global degree correction (binomial)",
               "rewiring z-score (max, cherry-picked)",
               "rewiring z-score (mean)",
               "both corrections"),
    n_nominated = c(length(conv), length(binom_hits), length(z_hits_max),
                    length(z_hits), length(combined)),
    auc_predictable_from_global_degree = c(
      predictability_auc(conv, genes, gd),
      predictability_auc(binom_hits, genes, gd),
      predictability_auc(z_hits_max, genes, gd),
      predictability_auc(z_hits, genes, gd),
      predictability_auc(combined, genes, gd)))

  write_csv(ev, P("tables", "hub_metric_comparison.csv"))

  log_msg("=================================================")
  log_msg("METRIC COMPARISON")
  log_msg("  AUC = how well global interactome degree predicts which genes")
  log_msg("  the method nominates. LOWER IS BETTER: it means the method is")
  log_msg("  reading the network rather than the annotation record.")
  log_msg("  0.5 would mean no relationship to how well studied a gene is.")
  log_msg("=================================================")
  print(as.data.frame(ev %>%
    mutate(auc_predictable_from_global_degree =
             round(auc_predictable_from_global_degree, 3))))

  base <- ev$auc_predictable_from_global_degree[1]
  ev_ok <- ev[!is.na(ev$auc_predictable_from_global_degree), ]
  best <- ev_ok[which.min(ev_ok$auc_predictable_from_global_degree), ]

  log_msg("Best method: ", best$method, " at AUC ",
          round(best$auc_predictable_from_global_degree, 3),
          " versus ", round(base, 3), " for raw centrality")

  if (best$auc_predictable_from_global_degree < 0.7) {
    log_msg("VERDICT: the correction substantially removes the dependence ",
            "on annotation density. Recommend it as a replacement metric.")
  } else if (best$auc_predictable_from_global_degree < base - 0.1) {
    log_msg("VERDICT: the correction reduces but does not remove the ",
            "dependence. Report it as a partial improvement and recommend ",
            "the diagnostic alongside it.")
  } else {
    log_msg("VERDICT: no correction meaningfully improves on raw centrality. ",
            "The recommendation becomes: run the diagnostic and treat hub ",
            "lists from this workflow with caution. Report this plainly.")
  }

  # --- Do the corrections agree with each other? -------------------------
  if (length(z_hits) && length(binom_hits)) {
    j <- length(intersect(z_hits, binom_hits)) /
         length(union(z_hits, binom_hits))
    log_msg("Jaccard between the two corrections: ", round(j, 3))
    log_msg("  They use different information, the global annotation record ",
            "versus the observed network topology, so agreement is a ",
            "meaningful check rather than a foregone conclusion.")
  }

  # --- Degree floor sensitivity -----------------------------------------
  # A minimum degree is imposed and the evaluation repeated, since a single
  # edge in a sparse network is a weak observation.
  log_msg("=== DEGREE FLOOR SENSITIVITY ===")
  deg <- setNames(degree(g), V(g)$name)
  z$degree_obs <- unname(deg[match(z$gene, names(deg))])

  floors <- c(1, 2, 3, 5)
  fl <- map_dfr(floors, function(mind) {
    elig <- z$gene[!is.na(z$degree_obs) & z$degree_obs >= mind]
    if (length(elig) < 20) return(tibble())
    zz <- z %>% filter(gene %in% elig)
    hits <- zz$gene[!is.na(zz$mean_z) & zz$mean_z > Z_THRESHOLD]
    # relax the threshold if too few clear it at this floor
    if (length(hits) < 3) {
      hits <- zz %>% arrange(desc(mean_z)) %>% head(7) %>% pull(gene)
      note <- "top 7 by mean z"
    } else note <- paste0("mean z > ", Z_THRESHOLD)
    tibble(min_degree = mind,
           n_eligible = length(elig),
           n_nominated = length(hits),
           rule = note,
           median_degree_nominated = median(zz$degree_obs[zz$gene %in% hits],
                                            na.rm = TRUE),
           auc = predictability_auc(hits, elig, gd),
           genes = paste(head(hits, 10), collapse = "; "))
  })

  write_csv(fl, P("tables", "rewiring_zscore_degree_floor.csv"))
  print(as.data.frame(fl %>% select(-genes) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))))
  for (i in seq_len(nrow(fl)))
    log_msg("  min degree ", fl$min_degree[i], ": ", fl$genes[i])

  if (nrow(fl) >= 2) {
    hi <- fl %>% filter(min_degree >= 3)
    if (nrow(hi) && all(hi$auc < 0.7, na.rm = TRUE)) {
      log_msg("The correction holds among genes with at least 3 edges. ",
              "It is not simply favouring sparsely connected genes.")
    } else if (nrow(hi)) {
      log_msg("AUC rises once low-degree genes are excluded (",
              paste(round(hi$auc, 3), collapse = ", "),
              "). The metric is partly exploiting sparsity, and that ",
              "limitation must be reported.")
    }
  }

  # Degree distribution of each method's nominations, for the same reason
  dd <- tibble(
    method = c("conventional", "binomial", "z-score (mean)"),
    median_degree = c(
      median(deg[names(deg) %in% conv], na.rm = TRUE),
      median(deg[names(deg) %in% binom_hits], na.rm = TRUE),
      median(deg[names(deg) %in% z_hits], na.rm = TRUE)),
    range_degree = c(
      paste(range(deg[names(deg) %in% conv]), collapse = "-"),
      paste(range(deg[names(deg) %in% binom_hits]), collapse = "-"),
      paste(range(deg[names(deg) %in% z_hits]), collapse = "-")))
  write_csv(dd, P("tables", "nomination_degree_profile.csv"))
  log_msg("Degree profile of each method's nominations:")
  print(as.data.frame(dd))

  saveRDS(list(z = z, evaluation = ev, degree_floor = fl),
          P("rds", "rewiring_zscore.rds"))
  log_msg("13_rewiring_zscore complete.")
  invisible(ev)
}
