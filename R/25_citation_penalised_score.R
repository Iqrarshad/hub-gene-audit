# Hub scoring that penalises study intensity rather than degree.
#
# Degree is part study intensity and part interaction propensity;
# penalising degree removes both. Publication count is measured
# independently of STRING, so the attention component can be removed alone.
#
# Three variants are tested: node-level residual adjustment, edge-level
# down-weighting of edges between heavily studied genes, and pure residual
# ranking. Each is evaluated on both bias measures and on DepMap
# selectivity, so a score cannot pass by reducing bias while finding
# nothing.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; distinct <- dplyr::distinct
arrange <- dplyr::arrange; summarise <- dplyr::summarise

N_NOMINATE <- 8

mcc_score <- function(g) {
  cl <- max_cliques(g, min = 1)
  s <- setNames(numeric(vcount(g)), V(g)$name)
  for (c in cl) { k <- length(c)
    if (k > 1) for (v in names(c)) s[v] <- s[v] + factorial(k - 1) }
  s[s == 0] <- 1; s
}

centrality_consensus <- function(g, weights = NULL, top_k = 10,
                                 min_methods = 4) {
  if (ecount(g) == 0) return(list(tab = NULL, hubs = character(0)))
  vs <- V(g)$name
  w <- if (is.null(weights)) NULL else weights
  tab <- tibble(
    gene = vs,
    degree = if (is.null(w)) degree(g) else strength(g, weights = w),
    betweenness = betweenness(g, weights = if (is.null(w)) NULL else 1 / w,
                              normalized = TRUE),
    closeness = closeness(g, weights = if (is.null(w)) NULL else 1 / w,
                          normalized = TRUE),
    eigenvector = eigen_centrality(g, weights = w)$vector,
    pagerank = page_rank(g, weights = w)$vector,
    mcc = mcc_score(g)[vs])
  mets <- setdiff(names(tab), "gene")
  top <- vapply(mets, function(m)
    as.integer(tab$gene %in% tab$gene[head(order(tab[[m]],
                                                 decreasing = TRUE), top_k)]),
    integer(nrow(tab)))
  tab$n_methods <- rowSums(top)
  list(tab = tab, hubs = tab$gene[tab$n_methods >= min_methods])
}

auc_of <- function(nominated, all_genes, predictor) {
  d <- tibble(gene = all_genes,
              y = as.integer(all_genes %in% nominated),
              x = predictor[match(all_genes, names(predictor))]) %>%
    filter(!is.na(x))
  if (sum(d$y) < 3 || sum(d$y) == nrow(d)) return(NA_real_)
  fit <- glm(y ~ x, data = d, family = binomial())
  as.numeric(pROC::auc(pROC::roc(d$y, predict(fit, type = "response"),
                                 quiet = TRUE)))
}

main_25 <- function() {
  lb <- readRDS(P("rds", "literature_bias.rds"))
  pap_all <- lb$data %>% select(gene, degree, n_papers, degree_resid)

  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared; if (is.data.frame(genes)) genes <- genes$gene
  net <- readRDS(P("rds", "network.rds"))
  conv <- net$consensus$gene
  score <- as.integer(net$chosen)

  ef <- file.path(CACHE_DIR, "string_density", paste0("t", score, ".rds"))
  edges <- if (file.exists(ef)) readRDS(ef) else NULL
  if (is.null(edges) || !nrow(edges))
    edges <- readRDS(P("rds", paste0("string_", score, ".rds")))
  edges <- edges %>% filter(from %in% genes, to %in% genes)
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = genes)
  log_msg("Network: ", vcount(g), " nodes, ", ecount(g), " edges")

  # Papers per gene, with a floor so unstudied genes are not infinite
  papers <- setNames(pap_all$n_papers, pap_all$gene)
  deg_glob <- setNames(pap_all$degree, pap_all$gene)
  resid <- setNames(pap_all$degree_resid, pap_all$gene)

  miss <- setdiff(genes, names(papers))
  if (length(miss)) {
    log_msg(length(miss), " genes lack publication counts; assigned the ",
            "10th percentile so they are neither rewarded nor excluded")
    fill <- quantile(pap_all$n_papers, 0.10, na.rm = TRUE)
    papers <- c(papers, setNames(rep(fill, length(miss)), miss))
  }

  # --- Variant A: node-level residual adjustment -------------------------
  cc <- centrality_consensus(g)
  tabA <- cc$tab %>%
    mutate(resid = resid[match(gene, names(resid))],
           resid = ifelse(is.na(resid), 0, resid),
           # rank-combine consensus count with the residual
           score_A = rank(n_methods) + rank(resid))
  hubsA <- tabA %>% arrange(desc(score_A)) %>% head(N_NOMINATE) %>% pull(gene)

  # --- Variant B: edge-level down-weighting ------------------------------
  # An edge between two heavily studied genes carries less evidential weight
  # than an edge between two rarely studied ones, because the former is more
  # likely to have been reported through attention alone.
  el <- as_data_frame(g, what = "edges")
  if (nrow(el)) {
    lp <- log10(papers + 1)
    w <- 1 / (1 + lp[el$from] * lp[el$to])
    w[!is.finite(w)] <- min(w[is.finite(w)], na.rm = TRUE)
    ccB <- centrality_consensus(g, weights = as.numeric(w))
    hubsB <- ccB$tab %>% arrange(desc(n_methods),
                                 desc(degree)) %>%
      head(N_NOMINATE) %>% pull(gene)
  } else hubsB <- character(0)

  # --- Variant C: pure residual ranking ----------------------------------
  hubsC <- tibble(gene = genes,
                  r = resid[match(genes, names(resid))]) %>%
    filter(!is.na(r)) %>% arrange(desc(r)) %>% head(N_NOMINATE) %>%
    pull(gene)

  sets <- list(`cytoHubba (conventional)` = conv,
               `A: consensus + residual`  = hubsA,
               `B: citation-weighted edges` = hubsB,
               `C: degree residual only`  = hubsC)
  for (nm in names(sets))
    log_msg(nm, ": ", paste(sets[[nm]], collapse = ", "))

  # --- Evaluation --------------------------------------------------------
  pool <- intersect(genes, names(deg_glob))
  ev <- imap_dfr(sets, function(h, nm) {
    tibble(method = nm, n = length(h),
           auc_degree_bias = auc_of(h, pool, log10(deg_glob + 1)),
           auc_citation_bias = auc_of(h, pool, log10(papers[pool] + 1)),
           median_papers = median(papers[h], na.rm = TRUE),
           median_degree = median(deg_glob[h], na.rm = TRUE))
  })

  # --- DepMap ------------------------------------------------------------
  self <- P("tables", "depmap_glioma_selectivity.csv")
  if (file.exists(self)) {
    sel <- read_csv(self, show_col_types = FALSE)
    dep <- imap_dfr(sets, function(h, nm) {
      s <- sel %>% mutate(in_set = gene %in% h)
      if (sum(s$in_set) < 3) return(tibble(method = nm,
                                           auc_selectivity = NA_real_,
                                           pct_pan_essential = NA_real_))
      tibble(method = nm,
             auc_selectivity = as.numeric(pROC::auc(
               pROC::roc(s$in_set, s$selectivity, quiet = TRUE))),
             pct_pan_essential = 100 * mean(s$pan_essential[s$in_set]))
    })
    ev <- ev %>% left_join(dep, by = "method")
  }

  write_csv(ev, P("tables", "citation_penalised_comparison.csv"))

  log_msg("=================================================")
  log_msg("CITATION-PENALISED SCORING")
  log_msg("  auc_degree_bias / auc_citation_bias: LOWER better, 0.5 ideal")
  log_msg("  auc_selectivity: HIGHER better, guards against finding nothing")
  log_msg("  benchmarks: cytoHubba 0.994 bias, binomial 0.825,")
  log_msg("              co-expression 0.576")
  log_msg("=================================================")
  print(as.data.frame(ev %>% mutate(across(where(is.numeric),
                                           ~ round(.x, 3)))))

  # --- Verdict -----------------------------------------------------------
  base <- ev %>% filter(grepl("cytoHubba", method))
  alts <- ev %>% filter(!grepl("cytoHubba", method))
  best <- alts[which.min(alts$auc_citation_bias), ]

  log_msg("---")
  log_msg("Best alternative: ", best$method,
          " (citation bias ", round(best$auc_citation_bias, 3),
          " against ", round(base$auc_citation_bias, 3), " for cytoHubba)")

  if (!is.na(best$auc_citation_bias) && best$auc_citation_bias < 0.65 &&
      !is.na(best$auc_selectivity) && best$auc_selectivity > 0.6) {
    log_msg("VERDICT: the citation-penalised score removes most of the ",
            "attention bias AND retains glioma-selective dependency ",
            "signal. This is the first correction to satisfy both. Worth ",
            "developing into a method.")
  } else if (!is.na(best$auc_citation_bias) &&
             best$auc_citation_bias < 0.65) {
    log_msg("VERDICT: the score removes attention bias but shows no ",
            "glioma-selective dependency (AUC ",
            round(best$auc_selectivity, 3),
            "). This repeats the rewiring z-score outcome: a less biased ",
            "score that finds nothing. Not a method.")
  } else {
    log_msg("VERDICT: the citation-penalised score does not meaningfully ",
            "reduce bias. Seven corrections have now failed. The paper ",
            "reports the problem, its mechanism (64 percent of degree is ",
            "publication count), and that no simple reweighting resolves it.")
  }

  saveRDS(ev, P("rds", "citation_penalised.rds"))
  log_msg("25_citation_penalised_score complete.")
  invisible(ev)
}
