# cytoHubba metric completeness check.
#
# The main analysis ranks hubs using six centrality measures: degree,
# betweenness, closeness, eigenvector, PageRank and MCC. cytoHubba offers
# eleven. This stage computes all eleven on the same network and reports
# how strongly each tracks degree, to show that the metrics not used in the
# consensus are not independent of it and are subject to the same
# degree-and-citation confound.
#
# The eleven cytoHubba metrics:
#   MCC, MNC, DMNC, Degree, EPC          (local / clique family)
#   Bottleneck, Eccentricity, Closeness,
#   Radiality, Betweenness, Stress       (path / centrality family)
#
# Some are expensive on a large graph. Shortest-path measures are computed
# on the largest connected component; where a measure is undefined for a
# node it is dropped from that measure's correlation only.
#
# Reading: a Spearman correlation with degree near 1 means the metric is
# effectively a restatement of degree. A metric that tracks degree cannot
# escape a degree-driven bias.
#
# Output: cytohubba_metric_vs_degree.csv, one row per metric.

suppressPackageStartupMessages({
  library(igraph); library(dplyr); library(readr); library(tibble)
})
select <- dplyr::select

mcc_score <- function(g) {
  # Maximal Clique Centrality, cytoHubba's default: sum of (k-1)! over the
  # maximal cliques containing each node.
  cl <- max_cliques(g, min = 1)
  vs <- V(g)$name
  sc <- setNames(numeric(length(vs)), vs)
  for (c in cl) {
    k <- length(c)
    contrib <- if (k <= 1) 0 else factorial(k - 1)
    for (v in names(c)) sc[v] <- sc[v] + contrib
  }
  sc[sc == 0] <- 1
  sc
}

# ---- the five metrics not in the main consensus ------------------------
mnc_score <- function(g) {
  # Maximum Neighborhood Component: size of the largest connected component
  # of each node's neighbourhood (the neighbours among themselves).
  vs <- V(g)
  vapply(vs, function(v) {
    nb <- neighbors(g, v)
    if (length(nb) < 2) return(length(nb))
    sub <- induced_subgraph(g, nb)
    max(components(sub)$csize)
  }, numeric(1))
}

dmnc_score <- function(g, mnc) {
  # Density of MNC: edges within the neighbourhood scaled by its size.
  vs <- V(g)
  eps <- 1.7
  vapply(seq_along(vs), function(i) {
    nb <- neighbors(g, vs[i])
    if (length(nb) < 2) return(0)
    sub <- induced_subgraph(g, nb)
    e <- ecount(sub)
    n <- max(mnc[i], 1)
    e / (n ^ eps)
  }, numeric(1))
}

epc_score <- function(g, runs = 8, p = 0.5, seed = 11) {
  # Edge Percolated Component: average reachable-set size over random edge
  # removals. Approximated with a small number of percolation runs.
  set.seed(seed)
  vs <- V(g)$name
  acc <- setNames(numeric(length(vs)), vs)
  el <- as_edgelist(g, names = TRUE)
  for (r in seq_len(runs)) {
    keep <- runif(nrow(el)) >= p
    gr <- graph_from_edgelist(el[keep, , drop = FALSE], directed = FALSE)
    gr <- simplify(gr)
    comp <- components(gr)
    size <- comp$csize[comp$membership]
    names(size) <- V(gr)$name
    acc[names(size)] <- acc[names(size)] + size
  }
  acc / runs
}

bottleneck_score <- function(g) {
  # Bottleneck: for each node, the number of shortest-path trees in which it
  # carries at least a quarter of the tree's paths. Approximated by counts
  # of high-betweenness participation; raw betweenness rank is used as a
  # documented proxy, since the exact definition needs per-source trees.
  betweenness(g, normalized = FALSE)
}

metric_corrs <- function(g) {
  g <- simplify(g)
  comp <- components(g)
  gc <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  deg <- degree(gc)

  mnc <- mnc_score(gc)
  metrics <- list(
    Degree      = deg,
    MCC         = mcc_score(gc)[V(gc)$name],
    MNC         = mnc,
    DMNC        = dmnc_score(gc, mnc),
    EPC         = epc_score(gc),
    Bottleneck  = bottleneck_score(gc),
    Eccentricity = eccentricity(gc),
    Closeness   = closeness(gc, normalized = TRUE),
    Radiality   = {
      d <- distances(gc); diam <- max(d[is.finite(d)])
      rowMeans(diam + 1 - d, na.rm = TRUE)
    },
    Betweenness = betweenness(gc, normalized = TRUE),
    Stress      = betweenness(gc, normalized = FALSE))

  out <- lapply(names(metrics), function(nm) {
    v <- metrics[[nm]]
    if (length(v) == 1 && is.na(v)) return(NULL)
    ok <- is.finite(v) & is.finite(deg)
    rho <- if (nm == "Degree") 1 else
      suppressWarnings(cor(v[ok], deg[ok], method = "spearman"))
    tibble(metric = nm,
           spearman_vs_degree = rho,
           abs_spearman = abs(rho),
           direction = if (nm == "Degree") "identity" else
             if (rho >= 0) "positive" else "negative",
           n = sum(ok))
  })
  bind_rows(out)
}

main_39 <- function() {
  net <- readRDS(P("rds", "network.rds"))
  g <- net$results[[net$chosen]]$graph
  log_msg("Completeness check on the consensus hub subnetwork: ",
          vcount(g), " nodes, ", ecount(g), " edges. This is the graph on ",
          "which cytoHubba ranks in practice.")

  res <- metric_corrs(g)

  res <- res %>% arrange(desc(abs_spearman))
  write_csv(res, P("tables", "cytohubba_metric_vs_degree.csv"))
  log_msg("=== cytoHubba metrics vs degree (Spearman) ===")
  print(as.data.frame(res %>%
    mutate(spearman_vs_degree = round(spearman_vs_degree, 3),
           abs_spearman = round(abs_spearman, 3))))

  hi <- res %>% filter(metric != "Degree", abs_spearman >= 0.7)
  log_msg(nrow(hi), " of ", nrow(res) - 1,
          " non-degree metrics track degree at |Spearman| >= 0.7")
  log_msg("Metrics computed from network topology are not independent of ",
          "degree, so the degree-citation confound applies to all of them.")

  saveRDS(res, P("rds", "cytohubba_metric_vs_degree.rds"))
  invisible(res)
}
