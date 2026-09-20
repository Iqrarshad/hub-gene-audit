# STRING network, topology diagnostics and hub selection.
#
# Connectivity diagnostics are computed before ranking. If the largest
# connected component holds under 30% of nodes, centrality measures are not
# interpretable; the threshold is then relaxed and both results reported.
#
# Six distinct centrality measures are used. Earlier versions also carried
# stress and radiality, which as implemented were rank-duplicates of
# betweenness and closeness and inflated the consensus count.

suppressPackageStartupMessages({
  library(igraph); library(dplyr); library(readr); library(tibble)
  library(httr); library(purrr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

STRING_API <- "https://string-db.org/api"

fetch_string_edges <- function(genes, required_score = 700, species = 9606) {
  cache <- P("rds", paste0("string_", required_score, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))

  log_msg("Querying STRING for ", length(genes), " genes at score >= ",
          required_score)

  # STRING caps identifiers per request; chunk conservatively
  chunks <- split(genes, ceiling(seq_along(genes) / 400))
  edges <- map_dfr(seq_along(chunks), function(i) {
    log_msg("  chunk ", i, " of ", length(chunks))
    r <- httr::POST(
      paste0(STRING_API, "/tsv/network"),
      body = list(identifiers    = paste(chunks[[i]], collapse = "%0d"),
                  species        = species,
                  required_score = required_score,
                  caller_identity = "hub_gene_audit"),
      encode = "form")
    if (httr::status_code(r) != 200) {
      warning("STRING chunk ", i, " returned ", httr::status_code(r))
      return(tibble())
    }
    txt <- httr::content(r, as = "text", encoding = "UTF-8")
    read_tsv(txt, show_col_types = FALSE)
  })

  if (nrow(edges) == 0) stop("STRING returned no interactions.")
  edges <- edges %>%
    select(from = preferredName_A, to = preferredName_B, score) %>%
    filter(from != to) %>% distinct()

  saveRDS(edges, cache)
  edges
}

network_diagnostics <- function(g, label) {
  comp <- components(g)
  lcc  <- max(comp$csize)
  d <- tibble(
    label            = label,
    n_nodes          = vcount(g),
    n_edges          = ecount(g),
    density          = edge_density(g),
    n_components     = comp$no,
    largest_component = lcc,
    lcc_fraction     = lcc / vcount(g),
    n_isolated       = sum(degree(g) == 0),
    mean_degree      = mean(degree(g)),
    transitivity     = transitivity(g, type = "global")
  )
  log_msg(label, ": ", d$n_nodes, " nodes, ", d$n_edges, " edges, density ",
          signif(d$density, 3), ", LCC ", lcc, " (",
          round(100 * d$lcc_fraction, 1), "%), isolated ", d$n_isolated)
  d
}

# --- cytoHubba-equivalent rankings --------------------------------------
# MCC (Maximal Clique Centrality) is cytoHubba's default and the one that
# most influenced the original hub list, so it is implemented explicitly.
mcc_score <- function(g) {
  cl <- max_cliques(g, min = 1)
  vs <- V(g)$name
  sc <- setNames(numeric(length(vs)), vs)
  for (c in cl) {
    k <- length(c)
    contrib <- if (k <= 1) 0 else factorial(k - 1)
    for (v in names(c)) sc[v] <- sc[v] + contrib
  }
  sc[sc == 0] <- 1   # cytoHubba convention for nodes in no clique
  sc
}

rank_hubs <- function(g) {
  vs <- V(g)$name
  m <- tibble(
    gene        = vs,
    degree      = degree(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness   = closeness(g, normalized = TRUE),
    eigenvector = eigen_centrality(g)$vector,
    pagerank    = page_rank(g)$vector,
    mcc         = mcc_score(g)[vs],
    clustering  = transitivity(g, type = "local", isolates = "zero")
  )

  # Six genuinely distinct measures. Earlier versions of this script also
  # Six distinct measures.
  methods <- c("degree", "betweenness", "closeness", "eigenvector",
               "pagerank", "mcc")

  # Top-10 membership per method, then consensus count
  top10 <- vapply(methods, function(mth) {
    ord <- order(m[[mth]], decreasing = TRUE)
    as.integer(m$gene %in% m$gene[head(ord, 10)])
  }, integer(nrow(m)))

  m$n_methods_top10 <- rowSums(top10)
  m %>% arrange(desc(n_methods_top10), desc(degree))
}

main_06 <- function() {
  inter  <- readRDS(P("rds", "intersect.rds"))
  genes  <- inter$shared$gene
  log_msg("Building PPI network from ", length(genes), " shared DEGs")

  results <- list()
  chosen  <- NULL

  for (sc in c(THRESH$string_score, THRESH$string_relax)) {
    edges <- tryCatch(fetch_string_edges(genes, sc),
                      error = function(e) { warning(conditionMessage(e)); NULL })
    if (is.null(edges)) next

    g <- graph_from_data_frame(edges, directed = FALSE)
    # Include DEGs with no interactions so the diagnostics are honest
    missing <- setdiff(genes, V(g)$name)
    if (length(missing) > 0) g <- add_vertices(g, length(missing), name = missing)

    diag <- network_diagnostics(g, paste0("STRING_", sc))
    hubs <- rank_hubs(g)

    write_csv(diag,  P("tables", paste0("network_diagnostics_", sc, ".csv")))
    write_csv(hubs,  P("tables", paste0("hub_ranking_", sc, ".csv")))
    write_csv(edges, P("tables", paste0("network_edges_", sc, ".csv")))

    results[[as.character(sc)]] <- list(graph = g, diag = diag, hubs = hubs)

    if (is.null(chosen) && diag$lcc_fraction >= THRESH$network_lcc_min) {
      chosen <- as.character(sc)
      log_msg("Score ", sc, " passes the LCC criterion (",
              round(100 * diag$lcc_fraction, 1), "% >= ",
              100 * THRESH$network_lcc_min, "%). Using as primary.")
    }
  }

  if (is.null(chosen)) {
    log_msg("WARNING: no threshold produces a connected component holding ",
            100 * THRESH$network_lcc_min, "% of nodes. Centrality-based hub ",
            "selection is not defensible on this network. Report the ",
            "diagnostics and reframe hub selection.")
    chosen <- names(results)[1]
  }

  hubs <- results[[chosen]]$hubs
  # Top-10 by at least four of six measures, a two-thirds majority.
  consensus <- hubs %>% filter(n_methods_top10 >= 4)
  write_csv(consensus, P("tables", "hub_genes_consensus.csv"))

  # Agreement with the prior hub list
  agree <- tibble(
    gene       = HUB_GENES_PRIOR,
    recovered  = HUB_GENES_PRIOR %in% consensus$gene,
    n_methods  = hubs$n_methods_top10[match(HUB_GENES_PRIOR, hubs$gene)],
    degree     = hubs$degree[match(HUB_GENES_PRIOR, hubs$gene)]
  )
  write_csv(agree, P("tables", "hub_agreement_with_prior.csv"))
  log_msg("Prior hub genes recovered: ", sum(agree$recovered, na.rm = TRUE),
          " of ", length(HUB_GENES_PRIOR))

  saveRDS(list(results = results, chosen = chosen, consensus = consensus,
               agreement = agree), P("rds", "network.rds"))
  log_msg("06_network complete.")
  invisible(consensus)
}
