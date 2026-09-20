# Three diagnostics for centrality-based hub selection.
#
# 1. Degree normalisation: local degree against the gene's total
#    interactome degree, which is the direct correction for study bias.
# 2. Physical-evidence-only network, since the text-mining channel scores
#    co-mention in abstracts.
# 3. Correlation among the centrality measures, so the consensus rule can
#    be judged rather than assumed.
#
# Adding more centrality measures does not help: all are functions of the
# same edge set and inherit any bias in it.

suppressPackageStartupMessages({
  library(igraph); library(dplyr); library(readr); library(tibble)
  library(httr); library(purrr); library(tidyr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

API_DELAY <- 0.5

# --- Global interactome degree ------------------------------------------
# STRING's interaction_partners endpoint returns a gene's partners across
# the whole interactome, which is the denominator we need.
global_degree <- function(gene, score) {
  cache <- file.path(CACHE_DIR, "string_global")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache, paste0(gene, "_", score, ".rds"))
  if (file.exists(f)) return(readRDS(f))

  r <- tryCatch(httr::GET("https://string-db.org/api/tsv/interaction_partners",
                          query = list(identifiers = gene, species = 9606,
                                       required_score = score, limit = 5000,
                                       caller_identity = "glioma_centrality")),
                error = function(e) NULL)
  n <- NA_integer_
  if (!is.null(r) && httr::status_code(r) == 200) {
    txt <- httr::content(r, as = "text", encoding = "UTF-8")
    d <- tryCatch(readr::read_tsv(I(txt), show_col_types = FALSE),
                  error = function(e) NULL)
    if (!is.null(d)) n <- nrow(d)
  }
  saveRDS(n, f)
  n
}

# --- Network builder with channel control -------------------------------
build_network <- function(genes, score, physical = FALSE) {
  tag <- paste0(if (physical) "phys" else "full", "_", score)
  f <- P("rds", paste0("string_", tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))

  q <- list(identifiers = paste(genes, collapse = "%0d"),
            species = 9606, required_score = score,
            caller_identity = "glioma_centrality")
  if (physical) q$network_type <- "physical"

  r <- tryCatch(httr::POST("https://string-db.org/api/tsv/network",
                           body = q, encode = "form"),
                error = function(e) NULL)
  if (is.null(r) || httr::status_code(r) != 200) return(NULL)
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  e <- tryCatch(readr::read_tsv(I(txt), show_col_types = FALSE),
                error = function(e) NULL)
  if (is.null(e) || nrow(e) == 0) return(NULL)
  e <- e %>% select(from = preferredName_A, to = preferredName_B) %>%
    filter(from != to) %>% distinct()
  saveRDS(e, f)
  e
}

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

centrality_table <- function(edges, all_genes) {
  g <- graph_from_data_frame(edges, directed = FALSE)
  miss <- setdiff(all_genes, V(g)$name)
  if (length(miss)) g <- add_vertices(g, length(miss), name = miss)
  vs <- V(g)$name
  tibble(gene = vs,
         degree = degree(g),
         betweenness = betweenness(g, normalized = TRUE),
         closeness = closeness(g, normalized = TRUE),
         eigenvector = eigen_centrality(g)$vector,
         pagerank = page_rank(g)$vector,
         mcc = mcc_score(g)[vs])
}

consensus_from <- function(tab) {
  methods <- c("degree", "betweenness", "closeness", "eigenvector",
               "pagerank", "mcc")
  top10 <- vapply(methods, function(m)
    as.integer(tab$gene %in% tab$gene[head(order(tab[[m]],
                                                 decreasing = TRUE), 10)]),
    integer(nrow(tab)))
  tab$n_methods_top10 <- rowSums(top10)
  tab
}

main_08 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene
  net   <- readRDS(P("rds", "network.rds"))
  score <- as.integer(net$chosen)
  log_msg("Genes: ", length(genes), "; STRING score ", score)

  # --- 1. Full network, with metric redundancy ---------------------------
  ed_full <- build_network(genes, score, physical = FALSE)
  if (is.null(ed_full)) stop("STRING full network query failed.")
  tab_full <- consensus_from(centrality_table(ed_full, genes))

  # Six distinct measures. "stress" and "radiality" were removed: as
  # implemented they were rank-identical to betweenness and closeness, and
  # their presence inflated apparent redundancy.
  mets <- c("degree", "betweenness", "closeness", "eigenvector",
            "pagerank", "mcc")
  cm <- cor(tab_full[, mets], method = "spearman", use = "complete.obs")
  write_csv(as.data.frame(cm) %>% rownames_to_column("metric"),
            P("tables", "centrality_metric_correlations.csv"))
  log_msg("Spearman correlations among the eight centrality metrics:")
  print(round(cm, 2))
  offdiag <- cm[upper.tri(cm)]
  log_msg("Median pairwise correlation: ", round(median(offdiag), 3),
          "; proportion above 0.85: ",
          round(mean(offdiag > 0.85), 2))
  if (median(offdiag) > 0.7) {
    log_msg("The metrics are largely redundant. A consensus rule across ",
            "them does not provide independent support.")
  }

  # --- 2. Physical-only network ------------------------------------------
  ed_phys <- build_network(genes, score, physical = TRUE)
  if (!is.null(ed_phys)) {
    tab_phys <- consensus_from(centrality_table(ed_phys, genes))
    hubs_phys <- tab_phys$gene[tab_phys$n_methods_top10 >= 6]
    log_msg("Physical-only network: ", nrow(ed_phys), " edges (full: ",
            nrow(ed_full), ")")
    log_msg("Physical-only hubs: ",
            if (length(hubs_phys)) paste(hubs_phys, collapse = ", ") else "none")
    write_csv(tab_phys, P("tables", "centrality_physical_only.csv"))

    hubs_full <- tab_full$gene[tab_full$n_methods_top10 >= 6]
    log_msg("Overlap with full-network hubs: ",
            length(intersect(hubs_full, hubs_phys)), " of ", length(hubs_full))
  } else {
    log_msg("Physical-only query failed; skipping.")
    tab_phys <- NULL
  }

  # --- 3. Degree normalisation -------------------------------------------
  # Query every gene in the set, not a high-degree subset. Restricting to
  # the top 30 by raw degree biased the comparison field toward exactly the
  # genes the conventional method favours, which understated the effect.
  cand <- union(tab_full$gene, net$consensus$gene)
  log_msg("Querying global interactome degree for all ", length(cand),
          " genes. Cached, so reruns are instant.")

  gd <- map_dfr(seq_along(cand), function(i) {
    g <- cand[i]
    n <- global_degree(g, score)
    if (i %% 25 == 0) log_msg("  ", i, "/", length(cand))
    Sys.sleep(API_DELAY)
    tibble(gene = g, global_degree = n)
  })
  log_msg("Global degree retrieved for ",
          sum(!is.na(gd$global_degree)), " of ", nrow(gd), " genes")

  norm <- tab_full %>% inner_join(gd, by = "gene") %>%
    filter(!is.na(global_degree), global_degree > 0) %>%
    mutate(
      local_degree = degree,
      # Proportion of a gene's whole-interactome partners that fall inside
      # this gene set. High values mean the gene is specifically connected
      # to these genes rather than to everything.
      capture_ratio = local_degree / pmax(global_degree, 1),
      # Expected local degree if this gene's partners were spread uniformly
      expected_degree = global_degree * (length(genes) / 19000),
      degree_enrichment = local_degree / pmax(expected_degree, 0.01)) %>%
    select(gene, local_degree, global_degree, expected_degree,
           degree_enrichment, capture_ratio, n_methods_top10) %>%
    arrange(desc(degree_enrichment))

  write_csv(norm, P("tables", "centrality_degree_normalised.csv"))

  log_msg("=================================================")
  log_msg("DEGREE-NORMALISED CENTRALITY")
  log_msg("  degree_enrichment = observed local degree / expected given")
  log_msg("  the gene's total STRING interactome degree")
  log_msg("=================================================")
  print(as.data.frame(head(norm, 25)))

  orig <- norm %>% filter(gene %in% net$consensus$gene) %>%
    arrange(desc(degree_enrichment))
  log_msg("The reported hubs under degree normalisation:")
  print(as.data.frame(orig))

  demoted <- orig$gene[orig$degree_enrichment < 1]
  if (length(demoted)) {
    log_msg("Hubs with fewer local connections than expected from their ",
            "global connectivity: ", paste(demoted, collapse = ", "),
            ". Their apparent centrality reflects interactome position.")
  }

  saveRDS(list(full = tab_full, physical = tab_phys, normalised = norm,
               metric_cor = cm), P("rds", "centrality_diagnostics.rds"))
  log_msg("08_centrality_diagnostics complete.")
  invisible(norm)
}
