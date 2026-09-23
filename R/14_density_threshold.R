# Does network density determine whether hub inference is possible?
#
# Density is varied two ways, by STRING confidence threshold and by gene
# set size, and hub list stability is measured at each level.
#
# Stability is the Jaccard between the hub list from a random half of the
# genes and the hub list from the full set. Comparing two disjoint halves
# to each other gives zero overlap by construction.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(httr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

N_SPLITS <- 20        # random half-splits per condition
API_DELAY <- 0.4

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

consensus_hubs <- function(g, top_k = 10, min_methods = 4) {
  if (ecount(g) == 0) return(character(0))
  vs <- V(g)$name
  tab <- tibble(
    gene = vs, degree = degree(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness = closeness(g, normalized = TRUE),
    eigenvector = eigen_centrality(g)$vector,
    pagerank = page_rank(g)$vector,
    mcc = mcc_score(g)[vs])
  mets <- setdiff(names(tab), "gene")
  top <- vapply(mets, function(m)
    as.integer(tab$gene %in% tab$gene[head(order(tab[[m]],
                                                 decreasing = TRUE), top_k)]),
    integer(nrow(tab)))
  tab$gene[rowSums(top) >= min_methods]
}

string_net <- function(genes, score, tag) {
  cd <- file.path(CACHE_DIR, "string_density")
  dir.create(cd, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cd, paste0(tag, ".rds"))
  if (file.exists(f)) {
    cached <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(cached) && nrow(cached)) return(cached)
  }

  # Reuse the stage-06 edge cache when it already covers this score. By the time
  # this stage runs, stages 06-13 have made many STRING calls, so a fresh POST is
  # often rate-limited; reusing the cache avoids the call entirely for 400/700.
  rf <- P("rds", paste0("string_", score, ".rds"))
  if (file.exists(rf)) {
    e0 <- tryCatch(readRDS(rf), error = function(err) NULL)
    if (!is.null(e0) && nrow(e0)) {
      e0 <- e0 %>% select(from, to) %>% filter(from != to) %>% distinct()
      saveRDS(e0, f); return(e0)
    }
  }

  fetch_once <- function() {
    r <- tryCatch(httr::POST("https://string-db.org/api/tsv/network",
          body = list(identifiers = paste(genes, collapse = "%0d"),
                      species = 9606, required_score = score,
                      caller_identity = "hub_gene_audit"), encode = "form"),
          error = function(e) NULL)
    if (is.null(r) || httr::status_code(r) != 200) return(NULL)
    d <- tryCatch(readr::read_tsv(I(httr::content(r, as = "text",
                                                  encoding = "UTF-8")),
                                  show_col_types = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(NULL)
    d %>% select(from = preferredName_A, to = preferredName_B) %>%
      filter(from != to) %>% distinct()
  }

  # Retry with exponential backoff to wait out STRING rate limiting.
  e <- NULL
  for (attempt in seq_len(6)) {
    e <- fetch_once()
    if (!is.null(e) && nrow(e)) break
    wait <- min(60, 5 * 2^(attempt - 1))   # 5, 10, 20, 40, 60, 60 s
    log_msg("  STRING score ", score, ": no edges (attempt ", attempt,
            "), backing off ", wait, "s")
    Sys.sleep(wait)
  }
  if (!is.null(e) && nrow(e)) saveRDS(e, f)
  e
}

build_graph <- function(edges, genes) {
  if (is.null(edges) || nrow(edges) == 0) {
    g <- make_empty_graph(n = 0, directed = FALSE)
    g <- add_vertices(g, length(genes), name = genes)
    return(g)
  }
  edges <- edges %>% filter(from %in% genes, to %in% genes)
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = genes)
  g
}

net_stats <- function(g) {
  comp <- components(g)
  tibble(n_nodes = vcount(g), n_edges = ecount(g),
         mean_degree = if (vcount(g)) mean(degree(g)) else 0,
         density = edge_density(g),
         pct_isolated = 100 * mean(degree(g) == 0),
         lcc_frac = if (vcount(g)) max(comp$csize) / vcount(g) else 0)
}

# Stability: Jaccard between hub lists from random halves
stability <- function(edges, genes, n_splits = N_SPLITS) {
  js <- numeric(0); sizes <- numeric(0)
  for (i in seq_len(n_splits)) {
    set.seed(SEED + i)
    half <- sample(genes, floor(length(genes) / 2))
    other <- setdiff(genes, half)
    ga <- build_graph(edges, half); gb <- build_graph(edges, other)
    ha <- consensus_hubs(ga); hb <- consensus_hubs(gb)
    sizes <- c(sizes, length(ha), length(hb))
    # Overlap between hub sets from disjoint halves is necessarily zero by
    # construction, so instead compare each half against the FULL network's
    # hub list: does a half recover what the whole finds?
    hf <- consensus_hubs(build_graph(edges, genes))
    if (!length(hf)) next
    ja <- length(intersect(ha, hf)) / length(union(ha, hf))
    jb <- length(intersect(hb, hf)) / length(union(hb, hf))
    js <- c(js, ja, jb)
  }
  tibble(stability = if (length(js)) mean(js, na.rm = TRUE) else NA_real_,
         stability_sd = if (length(js)) sd(js, na.rm = TRUE) else NA_real_,
         mean_hubs_per_half = mean(sizes))
}

main_14 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes148 <- inter$shared$gene

  # Larger gene sets for the size sweep, drawn from the conventional
  # intersection so the biology is comparable
  conv <- inter$conventional
  big_pool <- if (!is.null(conv) && "gene" %in% names(conv)) conv$gene else
    genes148
  log_msg("Gene pools: composition-robust ", length(genes148),
          ", conventional ", length(big_pool))

  results <- list()

  # ---- Sweep 1: STRING confidence threshold ----------------------------
  log_msg("=== SWEEP 1: STRING confidence threshold ===")
  for (sc in c(150, 250, 400, 700, 900)) {
    e <- string_net(genes148, sc, paste0("t", sc))
    g <- build_graph(e, genes148)
    st <- net_stats(g)
    stb <- stability(e, genes148)
    h <- consensus_hubs(g)
    results[[length(results) + 1]] <- bind_cols(
      tibble(sweep = "confidence", condition = paste0("score ", sc),
             n_genes = length(genes148)),
      st, stb,
      tibble(n_hubs = length(h),
             hubs = paste(head(h, 8), collapse = "; ")))
    log_msg("  score ", sc, ": ", st$n_edges, " edges, mean degree ",
            round(st$mean_degree, 2), ", stability ",
            round(stb$stability, 3))
    Sys.sleep(API_DELAY)
  }

  # ---- Sweep 2: gene set size at fixed threshold -----------------------
  log_msg("=== SWEEP 2: gene set size at score 400 ===")
  for (n in c(50, 100, 148, 300, 600)) {
    if (n > length(big_pool)) next
    set.seed(SEED)
    gs <- if (n <= length(genes148)) head(genes148, n) else
      unique(c(genes148, sample(setdiff(big_pool, genes148),
                                n - length(genes148))))
    e <- string_net(gs, 400, paste0("n", n))
    g <- build_graph(e, gs)
    st <- net_stats(g)
    stb <- stability(e, gs)
    h <- consensus_hubs(g)
    results[[length(results) + 1]] <- bind_cols(
      tibble(sweep = "set size", condition = paste0("n = ", n),
             n_genes = length(gs)),
      st, stb,
      tibble(n_hubs = length(h),
             hubs = paste(head(h, 8), collapse = "; ")))
    log_msg("  n = ", n, ": ", st$n_edges, " edges, mean degree ",
            round(st$mean_degree, 2), ", stability ",
            round(stb$stability, 3))
    Sys.sleep(API_DELAY)
  }

  out <- bind_rows(results)
  write_csv(out, P("tables", "density_threshold_sweep.csv"))

  log_msg("=================================================")
  log_msg("DENSITY AND HUB STABILITY")
  log_msg("  stability = mean Jaccard between the hub list from a random")
  log_msg("  half of the genes and the hub list from the full set.")
  log_msg("  Low stability means hub identity is not reproducible.")
  log_msg("=================================================")
  print(as.data.frame(out %>% select(sweep, condition, n_nodes, n_edges,
                                     mean_degree, pct_isolated, lcc_frac,
                                     n_hubs, stability) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))))

  # ---- Where is the boundary? -------------------------------------------
  d <- out %>% filter(!is.na(stability), n_edges > 0)
  if (nrow(d) >= 4) {
    rho <- suppressWarnings(cor(d$mean_degree, d$stability,
                                method = "spearman"))
    log_msg("Spearman, mean degree vs hub stability: ", round(rho, 3))

    lowd  <- d %>% filter(mean_degree < 2)
    highd <- d %>% filter(mean_degree >= 3)
    if (nrow(lowd) && nrow(highd)) {
      log_msg("Mean stability at mean degree < 2 : ",
              round(mean(lowd$stability, na.rm = TRUE), 3),
              " (", nrow(lowd), " conditions)")
      log_msg("Mean stability at mean degree >= 3: ",
              round(mean(highd$stability, na.rm = TRUE), 3),
              " (", nrow(highd), " conditions)")
    }

    # First condition where stability exceeds 0.5
    ok <- d %>% arrange(mean_degree) %>% filter(stability >= 0.5)
    if (nrow(ok)) {
      log_msg("Stability first reaches 0.5 at mean degree ",
              round(ok$mean_degree[1], 2), " (", ok$condition[1], ")")
      log_msg("PRACTICAL GUIDELINE: report mean degree and isolated-node ",
              "fraction before ranking hubs. Below mean degree ",
              round(ok$mean_degree[1], 1),
              ", hub identity is not reproducible and should not be reported.")
    } else {
      log_msg("Stability never reaches 0.5 in any condition tested. Hub ",
              "identity is unreliable across the whole range examined, ",
              "which is a stronger statement and must be reported as such.")
    }
  }

  saveRDS(out, P("rds", "density_threshold.rds"))
  log_msg("14_density_threshold complete.")
  invisible(out)
}
