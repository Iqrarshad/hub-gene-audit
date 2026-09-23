# Null model for hub selection, using a spike-in design.
#
# Each observed hub is placed into an otherwise random background of the
# same size. If a gene becomes a hub against random company, its hub status
# reflects its position in the interactome rather than the gene set.
#
# Drawing random sets and asking how often they contain a given gene does
# not work at this scale: with 148 genes from a 32,000-gene universe each
# gene appears in under 1% of draws.

suppressPackageStartupMessages({
  library(igraph); library(dplyr); library(readr); library(tibble)
  library(httr); library(purrr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

N_RANDOM <- 50          # STRING API calls; keep modest and polite
API_DELAY <- 1.0        # seconds between calls

# --- STRING query with caching ------------------------------------------
string_network <- function(genes, score, tag) {
  cache_dir <- file.path(CACHE_DIR, "string_null")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cache_dir, paste0(tag, "_", score, ".rds"))
  if (file.exists(f)) return(readRDS(f))

  r <- tryCatch(httr::POST(
    "https://string-db.org/api/tsv/network",
    body = list(identifiers = paste(genes, collapse = "%0d"),
                species = 9606, required_score = score,
                caller_identity = "glioma_hub_null"),
    encode = "form"), error = function(e) NULL)

  if (is.null(r) || httr::status_code(r) != 200) {
    saveRDS(NULL, f); return(NULL)
  }
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  e <- tryCatch(readr::read_tsv(I(txt), show_col_types = FALSE),
                error = function(e) NULL)
  if (is.null(e) || nrow(e) == 0) { saveRDS(NULL, f); return(NULL) }

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

consensus_hubs <- function(edges, all_genes) {
  if (is.null(edges) || nrow(edges) == 0) return(character(0))
  g <- graph_from_data_frame(edges, directed = FALSE)
  miss <- setdiff(all_genes, V(g)$name)
  if (length(miss)) g <- add_vertices(g, length(miss), name = miss)

  vs <- V(g)$name
  m <- tibble(gene = vs,
              degree = degree(g),
              betweenness = betweenness(g, normalized = TRUE),
              closeness = closeness(g, normalized = TRUE),
              eigenvector = eigen_centrality(g)$vector,
              pagerank = page_rank(g)$vector,
              mcc = mcc_score(g)[vs])

  methods <- c("degree", "betweenness", "closeness", "eigenvector",
               "pagerank", "mcc")
  top10 <- vapply(methods, function(mth) {
    as.integer(m$gene %in% m$gene[head(order(m[[mth]], decreasing = TRUE), 10)])
  }, integer(nrow(m)))
  m$n <- rowSums(top10)
  m$gene[m$n >= 4]
}

main_07 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  observed_genes <- inter$shared$gene
  n_set <- length(observed_genes)

  net <- readRDS(P("rds", "network.rds"))
  score <- as.integer(net$chosen)
  obs_hubs <- net$consensus$gene
  log_msg("Observed hubs at score ", score, ": ",
          paste(obs_hubs, collapse = ", "))

  micro <- readRDS(P("rds", "deg_microarray.rds"))
  rna   <- readRDS(P("rds", "deg_rnaseq.rds"))
  universe <- unique(c(unlist(lapply(micro, function(x)
                         if (!is.null(x)) x$gene)),
                       rna$LGG_vs_Normal$gene))
  universe <- universe[!is.na(universe) & universe != ""]
  log_msg("Universe: ", length(universe), " genes; set size ", n_set)

  # --- SPIKE-IN DESIGN ---------------------------------------------------
  # Drawing random sets and asking how often they happen to contain a given
  # gene fails: at 148 genes from a 32,792-gene universe, each gene appears
  # in roughly 0.5% of draws, so the conditional frequency is undefined.
  #
  # Instead each observed hub is placed INTO an otherwise random background
  # of the same size. If the gene becomes a consensus hub against random
  # company, its hub status reflects its position in the interactome rather
  # than anything about this gene set.
  N_BG <- 25

  res <- map_dfr(obs_hubs, function(g) {
    log_msg("Spike-in test: ", g)
    hits <- 0L; valid <- 0L
    for (i in seq_len(N_BG)) {
      set.seed(SEED + 1000 * match(g, obs_hubs) + i)
      bg <- sample(setdiff(universe, g), n_set - 1)
      gs <- c(g, bg)
      e <- string_network(gs, score, paste0("spike_", g, "_", i))
      if (is.null(e)) next
      valid <- valid + 1L
      if (g %in% consensus_hubs(e, gs)) hits <- hits + 1L
      Sys.sleep(API_DELAY)
    }
    log_msg("  ", g, ": hub in ", hits, " of ", valid, " random backgrounds")
    tibble(gene = g, n_backgrounds = valid, n_times_hub = hits,
           null_hub_frequency = if (valid > 0) hits / valid else NA_real_)
  })

  res <- res %>%
    mutate(
      condition_specific = !is.na(null_hub_frequency) &
                           null_hub_frequency < 0.10,
      interpretation = case_when(
        is.na(null_hub_frequency)   ~ "test failed",
        null_hub_frequency >= 0.50  ~ "interactome hub, not a finding",
        null_hub_frequency >= 0.10  ~ "partly explained by study bias",
        TRUE                        ~ "condition-specific")) %>%
    arrange(desc(null_hub_frequency))

  write_csv(res, P("tables", "hub_null_spikein.csv"))

  log_msg("=================================================")
  log_msg("HUB NULL MODEL, SPIKE-IN DESIGN")
  log_msg("  null_hub_frequency = proportion of RANDOM 148-gene backgrounds")
  log_msg("  in which this gene still emerges as a consensus hub")
  log_msg("=================================================")
  print(as.data.frame(res))

  n_spec <- sum(res$condition_specific, na.rm = TRUE)
  log_msg(n_spec, " of ", nrow(res), " hubs are condition-specific")
  if (n_spec < nrow(res)) {
    log_msg("Explained by interactome study bias: ",
            paste(res$gene[!res$condition_specific], collapse = ", "))
  }

  saveRDS(list(observed = obs_hubs, spikein = res), P("rds", "hub_null.rds"))
  log_msg("07_hub_null complete.")
  invisible(res)
}
