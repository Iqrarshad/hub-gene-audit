# Five tests quantifying hub selection bias.
#
# 1. Logistic prediction of hub status from global interactome degree,
#    reported as AUC.
# 2. Rank agreement restricted to genes with any centrality rank, since the
#    all-gene coefficient is dominated by ties at zero.
# 3. Top-k overlap between the conventional and corrected rankings.
# 4. Recurrence: hub selection run separately on each dataset's DEG set.
# 5. Degree-preserving rewiring: a gene that remains a hub after rewiring
#    owes its centrality to degree alone.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC); library(httr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

N_REWIRE <- 200

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

centrality_tab <- function(g) {
  vs <- V(g)$name
  tibble(gene = vs,
         degree = degree(g),
         betweenness = betweenness(g, normalized = TRUE),
         closeness = closeness(g, normalized = TRUE),
         eigenvector = eigen_centrality(g)$vector,
         pagerank = page_rank(g)$vector,
         mcc = mcc_score(g)[vs])
}

consensus_of <- function(tab, min_methods = 4) {
  mets <- c("degree", "betweenness", "closeness", "eigenvector",
            "pagerank", "mcc")
  top10 <- vapply(mets, function(m)
    as.integer(tab$gene %in% tab$gene[head(order(tab[[m]],
                                                 decreasing = TRUE), 10)]),
    integer(nrow(tab)))
  tab$n_methods <- rowSums(top10)
  tab$gene[tab$n_methods >= min_methods]
}

string_edges <- function(genes, score, tag) {
  cd <- file.path(CACHE_DIR, "string_bias"); dir.create(cd, showWarnings = FALSE,
                                                        recursive = TRUE)
  f <- file.path(cd, paste0(tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- tryCatch(httr::POST("https://string-db.org/api/tsv/network",
        body = list(identifiers = paste(genes, collapse = "%0d"),
                    species = 9606, required_score = score,
                    caller_identity = "glioma_hub_bias"), encode = "form"),
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

main_12 <- function() {
  net  <- readRDS(P("rds", "network.rds"))
  score <- as.integer(net$chosen)
  conv_hubs <- net$consensus$gene
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene

  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  sc <- read_csv(P("tables", "specific_connectivity_binomial.csv"),
                 show_col_types = FALSE)

  # ---- TEST 1: is hub status predictable from global degree alone? ------
  d <- norm %>%
    mutate(is_hub = as.integer(gene %in% conv_hubs),
           log_global = log10(global_degree + 1)) %>%
    filter(!is.na(global_degree))

  log_msg("=== TEST 1: predictability from global interactome degree ===")
  log_msg("Genes: ", nrow(d), "; conventional hubs among them: ",
          sum(d$is_hub))

  if (sum(d$is_hub) >= 3 && sum(d$is_hub) < nrow(d)) {
    fit <- glm(is_hub ~ log_global, data = d, family = binomial())
    pred <- predict(fit, type = "response")
    roc1 <- pROC::roc(d$is_hub, pred, quiet = TRUE)
    auc1 <- as.numeric(pROC::auc(roc1))
    ci1 <- as.numeric(pROC::ci.auc(roc1))
    co <- summary(fit)$coefficients

    log_msg("AUC of global degree predicting conventional hub status: ",
            round(auc1, 3), " (95% CI ", round(ci1[1], 3), " to ",
            round(ci1[3], 3), ")")
    log_msg("  log10 global degree: OR per unit = ",
            round(exp(co["log_global", "Estimate"]), 2),
            ", p = ", signif(co["log_global", "Pr(>|z|)"], 3))

    # same test for the binomial-significant set
    d2 <- d %>% mutate(is_sig = as.integer(gene %in%
                       sc$gene[sc$significant]))
    auc2 <- NA_real_
    if (sum(d2$is_sig) >= 3) {
      fit2 <- glm(is_sig ~ log_global, data = d2, family = binomial())
      auc2 <- as.numeric(pROC::auc(pROC::roc(d2$is_sig,
                predict(fit2, type = "response"), quiet = TRUE)))
      log_msg("AUC of global degree predicting BINOMIAL-significant status: ",
              round(auc2, 3))
      log_msg("  The corrected method should be far less predictable from ",
              "global degree than the conventional one.")
    }
    write_csv(tibble(target = c("conventional hub", "binomial significant"),
                     auc = c(auc1, auc2),
                     ci_lo = c(ci1[1], NA), ci_hi = c(ci1[3], NA)),
              P("tables", "hubbias_predictability.csv"))
  } else {
    log_msg("Too few hubs for a logistic fit.")
  }

  # ---- TEST 2: agreement among ranked genes only ------------------------
  log_msg("=== TEST 2: agreement excluding the zero mass ===")
  j <- norm %>% inner_join(sc %>% select(gene, p_binom), by = "gene")
  all_rho <- suppressWarnings(cor(j$n_methods_top10,
                                  -log10(j$p_binom + 1e-300),
                                  method = "spearman"))
  ranked <- j %>% filter(n_methods_top10 > 0)
  rk_rho <- if (nrow(ranked) > 5)
    suppressWarnings(cor(ranked$n_methods_top10,
                         -log10(ranked$p_binom + 1e-300),
                         method = "spearman")) else NA_real_
  log_msg("Spearman across all ", nrow(j), " genes: ", round(all_rho, 3))
  log_msg("Spearman across the ", nrow(ranked),
          " genes with any centrality rank: ", round(rk_rho, 3))
  log_msg("Genes with consensus count zero: ",
          sum(j$n_methods_top10 == 0), " of ", nrow(j),
          ". The all-gene coefficient is dominated by these ties.")

  # ---- TEST 3: top-k overlap curve --------------------------------------
  log_msg("=== TEST 3: top-k overlap ===")
  conv_rank <- norm %>% arrange(desc(n_methods_top10), desc(local_degree)) %>%
    pull(gene)
  spec_rank <- sc %>% arrange(p_binom) %>% pull(gene)
  ks <- c(5, 10, 20, 50)
  ov <- map_dfr(ks, function(k) {
    a <- head(conv_rank, k); b <- head(spec_rank, k)
    tibble(k = k, overlap = length(intersect(a, b)),
           jaccard = length(intersect(a, b)) / length(union(a, b)),
           expected_by_chance = round(k * k / nrow(norm), 2))
  })
  write_csv(ov, P("tables", "hubbias_topk_overlap.csv"))
  print(as.data.frame(ov))

  # ---- TEST 4: recurrence across different input gene sets --------------
  log_msg("=== TEST 4: do the same hubs appear from different inputs? ===")
  micro <- readRDS(P("rds", "deg_microarray.rds"))
  rna   <- readRDS(P("rds", "deg_rnaseq.rds"))

  input_sets <- list()
  for (nm in names(micro)) {
    x <- micro[[nm]]; if (is.null(x)) next
    g <- if ("composition_robust" %in% names(x) &&
             any(x$composition_robust, na.rm = TRUE))
      x$gene[which(x$composition_robust)] else
      x$gene[which(x$sig_raw)]
    if (length(g) >= 60) input_sets[[nm]] <- head(g, 148)
  }
  for (nm in names(rna)) {
    x <- rna[[nm]]
    g <- x$gene[which(x$composition_robust)]
    if (length(g) >= 60) input_sets[[paste0("GSE147352_", nm)]] <- head(g, 148)
  }

  log_msg("Independent input sets: ", length(input_sets))
  rec <- list()
  for (nm in names(input_sets)) {
    e <- string_edges(input_sets[[nm]], score, paste0("input_", make.names(nm)))
    if (is.null(e)) { log_msg("  ", nm, ": STRING returned nothing"); next }
    g <- graph_from_data_frame(e, directed = FALSE)
    miss <- setdiff(input_sets[[nm]], V(g)$name)
    if (length(miss)) g <- add_vertices(g, length(miss), name = miss)
    h <- consensus_of(centrality_tab(g))
    log_msg("  ", nm, " -> ", paste(head(h, 10), collapse = ", "))
    rec[[nm]] <- h
    Sys.sleep(0.5)
  }

  if (length(rec) >= 2) {
    tabf <- table(unlist(rec))
    recur <- tibble(gene = names(tabf), n_sets = as.integer(tabf),
                    pct_sets = round(100 * as.integer(tabf) / length(rec), 1)) %>%
      arrange(desc(n_sets))
    write_csv(recur, P("tables", "hubbias_recurrence.csv"))
    log_msg("Genes appearing as hubs across independent input sets:")
    print(as.data.frame(head(recur, 15)))

    # Jaccard between the hub sets from different inputs
    combs <- combn(names(rec), 2, simplify = FALSE)
    jac <- map_dbl(combs, function(p) {
      a <- rec[[p[1]]]; b <- rec[[p[2]]]
      if (!length(union(a, b))) return(NA_real_)
      length(intersect(a, b)) / length(union(a, b))
    })
    log_msg("Median Jaccard between hub sets from DIFFERENT gene lists: ",
            round(median(jac, na.rm = TRUE), 3))
    log_msg("  High values mean hub identity is insensitive to the input.")
  }

  # ---- TEST 5: degree-preserving rewiring -------------------------------
  log_msg("=== TEST 5: degree-preserving rewiring ===")
  e0 <- string_edges(genes, score, "observed")
  if (!is.null(e0)) {
    g0 <- graph_from_data_frame(e0, directed = FALSE)
    miss <- setdiff(genes, V(g0)$name)
    if (length(miss)) g0 <- add_vertices(g0, length(miss), name = miss)
    obs_h <- consensus_of(centrality_tab(g0))

    hits <- setNames(integer(length(obs_h)), obs_h)
    for (i in seq_len(N_REWIRE)) {
      gr <- rewire(g0, keeping_degseq(niter = ecount(g0) * 10))
      h <- consensus_of(centrality_tab(gr))
      for (x in intersect(obs_h, h)) hits[x] <- hits[x] + 1L
    }
    rw <- tibble(gene = names(hits), n_rewired_hub = as.integer(hits),
                 pct = round(100 * as.integer(hits) / N_REWIRE, 1)) %>%
      arrange(desc(pct))
    write_csv(rw, P("tables", "hubbias_rewiring.csv"))
    log_msg("Hub status retained after degree-preserving rewiring:")
    print(as.data.frame(rw))
    log_msg("  A gene that remains a hub after rewiring owes its centrality ",
            "to its degree alone, not to which genes it connects to.")
  }

  log_msg("12_hub_bias_quantified complete.")
}
