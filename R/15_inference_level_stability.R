# Stability of inference at three levels of abstraction, using the same
# split-half design: gene-level hubs, pathway enrichment, and cross-cohort
# expression replication.
#
# If the failure is specific to gene-level topological inference, pathway
# and expression-level inference on the same gene sets should be more
# stable. If all levels are comparable, the instability is a property of
# the gene sets rather than of the topological step.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(clusterProfiler); library(org.Hs.eg.db)
  library(igraph); library(httr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

N_SPLITS <- 20

jaccard <- function(a, b) {
  if (!length(union(a, b))) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# --- Level 1: gene-level hubs (recomputed here for a like-for-like base) -
hub_stability <- function(edges, genes, n_splits = N_SPLITS) {
  mcc <- function(g) {
    cl <- max_cliques(g, min = 1)
    s <- setNames(numeric(vcount(g)), V(g)$name)
    for (c in cl) { k <- length(c)
      if (k > 1) for (v in names(c)) s[v] <- s[v] + factorial(k - 1) }
    s[s == 0] <- 1; s
  }
  hubs_of <- function(gs) {
    e <- edges %>% filter(from %in% gs, to %in% gs)
    g <- graph_from_data_frame(e, directed = FALSE, vertices = gs)
    if (ecount(g) == 0) return(character(0))
    vs <- V(g)$name
    tb <- tibble(gene = vs, degree = degree(g),
                 betweenness = betweenness(g, normalized = TRUE),
                 closeness = closeness(g, normalized = TRUE),
                 eigenvector = eigen_centrality(g)$vector,
                 pagerank = page_rank(g)$vector, mcc = mcc(g)[vs])
    mets <- setdiff(names(tb), "gene")
    top <- vapply(mets, function(m)
      as.integer(tb$gene %in% tb$gene[head(order(tb[[m]], decreasing = TRUE),
                                           10)]), integer(nrow(tb)))
    tb$gene[rowSums(top) >= 4]
  }
  full <- hubs_of(genes)
  js <- map_dbl(seq_len(n_splits), function(i) {
    set.seed(SEED + i)
    jaccard(hubs_of(sample(genes, floor(length(genes) / 2))), full)
  })
  tibble(level = "gene-level hubs (topology)", n_full = length(full),
         stability = mean(js, na.rm = TRUE), sd = sd(js, na.rm = TRUE))
}

# --- Level 2: pathway enrichment ----------------------------------------
pathway_stability <- function(genes, universe, n_splits = N_SPLITS) {
  to_entrez <- function(g) {
    suppressWarnings(suppressMessages(
      bitr(g, "SYMBOL", "ENTREZID", org.Hs.eg.db)))$ENTREZID
  }
  uni <- to_entrez(universe)

  enrich_of <- function(gs) {
    ez <- to_entrez(gs)
    if (length(ez) < 10) return(character(0))
    k <- tryCatch(enrichKEGG(ez, universe = uni, organism = "hsa",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                             pAdjustMethod = "BH"),
                  error = function(e) NULL)
    if (is.null(k)) return(character(0))
    d <- as.data.frame(k)
    if (!nrow(d)) return(character(0))
    d$ID
  }

  full <- enrich_of(genes)
  log_msg("  full gene set: ", length(full), " KEGG pathways")
  if (!length(full)) {
    return(tibble(level = "pathway enrichment (KEGG)", n_full = 0,
                  stability = NA_real_, sd = NA_real_))
  }
  js <- map_dbl(seq_len(n_splits), function(i) {
    set.seed(SEED + i)
    jaccard(enrich_of(sample(genes, floor(length(genes) / 2))), full)
  })
  tibble(level = "pathway enrichment (KEGG)", n_full = length(full),
         stability = mean(js, na.rm = TRUE), sd = sd(js, na.rm = TRUE))
}

# --- Level 3: cross-cohort expression replication -----------------------
# Not a gene-set split: this asks whether the SUPPORTED associations from
# one cohort are recovered in another. The comparable quantity is the
# proportion of findings that replicate, which is the reproducibility of
# expression-level inference.
expression_stability <- function() {
  f <- P("tables", "validation_correlations_all.csv")
  if (!file.exists(f)) return(tibble(level = "expression replication",
                                     n_full = NA, stability = NA_real_,
                                     sd = NA_real_))
  val <- read_csv(f, show_col_types = FALSE) %>%
    filter(method == DECONV_PRIMARY)

  cohorts <- unique(val$cohort)
  if (length(cohorts) < 2) return(tibble(level = "expression replication",
                                         n_full = NA, stability = NA_real_,
                                         sd = NA_real_))

  sig_of <- function(co) {
    val %>% filter(cohort == co, !is.na(q), q < THRESH$corr_qval,
                   abs(rho) >= THRESH$corr_min_abs_r) %>%
      mutate(key = paste(grade, gene, cell, sep = "|")) %>% pull(key)
  }
  sets <- map(cohorts, sig_of); names(sets) <- cohorts
  sets <- sets[lengths(sets) > 0]
  if (length(sets) < 2) return(tibble(level = "expression replication",
                                      n_full = NA, stability = NA_real_,
                                      sd = NA_real_))

  combs <- combn(names(sets), 2, simplify = FALSE)
  js <- map_dbl(combs, function(p) jaccard(sets[[p[1]]], sets[[p[2]]]))
  log_msg("  pairwise cohort Jaccard: ",
          paste(round(js, 3), collapse = ", "))
  tibble(level = "expression replication (cross-cohort)",
         n_full = round(mean(lengths(sets))),
         stability = mean(js, na.rm = TRUE), sd = sd(js, na.rm = TRUE))
}

main_15 <- function() {
  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared$gene
  net <- readRDS(P("rds", "network.rds"))
  score <- as.integer(net$chosen)

  # Edges for the observed network
  ef <- file.path(CACHE_DIR, "string_density", paste0("t", score, ".rds"))
  edges <- if (file.exists(ef)) readRDS(ef) else NULL
  if (is.null(edges) || !nrow(edges))
    edges <- readRDS(P("rds", paste0("string_", score, ".rds")))
  if (is.null(edges) || !nrow(edges))
    stop("No cached STRING edges; run the network stage (06) first.")

  log_msg("Gene set: ", length(genes), "; edges: ", nrow(edges))

  log_msg("--- level 1: gene-level hubs ---")
  l1 <- hub_stability(edges, genes)

  log_msg("--- level 2: pathway enrichment ---")
  micro <- readRDS(P("rds", "deg_microarray.rds"))
  rna   <- readRDS(P("rds", "deg_rnaseq.rds"))
  universe <- unique(c(unlist(lapply(micro, function(x)
                        if (!is.null(x)) x$gene)),
                      rna$LGG_vs_Normal$gene))
  universe <- universe[!is.na(universe) & universe != ""]
  l2 <- pathway_stability(genes, universe)

  # Also the grade contrast, which is where the immune signal actually was
  hgg <- rna$HGG_vs_LGG
  l2b <- NULL
  if (!is.null(hgg)) {
    gg <- hgg$gene[which(hgg$composition_robust)]
    if (length(gg) > 200) {
      log_msg("--- level 2b: pathway enrichment, grade contrast (",
              length(gg), " genes) ---")
      l2b <- pathway_stability(gg, universe) %>%
        mutate(level = "pathway enrichment (grade contrast)")
    }
  }

  log_msg("--- level 3: cross-cohort expression replication ---")
  l3 <- expression_stability()

  out <- bind_rows(l1, l2, l2b, l3)
  write_csv(out, P("tables", "inference_level_stability.csv"))

  log_msg("=================================================")
  log_msg("STABILITY BY LEVEL OF INFERENCE")
  log_msg("  Jaccard against the full-data result after removing half the")
  log_msg("  input genes, except for expression replication which is the")
  log_msg("  agreement between independent patient cohorts.")
  log_msg("=================================================")
  print(as.data.frame(out %>% mutate(across(where(is.numeric),
                                            ~ round(.x, 3)))))

  g1 <- out$stability[out$level == "gene-level hubs (topology)"][1]
  gp <- out$stability[grepl("^pathway", out$level)]
  gp <- gp[!is.na(gp)]

  if (length(gp) && !is.na(g1)) {
    if (max(gp) > g1 + 0.2) {
      log_msg("CONCLUSION: pathway-level inference is substantially more ",
              "stable than gene-level hub selection (", round(max(gp), 3),
              " versus ", round(g1, 3), "). The recommendation is to report ",
              "enriched pathways rather than ranked hub genes, and to ",
              "validate individual genes by cross-cohort expression rather ",
              "than by network position.")
    } else if (max(gp) > g1) {
      log_msg("CONCLUSION: pathway inference is somewhat more stable (",
              round(max(gp), 3), " versus ", round(g1, 3),
              ") but not dramatically so. Report both figures and temper ",
              "the recommendation accordingly.")
    } else {
      log_msg("CONCLUSION: pathway inference is NOT more stable than gene ",
              "level. The instability is a property of the gene sets ",
              "themselves, not of the topological step. This is a broader ",
              "and less comfortable finding and must be reported as such.")
    }
  }

  saveRDS(out, P("rds", "inference_level_stability.rds"))
  log_msg("15_inference_level_stability complete.")
  invisible(out)
}
