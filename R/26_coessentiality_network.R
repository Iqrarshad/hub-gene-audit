# Network built from CRISPR gene effect rather than from literature.
#
# Two genes whose knockout effects correlate across cell lines are
# functionally related, and no publication record influences that
# measurement.
#
# The network is built twice, once on all lines and once with CNS lines
# held out, because the evaluation criterion is also derived from DepMap.
# The held-out version is the one to interpret.
#
# Requires CRISPRGeneEffect.csv and Model.csv from depmap.org.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(igraph); library(pROC); library(data.table)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; distinct <- dplyr::distinct
arrange <- dplyr::arrange; summarise <- dplyr::summarise

N_NOMINATE  <- 8
EDGE_BUDGET <- 139

mcc_score <- function(g) {
  cl <- max_cliques(g, min = 1)
  s <- setNames(numeric(vcount(g)), V(g)$name)
  for (c in cl) { k <- length(c)
    if (k > 1) for (v in names(c)) s[v] <- s[v] + factorial(k - 1) }
  s[s == 0] <- 1; s
}

consensus_hubs <- function(g, top_k = 10, min_methods = 4) {
  if (ecount(g) == 0) return(character(0))
  vs <- V(g)$name
  tab <- tibble(gene = vs, degree = degree(g),
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
  tab$n <- rowSums(top)
  tab$gene[tab$n >= min_methods]
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

main_26 <- function() {
  # --- DepMap ------------------------------------------------------------
  ge <- find_local("CRISPRGeneEffect")
  md <- find_local("Model")
  if (is.na(ge) || is.na(md))
    stop("CRISPRGeneEffect.csv and Model.csv required in ", DATA_DIR)

  log_msg("Reading DepMap gene effect matrix ...")
  eff <- data.table::fread(ge, data.table = FALSE, check.names = FALSE)
  rownames(eff) <- eff[[1]]; eff <- eff[, -1, drop = FALSE]
  colnames(eff) <- sub(" .*$", "", colnames(eff))
  mod <- data.table::fread(md, data.table = FALSE, check.names = FALSE)

  idc <- grep("ModelID|DepMap_ID", names(mod), value = TRUE)[1]
  linc <- grep("OncotreeLineage|lineage|primary_disease", names(mod),
               value = TRUE)[1]
  lin <- setNames(as.character(mod[[linc]]), mod[[idc]])
  cls <- rownames(eff)
  is_cns <- grepl("CNS|Brain|Glio", lin[cls], ignore.case = TRUE)
  is_cns[is.na(is_cns)] <- FALSE
  log_msg("Cell lines: ", length(cls), "; CNS/brain: ", sum(is_cns))

  inter <- readRDS(P("rds", "intersect.rds"))
  genes <- inter$shared; if (is.data.frame(genes)) genes <- genes$gene
  g_in <- intersect(genes, colnames(eff))
  log_msg("Genes present in DepMap: ", length(g_in), " of ", length(genes))
  if (length(g_in) < 50) stop("Too few genes in DepMap for a network.")

  E <- as.matrix(eff[, g_in, drop = FALSE])
  storage.mode(E) <- "numeric"

  # --- Co-essentiality network -------------------------------------------
  # Held-out construction: CNS lines excluded, so the selectivity evaluation
  # is not computed on the same cell lines that built the network.
  build_net <- function(rows, label) {
    M <- E[rows, , drop = FALSE]
    keep <- apply(M, 2, function(v) sum(!is.na(v)) > 50 &&
                    sd(v, na.rm = TRUE) > 0)
    M <- M[, keep, drop = FALSE]
    C <- suppressWarnings(cor(M, use = "pairwise.complete.obs",
                              method = "pearson"))
    C[!is.finite(C)] <- 0; diag(C) <- 0
    ut <- which(upper.tri(C), arr.ind = TRUE)
    v <- abs(C[ut])
    b <- min(EDGE_BUDGET, length(v))
    cut <- sort(v, decreasing = TRUE)[b]
    sel <- ut[v >= cut, , drop = FALSE]
    el <- tibble(from = colnames(C)[sel[, 1]], to = colnames(C)[sel[, 2]])
    gg <- graph_from_data_frame(el, directed = FALSE,
                                vertices = colnames(C))
    log_msg(label, ": ", vcount(gg), " nodes, ", ecount(gg),
            " edges at |r| >= ", round(cut, 3))
    gg
  }

  g_all <- build_net(seq_len(nrow(E)), "co-essentiality, all lines")
  g_ho  <- build_net(which(!is_cns), "co-essentiality, CNS held out")

  hubs_all <- consensus_hubs(g_all)
  hubs_ho  <- consensus_hubs(g_ho)
  log_msg("Co-essentiality hubs (all lines): ",
          paste(hubs_all, collapse = ", "))
  log_msg("Co-essentiality hubs (CNS held out): ",
          paste(hubs_ho, collapse = ", "))

  # --- Comparators -------------------------------------------------------
  net <- readRDS(P("rds", "network.rds"))
  conv <- net$consensus$gene
  cp_f <- P("rds", "citation_penalised.rds")
  best_prev <- NULL
  if (file.exists(cp_f)) {
    ev <- readRDS(cp_f)
    log_msg("Previous best correction: ",
            ev$method[which.min(ev$auc_citation_bias)])
  }

  lb <- readRDS(P("rds", "literature_bias.rds"))
  papers <- setNames(lb$data$n_papers, lb$data$gene)
  degg   <- setNames(lb$data$degree, lb$data$gene)

  sel <- read_csv(P("tables", "depmap_glioma_selectivity.csv"),
                  show_col_types = FALSE)

  sets <- list(`cytoHubba on STRING` = conv,
               `co-essentiality (all lines)` = hubs_all,
               `co-essentiality (CNS held out)` = hubs_ho)

  pool <- intersect(g_in, names(degg))
  ev <- imap_dfr(sets, function(h, nm) {
    s <- sel %>% mutate(in_set = gene %in% h)
    auc_sel <- if (sum(s$in_set) >= 3)
      as.numeric(pROC::auc(pROC::roc(s$in_set, s$selectivity,
                                     quiet = TRUE))) else NA_real_
    tibble(method = nm, n = length(h),
           auc_degree_bias = auc_of(h, pool, log10(degg + 1)),
           auc_citation_bias = auc_of(h, pool, log10(papers[pool] + 1)),
           median_papers = median(papers[h], na.rm = TRUE),
           auc_selectivity = auc_sel,
           pct_pan_essential = 100 * mean(s$pan_essential[s$in_set],
                                          na.rm = TRUE),
           overlap_cytoHubba = length(intersect(h, conv)))
  })

  write_csv(ev, P("tables", "coessentiality_comparison.csv"))

  log_msg("=================================================")
  log_msg("CO-ESSENTIALITY NETWORK VERSUS STRING")
  log_msg("  bias columns LOWER better, 0.5 ideal")
  log_msg("  selectivity HIGHER better")
  log_msg("  benchmarks: cytoHubba 0.948 citation bias / 0.541 selectivity")
  log_msg("              best correction so far 0.745 / 0.668")
  log_msg("=================================================")
  print(as.data.frame(ev %>% mutate(across(where(is.numeric),
                                           ~ round(.x, 3)))))

  # --- Verdict -----------------------------------------------------------
  ho <- ev %>% filter(grepl("held out", method))
  log_msg("---")
  log_msg("The CNS-held-out row is the one to believe, since its network ",
          "was built without the cell lines used to score selectivity.")

  if (nrow(ho) && !is.na(ho$auc_citation_bias) &&
      ho$auc_citation_bias < 0.65 && !is.na(ho$auc_selectivity) &&
      ho$auc_selectivity > 0.6) {
    log_msg("VERDICT: an experimentally derived network reduces attention ",
            "bias to ", round(ho$auc_citation_bias, 3),
            " while retaining glioma-selective dependency (",
            round(ho$auc_selectivity, 3),
            "). The substrate, not the scoring function, was the binding ",
            "constraint. This is a method.")
  } else if (nrow(ho) && !is.na(ho$auc_citation_bias) &&
             ho$auc_citation_bias < 0.65) {
    log_msg("VERDICT: bias is removed (", round(ho$auc_citation_bias, 3),
            ") but selectivity is ", round(ho$auc_selectivity, 3),
            ". The same outcome as the rewiring z-score and the citation ",
            "penalty: less biased, finds nothing. Eight attempts, none ",
            "successful. Write the diagnostics.")
  } else {
    log_msg("VERDICT: the co-essentiality network does not escape the bias ",
            "(citation AUC ", round(ho$auc_citation_bias, 3),
            "). Plausible reason: CRISPR libraries and cell line panels are ",
            "themselves shaped by prior interest. Eight attempts, none ",
            "successful. Write the diagnostics.")
  }

  saveRDS(list(evaluation = ev, hubs_all = hubs_all, hubs_ho = hubs_ho),
          P("rds", "coessentiality.rds"))
  log_msg("26_coessentiality_network complete.")
  invisible(ev)
}
