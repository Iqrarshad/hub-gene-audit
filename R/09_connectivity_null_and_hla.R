# Two analyses.
#
# A. Within-network null for degree enrichment. Comparing against random
#    gene backgrounds tests whether the gene set is functionally coherent,
#    which it is by construction. The comparison is therefore made against
#    the other genes in the same network.
#
# B. Direction of the antigen presentation genes in the grade contrast.
#    Enrichment alone does not give direction, and the interpretation
#    reverses depending on it: down in HGG indicates loss of antigen
#    presentation, up indicates inflammatory infiltration.

suppressPackageStartupMessages({
  library(igraph); library(dplyr); library(readr); library(tibble)
  library(httr); library(purrr); library(tidyr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

N_BG_ENRICH <- 25
API_DELAY   <- 0.5

# MHC class I, class II, and the processing/loading machinery
HLA_CLASS_I  <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M")
HLA_CLASS_II <- c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1",
                  "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CD74")
APM <- c("TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "CALR", "CANX",
         "PDIA3", "ERAP1", "ERAP2", "NLRC5", "CIITA")

# ---------------------------------------------------------------- Part A
string_net <- function(genes, score, tag) {
  cd <- file.path(CACHE_DIR, "string_enrich_null")
  dir.create(cd, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(cd, paste0(tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- tryCatch(httr::POST("https://string-db.org/api/tsv/network",
        body = list(identifiers = paste(genes, collapse = "%0d"),
                    species = 9606, required_score = score,
                    caller_identity = "glioma_enrich_null"),
        encode = "form"), error = function(e) NULL)
  e <- NULL
  if (!is.null(r) && httr::status_code(r) == 200) {
    txt <- httr::content(r, as = "text", encoding = "UTF-8")
    d <- tryCatch(readr::read_tsv(I(txt), show_col_types = FALSE),
                  error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0)
      e <- d %>% select(from = preferredName_A, to = preferredName_B) %>%
        filter(from != to) %>% distinct()
  }
  saveRDS(e, f); e
}

local_degree_of <- function(edges, gene) {
  if (is.null(edges) || nrow(edges) == 0) return(0L)
  sum(edges$from == gene) + sum(edges$to == gene)
}

# The first version of this test compared each gene's degree enrichment
# against random gene backgrounds. That asks whether the 148-gene SET is
# functionally coherent, which it is by construction, so 17 of 18 genes
# passed including TP53 at 1.36 and MYC at 1.87, and every p-value pinned to
# the 1/26 floor. The test answered a different question than intended.
#
# The right comparison is WITHIN the observed network: is this gene more
# specifically connected to the set than the other genes in the same set?
# The null distribution is the degree enrichment of all 148 genes, and the
# reference point is a gene's rank within it. No API calls, no permutation
# floor, and it directly addresses "is DLL3 unusual here".
enrichment_null <- function() {
  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  net  <- readRDS(P("rds", "network.rds"))

  # Genes with no STRING partners at all carry no information either way
  d <- norm %>% filter(!is.na(degree_enrichment), global_degree > 0)
  n <- nrow(d)
  log_msg("Within-network null across ", n, " genes with global degree data")
  if (n < 100) {
    log_msg("NOTE: the comparison field is only ", n, " genes. If the ",
            "normalised-degree table was built on a restricted gene field ",
            "rather than the full interactome, rebuild it (stage 08) first; ",
            "a field restricted to high-degree genes understates the ",
            "separation between conventional and normalised rankings.")
  }

  ecdf_e <- ecdf(d$degree_enrichment)

  res <- d %>%
    mutate(
      rank = rank(-degree_enrichment),
      percentile = 100 * ecdf_e(degree_enrichment),
      # One-sided empirical p: proportion of genes in the same network with
      # equal or greater enrichment.
      empirical_p = vapply(degree_enrichment,
                           function(x) mean(d$degree_enrichment >= x),
                           numeric(1)),
      fold_vs_median = degree_enrichment /
                       median(d$degree_enrichment, na.rm = TRUE),
      was_conventional_hub = gene %in% net$consensus$gene) %>%
    mutate(q = p.adjust(empirical_p, method = "BH"),
           top_decile = percentile >= 90) %>%
    arrange(desc(degree_enrichment))

  write_csv(res, P("tables", "degree_enrichment_within_network.csv"))

  log_msg("=== DEGREE ENRICHMENT, WITHIN-NETWORK NULL ===")
  log_msg("median enrichment across the network: ",
          round(median(d$degree_enrichment, na.rm = TRUE), 3))
  print(as.data.frame(res %>%
    select(gene, local_degree, global_degree, degree_enrichment,
           fold_vs_median, percentile, was_conventional_hub) %>%
    head(15)))

  log_msg("Conventional hubs, positioned within the same distribution:")
  print(as.data.frame(res %>% filter(was_conventional_hub) %>%
    select(gene, degree_enrichment, fold_vs_median, percentile, rank)))

  top <- res %>% filter(top_decile) %>% pull(gene)
  conv <- res %>% filter(was_conventional_hub) %>% pull(gene)
  log_msg("Top decile by specific connectivity: ",
          paste(head(top, 12), collapse = ", "))
  log_msg("Overlap with conventional hubs: ",
          length(intersect(top, conv)), " of ", length(conv))

  res
}

# ---------------------------------------------------------------- Part B
hla_direction <- function() {
  rna <- readRDS(P("rds", "deg_rnaseq.rds"))
  d <- rna$HGG_vs_LGG
  if (is.null(d)) { log_msg("No HGG_vs_LGG table."); return(NULL) }

  panels <- list(`MHC class I` = HLA_CLASS_I,
                 `MHC class II` = HLA_CLASS_II,
                 `Processing machinery` = APM)

  out <- imap_dfr(panels, function(genes, nm) {
    d %>% filter(gene %in% genes) %>%
      transmute(panel = nm, gene,
                logFC_raw, q_raw, logFC_adj, q_adj,
                composition_robust,
                direction = ifelse(logFC_adj > 0, "UP in HGG", "DOWN in HGG"))
  }) %>% arrange(panel, logFC_adj)

  write_csv(out, P("tables", "hla_direction_HGG_vs_LGG.csv"))

  log_msg("=== ANTIGEN PRESENTATION GENES, HGG vs LGG ===")
  print(as.data.frame(out))

  summ <- out %>% filter(!is.na(q_adj), q_adj < 0.05) %>%
    group_by(panel) %>%
    summarise(n_significant = n(),
              n_up = sum(logFC_adj > 0), n_down = sum(logFC_adj < 0),
              median_logFC = round(median(logFC_adj), 3), .groups = "drop")
  write_csv(summ, P("tables", "hla_direction_summary.csv"))
  log_msg("Summary (adjusted model, q < 0.05):")
  print(as.data.frame(summ))

  tot_up <- sum(summ$n_up); tot_dn <- sum(summ$n_down)
  if (tot_up + tot_dn == 0) {
    log_msg("No antigen presentation gene is significant after composition ",
            "adjustment. The KEGG enrichment may itself be composition ",
            "driven and must be reported with that caveat.")
  } else if (tot_dn > tot_up) {
    log_msg("INTERPRETATION: predominantly DOWN in HGG (", tot_dn, " down vs ",
            tot_up, " up). Consistent with loss of antigen presentation ",
            "during progression, i.e. immune escape.")
  } else {
    log_msg("INTERPRETATION: predominantly UP in HGG (", tot_up, " up vs ",
            tot_dn, " down). This is inflammatory infiltration, NOT immune ",
            "escape. The immune escape framing must be dropped.")
  }
  out
}

main_09 <- function() {
  log_msg("--- Part A: degree enrichment null ---")
  tryCatch(enrichment_null(),
           error = function(e) log_msg("Part A failed: ", conditionMessage(e)))
  log_msg("--- Part B: antigen presentation direction ---")
  tryCatch(hla_direction(),
           error = function(e) log_msg("Part B failed: ", conditionMessage(e)))
  log_msg("09_connectivity_null_and_hla complete.")
}
