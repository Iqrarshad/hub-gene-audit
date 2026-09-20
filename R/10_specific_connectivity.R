# Binomial test for specific connectivity.
#
# A gene has K partners across the interactome; the gene set covers a
# fraction p of annotated proteins. Under the null, k ~ Binomial(K, p),
# and the one-sided p-value is P(X >= k), BH-corrected.
#
# This replaces an observed-over-expected ratio, which rewards low
# denominators: a single edge on a gene with nine known partners gives a
# ratio above 14 without constituting evidence.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr); library(tidyr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# Number of protein-coding genes with STRING annotation. Used as the
# denominator for the background rate. STRING v11/v12 covers roughly 19,500
# human proteins.
STRING_PROTEOME <- 19500

main_10 <- function() {
  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  inter <- readRDS(P("rds", "intersect.rds"))
  net   <- readRDS(P("rds", "network.rds"))

  n_set <- length(inter$shared$gene)
  p_bg  <- n_set / STRING_PROTEOME
  log_msg("Gene set: ", n_set, " genes; background rate p = ",
          signif(p_bg, 4))

  d <- norm %>%
    filter(!is.na(global_degree), global_degree > 0) %>%
    mutate(
      k = local_degree,
      K = global_degree,
      expected = K * p_bg,
      obs_exp = k / pmax(expected, 1e-9),
      # P(X >= k) under Binomial(K, p_bg)
      p_binom = pbinom(k - 1, size = K, prob = p_bg, lower.tail = FALSE))

  d <- d %>%
    mutate(q_binom = p.adjust(p_binom, method = "BH"),
           significant = q_binom < 0.05,
           was_conventional_hub = gene %in% net$consensus$gene) %>%
    arrange(p_binom, desc(obs_exp))

  write_csv(d %>% select(gene, k, K, expected, obs_exp, p_binom, q_binom,
                         significant, was_conventional_hub,
                         n_methods_top10),
            P("tables", "specific_connectivity_binomial.csv"))

  sig <- d %>% filter(significant)

  log_msg("=================================================")
  log_msg("SPECIFIC CONNECTIVITY, BINOMIAL TEST")
  log_msg("  k = partners inside the gene set")
  log_msg("  K = partners across the whole interactome")
  log_msg("=================================================")
  log_msg(nrow(sig), " of ", nrow(d), " genes significant at BH q < 0.05")

  print(as.data.frame(sig %>%
    select(gene, k, K, expected, obs_exp, p_binom, q_binom,
           was_conventional_hub) %>%
    mutate(across(c(expected, obs_exp), ~ round(.x, 2)),
           across(c(p_binom, q_binom), ~ signif(.x, 3))) %>%
    head(30)))

  # --- What the ratio did that the test does not -------------------------
  ratio_top <- d %>% arrange(desc(obs_exp)) %>% head(15) %>% pull(gene)
  test_top  <- sig %>% head(15) %>% pull(gene)
  dropped <- setdiff(ratio_top, test_top)

  if (length(dropped)) {
    log_msg("Genes in the ratio's top 15 that the binomial test rejects:")
    print(as.data.frame(d %>% filter(gene %in% dropped) %>%
      select(gene, k, K, obs_exp, q_binom) %>%
      mutate(obs_exp = round(obs_exp, 1), q_binom = signif(q_binom, 3))))
    log_msg("These entered the ratio ranking on 1-2 edges. Low denominators ",
            "inflate a ratio but cannot reach significance.")
  }

  # --- Conventional hubs, positioned ------------------------------------
  conv <- d %>% filter(was_conventional_hub) %>%
    select(gene, k, K, expected, obs_exp, p_binom, q_binom, significant) %>%
    mutate(across(c(expected, obs_exp), ~ round(.x, 2)),
           across(c(p_binom, q_binom), ~ signif(.x, 3)))
  log_msg("Conventional hubs under the binomial test:")
  print(as.data.frame(conv))

  n_conv_sig <- sum(conv$significant, na.rm = TRUE)
  log_msg(n_conv_sig, " of ", nrow(conv), " conventional hubs are ",
          "significantly specifically connected")

  # --- Agreement between the two approaches ------------------------------
  overlap <- length(intersect(sig$gene, net$consensus$gene))
  log_msg("Overlap between binomial-significant genes and conventional ",
          "hubs: ", overlap, " of ", length(net$consensus$gene))

  # Rank correlation between conventional consensus and specific connectivity
  rc <- suppressWarnings(cor(d$n_methods_top10, -log10(d$p_binom + 1e-300),
                             method = "spearman", use = "complete.obs"))
  log_msg("Spearman correlation between conventional consensus count and ",
          "specific connectivity evidence: ", round(rc, 3))
  if (!is.na(rc) && abs(rc) < 0.2) {
    log_msg("The two rankings are close to independent: conventional ",
            "centrality and specific connectivity measure different things.")
  }

  saveRDS(d, P("rds", "specific_connectivity.rds"))
  log_msg("10_specific_connectivity complete.")
  invisible(d)
}
