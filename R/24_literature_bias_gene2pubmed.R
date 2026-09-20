# Decomposition of interactome degree into study intensity and residual.
#
#   log(degree) = a + b * log(papers) + residual
#
# NCBI gene2pubmed gives publication count per gene, independent of STRING.
# The fitted part is degree attributable to study intensity; the residual
# is degree not explained by recorded attention.
#
# The decisive test is which component predicts hub membership. The residual
# is not clean biology: it is degree unexplained by recorded attention, and
# curation gaps contribute to it.
#
# Requires gene2pubmed.gz from ftp.ncbi.nlm.nih.gov/gene/DATA/.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(data.table); library(pROC)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; distinct <- dplyr::distinct
summarise <- dplyr::summarise; arrange <- dplyr::arrange

# load_gene2pubmed() lives in the shared helper so the propagation-gates
# stage can build the same cache without depending on this stage having run
# first.
source("R/00f_citations.R")

main_24 <- function() {
  pap <- load_gene2pubmed()
  if (is.null(pap)) return(invisible(NULL))

  # Interactome-wide degree, produced by the propagation-gates stage as an
  # unweighted count of STRING partners at combined score >= 400. Two
  # variants: degree_all uses all channels including text-mining, which is the
  # score that builds the hub network (headline); degree_notm excludes
  # text-mining (robustness check).
  bg_f <- P("rds", "background_global_degrees.rds")
  if (!file.exists(bg_f))
    stop("background_global_degrees.rds not found. Run the propagation-gates ",
         "stage (27) first; it writes the interactome-wide degree table.")
  bg <- as_tibble(readRDS(bg_f))

  # Normalised degrees are retained only for the Q3 candidate filter below.
  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE) %>%
    select(gene, degree = global_degree)

  # --- Q1: how much of degree does citation explain? --------------------
  run_q1 <- function(degcol, label, csv) {
    dd <- tibble(gene = bg$gene, degree = bg[[degcol]]) %>%
      filter(!is.na(degree), degree > 0) %>%
      inner_join(pap, by = "gene") %>%
      mutate(ld = log10(degree + 1), lp = log10(n_papers + 1))
    if (nrow(dd) < 100) stop("Too few genes for the ", label, " regression.")
    fit <- lm(ld ~ lp, data = dd)
    dd$degree_resid <- residuals(fit)
    dd$degree_from_papers <- fitted(fit)
    r2 <- summary(fit)$r.squared
    rho <- cor(dd$ld, dd$lp, method = "spearman")
    log_msg(label, ": ", nrow(dd), " genes; Spearman ", round(rho, 3),
            "; R squared ", round(r2, 3), " -> ", round(100 * r2, 1),
            " percent of degree variance explained by publication count")
    write_csv(dd, P("tables", csv))
    list(d = dd, r2 = r2, rho = rho)
  }

  log_msg("=================================================")
  log_msg("Q1: DOES STUDY INTENSITY EXPLAIN INTERACTOME DEGREE?")
  log_msg("=================================================")
  all_q1  <- run_q1("degree_all",
                    "HEADLINE all channels (text-mining included)",
                    "degree_vs_papers.csv")
  notm_q1 <- run_q1("degree_notm",
                    "ROBUSTNESS text-mining excluded",
                    "degree_vs_papers_notextmining.csv")

  # The all-channels degree drives the downstream hub-component analysis, so
  # the degree that is shown to be citation-driven is the same degree that
  # selects the hubs.
  d   <- all_q1$d
  r2  <- all_q1$r2
  rho <- all_q1$rho
  log_msg("Genes with both degree and publication count: ", nrow(d))

  # --- Q2: what predicts hub membership? --------------------------------
  net <- readRDS(P("rds", "network.rds"))
  own_hubs <- net$consensus$gene
  lit_src <- if (file.exists("published_hub_lists.csv"))
    "published_hub_lists.csv" else P("tables", "published_hub_lists.csv")
  if (!file.exists(lit_src))
    stop("published_hub_lists.csv not found in the project root or tables. ",
         "It ships with the repository.")
  lit <- read_csv(lit_src, show_col_types = FALSE)
  lit_hubs <- lit %>% mutate(gene = strsplit(hub_genes, ";")) %>%
    unnest(gene) %>% mutate(gene = toupper(trimws(gene))) %>%
    filter(grepl("^[A-Z][A-Z0-9._-]{0,20}$", gene)) %>%
    distinct(gene) %>% pull(gene)

  auc_of <- function(y, x) {
    ok <- !is.na(x) & !is.na(y)
    if (sum(y[ok]) < 3) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(y[ok], x[ok], quiet = TRUE)))
  }

  q2 <- map_dfr(list(`our cytoHubba hubs` = own_hubs,
                     `published hub genes` = lit_hubs), function(h) {
    dd <- d %>% mutate(y = as.integer(gene %in% h))
    if (sum(dd$y) < 3) return(tibble())
    tibble(n_hubs = sum(dd$y),
           auc_degree = auc_of(dd$y, dd$ld),
           auc_papers = auc_of(dd$y, dd$lp),
           auc_degree_residual = auc_of(dd$y, dd$degree_resid),
           auc_degree_from_papers = auc_of(dd$y, dd$degree_from_papers))
  }, .id = "hub_set")

  write_csv(q2, P("tables", "hub_prediction_components.csv"))

  log_msg("=================================================")
  log_msg("Q2: WHICH COMPONENT PREDICTS HUB MEMBERSHIP?")
  log_msg("  auc_papers            study intensity alone")
  log_msg("  auc_degree_from_papers the citation-explained part of degree")
  log_msg("  auc_degree_residual   degree NOT explained by citations")
  log_msg("=================================================")
  print(as.data.frame(q2 %>% mutate(across(where(is.numeric),
                                           ~ round(.x, 3)))))

  # --- Q3: would a citation-penalised score differ? ---------------------
  # Rank genes by degree residual instead of raw degree, and compare the
  # nominations. The degree-penalised scores in stages 19 to 35 nominated
  # low-degree genes; a citation-penalised score should nominate genes that
  # are well connected relative to how little they have been studied.
  cand <- d %>% filter(gene %in% norm$gene)
  if (nrow(cand) > 20) {
    top_resid <- cand %>% arrange(desc(degree_resid)) %>% head(10)
    top_deg   <- cand %>% arrange(desc(degree)) %>% head(10)
    log_msg("Top 10 by raw degree:      ",
            paste(top_deg$gene, collapse = ", "))
    log_msg("Top 10 by degree residual: ",
            paste(top_resid$gene, collapse = ", "))
    log_msg("  residual top 10 median papers ",
            round(median(top_resid$n_papers)),
            " vs raw-degree top 10 median papers ",
            round(median(top_deg$n_papers)))
    write_csv(cand %>% arrange(desc(degree_resid)),
              P("tables", "citation_penalised_ranking.csv"))
  }

  # --- Verdict -----------------------------------------------------------
  lit_row <- q2 %>% filter(hub_set == "published hub genes")
  log_msg("---")
  if (nrow(lit_row)) {
    ap <- lit_row$auc_papers; ar <- lit_row$auc_degree_residual
    ad <- lit_row$auc_degree
    log_msg("For published hub genes: degree AUC ", round(ad, 3),
            ", citations AUC ", round(ap, 3),
            ", degree-residual AUC ", round(ar, 3))
    if (!is.na(ap) && !is.na(ar) && ap > ar + 0.1) {
      log_msg("VERDICT: publication count predicts hub membership better ",
              "than the part of degree unexplained by it. The bias claim is ",
              "measured rather than assumed, and a citation-penalised score ",
              "is worth building.")
    } else if (!is.na(ar) && !is.na(ap) && ar > ap + 0.1) {
      log_msg("VERDICT: the residual predicts hub membership better than ",
              "citations. Hub genes are genuinely well connected beyond ",
              "what study intensity explains. The literature-bias framing ",
              "is NOT supported and must be dropped; the claim stays as ",
              "'degree-driven', which the data do support.")
    } else {
      log_msg("VERDICT: citations and the residual predict hub membership ",
              "comparably (", round(ap, 3), " vs ", round(ar, 3),
              "). Study intensity and genuine connectivity cannot be ",
              "separated in this data. Report both and claim neither ",
              "exclusively.")
    }
  }

  saveRDS(list(data = d, q2 = q2, r2 = r2, rho = rho,
               r2_notm = notm_q1$r2, rho_notm = notm_q1$rho,
               data_notm = notm_q1$d),
          P("rds", "literature_bias.rds"))
  log_msg("24_literature_bias_gene2pubmed complete.")
  invisible(q2)
}
