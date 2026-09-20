# Gate checks for attention-penalised network propagation, at
# transcriptome scale.
#
# Kill criteria, fixed before the data were examined:
#   G1  mean weighted degree >= 2.0      graph not shattered
#   G2  R2(E, C) < 0.01                  input not already citation-driven
#   G3a rho(x*,C) - rho(k,C) > 0.05      propagation amplifies the bias
#                                        rather than restating the
#                                        degree-citation correlation
#   G3b rho(x*,E) >= 0.20                signal retained from the input
#   G4  rho(Tau,C) >= -0.30              tissue specificity validation is
#                                        not circular
#
# Text-mining and co-occurrence channels are excluded: both score
# co-mention in the literature and are proxies for study intensity.
#
# Requires the STRING full links and protein.info files, gene2pubmed, and
# the Human Protein Atlas table.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(dplyr); library(readr)
  library(tibble)
})

# load_gene2pubmed() builds and caches papers_per_gene.rds. Sourced here so
# this stage can build it directly when it runs before the degree-decomposition stage (24).
source("R/00f_citations.R")

MIN_DEGREE   <- 2.0
MAX_R2_EC    <- 0.01
MIN_RHO_GAIN <- 0.05
MIN_RHO_XE   <- 0.20
MIN_RHO_TAUC <- -0.30
STRING_CUT   <- 400
RESTART      <- 0.7

# ---- 1. STRING network, text-mining excluded ---------------------------
load_string <- function() {
  f <- find_local("9606.protein.links.full")
  if (is.na(f)) stop("Download 9606.protein.links.full.v12.0.txt.gz")
  log_msg("Reading STRING links (large) ...")
  d <- fread(f, data.table = TRUE, showProgress = FALSE)

  ch <- intersect(names(d), c("experiments", "experiments_transferred",
                              "database", "database_transferred",
                              "coexpression", "coexpression_transferred",
                              "neighborhood", "neighborhood_transferred",
                              "fusion"))
  log_msg("Channels used (text-mining excluded score): ",
          paste(ch, collapse = ", "))
  log_msg("Channels EXCLUDED: textmining, textmining_transferred, cooccurence")
  stopifnot(length(ch) >= 3)

  # Text-mining-excluded score: STRING combines channels as 1 - prod(1 - p_i)
  d[, comb := 1000 * (1 - Reduce(`*`, lapply(.SD, function(v) 1 - v / 1000))),
    .SDcols = ch]

  # All-channels score. STRING's own combined_score column includes the
  # text-mining channel and matches the combined score the STRING API returns
  # to the network stage (06), so degree computed from it is consistent with the network
  # that selects the hubs. If the column is absent, reconstruct it over every
  # numeric channel present, text-mining included.
  if (!"combined_score" %in% names(d)) {
    ch_all <- setdiff(names(d), c("protein1", "protein2", "comb"))
    ch_all <- ch_all[vapply(d[, ..ch_all], is.numeric, logical(1))]
    d[, combined_score := 1000 *
        (1 - Reduce(`*`, lapply(.SD, function(v) 1 - v / 1000))),
      .SDcols = ch_all]
    log_msg("combined_score absent; reconstructed all-channels score over: ",
            paste(ch_all, collapse = ", "))
  }

  # Keep only edges passing either threshold, so the symbol map is applied to
  # the union of the two networks rather than the full 13M-row file.
  keep <- d[comb >= STRING_CUT | combined_score >= STRING_CUT,
            .(protein1, protein2, comb, combined_score)]
  log_msg("Edges passing either threshold (score >= ", STRING_CUT, "): ",
          nrow(keep))

  # --- protein id to gene symbol (map built once, used for both outputs) --
  # The aliases file carries many alias TYPES per protein: transcript names
  # (ARF5-201), free-text descriptions, numeric IDs. Taking the first match
  # per protein returns junk for most entries, which collapsed the usable
  # protein.info gives exactly one preferred name per
  inf <- find_local("9606.protein.info")
  if (!is.na(inf)) {
    pi_ <- fread(inf, data.table = TRUE)
    pcol <- grep("protein_external_id|string_protein_id", names(pi_),
                 value = TRUE)[1]
    if (is.na(pcol)) pcol <- names(pi_)[1]
    ncol_ <- grep("preferred_name", names(pi_), value = TRUE)[1]
    if (is.na(ncol_)) stop("No preferred_name column in protein.info")
    map <- setNames(as.character(pi_[[ncol_]]), as.character(pi_[[pcol]]))
    log_msg("Mapping from protein.info: ", length(map), " proteins")
  } else {
    al <- find_local("9606.protein.aliases")
    if (is.na(al))
      stop("Download 9606.protein.info.v12.0.txt.gz (preferred), ",
           "or 9606.protein.aliases.v12.0.txt.gz")
    a <- fread(al, data.table = TRUE)
    setnames(a, 1:3, c("protein", "alias", "source"))
    a <- a[source %in% c("Ensembl_HGNC_symbol", "BioMart_HUGO",
                         "Ensembl_gene_name")]
    a <- a[grepl("^[A-Za-z][A-Za-z0-9._-]{0,20}$", alias)]
    a <- unique(a, by = "protein")
    map <- setNames(a$alias, a$protein)
    log_msg("Mapping from aliases fallback: ", length(map), " proteins")
  }

  keep[, g1 := map[protein1]][, g2 := map[protein2]]
  keep <- keep[!is.na(g1) & !is.na(g2) & g1 != g2]
  keep[, `:=`(a = pmin(g1, g2), b = pmax(g1, g2))]

  # --- interactome-wide unweighted degree (distinct partners per gene) ----
  # Two variants: all channels (headline, matches the hub network) and
  # text-mining excluded (robustness). Written for the degree-decomposition stage (24) regression.
  deg_count <- function(edges) {
    e <- unique(edges[, .(a, b)])
    rbindlist(list(e[, .(gene = a)], e[, .(gene = b)]))[, .(degree = .N),
                                                        by = gene]
  }
  d_all  <- deg_count(keep[combined_score >= STRING_CUT])
  d_notm <- deg_count(keep[comb          >= STRING_CUT])
  bg <- merge(setnames(copy(d_all),  "degree", "degree_all"),
              setnames(copy(d_notm), "degree", "degree_notm"),
              by = "gene", all = TRUE)
  bg[is.na(degree_all),  degree_all  := 0L]
  bg[is.na(degree_notm), degree_notm := 0L]
  saveRDS(as_tibble(bg), P("rds", "background_global_degrees.rds"))
  log_msg("Global degrees written: ", nrow(bg),
          " genes (all-channels and text-mining-excluded, unweighted, ",
          "score >= ", STRING_CUT, ")")

  # --- propagation edge table (text-mining-excluded, weighted) unchanged --
  ed <- unique(keep[comb >= STRING_CUT, .(a, b, comb)])
  log_msg("Symbol-mapped unique edges (propagation network): ", nrow(ed))

  # Fail loudly rather than silently proceeding on a broken mapping
  ref <- c("TP53", "MYC", "EGFR", "CDK1", "ACTB", "GAPDH", "SOX2", "PTEN")
  found <- sum(ref %in% unique(c(ed$a, ed$b)))
  log_msg("Symbol sanity: ", found, "/", length(ref),
          " reference genes present")
  if (found < 6)
    stop("Protein-to-symbol mapping looks wrong. Sample symbols: ",
         paste(head(unique(ed$a), 8), collapse = ", "))
  ed
}

# ---- 2. Citations -------------------------------------------------------
load_citations <- function() {
  pap <- load_gene2pubmed()
  if (is.null(pap))
    stop("gene2pubmed.gz required in DATA_DIR to build papers_per_gene.rds. ",
         "Download from ftp.ncbi.nlm.nih.gov/gene/DATA/ and place it in ",
         "DATA_DIR, then rerun.")
  pap
}

# ---- 3. Transcriptome-wide differential expression ---------------------
load_full_de <- function() {
  rna <- readRDS(P("rds", "deg_rnaseq.rds"))
  d <- rna$LGG_vs_Normal
  if (is.null(d)) stop("No LGG_vs_Normal table in deg_rnaseq.rds")
  log_msg("DE table rows: ", nrow(d))
  if (nrow(d) < 8000)
    stop("DE table has only ", nrow(d), " genes. Gate 2 requires the FULL ",
         "tested transcriptome, not a significance-filtered set.")
  col <- if ("logFC_adj" %in% names(d) && !all(is.na(d$logFC_adj)))
    "logFC_adj" else "logFC_raw"
  log_msg("Using effect column: ", col)
  tibble(gene = d$gene, E_raw = abs(d[[col]])) %>%
    filter(is.finite(E_raw)) %>%
    group_by(gene) %>% summarise(E_raw = max(E_raw), .groups = "drop")
}

# ---- 4. Tissue specificity ---------------------------------------------
load_tau <- function() {
  f <- find_local("proteinatlas")
  if (is.na(f)) { log_msg("HPA file absent; Gate 4 skipped."); return(NULL) }
  h <- fread(f, data.table = FALSE)
  gcol <- grep("^Gene$", names(h), value = TRUE)[1]
  scol <- grep("RNA tissue specificity|specificity|Tau", names(h),
               ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gcol) || is.na(scol)) {
    log_msg("HPA columns not identified; Gate 4 skipped."); return(NULL)
  }
  log_msg("HPA specificity column: ", scol)
  v <- h[[scol]]
  if (!is.numeric(v)) {
    lv <- c("Low tissue specificity" = 0, "Tissue enhanced" = 0.5,
            "Group enriched" = 0.75, "Tissue enriched" = 1)
    v <- unname(lv[as.character(v)])
    log_msg("Specificity is categorical; mapped to an ordinal scale. ",
            "This is coarser than a continuous Tau and should be replaced ",
            "with GTEx-derived Tau before publication.")
  }
  tibble(gene = as.character(h[[gcol]]), tau = v) %>% filter(!is.na(tau))
}

# ---- Random walk with restart ------------------------------------------
rwr <- function(M, E, r = RESTART, tol = 1e-6, maxit = 500) {
  x <- E
  for (i in seq_len(maxit)) {
    xn <- (1 - r) * as.vector(Matrix::crossprod(M, x)) + r * E
    if (sum(abs(xn - x)) < tol) {
      log_msg("RWR converged in ", i, " iterations"); return(xn)
    }
    x <- xn
  }
  log_msg("RWR hit the iteration cap without converging"); x
}

main_27 <- function() {
  ed  <- load_string()
  cit <- load_citations()
  de  <- load_full_de()
  tau <- load_tau()

  net_genes <- unique(c(ed$a, ed$b))
  log_msg("Genes in network: ", length(net_genes),
          "; with citations: ", nrow(cit), "; with DE: ", nrow(de))
  log_msg("Pairwise overlaps -- network/DE: ",
          length(intersect(net_genes, de$gene)),
          "; network/citations: ", length(intersect(net_genes, cit$gene)),
          "; DE/citations: ", length(intersect(de$gene, cit$gene)))

  genes <- sort(Reduce(intersect, list(net_genes, cit$gene, de$gene)))
  log_msg("Genes with all three: ", length(genes))
  if (length(genes) < 5000)
    stop("Only ", length(genes), " genes have network, citation and DE data. ",
         "Check the identifier formats before proceeding.")

  idx <- setNames(seq_along(genes), genes)
  e2  <- ed[a %in% genes & b %in% genes]
  A <- sparseMatrix(i = idx[e2$a], j = idx[e2$b], x = e2$comb / 1000,
                    dims = c(length(genes), length(genes)), symmetric = TRUE)
  A <- as(A, "generalMatrix")

  k <- Matrix::rowSums(A)
  kbar <- mean(k)

  E <- de$E_raw[match(genes, de$gene)]; E[is.na(E)] <- 0
  if (sum(E) == 0) stop("Effect vector is all zero.")
  E <- E / sum(E)
  C <- log2(1 + cit$n_papers[match(genes, cit$gene)]); C[is.na(C)] <- 0

  # ---- GATE 1: sparsity -------------------------------------------------
  log_msg("=================== GATE 1: sparsity ===================")
  log_msg("mean weighted degree = ", round(kbar, 3))
  log_msg("isolated nodes = ", sum(k == 0), " (",
          round(100 * mean(k == 0), 1), "%)")
  g1 <- kbar >= MIN_DEGREE
  log_msg("GATE 1 ", ifelse(g1, "PASS", "FAIL: graph shattered"))

  # ---- GATE 2: input cleanliness ---------------------------------------
  log_msg("============= GATE 2: input cleanliness ===============")
  fit <- lm(E ~ C)
  r2 <- summary(fit)$r.squared
  log_msg("R2(E, C) = ", signif(r2, 5))
  log_msg("Spearman(E, C) = ", round(cor(E, C, method = "spearman"), 4))
  g2 <- r2 < MAX_R2_EC
  log_msg("GATE 2 ", ifelse(g2, "PASS",
          "FAIL: expression signal is itself citation-driven"))

  # ---- GATE 3: corruption baseline -------------------------------------
  log_msg("============ GATE 3: corruption baseline ==============")
  rs <- Matrix::rowSums(A); rs[rs == 0] <- 1
  M <- A / rs
  x0 <- rwr(M, E)

  rho_xC <- cor(x0, C, method = "spearman")
  rho_kC <- cor(k,  C, method = "spearman")
  rho_xE <- cor(x0, E, method = "spearman")
  gain <- rho_xC - rho_kC

  log_msg("rho(x*, C) = ", round(rho_xC, 4))
  log_msg("rho(k,  C) = ", round(rho_kC, 4), "   <- the null to beat")
  log_msg("gain       = ", round(gain, 4),
          "   (propagation must AMPLIFY the bias, not merely restate it)")
  log_msg("rho(x*, E) = ", round(rho_xE, 4),
          "   (signal retained from the restart vector)")
  g3a <- gain > MIN_RHO_GAIN
  g3b <- rho_xE >= MIN_RHO_XE
  log_msg("GATE 3a amplification  ", ifelse(g3a, "PASS",
          "FAIL: propagation only restates the degree-citation correlation"))
  log_msg("GATE 3b signal retained ", ifelse(g3b, "PASS",
          "FAIL: the network erases the transcriptomic signal"))

  # ---- GATE 4: validation independence ---------------------------------
  g4 <- NA; rho_tc <- NA_real_
  if (!is.null(tau)) {
    log_msg("========= GATE 4: validation independence =============")
    tv <- tau$tau[match(genes, tau$gene)]
    ok <- !is.na(tv)
    rho_tc <- cor(tv[ok], C[ok], method = "spearman")
    log_msg("rho(Tau, C) = ", round(rho_tc, 4), " on ", sum(ok), " genes")
    g4 <- rho_tc >= MIN_RHO_TAUC
    log_msg("GATE 4 ", ifelse(g4, "PASS",
            "FAIL: Tau validation requires a degree-and-citation-matched null"))
  }

  res <- tibble(
    gate = c("1 sparsity", "2 input clean", "3a amplification",
             "3b signal retained", "4 Tau independence"),
    value = c(kbar, r2, gain, rho_xE, rho_tc),
    threshold = c(MIN_DEGREE, MAX_R2_EC, MIN_RHO_GAIN, MIN_RHO_XE,
                  MIN_RHO_TAUC),
    direction = c(">=", "<", ">", ">=", ">="),
    pass = c(g1, g2, g3a, g3b, g4))
  write_csv(res, P("tables", "propagation_gates.csv"))

  log_msg("======================= SUMMARY =======================")
  print(as.data.frame(res %>% mutate(across(where(is.numeric),
                                            ~ signif(.x, 4)))))

  log_msg("---")
  if (all(res$pass[1:4], na.rm = TRUE)) {
    log_msg("GATES 1-3 CLEARED. Propagation phase is warranted.")
    if (!isTRUE(g4))
      log_msg("Gate 4 failed: Tau and citations are anti-correlated, so the ",
              "Tau validation must be run against a degree-and-citation-",
              "matched null rather than against the unpenalised set.")
  } else {
    log_msg("KILL CRITERION MET. Abandon attention-penalised propagation ",
            "as specified.")
    log_msg("Failing gates: ",
            paste(res$gate[!res$pass & !is.na(res$pass)], collapse = ", "))
  }

  saveRDS(list(genes = genes, A = A, E = E, C = C, k = k, x0 = x0,
               tau = tau, gates = res), P("rds", "propagation_gates.rds"))
  invisible(res)
}