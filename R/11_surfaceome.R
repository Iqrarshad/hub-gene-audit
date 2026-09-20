# Surfaceome composition of the hub nominations.
#
# Conventional centrality nominates intracellular and secreted proteins,
# which are not addressable by antibody or CAR. This tests whether the
# corrected ranking is enriched for surface proteins relative to the
# conventional one.
#
# Reference is SURFY (Bausch-Fluck 2018 PNAS). A table with fewer than 500
# genes is rejected rather than used.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# Minimal fallback: well-established surface proteins relevant to glioma and
# immuno-oncology. Used only if the SURFY table is unavailable.
FALLBACK_SURFACE <- c(
  "DLL3","PTPRZ1","EGFR","CD44","IL13RA2","CD276","CSPG4","GPNMB","MET",
  "PDGFRA","ERBB2","NCAM1","L1CAM","CD70","EPHA2","EPHA3","ITGA6","ITGB1",
  "CD9","CD81","CD151","TFRC","LRP1","NOTCH1","NOTCH2","JAG1","DLL1","DLL4",
  "CDH2","NRCAM","CNTN1","THY1","PROM1","ATP1B2","SLC1A3","AQP4","GJA1",
  "MOG","MAG","PLP1","CD40","CD47","CD24","TNC","VCAN","SDC1","GPC1","GPC3",
  "ANTXR1","MSLN","FOLR1","TROP2","TACSTD2","MUC1","EPCAM","B2M","HLA-A",
  "HLA-B","HLA-C","HLA-DRA","CD74","PTPRC","ITGAM","CSF1R","P2RY12","TMEM119",
  "CX3CR1","TREM2","SIRPA","LILRB4","VSIR","HAVCR2","LAG3","TIGIT","PDCD1",
  "CD274","PDCD1LG2","CTLA4","TNFRSF9","ICOS","CD28","IL7R","IL2RA","KIT",
  "FLT1","KDR","TEK","PECAM1","CDH5","ENG","MCAM","PDGFRB","ACTA2","RGS5")

# Reads the SURFY master table (Bausch-Fluck 2018 PNAS, table S3). The file
# has two sheets; "in silico surfaceome only" holds the 2,799 predicted
# surface proteins. The header sits on row 2, since row 1 is a title.
#
# A reference with fewer than 500 genes is REJECTED rather than used. An
# A reference with fewer than 500 genes is rejected rather than used.
# failing.
load_surfaceome <- function() {
  cand <- list.files(DATA_DIR, pattern = "surfaceome|surfy|table_S3",
                     recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  cand <- c(cand, list.files(dirname(DATA_DIR), pattern = "surfaceome|surfy",
                             recursive = TRUE, full.names = TRUE,
                             ignore.case = TRUE))
  cand <- unique(cand[grepl("xlsx|xls|csv|tsv|txt", cand, ignore.case = TRUE)])

  for (f in cand) {
    g <- try_read_surfaceome(f)
    if (!is.null(g) && length(g) >= 500) {
      log_msg("Surfaceome reference: ", basename(f), " (", length(g),
              " surface proteins)")
      return(list(genes = g, source = basename(f), n_total = length(g),
                  fallback = FALSE))
    }
    if (!is.null(g)) {
      log_msg("Rejected ", basename(f), ": only ", length(g),
              " genes, implausible for a surfaceome reference")
    }
  }

  log_msg("No valid surfaceome table found. Falling back to a curated list ",
          "of ", length(FALLBACK_SURFACE), " established surface proteins. ",
          "This is weaker evidence and must be stated in the write-up.")
  list(genes = FALLBACK_SURFACE, source = "curated fallback",
       n_total = length(FALLBACK_SURFACE), fallback = TRUE)
}

try_read_surfaceome <- function(f) {
  ext <- tolower(tools::file_ext(f))
  d <- NULL

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      log_msg("readxl not installed; cannot read ", basename(f),
              '. Install with install.packages("readxl")')
      return(NULL)
    }
    sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) NULL)
    if (is.null(sheets)) return(NULL)
    # Prefer the in silico surfaceome sheet, else the master table
    pick <- grep("in silico|surfaceome only", sheets, ignore.case = TRUE,
                 value = TRUE)
    if (!length(pick)) pick <- sheets[1]
    # Header is on row 2 of this file, so skip the title row
    for (skip in c(1, 0)) {
      d <- tryCatch(readxl::read_excel(f, sheet = pick[1], skip = skip),
                    error = function(e) NULL)
      if (!is.null(d) && any(grepl("gene", names(d), ignore.case = TRUE))) break
    }
  } else {
    d <- tryCatch(data.table::fread(f, data.table = FALSE),
                  error = function(e) NULL)
  }
  if (is.null(d) || !nrow(d)) return(NULL)

  gcol <- grep("uniprot gene|^gene$|gene.?name|gene.?symbol|^symbol$",
               names(d), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gcol)) return(NULL)

  g <- toupper(trimws(as.character(d[[gcol]])))

  # If a label column exists, keep only entries marked as surface
  lcol <- grep("surfaceome label|^label$", names(d), ignore.case = TRUE,
               value = TRUE)[1]
  if (!is.na(lcol)) {
    keep <- grepl("^surface$", trimws(as.character(d[[lcol]])),
                  ignore.case = TRUE)
    if (sum(keep) >= 500) g <- g[keep]
  }

  g <- unique(g[g != "" & !is.na(g) & g != "NA"])
  g
}

main_11 <- function() {
  surf <- load_surfaceome()
  log_msg("Surfaceome reference: ", surf$source, " (", surf$n_total, " genes)")

  inter <- readRDS(P("rds", "intersect.rds"))
  genes148 <- inter$shared$gene

  net  <- readRDS(P("rds", "network.rds"))
  conv_hubs <- net$consensus$gene

  sc_f <- P("tables", "specific_connectivity_binomial.csv")
  if (!file.exists(sc_f)) stop("Run the specific-connectivity stage (10) first.")
  sc <- read_csv(sc_f, show_col_types = FALSE)
  corrected_hits <- sc %>% filter(significant) %>% pull(gene)

  is_surf <- function(g) toupper(g) %in% surf$genes

  # --- Composition of each nomination set --------------------------------
  sets <- list(
    `all 148 genes`         = genes148,
    `conventional hubs`     = conv_hubs,
    `corrected (binomial)`  = corrected_hits)

  tab <- imap_dfr(sets, function(g, nm) {
    tibble(set = nm, n = length(g),
           n_surface = sum(is_surf(g)),
           pct_surface = round(100 * sum(is_surf(g)) / max(length(g), 1), 1),
           surface_genes = paste(sort(g[is_surf(g)]), collapse = "; "))
  })
  write_csv(tab, P("tables", "surfaceome_composition.csv"))
  log_msg("Surfaceome composition of each nomination set:")
  print(as.data.frame(tab %>% select(-surface_genes)))
  for (i in seq_len(nrow(tab))) {
    if (tab$n_surface[i] > 0)
      log_msg("  ", tab$set[i], ": ", tab$surface_genes[i])
  }

  # --- Fisher exact, each set against the 148 background -----------------
  bg_surf <- sum(is_surf(genes148))
  bg_tot  <- length(genes148)

  ft <- imap_dfr(sets[-1], function(g, nm) {
    a <- sum(is_surf(g)); b <- length(g) - a
    c_ <- bg_surf - a;    d_ <- (bg_tot - length(g)) - c_
    m <- matrix(c(a, b, max(c_, 0), max(d_, 0)), nrow = 2)
    f <- fisher.test(m)
    tibble(set = nm, n_surface = a, n_total = length(g),
           odds_ratio = unname(f$estimate), p = f$p.value)
  })
  write_csv(ft, P("tables", "surfaceome_fisher.csv"))
  log_msg("Background: ", bg_surf, " of ", bg_tot, " genes in the ",
          "composition-robust set are surface proteins (",
          round(100 * bg_surf / bg_tot, 1), " percent)")
  log_msg("Enrichment against that background:")
  print(as.data.frame(ft %>% mutate(across(where(is.numeric),
                                           ~ signif(.x, 3)))))
  log_msg("NOTE: the nomination sets contain 4 and 8 genes. Fisher exact on ",
          "counts this small has very wide intervals. Treat the descriptive ",
          "contrast as the result and the p-value as illustrative only.")

  # --- Rank-based test across all 148 ------------------------------------
  # Do surfaceome genes rank higher by specific connectivity than by
  # conventional consensus?
  d <- sc %>% mutate(surface = is_surf(gene))
  if (sum(d$surface) >= 3) {
    w_spec <- wilcox.test(-log10(p_binom + 1e-300) ~ surface, data = d)
    w_conv <- wilcox.test(n_methods_top10 ~ surface, data = d)
    res <- tibble(
      ranking = c("specific connectivity (binomial)",
                  "conventional consensus count"),
      median_surface = c(
        median(-log10(d$p_binom[d$surface] + 1e-300)),
        median(d$n_methods_top10[d$surface])),
      median_nonsurface = c(
        median(-log10(d$p_binom[!d$surface] + 1e-300)),
        median(d$n_methods_top10[!d$surface])),
      p = c(w_spec$p.value, w_conv$p.value))
    write_csv(res, P("tables", "surfaceome_rank_test.csv"))
    log_msg("Do surface proteins rank higher under each approach?")
    print(as.data.frame(res %>% mutate(across(where(is.numeric),
                                              ~ signif(.x, 3)))))
  } else {
    log_msg("Fewer than 3 surface genes among the tested set; rank test ",
            "not informative.")
  }

  # --- Broader ranking comparison ----------------------------------------
  # Only 4 and 8 genes reach significance, which is too few to test. The
  # top 20 by each ranking gives a slightly larger comparison, still small
  # but more informative than the significant sets alone.
  top_n <- 20
  top_spec <- sc %>% arrange(p_binom) %>% head(top_n) %>% pull(gene)
  hub_tab_f <- P("tables", paste0("hub_ranking_", net$chosen, ".csv"))
  if (file.exists(hub_tab_f)) {
    ht <- read_csv(hub_tab_f, show_col_types = FALSE)
    top_conv <- ht %>% arrange(desc(n_methods_top10), desc(degree)) %>%
      head(top_n) %>% pull(gene)

    cmp2 <- tibble(
      ranking = c("conventional centrality", "specific connectivity"),
      n = c(length(top_conv), length(top_spec)),
      n_surface = c(sum(is_surf(top_conv)), sum(is_surf(top_spec))),
      surface_genes = c(paste(sort(top_conv[is_surf(top_conv)]), collapse = "; "),
                        paste(sort(top_spec[is_surf(top_spec)]), collapse = "; ")))
    write_csv(cmp2, P("tables", "surfaceome_top20_comparison.csv"))
    log_msg("Top ", top_n, " by each ranking:")
    print(as.data.frame(cmp2))

    m2 <- matrix(c(sum(is_surf(top_spec)), top_n - sum(is_surf(top_spec)),
                   sum(is_surf(top_conv)), top_n - sum(is_surf(top_conv))),
                 nrow = 2)
    f2 <- fisher.test(m2)
    log_msg("Fisher exact, top ", top_n, " specific vs top ", top_n,
            " conventional: p = ", signif(f2$p.value, 3),
            ", OR = ", signif(unname(f2$estimate), 3))
  }

  # --- Targetability summary --------------------------------------------
  log_msg("---")
  log_msg("Conventional hubs: ", paste(conv_hubs, collapse = ", "))
  log_msg("  surface proteins among them: ",
          sum(is_surf(conv_hubs)), " of ", length(conv_hubs))
  log_msg("Corrected hits: ", paste(corrected_hits, collapse = ", "))
  log_msg("  surface proteins among them: ",
          sum(is_surf(corrected_hits)), " of ", length(corrected_hits))
  log_msg("If the corrected set contains surface proteins and the ",
          "conventional set does not, the bias has a therapeutic cost: it ",
          "systematically favours undruggable intracellular proteins.")

  saveRDS(list(composition = tab, fisher = ft), P("rds", "surfaceome.rds"))
  log_msg("11_surfaceome complete.")
  invisible(tab)
}
