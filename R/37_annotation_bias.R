# Annotation bias, the enrichment-side analogue of the degree bias.
#
# The degree-decomposition stage regresses interactome degree on publication
# count and finds most of the variance explained by study attention. The same
# attention shapes the annotation record: a gene that has been studied more
# accumulates more Gene Ontology and pathway terms, independently of its
# biology. This stage measures that dependence and asks whether it carries
# through to enrichment.
#
# Two quantities are reported.
#
# First, per-gene annotation count against publication count. For each gene
# the number of GO biological-process terms and KEGG pathways it is annotated
# to is counted, and log10(1 + terms) is regressed on log10(1 + papers),
# mirroring the degree decomposition. The variance explained bounds how much
# of the annotation record tracks attention rather than biology.
#
# Second, the enrichment analogue of the degree-predicts-hubs test. For the
# consensus hub set, every KEGG pathway the hubs touch is tested, and a
# pathway's prior attention, the mean publication count of its member genes,
# is used to predict whether that pathway is called enriched. An area under
# the curve near one means the enriched terms are the well-studied ones, so
# enrichment restates prior attention in the same way hub nomination does.
# This second quantity is bounded rather than decisive: the consensus hub set
# is small, so the number of tested pathways is limited, and the figure is
# reported with its denominator.
#
# gene2pubmed and the KEGG annotation of org.Hs.eg.db are optional. If either
# is absent the stage records that and returns without error.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(AnnotationDbi); library(org.Hs.eg.db); library(pROC)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00f_citations.R")

to_entrez <- function(symbols) {
  suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
    columns = "ENTREZID"))
}

# GO biological-process term count per Entrez gene.
go_bp_counts <- function(entrez) {
  m <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(entrez), keytype = "ENTREZID",
    columns = c("GO", "ONTOLOGY")))
  m <- m[!is.na(m$GO) & m$ONTOLOGY == "BP", , drop = FALSE]
  tapply(m$GO, m$ENTREZID, function(x) length(unique(x)))
}

# KEGG pathway count per Entrez gene, from the legacy mapping. Returns NULL if
# the mapping is not present in this org.Hs.eg.db build.
kegg_counts <- function() {
  if (!exists("org.Hs.egPATH")) return(NULL)
  lst <- tryCatch(as.list(org.Hs.egPATH), error = function(e) NULL)
  if (is.null(lst)) return(NULL)
  vapply(lst, function(x) if (all(is.na(x))) 0L else length(x), integer(1))
}

# Mean prior attention of a KEGG pathway, from the publication counts of its
# member genes. path ids are given without the hsa prefix.
pathway_attention <- function(path_ids, papers_by_entrez) {
  if (!exists("org.Hs.egPATH2EG")) return(rep(NA_real_, length(path_ids)))
  p2e <- tryCatch(as.list(org.Hs.egPATH2EG), error = function(e) NULL)
  if (is.null(p2e)) return(rep(NA_real_, length(path_ids)))
  vapply(path_ids, function(pid) {
    eg <- p2e[[pid]]
    if (is.null(eg)) return(NA_real_)
    v <- papers_by_entrez[eg]
    v <- v[!is.na(v)]
    if (!length(v)) return(NA_real_)
    mean(log10(1 + v))
  }, numeric(1))
}

main_37 <- function() {
  pap <- load_gene2pubmed()
  if (is.null(pap)) {
    log_msg("gene2pubmed not available; skipping the annotation-bias stage.")
    return(invisible(NULL))
  }

  # ---- Part 1: annotation count against publication count ---------------
  map <- to_entrez(pap$gene)
  map <- map[!is.na(map$ENTREZID), , drop = FALSE]
  go_ct <- go_bp_counts(map$ENTREZID)
  kg_ct <- kegg_counts()
  if (is.null(kg_ct))
    log_msg("KEGG per-gene mapping (org.Hs.egPATH) absent; using GO only.")

  d <- tibble(gene = map$SYMBOL, ENTREZID = map$ENTREZID) %>%
    mutate(
      n_go   = as.integer(go_ct[ENTREZID]),
      n_kegg = if (is.null(kg_ct)) 0L else as.integer(kg_ct[ENTREZID])) %>%
    mutate(across(c(n_go, n_kegg), ~ ifelse(is.na(.x), 0L, .x)),
           n_terms = n_go + n_kegg) %>%
    inner_join(pap, by = "gene") %>%
    filter(n_terms > 0)

  fit <- lm(log10(1 + n_terms) ~ log10(1 + n_papers), data = d)
  r2  <- summary(fit)$r.squared
  rho <- suppressWarnings(
    cor(log10(1 + d$n_terms), log10(1 + d$n_papers), method = "spearman"))
  log_msg("Annotation count against publication count: R2 ", round(r2, 3),
          ", Spearman ", round(rho, 3), " across ", nrow(d), " genes")

  write_csv(d %>% select(gene, n_papers, n_go, n_kegg, n_terms),
            P("tables", "annotation_bias_decomposition.csv"))

  # ---- Part 2: enrichment analogue on the consensus hub set -------------
  auc_enr <- NA_real_; n_enr <- NA_integer_; n_tested <- NA_integer_
  net <- readRDS(P("rds", "network.rds"))
  hubs <- net$consensus$gene

  if (requireNamespace("clusterProfiler", quietly = TRUE) && !is.null(kg_ct)) {
    inter <- readRDS(P("rds", "intersect.rds"))
    micro <- readRDS(P("rds", "deg_microarray.rds"))
    rna   <- readRDS(P("rds", "deg_rnaseq.rds"))
    universe <- unique(c(unlist(lapply(inter$sets, function(s) s$gene)),
                         unlist(lapply(micro, function(x)
                           if (!is.null(x)) x$gene)),
                         rna$LGG_vs_Normal$gene))
    universe <- universe[!is.na(universe) & universe != ""]

    hub_eg <- to_entrez(hubs)$ENTREZID
    uni_eg <- to_entrez(universe)$ENTREZID
    ek <- tryCatch(clusterProfiler::enrichKEGG(
            gene = hub_eg[!is.na(hub_eg)],
            universe = uni_eg[!is.na(uni_eg)], organism = "hsa",
            pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH"),
          error = function(e) NULL)
    dk <- if (is.null(ek)) NULL else as.data.frame(ek)

    if (!is.null(dk) && nrow(dk) >= 5) {
      papers_by_eg <- setNames(
        pap$n_papers[match(map$SYMBOL, pap$gene)], map$ENTREZID)
      papers_by_eg <- papers_by_eg[
        !is.na(names(papers_by_eg)) & !is.na(papers_by_eg)]

      pid <- sub("^hsa", "", dk$ID)
      dk$attention <- pathway_attention(pid, papers_by_eg)
      dk$enriched  <- as.integer(dk$qvalue < 0.05)
      dk <- dk[is.finite(dk$attention), , drop = FALSE]

      n_tested <- nrow(dk); n_enr <- sum(dk$enriched)
      if (n_enr >= 3 && n_enr < n_tested) {
        auc_enr <- as.numeric(pROC::auc(pROC::roc(
          dk$enriched, dk$attention, quiet = TRUE)))
        log_msg("Enrichment analogue: prior attention predicts which of ",
                n_tested, " tested pathways are enriched at AUC ",
                round(auc_enr, 3), " (", n_enr, " enriched)")
      } else {
        log_msg("Enrichment analogue underpowered: ", n_enr, " enriched of ",
                n_tested, " tested. AUC not reported; the consensus hub set ",
                "is too small to support the term-level test.")
      }
      write_csv(dk %>% select(ID, Description, qvalue, attention, enriched),
                P("tables", "annotation_bias_enrichment.csv"))
    } else {
      log_msg("Too few KEGG pathways tested for the consensus hubs; the ",
              "enrichment analogue is not computed.")
    }
  } else {
    log_msg("clusterProfiler or KEGG mapping absent; the enrichment analogue ",
            "is skipped, and only the annotation-count decomposition is ",
            "reported.")
  }

  write_csv(tibble(
    n_genes = nrow(d), r2_terms_on_papers = round(r2, 4),
    spearman_terms_papers = round(rho, 4),
    enrichment_analogue_auc = round(auc_enr, 4),
    n_pathways_enriched = n_enr, n_pathways_tested = n_tested),
    P("tables", "annotation_bias_summary.csv"))

  # ---- Supplementary figure ---------------------------------------------
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    sd_ <- if (nrow(d) > 6000) d[sample(nrow(d), 6000), ] else d
    p <- ggplot2::ggplot(sd_,
        ggplot2::aes(log10(1 + n_papers), log10(1 + n_terms))) +
      ggplot2::geom_point(alpha = 0.15, size = 0.6,
                          colour = FIG$okabe_ito[5]) +
      ggplot2::geom_smooth(method = "lm", se = FALSE,
                           colour = FIG$okabe_ito[6], linewidth = 0.9) +
      ggplot2::annotate("text", x = min(log10(1 + sd_$n_papers)),
        y = max(log10(1 + sd_$n_terms)), hjust = 0, vjust = 1,
        label = paste0("R2 = ", round(r2, 2),
                       "\nSpearman = ", round(rho, 2)),
        size = 3.4, family = FIG$font) +
      ggplot2::labs(
        title = "Annotation density tracks publication count",
        x = "log10(1 + publications)",
        y = "log10(1 + GO BP and KEGG terms)") +
      ggplot2::theme_classic(base_size = 12, base_family = FIG$font) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,
                                           size = 13),
        axis.text = ggplot2::element_text(colour = "black"))

    for (fmt in FIG$formats) {
      f <- P("figures", paste0("supp_fig_13_annotation_bias.", fmt))
      tryCatch({
        if (fmt == "tiff")
          ggplot2::ggsave(f, p, width = 6, height = 4.6, dpi = FIG$dpi,
                          device = "tiff", compression = "lzw")
        else
          ggplot2::ggsave(f, p, width = 6, height = 4.6, device = cairo_pdf)
      }, error = function(e)
        log_msg("  ", fmt, " failed for supp_fig_13: ", conditionMessage(e)))
    }
    log_msg("Figure written: supp_fig_13_annotation_bias")
  }

  saveRDS(list(decomposition = d, r2 = r2, rho = rho,
               enrichment_auc = auc_enr), P("rds", "annotation_bias.rds"))
  log_msg("37_annotation_bias complete.")
  invisible(NULL)
}
