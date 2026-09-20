# KEGG and GO enrichment via clusterProfiler.
#
# clusterProfiler is used rather than Enrichr because it versions the
# pathway database, applies BH correction natively, and takes an explicit
# background universe.

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(dplyr); library(readr); library(ggplot2); library(tibble)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

to_entrez <- function(symbols) {
  suppressWarnings(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  )
}

run_enrichment <- function(genes, universe_genes, label) {
  map_g <- to_entrez(genes)
  map_u <- to_entrez(universe_genes)
  log_msg(label, ": ", nrow(map_g), " of ", length(genes),
          " genes mapped to Entrez")

  kegg <- enrichKEGG(gene          = map_g$ENTREZID,
                     universe      = map_u$ENTREZID,
                     organism      = "hsa",
                     pvalueCutoff  = THRESH$enrich_qval,
                     qvalueCutoff  = THRESH$enrich_qval,
                     pAdjustMethod = "BH")

  go <- enrichGO(gene          = map_g$ENTREZID,
                 universe      = map_u$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pvalueCutoff  = THRESH$enrich_qval,
                 qvalueCutoff  = THRESH$enrich_qval,
                 pAdjustMethod = "BH",
                 readable      = TRUE)

  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    k <- setReadable(kegg, org.Hs.eg.db, keyType = "ENTREZID")
    write_csv(as.data.frame(k), P("tables", paste0("kegg_", label, ".csv")))
    log_msg(label, ": ", nrow(as.data.frame(k)),
            " KEGG pathways at BH q < ", THRESH$enrich_qval)
  } else {
    log_msg(label, ": no KEGG pathways survive BH correction. ",
            "This is a result and must be reported as such.")
    write_csv(tibble(note = "No KEGG pathway significant after BH correction"),
              P("tables", paste0("kegg_", label, ".csv")))
  }

  if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
    write_csv(as.data.frame(go), P("tables", paste0("go_bp_", label, ".csv")))
    log_msg(label, ": ", nrow(as.data.frame(go)), " GO BP terms")
  }

  # Immune-relevant KEGG pathways, flagged rather than cherry-picked
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kd <- as.data.frame(kegg)
    immune_re <- paste("Th1 and Th2|Th17|T cell receptor|B cell receptor",
                       "Notch|Wnt|MAPK|VEGF|ErbB|Sphingolipid|Endocytosis",
                       "cytokine|chemokine|antigen|NF-kappa|JAK-STAT",
                       "PD-L1|checkpoint|Natural killer|Toll-like",
                       sep = "|")
    imm <- kd %>% filter(grepl(immune_re, Description, ignore.case = TRUE))
    write_csv(imm, P("tables", paste0("kegg_immune_", label, ".csv")))
    log_msg(label, ": ", nrow(imm), " immune-related pathways among these")
  }

  list(kegg = kegg, go = go)
}

main_05 <- function() {
  inter  <- readRDS(P("rds", "intersect.rds"))
  shared <- inter$shared$gene

  # Background universe: every gene tested in any dataset, not the genome
  universe <- unique(unlist(lapply(inter$sets, function(s) s$gene)))
  micro <- readRDS(P("rds", "deg_microarray.rds"))
  rna   <- readRDS(P("rds", "deg_rnaseq.rds"))
  universe <- unique(c(universe,
                       unlist(lapply(micro, function(x) if (!is.null(x)) x$gene)),
                       rna$LGG_vs_Normal$gene))
  universe <- universe[!is.na(universe) & universe != ""]

  log_msg("Enrichment universe: ", length(universe), " genes")
  res <- run_enrichment(shared, universe, "shared_degs")

  # Grade contrast enrichment
  hgg_lgg <- rna$HGG_vs_LGG
  hgg_lgg <- if ("composition_robust" %in% names(hgg_lgg))
    hgg_lgg %>% filter(composition_robust) else
    hgg_lgg %>% filter(abs(logFC_raw) >= THRESH$deg_logfc,
                       q_raw < THRESH$deg_adjp)
  if (nrow(hgg_lgg) > 10) {
    run_enrichment(hgg_lgg$gene, universe, "HGG_vs_LGG")
  }

  saveRDS(res, P("rds", "enrichment.rds"))
  log_msg("05_enrichment complete.")
  invisible(res)
}
