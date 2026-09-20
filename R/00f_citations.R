# Papers per gene from NCBI gene2pubmed.
#
# Shared by the degree-decomposition stage and the propagation-gates stage.
# Both need a symbol-to-publication-count table. Keeping the builder here
# removes the circular dependency that existed when only the decomposition
# stage defined it: the propagation-gates stage runs first in the pipeline
# order but reads the same table, so it now builds and caches it directly
# rather than stopping.
#
# Requires gene2pubmed.gz in DATA_DIR (ftp.ncbi.nlm.nih.gov/gene/DATA/).
# Returns a tibble with columns gene, n_papers, or NULL if the raw file is
# absent. The result is cached to rds/papers_per_gene.rds and reused.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(data.table)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

load_gene2pubmed <- function() {
  cache <- P("rds", "papers_per_gene.rds")
  if (file.exists(cache)) return(readRDS(cache))

  f <- find_local("gene2pubmed")
  if (is.na(f)) {
    log_msg("gene2pubmed not found in ", DATA_DIR)
    log_msg("Download: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz")
    log_msg("Place the file (gzipped is fine) in DATA_DIR and rerun.")
    return(NULL)
  }
  log_msg("Reading ", basename(f), " (large; allow a minute)")
  g2p <- data.table::fread(f, data.table = TRUE, showProgress = FALSE)
  names(g2p)[1:3] <- c("tax_id", "GeneID", "PubMed_ID")
  g2p <- g2p[tax_id == 9606]
  log_msg("Human gene-publication links: ", nrow(g2p))

  cnt <- g2p[, .(n_papers = uniqueN(PubMed_ID)), by = GeneID]
  log_msg("Genes with at least one publication: ", nrow(cnt))

  sym <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = as.character(cnt$GeneID),
    keytype = "ENTREZID", columns = "SYMBOL"))
  out <- as_tibble(cnt) %>%
    mutate(ENTREZID = as.character(GeneID)) %>%
    inner_join(sym, by = "ENTREZID") %>%
    filter(!is.na(SYMBOL)) %>%
    group_by(gene = SYMBOL) %>%
    summarise(n_papers = max(n_papers), .groups = "drop")

  log_msg("Genes mapped to symbols: ", nrow(out),
          "; median papers ", median(out$n_papers),
          "; max ", max(out$n_papers), " (",
          out$gene[which.max(out$n_papers)], ")")
  saveRDS(out, cache)
  out
}
