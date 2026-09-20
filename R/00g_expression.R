# Expression matrix loading and symbol verification.
#
# Extracted from the deconvolution stage so that the core RNA-seq DEG stage
# can load the expression matrix without depending on the immune analysis.
# Both the DEG stage and, in the immune build, the deconvolution stage source
# this file.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr)
})

source("R/00b_geneid.R")

load_expression_matrix <- function() {
  cache <- P("rds", "expr_symbols.rds")
  if (file.exists(cache)) return(readRDS(cache))

  g <- load_gse147352()
  mat <- normalise_to_symbols(g$expr, id_col = 1, label = "GSE147352")

  mx <- max(mat, na.rm = TRUE)
  if (mx < 30) {
    log_msg("Values look log-scaled (max ", round(mx, 2), "). Un-logging.")
    mat <- 2^mat - 1
  }
  mat[mat < 0] <- 0
  mat[is.na(mat)] <- 0

  saveRDS(mat, cache)
  mat
}

verify_symbols <- function(mat) {
  known <- c("ACTB", "GAPDH", "PTPRC", "CD3E", "CD68", "EGFR", "GFAP",
             "VIM", "CD44", "MBP", "SNAP25")
  found <- intersect(known, rownames(mat))
  log_msg("Symbol sanity check: ", length(found), " of ", length(known),
          " reference genes present (", paste(found, collapse = ", "), ")")
  if (length(found) < 5) {
    stop("Rownames do not look like HGNC symbols. First few rownames: ",
         paste(head(rownames(mat), 5), collapse = ", "))
  }
  invisible(TRUE)
}
