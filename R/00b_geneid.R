# Gene identifier detection and mapping to HGNC symbols.
# Handles Ensembl (with or without version suffix), Entrez, RefSeq and
# probe IDs. Duplicate symbols are collapsed to the highest mean.

suppressPackageStartupMessages({
  library(AnnotationDbi); library(org.Hs.eg.db)
})

# --- Detection -----------------------------------------------------------
detect_id_type <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & ids != ""]
  s <- head(ids, 2000)

  frac <- function(pat) mean(grepl(pat, s))

  if (frac("^ENSG[0-9]{11}") > 0.5)      return("ENSEMBL")
  if (frac("^ENST[0-9]{11}") > 0.5)      return("ENSEMBLTRANS")
  if (frac("^[0-9]+$") > 0.8)            return("ENTREZID")
  if (frac("^[0-9]+_at$|_at$") > 0.5)    return("PROBEID")
  if (frac("^N[MR]_[0-9]+") > 0.5)       return("REFSEQ")
  "SYMBOL"
}

# --- Mapping -------------------------------------------------------------
# Returns a character vector of symbols, same length as input, NA where
# no mapping exists.
map_to_symbol <- function(ids, id_type = NULL) {
  ids <- as.character(ids)
  if (is.null(id_type)) id_type <- detect_id_type(ids)

  if (id_type == "SYMBOL") return(ids)

  # Ensembl IDs frequently carry a version suffix that blocks lookup
  clean <- ids
  if (id_type %in% c("ENSEMBL", "ENSEMBLTRANS")) {
    clean <- sub("\\..*$", "", ids)
  }

  keytype <- switch(id_type,
                    ENSEMBL      = "ENSEMBL",
                    ENSEMBLTRANS = "ENSEMBLTRANS",
                    ENTREZID     = "ENTREZID",
                    REFSEQ       = "REFSEQ",
                    "SYMBOL")

  sym <- suppressMessages(tryCatch(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(clean),
                          column = "SYMBOL", keytype = keytype,
                          multiVals = "first"),
    error = function(e) {
      warning("Mapping failed for keytype ", keytype, ": ",
              conditionMessage(e)); NULL
    }))

  if (is.null(sym)) return(rep(NA_character_, length(ids)))
  unname(sym[clean])
}

# --- Matrix-level normalisation -----------------------------------------
# Takes a data.frame whose first column holds identifiers, returns a numeric
# matrix with unique HGNC symbols as rownames. Duplicate symbols are
# collapsed to the probe/transcript with the highest mean expression.
normalise_to_symbols <- function(df, id_col = 1, label = "matrix") {
  ids <- as.character(df[[id_col]])
  id_type <- detect_id_type(ids)
  log_msg(label, ": detected identifier type = ", id_type)

  sym <- map_to_symbol(ids, id_type)

  n_ok <- sum(!is.na(sym) & sym != "")
  log_msg(label, ": mapped ", n_ok, " of ", length(ids), " identifiers (",
          round(100 * n_ok / length(ids), 1), "%)")

  if (n_ok < 1000) {
    stop(label, ": only ", n_ok, " identifiers mapped to symbols. ",
         "Detected type was ", id_type, ". First few identifiers: ",
         paste(head(ids, 5), collapse = ", "),
         "\nMapping cannot proceed with this few genes.")
  }

  mat <- as.matrix(df[, -id_col, drop = FALSE])
  storage.mode(mat) <- "numeric"

  keep <- !is.na(sym) & sym != ""
  mat <- mat[keep, , drop = FALSE]
  sym <- sym[keep]

  # Collapse duplicates: keep the row with the highest mean
  rm_ <- rowMeans(mat, na.rm = TRUE)
  ord <- order(rm_, decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]
  sym <- sym[ord]
  dup <- duplicated(sym)
  mat <- mat[!dup, , drop = FALSE]
  rownames(mat) <- sym[!dup]

  log_msg(label, ": ", nrow(mat), " unique symbols x ", ncol(mat), " samples")
  mat
}
