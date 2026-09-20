# Loaders for the validation cohorts (TCGA-LGG, CGGA-325, CGGA-693).
#
# Deconvolution requires linear expression. cBioPortal z-score files are
# rejected rather than silently used, since z-scoring removes the
# between-gene magnitude information deconvolution depends on.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(tibble)
})

# AnnotationDbi masks several dplyr verbs. Bind them explicitly so this
# script is immune to package load order.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise

source("R/00b_geneid.R")

# --- Scale detection -----------------------------------------------------
# Returns "zscore", "log", or "linear".
detect_scale <- function(mat) {
  s <- sample(mat[is.finite(mat)], min(200000, sum(is.finite(mat))))
  neg  <- mean(s < 0)
  mx   <- max(s, na.rm = TRUE)
  if (neg > 0.15 && mx < 25) return("zscore")
  if (mx < 30) return("log")
  "linear"
}

prepare_for_deconv <- function(mat, label) {
  sc <- detect_scale(mat)
  log_msg(label, ": expression scale detected = ", sc)

  if (sc == "zscore") {
    stop(label, ": values look like z-scores (", round(100 * mean(mat < 0, na.rm = TRUE), 1),
         "% negative). Deconvolution requires linear expression. ",
         "Use an RSEM/FPKM/TPM file for this cohort, not the z-score file.")
  }
  if (sc == "log") {
    log_msg(label, ": un-logging for deconvolution")
    mat <- 2^mat - 1
  }
  mat[mat < 0] <- 0
  mat[is.na(mat)] <- 0
  mat
}

# --- Generic table reader ------------------------------------------------
read_matrix_file <- function(path, label) {
  log_msg(label, ": reading ", basename(path))
  df <- data.table::fread(path, data.table = FALSE, check.names = FALSE)

  # Drop a cBioPortal Entrez column if present
  drop <- intersect(names(df), c("Entrez_Gene_Id", "ENTREZ_GENE_ID",
                                 "Entrez_Gene_ID"))
  if (length(drop)) df <- df[, setdiff(names(df), drop), drop = FALSE]

  mat <- normalise_to_symbols(df, id_col = 1, label = label)
  mat
}

# --- TCGA-LGG ------------------------------------------------------------
load_tcga_lgg_matrix <- function() {
  cache <- P("rds", "cohort_TCGA_LGG.rds")
  if (file.exists(cache)) return(readRDS(cache))

  f <- find_local(LOCAL$tcga_lgg_expr)
  if (is.na(f)) { log_msg("TCGA-LGG expression not found; skipping."); return(NULL) }

  mat <- tryCatch(read_matrix_file(f, "TCGA-LGG"),
                  error = function(e) { warning(conditionMessage(e)); NULL })
  if (is.null(mat)) return(NULL)

  mat <- tryCatch(prepare_for_deconv(mat, "TCGA-LGG"),
                  error = function(e) {
                    log_msg("TCGA-LGG rejected: ", conditionMessage(e))
                    NULL })
  if (is.null(mat)) return(NULL)

  saveRDS(mat, cache)
  mat
}

# --- CGGA ----------------------------------------------------------------
# Two cohorts, RSEM gene-level files already in the data directory.
load_cgga_matrix <- function(which = c("325", "693")) {
  which <- match.arg(which)
  cache <- P("rds", paste0("cohort_CGGA_", which, ".rds"))
  if (file.exists(cache)) return(readRDS(cache))

  stem <- paste0("CGGA.mRNAseq_", which, ".RSEM-genes")
  f <- find_local(stem)
  if (is.na(f)) { log_msg("CGGA-", which, " not found; skipping."); return(NULL) }

  mat <- tryCatch(read_matrix_file(f, paste0("CGGA-", which)),
                  error = function(e) { warning(conditionMessage(e)); NULL })
  if (is.null(mat)) return(NULL)

  mat <- tryCatch(prepare_for_deconv(mat, paste0("CGGA-", which)),
                  error = function(e) {
                    log_msg("CGGA-", which, " rejected: ", conditionMessage(e))
                    NULL })
  if (is.null(mat)) return(NULL)

  saveRDS(mat, cache)
  mat
}

# --- CGGA clinical, for grade and IDH -----------------------------------
load_cgga_clinical <- function(which = c("325", "693")) {
  which <- match.arg(which)
  stem <- paste0("CGGA.mRNAseq_", which, "_clinical")
  f <- find_local(stem)
  if (is.na(f)) { log_msg("CGGA-", which, " clinical not found."); return(NULL) }

  cl <- data.table::fread(f, data.table = FALSE, check.names = FALSE)
  names(cl) <- toupper(gsub("[^A-Za-z0-9]", "_", names(cl)))

  id  <- names(cl)[1]
  grd <- grep("GRADE", names(cl), value = TRUE)[1]
  idh <- grep("IDH", names(cl), value = TRUE)[1]
  hist <- grep("HISTOLOGY", names(cl), value = TRUE)[1]

  if (is.na(grd)) { log_msg("CGGA-", which, ": no grade column."); return(NULL) }

  g <- toupper(trimws(as.character(cl[[grd]])))
  grade <- rep(NA_character_, length(g))
  grade[g %in% c("WHO II", "WHO_II", "II", "2")]  <- "LGG"
  grade[g %in% c("WHO III", "WHO_III", "III", "3")] <- "LGG"   # see note below
  grade[g %in% c("WHO IV", "WHO_IV", "IV", "4")]  <- "HGG"

  # NOTE: WHO grade III is treated as lower-grade here to match the
  # "lower grade glioma" definition used by GSE147352, which pools grade II
  # and III against grade IV. This is stated explicitly in the methods.

  tibble(
    sample = as.character(cl[[id]]),
    grade  = grade,
    IDH    = if (!is.na(idh)) as.character(cl[[idh]]) else NA_character_,
    histology = if (!is.na(hist)) as.character(cl[[hist]]) else NA_character_
  ) %>% filter(!is.na(grade))
}

# --- Registry ------------------------------------------------------------
# Each entry returns a list(matrix, meta) where meta has sample + grade.
get_validation_cohorts <- function() {
  out <- list()

  m <- load_tcga_lgg_matrix()
  if (!is.null(m)) {
    # TCGA-LGG is entirely lower-grade by construction
    out[["TCGA-LGG"]] <- list(
      matrix = m,
      meta = tibble(sample = colnames(m), grade = "LGG"))
  }

  for (w in c("325", "693")) {
    m <- load_cgga_matrix(w)
    if (is.null(m)) next
    cl <- load_cgga_clinical(w)
    if (is.null(cl)) {
      log_msg("CGGA-", w, ": no clinical, cannot assign grade. Skipping.")
      next
    }
    common <- intersect(colnames(m), cl$sample)
    log_msg("CGGA-", w, ": ", length(common), " samples with expression and grade")
    if (length(common) < 20) next
    out[[paste0("CGGA-", w)]] <- list(
      matrix = m[, common, drop = FALSE],
      meta   = cl %>% filter(sample %in% common))
  }

  if (length(out) == 0) {
    log_msg("No validation cohort could be loaded.")
  } else {
    for (nm in names(out)) {
      tb <- table(out[[nm]]$meta$grade)
      log_msg("Cohort ", nm, ": ",
              paste(names(tb), tb, sep = "=", collapse = ", "))
    }
  }
  out
}
