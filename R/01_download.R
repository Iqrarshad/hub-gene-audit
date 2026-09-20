# Loads local files, downloads what is missing, and builds the sample
# manifest with per-grade counts.
#
# Grade assignment fails loudly and writes unresolved samples to a log
# rather than proceeding with guessed labels.

suppressPackageStartupMessages({
  library(GEOquery); library(dplyr); library(readr)
  library(stringr); library(tibble); library(data.table)
})

# --- Generic delimited reader -------------------------------------------
read_any <- function(path) {
  stopifnot(file.exists(path))
  data.table::fread(path, data.table = FALSE, check.names = FALSE)
}

# --- GSE147352: local DESeq2-normalised counts ---------------------------
load_gse147352 <- function() {
  rds <- P("rds", "gse147352.rds")
  if (file.exists(rds)) return(readRDS(rds))

  f <- find_local(LOCAL$gse147352_counts)
  if (is.na(f)) {
    f <- find_local("GSE147352/GSE147352_DESeq_normalized_counts.csv.gz")
  }
  if (is.na(f)) stop("GSE147352 counts not found under ", D("GSE147352"))
  log_msg("Reading GSE147352 counts: ", f)

  expr <- read_any(f)
  names(expr)[1] <- "gene"

  # Phenotype from GEO (metadata only). Cached on its own so a slow or
  # dropped link does not force re-fetching the whole dataset, and retried
  # because this small call is the one most likely to stall.
  pd_rds <- P("rds", "gse147352_pdata.rds")
  if (file.exists(pd_rds)) {
    pd <- readRDS(pd_rds)
  } else {
    log_msg("Fetching GSE147352 phenotype from GEO ...")
    gse <- NULL
    for (i in 1:5) {
      gse <- tryCatch(
        getGEO("GSE147352", GSEMatrix = TRUE, getGPL = FALSE,
               destdir = CACHE_DIR)[[1]],
        error = function(e) {
          log_msg("  phenotype attempt ", i, " failed: ",
                  conditionMessage(e)); NULL })
      if (!is.null(gse)) break
      if (i < 5) { log_msg("  waiting ", 20 * i, "s"); Sys.sleep(20 * i) }
    }
    if (is.null(gse))
      stop("GSE147352 phenotype could not be fetched. Download the series ",
           "matrix into ", CACHE_DIR, " manually and rerun.")
    pd <- Biobase::pData(gse)
    saveRDS(pd, pd_rds)
  }

  out <- list(expr = expr, pdata = pd)
  saveRDS(out, rds)
  out
}

# --- Grade assignment ----------------------------------------------------
assign_grade <- function(pd) {
  cols <- grep("characteristics|title|source_name|description",
               names(pd), ignore.case = TRUE, value = TRUE)
  blob <- apply(pd[, cols, drop = FALSE], 1,
                function(r) tolower(paste(r, collapse = " | ")))

  grade <- rep(NA_character_, length(blob))
  grade[str_detect(blob, "grade\\s*iii|grade\\s*3|anaplastic")]        <- "HGG"
  grade[str_detect(blob, "lower.?grade|\\blgg\\b|grade\\s*ii\\b|grade\\s*2\\b|astrocytoma|oligodendroglioma|oligoastrocytoma")] <- "LGG"
  grade[str_detect(blob, "glioblastoma|\\bgbm\\b|grade\\s*iv|grade\\s*4")] <- "HGG"
  grade[str_detect(blob, "normal|non-?tumou?r|control|non-?neoplastic|epilepsy")] <- "Normal"

  if (anyNA(grade)) {
    unresolved <- data.frame(geo_accession = rownames(pd)[is.na(grade)],
                             text = blob[is.na(grade)])
    write_csv(unresolved, P("logs", "unresolved_grade.csv"))
    stop(sum(is.na(grade)), " samples could not be assigned a grade. ",
         "See logs/unresolved_grade.csv and extend the rules in ",
         "assign_grade(). Do NOT proceed with guessed labels.")
  }
  grade
}

# --- Manifest ------------------------------------------------------------
build_manifest <- function() {
  g  <- load_gse147352()
  pd <- g$pdata

  manifest <- tibble(
    dataset       = "GSE147352",
    geo_accession = rownames(pd),
    title         = as.character(pd$title),
    grade         = assign_grade(pd)
  )

  tab <- table(manifest$grade)
  log_msg("GSE147352 composition: ",
          paste(names(tab), tab, sep = "=", collapse = ", "))

  expected <- c(HGG = 85, LGG = 18, Normal = 15)   # from the GEO record
  for (k in names(expected)) {
    obs <- if (k %in% names(tab)) as.integer(tab[[k]]) else 0L
    if (obs != expected[[k]]) {
      warning(sprintf("%s: found %d, GEO record states %d. Reconcile.",
                      k, obs, expected[[k]]))
    }
  }

  write_csv(manifest, P("tables", "cohort_manifest.csv"))

  # Methods sentence generated from data, never typed by hand
  sent <- sprintf(paste0(
    "Immune deconvolution and correlation analyses were performed on ",
    "GSE147352, comprising %d low-grade glioma, %d high-grade glioma and ",
    "%d normal brain samples. Hub gene versus immune cell correlations were ",
    "computed separately within each grade (LGG n = %d; HGG n = %d) and were ",
    "not pooled. The difference in sample size between grades is noted, and ",
    "confidence intervals are reported for every correlation."),
    tab[["LGG"]], tab[["HGG"]], tab[["Normal"]], tab[["LGG"]], tab[["HGG"]])
  writeLines(sent, P("logs", "methods_sample_sizes.txt"))

  invisible(manifest)
}

# --- Microarrays not held locally ---------------------------------------
# getGEO with retries. On a slow or intermittent link a single attempt
# often stalls or drops; each retry waits longer before trying again. The
# global download timeout is set in 00_config.R.
get_geo_retry <- function(acc, tries = 5, base_wait = 20) {
  for (i in seq_len(tries)) {
    e <- tryCatch(
      getGEO(acc, GSEMatrix = TRUE, getGPL = TRUE, destdir = CACHE_DIR)[[1]],
      error = function(err) {
        log_msg("  attempt ", i, " for ", acc, " failed: ",
                conditionMessage(err))
        NULL
      })
    if (!is.null(e)) return(e)
    if (i < tries) {
      wait <- base_wait * i
      log_msg("  waiting ", wait, "s before retry")
      Sys.sleep(wait)
    }
  }
  stop("Could not download ", acc, " after ", tries, " attempts. ",
       "Download the series matrix manually into ", CACHE_DIR,
       " and rerun.")
}

download_missing <- function() {
  out <- list()
  for (acc in DOWNLOAD_IF_MISSING) {
    rds <- P("rds", paste0(acc, ".rds"))
    if (file.exists(rds)) { out[[acc]] <- readRDS(rds); next }
    log_msg("Downloading ", acc, " ...")
    e <- get_geo_retry(acc)
    saveRDS(e, rds); out[[acc]] <- e
  }
  invisible(out)
}

# --- Locally held microarrays -------------------------------------------
load_gse16011 <- function() {
  rds <- P("rds", "gse16011.rds")
  if (file.exists(rds)) return(readRDS(rds))

  mf <- find_local(LOCAL$gse16011_matrix)
  if (is.na(mf)) stop("GSE16011 series matrix not found in ", DATA_DIR)
  log_msg("Reading GSE16011: ", mf)

  e <- getGEO(filename = mf, getGPL = FALSE)

  ann_f <- find_local(LOCAL$gpl8542_annot)
  ann <- if (!is.na(ann_f)) read_any(ann_f) else NULL

  out <- list(eset = e, annot = ann)
  saveRDS(out, rds)
  out
}

load_gse108474 <- function() {
  rds <- P("rds", "gse108474.rds")
  if (file.exists(rds)) return(readRDS(rds))

  ef <- find_local(LOCAL$gse108474_expr)
  cf <- find_local(LOCAL$gse108474_clin)
  if (is.na(ef)) stop("GSE108474 expression not found in ", DATA_DIR)
  log_msg("Reading GSE108474: ", ef)

  out <- list(expr = read_any(ef),
              clin = if (!is.na(cf)) read_any(cf) else NULL)
  saveRDS(out, rds)
  out
}

main_01 <- function() {
  build_manifest()
  if (INCLUDE_SMALL_ARRAYS) download_missing()
  load_gse16011()
  load_gse108474()
  log_msg("01_download complete.")
}
