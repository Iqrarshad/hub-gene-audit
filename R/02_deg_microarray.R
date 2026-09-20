# Microarray differential expression with composition adjustment.
#
# GSE15824 (7 vs 3) and GSE21354 (10 vs 4) cannot support a covariate
# model; compositional_deg detects this and falls back to the unadjusted
# fit, flagging the dataset. GSE16011 has 276 tumours against 8 normals.

suppressPackageStartupMessages({
  library(limma); library(dplyr); library(readr); library(tibble)
  library(stringr); library(Biobase); library(data.table)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00e_compadj.R")

collapse_to_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; symbols <- symbols[ord]
  dup <- duplicated(symbols)
  mat <- mat[!dup, , drop = FALSE]
  rownames(mat) <- symbols[!dup]
  mat
}

infer_groups <- function(pd) {
  cols <- grep("characteristics|title|source_name|description",
               names(pd), ignore.case = TRUE, value = TRUE)
  blob <- apply(pd[, cols, drop = FALSE], 1,
                function(r) tolower(paste(r, collapse = " | ")))
  g <- rep(NA_character_, length(blob))
  g[str_detect(blob, "glioma|astrocytoma|oligodendroglioma|tumou?r|glioblastoma|gbm")] <- "Tumour"
  g[str_detect(blob, "normal|non-?tumou?r|control|non-?neoplastic|epilepsy|healthy")] <- "Normal"
  g
}

eset_to_matrix <- function(eset, label, annot = NULL) {
  mat <- exprs(eset)
  fd  <- fData(eset)

  sym_col <- grep("^gene.?symbol$|^symbol$|^gene_symbol$", names(fd),
                  ignore.case = TRUE, value = TRUE)
  if (length(sym_col) > 0 && nrow(fd) == nrow(mat)) {
    mat <- collapse_to_symbol(mat, as.character(fd[[sym_col[1]]]))
  } else if (!is.null(annot)) {
    id_col <- names(annot)[1]
    s_col  <- grep("symbol", names(annot), ignore.case = TRUE, value = TRUE)
    if (length(s_col) > 0) {
      idx <- match(rownames(mat), as.character(annot[[id_col]]))
      mat <- collapse_to_symbol(mat, as.character(annot[[s_col[1]]])[idx])
    }
  } else {
    warning(label, ": no symbol mapping available.")
  }
  mat
}

deg_from_eset <- function(eset, label, annot = NULL) {
  mat <- eset_to_matrix(eset, label, annot)
  grp <- infer_groups(pData(eset))
  if (all(is.na(grp))) {
    warning(label, ": no groups inferred."); return(NULL)
  }
  compositional_deg(mat, grp, label)
}

# REMBRANDT sample IDs need harmonising. Expression columns look like
# "00518392_U133P2"; clinical SUBJECT_ID is "518392". Strip the platform
# suffix and the leading zeros, then match.
normalise_rembrandt_id <- function(x) {
  x <- as.character(x)
  x <- sub("_.*$", "", x)          # drop _U133P2 and similar
  x <- sub("^0+", "", x)           # drop leading zeros
  trimws(x)
}

deg_gse108474 <- function() {
  obj  <- load_gse108474()
  expr <- obj$expr
  mat  <- as.matrix(expr[, -1, drop = FALSE])
  rownames(mat) <- as.character(expr[[1]])
  storage.mode(mat) <- "numeric"

  clin <- obj$clin
  if (is.null(clin)) { warning("GSE108474: no clinical."); return(NULL) }

  samp_norm <- normalise_rembrandt_id(colnames(mat))
  clin_norm <- normalise_rembrandt_id(clin$SUBJECT_ID)
  idx <- match(samp_norm, clin_norm)
  log_msg("GSE108474: matched ", sum(!is.na(idx)), " of ", ncol(mat),
          " expression columns to clinical records")
  if (sum(!is.na(idx)) < 50) {
    stop("GSE108474: ID harmonisation failed. Expression example: ",
         colnames(mat)[1], " -> ", samp_norm[1],
         "; clinical example: ", clin$SUBJECT_ID[1])
  }

  dt <- toupper(trimws(as.character(clin$DISEASE_TYPE)))[idx]
  grp <- rep(NA_character_, length(dt))
  grp[dt %in% c("ASTROCYTOMA", "GBM", "OLIGODENDROGLIOMA", "MIXED")] <- "Tumour"
  grp[dt == "NON_TUMOR"] <- "Normal"

  tb <- table(grp, useNA = "ifany")
  log_msg("GSE108474 groups: ",
          paste(names(tb), tb, sep = "=", collapse = ", "))

  # Probe IDs are Affymetrix; map to symbols via the GPL570 annotation if a
  # symbol column is not already present.
  if (grepl("_at$", rownames(mat)[1])) {
    ann_f <- find_local("GPL570_annot")
    if (!is.na(ann_f)) {
      ann <- data.table::fread(ann_f, data.table = FALSE, check.names = FALSE)
      id_c <- names(ann)[1]
      s_c  <- grep("symbol", names(ann), ignore.case = TRUE, value = TRUE)
      if (length(s_c) > 0) {
        i <- match(rownames(mat), as.character(ann[[id_c]]))
        mat <- collapse_to_symbol(mat, as.character(ann[[s_c[1]]])[i])
        log_msg("GSE108474: mapped to ", nrow(mat), " symbols via GPL570")
      }
    } else {
      warning("GSE108474: GPL570 annotation not found; probe IDs retained. ",
              "This dataset cannot join the symbol-level intersection.")
      return(NULL)
    }
  }

  compositional_deg(mat, grp, "GSE108474")
}

# Grade contrast within REMBRANDT, the composition-clean comparison
deg_gse108474_grade <- function() {
  obj  <- load_gse108474()
  expr <- obj$expr
  mat  <- as.matrix(expr[, -1, drop = FALSE])
  rownames(mat) <- as.character(expr[[1]])
  storage.mode(mat) <- "numeric"
  clin <- obj$clin
  if (is.null(clin)) return(NULL)

  idx <- match(normalise_rembrandt_id(colnames(mat)),
               normalise_rembrandt_id(clin$SUBJECT_ID))
  wg <- toupper(trimws(as.character(clin$WHO_GRADE)))[idx]
  grp <- rep(NA_character_, length(wg))
  grp[wg %in% c("II", "III")] <- "LGG"
  grp[wg == "IV"] <- "HGG"

  tb <- table(grp, useNA = "ifany")
  log_msg("GSE108474 grade groups: ",
          paste(names(tb), tb, sep = "=", collapse = ", "))
  if (sum(grp == "LGG", na.rm = TRUE) < 10 ||
      sum(grp == "HGG", na.rm = TRUE) < 10) return(NULL)

  if (grepl("_at$", rownames(mat)[1])) {
    ann_f <- find_local("GPL570_annot")
    if (is.na(ann_f)) return(NULL)
    ann <- data.table::fread(ann_f, data.table = FALSE, check.names = FALSE)
    s_c <- grep("symbol", names(ann), ignore.case = TRUE, value = TRUE)
    if (length(s_c) == 0) return(NULL)
    i <- match(rownames(mat), as.character(ann[[names(ann)[1]]]))
    mat <- collapse_to_symbol(mat, as.character(ann[[s_c[1]]])[i])
  }
  compositional_deg(mat, grp, "GSE108474_HGG_vs_LGG",
                    contrast = c("HGG", "LGG"))
}

deg_gse16011 <- function() {
  obj <- load_gse16011()
  deg_from_eset(obj$eset, "GSE16011", obj$annot)
}

main_02 <- function() {
  res <- list()

  res$GSE108474 <- tryCatch(deg_gse108474(),
    error = function(e) { warning("GSE108474: ", conditionMessage(e)); NULL })

  # Grade contrast, saved separately as an additional composition control
  gr <- tryCatch(deg_gse108474_grade(),
    error = function(e) { warning("GSE108474 grade: ", conditionMessage(e)); NULL })
  if (!is.null(gr)) {
    write_csv(gr, P("tables", "deg_GSE108474_HGG_vs_LGG_full.csv"))
    log_msg("GSE108474 grade contrast median attenuation: ",
            round(median(gr$attenuation[gr$sig_raw], na.rm = TRUE), 1), "%")
  }

  res$GSE16011 <- tryCatch(deg_gse16011(),
    error = function(e) { warning("GSE16011: ", conditionMessage(e)); NULL })

  small <- download_missing()
  for (acc in names(small)) {
    res[[acc]] <- tryCatch(deg_from_eset(small[[acc]], acc),
      error = function(e) { warning(acc, ": ", conditionMessage(e)); NULL })
  }

  res <- res[!vapply(res, is.null, logical(1))]
  if (length(res) == 0) stop("No microarray DEG succeeded.")

  all_res <- bind_rows(res)
  write_csv(all_res, P("tables", "deg_microarray_composition_adjusted.csv"))

  for (nm in names(res)) {
    write_csv(res[[nm]], P("tables", paste0("deg_", nm, "_full.csv")))
  }

  # Attenuation summary: how much of each dataset is composition driven
  summ <- all_res %>% group_by(dataset) %>%
    summarise(n_tested = n(),
              n_sig_raw = sum(sig_raw, na.rm = TRUE),
              n_robust  = sum(composition_robust, na.rm = TRUE),
              n_driven  = sum(composition_driven, na.rm = TRUE),
              pct_lost  = round(100 * sum(composition_driven, na.rm = TRUE) /
                                pmax(sum(sig_raw, na.rm = TRUE), 1), 1),
              median_attenuation = round(median(attenuation[sig_raw],
                                                na.rm = TRUE), 1),
              .groups = "drop")
  write_csv(summ, P("tables", "composition_adjustment_summary.csv"))

  log_msg("Composition adjustment summary:")
  print(as.data.frame(summ))

  # Where the original nine hub genes land
  hubs <- all_res %>% filter(gene %in% HUB_GENES_PRIOR) %>%
    select(dataset, gene, logFC_raw, q_raw, logFC_adj, q_adj,
           attenuation, composition_robust, composition_driven) %>%
    arrange(gene, dataset)
  write_csv(hubs, P("tables", "hub_genes_composition_check.csv"))
  log_msg("Original hub genes under composition adjustment:")
  print(as.data.frame(hubs))

  saveRDS(res, P("rds", "deg_microarray.rds"))
  log_msg("02_deg_microarray complete.")
  invisible(res)
}
