# GSE147352 differential expression with composition adjustment.
#
# The authors supply DESeq2-normalised counts rather than raw counts, so
# limma-trend on log2 values is used instead of the negative binomial model.
#
# The grade contrast (HGG vs LGG) adjusts for immune content rather than
# neuronal content: both arms are tumour tissue, so leukocyte infiltration
# is the relevant confounder.

suppressPackageStartupMessages({
  library(limma); library(dplyr); library(readr); library(tibble)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00g_expression.R")
source("R/00e_compadj.R")

map_columns <- function(mat, manifest) {
  cn <- colnames(mat)
  if (all(manifest$geo_accession %in% cn))
    return(unname(setNames(manifest$grade, manifest$geo_accession)[cn]))
  ttl <- trimws(manifest$title)
  if (sum(ttl %in% cn) > 0.8 * length(ttl))
    return(unname(setNames(manifest$grade, ttl)[cn]))
  stop("Could not map counts columns to manifest samples.")
}

main_03 <- function() {
  mat      <- load_expression_matrix()
  manifest <- read_csv(P("tables", "cohort_manifest.csv"), show_col_types = FALSE)
  grade    <- map_columns(mat, manifest)

  keep <- !is.na(grade)
  mat <- mat[, keep, drop = FALSE]; grade <- grade[keep]

  keep_g <- rowMeans(mat > 1) >= 0.2
  log_msg("Retaining ", sum(keep_g), " of ", nrow(mat), " genes")
  mat <- mat[keep_g, , drop = FALSE]

  contrasts <- list(
    LGG_vs_Normal = c("LGG", "Normal"),
    HGG_vs_Normal = c("HGG", "Normal"),
    HGG_vs_LGG    = c("HGG", "LGG"))

  # Tumour-vs-normal contrasts differ mainly in neuron and glia content.
  # The grade contrast is tumour against tumour, so neurons are not the
  # confound: leukocyte infiltration is. Immune content is therefore added
  # as a covariate there, and the unadjusted grade result is kept alongside
  # for comparison.
  adj_for <- list(
    LGG_vs_Normal = c("neuronal", "glial"),
    HGG_vs_Normal = c("neuronal", "glial"),
    HGG_vs_LGG    = c("neuronal", "glial", "immune"))

  out <- list()
  for (nm in names(contrasts)) {
    ct <- contrasts[[nm]]
    sub <- grade %in% ct
    r <- compositional_deg(mat[, sub, drop = FALSE], grade[sub],
                           paste0("GSE147352_", nm), contrast = ct,
                           adjust_for = adj_for[[nm]])
    if (is.null(r)) next
    write_csv(r, P("tables", paste0("deg_GSE147352_", nm, "_full.csv")))
    out[[nm]] <- r
  }

  if (length(out) == 0) stop("No RNA-seq contrast succeeded.")

  # The grade contrast is the clean one: both arms are tumour, so neuronal
  # dilution cannot manufacture the difference the way it does against
  # normal brain. Attenuation here should be far smaller.
  att <- bind_rows(out) %>% group_by(dataset) %>%
    summarise(median_attenuation = round(median(attenuation[sig_raw],
                                                na.rm = TRUE), 1),
              n_sig_raw = sum(sig_raw, na.rm = TRUE),
              n_robust = sum(composition_robust, na.rm = TRUE),
              .groups = "drop")
  write_csv(att, P("tables", "rnaseq_contrast_attenuation.csv"))
  log_msg("Attenuation by contrast (tumour-vs-normal should exceed grade contrast):")
  print(as.data.frame(att))

  saveRDS(out, P("rds", "deg_rnaseq.rds"))
  log_msg("03_deg_rnaseq complete.")
  invisible(out)
}
