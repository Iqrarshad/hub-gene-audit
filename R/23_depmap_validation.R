# Two tests.
#
# 1. Whether a degree-orthogonal residual of the expression signal is
#    applicable at all: it changes nothing unless effect size correlates
#    with degree in the first place.
# 2. DepMap CRISPR validation, measured as glioma-selective dependency
#    rather than raw essentiality. Pan-essential genes are essential
#    everywhere and are what a degree-biased method nominates, so raw
#    essentiality would reward the bias.
#
# Requires CRISPRGeneEffect.csv and Model.csv from depmap.org; test 1 runs
# without them.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(tidyr); library(data.table); library(pROC)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# --- TEST 1: DOS viability ----------------------------------------------
dos_viability <- function() {
  norm <- read_csv(P("tables", "centrality_degree_normalised.csv"),
                   show_col_types = FALSE)
  gd <- norm %>% filter(!is.na(global_degree), global_degree > 0) %>%
    select(gene, global_degree)

  # Effect sizes from every contrast available
  rna <- readRDS(P("rds", "deg_rnaseq.rds"))
  micro <- readRDS(P("rds", "deg_microarray.rds"))

  eff <- list()
  for (nm in names(rna)) {
    d <- rna[[nm]]
    if (is.null(d)) next
    eff[[paste0("GSE147352_", nm)]] <- d %>%
      transmute(gene, effect = abs(logFC_raw),
                effect_adj = abs(logFC_adj))
  }
  for (nm in names(micro)) {
    d <- micro[[nm]]
    if (is.null(d)) next
    eff[[nm]] <- d %>% transmute(gene, effect = abs(logFC_raw),
                                 effect_adj = abs(logFC_adj))
  }

  res <- imap_dfr(eff, function(d, nm) {
    j <- d %>% inner_join(gd, by = "gene") %>%
      filter(is.finite(effect), global_degree > 0)
    if (nrow(j) < 100) return(tibble())
    r_raw <- suppressWarnings(cor(j$effect, log10(j$global_degree + 1),
                                  method = "spearman", use = "complete.obs"))
    r_adj <- if (all(is.na(j$effect_adj))) NA_real_ else
      suppressWarnings(cor(j$effect_adj, log10(j$global_degree + 1),
                           method = "spearman", use = "complete.obs"))
    fit <- lm(effect ~ log10(global_degree + 1), data = j)
    tibble(contrast = nm, n_genes = nrow(j),
           rho_effect_vs_degree = r_raw,
           rho_adjusted_effect_vs_degree = r_adj,
           r_squared = summary(fit)$r.squared,
           slope_p = summary(fit)$coefficients[2, 4])
  })

  write_csv(res, P("tables", "dos_viability.csv"))

  log_msg("=================================================")
  log_msg("TEST 1: IS DOS APPLICABLE?")
  log_msg("  DOS removes the part of the expression effect explained by")
  log_msg("  degree. If effect size does not track degree, there is")
  log_msg("  nothing to remove and the metric does nothing.")
  log_msg("=================================================")
  print(as.data.frame(res %>% mutate(across(where(is.numeric),
                                            ~ signif(.x, 3)))))

  med_rho <- median(abs(res$rho_effect_vs_degree), na.rm = TRUE)
  med_r2  <- median(res$r_squared, na.rm = TRUE)
  log_msg("Median |rho| between effect size and log degree: ",
          round(med_rho, 3))
  log_msg("Median R squared: ", round(med_r2, 4),
          " (the proportion of effect size DOS would remove)")

  if (med_r2 < 0.01) {
    log_msg("VERDICT: degree explains under 1 percent of effect size ",
            "variance. DOS residuals would be almost identical to the ",
            "input. The metric is a no-op on this data and the second ",
            "paper has no basis.")
  } else if (med_r2 < 0.05) {
    log_msg("VERDICT: degree explains ", round(100 * med_r2, 1),
            " percent of effect size variance. DOS would change rankings ",
            "only marginally. Weak basis for a method paper.")
  } else {
    log_msg("VERDICT: degree explains ", round(100 * med_r2, 1),
            " percent of effect size variance. DOS has something real to ",
            "remove and is worth developing.")
  }
  list(results = res, median_r2 = med_r2, gd = gd)
}

# --- DOS score ------------------------------------------------------------
compute_dos <- function(deg_table, gd, effect_col = "logFC_raw") {
  j <- deg_table %>%
    transmute(gene, effect = abs(.data[[effect_col]])) %>%
    inner_join(gd, by = "gene") %>%
    filter(is.finite(effect), global_degree > 0)
  if (nrow(j) < 50) return(NULL)
  fit <- MASS::rlm(effect ~ log10(global_degree + 1), data = j,
                   maxit = 100)
  j$dos <- residuals(fit)
  j %>% arrange(desc(dos))
}

# --- TEST 2: DepMap -------------------------------------------------------
load_depmap <- function() {
  ge <- find_local("CRISPRGeneEffect")
  md <- find_local("Model")
  if (is.na(ge) || is.na(md)) {
    log_msg("DepMap files not found in ", DATA_DIR)
    log_msg("Download from https://depmap.org/portal/download/ :")
    log_msg("  CRISPRGeneEffect.csv   and   Model.csv")
    log_msg("Place both in DATA_DIR and rerun. Test 1 above is unaffected.")
    return(NULL)
  }
  log_msg("Reading ", basename(ge), " (large file, allow a minute)")
  eff <- data.table::fread(ge, data.table = FALSE, check.names = FALSE)
  mod <- data.table::fread(md, data.table = FALSE, check.names = FALSE)

  rownames(eff) <- eff[[1]]; eff <- eff[, -1, drop = FALSE]
  # Column names look like "SYMBOL (12345)"
  colnames(eff) <- sub(" .*$", "", colnames(eff))
  list(effect = eff, model = mod)
}

glioma_selectivity <- function(dm) {
  mod <- dm$model
  idc <- grep("ModelID|DepMap_ID", names(mod), value = TRUE)[1]
  linc <- grep("OncotreeLineage|lineage|primary_disease", names(mod),
               value = TRUE)[1]
  if (is.na(idc) || is.na(linc)) {
    log_msg("Could not identify lineage columns in Model.csv")
    return(NULL)
  }
  lin <- setNames(as.character(mod[[linc]]), mod[[idc]])
  cl <- rownames(dm$effect)
  is_cns <- grepl("CNS|Brain|Glio", lin[cl], ignore.case = TRUE)
  is_cns[is.na(is_cns)] <- FALSE
  log_msg("CNS/brain cell lines: ", sum(is_cns), " of ", length(cl))
  if (sum(is_cns) < 10) return(NULL)

  E <- as.matrix(dm$effect)
  storage.mode(E) <- "numeric"
  mean_cns <- colMeans(E[is_cns, , drop = FALSE], na.rm = TRUE)
  mean_oth <- colMeans(E[!is_cns, , drop = FALSE], na.rm = TRUE)

  tibble(gene = colnames(E),
         effect_cns = mean_cns, effect_other = mean_oth,
         # more negative in CNS than elsewhere means selective dependency
         selectivity = mean_oth - mean_cns,
         pan_essential = mean_oth < -0.5) %>%
    filter(is.finite(selectivity))
}

main_23 <- function() {
  v <- dos_viability()
  gd <- v$gd

  # --- DOS ranking, if viable -------------------------------------------
  rna <- readRDS(P("rds", "deg_rnaseq.rds"))
  base <- rna$LGG_vs_Normal
  dos <- if (!is.null(base)) compute_dos(base, gd) else NULL
  if (!is.null(dos)) {
    write_csv(dos, P("tables", "dos_scores.csv"))
    log_msg("Top 10 by DOS: ",
            paste(head(dos$gene, 10), collapse = ", "))
    log_msg("Top 10 by raw |logFC|: ",
            paste(head(dos %>% arrange(desc(effect)) %>% pull(gene), 10),
                  collapse = ", "))
    ov <- length(intersect(head(dos$gene, 50),
                           head(dos %>% arrange(desc(effect)) %>%
                                  pull(gene), 50)))
    log_msg("Overlap of top 50 by DOS and by raw effect: ", ov, "/50",
            ifelse(ov > 45,
                   ". Nearly identical, confirming the metric is close to a no-op.",
                   ". The ranking does change."))
  }

  # --- DepMap validation -------------------------------------------------
  dm <- load_depmap()
  if (is.null(dm)) {
    log_msg("38 complete: test 1 only.")
    return(invisible(v))
  }
  sel <- glioma_selectivity(dm)
  if (is.null(sel)) { log_msg("Selectivity could not be computed."); return(invisible(v)) }
  write_csv(sel, P("tables", "depmap_glioma_selectivity.csv"))

  net <- readRDS(P("rds", "network.rds"))
  inter <- readRDS(P("rds", "intersect.rds"))
  robust <- inter$shared; if (is.data.frame(robust)) robust <- robust$gene
  pb_f <- P("rds", "pipeline_benchmark.rds")
  modules <- if (file.exists(pb_f)) readRDS(pb_f)$modules else list()

  sets <- list(`cytoHubba hubs` = net$consensus$gene,
               `composition-robust set` = robust)
  if (length(modules))
    sets[["proposed modules"]] <- unlist(modules, use.names = FALSE)
  if (!is.null(dos)) sets[["DOS top 50"]] <- head(dos$gene, 50)

  # Is the set enriched for glioma-selective dependencies?
  bench <- imap_dfr(sets, function(g, nm) {
    s <- sel %>% mutate(in_set = gene %in% g)
    if (sum(s$in_set) < 3) return(tibble())
    wt <- wilcox.test(selectivity ~ in_set, data = s)
    auc <- as.numeric(pROC::auc(pROC::roc(s$in_set, s$selectivity,
                                          quiet = TRUE)))
    tibble(gene_set = nm, n_in_depmap = sum(s$in_set),
           median_selectivity_in = median(s$selectivity[s$in_set]),
           median_selectivity_out = median(s$selectivity[!s$in_set]),
           auc_selectivity = auc,
           p = wt$p.value,
           pct_pan_essential = 100 * mean(s$pan_essential[s$in_set]))
  })
  write_csv(bench, P("tables", "depmap_benchmark.csv"))

  log_msg("=================================================")
  log_msg("TEST 2: DEPMAP CRISPR VALIDATION")
  log_msg("  auc_selectivity: can the gene set be separated from the rest")
  log_msg("  by glioma-selective dependency? 0.5 means no.")
  log_msg("  pct_pan_essential guards against rewarding generic essentials.")
  log_msg("=================================================")
  print(as.data.frame(bench %>% mutate(across(where(is.numeric),
                                              ~ signif(.x, 3)))))

  ch <- bench %>% filter(gene_set == "cytoHubba hubs")
  if (nrow(ch) && ch$pct_pan_essential > 50)
    log_msg("Note: ", round(ch$pct_pan_essential), " percent of cytoHubba ",
            "hubs are pan-essential. They are essential everywhere, not ",
            "selectively in glioma, which is what a degree-biased method ",
            "would nominate.")

  best <- bench[which.max(bench$auc_selectivity), ]
  log_msg("Best by glioma-selective dependency: ", best$gene_set,
          " at AUC ", round(best$auc_selectivity, 3))

  saveRDS(list(viability = v, benchmark = bench, selectivity = sel),
          P("rds", "dos_depmap.rds"))
  log_msg("23_depmap_validation complete.")
  invisible(bench)
}
