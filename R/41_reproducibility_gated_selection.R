# Reproducibility-gated selection and the audit card.
#
# The audit shows that ranking genes by network position on a differential
# expression set is the failing step, and that corrections applied inside that
# framework either keep the attention bias or lose the signal. This stage
# tests a procedure that inverts the order of operations, and a card that lets
# any gene list be checked on its own data.
#
# The procedure gates first and ranks never. In each cohort that carries a
# grade contrast it computes a composition-adjusted high-grade against
# low-grade effect for every candidate gene, then keeps a gene only if the
# effect is present, of the same sign, and past a declared magnitude in every
# contrast-capable cohort. The network is not used to select. Study attention
# enters as a diagnostic, not as a ranking axis, because attention and effect
# size are entangled and penalising one removes the other.
#
# Cohorts carry a provenance label taken from the cohort name. Same-source
# cohorts agree for reasons of shared processing rather than shared biology,
# which is the pattern documented for the fixed network. The stage therefore
# reports whether the replicated set also replicates across at least two
# provenance groups, and marks the result provenance-limited when only one
# group carries a contrast.
#
# The audit card is three numbers with thresholds fixed here, before any
# result is seen: dependence of list membership on publication count
# (area under the curve, chance is 0.5), cross-cohort agreement of the
# per-cohort nominations (mean pairwise Jaccard), and the size of the
# replicated set. The card is computed for the conventional consensus hub list
# and for the gated list on the same data.
#
# A small replicated set is reported as a measured ceiling on reliable
# discovery in these cohorts, not as a failure of the procedure. All inputs
# are optional; a missing input is logged and the stage returns without error.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

source("R/00c_cohorts.R")
source("R/00e_compadj.R")
source("R/00f_citations.R")

# Second-provenance cohort. REMBRANDT (GSE108474) is a microarray series that
# carries grades II to IV in one matrix, from a different platform and
# consortium than CGGA, so it provides a genuine cross-provenance grade
# contrast. The loaders, ID harmonisation and grade coding are reused from the
# discovery stage so the labelling matches the rest of the study. The result
# is optional: any failure returns NULL and the stage falls back to the
# provenance-limited comparison.
load_rembrandt_cohort <- function() {
  ok <- tryCatch({
    source("R/01_download.R"); source("R/02_deg_microarray.R"); TRUE
  }, error = function(e) FALSE)
  if (!ok) { log_msg("REMBRANDT: helper sources unavailable; skipping."); return(NULL) }

  obj <- tryCatch(load_gse108474(), error = function(e) {
    log_msg("REMBRANDT: expression not loaded (", conditionMessage(e), ")")
    NULL })
  if (is.null(obj) || is.null(obj$clin)) {
    log_msg("REMBRANDT: no expression or clinical; skipping."); return(NULL)
  }

  expr <- obj$expr
  mat  <- as.matrix(expr[, -1, drop = FALSE])
  rownames(mat) <- as.character(expr[[1]])
  storage.mode(mat) <- "numeric"

  idx <- match(normalise_rembrandt_id(colnames(mat)),
               normalise_rembrandt_id(obj$clin$SUBJECT_ID))
  wg <- toupper(trimws(as.character(obj$clin$WHO_GRADE)))[idx]
  grade <- rep(NA_character_, length(wg))
  grade[wg %in% c("II", "III")] <- "LGG"
  grade[wg == "IV"] <- "HGG"
  if (sum(grade == "LGG", na.rm = TRUE) < MIN_PER_GRADE ||
      sum(grade == "HGG", na.rm = TRUE) < MIN_PER_GRADE) {
    log_msg("REMBRANDT: too few graded samples; skipping."); return(NULL)
  }

  if (grepl("_at$", rownames(mat)[1])) {
    ann_f <- find_local("GPL570_annot")
    if (is.na(ann_f)) {
      log_msg("REMBRANDT: GPL570 annotation absent, cannot map probes; skipping.")
      return(NULL)
    }
    ann <- data.table::fread(ann_f, data.table = FALSE, check.names = FALSE)
    s_c <- grep("symbol", names(ann), ignore.case = TRUE, value = TRUE)
    if (!length(s_c)) { log_msg("REMBRANDT: no symbol column; skipping."); return(NULL) }
    i <- match(rownames(mat), as.character(ann[[names(ann)[1]]]))
    mat <- collapse_to_symbol(mat, as.character(ann[[s_c[1]]])[i])
  }

  keep <- !is.na(grade)
  mat <- mat[, keep, drop = FALSE]
  meta <- tibble(sample = colnames(mat), grade = grade[keep])
  log_msg("REMBRANDT: ", ncol(mat), " graded samples (",
          sum(meta$grade == "LGG"), " LGG, ", sum(meta$grade == "HGG"),
          " HGG), ", nrow(mat), " symbols")
  list(matrix = mat, meta = meta)
}

# Third-provenance cohort. GSE16011 (Gravendeel) is a second microarray series
# on a different platform again (GPL8542), so it adds a further independent
# provenance group. The pipeline previously used it only for tumour against
# normal, so grade is read here from the phenotype table. Optional and
# fail-soft.
grade_from_pdata <- function(pd) {
  cols <- grep("characteristic|title|source|description|grade|who",
               names(pd), ignore.case = TRUE, value = TRUE)
  if (!length(cols)) return(rep(NA_character_, nrow(pd)))
  blob <- tolower(apply(pd[, cols, drop = FALSE], 1,
                        function(r) paste(r, collapse = " | ")))
  g <- rep(NA_character_, length(blob))
  lgg <- grepl("grade\\s*ii(i)?\\b|grade\\s*[23]\\b|who\\s*ii(i)?\\b", blob)
  hgg <- grepl("grade\\s*iv\\b|grade\\s*4\\b|who\\s*iv\\b|glioblastoma|\\bgbm\\b",
               blob)
  g[lgg & !hgg] <- "LGG"
  g[hgg] <- "HGG"
  g
}

load_gse16011_cohort <- function() {
  ok <- tryCatch({
    source("R/01_download.R"); source("R/02_deg_microarray.R"); TRUE
  }, error = function(e) FALSE)
  if (!ok) { log_msg("GSE16011: helper sources unavailable; skipping."); return(NULL) }

  obj <- tryCatch(load_gse16011(), error = function(e) {
    log_msg("GSE16011: not loaded (", conditionMessage(e), ")"); NULL })
  if (is.null(obj) || is.null(obj$eset)) {
    log_msg("GSE16011: no eset; skipping."); return(NULL)
  }

  mat <- tryCatch(eset_to_matrix(obj$eset, "GSE16011", obj$annot),
                  error = function(e) NULL)
  if (is.null(mat)) { log_msg("GSE16011: symbol mapping failed; skipping."); return(NULL) }

  grade <- grade_from_pdata(Biobase::pData(obj$eset))
  if (sum(grade == "LGG", na.rm = TRUE) < MIN_PER_GRADE ||
      sum(grade == "HGG", na.rm = TRUE) < MIN_PER_GRADE) {
    log_msg("GSE16011: too few graded samples (",
            sum(grade == "LGG", na.rm = TRUE), " LGG, ",
            sum(grade == "HGG", na.rm = TRUE), " HGG); skipping.")
    return(NULL)
  }
  keep <- !is.na(grade) & colnames(mat) %in% colnames(mat)
  mat <- mat[, keep, drop = FALSE]
  meta <- tibble(sample = colnames(mat), grade = grade[keep])
  log_msg("GSE16011: ", ncol(mat), " graded samples (",
          sum(meta$grade == "LGG"), " LGG, ", sum(meta$grade == "HGG"),
          " HGG), ", nrow(mat), " symbols")
  list(matrix = mat, meta = meta)
}

# Conventional STRING and cytoHubba nomination, run per cohort so its
# reproducibility is measured on the same footing as the gated procedure. For
# each cohort the STRING subnetwork induced by that cohort's grade-DEG genes is
# scored by the same six centrality measures used elsewhere, and hubs are genes
# selected by at least four of them.
conventional_per_cohort <- function(nom, string_edges, top_k = 10) {
  ok <- tryCatch({ source("R/17_coexpression_network.R"); TRUE },
                 error = function(e) FALSE)
  if (!ok || is.null(string_edges)) return(NULL)
  suppressPackageStartupMessages(library(igraph))
  hubs <- lapply(names(nom), function(nm) {
    g <- nom[[nm]]
    if (length(g) < 10) return(character(0))
    se <- string_edges[string_edges$from %in% g & string_edges$to %in% g, ,
                       drop = FALSE]
    gr <- igraph::graph_from_data_frame(se[, c("from", "to")],
                                        directed = FALSE, vertices = g)
    tryCatch(hub_select(gr, top_k = top_k)$hubs, error = function(e)
      tryCatch(hub_select(gr)$hubs, error = function(e2) character(0)))
  })
  names(hubs) <- names(nom)
  hubs
}

# ---- pre-declared criteria, fixed before results are seen --------------
MIN_PER_GRADE <- 20      # samples per grade for a usable contrast
Q_MIN         <- 0.05    # BH q for the per-cohort adjusted contrast
EFF_MIN       <- 0.25    # minimum |grade coefficient| on log2 expression
MIN_SET       <- 10      # replicated genes below which we report a ceiling
CARD_ATTN_MAX <- 0.60    # attention AUC at or below this passes the card
N_BOOT        <- 5000    # bootstrap resamples for BCa intervals

# Tie-corrected AUC as a Mann-Whitney statistic. Fast enough for the bootstrap
# loop, and matches the rank-based definition used elsewhere.
fast_auc <- function(y, x) {
  ok <- is.finite(x) & is.finite(y)
  y <- y[ok]; x <- x[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 < 1 || n0 < 1) return(NA_real_)
  r <- rank(x)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# BCa interval for the attention AUC of a list over a universe, by resampling
# genes with replacement. Acceleration is taken from the regression (empirical
# influence) estimate on the existing replicates, so no per-observation
# jackknife refit is needed. Falls back to the percentile interval if BCa fails.
bca_attention_ci <- function(list_genes, universe, papers, R = N_BOOT) {
  d <- data.frame(
    y = as.integer(universe %in% list_genes),
    x = log10(1 + papers[match(universe, names(papers))]))
  d <- d[is.finite(d$x), ]
  est <- fast_auc(d$y, d$x)
  if (is.na(est) || sum(d$y) < 3 || sum(d$y) >= nrow(d))
    return(c(est, NA_real_, NA_real_, NA_character_))
  if (!requireNamespace("boot", quietly = TRUE)) {
    stat <- replicate(R, {
      i <- sample.int(nrow(d), replace = TRUE)
      fast_auc(d$y[i], d$x[i]) })
    q <- quantile(stat, c(0.025, 0.975), na.rm = TRUE)
    return(c(est, q[[1]], q[[2]], "percentile"))
  }
  b <- boot::boot(d, function(dat, i) fast_auc(dat$y[i], dat$x[i]), R = R)
  ci <- tryCatch(boot::boot.ci(b, type = "bca", conf = 0.95),
                 error = function(e) NULL)
  if (is.null(ci) || is.null(ci$bca)) {
    q <- quantile(b$t, c(0.025, 0.975), na.rm = TRUE)
    return(c(est, q[[1]], q[[2]], "percentile"))
  }
  c(est, ci$bca[4], ci$bca[5], "bca")
}

# Apply the gate at given thresholds to the stored per-cohort contrasts.
# Returns retained gene sets for a chosen cohort subset.
apply_gate <- function(contrasts, beta_thr, q_thr, cohort_subset = NULL) {
  cs <- if (is.null(cohort_subset)) names(contrasts) else cohort_subset
  cs <- intersect(cs, names(contrasts))
  if (length(cs) < 1) return(character(0))
  tested <- Reduce(intersect, lapply(contrasts[cs], function(x) x$gene))
  keep <- vapply(tested, function(g) {
    rows <- lapply(contrasts[cs], function(x) x[x$gene == g, ])
    pass <- vapply(rows, function(r)
      isTRUE(!is.na(r$q) && r$q < q_thr && abs(r$beta) >= beta_thr), logical(1))
    signs <- vapply(rows, function(r) r$sign, numeric(1))
    all(pass) && length(unique(signs[!is.na(signs)])) == 1
  }, logical(1))
  tested[keep]
}

# Save a ggplot to TIFF (LZW), PDF, and PNG at the project figure standard.
save_fig41 <- function(p, name, w = 6.5, h = 4.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible())
  for (fmt in unique(c(FIG$formats, "png"))) {
    f <- P("figures", paste0(name, ".", fmt))
    tryCatch({
      if (fmt == "tiff")
        ggplot2::ggsave(f, p, width = w, height = h, dpi = FIG$dpi,
                        device = "tiff", compression = "lzw")
      else if (fmt == "pdf")
        ggplot2::ggsave(f, p, width = w, height = h, device = cairo_pdf)
      else
        ggplot2::ggsave(f, p, width = w, height = h, dpi = FIG$dpi,
                        device = "png")
    }, error = function(e)
      log_msg("  ", fmt, " failed for ", name, ": ", conditionMessage(e)))
  }
  log_msg("Figure written: ", name)
}

theme_ms41 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ggplot2::theme_classic(base_size = 12, base_family = FIG$font) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text = ggplot2::element_text(colour = "black"))
}

jaccard <- function(a, b) {
  if (!length(union(a, b))) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

attention_auc <- function(list_genes, universe, papers) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NA_real_)
  d <- tibble(gene = universe,
              y = as.integer(universe %in% list_genes),
              x = log10(1 + papers[match(universe, names(papers))])) %>%
    filter(!is.na(x))
  if (sum(d$y) < 3 || sum(d$y) >= nrow(d)) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(d$y, d$x, quiet = TRUE)))
}

# Composition-adjusted grade contrast in one cohort. Returns one row per gene
# with the grade coefficient, its BH q, and the sign, or NULL if the cohort
# does not carry both grades.
cohort_contrast <- function(mat, meta, genes, label,
                            use_panels = c("neuro", "glial", "immune")) {
  grade <- meta$grade[match(colnames(mat), meta$sample)]
  keep <- !is.na(grade)
  mat <- mat[, keep, drop = FALSE]; grade <- grade[keep]
  tab <- table(grade)
  if (!all(c("LGG", "HGG") %in% names(tab)) ||
      min(tab[c("LGG", "HGG")]) < MIN_PER_GRADE) {
    log_msg(label, ": no usable grade contrast (",
            paste(names(tab), tab, sep = "=", collapse = ", "), ")")
    return(NULL)
  }
  g <- intersect(genes, rownames(mat))
  hg <- as.integer(grade == "HGG")

  # RSEM counts are large and need log2; microarray intensities are already on
  # a log scale, so a second log2 would distort them. Detect and branch.
  already_log <- max(mat, na.rm = TRUE) < 50

  Z <- data.frame(hg = hg)
  if ("neuro" %in% use_panels)
    Z$neuro <- composition_score(mat, NEURONAL_PANEL, exclude = g,
                                 label = paste(label, "neuronal"))
  if ("glial" %in% use_panels)
    Z$glial <- composition_score(mat, GLIAL_PANEL, exclude = g,
                                 label = paste(label, "glial"))
  if ("immune" %in% use_panels)
    Z$immune <- composition_score(mat, IMMUNE_PANEL, exclude = g,
                                  label = paste(label, "immune"))
  keepz <- vapply(Z, function(v) all(is.finite(v)) && sd(v) > 0, logical(1))
  keepz[1] <- TRUE
  Z <- Z[, keepz, drop = FALSE]

  res <- map_dfr(g, function(gene) {
    y <- if (already_log) as.numeric(mat[gene, ]) else
      log2(as.numeric(mat[gene, ]) + 1)
    d <- cbind(y = y, Z)
    fit <- tryCatch(lm(y ~ ., data = d), error = function(e) NULL)
    if (is.null(fit) || !"hg" %in% rownames(summary(fit)$coefficients))
      return(tibble(gene = gene, beta = NA_real_, p = NA_real_))
    co <- summary(fit)$coefficients["hg", ]
    tibble(gene = gene, beta = co[["Estimate"]], p = co[["Pr(>|t|)"]])
  })
  res$q <- p.adjust(res$p, method = "BH")
  res$sign <- sign(res$beta)
  res$pass <- !is.na(res$q) & res$q < Q_MIN & abs(res$beta) >= EFF_MIN
  res$cohort <- label
  log_msg(label, ": ", sum(res$pass, na.rm = TRUE), " of ", nrow(res),
          " genes pass the adjusted contrast",
          if (length(use_panels) < 3)
            paste0(" [panels: ", paste(use_panels, collapse = "+"), "]") else "")
  res
}

provenance_of <- function(nm) sub("[-_].*$", "", nm)

main_41 <- function() {
  inter <- tryCatch(readRDS(P("rds", "intersect.rds")),
                    error = function(e) NULL)
  if (is.null(inter)) {
    log_msg("intersect.rds not found; run the intersection stage first. ",
            "Reproducibility-gated selection cannot proceed.")
    return(invisible(NULL))
  }
  genes <- inter$shared$gene

  cohorts <- get_validation_cohorts()
  rem <- load_rembrandt_cohort()
  if (!is.null(rem)) cohorts[["REMBRANDT"]] <- rem
  g16 <- load_gse16011_cohort()
  if (!is.null(g16)) cohorts[["GSE16011"]] <- g16

  contrasts <- list()
  for (nm in names(cohorts))
    contrasts[[nm]] <- cohort_contrast(cohorts[[nm]]$matrix,
                                       cohorts[[nm]]$meta, genes, nm)
  contrasts <- contrasts[!vapply(contrasts, is.null, logical(1))]

  if (length(contrasts) < 1) {
    log_msg("No cohort carries a usable grade contrast. The gate cannot be ",
            "applied on the available cohorts; stage stops here.")
    return(invisible(NULL))
  }
  prov <- provenance_of(names(contrasts))
  prov_testable <- length(unique(prov)) >= 2
  log_msg("Contrast-capable cohorts: ", paste(names(contrasts), collapse = ", "),
          " across provenance groups: ", paste(unique(prov), collapse = ", "))
  if (!prov_testable)
    log_msg("WARNING: only one provenance group carries a contrast, so ",
            "cross-provenance replication is not testable. The replicated set ",
            "below is within-provenance and its agreement is optimistic.")

  # per-cohort nomination sets, for the reproducibility axis of the card
  nom <- lapply(contrasts, function(x) x$gene[x$pass])
  names(nom) <- names(contrasts)

  # replication gate: pass, same sign, in every contrast-capable cohort
  tested <- Reduce(intersect, lapply(contrasts, function(x) x$gene))
  gate <- map_dfr(tested, function(gene) {
    rows <- lapply(contrasts, function(x) x[x$gene == gene, ])
    passes <- vapply(rows, function(r) isTRUE(r$pass), logical(1))
    signs  <- vapply(rows, function(r) r$sign, numeric(1))
    same_sign <- length(unique(signs[!is.na(signs)])) == 1
    n_pass <- sum(passes)
    prov_pass <- unique(prov[passes & same_sign])
    tibble(gene = gene, n_cohorts_pass = n_pass,
           consistent_sign = same_sign,
           replicated = n_pass == length(contrasts) && same_sign,
           cross_provenance = if (prov_testable)
             length(prov_pass) >= 2 else NA)
  })
  rg_list <- gate$gene[gate$replicated]
  cross_list <- if (prov_testable)
    gate$gene[gate$replicated & gate$cross_provenance] else character(0)

  write_csv(gate, P("tables", "rgselect_gate.csv"))
  write_csv(tibble(gene = rg_list), P("tables", "rgselect_genes.csv"))
  log_msg("Replicated set: ", length(rg_list), " genes",
          if (prov_testable) paste0("; cross-provenance: ",
                                    length(cross_list)) else "")

  # ---- audit card, conventional workflow versus the gated procedure -----
  pap <- load_gene2pubmed()
  papers <- if (is.null(pap)) NULL else setNames(pap$n_papers, pap$gene)
  net <- tryCatch(readRDS(P("rds", "network.rds")), error = function(e) NULL)
  conv_fixed <- if (is.null(net)) character(0) else net$consensus$gene

  # STRING edges, from the same cache the co-expression stage reads.
  string_edges <- NULL
  if (!is.null(net)) {
    score <- as.integer(net$chosen)
    sf <- file.path(CACHE_DIR, "string_density", paste0("t", score, ".rds"))
    string_edges <- tryCatch(
      if (file.exists(sf)) readRDS(sf)
      else if (file.exists(P("rds", paste0("string_", score, ".rds"))))
        readRDS(P("rds", paste0("string_", score, ".rds"))) else NULL,
      error = function(e) NULL)
  }

  # Conventional STRING and cytoHubba, run per cohort on each cohort's own
  # grade-DEG genes, so its reproducibility is data dependent and comparable.
  conv_hubs <- conventional_per_cohort(nom, string_edges)
  conv_repro <- NA_real_; conv_consensus <- character(0); conv_pc_attn <- NA_real_
  if (!is.null(conv_hubs) && length(conv_hubs) >= 2) {
    conv_repro <- mean(combn(names(conv_hubs), 2, function(p)
      jaccard(conv_hubs[[p[1]]], conv_hubs[[p[2]]])), na.rm = TRUE)
    tabh <- table(unlist(conv_hubs))
    conv_consensus <- names(tabh)[tabh >= ceiling(length(conv_hubs) / 2)]
    if (!is.null(papers))
      conv_pc_attn <- attention_auc(conv_consensus, genes, papers)
  }

  rg_repro <- if (length(nom) >= 2)
    mean(combn(names(nom), 2, function(p) jaccard(nom[[p[1]]], nom[[p[2]]])),
         na.rm = TRUE) else NA_real_
  rg_attn <- if (is.null(papers)) NA_real_ else attention_auc(rg_list, genes, papers)
  conv_fixed_attn <- if (is.null(papers) || !length(conv_fixed)) NA_real_ else
    attention_auc(conv_fixed, genes, papers)

  # BCa bootstrap intervals for the three attention AUCs (5000 resamples).
  ci_fixed <- if (is.null(papers) || !length(conv_fixed)) rep(NA, 4) else
    bca_attention_ci(conv_fixed, genes, papers)
  ci_pc    <- if (is.null(papers) || !length(conv_consensus)) rep(NA, 4) else
    bca_attention_ci(conv_consensus, genes, papers)
  ci_rg    <- if (is.null(papers)) rep(NA, 4) else
    bca_attention_ci(rg_list, genes, papers)

  card <- tibble(
    list = c("conventional consensus hubs (fixed network)",
             "conventional STRING + cytoHubba (per cohort)",
             "reproducibility-gated"),
    n_genes = c(length(conv_fixed), length(conv_consensus), length(rg_list)),
    attention_auc = c(conv_fixed_attn, conv_pc_attn, rg_attn),
    attention_ci_lo = as.numeric(c(ci_fixed[2], ci_pc[2], ci_rg[2])),
    attention_ci_hi = as.numeric(c(ci_fixed[3], ci_pc[3], ci_rg[3])),
    ci_method = c(ci_fixed[4], ci_pc[4], ci_rg[4]),
    cross_cohort_jaccard = c(NA_real_, conv_repro, rg_repro),
    reproducibility_note = c(
      "fixed network, identical across cohorts by construction",
      "per-cohort STRING subgraph on each cohort's grade-DEGs",
      "per-cohort nominations, data dependent"))
  write_csv(card, P("tables", "audit_card.csv"))

  log_msg("=================================================")
  log_msg("AUDIT CARD")
  log_msg("  attention_auc: publication count predicting list membership,")
  log_msg("  chance is 0.5. cross_cohort_jaccard: agreement of per-cohort")
  log_msg("  nominations. The conventional list is data-independent, so its")
  log_msg("  perfect agreement is not evidence about biology.")
  log_msg("=================================================")
  print(as.data.frame(card %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))))

  # ---- pre-declared verdict ---------------------------------------------
  crit <- tibble(
    criterion = c("replicated_set_size_ge_min", "gated_attention_auc_le_max",
                  "cross_provenance_nonempty"),
    threshold = c(paste0(">= ", MIN_SET), paste0("<= ", CARD_ATTN_MAX),
                  ">= 1 (if testable)"),
    observed = c(length(rg_list),
                 ifelse(is.na(rg_attn), NA, round(rg_attn, 3)),
                 if (prov_testable) length(cross_list) else NA),
    met = c(length(rg_list) >= MIN_SET,
            !is.na(rg_attn) && rg_attn <= CARD_ATTN_MAX,
            if (prov_testable) length(cross_list) >= 1 else NA))
  write_csv(crit, P("tables", "rgselect_criteria.csv"))

  size_ok <- length(rg_list) >= MIN_SET
  attn_ok <- !is.na(rg_attn) && rg_attn <= CARD_ATTN_MAX
  prov_ok <- if (prov_testable) length(cross_list) >= 1 else TRUE

  if (!size_ok) {
    log_msg("VERDICT: CEILING. The replicated set has ", length(rg_list),
            " genes, below the pre-declared floor of ", MIN_SET, ". This is a ",
            "measured limit on reliable cross-cohort discovery in these ",
            "cohorts, reported as such rather than as a failure of the ",
            "procedure.")
  } else if (!attn_ok) {
    log_msg("VERDICT: BIASED. The replicated set is non-trivial (",
            length(rg_list), " genes) but its attention AUC is ",
            round(rg_attn, 3), ", above ", CARD_ATTN_MAX, ". Gating on ",
            "replication did not remove the study-attention bias, and this ",
            "must be reported plainly.")
  } else if (!prov_ok) {
    log_msg("VERDICT: WITHIN-PROVENANCE ONLY. The set passes size and ",
            "attention criteria but does not replicate across provenance ",
            "groups, so its reproducibility may reflect shared processing.")
  } else {
    log_msg("VERDICT: PASS. The replicated set has ", length(rg_list),
            " genes, attention AUC ", round(rg_attn, 3), " at or below ",
            CARD_ATTN_MAX,
            if (prov_testable) paste0(", and replicates across provenance (",
                                      length(cross_list), " genes)") else
              ", provenance not testable here",
            ". The procedure clears its pre-declared bar on these data.")
  }

  # ---- sensitivity analyses ---------------------------------------------
  cgga <- names(contrasts)[prov == "CGGA"]
  attn_of <- function(set) {
    if (is.null(papers) || length(set) < 3) return(NA_real_)
    fast_auc(as.integer(genes %in% set),
             log10(1 + papers[match(genes, names(papers))]))
  }

  # (a) beta and q threshold sweep
  grid <- expand.grid(beta = c(0.20, 0.25, 0.30, 0.50),
                      q = c(0.01, 0.05, 0.10))
  sweep <- purrr::pmap_dfr(grid, function(beta, q) {
    within <- apply_gate(contrasts, beta, q, cohort_subset = cgga)
    cross  <- apply_gate(contrasts, beta, q)
    tibble(beta = beta, q = q, n_within_cgga = length(within),
           n_cross_source = length(cross),
           attn_auc_cross = round(attn_of(cross), 3))
  })
  write_csv(sweep, P("tables", "rgselect_sensitivity.csv"))
  log_msg("Threshold sweep written: rgselect_sensitivity.csv")

  # (b) leave-one-cohort-out
  loco <- purrr::map_dfr(names(contrasts), function(drop) {
    keep <- setdiff(names(contrasts), drop)
    set <- apply_gate(contrasts, EFF_MIN, Q_MIN, cohort_subset = keep)
    tibble(dropped = drop, n_retained = length(set),
           attn_auc = round(attn_of(set), 3),
           genes = paste(set, collapse = "; "))
  })
  write_csv(loco, P("tables", "rgselect_loco.csv"))
  log_msg("Leave-one-cohort-out written: rgselect_loco.csv")

  # (c) composition-panel sensitivity
  configs <- list("neuro+glial+immune" = c("neuro", "glial", "immune"),
                  "neuro+glial" = c("neuro", "glial"),
                  "neuro+immune" = c("neuro", "immune"),
                  "glial+immune" = c("glial", "immune"),
                  "none" = character(0))
  comp <- purrr::map_dfr(names(configs), function(cfg) {
    panels <- configs[[cfg]]
    ctr <- list()
    for (nm in names(cohorts)) {
      cc <- tryCatch(cohort_contrast(cohorts[[nm]]$matrix, cohorts[[nm]]$meta,
                                     genes, nm, use_panels = panels),
                     error = function(e) NULL)
      if (!is.null(cc)) ctr[[nm]] <- cc
    }
    if (length(ctr) < 2)
      return(tibble(config = cfg, n_cross_source = NA_integer_,
                    attn_auc = NA_real_, genes = NA_character_))
    set <- apply_gate(ctr, EFF_MIN, Q_MIN)
    tibble(config = cfg, n_cross_source = length(set),
           attn_auc = round(attn_of(set), 3),
           genes = paste(set, collapse = "; "))
  })
  write_csv(comp, P("tables", "rgselect_composition_sensitivity.csv"))
  log_msg("Composition sensitivity written: rgselect_composition_sensitivity.csv")

  # (d) k sweep for the per-cohort conventional reproducibility
  ksweep <- purrr::map_dfr(c(5, 10, 20), function(k) {
    ch <- conventional_per_cohort(nom, string_edges, top_k = k)
    if (is.null(ch) || length(ch) < 2)
      return(tibble(k = k, conv_repro = NA_real_, n_consensus = NA_integer_))
    jac <- mean(combn(names(ch), 2, function(p)
      jaccard(ch[[p[1]]], ch[[p[2]]])), na.rm = TRUE)
    tab <- table(unlist(ch)); cons <- names(tab)[tab >= ceiling(length(ch) / 2)]
    tibble(k = k, conv_repro = round(jac, 3), n_consensus = length(cons))
  })
  write_csv(ksweep, P("tables", "rgselect_ksweep.csv"))
  log_msg("k sweep written: rgselect_ksweep.csv")

  # ---- figures -----------------------------------------------------------
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    OK <- FIG$okabe_ito
    rd_tab <- function(name) {
      f <- P("tables", name)
      if (file.exists(f)) suppressMessages(readr::read_csv(f,
        show_col_types = FALSE)) else NULL
    }

    # Figure 1: cytoHubba metrics reduce to degree
    tryCatch({
      t1 <- rd_tab("cytohubba_metric_vs_degree.csv")
      if (!is.null(t1)) {
        t1 <- t1[order(-t1$abs_spearman), ]
        t1$metric <- factor(t1$metric, levels = rev(t1$metric))
        p <- ggplot(t1, aes(metric, abs_spearman)) +
          geom_col(width = 0.7, fill = OK[5]) +
          geom_hline(yintercept = 0.9, linetype = "dashed", colour = OK[6]) +
          coord_flip(ylim = c(0, 1)) +
          labs(title = "cytoHubba metrics rank by degree",
               x = NULL, y = "Absolute Spearman correlation with degree") +
          theme_ms41()
        save_fig41(p, "fig1_cytohubba_degree", 6.5, 4.2)
      }
    }, error = function(e) log_msg("fig1 failed: ", conditionMessage(e)))

    # Figure 2: degree tracks study attention
    tryCatch({
      d1 <- rd_tab("degree_vs_papers.csv")
      d2 <- rd_tab("degree_vs_papers_notextmining.csv")
      if (!is.null(d1)) {
        d1$channel <- "all channels"
        base <- d1
        if (!is.null(d2)) { d2$channel <- "text-mining removed"
          base <- rbind(d1[, c("ld", "lp", "channel")],
                        d2[, c("ld", "lp", "channel")]) }
        p <- ggplot(base, aes(lp, ld, colour = channel)) +
          geom_point(alpha = 0.12, size = 0.5) +
          geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
          scale_colour_manual(values = c("all channels" = OK[5],
                                         "text-mining removed" = OK[6])) +
          labs(title = "Interactome degree tracks publication count",
               x = "log10(1 + publications)", y = "log10(1 + degree)",
               colour = NULL) +
          theme_ms41() + theme(legend.position = c(0.75, 0.15))
        save_fig41(p, "fig2_degree_attention", 6, 4.6)
      }
    }, error = function(e) log_msg("fig2 failed: ", conditionMessage(e)))

    # Figure 3: reproducibility vanishes at the reported resolution
    tryCatch({
      t3 <- rd_tab("edge_vs_hub_agreement.csv")
      if (!is.null(t3)) {
        t3$level <- factor(t3$level, levels = t3$level[order(t3$agreement)])
        p <- ggplot(t3, aes(level, agreement)) +
          geom_col(width = 0.65, fill = OK[5]) +
          coord_flip(ylim = c(0, 1)) +
          labs(title = "Cross-cohort agreement falls at finer resolution",
               x = NULL, y = "Agreement between cohorts") +
          theme_ms41()
        save_fig41(p, "fig3_reproducibility", 6.5, 3.6)
      }
    }, error = function(e) log_msg("fig3 failed: ", conditionMessage(e)))

    # Figure 4: the trade-off (dependence versus reproducibility)
    tryCatch({
      sel <- rd_tab("Table3_selectors.csv")
      keep_lab <- c("Variance", "DE t-statistic", "Relevance (MI)", "mRMR")
      pts_sel <- NULL
      if (!is.null(sel)) {
        s <- sel[sel$Selector %in% keep_lab, ]
        if (nrow(s) > 0) pts_sel <- tibble(
          label = s$Selector, dependence = as.numeric(s$Citation_AUC),
          repro = as.numeric(s$Reproducibility),
          family = "network-free selector")
      }
      # single fixed reference, de-duplicated from the cytoHubba/STRING rows
      pts_fixed <- NULL
      if (!is.null(sel)) {
        fr <- sel[sel$Selector %in% c("cytoHubba", "STRING degree"), ]
        if (nrow(fr) > 0) pts_fixed <- tibble(
          label = "STRING + cytoHubba (fixed)",
          dependence = as.numeric(fr$Citation_AUC[1]),
          repro = as.numeric(fr$Reproducibility[1]),
          family = "network centrality")
      }
      extra <- tibble(
        label = c("cytoHubba (per cohort)", "reproducibility-gated"),
        dependence = c(conv_pc_attn, rg_attn),
        repro = c(conv_repro, rg_repro),
        family = c("network centrality", "reproducibility-gated"))
      pts <- bind_rows(pts_fixed, pts_sel, extra)
      pts <- pts[is.finite(pts$dependence) & is.finite(pts$repro), ]
      if (nrow(pts) > 0) {
        lab_layer <- if (requireNamespace("ggrepel", quietly = TRUE))
          ggrepel::geom_text_repel(aes(label = label), size = 2.9,
            family = FIG$font, box.padding = 0.5, point.padding = 0.3,
            min.segment.length = 0, seed = 42, max.overlaps = 20) else
          geom_text(aes(label = label), size = 2.6, vjust = -0.9,
                    family = FIG$font)
        p <- ggplot(pts, aes(dependence, repro)) +
          annotate("rect", xmin = 0, xmax = 0.6, ymin = 0.6, ymax = 1.02,
                   fill = OK[4], alpha = 0.15) +
          annotate("text", x = 0.04, y = 0.99, hjust = 0, size = 3,
                   family = FIG$font, fontface = "italic",
                   label = "reproducible and attention-neutral (empty)") +
          geom_point(aes(shape = family, colour = family), size = 3.4) +
          lab_layer +
          scale_colour_manual(values = c(
            "network centrality" = OK[6], "network-free selector" = OK[7],
            "reproducibility-gated" = OK[5])) +
          scale_shape_manual(values = c(
            "network centrality" = 17, "network-free selector" = 16,
            "reproducibility-gated" = 15)) +
          scale_x_continuous(limits = c(0, 1)) +
          scale_y_continuous(limits = c(-0.02, 1.05)) +
          labs(title = "No method is both reproducible and attention-neutral",
               x = "Dependence on study attention (AUC)",
               y = "Cross-cohort reproducibility", shape = NULL, colour = NULL) +
          theme_ms41() + theme(legend.position = "bottom")
        save_fig41(p, "fig4_tradeoff", 7, 5.6)
      }
    }, error = function(e) log_msg("fig4 failed: ", conditionMessage(e)))

    # Figure 5a: the diagnostic scorecard with bootstrap intervals
    tryCatch({
      pd <- card %>% filter(is.finite(attention_auc)) %>%
        mutate(short = c("fixed hubs", "per-cohort cytoHubba",
                         "gated")[seq_len(n())],
               grp = ifelse(short == "gated", "gated", "conventional"))
      p <- ggplot(pd, aes(reorder(short, attention_auc), attention_auc,
                          fill = grp)) +
        geom_col(width = 0.55) +
        geom_errorbar(aes(ymin = attention_ci_lo, ymax = attention_ci_hi),
                      width = 0.2, colour = "black") +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = OK[6]) +
        annotate("text", x = 0.62, y = 0.52, hjust = 0, size = 3,
                 family = FIG$font, label = "chance") +
        scale_fill_manual(values = c(conventional = OK[5], gated = OK[2]),
                          guide = "none") +
        coord_flip(ylim = c(0, 1)) +
        labs(title = "Attention dependence with 95% bootstrap intervals",
             x = NULL, y = "Attention AUC") +
        theme_ms41()
      save_fig41(p, "fig5a_scorecard", 6.5, 3.4)
    }, error = function(e) log_msg("fig5a failed: ", conditionMessage(e)))

    # Figure 5b: retained-set size across the threshold sweep
    tryCatch({
      p <- ggplot(sweep, aes(factor(q), factor(beta), fill = n_cross_source)) +
        geom_tile(colour = "white", linewidth = 1) +
        geom_text(aes(label = n_cross_source), size = 3.4,
                  family = FIG$font) +
        scale_fill_gradient(low = OK[3], high = OK[6], name = "genes") +
        labs(title = "Cross-source retained set across thresholds",
             x = "BH q threshold", y = "Effect-size threshold (log2)") +
        theme_ms41()
      save_fig41(p, "fig5b_sensitivity", 5.8, 4)
    }, error = function(e) log_msg("fig5b failed: ", conditionMessage(e)))

    # Supplementary: spike-in null recurrence
    tryCatch({
      sp <- rd_tab("hub_null_spikein.csv")
      if (!is.null(sp)) {
        sp <- sp[order(-sp$null_hub_frequency), ]
        sp$gene <- factor(sp$gene, levels = rev(sp$gene))
        p <- ggplot(sp, aes(gene, null_hub_frequency)) +
          geom_col(width = 0.7, fill = OK[6]) + coord_flip(ylim = c(0, 1)) +
          labs(title = "Consensus hubs recur in random gene sets",
               x = NULL, y = "Fraction of random backgrounds hub") +
          theme_ms41()
        save_fig41(p, "sfig_spikein_null", 6, 4.5)
      }
    }, error = function(e) log_msg("sfig spikein failed: ", conditionMessage(e)))

    # Supplementary: leave-one-cohort-out
    tryCatch({
      p <- ggplot(loco, aes(reorder(dropped, n_retained), n_retained)) +
        geom_col(width = 0.6, fill = OK[5]) + coord_flip() +
        labs(title = "Retained set with each cohort left out",
             x = "Cohort dropped", y = "Genes retained") +
        theme_ms41()
      save_fig41(p, "sfig_leave_one_cohort_out", 6, 3.4)
    }, error = function(e) log_msg("sfig loco failed: ", conditionMessage(e)))

    # Supplementary: composition-panel sensitivity
    tryCatch({
      cp <- comp[is.finite(comp$n_cross_source), ]
      if (nrow(cp) > 0) {
        p <- ggplot(cp, aes(reorder(config, n_cross_source), n_cross_source)) +
          geom_col(width = 0.6, fill = OK[5]) +
          geom_text(aes(label = ifelse(is.na(attn_auc), "",
                    paste0("AUC ", attn_auc))), hjust = -0.1, size = 2.8,
                    family = FIG$font) +
          coord_flip() +
          labs(title = "Retained set across composition adjustments",
               x = "Composition covariates", y = "Cross-source genes") +
          theme_ms41()
        save_fig41(p, "sfig_composition_sensitivity", 6.5, 4)
      }
    }, error = function(e) log_msg("sfig comp failed: ", conditionMessage(e)))

    # Supplementary: STRING confidence-threshold hub persistence
    tryCatch({
      dt <- rd_tab("density_threshold_sweep.csv")
      if (!is.null(dt)) {
        dt <- dt[dt$sweep == "confidence", ]
        p <- ggplot(dt, aes(reorder(condition, n_edges), n_hubs)) +
          geom_col(width = 0.6, fill = OK[5]) + coord_flip() +
          labs(title = "Hub count across STRING confidence thresholds",
               x = NULL, y = "Number of consensus hubs") +
          theme_ms41()
        save_fig41(p, "sfig_string_threshold", 6, 3.4)
      }
    }, error = function(e) log_msg("sfig threshold failed: ", conditionMessage(e)))
  }

  saveRDS(list(gate = gate, rg_list = rg_list, cross_list = cross_list,
               card = card, criteria = crit, sweep = sweep, loco = loco,
               composition = comp, ksweep = ksweep,
               provenance_testable = prov_testable),
          P("rds", "rgselect.rds"))
  log_msg("41_reproducibility_gated_selection complete.")
  invisible(card)
}
