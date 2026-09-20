# Cross-platform DEG intersection.
#
# Two intersections are produced: conventional (unadjusted) and
# composition-robust. The difference measures how much of a standard glioma
# DEG signature is tissue composition. Jaccard indices and direction
# consistency are reported.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr); library(tibble)
})

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

gather_sets <- function(which = c("robust", "raw")) {
  which <- match.arg(which)
  micro <- readRDS(P("rds", "deg_microarray.rds"))
  rna   <- readRDS(P("rds", "deg_rnaseq.rds"))

  sets <- list()
  for (nm in names(micro)) {
    d <- micro[[nm]]
    if (is.null(d)) next
    sel <- if (which == "robust" && any(!is.na(d$composition_robust)))
      d %>% filter(composition_robust) else d %>% filter(sig_raw)
    if (nrow(sel)) sets[[nm]] <- sel %>%
      transmute(gene, logFC = if (which == "robust" && !all(is.na(logFC_adj)))
                logFC_adj else logFC_raw)
  }
  d <- rna$LGG_vs_Normal
  if (!is.null(d)) {
    sel <- if (which == "robust") d %>% filter(composition_robust) else
      d %>% filter(sig_raw)
    if (nrow(sel)) sets[["GSE147352"]] <- sel %>%
      transmute(gene, logFC = if (which == "robust") logFC_adj else logFC_raw)
  }
  sets[vapply(sets, function(s) nrow(s) > 0, logical(1))]
}

jaccard_matrix <- function(sets) {
  ns <- names(sets)
  m <- matrix(NA_real_, length(ns), length(ns), dimnames = list(ns, ns))
  for (i in ns) for (j in ns) {
    a <- sets[[i]]$gene; b <- sets[[j]]$gene
    u <- length(union(a, b))
    m[i, j] <- if (u > 0) length(intersect(a, b)) / u else NA_real_
  }
  m
}

summarise_sets <- function(sets, label) {
  long <- imap_dfr(sets, function(s, nm)
    s %>% mutate(dataset = nm, direction = ifelse(logFC > 0, "up", "down")))

  summ <- long %>% group_by(gene) %>%
    summarise(n_datasets = n_distinct(dataset),
              n_up = sum(direction == "up"), n_down = sum(direction == "down"),
              datasets = paste(sort(unique(dataset)), collapse = ";"),
              mean_logFC = mean(logFC), .groups = "drop") %>%
    mutate(consistent = (n_up == 0) | (n_down == 0)) %>%
    arrange(desc(n_datasets), desc(abs(mean_logFC)))

  write_csv(summ, P("tables", paste0("deg_across_datasets_", label, ".csv")))

  shared <- summ %>% filter(n_datasets >= THRESH$deg_min_datasets, consistent)
  write_csv(shared, P("tables", paste0("shared_degs_", label, ".csv")))
  log_msg(label, ": ", nrow(shared), " genes in >= ",
          THRESH$deg_min_datasets, " datasets, direction-consistent")
  list(summary = summ, shared = shared)
}

main_04 <- function() {
  raw_sets <- gather_sets("raw")
  rob_sets <- gather_sets("robust")

  log_msg("Unadjusted DEG set sizes: ",
          paste(names(raw_sets), vapply(raw_sets, nrow, integer(1)),
                sep = "=", collapse = ", "))
  log_msg("Composition-robust set sizes: ",
          paste(names(rob_sets), vapply(rob_sets, nrow, integer(1)),
                sep = "=", collapse = ", "))

  jm <- jaccard_matrix(rob_sets)
  write_csv(as.data.frame(jm) %>% rownames_to_column("dataset"),
            P("tables", "jaccard_index.csv"))
  log_msg("Pairwise Jaccard indices (composition-robust):")
  print(round(jm, 3))

  raw <- summarise_sets(raw_sets, "conventional")
  rob <- summarise_sets(rob_sets, "composition_robust")

  lost <- setdiff(raw$shared$gene, rob$shared$gene)
  gained <- setdiff(rob$shared$gene, raw$shared$gene)

  cmp <- tibble(
    approach = c("conventional", "composition-robust"),
    n_shared = c(nrow(raw$shared), nrow(rob$shared)))
  write_csv(cmp, P("tables", "intersection_comparison.csv"))
  write_csv(tibble(gene = lost),
            P("tables", "degs_lost_to_composition_adjustment.csv"))

  log_msg("Conventional intersection: ", nrow(raw$shared), " genes")
  log_msg("Composition-robust intersection: ", nrow(rob$shared), " genes")
  log_msg(length(lost), " genes lost to composition adjustment, ",
          length(gained), " gained")

  hub_status <- tibble(
    gene = HUB_GENES_PRIOR,
    in_conventional = HUB_GENES_PRIOR %in% raw$shared$gene,
    in_robust = HUB_GENES_PRIOR %in% rob$shared$gene)
  write_csv(hub_status, P("tables", "hub_genes_intersection_status.csv"))
  log_msg("Original hub genes in each intersection:")
  print(as.data.frame(hub_status))

  saveRDS(list(sets = rob_sets, jaccard = jm, shared = rob$shared,
               conventional = raw$shared),
          P("rds", "intersect.rds"))
  log_msg("04_intersect complete.")
  invisible(rob$shared)
}
