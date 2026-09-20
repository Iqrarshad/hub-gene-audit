# Shared bootstrap helpers for the added interval, robustness, and
# corrections stages (42 to 44).
#
# Base R only. No dependency beyond what the pipeline already installs, so
# these stages run even if the tidyverse or pROC are unavailable. The AUC is
# the tie-corrected Mann-Whitney statistic, the same quantity stage 41 uses,
# and the interval is bias-corrected and accelerated (BCa) with a percentile
# fallback when the acceleration or bias term is not finite.
#
# The jackknife needed by BCa is computed in closed form from placement
# values, so a BCa interval is affordable even on the interactome-wide
# universe (about 19,000 genes) without an O(n^2) leave-one-out loop.

# Tie-corrected AUC = Pr(x_member > x_non-member), with half weight on ties.
fast_auc <- function(y, x) {
  y  <- as.integer(y)
  ok <- is.finite(x) & !is.na(y)
  x  <- x[ok]; y <- y[ok]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(x)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Counts of a reference vector below and equal to each value in v.
.placements <- function(v, ref) {
  s    <- sort(ref)
  n_le <- findInterval(v, s)                    # ref <= v
  n_lt <- findInterval(v, s, left.open = TRUE)  # ref <  v
  list(lt = n_lt, eq = n_le - n_lt)
}

# BCa interval for the AUC of predictor x separating the two classes in y.
# Resamples rows (genes) with replacement. Returns a named vector with the
# point estimate, the interval, the number of usable bootstrap replicates,
# and the method actually used ("bca" or "percentile").
bca_auc_ci <- function(y, x, R = 5000, conf = 0.95, seed = NULL) {
  y  <- as.integer(y)
  ok <- is.finite(x) & !is.na(y)
  x  <- x[ok]; y <- y[ok]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L); n <- n1 + n0
  if (n1 == 0L || n0 == 0L)
    return(c(estimate = NA_real_, lo = NA_real_, hi = NA_real_,
             R = 0, method = NA_real_))

  est <- fast_auc(y, x)
  if (!is.null(seed)) set.seed(seed)

  boot <- numeric(R)
  for (b in seq_len(R)) {
    i <- sample.int(n, n, replace = TRUE)
    boot[b] <- fast_auc(y[i], x[i])
  }
  boot <- boot[is.finite(boot)]
  a2   <- (1 - conf) / 2

  if (length(boot) < 50L) {
    out <- c(estimate = est, lo = NA_real_, hi = NA_real_, R = length(boot))
    attr(out, "method") <- "failed"
    return(out)
  }

  perc <- as.numeric(quantile(boot, c(a2, 1 - a2), names = FALSE,
                              na.rm = TRUE))
  lo <- perc[1]; hi <- perc[2]; meth <- "percentile"

  # Closed-form jackknife: positive placements among negatives and negative
  # placements among positives. Skipped when either class is a singleton.
  if (n1 >= 2L && n0 >= 2L) {
    pos <- x[y == 1L]; neg <- x[y == 0L]
    pp  <- .placements(pos, neg)
    p_i <- (pp$lt + 0.5 * pp$eq) / n0                 # positive beats negatives
    qq  <- .placements(neg, pos)
    s_j <- (n1 - qq$lt - 0.5 * qq$eq) / n1            # positives beat this negative
    jack <- c((n1 * est - p_i) / (n1 - 1),
              (n0 * est - s_j) / (n0 - 1))
    jbar <- mean(jack)
    den  <- 6 * (sum((jbar - jack)^2))^(3 / 2)
    acc  <- if (den == 0) NA_real_ else sum((jbar - jack)^3) / den
    prop <- mean(boot < est)
    z0   <- if (prop <= 0 || prop >= 1) NA_real_ else qnorm(prop)
    if (is.finite(z0) && is.finite(acc)) {
      zl <- qnorm(a2); zu <- qnorm(1 - a2)
      a1 <- pnorm(z0 + (z0 + zl) / (1 - acc * (z0 + zl)))
      a3 <- pnorm(z0 + (z0 + zu) / (1 - acc * (z0 + zu)))
      if (all(is.finite(c(a1, a3))) && a1 > 0 && a3 < 1 && a1 < a3) {
        g  <- as.numeric(quantile(boot, c(a1, a3), names = FALSE,
                                  na.rm = TRUE))
        lo <- g[1]; hi <- g[2]; meth <- "bca"
      }
    }
  }

  out <- c(estimate = est, lo = lo, hi = hi, R = length(boot))
  attr(out, "method") <- meth
  out
}

# --- Small base-R IO helpers, tolerant of missing files ------------------

# Read a table from RESULTS_DIR/tables, or NULL if absent.
read_table_safe <- function(file) {
  f <- P("tables", file)
  if (!file.exists(f)) return(NULL)
  tryCatch(utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
           error = function(e) NULL)
}

# One numeric cell from a keyed table: match key_val in key_col, return
# val_col as a number. NA when the file, columns, or row are absent.
pick_cell <- function(file, key_col, key_val, val_col) {
  d <- read_table_safe(file)
  if (is.null(d) || !all(c(key_col, val_col) %in% names(d)))
    return(NA_real_)
  i <- which(d[[key_col]] == key_val)
  if (!length(i)) return(NA_real_)
  suppressWarnings(as.numeric(d[[val_col]][i[1]]))
}

# One logical cell from a keyed table, returned as 1, 0, or NA. Handles the
# TRUE/FALSE strings that write.csv produces for logical columns.
pick_flag <- function(file, key_col, key_val, val_col) {
  d <- read_table_safe(file)
  if (is.null(d) || !all(c(key_col, val_col) %in% names(d)))
    return(NA_real_)
  i <- which(d[[key_col]] == key_val)
  if (!length(i)) return(NA_real_)
  v <- toupper(trimws(as.character(d[[val_col]][i[1]])))
  if (v %in% c("TRUE", "T", "1")) return(1)
  if (v %in% c("FALSE", "F", "0")) return(0)
  NA_real_
}

# Consensus hub genes, from the CSV if present, else the network RDS.
consensus_hub_genes <- function() {
  d <- read_table_safe("hub_genes_consensus.csv")
  if (!is.null(d)) {
    gcol <- intersect(c("gene", "Gene", "symbol", "SYMBOL"), names(d))
    if (length(gcol)) return(unique(as.character(d[[gcol[1]]])))
    return(unique(as.character(d[[1]])))
  }
  f <- P("rds", "network.rds")
  if (file.exists(f)) {
    net <- readRDS(f)
    if (!is.null(net$consensus$gene)) return(unique(as.character(net$consensus$gene)))
  }
  character(0)
}

# The 131 published hub genes from the shipped list.
published_hub_genes <- function() {
  src <- if (file.exists("published_hub_lists.csv"))
    "published_hub_lists.csv" else P("tables", "published_hub_lists.csv")
  if (!file.exists(src)) return(character(0))
  d <- utils::read.csv(src, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"hub_genes" %in% names(d)) return(character(0))
  g <- unlist(strsplit(paste(d$hub_genes, collapse = ";"), ";"))
  g <- toupper(trimws(g))
  g <- g[grepl("^[A-Z][A-Z0-9._-]{0,20}$", g)]
  unique(g)
}
