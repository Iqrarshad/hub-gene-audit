# Configuration: paths, thresholds, pre-declared criteria.
# No script hardcodes a path. See the resolution order below.

# --- Paths ---------------------------------------------------------------
# DATA_DIR holds the downloaded source files and is READ-ONLY to this
# pipeline. RESULTS_DIR receives everything generated. Nothing writes into
# DATA_DIR.
#
# Resolution order, first match wins:
#   1. --data-dir and --results-dir on the run_all.R command line, which
#      are passed through as environment variables
#   2. GLIOMA_DATA_DIR and GLIOMA_RESULTS_DIR in the environment
#   3. config.local.R in the project root, which is not tracked by git
#   4. data/ and results/ beside the project root
#
# Copy config.local.example.R to config.local.R to set machine-specific
# paths without touching a tracked file.

if (file.exists("config.local.R")) source("config.local.R")

.env_data <- Sys.getenv("GLIOMA_DATA_DIR", "")
.env_res  <- Sys.getenv("GLIOMA_RESULTS_DIR", "")
if (nzchar(.env_data)) DATA_DIR <- .env_data
if (nzchar(.env_res))  RESULTS_DIR <- .env_res
if (!exists("DATA_DIR"))    DATA_DIR    <- file.path(getwd(), "data")
if (!exists("RESULTS_DIR")) RESULTS_DIR <- file.path(getwd(), "results")

DATA_DIR    <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)
RESULTS_DIR <- normalizePath(RESULTS_DIR, winslash = "/", mustWork = FALSE)
CACHE_DIR   <- file.path(RESULTS_DIR, "cache")   # downloads land here

if (grepl("CHANGE_ME", DATA_DIR) || grepl("CHANGE_ME", RESULTS_DIR))
  stop("config.local.R still contains CHANGE_ME. Edit it with real paths, ",
       "or pass --data-dir and --results-dir instead.")
if (!dir.exists(DATA_DIR))
  warning("DATA_DIR does not exist: ", DATA_DIR,
          "\nSet it with --data-dir, GLIOMA_DATA_DIR, or config.local.R.",
          call. = FALSE)
if (normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE) ==
    normalizePath(RESULTS_DIR, winslash = "/", mustWork = FALSE))
  stop("DATA_DIR and RESULTS_DIR must differ. Results would overwrite the ",
       "source data.")

for (d in c(RESULTS_DIR, CACHE_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
for (s in c("tables", "figures", "rds", "logs")) {
  dir.create(file.path(RESULTS_DIR, s), showWarnings = FALSE,
             recursive = TRUE)
}

P <- function(...) file.path(RESULTS_DIR, ...)
D <- function(...) file.path(DATA_DIR, ...)
C <- function(...) file.path(CACHE_DIR, ...)

# --- Reproducibility -----------------------------------------------------
SEED <- 42
set.seed(SEED)
options(stringsAsFactors = FALSE, timeout = 7200)
# Slow or intermittent links: give libcurl a long connection timeout and
# allow low-speed transfers to continue rather than aborting.
options(HTTPUserAgent = "hub-gene-audit R pipeline")
Sys.setenv(
  R_DEFAULT_INTERNET_TIMEOUT = "7200",
  VROOM_CONNECTION_SIZE = "500000")
Sys.setenv(VROOM_CONNECTION_SIZE = 500000)

# --- Known local files ---------------------------------------------------
# Filenames as they appear in DATA_DIR. Extensions are resolved at runtime
# by find_local() because the directory listing does not show them.
LOCAL <- list(
  gse147352_counts   = "GSE147352/GSE147352_DESeq_normalized_counts.csv",
  gse16011_matrix    = "GSE16011_series_matrix",
  gpl8542_annot      = "GPL8542_entrez_symbol",
  gse108474_expr     = "GSE108474_REMBRANDT_GeneExpression",
  gse108474_clin     = "GSE108474_REMBRANDT_clinical.data",
  tcga_lgg_expr      = "TCGA_LGG_expression",
  lgg_clin           = "LGG_data_clinical_patient",
  lgg_mut            = "LGG_data_mutations",
  hgg_clin           = "HGG_data_clinical_patient",
  hgg_mut            = "HGG_data_mutations",
  hgg_expr           = "HGG_data_mrna_seq_v2_rsem_zscores_ref_diploid_samples",
  lgg_expr           = "LGG_data_mrna_seq_v2_rsem_zscores_ref_diploid_samples"
)

# Resolve a stem to an actual file regardless of extension.
find_local <- function(stem) {
  direct <- D(stem)
  if (file.exists(direct)) return(direct)
  dirn <- dirname(direct); base <- basename(direct)
  if (!dir.exists(dirn)) return(NA_character_)
  files <- list.files(dirn, full.names = FALSE)
  if (length(files) == 0) return(NA_character_)
  hit <- files[files == base]
  if (length(hit) == 0) hit <- files[startsWith(files, base)]
  if (length(hit) == 0) {
    lf <- tolower(files); lb <- tolower(base)
    hit <- files[lf == lb]
    if (length(hit) == 0) hit <- files[startsWith(lf, lb)]
  }
  if (length(hit) == 0) return(NA_character_)
  file.path(dirn, hit[order(nchar(hit))][1])
}

# --- Datasets ------------------------------------------------------------
# use_for controls which analyses each dataset feeds.
#   deg    : tumour vs normal differential expression
#   grade  : LGG vs HGG contrast
#   deconv : immune deconvolution
#   surv   : survival
DATASETS <- list(
  GSE147352 = list(type = "rnaseq", local = "gse147352_counts",
                   use_for = c("deg", "grade", "deconv")),
  GSE16011  = list(type = "microarray", local = "gse16011_matrix",
                   use_for = c("deg", "grade", "deconv")),
  GSE108474 = list(type = "microarray", local = "gse108474_expr",
                   use_for = c("deg")),
  GSE15824  = list(type = "microarray", local = NA, use_for = c("deg")),
  GSE21354  = list(type = "microarray", local = NA, use_for = c("deg"))
)

# Datasets not present locally are downloaded into CACHE_DIR by 01_download.R
DOWNLOAD_IF_MISSING <- c("GSE15824", "GSE21354")

# Set to FALSE to demote the two severely underpowered microarrays
# (GSE15824 = 7 vs 3; GSE21354 = 10 vs 4) to a supplementary sensitivity
# analysis rather than including them in the primary DEG intersection.
INCLUDE_SMALL_ARRAYS <- TRUE

# --- Statistical thresholds ---------------------------------------------
# PRE-DECLARED. Do not edit after the first run against real data.
THRESH <- list(
  deg_logfc       = 1.0,
  deg_adjp        = 0.05,
  deg_min_datasets = 3,
  enrich_qval     = 0.05,
  string_score    = 700,
  string_relax    = 400,
  network_lcc_min = 0.30,
  corr_qval       = 0.05,
  corr_boot_n     = 2000,
  corr_min_abs_r  = 0.30,
  surv_p          = 0.05
)

PASS_FAIL <- c(
  "HUB: hub genes reported only if the largest connected component of the",
  "  STRING graph at score >= 700 contains >= 30% of nodes. If it does not,",
  "  the threshold is relaxed to 400 and BOTH results are reported.",
  "CORR: a correlation is SUPPORTED only if BH q < 0.05 AND the bootstrap",
  "  95% CI excludes zero AND the sign is unchanged after adjusting for",
  "  ESTIMATE tumour purity AND |rho| >= 0.30. Failures are reported.",
  "ENRICH: BH q < 0.05 only. No nominal p-values in the results text.",
  "SURV: reported separately for IDH-mutant and IDH-wildtype. Never pooled.",
  "N: every per-grade sample size is generated from data, never typed.",
  "DEG: |log2FC| >= 1 and BH-adjusted p < 0.05, consistent across scripts."
)

# --- Hub genes from the prior analysis -----------------------------------
# Carried forward for continuity and comparison. 06_network.R re-derives
# hubs independently; disagreement between the two is reported.
HUB_GENES_PRIOR <- c("FYN", "SLC17A7", "SNCA", "SOX2", "VAMP2",
                     "CD44", "SNAP25", "EGFR", "SH3GL2")

# --- Deconvolution -------------------------------------------------------
DECONV_METHODS <- c("quantiseq", "epic", "mcp_counter", "xcell")
DECONV_PRIMARY <- "quantiseq"

# --- Figures -------------------------------------------------------------
FIG <- list(
  font = "Times New Roman", font_pdf = "Times", dpi = 1200,
  width = 7.0, height = 5.0,
  formats = c("tiff", "pdf"),
  # Okabe-Ito without yellow (#F0E442), which is illegible on white.
  okabe_ito = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#0072B2", "#D55E00", "#CC79A7")
)

# --- Logging -------------------------------------------------------------
log_msg <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  message(m)
  cat(m, "\n", file = P("logs", "run.log"), append = TRUE)
}

writeLines(PASS_FAIL, P("logs", "pre_declared_criteria.txt"))
log_msg("Config loaded. DATA=", DATA_DIR, " RESULTS=", RESULTS_DIR)
