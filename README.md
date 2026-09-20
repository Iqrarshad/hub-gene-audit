# hub-gene-audit

Code and analysis for a study of how hub genes are selected in cancer
transcriptomics. The short version of the finding: the genes returned by the
usual differential-expression to STRING to cytoHubba workflow are largely
predictable from how much each gene has already been studied, and none of
the corrections we tried removes that dependence without also removing the
biological signal. We also show that the cross-cohort reproducibility often
cited as evidence these hubs are meaningful is an artefact of the network
being fixed and data-independent.

The work began as a reanalysis of a specific glioma manuscript and grew into
a general audit of the workflow. The glioma data are the test bed; the
argument is about the method.

## What is here

The repository is a single R pipeline, run stage by stage. Each stage reads
cached inputs, does one thing, and writes named CSV and RDS outputs. Figures
are produced by separate stages that read those outputs and never recompute,
so the analysis and the plotting can be changed independently.

```
R/            numbered analysis stages plus helpers (00, 00b, 00c, 00e, 00f, 00g, 00h)
run_all.R     the stage runner
docs/         CHANGELOG with the full analysis history
published_hub_lists.csv   the one external supporting file
config.local.example.R    template for machine-specific paths
```

Stages are numbered sequentially, 01 to 44, in dependency order. Stages 42 to 44 are additive post-processing: they read the outputs of
01 to 41 and write new tables (bootstrap intervals, multi-seed
robustness, and a consolidated corrections table). They change no
earlier output.
`docs/CHANGELOG.md` records the analysis history, including the earlier
exploratory stages that were
replaced, so results from a superseded version are not cited by mistake.

## Requirements

R 4.3 or later. The pipeline depends on the usual Bioconductor and CRAN
stack for this kind of analysis (limma, edgeR, org.Hs.eg.db, AnnotationDbi,
igraph, RSpectra, Matrix, data.table, and the tidyverse). Install them with:

```r
source("install_dependencies.R")
```

This takes roughly half an hour on a clean machine.

## Data

The pipeline expects the source files under a directory you nominate. None
of them are redistributed here; they come from public repositories and, in
one case, a controlled-access portal. The datasets are listed in
`docs/CHANGELOG.md` under Inputs, with their accessions. In brief: five
discovery expression datasets (GSE147352, GSE16011, GSE108474, GSE15824,
GSE21354), four validation cohorts (TCGA-LGG, TCGA-GBM, CGGA-325, CGGA-693),
one single-cell reference (GSE84465), the STRING v12 human network, NCBI
gene2pubmed, and DepMap CRISPR essentiality.

## Paths

Nothing is hardcoded. The pipeline needs two directories, and they must be
different:

| Variable      | Contents                              | Access     |
|---------------|---------------------------------------|------------|
| `DATA_DIR`    | downloaded source files               | read only  |
| `RESULTS_DIR` | tables, figures, cached objects, logs | written to |

Set them in whichever way suits you; the first match wins.

```bash
# on the command line
Rscript run_all.R --data-dir /path/to/data --results-dir /path/to/out

# or in the environment
export GLIOMA_DATA_DIR=/path/to/data
export GLIOMA_RESULTS_DIR=/path/to/out
```

```r
# or in a local file, which git ignores
file.copy("config.local.example.R", "config.local.R")   # then edit it
```

With none of these set the pipeline falls back to `data/` and `results/`
beside the project root. Whatever it resolves to is printed on the first
line of every run, so a misconfigured path is obvious immediately.

## Running

A single stage, or several:

```bash
Rscript run_all.R 32
Rscript run_all.R 30 31 32
```

A named group:

```bash
Rscript run_all.R --group figures
```

Everything, in dependency order:

```bash
Rscript run_all.R
```

To discard all derived output and regenerate from scratch, leaving the
downloaded source data alone:

```bash
Rscript run_all.R --clean
```

`--list` prints the full stage table with groups.

Stages cache to `RESULTS_DIR/rds`, so a failed stage can be fixed and rerun
without redoing the earlier ones. A cold run takes a few hours, dominated by
the STRING and GEO queries in stages 12, 24 and 27; once cached, later runs
are quick.

A note on how the stages have been exercised: development was done by
sourcing individual stage files and calling their `main_NN()` functions
directly, which is how every reported number was produced. The end-to-end
`run_all.R` path is provided for reproduction and should be checked on your
own machine before you rely on a single full run.

## Stage groups

| Group       | Stages | Purpose                                              |
|-------------|--------|------------------------------------------------------|
| `core`      | 01-06  | data, DEG, enrichment, STRING network                |
| `validate`  | 07-10  | hub null models and centrality diagnostics           |
| `bias`      | 11-19, 36, 39 | hub-bias quantification and network diagnostics |
| `methods`   | 20-23  | attempted corrections, all reported as failures      |
| `lit`       | 24-28, 37 | literature-scale analysis and the propagation gates |
| `propagate` | 30-35  | propagation audit, information selection, baselines  |
| `figures`   | 29, 32, 38, 40 | manuscript figures and tables                |
| `recommend` | 41     | reproducibility-gated selection and audit card       |

## The main results, and where each comes from

Hub membership is predictable from a gene's total STRING degree at an AUC of
0.994; genes placed into random gene sets are still nominated as hubs
(stage 12). Publication count from gene2pubmed explains 54.5% of the variance
in interactome degree, or 28.4% with the text-mining channel removed
(stage 24). Eight attempted corrections, including
degree normalisation, degree-preserving nulls, cohort co-expression,
citation penalties and a CRISPR co-essentiality network, each either keep
the bias or lose the signal (stages 13, 20, 21, 22, 25, 26).

Network propagation on a text-mining-free STRING graph recovers degree
rather than differential expression: held-out recovery of a masked signal is
0.07 against the true values and 0.71 against degree, because degree
anti-predicts strong differential expression at an AUC of 0.42 while the
neighbourhood mean predicts it at 0.65 (stage 30). A calibration curve of
signals with graded graph smoothness places the achievable recovery at 0.26
for the observed autocorrelation, so the shortfall is specific and
measurable, not a general failure of the data (stage 30). Constraining the
diffusion operator to be orthogonal to degree or citation does not beat
simply adjusting the output for degree afterwards (stage 31).

Reproducibility, the property most often cited as validation, is the sharp
point. STRING hub lists replicate across cohorts almost perfectly, but that
is because the STRING degree vector is fixed and identical for every cohort.
When the network is built from each cohort's own expression, the full degree
ranking correlates only weakly across cohorts (Spearman 0.14 to 0.26) and
the top-k hub lists share no genes at all (stage 35 and the full-vector
diagnostic). High hub reproducibility measures insensitivity to the data,
not agreement about biology.

An information-theoretic selector that never touches the network, minimum
redundancy maximum relevance on the expression data, is less reproducible
than the workflow it would replace, and its reproducibility falls further as
the redundancy penalty is strengthened (stage 33). The ladder from raw
variance down to mRMR shows reproducibility decreasing as selection becomes
more elaborate (stage 34).

## Reproducing the figures

The figures for the audit are produced by stages 48 and 52 and depend only
on the CSV outputs of the analysis stages. If those outputs are present:

```bash
Rscript run_all.R 48 52 56
```

Figures are written as 1200 DPI TIFF (LZW), PDF and PNG, in Times New Roman,
using the Okabe-Ito colourblind-safe palette.

## Licence and citation

See `LICENSE`. If you use this code or build on the analysis, please cite the
associated paper; the citation will be added here on acceptance.
