# install_dependencies.R --------------------------------------------------

cran <- c("dplyr","tidyr","purrr","readr","tibble","stringr","data.table",
          "ggplot2","scales","boot","ppcor","igraph","httr","survival",
          "survminer","UpSetR","remotes","BiocManager","renv",
          "patchwork","RSpectra","pROC","WGCNA","dynamicTreeCut","ggrepel")

for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p); install.packages(p, repos = "https://cloud.r-project.org")
  }
}

# GO.db and impute are pulled in by WGCNA; listed so a fresh machine gets
# them from Bioconductor rather than failing at first use.
bioc <- c("GEOquery","limma","oligo","edgeR","Biobase","clusterProfiler",
          "org.Hs.eg.db","enrichplot","preprocessCore","AnnotationDbi",
          "GO.db","impute")
BiocManager::install(setdiff(bioc, rownames(installed.packages())),
                     ask = FALSE, update = FALSE)

# ESTIMATE is not on CRAN
if (!requireNamespace("estimate", quietly = TRUE)) {
  install.packages("estimate", repos = "http://r-forge.r-project.org",
                   dependencies = TRUE)
}

message("\n--- Verification ---")
need <- c(cran, bioc, "estimate")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("STILL MISSING: ", paste(missing, collapse = ", "))
} else {
  message("All dependencies present.")
}
