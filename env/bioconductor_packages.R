# Initial package manifest for the POINT-seq splicing pipeline.
# Run interactively or from Rscript to install missing dependencies.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "argparse",
  "yaml",
  "readxl",
  "data.table",
  "dplyr",
  "ggplot2",
  "knitr",
  "rmarkdown",
  "testthat"
)

bioc_packages <- c(
  "BiocParallel",
  "DESeq2",
  "DEXSeq",
  "edgeR",
  "GenomicAlignments",
  "GenomicFeatures",
  "GenomicRanges",
  "GenomeInfoDb",
  "IRanges",
  "Rsamtools",
  "rtracklayer",
  "trackViewer",
  "SummarizedExperiment"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    return(invisible(TRUE))
  }
  if (bioc) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing)
  }
}

install_if_missing(cran_packages, bioc = FALSE)
install_if_missing(bioc_packages, bioc = TRUE)
