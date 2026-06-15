#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Run POINT-seq pre-mRNA splicing analysis pipeline.")
parser$add_argument("--config", required = TRUE, help = "Path to experiment YAML config.")
args <- parser$parse_args()

source("R/config.R")
source("R/samples.R")
source("R/analysis_placeholders.R")

cfg <- read_experiment_config(args$config)
sample_table <- build_sample_table(cfg)

message("Loaded ", nrow(sample_table), " samples across ", length(unique(sample_table$condition)), " conditions.")
message("Samples: ", paste(sample_table$sample_id, collapse = ", "))
message("Pipeline scaffold is installed. Analysis modules will be implemented in upcoming phases.")
