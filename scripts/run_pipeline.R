#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Run POINT-seq pre-mRNA splicing analysis pipeline.")
parser$add_argument("--config", required = TRUE, help = "Path to experiment YAML config.")
args <- parser$parse_args()

source("R/config.R")
source("R/utils.R")
source("R/samples.R")
source("R/annotation.R")
source("R/counting.R")
source("R/differential.R")
source("R/replicate_comparisons.R")
source("R/browser_exports.R")
source("R/reporting.R")

cfg <- read_experiment_config(args$config)
sample_table <- build_sample_table(cfg)
cfg <- validate_input_files(sample_table, cfg)

message("Loaded ", nrow(sample_table), " samples across ", length(unique(sample_table$condition)), " conditions.")
message("Samples: ", paste(sample_table$sample_id, collapse = ", "))

annotation <- load_hg38_annotation(cfg)
gene_features <- make_gene_features(annotation)
intron_features <- make_intron_features(annotation)
exon_features <- make_exon_bin_features(annotation)

message("Counting gene, intron, and exon-bin features.")
gene_counts <- count_features(gene_features, sample_table, cfg, "gene")
intron_counts <- count_features(intron_features, sample_table, cfg, "intron")
exon_counts <- count_features(exon_features, sample_table, cfg, "exon_bin")
junction_counts <- count_splice_junctions(sample_table, cfg)

message("Writing merged-condition and replicate-pairwise count summaries.")
write_merged_condition_counts(gene_counts, feature_metadata(gene_features), sample_table, cfg, "gene")
write_merged_condition_counts(intron_counts, feature_metadata(intron_features), sample_table, cfg, "intron")
write_merged_condition_counts(exon_counts, feature_metadata(exon_features), sample_table, cfg, "exon_bin")
write_replicate_pairwise_counts(gene_counts, feature_metadata(gene_features), sample_table, cfg, "gene")
write_replicate_pairwise_counts(intron_counts, feature_metadata(intron_features), sample_table, cfg, "intron")
write_replicate_pairwise_counts(exon_counts, feature_metadata(exon_features), sample_table, cfg, "exon_bin")

message("Running differential analyses.")
gene_result <- run_differential_nascent_transcription(gene_counts, gene_features, sample_table, cfg)
normalized_intron <- run_intron_coverage_normalization(intron_counts, intron_features, gene_result, sample_table, cfg)
intron_retention <- run_differential_intron_retention(intron_counts, intron_features, gene_result, sample_table, cfg)
alternative_splicing <- run_alternative_splicing(exon_counts, exon_features, gene_result, sample_table, cfg)
splicing_efficiency <- run_splicing_efficiency(intron_counts, intron_features, gene_result, sample_table, cfg)

message("Exporting genome browser outputs.")
export_top_changed_regions(gene_result$table, cfg, "differential_nascent_transcription")
export_top_changed_regions(normalized_intron$table, cfg, "normalized_intron_coverage")
export_top_changed_regions(intron_retention$table, cfg, "differential_intron_retention")
export_top_changed_regions(alternative_splicing$table, cfg, "alternative_splicing_exon_usage")
export_top_changed_regions(splicing_efficiency$table, cfg, "differential_splicing_efficiency")

results <- list(
  gene = gene_result,
  normalized_intron = normalized_intron,
  intron_retention = intron_retention,
  alternative_splicing = alternative_splicing,
  splicing_efficiency = splicing_efficiency,
  junctions = junction_counts
)
saveRDS(results, file.path(project_output_dir(cfg), "pipeline_results.rds"))
export_bigwig_tracks(results, cfg)

message("Rendering PDF report.")
render_pdf_report(cfg, sample_table, results)

message("Pipeline complete. Results written to: ", project_output_dir(cfg))
