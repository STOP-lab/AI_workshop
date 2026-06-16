make_coldata <- function(sample_table) {
  data.frame(
    sample_id = sample_table$sample_id,
    condition = factor(sample_table$condition),
    condition_key = factor(sample_table$condition_key, levels = c("condition_a", "condition_b")),
    replicate_index = sample_table$replicate_index,
    row.names = sample_table$sample_id,
    stringsAsFactors = FALSE
  )
}

run_deseq2 <- function(counts, sample_table, metadata, cfg, label, normalization_factors = NULL) {
  require_package("DESeq2")
  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "tables", label))
  coldata <- make_coldata(sample_table)
  counts <- round(counts[, rownames(coldata), drop = FALSE])

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition_key
  )
  keep <- rowSums(DESeq2::counts(dds)) >= (cfg$analysis$min_feature_count %||% 10)
  dds <- dds[keep, ]
  if (nrow(dds) == 0) {
    stop("No features passed min_feature_count for analysis: ", label, call. = FALSE)
  }

  if (!is.null(normalization_factors)) {
    nf <- normalization_factors[rownames(dds), rownames(coldata), drop = FALSE]
    nf <- pmax(as.matrix(nf), 1e-6)
    nf <- nf / exp(rowMeans(log(nf)))
    DESeq2::normalizationFactors(dds) <- nf
  }

  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  result <- as.data.frame(DESeq2::results(dds, contrast = c("condition_key", "condition_b", "condition_a")))
  result$feature_id <- rownames(result)
  result <- merge(metadata, result, by = "feature_id", all.y = TRUE, sort = FALSE)
  result <- result[padj_or_pvalue_order(result), , drop = FALSE]
  write_tsv(result, file.path(out_dir, paste0(label, "_differential.tsv")))

  normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  normalized_counts$feature_id <- rownames(normalized_counts)
  write_tsv(normalized_counts, file.path(out_dir, paste0(label, "_normalized_counts.tsv")))

  list(dds = dds, table = result, normalized_counts = normalized_counts)
}

run_differential_nascent_transcription <- function(gene_counts, gene_features, sample_table, cfg) {
  metadata <- feature_metadata(gene_features)
  run_deseq2(gene_counts, sample_table, metadata, cfg, "differential_nascent_transcription")
}

build_gene_normalization_factors <- function(feature_metadata, gene_norm_counts, sample_table) {
  gene_matrix <- as.matrix(gene_norm_counts[, sample_table$sample_id, drop = FALSE])
  rownames(gene_matrix) <- gene_norm_counts$feature_id
  factors <- matrix(1, nrow = nrow(feature_metadata), ncol = nrow(sample_table))
  rownames(factors) <- feature_metadata$feature_id
  colnames(factors) <- sample_table$sample_id

  matched <- feature_metadata$gene_id %in% rownames(gene_matrix)
  factors[matched, ] <- gene_matrix[feature_metadata$gene_id[matched], , drop = FALSE] + 1
  factors
}

run_intron_coverage_normalization <- function(intron_counts, intron_features, gene_result, sample_table, cfg) {
  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "tables", "normalized_intron_coverage"))
  metadata <- feature_metadata(intron_features)
  gene_norm <- gene_result$normalized_counts
  gene_matrix <- as.matrix(gene_norm[, sample_table$sample_id, drop = FALSE])
  rownames(gene_matrix) <- gene_norm$feature_id

  intron_library_sizes <- colSums(intron_counts)
  intron_size_factors <- intron_library_sizes / median(intron_library_sizes[intron_library_sizes > 0])
  intron_size_factors[!is.finite(intron_size_factors) | intron_size_factors <= 0] <- 1
  intron_norm <- sweep(intron_counts, 2, intron_size_factors, "/")
  normalized <- metadata
  for (sample_id in sample_table$sample_id) {
    gene_values <- gene_matrix[metadata$gene_id, sample_id]
    normalized[[paste0(sample_id, "_intron_count")]] <- intron_counts[metadata$feature_id, sample_id]
    normalized[[paste0(sample_id, "_intron_per_gene")]] <- intron_norm[metadata$feature_id, sample_id] / (gene_values + 1)
  }

  conditions <- unique(sample_table$condition)
  for (condition in conditions) {
    samples <- sample_table$sample_id[sample_table$condition == condition]
    cols <- paste0(samples, "_intron_per_gene")
    normalized[[paste0(condition, "_mean_intron_per_gene")]] <- rowMeans(normalized[, cols, drop = FALSE], na.rm = TRUE)
  }

  if (length(conditions) == 2) {
    normalized$log2_fold_change <- safe_log2_ratio(
      normalized[[paste0(conditions[[2]], "_mean_intron_per_gene")]],
      normalized[[paste0(conditions[[1]], "_mean_intron_per_gene")]]
    )
  }

  write_tsv(normalized, file.path(out_dir, "normalized_intron_coverage.tsv"))
  list(table = normalized, normalized_matrix = intron_norm)
}

run_differential_intron_retention <- function(intron_counts, intron_features, gene_result, sample_table, cfg) {
  metadata <- feature_metadata(intron_features)
  factors <- build_gene_normalization_factors(metadata, gene_result$normalized_counts, sample_table)
  run_deseq2(intron_counts, sample_table, metadata, cfg, "differential_intron_retention", normalization_factors = factors)
}

run_alternative_splicing <- function(exon_counts, exon_features, gene_result, sample_table, cfg) {
  metadata <- feature_metadata(exon_features)
  factors <- build_gene_normalization_factors(metadata, gene_result$normalized_counts, sample_table)
  run_deseq2(exon_counts, sample_table, metadata, cfg, "alternative_splicing_exon_usage", normalization_factors = factors)
}

run_splicing_efficiency <- function(intron_counts, intron_features, gene_result, sample_table, cfg) {
  metadata <- feature_metadata(intron_features)
  factors <- build_gene_normalization_factors(metadata, gene_result$normalized_counts, sample_table)
  result <- run_deseq2(intron_counts, sample_table, metadata, cfg, "differential_splicing_efficiency", normalization_factors = factors)

  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "tables", "splicing_efficiency_metrics"))
  gene_norm <- gene_result$normalized_counts
  gene_matrix <- as.matrix(gene_norm[, sample_table$sample_id, drop = FALSE])
  rownames(gene_matrix) <- gene_norm$feature_id

  metrics <- metadata
  for (sample_id in sample_table$sample_id) {
    unspliced <- intron_counts[metadata$feature_id, sample_id]
    gene_signal <- gene_matrix[metadata$gene_id, sample_id]
    metrics[[paste0(sample_id, "_unspliced")]] <- unspliced
    metrics[[paste0(sample_id, "_spliced_proxy")]] <- pmax(gene_signal - unspliced, 0)
    metrics[[paste0(sample_id, "_spliced_unspliced_ratio")]] <- (metrics[[paste0(sample_id, "_spliced_proxy")]] + 1) / (unspliced + 1)
    metrics[[paste0(sample_id, "_splicing_index")]] <- log2(metrics[[paste0(sample_id, "_spliced_unspliced_ratio")]])
  }
  write_tsv(metrics, file.path(out_dir, "splicing_efficiency_metrics.tsv"))
  result$metrics <- metrics
  result
}
