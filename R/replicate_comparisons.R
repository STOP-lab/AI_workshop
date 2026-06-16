condition_sample_sets <- function(sample_table) {
  split(sample_table$sample_id, sample_table$condition)
}

write_merged_condition_counts <- function(counts, metadata, sample_table, cfg, label) {
  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "tables", "merged_condition_counts"))
  sets <- condition_sample_sets(sample_table)
  merged <- metadata
  for (condition in names(sets)) {
    merged[[condition]] <- rowSums(counts[metadata$feature_id, sets[[condition]], drop = FALSE])
  }
  if (length(sets) == 2) {
    merged$log2_fold_change <- safe_log2_ratio(merged[[names(sets)[[2]]]], merged[[names(sets)[[1]]]])
  }
  write_tsv(merged, file.path(out_dir, paste0(label, "_merged_condition_counts.tsv")))
  merged
}

write_replicate_pairwise_counts <- function(counts, metadata, sample_table, cfg, label) {
  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "tables", "replicate_pairwise"))
  condition_keys <- sort(unique(sample_table$condition_key))
  if (length(condition_keys) != 2) {
    warning("Skipping replicate pairwise output for ", label, "; expected exactly two condition keys.", call. = FALSE)
    return(invisible(NULL))
  }

  a <- sample_table[sample_table$condition_key == condition_keys[[1]], , drop = FALSE]
  b <- sample_table[sample_table$condition_key == condition_keys[[2]], , drop = FALSE]
  pair_indexes <- intersect(a$replicate_index, b$replicate_index)
  if (length(pair_indexes) == 0) {
    warning("Skipping replicate pairwise output for ", label, "; no shared replicate_index values.", call. = FALSE)
    return(invisible(NULL))
  }

  for (replicate_index in pair_indexes) {
    sample_a <- a$sample_id[a$replicate_index == replicate_index][[1]]
    sample_b <- b$sample_id[b$replicate_index == replicate_index][[1]]
    table <- metadata
    table[[sample_a]] <- counts[metadata$feature_id, sample_a]
    table[[sample_b]] <- counts[metadata$feature_id, sample_b]
    table$log2_fold_change <- safe_log2_ratio(table[[sample_b]], table[[sample_a]])
    table <- table[order(abs(table$log2_fold_change), decreasing = TRUE), , drop = FALSE]
    write_tsv(table, file.path(out_dir, paste0(label, "_replicate_", replicate_index, "_pairwise.tsv")))
  }
  invisible(TRUE)
}
