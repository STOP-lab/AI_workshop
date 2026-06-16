count_features <- function(features, sample_table, cfg, label) {
  require_package("GenomicAlignments")
  require_package("Rsamtools")
  require_package("SummarizedExperiment")

  ensure_dir(file.path(project_output_dir(cfg), "counts"))
  bam_files <- Rsamtools::BamFileList(sample_table$bam, yieldSize = cfg$analysis$yield_size %||% 2000000)
  names(bam_files) <- sample_table$sample_id

  args <- list(
    features = features,
    reads = bam_files,
    mode = cfg$analysis$count_mode %||% "Union",
    singleEnd = !is_paired_end(sample_table, cfg),
    ignore.strand = FALSE,
    fragments = is_paired_end(sample_table, cfg),
    inter.feature = TRUE
  )
  preprocess <- strand_preprocess_function(sample_table, cfg)
  if (!is.null(preprocess)) {
    args$preprocess.reads <- preprocess
  }
  counts <- do.call(GenomicAlignments::summarizeOverlaps, args)

  count_matrix <- SummarizedExperiment::assay(counts)
  colnames(count_matrix) <- sample_table$sample_id
  rownames(count_matrix) <- features$feature_id %||% names(features)

  write_tsv(
    data.frame(feature_id = rownames(count_matrix), count_matrix, check.names = FALSE),
    file.path(project_output_dir(cfg), "counts", paste0(label, "_counts.tsv"))
  )
  count_matrix
}

count_splice_junctions <- function(sample_table, cfg) {
  require_package("GenomicAlignments")
  require_package("Rsamtools")
  require_package("GenomicRanges")

  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "counts", "junctions"))
  all_junctions <- list()

  for (i in seq_len(nrow(sample_table))) {
    sample_id <- sample_table$sample_id[[i]]
    bam <- Rsamtools::BamFile(sample_table$bam[[i]], yieldSize = cfg$analysis$yield_size %||% 2000000)
    message("Counting splice junctions: ", sample_id)
    junctions <- tryCatch(
      GenomicAlignments::summarizeJunctions(bam, with.revmap = FALSE),
      error = function(e) {
        warning("Junction counting failed for ", sample_id, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(junctions) || length(junctions) == 0) {
      next
    }
    junctions$feature_id <- paste(
      as.character(GenomicRanges::seqnames(junctions)),
      GenomicRanges::start(junctions),
      GenomicRanges::end(junctions),
      as.character(GenomicRanges::strand(junctions)),
      sep = ":"
    )
    junctions$sample_id <- sample_id
    all_junctions[[sample_id]] <- junctions
    write_tsv(feature_metadata(junctions), file.path(out_dir, paste0(sample_id, "_junctions.tsv")))
  }

  all_junctions
}
