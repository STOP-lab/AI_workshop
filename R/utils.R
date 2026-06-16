`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required for this pipeline step.", call. = FALSE)
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

resolve_path <- function(path, base_dir = getwd(), must_work = FALSE) {
  if (is.null(path) || is.na(path) || path == "") {
    return(NA_character_)
  }
  candidates <- c(path, file.path(base_dir, path))
  for (candidate in unique(candidates)) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  if (must_work) {
    stop("Path does not exist: ", path, call. = FALSE)
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = FALSE)
}

project_output_dir <- function(cfg) {
  ensure_dir(resolve_path(cfg$project$output_dir %||% "results/run", cfg$.config_dir %||% getwd()))
}

write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  } else {
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  invisible(path)
}

strand_preprocess_function <- function(sample_table, cfg) {
  orientation <- unique(sample_table$strand_orientation %||% cfg$analysis$strand_orientation %||% "reverse")
  orientation <- orientation[!is.na(orientation)]
  if (length(orientation) > 1) {
    stop("All included samples must use the same strand_orientation for this pipeline version.", call. = FALSE)
  }
  if (length(orientation) == 0 || orientation == "forward") {
    return(NULL)
  }
  require_package("GenomicAlignments")
  GenomicAlignments::invertStrand
}

is_paired_end <- function(sample_table, cfg) {
  values <- unique(tolower(sample_table$paired_end %||% if (isTRUE(cfg$analysis$paired_end)) "yes" else "no"))
  values <- values[!is.na(values)]
  if (length(values) > 1) {
    stop("All included samples must use the same paired_end value for this pipeline version.", call. = FALSE)
  }
  identical(values, "yes")
}

validate_input_files <- function(sample_table, cfg) {
  if (!"bai" %in% names(sample_table)) {
    sample_table$bai <- NA_character_
  }
  missing_bam <- sample_table$bam[!file.exists(sample_table$bam)]
  if (length(missing_bam) > 0) {
    stop("Missing BAM files:\n", paste(missing_bam, collapse = "\n"), call. = FALSE)
  }

  inferred_bai <- ifelse(
    is.na(sample_table$bai) | sample_table$bai == "" | sample_table$bai == "NA",
    paste0(sample_table$bam, ".bai"),
    sample_table$bai
  )
  missing_bai <- inferred_bai[!file.exists(inferred_bai)]
  if (length(missing_bai) > 0) {
    stop("Missing BAM index files:\n", paste(missing_bai, collapse = "\n"), call. = FALSE)
  }

  annotation <- resolve_path(cfg$project$annotation_gtf, cfg$.config_dir %||% getwd(), must_work = TRUE)
  cfg$project$annotation_gtf <- annotation
  if (!is.null(cfg$project$chrom_sizes)) {
    cfg$project$chrom_sizes <- resolve_path(cfg$project$chrom_sizes, cfg$.config_dir %||% getwd(), must_work = FALSE)
  }
  cfg
}

safe_log2_ratio <- function(numerator, denominator, pseudocount = 1) {
  log2((numerator + pseudocount) / (denominator + pseudocount))
}

padj_or_pvalue_order <- function(x) {
  padj <- if ("padj" %in% names(x)) x$padj else NA_real_
  pvalue <- if ("pvalue" %in% names(x)) x$pvalue else NA_real_
  order(is.na(padj), padj, is.na(pvalue), pvalue)
}
