build_sample_table <- function(cfg) {
  if (!is.null(cfg$inputs$sample_sheet)) {
    return(read_sample_sheet(cfg))
  }

  rows <- list()
  for (condition_name in names(cfg$conditions)) {
    condition <- cfg$conditions[[condition_name]]
    condition_label <- condition$label %||% condition_name
    for (i in seq_along(condition$replicates)) {
      replicate <- condition$replicates[[i]]
      rows[[length(rows) + 1]] <- data.frame(
        sample_id = replicate$sample_id,
        condition = condition_label,
        condition_key = condition_name,
        replicate_index = i,
        bam = replicate$bam,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

read_sample_sheet <- function(cfg) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read Excel input templates.", call. = FALSE)
  }

  sample_sheet_path <- resolve_sample_sheet_path(cfg)
  sheet_name <- cfg$inputs$sample_sheet_name %||% "Samples"
  raw <- as.data.frame(readxl::read_excel(sample_sheet_path, sheet = sheet_name), stringsAsFactors = FALSE)
  names(raw) <- trimws(names(raw))
  raw[] <- lapply(raw, function(column) {
    if (is.character(column)) trimws(column) else column
  })
  raw <- drop_empty_rows(raw)

  required_columns <- c(
    "sample_id",
    "condition_key",
    "condition_label",
    "replicate_index",
    "bam_path",
    "stranded",
    "strand_orientation",
    "paired_end",
    "include"
  )
  missing_columns <- setdiff(required_columns, names(raw))
  if (length(missing_columns) > 0) {
    stop("Sample sheet is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  raw$include <- normalize_yes_no(raw$include, "include")
  raw <- raw[raw$include == "yes", , drop = FALSE]
  if (nrow(raw) == 0) {
    stop("Sample sheet contains no rows with include = yes.", call. = FALSE)
  }

  required_no_blank <- c("sample_id", "condition_key", "condition_label", "replicate_index", "bam_path")
  for (column in required_no_blank) {
    if (any(is.na(raw[[column]]) | trimws(as.character(raw[[column]])) == "")) {
      stop("Sample sheet column contains blank values: ", column, call. = FALSE)
    }
  }

  raw$condition_key <- tolower(trimws(as.character(raw$condition_key)))
  raw$strand_orientation <- tolower(trimws(as.character(raw$strand_orientation)))
  condition_keys <- sort(unique(raw$condition_key))
  if (length(condition_keys) != 2 || !all(condition_keys %in% c("condition_a", "condition_b"))) {
    stop("Sample sheet must include exactly condition_a and condition_b among included rows.", call. = FALSE)
  }

  raw$stranded <- normalize_yes_no(raw$stranded, "stranded")
  if (any(raw$stranded != "yes")) {
    stop("All included samples must have stranded = yes for this project.", call. = FALSE)
  }

  raw$paired_end <- normalize_yes_no(raw$paired_end, "paired_end")
  allowed_orientations <- c("forward", "reverse")
  if (any(!raw$strand_orientation %in% allowed_orientations)) {
    stop("strand_orientation must be one of: ", paste(allowed_orientations, collapse = ", "), call. = FALSE)
  }

  data.frame(
    sample_id = as.character(raw$sample_id),
    condition = as.character(raw$condition_label),
    condition_key = as.character(raw$condition_key),
    replicate_index = as.integer(raw$replicate_index),
    bam = as.character(raw$bam_path),
    bai = if ("bai_path" %in% names(raw)) as.character(raw$bai_path) else NA_character_,
    stranded = as.character(raw$stranded),
    strand_orientation = as.character(raw$strand_orientation),
    paired_end = as.character(raw$paired_end),
    notes = if ("notes" %in% names(raw)) as.character(raw$notes) else NA_character_,
    stringsAsFactors = FALSE
  )
}

resolve_sample_sheet_path <- function(cfg) {
  configured_path <- cfg$inputs$sample_sheet
  candidates <- c(
    configured_path,
    file.path(cfg$.config_dir %||% ".", configured_path)
  )
  for (candidate in unique(candidates)) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Sample sheet does not exist: ", configured_path, call. = FALSE)
}

drop_empty_rows <- function(x) {
  keep <- apply(x, 1, function(row) any(!is.na(row) & trimws(as.character(row)) != ""))
  x[keep, , drop = FALSE]
}

normalize_yes_no <- function(x, column) {
  value <- tolower(trimws(as.character(x)))
  allowed <- c("yes", "no")
  if (any(!value %in% allowed)) {
    stop(column, " must contain only yes/no values.", call. = FALSE)
  }
  value
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}
