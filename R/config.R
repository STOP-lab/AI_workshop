read_experiment_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read experiment configuration.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Config file does not exist: ", path, call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)
  cfg$.config_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  cfg$.config_dir <- dirname(cfg$.config_path)
  validate_experiment_config(cfg)
  cfg
}

validate_experiment_config <- function(cfg) {
  required_top_level <- c("project", "analysis")
  missing_top_level <- setdiff(required_top_level, names(cfg))
  if (length(missing_top_level) > 0) {
    stop("Config is missing required sections: ", paste(missing_top_level, collapse = ", "), call. = FALSE)
  }

  if (!is.null(cfg$inputs$sample_sheet)) {
    return(invisible(TRUE))
  }

  if (is.null(cfg$conditions)) {
    stop("Config requires either inputs.sample_sheet or inline conditions.", call. = FALSE)
  }

  condition_names <- names(cfg$conditions)
  if (length(condition_names) != 2) {
    stop("Exactly two conditions are required for this pipeline version.", call. = FALSE)
  }

  for (condition_name in condition_names) {
    condition <- cfg$conditions[[condition_name]]
    if (is.null(condition$replicates) || length(condition$replicates) == 0) {
      stop("Condition has no replicates: ", condition_name, call. = FALSE)
    }
    for (replicate in condition$replicates) {
      if (is.null(replicate$sample_id) || is.null(replicate$bam)) {
        stop("Each replicate requires sample_id and bam fields.", call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}
