read_experiment_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read experiment configuration.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Config file does not exist: ", path, call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)
  validate_experiment_config(cfg)
  cfg
}

validate_experiment_config <- function(cfg) {
  required_top_level <- c("project", "analysis", "conditions")
  missing_top_level <- setdiff(required_top_level, names(cfg))
  if (length(missing_top_level) > 0) {
    stop("Config is missing required sections: ", paste(missing_top_level, collapse = ", "), call. = FALSE)
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
