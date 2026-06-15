build_sample_table <- function(cfg) {
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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
