prepare_top_regions <- function(result_table, cfg, analysis_name) {
  n <- cfg$outputs$top_regions_n %||% 100
  required <- c("seqnames", "start", "end", "strand")
  missing <- setdiff(required, names(result_table))
  if (length(missing) > 0) {
    stop("Cannot export BED for ", analysis_name, "; missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if ("padj" %in% names(result_table) || "pvalue" %in% names(result_table)) {
    result_table <- result_table[padj_or_pvalue_order(result_table), , drop = FALSE]
  } else if ("log2_fold_change" %in% names(result_table)) {
    result_table <- result_table[order(abs(result_table$log2_fold_change), decreasing = TRUE), , drop = FALSE]
  } else if ("log2FoldChange" %in% names(result_table)) {
    result_table <- result_table[order(abs(result_table$log2FoldChange), decreasing = TRUE), , drop = FALSE]
  }
  result_table <- head(result_table, n)
  if ("padj" %in% names(result_table) || "pvalue" %in% names(result_table)) {
    p <- result_table$padj %||% result_table$pvalue %||% rep(1, nrow(result_table))
    min_positive <- suppressWarnings(min(p[p > 0], na.rm = TRUE))
    if (!is.finite(min_positive)) min_positive <- 1e-300
    p[is.na(p) | p <= 0] <- min_positive
    score <- pmin(1000, round(-log10(p) * 100))
  } else {
    effect <- abs(result_table$log2_fold_change %||% result_table$log2FoldChange %||% rep(0, nrow(result_table)))
    max_effect <- max(effect, na.rm = TRUE)
    if (!is.finite(max_effect) || max_effect == 0) max_effect <- 1
    score <- pmin(1000, round(effect / max_effect * 1000))
  }
  name <- result_table$feature_id %||% paste0(analysis_name, "_", seq_len(nrow(result_table)))

  data.frame(
    chrom = result_table$seqnames,
    start = pmax(as.integer(result_table$start) - 1L, 0L),
    end = as.integer(result_table$end),
    name = name,
    score = score,
    strand = result_table$strand,
    stringsAsFactors = FALSE
  )
}

export_top_changed_regions <- function(result_table, cfg, analysis_name) {
  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "browser", "top_regions"))
  regions <- prepare_top_regions(result_table, cfg, analysis_name)
  bed_path <- file.path(out_dir, paste0(analysis_name, ".top", cfg$outputs$top_regions_n %||% 100, ".bed"))
  write.table(regions, bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  score_column <- if ("log2FoldChange" %in% names(result_table)) "log2FoldChange" else if ("log2_fold_change" %in% names(result_table)) "log2_fold_change" else NA_character_
  if (!is.na(score_column)) {
    if ("padj" %in% names(result_table) || "pvalue" %in% names(result_table)) {
      ordered <- result_table[padj_or_pvalue_order(result_table), , drop = FALSE]
    } else {
      ordered <- result_table[order(abs(result_table[[score_column]]), decreasing = TRUE), , drop = FALSE]
    }
    ordered <- head(ordered, cfg$outputs$top_regions_n %||% 100)
    bedgraph <- data.frame(
      chrom = ordered$seqnames,
      start = pmax(as.integer(ordered$start) - 1L, 0L),
      end = as.integer(ordered$end),
      score = ordered[[score_column]],
      stringsAsFactors = FALSE
    )
    bedgraph <- bedgraph[stats::complete.cases(bedgraph), , drop = FALSE]
    write.table(
      bedgraph,
      file.path(out_dir, paste0(analysis_name, ".top", cfg$outputs$top_regions_n %||% 100, ".bedgraph")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  bed_path
}

export_bigwig_from_table <- function(result_table, cfg, track_name, score_column = "log2FoldChange") {
  require_package("GenomicRanges")
  require_package("GenomeInfoDb")
  require_package("IRanges")
  require_package("rtracklayer")

  if (!isTRUE(cfg$outputs$genome_browser_tracks$bigwig %||% TRUE)) {
    return(NA_character_)
  }
  if (!score_column %in% names(result_table)) {
    warning("Skipping BigWig track ", track_name, "; missing score column: ", score_column, call. = FALSE)
    return(NA_character_)
  }

  rows <- stats::complete.cases(result_table[, c("seqnames", "start", "end", score_column), drop = FALSE])
  result_table <- result_table[rows, , drop = FALSE]
  if (nrow(result_table) == 0) {
    warning("Skipping BigWig track ", track_name, "; no complete rows.", call. = FALSE)
    return(NA_character_)
  }

  gr <- GenomicRanges::GRanges(
    seqnames = result_table$seqnames,
    ranges = IRanges::IRanges(start = result_table$start, end = result_table$end),
    strand = result_table$strand %||% "*",
    score = result_table[[score_column]]
  )
  chrom_sizes <- read_chrom_sizes(cfg)
  if (!is.null(chrom_sizes)) {
    common <- intersect(names(chrom_sizes), GenomeInfoDb::seqlevels(gr))
    GenomeInfoDb::seqlengths(gr)[common] <- chrom_sizes[common]
  }
  gr <- sort(gr)

  out_dir <- ensure_dir(file.path(project_output_dir(cfg), "browser", "bigwig"))
  bedgraph_path <- file.path(out_dir, paste0(track_name, ".bedgraph"))
  bedgraph <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr) - 1L,
    end = GenomicRanges::end(gr),
    score = gr$score,
    stringsAsFactors = FALSE
  )
  write.table(bedgraph, bedgraph_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  bw_path <- file.path(out_dir, paste0(track_name, ".bw"))
  tryCatch(
    rtracklayer::export(gr, bw_path, format = "BigWig"),
    error = function(e) {
      warning(
        "BigWig export failed for ", track_name, ". BedGraph was written to ",
        bedgraph_path, ". Error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
  bw_path
}

read_chrom_sizes <- function(cfg) {
  if (is.null(cfg$project$chrom_sizes) || is.na(cfg$project$chrom_sizes) || cfg$project$chrom_sizes == "") {
    return(NULL)
  }
  path <- resolve_path(cfg$project$chrom_sizes, cfg$.config_dir %||% getwd(), must_work = FALSE)
  if (!file.exists(path)) {
    warning("Chromosome sizes file not found; BigWig export may fail: ", path, call. = FALSE)
    return(NULL)
  }
  sizes <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  stats::setNames(as.numeric(sizes[[2]]), sizes[[1]])
}

export_bigwig_tracks <- function(results, cfg) {
  tracks <- list()
  if (!is.null(results$alternative_splicing$table)) {
    tracks$alternative_splicing_event_signal <- export_bigwig_from_table(
      results$alternative_splicing$table,
      cfg,
      "alternative_splicing_event_signal",
      "log2FoldChange"
    )
  }
  if (!is.null(results$intron_retention$table)) {
    tracks$intron_retention_signal <- export_bigwig_from_table(
      results$intron_retention$table,
      cfg,
      "intron_retention_signal",
      "log2FoldChange"
    )
  }
  if (!is.null(results$splicing_efficiency$table)) {
    tracks$differential_splicing_efficiency <- export_bigwig_from_table(
      results$splicing_efficiency$table,
      cfg,
      "differential_splicing_efficiency",
      "log2FoldChange"
    )
  }
  tracks
}
