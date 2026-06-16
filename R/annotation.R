load_hg38_annotation <- function(cfg) {
  require_package("GenomicFeatures")
  require_package("GenomicRanges")
  require_package("rtracklayer")

  gtf <- resolve_path(cfg$project$annotation_gtf, cfg$.config_dir %||% getwd(), must_work = TRUE)
  message("Loading annotation: ", gtf)
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
  genes <- GenomicFeatures::genes(txdb)
  exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
  list(txdb = txdb, genes = genes, exons_by_gene = exons_by_gene)
}

make_gene_features <- function(annotation) {
  genes <- annotation$genes
  genes$feature_id <- names(genes)
  genes$gene_id <- names(genes)
  names(genes) <- genes$feature_id
  genes
}

make_exon_bin_features <- function(annotation) {
  require_package("GenomicRanges")
  exon_list <- GenomicRanges::reduce(annotation$exons_by_gene)
  exon_bins <- unlist(exon_list, use.names = FALSE)
  gene_ids <- rep(names(exon_list), lengths(exon_list))
  exon_bins$gene_id <- gene_ids
  exon_bins$feature_id <- paste0(gene_ids, ":exon_bin_", ave(seq_along(gene_ids), gene_ids, FUN = seq_along))
  names(exon_bins) <- exon_bins$feature_id
  exon_bins
}

make_intron_features <- function(annotation) {
  require_package("GenomicRanges")
  genes <- make_gene_features(annotation)
  exons <- GenomicRanges::reduce(annotation$exons_by_gene)
  shared <- intersect(names(genes), names(exons))
  intron_parts <- vector("list", length(shared))
  names(intron_parts) <- shared

  for (gene_id in shared) {
    gene <- genes[match(gene_id, names(genes))]
    introns <- GenomicRanges::setdiff(gene, exons[[gene_id]], ignore.strand = FALSE)
    if (length(introns) > 0) {
      introns$gene_id <- gene_id
      intron_parts[[gene_id]] <- introns
    }
  }

  intron_parts <- intron_parts[lengths(intron_parts) > 0]
  introns <- unlist(GenomicRanges::GRangesList(intron_parts), use.names = FALSE)
  gene_ids <- introns$gene_id
  introns$feature_id <- paste0(gene_ids, ":intron_", ave(seq_along(gene_ids), gene_ids, FUN = seq_along))
  names(introns) <- introns$feature_id
  introns
}

feature_metadata <- function(features) {
  feature_id <- features$feature_id
  if (is.null(feature_id)) feature_id <- names(features)
  if (is.null(feature_id)) feature_id <- paste0("feature_", seq_along(features))

  gene_id <- features$gene_id
  if (is.null(gene_id)) gene_id <- names(features)
  if (is.null(gene_id)) gene_id <- rep(NA_character_, length(features))

  data.frame(
    feature_id = as.character(feature_id),
    gene_id = as.character(gene_id),
    seqnames = as.character(GenomicRanges::seqnames(features)),
    start = GenomicRanges::start(features),
    end = GenomicRanges::end(features),
    width = GenomicRanges::width(features),
    strand = as.character(GenomicRanges::strand(features)),
    stringsAsFactors = FALSE
  )
}
