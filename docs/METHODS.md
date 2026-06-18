# Methods implemented in the first R pipeline

This document describes the current server-runnable implementation. The methods are intentionally modular so they can be replaced or extended after pilot POINT-seq data are inspected.

## Annotation

- Genome build: human hg38.
- Annotation input: GTF/GFF3.
- Gene features are imported with `GenomicFeatures`.
- Intron features are generated as gene ranges minus reduced exon ranges.
- Exon-bin features are reduced exonic intervals per gene.

## Counting

- BAM feature counts are generated with `GenomicAlignments::summarizeOverlaps`.
- Counting is strand-aware. Reverse-stranded libraries are handled with `GenomicAlignments::invertStrand`.
- Count matrices are generated for genes, intron intervals, and exon bins.
- Strand-split BAMs are supported by assigning multiple BAM rows to one `biological_sample_id`; raw BAM counts are written first, then counts are summed by biological sample for downstream analyses.
- Splice junction summaries are attempted per BAM with `GenomicAlignments::summarizeJunctions`.

## Differential nascent transcription

- Gene-level counts are tested with `DESeq2`.
- Output includes differential results and normalized counts.

## Intron coverage and intron retention

- Intron counts are normalized against matched gene-level nascent expression.
- Differential intron retention uses intron counts with gene-expression-derived normalization factors in `DESeq2`.
- This is a first-pass model designed to separate intronic signal from overall nascent transcription changes.

## Alternative splicing

- The first implementation uses exon-bin differential usage as a conservative alternative splicing proxy.
- Exon-bin counts are normalized against matched gene-level nascent expression and tested with `DESeq2`.
- Junction summaries are exported where possible and can be used later to add event-specific callers.

## Splicing efficiency

- Per-intron spliced/unspliced proxy metrics are calculated from gene-level normalized signal and intron counts.
- Differential splicing efficiency uses the same gene-normalized intron model as the intron retention step.

## Replicate and merged-condition outputs

- Full replicate-aware DESeq2 models are run across all included replicates.
- Merged-condition count summaries are written for genes, introns, and exon bins.
- Replicate-by-replicate pairwise fold-change tables are written by matching `replicate_index` between conditions.

## Browser outputs

- Top-region BED files are exported for each main analysis.
- bedGraph files are exported for effect-size tracks.
- BigWig export is attempted with `rtracklayer`; bedGraph files remain available if a BigWig conversion fails.

## PDF report

- A PDF report is rendered with `rmarkdown` after the pipeline finishes.
- The report includes input sample metadata, assigned count statistics, significant-result summaries, key top-result tables, and a list of genome browser output files.
