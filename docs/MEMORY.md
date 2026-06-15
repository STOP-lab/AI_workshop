# Project memory

## Stable decisions

- Project language is English.
- Pipeline should be reusable across multiple experiments.
- Input data are aligned POINT-seq BAM files stored and analyzed on a server.
- Analyses are always performed on human genome hg38.
- Each condition has a shared folder containing biological replicates, usually 2-3.
- Analyses must include replicate-by-replicate comparisons and merged-condition comparisons.
- Preferred implementation language is R with Bioconductor packages where possible.
- Git/GitHub should be used for version control.

## Required analysis modules

- Differential nascent transcription.
- Intron read coverage normalized to gene expression/nascent transcription.
- Differential intron retention.
- Alternative splicing between conditions.
- Spliced/unspliced ratio and splicing index.
- Differential splicing efficiency.

## Unknowns to resolve

- Annotation file source and version.
- POINT-seq strandedness and paired/single-end mode.
- Whether BAM merging is allowed/desirable on the server or count-level merging is preferred.
- Server scheduler environment, if any.
