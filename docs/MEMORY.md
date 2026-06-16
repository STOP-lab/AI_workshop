# Project memory

## Stable decisions

- Project language is English.
- Pipeline should be reusable across multiple experiments.
- Input data are aligned POINT-seq BAM files stored and analyzed on a server.
- Input BAM files are stranded.
- Analyses are always performed on human genome hg38.
- Each condition has a shared folder containing biological replicates, usually 2-3.
- Analyses must include replicate-by-replicate comparisons and merged-condition comparisons.
- Preferred implementation language is R with Bioconductor packages where possible.
- Git/GitHub should be used for version control.
- Final analyses will run on the server without Codex access.
- The server is expected to have R and all necessary libraries installed.
- Experiment sample metadata should be provided through an Excel workbook template when practical.
- The current implementation is server-runnable but should be validated on real or pilot POINT-seq BAMs before biological interpretation.

## Required analysis modules

- Differential nascent transcription.
- Intron read coverage normalized to gene expression/nascent transcription.
- Differential intron retention.
- Alternative splicing between conditions.
- Spliced/unspliced ratio and splicing index.
- Differential splicing efficiency.
- BED-like coordinate files with the 100 most significantly changed regions should be produced for each analysis module.
- BigWig tracks should be produced for alternative splicing events and other genome-browser-friendly splicing signals.

## Unknowns to resolve

- Annotation file source and version.
- POINT-seq strand orientation and paired/single-end mode.
- Whether BAM merging is allowed/desirable on the server or count-level merging is preferred.
- Exact BigWig signal definitions for alternative splicing and differential splicing efficiency.
- Server scheduler environment, if any.
