# Project journal

## 2026-06-15

- Created initial project scaffold for a reusable POINT-seq pre-mRNA splicing comparison pipeline.
- Captured project scope: two experimental conditions, 2-3 replicates per condition, server-side BAM inputs, replicate-level and merged-condition analyses.
- Selected R/Bioconductor as the preferred implementation direction.
- Added initial configuration template, repository layout, README, plan, memory, and placeholder scripts.
- Recorded human hg38 as the fixed genome build for analyses.
- Recorded that input BAM files should be treated as stranded; exact strand orientation remains configurable.
- Recorded that final analysis will run on the server without Codex, using installed R and R/Bioconductor libraries.
- Added a server runbook for clone/configure/run workflow.
- Added an Excel input template for sample/BAM metadata and wired the R scaffold to read it.
- Added requirements for top-100 BED-like coordinate files and BigWig genome browser tracks.

## Next entry should record

- Real server directory conventions.
- Example BAM names and annotation source.
- Confirmed strand orientation and library layout for POINT-seq datasets.
- Preferred GitHub repository URL.
