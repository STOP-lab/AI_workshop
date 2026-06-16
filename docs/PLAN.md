# Project plan

## Phase 1 - Foundation

- Done: Define experiment configuration schema.
- Done: Add Excel input workbook support for experiment sample metadata.
- Done: Build annotation import and validation utilities.
- Done: Validate BAM presence, indexes, chromosome naming, strandedness, and replicate metadata.
- Done: Establish output directory conventions.

## Phase 2 - Counting and coverage

- Done: Generate gene-level nascent transcription counts.
- Done: Generate intron-level coverage/count matrices.
- Done: Generate splice-junction or exon-boundary evidence tables where supported by POINT-seq data.
- Done: Implement replicate-level and merged-condition count generation.

## Phase 3 - Core differential analyses

- Done: Differential nascent transcription with `DESeq2`.
- Done: Intron coverage normalized to nascent gene expression for each condition.
- Done: Differential intron retention using gene-normalized intron signal and model-based testing.
- First implementation: Alternative splicing analysis using exon-bin differential usage; event-specific caller selection remains open after pilot data review.

## Phase 4 - Splicing efficiency metrics

- Done: Calculate spliced/unspliced proxy ratio and splicing index per intron.
- Done: Compare splicing efficiency between conditions.
- Done: Add replicate-by-replicate comparison mode.
- Done: Add merged-condition comparison mode.

## Phase 5 - Reporting and reuse

- Produce run-level HTML/Quarto report.
- Save machine-readable tables and diagnostic plots.
- Done: Export BED/BED-like coordinate files with the top 100 most significantly changed regions for each analysis module.
- Done: Export bedGraph files and attempt BigWig tracks for alternative splicing events and other browser-friendly splicing signals.
- Document assumptions, package versions, and server execution instructions.
- Ensure the pipeline can run from a cloned GitHub repository on the server without Codex.
- Add tests using small synthetic or downsampled BAM fixtures.

## Open method decisions

- Best statistical model for POINT-seq intron retention and splicing efficiency.
- Whether merged-condition analysis should merge BAMs or merge count matrices.
- Whether alternative splicing should rely on junction counts, exon/intron bins, or external event callers.
- Exact score definitions for BED top-region files and BigWig signal tracks per module.
- Organism-specific annotation cleanup rules.
