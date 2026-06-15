# Project plan

## Phase 1 - Foundation

- Define experiment configuration schema.
- Build annotation import and validation utilities.
- Validate BAM presence, indexes, chromosome naming, strandedness, and replicate metadata.
- Establish output directory conventions and logging.

## Phase 2 - Counting and coverage

- Generate gene-level nascent transcription counts.
- Generate intron-level coverage/count matrices.
- Generate splice-junction or exon-boundary evidence tables where supported by POINT-seq data.
- Implement replicate-level and merged-condition count generation.

## Phase 3 - Core differential analyses

- Differential nascent transcription with `DESeq2` or `edgeR`.
- Intron coverage normalized to nascent gene expression for each condition.
- Differential intron retention using normalized intron signal and model-based testing.
- Alternative splicing analysis using event or bin-level approaches selected after pilot data review.

## Phase 4 - Splicing efficiency metrics

- Calculate spliced reads, unspliced reads, spliced/unspliced ratio, and splicing index per gene/intron.
- Compare splicing efficiency between conditions.
- Add replicate-by-replicate comparison mode.
- Add merged-condition comparison mode.

## Phase 5 - Reporting and reuse

- Produce run-level HTML/Quarto report.
- Save machine-readable tables and diagnostic plots.
- Document assumptions, package versions, and server execution instructions.
- Add tests using small synthetic or downsampled BAM fixtures.

## Open method decisions

- Best statistical model for POINT-seq intron retention and splicing efficiency.
- Whether merged-condition analysis should merge BAMs or merge count matrices.
- Whether alternative splicing should rely on junction counts, exon/intron bins, or external event callers.
- Organism-specific annotation cleanup rules.
