# POINT-seq splicing analysis pipeline

Reusable R/Bioconductor-oriented pipeline for comparing pre-mRNA splicing patterns between two experimental conditions using aligned POINT-seq BAM files.

The pipeline is intended to run directly on the analysis server without Codex access. Codex is used only during development; server execution should rely on Git, R, installed R/Bioconductor libraries, and standard command-line tools.

## Goals

- Compare nascent transcription between two conditions.
- Quantify intron read coverage normalized to nascent gene expression.
- Test differential intron retention.
- Detect alternative splicing changes between conditions.
- Estimate spliced/unspliced ratios and splicing index per condition.
- Test differential splicing efficiency.
- Run both replicate-by-replicate comparisons and merged-condition comparisons.

## Expected inputs

- Coordinate-sorted and indexed BAM files for nascent transcriptomics data.
- Stranded BAM files. The pipeline will keep strand orientation configurable per experiment.
- One folder per condition, containing 2-3 biological replicates.
- Excel input workbook listing sample IDs, condition labels, replicate numbers, BAM paths, and optional BAI paths.
- Genome annotation in GTF/GFF3 format.
- Human hg38 genome annotation in GTF/GFF3 format.
- Reference genome metadata and chromosome naming compatible with hg38 BAM files.
- BAM index files (`.bai`) for efficient region-based access.

Example:

```text
server_project/
  condition_A/
    A_rep1.bam
    A_rep1.bam.bai
    A_rep2.bam
    A_rep2.bam.bai
  condition_B/
    B_rep1.bam
    B_rep1.bam.bai
    B_rep2.bam
    B_rep2.bam.bai
```

## Planned method stack

Primary implementation language: R.

Core Bioconductor packages to evaluate and use:

- `Rsubread` or `Rsamtools` for BAM-level counting.
- `GenomicFeatures`, `GenomicRanges`, `rtracklayer`, `IRanges` for annotation handling.
- `DESeq2` or `edgeR` for differential nascent transcription.
- `DEXSeq` for exon/intron-bin usage style testing where appropriate.
- `IsoformSwitchAnalyzeR`, `FRASER`, `limma`, or dedicated junction/count-based tools where suitable for alternative splicing and splicing efficiency modules.
- `BiocParallel` for scalable server execution.

The first implementation should keep package choices modular because POINT-seq nascent RNA has different assumptions than steady-state RNA-seq.

## Repository layout

```text
config/          Experiment configuration templates.
docs/            Project plan, journal, memory, and method notes.
env/             R/Bioconductor environment definitions.
R/               Reusable R functions.
scripts/         Command-line entry points.
workflow/        Workflow orchestration files.
tests/           Small tests and fixture documentation.
results/         Generated outputs, ignored by Git.
logs/            Runtime logs, ignored by Git.
```

## Quick start

1. Copy `config/experiment_template.yml` to a new experiment-specific file.
2. Copy `config/input_files_template.xlsx` to an experiment-specific workbook and fill in the `Samples` sheet.
3. Fill in annotation, output directory, and sample workbook path in the YAML config.
4. Confirm whether the stranded BAMs should be interpreted as forward/sense or reverse/antisense.
5. Install the R environment described in `env/bioconductor_packages.R`.
6. Run the planned pipeline entry point:

```bash
Rscript scripts/run_pipeline.R --config config/my_experiment.yml
```

The current repository is an initial scaffold. Implementation will be added module by module following `docs/PLAN.md`.

For server execution details, see `docs/SERVER_RUNBOOK.md`.

## Output design

Each run should produce:

- Per-replicate QC summaries.
- Merged-condition BAM/index or count summaries.
- Gene-level nascent transcription differential analysis.
- Intron coverage normalized by gene-level nascent transcription.
- Differential intron retention tables.
- Alternative splicing event tables.
- Spliced/unspliced ratio and splicing index tables.
- Differential splicing efficiency tables.
- HTML or Quarto report with plots and methods metadata.

## GitHub

This repository is ready to push to GitHub after choosing a remote repository name:

```bash
git remote add origin <github-repository-url>
git push -u origin main
```
