# Server runbook

This project is designed to run on the analysis server without access to Codex.

## Assumptions

- The server has R installed.
- Required CRAN and Bioconductor packages are installed before running the pipeline.
- BAM files and BAM indexes are available on the server filesystem.
- Human hg38 annotation files are available on the server filesystem.
- The repository is obtained from GitHub or copied to the server.

## First-time setup on the server

```bash
git clone <github-repository-url>
cd splicing_pipeline
```

If package installation is allowed on the server:

```bash
Rscript env/bioconductor_packages.R
```

If package installation is handled by an administrator, use `env/bioconductor_packages.R` as the package manifest.

## Per-experiment setup

Copy the template:

```bash
cp config/experiment_template.yml config/my_experiment.yml
```

Edit `config/my_experiment.yml` to point to:

- the experiment-specific Excel sample workbook,
- hg38 GTF/GFF3 annotation,
- hg38 chromosome sizes file,
- desired output directory,
- strand orientation,
- server resource settings.

Copy and edit the Excel input template:

```bash
cp config/input_files_template.xlsx config/my_experiment_input_files.xlsx
```

Fill in the `Samples` sheet with one included row per replicate. Use absolute server paths for BAM and BAI files.

## Run

```bash
Rscript scripts/run_pipeline.R --config config/my_experiment.yml
```

The run writes count matrices, differential tables, replicate-pairwise summaries, merged-condition summaries, BED/bedGraph browser files, attempted BigWig tracks, a PDF report, and `pipeline_results.rds` under the configured output directory.

For future workflow-managed runs, the optional Snakemake entry point is:

```bash
snakemake -s workflow/Snakefile --cores 8
```

## Practical checks before running

- Confirm BAMs are coordinate-sorted.
- Confirm every BAM has a `.bai` index.
- Confirm BAM chromosome names match the hg38 annotation style, for example `chr1` versus `1`.
- Confirm the hg38 chromosome sizes file uses the same chromosome naming style as the BAMs and annotation.
- Confirm the configured strand orientation matches the POINT-seq library preparation.
- Confirm output and log directories are writable.
- Confirm `DESeq2`, `GenomicAlignments`, `GenomicFeatures`, `Rsamtools`, `rtracklayer`, `readxl`, `rmarkdown`, and related dependencies load in R.
- Confirm the server has a LaTeX installation available for `rmarkdown::pdf_document`.
