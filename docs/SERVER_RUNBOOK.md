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

- condition folders,
- replicate BAM files,
- hg38 GTF/GFF3 annotation,
- desired output directory,
- strand orientation,
- server resource settings.

## Run

```bash
Rscript scripts/run_pipeline.R --config config/my_experiment.yml
```

For future workflow-managed runs, the optional Snakemake entry point is:

```bash
snakemake -s workflow/Snakefile --cores 8
```

## Practical checks before running

- Confirm BAMs are coordinate-sorted.
- Confirm every BAM has a `.bai` index.
- Confirm BAM chromosome names match the hg38 annotation style, for example `chr1` versus `1`.
- Confirm the configured strand orientation matches the POINT-seq library preparation.
- Confirm output and log directories are writable.
