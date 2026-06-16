# Server requirements

This pipeline is intended to run on the analysis server without Codex.

## System software

Required:

- Git, to clone/update the repository.
- R, preferably a recent release compatible with the installed Bioconductor version.
- Pandoc, required by `rmarkdown`.
- LaTeX distribution for PDF report rendering, for example TeX Live.

Strongly recommended:

- `samtools`, for checking, sorting, and indexing BAM files before running the pipeline.
- Shell utilities available on most Linux servers: `bash`, `sort`, `mkdir`, `cp`.

Optional:

- Snakemake, only if running the optional `workflow/Snakefile` layer.
- IGV, JBrowse, UCSC Genome Browser utilities, or another genome browser for inspecting BED/bedGraph/BigWig outputs.
- UCSC `bedGraphToBigWig`, only if you want command-line BigWig conversion outside R. The pipeline attempts BigWig export with `rtracklayer`.

## R packages from CRAN

- `argparse`
- `yaml`
- `readxl`
- `data.table`
- `dplyr`
- `ggplot2`
- `knitr`
- `rmarkdown`
- `testthat`

## R packages from Bioconductor

- `BiocManager`
- `BiocParallel`
- `DESeq2`
- `DEXSeq`
- `edgeR`
- `GenomicAlignments`
- `GenomicFeatures`
- `GenomicRanges`
- `GenomeInfoDb`
- `IRanges`
- `Rsamtools`
- `rtracklayer`
- `trackViewer`
- `SummarizedExperiment`

## Reference/input files

- Human hg38 GTF/GFF3 annotation file.
- hg38 chromosome sizes file using the same chromosome naming style as the BAMs and annotation.
- Coordinate-sorted stranded BAM files.
- Matching `.bai` index for every BAM file.
- Completed Excel input workbook based on `config/input_files_template.xlsx`.
- Experiment YAML config based on `config/experiment_template.yml`.

## Quick installation check in R

```r
packages <- c(
  "argparse", "yaml", "readxl", "data.table", "dplyr", "ggplot2",
  "knitr", "rmarkdown", "testthat", "BiocParallel", "DESeq2",
  "DEXSeq", "edgeR", "GenomicAlignments", "GenomicFeatures",
  "GenomicRanges", "GenomeInfoDb", "IRanges", "Rsamtools",
  "rtracklayer", "trackViewer", "SummarizedExperiment"
)

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
missing
```

If `missing` is empty, the R package set is available.

## PDF report check

```r
rmarkdown::pandoc_available()
tinytex::is_tinytex()
```

The second command is useful only if TinyTeX is used. On servers with TeX Live installed system-wide, check with:

```bash
pdflatex --version
```
