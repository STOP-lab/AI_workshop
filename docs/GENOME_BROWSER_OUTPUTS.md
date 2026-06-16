# Genome browser outputs

The pipeline should generate genome-browser-friendly outputs for each analysis where a genomic interval or signal can be assigned.

## Coordinate files

For each major analysis module, export the 100 most significantly changed regions between conditions.

Recommended files:

```text
results/<experiment>/browser/top_regions/differential_nascent_transcription.top100.bed
results/<experiment>/browser/top_regions/intron_retention.top100.bed
results/<experiment>/browser/top_regions/alternative_splicing.top100.bed
results/<experiment>/browser/top_regions/splicing_efficiency.top100.bed
```

Use BED6 as the minimum:

```text
chrom  start  end  name  score  strand
```

When event structure matters, prefer BED12. Alternative splicing events are especially good candidates for BED12 because blocks can represent exons or retained intronic segments.

## BigWig tracks

BigWig tracks should be generated for signal-style visualization:

- alternative splicing event support,
- intron retention signal,
- normalized intron coverage,
- spliced/unspliced signal,
- splicing efficiency signal,
- condition-specific nascent transcription coverage where useful.

Recommended files:

```text
results/<experiment>/browser/bigwig/control.normalized_intron_coverage.bw
results/<experiment>/browser/bigwig/treated.normalized_intron_coverage.bw
results/<experiment>/browser/bigwig/control.intron_retention_signal.bw
results/<experiment>/browser/bigwig/treated.intron_retention_signal.bw
results/<experiment>/browser/bigwig/delta_splicing_efficiency.bw
results/<experiment>/browser/bigwig/alternative_splicing_event_signal.bw
```

Generate bedGraph intermediates first, sort them, then convert to BigWig with hg38 chromosome sizes. The chromosome naming style must match the BAMs and annotation, for example `chr1` versus `1`.

## Other useful genome browser visualizations

The following outputs should be visually informative in a genome browser:

- Differential intron retention: highlights retained introns and condition-specific intronic signal.
- Normalized intron coverage: helps distinguish real intron retention from increased nascent transcription.
- Spliced/unspliced ratio and splicing index: useful as per-gene or per-intron signal tracks.
- Differential splicing efficiency: useful as delta tracks across introns or gene bodies.
- Differential nascent transcription: useful as gene-body coverage or log2 fold-change intervals.
- Splice junction support: useful as junction BED/BEDPE or sashimi-style input if junction extraction is added.

Genome browser tracks should always use hg38 coordinate conventions that match the BAMs and annotation.
