# Week 2 Report: LIN28A Motif Burden

## Goal

This week asks whether simple LIN28A sequence motifs explain the target groups defined in Week 1.
The scan uses DNA equivalents of the reported RNA motifs: `AAGNNG`, `AAGNG`, and `TGTG` for `UGUG`.

## Inputs

- Transcript annotation: `data/gencode.gtf`
- Transcript sequences: `data/gencode.vM27.transcripts.fa.gz`
- Week 1 grouped metrics: `w1/subdata/gene_metrics_grouped.tsv`

## Transcript Choice

One protein-coding transcript was selected per protein-coding gene.
The priority order was lowest transcript support level, APPRIS principal tag when available, then longest exon-composed transcript length.

## Results

- Genes with both Week 1 metrics and selected transcript motifs: 11,461
- Pearson correlation between `aag_like_per_kb` and log2 CLIP enrichment: 0.0331
- Pearson correlation between `aag_like_per_kb` and log2 ribosome-density change: -0.1715
- Median `aag_like_per_kb` in high-CLIP targets: 9.864
- Median `aag_like_per_kb` in functional targets: 9.202
- Median `aag_like_per_kb` in low-CLIP controls: 10.239

## Outputs

- Raw motif table: `subdata/transcript_motif_counts_raw.tsv`
- Merged motif and target table: `output/transcript_motif_counts.tsv`
- Boxplot: `output/w2-motif-burden-boxplot.png`
- Scatterplot: `output/w2-motif-vs-clip-scatter.png`

## Top CLIP-Enriched Genes

| gene_name | group | log2_clip | log2_rden | aag_like_per_kb | aagnng_per_kb | aagng_per_kb | tgtg_per_kb |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1f1 | high_clip_target | 8.039 | -0.785 | 37.787 | 31.039 | 6.748 | 1.350 |
| Ncbp3 | high_clip_target | 6.485 | -1.278 | 11.151 | 6.443 | 4.708 | 9.995 |
| Vmn2r120 | high_clip_target | 5.825 | -0.788 | 7.659 | 3.142 | 4.517 | 8.248 |
| Pgrmc1 | high_clip_target | 5.461 | -0.232 | 9.155 | 5.385 | 3.770 | 5.385 |
| Raet1d | functional_target | 5.336 | 1.383 | 14.949 | 8.655 | 6.294 | 3.934 |
| Cd59b | functional_target | 5.309 | 0.250 | 9.915 | 4.249 | 5.666 | 2.833 |
| Odad4 | high_clip_target | 5.151 | -1.635 | 20.536 | 12.054 | 8.482 | 0.893 |
| Dnajc25 | high_clip_target | 4.821 | -0.014 | 10.440 | 5.966 | 4.474 | 7.457 |
| Tmprss11a | high_clip_target | 4.763 | -1.148 | 12.044 | 3.764 | 8.280 | 4.516 |
| Ikbip | functional_target | 4.731 | 0.424 | 15.060 | 9.036 | 6.024 | 2.259 |

## Interpretation

A weak motif-vs-CLIP correlation would mean that the sequence motifs contribute to LIN28A binding but are not sufficient to explain target selection.
That is the outcome expected from Cho et al. 2012: LIN28A recognizes AAG-rich motifs, yet ER-proximal localization and transcript context help determine which mRNAs become strongly bound and translationally repressed.

## Caveats

- This scan counts linear sequence motifs only; it does not test the hairpin-loop structure reported for the main LIN28A motif.
- One transcript per gene hides isoform-specific motif differences.
- Motif counts are transcript-level, while Week 1 CLIP and Ribo-seq metrics are gene-level.
