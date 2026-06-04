# Week 3 Report: Do LIN28A Motifs Explain Binding and Repression?

## Goal

Week 3 combines the Week 1 target groups with the Week 2 transcript motif counts.
The central test is whether a simple linear count of AAG-like LIN28A motifs is enough to explain LIN28A CLIP enrichment and translational repression.

## Inputs

- Merged Week 2 motif table: `w2/output/transcript_motif_counts.tsv`
- Motif burden metric: `aag_like_per_kb = (AAGNNG + AAGNG motif counts) / transcript kb`

## Results

- Genes with motif burden, CLIP enrichment, and ribosome-density metrics: 11,461
- Spearman correlation between motif burden and log2 CLIP enrichment: 0.0143
- Spearman correlation between motif burden and log2 ribosome-density change: -0.2037
- Median motif burden in high-CLIP targets: 9.864 motifs/kb
- Median motif burden in functional targets: 9.202 motifs/kb
- Median motif burden in low-CLIP controls: 10.239 motifs/kb
- One-sided permutation test for high-CLIP targets having higher motif burden than low-CLIP controls: median difference = -0.375 motifs/kb, p = 0.9998

## Outputs

- Summary statistics: `subdata/w3-summary.json`
- Target/control scatterplot: `output/w3-target-control-scatter.png`
- Motif burden boxplot: `output/w3-motif-burden-boxplot.png`
- Motif-versus-response scatterplots: `output/w3-motif-vs-response.png`

## Top Functional Targets

| gene_name | group | log2_clip | log2_rden | aag_like_per_kb | tgtg_per_kb |
| --- | --- | --- | --- | --- | --- |
| Raet1d | functional_target | 5.336 | 1.383 | 14.949 | 3.934 |
| Vldlr | functional_target | 4.268 | 1.450 | 9.158 | 10.664 |
| Chac1 | functional_target | 4.161 | 1.451 | 9.774 | 6.720 |
| Vmn2r79 | functional_target | 4.206 | 1.403 | 9.110 | 4.758 |
| Cd59b | functional_target | 5.309 | 0.250 | 9.915 | 2.833 |
| Galnt3 | functional_target | 3.420 | 1.823 | 11.068 | 5.405 |
| Ulbp1 | functional_target | 3.471 | 1.721 | 8.126 | 8.650 |
| Ikbip | functional_target | 4.731 | 0.424 | 15.060 | 2.259 |
| Tmem260 | functional_target | 2.688 | 2.462 | 8.387 | 4.923 |
| Serpini1 | functional_target | 4.472 | 0.534 | 8.066 | 13.200 |

## Interpretation

The motif-burden signal is weak. AAG-like motif density is nearly uncorrelated with CLIP enrichment and is negatively associated with ribosome-density change.
High-CLIP targets also do not have higher median AAG-like motif burden than low-CLIP controls in this simplified transcript-level scan.
This supports the Cho et al. model more than a motif-only model: LIN28A recognizes AAG-rich sequence features, but target selection in mESCs also depends on context such as hairpin-loop presentation, transcript architecture, and the peri-ER localization that gives LIN28A preferential access to co-translationally targeted membrane and secretory mRNAs.

## Caveats

- The scan counts linear sequence motifs only and does not model the hairpin-loop context reported for the main LIN28A motif.
- One transcript was selected per gene, so isoform-specific motif differences are hidden.
- Motif counts are transcript-level, while CLIP and ribosome-density estimates are gene-level.
- The permutation test asks only whether high-CLIP targets have higher median motif burden than low-CLIP controls; it does not test motif position, structure, or accessibility.
