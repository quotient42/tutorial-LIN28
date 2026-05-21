# Week 1 Report: LIN28A Target Groups

## Goal

This week defines target and control gene groups for the LIN28A motif exploration project.
The analysis follows the tutorial formulas for LIN28A CLIP enrichment and ribosome-density change after Lin28a knockdown.

## Inputs

- Count matrix: `data/read-counts.txt`
- Gene annotation: `data/gencode.gtf`

## Filtering and Metrics

- `clip_enrichment = CLIP-35L33G / RNA-control`
- `rden_change = (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc)`
- Genes were kept only if the relevant denominators were nonzero and RNA counts were at least 10 in the three RNA libraries.
- Target/control groups were defined only among protein-coding genes, because Week 2 will scan mRNA transcript sequences.
- Final filtered protein-coding genes: 11,463
- Pearson correlation between log2 CLIP enrichment and log2 ribosome-density change: 0.3681

## Group Definitions

- High-CLIP cutoff: top 10% of filtered genes, log2 CLIP >= 2.400
- Low-CLIP control cutoff: bottom 50% of filtered genes, log2 CLIP <= 0.795
- Functional targets: high-CLIP genes with positive log2 ribosome-density change after Lin28a knockdown

## Group Counts

- Functional targets: 419
- High-CLIP targets without positive ribosome-density change: 728
- Low-CLIP controls: 5,733
- Other filtered genes: 4,583

## Outputs

- Cleaned metrics table: `subdata/gene_metrics_clean.tsv`
- Grouped metrics table: `subdata/gene_metrics_grouped.tsv`
- Gene annotation table: `subdata/gene_annotation.tsv`
- Target/control table: `output/target_gene_lists.tsv`
- Figure: `output/w1-target-groups.png`

## Top Functional Targets

| gene_name   | gene_id            |   log2_clip |   log2_rden |   target_score |
|:------------|:-------------------|------------:|------------:|---------------:|
| Raet1d      | ENSMUSG00000078452 |       5.336 |       1.383 |          6.719 |
| Vldlr       | ENSMUSG00000024924 |       4.268 |       1.45  |          5.719 |
| Chac1       | ENSMUSG00000027313 |       4.161 |       1.451 |          5.612 |
| Vmn2r79     | ENSMUSG00000090362 |       4.206 |       1.403 |          5.609 |
| Cd59b       | ENSMUSG00000068686 |       5.309 |       0.25  |          5.559 |
| Galnt3      | ENSMUSG00000026994 |       3.42  |       1.823 |          5.243 |
| Ulbp1       | ENSMUSG00000079685 |       3.471 |       1.721 |          5.192 |
| Ikbip       | ENSMUSG00000019975 |       4.731 |       0.424 |          5.155 |
| Tmem260     | ENSMUSG00000036339 |       2.688 |       2.462 |          5.15  |
| Serpini1    | ENSMUSG00000027834 |       4.472 |       0.534 |          5.007 |

## Interpretation

Genes with high LIN28A CLIP enrichment and positive ribosome-density change after Lin28a knockdown are plausible direct translational repression targets.
They are strong candidates for Week 2 motif counting because they combine evidence of LIN28A binding with a functional translation response.
Low-CLIP controls provide a comparison group for testing whether LIN28A recognition motifs are enriched among candidate targets.

## Caveats

- These are gene-level counts, so isoform-specific effects are not resolved.
- Ratio-based metrics are sensitive to low counts; the current filter is intentionally conservative for a first-pass target list.
- High CLIP enrichment does not prove direct functional repression unless it is paired with a translation response or further validation.
