# Final report plot data summary

The original cached TSV files used to build the figures were too large to upload with the report. They were reduced to the small tables below, keeping only the values actually needed by `final-report/code/figure1.py` and `final-report/code/figure2.py`.

## Reduced plot inputs

- `figure_gene_points.tsv`: per-gene group label, log2 CLIP enrichment, log2 ribosome density change, and AAG family motif density. Used for Figure 1A, Figure 1B, and Figure 2A.
- `figure_positional_region_fractions.tsv`: per-gene group label, CDS motif fraction, and 3UTR motif fraction. Used for Figure 2B.
- `figure_positional_metaplot_density.tsv`: 40-bin transcript-relative AAG family motif density by group. Used for Figure 2C.
- `w3-summary.json`: compact summary statistics used for annotations.
- `data_summary.json`: machine-readable description of this reduction.

## Original source sizes summarized

- `transcript_motif_counts.tsv`: 11,461 rows, 36 columns.
- `positional_context_gene_summary.tsv`: 11,461 rows, 28 columns.
- `positional_context_sites.tsv`: summarized into 40 bins for each plotted group.

The full TSV files remain derivable from the self-mission analysis pipeline, but are not required to regenerate the final report figures.
