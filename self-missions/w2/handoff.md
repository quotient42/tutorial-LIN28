# Week 2 Handoff

## Current State

Week 2 implements the LIN28A motif-count analysis from `w2-plan.md`.

The workflow is:

```bash
python3 self-missions/w2/code/build_motif_analysis.py
```

It runs:

1. `code/01_count_motifs.py`
2. `code/02_merge_plot_report.py`

## Inputs

- Week 1 grouped metrics: `self-missions/w1/subdata/gene_metrics_grouped.tsv`
- GENCODE annotation: `self-missions/data/gencode.gtf`
- Transcript FASTA cache: `self-missions/data/gencode.vM27.transcripts.fa.gz`

The FASTA is large and should not be staged.

## Generated Outputs

- `output/transcript_motif_counts.tsv`
- `output/w2-motif-burden-boxplot.png`
- `output/w2-motif-vs-clip-scatter.png`
- `output/report.md`
- `subdata/selected_transcripts.tsv`
- `subdata/transcript_motif_counts_raw.tsv`
- `subdata/motif_summary.json`

## Main Result To Check

The current report says motif burden is only weakly associated with LIN28A CLIP enrichment:

- Pearson r between `aag_like_per_kb` and log2 CLIP enrichment: `0.0331`
- Median `aag_like_per_kb`:
  - high-CLIP targets: `9.864`
  - functional targets: `9.202`
  - low-CLIP controls: `10.239`

This supports the intended interpretation: simple linear motif burden alone does not explain LIN28A target selection.

## Next Conversation Focus

Review and polish:

- Check whether `output/report.md` states the biological interpretation clearly enough.
- Inspect both figures for readability, labels, colors, and whether the plotted comparison is convincing.
- Decide whether the boxplot should compare only high-CLIP vs low-CLIP, or keep functional targets as a third group.
- Consider adding a short note that this scan ignores RNA secondary structure, even though Cho et al. emphasize a hairpin-loop motif context.
- Consider whether large generated TSVs should remain untracked and only the code/report/figures should be committed.

