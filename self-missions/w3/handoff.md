# Week 3 Handoff

## Current State

Week 3 implements the final interpretation pass from `w3-plan.md`.

The workflow is:

```bash
python3 self-missions/w3/code/build_w3_analysis.py
```

It reads the Week 2 merged motif table, calculates Spearman correlations and a permutation test, then writes the final report and figures.

## Inputs

- Week 2 merged motif/metric table: `self-missions/w2/output/transcript_motif_counts.tsv`
- Metric used for motif burden: `aag_like_per_kb`
- Target groups inherited from Week 1:
  - `functional_target`
  - `high_clip_target`
  - `low_clip_control`
  - `other`

## Generated Outputs

- `self-missions/w3/output/report.md`
- `self-missions/w3/output/w3-target-control-scatter.png`
- `self-missions/w3/output/w3-motif-burden-boxplot.png`
- `self-missions/w3/output/w3-motif-vs-response.png`
- `self-missions/w3/subdata/w3-summary.json`

## Main Result

The current report concludes that simple linear motif burden does not explain LIN28A binding or translational repression.

- Genes analyzed: `11,461`
- Spearman rho between `aag_like_per_kb` and log2 CLIP enrichment: `0.0143`
- Spearman rho between `aag_like_per_kb` and log2 ribosome-density change: `-0.2037`
- Median `aag_like_per_kb`:
  - high-CLIP targets: `9.864`
  - functional targets: `9.202`
  - low-CLIP controls: `10.239`
- One-sided permutation test for high-CLIP targets having higher motif burden than low-CLIP controls:
  - median difference: `-0.375` motifs/kb
  - p-value: `0.9998`

This supports the Cho et al. interpretation that LIN28A target selection is not explained by motif presence alone. Linear AAG-like motif density may be part of recognition, but hairpin-loop context, transcript architecture, and peri-ER access remain important.

## Figure Review

The current figures were visually checked after regeneration.

- `w3-target-control-scatter.png`: shows the Week 1 target/control grouping in CLIP-vs-ribosome-response space.
- `w3-motif-burden-boxplot.png`: shows high-CLIP and functional targets do not have elevated motif burden relative to low-CLIP controls.
- `w3-motif-vs-response.png`: shows weak motif-vs-CLIP association and a modest negative motif-vs-ribosome-response association.

## Commit Notes

Recent commit style is:

```text
feat: add LIN28A motif analysis workflow and outputs (self-mission, week2)
feat: add LIN28A target group outputs (self-mission, week1)
```

Suggested commit message for the current Week 3 update:

```text
feat: add LIN28A motif interpretation workflow and outputs (self-mission, week3)
```

Suggested files to stage:

```bash
git add self-missions/w3/w3-plan.md \
        self-missions/w3/handoff.md \
        self-missions/w3/code/build_w3_analysis.py \
        self-missions/w3/output/report.md \
        self-missions/w3/output/w3-target-control-scatter.png \
        self-missions/w3/output/w3-motif-burden-boxplot.png \
        self-missions/w3/output/w3-motif-vs-response.png \
        self-missions/w3/subdata/w3-summary.json
```

Do not stage the larger Week 1/2 generated TSVs unless intentionally committing prior outputs:

```text
self-missions/w1/output/target_gene_lists.tsv
self-missions/w2/output/transcript_motif_counts.tsv
self-missions/w2/subdata/selected_transcripts.tsv
self-missions/w2/subdata/transcript_motif_counts_raw.tsv
```

Also check `.gitignore` before committing; it is currently modified independently of Week 3.

## Next Conversation Focus

- Decide whether to keep all three Week 3 PNGs in the commit or commit only the report plus the two most explanatory figures.
- Decide whether Week 2 merged motif output should remain untracked even though Week 3 depends on it.
- If needed, add a small note to `README.md` or the self-mission overview explaining the w1 -> w2 -> w3 dependency chain.
