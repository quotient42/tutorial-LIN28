## Week 2: Count LIN28A Motifs

Goal: connect transcript sequence motifs to LIN28A target groups.

Tasks:

1. Obtain or generate transcript FASTA.
2. Choose one transcript per gene, for example the longest transcript or a transcript-support-level-1 transcript.
3. Count `AAGNNG`, `AAGNG`, and `UGUG` motifs in each transcript.
4. Normalize motif counts by transcript length.
5. Merge motif counts with the Week 1 target table.

Deliverables:

- `transcript_motif_counts.tsv`
- boxplot comparing motif burden between high-CLIP targets and controls
- scatterplot of motif burden versus log2 CLIP enrichment
