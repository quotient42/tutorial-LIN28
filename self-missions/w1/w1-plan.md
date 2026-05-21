## Week 1: Define LIN28A Target Groups

Goal: make a clean gene-level analysis table.

Tasks:

1. Load `read-counts.txt`.
2. Calculate:
   - `clip_enrichment = CLIP-35L33G / RNA-control`
   - `rden_change = (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc)`
3. Convert both to log2 scale.
4. Remove zero, infinite, and low-count rows.
5. Define groups:
   - high-CLIP targets: top 5-10% by log2 CLIP enrichment;
   - low-CLIP controls: low CLIP enrichment with enough RNA expression;
   - functional targets: high CLIP and positive ribosome-density change after Lin28a knockdown.

Deliverables:

- `target_gene_lists.tsv`
- Figure 4D-style scatterplot with target/control groups highlighted