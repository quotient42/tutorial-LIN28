## Project Question

Cho et al. 2012 shows that LIN28A binds RNA motifs such as `AAGNNG`, `AAGNG`, and `UGUG`, but also argues that motif presence alone does not fully explain LIN28A target selection. This project asks:

**Do LIN28A motif counts explain LIN28A CLIP enrichment and translational repression?**

## Main Data

Use existing tutorial data:

- `guided-missions/data/read-counts.txt`
- `guided-missions/data/gencode.gtf`
- optional localization annotation from `guided-missions/w1/code/draw_fig5bs6a.py`

Motif analysis also needs transcript sequences, which are not currently in `guided-missions/data`. Add or generate a transcript FASTA, preferably matching GENCODE mouse vM27.

## Motifs

Scan transcript sequences using DNA-style patterns:

| RNA motif | DNA scan pattern | meaning |
|---|---|---|
| `AAGNNG` | `AAG[ACGT]{2}G` | primary LIN28A motif |
| `AAGNG` | `AAG[ACGT]G` | secondary motif |
| `UGUG` | `TGTG` | less frequent motif |

Normalize motif counts by transcript length, for example motifs per kb.

## Expected Interpretation

Possible outcomes:

- If high-CLIP targets have more motifs, motif sequence likely contributes to LIN28A binding.
- If motif burden only weakly correlates with CLIP enrichment, motif presence alone is not enough.
- If integral membrane genes show high CLIP enrichment without much higher motif burden, this supports the Cho 2012 model that LIN28A target selection depends strongly on ER-proximal localization.