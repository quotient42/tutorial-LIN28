import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

cnts = pd.read_csv('../binfo1-work/read-counts.txt', sep='\t', comment='#', index_col=0)
# print(cnts.head())
cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
cnts['rden_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']) / (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam'])
# print(cnts.head())

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.scatter(np.log2(cnts['clip_enrichment']),
           np.log2(cnts['rden_change']))
fig.savefig("example.png")