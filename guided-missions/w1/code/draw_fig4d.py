import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

os.makedirs('../output', exist_ok=True)

cnts = pd.read_csv('../data/read-counts.txt', sep='\t', comment='#', index_col=0)

cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
cnts['rden_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']) / \
                     (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam'])

valid = cnts[(cnts['clip_enrichment'] > 0) & (cnts['rden_change'] > 0)].copy()
valid['log2_clip'] = np.log2(valid['clip_enrichment'])
valid['log2_rden'] = np.log2(valid['rden_change'])
valid = valid.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2_clip', 'log2_rden'])

r, _ = pearsonr(valid['log2_clip'], valid['log2_rden'])

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
ax.scatter(valid['log2_clip'], valid['log2_rden'],
           s=4, c='black', alpha=0.3, edgecolors='none')

ax.axhline(0, color='gray', lw=0.5)
ax.axvline(0, color='gray', lw=0.5)
ax.grid(True, which='major', linestyle=':', linewidth=0.5, color='gray', alpha=0.6)
ax.set_axisbelow(True)

ax.set_xlabel(r'LIN28A CLIP enrichment (log$_2$)')
ax.set_ylabel('Ribosome density change\nupon ' + r'$\it{Lin28a}$' + r' knockdown (log$_2$)')
ax.set_title('CLIP and ribosome footprinting upon\n' + r'$\it{Lin28a}$' + ' knockdown',
             fontsize=10, loc='left')

ax.text(0.97, 0.05, f'r = {r:.4f}', transform=ax.transAxes,
        fontsize=9, ha='right')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.tight_layout()
fig.savefig('../output/fig4d.png', dpi=200)
