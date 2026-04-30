import os
import ssl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

os.makedirs('../output', exist_ok=True)

ssl._create_default_https_context = ssl._create_unverified_context
mouselocal = pd.read_csv('https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt', sep='\t')

cnts = pd.read_csv('../data/read-counts.txt', sep='\t', comment='#', index_col=0)

cnts['clip_enrichment'] = cnts['CLIP-35L33G.bam'] / cnts['RNA-control.bam']
cnts['rden_change'] = (cnts['RPF-siLin28a.bam'] / cnts['RNA-siLin28a.bam']) / \
                     (cnts['RPF-siLuc.bam'] / cnts['RNA-siLuc.bam'])

cnts['gene_id'] = cnts.index.str.split('.').str[0]

merged = cnts.merge(mouselocal[['gene_id', 'type']], on='gene_id', how='inner')

valid = merged[(merged['clip_enrichment'] > 0) & (merged['rden_change'] > 0)].copy()
valid['log2_clip'] = np.log2(valid['clip_enrichment'])
valid['log2_rden'] = np.log2(valid['rden_change'])
valid = valid.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2_clip', 'log2_rden'])

palette = {
    'nucleus': '#3953a4',
    'integral membrane': '#ed1c24',
    'cytoplasm': '#39b54a',
}
order = ['nucleus', 'integral membrane', 'cytoplasm']
labels = {
    'nucleus': 'Nucleus',
    'integral membrane': 'Integral membrane',
    'cytoplasm': 'Cytoplasm',
}


def draw(ax, data, point_size, alpha, title):
    for term in order:
        sub = data[data['type'] == term]
        ax.scatter(sub['log2_clip'], sub['log2_rden'],
                   s=point_size, c=palette[term], alpha=alpha,
                   edgecolors='none', label=labels[term])
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.grid(True, which='major', linestyle=':', linewidth=0.5,
            color='gray', alpha=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel(r'LIN28A CLIP enrichment (log$_2$)')
    ax.set_ylabel('Ribosome density change\nupon ' + r'$\it{Lin28a}$' + r' knockdown (log$_2$)')
    ax.set_title(title, fontsize=10, loc='left')
    leg = ax.legend(loc='upper left', frameon=True, fontsize=8,
                    markerscale=1.5, handletextpad=0.3)
    leg.get_frame().set_edgecolor('black')
    leg.get_frame().set_linewidth(0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# Figure 5B: 200 random transcripts per term
rng = np.random.default_rng(0)
sampled = (valid.groupby('type', group_keys=False)
                .apply(lambda g: g.sample(n=min(200, len(g)), random_state=rng.integers(1 << 31))))

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
draw(ax, sampled, point_size=12, alpha=0.8, title='Protein localization')
fig.tight_layout()
fig.savefig('../output/fig5b.png', dpi=200)
plt.close(fig)

# Figure S6A: complete version (all transcripts)
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
draw(ax, valid, point_size=4, alpha=0.5, title='Linkage to localization')
fig.tight_layout()
fig.savefig('../output/figS6a.png', dpi=200)
plt.close(fig)
