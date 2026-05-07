import os
import logging
from collections import defaultdict
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

DATA_DIR = '../../data'
OUTPUT_DIR = '../output'
WINDOW = (-50, 50)

SAMPLES = [
    ('siLuc',    'fivepcounts-filtered-RPF-siLuc.txt'),
    ('siLin28a', 'fivepcounts-filtered-RPF-siLin28a.txt'),
]

os.makedirs(OUTPUT_DIR, exist_ok=True)


def accumulate_relative_counts(path, window):
    lo, hi = window
    counts = defaultdict(int)
    with open(path) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            five_p = int(fields[1])
            count = int(fields[3])
            start_codon = int(fields[8])
            rel = five_p - start_codon
            if lo <= rel <= hi:
                counts[rel] += count
    return counts


available = [(label, os.path.join(DATA_DIR, fname))
             for label, fname in SAMPLES
             if os.path.exists(os.path.join(DATA_DIR, fname))]

if not available:
    raise SystemExit('No fivepcounts-filtered-RPF-*.txt files found in ' + DATA_DIR)

profiles = [(label, accumulate_relative_counts(path, WINDOW))
            for label, path in available]

fig, axes = plt.subplots(len(profiles), 1, figsize=(6, 2.2 * len(profiles)),
                         sharex=True, sharey=True)
if len(profiles) == 1:
    axes = [axes]

xs = list(range(WINDOW[0], WINDOW[1] + 1))

for ax, (label, counts) in zip(axes, profiles):
    ys = [counts.get(x, 0) / 1000 for x in xs]
    ax.bar(xs, ys, width=0.6, color='black', edgecolor='none')
    ax.axvline(0, color='red', lw=1, label='start codon')
    ax.set_xlim(WINDOW[0] - 0.5, WINDOW[1] + 0.5)
    ax.set_ylabel(f'{label}\nRaw read count\n(x1000)', fontsize=7)
    ax.set_xticks(range(WINDOW[0], WINDOW[1] + 1, 10))
    ax.tick_params(axis='both', labelsize=7)
    ax.tick_params(axis='y', which='both', color='gray', width=0.3, length=0)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.grid(True, which='major', linestyle='-', linewidth=0.3, color='gray', alpha=0.5)
    ax.yaxis.grid(True, which='minor', linestyle='-', linewidth=0.3, color='gray', alpha=0.5)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.legend(loc='upper right', frameon=False, fontsize=7,
              handlelength=1.2, handletextpad=0.4, borderpad=0.3)

axes[-1].set_xlabel("Relative position to start codon of 5'-end of reads",
                    fontsize=7)

fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.suptitle('Ribosome footprint density near start codon',
             fontsize=8, y=0.95)
fig.savefig(os.path.join(OUTPUT_DIR, 'figS5a.png'), dpi=200)
