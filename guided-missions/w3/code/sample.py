"""
Mission 3 — CIMS (Crosslinking-Induced Mutation Sites) analysis for Mirlet7g.

각 position 별로
  1) pileup base 문자열에서 read마다 어떤 염기가 호출되었는지 카운트하고,
  2) Shannon entropy 를 계산한 뒤,
  3) UCSC Genome Browser 에 올릴 수 있는 bedgraph (4 columns) 로 저장합니다.
"""

import math
import re
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
SUBDATA_DIR = SCRIPT_DIR.parent / 'subdata'
OUTPUT_DIR = SCRIPT_DIR.parent / 'output'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PILEUP_PATH = SUBDATA_DIR / 'CLIP-let7g-gene.pileup'
BEDGRAPH_PATH = OUTPUT_DIR / 'CLIP-let7g-entropy.bedgraph'
COUNTS_PATH = OUTPUT_DIR / 'CLIP-let7g-basecounts.tsv'


# ---------------------------------------------------------------------------
# 1. pileup 로드 & 정리
# ---------------------------------------------------------------------------

pileup = pd.read_csv(
    PILEUP_PATH,
    sep='\t',
    names=['chrom', 'pos', '_ref', 'count', 'basereads', 'quals'],
)

# `^X` 다음에 오는 한 글자는 mapping quality 이므로 둘 다 통째로 제거한다.
# `+N...` / `-N...` 으로 표시되는 indel 도 제거한다. 그 외 영향 없는 marker (<,>,$,*,#) 도 제거.
indel_re = re.compile(r'[+-](\d+)')
read_start_re = re.compile(r'\^.')
skip_chars_re = re.compile(r'[<>$*#]')


def clean_basereads(s: str) -> str:
    s = read_start_re.sub('', s)

    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c in '+-':
            m = indel_re.match(s, i)
            length = int(m.group(1))
            i = m.end() + length
            continue
        out.append(c)
        i += 1
    cleaned = ''.join(out)
    return skip_chars_re.sub('', cleaned)


pileup['matches'] = pileup['basereads'].apply(clean_basereads)


# ---------------------------------------------------------------------------
# 2. position 별 base count
# ---------------------------------------------------------------------------
#
# reference base 가 'N' 으로 표시되어 있어서 `.`/`,` 가 의미를 갖지 않는다.
# 따라서 A/C/G/T 호출만 세면 충분하다.

BASES = ['A', 'C', 'G', 'T']


def base_counts(matches: str) -> dict:
    up = matches.upper()
    return {b: up.count(b) for b in BASES}


counts_df = pileup['matches'].apply(base_counts).apply(pd.Series)
counts_df['total'] = counts_df[BASES].sum(axis=1)


# ---------------------------------------------------------------------------
# 3. Shannon entropy
# ---------------------------------------------------------------------------

def shannon_entropy(row) -> float:
    total = row['total']
    if total == 0:
        return 0.0
    h = 0.0
    for b in BASES:
        p = row[b] / total
        if p > 0:
            h -= p * math.log2(p)
    return h


counts_df['entropy'] = counts_df.apply(shannon_entropy, axis=1)


# ---------------------------------------------------------------------------
# 4. 결과 저장
# ---------------------------------------------------------------------------

result = pd.concat(
    [pileup[['chrom', 'pos']].reset_index(drop=True), counts_df.reset_index(drop=True)],
    axis=1,
)

# 사람이 확인하기 좋은 base count + entropy TSV
result.to_csv(COUNTS_PATH, sep='\t', index=False)

# bedgraph: chrom, start(0-based), end, entropy
bedgraph = pd.DataFrame({
    'chrom': result['chrom'],
    'start': result['pos'] - 1,
    'end': result['pos'],
    'entropy': result['entropy'].round(6),
})

header = (
    'track type=bedGraph name="CLIP-let7g-entropy" '
    'description="Shannon entropy of CLIP-seq base calls at Mirlet7g" '
    'visibility=full color=40,40,200 altColor=200,40,40 '
    'autoScale=on viewLimits=0:2 priority=20\n'
)
with open(BEDGRAPH_PATH, 'w') as fh:
    fh.write(header)
    bedgraph.to_csv(fh, sep='\t', index=False, header=False)


# ---------------------------------------------------------------------------
# 5. 콘솔 요약
# ---------------------------------------------------------------------------

print(f'positions analyzed : {len(result)}')
print(f'mean total depth   : {result["total"].mean():.1f}')
print(f'mean entropy       : {result["entropy"].mean():.4f}')
print(f'max  entropy       : {result["entropy"].max():.4f} '
      f'@ chr9:{int(result.loc[result["entropy"].idxmax(), "pos"])}')
print()
print('top 10 positions by entropy:')
top = result.sort_values('entropy', ascending=False).head(10)
print(top[['chrom', 'pos', 'A', 'C', 'G', 'T', 'total', 'entropy']].to_string(index=False))
print()
print(f'wrote: {COUNTS_PATH.relative_to(SCRIPT_DIR.parent)}')
print(f'wrote: {BEDGRAPH_PATH.relative_to(SCRIPT_DIR.parent)}')
