"""
Mission 3 — CIMS (Crosslinking-Induced Mutation Sites) analysis.

Mirlet7g 미션과 동일한 파이프라인을 Mirlet7d, Mirlet7f-1 에 대해서도 수행한다.

각 gene 별로
  1) samtools view 로 해당 좌표 범위만 BAM 으로 솎아내고 index 를 만들고,
  2) samtools mpileup 으로 base 정보를 추출한 다음,
  3) gene 좌표 범위 안의 줄만 awk 로 추려내고,
  4) pileup base 문자열에서 read 마다 어떤 염기가 호출되었는지 카운트하고,
  5) Shannon entropy 를 계산해서
  6) 사람용 TSV 와 UCSC bedgraph (4 columns) 로 저장한다.

`seq` conda env 에 들어 있는 samtools 가 PATH 에서 잡혀야 한다.
    conda activate seq
    python code/run_cims.py
"""

import math
import re
import subprocess
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
W3_DIR = SCRIPT_DIR.parent
DATA_DIR = W3_DIR.parent / 'data'
SUBDATA_DIR = W3_DIR / 'subdata'
OUTPUT_DIR = W3_DIR / 'output'
SUBDATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SOURCE_BAM = DATA_DIR / 'CLIP-35L33G.bam'

# gencode.gtf 에서 찾은 mm39 좌표 (gene line 의 start/end, 1-based inclusive).
GENES = [
    {'name': 'Mirlet7d',   'chrom': 'chr13', 'start': 48689488, 'end': 48689590, 'tag': 'let7d'},
    {'name': 'Mirlet7f-1', 'chrom': 'chr13', 'start': 48691305, 'end': 48691393, 'tag': 'let7f1'},
]


# ---------------------------------------------------------------------------
# pileup base 문자열 정리
# ---------------------------------------------------------------------------
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


BASES = ['A', 'C', 'G', 'T']


def base_counts(matches: str) -> dict:
    up = matches.upper()
    return {b: up.count(b) for b in BASES}


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


def run(cmd, **kw) -> None:
    print('$', ' '.join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kw)


def process_gene(gene: dict) -> None:
    name = gene['name']
    chrom = gene['chrom']
    start = gene['start']
    end = gene['end']
    tag = gene['tag']
    region = f'{chrom}:{start}-{end}'

    bam = SUBDATA_DIR / f'CLIP-{tag}.bam'
    pileup_full = SUBDATA_DIR / f'CLIP-{tag}.pileup'
    pileup_gene = SUBDATA_DIR / f'CLIP-{tag}-gene.pileup'
    counts_tsv = OUTPUT_DIR / f'CLIP-{tag}-basecounts.tsv'
    bedgraph = OUTPUT_DIR / f'CLIP-{tag}-entropy.bedgraph'

    print(f'\n=== {name} ({region}) ===')

    # 1) 좌표 범위로 BAM 추출 + index.
    run(['samtools', 'view', '-b', '-o', str(bam), str(SOURCE_BAM), region])
    run(['samtools', 'index', str(bam)])
    n_reads = subprocess.check_output(
        ['samtools', 'view', '-c', str(bam)], text=True
    ).strip()
    print(f'reads in region    : {n_reads}')

    # 2) mpileup (reference 없이 Shannon entropy 만 볼 거라 -f 생략).
    with open(pileup_full, 'w') as fh:
        run(['samtools', 'mpileup', str(bam)], stdout=fh)

    # 3) gene 좌표 범위 안의 줄만 추리기. read 가 region 밖까지 펼쳐지면
    #    pileup 이 region 밖 좌표까지 나오기 때문에 awk 로 자른다.
    with open(pileup_gene, 'w') as fh:
        run(['awk',
             f'$2 >= {start} && $2 <= {end} {{ print $0; }}',
             str(pileup_full)],
            stdout=fh)

    # 4) pileup 로드 & base count.
    pileup = pd.read_csv(
        pileup_gene,
        sep='\t',
        names=['chrom', 'pos', '_ref', 'count', 'basereads', 'quals'],
    )
    if pileup.empty:
        print(f'WARNING: {name} pileup is empty — no reads in {region}')
        return

    pileup['matches'] = pileup['basereads'].apply(clean_basereads)

    counts_df = pileup['matches'].apply(base_counts).apply(pd.Series)
    counts_df['total'] = counts_df[BASES].sum(axis=1)
    counts_df['entropy'] = counts_df.apply(shannon_entropy, axis=1)

    result = pd.concat(
        [pileup[['chrom', 'pos']].reset_index(drop=True),
         counts_df.reset_index(drop=True)],
        axis=1,
    )

    # 5) 사람이 확인하기 좋은 base count + entropy TSV.
    result.to_csv(counts_tsv, sep='\t', index=False)

    # 6) bedgraph: chrom, start(0-based), end, entropy.
    bg = pd.DataFrame({
        'chrom': result['chrom'],
        'start': result['pos'] - 1,
        'end': result['pos'],
        'entropy': result['entropy'].round(6),
    })
    header = (
        f'track type=bedGraph name="CLIP-{tag}-entropy" '
        f'description="Shannon entropy of CLIP-seq base calls at {name}" '
        f'visibility=full color=40,40,200 altColor=200,40,40 '
        f'autoScale=on viewLimits=0:2 priority=20\n'
    )
    with open(bedgraph, 'w') as fh:
        fh.write(header)
        bg.to_csv(fh, sep='\t', index=False, header=False)

    # 7) 콘솔 요약.
    print(f'positions analyzed : {len(result)}')
    print(f'mean total depth   : {result["total"].mean():.1f}')
    print(f'mean entropy       : {result["entropy"].mean():.4f}')
    idxmax = result['entropy'].idxmax()
    print(f'max  entropy       : {result["entropy"].max():.4f} '
          f'@ {chrom}:{int(result.loc[idxmax, "pos"])}')
    print('top 10 positions by entropy:')
    top = result.sort_values('entropy', ascending=False).head(10)
    print(top[['chrom', 'pos', 'A', 'C', 'G', 'T', 'total', 'entropy']].to_string(index=False))
    print(f'wrote: {counts_tsv.relative_to(W3_DIR)}')
    print(f'wrote: {bedgraph.relative_to(W3_DIR)}')


def main() -> None:
    for gene in GENES:
        process_gene(gene)


if __name__ == '__main__':
    main()
