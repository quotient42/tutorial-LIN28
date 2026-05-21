#!/usr/bin/env python3
"""
Week 1 step 1: process raw count and annotation data.

Writes intermediate tables to ../subdata:
  - gene_annotation.tsv
  - gene_metrics_clean.tsv
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
W1_DIR = SCRIPT_DIR.parent
SELF_DIR = W1_DIR.parent
DATA_DIR = SELF_DIR / "data"
SUBDATA_DIR = W1_DIR / "subdata"

COUNTS_PATH = DATA_DIR / "read-counts.txt"
GTF_PATH = DATA_DIR / "gencode.gtf"
ANNOTATION_PATH = SUBDATA_DIR / "gene_annotation.tsv"
METRICS_PATH = SUBDATA_DIR / "gene_metrics_clean.tsv"

COUNT_COLUMNS = [
    "CLIP-35L33G.bam",
    "RNA-control.bam",
    "RNA-siLin28a.bam",
    "RNA-siLuc.bam",
    "RPF-siLin28a.bam",
    "RPF-siLuc.bam",
]

ATTR_RE = re.compile(r'(\S+) "([^"]+)"')


def stable_id(versioned_id: str) -> str:
    return versioned_id.split(".", 1)[0]


def parse_attrs(attr_text: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attr_text))


def build_gene_annotation(gtf_path: Path) -> pd.DataFrame:
    """Parse gene-level names/types from GTF and transcript counts per gene."""
    genes: dict[str, dict[str, object]] = {}
    transcript_counts: dict[str, int] = {}

    with gtf_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue

            sid = stable_id(gene_id)
            if feature == "gene":
                genes[sid] = {
                    "gene_id": sid,
                    "gene_id_versioned": gene_id,
                    "gene_name": attrs.get("gene_name", sid),
                    "gene_type": attrs.get("gene_type", ""),
                    "chrom": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                }
            elif feature == "transcript":
                transcript_counts[sid] = transcript_counts.get(sid, 0) + 1

    annotation = pd.DataFrame(genes.values())
    if annotation.empty:
        raise RuntimeError(f"No gene annotations parsed from {gtf_path}")

    annotation["n_transcripts"] = (
        annotation["gene_id"].map(transcript_counts).fillna(0).astype(int)
    )
    return annotation.sort_values("gene_id")


def calculate_metrics(counts_path: Path, annotation: pd.DataFrame) -> pd.DataFrame:
    counts = pd.read_csv(counts_path, sep="\t", comment="#")
    missing = sorted(set(COUNT_COLUMNS) - set(counts.columns))
    if missing:
        raise RuntimeError(f"Missing expected count columns: {missing}")

    counts["gene_id_versioned"] = counts["Geneid"]
    counts["gene_id"] = counts["Geneid"].map(stable_id)

    with np.errstate(divide="ignore", invalid="ignore"):
        counts["clip_enrichment"] = counts["CLIP-35L33G.bam"] / counts["RNA-control.bam"]
        counts["rden_change"] = (
            (counts["RPF-siLin28a.bam"] / counts["RNA-siLin28a.bam"])
            / (counts["RPF-siLuc.bam"] / counts["RNA-siLuc.bam"])
        )
        counts["log2_clip"] = np.log2(counts["clip_enrichment"])
        counts["log2_rden"] = np.log2(counts["rden_change"])

    counts["passes_count_filter"] = (
        (counts["CLIP-35L33G.bam"] > 0)
        & (counts["RNA-control.bam"] >= 10)
        & (counts["RNA-siLuc.bam"] >= 10)
        & (counts["RNA-siLin28a.bam"] >= 10)
        & (counts["RPF-siLuc.bam"] > 0)
        & (counts["RPF-siLin28a.bam"] > 0)
    )

    metrics = counts.replace([np.inf, -np.inf], np.nan)
    metrics = metrics.dropna(subset=["log2_clip", "log2_rden"])
    metrics = metrics[metrics["passes_count_filter"]].copy()

    metrics = metrics.merge(
        annotation[
            [
                "gene_id",
                "gene_name",
                "gene_type",
                "chrom",
                "start",
                "end",
                "strand",
                "n_transcripts",
            ]
        ],
        on="gene_id",
        how="left",
    )
    metrics["gene_name"] = metrics["gene_name"].fillna(metrics["gene_id"])
    metrics["gene_type"] = metrics["gene_type"].fillna("")

    # The motif project scans mRNA transcript sequences in Week 2.
    return metrics[metrics["gene_type"] == "protein_coding"].copy()


def main() -> None:
    SUBDATA_DIR.mkdir(parents=True, exist_ok=True)

    annotation = build_gene_annotation(GTF_PATH)
    annotation.to_csv(ANNOTATION_PATH, sep="\t", index=False)

    metrics = calculate_metrics(COUNTS_PATH, annotation)
    metrics.to_csv(METRICS_PATH, sep="\t", index=False)

    print(f"Wrote {ANNOTATION_PATH.relative_to(W1_DIR)} ({len(annotation):,} genes)")
    print(f"Wrote {METRICS_PATH.relative_to(W1_DIR)} ({len(metrics):,} protein-coding genes)")


if __name__ == "__main__":
    main()
