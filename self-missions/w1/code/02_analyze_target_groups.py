#!/usr/bin/env python3
"""
Week 1 step 2: define LIN28A target/control groups.

Reads:
  - ../subdata/gene_metrics_clean.tsv

Writes:
  - ../subdata/gene_metrics_grouped.tsv
  - ../subdata/target_group_summary.json
  - ../output/target_gene_lists.tsv
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
W1_DIR = SCRIPT_DIR.parent
SUBDATA_DIR = W1_DIR / "subdata"
OUTPUT_DIR = W1_DIR / "output"

METRICS_PATH = SUBDATA_DIR / "gene_metrics_clean.tsv"
GROUPED_METRICS_PATH = SUBDATA_DIR / "gene_metrics_grouped.tsv"
SUMMARY_PATH = SUBDATA_DIR / "target_group_summary.json"
TARGET_TABLE_PATH = OUTPUT_DIR / "target_gene_lists.tsv"

COUNT_COLUMNS = [
    "CLIP-35L33G.bam",
    "RNA-control.bam",
    "RNA-siLin28a.bam",
    "RNA-siLuc.bam",
    "RPF-siLin28a.bam",
    "RPF-siLuc.bam",
]


def assign_groups(metrics: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    high_clip_cutoff = float(metrics["log2_clip"].quantile(0.90))
    low_clip_cutoff = float(metrics["log2_clip"].quantile(0.50))

    grouped = metrics.copy()
    grouped["is_high_clip"] = grouped["log2_clip"] >= high_clip_cutoff
    grouped["is_low_clip_control"] = grouped["log2_clip"] <= low_clip_cutoff
    grouped["is_functional_target"] = grouped["is_high_clip"] & (grouped["log2_rden"] > 0)

    grouped["group"] = "other"
    grouped.loc[grouped["is_low_clip_control"], "group"] = "low_clip_control"
    grouped.loc[grouped["is_high_clip"], "group"] = "high_clip_target"
    grouped.loc[grouped["is_functional_target"], "group"] = "functional_target"

    grouped["target_score"] = grouped["log2_clip"] + grouped["log2_rden"].clip(lower=0)
    grouped = grouped.sort_values(
        ["group", "target_score", "log2_clip"],
        ascending=[True, False, False],
    )

    r = float(np.corrcoef(grouped["log2_clip"], grouped["log2_rden"])[0, 1])
    summary = {
        "n_filtered_protein_coding_genes": int(len(grouped)),
        "high_clip_quantile": 0.90,
        "low_clip_quantile": 0.50,
        "high_clip_cutoff_log2": high_clip_cutoff,
        "low_clip_cutoff_log2": low_clip_cutoff,
        "pearson_r_log2_clip_log2_rden": r,
        "group_counts": {
            key: int(value) for key, value in grouped["group"].value_counts().to_dict().items()
        },
    }
    return grouped, summary


def make_target_table(grouped: pd.DataFrame) -> pd.DataFrame:
    keep = grouped[grouped["group"] != "other"].copy()
    columns = [
        "gene_id",
        "gene_id_versioned",
        "gene_name",
        "gene_type",
        "group",
        "log2_clip",
        "log2_rden",
        "clip_enrichment",
        "rden_change",
        "target_score",
        "Length",
        *COUNT_COLUMNS,
        "chrom",
        "start",
        "end",
        "strand",
        "n_transcripts",
    ]
    return keep[columns].sort_values(["group", "target_score"], ascending=[True, False])


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SUBDATA_DIR.mkdir(parents=True, exist_ok=True)

    metrics = pd.read_csv(METRICS_PATH, sep="\t")
    grouped, summary = assign_groups(metrics)

    grouped.to_csv(GROUPED_METRICS_PATH, sep="\t", index=False)
    with SUMMARY_PATH.open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    target_table = make_target_table(grouped)
    target_table.to_csv(TARGET_TABLE_PATH, sep="\t", index=False)

    print(f"Wrote {GROUPED_METRICS_PATH.relative_to(W1_DIR)}")
    print(f"Wrote {SUMMARY_PATH.relative_to(W1_DIR)}")
    print(f"Wrote {TARGET_TABLE_PATH.relative_to(W1_DIR)} ({len(target_table):,} rows)")
    print(pd.Series(summary["group_counts"]).to_string())


if __name__ == "__main__":
    main()
