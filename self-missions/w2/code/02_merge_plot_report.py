#!/usr/bin/env python3
"""
Week 2 step 2: merge motif counts with Week 1 target groups and plot results.

Reads:
  - ../subdata/transcript_motif_counts_raw.tsv
  - ../../w1/subdata/gene_metrics_grouped.tsv

Writes:
  - ../output/transcript_motif_counts.tsv
  - ../output/w2-motif-burden-boxplot.png
  - ../output/w2-motif-vs-clip-scatter.png
  - ../output/report.md
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
W2_DIR = SCRIPT_DIR.parent
SELF_DIR = W2_DIR.parent
W1_DIR = SELF_DIR / "w1"
SUBDATA_DIR = W2_DIR / "subdata"
OUTPUT_DIR = W2_DIR / "output"

RAW_MOTIF_PATH = SUBDATA_DIR / "transcript_motif_counts_raw.tsv"
W1_GROUPED_PATH = W1_DIR / "subdata" / "gene_metrics_grouped.tsv"
MERGED_MOTIF_PATH = OUTPUT_DIR / "transcript_motif_counts.tsv"
SUMMARY_PATH = SUBDATA_DIR / "motif_summary.json"
BOXPLOT_PATH = OUTPUT_DIR / "w2-motif-burden-boxplot.png"
SCATTER_PATH = OUTPUT_DIR / "w2-motif-vs-clip-scatter.png"
REPORT_PATH = OUTPUT_DIR / "report.md"

MOTIF_METRIC = "aag_like_per_kb"
MOTIF_LABEL = "AAGNNG + AAGNG motifs per kb"


def pearsonr(x: pd.Series, y: pd.Series) -> float:
    clean = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) < 2:
        return float("nan")
    return float(np.corrcoef(clean.iloc[:, 0], clean.iloc[:, 1])[0, 1])


def merge_tables(metrics: pd.DataFrame, motifs: pd.DataFrame) -> pd.DataFrame:
    metric_columns = [
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
        "CLIP-35L33G.bam",
        "RNA-control.bam",
        "RNA-siLin28a.bam",
        "RNA-siLuc.bam",
        "RPF-siLin28a.bam",
        "RPF-siLuc.bam",
    ]
    merged = metrics[metric_columns].merge(
        motifs,
        on=["gene_id", "gene_id_versioned", "gene_name"],
        how="inner",
    )
    return merged.sort_values(["group", "log2_clip"], ascending=[True, False])


def summarize(merged: pd.DataFrame) -> dict[str, object]:
    high = merged[merged["group"].isin(["functional_target", "high_clip_target"])]
    low = merged[merged["group"] == "low_clip_control"]
    functional = merged[merged["group"] == "functional_target"]

    return {
        "n_genes_with_metrics_and_motifs": int(len(merged)),
        "motif_metric": MOTIF_METRIC,
        "pearson_r_motif_log2_clip": pearsonr(merged[MOTIF_METRIC], merged["log2_clip"]),
        "pearson_r_motif_log2_rden": pearsonr(merged[MOTIF_METRIC], merged["log2_rden"]),
        "median_motif_per_kb": {
            "high_clip_targets": float(high[MOTIF_METRIC].median()),
            "functional_targets": float(functional[MOTIF_METRIC].median()),
            "low_clip_controls": float(low[MOTIF_METRIC].median()),
            "all_filtered_genes": float(merged[MOTIF_METRIC].median()),
        },
        "group_counts": {
            key: int(value) for key, value in merged["group"].value_counts().to_dict().items()
        },
    }


def make_boxplot(merged: pd.DataFrame) -> Path:
    groups = [
        ("low_clip_control", "Low-CLIP\ncontrols", "#3b6ea8"),
        ("high_clip_target", "High-CLIP\ntargets", "#f0a202"),
        ("functional_target", "Functional\ntargets", "#d73027"),
    ]
    data = [merged.loc[merged["group"] == group, MOTIF_METRIC] for group, _, _ in groups]

    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    bp = ax.boxplot(
        data,
        patch_artist=True,
        tick_labels=[label for _, label, _ in groups],
        widths=0.55,
        showfliers=False,
        medianprops={"color": "black", "lw": 1.0},
        boxprops={"lw": 0.8},
        whiskerprops={"lw": 0.8},
        capprops={"lw": 0.8},
    )
    for patch, (_, _, color) in zip(bp["boxes"], groups):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)

    rng = np.random.default_rng(20260527)
    for i, (group, _, color) in enumerate(groups, start=1):
        values = merged.loc[merged["group"] == group, MOTIF_METRIC].to_numpy()
        if len(values) > 700:
            values = rng.choice(values, size=700, replace=False)
        x = rng.normal(i, 0.045, size=len(values))
        ax.scatter(x, values, s=5, alpha=0.14, c=color, edgecolors="none")

    ax.set_ylabel(MOTIF_LABEL)
    ax.set_title("LIN28A motif burden by target group")
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(BOXPLOT_PATH, dpi=220)
    plt.close(fig)
    return BOXPLOT_PATH


def make_scatter(merged: pd.DataFrame, summary: dict[str, object]) -> Path:
    palette = {
        "other": "#bdbdbd",
        "low_clip_control": "#3b6ea8",
        "high_clip_target": "#f0a202",
        "functional_target": "#d73027",
    }
    labels = {
        "other": "Other filtered genes",
        "low_clip_control": "Low-CLIP controls",
        "high_clip_target": "High-CLIP targets",
        "functional_target": "Functional targets",
    }
    order = ["other", "low_clip_control", "high_clip_target", "functional_target"]

    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for group in order:
        sub = merged[merged["group"] == group]
        ax.scatter(
            sub[MOTIF_METRIC],
            sub["log2_clip"],
            s=7 if group == "other" else 14,
            c=palette[group],
            alpha=0.22 if group == "other" else 0.50,
            edgecolors="none",
            label=f"{labels[group]} (n={len(sub)})",
        )

    ax.set_xlabel(MOTIF_LABEL)
    ax.set_ylabel(r"LIN28A CLIP enrichment (log$_2$)")
    ax.text(
        0.98,
        0.04,
        f"Pearson r = {summary['pearson_r_motif_log2_clip']:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
    )
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)
    ax.legend(loc="upper right", frameon=False, fontsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(SCATTER_PATH, dpi=220)
    plt.close(fig)
    return SCATTER_PATH


def markdown_table(df: pd.DataFrame) -> str:
    columns = list(df.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in columns) + " |")
    return "\n".join(lines)


def write_report(merged: pd.DataFrame, summary: dict[str, object]) -> Path:
    top_cols = [
        "gene_name",
        "group",
        "log2_clip",
        "log2_rden",
        "aag_like_per_kb",
        "aagnng_per_kb",
        "aagng_per_kb",
        "tgtg_per_kb",
    ]
    top = merged.nlargest(10, "log2_clip")[top_cols].copy()
    for col in ["log2_clip", "log2_rden", "aag_like_per_kb", "aagnng_per_kb", "aagng_per_kb", "tgtg_per_kb"]:
        top[col] = top[col].map(lambda x: f"{x:.3f}")

    med = summary["median_motif_per_kb"]
    lines = [
        "# Week 2 Report: LIN28A Motif Burden",
        "",
        "## Goal",
        "",
        "This week asks whether simple LIN28A sequence motifs explain the target groups defined in Week 1.",
        "The scan uses DNA equivalents of the reported RNA motifs: `AAGNNG`, `AAGNG`, and `TGTG` for `UGUG`.",
        "",
        "## Inputs",
        "",
        "- Transcript annotation: `data/gencode.gtf`",
        "- Transcript sequences: `data/gencode.vM27.transcripts.fa.gz`",
        "- Week 1 grouped metrics: `w1/subdata/gene_metrics_grouped.tsv`",
        "",
        "## Transcript Choice",
        "",
        "One protein-coding transcript was selected per protein-coding gene.",
        "The priority order was lowest transcript support level, APPRIS principal tag when available, then longest exon-composed transcript length.",
        "",
        "## Results",
        "",
        f"- Genes with both Week 1 metrics and selected transcript motifs: {summary['n_genes_with_metrics_and_motifs']:,}",
        f"- Pearson correlation between `{MOTIF_METRIC}` and log2 CLIP enrichment: {summary['pearson_r_motif_log2_clip']:.4f}",
        f"- Pearson correlation between `{MOTIF_METRIC}` and log2 ribosome-density change: {summary['pearson_r_motif_log2_rden']:.4f}",
        f"- Median `{MOTIF_METRIC}` in high-CLIP targets: {med['high_clip_targets']:.3f}",
        f"- Median `{MOTIF_METRIC}` in functional targets: {med['functional_targets']:.3f}",
        f"- Median `{MOTIF_METRIC}` in low-CLIP controls: {med['low_clip_controls']:.3f}",
        "",
        "## Outputs",
        "",
        "- Raw motif table: `subdata/transcript_motif_counts_raw.tsv`",
        "- Merged motif and target table: `output/transcript_motif_counts.tsv`",
        "- Boxplot: `output/w2-motif-burden-boxplot.png`",
        "- Scatterplot: `output/w2-motif-vs-clip-scatter.png`",
        "",
        "## Top CLIP-Enriched Genes",
        "",
        markdown_table(top),
        "",
        "## Interpretation",
        "",
        "A weak motif-vs-CLIP correlation would mean that the sequence motifs contribute to LIN28A binding but are not sufficient to explain target selection.",
        "That is the outcome expected from Cho et al. 2012: LIN28A recognizes AAG-rich motifs, yet ER-proximal localization and transcript context help determine which mRNAs become strongly bound and translationally repressed.",
        "",
        "## Caveats",
        "",
        "- This scan counts linear sequence motifs only; it does not test the hairpin-loop structure reported for the main LIN28A motif.",
        "- One transcript per gene hides isoform-specific motif differences.",
        "- Motif counts are transcript-level, while Week 1 CLIP and Ribo-seq metrics are gene-level.",
    ]

    REPORT_PATH.write_text("\n".join(lines) + "\n")
    return REPORT_PATH


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SUBDATA_DIR.mkdir(parents=True, exist_ok=True)

    motifs = pd.read_csv(RAW_MOTIF_PATH, sep="\t")
    metrics = pd.read_csv(W1_GROUPED_PATH, sep="\t")
    merged = merge_tables(metrics, motifs)
    merged.to_csv(MERGED_MOTIF_PATH, sep="\t", index=False)

    summary = summarize(merged)
    with SUMMARY_PATH.open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    boxplot_path = make_boxplot(merged)
    scatter_path = make_scatter(merged, summary)
    report_path = write_report(merged, summary)

    print(f"Wrote {MERGED_MOTIF_PATH.relative_to(W2_DIR)} ({len(merged):,} genes)")
    print(f"Wrote {SUMMARY_PATH.relative_to(W2_DIR)}")
    print(f"Wrote {boxplot_path.relative_to(W2_DIR)}")
    print(f"Wrote {scatter_path.relative_to(W2_DIR)}")
    print(f"Wrote {report_path.relative_to(W2_DIR)}")


if __name__ == "__main__":
    main()
