#!/usr/bin/env python3
"""
Week 1 step 3: visualize target groups and write the report.

Reads:
  - ../subdata/gene_metrics_grouped.tsv
  - ../subdata/target_group_summary.json
  - ../output/target_gene_lists.tsv

Writes:
  - ../output/w1-target-groups.png
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
W1_DIR = SCRIPT_DIR.parent
SELF_DIR = W1_DIR.parent
SUBDATA_DIR = W1_DIR / "subdata"
OUTPUT_DIR = W1_DIR / "output"

COUNTS_PATH = SELF_DIR / "data" / "read-counts.txt"
GTF_PATH = SELF_DIR / "data" / "gencode.gtf"
GROUPED_METRICS_PATH = SUBDATA_DIR / "gene_metrics_grouped.tsv"
SUMMARY_PATH = SUBDATA_DIR / "target_group_summary.json"
TARGET_TABLE_PATH = OUTPUT_DIR / "target_gene_lists.tsv"
PLOT_PATH = OUTPUT_DIR / "w1-target-groups.png"
REPORT_PATH = OUTPUT_DIR / "report.md"
LEGEND_FONTSIZE = 6
LABEL_FONTSIZE = 6


def load_summary() -> dict[str, object]:
    with SUMMARY_PATH.open() as fh:
        return json.load(fh)


def place_functional_area_labels(
    ax,
    labels: pd.DataFrame,
    high_clip_cutoff: float,
) -> None:
    """Place displayed labels in a non-overlapping column in functional-target space."""
    if labels.empty:
        return

    labels = labels.sort_values("log2_rden", ascending=False)
    label_x = min(8.0, max(labels["log2_clip"].max() + 0.35, high_clip_cutoff + 2.1))
    label_step = 0.33
    label_top = 4.15
    y_positions = label_top - label_step * np.arange(len(labels))

    for (_, row), label_y in zip(labels.iterrows(), y_positions):
        ax.annotate(
            str(row["gene_name"]),
            xy=(row["log2_clip"], row["log2_rden"]),
            xytext=(label_x, label_y),
            textcoords="data",
            fontsize=LABEL_FONTSIZE,
            ha="left",
            va="center",
            color="#7f0000",
            arrowprops={
                "arrowstyle": "-",
                "color": "#7f0000",
                "lw": 0.35,
                "alpha": 0.72,
                "shrinkA": 2,
                "shrinkB": 2,
                "relpos": (0.0, 0.5),
                "connectionstyle": "arc3,rad=0",
            },
        )


def make_plot(metrics: pd.DataFrame, summary: dict[str, object]) -> Path:
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

    fig, ax = plt.subplots(figsize=(5.0, 4.5))
    for group in order:
        sub = metrics[metrics["group"] == group]
        ax.scatter(
            sub["log2_clip"],
            sub["log2_rden"],
            s=7 if group == "other" else 16,
            c=palette[group],
            alpha=0.25 if group == "other" else 0.5,
            edgecolors="none",
            label=f"{labels[group]} (n={len(sub)})",
            zorder=1 if group == "other" else 2,
        )

    ax.axhline(0, color="black", lw=0.7, alpha=0.65)
    ax.axvline(
        summary["high_clip_cutoff_log2"],
        color="#d73027",
        lw=0.8,
        ls="--",
        alpha=0.8,
    )
    ax.axvline(
        summary["low_clip_cutoff_log2"],
        color="#3b6ea8",
        lw=0.8,
        ls="--",
        alpha=0.8,
    )
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)
    ax.set_axisbelow(True)
    ax.set_xlim(-10, 10)
    ax.set_ylim(-4.5, 4.5)
    ax.set_xlabel(r"LIN28A CLIP enrichment (log$_2$)")
    ax.set_ylabel("Ribosome density change after\n" + r"$\it{Lin28a}$ knockdown (log$_2$)")
    ax.text(
        0.98,
        0.04,
        f"Pearson r = {summary['pearson_r_log2_clip_log2_rden']:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=LEGEND_FONTSIZE,
    )

    top = metrics[metrics["group"] == "functional_target"].copy()
    place_functional_area_labels(
        ax,
        top.nlargest(6, "target_score"),
        summary["high_clip_cutoff_log2"],
    )

    ax.legend(loc="upper left", frameon=False, fontsize=LEGEND_FONTSIZE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(PLOT_PATH, dpi=220)
    plt.close(fig)
    return PLOT_PATH


def write_report(target_table: pd.DataFrame, summary: dict[str, object], plot_path: Path) -> Path:
    group_counts = summary["group_counts"]
    functional = target_table[target_table["group"] == "functional_target"]

    top_cols = ["gene_name", "gene_id", "log2_clip", "log2_rden", "target_score"]
    top_functional = functional.nlargest(10, "target_score")[top_cols].copy()
    for col in ["log2_clip", "log2_rden", "target_score"]:
        top_functional[col] = top_functional[col].map(lambda x: f"{x:.3f}")

    lines = [
        "# Week 1 Report: LIN28A Target Groups",
        "",
        "## Goal",
        "",
        "This week defines target and control gene groups for the LIN28A motif exploration project.",
        "The analysis follows the tutorial formulas for LIN28A CLIP enrichment and ribosome-density change after Lin28a knockdown.",
        "",
        "## Inputs",
        "",
        f"- Count matrix: `{COUNTS_PATH.relative_to(SELF_DIR)}`",
        f"- Gene annotation: `{GTF_PATH.relative_to(SELF_DIR)}`",
        "",
        "## Filtering and Metrics",
        "",
        "- `clip_enrichment = CLIP-35L33G / RNA-control`",
        "- `rden_change = (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc)`",
        "- Genes were kept only if the relevant denominators were nonzero and RNA counts were at least 10 in the three RNA libraries.",
        "- Target/control groups were defined only among protein-coding genes, because Week 2 will scan mRNA transcript sequences.",
        f"- Final filtered protein-coding genes: {summary['n_filtered_protein_coding_genes']:,}",
        f"- Pearson correlation between log2 CLIP enrichment and log2 ribosome-density change: {summary['pearson_r_log2_clip_log2_rden']:.4f}",
        "",
        "## Group Definitions",
        "",
        f"- High-CLIP cutoff: top 10% of filtered genes, log2 CLIP >= {summary['high_clip_cutoff_log2']:.3f}",
        f"- Low-CLIP control cutoff: bottom 50% of filtered genes, log2 CLIP <= {summary['low_clip_cutoff_log2']:.3f}",
        "- Functional targets: high-CLIP genes with positive log2 ribosome-density change after Lin28a knockdown",
        "",
        "## Group Counts",
        "",
        f"- Functional targets: {group_counts.get('functional_target', 0):,}",
        f"- High-CLIP targets without positive ribosome-density change: {group_counts.get('high_clip_target', 0):,}",
        f"- Low-CLIP controls: {group_counts.get('low_clip_control', 0):,}",
        f"- Other filtered genes: {group_counts.get('other', 0):,}",
        "",
        "## Outputs",
        "",
        "- Cleaned metrics table: `subdata/gene_metrics_clean.tsv`",
        "- Grouped metrics table: `subdata/gene_metrics_grouped.tsv`",
        "- Gene annotation table: `subdata/gene_annotation.tsv`",
        "- Target/control table: `output/target_gene_lists.tsv`",
        f"- Figure: `{plot_path.relative_to(W1_DIR)}`",
        "",
        "## Top Functional Targets",
        "",
        top_functional.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "Genes with high LIN28A CLIP enrichment and positive ribosome-density change after Lin28a knockdown are plausible direct translational repression targets.",
        "They are strong candidates for Week 2 motif counting because they combine evidence of LIN28A binding with a functional translation response.",
        "Low-CLIP controls provide a comparison group for testing whether LIN28A recognition motifs are enriched among candidate targets.",
        "",
        "## Caveats",
        "",
        "- These are gene-level counts, so isoform-specific effects are not resolved.",
        "- Ratio-based metrics are sensitive to low counts; the current filter is intentionally conservative for a first-pass target list.",
        "- High CLIP enrichment does not prove direct functional repression unless it is paired with a translation response or further validation.",
    ]

    REPORT_PATH.write_text("\n".join(lines) + "\n")
    return REPORT_PATH


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    metrics = pd.read_csv(GROUPED_METRICS_PATH, sep="\t")
    target_table = pd.read_csv(TARGET_TABLE_PATH, sep="\t")
    summary = load_summary()

    plot_path = make_plot(metrics, summary)
    report_path = write_report(target_table, summary, plot_path)

    print(f"Wrote {plot_path.relative_to(W1_DIR)}")
    print(f"Wrote {report_path.relative_to(W1_DIR)}")


if __name__ == "__main__":
    main()
