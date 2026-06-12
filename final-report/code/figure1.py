#!/usr/bin/env python3
"""Build Figure 1 for the final LIN28A report.

The script consolidates plotting code from the self-mission analyses and writes
the final two-panel figure used by the manuscript.
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
REPORT_DIR = SCRIPT_DIR.parent
DATA_DIR = REPORT_DIR / "data"

MOTIF_TABLE_PATH = DATA_DIR / "figure_gene_points.tsv"
SUMMARY_PATH = DATA_DIR / "w3-summary.json"
FIGURE_PATH = REPORT_DIR / "figures" / "results-figure-1-targets-and-motif-burden.png"
SVG_PATH = FIGURE_PATH.with_suffix(".svg")

MOTIF_METRIC = "aag_like_per_kb"
MOTIF_LABEL = "AAGNNG + AAGNG motifs per kb"
RNG_SEED = 20260612

FONT_SIZE = 7
LEGEND_SIZE = 6
PANEL_SIZE = 8

GROUP_PALETTE = {
    "other": "#bdbdbd",
    "low_clip_control": "#3b6ea8",
    "high_clip_target": "#f0a202",
    "functional_target": "#d73027",
}
GROUP_LABELS = {
    "other": "Other filtered genes",
    "low_clip_control": "Low-CLIP controls",
    "high_clip_target": "High-CLIP nonresponse",
    "functional_target": "Functional targets",
}
GROUP_ORDER = ["other", "low_clip_control", "high_clip_target", "functional_target"]


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.size": FONT_SIZE,
            "axes.labelsize": FONT_SIZE,
            "axes.titlesize": FONT_SIZE,
            "xtick.labelsize": FONT_SIZE,
            "ytick.labelsize": FONT_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "figure.dpi": 220,
            "savefig.dpi": 220,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
        }
    )


def setup_axes(ax: plt.Axes) -> None:
    ax.grid(True, linestyle=":", linewidth=0.35, alpha=0.45)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def clean_table(df: pd.DataFrame) -> pd.DataFrame:
    required = ["group", "log2_clip", "log2_rden", MOTIF_METRIC]
    return df.replace([np.inf, -np.inf], np.nan).dropna(subset=required).copy()


def add_panel_label(fig: plt.Figure, ax: plt.Axes, label: str) -> None:
    box = ax.get_position()
    fig.text(
        box.x0 - 0.025,
        box.y1 + 0.018,
        label,
        fontsize=PANEL_SIZE,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def plot_target_scatter(ax: plt.Axes, df: pd.DataFrame) -> None:
    for group in GROUP_ORDER:
        sub = df[df["group"] == group]
        ax.scatter(
            sub["log2_clip"],
            sub["log2_rden"],
            s=5 if group == "other" else 10,
            c=GROUP_PALETTE[group],
            alpha=0.18 if group == "other" else 0.48,
            edgecolors="none",
            label=f"{GROUP_LABELS[group]} (n={len(sub)})",
            zorder=1 if group == "other" else 2,
        )

    ax.axhline(0, color="black", lw=0.55, alpha=0.65)
    ax.axvline(0, color="black", lw=0.55, alpha=0.65)
    ax.set_xlim(-10, 10)
    ax.set_ylim(-4.8, 4.8)
    ax.set_xlabel(r"LIN28A CLIP enrichment (log$_2$)")
    ax.set_ylabel("Ribosome density change after\n" + r"$\it{Lin28a}$ knockdown (log$_2$)")
    ax.legend(loc="upper left", frameon=False, markerscale=0.9, handletextpad=0.4)
    setup_axes(ax)


def plot_motif_boxplot(ax: plt.Axes, df: pd.DataFrame, summary: dict[str, object]) -> None:
    groups = [
        ("low_clip_control", "Low-CLIP\ncontrols"),
        ("high_clip_target", "High-CLIP\nnonresponse"),
        ("functional_target", "Functional\ntargets"),
    ]
    data = [df.loc[df["group"] == group, MOTIF_METRIC] for group, _ in groups]
    bp = ax.boxplot(
        data,
        patch_artist=True,
        tick_labels=[label for _, label in groups],
        widths=0.55,
        showfliers=False,
        medianprops={"color": "black", "lw": 0.7},
        boxprops={"lw": 0.6},
        whiskerprops={"lw": 0.6},
        capprops={"lw": 0.6},
    )
    for patch, (group, _) in zip(bp["boxes"], groups):
        patch.set_facecolor(GROUP_PALETTE[group])
        patch.set_alpha(0.55)

    rng = np.random.default_rng(RNG_SEED)
    for i, (group, _) in enumerate(groups, start=1):
        values = df.loc[df["group"] == group, MOTIF_METRIC].to_numpy()
        if len(values) > 700:
            values = rng.choice(values, size=700, replace=False)
        x = rng.normal(i, 0.045, size=len(values))
        ax.scatter(
            x,
            values,
            s=3,
            alpha=0.12,
            c=GROUP_PALETTE[group],
            edgecolors="none",
            zorder=0,
        )

    test = summary["median_test_high_clip_vs_low_clip"]
    ax.text(
        0.98,
        0.96,
        "High vs control\n"
        + f"median diff = {test['observed_median_difference_high_minus_control']:.3f}\n"
        + f"p = {test['p_value_high_greater_than_control']:.4f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
    )
    ax.set_ylabel(MOTIF_LABEL)
    setup_axes(ax)


def main() -> None:
    set_style()
    df = clean_table(pd.read_csv(MOTIF_TABLE_PATH, sep="\t"))
    with SUMMARY_PATH.open() as fh:
        summary = json.load(fh)

    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.2), constrained_layout=False)
    plot_target_scatter(axes[0], df)
    plot_motif_boxplot(axes[1], df, summary)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.18, top=0.92, wspace=0.34)
    add_panel_label(fig, axes[0], "A")
    add_panel_label(fig, axes[1], "B")
    fig.savefig(FIGURE_PATH)
    fig.savefig(SVG_PATH)
    plt.close(fig)
    print(f"Wrote {FIGURE_PATH}")
    print(f"Wrote {SVG_PATH}")


if __name__ == "__main__":
    main()
