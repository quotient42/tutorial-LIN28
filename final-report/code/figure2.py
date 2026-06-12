#!/usr/bin/env python3
"""Build Figure 2 for the final LIN28A report.

The script consolidates motif-response and positional-context plotting code from
the self-mission analyses and writes the final three-panel manuscript figure.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec


SCRIPT_DIR = Path(__file__).resolve().parent
REPORT_DIR = SCRIPT_DIR.parent
DATA_DIR = REPORT_DIR / "data"

MOTIF_TABLE_PATH = DATA_DIR / "figure_gene_points.tsv"
W3_SUMMARY_PATH = DATA_DIR / "w3-summary.json"
POSITIONAL_SUMMARY_PATH = DATA_DIR / "figure_positional_region_fractions.tsv"
POSITIONAL_METAPLOT_PATH = DATA_DIR / "figure_positional_metaplot_density.tsv"
FIGURE_PATH = REPORT_DIR / "figures" / "results-figure-2-motif-response-and-position.png"
SVG_PATH = FIGURE_PATH.with_suffix(".svg")

MOTIF_METRIC = "aag_like_per_kb"
MOTIF_LABEL = "AAGNNG + AAGNG motifs per kb"
RNG_SEED = 20260612

FONT_SIZE = 7
LEGEND_SIZE = 6
PANEL_SIZE = 8

ALL_GROUP_ORDER = ["other", "low_clip_control", "high_clip_target", "functional_target"]
DISPLAY_GROUP_ORDER = ["low_clip_control", "high_clip_target", "functional_target"]
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


def add_panel_label(fig: plt.Figure, ax: plt.Axes, label: str, x_offset: float = 0.025) -> None:
    box = ax.get_position()
    fig.text(
        box.x0 - x_offset,
        box.y1 + 0.018,
        label,
        fontsize=PANEL_SIZE,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def pearsonr(x: pd.Series, y: pd.Series) -> float:
    clean = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) < 2:
        return float("nan")
    return float(np.corrcoef(clean.iloc[:, 0], clean.iloc[:, 1])[0, 1])


def spearmanr(x: pd.Series, y: pd.Series) -> float:
    clean = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) < 2:
        return float("nan")
    ranks = clean.rank(method="average")
    return pearsonr(ranks.iloc[:, 0], ranks.iloc[:, 1])


def clean_motif_table(df: pd.DataFrame) -> pd.DataFrame:
    required = ["group", MOTIF_METRIC, "log2_clip", "log2_rden"]
    return df.replace([np.inf, -np.inf], np.nan).dropna(subset=required).copy()


def scatter_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    y_col: str,
    y_label: str,
    rho: float,
) -> None:
    for group in ALL_GROUP_ORDER:
        sub = df[df["group"] == group]
        ax.scatter(
            sub[MOTIF_METRIC],
            sub[y_col],
            s=4 if group == "other" else 8,
            c=GROUP_PALETTE[group],
            alpha=0.16 if group == "other" else 0.44,
            edgecolors="none",
            label=GROUP_LABELS[group],
        )
    ax.axhline(0, color="black", lw=0.55, alpha=0.65)
    ax.set_xlabel(MOTIF_LABEL)
    ax.set_ylabel(y_label)
    ax.text(
        0.98,
        0.04,
        f"Spearman rho = {rho:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )
    setup_axes(ax)


def plot_response_panel(fig: plt.Figure, outer_spec, df: pd.DataFrame, summary: dict[str, object]) -> list[plt.Axes]:
    subgrid = outer_spec.subgridspec(1, 2, wspace=0.28)
    axes = [fig.add_subplot(subgrid[0, 0]), fig.add_subplot(subgrid[0, 1])]
    scatter_panel(
        axes[0],
        df,
        "log2_clip",
        r"LIN28A CLIP enrichment (log$_2$)",
        summary["spearman_r_motif_log2_clip"],
    )
    scatter_panel(
        axes[1],
        df,
        "log2_rden",
        "Ribosome density change after\n" + r"$\it{Lin28a}$ knockdown (log$_2$)",
        summary["spearman_r_motif_log2_rden"],
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.982),
        ncol=4,
        frameon=False,
        handletextpad=0.4,
        columnspacing=1.2,
    )
    return axes


def plot_region_boxplots(fig: plt.Figure, outer_spec, positional: pd.DataFrame) -> list[plt.Axes]:
    subgrid = outer_spec.subgridspec(1, 2, wspace=0.16)
    axes = [fig.add_subplot(subgrid[0, 0]), fig.add_subplot(subgrid[0, 1])]
    features = [
        ("fraction_aag_like_cds", "CDS fraction"),
        ("fraction_aag_like_3utr", "3UTR fraction"),
    ]
    labels = ["Low-CLIP\ncontrols", "High-CLIP\nnonresponse", "Functional\ntargets"]
    for ax, (feature, feature_label) in zip(axes, features):
        data = [
            positional.loc[positional["group"] == group, feature].dropna()
            for group in DISPLAY_GROUP_ORDER
        ]
        bp = ax.boxplot(
            data,
            patch_artist=True,
            tick_labels=labels,
            showfliers=False,
            medianprops={"color": "black", "lw": 0.7},
            boxprops={"lw": 0.6},
            whiskerprops={"lw": 0.6},
            capprops={"lw": 0.6},
        )
        for patch, group in zip(bp["boxes"], DISPLAY_GROUP_ORDER):
            patch.set_facecolor(GROUP_PALETTE[group])
            patch.set_alpha(0.55)
        ax.set_xlabel(feature_label)
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel("Fraction of AAG family sites" if feature == "fraction_aag_like_cds" else "")
        setup_axes(ax)
    return axes


def plot_metaplot(ax: plt.Axes, metaplot: pd.DataFrame) -> None:
    for group in DISPLAY_GROUP_ORDER:
        sub = metaplot[metaplot["group"] == group]
        ax.plot(
            sub["bin_center"],
            sub["motifs_per_gene_per_bin"],
            color=GROUP_PALETTE[group],
            lw=0.95,
            label=GROUP_LABELS[group],
        )
    ax.set_xlabel("Relative position in selected transcript")
    ax.set_ylabel("AAG family motifs per gene per bin")
    ax.legend(frameon=False, loc="upper right")
    setup_axes(ax)


def main() -> None:
    set_style()
    motif_df = clean_motif_table(pd.read_csv(MOTIF_TABLE_PATH, sep="\t"))
    positional = pd.read_csv(POSITIONAL_SUMMARY_PATH, sep="\t").replace([np.inf, -np.inf], np.nan)
    metaplot = pd.read_csv(POSITIONAL_METAPLOT_PATH, sep="\t")
    with W3_SUMMARY_PATH.open() as fh:
        summary = json.load(fh)

    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.1, 5.7), constrained_layout=False)
    grid = GridSpec(
        2,
        2,
        figure=fig,
        height_ratios=[1.0, 0.86],
        width_ratios=[2.0, 1.0],
        left=0.075,
        right=0.985,
        bottom=0.095,
        top=0.89,
        hspace=0.42,
        wspace=0.30,
    )
    axes_a = plot_response_panel(fig, grid[0, :], motif_df, summary)
    axes_b = plot_region_boxplots(fig, grid[1, 0], positional)
    ax_c = fig.add_subplot(grid[1, 1])
    plot_metaplot(ax_c, metaplot)

    add_panel_label(fig, axes_a[0], "A", x_offset=0.04)
    add_panel_label(fig, axes_b[0], "B", x_offset=0.04)
    add_panel_label(fig, ax_c, "C", x_offset=0.04)
    fig.savefig(FIGURE_PATH)
    fig.savefig(SVG_PATH)
    plt.close(fig)
    print(f"Wrote {FIGURE_PATH}")
    print(f"Wrote {SVG_PATH}")


if __name__ == "__main__":
    main()
