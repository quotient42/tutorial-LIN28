#!/usr/bin/env python3
"""
Week 3: interpret whether linear LIN28A motif burden explains binding and
translation response.

Reads:
  - ../../w2/output/transcript_motif_counts.tsv

Writes:
  - ../subdata/w3-summary.json
  - ../output/w3-target-control-scatter.png
  - ../output/w3-motif-burden-boxplot.png
  - ../output/w3-motif-vs-response.png
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
W3_DIR = SCRIPT_DIR.parent
SELF_DIR = W3_DIR.parent
W2_DIR = SELF_DIR / "w2"
SUBDATA_DIR = W3_DIR / "subdata"
OUTPUT_DIR = W3_DIR / "output"

MOTIF_TABLE_PATH = W2_DIR / "output" / "transcript_motif_counts.tsv"
SUMMARY_PATH = SUBDATA_DIR / "w3-summary.json"
TARGET_SCATTER_PATH = OUTPUT_DIR / "w3-target-control-scatter.png"
BOXPLOT_PATH = OUTPUT_DIR / "w3-motif-burden-boxplot.png"
RESPONSE_PATH = OUTPUT_DIR / "w3-motif-vs-response.png"
REPORT_PATH = OUTPUT_DIR / "report.md"

MOTIF_METRIC = "aag_like_per_kb"
MOTIF_LABEL = "AAGNNG + AAGNG motifs per kb"
RNG_SEED = 20260604

GROUP_PALETTE = {
    "other": "#bdbdbd",
    "low_clip_control": "#3b6ea8",
    "high_clip_target": "#f0a202",
    "functional_target": "#d73027",
}
GROUP_LABELS = {
    "other": "Other filtered genes",
    "low_clip_control": "Low-CLIP controls",
    "high_clip_target": "High-CLIP targets",
    "functional_target": "Functional targets",
}
GROUP_ORDER = ["other", "low_clip_control", "high_clip_target", "functional_target"]


def clean_analysis_table(df: pd.DataFrame) -> pd.DataFrame:
    required = ["group", MOTIF_METRIC, "log2_clip", "log2_rden"]
    clean = df.replace([np.inf, -np.inf], np.nan).dropna(subset=required).copy()
    return clean.sort_values(["group", "log2_clip"], ascending=[True, False])


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


def permutation_median_test(
    high_values: np.ndarray,
    control_values: np.ndarray,
    n_permutations: int = 20000,
) -> dict[str, float]:
    """One-sided test for whether high-CLIP targets have higher motif burden."""
    high_values = np.asarray(high_values, dtype=float)
    control_values = np.asarray(control_values, dtype=float)
    observed = float(np.median(high_values) - np.median(control_values))

    pooled = np.concatenate([high_values, control_values])
    n_high = len(high_values)
    rng = np.random.default_rng(RNG_SEED)
    hits = 0
    for _ in range(n_permutations):
        permuted = rng.permutation(pooled)
        diff = float(np.median(permuted[:n_high]) - np.median(permuted[n_high:]))
        if diff >= observed:
            hits += 1

    p_greater = (hits + 1) / (n_permutations + 1)
    return {
        "observed_median_difference_high_minus_control": observed,
        "p_value_high_greater_than_control": float(p_greater),
        "n_permutations": int(n_permutations),
    }


def summarize(df: pd.DataFrame) -> dict[str, object]:
    high = df[df["group"].isin(["functional_target", "high_clip_target"])]
    functional = df[df["group"] == "functional_target"]
    controls = df[df["group"] == "low_clip_control"]

    test = permutation_median_test(high[MOTIF_METRIC].to_numpy(), controls[MOTIF_METRIC].to_numpy())

    return {
        "n_genes": int(len(df)),
        "motif_metric": MOTIF_METRIC,
        "group_counts": {key: int(value) for key, value in df["group"].value_counts().to_dict().items()},
        "spearman_r_motif_log2_clip": spearmanr(df[MOTIF_METRIC], df["log2_clip"]),
        "spearman_r_motif_log2_rden": spearmanr(df[MOTIF_METRIC], df["log2_rden"]),
        "pearson_r_motif_log2_clip": pearsonr(df[MOTIF_METRIC], df["log2_clip"]),
        "pearson_r_motif_log2_rden": pearsonr(df[MOTIF_METRIC], df["log2_rden"]),
        "median_motif_per_kb": {
            "high_clip_targets": float(high[MOTIF_METRIC].median()),
            "functional_targets": float(functional[MOTIF_METRIC].median()),
            "low_clip_controls": float(controls[MOTIF_METRIC].median()),
            "all_filtered_genes": float(df[MOTIF_METRIC].median()),
        },
        "median_test_high_clip_vs_low_clip": test,
    }


def make_target_scatter(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(5.0, 4.5))
    for group in GROUP_ORDER:
        sub = df[df["group"] == group]
        ax.scatter(
            sub["log2_clip"],
            sub["log2_rden"],
            s=7 if group == "other" else 16,
            c=GROUP_PALETTE[group],
            alpha=0.22 if group == "other" else 0.50,
            edgecolors="none",
            label=f"{GROUP_LABELS[group]} (n={len(sub)})",
            zorder=1 if group == "other" else 2,
        )

    ax.axhline(0, color="black", lw=0.7, alpha=0.65)
    ax.axvline(0, color="black", lw=0.7, alpha=0.65)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)
    ax.set_axisbelow(True)
    ax.set_xlim(-10, 10)
    ax.set_ylim(-4.8, 4.8)
    ax.set_xlabel(r"LIN28A CLIP enrichment (log$_2$)")
    ax.set_ylabel("Ribosome density change after\n" + r"$\it{Lin28a}$ knockdown (log$_2$)")
    ax.set_title("Target and control groups")
    ax.legend(loc="upper left", frameon=False, fontsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(TARGET_SCATTER_PATH, dpi=220)
    plt.close(fig)
    return TARGET_SCATTER_PATH


def make_boxplot(df: pd.DataFrame, summary: dict[str, object]) -> Path:
    groups = [
        ("low_clip_control", "Low-CLIP\ncontrols", GROUP_PALETTE["low_clip_control"]),
        ("high_clip_target", "High-CLIP\ntargets", GROUP_PALETTE["high_clip_target"]),
        ("functional_target", "Functional\ntargets", GROUP_PALETTE["functional_target"]),
    ]
    data = [df.loc[df["group"] == group, MOTIF_METRIC] for group, _, _ in groups]

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

    rng = np.random.default_rng(RNG_SEED)
    for i, (group, _, color) in enumerate(groups, start=1):
        values = df.loc[df["group"] == group, MOTIF_METRIC].to_numpy()
        if len(values) > 700:
            values = rng.choice(values, size=700, replace=False)
        x = rng.normal(i, 0.045, size=len(values))
        ax.scatter(x, values, s=5, alpha=0.14, c=color, edgecolors="none")

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
        fontsize=7,
    )
    ax.set_ylabel(MOTIF_LABEL)
    ax.set_title("Motif burden by target group")
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(BOXPLOT_PATH, dpi=220)
    plt.close(fig)
    return BOXPLOT_PATH


def scatter_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    y_col: str,
    y_label: str,
    rho: float,
) -> None:
    for group in GROUP_ORDER:
        sub = df[df["group"] == group]
        ax.scatter(
            sub[MOTIF_METRIC],
            sub[y_col],
            s=6 if group == "other" else 13,
            c=GROUP_PALETTE[group],
            alpha=0.18 if group == "other" else 0.46,
            edgecolors="none",
            label=GROUP_LABELS[group],
        )
    ax.axhline(0, color="black", lw=0.6, alpha=0.55)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)
    ax.set_xlabel(MOTIF_LABEL)
    ax.set_ylabel(y_label)
    ax.text(
        0.98,
        0.04,
        f"Spearman rho = {rho:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_response_plot(df: pd.DataFrame, summary: dict[str, object]) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), sharex=True)
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
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=4,
        frameon=False,
        fontsize=7,
    )
    fig.suptitle("Motif burden versus binding and translation response", y=0.98, fontsize=10)
    fig.tight_layout(rect=(0, 0.10, 1, 0.92))
    fig.savefig(RESPONSE_PATH, dpi=220)
    plt.close(fig)
    return RESPONSE_PATH


def markdown_table(df: pd.DataFrame) -> str:
    columns = list(df.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in columns) + " |")
    return "\n".join(lines)


def write_report(df: pd.DataFrame, summary: dict[str, object]) -> Path:
    med = summary["median_motif_per_kb"]
    test = summary["median_test_high_clip_vs_low_clip"]

    top_cols = ["gene_name", "group", "log2_clip", "log2_rden", MOTIF_METRIC, "tgtg_per_kb"]
    top_functional = df[df["group"] == "functional_target"].nlargest(10, "target_score")[top_cols].copy()
    for col in ["log2_clip", "log2_rden", MOTIF_METRIC, "tgtg_per_kb"]:
        top_functional[col] = top_functional[col].map(lambda x: f"{x:.3f}")

    lines = [
        "# Week 3 Report: Do LIN28A Motifs Explain Binding and Repression?",
        "",
        "## Goal",
        "",
        "Week 3 combines the Week 1 target groups with the Week 2 transcript motif counts.",
        "The central test is whether a simple linear count of AAG-like LIN28A motifs is enough to explain LIN28A CLIP enrichment and translational repression.",
        "",
        "## Inputs",
        "",
        "- Merged Week 2 motif table: `w2/output/transcript_motif_counts.tsv`",
        "- Motif burden metric: `aag_like_per_kb = (AAGNNG + AAGNG motif counts) / transcript kb`",
        "",
        "## Results",
        "",
        f"- Genes with motif burden, CLIP enrichment, and ribosome-density metrics: {summary['n_genes']:,}",
        f"- Spearman correlation between motif burden and log2 CLIP enrichment: {summary['spearman_r_motif_log2_clip']:.4f}",
        f"- Spearman correlation between motif burden and log2 ribosome-density change: {summary['spearman_r_motif_log2_rden']:.4f}",
        f"- Median motif burden in high-CLIP targets: {med['high_clip_targets']:.3f} motifs/kb",
        f"- Median motif burden in functional targets: {med['functional_targets']:.3f} motifs/kb",
        f"- Median motif burden in low-CLIP controls: {med['low_clip_controls']:.3f} motifs/kb",
        f"- One-sided permutation test for high-CLIP targets having higher motif burden than low-CLIP controls: median difference = {test['observed_median_difference_high_minus_control']:.3f} motifs/kb, p = {test['p_value_high_greater_than_control']:.4f}",
        "",
        "## Outputs",
        "",
        "- Summary statistics: `subdata/w3-summary.json`",
        "- Target/control scatterplot: `output/w3-target-control-scatter.png`",
        "- Motif burden boxplot: `output/w3-motif-burden-boxplot.png`",
        "- Motif-versus-response scatterplots: `output/w3-motif-vs-response.png`",
        "",
        "## Top Functional Targets",
        "",
        markdown_table(top_functional),
        "",
        "## Interpretation",
        "",
        "The motif-burden signal is weak. AAG-like motif density is nearly uncorrelated with CLIP enrichment and is negatively associated with ribosome-density change.",
        "High-CLIP targets also do not have higher median AAG-like motif burden than low-CLIP controls in this simplified transcript-level scan.",
        "This supports the Cho et al. model more than a motif-only model: LIN28A recognizes AAG-rich sequence features, but target selection in mESCs also depends on context such as hairpin-loop presentation, transcript architecture, and the peri-ER localization that gives LIN28A preferential access to co-translationally targeted membrane and secretory mRNAs.",
        "",
        "## Caveats",
        "",
        "- The scan counts linear sequence motifs only and does not model the hairpin-loop context reported for the main LIN28A motif.",
        "- One transcript was selected per gene, so isoform-specific motif differences are hidden.",
        "- Motif counts are transcript-level, while CLIP and ribosome-density estimates are gene-level.",
        "- The permutation test asks only whether high-CLIP targets have higher median motif burden than low-CLIP controls; it does not test motif position, structure, or accessibility.",
    ]

    REPORT_PATH.write_text("\n".join(lines) + "\n")
    return REPORT_PATH


def main() -> None:
    SUBDATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = clean_analysis_table(pd.read_csv(MOTIF_TABLE_PATH, sep="\t"))
    summary = summarize(df)
    with SUMMARY_PATH.open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    target_scatter = make_target_scatter(df)
    boxplot = make_boxplot(df, summary)
    response_plot = make_response_plot(df, summary)
    report = write_report(df, summary)

    print(f"Wrote {SUMMARY_PATH.relative_to(W3_DIR)}")
    print(f"Wrote {target_scatter.relative_to(W3_DIR)}")
    print(f"Wrote {boxplot.relative_to(W3_DIR)}")
    print(f"Wrote {response_plot.relative_to(W3_DIR)}")
    print(f"Wrote {report.relative_to(W3_DIR)}")


if __name__ == "__main__":
    main()
