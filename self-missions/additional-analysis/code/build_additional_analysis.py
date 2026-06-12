#!/usr/bin/env python3
"""
Additional LIN28A motif-context analyses.

Runs positional context on all selected transcripts and structural-context
folding on a balanced MFE-window pilot.

Reads:
  - ../../w2/output/transcript_motif_counts.tsv
  - ../../w2/subdata/selected_transcripts.tsv
  - ../../data/gencode.gtf
  - ../../data/gencode.vM27.transcripts.fa.gz
  - ../../../tmp/RNAstructure/exe/Fold

Writes:
  - ../subdata/motif_occurrences.tsv
  - ../subdata/transcript_region_annotations.tsv
  - ../subdata/positional_context_sites.tsv
  - ../subdata/positional_context_gene_summary.tsv
  - ../subdata/structural_context_sites.tsv
  - ../subdata/structural_context_gene_summary.tsv
  - ../subdata/additional_analysis_summary.json
  - ../output/*.png
  - ../output/report.md
"""

from __future__ import annotations

import gzip
import json
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
AA_DIR = SCRIPT_DIR.parent
SELF_DIR = AA_DIR.parent
PROJECT_DIR = SELF_DIR.parent
DATA_DIR = SELF_DIR / "data"
W2_DIR = SELF_DIR / "w2"
SUBDATA_DIR = AA_DIR / "subdata"
OUTPUT_DIR = AA_DIR / "output"

GTF_PATH = DATA_DIR / "gencode.gtf"
FASTA_PATH = DATA_DIR / "gencode.vM27.transcripts.fa.gz"
MOTIF_TABLE_PATH = W2_DIR / "output" / "transcript_motif_counts.tsv"
SELECTED_TRANSCRIPTS_PATH = W2_DIR / "subdata" / "selected_transcripts.tsv"

FOLD_BIN = PROJECT_DIR / "tmp" / "RNAstructure" / "exe" / "Fold"
RNASTRUCTURE_DATAPATH = PROJECT_DIR / "tmp" / "RNAstructure" / "data_tables"

MOTIF_OCCURRENCES_PATH = SUBDATA_DIR / "motif_occurrences.tsv"
REGION_ANNOTATIONS_PATH = SUBDATA_DIR / "transcript_region_annotations.tsv"
POSITIONAL_SITES_PATH = SUBDATA_DIR / "positional_context_sites.tsv"
POSITIONAL_GENE_SUMMARY_PATH = SUBDATA_DIR / "positional_context_gene_summary.tsv"
STRUCTURAL_WINDOWS_FASTA_PATH = SUBDATA_DIR / "motif_windows_pilot_61nt.fa"
STRUCTURAL_SITES_PATH = SUBDATA_DIR / "structural_context_sites.tsv"
STRUCTURAL_GENE_SUMMARY_PATH = SUBDATA_DIR / "structural_context_gene_summary.tsv"
SUMMARY_PATH = SUBDATA_DIR / "additional_analysis_summary.json"
REPORT_PATH = OUTPUT_DIR / "report.md"

POSITION_METAPLOT_PATH = OUTPUT_DIR / "positional-context-metaplot.png"
POSITION_REGION_BOXPLOT_PATH = OUTPUT_DIR / "positional-region-density-boxplot.png"
POSITION_VS_CLIP_PATH = OUTPUT_DIR / "positional-context-vs-clip.png"
POSITION_VS_RDEN_PATH = OUTPUT_DIR / "positional-context-vs-rden.png"
STRUCT_UNPAIRED_BOXPLOT_PATH = OUTPUT_DIR / "structural-context-unpaired-boxplot.png"
STRUCT_VS_CLIP_PATH = OUTPUT_DIR / "structural-context-vs-clip.png"
STRUCT_VS_RDEN_PATH = OUTPUT_DIR / "structural-context-vs-rden.png"

ATTR_RE = re.compile(r'(\S+) "([^"]+)"')
ENERGY_RE = re.compile(r"ENERGY\s*=\s*([-+]?\d+(?:\.\d+)?)")
RNG_SEED = 20260612

GROUP_ORDER = ["low_clip_control", "high_clip_target", "functional_target"]
GROUP_LABELS = {
    "low_clip_control": "Low-CLIP controls",
    "high_clip_target": "High-CLIP nonresponse",
    "functional_target": "Functional targets",
}
GROUP_PALETTE = {
    "low_clip_control": "#3b6ea8",
    "high_clip_target": "#f0a202",
    "functional_target": "#d73027",
}


def stable_id(versioned_id: str) -> str:
    return str(versioned_id).split(".", 1)[0]


def parse_attrs(attr_text: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attr_text))


def fasta_records(fasta_path: Path):
    opener = gzip.open if fasta_path.suffix == ".gz" else open
    with opener(fasta_path, "rt") as fh:
        header = None
        chunks: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks).upper()


def transcript_id_from_header(header: str) -> str:
    return stable_id(header.split(None, 1)[0])


def load_selected_sequences(selected_ids: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for header, seq in fasta_records(FASTA_PATH):
        transcript_id = transcript_id_from_header(header)
        if transcript_id in selected_ids:
            sequences[transcript_id] = seq.replace("U", "T")
    missing = selected_ids - set(sequences)
    if missing:
        preview = ", ".join(sorted(missing)[:5])
        raise RuntimeError(f"Missing {len(missing):,} selected transcript sequences: {preview}")
    return sequences


def find_overlapping(pattern: str, seq: str) -> list[tuple[int, int, str]]:
    compiled = re.compile(f"(?=({pattern}))")
    rows = []
    for match in compiled.finditer(seq):
        motif = match.group(1)
        start = match.start(1) + 1
        end = start + len(motif) - 1
        rows.append((start, end, motif))
    return rows


def nonredundant_aag_like_sites(seq: str) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    for start, end, motif in find_overlapping(r"AAG[ACGT]{2}G", seq):
        candidates.append(
            {
                "motif_pattern": "AAGNNG",
                "motif_start_1based": start,
                "motif_end_1based": end,
                "motif_sequence_dna": motif,
            }
        )
    for start, end, motif in find_overlapping(r"AAG[ACGT]G", seq):
        candidates.append(
            {
                "motif_pattern": "AAGNG",
                "motif_start_1based": start,
                "motif_end_1based": end,
                "motif_sequence_dna": motif,
            }
        )

    candidates.sort(
        key=lambda row: (
            int(row["motif_start_1based"]),
            -(int(row["motif_end_1based"]) - int(row["motif_start_1based"]) + 1),
            str(row["motif_pattern"]),
        )
    )
    kept: list[dict[str, object]] = []
    occupied: set[int] = set()
    for row in candidates:
        positions = set(range(int(row["motif_start_1based"]), int(row["motif_end_1based"]) + 1))
        if positions & occupied:
            continue
        kept.append(row)
        occupied.update(positions)
    return kept


def build_motif_occurrences(metrics: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    rows = []
    for record in metrics.to_dict("records"):
        transcript_id = record["transcript_id"]
        seq = sequences[transcript_id]
        length = len(seq)
        site_index = 0

        for site in nonredundant_aag_like_sites(seq):
            site_index += 1
            start = int(site["motif_start_1based"])
            end = int(site["motif_end_1based"])
            rows.append(
                {
                    "gene_id": record["gene_id"],
                    "gene_name": record["gene_name"],
                    "group": record["group"],
                    "log2_clip": record["log2_clip"],
                    "log2_rden": record["log2_rden"],
                    "transcript_id": transcript_id,
                    "transcript_id_versioned": record["transcript_id_versioned"],
                    "motif_id": f"{transcript_id}:aag_like:{site_index}",
                    "motif_family": "aag_like",
                    "motif_pattern": site["motif_pattern"],
                    "motif_start_1based": start,
                    "motif_end_1based": end,
                    "motif_center_1based": (start + end) / 2,
                    "motif_sequence_rna": str(site["motif_sequence_dna"]).replace("T", "U"),
                    "transcript_length": length,
                    "relative_position": ((start + end) / 2) / length,
                    "distance_to_transcript_5p": start - 1,
                    "distance_to_transcript_3p": length - end,
                }
            )

        for idx, (start, end, motif) in enumerate(find_overlapping(r"TGTG", seq), start=1):
            rows.append(
                {
                    "gene_id": record["gene_id"],
                    "gene_name": record["gene_name"],
                    "group": record["group"],
                    "log2_clip": record["log2_clip"],
                    "log2_rden": record["log2_rden"],
                    "transcript_id": transcript_id,
                    "transcript_id_versioned": record["transcript_id_versioned"],
                    "motif_id": f"{transcript_id}:tgtg:{idx}",
                    "motif_family": "tgtg",
                    "motif_pattern": "TGTG",
                    "motif_start_1based": start,
                    "motif_end_1based": end,
                    "motif_center_1based": (start + end) / 2,
                    "motif_sequence_rna": motif.replace("T", "U"),
                    "transcript_length": length,
                    "relative_position": ((start + end) / 2) / length,
                    "distance_to_transcript_5p": start - 1,
                    "distance_to_transcript_3p": length - end,
                }
            )

    occurrences = pd.DataFrame(rows)
    return occurrences.sort_values(["gene_id", "transcript_id", "motif_start_1based", "motif_family"])


def parse_gtf_for_regions(selected: pd.DataFrame) -> pd.DataFrame:
    selected_ids = set(selected["transcript_id"])
    exons: dict[str, list[tuple[int, int, str]]] = {tid: [] for tid in selected_ids}
    cds_segments: dict[str, list[tuple[int, int, str]]] = {tid: [] for tid in selected_ids}

    with GTF_PATH.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in {"exon", "CDS"}:
                continue
            attrs = parse_attrs(fields[8])
            transcript_id = stable_id(attrs.get("transcript_id", ""))
            if transcript_id not in selected_ids:
                continue
            segment = (int(fields[3]), int(fields[4]), fields[6])
            if feature == "exon":
                exons[transcript_id].append(segment)
            elif feature == "CDS":
                cds_segments[transcript_id].append(segment)

    rows = []
    selected_by_tid = selected.set_index("transcript_id")
    for transcript_id in selected_ids:
        transcript_row = selected_by_tid.loc[transcript_id]
        transcript_length = int(transcript_row["transcript_length"])
        transcript_exons = exons.get(transcript_id, [])
        transcript_cds = cds_segments.get(transcript_id, [])
        strand = str(transcript_row["strand"])

        cds_coords: list[int] = []
        if transcript_exons and transcript_cds:
            ordered_exons = sorted(transcript_exons, key=lambda x: x[0], reverse=(strand == "-"))
            offset = 0
            for exon_start, exon_end, _ in ordered_exons:
                exon_len = exon_end - exon_start + 1
                for cds_start, cds_end, _ in transcript_cds:
                    overlap_start = max(exon_start, cds_start)
                    overlap_end = min(exon_end, cds_end)
                    if overlap_start > overlap_end:
                        continue
                    if strand == "-":
                        t_start = offset + (exon_end - overlap_end) + 1
                        t_end = offset + (exon_end - overlap_start) + 1
                    else:
                        t_start = offset + (overlap_start - exon_start) + 1
                        t_end = offset + (overlap_end - exon_start) + 1
                    cds_coords.extend([t_start, t_end])
                offset += exon_len

        if cds_coords:
            cds_start = int(min(cds_coords))
            cds_end = int(max(cds_coords))
            annotation_status = "cds_annotated"
        else:
            cds_start = math.nan
            cds_end = math.nan
            annotation_status = "missing_cds"

        rows.append(
            {
                "transcript_id": transcript_id,
                "transcript_id_versioned": transcript_row["transcript_id_versioned"],
                "gene_id": transcript_row["gene_id"],
                "gene_name": transcript_row["gene_name"],
                "transcript_length": transcript_length,
                "cds_start_1based": cds_start,
                "cds_end_1based": cds_end,
                "utr5_length": int(cds_start - 1) if not math.isnan(cds_start) else math.nan,
                "cds_length": int(cds_end - cds_start + 1) if not math.isnan(cds_start) else math.nan,
                "utr3_length": int(transcript_length - cds_end) if not math.isnan(cds_start) else math.nan,
                "annotation_status": annotation_status,
            }
        )

    return pd.DataFrame(rows).sort_values("gene_id")


def region_for_position(center: float, cds_start: float, cds_end: float) -> str:
    if pd.isna(cds_start) or pd.isna(cds_end):
        return "unknown_or_unannotated"
    if center < cds_start:
        return "5UTR"
    if center > cds_end:
        return "3UTR"
    return "CDS"


def annotate_positional_sites(occurrences: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    merged = occurrences.merge(
        regions[
            [
                "transcript_id",
                "cds_start_1based",
                "cds_end_1based",
                "utr5_length",
                "cds_length",
                "utr3_length",
                "annotation_status",
            ]
        ],
        on="transcript_id",
        how="left",
    )
    merged["region"] = [
        region_for_position(center, cds_start, cds_end)
        for center, cds_start, cds_end in zip(
            merged["motif_center_1based"], merged["cds_start_1based"], merged["cds_end_1based"]
        )
    ]
    merged["distance_to_cds_start"] = (merged["motif_center_1based"] - merged["cds_start_1based"]).abs()
    merged["distance_to_cds_stop"] = (merged["motif_center_1based"] - merged["cds_end_1based"]).abs()
    merged["distance_to_nearest_cds_boundary"] = merged[["distance_to_cds_start", "distance_to_cds_stop"]].min(axis=1)
    merged["in_start_proximal_window"] = merged["distance_to_cds_start"] <= 100
    merged["in_stop_proximal_window"] = merged["distance_to_cds_stop"] <= 100

    nearest_distances = pd.Series(np.nan, index=merged.index, dtype=float)
    cluster_50 = pd.Series(0, index=merged.index, dtype=int)
    cluster_100 = pd.Series(0, index=merged.index, dtype=int)
    for _, group in merged.groupby(["transcript_id", "motif_family"], sort=False):
        centers = group["motif_center_1based"].to_numpy(dtype=float)
        for idx, center in zip(group.index, centers):
            diffs = np.abs(centers - center)
            nonzero = diffs[diffs > 0]
            nearest_distances.loc[idx] = float(nonzero.min()) if len(nonzero) else math.nan
            cluster_50.loc[idx] = int(np.sum(diffs <= 50))
            cluster_100.loc[idx] = int(np.sum(diffs <= 100))

    merged["distance_to_nearest_same_family_motif"] = nearest_distances
    merged["cluster_size_50nt"] = cluster_50
    merged["cluster_size_100nt"] = cluster_100
    return merged.sort_values(["gene_id", "transcript_id", "motif_start_1based", "motif_family"])


def safe_divide(numerator: float, denominator: float) -> float:
    if denominator is None or pd.isna(denominator) or denominator == 0:
        return math.nan
    return float(numerator) / float(denominator)


def build_positional_gene_summary(positional_sites: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    base_cols = [
        "gene_id",
        "gene_name",
        "group",
        "log2_clip",
        "log2_rden",
        "transcript_id",
        "transcript_length",
        "aag_like_per_kb",
    ]
    base = metrics[base_cols].copy()
    aag_sites = positional_sites[positional_sites["motif_family"] == "aag_like"].copy()
    summary = base.set_index("gene_id", drop=False)

    n_sites = aag_sites.groupby("gene_id").size().rename("n_aag_like_sites")
    summary = summary.join(n_sites, how="left")
    summary["n_aag_like_sites"] = summary["n_aag_like_sites"].fillna(0).astype(int)

    region_counts = (
        aag_sites.groupby(["gene_id", "region"]).size().unstack(fill_value=0)
        if not aag_sites.empty
        else pd.DataFrame()
    )
    for region, col in [
        ("5UTR", "aag_like_5utr_count"),
        ("CDS", "aag_like_cds_count"),
        ("3UTR", "aag_like_3utr_count"),
    ]:
        summary[col] = region_counts[region] if region in region_counts.columns else 0
        summary[col] = summary[col].fillna(0).astype(int)

    for region_col, fraction_col in [
        ("aag_like_5utr_count", "fraction_aag_like_5utr"),
        ("aag_like_cds_count", "fraction_aag_like_cds"),
        ("aag_like_3utr_count", "fraction_aag_like_3utr"),
    ]:
        summary[fraction_col] = summary[region_col] / summary["n_aag_like_sites"].replace(0, np.nan)

    grouped = aag_sites.groupby("gene_id")
    aggregates = grouped.agg(
        n_aag_like_start_proximal=("in_start_proximal_window", "sum"),
        n_aag_like_stop_proximal=("in_stop_proximal_window", "sum"),
        mean_relative_position=("relative_position", "mean"),
        median_relative_position=("relative_position", "median"),
        max_cluster_size_50nt=("cluster_size_50nt", "max"),
        max_cluster_size_100nt=("cluster_size_100nt", "max"),
        mean_nearest_neighbor_distance=("distance_to_nearest_same_family_motif", "mean"),
        utr5_length=("utr5_length", "first"),
        cds_length=("cds_length", "first"),
        utr3_length=("utr3_length", "first"),
    )
    summary = summary.join(aggregates, how="left")

    for col in [
        "n_aag_like_start_proximal",
        "n_aag_like_stop_proximal",
        "max_cluster_size_50nt",
        "max_cluster_size_100nt",
    ]:
        summary[col] = summary[col].fillna(0).astype(int)

    summary["aag_like_5utr_per_kb"] = summary["aag_like_5utr_count"] / (summary["utr5_length"] / 1000)
    summary["aag_like_cds_per_kb"] = summary["aag_like_cds_count"] / (summary["cds_length"] / 1000)
    summary["aag_like_3utr_per_kb"] = summary["aag_like_3utr_count"] / (summary["utr3_length"] / 1000)
    summary["start_proximal_density"] = summary["n_aag_like_start_proximal"] / 0.2
    summary["stop_proximal_density"] = summary["n_aag_like_stop_proximal"] / 0.2
    summary["cluster_score_100nt"] = summary["max_cluster_size_100nt"]

    drop_cols = ["utr5_length", "cds_length", "utr3_length"]
    return summary.reset_index(drop=True).drop(columns=drop_cols).sort_values("gene_id")


def pearsonr(x: pd.Series, y: pd.Series) -> float:
    clean = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) < 2:
        return math.nan
    return float(np.corrcoef(clean.iloc[:, 0], clean.iloc[:, 1])[0, 1])


def spearmanr(x: pd.Series, y: pd.Series) -> float:
    clean = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) < 2:
        return math.nan
    return pearsonr(clean.iloc[:, 0].rank(method="average"), clean.iloc[:, 1].rank(method="average"))


def permutation_median_greater(
    target_values: np.ndarray,
    control_values: np.ndarray,
    n_permutations: int = 10000,
) -> dict[str, float]:
    target_values = np.asarray(target_values, dtype=float)
    control_values = np.asarray(control_values, dtype=float)
    target_values = target_values[np.isfinite(target_values)]
    control_values = control_values[np.isfinite(control_values)]
    observed = float(np.median(target_values) - np.median(control_values))
    pooled = np.concatenate([target_values, control_values])
    n_target = len(target_values)
    rng = np.random.default_rng(RNG_SEED)
    hits = 0
    for _ in range(n_permutations):
        permuted = rng.permutation(pooled)
        diff = float(np.median(permuted[:n_target]) - np.median(permuted[n_target:]))
        if diff >= observed:
            hits += 1
    return {
        "observed_median_difference": observed,
        "p_value_target_greater": float((hits + 1) / (n_permutations + 1)),
        "n_permutations": int(n_permutations),
        "n_target": int(n_target),
        "n_control": int(len(control_values)),
    }


def summarize_feature_by_group(df: pd.DataFrame, feature: str) -> dict[str, object]:
    controls = df[df["group"] == "low_clip_control"][feature].to_numpy(dtype=float)
    functional = df[df["group"] == "functional_target"][feature].to_numpy(dtype=float)
    high_nonresponse = df[df["group"] == "high_clip_target"][feature].to_numpy(dtype=float)
    high_all = df[df["group"].isin(["functional_target", "high_clip_target"])][feature].to_numpy(dtype=float)
    return {
        "median_by_group": {
            group: float(df.loc[df["group"] == group, feature].median())
            for group in GROUP_ORDER
            if group in set(df["group"])
        },
        "spearman_log2_clip": spearmanr(df[feature], df["log2_clip"]),
        "spearman_log2_rden": spearmanr(df[feature], df["log2_rden"]),
        "functional_vs_low_clip": permutation_median_greater(functional, controls),
        "all_high_clip_vs_low_clip": permutation_median_greater(high_all, controls),
        "functional_vs_high_clip_nonresponse": permutation_median_greater(functional, high_nonresponse),
    }


def select_structural_pilot_sites(positional_sites: pd.DataFrame) -> pd.DataFrame:
    aag = positional_sites[positional_sites["motif_family"] == "aag_like"].copy()
    rng = np.random.default_rng(RNG_SEED)
    selected_gene_ids: list[str] = []
    for group in GROUP_ORDER:
        genes = sorted(aag.loc[aag["group"] == group, "gene_id"].unique())
        if len(genes) > 50:
            genes = sorted(rng.choice(genes, size=50, replace=False).tolist())
        selected_gene_ids.extend(genes)

    pilot = aag[aag["gene_id"].isin(selected_gene_ids)].copy()
    sampled_rows = []
    for _, group_df in pilot.groupby("gene_id", sort=False):
        if len(group_df) > 3:
            sampled_rows.append(group_df.sample(n=3, random_state=RNG_SEED))
        else:
            sampled_rows.append(group_df)
    return pd.concat(sampled_rows, ignore_index=True).sort_values(["group", "gene_id", "motif_start_1based"])


def extract_window(seq: str, start: int, end: int, flank: int = 30) -> dict[str, object]:
    left = max(1, start - flank)
    right = min(len(seq), end + flank)
    window = seq[left - 1 : right].replace("T", "U")
    return {
        "window_start_1based": left,
        "window_end_1based": right,
        "window_sequence": window,
        "motif_window_start_1based": start - left + 1,
        "motif_window_end_1based": end - left + 1,
    }


def write_windows_fasta(pilot_sites: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    rows = []
    with STRUCTURAL_WINDOWS_FASTA_PATH.open("w") as fh:
        for record in pilot_sites.to_dict("records"):
            window = extract_window(
                sequences[record["transcript_id"]],
                int(record["motif_start_1based"]),
                int(record["motif_end_1based"]),
                flank=30,
            )
            row = {**record, **window}
            rows.append(row)
            fh.write(f">{row['motif_id']}\n{row['window_sequence']}\n")
    return pd.DataFrame(rows)


def parse_rnastructure_dbn(dbn_path: Path) -> tuple[float, str, str]:
    lines = [line.strip() for line in dbn_path.read_text().splitlines() if line.strip()]
    if len(lines) < 3:
        raise RuntimeError(f"Unexpected DBN output in {dbn_path}")
    energy_match = ENERGY_RE.search(lines[0])
    energy = float(energy_match.group(1)) if energy_match else math.nan
    return energy, lines[1], lines[2]


def fold_one_window(motif_id: str, sequence: str, tmp_dir: Path, env: dict[str, str]) -> tuple[float, str, str, str]:
    safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", motif_id)
    fasta_path = tmp_dir / f"{safe_name}.fa"
    dbn_path = tmp_dir / f"{safe_name}.dbn"
    fasta_path.write_text(f">{motif_id}\n{sequence}\n")
    cmd = [str(FOLD_BIN), str(fasta_path), str(dbn_path), "--MFE", "--bracket", "--quiet"]
    completed = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=30)
    if completed.returncode != 0:
        message = (completed.stderr or completed.stdout).strip().splitlines()
        return math.nan, "", "", message[-1] if message else f"Fold failed with code {completed.returncode}"
    energy, folded_seq, structure = parse_rnastructure_dbn(dbn_path)
    return energy, folded_seq, structure, ""


def contiguous_unpaired_run(structure: str, center_index_0based: int) -> int:
    if not structure or center_index_0based < 0 or center_index_0based >= len(structure):
        return 0
    if structure[center_index_0based] != ".":
        return 0
    left = center_index_0based
    while left - 1 >= 0 and structure[left - 1] == ".":
        left -= 1
    right = center_index_0based
    while right + 1 < len(structure) and structure[right + 1] == ".":
        right += 1
    return right - left + 1


def has_paired_base(structure: str, start_1based: int, end_1based: int) -> bool:
    start = max(1, start_1based)
    end = min(len(structure), end_1based)
    if start > end:
        return False
    return any(char in "()" for char in structure[start - 1 : end])


def structural_features_for_site(row: dict[str, object]) -> dict[str, object]:
    structure = str(row["dot_bracket"])
    motif_start = int(row["motif_window_start_1based"])
    motif_end = int(row["motif_window_end_1based"])
    motif_chars = structure[motif_start - 1 : motif_end]
    unpaired_count = motif_chars.count(".")
    motif_len = max(1, motif_end - motif_start + 1)
    center = int(round((motif_start + motif_end) / 2)) - 1
    left_paired = has_paired_base(structure, motif_start - 10, motif_start - 1)
    right_paired = has_paired_base(structure, motif_end + 1, motif_end + 10)
    return {
        "window_length": len(str(row["window_sequence"])),
        "mfe_per_nt": safe_divide(float(row["mfe_kcal_mol"]), len(str(row["window_sequence"]))),
        "motif_unpaired_count": unpaired_count,
        "motif_unpaired_fraction": unpaired_count / motif_len,
        "motif_all_unpaired": unpaired_count == motif_len,
        "motif_any_paired": unpaired_count < motif_len,
        "motif_loop_run_length": contiguous_unpaired_run(structure, center),
        "paired_bases_left_flank_10nt": left_paired,
        "paired_bases_right_flank_10nt": right_paired,
        "motif_in_terminal_loop_candidate": (unpaired_count == motif_len and left_paired and right_paired),
        "gc_fraction_window": (str(row["window_sequence"]).count("G") + str(row["window_sequence"]).count("C"))
        / max(1, len(str(row["window_sequence"]))),
    }


def run_structural_pilot(positional_sites: pd.DataFrame, sequences: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not FOLD_BIN.exists():
        raise FileNotFoundError(f"RNAstructure Fold executable not found: {FOLD_BIN}")
    if not RNASTRUCTURE_DATAPATH.exists():
        raise FileNotFoundError(f"RNAstructure DATAPATH not found: {RNASTRUCTURE_DATAPATH}")

    pilot_sites = select_structural_pilot_sites(positional_sites)
    windows = write_windows_fasta(pilot_sites, sequences)
    env = os.environ.copy()
    env["DATAPATH"] = str(RNASTRUCTURE_DATAPATH)

    folded_rows = []
    with tempfile.TemporaryDirectory(prefix="lin28_struct_", dir="/tmp") as tmp_name:
        tmp_dir = Path(tmp_name)
        for idx, row in enumerate(windows.to_dict("records"), start=1):
            energy, folded_seq, structure, error = fold_one_window(
                str(row["motif_id"]), str(row["window_sequence"]), tmp_dir, env
            )
            if idx % 25 == 0:
                print(f"Folded {idx:,}/{len(windows):,} motif windows", flush=True)
            folded = {
                **row,
                "mfe_kcal_mol": energy,
                "folded_sequence": folded_seq,
                "dot_bracket": structure,
                "fold_error": error,
            }
            if not error:
                folded.update(structural_features_for_site(folded))
            folded_rows.append(folded)

    sites = pd.DataFrame(folded_rows)
    successful = sites[sites["fold_error"].fillna("") == ""].copy()
    rows = []
    for (gene_id, gene_name, group, transcript_id), group_df in successful.groupby(
        ["gene_id", "gene_name", "group", "transcript_id"], sort=False
    ):
        rows.append(
            {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "group": group,
                "transcript_id": transcript_id,
                "log2_clip": float(group_df["log2_clip"].iloc[0]),
                "log2_rden": float(group_df["log2_rden"].iloc[0]),
                "n_folded_aag_like_sites": int(len(group_df)),
                "mean_motif_unpaired_fraction": float(group_df["motif_unpaired_fraction"].mean()),
                "max_motif_unpaired_fraction": float(group_df["motif_unpaired_fraction"].max()),
                "fraction_sites_all_unpaired": float(group_df["motif_all_unpaired"].mean()),
                "fraction_sites_terminal_loop_candidate": float(group_df["motif_in_terminal_loop_candidate"].mean()),
                "best_site_unpaired_fraction": float(group_df["motif_unpaired_fraction"].max()),
                "best_site_mfe_per_nt": float(group_df["mfe_per_nt"].min()),
                "mean_mfe_per_nt": float(group_df["mfe_per_nt"].mean()),
            }
        )
    gene_summary = pd.DataFrame(rows).sort_values("gene_id")
    return sites, gene_summary


def setup_plot(ax: plt.Axes) -> None:
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_metaplot(positional_sites: pd.DataFrame) -> None:
    aag = positional_sites[positional_sites["motif_family"] == "aag_like"].copy()
    bins = np.linspace(0, 1, 41)
    centers = (bins[:-1] + bins[1:]) / 2
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for group in GROUP_ORDER:
        sub = aag[aag["group"] == group]
        counts, _ = np.histogram(sub["relative_position"], bins=bins)
        n_genes = max(1, sub["gene_id"].nunique())
        density = counts / n_genes
        ax.plot(centers, density, color=GROUP_PALETTE[group], lw=1.8, label=GROUP_LABELS[group])
    ax.set_xlabel("Relative position in selected transcript")
    ax.set_ylabel("AAG-like motifs per gene per bin")
    ax.set_title("Transcript-relative motif density")
    ax.legend(frameon=False, fontsize=7)
    setup_plot(ax)
    fig.tight_layout()
    fig.savefig(POSITION_METAPLOT_PATH, dpi=220)
    plt.close(fig)


def make_region_boxplot(gene_summary: pd.DataFrame) -> None:
    features = [
        ("aag_like_5utr_per_kb", "5UTR"),
        ("aag_like_cds_per_kb", "CDS"),
        ("aag_like_3utr_per_kb", "3UTR"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.7), sharey=False)
    for ax, (feature, label) in zip(axes, features):
        data = [gene_summary.loc[gene_summary["group"] == group, feature].dropna() for group in GROUP_ORDER]
        bp = ax.boxplot(
            data,
            patch_artist=True,
            tick_labels=["Low", "High", "Func"],
            showfliers=False,
            medianprops={"color": "black", "lw": 1.0},
        )
        for patch, group in zip(bp["boxes"], GROUP_ORDER):
            patch.set_facecolor(GROUP_PALETTE[group])
            patch.set_alpha(0.55)
        ax.set_title(label)
        ax.set_ylabel("AAG-like motifs per kb" if feature == "aag_like_5utr_per_kb" else "")
        setup_plot(ax)
    fig.suptitle("Regional motif density by group", y=1.02)
    fig.tight_layout()
    fig.savefig(POSITION_REGION_BOXPLOT_PATH, dpi=220)
    plt.close(fig)


def make_feature_scatter(df: pd.DataFrame, feature: str, y_col: str, y_label: str, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    for group in GROUP_ORDER:
        sub = df[df["group"] == group]
        ax.scatter(
            sub[feature],
            sub[y_col],
            s=13,
            c=GROUP_PALETTE[group],
            alpha=0.45,
            edgecolors="none",
            label=GROUP_LABELS[group],
        )
    ax.axhline(0, color="black", lw=0.6, alpha=0.55)
    rho = spearmanr(df[feature], df[y_col])
    ax.text(0.98, 0.04, f"Spearman rho = {rho:.3f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
    ax.set_xlabel(feature.replace("_", " "))
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7)
    setup_plot(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def make_structural_boxplot(gene_summary: pd.DataFrame) -> None:
    feature = "fraction_sites_all_unpaired"
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    data = [gene_summary.loc[gene_summary["group"] == group, feature].dropna() for group in GROUP_ORDER]
    bp = ax.boxplot(
        data,
        patch_artist=True,
        tick_labels=["Low-CLIP", "High-CLIP", "Functional"],
        showfliers=False,
        medianprops={"color": "black", "lw": 1.0},
    )
    for patch, group in zip(bp["boxes"], GROUP_ORDER):
        patch.set_facecolor(GROUP_PALETTE[group])
        patch.set_alpha(0.55)
    ax.set_ylabel("Fraction of folded sites all-unpaired")
    ax.set_title("Predicted local motif accessibility")
    setup_plot(ax)
    fig.tight_layout()
    fig.savefig(STRUCT_UNPAIRED_BOXPLOT_PATH, dpi=220)
    plt.close(fig)


def write_report(summary: dict[str, object]) -> None:
    pos = summary["positional_features"]
    struct = summary["structural_features"]
    text = f"""# Additional Analysis Report: Motif Position and Predicted Structure

## Goal

This additional analysis asks whether the weak relationship between transcript-wide LIN28A motif burden and LIN28A targeting can be refined by considering motif position and predicted local RNA structure.

## Inputs

- Week 2 motif table: `self-missions/w2/output/transcript_motif_counts.tsv`
- Selected transcripts: `self-missions/w2/subdata/selected_transcripts.tsv`
- GENCODE annotation: `self-missions/data/gencode.gtf`
- Transcript FASTA: `self-missions/data/gencode.vM27.transcripts.fa.gz`
- Structural folding tool: `tmp/RNAstructure/exe/Fold`

## Positional Context

- Motif occurrence rows: {summary['n_motif_occurrences']:,}
- AAG-like motif occurrence rows: {summary['n_aag_like_occurrences']:,}
- Genes summarized: {summary['n_positional_genes']:,}
- Spearman rho, CDS motif fraction vs log2 CLIP: {pos['fraction_aag_like_cds']['spearman_log2_clip']:.4f}
- Spearman rho, CDS motif fraction vs log2 ribosome-density change: {pos['fraction_aag_like_cds']['spearman_log2_rden']:.4f}
- Functional-vs-control median difference for CDS motif fraction: {pos['fraction_aag_like_cds']['functional_vs_low_clip']['observed_median_difference']:.4f}
- Functional-vs-control p-value for higher CDS motif fraction: {pos['fraction_aag_like_cds']['functional_vs_low_clip']['p_value_target_greater']:.4f}

## Structural Context

Structural analysis was run as a balanced local-MFE pilot, not a full-transcript or all-site fold.

- Folded site rows attempted: {summary['n_structural_sites_attempted']:,}
- Folded site rows successful: {summary['n_structural_sites_successful']:,}
- Genes summarized: {summary['n_structural_genes']:,}
- Spearman rho, fraction all-unpaired sites vs log2 CLIP: {struct['fraction_sites_all_unpaired']['spearman_log2_clip']:.4f}
- Spearman rho, fraction all-unpaired sites vs log2 ribosome-density change: {struct['fraction_sites_all_unpaired']['spearman_log2_rden']:.4f}
- Functional-vs-control median difference for fraction all-unpaired sites: {struct['fraction_sites_all_unpaired']['functional_vs_low_clip']['observed_median_difference']:.4f}
- Functional-vs-control p-value for higher fraction all-unpaired sites: {struct['fraction_sites_all_unpaired']['functional_vs_low_clip']['p_value_target_greater']:.4f}

## Outputs

- `subdata/motif_occurrences.tsv`
- `subdata/transcript_region_annotations.tsv`
- `subdata/positional_context_sites.tsv`
- `subdata/positional_context_gene_summary.tsv`
- `subdata/motif_windows_pilot_61nt.fa`
- `subdata/structural_context_sites.tsv`
- `subdata/structural_context_gene_summary.tsv`
- `subdata/additional_analysis_summary.json`
- `output/positional-context-metaplot.png`
- `output/positional-region-density-boxplot.png`
- `output/positional-context-vs-clip.png`
- `output/positional-context-vs-rden.png`
- `output/structural-context-unpaired-boxplot.png`
- `output/structural-context-vs-clip.png`
- `output/structural-context-vs-rden.png`

## Interpretation Boundary

These analyses inspect positional and predicted structural context as possible refinements of the motif-burden model. They do not directly measure in vivo RNA structure, motif accessibility, ribosome position, ER localization, or LIN28A occupancy at individual motif sites.
"""
    REPORT_PATH.write_text(text)


def main() -> None:
    SUBDATA_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading Week 2 motif table and selected transcript sequences...", flush=True)
    metrics = pd.read_csv(MOTIF_TABLE_PATH, sep="\t")
    metrics = metrics.replace([np.inf, -np.inf], np.nan).dropna(subset=["group", "log2_clip", "log2_rden"])
    selected = pd.read_csv(SELECTED_TRANSCRIPTS_PATH, sep="\t")
    selected_ids = set(metrics["transcript_id"])
    sequences = load_selected_sequences(selected_ids)

    if MOTIF_OCCURRENCES_PATH.exists():
        print("Loading cached motif occurrence table...", flush=True)
        occurrences = pd.read_csv(MOTIF_OCCURRENCES_PATH, sep="\t")
    else:
        print("Building motif occurrence table...", flush=True)
        occurrences = build_motif_occurrences(metrics, sequences)
        occurrences.to_csv(MOTIF_OCCURRENCES_PATH, sep="\t", index=False)

    if REGION_ANNOTATIONS_PATH.exists() and POSITIONAL_SITES_PATH.exists() and POSITIONAL_GENE_SUMMARY_PATH.exists():
        print("Loading cached positional context tables...", flush=True)
        positional_sites = pd.read_csv(POSITIONAL_SITES_PATH, sep="\t")
        positional_gene_summary = pd.read_csv(POSITIONAL_GENE_SUMMARY_PATH, sep="\t")
    else:
        print("Annotating transcript regions and motif positions...", flush=True)
        regions = parse_gtf_for_regions(selected[selected["transcript_id"].isin(selected_ids)].copy())
        regions.to_csv(REGION_ANNOTATIONS_PATH, sep="\t", index=False)
        positional_sites = annotate_positional_sites(occurrences, regions)
        positional_sites.to_csv(POSITIONAL_SITES_PATH, sep="\t", index=False)
        positional_gene_summary = build_positional_gene_summary(positional_sites, metrics)
        positional_gene_summary.to_csv(POSITIONAL_GENE_SUMMARY_PATH, sep="\t", index=False)

    print("Running structural-context pilot with RNAstructure Fold...", flush=True)
    structural_sites, structural_gene_summary = run_structural_pilot(positional_sites, sequences)
    structural_sites.to_csv(STRUCTURAL_SITES_PATH, sep="\t", index=False)
    structural_gene_summary.to_csv(STRUCTURAL_GENE_SUMMARY_PATH, sep="\t", index=False)

    print("Making figures...", flush=True)
    make_metaplot(positional_sites)
    make_region_boxplot(positional_gene_summary)
    make_feature_scatter(
        positional_gene_summary,
        "fraction_aag_like_cds",
        "log2_clip",
        r"LIN28A CLIP enrichment (log$_2$)",
        POSITION_VS_CLIP_PATH,
        "CDS motif fraction vs CLIP enrichment",
    )
    make_feature_scatter(
        positional_gene_summary,
        "fraction_aag_like_cds",
        "log2_rden",
        "Ribosome density change after\n" + r"$\it{Lin28a}$ knockdown (log$_2$)",
        POSITION_VS_RDEN_PATH,
        "CDS motif fraction vs translation response",
    )
    make_structural_boxplot(structural_gene_summary)
    make_feature_scatter(
        structural_gene_summary,
        "fraction_sites_all_unpaired",
        "log2_clip",
        r"LIN28A CLIP enrichment (log$_2$)",
        STRUCT_VS_CLIP_PATH,
        "Predicted accessibility vs CLIP enrichment",
    )
    make_feature_scatter(
        structural_gene_summary,
        "fraction_sites_all_unpaired",
        "log2_rden",
        "Ribosome density change after\n" + r"$\it{Lin28a}$ knockdown (log$_2$)",
        STRUCT_VS_RDEN_PATH,
        "Predicted accessibility vs translation response",
    )

    summary = {
        "n_motif_occurrences": int(len(occurrences)),
        "n_aag_like_occurrences": int((occurrences["motif_family"] == "aag_like").sum()),
        "n_positional_genes": int(len(positional_gene_summary)),
        "n_structural_sites_attempted": int(len(structural_sites)),
        "n_structural_sites_successful": int((structural_sites["fold_error"].fillna("") == "").sum()),
        "n_structural_genes": int(len(structural_gene_summary)),
        "positional_features": {
            feature: summarize_feature_by_group(positional_gene_summary, feature)
            for feature in [
                "fraction_aag_like_cds",
                "fraction_aag_like_3utr",
                "start_proximal_density",
                "stop_proximal_density",
                "cluster_score_100nt",
            ]
        },
        "structural_features": {
            feature: summarize_feature_by_group(structural_gene_summary, feature)
            for feature in [
                "fraction_sites_all_unpaired",
                "fraction_sites_terminal_loop_candidate",
                "mean_motif_unpaired_fraction",
            ]
        },
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True))
    write_report(summary)

    print(f"Wrote {MOTIF_OCCURRENCES_PATH.relative_to(AA_DIR)}")
    print(f"Wrote {POSITIONAL_GENE_SUMMARY_PATH.relative_to(AA_DIR)}")
    print(f"Wrote {STRUCTURAL_GENE_SUMMARY_PATH.relative_to(AA_DIR)}")
    print(f"Wrote {REPORT_PATH.relative_to(AA_DIR)}")


if __name__ == "__main__":
    main()
