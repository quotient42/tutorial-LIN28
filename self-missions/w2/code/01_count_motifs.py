#!/usr/bin/env python3
"""
Week 2 step 1: choose one transcript per gene and count LIN28A motifs.

Reads:
  - ../../data/gencode.gtf
  - ../../data/gencode.vM27.transcripts.fa.gz

Writes:
  - ../subdata/selected_transcripts.tsv
  - ../subdata/transcript_motif_counts_raw.tsv
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
W2_DIR = SCRIPT_DIR.parent
SELF_DIR = W2_DIR.parent
DATA_DIR = SELF_DIR / "data"
SUBDATA_DIR = W2_DIR / "subdata"

GTF_PATH = DATA_DIR / "gencode.gtf"
FASTA_PATH = DATA_DIR / "gencode.vM27.transcripts.fa.gz"
SELECTED_TRANSCRIPTS_PATH = SUBDATA_DIR / "selected_transcripts.tsv"
MOTIF_COUNTS_PATH = SUBDATA_DIR / "transcript_motif_counts_raw.tsv"

ATTR_RE = re.compile(r'(\S+) "([^"]+)"')
MOTIFS = {
    "aagnng": re.compile(r"(?=(AAG[ACGT]{2}G))"),
    "aagng": re.compile(r"(?=(AAG[ACGT]G))"),
    "tgtg": re.compile(r"(?=(TGTG))"),
}


def stable_id(versioned_id: str) -> str:
    return versioned_id.split(".", 1)[0]


def parse_attrs(attr_text: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attr_text))


def transcript_sort_key(row: dict[str, object]) -> tuple[int, int, int, str]:
    tsl = row["transcript_support_level"]
    try:
        tsl_rank = int(tsl)
    except (TypeError, ValueError):
        tsl_rank = 99

    appris = str(row["appris_tag"])
    if appris.startswith("appris_principal"):
        appris_rank = 0
    elif appris:
        appris_rank = 1
    else:
        appris_rank = 2

    return (tsl_rank, appris_rank, -int(row["transcript_length"]), str(row["transcript_id"]))


def parse_transcripts_from_gtf(gtf_path: Path) -> pd.DataFrame:
    transcripts: dict[str, dict[str, object]] = {}
    exon_lengths: dict[str, int] = {}
    with gtf_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]

            attrs = parse_attrs(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            stable_transcript_id = stable_id(transcript_id)
            if feature == "exon":
                exon_lengths[stable_transcript_id] = exon_lengths.get(stable_transcript_id, 0) + (
                    int(fields[4]) - int(fields[3]) + 1
                )
                continue

            if feature != "transcript":
                continue
            if attrs.get("gene_type") != "protein_coding":
                continue
            if attrs.get("transcript_type") != "protein_coding":
                continue

            tags = re.findall(r'tag "([^"]+)"', fields[8])
            appris_tags = [tag for tag in tags if tag.startswith("appris")]
            transcripts[stable_transcript_id] = {
                "gene_id": stable_id(attrs["gene_id"]),
                "gene_id_versioned": attrs["gene_id"],
                "gene_name": attrs.get("gene_name", stable_id(attrs["gene_id"])),
                "transcript_id": stable_transcript_id,
                "transcript_id_versioned": transcript_id,
                "transcript_name": attrs.get("transcript_name", ""),
                "transcript_support_level": attrs.get("transcript_support_level", "NA"),
                "appris_tag": appris_tags[0] if appris_tags else "",
                "transcript_length": 0,
                "chrom": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            }

    for transcript_id, row in transcripts.items():
        row["transcript_length"] = exon_lengths.get(transcript_id, 0)

    transcript_table = pd.DataFrame(transcripts.values())
    if transcript_table.empty:
        raise RuntimeError(f"No protein-coding transcripts parsed from {gtf_path}")
    return transcript_table


def choose_transcripts(transcripts: pd.DataFrame) -> pd.DataFrame:
    chosen_rows = []
    for _, group in transcripts.groupby("gene_id", sort=False):
        ordered = sorted(group.to_dict("records"), key=transcript_sort_key)
        chosen_rows.append(ordered[0])
    return pd.DataFrame(chosen_rows).sort_values("gene_id")


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
    first = header.split(None, 1)[0]
    return stable_id(first)


def load_selected_sequences(fasta_path: Path, selected_ids: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for header, seq in fasta_records(fasta_path):
        transcript_id = transcript_id_from_header(header)
        if transcript_id in selected_ids:
            sequences[transcript_id] = seq
    missing = selected_ids - set(sequences)
    if missing:
        preview = ", ".join(sorted(missing)[:5])
        raise RuntimeError(
            f"Missing {len(missing):,} selected transcript sequences from {fasta_path}: {preview}"
        )
    return sequences


def count_overlapping(pattern: re.Pattern[str], seq: str) -> int:
    return sum(1 for _ in pattern.finditer(seq))


def count_motifs(selected: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    rows = []
    for row in selected.to_dict("records"):
        seq = sequences[row["transcript_id"]]
        length = len(seq)
        counts = {name: count_overlapping(pattern, seq) for name, pattern in MOTIFS.items()}
        total_aag_like = counts["aagnng"] + counts["aagng"]
        rows.append(
            {
                **row,
                "sequence_length": length,
                "aagnng_count": counts["aagnng"],
                "aagng_count": counts["aagng"],
                "tgtg_count": counts["tgtg"],
                "aag_like_count": total_aag_like,
                "aagnng_per_kb": counts["aagnng"] / length * 1000,
                "aagng_per_kb": counts["aagng"] / length * 1000,
                "tgtg_per_kb": counts["tgtg"] / length * 1000,
                "aag_like_per_kb": total_aag_like / length * 1000,
            }
        )
    return pd.DataFrame(rows).sort_values("gene_id")


def main() -> None:
    if not FASTA_PATH.exists():
        raise FileNotFoundError(
            f"{FASTA_PATH} does not exist. Download GENCODE mouse vM27 transcript FASTA first:\n"
            f"  curl -L -o {FASTA_PATH} "
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/"
            "gencode.vM27.transcripts.fa.gz"
        )

    SUBDATA_DIR.mkdir(parents=True, exist_ok=True)
    transcripts = parse_transcripts_from_gtf(GTF_PATH)
    selected = choose_transcripts(transcripts)
    selected.to_csv(SELECTED_TRANSCRIPTS_PATH, sep="\t", index=False)

    sequences = load_selected_sequences(FASTA_PATH, set(selected["transcript_id"]))
    motif_counts = count_motifs(selected, sequences)
    motif_counts.to_csv(MOTIF_COUNTS_PATH, sep="\t", index=False)

    print(f"Wrote {SELECTED_TRANSCRIPTS_PATH.relative_to(W2_DIR)} ({len(selected):,} genes)")
    print(f"Wrote {MOTIF_COUNTS_PATH.relative_to(W2_DIR)}")


if __name__ == "__main__":
    main()
