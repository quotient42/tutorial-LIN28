#!/usr/bin/env python3
"""
Run the Week 2 LIN28A motif analysis workflow.

This wrapper runs:
  1. 01_count_motifs.py
  2. 02_merge_plot_report.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
STEPS = [
    "01_count_motifs.py",
    "02_merge_plot_report.py",
]


def main() -> None:
    for step in STEPS:
        print(f"\n=== {step} ===", flush=True)
        subprocess.run([sys.executable, str(SCRIPT_DIR / step)], check=True)


if __name__ == "__main__":
    main()
