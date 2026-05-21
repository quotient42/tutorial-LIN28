#!/usr/bin/env python3
"""
Run the Week 1 LIN28A target-group workflow.

This wrapper preserves the original command:
  python code/build_target_groups.py

The workflow is split into:
  1. 01_process_data.py
  2. 02_analyze_target_groups.py
  3. 03_visualize_and_report.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
STEPS = [
    "01_process_data.py",
    "02_analyze_target_groups.py",
    "03_visualize_and_report.py",
]


def main() -> None:
    for step in STEPS:
        print(f"\n=== {step} ===", flush=True)
        subprocess.run([sys.executable, str(SCRIPT_DIR / step)], check=True)


if __name__ == "__main__":
    main()
