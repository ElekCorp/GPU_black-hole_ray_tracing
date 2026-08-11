#!/usr/bin/env python3
"""Mean/std per (variant, precision, region) from run_bench.sh's results.csv.

Usage: python3 benchmarks/summarize.py [path/to/results.csv]
"""
import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results" / "results.csv"

groups = defaultdict(list)
with open(path) as f:
    for row in csv.DictReader(f):
        key = (row["variant"], row["precision"], row["region"])
        try:
            groups[key].append(float(row["seconds"]))
        except ValueError:
            pass

VARIANT_LABEL = {
    "current": "current (dopri54 + christoffel_general/static)",
    "old87": "current dopri54 + old 84-temp christoffel",
    "rk6": "RK6 (fixed-step) + current christoffel",
    "old87_rk6": "RK6 (fixed-step) + old 84-temp christoffel",
}

print(f"{'variant':45s} {'prec':7s} {'region':7s} {'n':>3s} {'mean(s)':>9s} {'std(s)':>9s} {'cv%':>6s}")
for key in sorted(groups):
    vals = groups[key]
    n = len(vals)
    mean = st.mean(vals)
    std = st.stdev(vals) if n > 1 else 0.0
    cv = 100 * std / mean if mean else 0.0
    variant, prec, region = key
    print(f"{VARIANT_LABEL[variant]:45s} {prec:7s} {region:7s} {n:>3d} {mean:9.4f} {std:9.4f} {cv:6.1f}")
