#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Summarize formal cNMF K-selection statistics for GSE274498.

Input:
    Frozen cNMF k_selection_stats.df.npz

Outputs:
    formal_K5to20_selection_stats.csv

The formal analysis evaluated K=5-20.
K=8 was selected because it had the highest consensus stability
while retaining a comparatively parsimonious solution.
"""

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "source_data"
    / "cnmf"
    / "GSE274498_allMG_formal.k_selection_stats.df.npz"
)

OUT_DIR = PROJECT_ROOT / "results" / "cnmf"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = (
    OUT_DIR
    / "formal_K5to20_selection_stats.csv"
)


print("Reading frozen cNMF K-selection statistics...")

z = np.load(
    INPUT_FILE,
    allow_pickle=True,
)

required_keys = {
    "data",
    "index",
    "columns",
}

if not required_keys.issubset(z.files):
    raise ValueError(
        f"Unexpected NPZ structure. Keys found: {z.files}"
    )


stats = pd.DataFrame(
    z["data"],
    index=z["index"],
    columns=z["columns"],
)

required_columns = {
    "k",
    "local_density_threshold",
    "silhouette",
    "prediction_error",
}

missing = required_columns.difference(stats.columns)

if missing:
    raise ValueError(
        f"Missing K-selection columns: {sorted(missing)}"
    )


# Convert to numeric explicitly.
for column in required_columns:
    stats[column] = pd.to_numeric(
        stats[column],
        errors="raise",
    )


stats = (
    stats
    .sort_values("k")
    .reset_index(drop=True)
)


# Verify formal K range.
expected_k = list(range(5, 21))
observed_k = stats["k"].astype(int).tolist()

if observed_k != expected_k:
    raise ValueError(
        f"Unexpected K values: {observed_k}"
    )


# Identify highest consensus stability.
best_row = stats.loc[
    stats["silhouette"].idxmax()
]

best_k = int(best_row["k"])
best_stability = float(
    best_row["silhouette"]
)


stats.to_csv(
    OUTPUT_FILE,
    index=False,
)


print()
print(
    stats[
        [
            "k",
            "silhouette",
            "prediction_error",
        ]
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.6f}",
    )
)

print()
print(
    "HIGHEST STABILITY K:",
    best_k,
)

print(
    "HIGHEST STABILITY:",
    f"{best_stability:.6f}",
)

if best_k != 8:
    raise ValueError(
        f"Expected formal selected K=8, but highest stability was K={best_k}"
    )

print()
print("K=8 stability check: PASS")
print("Saved:", OUTPUT_FILE)
print("DONE")