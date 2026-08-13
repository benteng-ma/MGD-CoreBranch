#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reproduce formal K=8 cNMF usage summaries for GSE274498.

Inputs
------
1. Frozen K=8 consensus usage matrix
2. Frozen refined MG-state metadata

Analysis definition
-------------------
Full K=8 normalization:
    For each cell, divide P1-P8 usage by the sum of P1-P8 usage.

P2-excluded sensitivity:
    Remove P2, then renormalize P1/P3-P8 within each cell so that
    the remaining seven program usages sum to 1.

Outputs
-------
- formal_K8_mean_usage_by_MG_state.csv
- formal_K8_mean_usage_by_MG_state_sample.csv
- formal_K8_mean_usage_by_MG_state_sample_P2excluded.csv
- formal_K8_age_direction_sensitivity_summary.csv
- formal_K8_P7_P8_library_age_summary.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd


# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

USAGE_FILE = (
    PROJECT_ROOT
    / "source_data"
    / "cnmf"
    / "GSE274498_allMG_formal.usages.k_8.dt_0_1.consensus.txt"
)

METADATA_FILE = (
    PROJECT_ROOT
    / "source_data"
    / "GSE274498"
    / "05_cell_metadata_refined_MG_states.csv"
)

OUT_DIR = PROJECT_ROOT / "results" / "cnmf"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# 1. Read frozen inputs
# ------------------------------------------------------------

print("Reading frozen K=8 usage matrix...")

usage = pd.read_csv(
    USAGE_FILE,
    sep="\t",
    index_col=0,
)

print("Reading frozen refined MG metadata...")

meta = pd.read_csv(
    METADATA_FILE
)


programs = [
    "1", "2", "3", "4",
    "5", "6", "7", "8",
]

programs_p2excluded = [
    "1", "3", "4",
    "5", "6", "7", "8",
]


# ------------------------------------------------------------
# 2. Validate input structure
# ------------------------------------------------------------

missing_programs = [
    p for p in programs
    if p not in usage.columns
]

if missing_programs:
    raise ValueError(
        f"Missing K8 usage columns: {missing_programs}"
    )


required_meta = {
    "sample",
    "barcode",
    "age",
    "MG_state",
}

missing_meta = required_meta.difference(
    meta.columns
)

if missing_meta:
    raise ValueError(
        f"Missing metadata columns: {sorted(missing_meta)}"
    )


# cNMF cell IDs are sample + "__" + 10x barcode.
meta["cell_id"] = (
    meta["sample"].astype(str)
    + "__"
    + meta["barcode"].astype(str)
)

if meta["cell_id"].duplicated().any():
    raise ValueError(
        "Duplicated reconstructed cell IDs found in metadata."
    )


print()
print("USAGE CELLS:", len(usage))
print("METADATA CELLS:", len(meta))


if set(usage.index) != set(meta["cell_id"]):
    raise ValueError(
        "Usage and metadata cell sets do not match."
    )


# Reorder metadata to exactly match cNMF usage rows.
meta = (
    meta
    .set_index("cell_id")
    .loc[usage.index]
    .copy()
)


if usage.index.tolist() != meta.index.tolist():
    raise ValueError(
        "Usage and metadata order could not be aligned."
    )


print("CELL-ID MATCH: PASS")


# ------------------------------------------------------------
# 3. Full P1-P8 cell-level normalization
# ------------------------------------------------------------

usage8 = (
    usage[programs]
    .apply(pd.to_numeric, errors="raise")
    .copy()
)

row_sums8 = usage8.sum(axis=1)

if (row_sums8 <= 0).any():
    raise ValueError(
        "At least one cell has non-positive total K8 usage."
    )


usage8_norm = usage8.div(
    row_sums8,
    axis=0,
)


max_sum_error8 = float(
    np.abs(
        usage8_norm.sum(axis=1) - 1
    ).max()
)

print(
    "FULL K8 MAX CELL-SUM ERROR:",
    max_sum_error8,
)


# ------------------------------------------------------------
# 4. Mean normalized usage by MG state
# ------------------------------------------------------------

state_df = usage8_norm.copy()

state_df["MG_state"] = (
    meta["MG_state"].astype(str).values
)

mean_by_state = (
    state_df
    .groupby(
        "MG_state",
        sort=True,
    )[programs]
    .mean()
    .reset_index()
)

mean_by_state.to_csv(
    OUT_DIR
    / "formal_K8_mean_usage_by_MG_state.csv",
    index=False,
)


# ------------------------------------------------------------
# 5. Mean normalized usage by MG state x sample
# ------------------------------------------------------------

state_sample_df = usage8_norm.copy()

state_sample_df["MG_state"] = (
    meta["MG_state"].astype(str).values
)

state_sample_df["sample"] = (
    meta["sample"].astype(str).values
)

mean_by_state_sample = (
    state_sample_df
    .groupby(
        ["MG_state", "sample"],
        sort=True,
    )[programs]
    .mean()
    .reset_index()
)

mean_by_state_sample.to_csv(
    OUT_DIR
    / "formal_K8_mean_usage_by_MG_state_sample.csv",
    index=False,
)


# ------------------------------------------------------------
# 6. P2-excluded cell-level sensitivity normalization
# ------------------------------------------------------------

usage7 = (
    usage[programs_p2excluded]
    .apply(pd.to_numeric, errors="raise")
    .copy()
)

row_sums7 = usage7.sum(axis=1)

if (row_sums7 <= 0).any():
    raise ValueError(
        "At least one cell has non-positive P2-excluded usage."
    )


usage7_norm = usage7.div(
    row_sums7,
    axis=0,
)


max_sum_error7 = float(
    np.abs(
        usage7_norm.sum(axis=1) - 1
    ).max()
)

print(
    "P2-EXCLUDED MAX CELL-SUM ERROR:",
    max_sum_error7,
)


p2_df = usage7_norm.copy()

p2_df["MG_state"] = (
    meta["MG_state"].astype(str).values
)

p2_df["sample"] = (
    meta["sample"].astype(str).values
)

mean_by_state_sample_p2 = (
    p2_df
    .groupby(
        ["MG_state", "sample"],
        sort=True,
    )[programs_p2excluded]
    .mean()
    .reset_index()
)

mean_by_state_sample_p2.to_csv(
    OUT_DIR
    / "formal_K8_mean_usage_by_MG_state_sample_P2excluded.csv",
    index=False,
)


# ------------------------------------------------------------
# 7. Attach age to library-level summaries
# ------------------------------------------------------------

sample_age = (
    meta[
        [
            "sample",
            "age",
        ]
    ]
    .drop_duplicates()
    .copy()
)

age_counts = (
    sample_age
    .groupby("sample")["age"]
    .nunique()
)

if (age_counts > 1).any():
    raise ValueError(
        "A sample maps to more than one age."
    )


full_library = (
    mean_by_state_sample
    .merge(
        sample_age,
        on="sample",
        how="left",
        validate="many_to_one",
    )
)

p2_library = (
    mean_by_state_sample_p2
    .merge(
        sample_age,
        on="sample",
        how="left",
        validate="many_to_one",
    )
)


# ------------------------------------------------------------
# 8. Reproduce formal age-direction sensitivity summary
# ------------------------------------------------------------

targets = [
    (
        "Acinar basal",
        "P4",
        "4",
    ),
    (
        "Orifice",
        "P6",
        "6",
    ),
    (
        "Differentiating meibocyte",
        "P7",
        "7",
    ),
    (
        "Differentiated meibocyte",
        "P8",
        "8",
    ),
]


summary_rows = []

for mg_state, program_name, program_col in targets:

    full_sub = full_library[
        full_library["MG_state"] == mg_state
    ].copy()

    p2_sub = p2_library[
        p2_library["MG_state"] == mg_state
    ].copy()

    full_age = (
        full_sub
        .groupby("age")[program_col]
        .mean()
    )

    p2_age = (
        p2_sub
        .groupby("age")[program_col]
        .mean()
    )

    required_ages = {
        "8W",
        "21M",
    }

    if not required_ages.issubset(
        set(full_age.index.astype(str))
    ):
        raise ValueError(
            f"Missing age group for {mg_state} {program_name}"
        )

    if not required_ages.issubset(
        set(p2_age.index.astype(str))
    ):
        raise ValueError(
            f"Missing P2-excluded age group for {mg_state} {program_name}"
        )

    full_8w = float(
        full_age.loc["8W"]
    )

    full_21m = float(
        full_age.loc["21M"]
    )

    p2_8w = float(
        p2_age.loc["8W"]
    )

    p2_21m = float(
        p2_age.loc["21M"]
    )

    summary_rows.append(
        {
            "MG_state": mg_state,
            "Program": program_name,
            "Full_8W_mean": full_8w,
            "Full_21M_mean": full_21m,
            "Full_delta_21M_minus_8W": (
                full_21m - full_8w
            ),
            "P2excluded_8W_mean": p2_8w,
            "P2excluded_21M_mean": p2_21m,
            "P2excluded_delta_21M_minus_8W": (
                p2_21m - p2_8w
            ),
        }
    )


age_summary = pd.DataFrame(
    summary_rows
)

age_summary.to_csv(
    OUT_DIR
    / "formal_K8_age_direction_sensitivity_summary.csv",
    index=False,
)


# ------------------------------------------------------------
# 9. P7/P8 library-level source table for main figure
# ------------------------------------------------------------

p7 = (
    full_library[
        full_library["MG_state"]
        == "Differentiating meibocyte"
    ][
        [
            "sample",
            "age",
            "7",
        ]
    ]
    .rename(
        columns={
            "7": "normalized_P7_usage"
        }
    )
)

p7["Program"] = "P7"
p7["MG_state"] = (
    "Differentiating meibocyte"
)


p8 = (
    full_library[
        full_library["MG_state"]
        == "Differentiated meibocyte"
    ][
        [
            "sample",
            "age",
            "8",
        ]
    ]
    .rename(
        columns={
            "8": "normalized_P8_usage"
        }
    )
)

p8["Program"] = "P8"
p8["MG_state"] = (
    "Differentiated meibocyte"
)


p7p8 = pd.concat(
    [
        p7,
        p8,
    ],
    ignore_index=True,
)

p7p8 = p7p8[
    [
        "Program",
        "MG_state",
        "sample",
        "age",
        "normalized_P7_usage",
        "normalized_P8_usage",
    ]
]


p7p8.to_csv(
    OUT_DIR
    / "formal_K8_P7_P8_library_age_summary.csv",
    index=False,
)


# ------------------------------------------------------------
# 10. Console validation report
# ------------------------------------------------------------

print()
print("AGE-DIRECTION SUMMARY:")
print(
    age_summary.to_string(
        index=False
    )
)

print()

p7_row = age_summary[
    age_summary["Program"] == "P7"
].iloc[0]

p8_row = age_summary[
    age_summary["Program"] == "P8"
].iloc[0]

print(
    "P7 FULL delta 21M-8W:",
    p7_row["Full_delta_21M_minus_8W"],
)

print(
    "P7 P2-EXCLUDED delta 21M-8W:",
    p7_row["P2excluded_delta_21M_minus_8W"],
)

print(
    "P8 FULL delta 21M-8W:",
    p8_row["Full_delta_21M_minus_8W"],
)

print(
    "P8 P2-EXCLUDED delta 21M-8W:",
    p8_row["P2excluded_delta_21M_minus_8W"],
)


if (
    p7_row["Full_delta_21M_minus_8W"] <= 0
    or
    p7_row["P2excluded_delta_21M_minus_8W"] <= 0
    or
    p8_row["Full_delta_21M_minus_8W"] <= 0
    or
    p8_row["P2excluded_delta_21M_minus_8W"] <= 0
):
    raise ValueError(
        "Expected positive aging direction for P7/P8 was not reproduced."
    )


print()
print("P7/P8 AGE-DIRECTION CHECK: PASS")
print("DONE")
print("Results saved to:", OUT_DIR)