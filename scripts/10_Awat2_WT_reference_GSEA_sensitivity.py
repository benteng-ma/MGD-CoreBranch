#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
10_Awat2_WT_reference_GSEA_sensitivity.py

Awat2 WT-reference sensitivity and WT internal negative-control GSEA.

Input:
    results/cross_model/Awat2_wholeMG_allGene_pseudobulk_CPM.csv
    source_data/GSE261036/Awat2_defense_gene_sets.gmt

Comparisons:
    KO vs WT1
    KO vs WT2
    KO vs WT3
    KO vs WTmean
    WT3 vs WT1
    WT3 vs WT2

Ranking statistic:
    log2((CPM_A + 1) / (CPM_B + 1))

Per-comparison expressed-gene filter:
    max(CPM_A, CPM_B) >= 1

Each comparison is rerun with seeds 1-10 using GSEApy prerank.

This script is a reproducibility/sensitivity analysis.
The frozen original manuscript GSEA result tables are retained separately
under results/awat2/.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

CPM_FILE = (
    PROJECT_ROOT
    / "results"
    / "cross_model"
    / "Awat2_wholeMG_allGene_pseudobulk_CPM.csv"
)

GMT_FILE = (
    PROJECT_ROOT
    / "source_data"
    / "GSE261036"
    / "Awat2_defense_gene_sets.gmt"
)

OUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "awat2"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# GSEA parameters
# ============================================================

PERMUTATIONS = 1000

MIN_SIZE = 10

MAX_SIZE = 2000


# ============================================================
# Comparisons
# ============================================================

COMPARISONS = [
    (
        "KO_vs_WT1",
        "KO",
        "WT1",
    ),
    (
        "KO_vs_WT2",
        "KO",
        "WT2",
    ),
    (
        "KO_vs_WT3",
        "KO",
        "WT3",
    ),
    (
        "KO_vs_WTmean",
        "KO",
        "WTmean",
    ),
    (
        "WT3_vs_WT1",
        "WT3",
        "WT1",
    ),
    (
        "WT3_vs_WT2",
        "WT3",
        "WT2",
    ),
]


# ============================================================
# Main
# ============================================================

def main():

    print(
        "============================================"
    )

    print(
        "AWAT2 WT-REFERENCE GSEA SENSITIVITY"
    )

    print(
        "============================================"
    )

    print()

    # --------------------------------------------------------
    # Check files
    # --------------------------------------------------------

    if not CPM_FILE.exists():

        raise FileNotFoundError(
            CPM_FILE
        )

    if not GMT_FILE.exists():

        raise FileNotFoundError(
            GMT_FILE
        )

    # --------------------------------------------------------
    # Read whole-MG pseudobulk CPM
    # --------------------------------------------------------

    cpm = pd.read_csv(
        CPM_FILE
    )

    required = {
        "gene",
        "WT1",
        "WT2",
        "WT3",
        "KO",
    }

    missing = (
        required
        - set(
            cpm.columns
        )
    )

    if missing:

        raise RuntimeError(
            "Missing CPM columns: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

    # WTmean is the arithmetic mean CPM
    # across the three WT libraries.

    cpm[
        "WTmean"
    ] = (
        cpm[
            [
                "WT1",
                "WT2",
                "WT3",
            ]
        ]
        .mean(
            axis=1
        )
    )

    print(
        "Whole-MG CPM genes:",
        len(cpm)
    )

    print()

    all_results = []

    rank_summary = []

    # --------------------------------------------------------
    # Process each comparison
    # --------------------------------------------------------

    for (
        comparison,
        numerator,
        denominator,
    ) in COMPARISONS:

        print(
            "--------------------------------------------"
        )

        print(
            "COMPARISON:",
            comparison
        )

        # ----------------------------------------------------
        # Per-comparison expressed-gene filter
        # ----------------------------------------------------

        expressed = np.maximum(
            cpm[
                numerator
            ].to_numpy(
                dtype=float
            ),
            cpm[
                denominator
            ].to_numpy(
                dtype=float
            ),
        ) >= 1.0

        ranked = cpm.loc[
            expressed,
            [
                "gene",
                numerator,
                denominator,
            ]
        ].copy()

        # ----------------------------------------------------
        # Ranking statistic
        # ----------------------------------------------------

        ranked[
            "rank_stat"
        ] = np.log2(
            (
                ranked[
                    numerator
                ].to_numpy(
                    dtype=float
                )
                + 1.0
            )
            /
            (
                ranked[
                    denominator
                ].to_numpy(
                    dtype=float
                )
                + 1.0
            )
        )

        # Mouse symbols are converted to uppercase
        # to match the frozen GO GMT gene symbols.

        ranked[
            "gene_GSEA"
        ] = (
            ranked[
                "gene"
            ]
            .astype(str)
            .str.upper()
        )

        # Remove missing/empty symbols.

        ranked = ranked.loc[
            ranked[
                "gene_GSEA"
            ].ne("")
        ].copy()

        ranked = ranked.loc[
            ranked[
                "gene_GSEA"
            ].ne(
                "NAN"
            )
        ].copy()

        # If uppercasing produces duplicate symbols,
        # keep the row with the largest absolute ranking statistic.

        ranked[
            "abs_rank"
        ] = ranked[
            "rank_stat"
        ].abs()

        ranked = (
            ranked
            .sort_values(
                [
                    "gene_GSEA",
                    "abs_rank",
                ],
                ascending=[
                    True,
                    False,
                ],
                kind="mergesort",
            )
            .drop_duplicates(
                "gene_GSEA",
                keep="first",
            )
        )

        ranked = (
            ranked
            .sort_values(
                "rank_stat",
                ascending=False,
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )

        rank_n = len(
            ranked
        )

        print(
            "rank genes:",
            rank_n
        )

        rank_summary.append(
            {
                "Comparison":
                    comparison,

                "Numerator":
                    numerator,

                "Denominator":
                    denominator,

                "Rank_genes":
                    rank_n,
            }
        )

        # ----------------------------------------------------
        # 10-seed GSEA
        # ----------------------------------------------------

        for seed in range(
            1,
            11
        ):

            rnk = ranked[
                [
                    "gene_GSEA",
                    "rank_stat",
                ]
            ].copy()

            # Tiny seed-specific jitter is used only
            # to resolve exact ranking ties.

            rng = np.random.default_rng(
                seed
            )

            rnk[
                "rank_stat"
            ] = (
                rnk[
                    "rank_stat"
                ].to_numpy(
                    dtype=float
                )
                +
                rng.normal(
                    loc=0.0,
                    scale=1e-12,
                    size=len(
                        rnk
                    ),
                )
            )

            rnk = (
                rnk
                .sort_values(
                    "rank_stat",
                    ascending=False,
                    kind="mergesort",
                )
            )

            pre = gp.prerank(
                rnk=rnk,
                gene_sets=str(
                    GMT_FILE
                ),
                min_size=MIN_SIZE,
                max_size=MAX_SIZE,
                permutation_num=PERMUTATIONS,
                weight=1.0,
                ascending=False,
                threads=1,
                outdir=None,
                no_plot=True,
                seed=seed,
                verbose=False,
            )

            result = (
                pre.res2d
                .reset_index(
                    drop=True
                )
            )

            for _, row in result.iterrows():

                all_results.append(
                    {
                        "Comparison":
                            comparison,

                        "Seed":
                            seed,

                        "Rank_genes":
                            rank_n,

                        "Term":
                            row[
                                "Term"
                            ],

                        "ES":
                            row[
                                "ES"
                            ],

                        "NES":
                            row[
                                "NES"
                            ],

                        "NOM_p":
                            row[
                                "NOM p-val"
                            ],

                        "FDR_q":
                            row[
                                "FDR q-val"
                            ],
                    }
                )

            print(
                "seed",
                seed,
                "done"
            )

        print()

    # ========================================================
    # Save rank-gene counts
    # ========================================================

    rank_summary_df = pd.DataFrame(
        rank_summary
    )

    rank_summary_file = (
        OUT_DIR
        / "Awat2_WT_reference_GSEA_rank_gene_counts.csv"
    )

    rank_summary_df.to_csv(
        rank_summary_file,
        index=False,
    )

    # ========================================================
    # Save all seed-level GSEA results
    # ========================================================

    all_df = pd.DataFrame(
        all_results
    )

    for column in [
        "ES",
        "NES",
        "NOM_p",
        "FDR_q",
    ]:

        all_df[
            column
        ] = pd.to_numeric(
            all_df[
                column
            ],
            errors="coerce",
        )

    all_file = (
        OUT_DIR
        / "Awat2_WT_reference_GSEA_10seed.csv"
    )

    all_df.to_csv(
        all_file,
        index=False,
    )

    # ========================================================
    # 10-seed summary
    # ========================================================

    summary_rows = []

    for (
        comparison,
        term,
    ), group in all_df.groupby(
        [
            "Comparison",
            "Term",
        ],
        sort=False,
    ):

        summary_rows.append(
            {
                "Comparison":
                    comparison,

                "Term":
                    term,

                "Seeds_n":
                    len(
                        group
                    ),

                "Rank_genes":
                    int(
                        group[
                            "Rank_genes"
                        ].iloc[0]
                    ),

                "NES_min":
                    group[
                        "NES"
                    ].min(),

                "NES_max":
                    group[
                        "NES"
                    ].max(),

                "NOM_p_min":
                    group[
                        "NOM_p"
                    ].min(),

                "NOM_p_max":
                    group[
                        "NOM_p"
                    ].max(),

                "FDR_q_min":
                    group[
                        "FDR_q"
                    ].min(),

                "FDR_q_max":
                    group[
                        "FDR_q"
                    ].max(),

                "FDR_lt_0.05_all_10_seeds":
                    bool(
                        (
                            group[
                                "FDR_q"
                            ]
                            < 0.05
                        ).all()
                    ),
            }
        )

    summary = pd.DataFrame(
        summary_rows
    )

    summary_file = (
        OUT_DIR
        / "Awat2_WT_reference_GSEA_10seed_summary.csv"
    )

    summary.to_csv(
        summary_file,
        index=False,
    )

    # ========================================================
    # Print final summary
    # ========================================================

    print(
        "============================================"
    )

    print(
        "RANK-GENE COUNTS"
    )

    print(
        "============================================"
    )

    print()

    print(
        rank_summary_df.to_string(
            index=False
        )
    )

    print()

    print(
        "============================================"
    )

    print(
        "10-SEED GSEA SUMMARY"
    )

    print(
        "============================================"
    )

    print()

    print(
        summary.to_string(
            index=False
        )
    )

    print()

    print(
        "Saved:"
    )

    print(
        rank_summary_file
    )

    print(
        all_file
    )

    print(
        summary_file
    )

    print()

    print(
        "INTERPRETATION NOTE:"
    )

    print(
        "The frozen original manuscript GSEA tables "
        "remain the authoritative record of the "
        "original analysis. This 10-seed rerun "
        "evaluates robustness to WT reference choice "
        "and ranking ties using the frozen release CPM."
    )

    print()

    print(
        "DONE"
    )


if __name__ == "__main__":
    main()