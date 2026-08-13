#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
09_GSE17822_GSEA_zeroVariance_sensitivity.py

Sensitivity analysis for human GSE17822 GSEA.

For both representative-probe rules:
1. Remove genes whose selected probe has exactly zero variance
   across all 12 samples.
2. Rank remaining genes by the age/sex-adjusted MGD t statistic.
3. Add an extremely small seed-specific jitter (1e-12) only to
   break exact ranking ties.
4. Run GSEApy prerank for seeds 1-10.

This is a sensitivity analysis and does NOT replace the frozen
manuscript GSEA outputs.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp


PROJECT_ROOT = Path(__file__).resolve().parents[1]

SERIES_MATRIX = (
    PROJECT_ROOT
    / "source_data"
    / "GSE17822"
    / "GSE17822_series_matrix.txt.gz"
)

GMT_FILE = (
    PROJECT_ROOT
    / "source_data"
    / "GSE17822"
    / "gene_sets.gmt"
)

RESULT_DIR = (
    PROJECT_ROOT
    / "results"
    / "human"
)

RESULT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


RULES = {
    "representativeProbe":
        RESULT_DIR
        / "GSE17822_gene_DE_ageSex_adjusted_representativeProbe.csv",

    "maxLog2VarianceProbe":
        RESULT_DIR
        / "GSE17822_gene_DE_ageSex_adjusted_maxLog2VarianceProbe.csv",
}


PERMUTATIONS = 1000
MIN_SIZE = 10
MAX_SIZE = 2000


def main():

    print(
        "============================================"
    )

    print(
        "GSE17822 ZERO-VARIANCE FILTERED GSEA"
    )

    print(
        "============================================"
    )

    # --------------------------------------------------------
    # Read expression matrix and calculate probe variance
    # --------------------------------------------------------

    expr = pd.read_csv(
        SERIES_MATRIX,
        sep="\t",
        comment="!",
        index_col=0,
        compression="gzip",
    )

    expr.index = (
        expr.index
        .astype(str)
        .str.replace(
            '"',
            "",
            regex=False
        )
    )

    expr = expr.astype(float)

    log2_expr = np.log2(expr)

    probe_variance = (
        log2_expr.var(
            axis=1,
            ddof=1
        )
    )

    print()

    print(
        "Total probes:",
        len(probe_variance)
    )

    print(
        "Zero-variance probes:",
        int(
            probe_variance.eq(0).sum()
        )
    )

    all_results = []

    # --------------------------------------------------------
    # Analyze both representative-probe rules
    # --------------------------------------------------------

    for rule_name, input_file in RULES.items():

        print()

        print(
            "--------------------------------------------"
        )

        print(
            "RULE:",
            rule_name
        )

        gene_de = pd.read_csv(
            input_file
        )

        gene_de[
            "probe_variance"
        ] = gene_de[
            "probe_id"
        ].map(
            probe_variance
        )

        total_genes = len(
            gene_de
        )

        zero_variance_n = int(
            gene_de[
                "probe_variance"
            ].eq(0).sum()
        )

        ranked = gene_de.loc[
            gene_de[
                "probe_variance"
            ].gt(0),
            [
                "Symbol",
                "t",
            ]
        ].copy()

        ranked = ranked.dropna()

        ranked[
            "t"
        ] = pd.to_numeric(
            ranked["t"],
            errors="coerce"
        )

        ranked = ranked.dropna()

        print(
            "Total genes:",
            total_genes
        )

        print(
            "Zero-variance genes removed:",
            zero_variance_n
        )

        print(
            "Genes retained for GSEA:",
            len(ranked)
        )

        # ----------------------------------------------------
        # 10-seed sensitivity
        # ----------------------------------------------------

        for seed in range(
            1,
            11
        ):

            rnk = ranked.copy()

            rng = np.random.default_rng(
                seed
            )

            # Tiny jitter only breaks exact ties.
            rnk[
                "rank_stat"
            ] = (
                rnk["t"].to_numpy(
                    dtype=float
                )
                +
                rng.normal(
                    loc=0.0,
                    scale=1e-12,
                    size=len(rnk)
                )
            )

            rnk = (
                rnk[
                    [
                        "Symbol",
                        "rank_stat",
                    ]
                ]
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
                        "Probe_rule":
                            rule_name,

                        "Seed":
                            seed,

                        "Total_gene_level_rows":
                            total_genes,

                        "Zero_variance_removed":
                            zero_variance_n,

                        "Genes_ranked":
                            len(ranked),

                        "Term":
                            row["Term"],

                        "ES":
                            row["ES"],

                        "NES":
                            row["NES"],

                        "NOM_p":
                            row["NOM p-val"],

                        "FDR_q":
                            row["FDR q-val"],
                    }
                )

            print(
                "seed",
                seed,
                "done"
            )

    # --------------------------------------------------------
    # Save all seed-level results
    # --------------------------------------------------------

    all_df = pd.DataFrame(
        all_results
    )

    numeric_columns = [
        "ES",
        "NES",
        "NOM_p",
        "FDR_q",
    ]

    for column in numeric_columns:

        all_df[
            column
        ] = pd.to_numeric(
            all_df[
                column
            ],
            errors="coerce",
        )

    all_output = (
        RESULT_DIR
        / "GSE17822_GSEA_zeroVarianceFiltered_10seed.csv"
    )

    all_df.to_csv(
        all_output,
        index=False,
    )

    # --------------------------------------------------------
    # Summary across 10 seeds
    # --------------------------------------------------------

    summary_rows = []

    for (
        rule,
        term
    ), group in all_df.groupby(
        [
            "Probe_rule",
            "Term",
        ],
        sort=False,
    ):

        summary_rows.append(
            {
                "Probe_rule":
                    rule,

                "Term":
                    term,

                "Seeds_n":
                    len(group),

                "Genes_ranked":
                    int(
                        group[
                            "Genes_ranked"
                        ].iloc[0]
                    ),

                "Zero_variance_removed":
                    int(
                        group[
                            "Zero_variance_removed"
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

                "NOM_p_max":
                    group[
                        "NOM_p"
                    ].max(),

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

    summary_output = (
        RESULT_DIR
        / "GSE17822_GSEA_zeroVarianceFiltered_summary.csv"
    )

    summary.to_csv(
        summary_output,
        index=False,
    )

    print()

    print(
        "============================================"
    )

    print(
        "10-SEED SUMMARY"
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
        all_output
    )

    print(
        summary_output
    )

    print()

    print(
        "DONE"
    )


if __name__ == "__main__":
    main()