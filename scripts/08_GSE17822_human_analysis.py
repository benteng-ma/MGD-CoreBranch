#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
08_GSE17822_human_analysis.py

Reproduce the validated core GSE17822 human analysis:

1. Read GEO series-matrix expression and sample metadata
2. log2-transform the processed intensities
3. Probe-level OLS:
       log2(expression) ~ MGD + age + Male
4. BH-FDR across all probes
5. Read GPL6947 probe-to-gene annotation from family SOFT
6. Select one representative probe per gene by highest mean log2 expression
7. Sensitivity selection: highest log2-expression variance probe per gene
8. Reconstruct the fixed 16-gene human DefenseScore
9. Fit:
       DefenseScore ~ MGD + age + Male

Run:
    python scripts/08_GSE17822_human_analysis.py
"""

from pathlib import Path
import csv
import gzip

import numpy as np
import pandas as pd
import statsmodels.api as sm

from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_DIR = PROJECT_ROOT / "source_data" / "GSE17822"
OUTPUT_DIR = PROJECT_ROOT / "results" / "human"

SERIES_MATRIX = INPUT_DIR / "GSE17822_series_matrix.txt.gz"
FAMILY_SOFT = INPUT_DIR / "GSE17822_family.soft.gz"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Fixed mouse-derived defense candidate set
# ============================================================

DEFENSE_CANDIDATES = [
    "CXCL5",
    "LCN2",
    "IL1RN",
    "IFITM1",
    "IFITM3",
    "CXCL1",
    "CXCL2",
    "CCL20",
    "SLPI",
    "S100A8",
    "IFI47",
    "ZBP1",
    "ISG15",
    "RTP4",
    "CXCL10",
    "PLAC8",
    "WDFY1",
]


# ============================================================
# Read sample metadata from GEO series matrix
# ============================================================

def read_sample_metadata(path):

    accessions = None
    titles = None
    genders = None
    ages = None

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        for raw in fh:

            line = raw.rstrip("\n")

            if line.startswith("!Sample_geo_accession"):

                fields = next(
                    csv.reader(
                        [line],
                        delimiter="\t"
                    )
                )[1:]

                accessions = [
                    x.strip().strip('"')
                    for x in fields
                ]

            elif line.startswith("!Sample_title"):

                fields = next(
                    csv.reader(
                        [line],
                        delimiter="\t"
                    )
                )[1:]

                titles = [
                    x.strip().strip('"')
                    for x in fields
                ]

            elif line.startswith(
                "!Sample_characteristics_ch1"
            ):

                fields = next(
                    csv.reader(
                        [line],
                        delimiter="\t"
                    )
                )[1:]

                fields = [
                    x.strip().strip('"')
                    for x in fields
                ]

                if (
                    fields
                    and all(
                        x.lower().startswith("gender:")
                        for x in fields
                    )
                ):

                    genders = [
                        x.split(":", 1)[1].strip()
                        for x in fields
                    ]

                elif (
                    fields
                    and all(
                        x.lower().startswith("age:")
                        for x in fields
                    )
                ):

                    ages = [
                        int(
                            float(
                                x.split(":", 1)[1].strip()
                            )
                        )
                        for x in fields
                    ]

    if accessions is None:
        raise RuntimeError(
            "Sample accessions not found."
        )

    if titles is None:
        raise RuntimeError(
            "Sample titles not found."
        )

    if genders is None:
        raise RuntimeError(
            "Sample genders not found."
        )

    if ages is None:
        raise RuntimeError(
            "Sample ages not found."
        )

    mgd = []

    for title in titles:

        low = title.lower()

        if "normal" in low:
            mgd.append(0)

        elif "diseased" in low:
            mgd.append(1)

        else:
            raise RuntimeError(
                "Cannot determine disease state "
                f"from title: {title}"
            )

    male = []

    for gender in genders:

        low = gender.lower()

        if low == "male":
            male.append(1)

        elif low == "female":
            male.append(0)

        else:
            raise RuntimeError(
                f"Unexpected gender: {gender}"
            )

    metadata = pd.DataFrame(
        {
            "sample": accessions,
            "title": titles,
            "MGD": mgd,
            "age": ages,
            "Male": male,
        }
    )

    return metadata


# ============================================================
# Read expression matrix
# ============================================================

def read_expression(path):

    expr = pd.read_csv(
        path,
        sep="\t",
        comment="!",
        index_col=0,
        compression="gzip",
    )

    expr.index = (
        expr.index
        .astype(str)
        .str.replace('"', "", regex=False)
    )

    expr.columns = (
        expr.columns
        .astype(str)
        .str.replace('"', "", regex=False)
    )

    expr = expr.apply(
        pd.to_numeric,
        errors="raise"
    )

    if (expr <= 0).any().any():

        raise RuntimeError(
            "Non-positive expression values found."
        )

    log2_expr = np.log2(
        expr.astype(float)
    )

    log2_expr.index.name = "probe_id"

    return log2_expr


# ============================================================
# Probe-level OLS
# ============================================================

def run_probe_de(log2_expr, metadata):

    sample_order = metadata[
        "sample"
    ].tolist()

    if set(sample_order) != set(
        log2_expr.columns
    ):

        raise RuntimeError(
            "Expression samples do not match metadata."
        )

    log2_expr = log2_expr.loc[
        :,
        sample_order
    ]

    design = sm.add_constant(
        metadata[
            [
                "MGD",
                "age",
                "Male",
            ]
        ].astype(float),
        has_constant="add",
    )

    X = design.to_numpy(
        dtype=float
    )

    Y = log2_expr.to_numpy(
        dtype=float
    ).T

    # Ordinary least squares for all probes
    XtX_inv = np.linalg.inv(
        X.T @ X
    )

    beta = (
        XtX_inv
        @ X.T
        @ Y
    )

    fitted = X @ beta

    residuals = (
        Y - fitted
    )

    n = X.shape[0]
    p = X.shape[1]

    df_resid = n - p

    sigma2 = (
        residuals ** 2
    ).sum(axis=0) / df_resid

    mgd_index = list(
        design.columns
    ).index("MGD")

    mgd_beta = beta[
        mgd_index,
        :
    ]

    mgd_se = np.sqrt(
        sigma2
        * XtX_inv[
            mgd_index,
            mgd_index
        ]
    )

    mgd_t = (
        mgd_beta
        / mgd_se
    )

    mgd_p = (
        2
        * student_t.sf(
            np.abs(mgd_t),
            df=df_resid
        )
    )

    mgd_fdr = multipletests(
        mgd_p,
        method="fdr_bh"
    )[1]

    normal_samples = metadata.loc[
        metadata["MGD"] == 0,
        "sample"
    ].tolist()

    mgd_samples = metadata.loc[
        metadata["MGD"] == 1,
        "sample"
    ].tolist()

    result = pd.DataFrame(
        {
            "probe_id":
                log2_expr.index,

            "MGD_coef_ageSexAdj":
                mgd_beta,

            "SE":
                mgd_se,

            "t":
                mgd_t,

            "pvalue":
                mgd_p,

            "FDR_BH":
                mgd_fdr,

            "Normal_mean_log2":
                log2_expr[
                    normal_samples
                ].mean(
                    axis=1
                ).to_numpy(),

            "MGD_mean_log2":
                log2_expr[
                    mgd_samples
                ].mean(
                    axis=1
                ).to_numpy(),
        }
    )

    return result


# ============================================================
# Read GPL6947 annotation from family SOFT
# ============================================================

def read_platform_annotation(path):

    in_table = False
    header = None
    rows = []

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        for raw in fh:

            line = raw.rstrip("\n")

            if line.startswith(
                "!platform_table_begin"
            ):

                in_table = True
                header = None
                continue

            if line.startswith(
                "!platform_table_end"
            ):

                break

            if not in_table:
                continue

            fields = line.split("\t")

            if header is None:

                header = fields
                continue

            if not fields:
                continue

            if not fields[0].startswith(
                "ILMN_"
            ):
                continue

            if len(fields) < len(header):

                fields = fields + (
                    [""] *
                    (
                        len(header)
                        - len(fields)
                    )
                )

            row = dict(
                zip(
                    header,
                    fields
                )
            )

            rows.append(
                (
                    row.get(
                        "ID",
                        ""
                    ),
                    row.get(
                        "Symbol",
                        ""
                    ),
                )
            )

    annotation = pd.DataFrame(
        rows,
        columns=[
            "probe_id",
            "Symbol",
        ]
    )

    annotation[
        "probe_id"
    ] = (
        annotation[
            "probe_id"
        ]
        .astype(str)
        .str.strip()
    )

    annotation[
        "Symbol"
    ] = (
        annotation[
            "Symbol"
        ]
        .astype(str)
        .str.strip()
    )

    annotation.loc[
        annotation[
            "Symbol"
        ].isin(
            [
                "",
                "---",
                "NA",
                "nan",
            ]
        ),
        "Symbol"
    ] = np.nan

    annotation = (
        annotation
        .drop_duplicates(
            "probe_id",
            keep="first"
        )
    )

    return annotation


# ============================================================
# Representative probe selection
# ============================================================

def select_representative_probes(
    annotated_de,
    log2_expr
):

    table = annotated_de.copy()

    table = table[
        table["Symbol"].notna()
    ].copy()

    # With six Normal and six MGD samples,
    # this equals the mean across all 12 samples.
    table[
        "All_mean_log2"
    ] = (
        table[
            "Normal_mean_log2"
        ]
        +
        table[
            "MGD_mean_log2"
        ]
    ) / 2.0

    variance = (
        log2_expr.var(
            axis=1,
            ddof=1
        )
        .rename(
            "log2_variance"
        )
    )

    table = table.merge(
        variance,
        left_on="probe_id",
        right_index=True,
        how="left",
    )

    table[
        "_original_order"
    ] = np.arange(
        len(table)
    )

    # ----------------------------------
    # Rule 1:
    # highest mean log2 expression
    # ----------------------------------

    mean_probe = (
        table
        .sort_values(
            [
                "Symbol",
                "All_mean_log2",
                "_original_order",
            ],
            ascending=[
                True,
                False,
                True,
            ],
            kind="mergesort",
        )
        .drop_duplicates(
            "Symbol",
            keep="first"
        )
        .sort_values(
            "Symbol",
            kind="mergesort"
        )
        .copy()
    )

    # ----------------------------------
    # Rule 2:
    # highest log2 variance
    # ----------------------------------

    maxvar_probe = (
        table
        .sort_values(
            [
                "Symbol",
                "log2_variance",
                "_original_order",
            ],
            ascending=[
                True,
                False,
                True,
            ],
            kind="mergesort",
        )
        .drop_duplicates(
            "Symbol",
            keep="first"
        )
        .sort_values(
            "Symbol",
            kind="mergesort"
        )
        .copy()
    )

    output_columns = [
        "probe_id",
        "MGD_coef_ageSexAdj",
        "SE",
        "t",
        "pvalue",
        "FDR_BH",
        "Normal_mean_log2",
        "MGD_mean_log2",
        "Symbol",
        "All_mean_log2",
    ]

    mean_output = (
        mean_probe[
            output_columns
        ]
        .reset_index(
            drop=True
        )
    )

    maxvar_output = (
        maxvar_probe[
            output_columns
        ]
        .reset_index(
            drop=True
        )
    )

    comparison = (
        mean_probe[
            [
                "Symbol",
                "probe_id",
                "All_mean_log2",
            ]
        ]
        .rename(
            columns={
                "probe_id":
                    "meanProbe"
            }
        )
    )

    comparison = comparison.merge(
        maxvar_probe[
            [
                "Symbol",
                "probe_id",
                "log2_variance",
            ]
        ].rename(
            columns={
                "probe_id":
                    "maxVar_probe"
            }
        ),
        on="Symbol",
        how="inner",
    )

    comparison[
        "same_probe"
    ] = (
        comparison[
            "meanProbe"
        ]
        ==
        comparison[
            "maxVar_probe"
        ]
    )

    comparison = comparison[
        [
            "Symbol",
            "meanProbe",
            "All_mean_log2",
            "maxVar_probe",
            "log2_variance",
            "same_probe",
        ]
    ]

    return (
        mean_output,
        maxvar_output,
        comparison,
    )


# ============================================================
# DefenseScore
# ============================================================

def build_defense_score(
    representative_de,
    log2_expr,
    metadata
):

    rep = (
        representative_de
        .set_index(
            "Symbol"
        )
    )

    present = [
        gene
        for gene in DEFENSE_CANDIDATES
        if gene in rep.index
    ]

    missing = [
        gene
        for gene in DEFENSE_CANDIDATES
        if gene not in rep.index
    ]

    selected = (
        rep.loc[
            present,
            [
                "probe_id"
            ]
        ]
        .reset_index()
    )

    expression = (
        log2_expr.loc[
            selected[
                "probe_id"
            ].tolist()
        ]
        .copy()
    )

    expression.index = (
        selected[
            "Symbol"
        ].tolist()
    )

    # Important:
    # sample SD = ddof=1
    gene_mean = expression.mean(
        axis=1
    )

    gene_sd = expression.std(
        axis=1,
        ddof=1
    )

    zscore = (
        expression
        .sub(
            gene_mean,
            axis=0
        )
        .div(
            gene_sd,
            axis=0
        )
    )

    defense_score = (
        zscore.mean(
            axis=0
        )
    )

    score_table = metadata[
        [
            "sample",
            "MGD",
            "age",
            "Male",
        ]
    ].copy()

    score_table[
        "DefenseScore"
    ] = (
        defense_score
        .reindex(
            score_table[
                "sample"
            ]
        )
        .to_numpy()
    )

    X = sm.add_constant(
        score_table[
            [
                "MGD",
                "age",
                "Male",
            ]
        ].astype(float),
        has_constant="add",
    )

    model = sm.OLS(
        score_table[
            "DefenseScore"
        ],
        X,
    ).fit()

    ci = model.conf_int().loc[
        "MGD"
    ]

    summary = pd.DataFrame(
        [
            {
                "candidate_genes_n":
                    len(
                        DEFENSE_CANDIDATES
                    ),

                "measurable_genes_n":
                    len(present),

                "measurable_genes":
                    ";".join(
                        present
                    ),

                "missing_genes":
                    ";".join(
                        missing
                    ),

                "MGD_beta_ageSexAdj":
                    model.params[
                        "MGD"
                    ],

                "SE":
                    model.bse[
                        "MGD"
                    ],

                "t":
                    model.tvalues[
                        "MGD"
                    ],

                "pvalue":
                    model.pvalues[
                        "MGD"
                    ],

                "CI95_low":
                    ci.iloc[0],

                "CI95_high":
                    ci.iloc[1],
            }
        ]
    )

    return (
        score_table,
        summary,
        selected,
    )


# ============================================================
# Main
# ============================================================

def main():

    print(
        "========================================"
    )

    print(
        "GSE17822 HUMAN ANALYSIS"
    )

    print(
        "========================================"
    )

    print()

    if not SERIES_MATRIX.exists():

        raise FileNotFoundError(
            SERIES_MATRIX
        )

    if not FAMILY_SOFT.exists():

        raise FileNotFoundError(
            FAMILY_SOFT
        )

    # ----------------------------------
    # Metadata
    # ----------------------------------

    print(
        "1. Reading sample metadata..."
    )

    metadata = read_sample_metadata(
        SERIES_MATRIX
    )

    metadata.to_csv(
        OUTPUT_DIR
        / "GSE17822_sample_metadata.csv",
        index=False,
    )

    print(
        "   samples:",
        len(metadata)
    )

    print(
        "   Normal:",
        int(
            (
                metadata[
                    "MGD"
                ]
                == 0
            ).sum()
        )
    )

    print(
        "   MGD:",
        int(
            (
                metadata[
                    "MGD"
                ]
                == 1
            ).sum()
        )
    )

    print()

    # ----------------------------------
    # Expression
    # ----------------------------------

    print(
        "2. Reading series-matrix expression..."
    )

    log2_expr = read_expression(
        SERIES_MATRIX
    )

    print(
        "   expression shape:",
        log2_expr.shape
    )

    # ----------------------------------
    # Probe-level DE
    # ----------------------------------

    print()

    print(
        "3. Running probe-level "
        "age/sex-adjusted OLS..."
    )

    probe_de = run_probe_de(
        log2_expr,
        metadata
    )

    probe_de.to_csv(
        OUTPUT_DIR
        / "GSE17822_probe_DE_ageSex_adjusted.csv",
        index=False,
    )

    print(
        "   probe rows:",
        len(probe_de)
    )

    # ----------------------------------
    # Annotation
    # ----------------------------------

    print()

    print(
        "4. Reading GPL6947 annotation..."
    )

    annotation = (
        read_platform_annotation(
            FAMILY_SOFT
        )
    )

    annotated = (
        probe_de.merge(
            annotation,
            on="probe_id",
            how="left",
        )
    )

    annotated.to_csv(
        OUTPUT_DIR
        / "GSE17822_probe_DE_ageSex_adjusted_annotated.csv",
        index=False,
    )

    # ----------------------------------
    # Representative probes
    # ----------------------------------

    print()

    print(
        "5. Selecting representative probes..."
    )

    (
        representative,
        maxvar,
        comparison,
    ) = select_representative_probes(
        annotated,
        log2_expr,
    )

    representative.to_csv(
        OUTPUT_DIR
        / "GSE17822_gene_DE_ageSex_adjusted_representativeProbe.csv",
        index=False,
    )

    maxvar.to_csv(
        OUTPUT_DIR
        / "GSE17822_gene_DE_ageSex_adjusted_maxLog2VarianceProbe.csv",
        index=False,
    )

    comparison.to_csv(
        OUTPUT_DIR
        / "GSE17822_meanProbe_vs_maxLog2VarianceProbe.csv",
        index=False,
    )

    print(
        "   representative genes:",
        len(representative)
    )

    same_n = int(
        comparison[
            "same_probe"
        ].sum()
    )

    total_n = len(
        comparison
    )

    print(
        "   same probe under both rules:",
        f"{same_n}/{total_n}",
        f"({same_n / total_n:.2%})",
    )

    # ----------------------------------
    # DefenseScore
    # ----------------------------------

    print()

    print(
        "6. Reconstructing DefenseScore..."
    )

    (
        defense_score,
        defense_summary,
        defense_probes,
    ) = build_defense_score(
        representative,
        log2_expr,
        metadata,
    )

    defense_score.to_csv(
        OUTPUT_DIR
        / "GSE17822_mouseDefense_moduleScore.csv",
        index=False,
    )

    defense_summary.to_csv(
        OUTPUT_DIR
        / "GSE17822_mouseDefense_moduleScore_summary.csv",
        index=False,
    )

    defense_probes.to_csv(
        OUTPUT_DIR
        / "GSE17822_mouseDefense_selected_probes.csv",
        index=False,
    )

    # ----------------------------------
    # Validation output
    # ----------------------------------

    print()

    print(
        "========================================"
    )

    print(
        "VALIDATION"
    )

    print(
        "========================================"
    )

    first = probe_de.iloc[0]

    print()

    print(
        "First probe:"
    )

    print(
        "probe =",
        first[
            "probe_id"
        ]
    )

    print(
        "MGD coef =",
        first[
            "MGD_coef_ageSexAdj"
        ]
    )

    print(
        "SE =",
        first[
            "SE"
        ]
    )

    print(
        "t =",
        first[
            "t"
        ]
    )

    print(
        "p =",
        first[
            "pvalue"
        ]
    )

    print()

    s = defense_summary.iloc[0]

    print(
        "Defense module:"
    )

    print(
        "candidate genes =",
        int(
            s[
                "candidate_genes_n"
            ]
        )
    )

    print(
        "measurable genes =",
        int(
            s[
                "measurable_genes_n"
            ]
        )
    )

    print(
        "missing genes =",
        s[
            "missing_genes"
        ]
    )

    print(
        "MGD beta =",
        s[
            "MGD_beta_ageSexAdj"
        ]
    )

    print(
        "SE =",
        s[
            "SE"
        ]
    )

    print(
        "t =",
        s[
            "t"
        ]
    )

    print(
        "p =",
        s[
            "pvalue"
        ]
    )

    print(
        "95% CI =",
        (
            s[
                "CI95_low"
            ],
            s[
                "CI95_high"
            ],
        )
    )

    print()

    print(
        "IMPORTANT:"
    )

    print(
        "GSEA is intentionally not included "
        "in this script yet because the "
        "original GSEApy log did not record "
        "the exact permutation number and "
        "random seed. These parameters will "
        "be added only after they are verified."
    )

    print()

    print(
        "DONE"
    )

    print(
        "Outputs:",
        OUTPUT_DIR
    )


if __name__ == "__main__":
    main()