#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
11_Awat2_Leiden_resolution_sensitivity.py

Reproduce Awat2 MG Leiden-resolution sensitivity using the frozen
GSE261036_MG_Leiden_test.h5ad object.

Evaluates:
    Leiden resolution 0.2
    Leiden resolution 0.4
    Leiden resolution 0.6

Outputs:
1. cluster counts and ARI
2. r0.6 vs r0.4 cross-tab
3. major-state composition by sample and resolution
4. KO direction relative to the WT range

No clustering is rerun; the frozen Leiden labels are used directly.
"""

from pathlib import Path

import anndata as ad
import pandas as pd

from sklearn.metrics import adjusted_rand_score


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "source_data"
    / "GSE261036"
    / "GSE261036_MG_Leiden_test.h5ad"
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
# Fixed sample order
# ============================================================

SAMPLES = [
    "WT1",
    "WT2",
    "WT3",
    "KO",
]


# ============================================================
# Major-state mappings
# ============================================================

STATE_MAP_R02 = {
    "0": "Defense/stress duct-like",
    "1": "Early/differentiating meibocyte",
    "2": "Differentiated duct",
    "3": "Lipogenic/differentiated meibocyte",
    "4": "Cycling MG",
    "5": "Basal duct",
}


STATE_MAP_R04 = {
    "0": "Defense/stress duct-like",
    "1": "Early/differentiating meibocyte",
    "2": "Differentiated duct",
    "3": "Cycling MG",
    "4": "Lipogenic/differentiated meibocyte",
    "5": "Basal duct",
}


# r0.6 clusters are collapsed back to the r0.4 major-state framework.
# The mapping is based on the dominant overlap in the frozen
# r0.6 × r0.4 cross-tab.

STATE_MAP_R06 = {
    "0": "Defense/stress duct-like",
    "1": "Defense/stress duct-like",
    "2": "Differentiated duct",
    "3": "Cycling MG",
    "4": "Basal duct",
    "5": "Cycling MG",
    "6": "Early/differentiating meibocyte",
    "7": "Lipogenic/differentiated meibocyte",
    "8": "Lipogenic/differentiated meibocyte",
}


STATE_ORDER = [
    "Basal duct",
    "Cycling MG",
    "Defense/stress duct-like",
    "Differentiated duct",
    "Early/differentiating meibocyte",
    "Lipogenic/differentiated meibocyte",
]


RESOLUTIONS = [
    (
        "r0.2",
        "MG_leiden_r02",
        STATE_MAP_R02,
    ),
    (
        "r0.4",
        "MG_leiden_r04",
        STATE_MAP_R04,
    ),
    (
        "r0.6",
        "MG_leiden_r06",
        STATE_MAP_R06,
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
        "AWAT2 LEIDEN RESOLUTION SENSITIVITY"
    )

    print(
        "============================================"
    )

    print()

    if not INPUT_FILE.exists():

        raise FileNotFoundError(
            INPUT_FILE
        )

    # --------------------------------------------------------
    # Read frozen object
    # --------------------------------------------------------

    a = ad.read_h5ad(
        INPUT_FILE
    )

    print(
        "Frozen MG object:",
        a.shape
    )

    print()

    required_columns = [
        "sample",
        "MG_leiden_r02",
        "MG_leiden_r04",
        "MG_leiden_r06",
    ]

    missing = [
        c
        for c in required_columns
        if c not in a.obs.columns
    ]

    if missing:

        raise RuntimeError(
            "Missing obs columns: "
            + ", ".join(
                missing
            )
        )

    obs = a.obs[
        required_columns
    ].copy()

    for c in [
        "sample",
        "MG_leiden_r02",
        "MG_leiden_r04",
        "MG_leiden_r06",
    ]:

        obs[c] = (
            obs[c]
            .astype(str)
        )

    # ========================================================
    # 1. Cluster counts and ARI
    # ========================================================

    r02 = obs[
        "MG_leiden_r02"
    ]

    r04 = obs[
        "MG_leiden_r04"
    ]

    r06 = obs[
        "MG_leiden_r06"
    ]

    ari_02_04 = adjusted_rand_score(
        r02,
        r04
    )

    ari_04_06 = adjusted_rand_score(
        r04,
        r06
    )

    ari_02_06 = adjusted_rand_score(
        r02,
        r06
    )

    cluster_summary = pd.DataFrame(
        [
            {
                "Comparison":
                    "r0.2_vs_r0.4",

                "Clusters_A":
                    r02.nunique(),

                "Clusters_B":
                    r04.nunique(),

                "ARI":
                    ari_02_04,
            },
            {
                "Comparison":
                    "r0.4_vs_r0.6",

                "Clusters_A":
                    r04.nunique(),

                "Clusters_B":
                    r06.nunique(),

                "ARI":
                    ari_04_06,
            },
            {
                "Comparison":
                    "r0.2_vs_r0.6",

                "Clusters_A":
                    r02.nunique(),

                "Clusters_B":
                    r06.nunique(),

                "ARI":
                    ari_02_06,
            },
        ]
    )

    cluster_summary_file = (
        OUT_DIR
        / "Awat2_Leiden_resolution_cluster_counts_ARI.csv"
    )

    cluster_summary.to_csv(
        cluster_summary_file,
        index=False,
    )

    print(
        "Cluster counts:"
    )

    print(
        "r0.2 =",
        r02.nunique()
    )

    print(
        "r0.4 =",
        r04.nunique()
    )

    print(
        "r0.6 =",
        r06.nunique()
    )

    print()

    print(
        "ARI r0.2 vs r0.4 =",
        ari_02_04
    )

    print(
        "ARI r0.4 vs r0.6 =",
        ari_04_06
    )

    print(
        "ARI r0.2 vs r0.6 =",
        ari_02_06
    )

    print()

    # ========================================================
    # 2. r0.6 versus r0.4 cross-tab
    # ========================================================

    cross = pd.crosstab(
        r06,
        r04
    )

    cross.index.name = (
        "MG_leiden_r06"
    )

    cross.columns = [
        "r04_cluster_"
        + str(c)
        for c in cross.columns
    ]

    cross_file = (
        OUT_DIR
        / "Awat2_Leiden_r06_vs_r04_crosstab.csv"
    )

    cross.to_csv(
        cross_file
    )

    # ========================================================
    # 3. Composition by resolution / sample / major state
    # ========================================================

    composition_rows = []

    for (
        resolution,
        cluster_column,
        state_map,
    ) in RESOLUTIONS:

        temp = obs[
            [
                "sample",
                cluster_column,
            ]
        ].copy()

        temp[
            "MG_state"
        ] = (
            temp[
                cluster_column
            ]
            .map(
                state_map
            )
        )

        if temp[
            "MG_state"
        ].isna().any():

            unknown = sorted(
                temp.loc[
                    temp[
                        "MG_state"
                    ].isna(),
                    cluster_column
                ].unique()
            )

            raise RuntimeError(
                f"{resolution}: unmapped clusters: "
                + ", ".join(
                    unknown
                )
            )

        for sample in SAMPLES:

            sample_df = temp.loc[
                temp[
                    "sample"
                ].eq(
                    sample
                )
            ]

            total = len(
                sample_df
            )

            if total == 0:

                raise RuntimeError(
                    f"No cells found for sample {sample}"
                )

            counts = (
                sample_df[
                    "MG_state"
                ]
                .value_counts()
            )

            for state in STATE_ORDER:

                n = int(
                    counts.get(
                        state,
                        0
                    )
                )

                pct = (
                    100.0
                    * n
                    / total
                )

                composition_rows.append(
                    {
                        "Resolution":
                            resolution,

                        "Sample":
                            sample,

                        "MG_state":
                            state,

                        "Cells_n":
                            n,

                        "Sample_total_MG_cells":
                            total,

                        "Percent":
                            pct,
                    }
                )

    composition = pd.DataFrame(
        composition_rows
    )

    composition_file = (
        OUT_DIR
        / "Awat2_Leiden_resolution_state_composition.csv"
    )

    composition.to_csv(
        composition_file,
        index=False,
    )

    # ========================================================
    # 4. KO direction relative to WT range
    # ========================================================

    direction_rows = []

    for resolution in [
        "r0.2",
        "r0.4",
        "r0.6",
    ]:

        sub = composition.loc[
            composition[
                "Resolution"
            ].eq(
                resolution
            )
        ]

        for state in STATE_ORDER:

            s = (
                sub.loc[
                    sub[
                        "MG_state"
                    ].eq(
                        state
                    )
                ]
                .set_index(
                    "Sample"
                )[
                    "Percent"
                ]
            )

            wt_values = [
                float(
                    s.loc[
                        "WT1"
                    ]
                ),
                float(
                    s.loc[
                        "WT2"
                    ]
                ),
                float(
                    s.loc[
                        "WT3"
                    ]
                ),
            ]

            ko = float(
                s.loc[
                    "KO"
                ]
            )

            wt_min = min(
                wt_values
            )

            wt_max = max(
                wt_values
            )

            if ko > wt_max:

                direction = (
                    "above_WT_range"
                )

            elif ko < wt_min:

                direction = (
                    "below_WT_range"
                )

            else:

                direction = (
                    "within_WT_range"
                )

            direction_rows.append(
                {
                    "Resolution":
                        resolution,

                    "MG_state":
                        state,

                    "WT1_percent":
                        wt_values[0],

                    "WT2_percent":
                        wt_values[1],

                    "WT3_percent":
                        wt_values[2],

                    "KO_percent":
                        ko,

                    "WT_min_percent":
                        wt_min,

                    "WT_max_percent":
                        wt_max,

                    "KO_vs_WT_range":
                        direction,
                }
            )

    direction = pd.DataFrame(
        direction_rows
    )

    direction_file = (
        OUT_DIR
        / "Awat2_Leiden_resolution_KO_direction_summary.csv"
    )

    direction.to_csv(
        direction_file,
        index=False,
    )

    # ========================================================
    # Print compact validation table
    # ========================================================

    print(
        "============================================"
    )

    print(
        "KO DIRECTION RELATIVE TO WT RANGE"
    )

    print(
        "============================================"
    )

    print()

    compact = direction[
        [
            "Resolution",
            "MG_state",
            "KO_percent",
            "KO_vs_WT_range",
        ]
    ]

    print(
        compact.to_string(
            index=False
        )
    )

    print()

    print(
        "Saved:"
    )

    print(
        cluster_summary_file
    )

    print(
        cross_file
    )

    print(
        composition_file
    )

    print(
        direction_file
    )

    print()

    print(
        "INTERPRETATION NOTE:"
    )

    print(
        "Resolution 0.6 primarily subdivides "
        "existing major MG populations rather "
        "than defining new major MG identities. "
        "Resolution 0.4 is retained as an "
        "intermediate representative granularity."
    )

    print()

    print(
        "DONE"
    )


if __name__ == "__main__":
    main()