#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Reproduce Awat2 MG state composition and extended defense-panel pseudobulk."""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "source_data" / "GSE261036"
OUT_DIR = PROJECT_ROOT / "results" / "awat2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RAW_FILE = DATA_DIR / "GSE261036_MG_highconfidence_rawcounts.h5ad"
LABEL_FILE = DATA_DIR / "GSE261036_MG_Leiden_test.h5ad"

SAMPLES = ["WT1", "WT2", "WT3", "KO"]

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

EXTENDED10 = [
    "Cxcl5", "Cxcl1", "Cxcl2", "Ccl20", "Slpi",
    "S100a8", "Ifitm1", "Ifitm3", "Il1rn", "Lcn2",
]

print("Reading frozen Awat2 objects...")
raw = sc.read_h5ad(RAW_FILE)
lab = sc.read_h5ad(LABEL_FILE)

if set(raw.obs_names) != set(lab.obs_names):
    raise ValueError("Raw-count and Leiden objects do not contain the same cells.")

lab = lab[raw.obs_names, :].copy()
raw.obs["MG_leiden_r02"] = lab.obs["MG_leiden_r02"].astype(str).values
raw.obs["MG_leiden_r04"] = lab.obs["MG_leiden_r04"].astype(str).values

if "gene_symbol" not in raw.var.columns:
    raise ValueError("Expected raw.var['gene_symbol'] was not found.")

# 1) Composition sensitivity at r=0.2 and r=0.4.
rows = []
for resolution, column, state_map in [
    ("r0.2", "MG_leiden_r02", STATE_MAP_R02),
    ("r0.4", "MG_leiden_r04", STATE_MAP_R04),
]:
    d = raw.obs[["sample", column]].copy()
    d["state"] = d[column].astype(str).map(state_map)
    counts = pd.crosstab(d["state"], d["sample"])
    pct = counts.div(counts.sum(axis=0), axis=1) * 100

    for state in state_map.values():
        for sample in SAMPLES:
            rows.append({
                "resolution": resolution,
                "state": state,
                "sample": sample,
                "n_cells": int(counts.loc[state, sample]),
                "percent_of_sample": float(pct.loc[state, sample]),
            })

composition = pd.DataFrame(rows)
composition.to_csv(
    OUT_DIR / "Awat2_MG_state_composition_r02_r04.csv",
    index=False,
)

# 2) r=0.4 state x sample raw-count pseudobulk CPM for 10 representative genes.
gene_symbols = raw.var["gene_symbol"].astype(str).to_numpy()
symbol_to_index = {}
for i, symbol in enumerate(gene_symbols):
    if symbol not in symbol_to_index:
        symbol_to_index[symbol] = i

missing = [g for g in EXTENDED10 if g not in symbol_to_index]
if missing:
    raise ValueError(f"Missing target genes: {missing}")

rows = []
for cluster, state in STATE_MAP_R04.items():
    for sample in SAMPLES:
        mask = (
            (raw.obs["MG_leiden_r04"].astype(str) == cluster)
            & (raw.obs["sample"].astype(str) == sample)
        ).to_numpy()

        n_cells = int(mask.sum())
        library_counts = float(raw.X[mask, :].sum())

        for gene in EXTENDED10:
            idx = symbol_to_index[gene]
            gene_counts = float(raw.X[mask, idx].sum())
            cpm = gene_counts / library_counts * 1e6

            rows.append({
                "cluster": cluster,
                "state": state,
                "sample": sample,
                "gene": gene,
                "n_cells": n_cells,
                "library_counts": library_counts,
                "gene_counts": gene_counts,
                "CPM": cpm,
            })

state_cpm = pd.DataFrame(rows)
state_cpm.to_csv(
    OUT_DIR / "Awat2_MG_extended10_state_pseudobulk_CPM.csv",
    index=False,
)

wide = (
    state_cpm
    .pivot_table(index=["cluster", "state", "gene"], columns="sample", values="CPM")
    .reset_index()
)
wide["KO_above_allWT"] = (
    wide["KO"] > wide[["WT1", "WT2", "WT3"]].max(axis=1)
)
wide.to_csv(
    OUT_DIR / "Awat2_MG_extended10_KO_above_allWT.csv",
    index=False,
)

print()
print("DONE")
print("Composition file:",
      OUT_DIR / "Awat2_MG_state_composition_r02_r04.csv")
print("Extended-panel result:",
      int(wide["KO_above_allWT"].sum()), "/", len(wide),
      "state x gene comparisons with KO > all three WT libraries.")
