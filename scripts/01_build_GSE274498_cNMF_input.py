#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Build official cNMF input files for GSE274498 from the four raw 10x matrices
and the refined MG-state metadata.

Output:
  1) GSE274498_MG_all_rawcounts.h5ad
  2) GSE274498_acinar_meibocyte_rawcounts.h5ad
  3) input_build_summary.tsv

The AnnData .X matrix is raw integer counts (cells x genes).
Duplicate gene symbols are collapsed by summing counts across duplicated features.
Only cells present in the supplied refined-MG metadata are retained.
"""

from pathlib import Path
import argparse
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import anndata as ad

SAMPLES = [
    ("8W_R1", "8W",  "GSM8450395_8_Week_R1"),
    ("8W_R2", "8W",  "GSM8450396_8_Week_R2"),
    ("21M_R1","21M", "GSM8450397_21_Month_R1"),
    ("21M_R2","21M", "GSM8450398_21_Month_R2"),
]

ACINAR_MEIB_STATES = [
    "Acinar basal",
    "Differentiating meibocyte",
    "Differentiated meibocyte",
]

def find_one(root: Path, filename: str) -> Path:
    hits = list(root.rglob(filename))
    if len(hits) == 0:
        raise FileNotFoundError(
            f"Cannot find {filename} under {root}\n"
            "Please check --raw-dir and confirm the 12 GSE274498 files are present."
        )
    if len(hits) > 1:
        raise RuntimeError(f"Found more than one copy of {filename}: {hits}")
    return hits[0]

def read_10x_sample(raw_dir: Path, prefix: str):
    mtx_p = find_one(raw_dir, prefix + "_matrix.mtx.gz")
    feat_p = find_one(raw_dir, prefix + "_features.tsv.gz")
    bc_p = find_one(raw_dir, prefix + "_barcodes.tsv.gz")

    # Matrix Market is genes x cells
    X = mmread(mtx_p).tocsr()
    feat = pd.read_csv(feat_p, sep="\t", header=None, compression="gzip")
    bcs = pd.read_csv(bc_p, sep="\t", header=None, compression="gzip")[0].astype(str).to_numpy()

    if X.shape[0] != feat.shape[0] or X.shape[1] != len(bcs):
        raise ValueError(
            f"Dimension mismatch for {prefix}: matrix={X.shape}, "
            f"features={feat.shape[0]}, barcodes={len(bcs)}"
        )

    feat = feat.iloc[:, :3].copy()
    feat.columns = ["gene_id", "gene_symbol", "feature_type"]
    feat["gene_id"] = feat["gene_id"].astype(str)
    feat["gene_symbol"] = feat["gene_symbol"].astype(str)
    return X, feat, bcs

def collapse_duplicate_symbols(X_genes_by_cells, feat):
    """
    Collapse duplicated gene symbols by summing their raw counts.
    Returns unique_symbol x cells sparse matrix and a var DataFrame.
    """
    symbols = feat["gene_symbol"].astype(str).to_numpy()

    # Keep stable first-occurrence ordering.
    symbol_to_new = {}
    unique_symbols = []
    group_idx = np.empty(len(symbols), dtype=np.int64)
    for i, s in enumerate(symbols):
        if s not in symbol_to_new:
            symbol_to_new[s] = len(unique_symbols)
            unique_symbols.append(s)
        group_idx[i] = symbol_to_new[s]

    # Aggregation matrix: unique genes x original features.
    A = sp.csr_matrix(
        (np.ones(len(symbols), dtype=np.float32),
         (group_idx, np.arange(len(symbols), dtype=np.int64))),
        shape=(len(unique_symbols), len(symbols)),
    )
    X2 = (A @ X_genes_by_cells).tocsr()

    # Store all source Ensembl IDs for auditability.
    ids_by_symbol = feat.groupby("gene_symbol", sort=False)["gene_id"].apply(lambda x: "|".join(x.astype(str)))
    var = pd.DataFrame(index=pd.Index(unique_symbols, name="gene_symbol"))
    var["source_gene_ids"] = ids_by_symbol.reindex(unique_symbols).to_numpy()
    return X2, var

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True,
                    help="Folder containing (or containing subfolders with) the 12 raw GSE274498 10x files.")
    ap.add_argument("--metadata", required=True,
                    help="Path to 05_cell_metadata_refined_MG_states.csv")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for cNMF input files.")
    args = ap.parse_args()

    raw_dir = Path(args.raw_dir)
    meta_path = Path(args.metadata)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(meta_path)
    required = {"sample","age","barcode","MG_state"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"Metadata is missing required columns: {sorted(missing)}")

    all_adatas = []
    summary = []

    reference_symbols = None

    for sample, age, prefix in SAMPLES:
        print(f"[1/4] Reading {sample}: {prefix}")
        X, feat, bcs = read_10x_sample(raw_dir, prefix)

        # Collapse duplicate gene symbols consistently within each library.
        Xc, var = collapse_duplicate_symbols(X, feat)
        symbols = var.index.astype(str).to_numpy()

        if reference_symbols is None:
            reference_symbols = symbols
        else:
            if not np.array_equal(reference_symbols, symbols):
                raise ValueError(f"Gene symbol order differs in {sample}; stop for manual inspection.")

        bc_to_col = {bc: i for i, bc in enumerate(bcs)}
        md = meta.loc[meta["sample"].astype(str) == sample].copy()

        missing_bcs = [bc for bc in md["barcode"].astype(str) if bc not in bc_to_col]
        if missing_bcs:
            raise ValueError(
                f"{sample}: {len(missing_bcs)} metadata barcodes are absent from raw matrix. "
                f"Example: {missing_bcs[:5]}"
            )

        cols = np.array([bc_to_col[bc] for bc in md["barcode"].astype(str)], dtype=np.int64)
        # Convert to cells x genes.
        Xsub = Xc[:, cols].T.tocsr()

        # Unique cell IDs across libraries.
        unique_cell_ids = [f"{sample}__{bc}" for bc in md["barcode"].astype(str)]
        obs = md.copy()
        obs.index = pd.Index(unique_cell_ids, name="cell_id")
        obs["sample"] = sample
        obs["age"] = age
        obs["original_barcode"] = obs["barcode"].astype(str)

        adata = ad.AnnData(X=Xsub, obs=obs, var=var.copy())
        all_adatas.append(adata)
        summary.append({
            "sample": sample,
            "age": age,
            "retained_MG_cells": adata.n_obs,
            "genes_before_zero_filter": adata.n_vars,
        })

    print("[2/4] Concatenating four libraries")
    adata = ad.concat(all_adatas, axis=0, join="inner", merge="same", index_unique=None)

    # Ensure integer raw counts.
    adata.X = adata.X.tocsr()
    if not np.allclose(adata.X.data, np.round(adata.X.data)):
        raise ValueError("The input matrix does not look like integer raw counts.")

    # Remove genes with zero total counts across selected MG cells.
    gene_total = np.asarray(adata.X.sum(axis=0)).ravel()
    keep_gene = gene_total > 0
    adata = adata[:, keep_gene].copy()

    # Remove any zero-count cells (should not occur).
    cell_total = np.asarray(adata.X.sum(axis=1)).ravel()
    if np.any(cell_total == 0):
        adata = adata[cell_total > 0, :].copy()

    # Add simple audit fields.
    adata.uns["dataset"] = "GSE274498"
    adata.uns["input_type"] = "raw integer counts"
    adata.uns["cell_selection"] = "refined high-confidence MG epithelial states from Stage6B reconstruction"
    adata.uns["duplicate_gene_symbols"] = "collapsed by summing raw counts"
    adata.uns["statistical_unit_note"] = (
        "Each library pools four mice; nuclei/cells are not independent biological replicates."
    )

    all_out = out_dir / "GSE274498_MG_all_rawcounts.h5ad"
    print(f"[3/4] Writing all-MG input: {all_out}")
    adata.write_h5ad(all_out, compression="gzip")

    ac = adata[adata.obs["MG_state"].isin(ACINAR_MEIB_STATES), :].copy()
    gene_total_ac = np.asarray(ac.X.sum(axis=0)).ravel()
    ac = ac[:, gene_total_ac > 0].copy()

    ac_out = out_dir / "GSE274498_acinar_meibocyte_rawcounts.h5ad"
    print(f"[4/4] Writing acinar-meibocyte sensitivity input: {ac_out}")
    ac.write_h5ad(ac_out, compression="gzip")

    summary_df = pd.DataFrame(summary)
    summary_df.loc[len(summary_df)] = {
        "sample":"ALL_MG",
        "age":"8W+21M",
        "retained_MG_cells":adata.n_obs,
        "genes_before_zero_filter":adata.n_vars,
    }
    summary_df.loc[len(summary_df)] = {
        "sample":"ACINAR_MEIB",
        "age":"8W+21M",
        "retained_MG_cells":ac.n_obs,
        "genes_before_zero_filter":ac.n_vars,
    }
    summary_df.to_csv(out_dir / "input_build_summary.tsv", sep="\t", index=False)

    print("\nSUCCESS")
    print(f"All MG:           {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"Acinar/meibocyte: {ac.n_obs} cells x {ac.n_vars} genes")
    print("Files written to:", out_dir)
    print("\nExpected refined all-MG cell count is approximately 10,307 based on the supplied metadata.")

if __name__ == "__main__":
    main()
