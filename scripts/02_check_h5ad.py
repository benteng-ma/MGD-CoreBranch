#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
import sys
import anndata as ad
import numpy as np

if len(sys.argv) != 2:
    print("Usage: python 02_check_h5ad.py data/GSE274498_MG_all_rawcounts.h5ad")
    sys.exit(1)

p = Path(sys.argv[1])
adata = ad.read_h5ad(p)
print("FILE:", p)
print("SHAPE:", adata.shape)
print("X type:", type(adata.X))
print("obs columns:", list(adata.obs.columns))
print("MG states:")
print(adata.obs["MG_state"].value_counts())
print("Samples:")
print(adata.obs["sample"].value_counts())
print("Ages:")
print(adata.obs["age"].value_counts())
mins = adata.X.data.min() if adata.X.nnz else 0
print("Minimum non-zero count:", mins)
print("Raw-integer check:", bool(np.allclose(adata.X.data, np.round(adata.X.data))))
print("Zero-count cells:", int((np.asarray(adata.X.sum(axis=1)).ravel() == 0).sum()))
print("Zero-count genes:", int((np.asarray(adata.X.sum(axis=0)).ravel() == 0).sum()))
