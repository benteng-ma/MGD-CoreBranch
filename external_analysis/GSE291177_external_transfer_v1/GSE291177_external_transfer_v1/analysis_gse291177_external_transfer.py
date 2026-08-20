#!/usr/bin/env python3
"""
Reproduce the GSE291177 external-transfer analysis.

Expected inputs in the working directory or supplied through --tpm and --xml:
  GSE291177_tpm_counts.txt.gz
  GSE291177_family.xml.tgz

The primary gene panel predates this external dataset:
Cxcl5, Cxcl1, Cxcl2, Ccl20, Slpi, S100a8, Ifitm1, Ifitm3, Il1rn, Lcn2
"""
import argparse, itertools, re, tarfile, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import matplotlib.pyplot as plt

PANEL = ["Cxcl5","Cxcl1","Cxcl2","Ccl20","Slpi","S100a8","Ifitm1","Ifitm3","Il1rn","Lcn2"]

def exact_perm_diff(scores, n_treat, treat_indices):
    vals=np.asarray(scores,float)
    ti=set(treat_indices)
    obs=vals[list(ti)].mean()-vals[[i for i in range(len(vals)) if i not in ti]].mean()
    ds=[]
    for comb in itertools.combinations(range(len(vals)),n_treat):
        s=set(comb)
        ds.append(vals[list(s)].mean()-vals[[i for i in range(len(vals)) if i not in s]].mean())
    ds=np.asarray(ds)
    return obs, np.mean(ds>=obs-1e-12), np.mean(np.abs(ds)>=abs(obs)-1e-12)

def module(expr, genes, cols, refcols):
    present=[g for g in genes if g in expr.index]
    x=np.log2(expr.loc[present,cols]+1)
    r=np.log2(expr.loc[present,refcols]+1)
    z=x.sub(r.mean(axis=1),axis=0).div(r.std(axis=1,ddof=1).replace(0,np.nan),axis=0)
    return present,x,z,z.mean(axis=0)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--tpm",default="GSE291177_tpm_counts.txt.gz")
    ap.add_argument("--outdir",default="GSE291177_external_transfer_results")
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(exist_ok=True)
    d=pd.read_csv(args.tpm,sep="\t",compression="gzip",decimal=",")
    cols=list(d.columns[5:])
    expr=d.set_index("Gene")[cols].astype(float)
    C=[c for c in cols if c.startswith("C_WT_")]
    U=[c for c in cols if c.startswith("U_WT_")]
    primary=C+U
    present,logx,z,score=module(expr,PANEL,primary,primary)
    obs,p1,p2=exact_perm_diff(score.values,len(U),range(len(C),len(primary)))
    pd.DataFrame({"sample":primary,"group":["Control_WT"]*len(C)+["UVB_WT"]*len(U),"module_score":score.values}).to_csv(out/"primary_module_scores.csv",index=False)
    print("Present genes:",present)
    print("Missing genes:",[g for g in PANEL if g not in expr.index])
    print("UVB-Control module difference:",obs)
    print("Exact two-sided permutation P:",p2)

if __name__=="__main__":
    main()
