# GSE291177 external-transfer analysis specification

## Scientific role
Independent cross-species external transfer of an inflammatory-defense **program**, not replication of a universal fixed-gene signature.

## What was fixed before GSE291177 was examined
The 10-gene Awat2 inflammatory-defense panel already existed in the MGD manuscript:
Cxcl5, Cxcl1, Cxcl2, Ccl20, Slpi, S100a8, Ifitm1, Ifitm3, Il1rn, Lcn2.

## What is fixed from this analysis version forward
- Primary contrast: UVB WT vs untreated WT.
- Input: deposited GSE291177 TPM matrix.
- Transformation: log2(TPM+1).
- Missing genes: reported as not represented; no substitution.
- Module score: gene-wise z scores across the 7 WT libraries, then unweighted mean.
- Inferential unit: deposited RNA-seq library.
- Primary test: exact two-sided permutation of the sample-level module score, all 35 allocations.
- Sensitivities: all-17 z reference, percentile-rank score, leave-one-gene-out.
- P1.hMR analysis: secondary only.
- PCA: QC only.

## Important boundary
This GSE291177 analysis was added after the dataset became available and is not prospectively registered. Only the gene panel predates examination of GSE291177. Do not describe the whole analysis as preregistered or prospective.

## Primary result
- Genes measurable: 9/10 (Cxcl1, Cxcl2, Ccl20, Slpi, S100a8, Ifitm1, Ifitm3, Il1rn, Lcn2)
- Not represented: Cxcl5
- Gene directions: 9/9 higher in UVB WT
- Module difference: 1.257703
- Exact two-sided permutation P: 0.028571
- All three UVB-WT library scores exceed all four control-WT library scores: True

## Secondary P1.hMR result
- Module difference: 1.364774
- Exact two-sided permutation P: 0.007937
- UVB×genotype interaction P on a common all-17 score scale: 0.415871

## Experimental-unit caveat
The Nature Communications article reports 6 WT eyes for the UVB RNA-seq experiment, while the GEO deposit contains 3 UVB-WT RNA-seq libraries. The GEO metadata do not explicitly state the eye-to-library pooling scheme. Accordingly, this package treats each deposited RNA-seq library as the inferential unit and does not count eyes as independent replicates.
