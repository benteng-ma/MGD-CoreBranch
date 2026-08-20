# GSE291177 external-transfer package

This folder contains a reproducible external-transfer analysis for the MGD Core–Branch project.

Recommended manuscript interpretation:
**independent cross-species transfer of a predefined inflammatory-defense program**.

Do not describe this as:
- replication of a universal six-gene signature;
- causal validation of MGD;
- proof that MR amplifies the module;
- a DESeq2 reanalysis of GSE291177.

The deposited input used here is TPM-normalized; therefore the primary inference is a sample-level module-score permutation test.

Key files:
- `FROZEN_ANALYSIS_SPEC.md` — frozen analysis boundaries
- `methods_ready_text.md` — draft Methods paragraph
- `results_ready_text.md` — draft Results paragraph
- `sample_metadata.csv` — GEO GSM/library mapping reconstructed from MINiML XML
- `primary_module_scores.csv`
- `primary_panel_gene_effects.csv`
- `sensitivity_summary.csv`
- `leave_one_gene_out_sensitivity.csv`
- `secondary_P1hMR_module_scores.csv`
- `secondary_two_factor_module_model.csv`
- `pca_qc_coordinates.csv`
- three standalone 300-dpi figures

Important unit rule:
- Use deposited RNA-seq libraries as the inferential units (WT primary comparison: 4 control libraries vs 3 UVB libraries).
- Do not reinterpret the 6 UVB WT eyes reported in the source article as n=6 independent expression replicates.
