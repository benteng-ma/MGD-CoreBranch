# MGD Core–Branch Transcriptomic Integration

Reproducibility code and frozen analysis inputs for the study:

**Cross-etiology transcriptomic integration reveals an inflammatory-defense core and distinct epithelial remodeling branches in meibomian gland dysfunction**

## Current release scope

This repository contains reproducibility materials for the analyses directly reprocessed or recalculated in this study, including:

- GSE274498 aging mouse MG formal cNMF analysis
- GSE261036 Awat2-KO MG-state, pseudobulk, WT-reference, and Leiden-resolution sensitivity analyses
- GSE166784 Hsd3b6-KO to Awat2-KO cross-model directional concordance
- GSE17822 human MGD age/sex-adjusted expression and pathway-level sensitivity analyses

The public GitHub repository contains code, small frozen inputs, and analysis outputs. Large processed AnnData files (`*.h5ad`) are intentionally excluded from GitHub and are packaged in the archival Zenodo reproducibility release.

GSE274497 spatial localization and GSE274496 GLI2 perturbational evidence are used in the manuscript as contextual evidence from the source studies and are not presented as analyses fully reproduced by this repository.

## Repository structure

```text
MGD_CoreBranch_CodeRelease/
├─ environment/
├─ scripts/
├─ source_data/
├─ results/
└─ docs/
## Software environment

The analysis environment is recorded in:

- `environment/environment.yml`
- `environment/requirements.txt`

Core versions include Python 3.10.20, cNMF 1.7.1, Scanpy 1.11.5, AnnData 0.11.4, NumPy 2.2.5, pandas 2.3.3, SciPy 1.15.3, and GSEApy 1.3.1.

The environment files were checked successfully with `pip check`, with no broken requirements detected.

To reconstruct the environment:

```bash
conda env create -f environment/environment.yml
conda activate mgd-corebranch
```

---

## 1. GSE274498 aging mouse MG and formal cNMF analysis

### Frozen input

The primary frozen MG raw-count input is:

`source_data/GSE274498/GSE274498_MG_all_rawcounts.h5ad`

Refined MG-state metadata are stored in:

`source_data/GSE274498/05_cell_metadata_refined_MG_states.csv`

The frozen cNMF input contains 10,307 MG cells/nuclei.

The upstream input-building script is:

`scripts/01_build_GSE274498_cNMF_input.py`

The script reconstructs the MG raw-count object from the public GSE274498 10x matrices and refined MG-state metadata.

### Formal cNMF configuration

The formal cNMF analysis evaluated:

- K = 5–20
- 100 NMF replicates per K
- 3,000 overdispersed genes
- master random seed = 14
- beta loss = frobenius
- initialization = random
- maximum NMF iterations = 1000

This corresponds to 1,600 NMF runs in total.

The reconstructed formal Windows runner is:

`scripts/05_run_GSE274498_cNMF_formal.bat`

The workflow is:

```text
prepare
→ factorize
→ combine
→ k_selection_plot
→ K=8 consensus
```

The final consensus solution used K=8 and a local-density threshold of 0.1.

The runner was reconstructed from the original formal cNMF parameter files and checked against the stored formal outputs. A full 1,600-run factorization was not repeated solely for release validation because of its computational cost.

### K-selection summary

Frozen K-selection statistics are stored in:

`source_data/cnmf/GSE274498_allMG_formal.k_selection_stats.df.npz`

Run:

```bash
python scripts/06_cNMF_K_selection_summary.py
```

Output:

`results/cnmf/formal_K5to20_selection_stats.csv`

K=8 showed the highest consensus stability across K=5–20:

```text
K=8 consensus stability = 0.965913
```

### K=8 program usage

Frozen K=8 consensus usage is stored in:

`source_data/cnmf/GSE274498_allMG_formal.usages.k_8.dt_0_1.consensus.txt`

Run:

```bash
python scripts/07_cNMF_K8_usage_summary.py
```

Generated outputs include:

```text
results/cnmf/formal_K8_mean_usage_by_MG_state.csv
results/cnmf/formal_K8_mean_usage_by_MG_state_sample.csv
results/cnmf/formal_K8_mean_usage_by_MG_state_sample_P2excluded.csv
results/cnmf/formal_K8_age_direction_sensitivity_summary.csv
results/cnmf/formal_K8_P7_P8_library_age_summary.csv
```

For the full normalization, P1–P8 usage is normalized within each cell to sum to one before aggregation.

For the P2-excluded sensitivity analysis, P2 is removed and the remaining seven programs are renormalized within each cell before aggregation.

The aging directions of P7 and P8 remained positive under both normalization strategies.

These analyses are interpreted at the pooled-library level rather than treating individual cells as independent biological replicates.

### K=8 program interpretation

The final program interpretation table is:

`source_data/cnmf/formal_K8_program_final_interpretation.csv`

Final labels are:

```text
P1  Ductal basal developmental/regulatory GEP
P2  Ribosomal/translation-associated GEP
P3  Stress-responsive transcription/signaling GEP
P4  Acinar basal/neural-guidance-RTK GEP
P5  Ductal suprabasal/ductular structural-signaling GEP
P6  Orifice keratinization/cornification GEP
P7  Differentiating-meibocyte sterol/lipid-biosynthesis GEP
P8  Differentiated-meibocyte specialized meibum-lipid GEP
```

P2 is treated as a ribosomal/translation-associated program rather than a primary aging mechanism.

P3 is treated as a stress-responsive transcription/signaling program and is not labeled as the cross-dataset inflammatory-defense core.

### Pathway enrichment

The formal enrichment script is:

`scripts/K8_formal_enrichment.py`

Input gene scores are stored in:

`source_data/cnmf/GSE274498_allMG_formal.gene_spectra_score.k_8.dt_0_1.txt`

The analysis used the top 100 genes from each K=8 program and queried:

```text
GO_Biological_Process_2026
Reactome_Pathways_2024
```

Because online enrichment libraries can change over time, the exact formal outputs used for interpretation are frozen under:

`source_data/cnmf/K8_enrichment_top100/`

The curated pathway table used for final interpretation is:

`source_data/cnmf/formal_K8_curated_pathways_with_stats.csv`

---

## 2. GSE261036 Awat2-KO MG analysis

### Frozen processed inputs

The Awat2 release uses:

```text
source_data/GSE261036/GSE261036_MG_highconfidence_rawcounts.h5ad
source_data/GSE261036/GSE261036_MG_Leiden_test.h5ad
```

The high-confidence MG object contains 8,903 MG cells.

The experiment contains three WT pooled libraries and one pooled KO library.

Therefore, individual cells are not treated as independent biological replicates, and replicate-level KO differential-expression significance is not inferred from the single KO library.

### MG-state composition and extended defense panel

Run:

```bash
python scripts/03_Awat2_state_composition_and_defense.py
```

Generated outputs:

```text
results/awat2/Awat2_MG_state_composition_r02_r04.csv
results/awat2/Awat2_MG_extended10_state_pseudobulk_CPM.csv
results/awat2/Awat2_MG_extended10_KO_above_allWT.csv
```

The r=0.4 MG states are:

```text
Defense/stress duct-like
Early/differentiating meibocyte
Differentiated duct
Cycling MG
Lipogenic/differentiated meibocyte
Basal duct
```

The extended 10-gene defense/stress panel is:

```text
Cxcl5
Cxcl1
Cxcl2
Ccl20
Slpi
S100a8
Ifitm1
Ifitm3
Il1rn
Lcn2
```

Across six MG states and ten genes, 45 of 60 state × gene comparisons showed Awat2-KO expression above all three WT libraries.

This analysis is used as descriptive robustness support rather than as an independent formal enrichment test.

---

## 3. GSE166784 Hsd3b6-KO to Awat2-KO cross-model concordance

### Frozen Hsd3b6 source

The author-provided supplementary differential-expression table is stored as:

`source_data/GSE166784/Supplementary_Data_1.xls`

Hsd3b6 FDR-upregulated genes are defined as:

```text
log2FoldChange > 0
adjusted P-value < 0.05
```

### Cross-model directional concordance

Run:

```bash
python scripts/04_Hsd3b6_Awat2_crossmodel_concordance.py
```

The script:

1. extracts Hsd3b6 FDR-upregulated genes from the author-supplied table;
2. constructs whole-MG Awat2 sample-level pseudobulk CPM values;
3. defines the expressed Awat2 background as genes reaching at least 1 CPM in at least one sample;
4. defines reference-robust Awat2 directional elevation as:

```text
KO > WT1
AND
KO > WT2
AND
KO > WT3
```

5. performs a one-sided Fisher exact test.

Validated release result:

```text
Hsd3b6 FDR-upregulated genes total: 18
Detected in Awat2 expressed background: 9
KO > all three WT libraries: 6 / 9

Awat2 expressed background: 12,926 genes
Background KO > all WT: 3,276 / 12,926

Fisher odds ratio: 5.9003
One-sided Fisher P: 0.010704
95% conditional odds-ratio CI: 1.259–36.476
```

The six directionally concordant genes are:

```text
Ceacam10
H2-Aa
H2-Ab1
H2-Eb1
Wdfy1
Zbp1
```

The three detected Hsd3b6 FDR-upregulated genes that did not exceed all three WT libraries are:

```text
Ifi47
Irgm2
Rtp4
```

This result is interpreted as reference-robust directional cross-model concordance, not as replicated Awat2 differential-expression significance.

Generated outputs:

```text
results/cross_model/Hsd3b6_FDR_UP_from_author_supplement.csv
results/cross_model/Awat2_wholeMG_allGene_pseudobulk_CPM.csv
results/cross_model/Hsd3b6_FDR_UP_detected_in_Awat2.csv
results/cross_model/Hsd3b6_Awat2_directional_concordance_summary.csv
```

---

## Recommended execution order

For direct reproduction from the frozen processed inputs:

```text
1. Create the software environment
2. Check frozen AnnData files with 02_check_h5ad.py
3. Run 06_cNMF_K_selection_summary.py
4. Run 07_cNMF_K8_usage_summary.py
5. Run 03_Awat2_state_composition_and_defense.py
6. Run 04_Hsd3b6_Awat2_crossmodel_concordance.py
```

The full formal cNMF factorization can be rerun with:

`scripts/05_run_GSE274498_cNMF_formal.bat`

The formal enrichment can be rerun with:

`scripts/K8_formal_enrichment.py`

However, frozen original enrichment outputs are retained because online enrichment databases may change over time.

---

## Statistical-unit considerations

The analyses preserve the experimental sampling structure.

For GSE274498, the biological unit is the pooled sequencing library rather than the individual nucleus.

For GSE261036, the Awat2-KO group contains a single pooled KO library. Therefore:

- cells are not treated as independent biological replicates;
- Awat2 KO-vs-WT findings are interpreted primarily as descriptive composition, pseudobulk direction, sensitivity analysis, or cross-model concordance;
- replicate-level KO differential-expression significance is not claimed.

The Hsd3b6-to-Awat2 comparison uses an independently defined Hsd3b6 FDR-upregulated gene set and evaluates whether those genes disproportionately show Awat2-KO expression greater than all three WT libraries.

---

## Data availability

The original datasets directly reanalyzed in this study are publicly available from NCBI GEO:

- GSE274498
- GSE261036
- GSE166784
- GSE17822

GSE274497 and GSE274496 are referenced in the manuscript for contextual spatial-localization and GLI2 perturbational evidence from the source studies.

The public GitHub repository contains analysis code, small frozen inputs, and derived outputs:

https://github.com/benteng-ma/MGD-CoreBranch

The following large processed AnnData files are intentionally excluded from GitHub via `.gitignore` and are packaged in the archival Zenodo reproducibility release:

- `source_data/GSE274498/GSE274498_MG_all_rawcounts.h5ad`
- `source_data/GSE261036/GSE261036_MG_highconfidence_rawcounts.h5ad`
- `source_data/GSE261036/GSE261036_MG_Leiden_test.h5ad`

Large public raw sequencing files are not duplicated when they are already available from GEO.

Zenodo DOI for the archival reproducibility release: https://doi.org/10.5281/zenodo.21914654.

## Reproducibility status

Current script status:

```text
01_build_GSE274498_cNMF_input.py
    release input-building script

02_check_h5ad.py
    generic AnnData integrity checker

03_Awat2_state_composition_and_defense.py
    execution-validated

04_Hsd3b6_Awat2_crossmodel_concordance.py
    execution-validated

05_run_GSE274498_cNMF_formal.bat
    parameter-verified reconstructed formal runner

06_cNMF_K_selection_summary.py
    execution-validated

07_cNMF_K8_usage_summary.py
    execution-validated

K8_formal_enrichment.py
    formal enrichment script with frozen original outputs
```

The formal cNMF runner was checked against the original parameter files, including the K range, replicate allocation, random-seed sequence, NMF configuration, 3,000-gene feature set, and consensus settings.

---

## Numerical reproducibility note

Minor floating-point differences may occur across operating systems, package builds, and sparse-matrix implementations.

During release validation, one near-tie involving `Fra10ac1` differed only at a very small CPM boundary between earlier derived output and direct raw-count reconstruction.

The release statistics are defined by direct recalculation from the frozen raw-count AnnData object.

This numerical edge case does not affect the Hsd3b6 6/9 directional-concordance result.

---

## Current scope and planned completion

This README documents the analyses currently implemented in the release repository.



The public GitHub repository is available at https://github.com/benteng-ma/MGD-CoreBranch, and the archival reproducibility release uses the Zenodo DOI https://doi.org/10.5281/zenodo.21914654.

---

## Licensing

Original software code authored for this release is licensed under the MIT License; see `LICENSE_CODE.txt`. Original documentation and other original author-generated material are licensed under CC BY 4.0 unless otherwise noted and only to the extent that the authors hold the relevant rights. Files under `source_data/`, including processed derivatives of public GEO datasets, are not relicensed by this release; reuse remains subject to the terms and citation requirements of the original data sources. See `LICENSES.md` for the complete scope.

## Citation

If this code or the associated data are used, please cite the accompanying manuscript and the relevant GEO datasets.

The public repository is https://github.com/benteng-ma/MGD-CoreBranch. Zenodo DOI for the archival reproducibility release: https://doi.org/10.5281/zenodo.21914654. The final manuscript citation will be added after publication.

---

## Contact

Correspondence regarding the study and analysis workflow should be directed to the corresponding author listed in the associated manuscript.



---

## 4. GSE17822 human MGD analysis

### Dataset and frozen inputs

The human validation analysis uses GEO accession:

`GSE17822`

Frozen inputs are stored under:

```text
source_data/GSE17822/GSE17822_series_matrix.txt.gz
source_data/GSE17822/GSE17822_family.soft.gz
source_data/GSE17822/gene_sets.gmt
```

The study contains 12 human eyelid tarsal-plate samples:

```text
6 Normal
6 MGD
```

The sample covariates were reconstructed directly from GEO metadata and include disease state, age, and sex.

The GEO series-matrix expression values were generated from BeadStudio cubic-spline-normalized data. According to the GEO processing record, negative values had been set to zero and 16 had been added to all intensities before deposition. The release analysis therefore applies log2 transformation to the deposited positive intensities without additional between-array normalization.

### Age- and sex-adjusted probe-level model

The formal human model is:

```text
log2(expression) ~ MGD + age + Male
```

where a positive MGD coefficient indicates higher expression in MGD relative to Normal.

The analysis contains:

```text
48,803 Illumina probes
12 samples
```

The release script is:

`scripts/08_GSE17822_human_analysis.py`

Generated probe-level outputs include:

```text
results/human/GSE17822_probe_DE_ageSex_adjusted.csv
results/human/GSE17822_probe_DE_ageSex_adjusted_annotated.csv
results/human/GSE17822_sample_metadata.csv
```

For the 42,709 probes with non-zero variance across the 12 samples, the reconstructed age/sex-adjusted OLS statistics reproduce the original analysis to numerical precision.

There are 6,094 probes with exactly zero variance across all 12 samples. For such probes, coefficient and standard-error estimates are both at machine-precision scale and resulting t statistics are numerically unstable and have no biological interpretation. These probes are explicitly addressed by the downstream zero-variance-filtered GSEA sensitivity analysis.

No human probe-level result is used as evidence of genome-wide FDR-significant differential expression.

### Gene-level representative-probe analysis

Probe annotations are reconstructed directly from the GPL6947 platform table contained in the GEO family SOFT file.

One representative probe is selected per gene symbol using the probe with the highest mean log2 expression across all 12 samples.

This produces:

```text
25,159 gene-level representative probes
```

A probe-selection sensitivity analysis instead chooses the probe with the highest across-sample log2-expression variance for each gene.

The two rules select the identical probe for:

```text
22,822 / 25,159 genes
90.71%
```

Outputs are:

```text
results/human/GSE17822_gene_DE_ageSex_adjusted_representativeProbe.csv
results/human/GSE17822_gene_DE_ageSex_adjusted_maxLog2VarianceProbe.csv
results/human/GSE17822_meanProbe_vs_maxLog2VarianceProbe.csv
```

### Fixed mouse-derived defense module

The prespecified mouse-derived defense candidate set contained 17 genes:

```text
CXCL5
LCN2
IL1RN
IFITM1
IFITM3
CXCL1
CXCL2
CCL20
SLPI
S100A8
IFI47
ZBP1
ISG15
RTP4
CXCL10
PLAC8
WDFY1
```

Sixteen were measurable in the human representative-probe table.

The missing gene was:

```text
IFI47
```

For each measurable gene, expression was standardized across the 12 human samples using the sample standard deviation (`ddof=1`). The 16 standardized values were then averaged within each sample to obtain the DefenseScore.

The age- and sex-adjusted model was:

```text
DefenseScore ~ MGD + age + Male
```

Validated result:

```text
MGD beta = 0.225521
SE = 0.404505
t = 0.557524
P = 0.592423
95% CI = -0.707269 to 1.158312
```

Thus, the fixed mouse-derived 16-gene module was not significantly elevated in human MGD.

Outputs are:

```text
results/human/GSE17822_mouseDefense_moduleScore.csv
results/human/GSE17822_mouseDefense_moduleScore_summary.csv
results/human/GSE17822_mouseDefense_selected_probes.csv
```

### Human pathway-level GSEA

Human pathway analysis was intentionally restricted to two predefined defense-related GO gene sets stored in:

`source_data/GSE17822/gene_sets.gmt`

The tested pathways are:

```text
Response to Type II Interferon (GO:0034341)
Response to Cytokine (GO:0034097)
```

Genes are ranked by the age- and sex-adjusted MGD t statistic from the gene-level model.

The frozen representative-probe manuscript analysis gave:

```text
Response to Type II Interferon
NES = 1.702072
nominal P = 0.001674
FDR = 0.001977

Response to Cytokine
NES = 0.935811
nominal P = 0.621069
FDR = 0.602620
```

The highest-log2-variance-probe sensitivity analysis gave:

```text
Response to Type II Interferon
NES = 1.648251
nominal P = 0.006531
FDR = 0.007787

Response to Cytokine
NES = 0.911992
nominal P = 0.662004
FDR = 0.645456
```

A 10-seed tie-order sensitivity analysis using the max-log2-variance representative probes showed that the Type II interferon pathway remained positively enriched across all seeds, whereas the Response to Cytokine pathway remained non-significant.

Frozen GSEA outputs are:

```text
results/human/GSE17822_GSEA_representativeProbe_frozen.csv
results/human/GSE17822_GSEA_maxLog2VarianceProbe.csv
results/human/GSE17822_GSEA_maxLog2VarianceProbe_tieSensitivity.csv
```

### Zero-variance-filtered GSEA sensitivity

Because 2,212 of the 25,159 representative probes have exactly zero variance across the 12 samples, an additional sensitivity analysis removes all zero-variance representative probes before preranked GSEA.

The sensitivity script is:

`scripts/09_GSE17822_GSEA_zeroVariance_sensitivity.py`

After filtering:

```text
22,947 genes retained for GSEA
2,212 zero-variance gene-level probes removed
```

For the standard highest-mean representative-probe rule, Type II interferon enrichment remained significant in all 10 seeds:

```text
NES range = 1.631316 to 1.684006
maximum nominal P = 0.006818
maximum FDR = 0.010178
10 / 10 seeds FDR < 0.05
```

For the highest-log2-variance representative-probe rule, Type II interferon enrichment also remained significant in all 10 seeds:

```text
NES range = 1.543347 to 1.608861
maximum nominal P = 0.016432
maximum FDR = 0.017115
10 / 10 seeds FDR < 0.05
```

In contrast, Response to Cytokine remained non-significant under both representative-probe rules.

Outputs are:

```text
results/human/GSE17822_GSEA_zeroVarianceFiltered_10seed.csv
results/human/GSE17822_GSEA_zeroVarianceFiltered_summary.csv
```

### Interpretation of the human evidence

The human dataset provides partial rather than complete cross-species support.

The fixed 16-gene mouse-derived defense module does not show significant overall elevation in human MGD after adjustment for age and sex.

In contrast, pathway-level enrichment for Response to Type II Interferon is reproducible across representative-probe definitions, random-seed sensitivity, and explicit removal of zero-variance probes.

Accordingly, the human result is interpreted as robust pathway-level support for a Type II interferon-related defense signal, not as replication of a fixed universal gene signature.

---

## Awat2 WT-reference and internal-WT GSEA sensitivity

Because the GSE261036 Awat2 experiment contains three WT pooled libraries but only one pooled KO library, defense-pathway enrichment was evaluated against each WT reference separately as well as against the mean of all three WT libraries.

The reproducibility script is:

`scripts/10_Awat2_WT_reference_GSEA_sensitivity.py`

Frozen defense-related gene sets are stored in:

`source_data/GSE261036/Awat2_defense_gene_sets.gmt`

The three predefined pathways are:

```text
Defense Response to Virus (GO:0051607)
Response to Cytokine (GO:0034097)
Response to Type II Interferon (GO:0034341)
```

Whole-MG pseudobulk CPM values are taken from:

`results/cross_model/Awat2_wholeMG_allGene_pseudobulk_CPM.csv`

For each comparison, genes are retained when:

```text
max(CPM_A, CPM_B) >= 1
```

and ranked by:

```text
log2((CPM_A + 1) / (CPM_B + 1))
```

The following comparisons were evaluated:

```text
KO vs WT1
KO vs WT2
KO vs WT3
KO vs mean(WT1, WT2, WT3)
WT3 vs WT1
WT3 vs WT2
```

Each comparison was rerun across 10 random seeds.

The matched-run KO-versus-WT3 comparison showed significant positive enrichment of all three defense pathways across all 10 seeds. In contrast, KO-versus-WT1, KO-versus-WT2, and KO-versus-WTmean did not show stable significant enrichment.

Most importantly, WT3 itself showed strongly negative Type II interferon enrichment relative to both WT1 and WT2 across all 10 seeds:

```text
WT3 vs WT1:
Type II interferon NES range = -1.9575 to -1.8046
10 / 10 seeds FDR < 0.05

WT3 vs WT2:
Type II interferon NES range = -2.0940 to -2.0306
10 / 10 seeds FDR < 0.05
```

Therefore, the significant KO-versus-WT3 defense-pathway signal is considered reference-sensitive and is not interpreted as evidence that sequencing-run or reference-library effects were excluded.

The matched-run comparison is retained as a sensitivity analysis rather than as the principal Awat2 inflammatory-defense result.

The primary cross-model evidence instead relies on the independently defined Hsd3b6 FDR-upregulated gene set and the reference-robust Awat2 criterion requiring KO expression to exceed all three WT libraries.

Outputs are:

```text
results/awat2/Awat2_KO_vs_eachWT_and_WTmean_GSEA.csv
results/awat2/Awat2_WT_internal_negativeControl_GSEA.csv
results/awat2/Awat2_WT_reference_GSEA_rank_gene_counts.csv
results/awat2/Awat2_WT_reference_GSEA_10seed.csv
results/awat2/Awat2_WT_reference_GSEA_10seed_summary.csv
```

---

## Awat2 Leiden-resolution sensitivity

MG-state composition was evaluated using frozen Leiden labels at resolutions 0.2, 0.4, and 0.6 stored in:

`source_data/GSE261036/GSE261036_MG_Leiden_test.h5ad`

The reproducibility script is:

`scripts/11_Awat2_Leiden_resolution_sensitivity.py`

The three resolutions produced:

```text
r=0.2: 6 clusters
r=0.4: 6 clusters
r=0.6: 9 clusters
```

Adjusted Rand indices were:

```text
r0.2 vs r0.4: 0.770644
r0.4 vs r0.6: 0.662800
r0.2 vs r0.6: 0.599348
```

For interpretation at r=0.6, the nine clusters were collapsed back to the six major MG-state framework according to their dominant overlap with the frozen r=0.4 state assignments.

The principal Awat2-KO compositional directions were preserved across all three resolutions:

```text
Basal duct:
KO below WT range at r0.2, r0.4, and r0.6

Cycling MG:
KO above WT range at r0.2, r0.4, and r0.6

Differentiated duct:
KO above WT range at r0.2, r0.4, and r0.6

Early/differentiating meibocyte:
KO below WT range at r0.2, r0.4, and r0.6

Lipogenic/differentiated meibocyte:
KO within WT range at r0.2, r0.4, and r0.6
```

The Defense/stress duct-like compartment did not show a stable compositional expansion: KO remained within the WT range at r0.2 and r0.6 and was below the WT range at r0.4.

Thus, increasing Leiden resolution primarily subdivides existing MG populations without altering the principal Awat2-KO epithelial-remodeling pattern. Resolution 0.4 is retained as an intermediate representative granularity for downstream state-level analyses.

Outputs are:

```text
results/awat2/Awat2_Leiden_resolution_cluster_counts_ARI.csv
results/awat2/Awat2_Leiden_r06_vs_r04_crosstab.csv
results/awat2/Awat2_Leiden_resolution_state_composition.csv
results/awat2/Awat2_Leiden_resolution_KO_direction_summary.csv
```
