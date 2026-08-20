# MGD Core--Branch Transcriptomic Integration

## README v1.1.0 release

Reproducibility materials for:

**Cross-etiology transcriptomic integration reveals an
inflammatory-defense core and distinct epithelial remodeling branches in
meibomian gland dysfunction**

## Current release scope

This release contains reproducibility materials for directly reprocessed
analyses and manuscript-matched external validation modules.

Included modules:

-   GSE274498 aging mouse MG formal cNMF analysis
-   GSE261036 Awat2-KO MG-state and sensitivity analyses
-   GSE166784 Hsd3b6-KO to Awat2-KO cross-model concordance
-   GSE17822 human MGD analysis
-   GSE291177 rat UVB external transfer validation
-   GSE111496 human MG epithelial inducibility analysis
-   GSE185613 exploratory audit

## Repository structure

``` text
MGD_CoreBranch_CodeRelease/
├─ environment/
├─ scripts/
├─ source_data/
├─ external_analysis/
├─ results/
├─ docs/
└─ README.md
```

## External analysis modules

### GSE291177 external transfer

Path:

`external_analysis/GSE291177_external_transfer_v1/`

This module contains the independent cross-species transfer validation
of the frozen inflammatory-defense panel.

### GSE111496 human MG epithelial analysis

Path:

`external_analysis/MGD_external_tests_GSE111496_GSE185613/`

This module evaluates human MG epithelial inducibility and is
interpreted as pathway-level support rather than universal fixed-gene
replication.

### GSE185613 audit

The exploratory audit package is retained for transparency. The dataset
was excluded from primary validation because of biological design
limitations rather than statistical non-significance.

## Data availability

GitHub:

https://github.com/benteng-ma/MGD-CoreBranch

Zenodo version DOI:

To be updated after v1.1.0 archival release.

## Release history

### v1.1.0

-   Added GSE291177 external transfer module.
-   Added GSE111496 human epithelial validation module.
-   Added GSE185613 exploratory audit.
-   Updated reproducibility documentation.
-   Prepared Zenodo archival release.

## License

Original code is released under MIT License unless otherwise stated.
Dataset-derived files remain subject to original data source terms.
