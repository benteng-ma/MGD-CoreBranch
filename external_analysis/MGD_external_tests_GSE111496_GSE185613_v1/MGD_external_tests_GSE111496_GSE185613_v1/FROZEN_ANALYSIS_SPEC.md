# Frozen external-test specification: GSE111496 and GSE185613

## Shared principle
The 10-gene Awat2 inflammatory-defense panel predates these analyses:
Cxcl5, Cxcl1, Cxcl2, Ccl20, Slpi, S100a8, Ifitm1, Ifitm3, Il1rn, Lcn2.

These two external tests were added after the original MGD manuscript had been developed and were not prospectively registered. The gene panel is pre-existing; the exact external-analysis rules below are frozen from this version forward.

## GSE111496 — human immortalized meibomian gland epithelial cells
Scientific role: human MG-epithelial inflammatory inducibility / boundary test, not MGD disease replication.

- Primary cell type: HMGEC only.
- Design: 2×2 DHT × LPS/LBP, 3 arrays per cell/treatment condition.
- Mouse-to-human direct-symbol mapping:
  CXCL5, CXCL1, CXCL2, CCL20, SLPI, S100A8, IFITM1, IFITM3, IL1RN, LCN2.
- Deposited values: author-normalized Illumina intensities; log2 transformed for this analysis.
- Multiple probes: highest mean normalized intensity across all 12 HMGEC arrays, a treatment-effect-independent rule.
- Module: gene-wise z score across the 12 HMGEC arrays, then unweighted mean.
- Primary statistic: average LPS/LBP-minus-vehicle module difference across the DHT=0 and DHT=1 strata.
- Exact test: LPS labels permuted within each DHT stratum; 20×20=400 possible allocations.
- Sensitivities: leave-one-gene-out and maximum-variance probe rule.
- DHT×LPS interaction: secondary linear-model context only.

Primary result:
- targeted module effect = 0.437497
- restricted exact two-sided P = 0.015000
- DHT=0 stratum effect = 0.329497
- DHT=1 stratum effect = 0.545498
- interaction P = 0.407286
- leave-one-gene-out P range = 0.005000 to 0.045000
- maximum-variance probe-rule sensitivity P = 0.025000

Interpretation boundary:
This supports targeted inflammatory-defense inducibility in human MG epithelial cells under LPS/LBP. It does not establish MGD disease replication, IFN-gamma causality, or a universal fixed-gene signature.

## GSE185613 — HFD/circadian mouse MG-eyelid dataset
Scientific role: metabolic-context boundary test.

- 12 samples: 2 control + 2 HFD at each of 9h, 15h and 21h.
- The deposited expression file contains non-integer normalized values; GEO metadata state that FPKM standardization was used.
- All 10 fixed mouse-panel genes were directly represented.
- Transform: log2(expression+1).
- Module: gene-wise z score across all 12 samples, then unweighted mean.
- Primary effect: HFD coefficient adjusted for time-of-day.
- Exact test: HFD/control labels permuted only within each four-sample time stratum; 6^3=216 possible allocations.
- Sensitivities: leave-one-gene-out and percentile-rank score.
- Important tissue caveat: GEO overall design states that upper and lower eyelids were removed in whole for RNA extraction; this is therefore less MG-specific than purified/acinar datasets.
- Metadata duration caveat: GEO series-level text states 4 months of diet, while sample-level records in the uploaded SOFT file state 3 months. Do not silently resolve this inconsistency.

Primary result:
- time-adjusted HFD module effect = -0.563561
- restricted exact two-sided P = 0.296296
- time-specific effects: 9h -1.023, 15h 0.175, 21h -0.842
- percentile-rank sensitivity effect = -0.166667, P = 0.287037
- all leave-one-gene-out effects remained negative; none reached two-sided P<0.05.

Interpretation boundary:
GSE185613 does not support transfer of this fixed 10-gene defense panel. The result should be treated as a negative boundary result rather than rescued by post-hoc gene replacement.
