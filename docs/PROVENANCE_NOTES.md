# Provenance notes

## Figure 4 PC1 support-window rule

The frozen model-specification analysis uses:

- spline degrees of freedom: 3, 4, 5;
- within-disease PC1 trimming: 0%, 2.5%, 5% from each tail;
- for each trimming level, the LN and anti-GBM disease-specific PC1 ranges are computed separately;
- the intersection of the two ranges defines the common support used for refitting.

This yields 3 x 3 = 9 model specifications.

An earlier audit table (`reference_results/reproduction_provenance_notes.csv`) recorded this rule as
`DOCUMENTATION PENDING` because the original source cell had not yet been located at that time.
The supplied spatial-analysis notebook now contains the rule explicitly.

## Pathway libraries

Panel-aware pathway enrichment depends on external GO / Reactome / Enrichr library state.
External libraries can change, so frozen enrichment outputs should be archived with the release.

## CR2-HTZ

The exact fixed R code used for the final CR2-HTZ sensitivity analysis was recovered from the archived analysis record
and is included as `scripts/12_CR2_small_cluster_sensitivity.R`. The frozen outputs remain under
`reference_results/CR2_small_cluster_results.csv` and `reference_results/CR2_small_cluster_endpoint_audit.csv`.
