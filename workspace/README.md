# Runtime workspace

Large data and generated intermediate files are intentionally excluded from Git.

At minimum, place the processed public source file here:

```text
GSE294965_processed_data.h5ad
```

The source dataset is available under **GEO GSE294965**.

The analysis notebooks generate/use additional intermediates in this directory, including:

- `roi_782_PC1_primary.h5ad`
- `roi782_cell_metadata.pkl`
- `figure4_driver_neighbor_data_smoothed.csv`
- `figure4_driver_results_smoothed.csv`
- `figure4_smoothed_k_sensitivity.csv`
- `figure4_smoothed_GBM_LOO.csv`
- `figure4_final_model_specification_robustness.csv`
- `figure5_MAC_FIB_roi_mean_counts.csv`
- `figure5_MAC_FIB_trajectory_modules.csv`
- `figure5_frozen_module_scores.csv`
- `figure5_maintext_pathway_enrichment.csv`

Do not commit the large Xenium source file to GitHub.
