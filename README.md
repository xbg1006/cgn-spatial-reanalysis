# Shared remodeling in crescentic glomerulonephritis — analysis code

This repository contains the frozen computational workflow supporting the manuscript:

**Shared remodeling coexists with disease-associated cell states in crescentic glomerulonephritis**

The study is a secondary analysis of the publicly available human kidney Xenium dataset **GEO GSE294965**.  
Large source data are not redistributed in this repository.

## What this repository contains

- manuscript-relevant analysis notebooks with outputs removed;
- relative-path configuration for use on other computers;
- frozen reproduction/sensitivity result tables used for audit;
- final figure-generation code;
- documentation of known reproducibility boundaries.

The repository intentionally excludes exploratory notebooks that were not part of the frozen manuscript workflow.

## Data required

Place the processed Xenium source file in the working directory as:

```text
GSE294965_processed_data.h5ad
```

By default the working directory is:

```text
<repository>/workspace/
```

You can use a different directory by setting the environment variable `CGN_PROJECT_DIR`.

### Windows PowerShell

```powershell
$env:CGN_PROJECT_DIR="D:\path\to\CGN_workspace"
jupyter lab
```

### macOS/Linux

```bash
export CGN_PROJECT_DIR=/path/to/CGN_workspace
jupyter lab
```

Intermediate files such as `roi_782_PC1_primary.h5ad`,
`roi782_cell_metadata.pkl`, molecular tables, and spatial tables are generated
or consumed by the ordered notebooks below.

## Recommended run order

1. `notebooks/01_prepare_primary_roi_pc1.ipynb`  
   Reconstructs the primary 782-ROI crescent-associated PC1 object from the processed Xenium data.

2. `notebooks/02_spatial_neighborhood_analysis.ipynb`  
   Builds ROI-associated cell metadata, macrophage/stromal neighborhood metrics, the smoothed primary spatial endpoint,
   k sensitivity, anti-GBM leave-one-patient-out sensitivity, slide sensitivity inputs, and the 3 x 3 spline/support-window analysis.

3. `notebooks/03_composition_analysis.ipynb`  
   Reproduces the primary major-cell-composition analysis along the common crescent-associated PC1 support.

4. `notebooks/04_molecular_analysis.ipynb`  
   Generates macrophage/fibroblast ROI expression summaries, the formal gene-level screen, frozen trajectory modules,
   panel-aware enrichment outputs, and slide-adjusted module sensitivity analyses.

5. `notebooks/05_nearest_distance_sensitivity.ipynb`  
   Implements the orthogonal nearest-distance spatial sensitivity analysis.

6. `notebooks/06_reproduction_audit.ipynb`  
   Reproduces the core statistics, exact formal gene-testing family, model-specification checks,
   restricted wild-cluster bootstrap sensitivity, CLR compositional sensitivity, and master PASS/CHECK audit tables.

7. `notebooks/07_final_figures.ipynb`  
   Generates the frozen submission versions of Figures 1–4 from saved analysis outputs.

## Key statistical design

The manuscript workflow uses:

- common-support restriction along the crescent-associated PC1 axis;
- patient-balanced ROI weights;
- patient-clustered covariance;
- cubic spline models;
- global FDR control for the formal 955-gene MAC/FIB testing family;
- centered log-ratio sensitivity for cell composition;
- slide-adjusted and disease-overlapping-slide sensitivity analyses;
- leave-one-anti-GBM-patient-out checks;
- multiple spline/support-window specifications;
- small-cluster sensitivity analyses.

## CR2-HTZ small-cluster sensitivity

The exact R implementation used for the final four CR2-HTZ sensitivity endpoints is included as:

```text
scripts/12_CR2_small_cluster_sensitivity.R
```

It uses patient-balanced weighted linear models, patient-clustered CR2 covariance from `clubSandwich`,
and the HTZ small-sample joint Wald test for the three disease-by-spline interaction coefficients.
It regenerates:

- `CR2_small_cluster_results.csv`
- `CR2_small_cluster_endpoint_audit.csv`

The final frozen numerical tables are retained under `reference_results/` for comparison.

## Pathway-enrichment reproducibility

The pathway analysis used online pathway libraries. These libraries can change over time.
For that reason, the manuscript treated the pathway-enrichment state as a frozen output.
See `docs/PROVENANCE_NOTES.md`.

## Reference audit outputs

`reference_results/` contains frozen audit tables, including:

- core Figure 1–4 reproduction summaries;
- exact formal 955-gene testing-family checks;
- CLR composition sensitivity;
- molecular module slide sensitivity;
- spatial leave-one-out / slide / support-window sensitivity;
- wild-cluster bootstrap interpretation;
- final CR2-HTZ result table.

These files are included for transparency and comparison. They are not substitutes for the source data.

## Software

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate cgn-spatial
```

Exact package versions from the original workstation were not preserved in the public repository.
The environment file therefore specifies compatible package families rather than claiming an exact lockfile.

## Repository status

Public release accompanying manuscript submission.

The manuscript-relevant analysis notebooks, frozen reference outputs,
software environment files, citation metadata and MIT License are included.

A persistent archived release is available through Zenodo:

https://doi.

## Data availability

The source spatial transcriptomic dataset is publicly available through
GEO under accession **GSE294965**.

Analysis code and frozen reference outputs are available at:

https://github.com/xbg1006/cgn-spatial-reanalysis

## Citation

Citation metadata are provided in `CITATION.cff`.

If you use this code, please cite the associated manuscript.
Archived release DOI: 
