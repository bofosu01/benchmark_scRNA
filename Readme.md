# Benchmarking scRNA-seq DE Methods Under Sample-Size Imbalance

Benchmarking four differential-expression (DE) methods on monocytes from the COVID-19 PBMC dataset GSE158055, using a controlled gradient of biological sample-size imbalance (1:1:1 → 1:1:10) in a three-group design (control / mild-moderate / severe-critical).

## Project status

- [x] Data downloaded (GSE158055 counts, barcodes, features, annotation, metadata)
- [x] Dataset characterization (batches, groups, samples per patient)
- [x] Package installation
- [x] Real-data loading with streaming mtx reader (in progress)
- [ ] Preprocessing (normalization + cell-count downsampling)
- [ ] Method wrappers (Wilcoxon, MAST ± donor RE, NEBULA)
- [ ] Splatter simulation
- [ ] Benchmark driver
- [ ] Metrics (FDR, TPR, Jaccard, runtime)

## Methods compared

| Method                     | Package                | Models donor dependence |
|----------------------------|------------------------|-------------------------|
| Wilcoxon rank-sum          | base R (`wilcox.test`) | No                      |
| MAST (no random effect)    | `MAST`                 | No                      |
| MAST (donor random effect) | `MAST` + `glmer`       | Yes                     |
| NEBULA                     | `nebula`               | Yes                     |

## Data

- **Source**: GSE158055 (Ren et al., 2021)
- **Cell type**: monocytes (`majorType == "Mono"`)
- **Core batches**: Batch02, Batch03, Batch09, Batch15 (the only batches with all three disease groups)
- **Samples**: one per donor (`sampleID`)
- **Group label**: from `samples_per_patient.csv` (`control`, `mild/moderate`, `severe/critical`)

Raw input files are **not tracked** in this repo (the mtx is 7.6 GB). Download them from GEO into `dataset/`:

## Setup

### R dependencies

`Rscript scripts/install_packages.R`

CRAN: ggplot2, dplyr, tidyr, readr, purrr, stringr, Matrix, data.table, BiocManager, readxl, remotes

Bioconductor (devel): MAST, Splatter, SingleCellExperiment, scater

GitHub: nebula (via remotes::install_github("lhe17/nebula"))

### Running the pipeline

# 1. Load monocyte subset from GSE158055

nohup Rscript scripts/load_real_data.R \> nohup.out 2\>&1 &
