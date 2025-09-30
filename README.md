# Multi-Disease-Biomarker-Discovery-for-Drug-Repurposing-
R scripts — Thesis “Multi-Disease Biomarker Discovery for Drug Repurposing in Glioma”. These scripts comprise the codebase used to run all analyses in the thesis.

## Folder structure
- `scripts/` — R code (utils, visualization, analysis steps `00_`→`05_`)
- `outputs/`, `results/` — generated artifacts
- `NetworkAnalysisPython/` — auxiliary Python scripts for tirnagle count
- `data/` — (empty in git) see `data/README.md` to fetch inputs

## Quickstart
```r
# R env (optional but recommended)
install.packages("renv")
renv::init(); renv::restore()   # after lockfile exists

# Run the pipeline (example)
source("scripts/00_preprocessing.R")
source("scripts/03_SLR_GeNeIV_multiDistance_AGvsO.R")
source("scripts/04_SLR_twiner_AGvsO.R")
source("scripts/05_Gene_Genes.R")
