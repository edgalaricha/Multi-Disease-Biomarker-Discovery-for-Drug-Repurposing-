# Multi-Disease-Biomarker-Discovery-for-Drug-Repurposing-
R scripts — Thesis “Multi-Disease Biomarker Discovery for Drug Repurposing in Glioma”. These scripts comprise the codebase used to run all analyses in the thesis.

## Folder structure
- `scripts/` — R code (utils, visualization, analysis steps `00_`→`05_`)
- `outputs/`, `results/` — generated artifacts
- `NetworkAnalysisPython/` — auxiliary Python scripts for tirnagle count

The datasets used in this project are available in the GitHub Release “Glioma Data” (see the repository’s Releases tab).
The scripts/ directory contains the R scripts for each methodology presented in the thesis; helper functions and plotting utilities are in scripts/utils/. Start by running scripts/00_preprocessing.R (this prepares all downstream inputs). The 03_* and 04_* scripts require the network outputs produced by the Python workflow in NetworkAnalysisPython/—generate those files first before running these steps.
