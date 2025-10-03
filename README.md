# Multi-Disease-Biomarker-Discovery-for-Drug-Repurposing-
R scripts — Thesis “Multi-Disease Biomarker Discovery for Drug Repurposing in Glioma”. These scripts comprise the codebase used to run analyses in the thesis.

## Folder structure

- `scripts/` — Analysis pipeline (R)
  - `00_preprocessing.R` — Normalization and Joint Graphical Lasso
  - `01_SLR_edge_info_AGvsO.R` — Sparse Logistic Regression (SLR) with edge info penalization.
  - `02_SLR_GeNeIV_AGvsO.R` — GeNeIV-single-distance based SLR.
  - `03_SLR_GeNeIV_multiDistance_AGvsO.R` — GeNeIV-multi-distance based SLR.
  - `04_SLR_twiner_AGvsO.R` — Twiner based SLR.
  - `05_gene_selection_AGvsO.R` — Consolidates selection across models, Differentially Gene Expression and Enrichment Analysis.
  - `06_drugs_search_AGvsO.R` — Maps selected genes to DGIdb.

- `scripts/utils/` — Documented helpers for analysis
  -`functions.R` - All reusable functions to run analyses end-to-end: preprocessing, model fitting, aggregation...
  -`visualization.R` - Plotting utilities only.
- `outputs/` —
  - `graph_adjmatrices/` — Adjacency matrices for each disease-specific network.
  - `AGvsO_selected_genes.txt` — Final list of genes selected as biomarkers.
  - `AGvsO_drug_interaction_results.tsv` — DGIdb interaction table obtained by querying the selected genes.
  - `JGL-lam10.9-lam20.001.RData` — Saved Joint Graphical Lasso results produced in `00_preprocessing.R`.
- `NetworkAnalysisPython/` — Auxiliary Python for triangle counts.

## Information 
- **Working directory:** Run everything from the **repository root** (all script paths are relative to this folder).
- **Data:** The datasets used in this project are available in the GitHub Release **“Glioma Data”** (see the repository’s *Releases* tab).
- **Entry point:** `scripts/00_preprocessing.R` **must be run first** to prepare all downstream inputs.
- **Fast start:** `outputs/` already includes a **JGL result** so you don’t need to re-run network inference unless you want to retune λ.
- **Network features:** The **triangle-count/edge metrics** produced by `NetworkAnalysisPython/` are **required** by the `03_*` and `04_*` scripts—generate them before running those steps.


