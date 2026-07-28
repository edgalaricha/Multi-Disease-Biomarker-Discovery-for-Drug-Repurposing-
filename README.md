# GeNeIV: A Gene Network Importance Vector for Network-Guided Penalization for multi-disease biomarker discovery in glioma.
R scripts — These scripts comprise the codebase used to run analyses in the article.

## Folder structure

- `scripts/` — Analysis pipeline (R)
  - `00_preprocessing.R` — UMAP, Normalization and Joint Graphical Lasso
  - `01_SLR_edge_info_AGvsO.R` — Sparse Logistic Regression (SLR) with binary edge info penalization.
  - `02_SLR_GeNeIV_AGvsO.R` — Triangles Calculation and GeNeIV-single-distance based SLR.
  - `03_SLR_GeNeIV_multiDistance_AGvsO.R` — GeNeIV-multi-distance based SLR.
  - `04_SLR_twiner_AGvsO.R` — Twiner based SLR.
  - `05_gene_selection_AGvsO.R` — Consolidates selection across models, Differentially Gene Expression and Enrichment Analysis.
  - `06_drugs_search_AGvsO.R` — Maps selected genes to DGIdb.

- `scripts/utils/` — Documented helpers for analysis
  - `functions.R` - All reusable functions to run analyses end-to-end: preprocessing, model fitting, aggregation...
  - `visualization.R` - Plotting utilities only.
- `outputs/` —
  - `Supplementary_Tables.xlsx/` — Adjacency matrices for each disease-specific network, Edges Info, Triangle Counts, SLR results, DEG results.
  - `AGvsO_drug_interaction_results.tsv` — DGIdb interaction table obtained by querying the selected genes.
  - `JGL-lam10.9-lam20.001.RData` — Saved Joint Graphical Lasso results produced in `00_preprocessing.R`.

## Information 
- **Data:** The datasets used in this project are available in the GitHub Release **“Glioma Data”** (see the repository’s *Releases* tab).
- **Entry point:** `scripts/00_preprocessing.R` **must be run first** to prepare all downstream inputs.
- **Fast start:** `outputs/` already includes a **JGL result** so you don’t need to re-run network inference unless you want to retune λ.
