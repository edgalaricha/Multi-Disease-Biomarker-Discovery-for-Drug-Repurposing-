# ========================================================
# Script: 04_Twiner_AGvsO.R
# Purpose:
#   Run a "Twiner" Sparse Logistic Regression (SLR) and a
#   baseline Elastic Net using the SAME 561 genes used in the
#   distance-based models (i.e., genes that are not isolated
#   simultaneously in Astro and GBM). This keeps results
#   directly comparable across methods.
#
#   Twiner idea:
#     - Build correlation profiles (gene x gene) for Astro and GBM.
#     - For each gene, compute the angular distance between the two
#       correlation vectors:  ang = acos(cosine(Astro, GBM))/pi ∈ [0,1].
#     - Use ang as a penalty weight (normalized) → lower ang = lower penalty.
#
#   Steps:
#     1) Load data & utilities
#     2) Fix the gene universe (561 comparable genes)
#     3) Log2 + scale per gene
#     4) Correlation matrices (Astro, GBM)
#     5) Build SLR data (Oligo vs Astro+GBM)
#     6) Angular weights (normalize) + add Elastic Net baseline
#     7) Run bootstrapped SLR (twiner + elastic_net)
#     8) Diagnostics (PR AUC; selection frequencies; Venn)
#     9) Ridge-on-selected genes + coeffs vs weights
#    10) Weight distribution with vertical lines
#    11) Expression plots + Wilcoxon
#    12) Heatmap of selected genes
#    13) Save artifacts
# ========================================================

# ---------------------------
# 1) Load data & utilities
# ---------------------------
set.seed(138)
source("scripts/utils/functions.R")
source("scripts/utils/visualization.R")

# Provides: astro, gbm, oligo, common_genes_after_filter, networkwt_vector_raw_df, etc.
load("outputs/00_preprocessing.RData")

# Contains distance objects and the filtered network vector (used to define the 561 genes)
load("outputs/02_SLR_GeNeIV_AGvsO.RData")

# ---------------------------
# 2) Fix the gene universe (561 comparable genes)
#     - Same genes used across other models (not isolated in Astro & GBM)
# ---------------------------
genes_subset <- rownames(networkwt_vector_raw_df)  # 561 genes

xastro_less <- astro[, genes_subset, drop = FALSE]
xgbm_less   <- gbm[,   genes_subset, drop = FALSE]
xoligo_less <- oligo[, genes_subset, drop = FALSE]

# ---------------------------
# 3) Log2 + scale (per gene)
# ---------------------------
astro_norm <- scale(log2(xastro_less + 1))
gbm_norm   <- scale(log2(xgbm_less + 1))
# (Oligo not needed for Twiner weights; only for SLR later)

# ---------------------------
# 4) Correlation matrices (genes x genes)
#    Each column = correlation profile of a gene to all other genes
# ---------------------------
xastro_cor <- as.data.frame(cor(astro_norm, method = "pearson"))
xgbm_cor   <- as.data.frame(cor(gbm_norm,   method = "pearson"))

# ---------------------------
# 5) Build SLR data (Oligo vs Astro+GBM)
# ---------------------------
twiner_SLR_data <- pre_SLR_data(astro, gbm, oligo, genes_subset, label0 = "oligo")
Xdata <- twiner_SLR_data$Xdata
Ydata <- twiner_SLR_data$Ydata

# ---------------------------
# 6) Angular weights (normalize) + Elastic Net baseline
#    ang = acos(cosine(Astro_profile, GBM_profile))/pi in [0,1]
# ---------------------------
ang_weight <- numeric(ncol(Xdata))
for (i in seq_len(ncol(Xdata))) {
  ang_weight[i] <- acos(cosine(xastro_cor[, i], xgbm_cor[, i])) / pi
}

# Keep all 561 genes (no extra thresholding) to match other models
Xdata_less <- Xdata

# Normalize weights to [0,1]
twiner_weights <- ang_weight / max(ang_weight)

penalty_list <- list(
  twiner      = twiner_weights,
  elastic_net = rep(1, ncol(Xdata_less)) # baseline (no custom penalty)
)

# ---------------------------
# 7) Run bootstrapped SLR (twiner + elastic_net)
# ---------------------------
twiner_SLR <- run_SLR_penalized(Xdata_less, Ydata, penalty_list, times_boot = 100)

# ---------------------------
# 8) Diagnostics: PR AUC (TEST), selection frequency, Venn
# ---------------------------

# PR AUC distributions on TEST sets
plot_test_auc_boxplots(
  models = list(
    "Twiner"      = twiner_SLR$twiner,
    "Elastic Net" = twiner_SLR$elastic_net
  ),
  color_set = c("#66C2A5", "#E78AC3")
)

# Summaries + consistently selected genes (≥75%)
twiner_SLR_results <- summarize_SLR_results(twiner_SLR)

genes_75 <- consistently_selected_genes(twiner_SLR, threshold = 75)
twiner_genes_75      <- genes_75$twiner
elastic_net_genes_75 <- genes_75$elastic_net

# Selection counts for both models (restricted to union of selected genes)
genes_to_plot <- unique(c(elastic_net_genes_75, twiner_genes_75))

plot_selection_counts_multi(
  models_list = list(
    Elastic = twiner_SLR$elastic_net,
    Twiner  = twiner_SLR$twiner
  ),
  genes  = genes_to_plot,
  order  = "asc",
  title  = NULL,
  palette = c("#66C2A5", "#E78AC3")
)

# Venn (Twiner vs Elastic Net)
plot_venn_consistently_selected_2(
  selected_genes_list = list(elastic_net_genes_75, twiner_genes_75),
  category_names = c("Elastic Net", "Twiner"),
  plot_title = ""
)

# ---------------------------
# 9) Ridge on selected genes + coeffs vs weights
# ---------------------------

# Build a small weights data.frame for convenience
weights_df <- data.frame(gene = genes_subset, weight = twiner_weights)
rownames(weights_df) <- weights_df$gene
weights_df$gene <- NULL

# Ridge using Twiner selected genes
reg_data <- pre_SLR_data(astro, gbm, oligo, genes = twiner_genes_75, label0 = "oligo")
twiner_ridge_model <- cv.glmnet(x = reg_data$Xdata, y = reg_data$Ydata,
                                   family = "binomial",   
                                   alpha = 0,             
                                   nfolds = 10,
                                   type.measure = "mse")

# Coefficients 
twiner_ridge_model_beta <- as.matrix(coef(twiner_ridge_model, s = twiner_ridge_model$lambda.min ))


# Coefficients vs Twiner weights
best_twiner_genes_weights <- weights_df[twiner_genes_75, , drop = FALSE]

twiner_coeff_vs_weight <- data.frame(
  Gene        = twiner_genes_75,
  Weight      = as.numeric(best_twiner_genes_weights[, 1]),
  Coefficient = twiner_ridge_model_beta[-1]
)

plot_weights_vs_coefficients(
  df_plot = twiner_coeff_vs_weight,
  elastic_genes = elastic_net_genes_75,
  set_color         = "#E78AC3"
)

# ---------------------------
# 10) Weight distribution + vertical lines (selected genes)
# ---------------------------
plot_weight_hist_with_gene_lines(
  gene_scores_df   = weights_df,           # contains column 'weight'
  genes_lines      = twiner_genes_75,
  other_set_genes  = elastic_net_genes_75,
  weight_col       = "weight",
  binwidth         = 0.02,
  title            = NULL,
  x_label          = "Twiner Weights",
  show_legend      = FALSE,
  set_a_name       = "Only Twiner",
  set_b_name       = "Also in Elastic Net",
  base_line_color  = "#E78AC3"
)

# ---------------------------
# 11) Expression plots + Wilcoxon (selected by Twiner)
# ---------------------------
twiner_long_df <- prepare_long_expression_df(astro, gbm, oligo, twiner_genes_75)
plot_gene_expression_by_disease(twiner_long_df)
print(wilcox_test_one_vs_others(twiner_long_df))

# ---------------------------
# 12) Heatmap of selected genes (z-score per gene; GBM, Astro, Oligo order)
# ---------------------------
X_list <- build_heatmap_matrix(astro, gbm, oligo, genes = twiner_genes_75, zscore_rows = TRUE)

plot_gene_heatmap(
  X_list,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  use_brewer     = TRUE,
  brewer_name    = "BrBG",
  brewer_direction = 1
)

# ---------------------------
# 13) Save artifacts
# ---------------------------
save(
  twiner_SLR,
  twiner_SLR_results,
  twiner_genes_75,
  file = "outputs/04_SLR_twiner_AGvsO.RData"
)

save.image(file = "results/04_SLR_twiner_AGvsO.RData")
