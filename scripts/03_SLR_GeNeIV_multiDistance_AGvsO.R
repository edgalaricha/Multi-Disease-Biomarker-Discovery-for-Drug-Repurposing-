# ========================================================
# Script: 03_GeNeIV_multiDistance_AGvsO.R
# Approach: GeNeIV (Gene Network Identity Vector)
# Purpose:
#   Perform network-informed Sparse Logistic Regression (SLR) using gene-level 
#   network measures from glioblastoma, astrocytoma, and oligodendroglioma graphs.
#
#   The analysis uses custom score–based penalty factors (linear, ratio, and log)
#   derived from distance metrics between gene-level network vectors.
#
#   Key steps:
#     - Load precomputed distance matrices (Euclidean, Manhattan, Bray-Curtis, Canberra, Ratio).
#     - Add custom scores (linear, ratio, log) to each distance matrix.
#     - Build penalty vectors from custom scores.
#     - Run bootstrapped SLR (Elastic Net) for each configuration.
#     - Summarize and combine results across all models.
# ========================================================


# ---------------------------
# 1) Reproducibility + inputs
#    - Loads preprocessing objects and distance matrices for AGvsO
#    - Utility functions for SLR and visualization
# ---------------------------
set.seed(138)
load("outputs/00_preprocessing.RData")
load("outputs/02_SLR_GeNeIV_AGvsO.RData")
source("scripts/utils/functions.R")
source("scripts/utils/visualization.R")

# ========================================================
# 2) PREPARE DATA FOR SLR
#    - Build X/Y using genes present in score dataframes
#    - Classification: Oligo vs (Astro + GBM)
# ========================================================

networkwt_distances_data <- pre_SLR_data(
  astro, gbm, oligo,
  genes = rownames(distancewt_euclidean_df),
  label0 = "oligo"
)
Xdata <- networkwt_distances_data$Xdata
Ydata <- networkwt_distances_data$Ydata

# ========================================================
# 3) BASELINE ELASTIC NET (NO CUSTOM PENALTY)
#    - Bootstrapped SLR with uniform penalty factor
#    - Serves as a reference model against distance-informed penalties
# =======================================================

# Elastic baseline
default_penalty_vector <- rep(1, nrow(distancewt_euclidean_df))
elastic_net_penalty <- list(elastic_net = default_penalty_vector)
elastic_net_SLR <- run_SLR_penalized(Xdata, Ydata, elastic_net_penalty, times_boot = 100)

# ========================================================
# 4) ADD CUSTOM SCORES TO DISTANCE MATRICES
#    - GeNeIV parameter tuning:
#      For each distance (Euclidean/Canberra/Bray/Manhattan), try multiple
#      custom-score equations (linear / log / ratio) and pick best
#      λ/α settings via helper routines.
# ========================================================

# ---------------------------
# Euclidean distance → tune penalties → collect best
# ---------------------------
networkwt_euclidean_SLR_all <- run_parameter_tuning_SLR(
  distance_df  = distancewt_euclidean_df,
  primary      = 1,                     
  distance_name= "euclidean",
  Xdata = Xdata, Ydata = Ydata)

networkwt_euclidean_best_measures <- pick_best_model_per_equation(networkwt_euclidean_SLR_all$summary)
network_euclidean_SLR <- collect_SLR_details(networkwt_euclidean_best_measures, networkwt_euclidean_SLR_all) 

# ---------------------------
# Canberra distance → tune penalties → collect best
# ---------------------------
networkwt_canberra_SLR_all <- run_parameter_tuning_SLR(
  distance_df  = distancewt_canberra_df,
  primary      = 1,                         
  distance_name= "canberra",
  Xdata = Xdata, Ydata = Ydata)

networkwt_canberra_best_measures <- pick_best_model_per_equation(networkwt_canberra_SLR_all$summary)
network_canberra_SLR <- collect_SLR_details(networkwt_canberra_best_measures, networkwt_canberra_SLR_all) 

# ---------------------------
# Bray–Curtis distance → tune penalties → collect best
# ---------------------------
networkwt_bray_SLR_all <- run_parameter_tuning_SLR(
  distance_df  = distancewt_bray_df,
  primary      = 1,                         
  distance_name= "bray",
  Xdata = Xdata, Ydata = Ydata)

networkwt_bray_best_measures <- pick_best_model_per_equation(networkwt_bray_SLR_all$summary)
network_bray_SLR <- collect_SLR_details(networkwt_bray_best_measures, networkwt_bray_SLR_all) 

# ---------------------------
# Manhattan distance → tune penalties → collect best
# ---------------------------
networkwt_manhattan_SLR_all <- run_parameter_tuning_SLR(
  distance_df  = distancewt_manhattan_df,
  primary      = 1,                         
  distance_name= "manhattan",
  Xdata = Xdata, Ydata = Ydata)

networkwt_manhattan_best_measures <- pick_best_model_per_equation(networkwt_manhattan_SLR_all$summary)
network_manhattan_SLR <- collect_SLR_details(networkwt_manhattan_best_measures, networkwt_manhattan_SLR_all) 

# ---------------------------
# Aggregate best summaries across distances + baseline
# ---------------------------
summary_best_models_final <- rbind(summarize_SLR_results(elastic_net_SLR), 
                                   networkwt_bray_best_measures,
                                   networkwt_canberra_best_measures,
                                   networkwt_euclidean_best_measures,
                                   networkwt_manhattan_best_measures)

# ========================================================
# 5) TUNING DIAGNOSTICS
#    - Plot performance (e.g., Median PR AUC Test) per equation across distances
#    - Visual check for which score type works best
# ========================================================
all_models_summary <- rbind(networkwt_euclidean_SLR_all$summary,
                            networkwt_manhattan_SLR_all$summary,
                            networkwt_bray_SLR_all$summary,
                            networkwt_canberra_SLR_all$summary)

tuning_plots <- plot_tuning_by_equation(all_models_summary, metric = "Median_PR_AUC_Test")

# Show tuning panels
tuning_plots$linear
tuning_plots$log
tuning_plots$ratio



# ========================================================
# 6) KEEP BEST MODEL OBJECTS 
#    - These lists store the selected model 
# ========================================================
best_manhattan_SLR_model   <- network_manhattan_SLR
best_euclidean_SLR_model <- network_euclidean_SLR
best_elastic_SLR_model     <- elastic_net_SLR


# ========================================================
# 7) QUICK AUC COMPARISON (BEST MODELS)
#    - PR AUC distributions across selected best models (box + jitter)
# ========================================================
plot_test_auc_boxplots(list("Manhattan" = best_manhattan_SLR_model$linear,
                            "Euclidean" = best_euclidean_SLR_model$log,
                            "Elastic Net" = best_elastic_SLR_model$elastic_net))

# ========================================================
# 8) CONSISTENTLY SELECTED GENES (≥ 75 bootstraps)
#    - Recomputed explicitly for clarity (same as above)
# ========================================================

best_manhattan_genes_75 <- consistently_selected_genes(best_manhattan_SLR_model, threshold = 75)$linear
best_euclidean_genes_75 <- consistently_selected_genes(best_euclidean_SLR_model)$log
best_elastic_genes_75   <- consistently_selected_genes(best_elastic_SLR_model, threshold = 75)$elastic_net

# ========================================================
# 9) SELECTION COUNTS ACROSS MODELS (barplots per gene)
#    - Compare how many bootstraps selected each gene in each model
# ========================================================

genes_to_plot <- unique(c(best_manhattan_genes_75, best_euclidean_genes_75, best_elastic_genes_75))

plot_selection_counts_multi(
  models_list = list(
    Euclidean = best_euclidean_SLR_model$log,
    Manhattan = best_manhattan_SLR_model$linear,
    Elastic   = best_elastic_SLR_model$elastic_net
  ),
  genes = genes_to_plot,
  order = "asc",
  title = "Genes Selected Counts Across Models"
)

# ========================================================
# 10) OVERLAP (Venn diagram)
#    - Visual overlap among stable gene sets (Euclidean / Manhattan / Elastic)
# ========================================================

plot_venn_consistently_selected_3(
  selected_genes_list = list(
    best_euclidean_genes_75,
    best_manhattan_genes_75,
    best_elastic_genes_75
  ),
  category_names = c("Euclidean", "Manhattan", "Elastic Net"),
  plot_title = "Venn Diagram of Genes Selected in >75 Bootstraps"
)

# ========================================================
# 11) EXTRACT NETWORK VECTORS FOR SELECTED GENES
#    - Convenience data.frames for downstream summaries (if needed)
# ========================================================

best_manhattan_genes_df <- networkwt_vector_raw_df[rownames(networkwt_vector_raw_df) %in% best_manhattan_genes_75, ]
best_euclidean_genes_df <- networkwt_vector_raw_df[rownames(networkwt_vector_raw_df) %in% best_euclidean_genes_75, ]
best_elastic_genes_df   <- networkwt_vector_raw_df[rownames(networkwt_vector_raw_df) %in% best_elastic_genes_75, ]

# ========================================================
# 12) RIDGE (logistic) USING THOSE SELECTED GENE SETS
#    - Evaluate ridge GLM on each stable set to inspect coefficients later
# ========================================================

reg_data <- pre_SLR_data(astro, gbm, oligo, genes = best_manhattan_genes_75, label0 = "oligo")
manhattan_ridge_model <- cv.glmnet(x = reg_data$Xdata, y = reg_data$Ydata,
                                     family = "binomial",   
                                     alpha = 0,             
                                     nfolds = 10,
                                     type.measure = "mse")

# Coefficients 
manhattan_ridge_model_beta <- as.matrix(coef(manhattan_ridge_model, s = manhattan_ridge_model$lambda.min ))



reg_data <- pre_SLR_data(astro, gbm, oligo, genes = best_euclidean_genes_75, label0 = "oligo")
euclidean_ridge_model <- cv.glmnet(x = reg_data$Xdata, y = reg_data$Ydata,
                                   family = "binomial",   
                                   alpha = 0,             
                                   nfolds = 10,
                                   type.measure = "mse")

euclidean_ridge_model_beta <- as.matrix(coef(euclidean_ridge_model, s = euclidean_ridge_model$lambda.min ))


# ========================================================
# 13) COEFFICIENTS vs WEIGHTS (scatter)
#    - Attach the best (winning) penalty columns to distance dataframes
#      and align coefficients vs corresponding weights
# ========================================================

# Inject best penalty columns (per-equation winners) into distance tables
distancewt_manhattan_df <- add_best_penalties_to_df(networkwt_manhattan_SLR_all, networkwt_manhattan_best_measures,
                                                    distancewt_manhattan_df)

distancewt_euclidean_df <- add_best_penalties_to_df(networkwt_euclidean_SLR_all, networkwt_euclidean_best_measures,
                                                    distancewt_euclidean_df)

# Build coefficient-vs-weight frames for plotting
manhattan_coeff_vs_weight <- data.frame(
  Gene        = best_manhattan_genes_75,
  Weight      = distancewt_manhattan_df[best_manhattan_genes_75, "linear"],
  Coefficient = manhattan_ridge_model_beta[-1]
)

euclidean_coeff_vs_weight <- data.frame(
  Gene        = best_euclidean_genes_75,
  Weight      = distancewt_euclidean_df[best_euclidean_genes_75, "log"],
  Coefficient = euclidean_ridge_model_beta[-1]
)

# Scatterplots (highlight overlap with Elastic Net genes)
plot_weights_vs_coefficients(manhattan_coeff_vs_weight, best_elastic_genes_75, show_legend = FALSE)
plot_weights_vs_coefficients(euclidean_coeff_vs_weight, best_elastic_genes_75, set_color = "steelblue")

# ========================================================
# 14) WEIGHT DISTRIBUTIONS + VERTICAL LINES FOR SELECTED GENES
#    - Histogram + density of all weights with dashed lines for selected genes
# ========================================================

plot_weight_hist_with_gene_lines(
  gene_scores_df     = distancewt_manhattan_df,
  genes_lines        = best_manhattan_genes_75,
  other_set_genes    = best_elastic_genes_75,
  weight_col         = "linear",
  binwidth           = 0.02,
  title              = NULL,
  x_label            = "GeNeIV Manhattan-Linear Weights",
  show_legend        = FALSE,
  set_a_name         = "Only Manhattan",
  set_b_name         = "Also in Elastic Net",

)

plot_weight_hist_with_gene_lines(
  gene_scores_df     = distancewt_euclidean_df,
  genes_lines        = best_euclidean_genes_75,
  other_set_genes    = best_elastic_genes_75,
  weight_col         = "log",
  binwidth           = 0.02,
  title               = NULL,
  x_label            = "GeNeIV Euclidean-Log Weights",
  show_legend        = FALSE,
  base_line_color    = "#8DA0CB",
  set_a_name         = "Only Euclidean",
  set_b_name         = "Also in Elastic Net"
)

# ========================================================
# 15) EXPRESSION PLOTS + WILCOXON (by gene set)
#    - For each stable set, visualize expression by disease and test contrasts
# ========================================================

manhattan_long_df <- prepare_long_expression_df(astro, gbm, oligo, best_manhattan_genes_75)
plot_gene_expression_by_disease(manhattan_long_df)
print(wilcox_test_one_vs_others(manhattan_long_df))

euclidean_long_df <- prepare_long_expression_df(astro, gbm, oligo, best_euclidean_genes_75)
plot_gene_expression_by_disease(euclidean_long_df)
print(wilcox_test_one_vs_others(euclidean_long_df))

# ========================================================
# 16) CENTRALITY MEASURE HEATMAP (Top-20%)
#    - Build top-k% indicator per metric/disease and plot for union of genes
# ========================================================

centrality_top <- build_top_indicator_matrix(
  networkwt_vector_raw_df,
  diseases = c("Astro","GBM","Oligo"),
  metrics  = c("degree","betweenness","harmonic","eigen","triangles"),
  top_frac = 0.2
)

genes <- unique(c(best_euclidean_genes_75, best_manhattan_genes_75))

plot_top_heatmap(
  centrality_top$binary,
  genes_to_plot  = genes,
  highlight_genes = best_elastic_genes_75,
  diseases_order = c("GBM","Astro","Oligo"),
  metrics_order  = c("degree","betweenness","harmonic","eigen","triangles"),
  sort_by_hits   = FALSE
)

# ========================================================
# 17) SUBNETWORK ANALYSIS 
#    - Convert precision matrices to graphs
#    - Extract connected components (modules)
#    - Filter modules containing our genes
#    - Summarize presence by module and plot
# ========================================================

# Convert precision matrices to graphs
astro_graph <- graph_from_adjacency_matrix(Theta_a, mode = "undirected", diag = FALSE)
gbm_graph <- graph_from_adjacency_matrix(Theta_g,mode = "undirected", diag = FALSE)
oligo_graph <- graph_from_adjacency_matrix(Theta_o,mode = "undirected", diag = FALSE)

astro_subnetworks <- components(astro_graph)
gbm_subnetworks <- components(gbm_graph)
oligo_subnetworks <- components(oligo_graph)

astro_subnetworks <- extract_genes_from_subnetworks(astro_subnetworks)
gbm_subnetworks <- extract_genes_from_subnetworks(gbm_subnetworks)
oligo_subnetworks <- extract_genes_from_subnetworks(oligo_subnetworks)

# For each network, locate our genes of interest
astro_filtered <- filter_subnetworks(astro_subnetworks, genes_to_plot)
gbm_filtered <- filter_subnetworks(gbm_subnetworks, genes_to_plot)
oligo_filtered <- filter_subnetworks(oligo_subnetworks, genes_to_plot)

astro_filtered <- extract_module_graphs(astro_graph, astro_filtered)
gbm_filtered <- extract_module_graphs(gbm_graph, gbm_filtered)
oligo_filtered <- extract_module_graphs(oligo_graph, oligo_filtered)

# Map each gene to its module across diseases
selected_genes_subnetworks <- find_modules_for_genes(genes_to_plot, 
                                                     astro_modules = astro_filtered,
                                                     gbm_modules = gbm_filtered,
                                                     oligo_modules = oligo_filtered)

# Build presence-by-module summary (per disease)
genes_subnetworks_size_df <- build_presence_by_module(
  genes = genes,
  astro_modules = astro_filtered,
  gbm_modules   = gbm_filtered,
  oligo_modules = oligo_filtered
)

# Visualize module presence; optionally highlight Elastic Net genes
plot_presence_by_module(genes_subnetworks_size_df, highlight_genes = best_elastic_genes_75,
                        highlight_bold = TRUE)

# ========================================================
# 18) SAVE ARTIFACTS (AGvsO – GeNeIV)
#    - Save best model objects and summaries
#    - Save full workspace image for reproducibility
# ========================================================

save(best_manhattan_SLR_model,
     best_euclidean_SLR_model,
     best_elastic_SLR_model,
     best_elastic_genes_75, best_euclidean_genes_75, best_manhattan_genes_75,
     summary_best_models_final,
     file = "outputs/03_SLR_GeNeIV_multiDistance_AGvsO.RData")

save.image(file = "results/03_SLR_GeNeIV_multiDistance_AGvsO.Rdata")


