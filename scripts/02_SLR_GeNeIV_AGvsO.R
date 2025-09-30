# ========================================================
# Script: 02_GeNeIV_AGvsO.R
# Purpose:
#   Perform a network-informed Sparse Logistic Regression (SLR) analysis using gene-level 
#   network measures from glioblastoma, astrocytoma, and oligodendroglioma graphs.
#
#   The analysis is performed in two stages:
#     1. Without triangle counts:
#        - Build vectors of centrality measures (Degree, Betweenness, Closeness, Eigencentrality).
#        - Calculate distances between gene vectors across disease networks.
#        - Use these distances as weights in bootstrapped SLR models with Elastic Net.
#     2. With triangle counts:
#        - Extend the vectors with triangle counts.
#        - Repeat the process using the updated network features.
#
#   Key steps:
#     - Build igraph objects from precision matrices (JGL results).
#     - Extract network vectors per gene for each disease.
#     - Compute multiple distance metrics between disease-specific vectors.
#     - Use these distances as penalties in bootstrapped SLR.
#     - Extract consistently selected genes and compare results between both approaches.
# ========================================================

set.seed(138)
load("outputs/00_preprocessing.RData")
source("scripts/utils/functions.R")
source("scripts/utils/visualization.R")


# --------------------------------------------------------
# 1. Graph Construction and Basic Network Exploration
# --------------------------------------------------------

# Create igraph objects from precision matrices (Joint Graphical Lasso)
astro_graph <- graph_from_adjacency_matrix(Theta_a, mode = "undirected", diag = FALSE)
gbm_graph <- graph_from_adjacency_matrix(Theta_g, mode = "undirected", diag = FALSE)
oligo_graph <- graph_from_adjacency_matrix(Theta_o, mode = "undirected", diag = FALSE)

# Quick visualization of each network
plot_graph_network(astro_graph, "Astrocytoma")
plot_graph_network(gbm_graph, "Glioblastoma")
plot_graph_network(oligo_graph, "Oligodendroglioma")

# Compare basic topology across networks (number of components and sizes)
astro_modules <- components(astro_graph)
gbm_modules <- components(gbm_graph)
oligo_modules <- components(oligo_graph)

astro_modules$no # number of components
gbm_modules$no
oligo_modules$no

table(astro_modules$csize) # size distribution of components
table(gbm_modules$csize)
table(oligo_modules$csize)


# --------------------------------------------------------
# 2. Sparse Logistic Regression (SLR) without triangle information
# --------------------------------------------------------
# Build gene-level vectors of (degree, betweenness, closeness, eigencentrality) 
# for each disease and compute inter-disease distances.
# --------------------------------------------------------

## 2.1 Build gene network vectors
network_vector_raw_df <- build_network_df_from_three_graphs(astro_graph,
                                                            gbm_graph,
                                                            oligo_graph,
                                                            triangle_df = NULL)

# Remove genes with zero centrality values in both Astro and GBM (isolated nodes)
network_vector_raw_df <- remove_zero_vectors(network_vector_raw_df, diseases = c("Astro", "GBM")) # Exptected. 561 genes 

# Normalize each metric to ]0,1] range
network_vector_normalized_df <- normalize_network_df(network_vector_raw_df)


## 2.2 Compute pairwise distances between disease-specific gene vectors
# Note: Normalized vectors for scale-sensitive distances; raw for scale-invariant.
distance_euclidean_df <- compare_centrality_vectors(network_vector_normalized_df) #sem NA
distance_manhattan_df <- compare_centrality_vectors(network_vector_normalized_df, method = "manhattan") #sem NA
distance_bray_df <- compare_centrality_vectors(network_vector_normalized_df, method = "bray-curtis") #sem NA
distance_canberra_df <- compare_centrality_vectors(network_vector_normalized_df, method = "canberra") #sem NA


plot_distances_distribution(distance_euclidean_df, x_col_name = "Astro_vs_GBM", x_label = "Euclidean Distance")
plot_distances_distribution(distance_manhattan_df, x_col_name = "Astro_vs_GBM", x_label = "Manhattan Distance")
plot_distances_distribution(distance_bray_df, x_col_name = "Astro_vs_GBM", x_label = "Bray Distance")
plot_distances_distribution(distance_canberra_df, x_col_name = "Astro_vs_GBM", x_label = "Canberra Distance")



## 2.3 Convert distances into normalized weights for SLR 
distance_euclidean_weights <- get_normalized_weights(distance_euclidean_df) #561 genes
distance_manhattan_weights <- get_normalized_weights(distance_manhattan_df) #561 genes
distance_bray_weights <- get_normalized_weights(distance_bray_df) #561 genes
distance_canberra_weights <- get_normalized_weights(distance_canberra_df) #561 genes


# Combine all penalty vectors in a list
default_penalty_vector <- rep(1, nrow(distance_euclidean_df))
network_distances_penalty <- list(elastic_net = default_penalty_vector,
                                  euclidean = distance_euclidean_weights$Astro_vs_GBM,
                                  manhattan = distance_manhattan_weights$Astro_vs_GBM,
                                  bray_curtis = distance_bray_weights$Astro_vs_GBM,
                                  canberra = distance_canberra_weights$Astro_vs_GBM)





## 2.4 Prepare input for SLR
network_distances_data <- pre_SLR_data(astro, gbm, oligo,
                                       genes = rownames(distance_euclidean_weights),
                                       label0 = "oligo")
Xdata <- network_distances_data$Xdata
Ydata <- network_distances_data$Ydata

## 2.5 Run bootstrapped Sparse Logistic Regression
network_distances_SLR <- run_SLR_penalized(Xdata, Ydata, network_distances_penalty, times_boot = 100)

# Summarize results
network_distances_SLR_results <- summarize_SLR_results(network_distances_SLR)
network_distances_SLR_results


# --------------------------------------------------------
# 3. Sparse Logistic Regression with triangle counts
# --------------------------------------------------------

# Read precomputed triangle counts per gene
triangles_count <- read.csv("NetworkAnalysisPython/triangle_counts.csv", 
                            row.names = 1)



## 3.1 Build extended vectors including triangle counts
networkwt_vector_raw_df <- build_network_df_from_three_graphs(astro_graph,
                                                            gbm_graph,
                                                            oligo_graph,
                                                            triangle_df = triangles_count)

# Remove isolated genes
networkwt_vector_raw_df <- remove_zero_vectors(networkwt_vector_raw_df, diseases = c("Astro", "GBM")) #561 genes

# Normalize metrics
networkwt_vector_normalized_df <- normalize_network_df(networkwt_vector_raw_df) #confirmar


## 3.2 Compute distances between extended vectors
distancewt_euclidean_df <- compare_centrality_vectors(networkwt_vector_normalized_df) #sem NA
distancewt_manhattan_df <- compare_centrality_vectors(networkwt_vector_normalized_df, method = "manhattan") #sem NA
distancewt_bray_df <- compare_centrality_vectors(networkwt_vector_normalized_df, method = "bray-curtis") #sem NA
distancewt_canberra_df <- compare_centrality_vectors(networkwt_vector_normalized_df, method = "canberra") #sem NA


## 3.3 Convert to normalized weights
distancewt_euclidean_weights <- get_normalized_weights(distancewt_euclidean_df) #561 genes
distancewt_manhattan_weights <- get_normalized_weights(distancewt_manhattan_df) #561 genes
distancewt_bray_weights <- get_normalized_weights(distancewt_bray_df) #561 genes
distancewt_canberra_weights <- get_normalized_weights(distancewt_canberra_df) #561 genes


# Penalty vectors
default_penalty_vector <- rep(1, nrow(distancewt_euclidean_df))
networkwt_distances_penalty <- list(elastic_net_wt = default_penalty_vector,
                                    euclidean_wt = distancewt_euclidean_weights$Astro_vs_GBM,
                                    manhattan_wt = distancewt_manhattan_weights$Astro_vs_GBM,
                                    bray_curtis_wt = distancewt_bray_weights$Astro_vs_GBM,
                                    canberra_wt = distancewt_canberra_weights$Astro_vs_GBM)

## Optional: visualize distributions of distances
plot_distances_distribution(distancewt_euclidean_df, x_col_name = "Astro_vs_GBM", x_label = "Euclidean Distance")
plot_distances_distribution(distancewt_manhattan_df, x_col_name = "Astro_vs_GBM", x_label = "Manhattan Distance")
plot_distances_distribution(distancewt_bray_df, x_col_name = "Astro_vs_GBM", x_label = "Bray Distance")
plot_distances_distribution(distancewt_canberra_df, x_col_name = "Astro_vs_GBM", x_label = "Canberra Distance")


## 3.4 Prepare data for SLR
networkwt_distances_data <- pre_SLR_data(astro, gbm, oligo,
                                       genes = rownames(distance_bray_weights),
                                       label0 = "oligo")

Xdata <- networkwt_distances_data$Xdata
Ydata <- networkwt_distances_data$Ydata


## 3.5 Run bootstrapped SLR
networkwt_distances_SLR <- run_SLR_penalized(Xdata, Ydata, networkwt_distances_penalty, times_boot = 100)

# Summarize results
networkwt_distances_SLR_results <- summarize_SLR_results(networkwt_distances_SLR)
networkwt_distances_SLR_results

# --------------------------------------------------------
# 4. Compare Results
# --------------------------------------------------------

network_summary_full_results <- rbind(network_distances_SLR_results, networkwt_distances_SLR_results)
network_summary_full_results

# ========================================================
# Saving Results
# ========================================================

save(astro, gbm, oligo,
     networkwt_vector_raw_df,
     distancewt_euclidean_df, distancewt_manhattan_df, distancewt_bray_df, 
     distancewt_canberra_df, 
     file = "outputs/02_SLR_GeNeIV_AGvsO.RData")

save.image(file = "results/02_SLR_GeNeIV_AGvsO.RData")



