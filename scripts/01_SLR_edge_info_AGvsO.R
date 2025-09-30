# ========================================================
# Script: 01_SLR_edge_info_AGvsO.R
# Purpose:
#   Perform Sparse Logistic Regression (SLR) using gene sets derived from 
#   network analysis of glioblastoma, astrocytoma, and oligodendroglioma 
#   connections. Three strategies are tested:
#     1. All genes from JGL with binary penalization (lower penalties for genes
#        in common edges between glioblastoma and astrocytoma).
#     2. Only genes from common edges between glioblastoma and astrocytoma.
#     3. Genes from common edges plus exclusive edges from oligodendroglioma.
#
#   Key steps:
#     - Build adjacency matrices for common/exclusive edges.
#     - Define penalty vectors for Elastic Net models.
#     - Run bootstrapped SLR using different penalties.
#     - Extract consistently selected genes across bootstraps.
#     - Visualize gene overlaps and expression distributions.
#     - Perform Wilcoxon rank-sum tests to compare expression between groups.
# ========================================================

set.seed(38)
source("scripts/utils/functions.R")
source("scripts/utils/visualization.R")
load("outputs/00_preprocessing.RData")

# --------------------------------------------------------
# 1. Identify Genes in Common Edges Between Glioblastoma and Astrocytoma
# --------------------------------------------------------

# Build binary adjacency matrix for common edges (-2 in Ex_Oligo)
common_pairs_ag_matrix <- matrix(0, nrow = nrow(Ex_Oligo), ncol = ncol(Ex_Oligo))
rownames(common_pairs_ag_matrix) <- rownames(Ex_Oligo)
colnames(common_pairs_ag_matrix) <- colnames(Ex_Oligo)
common_pairs_ag_matrix[Ex_Oligo == -2] <- 1

# Create graph and extract edge list
common_pairs_ag_graph <- graph_from_adjacency_matrix(common_pairs_ag_matrix, mode = "undirected")
common_pairs_ag_edge_list <- igraph::as_data_frame(common_pairs_ag_graph, what = "edges")

# Genes in common edges
common_pairs_ag_genes <- unique(c(common_pairs_ag_edge_list[[1]], common_pairs_ag_edge_list[[2]]))
all_genes <- rownames(common_pairs_ag_matrix)  # all genes from JGL

# --------------------------------------------------------
# 2. SLR with All JGL Genes (Binary Penalization)
# --------------------------------------------------------

# Define penalty vectors (lower penalties for genes in common edges)
penalty1 <- ifelse(all_genes %in% common_pairs_ag_genes, 0.2, 1)
penalty2 <- ifelse(all_genes %in% common_pairs_ag_genes, 0.5, 1.5)
penalty3 <- ifelse(all_genes %in% common_pairs_ag_genes, 1, 1)
common_ag_binary_penalty_list <- list("0.2/1" = penalty1, "0.5/1.5" = penalty2, "1/1" = penalty3)

## 2.1 Prepare data for SLR
common_ag_binary_data <- pre_SLR_data(astro, gbm, oligo, genes = all_genes, label0 = "oligo")
Xdata <- common_ag_binary_data$Xdata
Ydata <- common_ag_binary_data$Ydata

## 2.2 Run SLR
common_pairs_ag_SLR <- run_SLR_penalized(Xdata, Ydata, common_ag_binary_penalty_list, times_boot = 100)
common_pairs_ag_SLR_results <- summarize_SLR_results(common_pairs_ag_SLR)

## 2.3 Extract consistently selected genes
common_ag_SLR_genes_100 <- consistently_selected_genes(common_pairs_ag_SLR, threshold = 100)
common_ag_SLR_genes_75  <- consistently_selected_genes(common_pairs_ag_SLR, threshold = 75)

# Visualize gene overlap (Venn diagram)
plot_venn_consistently_selected( 
  selected_genes_list = list(
    common_ag_SLR_genes_75$`0.2/1`,
    common_ag_SLR_genes_75$`0.5/1.5`,
    common_ag_SLR_genes_75$`1/1`
  ),
  category_names = c("Penalty 0.2/1", "Penalty 0.5/1.5", "Penalty 1/1"),
  plot_title = "Consistently Selected Genes - >75 Bootstraps"
)

## 2.4 Compare expression of selected genes across diseases
for (penalty_name in names(common_ag_SLR_genes_75)) {
  long_df <- prepare_long_expression_df(astro, gbm, oligo, selected_genes = common_ag_SLR_genes_75[[penalty_name]])
  plot_gene_expression_by_disease(long_df)
  print(wilcox_test_one_vs_others(long_df))
}

# --------------------------------------------------------
# 3. SLR Using Only Genes from Common Edges
# --------------------------------------------------------

## 3.1 Prepare data
only_common_ag_data <- pre_SLR_data(astro, gbm, oligo, genes = common_pairs_ag_genes, label0 = "oligo")
Xdata <- only_common_ag_data$Xdata
Ydata <- only_common_ag_data$Ydata

## 3.2 Run SLR with default penalty
default_penalty <- list(only_common_pairs = rep(1, ncol(Xdata)))
only_common_pairs_ag_SLR <- run_SLR_penalized(Xdata, Ydata, default_penalty, times_boot = 100)
only_common_pairs_ag_SLR_results <- summarize_SLR_results(only_common_pairs_ag_SLR)

## 3.3 Extract consistently selected genes
only_common_ag_SLR_genes_100 <- consistently_selected_genes(only_common_pairs_ag_SLR, threshold = 100)
only_common_ag_SLR_genes_75  <- consistently_selected_genes(only_common_pairs_ag_SLR, threshold = 75)

# Compare expression
long_df <- prepare_long_expression_df(astro, gbm, oligo, selected_genes = only_common_ag_SLR_genes_75$only_common_pairs)
plot_gene_expression_by_disease(long_df)
only_common_ag_SLR_genes_wilcox_results <- wilcox_test_one_vs_others(long_df)

# --------------------------------------------------------
# 4. SLR Using Common Edges + Exclusive Edges (Oligodendroglioma)
# --------------------------------------------------------

## 4.1 Identify exclusive edges (1 in Ex_Oligo)
exclusive_o_matrix <- matrix(0, nrow = nrow(Ex_Oligo), ncol = ncol(Ex_Oligo))
exclusive_o_matrix[Ex_Oligo == 1] <- 1
rownames(exclusive_o_matrix) <- rownames(Ex_Oligo)
colnames(exclusive_o_matrix) <- colnames(Ex_Oligo)

# Build graph and get exclusive genes
exclusive_o_graph <- graph_from_adjacency_matrix(exclusive_o_matrix, mode = "undirected")
exclusive_o_edge_list <- igraph::as_data_frame(exclusive_o_graph, what = "edges")
exclusive_o_genes <- unique(c(exclusive_o_edge_list[[1]], exclusive_o_edge_list[[2]]))

# Union of common and exclusive genes
common_ag_plus_exclusive_o_genes <- union(exclusive_o_genes, common_pairs_ag_genes)

## 4.2 Prepare data
common_ag_plus_exclusive_o_data <- pre_SLR_data(
  astro, gbm, oligo, genes = common_ag_plus_exclusive_o_genes, label0 = "oligo"
)
Xdata <- common_ag_plus_exclusive_o_data$Xdata
Ydata <- common_ag_plus_exclusive_o_data$Ydata

## 4.3 Run SLR with default penalty
default_penalty <- list(common_plus_exclusive = rep(1, ncol(Xdata)))
common_ag_plus_exclusive_o_SLR <- run_SLR_penalized(Xdata, Ydata, default_penalty, times_boot = 100)
common_ag_plus_exclusive_o_SLR_results <- summarize_SLR_results(common_ag_plus_exclusive_o_SLR)

## 4.4 Extract consistently selected genes
common_ag_plus_exclusive_o_SLR_genes_100 <- consistently_selected_genes(common_ag_plus_exclusive_o_SLR, threshold = 100)
common_ag_plus_exclusive_o_SLR_genes_75  <- consistently_selected_genes(common_ag_plus_exclusive_o_SLR, threshold = 75)

# Compare expression
long_df <- prepare_long_expression_df(astro, gbm, oligo, selected_genes = common_ag_plus_exclusive_o_SLR_genes_75$common_plus_exclusive)
plot_gene_expression_by_disease(long_df)
common_ag_plus_exclusive_o_SLR_genes_wilcox_results <- wilcox_test_one_vs_others(long_df)


# ========================================================
# Saving Results
# ========================================================

save.image(file = "results/01_SLR_edge_info_AGvsO")


