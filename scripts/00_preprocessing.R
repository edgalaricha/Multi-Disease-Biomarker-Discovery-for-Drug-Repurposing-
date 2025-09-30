# ========================================================
# Script: 01_preprocessing_JGL.R
# Purpose: preprocessing of data and fit the Joint Graphical Lasso (JGL) model
# Input: 
# Output: Inital Transcript Data, Binary Theta Matrix for each disease
# ========================================================


#setwd to where file is located
source("scripts/utils/functions.R")
source("scripts/utils/visualization.R")
load("data/glioma-RNASeq-2021-classification.RData")


# --------------------------------------------------------
# 1. Preprocessing
# --------------------------------------------------------

#Remove columns with sd == 0
astro_data <- remove_zero_sd(astro)
gbm_data <- remove_zero_sd(gbm)
oligo_data <- remove_zero_sd(oligo)


#Normalizing and filtering normally distributed variables
astro_filtered <- normalize_and_filter(astro_data)
gbm_filtered <- normalize_and_filter(gbm_data)
oligo_filtered <- normalize_and_filter(oligo_data)


#Leave each dataframe with only the genes in common
common_genes_after_filter <- Reduce(intersect, list(colnames(astro_filtered), colnames(gbm_filtered), colnames(oligo_filtered)))
astro_filtered <- astro_filtered[, common_genes_after_filter, drop = FALSE]
oligo_filtered <- oligo_filtered[, common_genes_after_filter, drop = FALSE]
gbm_filtered <- gbm_filtered[, common_genes_after_filter, drop = FALSE]

dim(astro_filtered) #254 16453
dim(gbm_filtered)   #199 16453
dim(oligo_filtered) #166 16453


# --------------------------------------------------------
# 2. Fit Joint Graphical Lasso
# --------------------------------------------------------
a <- list(astro_filtered, oligo_filtered, gbm_filtered)

#Run once and save results due to running time
#fgl.results = JGL(Y=a, penalty="fused",lambda1=0.9,lambda2=0.001) 
#save(fgl.results, file="outputs/JGL-lam10.9-lam20.001.RData") 

load("outputs/JGL-lam10.9-lam20.001.RData")
JGL_best_model <- fgl.results
plot.jgl(JGL_best_model,haslegend = TRUE)

#Network characteristics
JGL_best_model


# --------------------------------------------------------
# 3. Handling Precision Matrices for each Glioma
# --------------------------------------------------------

#Separation for each Glioma
Theta_a = fgl.results$theta[[1]]
Theta_o = fgl.results$theta[[2]]
Theta_g = fgl.results$theta[[3]]

#Removing the diagonal (auto connections)
Theta_a = Theta_a - diag(dim(Theta_a)[1]) * Theta_a
Theta_o = Theta_o - diag(dim(Theta_o)[1]) * Theta_o
Theta_g = Theta_g - diag(dim(Theta_g)[1]) * Theta_g

#Convert matrices to binary
Theta_a[Theta_a != 0] <- 1
Theta_o[Theta_o != 0] <- 1
Theta_g[Theta_g != 0] <- 1

# Function to Count isolated and non-isolated nodes
non_isolated_count <- function(g) sum(igraph::degree(g, mode = "all") > 0)
isolated_vertices <- function(g) sum(igraph::degree(g, mode = "all") == 0)

# Creating graph objets and counting nodes
astro_graph <- graph_from_adjacency_matrix(Theta_a, mode = "undirected")
non_isolated_count(astro_graph)
isolated_vertices(astro_graph)

oligo_graph <- graph_from_adjacency_matrix(Theta_o, mode = "undirected")
non_isolated_count(oligo_graph)
isolated_vertices(oligo_graph)

gbm_graph <- graph_from_adjacency_matrix(Theta_g, mode = "undirected")
non_isolated_count(gbm_graph)
isolated_vertices(gbm_graph)


# Quick visualization of each network
plot_graph_network(astro_graph, "Astrocytoma")
plot_graph_network(gbm_graph, "Glioblastoma")
plot_graph_network(oligo_graph, "Oligodendroglioma")


#Exporting Adjacency Matrice
export_to_csv(Theta_a, "outputs/graph_adjmatrices//binary_adjmatrix_astro.csv")
export_to_csv(Theta_g, "outputs/graph_adjmatrices//binary_adjmatrix_gbm.csv")
export_to_csv(Theta_o, "outputs/graph_adjmatrices//binary_adjmatrix_oligo.csv")


# --------------------------------------------------------
# 4. Identifying Connections for each type of Glioma:
#   0: No connection or connection shared with another type of glioma.
#   1: The connection only exists in the type of glioma analyzed.
#  -1: The connection only exists in another type of glioma.
#  -2: The connection is shared by the other two types of glioma, but not in the type analyzed.
# --------------------------------------------------------

Ex_Gbm = Theta_g - Theta_a - Theta_o
Ex_Astro = Theta_a - Theta_g - Theta_o
Ex_Oligo = Theta_o - Theta_a - Theta_g



# --------------------------------------------------------
# 5. Plotting Edge Networks for each glioma type
# --------------------------------------------------------

resA <- plot_overlap_edges(astro_graph, gbm_graph, oligo_graph, target = "A",
                           title = NULL)


resG <- plot_overlap_edges(astro_graph, gbm_graph, oligo_graph, target = "G",
                           title = NULL)


resO <- plot_overlap_edges(astro_graph, gbm_graph, oligo_graph, target = "O",
                           title = NULL)


# ========================================================
# Saving Results
# ========================================================

save.image(file = "results/00_preprocessing_full_results.RData")

save(common_genes_after_filter, Ex_Gbm, Ex_Oligo, Ex_Astro, 
     astro, gbm, oligo, 
     Theta_a, Theta_g, Theta_o,
     file = "outputs/00_preprocessing.RData")
