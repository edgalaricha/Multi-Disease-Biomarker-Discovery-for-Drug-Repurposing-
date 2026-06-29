# ========================================================
# Script: 00_preprocessing.R
# Purpose: preprocessing of data and fit the Joint Graphical Lasso (JGL) model
# Input: 
# Output: Inital Transcript Data, Binary Theta Matrix for each disease
# ========================================================

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

# ----------------------------------------------------------
# 2. UMAP
# ----------------------------------------------------------

#Ensure numeric matrices (keep rownames as sample IDs)
to_mat <- function(df) {
  m <- as.matrix(data.matrix(df))  # coerce to numeric, preserving colnames
  # Drop columns that are all-NA (just in case)
  m[, colSums(!is.na(m)) > 0, drop = FALSE]
}

astro_m <- to_mat(astro_data)
gbm_m   <- to_mat(gbm_data)
oligo_m <- to_mat(oligo_data)

# 1) Align genes (columns) across datasets and stack samples (rows)
common_genes <- Reduce(intersect, list(colnames(astro_m), colnames(gbm_m), colnames(oligo_m)))
# filter matrices
astro_m <- astro_m[, common_genes, drop = FALSE]
gbm_m   <- gbm_m[,   common_genes, drop = FALSE]
oligo_m <- oligo_m[, common_genes, drop = FALSE]

expr <- rbind(astro_m, gbm_m, oligo_m)  # samples x genes

type <- factor(c(
  rep("Astrocytoma",       nrow(astro_m)),
  rep("Glioblastoma",      nrow(gbm_m)),
  rep("Oligodendroglioma", nrow(oligo_m))
), levels = c("Astrocytoma","Glioblastoma","Oligodendroglioma"))

# Transform: TPM -> log1p (if already log-transformed, set use_log = TRUE)
use_log <- FALSE
norm_expr <- if (use_log) expr else log1p(expr)  # samples x genes


# PCA (on samples x genes)
pca <- prcomp(norm_expr, center = TRUE, scale. = TRUE)
npc <- min(50, ncol(pca$x))
Xpca <- pca$x[, seq_len(npc), drop = FALSE]  # samples x PCs

# UMAP on PCs
set.seed(123)
emb <- uwot::umap(
  Xpca,
  n_neighbors  = 20,
  min_dist     = 0.15,
  metric       = "cosine",
  n_components = 2,
  init         = "spectral",
  fast_sgd     = TRUE
)

# Plot
UMAP_plot <- data.frame(UMAP1 = emb[,1], UMAP2 = emb[,2], subtype = type)
# define named colors for each subtype
subtype_cols <- c(
  "Astrocytoma"       = "#DEA58C",  
  "Glioblastoma"      = "#A58CDE",  # green
  "Oligodendroglioma" = "#8CDEA5"   # blue
)

library(ggplot2)

UMAP_plot <- ggplot(UMAP_plot, aes(UMAP1, UMAP2, color = subtype)) +
  geom_point(size = 0.8, alpha = 0.9) +
  scale_color_manual(values = subtype_cols, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = 0.02)) +
  labs(x = "UMAP1", y = "UMAP2", title = NULL, color = "Subtype") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor     = element_blank(),
    plot.margin          = margin(3, 3, 3, 3, unit = "mm"),
    legend.position      = c(0.9, 0.01),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill = scales::alpha("white", 0.75), color = NA),
    legend.key           = element_blank(),
    legend.title         = element_text(size = 9),
    legend.text          = element_text(size = 8),
    # afastar títulos dos eixos:
    axis.title.x         = element_text(margin = margin(t = 6, unit = "mm")),
    axis.title.y         = element_text(margin = margin(r = 6, unit = "mm")),
    # (opcional) afastar também os rótulos dos ticks
    axis.text.x          = element_text(margin = margin(t = 2, unit = "mm")),
    axis.text.y          = element_text(margin = margin(r = 2, unit = "mm"))
  ) +
  guides(color = guide_legend(
    title.position = "top",
    override.aes   = list(size = 3, alpha = 1),
    byrow = TRUE
  ))

UMAP_plot

# ----------------------------------------------------------
# 3. Preparing for JGL
# -----------------------------------------------------------

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
