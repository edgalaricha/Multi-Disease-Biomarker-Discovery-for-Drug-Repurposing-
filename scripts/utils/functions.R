library(huge)
library(tseries)
library(glmnet)
library(PRROC)
library(reshape2)
library(igraph)
library(gprofiler2)
library("propagate")
library("lsa")
library(JGL)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tibble)
library(edgeR)

#Function to remove columns with sd == 0
remove_zero_sd <- function(df) {
  sd_values <- sapply(seq(ncol(df)), function(ix) {sd(df[,ix])})
  df_filtered <- df[, sd_values != 0]
  return(df_filtered)
}

# Function for normalizing and filtering normally distributed variables
normalize_and_filter <- function(df) {
  df_npn <- huge.npn(as.matrix(df))
  normal_vars <- apply(df_npn, 2, function(x) {
    p_value <- jarque.bera.test(x)$p.value
    return(p_value > 0.05)
  })
  return(as.data.frame(df_npn[, normal_vars, drop = FALSE]))
}

# Function to export a matrix or data.frame to CSV
export_to_csv <- function(mat, filepath = "my_matrix.csv") {
  # Ensure the input is a data.frame
  if (!is.data.frame(mat)) {
    mat <- as.data.frame(mat)
  }
  
  # Export the file to the given directory (including rownames as gene identifiers)
  write.csv(mat, file = filepath, row.names = TRUE)
  
  message("File successfully exported to: ", filepath)
}


# ========================================================
# Function: pre_SLR_data
# Purpose: 
#   Prepare data for sparse logistic regression:
#     1. Subset only the genes present in all datasets
#     2. Combine samples from astrocytoma, GBM, and oligodendroglioma
#     3. Apply log2 transformation and scaling
#     4. Create a binary label vector for classification
#
# Inputs:
#   astro  - Data frame/matrix of astrocytoma samples
#   gbm    - Data frame/matrix of GBM samples
#   oligo  - Data frame/matrix of oligodendroglioma samples
#   genes  - Vector of gene names OR matrix/data frame (column names used for subsetting)
#   label1 - String indicating which dataset should be labeled as 1 
#            (must be "astro", "gbm", or "oligo")
#
# Outputs:
#   A list containing:
#     Xdata - Combined, log2-transformed, and scaled expression matrix
#     Ydata - Binary label vector (1 for selected dataset, 0 for the others)
# ========================================================
pre_SLR_data <- function(astro, gbm, oligo, genes, label0 = "oligo") {
  # 0. Get gene names depending on the type of 'genes'
  if (is.vector(genes)) {
    gene_names <- genes
  } else if (!is.null(colnames(genes))) {
    gene_names <- colnames(genes)
  } else {
    stop("The 'genes' input must be either a character vector or a matrix/data frame with column names.")
  }
  
  # 1. Keep only genes present in all datasets
  common_cols <- intersect(gene_names, Reduce(intersect, list(colnames(astro), colnames(gbm), colnames(oligo))))
  if (length(common_cols) == 0) stop("No common genes found between 'genes' and datasets.")
  
  # Subset the data
  astro_x <- astro[, common_cols, drop = FALSE]
  gbm_x   <- gbm[, common_cols, drop = FALSE]
  oligo_x <- oligo[, common_cols, drop = FALSE]
  
  # 2. Combine all samples into one dataset
  Xdata <- rbind(astro_x, gbm_x, oligo_x)
  
  # 3. Log2-transform and scale
  Xdata_sc <- scale(log2(Xdata + 1))
  
  # 4. Get sample sizes
  n_astro <- nrow(astro_x)
  n_gbm   <- nrow(gbm_x)
  n_oligo <- nrow(oligo_x)
  
  # 5. Create binary label vector
  if (label0 == "astro") {
    Ylabel <- c(rep(0, n_astro), rep(1, n_gbm + n_oligo))
  } else if (label0 == "gbm") {
    Ylabel <- c(rep(1, n_astro), rep(0, n_gbm), rep(1, n_oligo))
  } else if (label0 == "oligo") {
    Ylabel <- c(rep(1, n_astro + n_gbm), rep(0, n_oligo))
  } else {
    stop("label0 must be 'astro', 'gbm', or 'oligo'.")
  }
  
  message("Selected ", length(common_cols), " common genes for modeling.")
  
  # 6. Return as list
  return(list(
    Xdata = Xdata_sc,
    Ydata = Ylabel
  ))
}


# ========================================================
# Function: run_SÇR_penalized
# Purpose:
#   Perform sparse logistic regression with bootstrapping using Elastic Net regularization.
#   For each bootstrap iteration, it splits the data into train/test sets, fits a cross-validated
#   Elastic Net model (cv.glmnet) with different penalty.factor vectors, and records variable selection
#   metrics, classification errors, MSE, and PR AUC for train/test sets.
#
# Inputs:
#   Xdata                - Data matrix (samples x features)
#   Ydata                - Binary response vector (0/1)
#   penalty_factors_list - Named list of penalty.factor vectors for different feature weighting schemes
#   times_boot           - Number of bootstrap iterations (default = 100)
#
# Outputs:
#   A named list containing, for each penalty scheme:
#     nvar_selected       - Number of selected variables in each bootstrap
#     miscl_train/test    - Misclassification errors for train/test sets
#     mse_train/test      - Mean squared error for train/test sets
#     pr_auc_train/test   - Precision-Recall AUC for train/test sets
#     var_selected_idx    - Indices of selected variables per bootstrap
#     var_selected_names  - Names of selected variables per bootstrap
# ========================================================
run_SLR_penalized <- function(Xdata, Ydata, penalty_factors_list, times_boot = 100) {
  
  # Initialize lists to store results
  results_list <- list()
  
  # Create ID matrix for test sets
  new_test_id_matrix <- matrix(0, round(nrow(Xdata) * 0.25), times_boot)
  
  for (i in 1:times_boot) {
    new_test_id_matrix[, i] <- sample(1:nrow(Xdata), round(nrow(Xdata) * 0.25), replace = FALSE)
  }
  
  # Create folds for Cross Validation
  my_foldid <- sample(1:10, size = nrow(Xdata[-new_test_id_matrix[, 1], ]), replace = TRUE)
  
  # Create lists to store the same train/test sets for each bootstrapping
  train_test_splits <- vector("list", times_boot)
  
  for (i in 1:times_boot) {
    train_test_splits[[i]] <- list(
      xtrain = Xdata[-new_test_id_matrix[, i], ],
      ytrain = Ydata[-new_test_id_matrix[, i]],
      xtest  = Xdata[new_test_id_matrix[, i], ],
      ytest  = Ydata[new_test_id_matrix[, i]]
    )
  }
  
  for (penalty in names(penalty_factors_list)) {
    # Obtain the penalty vector corresponding to the penalty
    penalty_factor <- penalty_factors_list[[penalty]]
    
    # Initialize result matrices for this penalty
    nvar_selected <- numeric(times_boot)
    miscl_train   <- numeric(times_boot)
    miscl_test    <- numeric(times_boot)
    mse_train     <- numeric(times_boot)
    mse_test      <- numeric(times_boot)
    pr_auc_train  <- numeric(times_boot)
    pr_auc_test   <- numeric(times_boot)
    var_selected_idx_list   <- vector("list", times_boot)
    var_selected_names_list <- vector("list", times_boot)
    
    for (i in 1:times_boot) {
      # Use the training and test sets already defined
      new_xtrain <- train_test_splits[[i]]$xtrain
      new_ytrain <- train_test_splits[[i]]$ytrain
      new_xtest  <- train_test_splits[[i]]$xtest
      new_ytest  <- train_test_splits[[i]]$ytest
      
      # ---- PRINT: penalty, boot, p, classes in train ----
      cls <- table(new_ytrain)
      cls_str <- paste(paste(names(cls), cls, sep = "="), collapse = " ")
      cat(sprintf("[Pen=%s] Boot %d | p=%d | classes train: %s\n",
                  penalty, i, ncol(new_xtrain), cls_str))
      
      # Train Elastic Net model with specific penalty 
      fit_EN_cv <- cv.glmnet(new_xtrain, as.factor(new_ytrain), 
                             family = "binomial", nfolds = 10, alpha = 0.9, 
                             foldid = my_foldid, penalty.factor = penalty_factor, 
                             type.measure = "mse")
      
      # Select variables
      best_lambda_idx <- which.min(fit_EN_cv$cvm)
      var_selected_idx <- which(fit_EN_cv$glmnet.fit$beta[, best_lambda_idx] != 0)
      var_selected_names <- colnames(new_xtrain)[var_selected_idx]
      
      nvar_selected[i] <- length(var_selected_idx)
      var_selected_idx_list[[i]] <- var_selected_idx
      var_selected_names_list[[i]] <- var_selected_names
      
      # Make predictions in training
      pred_train <- predict(fit_EN_cv, new_xtrain, type = "response")
      miscl_train[i] <- sum(new_ytrain != round(pred_train))
      mse_train[i]   <- mean((new_ytrain - pred_train)^2)
      pr_auc_train[i] <- pr.curve(scores.class0 = as.vector(round(pred_train)), 
                                  weights.class0 = new_ytrain)$auc.integral
      
      # Make predictions in test
      pred_test <- predict(fit_EN_cv, new_xtest, type = "response")
      miscl_test[i] <- sum(new_ytest != round(pred_test))
      mse_test[i]   <- mean((new_ytest - pred_test)^2)
      pr_auc_test[i] <- pr.curve(scores.class0 = as.vector(round(pred_test)), 
                                 weights.class0 = new_ytest)$auc.integral
      
      # ---- PRINT: MSE and PR-AUC for this boot ----
      cat(sprintf("  -> MSE train=%.4f test=%.4f | PR-AUC train=%.4f test=%.4f\n",
                  mse_train[i], mse_test[i], pr_auc_train[i], pr_auc_test[i]))
    }
    
    # Save the results of this penalty
    results_list[[penalty]] <- list(
      nvar_selected     = nvar_selected,
      miscl_train       = miscl_train,
      miscl_test        = miscl_test,
      mse_train         = mse_train,
      mse_test          = mse_test,
      pr_auc_train      = pr_auc_train,
      pr_auc_test       = pr_auc_test,
      var_selected_idx  = var_selected_idx_list,
      var_selected_names= var_selected_names_list
    )
  }
  
  return(results_list)
}


# ========================================================
# Function: summarize_SLR_results
# Purpose:
#   Summarize the results of sparse logistic regression bootstrapping.
#   Computes the median values for key performance metrics across bootstrap iterations
#   for each penalty scheme, returning them in a compact data frame.
#
# Inputs:
#   results_list - A named list of results from run_sparse_logistic_regression(),
#                  where each element contains performance metrics for a penalty scheme.
#
# Outputs:
#   A data frame with one row per penalty scheme, including:
#     Penalty               - Name of the penalty scheme
#     Median_Num_Vars       - Median number of selected variables
#     Median_Misclass_Train - Median misclassification error in training sets
#     Median_Misclass_Test  - Median misclassification error in test sets
#     Median_MSE_Train      - Median mean squared error in training sets
#     Median_MSE_Test       - Median mean squared error in test sets
#     Median_PR_AUC_Train   - Median precision-recall AUC in training sets
#     Median_PR_AUC_Test    - Median precision-recall AUC in test sets
# ========================================================
summarize_SLR_results <- function(results_list) {
  
  summary_df <- data.frame(
    Penalty = names(results_list),  
    Median_Num_Vars = sapply(results_list, function(x) median(x$nvar_selected)), 
    Median_Misclass_Train = sapply(results_list, function(x) median(x$miscl_train)),
    Median_Misclass_Test = sapply(results_list, function(x) median(x$miscl_test)),
    Median_MSE_Train = sapply(results_list, function(x) median(x$mse_train)),
    Median_MSE_Test = sapply(results_list, function(x) median(x$mse_test)),
    Median_PR_AUC_Train = sapply(results_list, function(x) median(x$pr_auc_train)),
    Median_PR_AUC_Test = sapply(results_list, function(x) median(x$pr_auc_test))
  )
  
  return(summary_df)
}


# ========================================================
# Function: consistently_selected_genes
# Purpose:
#   Identify genes that are consistently selected across bootstrap iterations
#   in sparse logistic regression models for each penalty scheme.
#
# Inputs:
#   results_list - A named list of results from run_sparse_logistic_regression(),
#                  where each element contains variable selection results for a penalty scheme.
#   threshold    - Minimum number of bootstrap iterations (frequency) in which a gene 
#                  must appear to be considered consistently selected (default = 75).
#
# Outputs:
#   A named list where each element corresponds to a penalty scheme and contains:
#     A character vector of gene names selected in at least 'threshold' bootstrap iterations.
# ========================================================
consistently_selected_genes <- function(results_list, threshold = 75) {
  
  # Create a dataframe to store the genes and their frequencies for each penalty
  gene_frequency_list <- list()
  
  for (penalty in names(results_list)) {
    # Obtain the list of genes selected along the bootstraps
    all_selected_genes <- unlist(results_list[[penalty]]$var_selected_names)
    
    # Count the frequency of each gene
    gene_counts <- table(all_selected_genes)
    
    # Filter only the genes that appeared >= threshold in all bootstraps
    consistently_selected_genes <- names(gene_counts[gene_counts >= threshold])
    
    # Store the result in the list
    gene_frequency_list[[penalty]] <- consistently_selected_genes
  }
  
  return(gene_frequency_list)
}


# ========================================================
# Function: prepare_long_expression_df
# Purpose:
#   Prepare a long-format data frame for ggplot visualization.
#   Combines astrocytoma, GBM, and oligodendroglioma samples,
#   subsets by selected genes, log2-transforms, scales, and reshapes to long format.
#
# Inputs:
#   astro        - Data frame/matrix of astrocytoma samples
#   gbm          - Data frame/matrix of GBM samples
#   oligo        - Data frame/matrix of oligodendroglioma samples
#   selected_genes - Vector of selected genes to subset
#   labels       - Vector of class labels for each dataset 
#                  (default: c("Astro.", "Gbm.", "Oligo."))
#
# Outputs:
#   A data frame in long format with columns:
#     Disease    - Cancer type (Astro., Gbm., Oligo.)
#     Gene       - Gene name
#     Expression - Scaled log2 expression
# ========================================================
prepare_long_expression_df <- function(astro, gbm, oligo, selected_genes,
                                       labels = c("Astro.", "Gbm.", "Oligo.")) {
  # 1. Keep only genes present in all datasets
  common_genes <- intersect(selected_genes, Reduce(intersect, list(colnames(astro), colnames(gbm), colnames(oligo))))
  if (length(common_genes) == 0) stop("No common genes found between selected_genes and datasets.")
  
  # 2. Subset the datasets
  astro_x <- astro[, common_genes, drop = FALSE]
  gbm_x   <- gbm[, common_genes, drop = FALSE]
  oligo_x <- oligo[, common_genes, drop = FALSE]
  
  # 3. Combine datasets
  Xdata <- rbind(astro_x, gbm_x, oligo_x)
  
  # 4. Log2-transform and scale
  Xdata_sc <- scale(log2(Xdata + 1))
  Xdata_df <- as.data.frame(Xdata_sc)
  
  # 5. Create Disease label vector
  cancer_type <- c(rep(labels[1], nrow(astro_x)),
                   rep(labels[2], nrow(gbm_x)),
                   rep(labels[3], nrow(oligo_x)))
  
  Xdata_df$Disease <- cancer_type
  
  # 6. Reshape to long format
  long_df <- melt(Xdata_df, id.vars = "Disease", variable.name = "Gene", value.name = "Expression")
  
  return(long_df)
}


# ========================================================
# Function: wilcox_test_one_vs_others
# Purpose:
#   Perform Wilcoxon rank-sum tests for each gene, comparing one chosen group
#   against the combined samples of the other groups.
#
# Inputs:
#   data       - Long-format data frame with columns: "Gene", "Disease", "Expression"
#   solo_group - The group to be compared against all others (default = "Oligo.")
#
# Outputs:
#   A data frame with:
#     Gene         - Gene name
#     p_value      - Raw Wilcoxon p-value
#     adj_p_value  - BH-adjusted p-value
#     significant  - Logical indicator for adjusted p-value < 0.05
# ========================================================
wilcox_test_one_vs_others <- function(data, solo_group = "Oligo.") {
  # Get list of unique genes
  genes <- unique(data[["Gene"]])
  
  # Prepare empty result frame
  results <- data.frame(Gene = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through genes
  for (g in genes) {
    subset_data <- data[data[["Gene"]] == g, ]
    
    group1 <- subset_data[subset_data[["Disease"]] == solo_group, "Expression"]
    group2 <- subset_data[!subset_data[["Disease"]] %in% solo_group, "Expression"]
    
    if (length(group1) > 0 && length(group2) > 0) {
      test <- wilcox.test(group1, group2)
      results <- rbind(results, data.frame(Gene = g, p_value = test$p.value))
    }
  }
  
  # Adjust p-values and add significance flag
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  results$significant <- results$adj_p_value < 0.05
  
  return(results)
}


# ========================================================
# Function: get_network_vector
# Purpose:
#   Compute normalized network centralities (degree, betweenness, harmonic,
#   eigenvector) for a set of genes in an igraph object, plus triangles if provided.
#
# Inputs:
#   graph        - igraph object
#   genes        - character vector of gene names to extract (must match vertex names)
#   triangle_vec - optional named numeric vector of triangle counts keyed by gene; default NULL
#
# Outputs:
#   Named list where each element is a named numeric vector for one gene with:
#     degree, betweenness, harmonic, eigen[, triangles]
#   Notes:
#     - Degree/betweenness/harmonic are computed with normalized = TRUE.
#     - Harmonic NaNs are coerced to 0.
#     - If a gene is absent from the graph, the function returns a NaN vector of length 4 or 5.
#     - If triangle_vec is given but missing for a gene, 'triangles' will be NA.
# ========================================================
get_network_vector <- function(graph, genes, triangle_vec = NULL) {
  deg <- igraph::degree(graph, normalized = TRUE)
  btw <- betweenness(graph, normalized = TRUE)
  hmc <- harmonic_centrality(graph, normalized = TRUE)
  eig <- eigen_centrality(graph)$vector
  
  hmc[is.nan(hmc)] <- 0
  
  # Create named list of centrality vectors
  centrality_list <- lapply(genes, function(g) {
    if (g %in% names(deg)) {
      if (!is.null(triangle_vec)) {
        c(
          degree = deg[g],
          betweenness = btw[g],
          harmonic = hmc[g],
          eigen = eig[g],
          triangles = triangle_vec[g]
        )
      } else {
        # Sem triangles
        c(
          degree = deg[g],
          betweenness = btw[g],
          harmonic = hmc[g],
          eigen = eig[g]
        )
      }
    } else {
      rep(NaN, ifelse(is.null(triangle_vec), 4, 5))
    }
  })
  
  names(centrality_list) <- genes
  return(centrality_list)
}



# ========================================================
# Function: build_network_df_from_three_graphs
# Purpose:
#   Build a data frame containing network measures for all unique genes
#   across three graphs (Astrocytoma, GBM, and Oligodendroglioma).
#   Adds the "triangles" metric only if a triangle_df is provided.
#
# Inputs:
#   graph_astro - igraph object for astrocytoma network
#   graph_gbm   - igraph object for GBM network
#   graph_oligo - igraph object for oligodendroglioma network
#   triangle_df - Optional data frame of triangle counts (genes x conditions)
#
# Outputs:
#   A data frame with:
#     Rows: genes (union of all nodes)
#     Columns: list-columns Astro, GBM, Oligo (each a named vector of measures)
# ========================================================
build_network_df_from_three_graphs <- function(graph_astro, graph_gbm, graph_oligo,
                                               triangle_df = NULL) {
  all_genes <- Reduce(union, list(
    V(graph_astro)$name,
    V(graph_gbm)$name,
    V(graph_oligo)$name
  ))
  
  # Function to retrieve triangle counts for a given disease
  get_tri_vec <- function(disease) {
    if (!is.null(triangle_df) && disease %in% colnames(triangle_df)) {
      vec <- triangle_df[, disease]
      names(vec) <- rownames(triangle_df)
      vec[setdiff(all_genes, names(vec))] <- 0
      vec <- vec[all_genes]  # ensure ordering
      vec[is.na(vec)] <- 0   # convert any remaining NAs to 0
      return(vec)
    } else {
      return(NULL)  # <-- Now properly returns NULL if no triangle data
    }
  }
  
  # Get triangle vectors only if available
  tri_astro <- get_tri_vec("Astro")
  tri_gbm   <- get_tri_vec("GBM")
  tri_oligo <- get_tri_vec("Oligo")
  
  # Get centrality/network measures
  astro_cent <- get_network_vector(graph_astro, all_genes, tri_astro)
  gbm_cent   <- get_network_vector(graph_gbm, all_genes, tri_gbm)
  oligo_cent <- get_network_vector(graph_oligo, all_genes, tri_oligo)
  
  # Build final data frame
  df <- data.frame(
    Astro = I(vector("list", length(all_genes))),
    GBM   = I(vector("list", length(all_genes))),
    Oligo = I(vector("list", length(all_genes))),
    row.names = all_genes
  )
  
  for (i in seq_along(all_genes)) {
    gene <- all_genes[i]
    df$Astro[[i]] <- astro_cent[[gene]]
    df$GBM[[i]]   <- gbm_cent[[gene]]
    df$Oligo[[i]] <- oligo_cent[[gene]]
  }
  
  return(df)
}




# ========================================================
# Function: remove_zero_vectors
# Purpose:
#   Filter out genes whose centrality vectors are entirely zeros for two
#   specified networks (given by 'diseases'). Useful to drop nodes with no
#   connectivity/centrality contribution in both.
#
# Inputs:
#   centrality_df - data frame from build_network_df_from_three_graphs() (or similar)
#                   with list-columns of named numeric vectors per gene.
#   diseases      - character vector of length 2 naming the two list-columns to check
#                   (e.g., c("Astro","GBM") or c("Astro","Oligo")).
#
# Outputs:
#   A filtered data frame excluding rows where both selected vectors are all zeros.
# Notes:
#   - NA/NaN entries are coerced to 0 before the check.
#   - Genes absent from a graph (NaN vectors) will be considered zeros.
#   - Columns in 'diseases' must exist in 'centrality_df'.
# ========================================================

remove_zero_vectors <- function(centrality_df, diseases) {
  is_zero_vector <- function(vec) {
    vals <- as.numeric(unname(vec))  
    vals[is.na(vals) | is.nan(vals)] <- 0
    all(vals == 0)
  }
  
  zero_a <- sapply(centrality_df[[diseases[1]]], is_zero_vector)
  zero_b <- sapply(centrality_df[[diseases[2]]], is_zero_vector)
  
  filtered_df <- centrality_df[!(zero_a & zero_b), ]
  return(filtered_df)
}



# ========================================================
# Function: normalize_network_df
# Purpose:
#   Apply min-max normalization by metric across all genes,
#   for each network (Astro, GBM, Oligo) in a centrality/network data frame.
#
# Inputs:
#   network_df - Data frame produced by build_network_df_from_three_graphs(),
#                where each cell is a named vector of metrics (degree, betweenness, etc.)
#
# Outputs:
#   A normalized data frame of the same structure, with values scaled to [0,1] per metric and per network.
# ========================================================
normalize_network_df <- function(network_df) {
  clean_names <- function(nms) sub("\\..*$", "", nms)  # keep only part before '.'
  metrics <- clean_names(names(network_df[[1]][[1]]))
  networks <- colnames(network_df)
  
  normalized_df <- network_df
  
  for (net in networks) {
    mat <- do.call(rbind, network_df[[net]])
    
    # Normalize each metric
    for (m in metrics) {
      col_idx <- grep(paste0("^", m, "\\."), colnames(mat))  # columns matching metric
      vals <- mat[, col_idx]
      rng <- range(vals, na.rm = TRUE)
      if (diff(rng) == 0) {
        mat[, col_idx] <- 0
      } else {
        mat[, col_idx] <- (vals - rng[1]) / diff(rng)
      }
    }
    
    # Rebuild proper numeric named vectors
    new_list <- lapply(seq_len(nrow(mat)), function(i) {
      vals <- as.numeric(mat[i, ])
      names(vals) <- metrics
      vals
    })
    names(new_list) <- rownames(mat)
    normalized_df[[net]] <- new_list
  }
  
  return(normalized_df)
}


# ========================================================
# Function: compute_distance
# Purpose:
#   Compute a distance between two centrality vectors using:
#     - "euclidean", "manhattan", "bray-curtis", "canberra", or "ratio".
#   Missing entries (NA/NaN) are ignored by comparing only common valid positions.
#
# Inputs:
#   v1, v2 - numeric vectors of equal-length centrality values
#   method - one of: "euclidean", "bray-curtis", "canberra", "manhattan", "ratio"
#
# Outputs:
#   A numeric scalar with the distance, or NA if no valid positions.
# Notes:
#   - "bray-curtis" assumes non-negative inputs; if negatives are possible,
#     use denom = sum(abs(x) + abs(y)).
#   - "canberra" skips terms with zero denominator.
#   - "ratio" is scale-invariant per component; non-positive pairs are penalized as 0.
#   - "cosine" and "log-ratio" are mentioned in earlier drafts but are NOT implemented here.
#   - Infinite values are not filtered; if they can occur, pre-clean v1/v2 upstream.
# ========================================================
compute_distance <- function(v1, v2, method = "euclidean") {
  # Identify valid (non-NA, non-NaN) positions in both vectors
  valid <- !(is.na(v1) | is.nan(v1) | is.na(v2) | is.nan(v2))
  
  # If no valid positions exist, return NA
  if (sum(valid) == 0) return(NA)
  
  # Filter vectors to valid entries only
  x <- v1[valid]
  y <- v2[valid]
  
  if (method == "euclidean") {
    return(sqrt(sum((x - y)^2)))
    
  } else if (method == "bray-curtis") {
    num <- sum(abs(x - y))
    denom <- sum(x + y)
    if (denom == 0) return(0) #if both are isolated, distance is 0
    return(num / denom) # distancia
    
  } else if (method == "canberra") {
    denom <- abs(x) + abs(y)
    denom[denom == 0] <- NA  # avoid division by zero
    return(sum(abs(x - y) / denom, na.rm = TRUE)) #distancia proxima de zero- vetores parecidos
    
  } else if (method == "manhattan") {
    return(sum(abs(x - y)))
    
  } else if (method == "ratio") {
    ratios <- numeric(length(x))
    for (i in seq_along(x)) {
      xi <- x[i]
      yi <- y[i]
      if (xi <= 0 || yi <= 0) {
        ratios[i] <- 0  # penaliza semelhanC'a
      } else {
        ratios[i] <- min(xi / yi, yi / xi)
      }
    }
    return(1 - mean(ratios))
    
  } else {
    stop("Unsupported distance method: ", method)
  }
}

# ========================================================
# Function: compare_centrality_vectors
# Purpose:
#   Compute pairwise distances between the centrality vectors of the same gene
#   across three networks (Astrocytoma, GBM, Oligodendroglioma), using a chosen metric.
#
# Inputs:
#   centrality_df - data frame from build_network_df_from_three_graphs(), where each
#                   cell is a named numeric vector of centrality metrics for that gene/network.
#   method        - distance metric passed to compute_distance()
#                   (e.g., "euclidean", "manhattan", "bray-curtis", "canberra", "ratio").
#
# Outputs:
#   A data frame with:
#     Rows: genes (matching rownames(centrality_df))
#     Cols: Astro_vs_GBM, Astro_vs_Oligo, GBM_vs_Oligo (pairwise distances per gene)
# Notes:
#   - If a pair has no common valid positions (all NA/NaN), the distance is NA.
#   - Requires columns named exactly "Astro", "GBM", "Oligo" in centrality_df.
# ========================================================
compare_centrality_vectors <- function(centrality_df, method = "euclidean") {
  genes <- rownames(centrality_df)
  
  # Initialize data frame to store pairwise distances
  distances <- data.frame(
    Astro_vs_GBM   = rep(NA_real_, length(genes)),
    Astro_vs_Oligo = rep(NA_real_, length(genes)),
    GBM_vs_Oligo   = rep(NA_real_, length(genes))
  )
  rownames(distances) <- genes
  
  # Loop through each gene and compute distances between its centrality vectors
  for (i in seq_along(genes)) {
    v_astro <- centrality_df$Astro[[i]]
    v_gbm   <- centrality_df$GBM[[i]]
    v_oligo <- centrality_df$Oligo[[i]]
    
    distances$Astro_vs_GBM[i]   <- compute_distance(v_astro, v_gbm, method)
    distances$Astro_vs_Oligo[i] <- compute_distance(v_astro, v_oligo, method)
    distances$GBM_vs_Oligo[i]   <- compute_distance(v_gbm, v_oligo, method)
  }
  
  return(distances)
}


# ========================================================
# Function: get_normalized_weights
# Purpose:
#   Min–max normalize one distance column in a centrality-distance data frame.
#   Optionally choose the column by name or index. Returned weights remain
#   aligned with the original row order (one weight per gene).
#
# Inputs:
#   distances_df - data frame from compare_centrality_vectors(), with pairwise
#                  distances per gene (e.g., Astro_vs_GBM, Astro_vs_Oligo, GBM_vs_Oligo).
#   column       - column name or 1-based index to normalize (default = 1)
#
# Outputs:
#   A data frame with a single column named after the chosen distance column,
#   containing min–max normalized weights in [0, 1] (NAs preserved).
# Notes:
#   - If all finite values are equal, all weights are set to 0.
#   - Exact zeros are replaced by the smallest non-zero weight (if any),
#     to avoid zero-penalty features downstream.
#   - Row names are preserved from distances_df.
# ========================================================

get_normalized_weights <- function(distances_df, column = 1) {
  # Resolve column
  if (is.character(column)) {
    if (!column %in% names(distances_df)) {
      stop(sprintf("Column '%s' not found in distances_df.", column))
    }
    values <- distances_df[[column]]
    col_name <- column
  } else if (is.numeric(column)) {
    if (column < 1 || column > ncol(distances_df)) {
      stop("Column index out of bounds.")
    }
    values <- distances_df[[column]]
    col_name <- names(distances_df)[column]
  } else {
    stop("Argument 'column' must be a column name or index.")
  }
  
  # Min–max normalization [0,1]
  rng <- range(values, na.rm = TRUE)
  if (diff(rng) == 0) {
    weights <- rep(0, length(values))   # all equal -> 0
  } else {
    weights <- (values - rng[1]) / diff(rng)
  }
  
  # Replace exact zeros by the smallest non-zero weight (if any)
  min_nonzero <- suppressWarnings(min(weights[weights > 0], na.rm = TRUE))
  if (is.finite(min_nonzero)) {
    weights[weights == 0 & !is.na(weights)] <- min_nonzero
  }
  
  # Build a one-column data frame with normalized weights
  out_df <- data.frame(weights)
  colnames(out_df) <- col_name
  rownames(out_df) <- rownames(distances_df)
  
  return(out_df)
}





# ========================================================
# Function: add_custom_scores
# Purpose:
#   From a 3-column distance data frame, pick one “primary” column and
#   compute three custom scores that favor high primary vs. low others:
#     • linear:  x − λ(y + z)
#     • ratio :  x / ( (y + z)^α + ε )
#     • log   :  log1p(x) − λ·log1p(y + z)
#   Each score is then min–max scaled to [0,1] and exact zeros are lifted
#   to the smallest non-zero value to avoid zero weights downstream.
#
# Inputs:
#   df      - data frame with exactly 3 numeric distance columns
#             (e.g., Astro_vs_GBM, Astro_vs_Oligo, GBM_vs_Oligo)
#   primary - column name or 1-based index to treat as x (default = 1)
#   lambda  - non-primary penalty weight (default = 0.5)
#   alpha   - exponent for the ratio score (default = 1)
#
# Outputs:
#   The same data frame with three added columns:
#     custom_score_linear_<primary>,
#     custom_score_ratio_<primary>,
#     custom_score_log_<primary>.
# Notes:
#   - Row order and row names are preserved.
#   - Min–max is applied separately per score.
#   - ε = 1e−6 guards against division by zero in the ratio score.
# ========================================================

add_custom_scores <- function(df, primary = 1, lambda = 0.5, alpha = 1) {
  # Preserve row names to restore later
  row_names <- rownames(df)
  
  # Helper: min–max scale to [0,1]
  normalize_minmax <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  }
  
  # Resolve primary column (x)
  if (is.character(primary)) {
    if (!primary %in% colnames(df)) {
      stop(sprintf("Column '%s' not found in dataframe.", primary))
    }
    x <- df[[primary]]
    col_x <- primary
  } else if (is.numeric(primary)) {
    if (primary < 1 || primary > ncol(df)) stop("Column index out of bounds.")
    x <- df[[primary]]
    col_x <- colnames(df)[primary]
  } else {
    stop("Argument 'primary' must be a column name or index.")
  }
  
  # The remaining two columns become y and z (order given by setdiff)
  other_cols <- setdiff(colnames(df), col_x)
  y <- df[[other_cols[1]]]
  z <- df[[other_cols[2]]]
  
  # Small constant for numerical safety in ratios
  eps <- 1e-6
  
  # ----------- Custom scores ----------- #
  s_linear <- x - lambda * (y + z)
  s_ratio  <- x / (((y + z)^alpha) + eps)
  s_log    <- log1p(x) - lambda * log1p(y + z)
  
  # ----------- Scale to [0,1] and avoid exact zeros ----------- #
  fix_zeros <- function(v) {
    v_norm <- normalize_minmax(v)
    min_nz <- suppressWarnings(min(v_norm[v_norm > 0], na.rm = TRUE))
    if (is.finite(min_nz)) {
      v_norm[v_norm == 0 & !is.na(v_norm)] <- min_nz
    }
    v_norm
  }
  
  s_linear <- fix_zeros(s_linear)
  s_ratio  <- fix_zeros(s_ratio)
  s_log    <- fix_zeros(s_log)
  
  # Attach score columns, names carry the chosen primary
  df[[paste0("custom_score_linear_", col_x)]] <- s_linear
  df[[paste0("custom_score_ratio_",  col_x)]] <- s_ratio
  df[[paste0("custom_score_log_",    col_x)]] <- s_log
  
  # Restore row names and return
  rownames(df) <- row_names
  df
}

# ========================================================
# Function: tuning_custom_parameters
# Purpose:
#   Generate a set of candidate penalty vectors by scanning λ (for linear/log)
#   and α (for ratio) around a given primary distance column. Useful for
#   hyper-parameter sweeps when building penalty.factor variants.
#
# Inputs:
#   distance_df  - 3-column distance data frame (as above)
#   primary      - primary column (name or index) passed to add_custom_scores()
#   distance_name- short tag used to name the returned vectors
#   lambdas      - numeric vector of λ values to scan (default c(0.20,0.30,0.40,0.5))
#   alphas       - numeric vector of α values to scan (default c(0.50,0.60,0.70,0.8))
#
# Outputs:
#   A named list of numeric vectors (one per setting), with keys like:
#     "<distance_name>_linear_l0.30", "<distance_name>_log_l0.30",
#     "<distance_name>_ratio_a0.70"
# Notes:
#   - For λ scans, α is fixed at alphas[1]; for α scans, λ is fixed at lambdas[1].
#   - Vectors are taken directly from the normalized score columns created by add_custom_scores().
# ========================================================

tuning_custom_parameters <- function(distance_df, primary, distance_name,
                                     lambdas = c(0.20, 0.30, 0.40, 0.5),
                                     alphas  = c(0.50, 0.60, 0.70, 0.8)) {
  pens <- list()
  
  # --- scan lambdas: produce linear/log once per lambda ---
  for (l in lambdas) {
    tmp <- add_custom_scores(distance_df, primary = primary, lambda = l, alpha = alphas[1])
    key_l <- sprintf("%s_linear_l%.2f", distance_name, l)
    key_g <- sprintf("%s_log_l%.2f",    distance_name, l)
    pens[[key_l]] <- tmp[[paste0("custom_score_linear_", colnames(distance_df)[primary])]]
    pens[[key_g]] <- tmp[[paste0("custom_score_log_",    colnames(distance_df)[primary])]]
  }
  
  # --- scan alphas: produce ratio once per alpha ---
  for (a in alphas) {
    tmp <- add_custom_scores(distance_df, primary = primary, lambda = lambdas[1], alpha = a)
    key_r <- sprintf("%s_ratio_a%.2f", distance_name, a)
    pens[[key_r]] <- tmp[[paste0("custom_score_ratio_", colnames(distance_df)[primary])]]
  }
  
  pens
}


# ========================================================
# Function: run_parameter_tuning_SLR
# Purpose:
#   Build multiple penalty vectors by scanning λ (linear/log) and α (ratio),
#   run bootstrapped sparse logistic regression for each, and summarize results.
#
# Inputs:
#   distance_df  - 3-column distance data frame (pairwise network distances)
#   primary      - primary column (name or index) to anchor custom scores
#   distance_name- short tag used to label penalty vectors
#   Xdata, Ydata - design matrix (n x p) and binary response (length n)
#   lambdas      - λ values for linear/log scores (default c(0.20,0.30,0.40,0.5))
#   alphas       - α values for ratio score (default c(0.50,0.60,0.70,0.8))
#   times_boot   - number of bootstrap iterations (default 100)
#
# Outputs:
#   A list with:
#     penalties - named list of penalty vectors (from tuning_custom_parameters)
#     slr       - raw results from run_SLR_penalized()
#     summary   - data frame from summarize_SLR_results(slr)
# Notes:
#   - Requires tuning_custom_parameters(), run_SLR_penalized(), and summarize_SLR_results().
#   - Penalty names encode the equation type (linear/log/ratio) and the scanned hyperparameter.
# ========================================================
run_parameter_tuning_SLR <- function(distance_df, primary, distance_name,
                                     Xdata, Ydata,
                                     lambdas = c(0.20,0.30,0.40, 0.5),
                                     alphas  = c(0.50,0.60,0.70, 0.8),
                                     times_boot = 100) {
  pens <- tuning_custom_parameters(distance_df, primary, distance_name, lambdas, alphas)
  
  # Prepare penalty set and run bootstrapped SLR for each penalty vector
  # Keys are preserved to allow grouping/summarization by equation later.
  network_penalty <- pens
  network_SLR     <- run_SLR_penalized(Xdata, Ydata, network_penalty, times_boot = times_boot)
  network_sum     <- summarize_SLR_results(network_SLR)
  
  # Pack artifacts for downstream selection/reporting
  list(penalties = pens, slr = network_SLR, summary = network_sum)
}

# ========================================================
# Function: pick_best_model_per_equation
# Purpose:
#   From a summary table of SLR runs, pick the best row per equation family
#   (linear/log/ratio) according to a chosen metric.
#
# Inputs:
#   SLR_summary - data frame with at least columns 'Penalty' and the target metric
#   metric      - column name to maximize (default "Median_PR_AUC_Test")
#
# Outputs:
#   A data frame with up to three rows (one per equation), containing the
#   best-scoring configuration for each family. Returns an empty df if none.
# Notes:
#   - Equation family is inferred from 'Penalty' via "_linear_", "_log_", "_ratio_".
#   - If ties exist, the first max row is returned (base R behavior of which.max).
#   - Rows with all-NA metric are skipped for that family.
# ========================================================
pick_best_model_per_equation <- function(SLR_summary,
                                         metric = "Median_PR_AUC_Test") {
  stopifnot(is.data.frame(SLR_summary))
  if (!"Penalty" %in% names(SLR_summary)) {
    stop("SLR_summary must contain a 'Penalty' column.")
  }
  if (!metric %in% names(SLR_summary)) {
    stop("Requested metric '", metric, "' not found in SLR_summary.")
  }
  
  # Robustly derive equation tag from the penalty name
  pen <- SLR_summary$Penalty
  eq <- ifelse(grepl("_linear_", pen), "linear",
               ifelse(grepl("_log_",    pen), "log",
                      ifelse(grepl("_ratio_",  pen), "ratio", NA_character_)))
  SLR_summary$equation <- eq
  
  # Split per equation and pick the row with the highest metric
  split_tbl <- split(SLR_summary, SLR_summary$equation, drop = TRUE)
  best_rows <- lapply(split_tbl, function(dd) {
    if (nrow(dd) == 0 || all(is.na(dd[[metric]]))) return(NULL)
    dd[which.max(dd[[metric]]), , drop = FALSE]
  })
  best <- do.call(rbind, best_rows)
  
  if (is.null(best) || nrow(best) == 0) {
    return(SLR_summary[0, , drop = FALSE])  # empty with same columns
  }
  
  # Remove helper column; return compact table
  best$equation <- NULL
  rownames(best) <- NULL
  best
}


# ========================================================
# Function: collect_SLR_details
# Purpose:
#   For each requested equation family (linear/log/ratio), retrieve the
#   corresponding raw SLR result object from SLR_results$slr using the
#   best penalty key taken from 'best_models'.
#
# Inputs:
#   best_models - data.frame with a 'Penalty' column and optionally 'equation';
#                 each row represents a selected (best) configuration.
#   SLR_results - list returned by run_parameter_tuning_SLR(), must contain $slr
#                 (a named list keyed by penalty names).
#   equations   - character vector of equation tags to collect (default: c("linear","log","ratio"))
#
# Outputs:
#   A named list with one element per equation; each element is the SLR result
#   object corresponding to the best penalty for that family. Missing entries
#   remain NULL and trigger a warning.
# Notes:
#   - If 'equation' column is absent in best_models, it infers the family from
#     the 'Penalty' name using the pattern "_<eq>_" (case-insensitive).
#   - Emits warnings when a best row or a corresponding SLR object is not found.
# ========================================================
collect_SLR_details <- function(best_models, SLR_results,
                                equations = c("linear", "log", "ratio")) {
  # --- sanity checks ---
  if (!is.data.frame(best_models) || !"Penalty" %in% names(best_models)) {
    stop("best_models must be a data.frame containing a 'Penalty' column.")
  }
  if (!is.list(SLR_results) || is.null(SLR_results$slr)) {
    stop("SLR_results must be a list that contains an 'slr' list.")
  }
  
  slr_pool <- SLR_results$slr
  out <- setNames(vector("list", length(equations)), equations)
  
  # helper: pick the row for a given equation (uses 'equation' if present; otherwise regex on Penalty)
  pick_row_for_eq <- function(eq) {
    if ("equation" %in% names(best_models)) {
      row <- best_models[tolower(best_models$equation) == tolower(eq), , drop = FALSE]
    } else {
      row <- best_models[grepl(paste0("_", eq, "_"), best_models$Penalty, ignore.case = TRUE),
                         , drop = FALSE]
    }
    row
  }
  
  # Iterate families and fetch the matching SLR object by penalty key
  for (eq in equations) {
    row_eq <- pick_row_for_eq(eq)
    if (nrow(row_eq) == 0) {
      warning("No best-model row found for equation '", eq, "'.")
      next
    }
    key <- row_eq$Penalty[1]
    if (is.null(slr_pool[[key]])) {
      warning("SLR object '", key, "' not found in SLR_results$slr.")
      next
    }
    out[[eq]] <- slr_pool[[key]]
  }
  
  out
}




# ========================================================
# Function: add_best_penalties_to_df
# Purpose:
#   Attach the best penalty vectors (one per equation family) to a distance
#   data frame, aligning by gene row names. Each added column is named after
#   the equation tag (e.g., "linear","log","ratio").
#
# Inputs:
#   SLR_results - list from run_parameter_tuning_SLR(); must contain $penalties
#                 (a named list keyed by penalty names).
#   best_models - data.frame with a 'Penalty' column listing selected winners.
#   distance_df - data frame whose rownames are gene IDs used to align penalties.
#   equations   - character vector of equation tags to add (default c("linear","log","ratio")).
#
# Outputs:
#   The same distance_df, augmented with one numeric column per equation
#   containing the aligned penalty weights (NAs may appear for missing genes).
# Notes:
#   - The best penalty key for each equation is detected via pattern "_<eq>_".
#   - If multiple rows match in best_models, the first is used.
#   - Alignment prefers names: penalties are reordered by rownames(distance_df).
#     If no names but lengths match, positional alignment is used.
#   - If a target column already exists in distance_df, it will be overwritten.
#   - Missing entries after alignment produce NA values.
# ========================================================
add_best_penalties_to_df <- function(SLR_results, best_models, distance_df,
                                     equations = c("linear","log","ratio")) {
  # --- basic checks ---
  if (is.null(SLR_results$penalties) || !is.list(SLR_results$penalties)) {
    stop("SLR_results must have a named list `penalties`.")
  }
  if (!("Penalty" %in% names(best_models))) {
    stop("best_models must have a column named `Penalty`.")
  }
  if (is.null(rownames(distance_df))) {
    stop("distance_df must have rownames (gene ids) so penalties can be aligned.")
  }
  
  # Helper: fetch the method key for each equation from best_models
  pick_key_for_eq <- function(eq) {
    rows <- best_models[grepl(paste0("_", eq, "_"), best_models$Penalty), , drop = FALSE]
    if (nrow(rows) == 0) return(NA_character_)
    rows$Penalty[1]
  }
  
  # For each equation, pull the penalty vector and align to distance_df rows
  for (eq in equations) {
    key <- pick_key_for_eq(eq)
    if (is.na(key)) {
      warning("No best model found for equation '", eq, "'. Skipping.")
      next
    }
    pen <- SLR_results$penalties[[key]]
    if (is.null(pen)) {
      warning("Penalty vector '", key, "' not found in SLR_results$penalties. Skipping.")
      next
    }
    
    # Align by names when possible; otherwise use positional if lengths match
    if (!is.null(names(pen)) && length(intersect(names(pen), rownames(distance_df))) > 0) {
      # Reorder by the dataframe's rownames; may introduce NAs if names missing
      pen_aligned <- pen[rownames(distance_df)]
    } else if (length(pen) == nrow(distance_df)) {
      pen_aligned <- pen
      names(pen_aligned) <- rownames(distance_df)
    } else {
      stop(
        sprintf("Cannot align penalties for '%s': length %d vs nrow(distance_df) %d, and no names.",
                key, length(pen), nrow(distance_df))
      )
    }
    
    distance_df[[eq]] <- as.numeric(pen_aligned)
  }
  
  distance_df
}



# --------------------------------------------------------
# extract_genes_from_subnetworks
# Purpose: Convert a community 'membership' vector into a list of genes per module.
# Input : components (e.g., igraph::cluster_* result with $membership)
# Output: named list {module_id -> character vector of genes}
# --------------------------------------------------------
extract_genes_from_subnetworks <- function(components) {
  modules_list <- split(names(components$membership), components$membership)
  return(modules_list)
}

# --------------------------------------------------------
# filter_subnetworks
# Purpose: Keep only modules that contain at least one of 'selected_genes'.
# Inputs : subnetworks (list of character vectors), selected_genes (character)
# Output : filtered list of modules
# --------------------------------------------------------
filter_subnetworks <- function(subnetworks, selected_genes) {
  subnetworks[sapply(subnetworks, function(module) {
    any(module %in% selected_genes)   # TRUE if the module has ≥1 selected gene
  })]
}

# --------------------------------------------------------
# extract_module_graphs
# Purpose: Build induced subgraphs per module and return basic artifacts.
# Inputs : graph (igraph), module_list (list of gene vectors)
# Output : list of {genes, graph, adj_matrix} per module (named "1","2",...)
# --------------------------------------------------------
extract_module_graphs <- function(graph, module_list) {
  out <- lapply(module_list, function(genes) {
    genes <- intersect(genes, igraph::V(graph)$name)  # keep only existing vertices
    subg  <- igraph::induced_subgraph(graph, vids = genes)
    list(
      genes      = genes,
      graph      = subg,
      adj_matrix = igraph::as_adjacency_matrix(subg, sparse = FALSE)
    )
  })
  names(out) <- seq_along(out)  # simple sequential labels
  out
}

# --------------------------------------------------------
# find_modules_for_genes
# Purpose: For each gene, fetch the first module containing it in Astro/GBM/Oligo.
# Inputs : genes (character), *_modules (lists from extract_module_graphs)
# Output : named list per gene with $Astro/$GBM/$Oligo (or NULL if absent)
# --------------------------------------------------------
find_modules_for_genes <- function(genes, astro_modules, gbm_modules, oligo_modules) {
  pick <- function(gene, mods) {
    idx <- which(vapply(mods, function(m) gene %in% m$genes, logical(1)))
    if (length(idx) == 0) return(NULL)
    m <- mods[[idx[1]]]
    m$module_index <- idx[1]
    m
  }
  out <- lapply(genes, function(g) {
    list(
      Astro = pick(g, astro_modules),
      GBM   = pick(g, gbm_modules),
      Oligo = pick(g, oligo_modules)
    )
  })
  names(out) <- genes
  out
}

# --------------------------------------------------------
# get_module_gene_lists
# Purpose: Reduce a module object list to a {Mi -> genes} named list.
# Input : mod_list (like output from extract_module_graphs)
# Output: list named M1, M2, ... with gene symbol vectors
# --------------------------------------------------------
get_module_gene_lists <- function(mod_list) {
  # return list: names = M1, M2, ..., values = gene vectors
  if (length(mod_list) == 0) return(list())
  out <- lapply(seq_along(mod_list), function(i) mod_list[[i]]$genes)
  names(out) <- paste0("M", seq_along(out))
  out
}


# --------------------------------------------------------
# save_genes_txt
# Purpose: Save a character vector of gene symbols to a single-line, tab-separated TXT.
# Input : genes (character or coercible), file_path (output path)
# Output: writes a file; returns invisibly via message with absolute path
# Notes :
#   - Uses a single row (t()) so the output is one line with tabs.
#   - Does not create parent folders; ensure dirname(file_path) exists.
# --------------------------------------------------------
save_genes_txt <- function(genes, file_path = "genes_list.txt") {
  # Ensure it's character
  genes <- as.character(genes)
  
  # Write as one line, tab-separated, without quotes
  write.table(
    t(genes), 
    file = file_path, 
    quote = FALSE, 
    sep = "\t", 
    row.names = FALSE, 
    col.names = FALSE
  )
  
  message("File saved at: ", normalizePath(file_path))
}

# --------------------------------------------------------
# symbols_to_entrez_robust
# Purpose: Clean/standardize gene symbols and map them to ENTREZ IDs.
# Inputs :
#   x       - character vector of gene symbols
#   species - passed to HGNChelper::checkGeneSymbols() (e.g., "human")
# Outputs:
#   list(mapped = data.frame(SYMBOL, ENTREZID),
#        not_mapped = character vector of symbols without mapping)
# Notes  :
#   - Normalization pipeline:
#       • trim/unique → HGNChelper::checkGeneSymbols() for up-to-date symbols
#       • fallback to original when no suggestion
#       • uppercase; replace '.' with '-' (rare cases like "ZNF585B.1")
#       • strip whitespace
#   - Mapping via clusterProfiler::bitr(..., OrgDb = org.Hs.eg.db).
#   - Requires packages: HGNChelper, clusterProfiler, org.Hs.eg.db.
# --------------------------------------------------------
symbols_to_entrez_robust <- function(x, species = "human") {
  stopifnot(is.character(x))
  x <- unique(trimws(x))
  x <- x[nzchar(x)]
  
  # Normalize: fix outdated symbols and common issues
  ch <- HGNChelper::checkGeneSymbols(x, species = species)
  fixed <- ch$Suggested                     # corrections when available
  fixed[is.na(fixed)] <- ch$x[is.na(fixed)] # keep original if no suggestion
  fixed <- toupper(fixed)
  fixed <- gsub("\\.", "-", fixed)          # e.g., "ZNF585B.1" -> "ZNF585B-1" (rare)
  fixed <- gsub("\\s+", "", fixed)
  
  # Map via org.Hs.eg.db
  map <- suppressMessages(
    clusterProfiler::bitr(fixed,
                          fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
  )
  map <- unique(map)
  # Drop empty symbols and missing ENTREZIDs
  map <- map[!is.na(map$SYMBOL) & !is.na(map$ENTREZID), , drop = FALSE]
  
  not_mapped <- setdiff(unique(fixed), unique(map$SYMBOL))
  
  list(
    mapped = map,               # data.frame with SYMBOL and ENTREZID
    not_mapped = not_mapped     # character vector of symbols not mapped
  )
}


# ========================================================
# Function: enrich_go_by_disease
# Purpose:
#   Run GO ORA per disease-specific gene set with a shared background.
#
# Inputs:
#   disease_gene_sets - named list {disease_label -> character vector of SYMBOLs}
#   universe_symbols  - optional character vector of SYMBOLs for the background;
#                       if NULL/empty, use full OrgDb ENTREZ universe
#   ont               - GO ontology: "BP", "MF", or "CC" (default "BP")
#   min_size          - minimum mapped gene set size (ENTREZ) to test (default 5)
#   p_adj_cut         - FDR threshold on p.adjust (default 0.05)
#   verbose           - print progress messages (default TRUE)
#
# Output:
#   Tibble/data.frame of enriched GO terms across diseases (readable SYMBOLs),
#   annotated with Disease, N_input, N_mapped. Returns empty tibble if none.
# Notes:
#   - SYMBOL→ENTREZ mapping via symbols_to_entrez_robust(); universe is in ENTREZ.
#   - Filtering: intersect with universe BEFORE size check; FDR applied after enrichGO.
#   - readable=TRUE converts core fields back to SYMBOLs for reporting.
# ========================================================
enrich_go_by_disease <- function(disease_gene_sets,
                                 universe_symbols = NULL,
                                 ont = "BP",
                                 min_size = 5,
                                 p_adj_cut = 0.05,
                                 verbose = TRUE) {
  # --- 0) Build universe in ENTREZ space (optional) ---
  if (is.null(universe_symbols) || length(universe_symbols) == 0) {
    # Use full OrgDb as background
    universe_entrez <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")
    if (verbose) {
      message(sprintf("Universe: full OrgDb — %d ENTREZ IDs.", length(universe_entrez)))
    }
  } else {
    # Ensure character (avoid factors)
    universe_symbols <- as.character(universe_symbols)
    conv_u <- symbols_to_entrez_robust(unique(universe_symbols))
    universe_entrez <- unique(conv_u$mapped$ENTREZID)
    if (verbose) {
      message(sprintf("Universe: %d symbols -> %d ENTREZ (not mapped: %d).",
                      length(unique(universe_symbols)),
                      length(universe_entrez),
                      length(conv_u$not_mapped)))
    }
  }
  
  # --- 1) Iterate diseases ---
  out_list <- lapply(names(disease_gene_sets), function(d) {
    # Coerce to character, drop duplicates
    g_sym <- unique(as.character(disease_gene_sets[[d]]))
    conv  <- symbols_to_entrez_robust(g_sym)
    g_ent <- unique(conv$mapped$ENTREZID)
    
    if (verbose) {
      message(sprintf("[%s] List: %d symbols -> %d ENTREZ (not mapped: %d).",
                      d, length(g_sym), length(g_ent), length(conv$not_mapped)))
    }
    
    # Intersect with the universe BEFORE size check
    g_ent <- intersect(g_ent, universe_entrez)
    if (verbose) {
      message(sprintf("[%s] After intersect with universe: %d ENTREZ.", d, length(g_ent)))
    }
    
    if (length(g_ent) < min_size) {
      if (verbose) message(sprintf("[%s] < min_size (%d) — skipped.", d, min_size))
      return(NULL)
    }
    
    # --- 2) GO enrichment ---
    ego <- clusterProfiler::enrichGO(
      gene          = g_ent,
      universe      = universe_entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      readable      = TRUE
    )
    
    if (is.null(ego)) {
      if (verbose) message(sprintf("[%s] enrichGO returned NULL.", d))
      return(NULL)
    }
    
    df <- as.data.frame(ego)
    if (!nrow(df)) {
      if (verbose) message(sprintf("[%s] enrichGO returned 0 rows.", d))
      return(NULL)
    }
    
    # --- 3) FDR filter ---
    df <- subset(df, p.adjust <= p_adj_cut)
    if (!nrow(df)) {
      if (verbose) message(sprintf("[%s] No terms with FDR <= %.3f.", d, p_adj_cut))
      return(NULL)
    }
    
    # --- 4) Extra annotations ---
    df$Disease  <- d
    df$N_input  <- length(g_sym)
    df$N_mapped <- length(g_ent)
    df
  })
  
  # --- 5) Bind or return empty tibble ---
  res <- dplyr::bind_rows(out_list)
  if (is.null(res) || !nrow(res)) {
    if (verbose) message("No enriched GO terms across diseases (after filtering).")
    return(tibble::tibble())
  }
  res
}


# ========================================================
# Function: enrich_kegg_by_disease
# Purpose:
#   Run KEGG ORA per disease-specific gene set, optionally with a custom
#   ENTREZ universe (if supported by your clusterProfiler::enrichKEGG).
#
# Inputs:
#   disease_gene_sets - named list {disease_label -> character vector of SYMBOLs}
#   universe_symbols  - optional SYMBOL vector for background; if NULL/empty,
#                       use KEGG's default organism background
#   organism          - KEGG organism code (default "hsa")
#   min_size          - minimum mapped gene set size (ENTREZ) to test (default 5)
#   p_adj_cut         - FDR threshold on p.adjust (default 0.05)
#   verbose           - print progress messages (default TRUE)
#
# Output:
#   Tibble/data.frame of enriched KEGG pathways across diseases (readable SYMBOLs),
#   annotated with Disease, N_input, N_mapped. Returns empty tibble if none.
# Notes:
#   - SYMBOL→ENTREZ via symbols_to_entrez_robust(); readable via setReadable().
#   - Some clusterProfiler versions don’t support 'universe' in enrichKEGG; this is detected.
#   - Intersects with universe (if provided) BEFORE size check; FDR applied after enrichment.
# ========================================================
enrich_kegg_by_disease <- function(disease_gene_sets,
                                   universe_symbols = NULL,
                                   organism = "hsa",
                                   min_size = 5,
                                   p_adj_cut = 0.05,
                                   verbose = TRUE) {
  # --- 0) Optional universe (ENTREZ) ---
  use_universe <- !is.null(universe_symbols) && length(universe_symbols) > 0
  if (use_universe) {
    conv_u <- symbols_to_entrez_robust(unique(universe_symbols))
    universe_entrez <- unique(conv_u$mapped$ENTREZID)
    if (verbose) {
      message(sprintf("Universe: %d symbols -> %d ENTREZ (not mapped: %d).",
                      length(unique(universe_symbols)),
                      length(universe_entrez),
                      length(conv_u$not_mapped)))
    }
  } else {
    universe_entrez <- NULL
    if (verbose) {
      message("Universe: using KEGG's default background for the selected organism.")
    }
  }
  
  # Does this clusterProfiler version support 'universe' in enrichKEGG?
  supports_universe <- "universe" %in% names(formals(clusterProfiler::enrichKEGG))
  if (use_universe && !supports_universe && verbose) {
    message("Note: your enrichKEGG() may not support 'universe'; falling back to KEGG default background.")
  }
  
  # --- 1) Iterate diseases ---
  out_list <- lapply(names(disease_gene_sets), function(d) {
    g_sym <- unique(disease_gene_sets[[d]])
    conv  <- symbols_to_entrez_robust(g_sym)
    g_ent <- unique(conv$mapped$ENTREZID)
    
    if (verbose) {
      message(sprintf("[%s] List: %d symbols -> %d ENTREZ (not mapped: %d).",
                      d, length(g_sym), length(g_ent), length(conv$not_mapped)))
    }
    
    # Intersect with universe if provided
    if (use_universe) {
      g_ent <- intersect(g_ent, universe_entrez)
      if (verbose) message(sprintf("[%s] After intersect with universe: %d ENTREZ.", d, length(g_ent)))
    }
    
    if (length(g_ent) < min_size) {
      if (verbose) message(sprintf("[%s] < min_size (%d) — skipped.", d, min_size))
      return(NULL)
    }
    
    # --- 2) KEGG enrichment ---
    ek <- try({
      if (use_universe && supports_universe) {
        clusterProfiler::enrichKEGG(
          gene          = g_ent,
          universe      = universe_entrez,   # only if supported
          organism      = organism,
          pvalueCutoff  = 1,
          pAdjustMethod = "BH",
          qvalueCutoff  = 1
        )
      } else {
        clusterProfiler::enrichKEGG(
          gene          = g_ent,
          organism      = organism,
          pvalueCutoff  = 1,
          pAdjustMethod = "BH",
          qvalueCutoff  = 1
        )
      }
    }, silent = TRUE)
    
    if (inherits(ek, "try-error") || is.null(ek)) {
      if (verbose) message(sprintf("[%s] enrichKEGG returned NULL/error.", d))
      return(NULL)
    }
    
    df_raw <- as.data.frame(ek)
    if (!nrow(df_raw)) {
      if (verbose) message(sprintf("[%s] enrichKEGG returned 0 rows.", d))
      return(NULL)
    }
    
    # --- 3) Make gene IDs readable (ENTREZ -> SYMBOL)
    ek_read <- clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    df <- as.data.frame(ek_read)
    
    # --- 4) FDR filter ---
    df <- subset(df, p.adjust <= p_adj_cut)
    if (!nrow(df)) {
      if (verbose) message(sprintf("[%s] No KEGG pathways with FDR <= %.3f.", d, p_adj_cut))
      return(NULL)
    }
    
    # --- 5) Extra annotations ---
    df$Disease  <- d
    df$N_input  <- length(g_sym)
    df$N_mapped <- length(g_ent)
    df
  })
  
  res <- dplyr::bind_rows(out_list)
  if (is.null(res) || !nrow(res)) {
    if (verbose) message("No enriched KEGG terms across diseases (after filtering).")
    return(tibble::tibble())
  }
  res
}



# ========================================================
# Function: add_deg_stars
# Purpose:
#   Append significance stars based on adjusted p-values, add a gene label
#   with stars, and return the data frame ordered by logFC.
#
# Inputs:
#   deg_df    - data frame with at least p_col and logfc_col
#   gene_col  - column name with gene IDs; if NULL, uses row names
#   p_col     - adjusted p-value column (default "adj.P.Val")
#   logfc_col - log2 fold-change column (default "logFC")
#   cuts      - ascending p-value thresholds (defaults c(1e-4, 1e-3, 5e-2))
#   symbols   - star symbols corresponding to cuts (same length as cuts)
#   decreasing- sort order for logFC (default TRUE)
#
# Outputs:
#   The input data frame with two added columns:
#     - stars: "***", "**", "*", or "" according to cuts
#     - gene_label: "<gene><space><stars>" when stars present
#   Rows are ordered by logFC (NA rows dropped by na.last = NA).
# Notes:
#   - Thresholds are inclusive (<=).
#   - No explicit check that length(cuts) == length(symbols).
# ========================================================
add_deg_stars <- function(deg_df,
                          gene_col   = NULL,          # NULL -> uses rownames
                          p_col      = "adj.P.Val",
                          logfc_col  = "logFC",
                          cuts       = c(1e-4, 1e-3, 5e-2),   # thresholds
                          symbols    = c("***", "**", "*"),
                          decreasing = TRUE) {
  
  stopifnot(all(c(p_col, logfc_col) %in% colnames(deg_df)))
  df <- deg_df
  
  # vector of adjusted p-values
  p <- as.numeric(df[[p_col]])
  
  # stars mapping: <=1e-4 '***', <=1e-3 '**', <=5e-2 '*', else ""
  stars <- ifelse(p <= cuts[1], symbols[1],
                  ifelse(p <= cuts[2], symbols[2],
                         ifelse(p <= cuts[3], symbols[3], "")))
  
  df$stars <- stars
  
  # optional: label with stars (handy for heatmaps/tables)
  genes <- if (is.null(gene_col)) rownames(df) else as.character(df[[gene_col]])
  df$gene_label <- paste0(genes, ifelse(stars == "", "", paste0(" ", stars)))
  
  # sort by logFC
  o <- order(df[[logfc_col]], decreasing = decreasing, na.last = NA)
  df[o, , drop = FALSE]
}








