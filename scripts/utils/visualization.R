
# ========================================================
# Statement on transparency and reproducibility
#
# As advised during the Master’s coursework, the code that generates every
# figure is provided to promote reproducibility and prevent any inadvertent
# manipulation of results. This file consolidates the visualization routines
# used to produce all images in the manuscript/analysis.
#
# To enhance clarity and visual consistency, I used generative AI tools —
# specifically ChatGPT (https://chatgpt.com/) — to assist with aesthetic refinements 
# of plotting functions. All AI-assisted edits were reviewed and validated by the author; 
# no underlying statistical analyses or results were modified by these changes. Methodological
# choices, interpretation, and validation remain entirely my responsibility.


library(VennDiagram)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(ggVennDiagram)
library(Cairo)
library(pheatmap)
library(ggrepel)
library(reshape2)
library(dplyr)
library(ggtext)
library(tidytext)
library(scales)
library(networkD3)
library(htmlwidgets)
library(circlize)
library(dplyr)
library(scales)
library(grid)         
library(ComplexHeatmap) 


# Plot: 2-set Venn diagram for two gene sets (with custom colors/alpha).
# Shows counts and category labels; optional title at the top.
plot_venn_consistently_selected_2 <- function(
    selected_genes_list,
    category_names,
    plot_title = "",
    set_colors = c("#66C2A5", "#E78AC3"),
    alpha_fill = 0.35,
    count_cex = 1.3,
    cat_cex = 1
) {
  stopifnot(length(selected_genes_list) == 2L, length(category_names) == 2L)
  if (length(set_colors) < 2L) stop("Provide at least 2 colors in set_colors.")
  names(selected_genes_list) <- category_names
  
  # --- Silence futile.logger globally + common VennDiagram loggers ---
  old_root <- futile.logger::flog.threshold()
  old_vd   <- futile.logger::flog.threshold(name = "VennDiagram")
  old_vdl  <- futile.logger::flog.threshold(name = "VennDiagramLogger")
  on.exit({
    futile.logger::flog.threshold(old_root)
    futile.logger::flog.threshold(old_vd,  name = "VennDiagram")
    futile.logger::flog.threshold(old_vdl, name = "VennDiagramLogger")
  }, add = TRUE)
  futile.logger::flog.threshold(futile.logger::FATAL)  # root
  futile.logger::flog.threshold(futile.logger::FATAL, name = "VennDiagram")
  futile.logger::flog.threshold(futile.logger::FATAL, name = "VennDiagramLogger")
  
  # --- Build grob quietly (hide output/messages/warnings) ---
  venn_grob <- NULL
  suppressWarnings(suppressMessages({
    invisible(capture.output({
      venn_grob <- VennDiagram::venn.diagram(
        x        = selected_genes_list,
        filename = NULL,
        output   = FALSE,
        disable.logging = TRUE,
        scaled   = FALSE,
        fill     = set_colors[1:2],
        alpha    = alpha_fill,
        col      = set_colors[1:2],
        lwd      = 2,
        cex      = count_cex,
        fontface = "plain",
        cat.col      = set_colors[1:2],
        cat.cex      = cat_cex,
        cat.fontface = "plain",
        cat.pos  = c(180, 0),
        cat.dist = c(0.06, 0.06),
        cat.just = list(c(1, 0.5), c(0, 0.5)),
        margin = 0.05
      )
    }))
  }))
  
  # --- Draw ---
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(0.95, "snpc"),
                                    height = grid::unit(0.95, "snpc")))
  grid::grid.draw(venn_grob)
  grid::popViewport()
  
  if (nzchar(plot_title)) {
    grid::grid.text(plot_title, x = 0.5, y = 0.97,
                    gp = grid::gpar(fontsize = 9, fontface = "plain"))
  }
  invisible(NULL)
}






# Plot: 3-set Venn diagram for three gene sets (custom colors/alpha).
# Shows counts and category labels; optional title at the top.
plot_venn_consistently_selected_3 <- function(
    selected_genes_list,
    category_names,
    plot_title = "",
    set_colors = c("#8DA0CB", "#FC8D62", "#66C2A5"),
    alpha_fill = 0.35,
    count_cex = 1.3,
    cat_cex = 1
) {
  stopifnot(length(selected_genes_list) == length(category_names))
  names(selected_genes_list) <- category_names
  
  # ---- silence VennDiagram logger (futile.logger) ----
  old_thr <- futile.logger::flog.threshold()
  on.exit(futile.logger::flog.threshold(old_thr), add = TRUE)
  futile.logger::flog.threshold(futile.logger::FATAL, name = "VennDiagram")
  
  # Capture any output/messages produced when creating the grob
  suppressWarnings(suppressMessages({
    invisible(capture.output({
      venn_grob <- VennDiagram::venn.diagram(
        x        = selected_genes_list,
        filename = NULL,
        output   = FALSE,
        disable.logging = TRUE,
        scaled   = FALSE,
        fill     = set_colors[seq_along(selected_genes_list)],
        alpha    = alpha_fill,
        col      = set_colors[seq_along(selected_genes_list)],
        lwd      = 2,
        cex      = count_cex,
        fontface = "plain",
        cat.col  = set_colors[seq_along(selected_genes_list)],
        cat.cex  = cat_cex,
        cat.fontface = "plain",
        margin = 0.05,
        cat.pos  = c(-10, 10, 180),
        cat.dist = c(0.04, 0.04, 0.04)
      )
    }))
  }))
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(0.95, "snpc"),
                                    height = grid::unit(0.95, "snpc")))
  grid::grid.draw(venn_grob)
  grid::popViewport()
  
  if (nzchar(plot_title)) {
    grid::grid.text(plot_title, x = 0.5, y = 0.97,
                    gp = grid::gpar(fontsize = 12, fontface = "plain"))
  }
  
  invisible(NULL)
}


# Plot: 4-set Venn diagram for four gene sets (custom colors/alpha).
# Shows counts and category labels; optional title at the top.
plot_venn_consistently_selected_4 <- function(
    selected_genes_list,
    category_names,
    plot_title = "",
    set_colors = c("#66C2A5", "#E78AC3", "#FC8D62", "#8DA0CB"),
    alpha_fill = 0.35,
    count_cex = 1.3,
    cat_cex = 1.2
) {
  # --- Safety checks ---
  if (length(selected_genes_list) != 4L)
    stop("This function expects exactly 4 sets (length(selected_genes_list) == 4).")
  if (length(category_names) != 4L)
    stop("Provide 4 category_names (one for each set).")
  if (length(set_colors) < 4L)
    stop("Provide at least 4 colors in set_colors.")
  names(selected_genes_list) <- category_names
  
  # ---- silence VennDiagram (futile.logger) ONLY within this function ----
  old_thr <- futile.logger::flog.threshold()
  old_app <- futile.logger::flog.appender()  # store current appender
  on.exit({
    futile.logger::flog.threshold(old_thr)
    futile.logger::flog.appender(old_app)
  }, add = TRUE)
  futile.logger::flog.threshold(futile.logger::FATAL, name = "VennDiagram")
  futile.logger::flog.appender(futile.logger::appender.console(), name = "VennDiagram")
  
  # --- Build Venn grob (without writing files) ---
  venn_grob <- NULL
  suppressWarnings(suppressMessages({
    invisible(capture.output({
      venn_grob <- VennDiagram::venn.diagram(
        x        = selected_genes_list,
        filename = NULL,     # do not write image
        output   = FALSE,    # do not write sidecar
        disable.logging = TRUE,  # ESSENTIAL to avoid creating the .log
        
        scaled   = FALSE,
        fill     = set_colors[seq_len(4)],
        alpha    = alpha_fill,
        col      = set_colors[seq_len(4)],
        lwd      = 2,
        cex      = count_cex,
        fontface = "plain",
        cat.col       = set_colors[seq_len(4)],
        cat.cex       = cat_cex,
        cat.fontface  = "plain",
        cat.pos  = c(-15, 15, 0, 0),
        cat.dist = c(0.22, 0.22, 0.12, 0.12),
        cat.just = list(c(1, 1), c(0, 1), c(1, 0), c(0, 0)),
        margin = 0.05
      )
    }))
  }))
  
  # --- Draw ---
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(0.95, "snpc"),
                                    height = grid::unit(0.95, "snpc")))
  grid::grid.draw(venn_grob)
  grid::popViewport()
  
  if (nzchar(plot_title)) {
    grid::grid.text(plot_title, x = 0.5, y = 0.97,
                    gp = grid::gpar(fontsize = 14, fontface = "plain"))
  }
  
  invisible(NULL)
}




# Plot: per-gene expression across diseases (jitter + boxplot), faceted by Gene.
# Uses custom colors and enforces a specific x-axis order.
plot_gene_expression_by_disease <- function(
    df,
    disease_colors = NULL,
    disease_order  = c("Gbm.", "Astro.", "Oligo.")
) {
  # Default colors if not provided
  if (is.null(disease_colors)) {
    disease_colors <- c(
      "Gbm."   = "#A58CDE",
      "Astro." = "#DEA58C",
      "Oligo." = "#8CDEA5"
    )
  }
  
  # --- Enforce desired order on x-axis (disease) ---
  # Convert to factor with the requested order
  df$Disease <- factor(df$Disease, levels = disease_order)
  
  ggplot(df, aes(x = Disease, y = Expression, color = Disease)) +
    geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    facet_wrap(~ Gene, scales = "free_y") +
    scale_x_discrete(limits = disease_order, drop = FALSE) +  # keep order
    scale_color_manual(values = disease_colors) +
    labs(
      title = "",
      x = "Disease Type",
      y = "Gene Expression Level"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(size = 5.5),
      axis.title.x = element_text(size = 10, margin = margin(t = 14)),
      axis.title.y = element_text(size = 10, margin = margin(r = 14)),
      plot.title   = element_text(size = 11, margin = margin(b = 13), hjust = 0.5)
    )
}


# Plot: force-directed network (Fruchterman–Reingold) with minimalist styling.
# Unlabeled nodes, thin edges; takes an igraph object and optional title.
plot_graph_network <- function(graph, title = "Graph") {
  layout <- layout_with_fr(graph, niter = 1000, grid = "nogrid")
  
  plot(graph,
       layout = layout,
       vertex.size = 3,
       vertex.label = NA,
       edge.width = 0.5,
       edge.color = "black",
       vertex.color = "gold",
       main = title)
}


# Parse lambda/alpha numbers from penalty names; returns data.frame(lambda, alpha).
parse_lambda_alpha <- function(x) {
  # returns data.frame(lambda, alpha) parsed from penalty names
  get_num <- function(pattern) {
    m <- regexpr(pattern, x, perl = TRUE, ignore.case = TRUE)
    out <- rep(NA_real_, length(x))
    hit <- m > 0
    if (any(hit)) {
      s <- regmatches(x, m)
      out[hit] <- as.numeric(sub("^[la]([0-9.]+)$", "\\1",
                                 sub(".*_([la][0-9.]+).*", "\\1", s), perl = TRUE))
    }
    out
  }
  data.frame(
    lambda = get_num("_l[0-9.]+"),
    alpha  = get_num("_a[0-9.]+")
  )
}

# Infer equation family ("linear","log","ratio") from the Penalty field.
infer_equation <- function(df) {
  if ("equation" %in% names(df)) return(df$equation)
  out <- ifelse(grepl("_linear_", df$Penalty, ignore.case = TRUE), "linear",
                ifelse(grepl("_log_",    df$Penalty, ignore.case = TRUE), "log",
                       ifelse(grepl("_ratio_",  df$Penalty, ignore.case = TRUE), "ratio", NA)))
  factor(out, levels = c("linear","log","ratio"))
}

# Infer distance label from Penalty (token before the first underscore).
infer_distance <- function(df) {
  if ("distance" %in% names(df)) return(df$distance)
  # token before first underscore in Penalty
  d <- sub("_.*$", "", df$Penalty)
  d
}

# Plot: tuning lines per equation (linear/log/ratio): metric vs λ (or α) colored by distance.
# Returns a list of three ggplot objects: $linear, $log, $ratio.
plot_tuning_by_equation <- function(summary_tbl,
                                    metric = "Median_PR_AUC_Test",
                                    distances_order = c("euclidean","manhattan","canberra","bray"),
                                    distance_colors = NULL) {
  stopifnot(is.data.frame(summary_tbl), metric %in% names(summary_tbl), "Penalty" %in% names(summary_tbl))
  
  df <- summary_tbl
  df$equation <- infer_equation(df)
  df$distance <- infer_distance(df)
  pa <- parse_lambda_alpha(df$Penalty)
  df$lambda <- pa$lambda
  df$alpha  <- pa$alpha
  
  present <- intersect(distances_order, unique(df$distance))
  df$distance <- factor(df$distance, levels = present)
  
  df_lin <- subset(df, equation == "linear" & is.finite(lambda))
  df_log <- subset(df, equation == "log"    & is.finite(lambda))
  df_rat <- subset(df, equation == "ratio"  & is.finite(alpha))
  
  yax <- if (metric == "Median_PR_AUC_Test") "PR-AUC (test)" else metric
  base_theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey90"),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
      legend.position = "top",
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20))
    )
  
  color_scale <- if (!is.null(distance_colors)) {
    ggplot2::scale_color_manual(values = distance_colors, name = "Distance")
  } else {
    ggplot2::scale_color_brewer(palette = "Set1", name = "Distance")
  }
  
  p_linear <- ggplot2::ggplot(df_lin,
                              ggplot2::aes(x = lambda, y = .data[[metric]], color = distance, group = distance)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df_lin$lambda))) +
    ggplot2::labs(title = NULL, x = expression(lambda), y = yax) +
    color_scale + base_theme
  
  p_log <- ggplot2::ggplot(df_log,
                           ggplot2::aes(x = lambda, y = .data[[metric]], color = distance, group = distance)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df_log$lambda))) +
    ggplot2::labs(title = NULL, x = expression(lambda), y = yax) +
    color_scale + base_theme
  
  p_ratio <- ggplot2::ggplot(df_rat,
                             ggplot2::aes(x = alpha, y = .data[[metric]], color = distance, group = distance)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df_rat$alpha))) +
    ggplot2::labs(title = NULL, x = expression(alpha), y = yax) +
    color_scale + base_theme
  
  list(linear = p_linear, log = p_log, ratio = p_ratio)
}


# Plot: histogram + density overlay of a chosen distance column.
# Inputs: df (data.frame), x_col_name (string), axis labels and title.
plot_distances_distribution <- function(df, x_col_name, x_label = "Distance", title = "Distance Distribution") {
  ggplot(df, aes_string(x = x_col_name)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "white", alpha = 0.6) +
    geom_density(color = "red", size = 1.2) +
    labs(title = title,
         x = x_label,
         y = "Density") +
    theme_minimal()
}

# Plot: scatter of weights vs. coefficients with text labels; colors indicate Elastic Net overlap.
# Returns a ggplot; optional legend on top when show_legend = TRUE.
plot_weights_vs_coefficients <- function(
    df_plot,
    elastic_genes = NULL,
    set_color    = "#FC8D62",
    legend_title = "Selection status",
    show_legend  = FALSE
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install 'ggrepel' to use this function.")
  }
  
  df_plot$ElasticNet <- ifelse(df_plot$Gene %in% elastic_genes,
                               "Also by Elastic Net", "Only by this model")
  df_plot$ElasticNet <- factor(df_plot$ElasticNet,
                               levels = c("Only by this model", "Also by Elastic Net"))
  
  color_map <- c(
    "Only by this model"  = set_color,
    "Also by Elastic Net" = "#228B22"
  )
  
  p <- ggplot2::ggplot(
    df_plot, ggplot2::aes(x = Weight, y = Coefficient, label = Gene, color = ElasticNet)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 2) +
    ggrepel::geom_text_repel(size = 2.5, show.legend = FALSE, max.overlaps = Inf) +
    ggplot2::scale_color_manual(name = legend_title, values = color_map) +
    ggplot2::labs(x = "GeNeIV Euclidean-Log Weights", y = "Regression Coefficient") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (show_legend) "top" else "none",  # legend hidden by default
      legend.title    = ggplot2::element_text(size = 10),
      legend.text     = ggplot2::element_text(size = 11),
      axis.title.x    = ggplot2::element_text(size = 13, margin = ggplot2::margin(t = 14)),
      axis.title.y    = ggplot2::element_text(size = 13, margin = ggplot2::margin(r = 14)),
      axis.text.x     = ggplot2::element_text(size = 9),
      axis.text.y     = ggplot2::element_text(size = 9)
    ) +
    ggplot2::guides(color = if (show_legend) ggplot2::guide_legend() else "none")
  
  p  # return ggplot object
}




# Plot: bar chart of gene selection frequency (top N) for a given penalty (bootstraps).
plot_gene_selection_frequency <- function(results_list, penalty_name, top_n = 20) {
  # Extract all selected genes across bootstraps for given penalty
  all_genes <- unlist(results_list[[penalty_name]]$var_selected_names)
  
  # Count frequency
  gene_counts <- sort(table(all_genes), decreasing = TRUE)
  
  # Convert to data frame
  df_counts <- as.data.frame(gene_counts)
  colnames(df_counts) <- c("Gene", "Frequency")
  
  # Keep top_n most frequent
  df_counts <- head(df_counts, top_n)
  
  # Plot
  ggplot(df_counts, aes(x = reorder(Gene, Frequency), y = Frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = paste("Gene Selection Frequency -", penalty_name),
      x = "Gene",
      y = "Frequency"
    ) +
    theme_minimal(base_size = 12)
}



# Plot: PR-AUC (test) comparison across models — boxplot + jitter, ordered by median.
plot_test_auc_boxplots <- function(models, color_set = c("#FC8D62","#8DA0CB", "#66C2A5", "#CC79A7",
                                                         "#0072B2", "#D55E00", "#999999")) {
  # 1. Build dataframe
  df_list <- lapply(names(models), function(label) {
    item <- models[[label]]
    if (is.numeric(item)) {
      return(data.frame(Model = label, Value = as.numeric(item)))
    }
    if (is.list(item) && !is.null(item$pr_auc_test)) {
      return(data.frame(Model = label, Value = as.numeric(item$pr_auc_test)))
    }
    stop(sprintf("Invalid input for '%s'.", label))
  })
  df_auc <- do.call(rbind, df_list)
  
  # 2. Order by median
  med_auc <- aggregate(Value ~ Model, df_auc, median, na.rm = TRUE)
  ord <- med_auc$Model[order(-med_auc$Value)]
  df_auc$Model <- factor(df_auc$Model, levels = ord)
  
  # 3. Assign one color per model (internal palette)
  model_names <- levels(df_auc$Model)
  n_models <- length(model_names)
  default_colors <- color_set
  color_vector <- setNames(default_colors[seq_len(n_models)], model_names)
  
  # 4. Plot
  ggplot(df_auc, aes(x = Model, y = Value, color = Model, fill = Model)) +
    geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    scale_color_manual(values = color_vector) +
    scale_fill_manual(values = color_vector) +
    labs(x = NULL, y = "PR-AUC (Test)", title = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(size = 12,  hjust = 0.5),  # smaller X ticks
      axis.text.y  = element_text(size =9),                  # smaller Y ticks
      axis.title.x = element_text(size = 14, margin = margin(t = 14)),
      axis.title.y = element_text(size = 14, margin = margin(r = 14)),
      plot.margin  = margin(10, 16, 30, 16),
      plot.title   = element_text(hjust = 0.5, vjust = -10, face = "bold", size = 20)
    ) +
    coord_cartesian(clip = "off")
}


# Prep: build log2(+1) expression matrix + sample annotations for heatmap; optional row z-scoring.
build_heatmap_matrix <- function(astro, gbm, oligo, genes, zscore_rows = TRUE) {
  common_genes <- Reduce(intersect, list(genes, colnames(astro), colnames(gbm), colnames(oligo)))
  if (length(common_genes) == 0) stop("No selected genes found in all datasets.")
  
  A <- astro[, common_genes, drop = FALSE]
  G <- gbm[,   common_genes, drop = FALSE]
  O <- oligo[, common_genes, drop = FALSE]
  
  A <- log2(A + 1); G <- log2(G + 1); O <- log2(O + 1)
  
  X <- rbind(G, A, O)
  disease <- c(rep("GBM", nrow(G)), rep("Astro", nrow(A)), rep("Oligo", nrow(O)))
  col_annot <- data.frame(Disease = disease, row.names = rownames(X))
  X <- t(X)
  
  if (zscore_rows) {
    X <- t(scale(t(X)))
    X[is.na(X)] <- 0
  }
  list(matrix = X, ann_col = col_annot)
}




# Plot: gene expression heatmap with optional DEG star labels/order and disease annotations.
plot_gene_heatmap <- function(
    X_list,
    title = NULL,
    res_star = NULL,            # pre-ordered DEGs df; rownames = genes; column "stars"
    gene_order = NULL,          # fallback if res_star is NULL
    append_stars = TRUE,
    cluster_rows = TRUE,        # ignored if res_star/gene_order provided
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_disease_legend = TRUE,
    show_heat_legend = TRUE,
    # base palette
    colors = c("#5D50CE", "white", "#ED3431"),
    n_colors = 101,
    # backward compatibility:
    use_brewer = FALSE,
    brewer_name = "BrBG",
    brewer_direction = 1,       # 1 = normal, -1 = reversed
    # text
    fontsize_row = 8,
    fontsize_title = 12,
    title_margin_mm = 2
) {
  X   <- X_list$matrix
  ann <- X_list$ann_col
  
  # --- order by res_star (preferred) ---
  labels_row <- NULL
  if (!is.null(res_star)) {
    genes_star <- rownames(res_star)
    idx <- match(genes_star, rownames(X))
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      X <- X[idx, , drop = FALSE]
      cluster_rows <- FALSE
    }
    if (append_stars && "stars" %in% colnames(res_star) && show_rownames) {
      s <- res_star$stars; names(s) <- genes_star
      s_now <- s[rownames(X)]; s_now[is.na(s_now)] <- ""
      labels_row <- paste0(rownames(X), ifelse(s_now == "", "", paste0(" ", s_now)))
    }
  } else if (!is.null(gene_order)) {
    idx <- match(gene_order, rownames(X))
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      X <- X[idx, , drop = FALSE]
      cluster_rows <- FALSE
    }
  } else {
    X <- X[order(rownames(X)), , drop = FALSE]
  }
  
  # --- palette ---
  if (isTRUE(use_brewer)) {
    base <- RColorBrewer::brewer.pal(
      RColorBrewer::brewer.pal.info[brewer_name, "maxcolors"], brewer_name
    )
    if (brewer_direction == -1) base <- rev(base)
    pal <- grDevices::colorRampPalette(base)(n_colors)
  } else {
    pal <- grDevices::colorRampPalette(colors)(n_colors)
  }
  
  ann_colors <- list(Disease = c(GBM="#A58CDE", Astro="#DEA58C", Oligo="#8CDEA5"))
  
  hm <- pheatmap::pheatmap(
    mat = X,
    color = pal,
    annotation_col = ann,
    annotation_colors = ann_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = FALSE,
    labels_row = labels_row,
    fontsize_row = fontsize_row,
    border_color = NA,
    legend = show_heat_legend,
    annotation_legend = show_disease_legend,
    main = ""
  )
  
  grid::grid.text(
    label = title,
    y = grid::unit(1, "npc") - grid::unit(title_margin_mm, "mm"),
    gp = grid::gpar(fontsize = fontsize_title, fontface = "plain")
  )
  
  invisible(hm)
}




# Build: per-gene TOP indicator and score matrices across diseases × metrics (top_frac threshold).
build_top_indicator_matrix <- function(network_df,
                                       diseases = c("Astro","GBM","Oligo"),
                                       metrics = NULL,
                                       top_frac = 0.10) {
  genes <- rownames(network_df)
  
  strip_metric <- function(x) sub("\\..*$", "", x)
  if (is.null(metrics)) {
    metrics <- Reduce(
      union,
      lapply(diseases, function(d)
        unique(strip_metric(unique(unlist(lapply(network_df[[d]], names)))))) )
  }
  
  colnames_all <- as.vector(outer(diseases, metrics, paste, sep = "_"))
  score_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(colnames_all),
                      dimnames = list(genes, colnames_all))
  bin_mat   <- score_mat
  
  for (d in diseases) {
    for (m in metrics) {
      vec <- sapply(network_df[[d]], function(v) {
        nm <- strip_metric(names(v))
        idx <- match(m, nm)
        if (is.na(idx)) NA_real_ else as.numeric(v[idx])
      })
      names(vec) <- genes
      
      thr <- stats::quantile(vec, probs = 1 - top_frac, na.rm = TRUE, type = 7)
      colname <- paste(d, m, sep = "_")
      
      score_mat[, colname] <- vec
      # Mark TOP if: non-NA, >= threshold AND strictly > 0
      bin_mat[, colname]   <- ifelse(!is.na(vec) & vec >= thr & vec > 0, 1, 0)
    }
  }
  
  list(binary = bin_mat, scores = score_mat,
       metrics = metrics, diseases = diseases)
}



# Plot: binary “TOP” heatmap (genes × metrics), faceted by disease; supports sorting/highlighting.
plot_top_heatmap <- function(
    bin_mat,
    genes_to_plot = NULL,
    diseases_order = c("Astro","GBM","Oligo"),
    metrics_order  = c("degree","betweenness","harmonic","eigen","triangles"),
    sort_by_hits   = FALSE,
    plot_title     = NULL,
    disease_colors = c(Astro = "#DEA58C", GBM = "#A58CDE", Oligo = "#8CDEA5"),
    show_legend    = FALSE,
    legend_title   = "Top status",
    legend_labels  = list(not_top = "Not top", top_prefix = "Top in "),
    disease_labels = c(Astro = "Astrocytoma", GBM = "Glioblastoma", Oligo = "Oligodendroglioma"),
    metric_labels  = c(degree="D.", betweenness="B.", harmonic="H.", eigen="E.", triangles="T."),
    highlight_genes = NULL,
    strip_text_size = 14   # <-- new: facet strip font size (disease names)
) {
  stopifnot(is.matrix(bin_mat))
  if (!all(diseases_order %in% names(disease_colors))) {
    stop("disease_colors must be named for: ", paste(diseases_order, collapse=", "))
  }
  
  if (!is.null(genes_to_plot)) {
    genes_to_plot <- intersect(genes_to_plot, rownames(bin_mat))
    bin_mat <- bin_mat[genes_to_plot, , drop = FALSE]
  }
  if (isTRUE(sort_by_hits)) {
    ord <- order(rowSums(bin_mat, na.rm = TRUE), decreasing = TRUE)
    bin_mat <- bin_mat[ord, , drop = FALSE]
  }
  
  gene_order <- rev(sort(rownames(bin_mat)))
  bin_mat <- bin_mat[gene_order, , drop = FALSE]
  
  df_long <- reshape2::melt(bin_mat, varnames = c("Gene","Col"), value.name = "Top")
  sp <- do.call(rbind, strsplit(as.character(df_long$Col), "_", fixed = TRUE))
  df_long$Disease <- factor(sp[,1], levels = diseases_order)
  df_long$Metric  <- factor(sp[,2], levels = metrics_order)
  df_long$Top     <- as.integer(df_long$Top)
  
  df_long$FillKey <- ifelse(df_long$Top == 1,
                            paste0("Top_", as.character(df_long$Disease)),
                            "Not_Top")
  
  fill_values <- c(
    Not_Top = "white",
    setNames(disease_colors[as.character(diseases_order)],
             paste0("Top_", diseases_order))
  )
  fill_labels <- c(
    Not_Top = legend_labels$not_top,
    setNames(paste0(legend_labels$top_prefix, diseases_order),
             paste0("Top_", diseases_order))
  )
  
  dis_labeller <- ggplot2::as_labeller(disease_labels)
  x_labeler <- function(x) unname(ifelse(x %in% names(metric_labels), metric_labels[x], x))
  
  highlight_set <- intersect(rownames(bin_mat), highlight_genes %||% character(0))
  y_lab_fun <- function(y) {
    y <- as.character(y)
    is_hi <- y %in% highlight_set
    y[is_hi] <- paste0("**", y[is_hi], "**")
    y
  }
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = Metric, y = Gene, fill = FillKey)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::facet_grid(
      . ~ Disease, scales = "free_x", space = "free_x",
      labeller = ggplot2::labeller(Disease = dis_labeller)
    ) +
    ggplot2::scale_x_discrete(labels = x_labeler) +
    ggplot2::scale_y_discrete(labels = y_lab_fun) +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = names(fill_values),
      labels = fill_labels,
      guide  = if (isTRUE(show_legend)) ggplot2::guide_legend(override.aes = list(color = NA)) else "none",
      name   = legend_title
    ) +
    ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      strip.text.x    = ggplot2::element_text(
        face = "bold", size = strip_text_size,
        margin = ggplot2::margin(b = 12)   # add space between facet strip and grid
      ),
      axis.text.x     = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y     = ggtext::element_markdown(size = 9)
    )
}



# Build: per-gene TOP indicator and score matrices across diseases × metrics (top_frac threshold).
build_top_indicator_matrix <- function(network_df,
                                       diseases = c("Astro","GBM","Oligo"),
                                       metrics = NULL,
                                       top_frac = 0.10) {
  genes <- rownames(network_df)
  
  strip_metric <- function(x) sub("\\..*$", "", x)
  if (is.null(metrics)) {
    metrics <- Reduce(
      union,
      lapply(diseases, function(d)
        unique(strip_metric(unique(unlist(lapply(network_df[[d]], names)))))) )
  }
  
  colnames_all <- as.vector(outer(diseases, metrics, paste, sep = "_"))
  score_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(colnames_all),
                      dimnames = list(genes, colnames_all))
  bin_mat   <- score_mat
  
  for (d in diseases) {
    for (m in metrics) {
      vec <- sapply(network_df[[d]], function(v) {
        nm <- strip_metric(names(v))
        idx <- match(m, nm)
        if (is.na(idx)) NA_real_ else as.numeric(v[idx])
      })
      names(vec) <- genes
      
      thr <- stats::quantile(vec, probs = 1 - top_frac, na.rm = TRUE, type = 7)
      colname <- paste(d, m, sep = "_")
      
      score_mat[, colname] <- vec
      # Mark TOP if: non-NA, >= threshold AND strictly > 0
      bin_mat[, colname]   <- ifelse(!is.na(vec) & vec >= thr & vec > 0, 1, 0)
    }
  }
  
  list(binary = bin_mat, scores = score_mat,
       metrics = metrics, diseases = diseases)
}

# Plot: binary “TOP” heatmap (genes × metrics), faceted by disease; supports sorting/highlighting.
plot_top_heatmap <- function(
    bin_mat,
    genes_to_plot = NULL,
    diseases_order = c("Astro","GBM","Oligo"),
    metrics_order  = c("degree","betweenness","harmonic","eigen","triangles"),
    sort_by_hits   = FALSE,
    plot_title     = NULL,
    disease_colors = c(Astro = "#DEA58C", GBM = "#A58CDE", Oligo = "#8CDEA5"),
    show_legend    = FALSE,
    legend_title   = "Top status",
    legend_labels  = list(not_top = "Not top", top_prefix = "Top in "),
    disease_labels = c(Astro = "Astrocytoma", GBM = "Glioblastoma", Oligo = "Oligodendroglioma"),
    metric_labels  = c(degree="D.", betweenness="B.", harmonic="H.", eigen="E.", triangles="T."),
    highlight_genes = NULL,
    strip_text_size = 14   # <-- new: facet strip font size (disease names)
) {
  stopifnot(is.matrix(bin_mat))
  if (!all(diseases_order %in% names(disease_colors))) {
    stop("disease_colors must be named for: ", paste(diseases_order, collapse=", "))
  }
  
  if (!is.null(genes_to_plot)) {
    genes_to_plot <- intersect(genes_to_plot, rownames(bin_mat))
    bin_mat <- bin_mat[genes_to_plot, , drop = FALSE]
  }
  if (isTRUE(sort_by_hits)) {
    ord <- order(rowSums(bin_mat, na.rm = TRUE), decreasing = TRUE)
    bin_mat <- bin_mat[ord, , drop = FALSE]
  }
  
  gene_order <- rev(sort(rownames(bin_mat)))
  bin_mat <- bin_mat[gene_order, , drop = FALSE]
  
  df_long <- reshape2::melt(bin_mat, varnames = c("Gene","Col"), value.name = "Top")
  sp <- do.call(rbind, strsplit(as.character(df_long$Col), "_", fixed = TRUE))
  df_long$Disease <- factor(sp[,1], levels = diseases_order)
  df_long$Metric  <- factor(sp[,2], levels = metrics_order)
  df_long$Top     <- as.integer(df_long$Top)
  
  df_long$FillKey <- ifelse(df_long$Top == 1,
                            paste0("Top_", as.character(df_long$Disease)),
                            "Not_Top")
  
  fill_values <- c(
    Not_Top = "white",
    setNames(disease_colors[as.character(diseases_order)],
             paste0("Top_", diseases_order))
  )
  fill_labels <- c(
    Not_Top = legend_labels$not_top,
    setNames(paste0(legend_labels$top_prefix, diseases_order),
             paste0("Top_", diseases_order))
  )
  
  dis_labeller <- ggplot2::as_labeller(disease_labels)
  x_labeler <- function(x) unname(ifelse(x %in% names(metric_labels), metric_labels[x], x))
  
  highlight_set <- intersect(rownames(bin_mat), highlight_genes %||% character(0))
  y_lab_fun <- function(y) {
    y <- as.character(y)
    is_hi <- y %in% highlight_set
    y[is_hi] <- paste0("**", y[is_hi], "**")
    y
  }
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = Metric, y = Gene, fill = FillKey)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
    ggplot2::facet_grid(
      . ~ Disease, scales = "free_x", space = "free_x",
      labeller = ggplot2::labeller(Disease = dis_labeller)
    ) +
    ggplot2::scale_x_discrete(labels = x_labeler) +
    ggplot2::scale_y_discrete(labels = y_lab_fun) +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = names(fill_values),
      labels = fill_labels,
      guide  = if (isTRUE(show_legend)) ggplot2::guide_legend(override.aes = list(color = NA)) else "none",
      name   = legend_title
    ) +
    ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      strip.text.x    = ggplot2::element_text(
        face = "bold", size = strip_text_size,
        margin = ggplot2::margin(b = 12)   # add space between facet strip and grid
      ),
      axis.text.x     = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y     = ggtext::element_markdown(size = 9)
    )
}

# Plot: side-by-side bar chart of per-gene selection counts across models (top-N or provided list).
plot_selection_counts_multi <- function(models_list,
                                        genes = NULL,
                                        top_n = 20,
                                        order = c("desc","asc","input"),
                                        palette = c("#66C2A5", "#8DA0CB","#FC8D62", "#CC79A7"),
                                        title = "Gene selection counts (side-by-side)",
                                        xlab = NULL,
                                        ylab = "Selections (count)") {
  stopifnot(is.list(models_list), length(models_list) >= 1)
  if (is.null(names(models_list)) || any(names(models_list) == "")) {
    stop("Please provide a *named* list for models_list (names are used as legend labels).")
  }
  order <- match.arg(order)
  
  # 1) Build long data.frame: Gene, Count, Model
  dfs <- lapply(names(models_list), function(nm) {
    df <- selection_counts_one(models_list[[nm]])
    df$Model <- nm
    df
  })
  long_all <- do.call(rbind, dfs)
  
  # 2) Define the gene set to plot
  if (is.null(genes)) {
    # pick genes with highest max count across models
    agg <- aggregate(Count ~ Gene, long_all, max)
    keep <- head(agg[order(-agg$Count), "Gene"], top_n)
  } else {
    keep <- intersect(genes, unique(long_all$Gene))
  }
  long <- long_all[long_all$Gene %in% keep, , drop = FALSE]
  
  # 3) Fill missing (Model, Gene) combos with 0 counts (so bars appear)
  template <- expand.grid(Gene = keep, Model = names(models_list), stringsAsFactors = FALSE)
  long <- merge(template, long, by = c("Gene","Model"), all.x = TRUE)
  long$Count[is.na(long$Count)] <- 0L
  
  # 4) Order genes on x-axis
  if (!is.null(genes) && order == "input") {
    # keep user-specified order
    long$Gene <- factor(long$Gene, levels = genes[genes %in% keep])
  } else {
    # order by total (or ascending if requested)
    tot <- aggregate(Count ~ Gene, long, sum)
    lev <- if (order == "desc") {
      tot$Gene[order(-tot$Count)]
    } else {
      tot$Gene[order(tot$Count)]
    }
    long$Gene <- factor(long$Gene, levels = lev)
  }
  
  # 5) Side-by-side bars
  ggplot2::ggplot(long, ggplot2::aes(x = Gene, y = Count, fill = Model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_fill_manual(values = rep(palette, length.out = length(unique(long$Model)))) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = xlab, y = ylab, fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "plain", size = 10, hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 11)
    )
}

# Plot: histogram of weights with dashed vertical lines for selected genes (two groups/legend optional).
plot_weight_hist_with_gene_lines <- function(
    gene_scores_df, genes_lines, other_set_genes = NULL,
    weight_col = "custom_score_log", binwidth = 0.02,
    hist_fill = "grey", density_color = "grey",
    base_line_color = "#FC8D62", highlight_line_color = "#66C2A5",
    show_legend = FALSE,
    title = "Distance Distribution with Selected Genes",
    x_label = "Distance", y_label = "Count",
    set_a_name = "Only set A", set_b_name = "Also in set B",
    axis_text_size = 11,    # <-- new: tick label size
    axis_title_size = 13    # <-- new: axis title size
) {
  if (is.null(rownames(gene_scores_df)))
    stop("gene_scores_df must have rownames with gene names.")
  if (!weight_col %in% colnames(gene_scores_df))
    stop(sprintf("Column '%s' not found in gene_scores_df.", weight_col))
  
  all_weights <- as.numeric(gene_scores_df[[weight_col]])
  df_all <- data.frame(weight = all_weights)
  df_all <- df_all[is.finite(df_all$weight), , drop = FALSE]
  
  genes_lines <- intersect(genes_lines, rownames(gene_scores_df))
  if (length(genes_lines) == 0)
    stop("None of 'genes_lines' found in gene_scores_df rownames.")
  
  sel_w <- as.numeric(gene_scores_df[genes_lines, weight_col])
  names(sel_w) <- genes_lines
  sel_w <- sel_w[is.finite(sel_w)]
  genes_lines <- names(sel_w)
  
  group <- rep(set_a_name, length(sel_w))
  if (!is.null(other_set_genes)) {
    group[genes_lines %in% other_set_genes] <- set_b_name
  }
  df_lines <- data.frame(
    gene  = genes_lines,
    x     = sel_w,
    group = factor(group, levels = c(set_a_name, set_b_name))
  )
  
  ggplot(df_all, aes(x = weight)) +
    geom_histogram(aes(y = ..count..), binwidth = binwidth,
                   fill = hist_fill, alpha = 0.6) +
    geom_vline(data = df_lines,
               aes(xintercept = x, color = group),
               linetype = "dashed", linewidth = 0.9,
               show.legend = show_legend) +
    scale_color_manual(values = setNames(c(base_line_color, highlight_line_color),
                                         c(set_a_name, set_b_name))) +
    labs(title = title, x = x_label, y = y_label, color = NULL) +
    theme_minimal() +
    theme(
      legend.position      = if (show_legend) c(0.98, 0.98) else "none",
      legend.justification = c("right", "top"),
      legend.background    = element_rect(
        fill = grDevices::adjustcolor("white", alpha.f = 0.75), colour = NA
      ),
      legend.key.width     = unit(18, "pt"),
      legend.margin        = margin(3, 4, 3, 4),
      axis.text.x          = element_text(size = axis_text_size),
      axis.text.y          = element_text(size = axis_text_size),
      axis.title.x         = element_text(size = axis_title_size, margin = margin(t = 10)),
      axis.title.y         = element_text(size = axis_title_size, margin = margin(r = 12))
    ) +
    guides(color = if (show_legend) guide_legend(
      override.aes = list(linetype = "dashed", linewidth = 1.2)
    ) else "none") +
    coord_cartesian(clip = "off")
}

# Plot: binary selection heatmap (gene × model) for provided model gene sets; optional legend.
plot_selection_heatmap <- function(
    gene_list,
    model_gene_sets,
    model_colors = c(
      Elastic   = "#66C2A5",
      Manhattan = "#FC8D62",
      Euclidean = "#8DA0CB",
      Twiner    = "#CC79A7"
    ),
    title = NULL,
    show_legend = TRUE   # <- new parameter
) {
  # Ensure model order is preserved
  model_names <- names(model_gene_sets)
  
  # Validate colors
  if (is.null(model_colors)) {
    model_colors <- RColorBrewer::brewer.pal(n = length(model_names), name = "Set2")
    names(model_colors) <- model_names
  } else {
    if (length(model_colors) != length(model_names) || is.null(names(model_colors))) {
      stop("model_colors must be a named vector with the same names as model_gene_sets.")
    }
  }
  
  # Binary matrix (genes x models)
  binary_matrix <- sapply(model_gene_sets, function(gset) gene_list %in% gset)
  rownames(binary_matrix) <- gene_list
  
  # Long format
  df_long <- reshape2::melt(binary_matrix)
  colnames(df_long) <- c("Gene", "Model", "Selected")
  df_long$Selected <- as.logical(df_long$Selected)
  
  # alphabetical order for genes (A at top if you remove rev())
  sorted_genes <- sort(unique(gene_list))
  df_long$Gene  <- factor(df_long$Gene,  levels = rev(sorted_genes))
  df_long$Model <- factor(df_long$Model, levels = model_names)
  
  ggplot(df_long[df_long$Selected, ], aes(x = Model, y = Gene, fill = Model)) +
    geom_tile(color = "white", width = 0.85, height = 1) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    ) +
    guides(fill = if (show_legend) ggplot2::guide_legend() else "none") +
    coord_fixed(ratio = 0.85)
}

# Build: long table of gene presence by module size per disease (for presence heatmaps/plots).
build_presence_by_module <- function(genes,
                                     astro_modules, gbm_modules, oligo_modules) {
  
  # Map a gene to (module_id, module_size) for one disease
  map_gene_to_module <- function(mod_list, genes) {
    if (length(mod_list) == 0) {
      return(data.frame(Gene = character(0),
                        ModuleID = integer(0),
                        ModSize  = integer(0)))
    }
    mod_genes <- lapply(mod_list, `[[`, "genes")
    mod_size  <- vapply(mod_genes, length, integer(1))
    out <- lapply(genes, function(g) {
      hit <- which(vapply(mod_genes, function(v) g %in% v, logical(1)))
      if (length(hit) == 0) return(NULL)
      i <- hit[1]
      data.frame(Gene = g, ModuleID = i, ModSize = mod_size[i], stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
  }
  
  build_df_for_disease <- function(disease, mod_list) {
    gm <- map_gene_to_module(mod_list, genes)
    if (nrow(gm) == 0) {
      return(data.frame(Gene=character(), Disease=character(),
                        ModuleLabel=factor(), Present=integer(),
                        XKey=character(), XLabel=integer()))
    }
    # label = module size; 1 = isolated
    gm$ModuleLabel <- gm$ModSize
    
    # unique sizes (desc), with 1 forced to the end
    lev <- sort(unique(gm$ModuleLabel), decreasing = TRUE)
    if (1L %in% lev) lev <- c(setdiff(lev, 1L), 1L)
    
    # complete grid (genes x labels existing in this disease)
    grid <- expand.grid(Gene = unique(genes),
                        ModuleLabel = lev,
                        stringsAsFactors = FALSE)
    grid$Disease <- disease
    grid$Present <- 0L
    
    # mark presence
    key <- paste(gm$Gene, gm$ModuleLabel)
    grid$Present[paste(grid$Gene, grid$ModuleLabel) %in% key] <- 1L
    
    # factor with per-disease order
    grid$ModuleLabel <- factor(grid$ModuleLabel, levels = lev)
    
    # per-disease x-key (ensures correct order inside facet)
    pos_map <- setNames(seq_along(lev), as.character(lev))
    grid$XKey   <- paste(disease, pos_map[as.character(grid$ModuleLabel)], sep = "::")
    grid$XLabel <- as.integer(as.character(grid$ModuleLabel))  # what we show on axis
    
    grid
  }
  
  df <- rbind(
    build_df_for_disease("Glioblastoma",      gbm_modules),
    build_df_for_disease("Astrocytoma",       astro_modules),
    build_df_for_disease("Oligodendroglioma", oligo_modules)
  )
  
  df$Disease <- factor(df$Disease,
                       levels = c("Glioblastoma","Astrocytoma","Oligodendroglioma"))
  df
}

# Plot: presence heatmap of genes across module-size bins per disease; optional equal tile width & bold highlights.
plot_presence_by_module <- function(
    df_mod,
    disease_colors = c(
      "Glioblastoma"      = "#A58CDE",
      "Astrocytoma"       = "#DEA58C",
      "Oligodendroglioma" = "#8CDEA5"
    ),
    show_numbers   = FALSE,
    highlight_genes = NULL,
    highlight_bold  = FALSE,
    equal_tile_width = TRUE   # <-- NEW
) {
  d <- df_mod[df_mod$Present == 1, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  
  # order genes (reverse alphabetical on Y)
  gene_levels <- rev(sort(unique(d$Gene)))
  
  # map XKey -> XLabel per disease (to order by size desc.)
  if (!"XLabel" %in% names(df_mod)) stop("df_mod needs column 'XLabel'.")
  lab_map_df <- unique(df_mod[, c("Disease","XKey","XLabel")])
  lab_map_df$XLabel <- suppressWarnings(as.numeric(lab_map_df$XLabel))
  
  ord_by_dis <- lapply(split(lab_map_df, lab_map_df$Disease), function(dd) {
    dd <- dd[order(-dd$XLabel, dd$XKey), , drop = FALSE]
    dd$XKey
  })
  
  make_xid <- function(dis, key) paste(dis, key, sep = "__")
  
  # background grid (only keys present per disease)
  bg_list <- lapply(names(ord_by_dis), function(dis) {
    xlev <- ord_by_dis[[dis]]
    expand.grid(
      Disease = dis,
      XKey    = xlev,
      Gene    = gene_levels,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
  })
  bg <- do.call(rbind, bg_list)
  bg$XID <- make_xid(bg$Disease, bg$XKey)
  
  d$XID <- make_xid(d$Disease, d$XKey)
  
  xid_levels <- unlist(lapply(names(ord_by_dis), function(dis) {
    make_xid(dis, ord_by_dis[[dis]])
  }), use.names = FALSE)
  
  # factors/levels
  d$Disease  <- factor(d$Disease,  levels = names(ord_by_dis))
  d$Gene     <- factor(d$Gene,     levels = gene_levels)
  d$XID      <- factor(d$XID,      levels = xid_levels)
  bg$Disease <- factor(bg$Disease, levels = names(ord_by_dis))
  bg$Gene    <- factor(bg$Gene,    levels = gene_levels)
  bg$XID     <- factor(bg$XID,     levels = xid_levels)
  
  # x-axis labels (show XLabel)
  lab_map <- split(setNames(lab_map_df$XLabel, make_xid(lab_map_df$Disease, lab_map_df$XKey)),
                   lab_map_df$Disease)
  x_labeller <- function(keys) {
    dis <- as.character(unique(bg$Disease[bg$XID %in% keys]))[1]
    if (is.na(dis)) return(keys)
    mp <- lab_map[[dis]]
    unname(ifelse(keys %in% names(mp), mp[keys], keys))
  }
  
  y_labeller <- function(genes) {
    if (is.null(highlight_genes) || !isTRUE(highlight_bold)) return(genes)
    ifelse(genes %in% highlight_genes, paste0("<b>", genes, "</b>"), genes)
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = bg,
      ggplot2::aes(x = XID, y = Gene),
      width = 0.95, height = 0.95,
      fill = NA, color = "grey90", linewidth = 0.35
    ) +
    ggplot2::geom_tile(
      data = d,
      ggplot2::aes(x = XID, y = Gene, fill = Disease),
      width = 0.95, height = 0.95,
      color = NA
    ) +
    ggplot2::scale_fill_manual(values = disease_colors) +
    ggplot2::scale_x_discrete(labels = x_labeller) +
    ggplot2::scale_y_discrete(labels = y_labeller) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid      = ggplot2::element_blank(),
      panel.border    = ggplot2::element_blank(),
      strip.text      = ggplot2::element_text(face = "bold", size = 14),
      axis.text.y     = ggtext::element_markdown(size = 9),
      panel.spacing.x = grid::unit(6, "mm")
    )
  
  # facet: equalize tile width (panel width proportional to #columns)
  if (isTRUE(equal_tile_width)) {
    p <- p + ggplot2::facet_grid(~ Disease, scales = "free_x", space = "free_x")
  } else {
    # panels with same width (tiles may vary in width)
    p <- p + ggplot2::facet_wrap(~ Disease, nrow = 1, scales = "free_x")
  }
  
  if (isTRUE(show_numbers)) {
    p <- p + ggplot2::geom_text(
      data = d,
      ggplot2::aes(x = XID, y = Gene, label = XLabel),
      color = "black", size = 3
    )
  }
  
  p
}


# Plot: interactive gene–drug Sankey diagram (networkD3), optional approved-only view and top-K per gene.
sankey_gene_drug <- function(
    df_raw,
    top_k_per_gene = 5,                 # keep top-K drugs per gene (by score)
    width_range = c(2, 15),             # rescale link widths to this range
    font_family = "Arial",
    font_size = 9,
    # default palette; if approved_only=TRUE, we will use domain with Single/Multi automatically
    colour_scale = 'd3.scaleOrdinal()
                    .range(["#CCEDC5","#76B3E5","#C3DEF4",
                            "#FB73B9","#7375FB","#26828EFF",
                            "#31688EFF","#3E4A89FF","#482878FF",
                            "#440154FF"])',
    approved_only = FALSE               # if TRUE: keep only approved drugs and split Single vs Multi
) {
  # ---- 1) Normalize input and optional filter ----
  # Expected columns: gene, drug, regulatory.approval, interaction.score
  df <- df_raw %>%
    rename(score = interaction.score, approved = regulatory.approval) %>%
    mutate(
      gene     = as.character(gene),
      drug     = as.character(drug),
      approved = as.character(approved)
    ) %>%
    filter(!is.na(gene), !is.na(drug))
  
  if (approved_only) {
    df <- df %>%
      filter(tolower(approved) %in% c("true","approved","yes","1"))
  }
  
  df <- df %>%
    group_by(gene, drug, approved) %>%
    summarize(score = max(score, na.rm = TRUE), .groups = "drop")
  
  # Optional: reduce clutter — keep top-K drugs per gene
  if (!is.null(top_k_per_gene) && top_k_per_gene > 0) {
    df <- df %>%
      group_by(gene) %>%
      slice_max(order_by = score, n = top_k_per_gene, with_ties = FALSE) %>%
      ungroup()
  }
  
  # ---- 2) Build nodes with proper grouping ----
  genes <- df %>% distinct(gene) %>% mutate(type = "Gene")
  
  if (approved_only) {
    # Split approved drugs into Single vs Multi degree (after top-k/filtering)
    drug_deg <- df %>% count(drug, name = "deg")
    drugs <- df %>%
      distinct(drug, approved) %>%
      left_join(drug_deg, by = "drug") %>%
      mutate(type = ifelse(deg > 1, "Drug_Approved_Multi", "Drug_Approved_Single")) %>%
      select(-deg)
    
    # Fix colour scale domain explicitly for these groups + Gene
    colour_scale <- 'd3.scaleOrdinal()
                       .domain(["Gene","Drug_Approved_Single","Drug_Approved_Multi"])
                       .range(["#4E9EDF","#C3DEF4","#76B3E5"])'
  } else {
    # Regular grouping: approved vs other
    drugs <- df %>%
      distinct(drug, approved) %>%
      mutate(type = ifelse(tolower(approved) %in% c("true","approved","yes","1"),
                           "Drug_Approved", "Drug_Other"))
    # If you want explicit domain here, you can also set it; we keep user's palette by default
  }
  
  nodes <- bind_rows(
    genes %>% transmute(name = gene, group = type),
    drugs %>% transmute(name = drug, group = type)
  )
  
  # 0-based indices for networkD3
  id_map <- setNames(seq_len(nrow(nodes)) - 1L, nodes$name)
  
  # ---- 3) Links (color by target node group) ----
  links <- df %>%
    transmute(
      source  = id_map[gene],
      target  = id_map[drug],
      value   = pmax(score, 1e-12),
      tooltip = paste0(
        "Gene: ", gene,
        "<br>Drug: ", drug,
        "<br>Approved: ", approved,
        "<br>Interaction score: ", signif(score, 3)
      )
    )
  
  # Link group = target node group (so links match destination node color)
  target_groups <- nodes$group
  links$group <- target_groups[links$target + 1L]
  
  # ---- 4) Type-safety and rescaling ----
  links_df <- as.data.frame(links)
  nodes_df <- as.data.frame(nodes)
  links_df$source <- as.integer(links_df$source)
  links_df$target <- as.integer(links_df$target)
  links_df$value  <- scales::rescale(links_df$value, to = width_range)
  nodes_df$name   <- as.character(nodes_df$name)
  nodes_df$group  <- as.character(nodes_df$group)
  
  # ---- 5) Sankey ----
  w <- sankeyNetwork(
    Links = links_df, Nodes = nodes_df,
    Source = "source", Target = "target", Value = "value",
    NodeID = "name", NodeGroup = "group", LinkGroup = "group",
    fontSize = font_size, fontFamily = font_family,
    nodeWidth = 7, nodePadding = 8,
    sinksRight = FALSE, colourScale = colour_scale
  )
  
  # ---- 6) Tooltips for links ----
  w <- htmlwidgets::onRender(
    w,
    "
    function(el, x){
      d3.select(el).selectAll('.link').select('title').remove();
      d3.select(el).selectAll('.link')
        .append('title')
        .text(function(d){
          return (d.tooltip || (d.source.name + ' \u2192 ' + d.target.name + '\\nScore width: ' + d.value));
        });
    }"
  )
  
  # ---- 7) Font-size rules (Drugs smaller; very long names even smaller) ----
  w <- htmlwidgets::onRender(
    w,
    "
    function(el,x){
      d3.select(el).selectAll('.node text')
        .style('font-size', function(d){
          var L = (d.name || '').length;

          // DRUG nodes (both Single/Multi or Approved/Other)
          if (d.group && d.group.indexOf('Drug') === 0) {
            if (L > 28) return '8px';
            if (L > 18) return '9px';
            return '10px';
          }

          // GENE nodes
          if (d.group && d.group.indexOf('Gene') === 0) {
            if (L > 28) return '9px';
            if (L > 18) return '10px';
            return '11px';
          }

          return '11px';
        });
    }"
  )
  
  w
}


# Plot: chord diagram of gene–drug interactions (DGIdb); options for approved-only, top-K per gene, and link scaling.
chord_gene_drug_overview <- function(
    df_raw,
    approved_only      = FALSE,       # keep only approved drugs
    top_k_per_gene     = 0,           # 0 = keep all; otherwise keep top-K drugs by score
    min_score          = NULL,        # optional score filter
    link_width_range   = c(1, 6),     # rescale widths
    show_drug_labels   = FALSE,       # TRUE = show tiny drug labels; FALSE = hide
    drug_label_cex     = 0.4,        # size for drug labels (if shown)
    gene_label_cex     = 0.65,        # size for gene labels
    gene_color         = "#C1CBBE",   # sector color for genes
    approved_col       = "#76B3E5",   # link color (approved)
    other_col          = "#CCEDC5",   # link color (other)
    alpha_links        = 0.8,        # transparency for links
    gap_degree_genes   = 2,           # gap between gene sectors
    gap_degree_drugs   = 1,           # gap between drug sectors
    seed               = 1
){
  set.seed(seed)
  
  # --- 1) normalize + filters ---
  df <- df_raw %>%
    rename(score = interaction.score,
           approved = regulatory.approval) %>%
    mutate(
      gene = as.character(gene),
      drug = as.character(drug),
      approved = as.character(approved)
    ) %>%
    filter(!is.na(gene), !is.na(drug))
  
  if (approved_only) df <- df %>% filter(tolower(approved) %in% c("true","approved","yes","1"))
  if (!is.null(min_score)) df <- df %>% filter(score >= min_score)
  
  df <- df %>%
    group_by(gene, drug, approved) %>%
    summarize(score = max(score, na.rm = TRUE), .groups = "drop")
  
  if (top_k_per_gene > 0) {
    df <- df %>% group_by(gene) %>%
      slice_max(order_by = score, n = top_k_per_gene, with_ties = FALSE) %>%
      ungroup()
  }
  
  # --- 2) sets & ordering ---
  genes <- sort(unique(df$gene))
  drugs <- sort(unique(df$drug))
  
  # sector order: genes first (big block), then drugs
  sector_order <- c(genes, drugs)
  
  # sector colors: genes gray; drugs light gray (so links stand out)
  # sector colors: genes (fixed) + drugs by approval status
  drug_status <- df %>%
    mutate(approved_flag = tolower(approved) %in% c("true","approved","yes","1")) %>%
    distinct(drug, approved_flag)
  
  drug_sector_col <- setNames(
    ifelse(drug_status$approved_flag,
           scales::alpha(approved_col, 0.25),   # softer tone for approved
           scales::alpha(other_col,    0.25)),  # softer tone for other
    drug_status$drug
  )
  
  grid_col <- c(
    setNames(rep(gene_color, length(genes)), genes),
    drug_sector_col
  )
  
  
  # --- 3) link data & aesthetics ---
  df <- df %>% mutate(
    approved_flag = tolower(approved) %in% c("true","approved","yes","1"),
    link_col = ifelse(approved_flag,
                      scales::alpha(approved_col, alpha_links),
                      scales::alpha(other_col,    alpha_links))
  )
  
  links <- df %>% select(from = gene, to = drug, value = score, col = link_col)
  links$value <- scales::rescale(links$value, to = link_width_range,
                                 from = range(links$value, finite = TRUE))
  
  # --- 4) layout gaps: larger gap between the two blocks (genes vs drugs) ---
  gap_after <- c(
    rep(gap_degree_genes, length(genes) - 1),
    8,  # big gap between genes and drugs
    rep(gap_degree_drugs, length(drugs) - 1),
    8
  )
  
  circos.clear()
  circos.par(start.degree = 90, gap.after = gap_after, track.margin = c(0.004, 0.004))
  
  chordDiagram(
    x = links[, c("from","to","value")],
    order = sector_order,
    grid.col = grid_col,
    transparency = 0,                  # we control transparency via link colors
    col = links$col,                   # color per link
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = mm_h(5))
  )
  
  
  # --- 5) labels: large for genes, tiny/hidden for drugs ---
  circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    sector = CELL_META$sector.index
    xlim   = CELL_META$xlim
    ylim   = CELL_META$ylim
    is_gene = sector %in% genes
    
    if (is_gene) {
      cex_use <- gene_label_cex
      lab <- sector
    } else {
      if (!show_drug_labels) return()  # hide drug labels entirely
      cex_use <- drug_label_cex
      # shorten very long drug names
      lab <- if (nchar(sector) > 28) paste0(substr(sector, 1, 28), "…") else sector
    }
    
    circos.text(
      x = mean(xlim), y = ylim[1] + 0.1,
      labels = lab, facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.5), cex = cex_use
    )
  })
}

# Plot: single-disease lollipop chart for enrichment results (GO/KEGG); colored by −log10(FDR).
plot_go_lollipop_one <- function(df, disease,
                                 x_metric = c("GeneRatio","RichFactor"),
                                 top_n = 10,
                                 low_col = "#67a9cf",
                                 high_col = "#d7301f",
                                 point_stroke = 0.6,
                                 show_title = TRUE) {   # <- new parameter
  
  x_metric <- match.arg(x_metric)
  
  d <- df[df$Disease == disease, , drop = FALSE]
  if (!nrow(d)) stop("No rows for disease: ", disease)
  
  if (x_metric == "GeneRatio") {
    num <- as.numeric(sub("/.*", "", d$GeneRatio))
    den <- as.numeric(sub(".*/", "", d$GeneRatio))
    d$Xvalue <- num / den
    x_lab <- "GeneRatio"
  } else {
    d$Xvalue <- as.numeric(d$RichFactor)
    x_lab <- "RichFactor"
  }
  
  d$logFDR <- -log10(as.numeric(d$p.adjust))
  
  d_top <- d |>
    dplyr::slice_min(p.adjust, n = top_n) |>
    dplyr::arrange(Xvalue)
  
  d_top$Desc_ord <- reorder(d_top$Description, d_top$Xvalue)
  
  p <- ggplot2::ggplot(d_top) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = Xvalue, y = Desc_ord, yend = Desc_ord),
      colour = "grey55", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = Xvalue, y = Desc_ord, size = Count, fill = logFDR),
      shape = 21, colour = "black", stroke = point_stroke
    ) +
    ggplot2::scale_fill_gradient(
      name = expression(-log[10]~FDR),
      low  = low_col,
      high = high_col
    ) +
    ggplot2::labs(title = if (show_title) disease else NULL,  # <- conditional title
                  x = x_lab, y = NULL, size = "Count") +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.6),
      panel.grid.major.x = ggplot2::element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "plain"),
      legend.position = "right"
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, width = 30),
      expand = ggplot2::expansion(mult = c(0.03, 0.04))
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1),
      size = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0.005, 0.05)),
      breaks = function(x) { b <- pretty(x); b[b != 0] }
    )
  
  p
}

# Plot: volcano (log2FC vs −log10(FDR)) with thresholds and optional labeled genes.
make_volcano <- function(
    tt, title = "", lfc = 1, fdr = 0.05, top_labels = 15,
    up_col = "#018571", down_col = "#a6611a", ns_col = "grey70",
    label_genes = NULL,              # genes to highlight/label
    label_only_if_sig = FALSE,       # keep only if passes FDR/lfc
    bg_alpha_updown = 0.15,          # background (Up/Down) lighter
    bg_alpha_ns = 0.10,              # background (NS) very light
    sel_alpha = 0.90,                # selected points more saturated
    show_legend = TRUE               # show/hide color legend
) {
  df <- tt %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::mutate(
      neglog10FDR = -log10(adj.P.Val),
      sig = dplyr::case_when(
        adj.P.Val < fdr & logFC >=  lfc ~ "Up",
        adj.P.Val < fdr & logFC <= -lfc ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # genes to highlight
  label_df <- tibble::tibble()
  if (!is.null(label_genes)) {
    present <- intersect(unique(label_genes), df$Gene)
    label_df <- dplyr::filter(df, Gene %in% present, abs(logFC) >= lfc)
    if (label_only_if_sig) label_df <- dplyr::filter(label_df, sig != "NS")
  } else {
    label_df <- df %>%
      dplyr::filter(sig != "NS") %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = top_labels)
  }
  
  legend_pos <- if (show_legend) "right" else "none"
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = neglog10FDR, color = sig)) +
    ggplot2::geom_point(ggplot2::aes(alpha = sig), size = 1.3, show.legend = show_legend) +
    ggplot2::scale_color_manual(values = c(Down = down_col, Up = up_col, NS = ns_col)) +
    ggplot2::scale_alpha_manual(values = c(Down = bg_alpha_updown, Up = bg_alpha_updown, NS = bg_alpha_ns),
                                guide = "none") +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = 2) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = 2)
  
  if (nrow(label_df) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = label_df,
        inherit.aes = FALSE,
        ggplot2::aes(x = logFC, y = neglog10FDR, color = sig),
        size = 1.8, alpha = sel_alpha, show.legend = FALSE
      ) +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = logFC, y = neglog10FDR, label = Gene),
        color = "black", fontface = "bold",
        size = 3, max.overlaps = Inf, segment.alpha = 0.5
      )
  }
  
  p +
    ggplot2::labs(title = title, x = "Biological Variation (log2 Fold-Change)", y = "Statistical Significance (-log10(FDR)", 
                  color = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = legend_pos)
}

# Helper: build consistent edge keys from a graph (keeps direction if directed).
edge_keys <- function(g) {
  if (is.null(V(g)$name)) {
    stop("All graphs must have vertex names (V(g)$name).")
  }
  el <- igraph::as_edgelist(g, names = TRUE)
  if (nrow(el) == 0) return(character())
  if (igraph::is_directed(g)) {
    apply(el, 1, function(x) paste(x, collapse = "|"))  # keep direction
  } else {
    apply(el, 1, function(x) paste(sort(x), collapse = "|"))  # undirected
  }
}

# Helper: classify edges of a target graph by overlap with the other two (e.g., "A only", "Shared all").
edge_classes_for_target <- function(gA, gG, gO, target = c("A","G","O")) {
  target <- match.arg(target)
  # Collect key sets
  sets <- list(A = edge_keys(gA), G = edge_keys(gG), O = edge_keys(gO))
  gT <- switch(target, A = gA, G = gG, O = gO)
  kT <- edge_keys(gT)
  inA <- kT %in% sets$A
  inG <- kT %in% sets$G
  inO <- kT %in% sets$O
  
  cls <- switch(
    target,
    "A" = ifelse(inG & inO, "Shared all",
                 ifelse(inG, "Shared A-G",
                        ifelse(inO, "Shared A-O", "A only"))),
    "G" = ifelse(inA & inO, "Shared all",
                 ifelse(inA, "Shared G-A",
                        ifelse(inO, "Shared G-O", "G only"))),
    "O" = ifelse(inA & inG, "Shared all",
                 ifelse(inA, "Shared O-A",
                        ifelse(inG, "Shared O-G", "O only")))
  )
  
  # Consistent order for legend
  lvls <- switch(
    target,
    A = c("A only", "Shared A-G", "Shared A-O", "Shared all"),
    G = c("G only", "Shared G-A", "Shared G-O", "Shared all"),
    O = c("O only", "Shared O-A", "Shared O-G", "Shared all")
  )
  factor(cls, levels = lvls)
}

# Plot: edge-colored network for a chosen target graph (A/G/O) highlighting overlap classes.
plot_overlap_edges <- function(gA, gG, gO,
                               target = c("A","G","O"),
                               title = NULL,
                               vertex.size = 0.3,
                               vertex.label = NA,
                               vertex.color = "grey85",
                               layout = NULL,
                               edge_palette = NULL,   # named vector with the 4 classes
                               edge_width_map = NULL, # named vector widths per class
                               seed = 123,
                               fr_niter = 1000,
                               ...) {
  target <- match.arg(target)
  gT <- switch(target, A = gA, G = gG, O = gO)
  
  # Warn if directedness mismatches
  if (is.directed(gA) != is.directed(gG) || is.directed(gA) != is.directed(gO)) {
    warning("Graphs have mixed directedness; edge matching may be inconsistent.")
  }
  
  # Layout
  if (is.null(layout)) {
    set.seed(seed)
    layout <- igraph::layout_with_fr(gT, niter = fr_niter, grid = "nogrid")
  }
  
  # Classes (factor with target-specific levels)
  classes <- edge_classes_for_target(gA, gG, gO, target)
  
  # Default palettes keyed to EXACT class labels
  if (is.null(edge_palette)) {
    edge_palette <- switch(
      target,
      A = c("A only"="#DEA58C", "Shared A-G"="#E6E670", "Shared A-O"="#8CC5DE", "Shared all"="#D63A3A"),
      G = c("G only"="#A58CDE", "Shared G-A"="#E6E670", "Shared G-O"="#E670AB", "Shared all"="#D63A3A"),
      O = c("O only"="#8CDEA5", "Shared O-A"="#8CC5DE", "Shared O-G"="#E670AB", "Shared all"="#D63A3A")
    )
  }
  
  if (is.null(edge_width_map)) {
    edge_width_map <- setNames(c(3, 3, 3, 3), names(edge_palette)) # tweak if you want
  }
  
  # Map classes to aesthetics
  edge_cols <- unname(edge_palette[as.character(classes)])
  edge_w    <- unname(edge_width_map[as.character(classes)])
  
  # Title style (lighter)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(font.main = 1, cex.main = 1.0)
  
  graphics::plot(
    gT,
    layout = layout,
    vertex.size  = vertex.size,
    vertex.label = vertex.label,
    vertex.color = vertex.color,
    edge.width   = edge_w,
    edge.color   = edge_cols,
    main = title,
    ...
  )
  
  cls_tab <- table(classes)
  all_classes <- names(edge_palette)
  counts <- setNames(rep(0L, length(all_classes)), all_classes)
  counts[names(cls_tab)] <- as.integer(cls_tab)
  
  invisible(list(classes = classes, counts = counts, palette = edge_palette))
}

