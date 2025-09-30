# ========================================================
# Script: 05_Gene_Genes.R
# Purpose:
#   Compare and interpret gene sets consistently selected by multiple SLR-based
#   models in the AGvsO task (Oligo vs Astro+GBM):
#     - GeNeIV-derived models (best Manhattan / best Euclidean),
#     - Baseline Elastic Net,
#     - Twiner (angular-profile penalty).
#
#   What this script produces:
#     • PR AUC comparison across models (boxplot+jitter of bootstrap runs)
#     • Per-gene selection counts across models (side-by-side bars)
#     • 4-way Venn diagram of consistently selected gene sets (≥75% bootstraps)
#     • Expression heatmap across diseases for the union of genes
#     • Binary selection heatmap (gene × model) for a filtered gene set
#     • Drug–Gene exploration: chord diagram and Sankey diagrams (top-5, approved-only)
# ========================================================

set.seed(138)  # Reproducibility anchor for any stochastic steps (e.g., bootstrap order, plotting jitter)

# ---------------------------
# 1) Load required data, best models and functions
# ---------------------------
load("outputs/00_preprocessing.RData")
load("outputs/03_SLR_GeNeIV_multiDistance_AGvsO.RData")
load("outputs/04_SLR_twiner_AGvsO.RData")
source("scripts/utils/functions.R")
source("scripts/utils/visualization.R")

# ---------------------------
# 2) Global comparison across models
# ---------------------------
plot_test_auc_boxplots(list("Manhattan" = best_manhattan_SLR_model$linear,
                            "Euclidean" = best_euclidean_SLR_model$log,
                            "Elastic Net" = best_elastic_SLR_model$elastic_net,
                            "Twiner" = twiner_SLR$twiner))

# twiner genes
best_twiner_genes_75 <- twiner_genes_75

# Union of all gene sets across models
selected_genes <- unique(c(best_manhattan_genes_75, best_euclidean_genes_75, best_elastic_genes_75, 
                           best_twiner_genes_75))

# Per-gene selection frequency across models
plot_selection_counts_multi(
  models_list = list(
    Euclidean = best_euclidean_SLR_model$log,
    Manhattan = best_manhattan_SLR_model$linear,
    Elastic = best_elastic_SLR_model$elastic_net,
    Twiner = twiner_SLR$twiner
  ),
  genes = selected_genes,
  order = "asc",
  title = "Genes Selected Counts Across Models"
)

# Venn diagram of consistently selected sets (Elastic, Twiner, Manhattan, Euclidean)
plot_venn_consistently_selected_4(
  selected_genes_list = list(
    best_elastic_genes_75,
    best_twiner_genes_75,
    best_manhattan_genes_75,
    best_euclidean_genes_75
  ),
  category_names = c("Elastic Net","Twiner", "Manhattan", "Euclidean"),
  plot_title = "Consistently Selected Genes - >75 Bootstraps"
)

################## DEGS
# Differential Expression (AG vs Oligo) and functional interpretation.
#############################

load("data/glioma-raw-RNAseq-2021-classification.RData")

## ---------------------- Inputs ------------------------
## data.frames with RAW counts

data_raw_list <- list(
  Astrocytoma       = astro_raw_RNA,
  Glioblastoma      = gbm_raw_RNA,
  Oligodendroglioma = oligo_raw_RNA
)


# Make gene names unique 
to_matrix <- function(df) {
  m <- t(as.matrix(df))                   # genes x samples
  storage.mode(m) <- "numeric"
  if (anyDuplicated(rownames(m))) rownames(m) <- make.unique(rownames(m))
  m
}
raw_counts_list <- lapply(data_raw_list, to_matrix)

# Keep only genes present in all subtypes to ensure a fair comparison
common_genes <- Reduce(intersect, lapply(raw_counts_list, rownames))
raw_counts_list <- lapply(raw_counts_list, \(m) m[common_genes, , drop = FALSE])

# Concatenate all samples across subtypes into one matrix
raw_counts <- do.call(cbind, raw_counts_list)     # genes x samples


# ===================== Grouping AG =====================

sample_info <- data.frame(
  sample_id = colnames(raw_counts),
  subtype   = factor(rep(names(raw_counts_list), sapply(raw_counts_list, ncol)),
                     levels = c("Astrocytoma","Glioblastoma","Oligodendroglioma"))
)

# Collapse Astro + GBM into a single group
sample_info$group <- ifelse(sample_info$subtype == "Oligodendroglioma", "Oligo", "AG")
sample_info$group <- factor(sample_info$group, levels = c("AG","Oligo"))

# Two-level design matrix (no intercept) for AG and Oligo
design_AG <- model.matrix(~ 0 + group, data = sample_info)
colnames(design_AG) <- levels(sample_info$group)  

#edgeR preprocessing 
y <- DGEList(counts = raw_counts)

# Filter out lowly-expressed genes
keep <- filterByExpr(y, design = design_AG)
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalization
y <- calcNormFactors(y, method = "TMM")

# ===================== voom transform ==========================
# 'voom' models the mean–variance relationship and returns precision weights for limma
v <- voom(y, design_AG, plot = TRUE)

#fit linear model
fit_ag <- lmFit(v, design_AG)
fit_ag <- contrasts.fit(fit_ag, makeContrasts(AG_vs_Oligo = AG - Oligo, levels = design_AG))
fit_ag <- eBayes(fit_ag, robust = TRUE)   # robust=TRUE mitigates outlier influence

# ===================== Results =================================
res_AG_Oligo <- topTable(fit_ag, coef = "AG_vs_Oligo", number = Inf, sort.by = "P")

# Volcano plot 
AGvsO_volcano <- make_volcano(res_AG_Oligo, title = NULL, lfc = 1,
                              show_legend = FALSE,
                              sel_alpha = 2)
AGvsO_volcano

# Subset DE results to the union of model-selected genes
res_AG_oligo_sel <- res_AG_Oligo[rownames(res_AG_Oligo) %in% selected_genes, , drop = FALSE]

# Vector of selected gene symbols (FDR < 0.05) to label on volcano
selected_genes_pdaj <- rownames(res_AG_oligo_sel)[res_AG_oligo_sel$adj.P.Val < 0.05]

# Volcano highlighting only genes above |log2FC| ≥ 1.5 and FDR < 0.05 (labeled)
AGvsO_volcano_selected <- make_volcano(res_AG_Oligo, title = NULL, lfc = 1.5,
                                       label_genes = selected_genes_pdaj,
                                       show_legend = FALSE,
                                       sel_alpha = 2)

AGvsO_volcano_selected

# Order selected genes by effect size 
res_AG_oligo_sel <- res_AG_oligo_sel[order(res_AG_oligo_sel$logFC, decreasing = TRUE), , drop = FALSE]
deg_gene_order <- rownames(res_AG_oligo_sel)

# Append significance stars to each gene
res_AG_oligo_sel <- add_deg_stars(res_AG_oligo_sel)

#----------------------------------------------------------------------
# Enrichment Analysis
# Build a tidy DE table with SYMBOL column (requires dplyr/tibble/magrittr in session)
DEGs_AGvsO <- res_AG_Oligo %>%
  rownames_to_column("SYMBOL") %>%
  mutate(SYMBOL = as.character(SYMBOL))

# --- 1) Define significant, up and down -------------------------------------
alpha <- 0.05
lfc_up   <- 1.5   
lfc_down <- -1.5

sig_tbl <- DEGs_AGvsO %>% filter(adj.P.Val < alpha)

# Extract unique SYMBOLs for Up and Down sets
genes_up   <- sig_tbl %>% filter(logFC >  lfc_up)   %>% pull(SYMBOL) %>% unique()
genes_down <- sig_tbl %>% filter(logFC <  lfc_down) %>% pull(SYMBOL) %>% unique()

# --- 2) Build disease_gene_sets for ORA -------------------------------------
disease_sets <- list(
  Up   = genes_up,
  Down = genes_down
)

# --- 4) Perform GO ORA (BP) using your helper ---------------------------------------
DEGs_set <- c("Up","Down")

# Color palettes per direction for lollipop plots (low->high)
pal <- list(
  Up = c("#9CDE8C", "#8CDECE"),
  Down      = c("#8C9CDE", "#CE8CDE")
)

# Run GO BP enrichment (helper returns a named list per direction)
ego_disease_BP <- enrich_go_by_disease(disease_sets, ont = "BP", p_adj_cut = 0.05)

# Plot top-10 terms by GeneRatio for each direction (no titles for panel assembly later)
for (d in DEGs_set) {
  cols <- pal[[d]]
  p <- plot_go_lollipop_one(
    ego_disease_BP, d, x_metric = "GeneRatio", top_n = 10,
    low_col = cols[1], high_col = cols[2], show_title = FALSE)
  print(p)  
} 

# Repeat for GO Molecular Function
ego_disease_MF <- enrich_go_by_disease(disease_sets, ont = "MF", p_adj_cut = 0.05)

for (d in DEGs_set) {
  cols <- pal[[d]]
  p <- plot_go_lollipop_one(
    ego_disease_MF, d, x_metric = "GeneRatio", top_n = 10,
    low_col = cols[1], high_col = cols[2], show_title = FALSE)
  print(p)  
} 

# KEGG pathway ORA 
kegg_disease <- enrich_kegg_by_disease(disease_sets, p_adj_cut = 0.05)

for (d in DEGs_set) {
  cols <- pal[[d]]
  p <- plot_go_lollipop_one(
    kegg_disease, d, x_metric = "GeneRatio", top_n = 10,
    low_col = cols[1], high_col = cols[2], show_title = FALSE)
  print(p)  # mostra no device atual
}


################## GSEA
# Ranked enrichment using t from DEGs

DEGs_rank_df <- DEGs_AGvsO %>%
  transmute(SYMBOL = as.character(SYMBOL),
            score  = if ("t" %in% names(.)) t else sign(logFC) * (-log10(P.Value + 1e-300))) %>%
  filter(!is.na(SYMBOL), !is.na(score)) %>%
  group_by(SYMBOL) %>%                      # collapse duplicates if present
  summarise(score = max(score), .groups="drop") %>%
  arrange(desc(score))

# Named numeric vector for clusterProfiler's GSEA functions
geneList_SYM <- setNames(DEGs_rank_df$score, DEGs_rank_df$SYMBOL)

# GO:BP GSEA 
gsea_BP <- gseGO(
  geneList       = geneList_SYM,
  OrgDb          = org.Hs.eg.db,
  keyType        = "SYMBOL",     # use gene symbols directly
  ont            = "BP",         # "BP", "MF" or "CC"
  minGSSize      = 10,
  maxGSSize      = 500,
  pAdjustMethod  = "BH",
  pvalueCutoff   = 0.05,
  verbose        = FALSE
)

# Filter 
gsea_BP_res <- as.data.frame(gsea_BP) 
gsea_BP_filt <- gsea_BP_res %>%
  filter(!is.na(p.adjust), p.adjust <= 0.05,
         !is.na(NES), is.finite(NES),
         setSize >= 50)

BP_top_up <- gsea_BP_filt %>%
  filter(NES > 0) %>%
  arrange(desc(NES), p.adjust, desc(setSize)) %>%
  slice_head(n = 10) %>%
  mutate(direction = "AG_up")

BP_top_down <- gsea_BP_filt %>%
  filter(NES < 0) %>%
  arrange(NES, p.adjust, desc(setSize)) %>%
  slice_head(n = 10) %>%
  mutate(direction = "Oligo_up")

# GO:MF GSEA 
gsea_MF <- gseGO(
  geneList       = geneList_SYM,
  OrgDb          = org.Hs.eg.db,
  keyType        = "SYMBOL",     # usa SYMBOL diretamente
  ont            = "MF",         # "BP", "MF" ou "CC"
  minGSSize      = 10,
  maxGSSize      = 500,
  pAdjustMethod  = "BH",
  pvalueCutoff   = 0.05,
  verbose        = FALSE
)

gsea_MF_res <- as.data.frame(gsea_MF) 
gsea_MF_filt <- gsea_MF_res %>%
  filter(!is.na(p.adjust), p.adjust <= 0.05,
         !is.na(NES), is.finite(NES),
         setSize >= 50)

MF_top_up <- gsea_MF_filt %>%
  filter(NES > 0) %>%
  arrange(desc(NES), p.adjust, desc(setSize)) %>%
  slice_head(n = 10) %>%
  mutate(direction = "AG_up")

MF_top_down <- gsea_MF_filt %>%
  filter(NES < 0) %>%
  arrange(NES, p.adjust, desc(setSize)) %>%
  slice_head(n = 10) %>%
  mutate(direction = "Oligo_up")

################## KEGG
# KEGG GSEA requires ENTREZ IDs:

sym2ent <- bitr(names(geneList_SYM),
                fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

geneList_ENTREZ <- geneList_SYM[sym2ent$SYMBOL]
names(geneList_ENTREZ) <- sym2ent$ENTREZID
geneList_ENTREZ <- sort(geneList_ENTREZ, decreasing = TRUE)  

# 2) GSEA KEGG (Homo sapiens)
set.seed(42)
gsea_KEGG <- gseKEGG(
  geneList      = geneList_ENTREZ,
  organism      = "hsa",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

# Same filtering 
gsea_KEGG_res <- as.data.frame(gsea_KEGG) 
gsea_KEGG_filt <- gsea_KEGG_res %>%
  filter(!is.na(p.adjust), p.adjust <= 0.05,
         !is.na(NES), is.finite(NES),
         setSize >= 50)

KEGG_top_up <- gsea_KEGG_filt %>%
  filter(NES > 0) %>%
  arrange(desc(NES), p.adjust, desc(setSize)) %>%
  slice_head(n = 10) %>%
  mutate(direction = "AG_up")

KEGG_top_down <- gsea_KEGG_filt %>%
  filter(NES < 0) %>%
  arrange(NES, p.adjust, desc(setSize)) %>%
  slice_head(n = 10) %>%
  mutate(direction = "Oligo_up")


# ---------------------------
# 3) Gene filtering for downstream visualization
# ---------------------------

heatmap_matrix <- build_heatmap_matrix(astro, gbm, oligo,
                                       genes = selected_genes,
                                       zscore_rows = TRUE)


plot_gene_heatmap(heatmap_matrix,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  show_rownames = TRUE,
                  use_brewer = TRUE,
                  brewer_name = "BrBG",
                  brewer_direction = 1,
                  show_disease_legend = FALSE,
                  res_star = res_AG_oligo_sel)

# Binary selection (gene × model) heatmap f
plot_selection_heatmap(
  gene_list = selected_genes,
  model_gene_sets = list(
    Elastic = best_elastic_genes_75,
    Manhattan = best_manhattan_genes_75,
    Euclidean = best_euclidean_genes_75,
    Twiner = best_twiner_genes_75),
  show_legend = FALSE
  
)

# Persist the union list for downstream drug–gene exploration
save_genes_txt(selected_genes, "outputs/AGvsO_selected_genes.txt")



# ========================================================
# 5) Save workspace (AGvsO – Selected Genes)
# ========================================================
save.image(file = "results/05_gene_selection_AGvsO.RData")
