# ========================================================
# Script: 06_drugs_search_AGvsO.R
#    - Query DGIdb results for selected genes (precomputed TSV)
#    - Create: (i) chord diagram for all drug–gene links
#              (ii) Sankey (top-5 drugs per gene)
#              (iii) Sankey (approved drugs only)


# We queried the DGIdb (Drug–Gene Interaction database) to 
# identify existing drugs targeting the selected genes,
# and downloaded that same tsv file

# --- Load DGIdb results
source("scripts/utils/visualization.R")  

dbidb_data <- read.delim("outputs/AGvsO_drug_interaction_results.tsv", header = TRUE, sep = "\t")
dbidb_data <- dbidb_data %>%
  mutate(
    drug = sub("^CHEMBL:", "", drug)   # remove "CHEMBL:" prefix from drug identifier
  )
# Tip: If downstream needs CHEMBL IDs, keep an extra column before stripping the prefix.

# Count distinct approved drugs (as per DGIdb 'regulatory.approval' field)
n_aprovadas <- dbidb_data %>%
  dplyr::filter(regulatory.approval == "Approved") %>%
  dplyr::distinct(drug) %>%
  nrow()

# Vector of unique approved drug names (useful for reporting/tables)
drugs_aprovadas <- dbidb_data %>%
  dplyr::filter(regulatory.approval == "Approved") %>%
  dplyr::distinct(drug) %>%
  dplyr::pull(drug)



# --- Chord diagram (ALL drugs–genes). High-resolution PNG output ---
# Renders a global overview of gene–drug connectivity.
# Large canvas (3000×3000 @ 300 dpi) to keep labels readable.
png("images/AGvsO_Drugs_chord_big.png", width = 3000, height = 3000, res = 300)  # big canvas for better readability
chord_gene_drug_overview(dbidb_data, show_drug_labels = TRUE)
dev.off()

# --- Sankey diagram: top-5 drugs per gene (width scaled between 2 and 15) ---
# Produces an interactive HTML widget; width_range controls edge thickness scaling.
sankey_top5 <- sankey_gene_drug(dbidb_data, top_k_per_gene = 5, width_range = c(2, 15))
sankey_top5
htmlwidgets::saveWidget(sankey_top5, "images/AGvsO_sankey.html", selfcontained = TRUE)
browseURL("images/AGvsO_sankey.html")  # open locally after save

# --- Sankey diagram: approved drugs only (no top-k restriction) ---
# Filters to DGIdb entries flagged as 'Approved'; shows the full set per gene.
sankey_approved_drugs <- sankey_gene_drug(dbidb_data, approved_only = TRUE, top_k_per_gene = 0)
sankey_approved_drugs
htmlwidgets::saveWidget(sankey_approved_drugs, "images/AGvsO_sankey_approved.html", selfcontained = TRUE)
browseURL("images/AGvsO_sankey_approved.html")  
