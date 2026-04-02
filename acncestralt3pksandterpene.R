# ==============================================================================
# PIPELINE: EVOLUTIONARY DYNAMICS AND ASR OF BGCs IN VIRGIBACILLUS
# Description: This script processes genomic data, calibrates phylogenetic trees, 
# prepares presence/absence matrices, and performs Ancestral State Reconstruction 
# (ASR) for Biosynthetic Gene Clusters (Terpenes and T3PKS).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. INITIAL SETUP & LIBRARIES
# ------------------------------------------------------------------------------

# Load required packages (install if not present)
library(ape)
library(phytools)
library(tidyverse)

cat("✅ Setup complete. Libraries loaded.\n\n")

# ------------------------------------------------------------------------------
# 1. TREE CALIBRATION (ULTRAMETRIC CONVERSION)
# ------------------------------------------------------------------------------
cat("--- STEP 1: Calibrating Phylogenetic Tree ---\n")

# Read the maximum-likelihood species tree
tree <- read.tree("Species_Tree_OrthoFinder.nwk")

# Ensure the tree is dichotomous and rooted
if (!is.binary(tree)) {
  tree <- multi2di(tree)
}
if (!is.rooted(tree)) {
  tree <- midpoint.root(tree)
}

# Transform into an ultrametric tree using chronos (discrete model, max flexibility)
ultra_tree <- chronos(tree, lambda = 0, model = "discrete", 
                      control = chronos.control(iter.max = 50000, eval.max = 50000))

# CRITICAL FIX FOR DOWNSTREAM ANALYSIS (e.g., CAFE5 & phytools):
# Replace zero or negative branch lengths with a minimal positive constant
ultra_tree$edge.length[ultra_tree$edge.length <= 0] <- 1e-6

# Save the calibrated tree
write.tree(ultra_tree, file = "virgibacillus_ultrametric.nwk")

cat("✅ SUCCESS: Ultrametric tree saved as 'virgibacillus_ultrametric.nwk'\n\n")

# ------------------------------------------------------------------------------
# 2. BiG-SCAPE TO PRESENCE/ABSENCE MATRIX (CAFE5 FORMAT)
# ------------------------------------------------------------------------------
cat("--- STEP 2: Creating Presence/Absence Matrix ---\n")

# Define network file path
network_dir <- "bigscape_results/bigscape_output_final/network_files/2025-11-04_18-42-49_hybrids_glocal/mix/"
clustering_file <- paste0(network_dir, "mix_onlyvirgi_clustering.tsv")

# Read and format the clustering file
clusters <- read_tsv(clustering_file, show_col_types = FALSE) %>%
  rename(bgc = `#BGC Name`, family = `Family Number`) %>%
  mutate(gcf = str_extract(bgc, "GCF_[0-9]+\\.[0-9]+"))

# Pivot to CAFE5 format matrix
cafe_matrix <- clusters %>%
  select(gcf, family) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = gcf,
    values_from = value,
    values_fill = 0,
    values_fn = length
  ) %>%
  rename(`Family ID` = family) %>%
  mutate(Desc = "Virgibacillus_BGC") %>% 
  relocate(Desc, `Family ID`)

write_tsv(cafe_matrix, "cafe5_bgc_counts.tsv")
cat("✅ SUCCESS: BGC Matrix saved as 'cafe5_bgc_counts.tsv'\n\n")

# ------------------------------------------------------------------------------
# 3. TREE PRUNING (MATCHING TREE TO MATRIX)
# ------------------------------------------------------------------------------
cat("--- STEP 3: Pruning Tree to match Matrix ---\n")

matrix_names <- colnames(cafe_matrix)[-c(1, 2)]

# Clean tip labels to match GCF format
clean_names <- str_extract(ultra_tree$tip.label, "GCF_[0-9]+\\.[0-9]+")
ultra_tree$tip.label <- ifelse(is.na(clean_names), ultra_tree$tip.label, clean_names)

# Prune tree (remove outgroups or genomes not in the matrix)
tree_pruned <- keep.tip(ultra_tree, matrix_names)
write.tree(tree_pruned, "virgibacillus_ultrametric_pruned.nwk")

# Sanity check
if(length(setdiff(tree_pruned$tip.label, matrix_names)) == 0) {
  cat("✅ SUCCESS: Tree and Matrix names are 100% perfectly aligned.\n\n")
} else {
  cat("❌ WARNING: Name mismatch detected between tree and matrix.\n\n")
}

# ------------------------------------------------------------------------------
# 4. GENOMIC EXPANSION VS BIOSYNTHETIC DIVERSITY (CORRELATION)
# ------------------------------------------------------------------------------
cat("--- STEP 4: Correlation Analysis ---\n")
# NOTE: Make sure 'df_tamanhos' (genome sizes) and 'df_counts' (bgc counts) are loaded 
# in your environment before running this block.

# df_tamanhos$GCF_clean <- substr(df_tamanhos$GCF, 1, 13)
# df_counts$GCF_clean <- substr(df_counts$GCF, 1, 13)
# final_df <- merge(df_tamanhos, df_counts, by="GCF_clean")

# correlacao <- cor.test(final_df$Tamanho_Mb, final_df$Total_BGCs)
# p_val <- format(correlacao$p.value, digits = 4)
# r_val <- format(correlacao$estimate, digits = 3)

# plot_cor <- ggplot(final_df, aes(x=Tamanho_Mb, y=Total_BGCs)) +
#   geom_point(color="dodgerblue4", size=3, alpha=0.7) +
#   geom_smooth(method="lm", color="firebrick", fill="lightpink", se=TRUE) +
#   theme_minimal(base_size = 14) +
#   labs(
#     title = "Genomic Expansion vs. Biosynthetic Diversity",
#     subtitle = paste("Pearson R =", r_val, "| p-value =", p_val),
#     x = "Genome Size (Mb)",
#     y = "Total Number of BGCs (Terpenes + T3PKS)"
#   )

# ggsave("Correlation_Final_Virgibacillus.pdf", plot = plot_cor, width=8, height=6)
cat("✅ Script template for correlation ready (uncomment to execute).\n\n")

# ------------------------------------------------------------------------------
# 5. ANCESTRAL STATE RECONSTRUCTION (ASR) - TERPENES & T3PKS
# ------------------------------------------------------------------------------
cat("--- STEP 5: Ancestral State Reconstruction (ARD Model) ---\n")

# Reload clean tree and matrix
tree <- read.tree("virgibacillus_ultrametric_pruned.nwk")
bgc_matrix <- read.table("cafe5_bgc_counts.tsv", header=TRUE, sep="\t", check.names=FALSE)

# Prepare species mapping for plot labels
mapping <- read.csv("onlyvirgi_cleaned.csv", stringsAsFactors = FALSE)
mapping$gcf <- gsub(" ", "_", mapping$gcf) 
species_dict <- setNames(mapping$species, mapping$gcf)

# Create a dedicated plot tree with aligned species names [GCF \n Species]
tree_plot <- tree
tree_plot$tip.label <- sapply(tree$tip.label, function(gcf) {
  sp <- species_dict[gcf]
  if(!is.na(sp)) return(paste0(gcf, "\n[", sp, "]")) else return(gcf)
})

# --- TERPENES PIPELINE ---
terpene_ids <- c(22, 58, 180, 187, 246, 348, 405, 407, 440, 442, 478, 490, 501, 506, 620)
dir.create("ASR_Trees_Terpenes", showWarnings = FALSE)
col_terpenes <- c("grey85", "forestgreen") 

cat("Processing Terpene GCFs...\n")
for (id in terpene_ids) {
  row_data <- bgc_matrix[bgc_matrix[,2] == id, ]
  if(nrow(row_data) == 0) next 
  
  # Extract presence/absence (binary)
  core_vals <- as.numeric(unlist(row_data[1, 3:ncol(bgc_matrix)]))
  core_vals <- ifelse(is.na(core_vals) | core_vals == 0, 0, 1) 
  
  trait_vector <- as.factor(core_vals)
  names(trait_vector) <- colnames(bgc_matrix)[3:ncol(bgc_matrix)]
  trait_vector <- trait_vector[tree$tip.label] # Order alignment
  
  # Skip if no evolutionary variation
  if (length(levels(trait_vector)) < 2) next
  
  # ASR using Asymmetric Rates (ARD)
  fit_model <- fitMk(tree, trait_vector, model="ARD")
  anc_states <- ancr(fit_model)
  
  # Plot & Save
  pdf(paste0("ASR_Trees_Terpenes/Terpene_ID_", id, ".pdf"), width=7, height=7)
  plotTree(tree_plot, type="fan", fsize=0.3, lwd=0.8, offset=0.5)
  tiplabels(pie=to.matrix(trait_vector, c("0", "1")), piecol=col_terpenes, cex=0.25)
  nodelabels(pie=anc_states$ace, piecol=col_terpenes, cex=0.35)
  legend("topleft", legend="", title=paste("Terpene ID:", id), bty="n", cex=1.2)
  dev.off()
}

# --- T3PKS PIPELINE ---
t3pks_ids <- c(15, 29, 63, 70, 84, 150, 160, 176, 177, 178, 182, 408, 411, 426, 
               429, 437, 439, 441, 448, 467, 477, 480, 492, 502, 510, 511, 552, 581, 637)
dir.create("ASR_Trees_T3PKS", showWarnings = FALSE)
col_t3pks <- c("grey85", "darkorange")

cat("Processing T3PKS GCFs...\n")
for (id in t3pks_ids) {
  row_data <- bgc_matrix[bgc_matrix[,2] == id, ]
  if(nrow(row_data) == 0) next
  
  core_vals <- as.numeric(unlist(row_data[1, 3:ncol(bgc_matrix)]))
  core_vals <- ifelse(is.na(core_vals) | core_vals == 0, 0, 1)
  
  trait_vector <- as.factor(core_vals)
  names(trait_vector) <- colnames(bgc_matrix)[3:ncol(bgc_matrix)]
  trait_vector <- trait_vector[tree$tip.label] 
  
  if (length(levels(trait_vector)) < 2) next
  
  fit_model <- fitMk(tree, trait_vector, model="ARD")
  anc_states <- ancr(fit_model)
  
  pdf(paste0("ASR_Trees_T3PKS/T3PKS_ID_", id, ".pdf"), width=7, height=7)
  plotTree(tree_plot, type="fan", fsize=0.3, lwd=0.8, offset=0.5)
  tiplabels(pie=to.matrix(trait_vector, c("0", "1")), piecol=col_t3pks, cex=0.25)
  nodelabels(pie=anc_states$ace, piecol=col_t3pks, cex=0.35)
  legend("topleft", legend="", title=paste("T3PKS ID:", id), bty="n", cex=1.2)
  dev.off()
}

cat("🎉 DONE! All analyses completed and ASR plots generated successfully.\n")
