# *************************************************************
# Post hoc analyses of genomic accipiter data
# *************************************************************
# Script to create a tanglegram from mitochondrial and genomic phylogenetic trees
# This script reads pre-computed trees, matches sample names, and creates a tanglegram
# Author: Martin Kapun with the help of GitHub Copilot (Anthropic - Claude Sonnet 4.5)
# last changed: February 8, 2026

# Load required libraries
library(ape)
library(phytools)
library(dendextend)
library(phangorn)

# Set working directory
setwd("")

# Read the name mapping file (note: as samples have been named in different formats, this file contains the key to align all names to the same format (e.g., A.g.gentilis01))
cat("Reading name mapping file...\n")
name_mapping <- read.csv("data/names_mastersheet.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)

# Clean up the data - remove rows with missing values in key columns
name_mapping <- name_mapping[!is.na(name_mapping$Final_ID) & name_mapping$Final_ID != "", ]

cat("Total samples in mapping file:", nrow(name_mapping), "\n")

# Read the phylogenetic trees
cat("Reading mitochondrial tree (NEXUS format)...\n")
tree_mito <- read.nexus("data/tree_COM_FinalUsed_haplo_HKY.nwk")

cat("Reading genomic tree (Newick format)...\n")
tree_genomic <- read.tree("data/FINAL_52inds_manuel.nexus.varsites.phy.treefile")

# Print tree information before renaming
cat("\nMitochondrial tree:\n")
cat("  Number of tips:", Ntip(tree_mito), "\n")
cat("  Sample tip labels:", paste(tree_mito$tip.label[1:min(10, length(tree_mito$tip.label))], collapse = ", "), "...\n")

cat("\nGenomic tree:\n")
cat("  Number of tips:", Ntip(tree_genomic), "\n")
cat("  Sample tip labels:", paste(tree_genomic$tip.label[1:min(10, length(tree_genomic$tip.label))], collapse = ", "), "...\n")

# Create mapping dictionaries from genomic and mitochondrial names to final names
# Only keep samples that have entries in the relevant columns
genomic_to_final <- setNames(
    name_mapping$Final_ID[!is.na(name_mapping$FINAL_63inds_cleaned.recode.vcf) & name_mapping$FINAL_63inds_cleaned.recode.vcf != ""],
    name_mapping$FINAL_63inds_cleaned.recode.vcf[!is.na(name_mapping$FINAL_63inds_cleaned.recode.vcf) & name_mapping$FINAL_63inds_cleaned.recode.vcf != ""]
)

mito_to_final <- setNames(
    name_mapping$Final_ID[!is.na(name_mapping$tree_COM_FinalUsed_haplo_HKY.nwk) & name_mapping$tree_COM_FinalUsed_haplo_HKY.nwk != ""],
    name_mapping$tree_COM_FinalUsed_haplo_HKY.nwk[!is.na(name_mapping$tree_COM_FinalUsed_haplo_HKY.nwk) & name_mapping$tree_COM_FinalUsed_haplo_HKY.nwk != ""]
)

cat("\nSamples with genomic data:", length(genomic_to_final), "\n")
cat("Samples with mitochondrial data:", length(mito_to_final), "\n")

# Rename tree tips to final names
cat("\nRenaming tree tips to final IDs...\n")

# For genomic tree
genomic_tips_to_keep <- tree_genomic$tip.label %in% names(genomic_to_final)
tree_genomic_renamed <- tree_genomic
tree_genomic_renamed$tip.label[genomic_tips_to_keep] <- genomic_to_final[tree_genomic$tip.label[genomic_tips_to_keep]]

# For mitochondrial tree
mito_tips_to_keep <- tree_mito$tip.label %in% names(mito_to_final)
tree_mito_renamed <- tree_mito
tree_mito_renamed$tip.label[mito_tips_to_keep] <- mito_to_final[tree_mito$tip.label[mito_tips_to_keep]]

# Find common samples (using final IDs)
common_taxa <- intersect(tree_genomic_renamed$tip.label, tree_mito_renamed$tip.label)
cat("\nNumber of samples in both datasets:", length(common_taxa), "\n")
cat("Common samples:", paste(sort(common_taxa), collapse = ", "), "\n")

# Prune trees to only include common taxa
cat("\nPruning trees to common taxa...\n")
tree_genomic_pruned <- keep.tip(tree_genomic_renamed, common_taxa)
tree_mito_pruned <- keep.tip(tree_mito_renamed, common_taxa)

# Root both trees by A.nisus
cat("Rooting trees by A.nisus...\n")
if ("A.nisus" %in% tree_genomic_pruned$tip.label) {
    tree_genomic_pruned <- root(tree_genomic_pruned, outgroup = "A.nisus", resolve.root = TRUE)
    cat("  Genomic tree rooted by A.nisus\n")
} else {
    cat("  Warning: A.nisus not found in genomic tree, skipping rooting\n")
}

if ("A.nisus" %in% tree_mito_pruned$tip.label) {
    tree_mito_pruned <- root(tree_mito_pruned, outgroup = "A.nisus", resolve.root = TRUE)
    cat("  Mitochondrial tree rooted by A.nisus\n")
} else {
    cat("  Warning: A.nisus not found in mitochondrial tree, skipping rooting\n")
}

# Ladderize trees to improve visual layout (optional but helps with rooted trees)
tree_genomic_pruned <- ladderize(tree_genomic_pruned)
tree_mito_pruned <- ladderize(tree_mito_pruned)

# Convert rooted phylogenetic trees to dendrograms
cat("Converting rooted trees to dendrograms...\n")
# Make trees ultrametric using compute.brlen method
tree_genomic_ultra <- compute.brlen(tree_genomic_pruned, method = "Grafen")
tree_mito_ultra <- compute.brlen(tree_mito_pruned, method = "Grafen")
# Convert to dendrograms
dend_genomic <- as.dendrogram(as.hclust.phylo(tree_genomic_ultra))
dend_mito <- as.dendrogram(as.hclust.phylo(tree_mito_ultra))

# Create dendlist object
dend_list <- dendlist(dend_genomic, dend_mito)

# Untangle dendrograms to minimize crossing lines
cat("Rotating branches to minimize line overlaps...\n")
dend_list <- untangle(dend_list, method = "step2side")
dend_genomic <- dend_list[[1]]
dend_mito <- dend_list[[2]]

# Calculate entanglement (measure of how tangled the trees are)
cat("Calculating entanglement...\n")
entanglement_value <- tryCatch({
  entanglement(dend_list)
}, error = function(e) {
  cat("Warning: Could not calculate entanglement due to label issues\n")
  return(NA)
})
if(!is.na(entanglement_value)) {
  cat("Entanglement:", round(entanglement_value, 3), "(0=perfect alignment, 1=complete disorder)\n")
} else {
  cat("Entanglement: Not calculated\n")
}

# Create and save tanglegram as PNG
cat("\nCreating tanglegram plot...\n")
png("tanglegram.png", width = 5000, height = 2600, res = 150)
par(mar = c(6, 0, 2, 0), xaxt = "n")
tanglegram(
    dend_genomic, dend_mito,
    highlight_distinct_edges = TRUE,
    common_subtrees_color_branches = TRUE,
    main = "",
    main_left = "Genomic (ddRAD SNPs)",
    main_right = "Mitochondrial",
    cex_main_left = 5.0,
    cex_main_right = 5.0,
    lab.cex = 3.5,
    margin_inner = 38,
    lwd = 1.2,
    lty = 2, # dashed lines for connecting edges
    edge.lwd = 3.5,
    axes = FALSE
)
# Add legend - matching actual tanglegram colors
legend("bottom",
    legend = c("Common subtrees", "Distinct edges"),
    col = c("#ED0000FF", "#868686FF"),
    lty = c(2, 2),
    lwd = 4,
    cex = 2.8,
    bty = "n",
    horiz = TRUE,
    xpd = TRUE
)
dev.off()
cat("Tanglegram saved as: tanglegram.png\n")

# Create and save tanglegram as PDF
pdf("tanglegram.pdf", width = 44, height = 26)
par(mar = c(6, 0, 2, 0), xaxt = "n")
tanglegram(
    dend_genomic, dend_mito,
    highlight_distinct_edges = TRUE,
    common_subtrees_color_branches = TRUE,
    main = "",
    main_left = "Genomic (ddRAD SNPs)",
    main_right = "Mitochondrial",
    cex_main_left = 5.0,
    cex_main_right = 5.0,
    lab.cex = 3.5,
    margin_inner = 38,
    lwd = 1.2,
    lty = 2, # dashed lines for connecting edges
    edge.lwd = 3.5,
    axes = FALSE
)
# Add legend - matching actual tanglegram colors
legend("bottom",
    legend = c("Common subtrees", "Distinct edges"),
    col = c("#ED0000FF", "#868686FF"),
    lty = c(2, 2),
    lwd = 4,
    cex = 2.8,
    bty = "n",
    horiz = TRUE,
    xpd = TRUE
)
dev.off()
cat("Tanglegram saved as: tanglegram.pdf\n")
cat("Tanglegram saved as: tanglegram.pdf\n")

# Calculate and report Robinson-Foulds distance
cat("\nCalculating Robinson-Foulds distance...\n")
rf_dist <- RF.dist(tree_genomic_pruned, tree_mito_pruned)
cat("Robinson-Foulds distance:", rf_dist, "\n")
cat("(Lower values indicate more similar tree topologies)\n")

# Save statistics to text file
cat("\nSaving statistics to file...\n")
writeLines(c(
    "Robinson-Foulds Distance",
    "========================",
    paste("Value:", rf_dist),
    "",
    "Interpretation: Lower values indicate more similar tree topologies.",
    paste("For", length(common_taxa), "taxa, maximum possible RF =", 2 * (length(common_taxa) - 3)),
    paste("Observed RF represents", round(100 * rf_dist / (2 * (length(common_taxa) - 3)), 1), "% topological dissimilarity"),
    "",
    paste("Number of common taxa:", length(common_taxa)),
    "",
    "Additional Statistics:",
    paste("Entanglement:", round(entanglement_value, 4)),
    "(0 = perfect alignment, 1 = complete disorder)"
), "robinson_foulds_distance.txt")
cat("Robinson-Foulds distance saved to: robinson_foulds_distance.txt\n")

cat("\nDone! Check tanglegram.png and tanglegram.pdf\n")
