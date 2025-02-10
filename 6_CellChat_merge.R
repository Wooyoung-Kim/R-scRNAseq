# =============================================================================
# 1. Environment Setup and Library Loading
# =============================================================================
# (Optional) Remove and reinstall Miniconda (uncomment if needed)
# rminiconda::remove_miniconda(path = "/home/kwy7605/.rminiconda", name = "my_python")
# remotes::install_github("hafen/rminiconda", force = TRUE)
# rminiconda::install_miniconda(name = "my_python")
# py <- rminiconda::find_miniconda_python("my_python")

# Load the reticulate package to interface with Python.
library(reticulate)
# Specify the Python binary required for analysis.
reticulate::use_python('/home/kwy7605/.rminiconda/my_python/bin/python3.12', required = TRUE)
# Alternatively, use the python path found by rminiconda:
# reticulate::use_python(py, required = TRUE)

# (Optional) If necessary, disable SSL verification or install Python packages (e.g., umap-learn)
# system("conda config --set ssl_verify false")
# reticulate::py_install(packages = 'umap-learn')

# Set the library path for R packages to ensure compatibility with R version 4.3.
.libPaths("/usr/lib/R/library/4.3")

# Load required libraries for parallel computing, data manipulation, visualization,
# and cell???cell communication analysis.
library(future)
library(dplyr)
library(CellChat)
library(Seurat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(reticulate)
library(qs)

# =============================================================================
# 2. Load Previously Processed CellChat Objects
# =============================================================================
# Load the CellChat objects (stored in qs format) corresponding to different conditions.
NCD    <- qread("~/Exercise_Full/RDS/CellChat_NCD.qs")
NCD_TR <- qread("~/Exercise_Full/RDS/CellChat_NCD_TR.qs")
HFD    <- qread("~/Exercise_Full/RDS/CellChat_HFD.qs")
HFD_TR <- qread("~/Exercise_Full/RDS/CellChat_HFD_TR.qs")

# Update each CellChat object to the latest version (if any updates exist).
NCD    <- updateCellChat(NCD)
NCD_TR <- updateCellChat(NCD_TR)
HFD    <- updateCellChat(HFD)
HFD_TR <- updateCellChat(HFD_TR)

# =============================================================================
# 3. Define Common Cell Types and Subset the Data
# =============================================================================
# Extract unique cell-type identifiers ("celltype_2") from each dataset.
a <- NCD@meta$celltype_2   %>% unique()
b <- NCD_TR@meta$celltype_2 %>% unique()
c <- HFD@meta$celltype_2   %>% unique()
d <- HFD_TR@meta$celltype_2 %>% unique()

# Identify the common cell types across all datasets.
common_elements <- Reduce(intersect, list(a, b, c, d))

# Subset each CellChat object to include only the common cell types.
NCD    <- subsetCellChat(object = NCD,    idents.use = common_elements)
NCD_TR <- subsetCellChat(object = NCD_TR, idents.use = common_elements)
HFD    <- subsetCellChat(object = HFD,    idents.use = common_elements)
HFD_TR <- subsetCellChat(object = HFD_TR, idents.use = common_elements)

# =============================================================================
# 4. Merge Datasets for Comparative Analysis
# =============================================================================
# Create lists for merging datasets.
object.list_NH      <- list(NCD = NCD, HFD = HFD)
object.list_NH_TH   <- list(HFD = HFD, HFD_TR = HFD_TR)

# Merge CellChat objects; the 'add.names' argument retains the original dataset names.
cellchat_NH    <- mergeCellChat(object.list_NH,      add.names = names(object.list_NH))
cellchat_NH_TH <- mergeCellChat(object.list_NH_TH,   add.names = names(object.list_NH_TH))

# =============================================================================
# 5. Compare Global Interaction Patterns
# =============================================================================
# Compare overall cell???cell interactions between datasets.
gg1 <- compareInteractions(cellchat_NH,    show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat_NH_TH, show.legend = FALSE, group = c(1, 2))
# Display the side-by-side comparison plots.
gg1 + gg2

# =============================================================================
# 6. Heatmap Visualization of Communication Networks
# =============================================================================
# Define tissue-specific cell-type subsets based on metadata.
vWAT  <- cellchat_NH@meta %>% filter(tissue == "vWAT")  %>% pull(celltype_2) %>% unique()  # e.g., indexes 34:51
scWAT <- cellchat_NH@meta %>% filter(tissue == "scWAT") %>% pull(celltype_2) %>% unique()  # e.g., indexes 1:15
SkM   <- cellchat_NH@meta %>% filter(tissue == "SkM")   %>% pull(celltype_2) %>% unique()  # e.g., indexes 16:33
type  <- cellchat_NH@meta %>% pull(celltype_2) %>% unique()

# Generate heatmaps for comparing communication between conditions.
gg1 <- netVisual_heatmap(cellchat_NH, comparison = c(1,2))
# Crop the heatmap to focus on specific tissue comparisons (e.g., vWAT vs. scWAT)
gg1 <- gg1[34:51, 1:15]
gg2 <- netVisual_heatmap(cellchat_NH, comparison = c(1,2))
gg2 <- gg2[34:51, 1:15]
gg1 + gg2

# =============================================================================
# 7. Signaling Role and Network Similarity Analyses
# =============================================================================
# Compute the total number of interactions (links) in each dataset.
# Note: 'object.list' should refer to the merged list of CellChat objects.
num.link <- sapply(object.list_NH, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
# Determine the minimum and maximum number of links across datasets.
weight.MinMax <- c(min(num.link), max(num.link))

# Generate scatter plots of the signaling roles for each dataset,
# controlling the dot size using the computed weight range.
gg <- list()
for (i in 1:length(object.list_NH)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list_NH[[i]],
                                               title = names(object.list_NH)[i],
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

# Compute pairwise network similarity using the functional attributes.
cellchat <- computeNetSimilarityPairwise(cellchat_NH, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")  # Manifold learning
cellchat <- netClustering(cellchat, type = "functional")  # Clustering/classification
# Visualize the embedding in a 2D space.
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

# Repeat the analysis for the structural properties.
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural",
                            label.size = 3.5, top.label = TRUE)
# Zoom in on the pairwise structural embedding.
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
# Rank the similarity scores.
rankSimilarity(cellchat, type = "structural")

# =============================================================================
# 8. Ranking and Bubble Plot Visualizations
# =============================================================================
# Compare ranked networks between datasets (e.g., vWAT vs. SkM).
gg1 <- rankNet(cellchat_NH, mode = "comparison", stacked = TRUE, do.stat = TRUE,
               comparison = c(1, 2), sources.use = vWAT, targets.use = SkM[1])
gg2 <- rankNet(cellchat_NH_TH, mode = "comparison", stacked = TRUE, do.stat = TRUE,
               comparison = c(1, 2), sources.use = vWAT, targets.use = SkM[1])
gg1 + gg2

# Visualize cell???cell communication for selected signaling pathways via bubble plots.
netVisual_bubble(cellchat_NH, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(1, 2), angle.x = 45,
                 signaling = c("PERIOSTIN", "PLAU", "FASLG", "ANGPTL", "CXCL")) +
  theme(axis.text = element_text(size = 6))
netVisual_bubble(cellchat_NH_TH, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(1, 2), angle.x = 45,
                 signaling = c("PERIOSTIN", "PLAU", "FASLG", "ANGPTL", "CXCL")) +
  theme(axis.text = element_text(size = 6))

# Example: Bubble plot for a single signaling pathway ("PERIOSTIN").
netVisual_bubble(cellchat_NH, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(1, 2), angle.x = 45, signaling = "PERIOSTIN") +
  theme(axis.text = element_text(size = 6))
netVisual_bubble(cellchat_NH_TH, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(1, 2), angle.x = 45, signaling = "PERIOSTIN") +
  theme(axis.text = element_text(size = 6))

# Additional bubble plots for other pathways (e.g., "PLAU" and "FASLG")
netVisual_bubble(cellchat, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(1, 2), angle.x = 45, signaling = "PLAU")
netVisual_bubble(cellchat, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(3, 4), angle.x = 45, signaling = "PLAU")
netVisual_bubble(cellchat, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(1, 2), angle.x = 45, signaling = "FASLG")
netVisual_bubble(cellchat, sources.use = vWAT, targets.use = SkM[1],
                 comparison = c(3, 4), angle.x = 45, signaling = "FASLG")

# =============================================================================
# 9. Differential Expression Analysis and Ligand???Receptor Mapping
# =============================================================================
# Load the ComplexHeatmap package for advanced heatmap visualizations.
library(ComplexHeatmap)

# Define the positive datasets (datasets with positive fold change).
pos.dataset_HFD    <- "HFD"
pos.dataset_HFD_TR <- "HFD_TR"

# Assign feature names for storing differential expression results.
features.name_HFD    <- pos.dataset_HFD
features.name_HFD_TR <- pos.dataset_HFD_TR

# Perform differential expression analysis for cellchat_NH.
cellchat_NH <- identifyOverExpressedGenes(cellchat_NH,
                                          group.dataset = "datasets",
                                          pos.dataset = pos.dataset_HFD,
                                          features.name = features.name_HFD,
                                          only.pos = FALSE,
                                          thresh.pc = 0.1,
                                          thresh.fc = 0.1,
                                          thresh.p  = 1)
# Map the DEG results onto the inferred communications to manage ligand???receptor pairs.
net_HFD <- netMappingDEG(cellchat_NH, features.name = features.name_HFD)
# Extract ligand???receptor pairs with upregulated ligands in the dataset with higher expression.
net.up_HFD <- subsetCommunication(cellchat_NH, net = net_HFD, datasets = "HFD",
                                  ligand.logFC = 0.2, receptor.logFC = NULL)
# Extract pairs with downregulated ligands and receptors (i.e., upregulated in the alternative dataset).
net.down_HFD <- subsetCommunication(cellchat_NH, net = net_HFD, datasets = "NCD",
                                    ligand.logFC = -0.1, receptor.logFC = -0.1)

# Repeat differential expression analysis for the second merged object.
cellchat_NH_TH <- identifyOverExpressedGenes(cellchat_NH_TH,
                                             group.dataset = "datasets",
                                             pos.dataset = pos.dataset_HFD_TR,
                                             features.name = features.name_HFD_TR,
                                             only.pos = FALSE,
                                             thresh.pc = 0.1,
                                             thresh.fc = 0.1,
                                             thresh.p  = 1)
net_HFD_TR <- netMappingDEG(cellchat_NH_TH, features.name = features.name_HFD_TR)
net.up_HFD_TR <- subsetCommunication(cellchat_NH_TH, net = net_HFD_TR, datasets = "HFD_TR",
                                     ligand.logFC = 0.2, receptor.logFC = NULL)
net.down_HFD_TR <- subsetCommunication(cellchat_NH_TH, net = net_HFD_TR, datasets = "HFD",
                                       ligand.logFC = -0.1, receptor.logFC = -0.1)

# Extract gene subsets from the ligand???receptor pairs.
gene.up_HFD      <- extractGeneSubsetFromPair(net.up_HFD, cellchat_NH)
gene.down_HFD    <- extractGeneSubsetFromPair(net.down_HFD, cellchat_NH)
gene.up_HFD_TR   <- extractGeneSubsetFromPair(net.up_HFD_TR, cellchat_NH_TH)
gene.down_HFD_TR <- extractGeneSubsetFromPair(net.down_HFD_TR, cellchat_NH_TH)

# =============================================================================
# 10. Visualization of Differentially Regulated Ligand???Receptor Pairs
# =============================================================================
# Create a bubble plot for upregulated interactions.
pairLR.use.up_HFD <- net.up_HFD[, "interaction_name", drop = FALSE]
gg1 <- netVisual_bubble(cellchat_NH,
                        pairLR.use = pairLR.use.up_HFD,
                        sources.use = vWAT,
                        targets.use = SkM,
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = TRUE,
                        title.name = paste0("Up-regulated signaling in ", names(object.list_NH)[2]))

# Create a bubble plot for downregulated interactions.
pairLR.use.down_HFD <- net.down_HFD[, "interaction_name", drop = FALSE]
gg2 <- netVisual_bubble(cellchat_NH,
                        pairLR.use = pairLR.use.down_HFD,
                        sources.use = vWAT,
                        targets.use = SkM,
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = TRUE,
                        title.name = paste0("Down-regulated signaling in ", names(object.list_NH)[1]))
# Adjust text size in the plot.
gg1 + theme(axis.text = element_text(size = 5))

# =============================================================================
# 11. Chord Diagrams and Enrichment Score Computation
# =============================================================================
# Set up a 1x2 plotting layout.
par(mfrow = c(1,2), xpd = TRUE)
netVisual_chord_gene(object.list_NH[[2]], sources.use = vWAT, targets.use = SkM,
                     slot.name = 'net', net = net.up_HFD,
                     lab.cex = 0.8, small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ", names(object.list_NH)[2]))
netVisual_chord_gene(object.list_NH[[1]], sources.use = vWAT, targets.use = SkM,
                     slot.name = 'net', net = net.down_HFD,
                     lab.cex = 0.8, small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(object.list_NH)[2]))

# Compute enrichment scores for the up- and down-regulated networks.
computeEnrichmentScore(net.down_HFD, species = 'mouse')
computeEnrichmentScore(net.up_HFD,   species = 'mouse')
computeEnrichmentScore(net.down_HFD_TR, species = 'mouse')
computeEnrichmentScore(net.up_HFD_TR,   species = 'mouse')

# =============================================================================
# 12. Aggregated Network Visualization for Specific Pathways
# =============================================================================
# Specify signaling pathways of interest.
pathways.show <- c("PERIOSTIN")
# Get the maximum edge weight across datasets (for consistent scaling).
weight.max <- getMaxWeight(object.list_NH_TH, slot.name = c("netP"), attribute = pathways.show)

# Plot aggregated network visualizations in a circular layout.
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list_NH_TH)) {
  netVisual_aggregate(object.list_NH_TH[[i]],
                      signaling = pathways.show,
                      layout = "circle",
                      edge.weight.max = weight.max[1],
                      edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list_NH_TH)[i]),
                      sources.use = vWAT,
                      targets.use = SkM[1])
}

# Generate heatmaps for the selected pathway across datasets.
ht <- list()
for (i in 1:length(object.list_NH_TH)) {
  ht[[i]] <- netVisual_heatmap(object.list_NH_TH[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling", names(object.list_NH_TH)[i]),
                               font.size = 6)
}
ComplexHeatmap::draw(ht[[1]][34:51, 16:33] + ht[[2]][34:51, 16:33],
                     ht_gap = unit(0.5, "cm"))

# Repeat the heatmap visualization for the other merged object.
ht <- list()
for (i in 1:length(object.list_NH)) {
  ht[[i]] <- netVisual_heatmap(object.list_NH[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling", names(object.list_NH)[i]),
                               font.size = 6)
}
ComplexHeatmap::draw(ht[[1]][34:51, 16:33] + ht[[2]][34:51, 16:33],
                     ht_gap = unit(0.5, "cm"))

# =============================================================================
# 13. Gene Expression Visualization and Data Export
# =============================================================================
# Set factor levels for the dataset labels to ensure proper ordering in plots.
cellchat_NH@meta$datasets    <- factor(cellchat_NH@meta$datasets, levels = c("NCD", "HFD"))
cellchat_NH_TH@meta$datasets <- factor(cellchat_NH_TH@meta$datasets, levels = c("HFD", "HFD_TR"))

# Plot gene expression for the "PERIOSTIN" signaling pathway, split by datasets.
plotGeneExpression(cellchat_NH,    signaling = "PERIOSTIN", split.by = "datasets", colors.ggplot = TRUE)
plotGeneExpression(cellchat_NH_TH, signaling = "PERIOSTIN", split.by = "datasets", colors.ggplot = TRUE)

# Plot gene expression grouped by tissue.
plotGeneExpression(cellchat_NH,    signaling = "PERIOSTIN", split.by = "datasets", colors.ggplot = TRUE, group.by = "tissue")
plotGeneExpression(cellchat_NH_TH, signaling = "PERIOSTIN", split.by = "datasets", colors.ggplot = TRUE, group.by = "tissue")

# Export the differential network data to CSV files.
write.csv(net_HFD,    "~/Exercise_Full/RDS/HFD_Net.csv")
write.csv(net_HFD_TR, "~/Exercise_Full/RDS/HFD_TR_Net.csv")
