# =============================================================================
# 1. Set Up the Python Environment using rminiconda and reticulate
# =============================================================================

# If you need to remove or reinstall Miniconda, uncomment and run the following commands:
# rminiconda::remove_miniconda(path = "/home/kwy7605/.rminiconda", name = "my_python")
# remotes::install_github("hafen/rminiconda", force = TRUE)
# rminiconda::install_miniconda(name = "my_python")
# py <- rminiconda::find_miniconda_python("my_python")

# Load the reticulate package to enable integration with Python.
library(reticulate)

# Specify the Python binary via reticulate. Adjust the path if necessary.
reticulate::use_python('/home/kwy7605/.rminiconda/my_python/bin/python3.12', required = TRUE)
# Alternatively, use the Python path discovered by rminiconda:
# reticulate::use_python(py, required = TRUE)

# Optionally, disable SSL verification for conda configuration if needed.
# system("conda config --set ssl_verify false")

# If required, install the Python package "umap-learn" using reticulate:
# reticulate::py_install(packages = 'umap-learn')

# =============================================================================
# 2. Load Required Libraries and Set Global Options
# =============================================================================

# Ensure the R library path is compatible with R version 4.3.
.libPaths("/usr/lib/R/library/4.3")

# Load libraries for parallel processing, data manipulation, visualization,
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

# (Optional) Enable parallel processing using 16 workers.
# future::plan("multisession", workers = 16)
options(future.rng.onMisuse = "ignore")

# Increase the maximum allowed size for global objects (~100GB) to accommodate large datasets.
options(future.globals.maxSize = 100000 * 1024^2)

# =============================================================================
# 3. Data Loading and Preprocessing
# =============================================================================

# Load the tissue data stored in qs format.
tissue <- qread("~/Exercise_Full/RDS/tissue.qs")

# Extract the metadata from the Seurat object and print unique cell types.
meta <- tissue@meta.data
print(unique(meta$celltype))

# Preprocess the Seurat object using SCTransform and prepare it for marker detection.
tissue <- PrepSCTFindMarkers(tissue)

# =============================================================================
# 4. Define the CellChat Analysis Function
# =============================================================================

# This function performs cell???cell communication analysis using CellChat.
# Arguments:
#   x: Seurat object containing single-cell transcriptomic data.
#   y: Metadata corresponding to the Seurat object.
#   z: Grouping variable (default is "ident" if not specified).
chat <- function(x, y, z = "ident") {
  
  # Create a CellChat object using the specified grouping variable and assay ("SCT").
  cellchat <- createCellChat(object = x, meta = y, group.by = z, assay = "SCT")
  
  # Load the mouse CellChat database for analysis of mouse data.
  CellChatDB <- CellChatDB.mouse
  
  # Display the database categories and inspect the structure of interactions.
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  
  # Subset the database to include only secreted signaling pathways.
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  # To use the complete database, comment out the above line and uncomment the following:
  # CellChatDB.use <- CellChatDB
  
  # Assign the chosen database to the CellChat object.
  cellchat@DB <- CellChatDB.use
  
  # Subset the expression data to include only signaling genes, reducing computation.
  cellchat <- subsetData(cellchat)
  
  # Identify overexpressed genes and corresponding interactions.
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # (Optional) Project gene expression data onto a protein???protein interaction (PPI) network.
  cellchat <- projectData(cellchat, PPI.mouse)
  
  # Compute the communication probability. Setting raw.use = FALSE uses the projected data.
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
  
  # Filter out communications that are supported by fewer than 10 cells.
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Calculate the communication probability at the signaling pathway level.
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate the cell???cell communication network.
  cellchat <- aggregateNet(cellchat)
  
  # Compute network centrality scores to pinpoint key signaling roles.
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  # Identify signaling groups based on functional similarity.
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netEmbedding(cellchat, type = "functional")
  cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)
  
  # Identify signaling groups based on structural similarity.
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural")
  cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE)
  
  # Return the processed CellChat object.
  return(cellchat)
}

# =============================================================================
# 5. Modify the Seurat Object and Subset by Condition
# =============================================================================

# Set the default assay to "SCT" and ensure it is of the correct class.
DefaultAssay(tissue) <- "SCT"
tissue[["SCT"]] <- as(object = tissue[["SCT"]], Class = "Assay")

# Print the unique experimental conditions from the metadata.
print(unique(tissue@meta.data$Condition))

# Subset the tissue object based on experimental conditions.
NCD    <- subset(tissue, subset = Condition == "Chow_SD")
NCD_TR <- subset(tissue, subset = Condition == "Chow_TR")
HFD    <- subset(tissue, subset = Condition == "HFD_SD")
HFD_TR <- subset(tissue, subset = Condition == "HFD_TR")

# Remove the full tissue object from memory to free up resources.
rm(tissue)

# =============================================================================
# 6. CellChat Analysis with Immune-Type Grouping
# =============================================================================

# Create a new metadata column 'immune_type' by concatenating 'tissue' and 'Type'.
NCD@meta.data$immune_type    <- paste0(NCD@meta.data$tissue, "_", NCD@meta.data$Type)
NCD_TR@meta.data$immune_type <- paste0(NCD_TR@meta.data$tissue, "_", NCD_TR@meta.data$Type)
HFD@meta.data$immune_type    <- paste0(HFD@meta.data$tissue, "_", HFD@meta.data$Type)
HFD_TR@meta.data$immune_type <- paste0(HFD_TR@meta.data$tissue, "_", HFD_TR@meta.data$Type)

# Extract updated metadata for each subset.
NCD_meta    <- NCD@meta.data
NCD_TR_meta <- NCD_TR@meta.data
HFD_meta    <- HFD@meta.data
HFD_TR_meta <- HFD_TR@meta.data

# Perform CellChat analysis using 'immune_type' as the grouping variable.
NCD    <- chat(NCD, NCD_meta, "immune_type")
NCD_TR <- chat(NCD_TR, NCD_TR_meta, "immune_type")
HFD    <- chat(HFD, HFD_meta, "immune_type")
HFD_TR <- chat(HFD_TR, HFD_TR_meta, "immune_type")

# Save the resulting CellChat objects for the immune-type analysis.
qsave(NCD,    "~/Exercise_Full/RDS/CellChat_NCD_Immune.qs")
qsave(NCD_TR, "~/Exercise_Full/RDS/CellChat_NCD_TR_Immune.qs")
qsave(HFD,    "~/Exercise_Full/RDS/CellChat_HFD_Immune.qs")
qsave(HFD_TR, "~/Exercise_Full/RDS/CellChat_HFD_TR_Immune.qs")

# =============================================================================
# 7. CellChat Analysis with Cell-Type Grouping
# =============================================================================

# Create a new metadata column 'celltype_2' by concatenating 'tissue' and 'celltype'.
NCD@meta.data$celltype_2    <- paste0(NCD@meta.data$tissue, "_", NCD@meta.data$celltype)
NCD_TR@meta.data$celltype_2 <- paste0(NCD_TR@meta.data$tissue, "_", NCD_TR@meta.data$celltype)
HFD@meta.data$celltype_2    <- paste0(HFD@meta.data$tissue, "_", HFD@meta.data$celltype)
HFD_TR@meta.data$celltype_2 <- paste0(HFD_TR@meta.data$tissue, "_", HFD_TR@meta.data$celltype)

# Update the metadata after adding the new 'celltype_2' column.
NCD_meta    <- NCD@meta.data
NCD_TR_meta <- NCD_TR@meta.data
HFD_meta    <- HFD@meta.data
HFD_TR_meta <- HFD_TR@meta.data

# Perform CellChat analysis using the default grouping variable ("ident").
NCD    <- chat(NCD, NCD_meta)
NCD_TR <- chat(NCD_TR, NCD_TR_meta)
HFD    <- chat(HFD, HFD_meta)
HFD_TR <- chat(HFD_TR, HFD_TR_meta)

# Save the resulting CellChat objects for cell-type analysis.
qsave(NCD,    "~/Exercise_Full/RDS/CellChat_NCD.qs")
qsave(NCD_TR, "~/Exercise_Full/RDS/CellChat_NCD_TR.qs")
qsave(HFD,    "~/Exercise_Full/RDS/CellChat_HFD.qs")
qsave(HFD_TR, "~/Exercise_Full/RDS/CellChat_HFD_TR.qs")

# =============================================================================
# 8. Communication Pattern and Signaling Role Analysis
# =============================================================================

# ---- NCD Dataset: Outgoing Communication Patterns ----
# Identify the optimal number of patterns.
selectK(NCD, pattern = "outgoing")
npatterns <- 2  # Define the number of outgoing communication patterns.
dev.off()      # Close any open graphics devices.
NCD <- identifyCommunicationPatterns(NCD, pattern = "outgoing",
                                     k = npatterns, width = 5, height = 9)

# ---- NCD Dataset: Incoming Communication Patterns ----
selectK(NCD, pattern = "incoming")
npatterns <- 3  # Define the number of incoming communication patterns.
dev.off()
NCD <- identifyCommunicationPatterns(NCD, pattern = "incoming",
                                     k = npatterns, width = 5, height = 9)

# Aggregate the communication network for the NCD dataset.
NCD <- aggregateNet(NCD)

# ---- HFD Dataset: Outgoing Communication Patterns ----
selectK(HFD, pattern = "outgoing")
npatterns <- 2
dev.off()
HFD <- identifyCommunicationPatterns(HFD, pattern = "outgoing",
                                     k = npatterns, width = 5, height = 9)

# ---- HFD Dataset: Incoming Communication Patterns ----
selectK(HFD, pattern = "incoming")
npatterns <- 2
dev.off()
HFD <- identifyCommunicationPatterns(HFD, pattern = "incoming",
                                     k = npatterns, width = 5, height = 9)

# Aggregate the communication network for the HFD dataset.
HFD <- aggregateNet(HFD)

# ---- Signaling Role Analysis ----
# Generate scatter plots for the signaling roles in the HFD and HFD_TR networks.
gg1 <- netAnalysis_signalingRole_scatter(HFD)     # HFD network analysis.
gg2 <- netAnalysis_signalingRole_scatter(HFD_TR)    # HFD_TR network analysis.
combined_plot <- gg1 + gg2                          # Combine plots using patchwork.

# ---- UMAP Visualization ----
# Visualize the UMAP embedding split by Condition.
# (Ensure that the Seurat object 'SFTS' exists and contains a reduction named "umap.harmony".)
DimPlot(SFTS, split.by = "Condition", reduction = "umap.harmony")

# =============================================================================
# References
# =============================================================================
# For details on the CellChat methodology:
# Jin, S., Guerrero-Juarez, C. F., Zhang, L., Chang, I., Ramos, R., Kuan, C. H.,
# Myung, P., Plikus, M. V., & Su??rez-Fari??as, M. (2021).
# "CellChat: Inference and analysis of cell???cell communication using single-cell transcriptomic data."
# Nature Communications, 12, 1088.
# DOI: https://doi.org/10.1038/s41467-020-20260-8
#
# GitHub repository for rminiconda (Python environment management):
# https://github.com/hafen/rminiconda
