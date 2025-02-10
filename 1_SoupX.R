# Load necessary libraries for parallel processing, single-cell analysis,
# ambient RNA correction, droplet utilities, plotting, doublet detection, and dynamic reporting.
library(future)        # Enables parallel and distributed processing
options(future.globals.maxSize = 1024 * 1024 * 1024 * 20)  # Set max global object size to 20 GB

library(Seurat)        # For single-cell RNA-seq data processing and analysis (Stuart et al., 2019)
library(SoupX)         # For correction of ambient RNA contamination (Young & Behjati, 2020)
library(DropletUtils)  # For handling 10X Genomics droplet-based data
library(ggplot2)       # For data visualization
library(DoubletFinder) # For identification of doublets in scRNA-seq data (McGinnis et al., 2019)
library(knitr)         # For dynamic report generation

# Function to process 10X Genomics data and perform ambient RNA correction
process_10x_data <- function(input_path, output_path) {
  
  ## 1. Load Filtered and Raw Count Matrices
  # Read the filtered gene expression matrix from the HDF5 file (10X Genomics output)
  filt.matrix <- Read10X_h5(file.path(input_path, "outs/filtered_feature_bc_matrix.h5"))
  # Read the raw gene expression matrix (including potential ambient RNA contamination)
  raw.matrix <- Read10X_h5(file.path(input_path, "outs/raw_feature_bc_matrix.h5"))
  
  ## 2. Create a Seurat Object for Downstream Analysis
  # Instantiate a Seurat object using the filtered count matrix.
  seurat <- CreateSeuratObject(counts = filt.matrix)
  
  ## 3. Initialize SoupChannel for Ambient RNA Estimation
  # Create a SoupChannel object that takes both raw and filtered matrices.
  # This object will be used to estimate and correct for ambient RNA contamination.
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)
  
  ## 4. Preprocess the Data Using Seurat Workflow
  # Normalize and stabilize variance using SCTransform (Hafemeister & Satija, 2019)
  seurat <- SCTransform(seurat, verbose = FALSE)
  # Perform principal component analysis (PCA) for dimensionality reduction
  seurat <- RunPCA(seurat, verbose = FALSE)
  # Compute UMAP embeddings for visualization in two dimensions using the first 30 PCs
  seurat <- RunUMAP(seurat, dims = 1:30, verbose = FALSE)
  # Construct a k-nearest neighbor graph based on the reduced dimensions
  seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)
  # Identify clusters (i.e., putative cell populations) in the dataset
  seurat <- FindClusters(seurat, verbose = TRUE)
  
  ## 5. Extract Metadata and Dimensionality Reduction Coordinates
  # Retrieve cell-level metadata (including cluster assignments) from the Seurat object
  meta <- seurat@meta.data
  # Extract UMAP embeddings (used as a low-dimensional representation of the data)
  umap <- seurat@reductions$umap@cell.embeddings
  
  ## 6. Integrate Clustering and Dimensionality Information into SoupChannel
  # Update the SoupChannel object with cell cluster identities from Seurat
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  # Incorporate UMAP coordinates into the SoupChannel object for better ambient RNA estimation
  soup.channel <- setDR(soup.channel, umap)
  # Automatically estimate the level of ambient RNA contamination.
  # The parameter forceAccept = TRUE bypasses manual inspection.
  soup.channel <- autoEstCont(soup.channel, forceAccept = TRUE)
  
  ## 7. Adjust Count Matrix for Ambient RNA Contamination
  # Correct the count matrix for ambient RNA using SoupX, rounding values to integers
  adj.matrix <- adjustCounts(soup.channel, roundToInt = TRUE)
  
  ## 8. Write the Corrected Matrix in 10X Genomics Format
  # Save the contamination-corrected count matrix to the specified output directory
  DropletUtils:::write10xCounts(output_path, adj.matrix)
  
  ## (Optional) Print a message upon completion of processing
  # cat("Processing complete. Results saved to:", output_path, "\n")
}

## Process Multiple Datasets
# The following calls process individual datasets by specifying the input directory 
# (containing the 10X Genomics output) and the desired output directory for the corrected matrix.

# Process scWAT samples
process_10x_data("~/Exercise_Full/GSM5554964", "~/Exercise_Full/scWAT/GSM5554964_SoupX")
process_10x_data("~/Exercise_Full/GSM5554965", "~/Exercise_Full/scWAT/GSM5554965_SoupX")
process_10x_data("~/Exercise_Full/GSM5554966", "~/Exercise_Full/scWAT/GSM5554966_SoupX")
process_10x_data("~/Exercise_Full/GSM5554967", "~/Exercise_Full/scWAT/GSM5554967_SoupX")
process_10x_data("~/Exercise_Full/GSM5554968", "~/Exercise_Full/scWAT/GSM5554968_SoupX")
process_10x_data("~/Exercise_Full/GSM5554969", "~/Exercise_Full/scWAT/GSM5554969_SoupX")
process_10x_data("~/Exercise_Full/GSM5554970", "~/Exercise_Full/scWAT/GSM5554970_SoupX")
process_10x_data("~/Exercise_Full/GSM5554971", "~/Exercise_Full/scWAT/GSM5554971_SoupX")
process_10x_data("~/Exercise_Full/GSM5554972", "~/Exercise_Full/scWAT/GSM5554972_SoupX")
process_10x_data("~/Exercise_Full/GSM5554973", "~/Exercise_Full/scWAT/GSM5554973_SoupX")
process_10x_data("~/Exercise_Full/GSM5554974", "~/Exercise_Full/scWAT/GSM5554974_SoupX")
process_10x_data("~/Exercise_Full/GSM5554975", "~/Exercise_Full/scWAT/GSM5554975_SoupX")
process_10x_data("~/Exercise_Full/GSM5554976", "~/Exercise_Full/scWAT/GSM5554976_SoupX")
process_10x_data("~/Exercise_Full/GSM5554977", "~/Exercise_Full/scWAT/GSM5554977_SoupX")

# Process vWAT samples
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554978", "~/Exercise_Full/vWAT/GSM5554978_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554979", "~/Exercise_Full/vWAT/GSM5554979_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554980", "~/Exercise_Full/vWAT/GSM5554980_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554981", "~/Exercise_Full/vWAT/GSM5554981_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554982", "~/Exercise_Full/vWAT/GSM5554982_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554983", "~/Exercise_Full/vWAT/GSM5554983_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554984", "~/Exercise_Full/vWAT/GSM5554984_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554985", "~/Exercise_Full/vWAT/GSM5554985_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554986", "~/Exercise_Full/vWAT/GSM5554986_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554987", "~/Exercise_Full/vWAT/GSM5554987_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554988", "~/Exercise_Full/vWAT/GSM5554988_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554989", "~/Exercise_Full/vWAT/GSM5554989_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554990", "~/Exercise_Full/vWAT/GSM5554990_SoupX")
process_10x_data("/data2/kbsi/DATA/Wooyoung/vWAT/GSM5554991", "~/Exercise_Full/vWAT/GSM5554991_SoupX")

# Process SkM (skeletal muscle) samples
process_10x_data("/data1/kwy/SkM/GSM5554992", "~/Exercise_Full/SkM/GSM5554992_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5554993", "~/Exercise_Full/SkM/GSM5554993_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5554994", "~/Exercise_Full/SkM/GSM5554994_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5554995", "~/Exercise_Full/SkM/GSM5554995_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5554996", "~/Exercise_Full/SkM/GSM5554996_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5554997", "~/Exercise_Full/SkM/GSM5554997_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5554998", "~/Exercise_Full/SkM/GSM5554998_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5554999", "~/Exercise_Full/SkM/GSM5554999_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5555000", "~/Exercise_Full/SkM/GSM5555000_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5555001", "~/Exercise_Full/SkM/GSM5555001_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5555002", "~/Exercise_Full/SkM/GSM5555002_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5555003", "~/Exercise_Full/SkM/GSM5555003_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5555004", "~/Exercise_Full/SkM/GSM5555004_SoupX")
process_10x_data("/data1/kwy/SkM/GSM5555005", "~/Exercise_Full/SkM/GSM5555005_SoupX")
