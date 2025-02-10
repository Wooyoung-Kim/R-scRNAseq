# Single-Cell RNA-seq Analysis Pipeline

This repository contains a set of R scripts used for single-cell RNA sequencing (scRNA-seq) analysis, including preprocessing, doublet detection, annotation, data integration, and cell-cell communication analysis using CellChat.

## Table of Contents

- [Installation](#installation)
- [Scripts Overview](#scripts-overview)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [License](#license)

## Installation

To use this pipeline, install the required R packages before running the scripts:

```r
install.packages(c("Seurat", "tidyverse", "SoupX", "DropletUtils", "DoubletFinder", "scCustomize", "qs", "CellChat", "patchwork", "NMF", "ggalluvial", "reticulate"))
```

Ensure that Python dependencies are installed via `reticulate`:

```r
library(reticulate)
reticulate::py_install(packages = 'umap-learn')
```

## Scripts Overview

### 1. `1_SoupX.R`

- **Purpose**: Preprocess raw 10X data, remove background RNA contamination using SoupX, and save cleaned matrices.
- **Input**: `filtered_feature_bc_matrix.h5` and `raw_feature_bc_matrix.h5`
- **Output**: Processed 10X format data with background correction.

### 2. `2_doublet.R`

- **Purpose**: Detect and remove doublets from scRNA-seq datasets using DoubletFinder.
- **Input**: Preprocessed 10X count matrices.
- **Output**: Filtered Seurat objects with doublets removed.

### 3. `3_Annotation.R`

- **Purpose**: Annotate cell types using SciBet and ScType based on reference datasets.
- **Input**: Seurat objects.
- **Output**: Annotated Seurat objects with predicted cell types.

### 4. `4_Integration.R`

- **Purpose**: Merge datasets from different tissues and conditions, normalize data, and perform batch correction using Harmony integration.
- **Input**: Annotated Seurat objects from multiple sources.
- **Output**: Integrated Seurat object with batch correction applied.

### 5. `5_Cellchat.R`

- **Purpose**: Perform cell-cell communication analysis using CellChat.
- **Input**: Integrated Seurat object.
- **Output**: CellChat objects for different conditions.

### 6. `6_CellChat_merge.R`

- **Purpose**: Merge and compare CellChat results across different experimental conditions.
- **Input**: CellChat objects from different conditions.
- **Output**: Visualization and statistical comparisons of ligand-receptor interactions.

## Usage

Run the scripts sequentially in the given order. For example:

```r
source("1_SoupX.R")
source("2_doublet.R")
source("3_Annotation.R")
source("4_Integration.R")
source("5_Cellchat.R")
source("6_CellChat_merge.R")
```

Make sure input files are correctly placed in expected directories before running the scripts.

## Dependencies

This pipeline relies on the following R packages:

- `Seurat`
- `SoupX`
- `DoubletFinder`
- `scCustomize`
- `qs`
- `CellChat`
- `reticulate`
- `patchwork`
- `NMF`
- `ggalluvial`

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

Save this file as `README.md` for use in your GitHub repository.
