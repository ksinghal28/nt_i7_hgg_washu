# Single-Cell RNA-Seq Analysis: T-Cell Clustering and Annotation

This repository contains the analysis pipeline for single-cell RNA sequencing (scRNA-seq) data focusing on T-cell populations from a longitudinal IL-7 treatment study. The pipeline uses Seurat for preprocessing, clustering, and annotation, with SingleR for automated cell type identification.

## Overview

This analysis workflow processes scRNA-seq data through the following stages:

1. **Preprocessing**: Normalization, dimensionality reduction (PCA), and quality control
2. **Clustering & Annotation**: Cell clustering, UMAP visualization, and cell type annotation
3. **Visualization**: Generation of publication-quality figures

## Repository Structure

```
code_ocean_git/
├── README.md
├── 01_preprocessing.R          # Data normalization and PCA
├── 02_clustering_annotation.R  # Clustering, annotation, and DE analysis
└── 03_plots.R                 # Figure generation for publication
```

## Prerequisites

### Required R Packages

Install the following R packages before running the scripts (package version information at bottom of README):
```
Seurat, ggplot2, cowplot, dplyr, Matrix, viridis, gridExtra, stringr, celldex, SingleR, scater, ggrepel, EnhancedVolcano, pheatmap, RColorBrewer,scRepertoire,ggpubr,janitor,readxl,purrr
```

### Data Requirements

- **Input**: Seurat object (`.rds` file) containing single-cell expression data
- The object should contain:
  - Raw or filtered count matrix in the `RNA` assay
  - Cell metadata including patient IDs, timepoints, and sample information
  - Quality metrics (nCount_RNA, nFeature_RNA, percent_mt)

## Usage

### Step 1: Preprocessing (`01_preprocessing.R`)

This script performs:
- Log normalization
- Variable feature identification
- Cell cycle scoring and regression
- Principal component analysis (PCA)
- JackStraw analysis for significant PC identification
- Elbow plot generation

**Before running**, update the input file path:
```r
input_file <- "/path/to/your/tcell_subcluster.rds"
```

**Run the script:**
```r
source("01_preprocessing.R")
```

**Outputs:**
- Preprocessed Seurat object (`.rds` file)
- JackStraw plot (`.jpg`)
- Elbow plot (`.jpg`)
- Dimension heatmaps (displayed in R console)

### Step 2: Clustering and Annotation (`02_clustering_annotation.R`)

This script performs:
- Neighbor graph construction and clustering
- UMAP dimensionality reduction
- Quality control filtering 
  Note: This cluster list may be different due to the version differences. We removed clusters due to low nCount_RNA, nFeature_RNA, high percent_MT, or high PPBP expression based on violin plots and UMAPs. 
- SingleR cell type annotation
- Marker gene identification
- Doublet detection (identified a previously overlooked set of monocyte-T cell doublets, and labeled them as such)
- Differential expression analysis between timepoints
- Heatmap generation

**Before running**, update the input file path:
```r
input_file <- "/path/to/your/preprocessed_seurat_object.rds"
```

**Run the script:**
```r
source("02_clustering_annotation.R")
```

**Outputs:**
- Annotated Seurat object (`.rds` file)
- Marker gene lists (`.csv` files)
- SingleR annotation scores (`.csv` file)
- Differential expression results (`.csv` files)
- Volcano plots (displayed in R console)
- Marker gene heatmap (`.jpeg`)

### Step 3: Visualization (`03_plots.R`)

This script generates publication-quality figures including:
- UMAP plots colored by timepoint, patient, and cell type
- Feature plots for marker genes
- Proportion plots
- Clonal diversity analysis
- Additional custom visualizations

**Before running**, update the input file path:
```r
object <- readRDS(file = "/path/to/your/annotated_seurat_object.rds")
```

**Run the script:**
```r
source("03_plots.R")
```

**Outputs:**
- Publication figures (`.pdf` and `.png` files)
- Saved in the `FinalFigures/` directory (update path as needed)

## Key Parameters

### Preprocessing Parameters
- `nfeatures = 2000`: Number of variable features to identify
- `npcs = 50`: Number of principal components to compute
- `scale.factor = 10000`: Scaling factor for normalization

### Clustering Parameters
- `PC = 35`: Number of principal components for clustering
- `resolution = 0.85`: Initial clustering resolution (high for fine-grained clusters)
- `resolution = 0.5`: Final clustering resolution after filtering

### Filtering Parameters
The script removes clusters identified as:
- Non-T cells (based on CD3D expression)
- Platelets (PPBP expression)
- Low-quality cells
- Doublets (T cell-monocyte doublets based on LYZ, S100A8, S100A9)

**Note**: Adjust `clusters_to_remove` in `02_clustering_annotation.R` based on your specific data.

## Cell Type Annotations

The pipeline identifies the following T-cell populations:
- **CD4 memory**: Memory CD4+ T cells
- **CD4 naïve**: Naïve CD4+ T cells
- **CD8 memory/ TEMRA**: Memory and terminally differentiated CD8+ T cells
- **CD8 naïve/memory**: Transitional CD8+ T cell populations
- **CD8 naïve**: Naïve CD8+ T cells
- **Treg**: Regulatory T cells
- **gamma/delta**: γδ T cells
- **Proliferating T**: Actively dividing T cells
- **CD3+; CD56+**: NKT-like cells
- **T/monocyte doublets**: Cell doublets (filtered)
- **B cells**: B lymphocytes (contaminants, filtered)

## Differential Expression Analysis

The pipeline performs differential expression analysis comparing:
- Week 2 vs Week 1
- Week 4 vs Week 1
- Week 4 vs Week 2

## Minor figure labels note
The timepoints were updated from the code for plotting in the figures as follows-
- Week 1 in the code is Week 0 in the publication.
- Week 2 in the code is Week 1 in the publication.
- Week 4 in the code is Week 3 in the publication.
- Week 14 in the code is Week 12 in the publication.
- We plan on updating the code to reflect the correct timepoints as sooon as possible, and apologize for any confusion in the meantime!

For questions or issues, please contact Kartik Singhal (k.singhal@wustl.edu)

## Package versions

```
> sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Ventura 13.6.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-3          pheatmap_1.0.12             EnhancedVolcano_1.14.0      ggrepel_0.9.4              
 [5] scater_1.24.0               scuttle_1.8.4               SingleCellExperiment_1.20.1 SingleR_1.10.0             
 [9] celldex_1.6.0               SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
[13] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0        
[17] MatrixGenerics_1.10.0       matrixStats_1.0.0           stringr_1.5.0               gridExtra_2.3              
[21] viridis_0.6.4               viridisLite_0.4.2           Matrix_1.5-4.1              dplyr_1.1.4                
[25] cowplot_1.1.1               ggplot2_3.4.4               SeuratObject_4.1.3          Seurat_4.3.0.1
```
