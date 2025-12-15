# Single-Cell RNA-Seq Analysis: T-Cell Clustering and Annotation with ADT

This repository contains the analysis pipeline for single-cell RNA sequencing (scRNA-seq) data with antibody-derived tags (ADT) focusing on T-cell populations from a longitudinal IL-7 treatment study. The pipeline uses Seurat for preprocessing, clustering, and annotation, with SingleR for automated cell type identification and DSB for ADT denoising.

## Overview

This analysis workflow processes scRNA-seq data through the following stages:

0. **Get T cells**: Subset all PBMCs object to just T cells
1. **ADT Preprocessing**: Denoising antibody-derived tags using DSB normalization
2. **GEX Preprocessing**: Normalization, dimensionality reduction (PCA), and integration with ADT
3. **Clustering & Annotation**: Cell clustering, UMAP visualization, filtering, cell type annotation, and differential expression
4. **Visualization**: Generation of publication figures

## Repository Structure

```
nt_i7_hgg_washu/
├── README.md
├── data/
│   └── ADT_aggregation_csv.csv
└── scripts/
    ├── 00_Get_Tcells.R          # Placeholder for T-cell subsetting
    ├── 01_ADT_preprocessing.R   # ADT data processing and DSB normalization
    ├── 02_GEX_preprocessing.R   # GEX preprocessing and ADT integration
    ├── 03_clustering_annotation.R # Clustering, annotation, Filtering, and DE analysis
    └── 04_plots.R               # Figure generation for publication
```

## Prerequisites

### Required R Packages

Install the following R packages before running the scripts (package versions at bottom of README):
```
Seurat, ggplot2, cowplot, dplyr, Matrix, viridis, gridExtra, stringr, celldex, SingleR, scater, ggrepel, EnhancedVolcano, pheatmap, RColorBrewer, scRepertoire, ggpubr, janitor, readxl, purrr, dsb, patchwork
```

## Usage

### Step 0: Get T cells (`00_Get_Tcells.R`)

This script processes GEX data from aggregated CellRanger output:
- Reads aggregated GEX data from h5 file
- Performs log normalization and variable feature identification
- SingleR cell type annotation using MonacoImmuneData()
- Subsets object to T cells
- Update the input file path before running

**Outputs:**
- `tcell_seurat_object.rds`: T cells Seurat object

### Step 1: ADT Preprocessing (`01_ADT_preprocessing.R`)

This script processes antibody-derived tag (ADT) data from CellRanger outputs:
- Reads and aggregates ADT and GEX data across samples
- Performs DSB normalization to denoise ADT counts using isotype controls
- Generates quality control plots for library sizes
- Update the input file path before running

**Outputs:**
- `dsb_norm_prot_isowithquant.Rdata`: Denoised ADT matrix

### Step 2: GEX Preprocessing (`02_GEX_preprocessing.R`)

This script performs preprocessing on gene expression (GEX) data:
- Log normalization and variable feature identification
- Cell cycle scoring and regression
- Principal component analysis (PCA)
- Update the input file path before running

**Outputs:**
- `preprocessed_seurat.rds`: Preprocessed Seurat object

### Step 3: Clustering and Annotation (`03_clustering_annotation.R`)

This script performs clustering, celltype annotation, and differential expression:
- Neighbor graph construction and clustering with UMAP
- Quality control filtering (removes low-quality, non-T cell, and doublet clusters)
- SingleR cell type annotation using immune reference datasets
- Integration with denoised ADT data
- Manual annotation based on markers
- Marker gene identification
- Differential expression analysis between timepoints
- Update the input file path before running

**Outputs:**
- `annotated_seurat.rds`: Annotated Seurat object

### Step 4: Visualization (`04_plots.R`)

This script generates the plots that are in Figure 5 and Supplementary Figure 7.

## Minor figure labels note
The timepoints were updated from the code for plotting in the figures as follows-
- Week 1 in the code is Week 0 in the publication.
- Week 2 in the code is Week 1 in the publication.
- Week 4 in the code is Week 3 in the publication.
- Week 14 in the code is Week 12 in the publication.
- We plan on updating the code to reflect the correct timepoints as soon as possible, and apologize for any confusion in the meantime!

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
[29] dsb_0.3.0                   patchwork_1.1.2
```

