# ==============================================================================
# Single-Cell RNA-Seq Analysis: Preprocessing Pipeline
# ==============================================================================
# This script performs normalization, dimensionality reduction, and 
# quality control for single-cell RNA-seq data using Seurat
#
# Input: Seurat object (RDS file)
# Output: Preprocessed Seurat object with PCA, JackStraw, and Elbow plots
# ==============================================================================

# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)

# ==============================================================================
# Load Data
# ==============================================================================
cat("Loading Seurat object...\n")
object <- readRDS(file = '/path/to/tcell seurat object.rds')

# ==============================================================================
# Normalization and Feature Selection
# ==============================================================================
cat("Normalizing data...\n")
object <- NormalizeData(object = object, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000, 
                        assay = "RNA")

cat("Finding variable features...\n")
object <- FindVariableFeatures(object = object, 
                                selection.method = 'vst', 
                                nfeatures = 2000, 
                                mean.cutoff = c(0.1, 8), 
                                dispersion.cutoff = c(1, Inf), 
                                assay = "RNA")

# ==============================================================================
# Cell Cycle Scoring and Regression
# ==============================================================================
cat("Scoring cell cycle phases...\n")
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

object <- CellCycleScoring(object = object, 
                            s.features = s.genes, 
                            g2m.features = g2m.genes)

cat("Scaling data with cell cycle regression...\n")
object <- ScaleData(object = object, 
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    verbose = TRUE)

# ==============================================================================
# Principal Component Analysis
# ==============================================================================
cat("Running PCA...\n")
object <- RunPCA(object, npcs = 50, assay = "RNA")

# ==============================================================================
# Dimension Heatmaps
# ==============================================================================
cat("Generating dimension heatmaps...\n")
DimHeatmap(object, dims = 1:12, balanced = TRUE, cells = 500)
DimHeatmap(object, dims = 13:24, balanced = TRUE, cells = 500)
DimHeatmap(object, dims = 25:36, balanced = TRUE, cells = 500)

# ==============================================================================
# Elbow Plot
# ==============================================================================
cat("Generating elbow plot...\n")
elbow <- ElbowPlot(object, ndims = 50)
print(elbow)

# ==============================================================================
# Add ADT object's data
# ==============================================================================
cat("Reading in and adding ADT to seurat object...\n")
load('dsb_norm_prot_isowithquant.Rdata')

# Find common cells
common_cells <- intersect(colnames(object), colnames(dsb_norm_prot_iso_quant))

# Subset ADT objects to common cells
adt_data <- dsb_norm_prot_iso_quant[, common_cells]

# Now this will work
object[["ADT_denoised_iso_quant"]] <- CreateAssayObject(counts = adt_data)

# ==============================================================================
# Save Preprocessed Object
# ==============================================================================
cat("Saving preprocessed object...\n")
saveRDS(object, file = 'preprocessed_with_ADT_seurat.rds')
