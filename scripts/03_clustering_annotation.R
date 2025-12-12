# ==============================================================================
# Single-Cell RNA-Seq Analysis: Clustering and Annotation Pipeline
# ==============================================================================
# This script performs clustering, UMAP, cell type annotation using SingleR,
# marker identification, and differential expression analysis
#
# Input: Preprocessed Seurat object (from 01_preprocessing.R)
# Output: Annotated Seurat object with cell type labels and marker genes
# ==============================================================================

# Load required libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(viridis)
library(gridExtra)
library(stringr)
library(celldex)
library(SingleR)
library(scater)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# ==============================================================================
# Load Preprocessed Data
# ==============================================================================
cat("Loading preprocessed Seurat object...\n")
object <- readRDS(file = 'preprocessed_with_ADT_seurat.rds')

# ==============================================================================
# Clustering and UMAP
# ==============================================================================
cat("Building neighbor graph...\n")
object <- FindNeighbors(object, dims = 1:35)

cat("Finding clusters (resolution = 0.85)...\n")
object <- FindClusters(object, resolution = 0.85)

cat("Running UMAP...\n")
object <- RunUMAP(object, dims = 1:35)

Idents(object) <- object@meta.data$seurat_clusters

# Visualize clusters
cluster_plot <- DimPlot(object, split.by = "orig.ident", label = TRUE) + 
  ggtitle("res 0.85")
print(cluster_plot)

# ==============================================================================
# Quality Control and Marker Visualization
# ==============================================================================
cat("Generating QC plots...\n")
qc_features <- c("nCount_RNA", "nFeature_RNA", "percent_mt", 
                 "CD3D", "CD4", "CD8A", "IL2RA", "FOXP3", "IL7R", "SELL", 
                 "CCR7", "CD27", "TCF7", "PPBP", "MKI67", "LYZ", "S100A8", "S100A9",
                 "CD45RO-ADT", "CD45RA-ADT")

qc_plot <- FeaturePlot(object, features = qc_features, ncol = 5, 
                       cols = c("lightgrey", "red"))
print(qc_plot)
qc_vlnplot <- VlnPlot(object, features = qc_features, pt.size=0)
print(qc_vlnplot)
# ==============================================================================
# Filtering Low-Quality Clusters
# ==============================================================================
# Remove clusters identified as non-T cells, platelets, or low quality
# Based on marker expression analysis
cat("Filtering low-quality clusters...\n")
Idents(object) <- object@meta.data$seurat_clusters

# Clusters to remove (adjust based on your data)
clusters_to_remove <- c("12", "15", "17", "21", "24") 
# Note: This cluster list may be different due to the inherent stochastic nature of these algorithms.
# Here, we removed clusters due to low nCount_RNA and nFeature_RNA (12,17,24), high percent_MT (15), or high PPBP expression (21). 

object <- subset(object, idents = clusters_to_remove, invert = TRUE)

# Verify filtering
DimPlot(object, split.by = "orig.ident", label = TRUE) + 
  ggtitle("After filtering")

# ==============================================================================
# Find Variable Features, Scale and Cluster after Filtering
# ==============================================================================
cat("Re-processing filtered data...\n")
object <- FindVariableFeatures(object = object, 
                               selection.method = 'vst', 
                               nfeatures = 2000, 
                               mean.cutoff = c(0.1, 8), 
                               dispersion.cutoff = c(1, Inf), 
                               assay = "RNA")
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
cat("Scoring cell cycle phases...\n")
object <- CellCycleScoring(object = object, 
                           s.features = s.genes, 
                           g2m.features = g2m.genes)

cat("Scaling data with cell cycle regression...\n")
object <- ScaleData(object = object, 
                    vars.to.regress = c("S.Score", "G2M.Score"), 
                    verbose = TRUE)
cat("Running PCA...\n")
object <- RunPCA(object, npcs = 50, assay = "RNA")

cat("Building neighbor graph...\n")
object <- FindNeighbors(object, dims = 1:35)

cat("Finding clusters (resolution = 0.5)...\n")
object <- FindClusters(object, resolution = 0.5)

cat("Running UMAP...\n")
object <- RunUMAP(object, dims = 1:35)

Idents(object) <- object@meta.data$seurat_clusters

# ==============================================================================
# SingleR Cell Type Annotation
# ==============================================================================
cat("Annotating cell types using SingleR...\n")

# Load reference datasets
ref_monaco <- MonacoImmuneData()
ref_hpca <- HumanPrimaryCellAtlasData()

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(object, assay = "RNA")

# Run SingleR prediction
pred <- SingleR(test = sce, 
                ref = list(ref_monaco, ref_hpca), 
                labels = list(ref_monaco$label.fine, ref_hpca$label.fine), 
                method = "cluster", 
                clusters = sce@colData@listData[["seurat_clusters"]])

# Match references
matched <- matchReferences(ref_monaco, ref_hpca, 
                           ref_monaco$label.fine, ref_hpca$label.fine)
pheatmap::pheatmap(matched, col = viridis::plasma(100), fontsize = 5)

# Export correlation matrix
write.csv(matched, file = "matchedSingleRaggr.csv", row.names = TRUE)

# Plot score heatmap
plotScoreHeatmap(pred, clusters = rownames(pred), 
                 show_colnames = TRUE, fontsize = 6)

# Export scores
scores_matrix <- as.matrix(pred@listData[["scores"]])
write.csv(scores_matrix, file = "SingleRscoreaggr.csv")

# Add labels to Seurat object
object[["SingleR.cluster.labels.bulk"]] <- pred$labels[match(object[[]][["seurat_clusters"]], rownames(pred))]

Idents(object) <- "SingleR.cluster.labels.bulk"
DimPlot(object, label = TRUE, repel = TRUE) + 
  labs(title = "SingleR Bulk Idents")

# ==============================================================================
# Marker Gene Identification
# ==============================================================================
cat("Finding marker genes...\n")

# Markers by SingleR labels
Idents(object) <- "SingleR.cluster.labels.bulk"
markers_singler <- FindAllMarkers(object, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
write.csv(markers_singler, file = "markers_singler.csv", row.names = TRUE)

# Markers by Seurat clusters
Idents(object) <- "seurat_clusters"
markers_cluster <- FindAllMarkers(object, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
write.csv(markers_cluster, file = "markers_cluster.csv", row.names = TRUE)

# ==============================================================================
# Manual Cell Type Annotation
# ==============================================================================
# Based on SingleR results and marker expression, manually annotate clusters
Idents(object) <- "seurat_clusters"
object <- RenameIdents(object, 
                       "0" = "CD4 memory", 
                       "1" = "CD4 memory", 
                       "2" = "CD8 memory/ TEMRA", 
                       "3" = "CD8 naïve/memory",
                       "4" = "CD4 memory", 
                       "5" = "CD8 naïve/memory", 
                       "6" = "CD4 naïve", 
                       "7" = "CD8 naïve",
                       "8" = "Treg", 
                       "9" = "gamma/delta", 
                       "10" = "CD8 naïve/memory", 
                       "11" = "T/monocyte doublets",
                       "12" = "Proliferating T", 
                       "13" = "CD3+; CD56+", 
                       "14" = "B cells")

object[["cell_classification_final"]] <- Idents(object)

# ==============================================================================
# Further Doublet Identification 
# ==============================================================================
cat("Identifying T cell-monocyte doublets...\n")
mono <- WhichCells(object, 
                   expression = LYZ > 0.5 & S100A8 > 0.5 | 
                               LYZ > 0.5 & S100A9 > 0.5)
mono <- as.data.frame(mono)
rownames(mono) <- mono$mono
mono$Doublet <- "Doublet"
mono$mono <- NULL
object <- AddMetaData(object, metadata = mono, col.name = "Doublet")

object$doublet_cluster <- paste(object$cell_classification_final, 
                                sep = "_", object$Doublet)
Idents(object) <- "doublet_cluster"

# Rename doublet clusters
object <- RenameIdents(object, 
                       "CD4 memory_Doublet" = "T/monocyte doublets", 
                       "CD4 memory_NA" = "CD4 memory", 
                       "CD4 naïve_Doublet" = "T/monocyte doublets", 
                       "CD4 naïve_NA" = "CD4 naïve",
                       "CD8 memory/ TEMRA_Doublet" = "T/monocyte doublets", 
                       "CD8 memory/ TEMRA_NA" = "CD8 memory/ TEMRA",
                       "CD8 naïve/memory_Doublet" = "T/monocyte doublets", 
                       "CD8 naïve/memory_NA" = "CD8 naïve/memory",
                       "CD8 naïve_Doublet" = "T/monocyte doublets", 
                       "CD8 naïve_NA" = "CD8 naïve",
                       "gamma/delta_Doublet" = "T/monocyte doublets", 
                       "gamma/delta_NA" = "gamma/delta",
                       "Proliferating T_Doublet" = "T/monocyte doublets", 
                       "Proliferating T_NA" = "Proliferating T",
                       "T/monocyte doublets_Doublet" = "T/monocyte doublets", 
                       "T/monocyte doublets_NA" = "T/monocyte doublets",
                       "Treg_Doublet" = "T/monocyte doublets", 
                       "Treg_NA" = "Treg",
                       "CD3+; CD56+_Doublet" = "T/monocyte doublets", 
                       "CD3+; CD56+_NA" = "CD3+; CD56+",
                       "B cells_NA" = "B cells")

object[["cell_classification_final"]] <- Idents(object)

# ==============================================================================
# Patient-Specific Subsets
# ==============================================================================
cat("Creating patient-specific subsets...\n")
object_pt1012 <- subset(x = object, subset = patient == "pt1012")
object_pt1013 <- subset(x = object, subset = patient == "pt1013")
object_pt1014 <- subset(x = object, subset = patient == "pt1014")
object_pt1019 <- subset(x = object, subset = patient == "pt1019")

# ==============================================================================
# Differential Expression Analysis by Timepoint
# ==============================================================================
cat("Performing differential expression analysis...\n")

object_test <- object
object_test <- SetIdent(object_test, value = "timepoint")

# Find markers between timepoints
object_w1_w2_markers <- FindMarkers(object_test, 
                                     ident.1 = "week2", 
                                     ident.2 = "week1")
object_w1_w4_markers <- FindMarkers(object_test, 
                                     ident.1 = "week4", 
                                     ident.2 = "week1")
object_w2_w4_markers <- FindMarkers(object_test, 
                                     ident.1 = "week4", 
                                     ident.2 = "week2")

write.csv(object_w1_w2_markers, 
          file = "object_w1_w2_markers.csv", 
          row.names = TRUE)


# ==============================================================================
# Final Marker Heatmap
# ==============================================================================
cat("Generating final marker heatmap...\n")
markers_final <- FindAllMarkers(object, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

top10 <- markers_final %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


DoHeatmap(object, features = top10$gene) + 
  theme(plot.margin = unit(c(4, 4, 0, 0), "cm"))

# ==============================================================================
# Save Final Object
# ==============================================================================
cat("Saving final annotated object...\n")
saveRDS(object, file = output_file)

