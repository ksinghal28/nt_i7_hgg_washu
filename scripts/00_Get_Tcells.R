# ==============================================================================
# Single-Cell RNA-Seq Analysis: T-Cell Subsetting Pipeline
# ==============================================================================
# This script subsets T-cell populations from a larger single-cell RNA-seq dataset
#
# Input: Aggregated Cell Ranger output
# Output: Seurat object containing only T cells (tcell_seurat.rds)
# ==============================================================================

library(Seurat)
library(SingleR)
library(celldex)

# ==============================================================================
# Load CellRanger outputs, Create Seurat Object, Process
# ==============================================================================

dat_GEX <- Read10X_h5('../data/filtered_feature_bc_matrix.h5')

dat_GEX <- CreateSeuratObject(counts = dat_GEX[["Gene Expression"]], min.cells = 10)

# Calculate mitochondrial and murine gene percentages 
dat_GEX[["percent.mt"]] <- PercentageFeatureSet(dat_GEX, pattern = "^GRCh38-MT-")
dat_GEX[["percent.mouse"]] <- PercentageFeatureSet(dat_GEX, pattern = "^mm10")

# Filter cells based on QC metrics
dat_GEX <- subset(dat_GEX, subset = percent.mouse < 10 & percent.mt < 10)

# Remove murine genes from the Gene Expression assay
DefaultAssay(dat_GEX) <- "RNA"
human_genes <- rownames(dat_GEX[["RNA"]])[grepl("^GRCh38-", rownames(dat_GEX[["RNA"]]))]
dat_GEX <- subset(dat_GEX, features = human_genes)

# Standard Seurat workflow on Gene Expression data
DefaultAssay(dat_GEX) <- "RNA"
dat_GEX <- NormalizeData(dat_GEX, normalization.method = "LogNormalize", scale.factor = 10000)
dat_GEX <- FindVariableFeatures(dat_GEX, selection.method = "vst", nfeatures = 2000)
dat_GEX <- ScaleData(dat_GEX)
dat_GEX <- RunPCA(dat_GEX, npcs = 50)
dat_GEX <- FindNeighbors(dat_GEX, dims = 1:10)
dat_GEX <- FindClusters(dat_GEX, resolution = 0.5)
dat_GEX <- RunUMAP(dat_GEX, dims = 1:10)

# Clean up gene names - remove GRCh38- prefix
counts <- GetAssayData(dat_GEX, slot = "counts")
rownames(counts) <- sub("^GRCh38-", "", rownames(counts))

# Recreate the RNA assay with cleaned gene names
dat_GEX[["RNA"]] <- CreateAssayObject(counts = counts)

# ==============================================================================
# Subset T Cells
# ==============================================================================

# Annotate with SingleR using Monaco Immune Dataset
monaco_ref <- celldex::MonacoImmuneData()
sce <- as.SingleCellExperiment(dat_GEX, assay = "RNA")
predictions <- SingleR(test = sce, ref = monaco_ref, labels = monaco_ref$label.main)

# Add SingleR predictions to Seurat object
dat_GEX$SingleR_labels <- predictions$labels

# Check what cell types were identified
table(dat_GEX$SingleR_labels)

# Subset T cells based on SingleR annotations
tcell_labels <- c("CD4+ T cells", "CD8+ T cells", "T cells")

dat_tcells <- subset(dat_GEX, subset = SingleR_labels %in% tcell_labels)

# ==============================================================================
# Save T-Cell Object
# ==============================================================================

saveRDS(dat_tcells, file = 'tcell_seurat_object.rds')



