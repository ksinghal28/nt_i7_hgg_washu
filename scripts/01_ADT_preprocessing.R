# ==============================================================================
# Single-Cell RNA-Seq Analysis: ADT Preprocessing Pipeline
# ==============================================================================
# This script processes antibody-derived tag (ADT) data from CellRanger outputs,
# aggregates samples, and performs DSB normalization for denoising
#
# Input: CellRanger output directories for each sample, aggregation CSV
# Output: Denoised ADT matrix (dsb_norm_prot_isowithquant.Rdata)
# ==============================================================================

library(Seurat)
library(dsb)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(Matrix)

# ==============================================================================
# Data Directories and Paths
# ==============================================================================
# Base path to CellRanger output directory - UPDATE THIS PATH
cellranger_base <- "path/to/cellranger outs"

# Aggregation CSV file - UPDATE THIS PATH
aggregation_csv <- "path/to/ADT_aggregation_csv.csv"

# Construct paths to per_sample_outs directories
scrna.dir <- list()
scrna.dir["pt1012_week1"] <- file.path(cellranger_base, "pt1012", "pt1012_week1", "outs", "per_sample_outs", "pt1012_week1")
scrna.dir["pt1012_week2"] <- file.path(cellranger_base, "pt1012", "pt1012_week2", "outs", "per_sample_outs", "pt1012_week2")
scrna.dir["pt1012_week4"] <- file.path(cellranger_base, "pt1012", "pt1012_week4", "outs", "per_sample_outs", "pt1012_week4")
scrna.dir["pt1012_week14"] <- file.path(cellranger_base, "pt1012", "pt1012_week14", "outs", "per_sample_outs", "pt1012_week14")
scrna.dir["pt1013_week1"] <- file.path(cellranger_base, "pt1013", "pt1013_week1", "outs", "per_sample_outs", "pt1013_week1")
scrna.dir["pt1013_week2"] <- file.path(cellranger_base, "pt1013", "pt1013_week2", "outs", "per_sample_outs", "pt1013_week2")
scrna.dir["pt1013_week4"] <- file.path(cellranger_base, "pt1013", "pt1013_week4", "outs", "per_sample_outs", "pt1013_week4")
scrna.dir["pt1013_week14"] <- file.path(cellranger_base, "pt1013", "pt1013_week14", "outs", "per_sample_outs", "pt1013_week14")
scrna.dir["pt1014_week1"] <- file.path(cellranger_base, "pt1014", "pt1014_week1", "outs", "per_sample_outs", "pt1014_week1")
scrna.dir["pt1014_week2"] <- file.path(cellranger_base, "pt1014", "pt1014_week2", "outs", "per_sample_outs", "pt1014_week2")
scrna.dir["pt1014_week4"] <- file.path(cellranger_base, "pt1014", "pt1014_week4", "outs", "per_sample_outs", "pt1014_week4")
scrna.dir["pt1019_week1"] <- file.path(cellranger_base, "pt1019", "pt1019_week1", "outs", "per_sample_outs", "pt1019_week1")
scrna.dir["pt1019_week4"] <- file.path(cellranger_base, "pt1019", "pt1019_week4", "outs", "per_sample_outs", "pt1019_week4")
scrna.dir["pt1019_week14"] <- file.path(cellranger_base, "pt1019", "pt1019_week14", "outs", "per_sample_outs", "pt1019_week14")

isotypes <- c("IgG1_ADT",  "IgG2a_ADT", "IgG2b_ADT")

# Construct paths to per_sample_outs directories
scrna.dir <- list()
scrna.dir["pt1012_week1"] <- file.path(cellranger_base, "pt1012", "pt1012_week1", "outs", "per_sample_outs", "pt1012_week1")
scrna.dir["pt1012_week2"] <- file.path(cellranger_base, "pt1012", "pt1012_week2", "outs", "per_sample_outs", "pt1012_week2")
scrna.dir["pt1012_week4"] <- file.path(cellranger_base, "pt1012", "pt1012_week4", "outs", "per_sample_outs", "pt1012_week4")
scrna.dir["pt1012_week14"] <- file.path(cellranger_base, "pt1012", "pt1012_week14", "outs", "per_sample_outs", "pt1012_week14")
scrna.dir["pt1013_week1"] <- file.path(cellranger_base, "pt1013", "pt1013_week1", "outs", "per_sample_outs", "pt1013_week1")
scrna.dir["pt1013_week2"] <- file.path(cellranger_base, "pt1013", "pt1013_week2", "outs", "per_sample_outs", "pt1013_week2")
scrna.dir["pt1013_week4"] <- file.path(cellranger_base, "pt1013", "pt1013_week4", "outs", "per_sample_outs", "pt1013_week4")
scrna.dir["pt1013_week14"] <- file.path(cellranger_base, "pt1013", "pt1013_week14", "outs", "per_sample_outs", "pt1013_week14")
scrna.dir["pt1014_week1"] <- file.path(cellranger_base, "pt1014", "pt1014_week1", "outs", "per_sample_outs", "pt1014_week1")
scrna.dir["pt1014_week2"] <- file.path(cellranger_base, "pt1014", "pt1014_week2", "outs", "per_sample_outs", "pt1014_week2")
scrna.dir["pt1014_week4"] <- file.path(cellranger_base, "pt1014", "pt1014_week4", "outs", "per_sample_outs", "pt1014_week4")
scrna.dir["pt1019_week1"] <- file.path(cellranger_base, "pt1019", "pt1019_week1", "outs", "per_sample_outs", "pt1019_week1")
scrna.dir["pt1019_week4"] <- file.path(cellranger_base, "pt1019", "pt1019_week4", "outs", "per_sample_outs", "pt1019_week4")
scrna.dir["pt1019_week14"] <- file.path(cellranger_base, "pt1019", "pt1019_week14", "outs", "per_sample_outs", "pt1019_week14")

isotypes <- c("IgG1_ADT",  "IgG2a_ADT", "IgG2b_ADT") 
samples <- list()
sample.data <- list()
sample.rawdata <- list()

# ==============================================================================
# Load Sample Data
# ==============================================================================
for (i in 1:length(scrna.dir)) {
  sample.data[[i]] <- Read10X(paste0(scrna.dir[[i]],"/count/sample_feature_bc_matrix/"))
  rownames(sample.data[[i]]$`Antibody Capture`) <- gsub("TotalSeqC", replacement = "ADT", rownames(sample.data[[i]]$`Antibody Capture`) )
  sample.data[[i]][["Sample"]] = scrna.dir[i]
}
names(sample.data) <- names(scrna.dir)

print("separate list for reading in raw files since multi cellranger output")
# Construct paths to main outs directories for raw data
scrna.dir.raw <- list()
scrna.dir.raw["pt1012_week1"] <- file.path(cellranger_base, "pt1012", "pt1012_week1", "outs")
scrna.dir.raw["pt1012_week2"] <- file.path(cellranger_base, "pt1012", "pt1012_week2", "outs")
scrna.dir.raw["pt1012_week4"] <- file.path(cellranger_base, "pt1012", "pt1012_week4", "outs")
scrna.dir.raw["pt1012_week14"] <- file.path(cellranger_base, "pt1012", "pt1012_week14", "outs")
scrna.dir.raw["pt1013_week1"] <- file.path(cellranger_base, "pt1013", "pt1013_week1", "outs")
scrna.dir.raw["pt1013_week2"] <- file.path(cellranger_base, "pt1013", "pt1013_week2", "outs")
scrna.dir.raw["pt1013_week4"] <- file.path(cellranger_base, "pt1013", "pt1013_week4", "outs")
scrna.dir.raw["pt1013_week14"] <- file.path(cellranger_base, "pt1013", "pt1013_week14", "outs")
scrna.dir.raw["pt1014_week1"] <- file.path(cellranger_base, "pt1014", "pt1014_week1", "outs")
scrna.dir.raw["pt1014_week2"] <- file.path(cellranger_base, "pt1014", "pt1014_week2", "outs")
scrna.dir.raw["pt1014_week4"] <- file.path(cellranger_base, "pt1014", "pt1014_week4", "outs")
scrna.dir.raw["pt1019_week1"] <- file.path(cellranger_base, "pt1019", "pt1019_week1", "outs")
scrna.dir.raw["pt1019_week4"] <- file.path(cellranger_base, "pt1019", "pt1019_week4", "outs")
scrna.dir.raw["pt1019_week14"] <- file.path(cellranger_base, "pt1019", "pt1019_week14", "outs")

sample.raw <- list()
for (i in 1:length(scrna.dir.raw)) {
  sample.raw[[i]] <- Read10X(paste0(scrna.dir.raw[[i]],"/multi/count/raw_feature_bc_matrix/"))
  rownames(sample.raw[[i]]$`Antibody Capture`) <- gsub("TotalSeqC", replacement = "ADT", rownames(sample.raw[[i]]$`Antibody Capture`) )
}
names(sample.raw) <- names(scrna.dir)

# ==============================================================================
# Aggregate Samples
# ==============================================================================
print("making barcodes match for adding to aggregated data")
aggrfile <- read.csv(file = aggregation_csv)

sample.data.ordered <- sample.data[aggrfile$sample_id]
for (i in 1:length(sample.data.ordered)) {
  j<- names(sample.data.ordered[i])
  t <- aggrfile[aggrfile$sample_id == j,]
  colnames(sample.data.ordered[[i]]$`Gene Expression`) <- gsub("1$", replacement = t$Order, colnames(sample.data.ordered[[i]]$`Gene Expression`))
  colnames(sample.data.ordered[[i]]$`Antibody Capture`) <- gsub("1$", replacement = t$Order, colnames(sample.data.ordered[[i]]$`Antibody Capture`))
}
sample.raw.ordered <- sample.raw[aggrfile$sample_id]
for (i in 1:length(sample.raw.ordered)) {
  j<- names(sample.raw.ordered[i])
  t <- aggrfile[aggrfile$sample_id == j,]
  colnames(sample.raw.ordered[[i]]$`Gene Expression`) <- gsub("1$", replacement = t$Order, colnames(sample.raw.ordered[[i]]$`Gene Expression`))
  colnames(sample.raw.ordered[[i]]$`Antibody Capture`) <- gsub("1$", replacement = t$Order, colnames(sample.raw.ordered[[i]]$`Antibody Capture`))
}

print("initialize rna object")
rna <- sample.raw.ordered[[1]]$`Gene Expression`
for (i in 2:length(sample.raw.ordered)) {
  rna <- cbind(sample.raw.ordered[[i]]$`Gene Expression`, rna)
}
print("initialize adt object")
prot <- sample.raw.ordered[[1]]$`Antibody Capture`
for (i in 2:length(sample.raw.ordered)) {
  prot <- cbind(sample.raw.ordered[[i]]$`Antibody Capture`, prot)
}

# ==============================================================================
# Quality Control and Background Detection
# ==============================================================================
rna_size = log10(Matrix::colSums(rna))
prot_size = log10(Matrix::colSums(prot))
ngene = Matrix::colSums(rna > 0)

md = as.data.frame(cbind(rna_size, ngene, prot_size))

md$bc = rownames(md)
md$category <- 'neutral'
print("cellranger cells object")
forcells <- sample.data.ordered[[1]]$`Gene Expression`
for (i in 2:length(sample.data.ordered)) {
  forcells <- cbind(sample.data.ordered[[i]]$`Gene Expression`, forcells)
}

positive_cells = md[md$bc %in% colnames(forcells),]$bc # using cellranger called cells not Seurat to call positive cellss
md[md$bc %in% positive_cells,]$category <- "positive cells"
md = md[md$rna_size > 0 & md$prot_size > 0, ] # remove barcodes with no capture

p1 = ggplot(md, aes(x = rna_size)) + geom_histogram(fill = "dodgerblue") + ggtitle("RNA library size \n distribution")
p2 = ggplot(md, aes(x = prot_size)) + geom_density(fill = "firebrick2") + ggtitle("Protein library size \n distribution")
p3 <- cowplot::plot_grid(p1, p2, nrow = 1)
droplets <- ggplot(md, aes(x = prot_size, y = rna_size, colour=category)) + geom_point(size = 0.5) 

print(p3)
print(droplets)

droplets2 <- ggplot(md[md$ngene> 0 & md$prot_size > 0, ], aes(x = ngene, y = prot_size, colour=category)) + geom_point(size = 0.5) 
print(droplets2)
# define a vector of background / empty droplet barcodes based on protein library size and mRNA content
background_drops = md[md$prot_size <3.5 & md$rna_size <2 & md$prot_size >1.5, ]$bc

negative_mtx_rawprot = prot[ , background_drops] %>% as.matrix()

# define a vector of cell-containing droplet barcodes based on protein library size and mRNA content 
# use cell-ranger filtered cells as positives
#length(positive_cells)
cells_mtx_rawprot = prot[ , positive_cells] %>% as.matrix()

md[md$bc %in% background_drops,]$category <- "background drops"

#visualize assignments
droplets <- ggplot(md[md$prot_size> 0 & md$rna_size > 0, ], aes(x = prot_size, y = rna_size, colour=category)) + geom_point(size = 0.5) 
print(droplets)

# ==============================================================================
# DSB Normalization
# ==============================================================================
# remove ADTs with poor staining or not used
dsb_norm_prot_iso_quant = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot, # cell containing droplets
  empty_drop_matrix = negative_mtx_rawprot, # estimate ambient noise with the background drops 
  denoise.counts = TRUE, # model and remove each cell's technical component
  use.isotype.control = TRUE, # use isotype controls to define the technical component
  isotype.control.name.vec = isotypes,
  quantile.clipping = TRUE)

# ==============================================================================
# Save Denoised ADT Data
# ==============================================================================
save(dsb_norm_prot_iso_quant, file = "dsb_norm_prot_isowithquant.Rdata")
