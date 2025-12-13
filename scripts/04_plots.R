# ==============================================================================
# Single-Cell RNA-Seq Analysis: Visualization Pipeline
# ==============================================================================
# This script generates publication-quality figures for the scRNA-seq analysis,
# including UMAP plots, feature plots, proportion plots, and patient-specific visualizations
#
# Input: Annotated Seurat object (RDS file)
# Output: Various plots and figures saved as PDF/PNG files
# ==============================================================================

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(viridis)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(scRepertoire)
library(ggpubr)
library(janitor)
library(readxl)
library(purrr)

#cite RColorBrewer for color palette
col_pal <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#E5C494", "#BC80BD", "#FBB4AE")
col_pal_extra <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#E5C494", "#BC80BD", "#FBB4AE", "#E78AC3", "#66C2A5", "#8DA0CB", "#FFD92F", "#B3B3B3")
object <- readRDS(file = "annotated_seurat.rds")

# ==============================================================================
# Panel 5A - DimPlot split by time
# ==============================================================================
####plot all timepoints together and highlight 
#Get cell barcodes by timepoint
object <- SetIdent(object, value = "timepoint")
week1_bc <- WhichCells(object, idents = "week1")
week2_bc <- WhichCells(object, idents = "week2")
week4_bc <- WhichCells(object, idents = "week4")
week14_bc <- WhichCells(object, idents = "week14")

dim_w1 <- DimPlot(object, cells.highlight = list(week1_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal_extra[1])) +
  theme(plot.title = element_text(color = col_pal_extra[1], size = 40), axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + labs(title = "Week 1")
dim_w2 <- DimPlot(object, cells.highlight = list(week2_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal_extra[2])) +
  theme(plot.title = element_text(color = col_pal_extra[2], size = 40), axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + labs(title = "Week 2")
dim_w4 <- DimPlot(object, cells.highlight = list(week4_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal_extra[14])) +
  theme(plot.title = element_text(color = col_pal_extra[14], size = 40), axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + labs(title = "Week 4")
dim_w14 <- DimPlot(object, cells.highlight = list(week14_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal_extra[12])) +
  theme(plot.title = element_text(color = col_pal_extra[12], size = 40), axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + labs(title = "Week 14")

A_main <- plot_grid(dim_w1+NoLegend(), dim_w2+NoLegend(), dim_w4+NoLegend(), dim_w14+NoLegend(), ncol=4)

# ==============================================================================
# Supplementary Panel 7B - DimPlot with seurat cluster labels
# ==============================================================================
A2_main_dimplot <- DimPlot(object, group.by = 'seurat_clusters', label = TRUE, label.size = 10) + 
  scale_color_manual(values = col_pal_extra) +
  theme(plot.title = element_text(size = 40), axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))


# ==============================================================================
# Panel 5B - DimPlot with cell classification final
# ==============================================================================
A3_main_dimplot <- DimPlot(object, group.by = 'cell_classification_final') + 
  scale_color_manual(values = c("light grey", "dark grey", col_pal[3:11])) +
  theme(plot.title = element_text(size = 40), axis.title = element_text(size = 40), 
        axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))


# ==============================================================================
# Supplementary Panel 7A - DimPlot with seurat cluster labels
# ==============================================================================
object <- SetIdent(object, value = "patient")
object_pt1012 <- subset(object, idents = "pt1012")
object_pt1013 <- subset(object, idents = "pt1013")
object_pt1014 <- subset(object, idents = "pt1014")
object_pt1019 <- subset(object, idents = "pt1019")

dim_pt1012_w1 <- DimPlot(object_pt1012, cells.highlight = list(week1_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal[1])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1012_w2 <- DimPlot(object_pt1012, cells.highlight = list(week2_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 2"), values=c("light grey", col_pal[2])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1012_w4 <- DimPlot(object_pt1012, cells.highlight = list(week4_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 4"), values=c("light grey", col_pal_extra[14])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) 
dim_pt1012_w14 <- DimPlot(object_pt1012, cells.highlight = list(week14_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 14"), values=c("light grey", col_pal_extra[12])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))

dim_pt1013_w1 <- DimPlot(object_pt1013, cells.highlight = list(week1_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal[1])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1013_w2 <- DimPlot(object_pt1013, cells.highlight = list(week2_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 2"), values=c("light grey", col_pal[2])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) 
dim_pt1013_w4 <- DimPlot(object_pt1013, cells.highlight = list(week4_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 4"), values=c("light grey", col_pal_extra[14])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) 
dim_pt1013_w14 <- DimPlot(object_pt1013, cells.highlight = list(week14_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 14"), values=c("light grey", col_pal_extra[12])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))

dim_pt1014_w1 <- DimPlot(object_pt1014, cells.highlight = list(week1_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal[1])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1014_w2 <- DimPlot(object_pt1014, cells.highlight = list(week2_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 2"), values=c("light grey", col_pal[2])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1014_w4 <- DimPlot(object_pt1014, cells.highlight = list(week4_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 4"), values=c("light grey", col_pal_extra[14])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1014_w14 <- DimPlot(object_pt1014, cells.highlight = list(week14_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 14"), values=c("light grey", col_pal_extra[12])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))

dim_pt1019_w1 <- DimPlot(object_pt1019, cells.highlight = list(week1_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 1"), values=c("light grey", col_pal[1])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1019_w2 <- DimPlot(object_pt1019, cells.highlight = list(week2_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 2"), values=c("light grey", col_pal[2])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1019_w4 <- DimPlot(object_pt1019, cells.highlight = list(week4_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 4"), values=c("light grey", col_pal_extra[14])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dim_pt1019_w14 <- DimPlot(object_pt1019, cells.highlight = list(week14_bc), sizes.highlight = 1, pt.size = 0.05) + 
  scale_color_manual(labels=c("Other timepoints", "Week 14"), values=c("light grey", col_pal_extra[12])) +
  theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.title=element_text(size=40), legend.text=element_text(size=40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))

A_supp_pt1012 <- plot_grid(dim_pt1012_w1 + NoLegend(), dim_pt1012_w2 + NoLegend(), dim_pt1012_w4 + NoLegend(), dim_pt1012_w14 + NoLegend(), ncol = 4)
A_supp_pt1013 <- plot_grid(dim_pt1013_w1 + NoLegend(), dim_pt1013_w2 + NoLegend(), dim_pt1013_w4 + NoLegend(), dim_pt1013_w14 + NoLegend(), ncol = 4)
A_supp_pt1014 <- plot_grid(dim_pt1014_w1 + NoLegend(), dim_pt1014_w2 + NoLegend(), dim_pt1014_w4 + NoLegend(), ncol = 4)
A_supp_pt1019 <- plot_grid(dim_pt1019_w1 + NoLegend(), dim_pt1019_w4 + NoLegend(), dim_pt1019_w14 + NoLegend(), ncol = 4)

A_supp <- plot_grid(A_supp_pt1012, A_supp_pt1013, A_supp_pt1014, A_supp_pt1019, ncol = 1)

################# Panel 5G,I; Supplementary Panel 7F,H - DimPlot with expanded clonotypes ###################

#pt1012
s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w14_filtered_contig_annotations.csv")



#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)


virus_CDR3 <- list('CASSRDSSNQPQHF','CASSLDGGPYEQYF','CATSIGGNQPQHF','CASSGGINEQFF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSYGAGGYNEQFF','CASSYGAGGYNEQFF','CASSYGAGGYNEQFF','CASSLIINEQFF','CASSLIINEQFF','CASSLIINEQFF')
toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

contig_list_pt1012 <- list(s1, s2, s3, s4)
contig_list_pt1012_sub <- list(s1_sub, s2_sub, s3_sub, s4_sub)

combined <- combineTCR(contig_list_pt1012_sub, 
                       samples = c("weekA1", "weekB2", "weekC4", "weekD14"), 
                       ID = c("ID", "ID", "ID", "ID"), cells ="T-AB")

combined$weekA1_ID$barcode <- gsub("weekA1_ID_","", combined$weekA1_ID$barcode)
combined$weekB2_ID$barcode <- gsub("weekB2_ID_","", combined$weekB2_ID$barcode)
combined$weekC4_ID$barcode <- gsub("weekC4_ID_","", combined$weekC4_ID$barcode)
combined$weekD14_ID$barcode <- gsub("weekD14_ID_","", combined$weekD14_ID$barcode)

integrated_pt1012 <- combineExpression(combined, object_pt1012, 
                                       cloneCall="gene", proportion = FALSE, group.by = "sample",
                                       cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

integrated_pt1012 <- highlightClonotypes(integrated_pt1012, cloneCall= "aa", 
                                         sequence = c("CAASGGADGLTF_CASSSASSGEETQYF",
                                                      "CALSEAQAAGNKLTF_CASSDEDRVEKAFF",
                                                      "CAMSVLSGNTGKLIF_CASSVASLTVGTGELFF",
                                                      "CAVGGDDAGNMLTF_CASMPSTGPNEKLFF",
                                                      "CAVNNNFNKFYF_CASSFGGLDAEAFF",
                                                      "CAASDNTGNQFYF_CASSLAQEVNQPQHF",
                                                      "CAPSVGATNKLIF_CSVELPGDTIYF",
                                                      "CALRAPYGGSQGNLIF_CASSLPGHEKLFF",
                                                      "CAVPRFSDGQKLLF_CASSFGRNQPQHF"))

clonotype_table <- compareClonotypes(combined, numbers = 50, samples = c("weekA1_ID", "weekB2_ID", "weekC4_ID", "weekD14_ID"), 
                                     cloneCall="aa", graph = "alluvial", exportTable = TRUE)
#Print top 10 clonotypes easy (11 for pt1012 because 2 clonotypes have same proportion and only print 9)
head(unique(clonotype_table[order(-clonotype_table$Proportion),]$Clonotypes), 11)

#looking at cell types for specific clonotypes that expanded based on alluvial plots
df_clone_celltype <- data.frame(integrated_pt1012@meta.data$CTaa, integrated_pt1012@meta.data$cell_classification_final, integrated_pt1012@meta.data$timepoint, integrated_pt1012@meta.data$cloneType, integrated_pt1012@meta.data$barcode)
#merge with clonotype proportions table
df_clone_cell_prop <- merge(df_clone_celltype, clonotype_table, by.x="integrated_pt1012.meta.data.CTaa", by.y="Clonotypes")
#pick top clonotypes
df_clone_cell_prop_top <- df_clone_cell_prop[df_clone_cell_prop$integrated_pt1012.meta.data.CTaa == c(head(unique(df_clone_cell_prop[order(-df_clone_cell_prop$Proportion), ]$integrated_pt1012.meta.data.CTaa), 11)),]

CAMSVLSGNTGKLIF_CASSVASLTVGTGELFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CAMSVLSGNTGKLIF_CASSVASLTVGTGELFF", ]$integrated_pt1012.meta.data.barcode
CAVGGDDAGNMLTF_CASMPSTGPNEKLFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CAVGGDDAGNMLTF_CASMPSTGPNEKLFF", ]$integrated_pt1012.meta.data.barcode 
CALRAPYGGSQGNLIF_CASSLPGHEKLFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CALRAPYGGSQGNLIF_CASSLPGHEKLFF", ]$integrated_pt1012.meta.data.barcode
CAVNNNFNKFYF_CASSFGGLDAEAFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CAVNNNFNKFYF_CASSFGGLDAEAFF", ]$integrated_pt1012.meta.data.barcode
CAPSVGATNKLIF_CSVELPGDTIYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CAPSVGATNKLIF_CSVELPGDTIYF", ]$integrated_pt1012.meta.data.barcode
CVVSDRGSTLGRLYF_CASSEPPAGTGAEKLFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CVVSDRGSTLGRLYF_CASSEPPAGTGAEKLFF", ]$integrated_pt1012.meta.data.barcode
CAVPRFSDGQKLLF_CASSFGRNQPQHF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CAVPRFSDGQKLLF_CASSFGRNQPQHF", ]$integrated_pt1012.meta.data.barcode
CAASDNTGNQFYF_CASSLAQEVNQPQHF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CAASDNTGNQFYF_CASSLAQEVNQPQHF", ]$integrated_pt1012.meta.data.barcode
CAASGGADGLTF_CASSSASSGEETQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CAASGGADGLTF_CASSSASSGEETQYF", ]$integrated_pt1012.meta.data.barcode
CALSEAQAAGNKLTF_CASSDEDRVEKAFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1012.meta.data.CTaa == "CALSEAQAAGNKLTF_CASSDEDRVEKAFF", ]$integrated_pt1012.meta.data.barcode
# CAASDNTGNQFYF_CASSLAQEVNQPQHF        CAASGGADGLTF_CASSSASSGEETQYF      CALRAPYGGSQGNLIF_CASSLPGHEKLFF 
# CALSEAQAAGNKLTF_CASSDEDRVEKAFF   CAMSVLSGNTGKLIF_CASSVASLTVGTGELFF          CAPSVGATNKLIF_CSVELPGDTIYF 
# CAVGGDDAGNMLTF_CASMPSTGPNEKLFF CAVKGNYGGSQGNLIF_CASSGRTGGLSGANVLTF         CAVNNNFNKFYF_CASSFGGLDAEAFF 
# CAVPRFSDGQKLLF_CASSFGRNQPQHF   CVVSDRGSTLGRLYF_CASSEPPAGTGAEKLFF 

highlight_plot_1 <- DimPlot(integrated_pt1012, cells.highlight = list(CAASGGADGLTF_CASSSASSGEETQYF,
                                                                      CALSEAQAAGNKLTF_CASSDEDRVEKAFF,
                                                                      CAMSVLSGNTGKLIF_CASSVASLTVGTGELFF,
                                                  CAVGGDDAGNMLTF_CASMPSTGPNEKLFF,
                                                  CAVNNNFNKFYF_CASSFGGLDAEAFF,
                                                  CAASDNTGNQFYF_CASSLAQEVNQPQHF,
                                                  CAPSVGATNKLIF_CSVELPGDTIYF,
                                                  CALRAPYGGSQGNLIF_CASSLPGHEKLFF,
                                                  CAVPRFSDGQKLLF_CASSFGRNQPQHF), sizes.highlight = 5, pt.size = 0.5, 
        cols.highlight = col_pal[1])  +  ggtitle("pt1012 Top 10 expanded clonotypes") +
  scale_color_manual(labels=c("NA", "CAASGGADGLTF_CASSSASSGEETQYF",
                              "CALSEAQAAGNKLTF_CASSDEDRVEKAFF",
                              "CAMSVLSGNTGKLIF_CASSVASLTVGTGELFF",
                              "CAVGGDDAGNMLTF_CASMPSTGPNEKLFF",
                              "CAVNNNFNKFYF_CASSFGGLDAEAFF",
                              "CAASDNTGNQFYF_CASSLAQEVNQPQHF",
                              "CAPSVGATNKLIF_CSVELPGDTIYF",
                              "CALRAPYGGSQGNLIF_CASSLPGHEKLFF",
                              "CAVPRFSDGQKLLF_CASSFGRNQPQHF"), 
                     breaks = c("Unselected", "Group_1", "Group_2", "Group_3", "Group_4", "Group_5",
                                "Group_6", "Group_7", "Group_8", "Group_9"),
                     values=c("light grey", col_pal)) + theme(legend.title=element_text(size=40), 
          legend.text=element_text(size=40), plot.title = element_text(size = 40), axis.title = element_text(size = 40), axis.text = element_text(size = 40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
  # geom_segment(aes(x = -6.5, y = -2, xend = -4, yend = 5.25)) + geom_segment(aes(x = -4, y = 5.25, xend = 1.75, yend = 2.75)) + 
  # geom_segment(aes(x = 1.75, y = 2.75, xend = -3, yend = -4)) + geom_segment(aes(x = -3, y = -4, xend = -6.5, yend = -2))+ 
  # geom_segment(aes(x = -2.5, y = -4.5, xend = 3.25, yend = 3.75)) + geom_segment(aes(x = 3.25, y = 3.75, xend = 10, yend = 2.25)) + 
  # geom_segment(aes(x = 10, y = 2.25, xend = 7, yend = -8)) + geom_segment(aes(x = 7, y = -8, xend = -2.5, yend = -4.5))
  

#pt1013
s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w14_filtered_contig_annotations.csv")

#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSPGGQETQYF','CASGTGSYEQYF','CASGTGSYEQYF','CASSQIRGGYEQYF','CASSIVGQINTEAFF')
toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]


contig_list_pt1013 <- list(s1, s2, s3, s4)
contig_list_pt1013_sub <- list(s1_sub, s2_sub, s3_sub, s4_sub)


combined <- combineTCR(contig_list_pt1013_sub, 
                       samples = c("weekA1", "weekB2", "weekC4", "weekD14"), 
                       ID = c("ID", "ID", "ID", "ID"), cells ="T-AB")
#for list of clonotypes
clonotype_table <- compareClonotypes(combined, numbers = 50, samples = c("weekA1_ID", "weekB2_ID", "weekC4_ID", "weekD14_ID"), 
                                     cloneCall="aa", graph = "alluvial", exportTable = TRUE)
#Print top 10 clonotypes easy
head(unique(clonotype_table[order(-clonotype_table$Proportion),]$Clonotypes), 11)

combined$weekA1_ID$barcode <- gsub("weekA1_ID_","", combined$weekA1_ID$barcode)
combined$weekB2_ID$barcode <- gsub("weekB2_ID_","", combined$weekB2_ID$barcode)
combined$weekC4_ID$barcode <- gsub("weekC4_ID_","", combined$weekC4_ID$barcode)
combined$weekD14_ID$barcode <- gsub("weekD14_ID_","", combined$weekD14_ID$barcode)

integrated_pt1013 <- combineExpression(combined, object_pt1013, 
                                       cloneCall="gene", proportion = FALSE, group.by = "sample",
                                       cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

integrated_pt1013 <- highlightClonotypes(integrated_pt1013, cloneCall= "aa", 
                                         sequence = c("CACQDGGGADGLTF_CASSVKISTDTQYF",
                                                      "CAGQLGGGSNYKLTF_CASSERGVMVSSNQPQHF",
                                                      "CALDMRQGGKLIF_CASRLHTGTGTSGANVLTF",
                                                      "CALINDYKLSF_CASSQTLWGSGELTDTQYF",
                                                      "CALSEARNAGNMLTF_CASSLTGQGDYEQYF",
                                                      "CAMREGLIKAAGNKLTF_CASSVRSRGDDSPLHF",                                                      
                                                      "CAVSTSGGSYIPTF_CASSHPTSGRETQYF",
                                                      "CIVRVLPSNTGKLIF_CASSQDSVGRSSYEQYF",
                                                      "NA_CASSLDGVSTDTQYF",
                                                      "CAVSGIKAAGNKLTF_CASSVAALDQPQHF",
                                                      "CAVSGIKAAGNKLTF_CASSVAALDQPQHF",
                                                      "CAVTPGGGKLIF_CVVRGLTGGGNKLTF_CASSQSGGGGGAYEQYF"
))

#looking at cell types for specific clonotypes that expanded based on alluvial plots
df_clone_celltype <- data.frame(integrated_pt1013@meta.data$CTaa, integrated_pt1013@meta.data$cell_classification_final, integrated_pt1013@meta.data$timepoint, integrated_pt1013@meta.data$cloneType, integrated_pt1013@meta.data$barcode)
#merge with clonotype proportions table
df_clone_cell_prop <- merge(df_clone_celltype, clonotype_table, by.x="integrated_pt1013.meta.data.CTaa", by.y="Clonotypes")
#pick top clonotypes
df_clone_cell_prop_top <- df_clone_cell_prop[df_clone_cell_prop$integrated_pt1013.meta.data.CTaa == c(head(unique(df_clone_cell_prop[order(-df_clone_cell_prop$Proportion), ]$integrated_pt1013.meta.data.CTaa), 10)),]


CAVSGIKAAGNKLTF_CASSVAALDQPQHF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CAVSGIKAAGNKLTF_CASSVAALDQPQHF", ]$integrated_pt1013.meta.data.barcode
CAVSTSGGSYIPTF_CASSHPTSGRETQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CAVSTSGGSYIPTF_CASSHPTSGRETQYF", ]$integrated_pt1013.meta.data.barcode
CIVRVLPSNTGKLIF_CASSQDSVGRSSYEQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CIVRVLPSNTGKLIF_CASSQDSVGRSSYEQYF", ]$integrated_pt1013.meta.data.barcode
NA_CASSLDGVSTDTQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "NA_CASSLDGVSTDTQYF", ]$integrated_pt1013.meta.data.barcode
CAGQLGGGSNYKLTF_CASSERGVMVSSNQPQHF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CAGQLGGGSNYKLTF_CASSERGVMVSSNQPQHF", ]$integrated_pt1013.meta.data.barcode
CAVTPGGGKLIF_CVVRGLTGGGNKLTF_CASSQSGGGGGAYEQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CAVTPGGGKLIF;CVVRGLTGGGNKLTF_CASSQSGGGGGAYEQYF", ]$integrated_pt1013.meta.data.barcode
CAMREGLIKAAGNKLTF_CASSVRSRGDDSPLHF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CAMREGLIKAAGNKLTF_CASSVRSRGDDSPLHF", ]$integrated_pt1013.meta.data.barcode
CALINDYKLSF_CASSQTLWGSGELTDTQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CALINDYKLSF_CASSQTLWGSGELTDTQYF", ]$integrated_pt1013.meta.data.barcode
CACQDGGGADGLTF_CASSVKISTDTQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CACQDGGGADGLTF_CASSVKISTDTQYF", ]$integrated_pt1013.meta.data.barcode
CALDMRQGGKLIF_CASRLHTGTGTSGANVLTF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CALDMRQGGKLIF_CASRLHTGTGTSGANVLTF", ]$integrated_pt1013.meta.data.barcode
CALSEARNAGNMLTF_CASSLTGQGDYEQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1013.meta.data.CTaa == "CALSEARNAGNMLTF_CASSLTGQGDYEQYF", ]$integrated_pt1013.meta.data.barcode

highlight_plot_2 <- DimPlot(integrated_pt1013, cells.highlight = list(CACQDGGGADGLTF_CASSVKISTDTQYF,
                                                                      CAGQLGGGSNYKLTF_CASSERGVMVSSNQPQHF,
                                                                      CALDMRQGGKLIF_CASRLHTGTGTSGANVLTF,
                                                                      CALINDYKLSF_CASSQTLWGSGELTDTQYF,
                                                                      CALSEARNAGNMLTF_CASSLTGQGDYEQYF,
                                                                      CAMREGLIKAAGNKLTF_CASSVRSRGDDSPLHF,                                                      
                                                                      CAVSTSGGSYIPTF_CASSHPTSGRETQYF,
                                                                      CIVRVLPSNTGKLIF_CASSQDSVGRSSYEQYF,
                                                                      NA_CASSLDGVSTDTQYF,
                                                                      CAVSGIKAAGNKLTF_CASSVAALDQPQHF), sizes.highlight = 5, pt.size = 0.5, 
        cols.highlight = col_pal[1])  +  ggtitle("pt1013 Top 10 expanded clonotypes") +
  scale_color_manual(labels=c("NA", "CACQDGGGADGLTF_CASSVKISTDTQYF",
                              "CAGQLGGGSNYKLTF_CASSERGVMVSSNQPQHF",
                              "CALDMRQGGKLIF_CASRLHTGTGTSGANVLTF",
                              "CALINDYKLSF_CASSQTLWGSGELTDTQYF",
                              "CALSEARNAGNMLTF_CASSLTGQGDYEQYF",
                              "CAMREGLIKAAGNKLTF_CASSVRSRGDDSPLHF",                                                      
                              "CAVSTSGGSYIPTF_CASSHPTSGRETQYF",
                              "CIVRVLPSNTGKLIF_CASSQDSVGRSSYEQYF",
                              "NA_CASSLDGVSTDTQYF",
                              "CAVSGIKAAGNKLTF_CASSVAALDQPQHF"), 
                     breaks = c("Unselected", "Group_1", "Group_2", "Group_3", "Group_4", "Group_5",
                                "Group_6", "Group_7", "Group_8", "Group_9", "Group_10"),
                     values=c("light grey", col_pal)) + 
  theme(legend.title=element_text(size=40), legend.text=element_text(size=40), plot.title = element_text(size = 40), 
        axis.title = element_text(size = 40), axis.text = element_text(size = 40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
  # geom_segment(aes(x = -6.5, y = -2, xend = -4, yend = 5.25)) + geom_segment(aes(x = -4, y = 5.25, xend = 1.75, yend = 2.75)) + 
  # geom_segment(aes(x = 1.75, y = 2.75, xend = -3, yend = -4)) + geom_segment(aes(x = -3, y = -4, xend = -6.5, yend = -2))+ 
  # geom_segment(aes(x = -2.5, y = -4.5, xend = 3.25, yend = 3.75)) + geom_segment(aes(x = 3.25, y = 3.75, xend = 10, yend = 2.25)) + 
  # geom_segment(aes(x = 10, y = 2.25, xend = 7, yend = -8)) + geom_segment(aes(x = 7, y = -8, xend = -2.5, yend = -4.5))



#pt1014
s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w4_filtered_contig_annotations.csv")
#s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w14_filtered_contig_annotations.csv")


#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
# s4$is_cell = toupper(s4$is_cell)
# s4$high_confidence = toupper(s4$high_confidence)
# s4$full_length = toupper(s4$full_length)
# s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSPPGGSYEQYF','CAISESKGNYGYTF','CASSLTGTGSYNEQFF','CASSPPGQVNNEQFF','CAISTGNTEAFF','CAISTGNTEAFF')

toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

#toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
#s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]


contig_list_pt1014 <- list(s1, s2, s3)
contig_list_pt1014_sub <- list(s1_sub, s2_sub, s3_sub)

combined <- combineTCR(contig_list_pt1014_sub, 
                       samples = c("weekA1", "weekB2", "weekC4"), 
                       ID = c("ID", "ID", "ID"), cells ="T-AB")

#for list of clonotypes
clonotype_table <- compareClonotypes(combined, numbers = 50, samples = c("weekA1_ID", "weekB2_ID", "weekC4_ID"), 
                                     cloneCall="aa", graph = "alluvial", exportTable = TRUE)
#Print top 10 clonotypes easy
head(unique(clonotype_table[order(-clonotype_table$Proportion),]$Clonotypes), 11)


combined$weekA1_ID$barcode <- gsub("weekA1_ID_","", combined$weekA1_ID$barcode)
combined$weekB2_ID$barcode <- gsub("weekB2_ID_","", combined$weekB2_ID$barcode)
combined$weekC4_ID$barcode <- gsub("weekC4_ID_","", combined$weekC4_ID$barcode)
#combined$weekD14_ID$barcode <- gsub("weekD14_ID_","", combined$weekD14_ID$barcode)

integrated_pt1014 <- combineExpression(combined, object_pt1014, 
                                       cloneCall="gene", proportion = FALSE, group.by = "sample",
                                       cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

integrated_pt1014 <- highlightClonotypes(integrated_pt1014, cloneCall= "aa", 
                                         sequence = c("CAASVNDYKLSF_CASSLELAETGELFF",
                                                      "CAAYEGTDKLIF_CASRQPSSYEQYF",
                                                      "CAGRPSNNFGNEKLTF_CASSLGGRNIQYF",
                                                      "CALRSNRDDKIIF_CSARDGDRGSYEQYF",                                                      
                                                      "CALSEDNFNKFYF_CASSLDPGKGEPLHF",
                                                      "CARASGGSYIPTF_CASSDTGTGRDTQYF",
                                                      "CAVPGGGADGLTF_CATSWGRGEQFF",                                                      
                                                      "CAVPTGGADGLTF_CATSAGRGELFF",
                                                      "NA_CASSLELAETGELFF",
                                                      "NA_CASSLGGRNIQYF"
))

#looking at cell types for specific clonotypes that expanded based on alluvial plots
df_clone_celltype <- data.frame(integrated_pt1014@meta.data$CTaa, integrated_pt1014@meta.data$cell_classification_final, integrated_pt1014@meta.data$timepoint, integrated_pt1014@meta.data$cloneType, integrated_pt1014@meta.data$barcode)
#merge with clonotype proportions table
df_clone_cell_prop <- merge(df_clone_celltype, clonotype_table, by.x="integrated_pt1014.meta.data.CTaa", by.y="Clonotypes")
#pick top clonotypes
df_clone_cell_prop_top <- df_clone_cell_prop[df_clone_cell_prop$integrated_pt1014.meta.data.CTaa == c(head(unique(df_clone_cell_prop[order(-df_clone_cell_prop$Proportion), ]$integrated_pt1014.meta.data.CTaa), 10)),]

CALSEDNFNKFYF_CASSLDPGKGEPLHF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CALSEDNFNKFYF_CASSLDPGKGEPLHF", ]$integrated_pt1014.meta.data.barcode
CAASVNDYKLSF_CASSLELAETGELFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CAASVNDYKLSF_CASSLELAETGELFF", ]$integrated_pt1014.meta.data.barcode
CAGRPSNNFGNEKLTF_CASSLGGRNIQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CAGRPSNNFGNEKLTF_CASSLGGRNIQYF", ]$integrated_pt1014.meta.data.barcode
CAVPTGGADGLTF_CATSAGRGELFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CAVPTGGADGLTF_CATSAGRGELFF", ]$integrated_pt1014.meta.data.barcode
NA_CASSLELAETGELFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "NA_CASSLELAETGELFF", ]$integrated_pt1014.meta.data.barcode
NA_CASSLGGRNIQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "NA_CASSLGGRNIQYF", ]$integrated_pt1014.meta.data.barcode
CAVPGGGADGLTF_CATSWGRGEQFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CAVPGGGADGLTF_CATSWGRGEQFF", ]$integrated_pt1014.meta.data.barcode
CALRSNRDDKIIF_CSARDGDRGSYEQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CALRSNRDDKIIF_CSARDGDRGSYEQYF", ]$integrated_pt1014.meta.data.barcode
CAAYEGTDKLIF_CASRQPSSYEQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CAAYEGTDKLIF_CASRQPSSYEQYF", ]$integrated_pt1014.meta.data.barcode
CARASGGSYIPTF_CASSDTGTGRDTQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1014.meta.data.CTaa == "CARASGGSYIPTF_CASSDTGTGRDTQYF", ]$integrated_pt1014.meta.data.barcode


highlight_plot_3 <- DimPlot(integrated_pt1014, cells.highlight = list(CAASVNDYKLSF_CASSLELAETGELFF,
                                                                      CAAYEGTDKLIF_CASRQPSSYEQYF,
                                                                      CAGRPSNNFGNEKLTF_CASSLGGRNIQYF,
                                                                      CALRSNRDDKIIF_CSARDGDRGSYEQYF,                                                      
                                                                      CALSEDNFNKFYF_CASSLDPGKGEPLHF,
                                                                      CARASGGSYIPTF_CASSDTGTGRDTQYF,
                                                                      CAVPGGGADGLTF_CATSWGRGEQFF,                                                      
                                                                      CAVPTGGADGLTF_CATSAGRGELFF,
                                                                      NA_CASSLELAETGELFF,
                                                                      NA_CASSLGGRNIQYF), sizes.highlight = 5, pt.size = 0.5, 
        cols.highlight = col_pal[1])  +  ggtitle("pt1014 Top 10 expanded clonotypes") +
  scale_color_manual(labels=c("NA", "CAASVNDYKLSF_CASSLELAETGELFF",
                              "CAAYEGTDKLIF_CASRQPSSYEQYF",
                              "CAGRPSNNFGNEKLTF_CASSLGGRNIQYF",
                              "CALRSNRDDKIIF_CSARDGDRGSYEQYF",                                                      
                              "CALSEDNFNKFYF_CASSLDPGKGEPLHF",
                              "CARASGGSYIPTF_CASSDTGTGRDTQYF",
                              "CAVPGGGADGLTF_CATSWGRGEQFF",                                                      
                              "CAVPTGGADGLTF_CATSAGRGELFF",
                              "NA_CASSLELAETGELFF",
                              "NA_CASSLGGRNIQYF"), 
                     breaks = c("Unselected", "Group_1", "Group_2", "Group_3", "Group_4", "Group_5",
                                "Group_6", "Group_7", "Group_8", "Group_9", "Group_10"), values=c("light grey", col_pal)) + theme(legend.title=element_text(size=40), 
                                                                                                         legend.text=element_text(size=40), plot.title = element_text(size = 40), axis.title = element_text(size = 40), axis.text = element_text(size = 40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
  # geom_segment(aes(x = -6.5, y = -2, xend = -4, yend = 5.25)) + geom_segment(aes(x = -4, y = 5.25, xend = 1.75, yend = 2.75)) + 
  # geom_segment(aes(x = 1.75, y = 2.75, xend = -3, yend = -4)) + geom_segment(aes(x = -3, y = -4, xend = -6.5, yend = -2))+ 
  # geom_segment(aes(x = -2.5, y = -4.5, xend = 3.25, yend = 3.75)) + geom_segment(aes(x = 3.25, y = 3.75, xend = 10, yend = 2.25)) + 
  # geom_segment(aes(x = 10, y = 2.25, xend = 7, yend = -8)) + geom_segment(aes(x = 7, y = -8, xend = -2.5, yend = -4.5))

#pt1019

s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w1_filtered_contig_annotations.csv")
#s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w14_filtered_contig_annotations.csv")

#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
# s2$is_cell = toupper(s2$is_cell)
# s2$high_confidence = toupper(s2$high_confidence)
# s2$full_length = toupper(s2$full_length)
# s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSSRLAGSYNEQFF','CASSYTETQYF','CASSYTETQYF','CASSYTETQYF','CASSYTETQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF')

toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

#toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
#s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

contig_list_pt1019 <- list(s1, s3, s4)
contig_list_pt1019_sub <- list(s1_sub, s3_sub, s4_sub)

combined <- combineTCR(contig_list_pt1019_sub, 
                       samples = c("weekA1", "weekC4", "weekD14"), 
                       ID = c("ID", "ID", "ID"), cells ="T-AB")

#for list of clonotypes
clonotype_table <- compareClonotypes(combined, numbers = 50, samples = c("weekA1_ID", "weekC4_ID", "weekD14_ID"), 
                                     cloneCall="aa", graph = "alluvial", exportTable = TRUE)
#Print top 10 clonotypes easy
head(unique(clonotype_table[order(-clonotype_table$Proportion),]$Clonotypes), 11)


combined$weekA1_ID$barcode <- gsub("weekA1_ID_","", combined$weekA1_ID$barcode)
#combined$weekB2_ID$barcode <- gsub("weekB2_ID_","", combined$weekB2_ID$barcode)
combined$weekC4_ID$barcode <- gsub("weekC4_ID_","", combined$weekC4_ID$barcode)
combined$weekD14_ID$barcode <- gsub("weekD14_ID_","", combined$weekD14_ID$barcode)

integrated_pt1019 <- combineExpression(combined, object_pt1019, 
                                       cloneCall="gene", proportion = FALSE, group.by = "sample",
                                       cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

integrated_pt1019 <- highlightClonotypes(integrated_pt1019, cloneCall= "aa", 
                                         sequence = c("CAAPQPGYSTLTF_CASSLRGGSTEAFF",
                                                      "CAARLGGANNLFF_CASSRTGVVTEAFF",
                                                      "CAGPRRYNFNKFYF_CASSSPGLAGGPGNTGELFF",                                                      
                                                      "CAMSGRWAGGGNKLTF_CAWSPTESPEQFF",
                                                      "CAVADSGANSKLTF_CASSFLGTDTQYF",
                                                      "CAVRDRDYGGATNKLIF_CASSPQPRTGLDGYTF",
                                                      "CAYQGGKLIF_CASSSREAGTYNEQFF",
                                                      "CVVSELNNARLMF;CVFRRGQKLLF_CASSIVGVLSGELFF",
                                                      "NA_CASSPQPRTGLDGYTF",                                                      
                                                      "CAMRNRDDKIIF_CASSDWAQAVAAYEQYF"
))

#looking at cell types for specific clonotypes that expanded based on alluvial plots
df_clone_celltype <- data.frame(integrated_pt1019@meta.data$CTaa, integrated_pt1019@meta.data$cell_classification_final, integrated_pt1019@meta.data$timepoint, integrated_pt1019@meta.data$cloneType, integrated_pt1019@meta.data$barcode)
#merge with clonotype proportions table
df_clone_cell_prop <- merge(df_clone_celltype, clonotype_table, by.x="integrated_pt1019.meta.data.CTaa", by.y="Clonotypes")
#pick top clonotypes
df_clone_cell_prop_top <- df_clone_cell_prop[df_clone_cell_prop$integrated_pt1019.meta.data.CTaa == c(head(unique(df_clone_cell_prop[order(-df_clone_cell_prop$Proportion), ]$integrated_pt1019.meta.data.CTaa), 10)),]

CVVSELNNARLMF_CVFRRGQKLLF_CASSIVGVLSGELFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CVVSELNNARLMF;CVFRRGQKLLF_CASSIVGVLSGELFF", ]$integrated_pt1019.meta.data.barcode
CAMSGRWAGGGNKLTF_CAWSPTESPEQFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAMSGRWAGGGNKLTF_CAWSPTESPEQFF", ]$integrated_pt1019.meta.data.barcode
CAARLGGANNLFF_CASSRTGVVTEAFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAARLGGANNLFF_CASSRTGVVTEAFF", ]$integrated_pt1019.meta.data.barcode
CAVRDRDYGGATNKLIF_CASSPQPRTGLDGYTF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAVRDRDYGGATNKLIF_CASSPQPRTGLDGYTF", ]$integrated_pt1019.meta.data.barcode
CAGPRRYNFNKFYF_CASSSPGLAGGPGNTGELFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAGPRRYNFNKFYF_CASSSPGLAGGPGNTGELFF", ]$integrated_pt1019.meta.data.barcode
CAAPQPGYSTLTF_CASSLRGGSTEAFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAAPQPGYSTLTF_CASSLRGGSTEAFF", ]$integrated_pt1019.meta.data.barcode
CAYQGGKLIF_CASSSREAGTYNEQFF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAYQGGKLIF_CASSSREAGTYNEQFF", ]$integrated_pt1019.meta.data.barcode
CAMRNRDDKIIF_CASSDWAQAVAAYEQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAMRNRDDKIIF_CASSDWAQAVAAYEQYF", ]$integrated_pt1019.meta.data.barcode
CAVADSGANSKLTF_CASSFLGTDTQYF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "CAVADSGANSKLTF_CASSFLGTDTQYF", ]$integrated_pt1019.meta.data.barcode
NA_CASSPQPRTGLDGYTF <- df_clone_cell_prop_top[df_clone_cell_prop_top$integrated_pt1019.meta.data.CTaa == "NA_CASSPQPRTGLDGYTF", ]$integrated_pt1019.meta.data.barcode

highlight_plot_4 <- DimPlot(integrated_pt1019, cells.highlight = list(CAAPQPGYSTLTF_CASSLRGGSTEAFF,
                                                                      CAARLGGANNLFF_CASSRTGVVTEAFF,
                                                                      CAGPRRYNFNKFYF_CASSSPGLAGGPGNTGELFF,                                                      
                                                                      CAMSGRWAGGGNKLTF_CAWSPTESPEQFF,
                                                                      CAVADSGANSKLTF_CASSFLGTDTQYF,
                                                                      CAVRDRDYGGATNKLIF_CASSPQPRTGLDGYTF,
                                                                      CAYQGGKLIF_CASSSREAGTYNEQFF,
                                                                      CVVSELNNARLMF_CVFRRGQKLLF_CASSIVGVLSGELFF,
                                                                      NA_CASSPQPRTGLDGYTF,                                                      
                                                                      CAMRNRDDKIIF_CASSDWAQAVAAYEQYF), sizes.highlight = 5, pt.size = 0.5, 
        cols.highlight = col_pal[1])  +  ggtitle("pt1019 Top 10 expanded clonotypes") +
  scale_color_manual(labels=c("NA", "CAAPQPGYSTLTF_CASSLRGGSTEAFF",
                              "CAARLGGANNLFF_CASSRTGVVTEAFF",
                              "CAGPRRYNFNKFYF_CASSSPGLAGGPGNTGELFF",                                                      
                              "CAMSGRWAGGGNKLTF_CAWSPTESPEQFF",
                              "CAVADSGANSKLTF_CASSFLGTDTQYF",
                              "CAVRDRDYGGATNKLIF_CASSPQPRTGLDGYTF",
                              "CAYQGGKLIF_CASSSREAGTYNEQFF",
                              "CVVSELNNARLMF;CVFRRGQKLLF_CASSIVGVLSGELFF",
                              "NA_CASSPQPRTGLDGYTF",                                                      
                              "CAMRNRDDKIIF_CASSDWAQAVAAYEQYF"), 
                     breaks = c("Unselected", "Group_1", "Group_2", "Group_3", "Group_4", "Group_5",
                                "Group_6", "Group_7", "Group_8", "Group_9", "Group_10"),values=c("light grey", col_pal)) + theme(legend.title=element_text(size=40), 
                                                                                               legend.text=element_text(size=40), plot.title = element_text(size = 40), axis.title = element_text(size = 40), axis.text = element_text(size = 40)) + 
  guides(color = guide_legend(override.aes = list(size = 10)))

highlight_main <- ggarrange(highlight_plot_1+NoLegend()+ggtitle(""), highlight_plot_2+NoLegend()+ggtitle(""), highlight_plot_3+NoLegend()+ggtitle(""), highlight_plot_4+NoLegend()+ggtitle(""), ncol=2, nrow=2)

################# Panel 5F,H; Supplementary Panel 7E,G - Alluvial plots with expanded clonotypes ###################

#pt1012
s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w14_filtered_contig_annotations.csv")



#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)


virus_CDR3 <- list('CASSRDSSNQPQHF','CASSLDGGPYEQYF','CATSIGGNQPQHF','CASSGGINEQFF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSYGAGGYNEQFF','CASSYGAGGYNEQFF','CASSYGAGGYNEQFF','CASSLIINEQFF','CASSLIINEQFF','CASSLIINEQFF')
toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

object_pt1012 <- SetIdent(object_pt1012, value = "cell_classification_final")
object_pt1012 <- RenameIdents(object_pt1012, "B cells" = "Other", "T/monocyte doublets" = "Other", "CD3+; CD56+" = "Other", 
                       "CD4 memory" = "CD4", "CD4 naïve" = "CD4", "CD8 memory/ TEMRA" = "CD8", "CD8 naïve/memory" = "CD8", 
                       "CD8 naïve" = "CD8", "gamma/delta" = "Other", "Proliferating T" = "Other","Treg" = "Other")
object_pt1012[["CD4vsCD8.ids"]] <- Idents(object_pt1012)

object_pt1012_CD4 <- subset(x = object_pt1012, subset = CD4vsCD8.ids == "CD4")
object_pt1012_CD8 <- subset(x = object_pt1012, subset = CD4vsCD8.ids == "CD8")

pt1012_CD4_barcodes <- WhichCells(object_pt1012_CD4, idents = "CD4")
pt1012_CD8_barcodes <- WhichCells(object_pt1012_CD8, idents = "CD8")

s1_sub_pt1012_CD8 <- s1_sub[s1_sub$barcode %in% pt1012_CD8_barcodes,]
s2_sub_pt1012_CD8 <- s2_sub[s2_sub$barcode %in% pt1012_CD8_barcodes,]
s3_sub_pt1012_CD8 <- s3_sub[s3_sub$barcode %in% pt1012_CD8_barcodes,]
s4_sub_pt1012_CD8 <- s4_sub[s4_sub$barcode %in% pt1012_CD8_barcodes,]

contig_list_pt1012_CD8 <- list(s1_sub_pt1012_CD8, s2_sub_pt1012_CD8, s3_sub_pt1012_CD8, s4_sub_pt1012_CD8)

combined_pt1012_CD8 <- combineTCR(contig_list_pt1012_CD8, 
                                           samples = c("weekA1", "weekB2", "weekC4", "weekD14"), 
                                           ID = c("ID", "ID", "ID", "ID"), cells ="T-AB")

alluv_plot1 <- compareClonotypes(combined_pt1012_CD8, samples = c("weekA1_ID", "weekB2_ID", "weekC4_ID", "weekD14_ID"), 
                  cloneCall="aa", graph = "alluvial", clonotypes = c("CAMSVLSGNTGKLIF_CASSVASLTVGTGELFF",
                                                                     "CAVGGDDAGNMLTF_CASMPSTGPNEKLFF",
                                                                     "CALRAPYGGSQGNLIF_CASSLPGHEKLFF",
                                                                     "CAVNNNFNKFYF_CASSFGGLDAEAFF",
                                                                     "CAPSVGATNKLIF_CSVELPGDTIYF",
                                                                     "CVVSDRGSTLGRLYF_CASSEPPAGTGAEKLFF",
                                                                     "CAVPRFSDGQKLLF_CASSFGRNQPQHF",
                                                                     "CAASDNTGNQFYF_CASSLAQEVNQPQHF",
                                                                     "CAASGGADGLTF_CASSSASSGEETQYF",
                                                                     "CALSEAQAAGNKLTF_CASSDEDRVEKAFF")) + 
  ggtitle("Clonotype dynamics - pt1012 - CD8") + scale_x_discrete(labels=c("Week 1","Week 2","Week 4","Week 14")) +
  scale_fill_manual(values = col_pal) + theme(legend.title=element_text(size=50), legend.text=element_text(size=50), 
                                              plot.title = element_text(size = 50), axis.title = element_text(size = 50),
                                              axis.text.x = element_text(size = 40), axis.text.y = element_text(size = 40))


#pt1013

s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w14_filtered_contig_annotations.csv")

#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSPGGQETQYF','CASGTGSYEQYF','CASGTGSYEQYF','CASSQIRGGYEQYF','CASSIVGQINTEAFF')
toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

object_pt1013 <- SetIdent(object_pt1013, value = "cell_classification_final")
object_pt1013 <- RenameIdents(object_pt1013, "B cells" = "Other", "T/monocyte doublets" = "Other", "CD3+; CD56+" = "Other", 
                              "CD4 memory" = "CD4", "CD4 naïve" = "CD4", "CD8 memory/ TEMRA" = "CD8", "CD8 naïve/memory" = "CD8", 
                              "CD8 naïve" = "CD8", "gamma/delta" = "Other", "Proliferating T" = "Other","Treg" = "Other")
object_pt1013[["CD4vsCD8.ids"]] <- Idents(object_pt1013)

object_pt1013_CD4 <- subset(x = object_pt1013, subset = CD4vsCD8.ids == "CD4")
object_pt1013_CD8 <- subset(x = object_pt1013, subset = CD4vsCD8.ids == "CD8")

pt1013_CD4_barcodes <- WhichCells(object_pt1013_CD4, idents = "CD4")
pt1013_CD8_barcodes <- WhichCells(object_pt1013_CD8, idents = "CD8")

s1_sub_pt1013_CD8 <- s1_sub[s1_sub$barcode %in% pt1013_CD8_barcodes,]
s2_sub_pt1013_CD8 <- s2_sub[s2_sub$barcode %in% pt1013_CD8_barcodes,]
s3_sub_pt1013_CD8 <- s3_sub[s3_sub$barcode %in% pt1013_CD8_barcodes,]
s4_sub_pt1013_CD8 <- s4_sub[s4_sub$barcode %in% pt1013_CD8_barcodes,]

contig_list_pt1013_CD8 <- list(s1_sub_pt1013_CD8, s2_sub_pt1013_CD8, s3_sub_pt1013_CD8, s4_sub_pt1013_CD8)

combined_pt1013_CD8 <- combineTCR(contig_list_pt1013_CD8, 
                                  samples = c("weekA1", "weekB2", "weekC4", "weekD14"), 
                                  ID = c("ID", "ID", "ID", "ID"), cells ="T-AB")

alluv_plot2 <- compareClonotypes(combined_pt1013_CD8, samples = c("weekA1_ID", "weekB2_ID", "weekC4_ID", "weekD14_ID"), 
                  cloneCall="aa", graph = "alluvial", clonotypes = c("CAVSGIKAAGNKLTF_CASSVAALDQPQHF",
                                                                     "CAVSTSGGSYIPTF_CASSHPTSGRETQYF",
                                                                     "CIVRVLPSNTGKLIF_CASSQDSVGRSSYEQYF",
                                                                     "NA_CASSLDGVSTDTQYF",
                                                                     "CAGQLGGGSNYKLTF_CASSERGVMVSSNQPQHF",
                                                                     "CAMREGLIKAAGNKLTF_CASSVRSRGDDSPLHF",
                                                                     "CALINDYKLSF_CASSQTLWGSGELTDTQYF",
                                                                     "CACQDGGGADGLTF_CASSVKISTDTQYF",
                                                                     "CALDMRQGGKLIF_CASRLHTGTGTSGANVLTF",
                                                                     "CALSEARNAGNMLTF_CASSLTGQGDYEQYF")) + 
  ggtitle("Clonotype dynamics - pt1013 - CD8") + scale_x_discrete(labels=c("Week 1","Week 2","Week 4","Week 14")) +
  scale_fill_manual(values = col_pal) + theme(legend.title=element_text(size=50), legend.text=element_text(size=50), 
                                           plot.title = element_text(size = 50), axis.title = element_text(size = 50),
                                           axis.text.x = element_text(size = 40), axis.text.y = element_text(size = 40))


#pt1014

s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w4_filtered_contig_annotations.csv")
#s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w14_filtered_contig_annotations.csv")


#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
# s4$is_cell = toupper(s4$is_cell)
# s4$high_confidence = toupper(s4$high_confidence)
# s4$full_length = toupper(s4$full_length)
# s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSPPGGSYEQYF','CAISESKGNYGYTF','CASSLTGTGSYNEQFF','CASSPPGQVNNEQFF','CAISTGNTEAFF','CAISTGNTEAFF')

toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

#toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
#s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

object_pt1014 <- SetIdent(object_pt1014, value = "cell_classification_final")
object_pt1014 <- RenameIdents(object_pt1014, "B cells" = "Other", "T/monocyte doublets" = "Other", "CD3+; CD56+" = "Other", 
                              "CD4 memory" = "CD4", "CD4 naïve" = "CD4", "CD8 memory/ TEMRA" = "CD8", "CD8 naïve/memory" = "CD8", 
                              "CD8 naïve" = "CD8", "gamma/delta" = "Other", "Proliferating T" = "Other","Treg" = "Other")
object_pt1014[["CD4vsCD8.ids"]] <- Idents(object_pt1014)

object_pt1014_CD4 <- subset(x = object_pt1014, subset = CD4vsCD8.ids == "CD4")
object_pt1014_CD8 <- subset(x = object_pt1014, subset = CD4vsCD8.ids == "CD8")

pt1014_CD4_barcodes <- WhichCells(object_pt1014_CD4, idents = "CD4")
pt1014_CD8_barcodes <- WhichCells(object_pt1014_CD8, idents = "CD8")

s1_sub_pt1014_CD8 <- s1_sub[s1_sub$barcode %in% pt1014_CD8_barcodes,]
s2_sub_pt1014_CD8 <- s2_sub[s2_sub$barcode %in% pt1014_CD8_barcodes,]
s3_sub_pt1014_CD8 <- s3_sub[s3_sub$barcode %in% pt1014_CD8_barcodes,]
#s4_sub_pt1014_CD8 <- s4_sub[s4_sub$barcode %in% pt1014_CD8_barcodes,]

contig_list_pt1014_CD8 <- list(s1_sub_pt1014_CD8, s2_sub_pt1014_CD8, s3_sub_pt1014_CD8)

combined_pt1014_CD8 <- combineTCR(contig_list_pt1014_CD8, 
                                  samples = c("weekA1", "weekB2", "weekC4"), 
                                  ID = c("ID", "ID", "ID"), cells ="T-AB")

alluv_plot3 <- compareClonotypes(combined_pt1014_CD8, samples = c("weekA1_ID", "weekB2_ID", "weekC4_ID"), 
                  cloneCall="aa", graph = "alluvial", clonotypes = c("CALSEDNFNKFYF_CASSLDPGKGEPLHF",
                                                                     "CAASVNDYKLSF_CASSLELAETGELFF",
                                                                     "CAGRPSNNFGNEKLTF_CASSLGGRNIQYF",
                                                                     "CAVPTGGADGLTF_CATSAGRGELFF",
                                                                     "NA_CASSLELAETGELFF",
                                                                     "NA_CASSLGGRNIQYF",
                                                                     "CAVPGGGADGLTF_CATSWGRGEQFF",
                                                                     "CALRSNRDDKIIF_CSARDGDRGSYEQYF",
                                                                     "CAAYEGTDKLIF_CASRQPSSYEQYF",
                                                                     "CARASGGSYIPTF_CASSDTGTGRDTQYF")) + 
  ggtitle("Clonotype dynamics - pt1014 - CD8") + scale_x_discrete(labels=c("Week 1","Week 2","Week 4")) +
  scale_fill_manual(values = col_pal) + theme(legend.title=element_text(size=50), legend.text=element_text(size=50), 
                                              plot.title = element_text(size = 50), axis.title = element_text(size = 50),
                                              axis.text.x = element_text(size = 40), axis.text.y = element_text(size = 40))


#pt1019

s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w1_filtered_contig_annotations.csv")
#s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w14_filtered_contig_annotations.csv")

#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
# s2$is_cell = toupper(s2$is_cell)
# s2$high_confidence = toupper(s2$high_confidence)
# s2$full_length = toupper(s2$full_length)
# s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSSRLAGSYNEQFF','CASSYTETQYF','CASSYTETQYF','CASSYTETQYF','CASSYTETQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF')

toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

#toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
#s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

object_pt1019 <- SetIdent(object_pt1019, value = "cell_classification_final")
object_pt1019 <- RenameIdents(object_pt1019, "B cells" = "Other", "T/monocyte doublets" = "Other", "CD3+; CD56+" = "Other", 
                              "CD4 memory" = "CD4", "CD4 naïve" = "CD4", "CD8 memory/ TEMRA" = "CD8", "CD8 naïve/memory" = "CD8", 
                              "CD8 naïve" = "CD8", "gamma/delta" = "Other", "Proliferating T" = "Other","Treg" = "Other")
object_pt1019[["CD4vsCD8.ids"]] <- Idents(object_pt1019)

object_pt1019_CD4 <- subset(x = object_pt1019, subset = CD4vsCD8.ids == "CD4")
object_pt1019_CD8 <- subset(x = object_pt1019, subset = CD4vsCD8.ids == "CD8")

pt1019_CD4_barcodes <- WhichCells(object_pt1019_CD4, idents = "CD4")
pt1019_CD8_barcodes <- WhichCells(object_pt1019_CD8, idents = "CD8")

s1_sub_pt1019_CD8 <- s1_sub[s1_sub$barcode %in% pt1019_CD8_barcodes,]
#s2_sub_pt1019_CD8 <- s2_sub[s2_sub$barcode %in% pt1019_CD8_barcodes,]
s3_sub_pt1019_CD8 <- s3_sub[s3_sub$barcode %in% pt1019_CD8_barcodes,]
s4_sub_pt1019_CD8 <- s4_sub[s4_sub$barcode %in% pt1019_CD8_barcodes,]

contig_list_pt1019_CD8 <- list(s1_sub_pt1019_CD8, s3_sub_pt1019_CD8, s4_sub_pt1019_CD8)

combined_pt1019_CD8 <- combineTCR(contig_list_pt1019_CD8, 
                                  samples = c("weekA1", "weekC4", "weekD14"), 
                                  ID = c("ID", "ID", "ID"), cells ="T-AB")

alluv_plot4 <- compareClonotypes(combined_pt1019_CD8, samples = c("weekA1_ID", "weekC4_ID", "weekD14_ID"), 
                  cloneCall="aa", graph = "alluvial", clonotypes = c("NA", "CVVSELNNARLMF;CVFRRGQKLLF_CASSIVGVLSGELFF",
                                                                     "CAMSGRWAGGGNKLTF_CAWSPTESPEQFF",
                                                                     "CAARLGGANNLFF_CASSRTGVVTEAFF",
                                                                     "CAVRDRDYGGATNKLIF_CASSPQPRTGLDGYTF",
                                                                     "CAGPRRYNFNKFYF_CASSSPGLAGGPGNTGELFF",
                                                                     "CAAPQPGYSTLTF_CASSLRGGSTEAFF",
                                                                     "CAYQGGKLIF_CASSSREAGTYNEQFF",
                                                                     "CAMRNRDDKIIF_CASSDWAQAVAAYEQYF",
                                                                     "CAVADSGANSKLTF_CASSFLGTDTQYF",
                                                                     "NA_CASSPQPRTGLDGYTF")) + 
  ggtitle("Clonotype dynamics - pt1019 - CD8") + scale_x_discrete(labels=c("Week 1","Week 4","Week 14")) +
  scale_fill_manual(values = col_pal) + theme(legend.title=element_text(size=50), legend.text=element_text(size=50), 
                                              plot.title = element_text(size = 50), axis.title = element_text(size = 50),
                                              axis.text.x = element_text(size = 40), axis.text.y = element_text(size = 40))


Alluv_main <- ggarrange(alluv_plot1+NoLegend(), alluv_plot2+NoLegend(), alluv_plot3+NoLegend(), alluv_plot4+NoLegend(), ncol=2, nrow=2)

################# Panel 5C - Proportion of cell types (T cells) by time, split by patient ###################
#pt1012
cell_timepoint_pt1012 <- table(object_pt1012$cell_classification_final, object_pt1012$timepoint)
cell_timepoint_pt1012_sub <- cell_timepoint_pt1012[!(row.names(cell_timepoint_pt1012) %in% c('B cells', 'T/monocyte doublets')),]
prop_cell_timepoint_pt1012 <- data.frame(prop.table(cell_timepoint_pt1012_sub, margin = 2))

prop_plot1 <- ggplot(data=subset(prop_cell_timepoint_pt1012, !is.na(Freq)), aes(x=Var2, y=Freq, group=Var1)) + ggtitle("pt1012") +
  geom_line(aes(color=Var1),size = 3) + geom_point(aes(color=Var1), size=6) + xlab("Timepoint") + ylab("Percent") + 
  scale_y_continuous(labels = scales::percent) + theme_light() +
  scale_color_manual(values=c(col_pal[3],col_pal[4],col_pal[5],col_pal[6],col_pal[7],col_pal[8],col_pal[9],
                              col_pal[10],col_pal[11])) + theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title=element_text(size=30), legend.text=element_text(size=30), 
        plot.title = element_text(size = 30), axis.title = element_text(size = 30), axis.text = element_text(size = 20)) + 
  guides(color=guide_legend(title="Cell types"))

#pt1013
cell_timepoint_pt1013 <- table(object_pt1013$cell_classification_final, object_pt1013$timepoint)
cell_timepoint_pt1013_sub <- cell_timepoint_pt1013[!(row.names(cell_timepoint_pt1013) %in% c('B cells', 'T/monocyte doublets')),]
prop_cell_timepoint_pt1013 <- data.frame(prop.table(cell_timepoint_pt1013_sub, margin = 2))

prop_plot2 <- ggplot(data=subset(prop_cell_timepoint_pt1013, !is.na(Freq)), aes(x=Var2, y=Freq, group=Var1)) + ggtitle("pt1013") +
  geom_line(aes(color=Var1),size = 3) + geom_point(aes(color=Var1), size=6) + xlab("Timepoint") + ylab("Percent") + 
  scale_y_continuous(labels = scales::percent) + theme_light() +
  scale_color_manual(values=c(col_pal[3],col_pal[4],col_pal[5],col_pal[6],col_pal[7],col_pal[8],col_pal[9],
                              col_pal[10],col_pal[11])) + theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title=element_text(size=30), legend.text=element_text(size=30), 
        plot.title = element_text(size = 30), axis.title = element_text(size = 30), axis.text = element_text(size = 20)) + 
  guides(color=guide_legend(title="Cell types"))

#pt1014
cell_timepoint_pt1014 <- table(object_pt1014$cell_classification_final, object_pt1014$timepoint)
cell_timepoint_pt1014_sub <- cell_timepoint_pt1014[!(row.names(cell_timepoint_pt1014) %in% c('B cells', 'T/monocyte doublets')),]
prop_cell_timepoint_pt1014 <- data.frame(prop.table(cell_timepoint_pt1014_sub, margin = 2))

prop_plot3 <- ggplot(data=subset(prop_cell_timepoint_pt1014, !is.na(Freq)), aes(x=Var2, y=Freq, group=Var1)) + ggtitle("pt1014") +
  geom_line(aes(color=Var1),size = 3) + geom_point(aes(color=Var1), size=6) + xlab("Timepoint") + ylab("Percent") + 
  scale_y_continuous(labels = scales::percent) + theme_light() +
  scale_color_manual(values=c(col_pal[3],col_pal[4],col_pal[5],col_pal[6],col_pal[7],col_pal[8],col_pal[9],
                              col_pal[10],col_pal[11])) + theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title=element_text(size=30), legend.text=element_text(size=30), 
        plot.title = element_text(size = 30), axis.title = element_text(size = 30), axis.text = element_text(size = 20)) + 
  guides(color=guide_legend(title="Cell types"))

#pt1019
cell_timepoint_pt1019 <- table(object_pt1019$cell_classification_final, object_pt1019$timepoint)
cell_timepoint_pt1019_sub <- cell_timepoint_pt1019[!(row.names(cell_timepoint_pt1019) %in% c('B cells', 'T/monocyte doublets')),]
prop_cell_timepoint_pt1019 <- data.frame(prop.table(cell_timepoint_pt1019_sub, margin = 2))

prop_plot4 <- ggplot(data=subset(prop_cell_timepoint_pt1019, !is.na(Freq)), aes(x=Var2, y=Freq, group=Var1)) + ggtitle("pt1019") +
  geom_line(aes(color=Var1),size = 3) + geom_point(aes(color=Var1), size=6) + xlab("Timepoint") + ylab("Percent") + 
  scale_y_continuous(labels = scales::percent) + theme_light() +
  scale_color_manual(values=c(col_pal[3],col_pal[4],col_pal[5],col_pal[6],col_pal[7],col_pal[8],col_pal[9],
                              col_pal[10],col_pal[11])) + theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title=element_text(size=30), legend.text=element_text(size=30), 
        plot.title = element_text(size = 30), axis.title = element_text(size = 30), axis.text = element_text(size = 20)) + 
  guides(color=guide_legend(title="Cell types"))

main_prop_plot <- ggarrange(prop_plot1, prop_plot2, prop_plot3, prop_plot4, ncol=2, nrow=2, common.legend = TRUE, legend = "right")

################# Panel 5D,E - heatmaps split by celltype without week14 ###################
#removing week14
obj_rm_w14 <- object
obj_rm_w14 <- SetIdent(obj_rm_w14, value = "timepoint")
obj_rm_w14 <- subset(x = obj_rm_w14, idents = c("week1", "week2", "week4"))
obj_rm_w14 <- SetIdent(obj_rm_w14, value = "cell_classification_final")

adt_data <- GetAssayData(obj_rm_w14, slot = "data", assay = "ADT_denoised_iso_quant")
obj_rm_w14[["ADT_denoised_iso_quant"]] = CreateAssayObject(counts = adt_data)

genelist_main <- c("NR4A2", "JUN", "JUNB", "JUND", "FOSB", "FOS", 
              "IL7R", "CCR7", "SMAD7", 
              "LTB", "CISH", "TIGIT", "PRF1", "GZMA", "IFNG", "TNF", "GZMB", "CCL4", "CCL5", "XCL2", "CD74", "HLA-DRA", "HLA-DRB1")
genelist_adt <- c("CD45RO-ADT", "CD45RA-ADT", "CD4-ADT", "CD8-ADT", "CD56-ADT")
genelist_class <- c("CD4", "CISH", "CD8A", "CD8B", "GZMH", "NKG7", "IL7R", "GZMB", "GATA3", "CCR7", "CTLA4", "FOXP3",
                    "MKI67", "GNLY", "FCER1G", "CD79A", "JCHAIN", "IGKC", "CD3D", "LYZ", "S100A8", "S100A9")
# genelist_supp <- c("CD45RO-ADT", "CD45RA-ADT", "CD4-ADT", "CD8-ADT", "CD56-ADT",
#                    "CD4", "CISH", "CD8A", "CD8B", "GZMH", "NKG7", "IL7R", "GZMB", "GATA3", "CCR7", "CTLA4", "FOXP3",
#                    "MKI67", "GNLY", "FCER1G", "CD79A", "JCHAIN", "IGKC", "CD3D", "LYZ", "S100A8", "S100A9")

sample_order <- data.frame(Patient = c("1012", "1012", "1012", "1013", "1013", "1013", "1014", "1014", "1014", "1019", "1019"), 
                     Timepoint = c("Week 1", "Week 2", "Week 4", "Week 1", "Week 2", "Week 4", "Week 1", "Week 2", "Week 4", "Week 1", "Week 4"))

pheatmap_order <- c("pt1012_week1", "pt1013_week1", "pt1014_week1", "pt1019_week1", "pt1012_week2", "pt1013_week2", 
                    "pt1014_week2","pt1012_week4", "pt1013_week4", "pt1014_week4", "pt1019_week4")

ann_colors = list(Patient = c(`1012` = col_pal[1], `1013` = col_pal[2], `1014`= col_pal[3], `1019`= col_pal[4]), 
                  Timepoint = c(`Week 1` = col_pal[5], `Week 2` = col_pal[6], `Week 4` = col_pal[7]))

run_sig_features <- function(seurat_obj, cell_idents, obj_name, genelist_main) {
obj_subset <- subset(seurat_obj, idents = cell_idents)
obj_subset <- SetIdent(obj_subset, value = "timepoint")
obj_w1_w2 <- FindMarkers(obj_subset, ident.1 = "week2", ident.2 = "week1", min.pct = 0.25)
obj_w1_w4 <- FindMarkers(obj_subset, ident.1 = "week4", ident.2 = "week1", min.pct = 0.25)
obj_w2_w4 <- FindMarkers(obj_subset, ident.1 = "week4", ident.2 = "week2", min.pct = 0.25)

obj_w1_w2_vector <- rownames(obj_w1_w2[obj_w1_w2$p_val_adj < 0.05, ])
obj_w1_w4_vector <- rownames(obj_w1_w4[obj_w1_w4$p_val_adj < 0.05, ])
obj_w2_w4_vector <- rownames(obj_w2_w4[obj_w2_w4$p_val_adj < 0.05, ])
 
obj_siglist <- unique(append(obj_w1_w2_vector, c(obj_w1_w4_vector, obj_w2_w4_vector)))
temp_features <- intersect(obj_siglist, genelist_main)

temp_name <- paste(obj_name, "_features", sep = "")
return(temp_features)
}

CD4memory_features <- run_sig_features(obj_rm_w14, "CD4 memory", "CD4memory", genelist_main)

CD4naive_features <- run_sig_features(obj_rm_w14, "CD4 naïve", "CD4naive", genelist_main)

CD8memory_TEMRA_features <- run_sig_features(obj_rm_w14, "CD8 memory/ TEMRA", "CD8memory_TEMRA", genelist_main)

CD8naive_memory_features <- run_sig_features(obj_rm_w14, "CD8 naïve/memory", "CD8naive_memory", genelist_main)

CD8naive_features <- run_sig_features(seurat_obj = obj_rm_w14, cell_idents = "CD8 naïve", obj_name = "CD8naive", genelist_main)


##CD4 memory
obj_CD4mem <- subset(obj_rm_w14, idents = "CD4 memory")
obj_CD4mem <- SetIdent(obj_CD4mem, value = "sample")
#Setting up the matrix and its order
#CD4mem_main_matrix <- AverageExpression(obj_CD4mem, assays = "RNA", features = genelist_main)
CD4memory_features_reviewed <- c("CISH","IL7R","LTB","JUND","JUNB","CD74","JUN","SMAD7","FOS")

# CD4mem_main_matrix <- AverageExpression(obj_CD4mem, assays = "RNA", features = CD4memory_features)
CD4mem_main_matrix <- AverageExpression(obj_CD4mem, assays = "RNA", features = CD4memory_features_reviewed)
CD4mem_main_matrix <- CD4mem_main_matrix$RNA
row.names(sample_order) <- colnames(CD4mem_main_matrix)
CD4mem_main_matrix <- CD4mem_main_matrix[, pheatmap_order]

pheatmap(CD4mem_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD4 memory main heatmap \n by patient and timepoint")

#row clustering on
pheatmap(CD4mem_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = TRUE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD4 memory main heatmap \n by patient and timepoint")


#Setting up the matrix and its order
CD4mem_supp_matrix_RNA <- AverageExpression(obj_CD4mem, features = genelist_class, assays = "RNA")$RNA
CD4mem_supp_matrix_ADT <- AverageExpression(obj_CD4mem, features = genelist_adt, assays = "ADT_denoised_iso_quant", slot = "counts")$ADT_denoised_iso_quant

CD4mem_supp_matrix <- rbind(CD4mem_supp_matrix_ADT, CD4mem_supp_matrix_RNA)
row.names(sample_order) <- colnames(CD4mem_supp_matrix)
CD4mem_supp_matrix <- CD4mem_supp_matrix[, pheatmap_order]

pheatmap(CD4mem_supp_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE,
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order,
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD4 memory supplementary heatmap \n by patient and timepoint")

#Dotplot
DotPlot(obj_CD4mem, features = genelist_main, group.by = 'sample') + 
  coord_flip() + ggtitle("CD4 memory - main genes")

##CD4 naive

obj_CD4naive <- subset(obj_rm_w14, idents = "CD4 naïve")
obj_CD4naive <- SetIdent(obj_CD4naive, value = "sample")

#Setting up the matrix and its order
#CD4naive_main_matrix <- AverageExpression(obj_CD4naive, assays = "RNA", features = genelist_main)
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD4naive_features_reviewed <- c("LTB","IL7R","CISH","CCR7","JUNB","CD74","JUND","JUN","SMAD7","NR4A2")

# CD4naive_main_matrix <- AverageExpression(obj_CD4naive, assays = "RNA", features = CD4naive_features)
CD4naive_main_matrix <- AverageExpression(obj_CD4naive, assays = "RNA", features = CD4naive_features_reviewed)

CD4naive_main_matrix <- CD4naive_main_matrix$RNA
row.names(sample_order) <- colnames(CD4naive_main_matrix)
CD4naive_main_matrix <- CD4naive_main_matrix[, pheatmap_order]

pheatmap(CD4naive_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD4 naïve main heatmap \n by patient and timepoint")

#row clustering on
pheatmap(CD4naive_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = TRUE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD4 naïve main heatmap \n by patient and timepoint")
#Dotplot
DotPlot(obj_CD4naive, features = genelist_main, group.by = 'sample') + 
  coord_flip() + ggtitle("CD4 naive - main genes")

##CD8 memory and TEMRA

obj_CD8mem_TEMRA <- subset(obj_rm_w14, idents = "CD8 memory/ TEMRA")
obj_CD8mem_TEMRA <- SetIdent(obj_CD8mem_TEMRA, value = "sample")

#Setting up the matrix and its order
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD8mem_TEMRA_features_reviewed <- c("GZMA","PRF1","CISH","CCL4","JUND","GZMB","IL7R","JUN","LTB","XCL2",    
                                    "TIGIT","TNF","IFNG","SMAD7","CD74","HLA-DRA","FOSB")

# CD8mem_TEMRA_main_matrix <- AverageExpression(obj_CD8mem_TEMRA, assays = "RNA", features = genelist_main)
# CD8mem_TEMRA_main_matrix <- AverageExpression(obj_CD8mem_TEMRA, assays = "RNA", features = CD8memory_TEMRA_features)
CD8mem_TEMRA_main_matrix <- AverageExpression(obj_CD8mem_TEMRA, assays = "RNA", features = CD8mem_TEMRA_features_reviewed)

CD8mem_TEMRA_main_matrix <- CD8mem_TEMRA_main_matrix$RNA
row.names(sample_order) <- colnames(CD8mem_TEMRA_main_matrix)
CD8mem_TEMRA_main_matrix <- CD8mem_TEMRA_main_matrix[, pheatmap_order]

pheatmap(CD8mem_TEMRA_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 memory/ TEMRA main heatmap \n by patient and timepoint")

#row clustering on
pheatmap(CD8mem_TEMRA_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = TRUE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 memory/ TEMRA main heatmap \n by patient and timepoint")

#Dotplot
DotPlot(obj_CD8mem_TEMRA, features = genelist_main, group.by = 'sample') + 
  coord_flip() + ggtitle("CD8 memory/TEMRA - main genes")

# #Setting up the matrix and its order
# CD8mem_TEMRA_supp_matrix_RNA <- AverageExpression(obj_CD8mem_TEMRA, features = genelist_class, assays = "RNA")$RNA
# CD8mem_TEMRA_supp_matrix_ADT <- AverageExpression(obj_CD8mem_TEMRA, features = genelist_adt, assays = "ADT_denoised_iso_quant", slot = "counts")$ADT
# 
# CD8mem_TEMRA_supp_matrix <- rbind(CD8mem_TEMRA_supp_matrix_ADT, CD8mem_TEMRA_supp_matrix_RNA)
# row.names(sample_order) <- colnames(CD8mem_TEMRA_supp_matrix)
# CD8mem_TEMRA_supp_matrix <- CD8mem_TEMRA_supp_matrix[, pheatmap_order]
# 
# pheatmap(CD8mem_TEMRA_supp_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
#          clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
#          show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
#          annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 memory/ TEMRA supplementary heatmap \n by patient and timepoint")



##CD8 naive and memory

obj_CD8naive_mem <- subset(obj_rm_w14, idents = "CD8 naïve/memory")
obj_CD8naive_mem <- SetIdent(obj_CD8naive_mem, value = "sample")

#Setting up the matrix and its order
# CD8naive_mem_main_matrix <- AverageExpression(obj_CD8naive_mem, assays = "RNA", features = genelist_main)
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD8naive_mem_features_reviewed <- c("IL7R","JUNB","CD74","JUND","GZMA","PRF1","LTB","CISH","TIGIT","HLA-DRB1",
                                    "JUN","CCL4","HLA-DRA","FOS","NR4A2","CCR7","SMAD7")

# CD8naive_mem_main_matrix <- AverageExpression(obj_CD8naive_mem, assays = "RNA", features = CD8naive_memory_features)
CD8naive_mem_main_matrix <- AverageExpression(obj_CD8naive_mem, assays = "RNA", features = CD8naive_mem_features_reviewed)

CD8naive_mem_main_matrix <- CD8naive_mem_main_matrix$RNA
row.names(sample_order) <- colnames(CD8naive_mem_main_matrix)
CD8naive_mem_main_matrix <- CD8naive_mem_main_matrix[, pheatmap_order]

pheatmap(CD8naive_mem_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 naïve/memory main heatmap \n by patient and timepoint")

#row clustering on
pheatmap(CD8naive_mem_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = TRUE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 naïve/memory main heatmap \n by patient and timepoint")

#Dotplot
DotPlot(obj_CD8naive_mem, features = genelist_main, group.by = 'sample') + 
  coord_flip() + ggtitle("CD8 naive/memory - main genes")

##CD8 naive

obj_CD8naive <- subset(obj_rm_w14, idents = "CD8 naïve/memory")
obj_CD8naive <- SetIdent(obj_CD8naive, value = "sample")

#Setting up the matrix and its order
# CD8naive_main_matrix <- AverageExpression(obj_CD8naive, assays = "RNA", features = genelist_main)
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD8naive_features_reviewed <- c()

CD8naive_main_matrix <- AverageExpression(obj_CD8naive, assays = "RNA", features = CD8naive_features)
# CD8naive_main_matrix <- AverageExpression(obj_CD8naive, assays = "RNA", features = CD8naive_features_reviewed)


CD8naive_main_matrix <- CD8naive_main_matrix$RNA
row.names(sample_order) <- colnames(CD8naive_main_matrix)
CD8naive_main_matrix <- CD8naive_main_matrix[, pheatmap_order]

pheatmap(CD8naive_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 naïve main heatmap \n by patient and timepoint")

#row clustering on
pheatmap(CD8naive_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = TRUE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 naïve main heatmap \n by patient and timepoint")

#Dotplot
DotPlot(obj_CD8naive, features = genelist_main, group.by = 'sample') + 
  coord_flip() + ggtitle("CD8 naive - main genes")

# #Setting up the matrix and its order
# CD8naive_mem_supp_matrix_RNA <- AverageExpression(obj_CD8naive_mem, features = genelist_class, assays = "RNA")$RNA
# CD8naive_mem_supp_matrix_ADT <- AverageExpression(obj_CD8naive_mem, features = genelist_adt, assays = "ADT_denoised_iso_quant", slot = "counts")$ADT
# 
# CD8naive_mem_supp_matrix <- rbind(CD8naive_mem_supp_matrix_ADT, CD8naive_mem_supp_matrix_RNA)
# row.names(sample_order) <- colnames(CD8naive_mem_supp_matrix)
# CD8naive_mem_supp_matrix <- CD8naive_mem_supp_matrix[, pheatmap_order]
# 
# pheatmap(CD8naive_mem_supp_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
#          clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
#          show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
#          annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "CD8 naïve/memory supplementary heatmap \n by patient and timepoint")


#Heatmap Supplementary by celltype
object <- SetIdent(object, value = "cell_classification_final")

sample_order_supp <- data.frame(CellType = c("B cells", "T/monocyte doublets", "CD3+; CD56+", "CD4 memory", "CD4 naïve", "CD8 memory/ TEMRA", 
                                             "CD8 naïve/memory", "CD8 naïve", "gamma/delta", "Proliferating T", "Treg"))

pheatmap_order_supp <- c("CD4 memory", "CD4 naïve", "CD8 memory/ TEMRA", "CD8 naïve/memory", "CD8 naïve", "gamma/delta", 
                    "Proliferating T", "Treg", "T/monocyte doublets", "CD3+; CD56+", "B cells")  

ann_colors_supp = list(CellType = c(`B cells` = "light grey", `T/monocyte doublets` = "dark grey", `CD3+; CD56+`= col_pal[3], `CD4 memory`= col_pal[4],
                                    `CD4 naïve` = col_pal[5], `CD8 memory/ TEMRA` = col_pal[6], `CD8 naïve/memory`= col_pal[7], `CD8 naïve`= col_pal[8],
                                    `gamma/delta` = col_pal[9], `Proliferating T` = col_pal[10], `Treg`= col_pal[11]))

object_adt <- object

adt_data_whole <- GetAssayData(object_adt, slot = "data", assay = "ADT_denoised_iso_quant")
object_adt[["ADT_denoised_iso_quant"]] = CreateAssayObject(counts = adt_data_whole)

#Setting up the matrix and its order
object_adt_matrix_RNA <- AverageExpression(object_adt, features = genelist_class, assays = "RNA")$RNA
object_adt_matrix_ADT <- AverageExpression(object_adt, features = genelist_adt, assays = "ADT_denoised_iso_quant", slot = "counts")$ADT_denoised_iso_quant

object_adt_matrix <- rbind(object_adt_matrix_ADT,object_adt_matrix_RNA)

row.names(sample_order_supp) <- colnames(object_adt_matrix)
object_adt_matrix <- object_adt_matrix[, pheatmap_order_supp]

supp_celltypes_heatmap <- pheatmap(object_adt_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE,
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order_supp,
         labels_row = make_italic_names(object_adt_matrix, rownames, rownames(object_adt_matrix)),
         annotation_colors = ann_colors_supp, fontsize = 16, cellwidth = 30, cellheight = 15)
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/Heatmap_celltypes_supp_0713.pdf", supp_celltypes_heatmap, width = 20, height = 20, units = "in")
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/Heatmap_celltypes_supp_0713.png", supp_celltypes_heatmap, width = 10, height = 10, units = "in")

#row clustering on
pheatmap(object_adt_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = TRUE,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE,
         show_colnames = TRUE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order_supp,
         annotation_colors = ann_colors_supp, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "Cell classification heatmap")


##Heatmap split by patient 
object_pt1014 <- SetIdent(object_pt1014, value = "cell_classification_final")
obj <- subset(object_pt1014, idents = "CD8 naïve/memory")
obj <- SetIdent(obj, value = "timepoint")

sample_order <- data.frame(Timepoint = c("Week 1", "Week 2", "Week 4"))
ann_colors = list(Timepoint = c(`Week 1` = col_pal[5], `Week 2` = col_pal[6], `Week 4` = col_pal[7]))

#Setting up the matrix and its order
main_matrix <- AverageExpression(obj, assays = "RNA", features = genelist_main)
main_matrix <- main_matrix$RNA
row.names(sample_order) <- colnames(main_matrix)
#CD4mem_main_matrix <- CD4mem_main_matrix[, pheatmap_order]

pheatmap(main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "pt1014 CD8 naïve/memory main heatmap by timepoint")

##Make a function for it!!!
CD4memory_features <- c("CISH", "IL7R", "LTB" , "JUND", "JUNB", "CD74", "JUN", "SMAD7", 
                        "FOSB", "FOS")
CD4naive_features <- c("LTB", "IL7R", "CISH", "CCR7", "JUNB", "CD74", "JUND", 
                       "JUN", "SMAD7", "NR4A2")
CD8memory_TEMRA_features <- c("LTB", "XCL2", "TIGIT", "TNF", "IFNG", "JUNB", 
                              "SMAD7", "CD74", "HLA-DRA", "FOSB")
CD8naive_memory_features <- c("IL7R", "JUNB", "CD74", "JUND", "GZMA", "PRF1", "LTB", 
                              "CISH", "TIGIT", "HLA-DRB1", "FOSB", "JUN", "CCL4", 
                              "HLA-DRA", "FOS", "CCL5", "NR4A2", "CCR7", "SMAD7")

run_patient_pheatmap <- function(seurat_obj, patient, cell_idents, features_to_plot, plot_title) {
  temp_obj <- SetIdent(seurat_obj, value = "patient")
  temp_obj <- subset(temp_obj, idents = patient)
  
  temp_obj <- SetIdent(temp_obj, value = "cell_classification_final")
  temp_obj <- subset(temp_obj, idents = cell_idents)
  temp_obj <- SetIdent(temp_obj, value = "timepoint")
  
  if (patient == "pt1019") {
    sample_order <- data.frame(Timepoint = c("Week 1", "Week 4"))
    ann_colors = list(Timepoint = c(`Week 1` = col_pal[5], `Week 4` = col_pal[7]))
  } else {
  sample_order <- data.frame(Timepoint = c("Week 1", "Week 2", "Week 4"))
  ann_colors = list(Timepoint = c(`Week 1` = col_pal[5], `Week 2` = col_pal[6], `Week 4` = col_pal[7]))
  }
  #Setting up the matrix and its order
  temp_main_matrix <- AverageExpression(temp_obj, assays = "RNA", features = features_to_plot)
  temp_main_matrix <- temp_main_matrix$RNA
  row.names(sample_order) <- colnames(temp_main_matrix)
  
  temp_title <- paste(plot_title, patient, sep = "_")
  temp_plot <- pheatmap(temp_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
           clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
           show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
           annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = temp_title)
  
  ggsave(paste("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/pheatmaps/",temp_title, ".png", sep = ""), temp_plot, width = 7, height = 3, units = "in")
}

patient_list <- c("pt1012", "pt1013", "pt1014", "pt1019")

for (i in 1:length(patient_list)) {
  run_patient_pheatmap(obj_rm_w14, patient_list[i], "CD4 memory", CD4memory_features, "CD4_memory_main_heatmap")
  run_patient_pheatmap(obj_rm_w14, patient_list[i], "CD4 naïve", CD4naive_features, "CD4_naive_main_heatmap")
  run_patient_pheatmap(obj_rm_w14, patient_list[i], "CD8 memory/ TEMRA", CD8memory_TEMRA_features, "CD8_memory_TEMRA_main_heatmap")
  run_patient_pheatmap(obj_rm_w14, patient_list[i], "CD8 naïve/memory", CD8naive_memory_features, "CD8_naive_memory_main_heatmap")
  run_patient_pheatmap(obj_rm_w14, patient_list[i], "CD8 naïve", CD8naive_features, "CD8_naive_main_heatmap")
}

#Function to italicize gene names
make_italic_names <- function(mat, rc_fun, rc_names) {
  italic_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        italic_names[i] <<-
        bquote(italic(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  italic_names
}

##### Panel H (final) Heatmap by timepoint combining patient #####
obj_rm_w14 <- SetIdent(obj_rm_w14, value = "cell_classification_final")
obj <- subset(obj_rm_w14, idents = "CD8 naïve/memory")
obj <- SetIdent(obj, value = "timepoint")

sample_order <- data.frame(Timepoint = c("Week 1", "Week 2", "Week 4"))
ann_colors_supp = list(Timepoint = c(`Week 1` = col_pal_extra[1], `Week 2` = col_pal_extra[2], `Week 4` = col_pal_extra[14]))

# #Setting up the matrix and its order
# main_matrix <- AverageExpression(obj, assays = "RNA", features = genelist_main)
# main_matrix <- main_matrix$RNA
# row.names(sample_order) <- colnames(main_matrix)
# #CD4mem_main_matrix <- CD4mem_main_matrix[, pheatmap_order]
# 
# pheatmap(main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
#          clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = TRUE, show_rownames = TRUE, 
#          show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
#          annotation_colors = ann_colors, fontsize_row = 10, cellwidth = 20, cellheight = 9, main = "Combined CD8 naïve/memory main heatmap by timepoint")

##CD4 memory
obj_CD4mem <- subset(obj_rm_w14, idents = "CD4 memory")
obj_CD4mem <- SetIdent(obj_CD4mem, value = "timepoint")

#Setting up the matrix and its order
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD4memory_features_reviewed <- c("IL7R", "LTB", "CISH", "CD74", "JUND", "JUNB", "JUN", "FOS", "SMAD7")

CD4mem_main_matrix <- AverageExpression(obj_CD4mem, assays = "RNA", features = CD4memory_features_reviewed)
CD4mem_main_matrix <- CD4mem_main_matrix$RNA
row.names(sample_order) <- colnames(CD4mem_main_matrix)

CD4mem_main_heatmap <- pheatmap(CD4mem_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = FALSE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors_supp, fontsize = 16, cellwidth = 30, cellheight = 15,
         labels_row = make_italic_names(CD4mem_main_matrix, rownames, rownames(CD4mem_main_matrix)))
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD4mem_main_heatmap_0713.pdf", 
       CD4mem_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD4mem_main_heatmap_0713.png", 
       CD4mem_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)

##CD4 naive

obj_CD4naive <- subset(obj_rm_w14, idents = "CD4 naïve")
obj_CD4naive <- SetIdent(obj_CD4naive, value = "timepoint")

#Setting up the matrix and its order
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD4naive_features_reviewed <- c("IL7R", "LTB", "CISH", "CD74", "CCR7", "JUND", "JUNB", "JUN", "NR4A2", "SMAD7")

CD4naive_main_matrix <- AverageExpression(obj_CD4naive, assays = "RNA", features = CD4naive_features_reviewed)

CD4naive_main_matrix <- CD4naive_main_matrix$RNA
row.names(sample_order) <- colnames(CD4naive_main_matrix)

CD4naive_main_heatmap <- pheatmap(CD4naive_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = FALSE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors_supp, fontsize = 16, cellwidth = 30, cellheight = 15,
         labels_row = make_italic_names(CD4naive_main_matrix, rownames, rownames(CD4naive_main_matrix)))
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD4naive_main_heatmap_0713.pdf", 
       CD4naive_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD4naive_main_heatmap_0713.png", 
       CD4naive_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)

##CD8 memory and TEMRA

obj_CD8mem_TEMRA <- subset(obj_rm_w14, idents = "CD8 memory/ TEMRA")
obj_CD8mem_TEMRA <- SetIdent(obj_CD8mem_TEMRA, value = "timepoint")

#Setting up the matrix and its order
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD8mem_TEMRA_features_reviewed <- c("IL7R", "LTB", "CISH", "TIGIT", "GZMA", "GZMB", "PRF1", "IFNG", "TNF", "CCL4", "XCL2", 
                                    "CD74", "HLA-DRA", "JUND", "JUN", "FOSB", "SMAD7")

CD8mem_TEMRA_main_matrix <- AverageExpression(obj_CD8mem_TEMRA, assays = "RNA", features = CD8mem_TEMRA_features_reviewed)

CD8mem_TEMRA_main_matrix <- CD8mem_TEMRA_main_matrix$RNA
row.names(sample_order) <- colnames(CD8mem_TEMRA_main_matrix)

CD8mem_TEMRA_main_heatmap <- pheatmap(CD8mem_TEMRA_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = FALSE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors_supp, fontsize = 16, cellwidth = 30, cellheight = 15,
         labels_row = make_italic_names(CD8mem_TEMRA_main_matrix, rownames, rownames(CD8mem_TEMRA_main_matrix)))
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD8mem_TEMRA_main_heatmap_0713.pdf", 
       CD8mem_TEMRA_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD8mem_TEMRA_main_heatmap_0713.png", 
       CD8mem_TEMRA_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)

##CD8 naive and memory

obj_CD8naive_mem <- subset(obj_rm_w14, idents = "CD8 naïve/memory")
obj_CD8naive_mem <- SetIdent(obj_CD8naive_mem, value = "timepoint")

#Setting up the matrix and its order
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD8naive_mem_features_reviewed <- c("IL7R", "LTB", "CISH", "TIGIT", "GZMA", "PRF1", "CCL4", "CCR7", "CD74", 
                                    "HLA-DRA", "HLA-DRB1", "JUND", "JUNB", "JUN", "FOS", "NR4A2", "SMAD7")

CD8naive_mem_main_matrix <- AverageExpression(obj_CD8naive_mem, assays = "RNA", features = CD8naive_mem_features_reviewed)

CD8naive_mem_main_matrix <- CD8naive_mem_main_matrix$RNA
row.names(sample_order) <- colnames(CD8naive_mem_main_matrix)

CD8naive_mem_main_heatmap <- pheatmap(CD8naive_mem_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = FALSE, show_rownames = TRUE, 
         show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_col = sample_order, 
         annotation_colors = ann_colors_supp, fontsize = 16, cellwidth = 30, cellheight = 15,
         labels_row = make_italic_names(CD8naive_mem_main_matrix, rownames, rownames(CD8naive_mem_main_matrix)))
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD8naive_mem_main_heatmap_0713.pdf", 
       CD8naive_mem_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD8naive_mem_main_heatmap_0713.png", 
       CD8naive_mem_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)

##CD8 naive

obj_CD8naive <- subset(obj_rm_w14, idents = "CD8 naïve")
obj_CD8naive <- SetIdent(obj_CD8naive, value = "timepoint")

#Setting up the matrix and its order
#for final plot, we only need these features (manually reviewed by Jennifer Foltz)
CD8naive_features_reviewed <- c("IL7R", "LTB", "CISH", "CD74", "CCR7", "JUND", "JUNB", "JUN", "NR4A2", "SMAD7")

CD8naive_main_matrix <- AverageExpression(obj_CD8naive, assays = "RNA", features = CD8naive_features_reviewed)

CD8naive_main_matrix <- CD8naive_main_matrix$RNA
row.names(sample_order) <- colnames(CD8naive_main_matrix)

CD8naive_main_heatmap <- pheatmap(CD8naive_main_matrix, scale = "row", color = rev(brewer.pal(n = 9, name = "RdBu")), cluster_rows = FALSE, cluster_cols = FALSE, 
                                      clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", legend = FALSE, show_rownames = TRUE, 
                                      show_colnames = FALSE, annotation_names_col = TRUE, annotation_names_row = FALSE, annotation_col = sample_order, 
                                      annotation_colors = ann_colors_supp, fontsize = 16, cellwidth = 30, cellheight = 15,
                                      labels_row = make_italic_names(CD8naive_main_matrix, rownames, rownames(CD8naive_main_matrix)))
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD8naive_main_heatmap_0713.pdf", 
       CD8naive_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/CD8naive_main_heatmap_0713.png", 
       CD8naive_main_heatmap, width = 7, height = 5, units = "in", limitsize = FALSE)


#################Panel I - Clonal diversity###################

#pt1012
s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1012_w14_filtered_contig_annotations.csv")



#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)


virus_CDR3 <- list('CASSRDSSNQPQHF','CASSLDGGPYEQYF','CATSIGGNQPQHF','CASSGGINEQFF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSLGRSYEQYF','CASSYGAGGYNEQFF','CASSYGAGGYNEQFF','CASSYGAGGYNEQFF','CASSLIINEQFF','CASSLIINEQFF','CASSLIINEQFF')
toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

contig_list_pt1012 <- list(s1, s2, s3, s4)
contig_list_pt1012_sub <- list(s1_sub, s2_sub, s3_sub, s4_sub)

#pt1013
s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1013_w14_filtered_contig_annotations.csv")

#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSPGGQETQYF','CASGTGSYEQYF','CASGTGSYEQYF','CASSQIRGGYEQYF','CASSIVGQINTEAFF')
toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]


contig_list_pt1013 <- list(s1, s2, s3, s4)
contig_list_pt1013_sub <- list(s1_sub, s2_sub, s3_sub, s4_sub)


#pt1014
s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w1_filtered_contig_annotations.csv")
s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w4_filtered_contig_annotations.csv")
#s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1014_w14_filtered_contig_annotations.csv")


#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
s2$is_cell = toupper(s2$is_cell)
s2$high_confidence = toupper(s2$high_confidence)
s2$full_length = toupper(s2$full_length)
s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
# s4$is_cell = toupper(s4$is_cell)
# s4$high_confidence = toupper(s4$high_confidence)
# s4$full_length = toupper(s4$full_length)
# s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSPPGGSYEQYF','CAISESKGNYGYTF','CASSLTGTGSYNEQFF','CASSPPGQVNNEQFF','CAISTGNTEAFF','CAISTGNTEAFF')

toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

#toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
#s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]


contig_list_pt1014 <- list(s1, s2, s3)
contig_list_pt1014_sub <- list(s1_sub, s2_sub, s3_sub)

#pt1019

s1 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w1_filtered_contig_annotations.csv")
#s2 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w2_filtered_contig_annotations.csv")
s3 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w4_filtered_contig_annotations.csv")
s4 <- read.csv("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/Integrating with Seurat/pt1019_w14_filtered_contig_annotations.csv")

#might be important to change all boolean columns to upper case for later on, not sure
s1$is_cell = toupper(s1$is_cell)
s1$high_confidence = toupper(s1$high_confidence)
s1$full_length = toupper(s1$full_length)
s1$productive = toupper(s1$productive)
# s2$is_cell = toupper(s2$is_cell)
# s2$high_confidence = toupper(s2$high_confidence)
# s2$full_length = toupper(s2$full_length)
# s2$productive = toupper(s2$productive)
s3$is_cell = toupper(s3$is_cell)
s3$high_confidence = toupper(s3$high_confidence)
s3$full_length = toupper(s3$full_length)
s3$productive = toupper(s3$productive)
s4$is_cell = toupper(s4$is_cell)
s4$high_confidence = toupper(s4$high_confidence)
s4$full_length = toupper(s4$full_length)
s4$productive = toupper(s4$productive)

virus_CDR3 <- list('CASSSRLAGSYNEQFF','CASSYTETQYF','CASSYTETQYF','CASSYTETQYF','CASSYTETQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSVAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF','CASSLAGGSYEQYF')

toremove_s1 <- s1[s1$cdr3 %in% virus_CDR3,]
s1_sub <- s1[!s1$barcode %in% toremove_s1$barcode,]

#toremove_s2 <- s2[s2$cdr3 %in% virus_CDR3,]
#s2_sub <- s2[!s2$barcode %in% toremove_s2$barcode,]

toremove_s3 <- s3[s3$cdr3 %in% virus_CDR3,]
s3_sub <- s3[!s3$barcode %in% toremove_s3$barcode,]

toremove_s4 <- s4[s4$cdr3 %in% virus_CDR3,]
s4_sub <- s4[!s4$barcode %in% toremove_s4$barcode,]

contig_list_pt1019 <- list(s1, s3, s4)
contig_list_pt1019_sub <- list(s1_sub, s3_sub, s4_sub)


######## all together clonal diversity analysis ######

combined_pt1012 <- combineTCR(contig_list_pt1012_sub, 
                       samples = c("pt1012", "pt1012", "pt1012", "pt1012"), 
                       ID = c("week1", "week2", "week4", "week14"), cells ="T-AB")
combined_pt1013 <- combineTCR(contig_list_pt1013_sub, 
                       samples = c("pt1013", "pt1013", "pt1013", "pt1013"), 
                       ID = c("week1", "week2", "week4", "week14"), cells ="T-AB")

combined_pt1014 <- combineTCR(contig_list_pt1014_sub, 
                       samples = c("pt1014", "pt1014", "pt1014"), 
                       ID = c("week1", "week2", "week4"), cells ="T-AB")

combined_pt1019 <- combineTCR(contig_list_pt1019_sub, 
                       samples = c("pt1019", "pt1019", "pt1019"), 
                       ID = c("week1", "week4", "week14"), cells ="T-AB")

combined_total <- c(combined_pt1012, combined_pt1013, combined_pt1014, combined_pt1019)

clonal_diversity_combined <- clonalDiversity(combined_total, cloneCall = "gene", n.boots = 100, x.axis = "ID", group.by = "sample") + 
  scale_x_discrete(limits = c("week1", "week2", "week4", "week14"), labels = c("Week 1", "Week 2", "Week 4", "Week 14")) + 
  scale_color_manual(values = c(col_pal[1], col_pal[3], col_pal[4], col_pal[7])) +
  theme(text = element_text(size = 20))
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/clonal_diversity_combined_0713.png", 
       clonal_diversity_combined, width = 24, height = 6, units = "in", limitsize = FALSE)
ggsave("/Users/kartiksinghal/Desktop/Rotation3-FehnigerLab/IL-7 scRNA/FinalFigures/clonal_diversity_combined_0713.pdf", 
       clonal_diversity_combined, width = 24, height = 6, units = "in", limitsize = FALSE)







