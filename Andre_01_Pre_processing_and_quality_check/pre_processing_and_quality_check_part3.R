library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(DoubletFinder)
library(cowplot)

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/01_Pre_processing_and_quality_check/doublet_cleaning")

#Load filtered files to check the dimension

merged_filtered_list_sct_dbl <- readRDS("merged_filtered_list_sct_dbl.rds")

# Retrieve the column names of the RNA assay feature in each Seurat object
selected_list <- lapply(merged_filtered_list_sct_dbl, function(x) {
  selected_colname <- colnames(x@meta.data)[grep("^DF.classifications", colnames(x@meta.data))]
  return(selected_colname)
})

# Generate a DimPlot for each Seurat object in the list, grouping by a specific column name
for (i in 1:length(merged_filtered_list_sct_dbl)) {
  # Get the column name for this Seurat object
  colname <- selected_list[[i]][1]
  # Generate the DimPlot
  p <- DimPlot(merged_filtered_list_sct_dbl[[i]], reduction = "umap.rna", label = FALSE, pt.size = 0.1, group.by = colname) + 
    ggtitle(paste0("Doublets", i))
    # Save using ggsave
  ggsave(filename = paste0("UMAP_doublets_", i, ".png"), plot = p)
}

# Specify the names for each object
object_names <- c("P27_Biorep1", "P27_Biorep2", "P27_Biorep3", "P27_Biorep4", "P27_Biorep5", "P27_Biorep6",
                  "LDR0_Biorep1", "LDR0_Biorep2", "LDR0_Biorep3",
                  "LDR30_Biorep1", "LDR30_Biorep2", "LDR30_Biorep3",
                  "LDR2_Biorep1", "LDR2_Biorep2", "LDR2_Biorep3",
                  "LDR4_Biorep1", "LDR4_Biorep2", "LDR4_Biorep3",
                  "LDR6_Biorep1", "LDR6_Biorep2", "LDR6_Biorep3")

# Split the list into separate objects with specific names
for (i in seq_along(merged_filtered_list_sct_dbl)) {
  assign(object_names[i], merged_filtered_list_sct_dbl[[i]])
}

DefaultAssay(P27_Biorep1) <- "RNA"
DefaultAssay(P27_Biorep2) <- "RNA"
DefaultAssay(P27_Biorep3) <- "RNA"
DefaultAssay(P27_Biorep4) <- "RNA"
DefaultAssay(P27_Biorep5) <- "RNA"
DefaultAssay(P27_Biorep6) <- "RNA"
DefaultAssay(LDR0_Biorep1) <- "RNA"
DefaultAssay(LDR0_Biorep2) <- "RNA"
DefaultAssay(LDR0_Biorep3) <- "RNA"
DefaultAssay(LDR30_Biorep1) <- "RNA"
DefaultAssay(LDR30_Biorep2) <- "RNA"
DefaultAssay(LDR30_Biorep3) <- "RNA"
DefaultAssay(LDR2_Biorep1) <- "RNA"
DefaultAssay(LDR2_Biorep2) <- "RNA"
DefaultAssay(LDR2_Biorep3) <- "RNA"
DefaultAssay(LDR4_Biorep1) <- "RNA"
DefaultAssay(LDR4_Biorep2) <- "RNA"
DefaultAssay(LDR4_Biorep3) <- "RNA"
DefaultAssay(LDR6_Biorep1) <- "RNA"
DefaultAssay(LDR6_Biorep2) <- "RNA"
DefaultAssay(LDR6_Biorep3) <- "RNA"

P27_Biorep1_no_doublets <- subset(x = P27_Biorep1, subset = DF.classifications_0.25_0.005_250.5 == "Singlet")
P27_Biorep2_no_doublets <- subset(x = P27_Biorep2, subset = DF.classifications_0.25_0.005_665.7 == "Singlet")
P27_Biorep3_no_doublets <- subset(x = P27_Biorep3, subset = DF.classifications_0.25_0.05_403.2 == "Singlet")
P27_Biorep4_no_doublets <- subset(x = P27_Biorep4, subset = DF.classifications_0.25_0.2_496.425 == "Singlet")
P27_Biorep5_no_doublets <- subset(x = P27_Biorep5, subset = DF.classifications_0.25_0.21_538.425 == "Singlet")
P27_Biorep6_no_doublets <- subset(x = P27_Biorep6, subset = DF.classifications_0.25_0.22_425.475 == "Singlet")
LDR0_Biorep1_no_doublets <- subset(x = LDR0_Biorep1, subset = DF.classifications_0.25_0.005_519.75 == "Singlet")
LDR0_Biorep2_no_doublets <- subset(x = LDR0_Biorep2, subset = DF.classifications_0.25_0.27_517.275 == "Singlet")
LDR0_Biorep3_no_doublets <- subset(x = LDR0_Biorep3, subset = DF.classifications_0.25_0.005_411.9 == "Singlet")
LDR30_Biorep1_no_doublets <- subset(x = LDR30_Biorep1, subset = DF.classifications_0.25_0.17_420.075 == "Singlet")
LDR30_Biorep2_no_doublets <- subset(x = LDR30_Biorep2, subset = DF.classifications_0.25_0.17_410.625 == "Singlet")
LDR30_Biorep3_no_doublets <- subset(x = LDR30_Biorep3, subset = DF.classifications_0.25_0.1_341.175 == "Singlet")
LDR2_Biorep1_no_doublets <- subset(x = LDR2_Biorep1, subset = DF.classifications_0.25_0.005_609.15 == "Singlet")
LDR2_Biorep2_no_doublets <- subset(x = LDR2_Biorep2, subset = DF.classifications_0.25_0.005_618.75 == "Singlet")
LDR2_Biorep3_no_doublets <- subset(x = LDR2_Biorep3, subset = DF.classifications_0.25_0.04_668.4 == "Singlet")
LDR4_Biorep1_no_doublets <- subset(x = LDR4_Biorep1, subset = DF.classifications_0.25_0.25_340.8 == "Singlet")
LDR4_Biorep2_no_doublets <- subset(x = LDR4_Biorep2, subset = DF.classifications_0.25_0.04_494.1 == "Singlet")
LDR4_Biorep3_no_doublets <- subset(x = LDR4_Biorep3, subset = DF.classifications_0.25_0.22_618.075 == "Singlet")
LDR6_Biorep1_no_doublets <- subset(x = LDR6_Biorep1, subset = DF.classifications_0.25_0.17_542.025 == "Singlet")
LDR6_Biorep2_no_doublets <- subset(x = LDR6_Biorep2, subset = DF.classifications_0.25_0.01_523.425 == "Singlet")
LDR6_Biorep3_no_doublets <- subset(x = LDR6_Biorep3, subset = DF.classifications_0.25_0.005_612.525 == "Singlet")

merged_filtered_list_sct_nodbl <- list(P27_Biorep1_no_doublets, P27_Biorep2_no_doublets, P27_Biorep3_no_doublets, P27_Biorep4_no_doublets, P27_Biorep5_no_doublets, P27_Biorep6_no_doublets,
                                         LDR0_Biorep1_no_doublets, LDR0_Biorep2_no_doublets, LDR0_Biorep3_no_doublets,
                                         LDR30_Biorep1_no_doublets, LDR30_Biorep2_no_doublets, LDR30_Biorep3_no_doublets,
                                         LDR2_Biorep1_no_doublets, LDR2_Biorep2_no_doublets, LDR2_Biorep3_no_doublets,
                                         LDR4_Biorep1_no_doublets, LDR4_Biorep2_no_doublets, LDR4_Biorep3_no_doublets,
                                         LDR6_Biorep1_no_doublets, LDR6_Biorep2_no_doublets, LDR6_Biorep3_no_doublets)

names(merged_filtered_list_sct_nodbl) <- object_names

# Check the dimension of each object in the list
dimensions_dbl <- lapply(merged_filtered_list_sct_dbl, function(x) dim(x))
dimensions_no_dbl <- lapply(merged_filtered_list_sct_nodbl, function(x) dim(x))

# Apply the DietSeurat function to each object in the list
dieted_seurat_list <- lapply(X = merged_filtered_list_sct_nodbl, FUN = function(x) DietSeurat(x, counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA"))

# Remove the specified columns from each object in the list
column_numbers <- c(1,5:6,14:18)

dieted_seurat_list = lapply(dieted_seurat_list, function(x) {
  x@meta.data = x@meta.data[,-c(column_numbers)]
  x
})

rename_cols <- function(x, col_pos, new_name) {
  names(x@meta.data)[col_pos] <- new_name
  return(x)
}
dieted_seurat_list <- lapply(dieted_seurat_list, rename_cols, col_pos = 11, new_name = "DoubletFinder")

saveRDS(dieted_seurat_list, "dieted_seurat_list.rds")

save.image("pre_processing_and_quality_check_part3.RData")