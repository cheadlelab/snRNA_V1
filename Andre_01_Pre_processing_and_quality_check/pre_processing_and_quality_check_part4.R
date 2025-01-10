library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(dplyr)
library(stringr)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(celda)

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_final/01_Pre_processing_and_quality_check/Ambient_RNA_cleaning")

#Load cleaned Seurat objects
dieted_seurat_list  <- readRDS("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_final/01_Pre_processing_and_quality_check/doublet_cleaning/dieted_seurat_list.rds")

# Specify the names for each object
object_names_filtered <- c("P27_Biorep1_filtered", "P27_Biorep2_filtered", "P27_Biorep3_filtered", "P27_Biorep4_filtered", "P27_Biorep5_filtered", "P27_Biorep6_filtered",
                           "LDR0_Biorep4_filtered", "LDR0_Biorep5_filtered", "LDR0_Biorep6_filtered",
                           "LDR30_Biorep1_filtered", "LDR30_Biorep2_filtered", "LDR30_Biorep3_filtered",
                           "LDR2_Biorep1_filtered", "LDR2_Biorep2_filtered", "LDR2_Biorep3_filtered",
                           "LDR4_Biorep1_filtered", "LDR4_Biorep2_filtered", "LDR4_Biorep3_filtered",
                           "LDR6_Biorep1_filtered", "LDR6_Biorep2_filtered", "LDR6_Biorep3_filtered")

# Split the list into separate objects with specific names
for (i in seq_along(dieted_seurat_list)) {
  assign(object_names_filtered[i], dieted_seurat_list[[i]])
}

P27_1_filtered_counts <- GetAssayData(object = P27_Biorep1_filtered, slot = "counts")
P27_2_filtered_counts <- GetAssayData(object = P27_Biorep2_filtered, slot = "counts")
P27_3_filtered_counts <- GetAssayData(object = P27_Biorep3_filtered, slot = "counts")
P27_4_filtered_counts <- GetAssayData(object = P27_Biorep4_filtered, slot = "counts")
P27_5_filtered_counts <- GetAssayData(object = P27_Biorep5_filtered, slot = "counts")
P27_6_filtered_counts <- GetAssayData(object = P27_Biorep6_filtered, slot = "counts")
LDR0_4_filtered_counts <- GetAssayData(object = LDR0_Biorep4_filtered, slot = "counts")
LDR0_5_filtered_counts <- GetAssayData(object = LDR0_Biorep5_filtered, slot = "counts")
LDR0_6_filtered_counts <- GetAssayData(object = LDR0_Biorep6_filtered, slot = "counts")
LDR30_1_filtered_counts <- GetAssayData(object = LDR30_Biorep1_filtered, slot = "counts")
LDR30_2_filtered_counts <- GetAssayData(object = LDR30_Biorep2_filtered, slot = "counts")
LDR30_3_filtered_counts <- GetAssayData(object = LDR30_Biorep3_filtered, slot = "counts")
LDR2_1_filtered_counts <- GetAssayData(object = LDR2_Biorep1_filtered, slot = "counts")
LDR2_2_filtered_counts <- GetAssayData(object = LDR2_Biorep2_filtered, slot = "counts")
LDR2_3_filtered_counts <- GetAssayData(object = LDR2_Biorep3_filtered, slot = "counts")
LDR4_1_filtered_counts <- GetAssayData(object = LDR4_Biorep1_filtered, slot = "counts")
LDR4_2_filtered_counts <- GetAssayData(object = LDR4_Biorep2_filtered, slot = "counts")
LDR4_3_filtered_counts <- GetAssayData(object = LDR4_Biorep3_filtered, slot = "counts")
LDR6_1_filtered_counts <- GetAssayData(object = LDR6_Biorep1_filtered, slot = "counts")
LDR6_2_filtered_counts <- GetAssayData(object = LDR6_Biorep2_filtered, slot = "counts")
LDR6_3_filtered_counts <- GetAssayData(object = LDR6_Biorep3_filtered, slot = "counts")

P27_1_sce_filtered <- SingleCellExperiment(list(counts = P27_1_filtered_counts))
P27_2_sce_filtered <- SingleCellExperiment(list(counts = P27_2_filtered_counts))
P27_3_sce_filtered <- SingleCellExperiment(list(counts = P27_3_filtered_counts))
P27_4_sce_filtered <- SingleCellExperiment(list(counts = P27_4_filtered_counts))
P27_5_sce_filtered <- SingleCellExperiment(list(counts = P27_5_filtered_counts))
P27_6_sce_filtered <- SingleCellExperiment(list(counts = P27_6_filtered_counts))
LDR0_4_sce_filtered <- SingleCellExperiment(list(counts = LDR0_4_filtered_counts))
LDR0_5_sce_filtered <- SingleCellExperiment(list(counts = LDR0_5_filtered_counts))
LDR0_6_sce_filtered <- SingleCellExperiment(list(counts = LDR0_6_filtered_counts))
LDR30_1_sce_filtered <- SingleCellExperiment(list(counts = LDR30_1_filtered_counts))
LDR30_2_sce_filtered <- SingleCellExperiment(list(counts = LDR30_2_filtered_counts))
LDR30_3_sce_filtered <- SingleCellExperiment(list(counts = LDR30_3_filtered_counts))
LDR2_1_sce_filtered <- SingleCellExperiment(list(counts = LDR2_1_filtered_counts))
LDR2_2_sce_filtered <- SingleCellExperiment(list(counts = LDR2_2_filtered_counts))
LDR2_3_sce_filtered <- SingleCellExperiment(list(counts = LDR2_3_filtered_counts))
LDR4_1_sce_filtered <- SingleCellExperiment(list(counts = LDR4_1_filtered_counts))
LDR4_2_sce_filtered <- SingleCellExperiment(list(counts = LDR4_2_filtered_counts))
LDR4_3_sce_filtered <- SingleCellExperiment(list(counts = LDR4_3_filtered_counts))
LDR6_1_sce_filtered <- SingleCellExperiment(list(counts = LDR6_1_filtered_counts))
LDR6_2_sce_filtered <- SingleCellExperiment(list(counts = LDR6_2_filtered_counts))
LDR6_3_sce_filtered <- SingleCellExperiment(list(counts = LDR6_3_filtered_counts))

are_equal <- sapply(dieted_seurat_list, nrow) == nrow(dieted_seurat_list[[1]])
are_equal

keep_genes <- rownames(dieted_seurat_list$P27_Biorep1)
length(keep_genes)

# Use for FeatureMatrix and FindMarkers in Seurat -----------------------------------------------------------------------------------------------------------------------
library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 16000 * 1024^2)

samples = c("P27_1_raw_feature_bc_matrix", "P27_2_raw_feature_bc_matrix", "P27_3_raw_feature_bc_matrix", "P27_4_raw_feature_bc_matrix", "P27_5_raw_feature_bc_matrix","P27_6_raw_feature_bc_matrix",
            "LDR0_4_raw_feature_bc_matrix", "LDR0_5_raw_feature_bc_matrix", "LDR0_6_raw_feature_bc_matrix", "LDR30_1_raw_feature_bc_matrix", "LDR30_2_raw_feature_bc_matrix", "LDR30_3_raw_feature_bc_matrix",
            "LDR2_1_raw_feature_bc_matrix", "LDR2_2_raw_feature_bc_matrix", "LDR2_3_raw_feature_bc_matrix", "LDR4_1_raw_feature_bc_matrix", "LDR4_2_raw_feature_bc_matrix", "LDR4_3_raw_feature_bc_matrix",
            "LDR6_1_raw_feature_bc_matrix", "LDR6_2_raw_feature_bc_matrix", "LDR6_3_raw_feature_bc_matrix")

# Create a Seurat object for each sample
for (i in samples){
  seurat_data <- Read10X(data.dir = paste0("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/Data/", i))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = i)
  assign(i, seurat_obj)
}

# Create a merged Seurat object
All_combined <- merge(x = P27_1_raw_feature_bc_matrix, 
                      y = c(P27_2_raw_feature_bc_matrix,
                            P27_3_raw_feature_bc_matrix,
                            P27_4_raw_feature_bc_matrix,
                            P27_5_raw_feature_bc_matrix,
                            P27_6_raw_feature_bc_matrix,
                            LDR0_4_raw_feature_bc_matrix,
                            LDR0_5_raw_feature_bc_matrix,
                            LDR0_6_raw_feature_bc_matrix,
                            LDR30_1_raw_feature_bc_matrix,
                            LDR30_2_raw_feature_bc_matrix,
                            LDR30_3_raw_feature_bc_matrix,
                            LDR2_1_raw_feature_bc_matrix,
                            LDR2_2_raw_feature_bc_matrix,
                            LDR2_3_raw_feature_bc_matrix,
                            LDR4_1_raw_feature_bc_matrix,
                            LDR4_2_raw_feature_bc_matrix,
                            LDR4_3_raw_feature_bc_matrix,
                            LDR6_1_raw_feature_bc_matrix,
                            LDR6_2_raw_feature_bc_matrix,
                            LDR6_3_raw_feature_bc_matrix), 
                      add.cell.ids = c("P27_1","P27_2", "P27_3", "P27_4", "P27_5", "P27_6", "LDR0_4", "LDR0_5", "LDR0_6", "LDR30_1", "LDR30_2", "LDR30_3", "LDR2_1", "LDR2_2", "LDR2_3" , "LDR4_1", "LDR4_2", "LDR4_3", "LDR6_1", "LDR6_2", "LDR6_3"), 
                      project = "All_combined")

cells <- colnames(All_combined)
dim(All_combined)
length(cells)
cells <- gsub("_", "", cells)
length(cells)
All_combined <- RenameCells(object = All_combined, new.names = cells)
dim(All_combined)
head(All_combined@meta.data)

DefaultAssay(All_combined) <- "RNA"

# Create metadata dataframe
metadata <- All_combined@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
head(metadata)

# Create Sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^P271"))] <- "P27_Biorep1"
metadata$sample[which(str_detect(metadata$cells, "^P272"))] <- "P27_Biorep2"
metadata$sample[which(str_detect(metadata$cells, "^P273"))] <- "P27_Biorep3"
metadata$sample[which(str_detect(metadata$cells, "^P274"))] <- "P27_Biorep4"
metadata$sample[which(str_detect(metadata$cells, "^P275"))] <- "P27_Biorep5"
metadata$sample[which(str_detect(metadata$cells, "^P276"))] <- "P27_Biorep6"
metadata$sample[which(str_detect(metadata$cells, "^LDR04"))] <- "LDR0_Biorep4"
metadata$sample[which(str_detect(metadata$cells, "^LDR05"))] <- "LDR0_Biorep5"
metadata$sample[which(str_detect(metadata$cells, "^LDR06"))] <- "LDR0_Biorep6"
metadata$sample[which(str_detect(metadata$cells, "^LDR301"))] <- "LDR30_Biorep1"
metadata$sample[which(str_detect(metadata$cells, "^LDR302"))] <- "LDR30_Biorep2"
metadata$sample[which(str_detect(metadata$cells, "^LDR303"))] <- "LDR30_Biorep3"
metadata$sample[which(str_detect(metadata$cells, "^LDR21"))] <- "LDR2_Biorep1"
metadata$sample[which(str_detect(metadata$cells, "^LDR22"))] <- "LDR2_Biorep2"
metadata$sample[which(str_detect(metadata$cells, "^LDR23"))] <- "LDR2_Biorep3"
metadata$sample[which(str_detect(metadata$cells, "^LDR41"))] <- "LDR4_Biorep1"
metadata$sample[which(str_detect(metadata$cells, "^LDR42"))] <- "LDR4_Biorep2"
metadata$sample[which(str_detect(metadata$cells, "^LDR43"))] <- "LDR4_Biorep3"
metadata$sample[which(str_detect(metadata$cells, "^LDR61"))] <- "LDR6_Biorep1"
metadata$sample[which(str_detect(metadata$cells, "^LDR62"))] <- "LDR6_Biorep2"
metadata$sample[which(str_detect(metadata$cells, "^LDR63"))] <- "LDR6_Biorep3"

All_combined@meta.data <- metadata

#Split object to run doubletfinder
merged_raw_list <- SplitObject(All_combined, split.by = "sample")
merged_raw_list

# Specify the names for each object
object_names <- c("P27_Biorep1_raw", "P27_Biorep2_raw", "P27_Biorep3_raw", "P27_Biorep4_raw", "P27_Biorep5_raw", "P27_Biorep6_raw",
                  "LDR0_Biorep4_raw", "LDR0_Biorep5_raw", "LDR0_Biorep6_raw",
                  "LDR30_Biorep1_raw", "LDR30_Biorep2_raw", "LDR30_Biorep3_raw",
                  "LDR2_Biorep1_raw", "LDR2_Biorep2_raw", "LDR2_Biorep3_raw",
                  "LDR4_Biorep1_raw", "LDR4_Biorep2_raw", "LDR4_Biorep3_raw",
                  "LDR6_Biorep1_raw", "LDR6_Biorep2_raw", "LDR6_Biorep3_raw")

# Split the list into separate objects with specific names
for (i in seq_along(merged_raw_list)) {
  assign(object_names[i], merged_raw_list[[i]])
}

P27_1_raw_counts <- GetAssayData(object = P27_Biorep1_raw, slot = "counts")
P27_2_raw_counts <- GetAssayData(object = P27_Biorep2_raw, slot = "counts")
P27_3_raw_counts <- GetAssayData(object = P27_Biorep3_raw, slot = "counts")
P27_4_raw_counts <- GetAssayData(object = P27_Biorep4_raw, slot = "counts")
P27_5_raw_counts <- GetAssayData(object = P27_Biorep5_raw, slot = "counts")
P27_6_raw_counts <- GetAssayData(object = P27_Biorep6_raw, slot = "counts")
LDR0_4_raw_counts <- GetAssayData(object = LDR0_Biorep4_raw, slot = "counts")
LDR0_5_raw_counts <- GetAssayData(object = LDR0_Biorep5_raw, slot = "counts")
LDR0_6_raw_counts <- GetAssayData(object = LDR0_Biorep6_raw, slot = "counts")
LDR30_1_raw_counts <- GetAssayData(object = LDR30_Biorep1_raw, slot = "counts")
LDR30_2_raw_counts <- GetAssayData(object = LDR30_Biorep2_raw, slot = "counts")
LDR30_3_raw_counts <- GetAssayData(object = LDR30_Biorep3_raw, slot = "counts")
LDR2_1_raw_counts <- GetAssayData(object = LDR2_Biorep1_raw, slot = "counts")
LDR2_2_raw_counts <- GetAssayData(object = LDR2_Biorep2_raw, slot = "counts")
LDR2_3_raw_counts <- GetAssayData(object = LDR2_Biorep3_raw, slot = "counts")
LDR4_1_raw_counts <- GetAssayData(object = LDR4_Biorep1_raw, slot = "counts")
LDR4_2_raw_counts <- GetAssayData(object = LDR4_Biorep2_raw, slot = "counts")
LDR4_3_raw_counts <- GetAssayData(object = LDR4_Biorep3_raw, slot = "counts")
LDR6_1_raw_counts <- GetAssayData(object = LDR6_Biorep1_raw, slot = "counts")
LDR6_2_raw_counts <- GetAssayData(object = LDR6_Biorep2_raw, slot = "counts")
LDR6_3_raw_counts <- GetAssayData(object = LDR6_Biorep3_raw, slot = "counts")

P27_1_raw_counts <- P27_1_raw_counts[keep_genes, ]
P27_2_raw_counts <- P27_2_raw_counts[keep_genes, ]
P27_3_raw_counts <- P27_3_raw_counts[keep_genes, ]
P27_4_raw_counts <- P27_4_raw_counts[keep_genes, ]
P27_5_raw_counts <- P27_5_raw_counts[keep_genes, ]
P27_6_raw_counts <- P27_6_raw_counts[keep_genes, ]
LDR0_4_raw_counts <- LDR0_4_raw_counts[keep_genes, ]
LDR0_5_raw_counts <- LDR0_5_raw_counts[keep_genes, ]
LDR0_6_raw_counts <- LDR0_6_raw_counts[keep_genes, ]
LDR30_1_raw_counts <- LDR30_1_raw_counts[keep_genes, ]
LDR30_2_raw_counts <- LDR30_2_raw_counts[keep_genes, ]
LDR30_3_raw_counts <- LDR30_3_raw_counts[keep_genes, ]
LDR2_1_raw_counts <- LDR2_1_raw_counts[keep_genes, ]
LDR2_2_raw_counts <- LDR2_2_raw_counts[keep_genes, ]
LDR2_3_raw_counts <- LDR2_3_raw_counts[keep_genes, ]
LDR4_1_raw_counts <- LDR4_1_raw_counts[keep_genes, ]
LDR4_2_raw_counts <- LDR4_2_raw_counts[keep_genes, ]
LDR4_3_raw_counts <- LDR4_3_raw_counts[keep_genes, ]
LDR6_1_raw_counts <- LDR6_1_raw_counts[keep_genes, ]
LDR6_2_raw_counts <- LDR6_2_raw_counts[keep_genes, ]
LDR6_3_raw_counts <- LDR6_3_raw_counts[keep_genes, ]

P27_1_sce_raw <- SingleCellExperiment(list(counts = P27_1_raw_counts))
P27_2_sce_raw <- SingleCellExperiment(list(counts = P27_2_raw_counts))
P27_3_sce_raw <- SingleCellExperiment(list(counts = P27_3_raw_counts))
P27_4_sce_raw <- SingleCellExperiment(list(counts = P27_4_raw_counts))
P27_5_sce_raw <- SingleCellExperiment(list(counts = P27_5_raw_counts))
P27_6_sce_raw <- SingleCellExperiment(list(counts = P27_6_raw_counts))
LDR0_4_sce_raw <- SingleCellExperiment(list(counts = LDR0_4_raw_counts))
LDR0_5_sce_raw <- SingleCellExperiment(list(counts = LDR0_5_raw_counts))
LDR0_6_sce_raw <- SingleCellExperiment(list(counts = LDR0_6_raw_counts))
LDR30_1_sce_raw <- SingleCellExperiment(list(counts = LDR30_1_raw_counts))
LDR30_2_sce_raw <- SingleCellExperiment(list(counts = LDR30_2_raw_counts))
LDR30_3_sce_raw <- SingleCellExperiment(list(counts = LDR30_3_raw_counts))
LDR2_1_sce_raw <- SingleCellExperiment(list(counts = LDR2_1_raw_counts))
LDR2_2_sce_raw <- SingleCellExperiment(list(counts = LDR2_2_raw_counts))
LDR2_3_sce_raw <- SingleCellExperiment(list(counts = LDR2_3_raw_counts))
LDR4_1_sce_raw <- SingleCellExperiment(list(counts = LDR4_1_raw_counts))
LDR4_2_sce_raw <- SingleCellExperiment(list(counts = LDR4_2_raw_counts))
LDR4_3_sce_raw <- SingleCellExperiment(list(counts = LDR4_3_raw_counts))
LDR6_1_sce_raw <- SingleCellExperiment(list(counts = LDR6_1_raw_counts))
LDR6_2_sce_raw <- SingleCellExperiment(list(counts = LDR6_2_raw_counts))
LDR6_3_sce_raw <- SingleCellExperiment(list(counts = LDR6_3_raw_counts))

save.image("pre_processing_and_quality_check_part4.RData")

# Running decontX
P27_1_sce_x <- decontX(P27_1_sce_filtered, background = P27_1_sce_raw)
P27_2_sce_x <- decontX(P27_2_sce_filtered, background = P27_2_sce_raw)
P27_3_sce_x <- decontX(P27_3_sce_filtered, background = P27_3_sce_raw)
P27_4_sce_x <- decontX(P27_4_sce_filtered, background = P27_4_sce_raw)
P27_5_sce_x <- decontX(P27_5_sce_filtered, background = P27_5_sce_raw)
P27_6_sce_x <- decontX(P27_6_sce_filtered, background = P27_6_sce_raw)
LDR0_3_sce_x <- decontX(LDR0_3_sce_filtered, background = LDR0_3_sce_raw)
LDR0_4_sce_x <- decontX(LDR0_4_sce_filtered, background = LDR0_4_sce_raw)
LDR0_5_sce_x <- decontX(LDR0_5_sce_filtered, background = LDR0_5_sce_raw)
LDR30_1_sce_x <- decontX(LDR30_1_sce_filtered, background = LDR30_1_sce_raw)
LDR30_2_sce_x <- decontX(LDR30_2_sce_filtered, background = LDR30_2_sce_raw)
LDR30_3_sce_x <- decontX(LDR30_3_sce_filtered, background = LDR30_3_sce_raw)
LDR2_1_sce_x <- decontX(LDR2_1_sce_filtered, background = LDR2_1_sce_raw)
LDR2_2_sce_x <- decontX(LDR2_2_sce_filtered, background = LDR2_2_sce_raw)
LDR2_3_sce_x <- decontX(LDR2_3_sce_filtered, background = LDR2_3_sce_raw)
LDR4_1_sce_x <- decontX(LDR4_1_sce_filtered, background = LDR4_1_sce_raw)
LDR4_2_sce_x <- decontX(LDR4_2_sce_filtered, background = LDR4_2_sce_raw)
LDR4_3_sce_x <- decontX(LDR4_3_sce_filtered, background = LDR4_3_sce_raw)
LDR6_1_sce_x <- decontX(LDR6_1_sce_filtered, background = LDR6_1_sce_raw)
LDR6_2_sce_x <- decontX(LDR6_2_sce_filtered, background = LDR6_2_sce_raw)
LDR6_3_sce_x <- decontX(LDR6_3_sce_filtered, background = LDR6_3_sce_raw)

#Round counts to integer and create a Seurat assay from a SCE with decontX results
P27_Biorep1_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(P27_1_sce_x)))
P27_Biorep2_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(P27_2_sce_x)))
P27_Biorep3_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(P27_3_sce_x)))
P27_Biorep4_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(P27_4_sce_x)))
P27_Biorep5_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(P27_5_sce_x)))
P27_Biorep6_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(P27_6_sce_x)))
LDR0_Biorep4_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR0_4_sce_x)))
LDR0_Biorep5_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR0_5_sce_x)))
LDR0_Biorep6_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR0_6_sce_x)))
LDR30_Biorep1_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR30_1_sce_x)))
LDR30_Biorep2_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR30_2_sce_x)))
LDR30_Biorep3_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR30_3_sce_x)))
LDR2_Biorep1_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR2_1_sce_x)))
LDR2_Biorep2_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR2_2_sce_x)))
LDR2_Biorep3_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR2_3_sce_x)))
LDR4_Biorep1_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR4_1_sce_x)))
LDR4_Biorep2_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR4_2_sce_x)))
LDR4_Biorep3_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR4_3_sce_x)))
LDR6_Biorep1_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR6_1_sce_x)))
LDR6_Biorep2_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR6_2_sce_x)))
LDR6_Biorep3_filtered[["decontXcounts"]] <- CreateAssayObject(counts = round(decontXcounts(LDR6_3_sce_x)))







## Barplot of markers detected in cell clusters
#sce <- decontX(sce)
#seuratObj[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce))

#counts.raw <- Read10X("sample/outs/raw_feature_bc_matrix/")
#sce.raw <- SingleCellExperiment(list(counts = counts.raw))
#sce <- decontX(sce, background = sce.raw)

# Create a SingleCellExperiment object and run decontX
#sce <- SingleCellExperiment(list(counts = counts))
#sce <- decontX(sce)


#Note that the decontaminated matrix of decontX consists of floating point numbers and must be rounded to integers before adding it to a Seurat object. If you already have a Seurat object containing the counts matrix and would like to run decontX, you can retrieve the count matrix, create a SCE object, and run decontX, and then add it back to the Seurat object:
  
