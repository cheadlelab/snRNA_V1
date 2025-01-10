library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(dplyr)
library(stringr)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(celda)
library(scCustomize)

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/01_Pre_processing_and_quality_check/Ambient_RNA_cleaning")

#Load cleaned Seurat objects
dieted_seurat_list  <- readRDS("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/01_Pre_processing_and_quality_check/doublet_cleaning/dieted_seurat_list.rds")

merged_filtered <- Merge_Seurat_List(list_seurat = dieted_seurat_list, project = "merged_filtered")

table(merged_filtered$condition_simpl)
table(merged_filtered$condition)

merged_filtered_counts <- GetAssayData(object = merged_filtered, slot = "counts")

merged_filtered_sce <- SingleCellExperiment(list(counts = merged_filtered_counts))
nrow(merged_filtered_sce)

keep_genes <- rownames(merged_filtered_sce)
length(keep_genes)

# Load the raw datasets for ambient RNA cleaning -----------------------------------------------------------------------------------------------------------------------
library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 16000 * 1024^2)

samples = c("P27_1_raw_feature_bc_matrix", "P27_2_raw_feature_bc_matrix", "P27_3_raw_feature_bc_matrix", "P27_4_raw_feature_bc_matrix", "P27_5_raw_feature_bc_matrix","P27_6_raw_feature_bc_matrix",
            "LDR0_1_raw_feature_bc_matrix", "LDR0_2_raw_feature_bc_matrix", "LDR0_3_raw_feature_bc_matrix", "LDR30_1_raw_feature_bc_matrix", "LDR30_2_raw_feature_bc_matrix", "LDR30_3_raw_feature_bc_matrix",
            "LDR2_1_raw_feature_bc_matrix", "LDR2_2_raw_feature_bc_matrix", "LDR2_3_raw_feature_bc_matrix", "LDR4_1_raw_feature_bc_matrix", "LDR4_2_raw_feature_bc_matrix", "LDR4_3_raw_feature_bc_matrix",
            "LDR6_1_raw_feature_bc_matrix", "LDR6_2_raw_feature_bc_matrix", "LDR6_3_raw_feature_bc_matrix")

# Create a Seurat object for each sample
for (i in samples){
  seurat_data <- Read10X(data.dir = paste0("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/Data_Lin/", i))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = i)
  assign(i, seurat_obj)
}

# Create a merged Seurat object
merged_raw <- merge(x = P27_1_raw_feature_bc_matrix, 
                      y = c(P27_2_raw_feature_bc_matrix,
                            P27_3_raw_feature_bc_matrix,
                            P27_4_raw_feature_bc_matrix,
                            P27_5_raw_feature_bc_matrix,
                            P27_6_raw_feature_bc_matrix,
                            LDR0_1_raw_feature_bc_matrix,
                            LDR0_2_raw_feature_bc_matrix,
                            LDR0_3_raw_feature_bc_matrix,
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
                      add.cell.ids = c("P27_1","P27_2", "P27_3", "P27_4", "P27_5", "P27_6", "LDR0_1", "LDR0_2", "LDR0_3", "LDR30_1", "LDR30_2", "LDR30_3", "LDR2_1", "LDR2_2", "LDR2_3" , "LDR4_1", "LDR4_2", "LDR4_3", "LDR6_1", "LDR6_2", "LDR6_3"), 
                      project = "merged_raw")

cells <- colnames(merged_raw)
dim(merged_raw)
length(cells)
cells <- gsub("_", "", cells)
length(cells)
merged_raw <- RenameCells(object = merged_raw, new.names = cells)
dim(merged_raw)
head(merged_raw@meta.data)

DefaultAssay(merged_raw) <- "RNA"

# Create metadata dataframe
metadata <- merged_raw@meta.data

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
metadata$sample[which(str_detect(metadata$cells, "^LDR01"))] <- "LDR0_Biorep1"
metadata$sample[which(str_detect(metadata$cells, "^LDR02"))] <- "LDR0_Biorep2"
metadata$sample[which(str_detect(metadata$cells, "^LDR03"))] <- "LDR0_Biorep3"
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

merged_raw@meta.data <- metadata

merged_raw_counts <- GetAssayData(object = merged_raw, slot = "counts")
merged_raw_counts <- merged_raw_counts[keep_genes, ]

merged_raw_sce <- SingleCellExperiment(list(counts = merged_raw_counts))
nrow(merged_raw_sce)

merged_filtered_sce
merged_raw_sce

#Running decontX
merged_filtered_sce_x <- decontX(merged_filtered_sce, background = merged_raw_sce)

#Convert the decontXmatrix in integer and create a new assay in Seurat object
mat <- as.matrix(merged_filtered_sce_x@assays@data$decontXcounts)
storage.mode(mat) <- "integer"
merged_filtered[["decontXcounts"]] <- CreateAssayObject(counts = mat)

saveRDS(merged_filtered, "merged_filtered.rds")
saveRDS(merged_raw, "merged_raw.rds")
saveRDS(merged_filtered_sce, "merged_filtered_sce.rds")
saveRDS(merged_raw_sce, "merged_raw_sce.rds")
saveRDS(merged_filtered_sce_x, "merged_filtered_sce_x.rds")

merged_filtered_sce_x <- readRDS("merged_filtered_sce_x.rds")
merged_filtered <- readRDS("merged_filtered.rds")

# Plotting DecontX results
# Cluster labels on UMAP
umap <- reducedDim(merged_filtered_sce_x, "decontX_UMAP")
Amb_RNA_plot_cluster <- plotDimReduceCluster(x = merged_filtered_sce_x$decontX_clusters,
                                             dim1 = umap[, 1], dim2 = umap[, 2], labelClusters = TRUE)
ggsave(Amb_RNA_plot_cluster, filename = "Amb_RNA_plot_cluster.png", height = 7, width = 7)

#Percentage of contamination
perc_cont <- plotDecontXContamination(merged_filtered_sce_x)
ggsave(perc_cont, filename = "perc_cont.png", height = 7, width = 7)

#Plotting based in markers
merged_filtered_sce_x <- logNormCounts(merged_filtered_sce_x)
cluster_decontX <- plotDimReduceFeature(as.matrix(logcounts(merged_filtered_sce_x)),
                                        dim1 = umap[, 1],
                                        dim2 = umap[, 2],
                                        features = c("Slc17a7", "Gad1", "Tgfbr1", "Pdgfra", "Mobp", "Atp1a2"),
                                        exactMatch = TRUE)
ggsave(cluster_decontX, filename = "cluster_decontX.png", height = 20, width = 30)

## Barplot of markers detected in cell clusters
markers <- list(Exc_Neurons_Markers = c("Slc17a7", "Rora"),
                Inh_Neurons_Markers = c("Gad1", "Gad2"),
                Microglia_Markers = c("Tgfbr1", "Cx3cr1"),
                OPCs_Markers = c("Pdgfra", "Cspg4"),
                Astrocytes_Markers = c("Atp1a2", "Aqp4"),
                Oligodendrocytes_Markers = c("Mobp", "Mog"))

cellTypeMappings <- list('Excitatory Neurons' = c(4, 5, 6, 8, 10, 13), 'Inhibitory Neurons' = c(1, 11), Microglia = 7, 'Oligo/OPCs' = 2, Astrocytes = 3)

decontx_marker_per <- plotDecontXMarkerPercentage(merged_filtered_sce_x,
                            markers = markers,
                            groupClusters = cellTypeMappings,
                            assayName = "counts")
ggsave(decontx_marker_per, filename = "decontx_marker_per.png", height = 7, width = 7)

decontx_marker_per_clean <- plotDecontXMarkerPercentage(merged_filtered_sce_x,
                            markers = markers,
                            groupClusters = cellTypeMappings,
                            assayName = "decontXcounts")
ggsave(decontx_marker_per_clean, filename = "decontx_marker_per_clean.png", height = 7, width = 7)

decontx_marker_per_paired <- plotDecontXMarkerPercentage(merged_filtered_sce_x,
                            markers = markers,
                            groupClusters = cellTypeMappings,
                            assayName = c("counts", "decontXcounts"))
ggsave(decontx_marker_per_paired, filename = "decontx_marker_per_paired.png", height = 7, width = 14)

#Violin plot to check the gene distribution in the cells
Neuronal_Markers = c("Erbb4", "Nrxn3", "Lrrtm4", "Nrxn1", "Nrg3", "Csmd1", "Ptprd", "Kcnip4")
decontx_violin_plot_paired <- plotDecontXMarkerExpression(merged_filtered_sce_x,
                            markers = Neuronal_Markers,
                            groupClusters = cellTypeMappings,
                            ncol = 8)
ggsave(decontx_violin_plot_paired, filename = "decontx_violin_plot_paired.png", height = 14, width = 28)

merged_filtered_sce_x <- logNormCounts(merged_filtered_sce_x,
                     exprs_values = "decontXcounts",
                     name = "decontXlogcounts")

decontx_violin_plot_paired_lognorm <- plotDecontXMarkerExpression(merged_filtered_sce_x,
                            markers = Neuronal_Markers,
                            groupClusters = cellTypeMappings,
                            ncol = 8,
                            assayName = c("logcounts", "decontXlogcounts"))
ggsave(decontx_violin_plot_paired_lognorm, filename = "decontx_violin_plot_paired_lognorm.png", height = 14, width = 28)
