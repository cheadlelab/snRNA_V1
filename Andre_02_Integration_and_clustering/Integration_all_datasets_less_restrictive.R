library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(cowplot)

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/02_Integration_and_clustering/less_restrictive/")

merged_filtered <- readRDS("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/01_Pre_processing_and_quality_check/Ambient_RNA_cleaning/merged_filtered.rds")

table(merged_filtered$condition)
table(merged_filtered$condition_simpl)
table(merged_filtered$condition_simpl, merged_filtered$sample)

DefaultAssay(merged_filtered) <- "decontXcounts"

dim(merged_filtered)
merged_filtered <- subset(x = merged_filtered, 
              subset= (nCount_decontXcounts >= 150) & nFeature_decontXcounts >= 90)
dim(merged_filtered)

# Extract counts
decontXcounts <- GetAssayData(object = merged_filtered, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- decontXcounts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(decontXcounts) >= 10

# Only keeping those genes expressed in more than 10 cells
decontXfiltered_counts <- decontXcounts[keep_genes, ]

# Reassign to filtered Seurat object
merged_filtered[["decontXcountsclean"]] <- CreateAssayObject(counts = decontXfiltered_counts)
dim(merged_filtered)

DefaultAssay(merged_filtered) <- "decontXcountsclean"

merged_filtered$sample <- factor(merged_filtered$sample, levels = c("P27_Biorep1", "P27_Biorep2", "P27_Biorep3", "P27_Biorep4", "P27_Biorep5", "P27_Biorep6",
                                                              "LDR0_Biorep1", "LDR0_Biorep2", "LDR0_Biorep3",
                                                              "LDR30_Biorep1", "LDR30_Biorep2", "LDR30_Biorep3",
                                                              "LDR2_Biorep1", "LDR2_Biorep2", "LDR2_Biorep3",
                                                              "LDR4_Biorep1", "LDR4_Biorep2", "LDR4_Biorep3",
                                                              "LDR6_Biorep1", "LDR6_Biorep2", "LDR6_Biorep3"))

merged_filtered$batch <- factor(merged_filtered$batch, levels = c("batch1", "batch2", "batch3", "batch4", "batch5", "batch6"))

merged_filtered$condition <- factor(merged_filtered$condition, levels = c("ctrl", "dark reared", "dark reared + light stim30m", "dark reared + light stim2h", "dark reared + light stim4h", "dark reared + light stim6h"))

merged_filtered$condition_simpl <- factor(merged_filtered$condition_simpl, levels = c("NR", "LDR", "LDR30m", "LDR2h", "LDR4h", "LDR6h"))

merged_filtered_cleaned <- merged_filtered

saveRDS(merged_filtered_cleaned, "merged_filtered_cleaned.rds")

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/02_Integration_and_clustering/less_restrictive/no_batch_correction")

#SCT normalization
Allexp.sct_nobatch_corrected <- SCTransform(merged_filtered, verbose = FALSE, vars.to.regress = "percent.mt", variable.features.n = 3000, assay = "decontXcountsclean", new.assay.name = "sctdecontXcountsclean")

#Perform dimensionality reduction for clustering and check for the batch effects
Allexp.sct_nobatch_corrected <- RunPCA(Allexp.sct_nobatch_corrected, features = VariableFeatures(object = Allexp.sct_nobatch_corrected), reduction.name = "pca")
elbow_plot <- ElbowPlot(Allexp.sct_nobatch_corrected, ndims = 50)
ggsave(elbow_plot, filename = "elbow_plot_nobatch_corrected.png", height = 7, width = 7)

Allexp.sct_nobatch_corrected.40 <- FindNeighbors(Allexp.sct_nobatch_corrected, dims = 1:40)
Allexp.sct_nobatch_corrected.40nres <- FindClusters(Allexp.sct_nobatch_corrected.40)
Allexp.sct_nobatch_corrected.40.01 <- FindClusters(Allexp.sct_nobatch_corrected.40, resolution = 0.1)
Allexp.sct_nobatch_corrected.40.02 <- FindClusters(Allexp.sct_nobatch_corrected.40, resolution = 0.2)
Allexp.sct_nobatch_corrected.40.03 <- FindClusters(Allexp.sct_nobatch_corrected.40, resolution = 0.3)
Allexp.sct_nobatch_corrected.40.04 <- FindClusters(Allexp.sct_nobatch_corrected.40, resolution = 0.4)
Allexp.sct_nobatch_corrected.40.05 <- FindClusters(Allexp.sct_nobatch_corrected.40, resolution = 0.5)

Allexp.sct_nobatch_corrected.40nres <- RunUMAP(Allexp.sct_nobatch_corrected.40nres, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_nobatch_corrected.40.01 <- RunUMAP(Allexp.sct_nobatch_corrected.40.01, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_nobatch_corrected.40.02 <- RunUMAP(Allexp.sct_nobatch_corrected.40.02, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_nobatch_corrected.40.03 <- RunUMAP(Allexp.sct_nobatch_corrected.40.03, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_nobatch_corrected.40.04 <- RunUMAP(Allexp.sct_nobatch_corrected.40.04, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_nobatch_corrected.40.05 <- RunUMAP(Allexp.sct_nobatch_corrected.40.05, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')

UMAPAllexp.sct_nobatch_corrected.40nres1 <- DimPlot(Allexp.sct_nobatch_corrected.40nres, reduction = "UMAP", label = FALSE, split.by = "batch")
ggsave(UMAPAllexp.sct_nobatch_corrected.40nres1, filename = "UMAPAllexp.sct_nobatch_corrected.40nres1.png", height = 7, width = 23)

UMAPAllexp.sct_nobatch_corrected.40nres2 <- DimPlot(Allexp.sct_nobatch_corrected.40nres, reduction = "UMAP", label = FALSE, group.by = "batch")
ggsave(UMAPAllexp.sct_nobatch_corrected.40nres2, filename = "UMAPAllexp.sct_nobatch_corrected.40nres2.png", height = 7, width = 7)

UMAPAllexp.sct_nobatch_corrected.40nres3 <- DimPlot(Allexp.sct_nobatch_corrected.40nres, reduction = "UMAP", label = FALSE, split.by = "condition")
ggsave(UMAPAllexp.sct_nobatch_corrected.40nres3, filename = "UMAPAllexp.sct_nobatch_corrected.40nres3.png", height = 7, width = 23)

UMAPAllexp.sct_nobatch_corrected.40nres4 <- DimPlot(Allexp.sct_nobatch_corrected.40nres, reduction = "UMAP", label = FALSE, group.by = "condition")
ggsave(UMAPAllexp.sct_nobatch_corrected.40nres4, filename = "UMAPAllexp.sct_nobatch_corrected.40nres4.png", height = 7, width = 7)

UMAPAllexp.sct_nobatch_corrected.40.01 <- DimPlot(Allexp.sct_nobatch_corrected.40.01, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_nobatch_corrected.40.01, filename = "UMAPAllexp.sct_nobatch_corrected.40.01.png", height = 7, width = 7)
UMAPAllexp.sct_nobatch_corrected.40.02 <- DimPlot(Allexp.sct_nobatch_corrected.40.02, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_nobatch_corrected.40.02, filename = "UMAPAllexp.sct_nobatch_corrected.40.02.png", height = 7, width = 7)
UMAPAllexp.sct_nobatch_corrected.40.03 <- DimPlot(Allexp.sct_nobatch_corrected.40.03, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_nobatch_corrected.40.03, filename = "UMAPAllexp.sct_nobatch_corrected.40.03.png", height = 7, width = 7)
UMAPAllexp.sct_nobatch_corrected.40.04 <- DimPlot(Allexp.sct_nobatch_corrected.40.04, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_nobatch_corrected.40.04, filename = "UMAPAllexp.sct_nobatch_corrected.40.04.png", height = 7, width = 7)
UMAPAllexp.sct_nobatch_corrected.40.05 <- DimPlot(Allexp.sct_nobatch_corrected.40.05, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_nobatch_corrected.40.05, filename = "UMAPAllexp.sct_nobatch_corrected.40.05.png", height = 7, width = 7)

saveRDS(Allexp.sct_nobatch_corrected.40.01, "Allexp.sct_nobatch_corrected.40.01.rds")
saveRDS(Allexp.sct_nobatch_corrected.40.02, "Allexp.sct_nobatch_corrected.40.02.rds")
saveRDS(Allexp.sct_nobatch_corrected.40.03, "Allexp.sct_nobatch_corrected.40.03.rds")
saveRDS(Allexp.sct_nobatch_corrected.40.04, "Allexp.sct_nobatch_corrected.40.04.rds")
saveRDS(Allexp.sct_nobatch_corrected.40.05, "Allexp.sct_nobatch_corrected.40.05.rds")

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/02_Integration_and_clustering/less_restrictive/batch_correction")

#In case of batch effects, proceed with this code
#Splitting dataset, normalizing and scaling for RedDim
merged_filtered_list <- SplitObject(merged_filtered, split.by = "batch")
merged_filtered_list_sct <- lapply(X = merged_filtered_list, FUN = function(x) {
  x <- SCTransform(x, verbose = FALSE, return.only.var.genes = FALSE, vars.to.regress = "percent.mt", variable.features.n = 3000, assay = "decontXcountsclean", new.assay.name = "sctdecontXcountsclean")
  })
features <- SelectIntegrationFeatures(object.list = merged_filtered_list_sct, nfeatures = 3000)
merged_filtered_list_sct <- PrepSCTIntegration(object.list = merged_filtered_list_sct, anchor.features = features)
all.nuclei.anchors <- FindIntegrationAnchors(object.list = merged_filtered_list_sct, reference = 1, normalization.method = "SCT",
                                            anchor.features = features)
Allexp.sct_batch_corrected <- IntegrateData(anchorset = all.nuclei.anchors, normalization.method = "SCT")

saveRDS(Allexp.sct_batch_corrected, "Allexp.sct_batch_corrected.rds")

Allexp.sct_batch_corrected <- RunPCA(Allexp.sct_batch_corrected, features = VariableFeatures(object = Allexp.sct_batch_corrected), reduction.name = "pca")
elbow_plot <- ElbowPlot(Allexp.sct_batch_corrected, ndims = 50)
ggsave(elbow_plot, filename = "elbow_plot_batch_corrected.png", height = 7, width = 7)

Allexp.sct_batch_corrected.40 <- FindNeighbors(Allexp.sct_batch_corrected, dims = 1:40)
Allexp.sct_batch_corrected.40nres <- FindClusters(Allexp.sct_batch_corrected.40)
Allexp.sct_batch_corrected.40.01 <- FindClusters(Allexp.sct_batch_corrected.40, resolution = 0.1)
Allexp.sct_batch_corrected.40.02 <- FindClusters(Allexp.sct_batch_corrected.40, resolution = 0.2)
Allexp.sct_batch_corrected.40.03 <- FindClusters(Allexp.sct_batch_corrected.40, resolution = 0.3)
Allexp.sct_batch_corrected.40.04 <- FindClusters(Allexp.sct_batch_corrected.40, resolution = 0.4)
Allexp.sct_batch_corrected.40.05 <- FindClusters(Allexp.sct_batch_corrected.40, resolution = 0.5)

Allexp.sct_batch_corrected.40nres <- RunUMAP(Allexp.sct_batch_corrected.40nres, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_batch_corrected.40.01 <- RunUMAP(Allexp.sct_batch_corrected.40.01, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_batch_corrected.40.02 <- RunUMAP(Allexp.sct_batch_corrected.40.02, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_batch_corrected.40.03 <- RunUMAP(Allexp.sct_batch_corrected.40.03, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_batch_corrected.40.04 <- RunUMAP(Allexp.sct_batch_corrected.40.04, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')
Allexp.sct_batch_corrected.40.05 <- RunUMAP(Allexp.sct_batch_corrected.40.05, dims = 1:40, reduction.name = 'UMAP', reduction.key = 'rnaUMAP_')

UMAPAllexp.sct_batch_corrected.40nres1 <- DimPlot(Allexp.sct_batch_corrected.40nres, reduction = "UMAP", label = FALSE, split.by = "batch")
ggsave(UMAPAllexp.sct_batch_corrected.40nres1, filename = "UMAPAllexp.sct_batch_corrected.40nres1.png", height = 7, width = 23)

UMAPAllexp.sct_batch_corrected.40nres2 <- DimPlot(Allexp.sct_batch_corrected.40nres, reduction = "UMAP", label = FALSE, group.by = "batch")
ggsave(UMAPAllexp.sct_batch_corrected.40nres2, filename = "UMAPAllexp.sct_batch_corrected.40nres2.png", height = 7, width = 7)

UMAPAllexp.sct_batch_corrected.40nres3 <- DimPlot(Allexp.sct_batch_corrected.40nres, reduction = "UMAP", label = FALSE, split.by = "condition")
ggsave(UMAPAllexp.sct_batch_corrected.40nres3, filename = "UMAPAllexp.sct_batch_corrected.40nres3.png", height = 7, width = 23)

UMAPAllexp.sct_batch_corrected.40nres4 <- DimPlot(Allexp.sct_batch_corrected.40nres, reduction = "UMAP", label = FALSE, group.by = "condition")
ggsave(UMAPAllexp.sct_batch_corrected.40nres4, filename = "UMAPAllexp.sct_batch_corrected.40nres4.png", height = 7, width = 10)

UMAPAllexp.sct_batch_corrected.40.01 <- DimPlot(Allexp.sct_batch_corrected.40.01, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_batch_corrected.40.01, filename = "UMAPAllexp.sct_batch_corrected.40.01.png", height = 7, width = 7)
UMAPAllexp.sct_batch_corrected.40.02 <- DimPlot(Allexp.sct_batch_corrected.40.02, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_batch_corrected.40.02, filename = "UMAPAllexp.sct_batch_corrected.40.02.png", height = 7, width = 7)
UMAPAllexp.sct_batch_corrected.40.03 <- DimPlot(Allexp.sct_batch_corrected.40.03, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_batch_corrected.40.03, filename = "UMAPAllexp.sct_batch_corrected.40.03.png", height = 7, width = 7)
UMAPAllexp.sct_batch_corrected.40.04 <- DimPlot(Allexp.sct_batch_corrected.40.04, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_batch_corrected.40.04, filename = "UMAPAllexp.sct_batch_corrected.40.04.png", height = 7, width = 7)
UMAPAllexp.sct_batch_corrected.40.05 <- DimPlot(Allexp.sct_batch_corrected.40.05, reduction = "UMAP", label = TRUE)
ggsave(UMAPAllexp.sct_batch_corrected.40.05, filename = "UMAPAllexp.sct_batch_corrected.40.05.png", height = 7, width = 7)

saveRDS(Allexp.sct_batch_corrected.40.01, "Allexp.sct_batch_corrected.40.01.rds")
saveRDS(Allexp.sct_batch_corrected.40.02, "Allexp.sct_batch_corrected.40.02.rds")
saveRDS(Allexp.sct_batch_corrected.40.03, "Allexp.sct_batch_corrected.40.03.rds")
saveRDS(Allexp.sct_batch_corrected.40.04, "Allexp.sct_batch_corrected.40.04.rds")
saveRDS(Allexp.sct_batch_corrected.40.05, "Allexp.sct_batch_corrected.40.05.rds")
