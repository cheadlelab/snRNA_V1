library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(cowplot)

Neuron_subset <- readRDS("/grid/cheadle/home/machado/scRNAseq_V1/Andre-snRNAseq/R-analysis/03_Integration_and_clustering/Analysis_P27_LDR0_LDR2_LDR6/Sub-Cluster_Neurons/Neuron_subset.rds")

DefaultAssay(Neuron_subset) <- "RNA"

setwd("/grid/cheadle/home/machado/scRNAseq_V1/Andre-snRNAseq/R-analysis/03_Integration_and_clustering/Analysis_P27_LDR0_LDR2_LDR6/Sub-Cluster_Neurons/Not_batch_corrected")

#SCT normalization
Neuron_subset_nobatch_corrected <- SCTransform(Neuron_subset, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = 3000)

#Perform dimensionality reduction for clustering and check for the batch effects
Neuron_subset_nobatch_corrected <- RunPCA(Neuron_subset_nobatch_corrected, features = VariableFeatures(object = Neuron_subset_nobatch_corrected))
elbow_plot <- ElbowPlot(Neuron_subset_nobatch_corrected, ndims = 50)
ggsave(elbow_plot, filename = "elbow_plot_nobatch_corrected.png", height = 7, width = 7, type = "cairo")

Neuron_subset_nobatch_corrected.40 <- FindNeighbors(Neuron_subset_nobatch_corrected, dims = 1:40)
Neuron_subset_nobatch_corrected.40nres <- FindClusters(Neuron_subset_nobatch_corrected.40)
Neuron_subset_nobatch_corrected.40.01 <- FindClusters(Neuron_subset_nobatch_corrected.40, resolution = 0.1)
Neuron_subset_nobatch_corrected.40.02 <- FindClusters(Neuron_subset_nobatch_corrected.40, resolution = 0.2)
Neuron_subset_nobatch_corrected.40.03 <- FindClusters(Neuron_subset_nobatch_corrected.40, resolution = 0.3)
Neuron_subset_nobatch_corrected.40.04 <- FindClusters(Neuron_subset_nobatch_corrected.40, resolution = 0.4)
Neuron_subset_nobatch_corrected.40.05 <- FindClusters(Neuron_subset_nobatch_corrected.40, resolution = 0.5)

Neuron_subset_nobatch_corrected.40nres <- RunUMAP(Neuron_subset_nobatch_corrected.40nres, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_nobatch_corrected.40.01 <- RunUMAP(Neuron_subset_nobatch_corrected.40.01, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_nobatch_corrected.40.02 <- RunUMAP(Neuron_subset_nobatch_corrected.40.02, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_nobatch_corrected.40.03 <- RunUMAP(Neuron_subset_nobatch_corrected.40.03, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_nobatch_corrected.40.04 <- RunUMAP(Neuron_subset_nobatch_corrected.40.04, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_nobatch_corrected.40.05 <- RunUMAP(Neuron_subset_nobatch_corrected.40.05, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

UMAPNeuron_subset_nobatch_corrected.40nres1 <- DimPlot(Neuron_subset_nobatch_corrected.40nres, reduction = "umap.rna", label = TRUE, split.by = "batch")
ggsave(UMAPNeuron_subset_nobatch_corrected.40nres1, filename = "UMAPNeuron_subset_nobatch_corrected.40nres1.png", height = 7, width = 14, type = "cairo")

UMAPNeuron_subset_nobatch_corrected.40nres2 <- DimPlot(Neuron_subset_nobatch_corrected.40nres, reduction = "umap.rna", label = TRUE, group.by = "batch")
ggsave(UMAPNeuron_subset_nobatch_corrected.40nres2, filename = "UMAPNeuron_subset_nobatch_corrected.40nres2.png", height = 7, width = 7, type = "cairo")

UMAPNeuron_subset_nobatch_corrected.40nres3 <- DimPlot(Neuron_subset_nobatch_corrected.40nres, reduction = "umap.rna", label = TRUE, split.by = "stim")
ggsave(UMAPNeuron_subset_nobatch_corrected.40nres3, filename = "UMAPNeuron_subset_nobatch_corrected.40nres3.png", height = 7, width = 14, type = "cairo")

UMAPNeuron_subset_nobatch_corrected.40nres4 <- DimPlot(Neuron_subset_nobatch_corrected.40nres, reduction = "umap.rna", label = TRUE, group.by = "stim")
ggsave(UMAPNeuron_subset_nobatch_corrected.40nres4, filename = "UMAPNeuron_subset_nobatch_corrected.40nres4.png", height = 7, width = 7, type = "cairo")

UMAPNeuron_subset_nobatch_corrected.40.01 <- DimPlot(Neuron_subset_nobatch_corrected.40.01, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_nobatch_corrected.40.01, filename = "UMAPNeuron_subset_nobatch_corrected.40.01.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_nobatch_corrected.40.02 <- DimPlot(Neuron_subset_nobatch_corrected.40.02, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_nobatch_corrected.40.02, filename = "UMAPNeuron_subset_nobatch_corrected.40.02.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_nobatch_corrected.40.03 <- DimPlot(Neuron_subset_nobatch_corrected.40.03, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_nobatch_corrected.40.03, filename = "UMAPNeuron_subset_nobatch_corrected.40.03.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_nobatch_corrected.40.04 <- DimPlot(Neuron_subset_nobatch_corrected.40.04, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_nobatch_corrected.40.04, filename = "UMAPNeuron_subset_nobatch_corrected.40.04.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_nobatch_corrected.40.05 <- DimPlot(Neuron_subset_nobatch_corrected.40.05, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_nobatch_corrected.40.05, filename = "UMAPNeuron_subset_nobatch_corrected.40.05.png", height = 7, width = 7, type = "cairo")

saveRDS(Neuron_subset_nobatch_corrected.40.01, "Neuron_subset_nobatch_corrected.40.01.rds")
saveRDS(Neuron_subset_nobatch_corrected.40.02, "Neuron_subset_nobatch_corrected.40.02.rds")
saveRDS(Neuron_subset_nobatch_corrected.40.03, "Neuron_subset_nobatch_corrected.40.03.rds")
saveRDS(Neuron_subset_nobatch_corrected.40.04, "Neuron_subset_nobatch_corrected.40.04.rds")
saveRDS(Neuron_subset_nobatch_corrected.40.05, "Neuron_subset_nobatch_corrected.40.05.rds")

setwd("/grid/cheadle/home/machado/scRNAseq_V1/Andre-snRNAseq/R-analysis/03_Integration_and_clustering/Analysis_P27_LDR0_LDR2_LDR6/Sub-Cluster_Neurons/batch_corrected")

#In case of batch effects, proceed with this code
Neuron_subset_list <- SplitObject(Neuron_subset, split.by = "batch")
Neuron_subset_list <- lapply(X = Neuron_subset_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = Neuron_subset_list, nfeatures = 3000)
Neuron_subset_list <- PrepSCTIntegration(object.list = Neuron_subset_list, anchor.features = features)
neuron.anchors <- FindIntegrationAnchors(object.list = Neuron_subset_list, reference = 1, normalization.method = "SCT",
                                            anchor.features = features)
Neuron_subset_batch_corrected <- IntegrateData(anchorset = neuron.anchors, normalization.method = "SCT")

saveRDS(Neuron_subset_batch_corrected, "Neuron_subset_batch_corrected.rds")

Neuron_subset_batch_corrected <- RunPCA(Neuron_subset_batch_corrected, features = VariableFeatures(object = Neuron_subset_batch_corrected))
elbow_plot <- ElbowPlot(Neuron_subset_batch_corrected, ndims = 50)
ggsave(elbow_plot, filename = "elbow_plot_batch_corrected.png", height = 7, width = 7, type = "cairo")

Neuron_subset_batch_corrected.40 <- FindNeighbors(Neuron_subset_batch_corrected, dims = 1:40)
Neuron_subset_batch_corrected.40nres <- FindClusters(Neuron_subset_batch_corrected.40)
Neuron_subset_batch_corrected.40.01 <- FindClusters(Neuron_subset_batch_corrected.40, resolution = 0.1)
Neuron_subset_batch_corrected.40.02 <- FindClusters(Neuron_subset_batch_corrected.40, resolution = 0.2)
Neuron_subset_batch_corrected.40.03 <- FindClusters(Neuron_subset_batch_corrected.40, resolution = 0.3)
Neuron_subset_batch_corrected.40.04 <- FindClusters(Neuron_subset_batch_corrected.40, resolution = 0.4)
Neuron_subset_batch_corrected.40.05 <- FindClusters(Neuron_subset_batch_corrected.40, resolution = 0.5)

Neuron_subset_batch_corrected.40nres <- RunUMAP(Neuron_subset_batch_corrected.40nres, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_batch_corrected.40.01 <- RunUMAP(Neuron_subset_batch_corrected.40.01, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_batch_corrected.40.02 <- RunUMAP(Neuron_subset_batch_corrected.40.02, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_batch_corrected.40.03 <- RunUMAP(Neuron_subset_batch_corrected.40.03, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_batch_corrected.40.04 <- RunUMAP(Neuron_subset_batch_corrected.40.04, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
Neuron_subset_batch_corrected.40.05 <- RunUMAP(Neuron_subset_batch_corrected.40.05, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

UMAPNeuron_subset_batch_corrected.40nres1 <- DimPlot(Neuron_subset_batch_corrected.40nres, reduction = "umap.rna", label = TRUE, split.by = "batch")
ggsave(UMAPNeuron_subset_batch_corrected.40nres1, filename = "UMAPNeuron_subset_batch_corrected.40nres1.png", height = 7, width = 14, type = "cairo")

UMAPNeuron_subset_batch_corrected.40nres2 <- DimPlot(Neuron_subset_batch_corrected.40nres, reduction = "umap.rna", label = TRUE, group.by = "batch")
ggsave(UMAPNeuron_subset_batch_corrected.40nres2, filename = "UMAPNeuron_subset_batch_corrected.40nres2.png", height = 7, width = 7, type = "cairo")

UMAPNeuron_subset_batch_corrected.40nres3 <- DimPlot(Neuron_subset_batch_corrected.40nres, reduction = "umap.rna", label = TRUE, split.by = "stim")
ggsave(UMAPNeuron_subset_batch_corrected.40nres3, filename = "UMAPNeuron_subset_batch_corrected.40nres3.png", height = 7, width = 14, type = "cairo")

UMAPNeuron_subset_batch_corrected.40nres4 <- DimPlot(Neuron_subset_batch_corrected.40nres, reduction = "umap.rna", label = TRUE, group.by = "stim")
ggsave(UMAPNeuron_subset_batch_corrected.40nres4, filename = "UMAPNeuron_subset_batch_corrected.40nres4.png", height = 7, width = 7, type = "cairo")

UMAPNeuron_subset_batch_corrected.40.01 <- DimPlot(Neuron_subset_batch_corrected.40.01, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_batch_corrected.40.01, filename = "UMAPNeuron_subset_batch_corrected.40.01.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_batch_corrected.40.02 <- DimPlot(Neuron_subset_batch_corrected.40.02, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_batch_corrected.40.02, filename = "UMAPNeuron_subset_batch_corrected.40.02.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_batch_corrected.40.03 <- DimPlot(Neuron_subset_batch_corrected.40.03, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_batch_corrected.40.03, filename = "UMAPNeuron_subset_batch_corrected.40.03.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_batch_corrected.40.04 <- DimPlot(Neuron_subset_batch_corrected.40.04, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_batch_corrected.40.04, filename = "UMAPNeuron_subset_batch_corrected.40.04.png", height = 7, width = 7, type = "cairo")
UMAPNeuron_subset_batch_corrected.40.05 <- DimPlot(Neuron_subset_batch_corrected.40.05, reduction = "umap.rna", label = TRUE)
ggsave(UMAPNeuron_subset_batch_corrected.40.05, filename = "UMAPNeuron_subset_batch_corrected.40.05.png", height = 7, width = 7, type = "cairo")

saveRDS(Neuron_subset_batch_corrected.40.01, "Neuron_subset_batch_corrected.40.01.rds")
saveRDS(Neuron_subset_batch_corrected.40.02, "Neuron_subset_batch_corrected.40.02.rds")
saveRDS(Neuron_subset_batch_corrected.40.03, "Neuron_subset_batch_corrected.40.03.rds")
saveRDS(Neuron_subset_batch_corrected.40.04, "Neuron_subset_batch_corrected.40.04.rds")
saveRDS(Neuron_subset_batch_corrected.40.05, "Neuron_subset_batch_corrected.40.05.rds")








# To subset and remove single cluster and keep the remaining clusters for new analysis
#sub_obj <- subset(object = obj_name, idents = 1, invert = TRUE)
#table(object@meta.data$res.0.8, object@meta.data$orig.ident)