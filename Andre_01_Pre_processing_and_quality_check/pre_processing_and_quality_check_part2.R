library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(cowplot)
library(DoubletFinder)

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/01_Pre_processing_and_quality_check/doublet_cleaning")

merged_filtered <- readRDS("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/01_Pre_processing_and_quality_check/low_quality_cells_cleaning/All_combined_posfilt.rds")
merged_filtered

#Split object to run doubletfinder
merged_filtered_list <- SplitObject(merged_filtered, split.by = "sample")
merged_filtered_list

#SCT normalization and scaling and PCA
merged_filtered_list_sct <- lapply(X = merged_filtered_list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = 3000)
  x <- RunPCA(x)
  x <- FindNeighbors(x, dims = 1:40)
  x <- FindClusters(x, resolution = 0.1)
  x <- RunUMAP(x, dims = 1:40, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
})
merged_filtered_list_sct

saveRDS(merged_filtered_list_sct, "merged_filtered_list_sct.rds")

UmapPlots <- lapply(X = merged_filtered_list_sct, FUN = function(x) DimPlot(x, reduction = "umap.rna", label = TRUE))
lapply(names(UmapPlots), function(x) ggsave(filename=paste(x,"_UmapPlot.png",sep=""), height = 7, width = 7, type = "cairo", plot=UmapPlots[[x]]))  

gene_markers <- c("Cx3cr1", "Tmem119", "Gad1", "Slc17a7", "Pdgfra", "Cspg4", "Atp1a2", "Ndrg2", "Mobp", "Mog")
VlnPlots <- lapply(X = merged_filtered_list_sct, FUN = function(x) VlnPlot(x, features = gene_markers, slot = "counts", log = TRUE))
lapply(names(VlnPlots), function(x) ggsave(filename=paste(x,"_VlnPlot.png",sep=""), height = 7, width = 14, type = "cairo", plot=VlnPlots[[x]]))  

#Run doubletfinder ----------------------------------

sweep.stats.list <- list()
for (i in 1:length(merged_filtered_list_sct)) {
  seu_temp <- merged_filtered_list_sct[[i]]
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats.list[[i]] <- sweep.stats
}

sweep.stats.list <- as.list(sweep.stats.list)

names(sweep.stats.list) <- c("P27_1", "P27_2", "P27_3", "P27_4", "P27_5", "P27_6",
                      "LDR0_1", "LDR0_2", "LDR0_3",
                      "LDR30_1", "LDR30_2", "LDR30_3",
                      "LDR2_1", "LDR2_2", "LDR2_3",
                      "LDR4_1", "LDR4_2", "LDR4_3",
                      "LDR6_1", "LDR6_2", "LDR6_3")

sweep.stats.list.pK <- lapply(sweep.stats.list, FUN = find.pK)

P27_1_pK <- as.numeric(as.character((sweep.stats.list.pK$P27_1[which.max(sweep.stats.list.pK$P27_1$BCmetric),"pK"])))
P27_2_pK <- as.numeric(as.character((sweep.stats.list.pK$P27_2[which.max(sweep.stats.list.pK$P27_2$BCmetric),"pK"])))
P27_3_pK <- as.numeric(as.character((sweep.stats.list.pK$P27_3[which.max(sweep.stats.list.pK$P27_3$BCmetric),"pK"])))
P27_4_pK <- as.numeric(as.character((sweep.stats.list.pK$P27_4[which.max(sweep.stats.list.pK$P27_4$BCmetric),"pK"])))
P27_5_pK <- as.numeric(as.character((sweep.stats.list.pK$P27_5[which.max(sweep.stats.list.pK$P27_5$BCmetric),"pK"])))
P27_6_pK <- as.numeric(as.character((sweep.stats.list.pK$P27_6[which.max(sweep.stats.list.pK$P27_6$BCmetric),"pK"])))
LDR0_1_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR0_1[which.max(sweep.stats.list.pK$LDR0_1$BCmetric),"pK"])))
LDR0_2_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR0_2[which.max(sweep.stats.list.pK$LDR0_2$BCmetric),"pK"])))
LDR0_3_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR0_3[which.max(sweep.stats.list.pK$LDR0_3$BCmetric),"pK"])))
LDR30_1_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR30_1[which.max(sweep.stats.list.pK$LDR30_1$BCmetric),"pK"])))
LDR30_2_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR30_2[which.max(sweep.stats.list.pK$LDR30_2$BCmetric),"pK"])))
LDR30_3_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR30_3[which.max(sweep.stats.list.pK$LDR30_3$BCmetric),"pK"])))
LDR2_1_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR2_1[which.max(sweep.stats.list.pK$LDR2_1$BCmetric),"pK"])))
LDR2_2_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR2_2[which.max(sweep.stats.list.pK$LDR2_2$BCmetric),"pK"])))
LDR2_3_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR2_3[which.max(sweep.stats.list.pK$LDR2_3$BCmetric),"pK"])))
LDR4_1_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR4_1[which.max(sweep.stats.list.pK$LDR4_1$BCmetric),"pK"])))
LDR4_2_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR4_2[which.max(sweep.stats.list.pK$LDR4_2$BCmetric),"pK"])))
LDR4_3_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR4_3[which.max(sweep.stats.list.pK$LDR4_3$BCmetric),"pK"])))
LDR6_1_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR6_1[which.max(sweep.stats.list.pK$LDR6_1$BCmetric),"pK"])))
LDR6_2_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR6_2[which.max(sweep.stats.list.pK$LDR6_2$BCmetric),"pK"])))
LDR6_3_pK <- as.numeric(as.character((sweep.stats.list.pK$LDR6_3[which.max(sweep.stats.list.pK$LDR6_3$BCmetric),"pK"])))

pK.vec <- c(P27_1_pK, P27_2_pK, P27_3_pK, P27_4_pK, P27_5_pK, P27_6_pK,
                             LDR0_1_pK, LDR0_2_pK, LDR0_3_pK,
                             LDR30_1_pK, LDR30_2_pK, LDR30_3_pK,
                             LDR2_1_pK, LDR2_2_pK, LDR2_3_pK,
                             LDR4_1_pK, LDR4_2_pK , LDR4_3_pK ,
                             LDR6_1_pK , LDR6_2_pK , LDR6_3_pK)
pK.vec

for (i in 1:length(merged_filtered_list_sct)) {
  seu_temp <- merged_filtered_list_sct[[i]]
  nExp_poi <- 0.075*nrow(seu_temp@meta.data)
  seu_temp <- doubletFinder_v3(seu_temp, PCs = seu_temp@commands$RunUMAP.SCT.pca$dims, pN = 0.25, pK = pK.vec[i], nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  merged_filtered_list_sct[[i]] <- seu_temp
}

head(merged_filtered_list_sct$P27_Biorep1@meta.data)
head(merged_filtered_list_sct$LDR0_Biorep1@meta.data)
head(merged_filtered_list_sct$LDR30_Biorep1@meta.data)
head(merged_filtered_list_sct$LDR2_Biorep1@meta.data)
head(merged_filtered_list_sct$LDR4_Biorep1@meta.data)
head(merged_filtered_list_sct$LDR6_Biorep1@meta.data)

merged_filtered_list_sct_dbl <- merged_filtered_list_sct

saveRDS(merged_filtered_list_sct_dbl, "merged_filtered_list_sct_dbl.rds")
