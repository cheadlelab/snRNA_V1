library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)

setwd("/grid/cheadle/home/machado/scRNAseq_V1/Andre-snRNAseq/R-analysis/01_Data_preparation")

# Use for FeatureMatrix and FindMarkers in Seurat -----------------------------------------------------------------------------------------------------------------------
library(future)
plan("multicore", workers = 10)
options(future.globals.maxSize = 16000 * 1024^2)

# Load in individual feature matrix files for each samples -------------------------------------------------------------------------------------------
P27_1 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX03_p27/outs/filtered_feature_bc_matrix")
P27_2 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX05_p27/outs/filtered_feature_bc_matrix")
P27_3 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX06_P27/outs/filtered_feature_bc_matrix")
LDR0_1 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX03_LDR0/outs/filtered_feature_bc_matrix")
LDR0_2 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX05_LDR0/outs/filtered_feature_bc_matrix")
LDR0_3 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX06_LDR0/outs/filtered_feature_bc_matrix")
LDR2_1 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX03_LDR2/outs/filtered_feature_bc_matrix")
LDR2_2 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX05_LDR2/outs/filtered_feature_bc_matrix")
LDR2_3 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX06_LDR2/outs/filtered_feature_bc_matrix")
LDR6_1 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX03_LDR6/outs/filtered_feature_bc_matrix")
LDR6_2 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX05_LDR6/outs/filtered_feature_bc_matrix")
LDR6_3 <- Read10X(data.dir = "/grid/cheadle/home/machado/scRNAseq_V1/MIA-snRNAseq/count/Cheadle_AMX06_LDR6/outs/filtered_feature_bc_matrix")

# Create Seurat objects for Gene expression -------------------------------------------------------------------------------------------
P27_1_seurat <- CreateSeuratObject(counts = P27_1, assay = "RNA", project = "P27rep1")
P27_2_seurat <- CreateSeuratObject(counts = P27_2, assay = "RNA", project = "P27rep2")
P27_3_seurat <- CreateSeuratObject(counts = P27_3, assay = "RNA", project = "P27rep3")

LDR2_1_seurat <- CreateSeuratObject(counts = LDR2_1, assay = "RNA", project = "LDR2rep1")
LDR2_2_seurat <- CreateSeuratObject(counts = LDR2_2, assay = "RNA", project = "LDR2rep2")
LDR2_3_seurat <- CreateSeuratObject(counts = LDR2_3, assay = "RNA", project = "LDR2rep3")

LDR0_1_seurat <- CreateSeuratObject(counts = LDR0_1, assay = "RNA", project = "LDR0rep1")
LDR0_2_seurat <- CreateSeuratObject(counts = LDR0_2, assay = "RNA", project = "LDR0rep2")
LDR0_3_seurat <- CreateSeuratObject(counts = LDR0_3, assay = "RNA", project = "LDR0rep3")

LDR6_1_seurat <- CreateSeuratObject(counts = LDR6_1, assay = "RNA", project = "LDR6rep1")
LDR6_2_seurat <- CreateSeuratObject(counts = LDR6_2, assay = "RNA", project = "LDR6rep2")
LDR6_3_seurat <- CreateSeuratObject(counts = LDR6_3, assay = "RNA", project = "LDR6rep3")

#Split samples by group and add metadata

#P27_1
samplename = P27_1_seurat@meta.data$orig.ident
table(samplename)

stim = rep("CTRL",length(samplename))
batch = rep("batch1",length(samplename))

names(stim) = rownames(P27_1_seurat@meta.data)
names(batch) = rownames(P27_1_seurat@meta.data)

P27_1_seurat <- AddMetaData(
  object = P27_1_seurat,
  metadata = stim,
  col.name = "stim")
P27_1_seurat <- AddMetaData(
  object = P27_1_seurat,
  metadata = batch,
  col.name = "batch")

P27_1_seurat@meta.data
table(P27_1_seurat@meta.data$stim)
table(P27_1_seurat@meta.data$batch)
#P27_2
samplename = P27_2_seurat@meta.data$orig.ident
table(samplename)

stim = rep("CTRL",length(samplename))
batch = rep("batch2",length(samplename))

names(stim) = rownames(P27_2_seurat@meta.data)
names(batch) = rownames(P27_2_seurat@meta.data)

P27_2_seurat <- AddMetaData(
  object = P27_2_seurat,
  metadata = stim,
  col.name = "stim")
P27_2_seurat <- AddMetaData(
  object = P27_2_seurat,
  metadata = batch,
  col.name = "batch")

P27_2_seurat@meta.data
table(P27_2_seurat@meta.data$stim)
table(P27_2_seurat@meta.data$batch)
#P27_3
samplename = P27_3_seurat@meta.data$orig.ident
table(samplename)

stim = rep("CTRL",length(samplename))
batch = rep("batch3",length(samplename))

names(stim) = rownames(P27_3_seurat@meta.data)
names(batch) = rownames(P27_3_seurat@meta.data)

P27_3_seurat <- AddMetaData(
  object = P27_3_seurat,
  metadata = stim,
  col.name = "stim")
P27_3_seurat <- AddMetaData(
  object = P27_3_seurat,
  metadata = batch,
  col.name = "batch")

P27_3_seurat@meta.data
table(P27_3_seurat@meta.data$stim)
table(P27_3_seurat@meta.data$batch)


#LDR0
#LDR0_1
samplename = LDR0_1_seurat@meta.data$orig.ident
table(samplename)

stim = rep("NO_STIM",length(samplename))
batch = rep("batch1",length(samplename))

names(stim) = rownames(LDR0_1_seurat@meta.data)
names(batch) = rownames(LDR0_1_seurat@meta.data)

LDR0_1_seurat <- AddMetaData(
  object = LDR0_1_seurat,
  metadata = stim,
  col.name = "stim")
LDR0_1_seurat <- AddMetaData(
  object = LDR0_1_seurat,
  metadata = batch,
  col.name = "batch")

LDR0_1_seurat@meta.data
table(LDR0_1_seurat@meta.data$stim)
table(LDR0_1_seurat@meta.data$batch)
#LDR0_2
samplename = LDR0_2_seurat@meta.data$orig.ident
table(samplename)

stim = rep("NO_STIM",length(samplename))
batch = rep("batch2",length(samplename))

names(stim) = rownames(LDR0_2_seurat@meta.data)
names(batch) = rownames(LDR0_2_seurat@meta.data)

LDR0_2_seurat <- AddMetaData(
  object = LDR0_2_seurat,
  metadata = stim,
  col.name = "stim")
LDR0_2_seurat <- AddMetaData(
  object = LDR0_2_seurat,
  metadata = batch,
  col.name = "batch")

LDR0_2_seurat@meta.data
table(LDR0_2_seurat@meta.data$stim)
table(LDR0_2_seurat@meta.data$batch)
#LDR0_3
samplename = LDR0_3_seurat@meta.data$orig.ident
table(samplename)

stim = rep("NO_STIM",length(samplename))
batch = rep("batch3",length(samplename))

names(stim) = rownames(LDR0_3_seurat@meta.data)
names(batch) = rownames(LDR0_3_seurat@meta.data)

LDR0_3_seurat <- AddMetaData(
  object = LDR0_3_seurat,
  metadata = stim,
  col.name = "stim")
LDR0_3_seurat <- AddMetaData(
  object = LDR0_3_seurat,
  metadata = batch,
  col.name = "batch")

LDR0_3_seurat@meta.data
table(LDR0_3_seurat@meta.data$stim)
table(LDR0_3_seurat@meta.data$batch)

#LDR2
#LDR2_1
samplename = LDR2_1_seurat@meta.data$orig.ident
table(samplename)

stim = rep("STIM2h",length(samplename))
batch = rep("batch1",length(samplename))

names(stim) = rownames(LDR2_1_seurat@meta.data)
names(batch) = rownames(LDR2_1_seurat@meta.data)

LDR2_1_seurat <- AddMetaData(
  object = LDR2_1_seurat,
  metadata = stim,
  col.name = "stim")
LDR2_1_seurat <- AddMetaData(
  object = LDR2_1_seurat,
  metadata = batch,
  col.name = "batch")

LDR2_1_seurat@meta.data
table(LDR2_1_seurat@meta.data$stim)
table(LDR2_1_seurat@meta.data$batch)
#LDR2_2
samplename = LDR2_2_seurat@meta.data$orig.ident
table(samplename)

stim = rep("STIM2h",length(samplename))
batch = rep("batch2",length(samplename))

names(stim) = rownames(LDR2_2_seurat@meta.data)
names(batch) = rownames(LDR2_2_seurat@meta.data)

LDR2_2_seurat <- AddMetaData(
  object = LDR2_2_seurat,
  metadata = stim,
  col.name = "stim")
LDR2_2_seurat <- AddMetaData(
  object = LDR2_2_seurat,
  metadata = batch,
  col.name = "batch")

LDR2_2_seurat@meta.data
table(LDR2_2_seurat@meta.data$stim)
table(LDR2_2_seurat@meta.data$batch)
#LDR2_3
samplename = LDR2_3_seurat@meta.data$orig.ident
table(samplename)

stim = rep("STIM2h",length(samplename))
batch = rep("batch3",length(samplename))

names(stim) = rownames(LDR2_3_seurat@meta.data)
names(batch) = rownames(LDR2_3_seurat@meta.data)

LDR2_3_seurat <- AddMetaData(
  object = LDR2_3_seurat,
  metadata = stim,
  col.name = "stim")
LDR2_3_seurat <- AddMetaData(
  object = LDR2_3_seurat,
  metadata = batch,
  col.name = "batch")

LDR2_3_seurat@meta.data
table(LDR2_3_seurat@meta.data$stim)
table(LDR2_3_seurat@meta.data$batch)

#LDR6
#LDR6_1
samplename = LDR6_1_seurat@meta.data$orig.ident
table(samplename)

stim = rep("STIM6h",length(samplename))
batch = rep("batch1",length(samplename))

names(stim) = rownames(LDR6_1_seurat@meta.data)
names(batch) = rownames(LDR6_1_seurat@meta.data)

LDR6_1_seurat <- AddMetaData(
  object = LDR6_1_seurat,
  metadata = stim,
  col.name = "stim")
LDR6_1_seurat <- AddMetaData(
  object = LDR6_1_seurat,
  metadata = batch,
  col.name = "batch")

LDR6_1_seurat@meta.data
table(LDR6_1_seurat@meta.data$stim)
table(LDR6_1_seurat@meta.data$batch)
#LDR6_2
samplename = LDR6_2_seurat@meta.data$orig.ident
table(samplename)

stim = rep("STIM6h",length(samplename))
batch = rep("batch2",length(samplename))

names(stim) = rownames(LDR6_2_seurat@meta.data)
names(batch) = rownames(LDR6_2_seurat@meta.data)

LDR6_2_seurat <- AddMetaData(
  object = LDR6_2_seurat,
  metadata = stim,
  col.name = "stim")
LDR6_2_seurat <- AddMetaData(
  object = LDR6_2_seurat,
  metadata = batch,
  col.name = "batch")

LDR6_2_seurat@meta.data
table(LDR6_2_seurat@meta.data$stim)
table(LDR6_2_seurat@meta.data$batch)
#LDR6_3
samplename = LDR6_3_seurat@meta.data$orig.ident
table(samplename)

stim = rep("STIM6h",length(samplename))
batch = rep("batch3",length(samplename))

names(stim) = rownames(LDR6_3_seurat@meta.data)
names(batch) = rownames(LDR6_3_seurat@meta.data)

LDR6_3_seurat <- AddMetaData(
  object = LDR6_3_seurat,
  metadata = stim,
  col.name = "stim")
LDR6_3_seurat <- AddMetaData(
  object = LDR6_3_seurat,
  metadata = batch,
  col.name = "batch")

LDR6_3_seurat@meta.data
table(LDR6_3_seurat@meta.data$stim)
table(LDR6_3_seurat@meta.data$batch)

#Save the files for pre-processing ----------------------------------------------------------------------
saveRDS(P27_1_seurat, "P27_1_seurat_v1.rds")
saveRDS(P27_2_seurat, "P27_2_seurat_v1.rds")
saveRDS(P27_3_seurat, "P27_3_seurat_v1.rds")
saveRDS(LDR0_1_seurat, "LDR0_1_seurat_v1.rds")
saveRDS(LDR0_2_seurat, "LDR0_2_seurat_v1.rds")
saveRDS(LDR0_3_seurat, "LDR0_3_seurat_v1.rds")
saveRDS(LDR2_1_seurat, "LDR2_1_seurat_v1.rds")
saveRDS(LDR2_2_seurat, "LDR2_2_seurat_v1.rds")
saveRDS(LDR2_3_seurat, "LDR2_3_seurat_v1.rds")
saveRDS(LDR6_1_seurat, "LDR6_1_seurat_v1.rds")
saveRDS(LDR6_2_seurat, "LDR6_2_seurat_v1.rds")
saveRDS(LDR6_3_seurat, "LDR6_3_seurat_v1.rds")



# Merge files into single object
# merge_test <- merge(P27_1_seurat, y = c(P27_2_seurat, P27_3_seurat))
# total_data <- merge(merge_test, y= c(LDR0_1_seurat, LDR0_2_seurat, LDR0_3_seurat, LDR2_1_seurat, LDR2_2_seurat, LDR2_3_seurat, LDR6_1_seurat, LDR6_2_seurat, LDR6_3_seurat))
