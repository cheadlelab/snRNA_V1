# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(monocle3)
library(patchwork)

library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)

################################################################################################################################
###################################################### Seurat -> Cellchat ######################################################
################################################################################################################################
path_string <- "/media/linqianyu/MOSS/snRNAseq_all_nuclei_V1/"
seurat_object <- readRDS(paste0(path_string, "output/Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean2_corrected.rds"))


seurat_object@meta.data$spec_Celltype <- as.factor(seurat_object@meta.data$spec_Celltype)
Idents(seurat_object) <- seurat_object@meta.data$spec_Celltype
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(seurat_object)
meta <- as.data.frame(seurat_object@meta.data)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "spec_Celltype")

################################################################################################################################
################################################################################################################################
################################################################################################################################
# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

################################################################################################################################
################################################################################################################################
################################################################################################################################
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 12) # do parallel
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

################################################################################################################################
################################################################################################################################
################################################################################################################################
# save and load cellchat to save memory
# turn off Rstudio after saving
# restart Rstudio and reload cellchat
path_string <- "/media/linqianyu/MOSS/snRNAseq_all_nuclei_V1/"
save(cellchat, file = paste0(path_string, "output/step5_CellChat_tmp.RData"))

path_string <- "/media/linqianyu/MOSS/snRNAseq_all_nuclei_V1/"
load(paste0(path_string, "output/step5_CellChat_tmp.RData"))

future::plan("multisession", workers = 6) # do parallel
# set to 10 and computer will die

################################################################################################################################
################################################################################################################################
################################################################################################################################
# Compute the communication probability and infer cellular communication network
options(future.globals.maxSize = 1000 * 1024^2)  # Set limit to 1000 MiB
cellchat <- computeCommunProb(cellchat)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-1https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html1-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

save(cellchat, file = paste0(path_string, "output/step5_CellChat.RData"))

################################################################################################################################
################################################################################################################################
################################################################################################################################
path_string <- "/media/linqianyu/MOSS/snRNAseq_all_nuclei_V1/"
load(paste0(path_string, "output/step5_CellChat.RData"))














