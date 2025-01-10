## match the data from total RNA analysis
## quite complicated, if want to apply on other data, pls check and modify the code
rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(dplyr)
library(velocyto.R)
library(rhdf5)
library(celda)

##################################################################################################################
###################################################### Load ######################################################
##################################################################################################################

path_backup <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction/Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean2.rds"
path_string <- "/media/linqianyu/MOSS/snRNAseq_all_nuclei_V1/"
load(paste0(path_string, "output/step1_DecontX.RData"))
V1_obj <- readRDS(paste0(path_string, "output/Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean2.rds"))

## merge
for (name in names(seurat_obj)) {
  seurat_obj[[name]]@meta.data$orig.ident <- name
}

first_obj <- seurat_obj[[1]]
rest_objs <- seurat_obj[2:length(seurat_obj)]
merged_seurat_obj <- merge(first_obj, 
                           y = unlist(rest_objs, recursive = FALSE),
                           add.cell.ids = names(seurat_obj),
                           project = "V1",
                           merge.data = TRUE)

########################################################################################################################
###################################################### Match Name ######################################################
########################################################################################################################
seurat_obj <- merged_seurat_obj

## match
tmp <- colnames(V1_obj$RNA)
V1_obj@meta.data$simplify <- gsub(".{18}$", "", tmp)

# Generate a cross table
cross_tab <- table(V1_obj@meta.data$simplify, V1_obj@meta.data$orig.ident)
# Initialize flag for successful check
check_success <- TRUE
# Generate NameTransfer_Ref
unique_names <- names(table(seurat_obj@meta.data$orig.ident))
NameTransfer_Ref <- setNames(vector("list", length(unique_names)), unique_names)

# For each non-zero element in the table, check the conditions
for (i in seq_len(nrow(cross_tab))) {
  for (j in seq_len(ncol(cross_tab))) {
    if (cross_tab[i, j] != 0) {
      col_name <- colnames(cross_tab)[j]
      NameTransfer_Ref[col_name] <- rownames(cross_tab)[i]
      if ((sum(cross_tab[i, ]) - cross_tab[i,j] != 0) || (sum(cross_tab[, j]) - cross_tab[i,j] != 0)) {
        check_success <- FALSE
        break
      }
    }
  }
  if (!check_success) {
    break
  }
}

if (!check_success) {
  stop("The elements of 'simplify' and 'orig.ident' do not correspond to each other.")
} else {
  cat("Name check success.\n")
}

## rename seurat_obj@meta.data$cells
seurat_obj@meta.data$cells <- colnames(seurat_obj@assays$spliced)
seurat_obj@meta.data$cells <- gsub("x","-1",seurat_obj@meta.data$cells)
seurat_obj@meta.data$cells <- gsub(".*(?=.{18}$)", "", seurat_obj@meta.data$cells, perl = TRUE)
seurat_obj@meta.data$orig.ident2 <- NameTransfer_Ref[seurat_obj@meta.data$orig.ident]

# 1. Replace NULL with "null" in NameTransfer_Ref
NameTransfer_Ref[sapply(NameTransfer_Ref, is.null)] <- "null"
# 2. Use unlist
seurat_obj@meta.data$orig.ident2 <- unlist(NameTransfer_Ref[seurat_obj@meta.data$orig.ident])
# rename finished
seurat_obj@meta.data$cells <- paste0(seurat_obj@meta.data$orig.ident2, seurat_obj@meta.data$cells)


########################################################################################################################
###################################################### Match Cell ######################################################
########################################################################################################################
# Find the match indices
match_indices <- match(seurat_obj@meta.data$cells, V1_obj@meta.data$cells)

# Replace NA values with "null" or NA based on your preference
match_indices[is.na(match_indices)] <- "null"

# Store the indices in the 'match' column
seurat_obj@meta.data$match <- match_indices

###########################################################################################################################
###################################################### Match Reorder ######################################################
###########################################################################################################################
# Subset the metadata to remove rows where match is "null"
tmp_meta <- subset(seurat_obj@meta.data, match != "null")

# Obtain the sorted order of the rows based on the match values
sorted_cells <- rownames(tmp_meta)[order(as.numeric(tmp_meta$match))]

# Subset and sort the Seurat object based on the sorted order
matched_seurat <- subset(seurat_obj, cells = sorted_cells)

#######################################################################################################################
###################################################### Copy Info ######################################################
#######################################################################################################################

matched_seurat@assays$RNA <- matched_seurat@assays$spliced
DefaultAssay(matched_seurat) <- "RNA"
matched_seurat@meta.data$gen_CellType <- V1_obj@meta.data$gen_CellType
matched_seurat@meta.data$spec_Celltype <- V1_obj@meta.data$spec_Celltype
matched_seurat@meta.data$spec_Celltype_cluster <- V1_obj@meta.data$spec_Celltype_cluster

# just give the object of obsm umap for scvelo
matched_seurat <- FindVariableFeatures(matched_seurat, selection.method = "vst", nfeatures = 2000)
matched_seurat <- ScaleData(matched_seurat)
matched_seurat <- RunPCA(matched_seurat)
ElbowPlot(matched_seurat)
dims_parameter <- 15
matched_seurat <- RunUMAP(matched_seurat, dims = 1:dims_parameter)
matched_seurat <- FindNeighbors(matched_seurat, dims = 1:dims_parameter)
matched_seurat <- FindClusters(matched_seurat, resolution = 0.5)
DimPlot(matched_seurat, reduction = "umap",label = TRUE, raster=FALSE)
DimPlot(matched_seurat, reduction = "umap",label = TRUE, group.by = "gen_CellType", raster=FALSE)
DimPlot(matched_seurat, reduction = "umap",label = TRUE, group.by = "spec_Celltype", raster=FALSE)
# replace the umap coordination
matched_seurat@reductions$umap@cell.embeddings[] <- V1_obj@reductions$UMAP@cell.embeddings

# use raw RNA assay for scvelo, otherwise with error in scvelo
matched_seurat@assays$RNA <- matched_seurat@assays$spliced
DefaultAssay(matched_seurat) <- "RNA"

##################################################################################################################
###################################################### Save ######################################################
##################################################################################################################

save(matched_seurat, file = paste0(path_string, "output/step2_Match.RData"))
# copy and paste the file manually
# unknown bug, just save locally and copy and paste
SaveH5Seurat(matched_seurat, filename = "/home/linqianyu/Desktop/step2_Match.h5Seurat")
# SaveH5Seurat(matched_seurat, filename = paste0(path_string, "output/step2_Match.h5Seurat"))
Convert("/home/linqianyu/Desktop/step2_Match.h5Seurat", dest = "h5ad")
#Convert(paste0(path_string, "output/step2_Match.h5Seurat"), dest = "h5ad")










