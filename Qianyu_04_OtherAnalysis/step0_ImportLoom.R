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
###################################################### Path ######################################################
##################################################################################################################

exp1 <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home/qianyu/snRNAseq_all_nuclei_V1/cellranger_output/Exp1"
exp1_subdirectories <- list.dirs(path = exp1, full.names = TRUE, recursive = FALSE)
exp2 <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home/qianyu/snRNAseq_all_nuclei_V1/cellranger_output/Exp2"
exp2_subdirectories <- list.dirs(path = exp2, full.names = TRUE, recursive = FALSE)

append_path <- function(directory) {
  folder_name <- basename(directory)
  return(paste0(directory, "/velocyto/", folder_name, ".loom"))
}
exp1_modified <- sapply(exp1_subdirectories, append_path)
exp2_modified <- sapply(exp2_subdirectories, append_path)

data_path <- c(exp1_modified, exp2_modified)
rm(exp1, exp2, exp1_subdirectories, exp2_subdirectories, exp1_modified, exp2_modified)

##################################################################################################################
###################################################### Load ######################################################
##################################################################################################################

result_list <- lapply(data_path, read.loom.matrices)
names(result_list) <- tools::file_path_sans_ext(basename(data_path))
seurat_list <- lapply(result_list, as.Seurat)

update_metadata <- function(seurat_obj) {
  # Get the name of the Seurat object
  name <- deparse(substitute(seurat_obj))
  # Use strsplit to split the name
  parts <- unlist(strsplit(name, "_"))
  # Extract batch and exp from the split parts
  batch <- parts[length(parts) - 1]  # second to last part
  exp <- parts[length(parts)]         # last part
  # Update the Seurat object's metadata
  seurat_obj[["batch"]] <- batch
  seurat_obj[["exp"]] <- exp
  
  return(seurat_obj)
}

# Apply the function to each Seurat object in seurat_list
seurat_obj <- lapply(names(seurat_list), function(x) {
  obj <- seurat_list[[x]]
  parts <- unlist(strsplit(x, "_"))
  batch <- parts[length(parts) - 1]  # second to last part
  exp <- parts[length(parts)]         # last part
  obj[["batch"]] <- batch
  obj[["exp"]] <- exp
  return(obj)
})

# Use the original names for the updated list
names(seurat_obj) <- names(seurat_list)

# clean
rm(list = setdiff(ls(), "seurat_obj"))
gc()


##################################################################################################################
###################################################### Save ######################################################
##################################################################################################################

path_string <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home/qianyu/snRNAseq_all_nuclei_V1/"
save(seurat_obj, file = paste0(path_string, "output/step0_ImportlLoom.RData"))
