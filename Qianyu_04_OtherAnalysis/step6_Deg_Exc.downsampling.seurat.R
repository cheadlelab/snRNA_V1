rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(venn)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

################################################################################################################################
performGOAnalysis <- function(genes, direction, condition, path_string, prefix) {
  if (length(genes) > 0) {
    ego <- enrichGO(gene = genes,
                    OrgDb = org.Mm.eg.db, 
                    keyType = "SYMBOL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)
    if (!is.null(ego) && dim(ego)[1] > 0) {
      pdf_path <- paste0(path_string, "output/step6_DEGs_seurat/", prefix, ".GO.", condition, ".", direction, ".pdf")
      pdf(pdf_path, width = 10, height = 10)
      print(barplot(ego, showCategory = 20) + ggtitle(paste(prefix, "genes: ", direction, "regulated in ", condition)))
      dev.off()
    }
  }
}

################################################################################################################################
############################################################# Load #############################################################
################################################################################################################################
path_string <- "/home/qianyu/Desktop/grid/snRNAseq_all_nuclei_V1/"
script_path <- "/home/qianyu/Desktop/grid/snRNAseq_all_nuclei_V1/script/VolcanoPlot.R"
source(script_path)
V1 <- readRDS(paste0(path_string, "output/Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean2_corrected.rds"))
# V1@assays$RNA <- V1@assays$decontXcounts

set.seed(666)
# Create a combined grouping variable
V1@meta.data$group <- paste(V1$condition_simpl, V1$spec_Celltype, sep = "_")

# Downsample cells for each group (orig.ident + spec_Celltype)
cells_to_keep <- unlist(lapply(split(rownames(V1@meta.data), V1@meta.data$group), function(cells) {
  sample(cells, size = min(400, length(cells))) # Adjust the fraction as needed
}))

# Subset the Seurat object with the selected cells
V1_downsampled <- subset(V1, cells = cells_to_keep)
table(V1_downsampled@meta.data$condition_simpl, V1_downsampled@meta.data$spec_Celltype)
# Aggregate expression
pseudo_V1 <- AggregateExpression(
  V1_downsampled,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("condition_simpl", "spec_Celltype", "batch")
)

DEG <- list()


######################################################################################################################################
############################################################# DEG and GO #############################################################
######################################################################################################################################
# Loop through each neuron type
for (celltype in c("Exc L6IT", "Exc L6a", "Exc L5PT", "Exc L5IT", "Exc L5NP", "Exc L6b")) {
  obj <- subset(pseudo_V1, subset = spec_Celltype %in% celltype)
  Idents(obj) <- "condition_simpl"
  conditions <- setdiff(unique(obj@meta.data[["condition_simpl"]]), "LDR")
  
  # Determine the prefix based on cell type for naming consistency
  prefix <- substring(celltype, 5)
  DEG[[prefix]] <- list()
  
  deg_list <- list()
  go_list <- list()
  
  # DEGs analysis loop
  for (condition in conditions) {
    tmp.V1 <- subset(obj, subset = condition_simpl %in% c(condition, "LDR"))
    bulk.mono.de <- FindMarkers(object = tmp.V1, ident.1 = condition, ident.2 = "LDR", test.use = "DESeq2", logfc.threshold = log2(1.5))
    bulk.mono.de <- na.omit(bulk.mono.de)
    DEG[[prefix]][[condition]] <- bulk.mono.de[bulk.mono.de$p_val_adj < 0.05 & abs(bulk.mono.de$avg_log2FC) > log2(1.5),]
    # Save DEG CSV with appropriate prefix
    csv_path <- paste0(path_string, "output/step6_DEGs_seurat/", prefix, ".DEG.", condition, ".downsampling.csv")
    write.csv(bulk.mono.de[bulk.mono.de$p_val_adj < 0.05 & abs(bulk.mono.de$avg_log2FC) > log2(1.5),], csv_path, row.names = TRUE)
    
    # Prepare data for volcano plot
    bulk.mono.de$log2FoldChange <- bulk.mono.de$avg_log2FC
    bulk.mono.de$padj <- bulk.mono.de$p_val_adj
    bulk.mono.de$symbol <- rownames(bulk.mono.de)
    bulk.mono.de$Condition <- condition
    
    # Volcano plot PDF with appropriate prefix
    pdf_path <- paste0(path_string, "output/step6_DEGs_seurat/", prefix, ".Volcano.", condition, ".downsampling.pdf")
    cairo_pdf(pdf_path, width = 8, height = 8)
    print(VolcanoPlot(dif = bulk.mono.de, title = paste0(": ", condition, ".vs.LDR")))
    dev.off()
    
    deg_list[[condition]] <- bulk.mono.de
  }
  
  # GO analysis loop
  # Inside GO analysis loop
  for (condition in conditions) {
    dif <- deg_list[[condition]] # Use deg_list with dynamic referencing
    upregulated_genes <- dif$symbol[dif$padj < 0.05 & dif$log2FoldChange > log2(1.5)]
    downregulated_genes <- dif$symbol[dif$padj < 0.05 & dif$log2FoldChange < -log2(1.5)]
    
    # Perform GO analysis for upregulated genes
    performGOAnalysis(upregulated_genes, "up", condition, path_string, prefix)
    # Perform GO analysis for downregulated genes
    performGOAnalysis(downregulated_genes, "down", condition, path_string, prefix)
  }
}
