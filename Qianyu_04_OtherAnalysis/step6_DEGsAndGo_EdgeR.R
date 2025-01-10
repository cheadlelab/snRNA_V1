# DEGs Vocanno
# Venne
# GO/KEGG Enrichment

rm(list = ls(all = TRUE))
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(egg)
library(ggplot2)
library(monocle3)
library(tidyverse)
library(edgeR)

################################################################################################################################
############################################################# Load #############################################################
################################################################################################################################
path_string <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/04_Differential_expression_analysis"
save_path <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/qianyu/snRNAseq_all_nuclei_V1/figures/vocano_edgeR/"
script_path <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/qianyu/snRNAseq_all_nuclei_V1/script/VolcanoPlot.R"
source(script_path)

Allexp <- readRDS(paste0(path_string, "/no_batch_correction/Monocle3/Allexp.cds_monocle3_DEG.rds"))
Allexp <- Allexp[, Allexp@colData@listData$gen_CellType %in% c("GABAergic Neurons", "Glutamatergic Neurons")]
# Allexp <- Allexp[, Allexp@colData@listData$gen_CellType %in% c("Glia cells")]
# Allexp <- Allexp[, Allexp@colData@listData$spec_Celltype %in% c("Olg","Ast","Mg")]
# Allexp <- Allexp[, Allexp@colData@listData$spec_Celltype %in% c("Olg","Ast","Mg")]

# Results list of GO analysis
results_list <- list()


################################################################################################################################
########################################################## PreProcess ##########################################################
################################################################################################################################
# Get all unique spec_Celltype values
all_cell_types <- unique(Allexp@colData@listData$spec_Celltype)
# Initialize a list to store results for each cell type
all_results_list <- list()

for(cell_type in all_cell_types) {  # cell_type <- all_cell_types[1]
  # Extract cells of a specific cell type
  tmp <- Allexp[, Allexp@colData@listData$spec_Celltype == cell_type]
  # Extract the count matrix
  DGE_DATA_tmp <- tmp@assays@data$counts
  # Compute pseudobulk
  mm_tmp <- Matrix::sparse.model.matrix(~0 + tmp$sample)
  pseudobulk_tmp <- DGE_DATA_tmp %*% mm_tmp
  # Define the groups for edgeR
  bulk.labels_tmp = c(rep("P27", 6), rep("LDR0", 3), rep("LDR30m", 3),  rep("LDR2h", 3), rep("LDR4h", 3), rep("LDR6h", 3))
  # Factor and Change the order
  bulk.labels_tmp <- factor(bulk.labels_tmp, levels = c("LDR0", "LDR30m", "LDR2h", "LDR4h", "LDR6h", "P27"))
  # Create a DGEList for differential expression analysis
  dge.list_tmp <- DGEList(counts = pseudobulk_tmp, group = bulk.labels_tmp)
  # Filter out genes with low expression
  keep_tmp <- filterByExpr(dge.list_tmp)
  dge.list_tmp <- dge.list_tmp[keep_tmp, , keep.lib.sizes = FALSE]
  # Normalize the data
  dge.list_tmp <- calcNormFactors(dge.list_tmp)
  # Create the design matrix
  design_tmp = model.matrix(~bulk.labels_tmp)
  # Voom transformation considering quality weights
  v.dge.list_tmp <- voomWithQualityWeights(dge.list_tmp, design_tmp)
  # Fit a linear model
  fit.dge.list_tmp <- lmFit(v.dge.list_tmp)
  # Apply empirical Bayes moderation to the model
  fit.dge.list_tmp <- eBayes(fit.dge.list_tmp, robust=TRUE)
  
  ################################################################################################################################
  ############################################################# Loop #############################################################
  ################################################################################################################################
  # Loop for each condition and gather the results
  conditions <- c("LDR0", "LDR30m", "LDR2h", "LDR4h", "LDR6h", "P27")
  
  # Initialize a list to store the pairwise comparison results
  pairwise_results <- list()
  
  # Loop through each pair of conditions for pairwise comparison
  for (i in 1:(length(conditions) - 1)) {
    for (j in (i+1):length(conditions)) {
      coef_index <- which(levels(dge.list_tmp$samples$group) == conditions[j])  # Updated here
      res_name <- paste(conditions[i], "vs", conditions[j], sep=".")
      pairwise_results[[res_name]] <- topTable(fit.dge.list_tmp, sort.by="p", n=Inf, coef=coef_index)  # Updated here
      pairwise_results[[res_name]]$comparison <- res_name  # This column describes the comparison
      pairwise_results[[res_name]] <- tibble::rownames_to_column(pairwise_results[[res_name]], "gene")
      
      dif <- data.frame(
        symbol = pairwise_results[[res_name]]$gene,
        log2FoldChange = pairwise_results[[res_name]]$logFC,
        padj = pairwise_results[[res_name]]$P.Value
      )
      ################################################################################################################################
      # Generate the volcano plot
      p1 <- VolcanoPlot(dif, padj=0.05, title=res_name, label.max = 50)
      # Save the volcano plot
      file_name <- paste0(save_path, gsub("/", "-", cell_type), "_", res_name, ".pdf")
      ggsave(file_name, egg::set_panel_size(p1, width=unit(66, "mm"), height=unit(66, "mm")), width = 110, height = 110, units = 'mm', dpi = 300)
      
      ################################################################################################################################
      # GO analysis
      # Select upregulated and downregulated genes
      upregulated_genes <- dif$symbol[dif$padj < 0.05 & dif$log2FoldChange > 0]
      downregulated_genes <- dif$symbol[dif$padj < 0.05 & dif$log2FoldChange < 0]
      
      # Perform GO enrichment analysis for upregulated genes
      ego_up <- enrichGO(gene         = upregulated_genes,
                         OrgDb        = org.Mm.eg.db, 
                         keyType      = "SYMBOL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)
      
      # Perform GO enrichment analysis for downregulated genes
      ego_down <- enrichGO(gene         = downregulated_genes,
                           OrgDb        = org.Mm.eg.db, 
                           keyType      = "SYMBOL",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)
      
      # Save results
      results_name_up <- paste(cell_type, conditions[i], "vs", conditions[j], "upregulated", sep = "_")
      results_name_down <- paste(cell_type, conditions[i], "vs", conditions[j], "downregulated", sep = "_")
      
      results_list[[results_name_up]] <- ego_up
      results_list[[results_name_down]] <- ego_down
      
      
      
      ################################################################################################################################
      # save csv
      write.csv(pairwise_results[[res_name]], file = paste0("/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/qianyu/snRNAseq_all_nuclei_V1/output/step6_DEGs_edgeR/",
                                                            gsub("/", "-", cell_type),"_",res_name, ".csv"))
    }
  }
  
  # Combining all pairwise results into a single data frame
  results_df <- Reduce(dplyr::full_join, pairwise_results)
  all_results_list[[cell_type]] <- results_df
  
}






# A function to break the text at the space closest to the middle
break_long_text <- function(string, max_len = 30) {
  # If the string is shorter than max_len or has no spaces, return it as is
  if (nchar(string) <= max_len || !grepl(" ", string)) {
    return(string)
  }
  
  # Find the middle position
  mid_pos <- floor(nchar(string) / 2)
  
  # Find the nearest spaces before and after the midpoint
  spaces_before <- which(strsplit(string, '')[[1]][1:mid_pos] == ' ')
  spaces_after <- which(strsplit(string, '')[[1]][(mid_pos+1):nchar(string)] == ' ') + mid_pos
  
  # If no spaces were found before or after the midpoint, return the original string
  if (length(spaces_before) == 0 && length(spaces_after) == 0) {
    return(string)
  } else if (length(spaces_before) == 0) {
    position_to_break <- min(spaces_after)
  } else if (length(spaces_after) == 0) {
    position_to_break <- max(spaces_before)
  } else {
    space_before <- max(spaces_before)
    space_after <- min(spaces_after)
    
    # Choose the closest space to the midpoint
    if (abs(mid_pos - space_before) <= abs(mid_pos - space_after)) {
      position_to_break <- space_before
    } else {
      position_to_break <- space_after
    }
  }
  
  # Return the string with a newline at the chosen break point
  return(paste0(substr(string, 1, position_to_break - 1), "\n", substr(string, position_to_break + 1, nchar(string))))
}

############################################################################################################################################
############################################################# GO result output #############################################################
############################################################################################################################################
# Define the output path
output_path <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/qianyu/snRNAseq_all_nuclei_V1/figures/GO_edgeR/"

# Iterate over each cell type result
for (celltype in all_cell_types) {
  for (i in 1:(length(conditions) - 1)) {
    for (j in (i + 1):length(conditions)) {
      condition_1 <- conditions[i]
      condition_2 <- conditions[j]
      
      results_name_up <- paste(celltype, condition_1, "vs", condition_2, "upregulated", sep = "_")
      results_name_down <- paste(celltype, condition_1, "vs", condition_2, "downregulated", sep = "_")
      
      ego_up <- results_list[[results_name_up]]
      ego_down <- results_list[[results_name_down]]
      
      if (is.null(ego_up) || is.null(ego_down)) {
        next
      }
      
      up_data <- head(ego_up@result[, c("Description", "p.adjust")], 20)
      down_data <- head(ego_down@result[, c("Description", "p.adjust")], 20)
      
      up_data$neglog_pval <- -log10(up_data$"p.adjust")
      down_data$neglog_pval <- -log10(down_data$"p.adjust")
      
      all_data <- rbind(
        data.frame(Category = up_data$Description, Value = up_data$neglog_pval, Direction = "Up"),
        data.frame(Category = down_data$Description, Value = -down_data$neglog_pval, Direction = "Down")
      )
      
      plot_title <- paste(celltype, condition_1, "vs", condition_2)
      file_name <- paste0(output_path, gsub("/", "-", plot_title), ".pdf")
      
      pdf(file_name, width = 300/25.4, height = 150/25.4)
      
      
      # Adjust the bottom margin
      old_par <- par(mar = c(10, 4, 10, 2) + 0.1)
      
      all_data$Category <- sapply(all_data$Category, break_long_text, max_len = 30)
      
      ylim_adjusted <- c(min(all_data$Value) - 1.5, max(all_data$Value))
      midpoints <- barplot(all_data$Value, 
                           col = c(rep("red", nrow(up_data)), rep("blue", nrow(down_data))), 
                           ylim = ylim_adjusted,
                           cex.names = 0.6, 
                           border = "white",
                           axes = FALSE,
                           ylab = "-log10(p-value)")
      
      axis(2)
      
      # Staggered text labels with corresponding dashed lines
      for (idx in 1:length(midpoints)) {
        if (idx %% 2 == 1) {  # odd index
          text(midpoints[idx], ylim_adjusted[1] - 0.2, srt = 45, adj = 1, labels = all_data$Category[idx], xpd = TRUE, cex=0.6)
          segments(x0 = midpoints[idx], y0 = ylim_adjusted[1], x1 = midpoints[idx], y1 = 0, lty = 2)
        } else {  # even index
          text(midpoints[idx], ylim_adjusted[2] + 0.2, srt = 45, adj = 0, labels = all_data$Category[idx], xpd = TRUE, cex=0.6) # adj is set to 0 for left alignment
          segments(x0 = midpoints[idx], y0 = ylim_adjusted[2], x1 = midpoints[idx], y1 = 0, lty = 2)
        }
      }
      
      # Add the title at the bottom
      title(main = plot_title, line = -2, outer = TRUE)
      
      par(old_par)
      dev.off()
    }
  }
}
