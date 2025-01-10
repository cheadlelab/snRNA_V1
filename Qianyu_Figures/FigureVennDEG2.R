rm(list = ls(all = TRUE))
#Run in the R 4.1.2 version
library(ggplot2)
library(monocle3)
library(dplyr)
library(tidyverse)
library(edgeR)
library(venn)
library(DESeq2)

path <- "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/qianyu/snRNAseq_all_nuclei_V1/output/step6_DEGs_edgeR"
save_path = "/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/qianyu/snRNAseq_all_nuclei_V1/figures/Venn"

############################################################################################################################################
################################################################ L4 vs L2/3 ################################################################
############################################################################################################################################
#  (1) L4 vs L2/3;

timepoints <- c("P27", "LDR30m", "LDR2h", "LDR4h", "LDR6h")
for (i in timepoints) {
  tmp1 <- read.csv(paste0(path,"/Exc L4_LDR0.vs.",i,".csv"))
  tmp1_DEG <- tmp1[abs(tmp1$logFC) > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_UPR <- tmp1[tmp1$logFC > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_DWR <- tmp1[tmp1$logFC < -log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  
  tmp2 <- read.csv(paste0(path,"/Exc L2-3_LDR0.vs.",i,".csv"))
  tmp2_DEG <- tmp2[abs(tmp2$logFC) > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_UPR <- tmp2[tmp2$logFC > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_DWR <- tmp2[tmp2$logFC < -log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  
  deg.use.gene <- list(Exc_L4 = tmp1_DEG, Exc_L2_3 = tmp2_DEG)
  upr.use.gene <- list(Exc_L4 = tmp1_UPR, Exc_L2_3 = tmp2_UPR)
  dwr.use.gene <- list(Exc_L4 = tmp1_DWR, Exc_L2_3 = tmp2_DWR)
  
  p <- venn(deg.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".L234.DEG.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(deg.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(upr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".L234.UPR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(upr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(dwr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".L234.DWR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(dwr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
}


############################################################################################################################################
#################################################################### L5 ####################################################################
############################################################################################################################################
#  (2) All 3 L5 populations;
timepoints <- c("P27", "LDR30m", "LDR2h", "LDR4h", "LDR6h")
for (i in timepoints) {
  tmp1 <- read.csv(paste0(path,"/Exc L5IT_LDR0.vs.",i,".csv"))
  tmp1_DEG <- tmp1[abs(tmp1$logFC) > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_UPR <- tmp1[tmp1$logFC > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_DWR <- tmp1[tmp1$logFC < -log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  
  tmp2 <- read.csv(paste0(path,"/Exc L5NP_LDR0.vs.",i,".csv"))
  tmp2_DEG <- tmp2[abs(tmp2$logFC) > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_UPR <- tmp2[tmp2$logFC > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_DWR <- tmp2[tmp2$logFC < -log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  
  tmp3 <- read.csv(paste0(path,"/Exc L5PT_LDR0.vs.",i,".csv"))
  tmp3_DEG <- tmp3[abs(tmp3$logFC) > log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  tmp3_UPR <- tmp3[tmp3$logFC > log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  tmp3_DWR <- tmp3[tmp3$logFC < -log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  
  deg.use.gene <- list(Exc_L5IT = tmp1_DEG, Exc_L5NP = tmp2_DEG, Exc_L5PT = tmp3_DEG)
  upr.use.gene <- list(Exc_L5IT = tmp1_UPR, Exc_L5NP = tmp2_UPR, Exc_L5PT = tmp3_UPR)
  dwr.use.gene <- list(Exc_L5IT = tmp1_DWR, Exc_L5NP = tmp2_DWR, Exc_L5PT = tmp3_DWR)
  
  p <- venn(deg.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".L5.DEG.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(deg.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(upr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".L5.UPR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(upr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(dwr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".L5.DWR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(dwr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
}

#############################################################################################################################################
#################################################################### Inh ####################################################################
#############################################################################################################################################
#  (3) all 4 inhibitory populations
timepoints <- c("P27", "LDR30m", "LDR2h", "LDR4h", "LDR6h")
for (i in timepoints) {
  tmp1 <- read.csv(paste0(path,"/Inh Pvalb_LDR0.vs.",i,".csv"))
  tmp1_DEG <- tmp1[abs(tmp1$logFC) > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_UPR <- tmp1[tmp1$logFC > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_DWR <- tmp1[tmp1$logFC < -log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  
  tmp2 <- read.csv(paste0(path,"/Inh Grin3a_LDR0.vs.",i,".csv"))
  tmp2_DEG <- tmp2[abs(tmp2$logFC) > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_UPR <- tmp2[tmp2$logFC > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_DWR <- tmp2[tmp2$logFC < -log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  
  tmp3 <- read.csv(paste0(path,"/Inh Vip_LDR0.vs.",i,".csv"))
  tmp3_DEG <- tmp3[abs(tmp3$logFC) > log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  tmp3_UPR <- tmp3[tmp3$logFC > log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  tmp3_DWR <- tmp3[tmp3$logFC < -log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  
  tmp4 <- read.csv(paste0(path,"/Inh Npy_LDR0.vs.",i,".csv"))
  tmp4_DEG <- tmp4[abs(tmp4$logFC) > log2(1.5) & tmp4$P.Value < 0.05, ]$gene
  tmp4_UPR <- tmp4[tmp4$logFC > log2(1.5) & tmp4$P.Value < 0.05, ]$gene
  tmp4_DWR <- tmp4[tmp4$logFC < -log2(1.5) & tmp4$P.Value < 0.05, ]$gene
  
  deg.use.gene <- list(Inh_Pvalb = tmp1_DEG, Inh_Grin3a = tmp2_DEG, Inh_Vip = tmp3_DEG, Inh_Npy = tmp4_DEG)
  upr.use.gene <- list(Inh_Pvalb = tmp1_UPR, Inh_Grin3a = tmp2_UPR, Inh_Vip = tmp3_UPR, Inh_Npy = tmp4_UPR)
  dwr.use.gene <- list(Inh_Pvalb = tmp1_DWR, Inh_Grin3a = tmp2_DWR, Inh_Vip = tmp3_DWR, Inh_Npy = tmp4_DWR)
  
  p <- venn(deg.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".Inh.DEG.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(deg.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(upr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".Inh.UPR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(upr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(dwr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".Inh.DWR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(dwr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
}

###############################################################################################################################################
#################################################################### glial ####################################################################
###############################################################################################################################################
#  (4) all glial populations
timepoints <- c("P27", "LDR30m", "LDR2h", "LDR4h", "LDR6h")
for (i in timepoints) {
  tmp1 <- read.csv(paste0(path,"/Ast_LDR0.vs.",i,".csv"))
  tmp1_DEG <- tmp1[abs(tmp1$logFC) > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_UPR <- tmp1[tmp1$logFC > log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  tmp1_DWR <- tmp1[tmp1$logFC < -log2(1.5) & tmp1$P.Value < 0.05, ]$gene
  
  tmp2 <- read.csv(paste0(path,"/Olg_LDR0.vs.",i,".csv"))
  tmp2_DEG <- tmp2[abs(tmp2$logFC) > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_UPR <- tmp2[tmp2$logFC > log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  tmp2_DWR <- tmp2[tmp2$logFC < -log2(1.5) & tmp2$P.Value < 0.05, ]$gene
  
  tmp3 <- read.csv(paste0(path,"/Mg_LDR0.vs.",i,".csv"))
  tmp3_DEG <- tmp3[abs(tmp3$logFC) > log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  tmp3_UPR <- tmp3[tmp3$logFC > log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  tmp3_DWR <- tmp3[tmp3$logFC < -log2(1.5) & tmp3$P.Value < 0.05, ]$gene
  
  opcs_file_path <- paste0(path, "/OPCs_LDR0.vs.", i, ".csv")
  if (file.exists(opcs_file_path)) {
    tmp4 <- read.csv(opcs_file_path)
    tmp4_DEG <- tmp3[abs(tmp4$logFC) > log2(1.5) & tmp4$P.Value < 0.05, ]$gene
    tmp4_UPR <- tmp3[tmp4$logFC > log2(1.5) & tmp4$P.Value < 0.05, ]$gene
    tmp4_DWR <- tmp3[tmp4$logFC < -log2(1.5) & tmp4$P.Value < 0.05, ]$gene
    # Existing logic to filter DEGs
  } else {
    tmp4_DEG <- character(0)  # No data, set to empty character vector
    tmp4_UPR <- character(0)
    tmp4_DWR <- character(0)
  }

  deg.use.gene <- list(Ast = tmp1_DEG, Olg = tmp2_DEG, Mg = tmp3_DEG, OPCs = tmp4_DEG)
  upr.use.gene <- list(Ast = tmp1_UPR, Olg = tmp2_UPR, Mg = tmp3_UPR, OPCs = tmp4_UPR)
  dwr.use.gene <- list(Ast = tmp1_DWR, Olg = tmp2_DWR, Mg = tmp3_DWR, OPCs = tmp4_DWR)
  
  p <- venn(deg.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".Glial.DEG.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(deg.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(upr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".Glial.UPR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(upr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
  
  p <- venn(dwr.use.gene,  zcolor = "style")
  pdf(file = paste0(save_path, "/", i, ".Glial.DWR.pdf"), width = 8, height = 8) # Adjust width and height as needed
  venn(dwr.use.gene,  zcolor = "style", counts =p$counts, box = FALSE)
  dev.off() # Close the PDF device
}