rm(list = ls(all = TRUE))
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

path_string <- "/home/qianyu/Desktop/grid/snRNAseq_all_nuclei_V1/"

V1 <- readRDS(paste0(path_string, "output/Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean2_corrected.rds"))
V1 <- subset(V1, subset = spec_Celltype != "Endo")

p1 <- DimPlot(V1, reduction = "UMAP",label = FALSE, raster=FALSE, group.by="spec_Celltype")
p2 <- DimPlot(V1, reduction = "UMAP",label = FALSE, raster=FALSE, group.by="condition_simpl")
p3 <- DimPlot(V1, reduction = "UMAP",label = FALSE, raster=FALSE, group.by="batch")
p4 <- DimPlot(V1, reduction = "UMAP",label = TRUE, raster=FALSE, group.by="gen_CellType")

######################################################################################################
############################################### Fig 1A ###############################################
######################################################################################################
# save as pdf 10 inch squre for images
p1 + NoLegend() + NoAxes() + labs(title = "")

# save as pdf 5 inch squre only for legend
p1 + labs(title = "") + 
  theme(text = element_text(family = "Arial", size = 5),
        plot.title = element_text(size = 6),
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5))

tmp <- p1 + labs(title = "") + 
  theme(text = element_text(family = "Arial", size = 5),
        plot.title = element_text(size = 6),
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5))

ggsave("/home/qianyu/Desktop/a.svg", plot = tmp, width = 5, height = 3)
######################################################################################################
############################################### Fig 1B ###############################################
######################################################################################################
# save as pdf 10 inch squre for images
p2 + NoLegend() + NoAxes() + labs(title = "")

# save as pdf 5 inch squre only for legend
p2 + labs(title = "") + 
  theme(text = element_text(family = "Arial", size = 5),
        plot.title = element_text(size = 6),
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5))

tmp <- p2 + labs(title = "") + 
  theme(text = element_text(family = "Arial", size = 5),
        plot.title = element_text(size = 6),
        legend.text = element_text(size = 5), 
        legend.title = element_text(size = 5))

ggsave("C:\\Users\\lucas\\Desktop\\a.svg", plot = tmp, width = 5, height = 3)

######################################################################################################
############################################### Fig 1C ###############################################
######################################################################################################
Idents(V1) <- "spec_Celltype"
new_levels <- c("Exc L2/3", "Exc L4", "Exc L5IT", "Exc L5NP", "Exc L5PT", "Exc L6a", 
                "Exc L6b", "Exc L6IT", "Inh Grin3a", "Inh Npy", "Inh Pval", "Inh Vip", 
                "Ast", "Mg", "Olg", "OPCs")
V1@meta.data$spec_Celltype <- factor(V1@meta.data$spec_Celltype, levels = new_levels)
tmp <- V1
# tmp <- subset(V1, subset = gen_CellType == "Glutamatergic Neurons") # "Glutamatergic Neurons"/"GABAergic Neurons"/"Glia cells"/"Glutamatergic Neurons"/"Endo"    
tmp <- subset(tmp, cells = sample(colnames(tmp), size = round(length(colnames(tmp)) * 0.1)))
markers <- FindAllMarkers(tmp, min.pct = 0.25, logfc.threshold = 0.58)
top <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(tmp, features = top$gene, slot = "data") + NoLegend()

DEGs <- c(
  "Gpc6", "Lrrtm4", "Dgkb", # Exc L2/3
  "Rorb", "Zmat4", "Brinp3", # Exc L4
  "Il1rapl2", "Ptprm", "Galntl6", # Exc L5IT
  "Tshz2", "Nxph1", "Vwc2l", # Exc L5NP
  "Hcn1", "Robo1", "Sorcs1", # Exc L5PT
  "Zfpm2", "Foxp2", "Cdh18", # Exc L6a
  "Inpp4b", "Kcnab1", "Fbxl7", # Exc L6b
  "Zfp804b", "Galnt14", "Zfp804a", # Exc L6IT
  "Grin3a", # Inh Grin3a
  "Npy", "Pde11a", "Cacna2d1", # Inh Npy
  "Pvalb", "Dock4", "6330411D24Rik", # Inh Pval
  "Vip", "Sorcs3", "Adarb2", # Inh Vip
  "Slc1a2", "Ntm", # Ast
  "Inpp5d", "Zfhx3", "Tgfbr1", # Mg
  "Nkain2", "Mbp", "Slc24a2", # Olg
  "Lhfpl3", "Nxph1", "Pdgfra", # OPCs
  "Flt1" # Endo
)

genes <- c(
  "Lrrtm4", # Exc + Inh
  "Sv2b", # Exc
  "Gpc6", # Exc L2/3
  "Rorb", # Exc L4
  "Il1rapl2", # Exc L5IT
  "Tshz2", # Exc L5NP
  "Galnt14", # Exc L5PT
  "Foxp2", # Exc L6a
  "Inpp4b", # Exc L6b
  "Zfp804b", # Exc L6IT
  "Kcnmb2", # Inh
  "Grin3a", # Inh Grin3a
  "Pde11a", # Inh Npy
  "Cemip", # Inh Pvalb
  "Vip", # Inh Vip
  "Atp1a2", # Ast
  "Tgfbr1", # Mg
  "Mbp", # Olg
  "Pdgfra" # OPCs
)

ff <-
VlnPlot(
  V1,
  features = genes,
  stack = TRUE,
  pt.size=0,
  flip = TRUE
) + NoLegend()

ggsave(paste0("/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home/qianyu/snRNAseq_all_nuclei_V1/paper/Figure1.markers.svg"), plot = p, width = 10, height = 10)

######################################################################################################
############################################### Fig 1D ###############################################
######################################################################################################

library(scales)
vv <- V1
pt <- table(vv$spec_Celltype, vv$condition_simpl)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

p1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_classic() +
  #geom_col(position = "fill", width = 0.8) +
  geom_col(width = 0.8) +
  xlab("V1") +
  ylab("Count") +
  theme(legend.position="none") +
  theme(legend.title = element_blank())+
  theme(legend.position = "none",
        legend.title = element_blank(),
        text = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 5))

ggsave(paste0("/home/qianyu/Desktop/Allbatch.pdf"), plot = p1, width = 1.5, height = 1.5, dpi = 300)





batch_choose = "batch3" # 1 2 3
vv <- subset(V1, subset = batch == batch_choose)
pt <- table(vv$spec_Celltype, vv$condition_simpl)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

p1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = label_number_si(),
                     limits = c(0, 8000))+
  theme_classic() +
  #geom_col(position = "fill", width = 0.8) +
  geom_col(width = 0.8) +
  xlab(batch_choose) +
  ylab("Count") +
  theme(legend.position="none") +
  theme(legend.title = element_blank())+
  theme(legend.position = "none",
        legend.title = element_blank(),
        text = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 5))

ggsave(paste0("C:\\Users\\lucas\\Desktop\\",batch_choose,".svg"), plot = p1, width = 1.5, height = 1.5, dpi = 300)


# batch 4 5 6
batch_choose <- c("batch4", "batch5", "batch6")
vv <- subset(V1, subset = batch %in% batch_choose)
vv$batch <- factor(vv$batch, levels = batch_choose)
pt <- table(vv$spec_Celltype, vv$batch)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var2 <- paste("NR (", pt$Var2, ")", sep = "")

p1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  scale_y_continuous(expand = c(0, 0),
                     labels = label_number_si(),
                     limits = c(0, 8000))+
  theme_classic() +
  #geom_col(position = "fill", width = 0.8) +
  geom_col(width = 0.8) +
  xlab("batch 4/5/6") +
  ylab("Count") +
  theme(legend.position="none") +
  theme(legend.title = element_blank())+
  theme(legend.position = "none",
        legend.title = element_blank(),
        text = element_text(size = 5),
        axis.title = element_text(size = 5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 5))

ggsave(paste0("C:\\Users\\lucas\\Desktop\\batch456.svg"), plot = p1, width = 1, height = 1.5, dpi = 300)

