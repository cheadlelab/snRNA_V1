library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(cowplot)
library(purrr)
library(tibble)
library(viridis)

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction")

Allexp.sct_nobatch_corrected.40.05_norm <- readRDS("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction/Allexp.sct_nobatch_corrected.40.05_norm.rds")

MAST_res05_top20_filtered_conserved_markers_final <- read.csv("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction/MAST_res05_top20_filtered_conserved_markers_final.csv")

MAST_res05_top20_conserved_markers_final <- read.csv("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction/MAST_res05_top20_conserved_markers_final.csv")

MAST_res05_top3_genes <- read.csv("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction/MAST_res05_top3_genes.csv")

DefaultAssay(Allexp.sct_nobatch_corrected.40.05_norm) <- "decontXcountsclean"

#Annotating cluster cells
library(rcellmarker)

df1<- MAST_res05_top20_filtered_conserved_markers_final[, c(3,2)]
df2<- MAST_res05_top20_conserved_markers_final[, c(3,2)]
colnames(df1) <- c("gene", "Cluster")
colnames(df2) <- c("gene", "Cluster")

res1 <- cellMarker(df1,type='custom',species='mouse',keytype='SYMBOL')
res2 <- cellMarker(df1,type='custom',species='mouse', minSize = 1, keytype='SYMBOL')
res3 <- cellMarker(df2,type='custom',species='mouse', minSize = 1, keytype='SYMBOL')

write.csv(res1, "MAST_cluster_annotation_top20_filtered.csv")
write.csv(res2, "MAST_cluster_annotation_top20_filtered_minSize1.csv")
write.csv(res3, "MAST_cluster_annotation_top20.csv")

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction/feature_plots")

# set up list of canonical cell type markers
canonical_markers <- list(
  'Astrocyte' = c('Atp1a2', 'Slc1a3', 'Aqp4'),
  'Pan-neuronal' = c('Snap25', 'Syt1'),
  'Excitatory Neuron' = c('Slc17a7', 'Satb2'),
  'Inhibitory Neuron' = c('Gad1', 'Gad2', 'Pvalb', 'Vip'),
  'Microglia' = c('Csf1r', 'Cx3cr1', 'Tgfbr1'),
  'Oligodendrocyte' = c('Mobp', 'Mbp', 'Mog'),
  'OPCs' = c('Pdgfra', 'Cspg4'),
  'Endothelial cells' = c('Flt1', 'Id1', 'Xdh', 'Tbc1d4', 'Exosc7', 'Pecam1', 'Cdh5', 'Tie1'),
  'VLMC' = c('Ogn', 'Lum', 'Dcn', 'Bgn', 'Myl9', 'Crip1', 'Vtn')
)

# create feature plots, cutoff expression values for the 98th and 99th percentile
#panel.background=element_rect(fill = 'white', color = 'black')
umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)

plot_list <- FeaturePlot(
  Allexp.sct_nobatch_corrected.40.05_norm,
  features=unlist(canonical_markers),
  combine=FALSE, cols=viridis(256),
  max.cutoff='q98'
)

# apply theme to each feature plot
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme + NoLegend()
}
combined_plot <- CombinePlots(plot_list)
ggsave(combined_plot, filename = "basic_canonical_marker_featurePlot.png", height = 10, width = 10)

#plot Vlnplot
for(celltype in names(general_celltype_annotations)){
  
  print(celltype)
  cur_features <- general_celltype_annotations[[celltype]]
  
  # plot distributions for marker genes:
  p <- VlnPlot(
    Allexp.sct_nobatch_corrected.40.05_norm,
    group.by='seurat_clusters',
    features=cur_features,
    pt.size = 0, ncol=1
  )
  png(paste0('basic_canonical_marker_',celltype,'_vlnPlot.png'), width=10, height=3*length(cur_features), units='in', res=200, type = "cairo")
  print(p)
  dev.off()
  
}

#Module score
cluster_enriched_list <- list()
unique_clusters <- unique(MAST_res05_top20_conserved_markers_final$cluster_id)
for (i in unique_clusters) {
  cluster_enriched <- MAST_res05_top20_conserved_markers_final %>%
    filter(cluster_id == i) %>%
    arrange(-NR_avg_log2FC) %>%
    pull(gene) %>%
    .[1:20]
  # Store the top 20 genes in a named list
  cluster_enriched_list[[paste0("Cluster", i, "_top20_genes")]] <- cluster_enriched
}

# Remove NA values from the list
cluster_enriched_list <- lapply(cluster_enriched_list, na.omit)

Allexp.sct_nobatch_corrected.40.05_norm_MODULE <- AddModuleScore(Allexp.sct_nobatch_corrected.40.05_norm,
                                                        features = cluster_enriched_list,
                                                        name=c("cluster_enriched0", "cluster_enriched1", "cluster_enriched2", "cluster_enriched3", "cluster_enriched4", "cluster_enriched5", "cluster_enriched6", "cluster_enriched7", "cluster_enriched8", "cluster_enriched9",
                                                        "cluster_enriched10", "cluster_enriched11", "cluster_enriched12", "cluster_enriched13", "cluster_enriched14", "cluster_enriched15", "cluster_enriched17", "cluster_enriched18", "cluster_enriched19", "cluster_enriched20"
                                                        , "cluster_enriched21", "cluster_enriched22", "cluster_enriched23", "cluster_enriched24", "cluster_enriched25", "cluster_enriched26", "cluster_enriched28", "cluster_enriched29"
                                                        , "cluster_enriched30", "cluster_enriched31", "cluster_enriched32", "cluster_enriched33", "cluster_enriched34"))

#Changing the metadata column names
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[20] <- "Cluster_enriched0"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[21] <- "Cluster_enriched1"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[22] <- "Cluster_enriched2"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[23] <- "Cluster_enriched3"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[24] <- "Cluster_enriched4"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[25] <- "Cluster_enriched5"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[26] <- "Cluster_enriched6"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[27] <- "Cluster_enriched7"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[28] <- "Cluster_enriched8"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[29] <- "Cluster_enriched9"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[30] <- "Cluster_enriched10"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[31] <- "Cluster_enriched11"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[32] <- "Cluster_enriched12"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[33] <- "Cluster_enriched13"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[34] <- "Cluster_enriched14"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[35] <- "Cluster_enriched15"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[36] <- "Cluster_enriched17"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[37] <- "Cluster_enriched18"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[38] <- "Cluster_enriched19"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[39] <- "Cluster_enriched20"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[40] <- "Cluster_enriched21"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[41] <- "Cluster_enriched22"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[42] <- "Cluster_enriched23"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[43] <- "Cluster_enriched24"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[44] <- "Cluster_enriched25"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[45] <- "Cluster_enriched26"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[46] <- "Cluster_enriched28"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[47] <- "Cluster_enriched29"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[48] <- "Cluster_enriched30"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[49] <- "Cluster_enriched31"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[50] <- "Cluster_enriched32"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[51] <- "Cluster_enriched33"
colnames(Allexp.sct_nobatch_corrected.40.05_norm_MODULE@meta.data)[52] <- "Cluster_enriched34"

list_names <- list(c("Cluster_enriched0", "Cluster_enriched1", "Cluster_enriched2", "Cluster_enriched3", "Cluster_enriched4", "Cluster_enriched5", "Cluster_enriched6", "Cluster_enriched7", "Cluster_enriched8", "Cluster_enriched9",
                     "Cluster_enriched10", "Cluster_enriched11", "Cluster_enriched12", "Cluster_enriched13", "Cluster_enriched14", "Cluster_enriched15", "Cluster_enriched17", "Cluster_enriched18", "Cluster_enriched19", "Cluster_enriched20"
                     , "Cluster_enriched21", "Cluster_enriched22", "Cluster_enriched23", "Cluster_enriched24", "Cluster_enriched25", "Cluster_enriched26", "Cluster_enriched28", "Cluster_enriched29"
                     , "Cluster_enriched30", "Cluster_enriched31", "Cluster_enriched32", "Cluster_enriched33", "Cluster_enriched34"))

# Plot scores
Cluster_enriched_Plot_list <- FeaturePlot(Allexp.sct_nobatch_corrected.40.05_norm_MODULE,
                                          features = unlist(list_names), label = TRUE, repel = TRUE, cols=viridis(256)) 

ggsave(Cluster_enriched_Plot_list,
       filename = paste0("Cluster", "_enriched.png"),
       height = 30,
       width = 30, dpi = 600)

saveRDS(Allexp.sct_nobatch_corrected.40.05_norm_MODULE, "Allexp.sct_nobatch_corrected.40.05_norm_MODULE.rds")

# plot heatmap:
MAST_res05_top3_genes <- MAST_res05_top3_genes$x
DefaultAssay(Allexp.sct_nobatch_corrected.40.05_norm) <- "sctdecontXcountsclean"
heatmap_clusters <- DoHeatmap(Allexp.sct_nobatch_corrected.40.05_norm, group.by ="seurat_clusters", label=FALSE, features= MAST_res05_top3_genes) + NoLegend()
ggsave(heatmap_clusters, filename = "basic_canonical_marker_heatmap.png", height = 10, width = 20, dpi = 300)

#DotPlot with 1 representative gene per cluster
DefaultAssay(Allexp.sct_nobatch_corrected.40.05_norm) <- "decontXcountsclean"
markers.to.plot <- c("Nectin3", "Nell1", "Foxp2", "Atp1a2", "Sema3e", "Mobp", "Deptor", "Thsd7a", "Kcnmb2", "Gm2164", "Tgfbr1", "Grin3a", "Tshz2", "Pdgfra", "Adarb2", "Gad2", "Nxph1", "Bcas1", "Zfp536", "Thsd7b", "Tmem163", "Slc47a1", "Tead1", "Cped1", "Zfhx3", "Rxfp1", "Slc1a3", "Ebf1", "Mertk")
Cluster_genes <- DotPlot(Allexp.sct_nobatch_corrected.40.05_norm, features = markers.to.plot, cols = c("black", "red", "orange", "yellow", "green", "blue"), dot.scale = 10, split.by = "condition_simpl") #+ RotatedAxis()
ggsave(Cluster_genes, filename = "Cluster_genes_dotplot.png", height = 20, width = 30)

#Highlighting only one cluster of cells (In case you are not sure about the distribution of a cluster)
cluster22 <- WhichCells(Allexp.sct_nobatch_corrected.40.05_norm, idents = "22")
cluster22_plot <- DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, label=T, cells.highlight= cluster22, cols.highlight = c("darkblue", "darkred"), cols= "grey")
ggsave(cluster22_plot, filename = "cluster22_plot.png", height = 7, width = 10, type = "cairo")

cluster30 <- WhichCells(Allexp.sct_nobatch_corrected.40.05_norm, idents = "30")
cluster30_plot <- DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, label=T, cells.highlight= cluster30, cols.highlight = c("darkblue", "darkred"), cols= "grey")
ggsave(cluster30_plot, filename = "cluster30_plot.png", height = 7, width = 10, type = "cairo")

cluster31 <- WhichCells(Allexp.sct_nobatch_corrected.40.05_norm, idents = "31")
cluster31_plot <- DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, label=T, cells.highlight= cluster31, cols.highlight = c("darkblue", "darkred"), cols= "grey")
ggsave(cluster31_plot, filename = "cluster31_plot.png", height = 7, width = 10, type = "cairo")

cluster32 <- WhichCells(Allexp.sct_nobatch_corrected.40.05_norm, idents = "32")
cluster32_plot <- DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, label=T, cells.highlight= cluster32, cols.highlight = c("darkblue", "darkred"), cols= "grey")
ggsave(cluster32_plot, filename = "cluster32_plot.png", height = 7, width = 10, type = "cairo")

# set up list of general cell type markers
general_celltype_annotations <- list(
  '0' = 'Glutamatergic Neurons',
  '1' = 'Glutamatergic Neurons',
  '2' = 'Glutamatergic Neurons',
  '3' = 'Glia cells',
  '4' = 'Glutamatergic Neurons',
  '5' = 'Glutamatergic Neurons',
  '6' = 'Glutamatergic Neurons', 
  '7' = 'Glia cells',
  '8' = 'Glutamatergic Neurons',
  '9' = 'Glutamatergic Neurons',
  '10' = 'GABAergic Neurons',
  '11' = 'Glutamatergic Neurons',
  '12' = 'Glia cells',
  '13' = 'GABAergic Neurons',
  '14' = 'Glutamatergic Neurons',
  '15' = 'Glia cells',
  '16' = 'Glutamatergic Neurons',
  '17' = 'Mito-genes',
  '18' = 'GABAergic Neurons',
  '19' = 'GABAergic Neurons',
  '20' = 'Glutamatergic Neurons',
  '21' = 'Glia cells',
  '22' = 'Glutamatergic Neurons',
  '23' = 'Glutamatergic Neurons',
  '24' = 'Glutamatergic Neurons',
  '25' = 'Unk',
  '26' = 'Glutamatergic Neurons',
  '27' = 'Glutamatergic Neurons',
  '28' = 'Unk',
  '29' = 'Unk',
  '30' = 'Glutamatergic Neurons',
  '31' = 'Glutamatergic Neurons',
  '32' = 'Glutamatergic Neurons',
  '33' = 'Endo',
  '34' = 'Glia cells'
)

# add General CellType to seurat metadata
Allexp.sct_nobatch_corrected.40.05_norm@meta.data$gen_CellType <- unlist(general_celltype_annotations[Allexp.sct_nobatch_corrected.40.05_norm@meta.data$seurat_clusters])

png('umap_gen_celltypes.png', width=8, height=7, res=200, units='in', type = "cairo")
DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, reduction = "UMAP", group.by='gen_CellType') +
  umap_theme + ggtitle('UMAP colored by general cell type annotations')
dev.off()

#Finding specific genes by celltype

spec_markers <- list(
  'Sub_Inh' = c('Sst', 'Ndnfm', 'Vip', 'Sncg', 'Pvalb', 'Igtp', 'Kcnmb2', 'Grin3a', 'Kcnip1', 'Dlx6os1'),
  'Sub_exc' = c('Calb2', 'Enpp2', 'Otof', 'Ngb', 'Rorb','Cux2', 'Ccbe1', 'Mdga1', 'Stard8', 'Whrn', 'Rorb', 'Deptor',
                'Foxo1', 'Ptprm', 'Nxph1', 'Tshz2', 'Trhr', 'Slc17a8', 
                'Bcl6', 'Cdh13', 'Erg', 'Reln', 'Foxp2', 'Syt6', 'Zfp804b', 'Cdh9', 'Ctgf', 'Inpp4b', 'Svil', 'Tnmd', 'Serpinb11'))

# create feature plots, cutoff expression values for the 98th and 99th percentile
#panel.background=element_rect(fill = 'white', color = 'black')
umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)

plot_list <- FeaturePlot(
  Allexp.sct_nobatch_corrected.40.05_norm,
  features=unlist(spec_markers),
  combine=FALSE, cols=viridis(256),
  max.cutoff='q98'
)

# apply theme to each feature plot
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme + NoLegend()
}
combined_plot <- CombinePlots(plot_list)
ggsave(combined_plot, filename = "spec_marker_featurePlot.png", height = 10, width = 10)

# set up list of specific markers
specific_markers_annotations <- list(
  '0' = 'Exc L2/3',
  '1' = 'Exc L4',
  '2' = 'Exc L6a',
  '3' = 'Ast',
  '4' = 'Exc L6IT',
  '5' = 'Exc L4',
  '6' = 'Exc L2/3',
  '7' = 'Olg',
  '8' = 'Exc L5IT',
  '9' = 'Exc L4',
  '10' = 'Inh Pvalb',
  '11' = 'Exc L5PT',
  '12' = 'Mg',
  '13' = 'Inh Grin3a',
  '14' = 'Exc L5NP',
  '15' = 'OPCs',
  '16' = 'Exc L2/3',
  '17' = 'Mito-genes',
  '18' = 'Inh Vip',
  '19' = 'Inh Npy',
  '20' = 'Exc L5IT',
  '21' = 'Olg',
  '22' = 'Unk', #Badly distributed
  '23' = 'Exc L6a',
  '24' = 'Exc L6b',
  '25' = 'Dbl', #More than one marker (Ast, Neurons, OPCs)
  '26' = 'Exc L5NP',
  '27' = 'Exc L6a',
  '28' = 'Dbl', #More than one marker (Ast, Neurons, OPCs)
  '29' = 'Dbl', #More than one marker (Ast, Neurons, OPCs)
  '30' = 'Unk', #Badly distributed
  '31' = 'Exc L6IT',
  '32' = 'Dbl', #More than one marker (Ast, Neurons, OPCs)
  '33' = 'Endo',
  '34' = 'Ast'
)

Allexp.sct_nobatch_corrected.40.05_norm@meta.data$spec_Celltype <- unlist(specific_markers_annotations[Allexp.sct_nobatch_corrected.40.05_norm@meta.data$seurat_clusters])
Allexp.sct_nobatch_corrected.40.05_norm@meta.data$spec_Celltype_cluster <- paste0(Allexp.sct_nobatch_corrected.40.05_norm@meta.data$spec_Celltype, '-', Allexp.sct_nobatch_corrected.40.05_norm@meta.data$seurat_clusters)

png('umap_spec_celltype_clusters.png', width=8, height=8, res=200, units='in', type = "cairo")
DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, reduction = "UMAP", group.by='spec_Celltype_cluster', label=TRUE) +
  umap_theme + ggtitle('UMAP colored by spec cell type + cluster') + NoLegend()
dev.off()

png('umap_spec_celltype.png', width=8, height=8, res=200, units='in', type = "cairo")
DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, reduction = "UMAP", group.by='spec_Celltype', label=TRUE) +
  umap_theme + ggtitle('UMAP colored by spec cell type') + NoLegend()
dev.off()

Allexp.sct_nobatch_corrected.40.05_norm_labelled <- Allexp.sct_nobatch_corrected.40.05_norm

#subsetting the clusters for DEG analysis (remove the unknown, doublets, misc and Mito-genes clusters)
Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean <- subset(Allexp.sct_nobatch_corrected.40.05_norm_labelled, idents = c(17, 22, 25, 28, 29, 30, 32), invert = TRUE)

png('umap_gen_celltypes_clean.png', width=8, height=7, res=200, units='in', type = "cairo")
DimPlot(Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean, reduction = "UMAP", group.by='gen_CellType') +
  umap_theme + ggtitle('UMAP colored by general cell type annotations')
dev.off()

png('umap_spec_celltype_clusters_clean.png', width=8, height=8, res=200, units='in', type = "cairo")
DimPlot(Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean, reduction = "UMAP", group.by='spec_Celltype_cluster', label=TRUE) +
  umap_theme + ggtitle('UMAP colored by spec cell type + cluster') + NoLegend()
dev.off()

png('umap_spec_celltype_clean.png', width=8, height=8, res=200, units='in', type = "cairo")
DimPlot(Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean, reduction = "UMAP", group.by='spec_Celltype', label=TRUE, repel = TRUE) +
  umap_theme + ggtitle('UMAP colored by spec cell type') + NoLegend()
dev.off()

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/NGS/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction")

saveRDS(Allexp.sct_nobatch_corrected.40.05_norm_labelled, "Allexp.sct_nobatch_corrected.40.05_norm_labelled.rds")
saveRDS(Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean, "Allexp.sct_nobatch_corrected.40.05_norm_labelled_clean.rds")

#UMAPAllexp.sct_nobatch_corrected.40.05_norm_labelled <- DimPlot(Allexp.sct_nobatch_corrected.40.05_norm, reduction = "UMAP", label = F, repel = F)
#ggsave(UMAPAllexp.sct_nobatch_corrected.40.05_norm_labelled, filename = "UMAPAllexp.sct_nobatch_corrected.labelled.png", height = 7, width = 10, type = "cairo")

#Violin plot per cluster with top 10
#Cluster0_VlnPlot <- VlnPlot(Allexp.sct_nobatch_corrected.40.01_norm, features = c("Nell1", "Gm15398", "Kcnh5", "Zmat4", "Prr16"), slot = "data", log = F)
#ggsave(Cluster0_VlnPlot, filename = "Cluster0_VlnPlot.png", height = 7, width = 10, type = "cairo")

