library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(dplyr)
library(cowplot)
library(purrr)
library(tibble)
library(EnsDb.Mmusculus.v79)
library(AnnotationHub)
library(stringr)

setwd("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/03_Clusters_labeling/less_restrictive/no_batch_correction")

#Allexp.sct_nobatch_corrected.40.05 <- readRDS("/grid/cheadle/home/machado/Projects/Microglia_sensory_experience_genomics/snRNAseq_all_nuclei_V1/New_pipeline-snRNAseq/R-analysis_Lin_align/02_Integration_and_clustering/less_restrictive/no_batch_correction/Allexp.sct_nobatch_corrected.40.05.rds")
#Allexp.sct_nobatch_corrected.40.05

#dim(Allexp.sct_nobatch_corrected.40.05)
#table(Allexp.sct_nobatch_corrected.40.05$condition_simpl)
#table(Allexp.sct_nobatch_corrected.40.05$sample, Allexp.sct_nobatch_corrected.40.05$seurat_clusters)

# Select the RNA counts slot to be the default assay
#DefaultAssay(Allexp.sct_nobatch_corrected.40.05) <- "decontXcountsclean"

# Normalize RNA data for visualization purposes
#Allexp.sct_nobatch_corrected.40.05_norm <- NormalizeData(Allexp.sct_nobatch_corrected.40.05, verbose = FALSE)

#saveRDS(Allexp.sct_nobatch_corrected.40.05_norm, "Allexp.sct_nobatch_corrected.40.05_norm.rds")

Allexp.sct_nobatch_corrected.40.05_norm <- readRDS("Allexp.sct_nobatch_corrected.40.05_norm.rds")

DefaultAssay(Allexp.sct_nobatch_corrected.40.05_norm) <- "decontXcountsclean"

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

head(annotations)

#markers<- presto::wilcoxauc(Allexp.sct_batch_corrected.40.02_norm, 'seurat_clusters', assay = 'data')
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(Allexp.sct_nobatch_corrected.40.05_norm,
                       ident.1 = cluster,
                       grouping.var = "condition_simpl",
                       only.pos = TRUE,
                       min.pct = 0.25) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "gene_biotype", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:34), get_conserved)

# Extract top 10 markers per cluster
DFLT_res05_top20 <- conserved_markers %>% 
  mutate(avg_fc = (NR_avg_log2FC + LDR_avg_log2FC + LDR30m_avg_log2FC + LDR2h_avg_log2FC + LDR4h_avg_log2FC + LDR6h_avg_log2FC) /6) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)

# Extract top 3 markers per cluster for the HeATMAP
DFLT_res05_top3 <- conserved_markers %>% 
  mutate(avg_fc = (NR_avg_log2FC + LDR_avg_log2FC + LDR30m_avg_log2FC + LDR2h_avg_log2FC + LDR4h_avg_log2FC + LDR6h_avg_log2FC) /6) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 3, wt = avg_fc) %>%
  .$gene

write.csv(DFLT_res05_top20, file = "DFLT_res05_top20_conserved_markers_final.csv")

DFLT_res05_top20_filtered <- DFLT_res05_top20 %>% dplyr::filter(NR_p_val_adj < 0.01 & NR_pct.1 > 0.7 & NR_pct.2 < 0.3)

write.csv(DFLT_res05_top20_filtered, file = "DFLT_res05_top20_filtered_conserved_markers_final.csv")

write.csv(DFLT_res05_top3, file = "DFLT_res05_top3_genes.csv")
