rm(list = ls(all = TRUE))
# Load necessary libraries
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(egg)
library(ggplot2)

# Set the file path
path_string <- "/media/linqianyu/MOSS/snRNAseq_all_nuclei_V1/figures/pairwise2_L2346_velocity/"
data <- read.csv(paste0(path_string, "final_results.csv"))
data <- data %>% mutate(OriginalOrder = as.numeric(factor(Group, unique(Group))))

# Define cell types and conditions
cells <- c("ExcL23","ExcL4","ExcL6a","ExcL6IT")
cond <- c("induction","repression")

# Loop through each cell type and condition
for (cell in cells) {
  for (condition in cond) {
    # Filter data for the current cell type and condition
    ranked_data <- data %>%
      filter(grepl(cell, Group) & grepl(condition, Group)) %>%
      group_by(Group) %>%
      arrange(desc(Count), desc(PAdjust)) %>%
      ungroup() %>%
      arrange(OriginalOrder)
    
    ranked_data$Group <- gsub(paste0(cell, "_"), "", ranked_data$Group)
    ranked_data$Group <- gsub(paste0("_", condition), "", ranked_data$Group)
    
    # Check if there is data to plot
    if (nrow(ranked_data) > 0) {
      # Set the order of 'Pathway' based on its appearance in the dataset
      ranked_data <- ranked_data %>%
        mutate(Pathway = factor(Pathway, levels = unique(Pathway)))
      
      ranked_data$Group <- factor(ranked_data$Group, levels = unique(ranked_data$Group))
      
      # Convert GeneRatio to a continuous measure
      ranked_data <- ranked_data %>%
        mutate(GeneRatioValue =  as.numeric(sub("/.*$", "", GeneRatio))/as.numeric(sub("^.*/", "", GeneRatio)))
      
      # Calculate the plot dimensions
      num_pathways <- length(unique(ranked_data$Pathway))
      num_groups <- length(unique(ranked_data$Group))
      
      plot_width <- num_pathways / 8 # Reduced width
      plot_height <- num_groups / 1.2 # Reduced height
      
      # Create and save the plot
      p <- ggplot(ranked_data, aes(x = Pathway, y = Group, size = GeneRatioValue, color = PAdjust)) +
        geom_point(alpha = 0.6) +
        scale_size_continuous(range = c(1, 3)) +
        scale_color_gradient(low = "blue", high = "red", limits = c(0, 0.05)) + # Adjusted here
        theme_minimal() +
        theme(text = element_text(family = "Arial", size = 5),
              plot.title = element_text(size = 6),
              axis.text.x = element_text(angle = 30, hjust = 1),
              legend.text = element_text(size = 5), 
              legend.title = element_text(size = 5)
              
        ) +
        labs(title = paste(cell, condition),
             x = "Pathway",
             y = "Group",
             size = "GeneRatio",
             color = "P.Adjust") +
        guides(color = guide_colorbar(title = "P.Adjust"))
      
      
      # Uncomment the following line to save the plot
      ggsave(paste0("~/Desktop/bubble_plot_", cell, "_", condition, ".svg"), plot = p, width = plot_width, height = plot_height)
    }
  }
}
