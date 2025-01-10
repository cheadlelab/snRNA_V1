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
data <- data %>% group_by(Group) %>% top_n(n = -10, wt = PAdjust)

# Calculate global minimum and maximum of PAdjust (assuming you want to use -log10(PAdjust) for size)
global_min <- min(-log10(data$PAdjust), na.rm = TRUE)
global_max <- max(-log10(data$PAdjust), na.rm = TRUE)

# Define cell types and conditions
cells <- c("ExcL23","ExcL4","ExcL6a","ExcL6IT")
cond <- c("induction","repression")
# Define a consistent color palette for the groups
group_colors <- c("LDR_vs_LDR30m" = "blue", "LDR30m_vs_LDR2h" = "green", "LDR2h_vs_LDR4h" = "red", "LDR4h_vs_LDR6h" = "purple")

t <- 0
ww <- c(4.5, 3.8, 3.8, 5.5, 5, 6.5, 5.5, 4.5)
hh <- c(2, 1.8, 1.5, 2.2, 2, 2.5, 2, 2)


# Loop through each cell type and condition
for (cell in cells) {
  for (condition in cond) {
    # Filter data for the current cell type and condition
    ranked_data <- data %>%
      filter(grepl(cell, Group) & grepl(condition, Group)) %>%
      group_by(Group) %>%
      arrange(desc(PAdjust)) %>%
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
      
      
      # Calculate the plot dimensions
      num_pathways <- length(unique(ranked_data$Pathway))
      num_groups <- length(unique(ranked_data$Group))
      
      t <- t+1
      plot_width <- ww[t] # Reduced width
      plot_height <- hh[t] # Reduced height
      
      # Create and save the plot
      p <- ggplot(ranked_data, aes(x = Pathway, y = Group, size = -log10(PAdjust), color = Group)) +
        geom_point(alpha = 0.6) +
        scale_size_continuous(range = c(global_min, global_max)) +
        theme_minimal() +
        theme(text = element_text(family = "Arial", size = 5),
              plot.title = element_text(size = 6),
              axis.text.x = element_text(angle = 30, hjust = 1),
              axis.text.y = element_text(size = 5),
              legend.position = "none",  # Hide all legends
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        labs(title = paste(cell, condition))
      
      # Uncomment the following line to save the plot
      ggsave(paste0("~/Desktop/bubble_plot_", cell, "_", condition, ".svg"), plot = p, width = plot_width, height = plot_height)
    }
  }
}







p <- ggplot(ranked_data, aes(x = Pathway, y = Group, size = -log10(PAdjust), color = Group)) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(range = c(1,3), limits = c(global_min, global_max)) +  # Use global min and max here for size
  scale_color_manual(values = group_colors) + # Use the defined color palette
  theme_minimal() +
  theme(text = element_text(family = "Arial", size = 5),
        plot.title = element_text(size = 6),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(size = 5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = paste(cell, condition),
       size = "-log10(PAdjust)")

# Uncomment the following line to save the plot
ggsave(paste0("~/Desktop/bubble_plot_legend.svg"), plot = p, width = plot_width, height = plot_height)