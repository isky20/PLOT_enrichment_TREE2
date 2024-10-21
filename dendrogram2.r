library(ape)
library(dendextend)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(ggtreeExtra)

# Define the main function
create_dendrogram_plots <- function(file_path, task_name, width = 200, offset= 5, xlim_max= 20, height = 240) {
  # Load data
  data <- read.csv(file_path)
  
  # Create a subdataframe based on unique categories
  categories <- unique(data$category)
  sub_ref_list <- split(data, data$category)
  
  # Normalization function
  normalize <- function(data) {
    # Extract 'id' and 'jds' columns
    subdata_AV <- data[, grep("^Av_", colnames(data))]
    
    # Normalize row-wise
    normalized_df <- t(apply(subdata_AV, 1, function(row) {
      max_value <- max(row)  # Find the maximum value in the current row
      normalized_row <- (row / max_value) * 100  # Normalize each value in the row
      return(normalized_row)
    }))
    
    # Convert back to dataframe and combine with non-normalized columns
    normalized_df <- as.data.frame(normalized_df)
    therest <- data[, grep("^Av_", colnames(data), invert = TRUE)]
    return(cbind(therest, normalized_df))
  }
  
  # Normalize each category in sub_ref_list
  normalize_list <- lapply(sub_ref_list, normalize)
  
  # Function to create plots
  create_plot <- function(df, title) {
    av_columns <- grep("Av_", colnames(df), value = TRUE)
    all_genes <- unique(unlist(strsplit(gsub(" ", "", df$MergedGenes), ",")))
    
    # Create binary matrix
    binary_matrix <- matrix(0, nrow = nrow(df), ncol = length(all_genes), dimnames = list(df$description, all_genes))
    for (i in seq_along(df$MergedGenes)) {
      gene_list <- unlist(strsplit(gsub(" ", "", df$MergedGenes[i]), ","))
      binary_matrix[df$description[i], intersect(all_genes, gene_list)] <- 1
    }
    binary_position_matrix <- (binary_matrix == 1)
    
    # Compute Jaccard similarity matrix
    dist_matrix <- dist(binary_position_matrix, method = "binary")
    hclust_object <- hclust(dist_matrix, method = "ward.D2")
    dendrogram_object <- as.dendrogram(hclust_object)
    newick <- as.phylo(dendrogram_object)
    
    # Create initial ggtree plot
    p1 <- ggtree(newick, branch.length = "none", layout = "circular", open.angle = 15) + 
      geom_tiplab(align = TRUE, linesize = 0, size = 2.5, offset = offset) + 
      ggtitle(title) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      xlim(0, xlim_max)
    
    # Add columns to the plot
    for (case in av_columns) {
      p1 <- p1 + 
        geom_fruit(
          data = df,
          geom = geom_col,
          mapping = aes_string(y = "description", x = case, fill = shQuote(case))
        )
    }
    
    # Save the plot
    ggsave(paste0(title, "_", task_name, ".svg"), plot = p1, device = "svg", 
           width = width, height = height, units = "mm", dpi = 1000)
  }
  
  # Create plots for each normalized dataframe
  for (df in normalize_list) {
    create_plot(df, unique(df$category))
  }
}

# Example usage
create_dendrogram_plots("ATTR_AL_FULL_filter_adj.csv", "ATTR_AL_FULL_filter",
                        width = 200, offset= 6, xlim_max= 25, height = 300)
