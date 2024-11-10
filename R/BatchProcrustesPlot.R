#' Title
#' @name BatchProcrustesPlot
#' @import ggplot2
#' @import vegan
#' @import patchwork
#' @importFrom dplyr %>% filter
#' @importFrom Seurat FetchData
#'
#' @param seurat.obj1  the first Seurat object with completed dimensionality reduction analysis
#' @param seurat.obj2  the second Seurat object with completed dimensionality reduction analysis
#' @param DR.method  specify a common dimensional reduction method for two Seurat objects (e.g., umap, tsne, pca)
#'
#' @examples
#' brain.harmony$cell4link = brain.harmony$orig.ident
#' BatchProcrustesPlot(seurat.obj1 = brain.harmony, seurat.obj2 = brain.cca, DR.method = "pca")
#'
#' @export

# Load necessary libraries
library(ggplot2)
library(patchwork)
library(vegan)

BatchProcrustesPlot <- function(seurat.obj1, seurat.obj2, DR.method = "umap", dim1 = NULL, dim2 = NULL) {
  
  # Define dimension names based on the DR.method
  if (DR.method == "umap") {
    dim1 <- "UMAP_1"
    dim2 <- "UMAP_2"
  } else if (DR.method == "pca") {
    dim1 <- "PC_1"
    dim2 <- "PC_2"
  } else if (DR.method == "tsne") {
    dim1 <- "tSNE_1"
    dim2 <- "tSNE_2"
  } else {
    if (is.null(dim1) || is.null(dim2)) {
      embeddings <- suppressMessages(Embeddings(seurat.obj1[[DR.method]]))
      dim1 <- colnames(embeddings)[1]
      dim2 <- colnames(embeddings)[2]
      message("Using columns ", dim1, " and ", dim2, " as dimensions for ", DR.method)
    }
  }
  
  # Extract embeddings for the two Seurat objects
  DR.res1 <- suppressMessages(Embeddings(seurat.obj1[[DR.method]]))
  DR.res2 <- suppressMessages(Embeddings(seurat.obj2[[DR.method]]))
  
  # Find common cells and align embeddings
  common_cells <- intersect(rownames(DR.res1), rownames(DR.res2))
  DR.res1 <- DR.res1[common_cells, ]
  DR.res2 <- DR.res2[common_cells, ]

  # Create directories for saving plots if not exist
  dir.create("Procrustes_plots", showWarnings = FALSE)

  # Iterate over each unique cluster ID in seurat.obj1$cell4link
  for (cluster.i in unique(seurat.obj1$cell4link)) {
    cluster_cells <- which(seurat.obj1$cell4link[common_cells] == cluster.i)
    
    if (length(cluster_cells) == 0) {
      message("No cells found for Cluster ", cluster.i, " in both objects.")
      next
    }
    
    # Extract embeddings for the current cluster
    DR.res1_cluster <- DR.res1[cluster_cells, ]
    DR.res2_cluster <- DR.res2[cluster_cells, ]
    
    # Perform Procrustes analysis
    pro.s.e <- suppressMessages(procrustes(DR.res1_cluster, DR.res2_cluster, symmetric = TRUE))
    pro.s.e_t <- suppressMessages(protest(DR.res1_cluster, DR.res2_cluster, permutations = 999))
    
    # Output Procrustes analysis summary
    summary_pro <- summary(pro.s.e)
    # cat("Summary for Cluster", cluster.i, ":\n")
    # print(summary_pro)
    
    # Save Procrustes error plot
    pdf(paste0("Procrustes_plots/Procrustes_errors_cluster_", cluster.i, ".pdf"))
    plot(pro.s.e, kind = 2, main = paste("Procrustes errors for Cluster", cluster.i))
    dev.off()
    
    # Set M2 and p-value
    M2_value <- pro.s.e_t$ss
    p_value <- pro.s.e_t$signif
    
    # Obtain coordinates for x and y axes and rotated coordinates
    Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
    
    # Calculate x and y axis range
    x_min <- min(c(Pro_Y$X1, Pro_Y[[dim1]])) - 0.005
    x_max <- max(c(Pro_Y$X1, Pro_Y[[dim1]])) + 0.005
    y_min <- min(c(Pro_Y$X2, Pro_Y[[dim2]])) - 0.005
    y_max <- max(c(Pro_Y$X2, Pro_Y[[dim2]])) + 0.005

    Pro_Y[[dim1]] <- as.numeric(Pro_Y[[dim1]])
    Pro_Y[[dim2]] <- as.numeric(Pro_Y[[dim2]])
    Pro_Y$X1 <- as.numeric(Pro_Y$X1)
    Pro_Y$X2 <- as.numeric(Pro_Y$X2)
    
    # Plot 1: Transformed DR.res1 for the current cluster
    plot1 <- ggplot(data = Pro_Y, aes_string(x = dim1, y = dim2)) +
      geom_point(color = "#B0E0E6", size = 2, shape = 21, alpha = 1) +
      labs(x = "Dimension 1", y = "Dimension 2", title = paste("Transformed DR.res1 - Cluster", cluster.i)) +
      geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      theme_minimal() + xlim(x_min, x_max) + ylim(y_min, y_max)
    
    # Plot 2: Transformed DR.res2 for the current cluster
    plot2 <- ggplot(data = Pro_Y, aes(x = X1, y = X2)) +
      geom_point(color = "#D8BFD8", size = 2, shape = 16, alpha = 1) +
      labs(x = "Dimension 1", y = "Dimension 2", title = paste("Transformed DR.res2 - Cluster", cluster.i)) +
      geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      theme_minimal() + xlim(x_min, x_max) + ylim(y_min, y_max)
    
    # Plot 3: Combined DR and Procrustes transformation for the current cluster
    plot3 <- ggplot(data = Pro_Y) +
      geom_segment(aes(x = X1, y = X2, xend = (X1 + Pro_Y[[dim1]]) / 2, yend = (X2 + Pro_Y[[dim2]]) / 2),
                   arrow = arrow(length = unit(0, 'cm')), color = "#D8BFD8", size = 1) +
      geom_segment(aes(x = (X1 + Pro_Y[[dim1]]) / 2, y = (X2 + Pro_Y[[dim2]]) / 2, xend = Pro_Y[[dim1]], yend = Pro_Y[[dim2]]),
                   arrow = arrow(length = unit(0.2, 'cm')), color = "#B0E0E6", size = 1) +
      geom_point(aes(x = X1, y = X2), color = "#D8BFD8", size = 2, shape = 16, alpha = 0.5) +
      geom_point(aes_string(x = dim1, y = dim2), color = "#B0E0E6", size = 2, shape = 21, alpha = 0.5) +
      labs(x = 'Dimension 1', y = 'Dimension 2') +
      annotate('text', label = sprintf("Procrustes analysis:\nM2 = %.4f, p-value = %.3f", M2_value, p_value),
               x = x_min, y = y_max, size = 4, hjust = 0) +
      theme_minimal() + xlim(x_min, x_max) + ylim(y_min, y_max)
    
    # Combine the three plots for the current cluster
    combined_plot <- plot1 | plot2 | plot3
    
    # Save the combined Procrustes plot for the current cluster
    suppressMessages(ggsave(filename = paste0("Procrustes_plots/Procrustes_DR_cluster_", cluster.i, ".pdf"), 
           plot = combined_plot, width = 525, height = 175, units = "mm", limitsize = FALSE))
    
    cat("Cluster", cluster.i, "Procrustes analysis completed, plots successfully saved!\n")
  }
}
