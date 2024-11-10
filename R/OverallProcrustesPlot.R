#' Title
#' @name OverallProcrustesPlot
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
#' OverallProcrustesPlot(seurat.obj1 = brain.harmony, seurat.obj2 = brain.cca, DR.method = "pca")
#'
#' @export

# Load necessary libraries
library(ggplot2)
library(patchwork)
library(vegan)

# Define the Procrustes analysis function
OverallProcrustesPlot <- function(seurat.obj1, seurat.obj2, DR.method = "umap", dim1 = NULL, dim2 = NULL) {
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
  
  # Unique cluster IDs
  cell4link_clusters <- unique(seurat.obj1$cell4link)
  
  # Create directories to save plots
  dir.create("Procrustes_plots", showWarnings = FALSE)
  
  # Perform Procrustes analysis and plotting
  perform_procrustes_analysis <- function(DR.res1, DR.res2) {
    # Perform Procrustes analysis and permutation test
    pro.s.e <- suppressMessages(procrustes(DR.res1, DR.res2, symmetric = TRUE))
    pro.s.e_t <- suppressMessages(protest(DR.res1, DR.res2, permutations = 999))
    summary_pro <- summary(pro.s.e)
    # cat("Summary for all clusters", ":\n")
    # print(summary_pro)
    
    # Save Procrustes error plot
    pdf("Procrustes_plots/Procrustes_errors_all_clusters.pdf")
    plot(pro.s.e, kind = 2, main = "Procrustes errors for all clusters")
    dev.off()
    
    # Extract M2 and p-value
    M2_value <- pro.s.e_t$ss
    p_value <- pro.s.e_t$signif
    
    # Coordinates for plotting
    Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
    Pro_X <- data.frame(pro.s.e$rotation)
    x_min <- min(c(Pro_Y$X1, Pro_Y[[dim1]])) - 0.005
    x_max <- max(c(Pro_Y$X1, Pro_Y[[dim1]])) + 0.005
    y_min <- min(c(Pro_Y$X2, Pro_Y[[dim2]])) - 0.005
    y_max <- max(c(Pro_Y$X2, Pro_Y[[dim2]])) + 0.005

    message("Using dimensions: ", dim1, " and ", dim2)

    # Define individual plots
    plot1 <- ggplot(data = Pro_Y, aes_string(x = dim1, y = dim2)) +
      geom_point(color = "#B0E0E6", size = 2, shape = 21, alpha = 1) +
      labs(x = "Dimension 1", y = "Dimension 2", title = "Transformed DR.res1 - all clusters") +
      geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      theme_minimal() + xlim(x_min, x_max) + ylim(y_min, y_max)
    
    plot2 <- ggplot(data = Pro_Y, aes(x = X1, y = X2)) +
      geom_point(color = "#D8BFD8", size = 2, shape = 16, alpha = 1) +
      labs(x = "Dimension 1", y = "Dimension 2", title = "Transformed DR.res2 - all clusters") +
      geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
      theme_minimal() + xlim(x_min, x_max) + ylim(y_min, y_max)
    
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
    
    # Combine plots and save
    combined_plot <- plot1 | plot2 | plot3
    suppressMessages(ggsave(filename = "Procrustes_plots/Procrustes_plot_all_clusters.pdf", 
           plot = combined_plot, width = 525, height = 175, units = "mm", limitsize = FALSE))
  }
  
  # Run overall Procrustes analysis
  perform_procrustes_analysis(DR.res1, DR.res2)

  cat("Overall Procrustes analysis completed, plots successfully saved!\n")
  
 }
