#' Title
#' @import ggplot2
#' @importFrom dplyr %>% filter
#' @importFrom Seurat FetchData
#'
#' @param seurat.obj1  the first Seurat object with completed dimensionality reduction analysis
#' @param seurat.obj2  the second Seurat object with completed dimensionality reduction analysis
#' @param cell4link.order  specify the order of clusters for cell4link
#' @param cell4link.color  specify the color of clusters for cell4link
#' @param cellnot4link.color  specify the color of other cells
#' @param line.type  specify the type of lines for linking cells (refer to geom_line in ggplot2)
#' @param line.color  specify the color of lines for linking cells (refer to geom_line in ggplot2)
#' @param line.alpha  specify the transparency of lines for linking cells (refer to geom_line in ggplot2)
#' @param save.format  specify the format for saving the plot
#' @param DR.method  specify a common dimensional reduction method for two Seurat objects (e.g., umap, tsne, pca)
#' @param dim1  the name of Dimension 1 in the dimensional reduction method usually does not need to be specified
#' @param dim2  the name of Dimension 2 in the dimensional reduction method usually does not need to be specified
#'
#' @return return a list of plots for all clusters
#'
#' @examples
#' brain.harmony$cell4link = brain.harmony$orig.ident
#' BatchLinkPlot(
#'     seurat.obj1 = brain.harmony,
#'     seurat.obj2 = brain.cca,
#'     cell4link.order = c("cortex", "hippocampus"),
#'     save.format = "png"
#'     )
#'
#' @export

BatchLinkPlot <- function(seurat.obj1,
                          seurat.obj2,
                          cell4link.order = NULL,
                          cell4link.color = c("#BC3C29", "#A1433D", "#864B51", "#6B5365", "#505A79", "#35628D", "#1A6AA1", "#0072B5", "#2075A0", "#40788C",
                                              "#607B78", "#807E63", "#A0814F", "#C0833B", "#E18727", "#C5862C", "#A98632", "#8E8637", "#72853D", "#578542",
                                              "#3B8548", "#20854E", "#2C825C", "#39806A", "#457E78", "#527C86", "#5E7A94", "#6B78A2", "#7876B1", "#767BB0",
                                              "#7580AF", "#7485AF", "#728AAE", "#718FAE", "#7094AD", "#6F99AD", "#83A2A9", "#98ACA5", "#ACB5A1", "#C1BF9D",
                                              "#D5C899", "#EAD295", "#FFDC91", "#FCC791", "#FAB292", "#F79E93", "#F58994", "#F27595", "#F06096", "#EE4C97"),
                          cellnot4link.color = "#DCDCDC",
                          line.type = "solid",
                          line.color = '#A9A9A9',
                          line.alpha = 0.1,
                          save.format = "pdf",
                          DR.method = "umap",
                          dim1 = NULL,
                          dim2 = NULL) {

  # Define column names based on DR.method
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
      embeddings <- Embeddings(seurat.obj1[[DR.method]])
      dim1 <- colnames(embeddings)[1]
      dim2 <- colnames(embeddings)[2]
      message("Using columns ", dim1, " and ", dim2, " as dimensions for ", DR.method)
    }
  }

  # Check if the number of clusters in cell4link exceeds 50
  num_clusters <- length(unique(seurat.obj1$cell4link))
  if (num_clusters > 50) {
    stop("The number of clusters in cell4link exceeds 50. Please provide a custom color palette with at least ", num_clusters, " colors.")
  }

  # If cell4link.order is not NULL, set cell4link levels according to specified order
  if (!is.null(cell4link.order)) {
    seurat.obj1$cell4link <- factor(seurat.obj1$cell4link, levels = cell4link.order)
  }

  projection.long.total <- list()  # Store results for each cluster

  for (cluster.i in levels(seurat.obj1$cell4link)) {
    # Extract data based on the current cluster
    projection.ref <- seurat.obj1[, seurat.obj1$cell4link %in% cluster.i]

    if (ncol(projection.ref) == 0) {
      warning(paste("No cells found for cluster:", cluster.i))
      next  # Skip clusters with no cells
    }

    seurat.obj2@meta.data$projection.index <- "NO"
    seurat.obj2@meta.data$projection.index[colnames(seurat.obj2) %in% colnames(projection.ref)] <- colnames(projection.ref)

    seurat.obj2.projection <- FetchData(seurat.obj2, vars = c(dim1, dim2, "projection.index"))
    seurat.obj2.projection[[dim1]] <- (seurat.obj2.projection[[dim1]] + 30)

    seurat.obj2.projection.long <- data.frame(
      Name = row.names(seurat.obj2.projection),
      Dim1 = seurat.obj2.projection[[dim1]],
      Dim2 = seurat.obj2.projection[[dim2]],
      projection.index = seurat.obj2.projection$projection.index,
      method = 'seurat.obj2.method'
    )

    seurat.obj1@meta.data$projection.index <- "NO"
    seurat.obj1@meta.data$projection.index[seurat.obj1$cell4link %in% cluster.i] <- colnames(projection.ref)

    seurat.obj1.projection <- FetchData(seurat.obj1, vars = c(dim1, dim2, "projection.index"))
    seurat.obj1.projection.long <- data.frame(
      Name = row.names(seurat.obj1.projection),
      Dim1 = seurat.obj1.projection[[dim1]],
      Dim2 = seurat.obj1.projection[[dim2]],
      projection.index = seurat.obj1.projection$projection.index,
      method = 'seurat.obj1.method'
    )

    projection.long <- rbind(seurat.obj2.projection.long, seurat.obj1.projection.long)
    projection.long.yes <- projection.long %>% filter(projection.index != "NO")
    projection.long.yes.seurat.obj1.method <- projection.long.yes %>% filter(method == "seurat.obj1.method")
    projection.long.yes.seurat.obj1.method$INDEX <- as.factor("yes")
    projection.long.yes.seurat.obj2.method <- projection.long.yes %>% filter(method == "seurat.obj2.method")

    # Plot
    plot <- ggplot() +
      geom_point(data = projection.long, aes(x = Dim1, y = Dim2), color = cellnot4link.color) +
      geom_point(data = projection.long.yes.seurat.obj1.method, aes(x = Dim1, y = Dim2),
                 color = cell4link.color[which(cell4link.order == cluster.i)], shape = 19, size = 3) +
      geom_point(data = projection.long.yes.seurat.obj2.method, aes(x = Dim1, y = Dim2, color = method),
                 fill = cellnot4link.color, shape = 21, size = 3, alpha = 1, stroke = 0.5) +
      scale_color_manual(values = c(cell4link.color[which(cell4link.order == cluster.i)], cellnot4link.color)) +
      geom_line(data = projection.long.yes, aes(x = Dim1, y = Dim2, group = Name),
                linetype = line.type, color = line.color, alpha = line.alpha) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()
      )

    # Save or return plot for each cluster
    projection.long.total[[cluster.i]] <- plot

    # Choose file extension based on save format
    if (save.format == "pdf") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cell4link_cluster_",  cluster.i, ".pdf"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "png") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cell4link_cluster_",  cluster.i, ".png"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "jpeg") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cell4link_cluster_",  cluster.i, ".jpeg"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "tiff") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cell4link_cluster_",  cluster.i, ".tiff"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "bmp") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cell4link_cluster_",  cluster.i, ".bmp"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "svg") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cell4link_cluster_",  cluster.i, ".svg"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else {
      stop("Unsupported save.format. Please use 'pdf', 'png', 'jpeg', 'tiff', 'bmp', or 'svg'.")
    }

    # Display a message indicating the cluster plot has been completed
    message(paste("Cell4link cluster", cluster.i, "has been plotted and saved."))
  }

  return(projection.long.total)  # Return a list of plots for all clusters
}
