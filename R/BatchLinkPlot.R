#' Title
#' @import ggplot2
#' @importFrom dplyr %>% filter
#' @importFrom Seurat FetchData
#'
#' @param seurat.obj1 the first Seurat object with completed dimensionality reduction analysis
#' @param seurat.obj2 the second Seurat object with completed dimensionality reduction analysis
#' @param cell4plot.order  specify the order of clusters for cell4plot
#' @param cell4plot.color  specify the color of clusters for cell4plot
#' @param save.format  specify the format for saving the plot
#'
#' @return return a list of plots for all clusters
#'
#' @examples
#' brain.harmony$cell4plot = brain.harmony$orig.ident
#' BatchLinkPlot(
#'     seurat.obj1 = brain.harmony,
#'     seurat.obj2 = brain.cca,
#'     cell4plot.order = c("cortex", "hippocampus"),
#'     save.format = "png"
#'     )
#'
#' @export

BatchLinkPlot <- function(seurat.obj1,
                          seurat.obj2,
                          cell4plot.order = NULL,
                          cell4plot.color = c("#BC3C29", "#A1433D", "#864B51", "#6B5365", "#505A79", "#35628D", "#1A6AA1", "#0072B5", "#2075A0", "#40788C",
                                              "#607B78", "#807E63", "#A0814F", "#C0833B", "#E18727", "#C5862C", "#A98632", "#8E8637", "#72853D", "#578542",
                                              "#3B8548", "#20854E", "#2C825C", "#39806A", "#457E78", "#527C86", "#5E7A94", "#6B78A2", "#7876B1", "#767BB0",
                                              "#7580AF", "#7485AF", "#728AAE", "#718FAE", "#7094AD", "#6F99AD", "#83A2A9", "#98ACA5", "#ACB5A1", "#C1BF9D",
                                              "#D5C899", "#EAD295", "#FFDC91", "#FCC791", "#FAB292", "#F79E93", "#F58994", "#F27595", "#F06096", "#EE4C97"),
                          save.format = "pdf") {

  # Check if the number of clusters in cell4plot exceeds 50
  num_clusters <- length(unique(seurat.obj1$cell4plot))
  if (num_clusters > 50) {
    stop("The number of clusters in cell4plot exceeds 50. Please provide a custom color palette with at least ", num_clusters, " colors.")
  }

  # If cell4plot.order is not NULL, set cell4plot levels according to specified order
  if (!is.null(cell4plot.order)) {
    seurat.obj1$cell4plot <- factor(seurat.obj1$cell4plot, levels = cell4plot.order)
  }

  projection.long.total <- list()  # Store results for each cluster

  for (cluster.i in levels(seurat.obj1$cell4plot)) {
    # Extract data based on the current cluster
    projection.ref <- seurat.obj1[, seurat.obj1$cell4plot %in% cluster.i]

    if (ncol(projection.ref) == 0) {
      warning(paste("No cells found for cluster:", cluster.i))
      next  # Skip clusters with no cells
    }

    seurat.obj2@meta.data$projection.index <- "NO"
    seurat.obj2@meta.data$projection.index[colnames(seurat.obj2) %in% colnames(projection.ref)] <- colnames(projection.ref)

    seurat.obj2.projection <- FetchData(seurat.obj2, vars = c("UMAP_1", "UMAP_2", "projection.index"))
    seurat.obj2.projection$UMAP_1 <- (seurat.obj2.projection$UMAP_1 + 30)

    seurat.obj2.projection.long <- data.frame(
      Name = row.names(seurat.obj2.projection),
      UMAP_1 = seurat.obj2.projection$UMAP_1,
      UMAP_2 = seurat.obj2.projection$UMAP_2,
      projection.index = seurat.obj2.projection$projection.index,
      method = 'seurat.obj2.method'
    )

    seurat.obj1@meta.data$projection.index <- "NO"
    seurat.obj1@meta.data$projection.index[seurat.obj1$cell4plot %in% cluster.i] <- colnames(projection.ref)

    seurat.obj1.projection <- FetchData(seurat.obj1, vars = c("UMAP_1", "UMAP_2", "projection.index"))
    seurat.obj1.projection.long <- data.frame(
      Name = row.names(seurat.obj1.projection),
      UMAP_1 = seurat.obj1.projection$UMAP_1,
      UMAP_2 = seurat.obj1.projection$UMAP_2,
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
      geom_point(data = projection.long, aes(x = UMAP_1, y = UMAP_2), color = '#DCDCDC') +
      geom_point(data = projection.long.yes.seurat.obj1.method, aes(x = UMAP_1, y = UMAP_2),
                 color = cell4plot.color[which(cell4plot.order == cluster.i)], shape = 19, size = 3) +
      geom_point(data = projection.long.yes.seurat.obj2.method, aes(x = UMAP_1, y = UMAP_2, color = method),
                 fill = "#DCDCDC", shape = 21, size = 3, alpha = 1, stroke = 0.5) +
      scale_color_manual(values = c(cell4plot.color[which(cell4plot.order == cluster.i)], "#DCDCDC")) +
      geom_line(data = projection.long.yes, aes(x = UMAP_1, y = UMAP_2, group = Name),
                linetype = "solid", color = '#A9A9A9', alpha = 0.1) +
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
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cluster_",  cluster.i, "_theme.pdf"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "png") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cluster_",  cluster.i, "_theme.png"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "jpeg") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cluster_",  cluster.i, "_theme.jpeg"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "tiff") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cluster_",  cluster.i, "_theme.tiff"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "bmp") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cluster_",  cluster.i, "_theme.bmp"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else if (save.format == "svg") {
      ggsave(plot = plot, filename = paste0("BatchLinkPlot_cluster_",  cluster.i, "_theme.svg"), width = 850, height = 350, units = "mm", limitsize = FALSE)
    } else {
      stop("Unsupported save.format. Please use 'pdf', 'png', 'jpeg', 'tiff', 'bmp', or 'svg'.")
    }

    # Display a message indicating the cluster plot has been completed
    message(paste("Cell4plot cluster", cluster.i, "has been plotted and saved."))
  }

  return(projection.long.total)  # Return a list of plots for all clusters
}
