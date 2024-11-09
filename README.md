# cellLink

When analyzing single-cell transcriptomics data, we often perform various analyses on the same cells, which generates different dimensional reduction results. This makes it challenging to locate target cells across different reduction results and to compare them. 

***CellLink*** is a lightweight R package designed to help project shared cells (target cell groups) across different dimensional reduction results, providing a function for batch visualization. It identifies target cell groups in various dimensional reduction results and connects them with lines to facilitate comparison.

## Features

- Map cells between various dimensionality reduction plots (e.g., UMAP, t-SNE).
- Visually track specific cells across different analysis results.
- Specify cells for mapping based on metadata (e.g., `orig.ident` or `seurat_clusters`).
- Customizable color and order of cell clusters in plots.
- Save plots in various formats, including PDF, PNG, JPEG, and more.

## Installation

To install the development version of `cellLink` from GitHub, you can use the following commands:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install cellLink from GitHub
devtools::install_github("mushan-guan/cellLink")
```

## Usage
### Dependency
The cellLink package depends on the Seurat, ggplot2, and dplyr packages.
### Loading the Package
Load the cellLink package into your R session:
```r
library(cellLink)
```
### Example Workflow
This example demonstrates how to use the BatchLinkPlot function to map cells between two Seurat objects with different dimensionality reduction results.

Step 1: Prepare Seurat Objects
```r
# Assuming you have two Seurat objects (e.g., seurat.obj1 and seurat.obj2)
```
Step 2: Map Cells Between Dimensionality Reductions
Use the BatchLinkPlot function to map cells between the two Seurat objects. Customize the color palette, cluster order, and output format as needed.
```r
# choose cells for plot
seurat.obj1@meta.data$cell4plot = seurat.obj1@meta.data$cluster

# Map cells between the two Seurat objects (from seurat.obj1 to seurat.obj2)
mapped_plots <- BatchLinkPlot(
  seurat.obj1 = seurat.obj1,
  seurat.obj2 = seurat.obj2,
  cell4link.order = c("Cluster1", "Cluster2", "Cluster3"),  # Specify the order of clusters
  cell4link.color = c("#BC3C29", "#A1433D", "#864B51"),     # Customize the colors of clusters for linking
  cellnot4link.color = "#DCDCDC",                           # Customize the colors of other cells (not for linking)
  save.format = "pdf"                                       # Save plots as PDF files
)
```
Step 3: Plot Saving
The BatchLinkPlot function will save the generated plots to your working directory in the specified format (e.g., PDF). You can also access individual plots in the mapped_plots list for further customization or visualization.
```r
mapped_plots
# $`Cluster1`

# $`Cluster2`

# $`Cluster3`

class(mapped_plots$`Cluster1`)
# [1] "gg"     "ggplot"
```

## Parameters
 - seurat.obj1: First Seurat object, typically with an initial dimensionality reduction.
 - seurat.obj2: Second Seurat object with an alternative dimensionality reduction.
 - cell4link.order: Order of clusters to be displayed in the plot.
 - cell4link.color: Color palette for clusters, specified as a vector of color codes ***(A default palette of 50 colors is provided)***.
 - cellnot4link.color: Color of other cells ***(Default: #DCDCDC)***.
 - line.type: Type of lines for linking cells (refer to geom_line in ggplot2) ***(Default: solid)***.
 - line.color: Color of lines for linking cells (refer to geom_line in ggplot2) ***(Default: #A9A9A9)***.
 - line.alpha: Transparency of lines for linking cells (refer to geom_line in ggplot2) ***(Default: 0.1)***.
 - save.format: Format for saving plots (e.g., "pdf", "png", "jpeg", etc.) ***(Default: pdf)***.
 - DR.method: specify a common dimensional reduction method for two Seurat objects (e.g., umap, tsne, pca) ***(Default: umap)***.
 - dim1: the name of Dimension 1 in the dimensional reduction method (usually does not need to be specified)
 - dim2: the name of Dimension 2 in the dimensional reduction method (usually does not need to be specified)

## Notes
 - Ensure the cells used for plotting are in the intersection of both Seurat objects, with consistent barcodes across objects.
 - Adjust cell4link.order and cell4link.color to customize the plot order and color scheme.

# License
This package is licensed under the MIT License. See the LICENSE file for more details.

# Acknowledgments
The cellLink package relies on the following packages for data manipulation and visualization:

 - Seurat for handling single-cell data.
 - ggplot2 for plotting and visualizations.
 - dplyr for data manipulation.
 - Special thanks to the contributors and users who provided valuable feedback to improve the package.

# Contact
For help, please contact qibiaoguan@163.com.
