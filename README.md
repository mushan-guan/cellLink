# cellLink

**cellLink** is an lightweight R package designed to help map cells between different dimensionality reduction plots generated by Seurat. This visualization tool allows you to track specific cells across multiple dimensionality reduction results and connect them with lines for easy comparison.

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
# Map cells between the two Seurat objects
mapped_plots <- BatchLinkPlot(
  seurat.obj1 = seurat.obj1,
  seurat.obj2 = seurat.obj2,
  cell4plot.order = c("Cluster1", "Cluster2", "Cluster3"),  # Specify the order of clusters
  cell4plot.color = c("#BC3C29", "#A1433D", "#864B51"),     # Customize cluster colors
  save.format = "pdf"                                       # Save plots as PDF files
)
```
Step 3: Plot Saving
The BatchLinkPlot function will save the generated plots to your working directory in the specified format (e.g., PDF). You can also access individual plots in the mapped_plots list for further customization or visualization.

## Parameters
 - seurat.obj1: First Seurat object, typically with an initial dimensionality reduction.
 - seurat.obj2: Second Seurat object with an alternative dimensionality reduction.
 - cell4plot.order: Order of clusters to be displayed in the plot.
 - cell4plot.color: Color palette for clusters, specified as a vector of color codes.
 - save.format: Format for saving plots (e.g., "pdf", "png", "jpeg", etc.).

## Notes
 - Ensure the cells used for plotting are in the intersection of both Seurat objects, with consistent barcodes across objects.
 - Adjust cell4plot.order and cell4plot.color to customize the plot order and color scheme.

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