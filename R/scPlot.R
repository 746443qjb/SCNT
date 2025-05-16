#' Create Customized Single-Cell Plots with Cluster Labels
#'
#' This function creates publication-quality scatter plots for single-cell data,
#' displaying dimensionality reduction results (e.g., UMAP, t-SNE, PCA) with cluster
#' labels at the median position of each group. The plot includes customizable
#' point size, label size, colors, and legend layout.
#'
#' @param sce A Seurat object containing dimensionality reduction results.
#' @param group_by Character string specifying the metadata column to use for coloring points.
#'   Default is "orig.ident".
#' @param reduction Character string specifying the dimensionality reduction to use.
#'   Default is "umap".
#' @param point_size Numeric value for the size of points in the plot. Default is 1.2.
#' @param label_size Numeric value for the size of cluster labels. Default is 5.
#' @param colors Optional vector of colors to use for clusters. If NULL, default color
#'   palette will be used.
#' @param legend_cols Number of columns to use in the legend. Default is 2.
#' @param raster_dpi Resolution (dots per inch) for rasterized points. Higher values
#'   give better quality but larger file size. Default is 600.
#'
#' @return A ggplot object representing the dimensionality reduction plot with
#'   clusters colored according to the specified grouping variable and labeled
#'   at the median position of each cluster.
#' @export
scPlot <- function(sce, group_by = "orig.ident", reduction = "umap", point_size = 1.2,
                   label_size = 5, colors = NULL, legend_cols = 2, raster_dpi = 600) {
  # Load required packages
  require(ggplot2)
  require(dplyr)
  require(ggrastr)
  require(rlang)


  metaData <- sce@meta.data
  dim_data <- Embeddings(sce, reduction = reduction)


  dim_names <- paste0(toupper(reduction), c(" 1", " 2"))
  colnames(dim_data) <- c("dim_1", "dim_2")

  metaData <- metaData %>%
    cbind(dim_1 = dim_data[, 1], dim_2 = dim_data[, 2])

  plotData <- metaData %>%
    dplyr::select(dim_1, dim_2, !!sym(group_by)) %>%
    dplyr::mutate(cellClusters = !!sym(group_by))

  IdentLabel <- plotData %>%
    group_by(cellClusters) %>%
    summarise(
      dim_1 = median(dim_1),
      dim_2 = median(dim_2)
    )
  IdentLabel$Label <- IdentLabel$cellClusters

  if (is.null(colors)) {
    default_colors <- c(
      "#FB9A99", "#377EB8", "#4DAF4A", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", "#F5D2A8",
      "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#E7298A", "#E41A1C", "#D2EBC8", "#00CDD1",
      "#5E4FA2", "#8CA77B", "#8DD3C7", "#7DBFC7", "#B3DE69", "#999999", "#FCED82",
      "#F5CFE4", "#F5D2A8", "#B383B9", "#EE934E", "#BBDD78"
    )
    n_clusters <- nrow(IdentLabel)
    if (n_clusters <= length(default_colors)) {
      cell_colors <- default_colors[1:n_clusters]
    } else {
      cell_colors <- colorRampPalette(default_colors)(n_clusters)
    }
  } else {
    n_clusters <- nrow(IdentLabel)
    if (length(colors) < n_clusters) {
      warning("提供的颜色数量不足，将循环使用")
      cell_colors <- rep(colors, length.out = n_clusters)
    } else {
      cell_colors <- colors[1:n_clusters]
    }
  }

  sorted_clusters <- sort(unique(IdentLabel$cellClusters))
  names(cell_colors) <- sorted_clusters


  p <- ggplot(plotData, aes(x = dim_1, y = dim_2, colour = cellClusters)) +
    geom_point_rast(size = point_size, stroke = 0, shape = 16, raster.dpi = raster_dpi) +
    geom_point(data = IdentLabel, aes(color = cellClusters),
               size = 15, alpha = 0.6, shape = 16) +
    geom_text(data = IdentLabel, aes(label = Label), colour = "black", size = label_size) +
    scale_colour_manual(values = cell_colors) +
    labs(x = dim_names[1], y = dim_names[2], color = group_by) +
    theme(
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "right",
      legend.spacing.x = unit(0, "cm")
    ) +
    guides(
      color = guide_legend(
        ncol = legend_cols,
        override.aes = list(size = 10, stroke = 0)
      )
    )

  return(p)
}
# 示例用法：
# p1 <- scPlot(seurat_obj, group_by = "orig.ident", reduction = "umap", point_size = 1.2)
