#' Create Spatial Tissue Plots for Single-Cell Spatial Transcriptomics Data
#'
#' This function generates visualization plots for spatial transcriptomics data, displaying
#' cells/spots on the tissue image with coloring based on a specified grouping variable.
#' It supports custom coordinate ranges, image display options, and appearance customization.
#'
#' @param object A Seurat object containing spatial transcriptomics data with tissue coordinates
#'   and image information.
#' @param group Character string specifying which metadata column to use for coloring spots.
#' @param color_list Named vector of colors where names correspond to group levels.
#' @param point_size Numeric value for the size of points in the plot. Default is 0.5.
#' @param title Character string for the plot title. Default is empty.
#' @param x_label Character string for x-axis label. Default is empty.
#' @param y_label Character string for y-axis label. Default is empty.
#' @param image_scale Character string specifying which image in the Seurat object to use.
#'   This should match a name in object@images. Default is "slice1.008um".
#' @param alpha Numeric value between 0 and 1 for point transparency. Default is 0.5.
#' @param x_range Numeric vector of length 2 specifying the x-axis range to display.
#'   If NULL, shows the entire range. Default is NULL.
#' @param y_range Numeric vector of length 2 specifying the y-axis range to display.
#'   If NULL, shows the entire range. Default is NULL.
#' @param max Numeric value specifying the maximum coordinate value (image height).
#'   If NULL, calculated automatically from data. Default is NULL.
#' @param show_image Logical indicating whether to display the background tissue image.
#'   Default is TRUE.
#'
#' @return A ggplot object representing the spatial plot with cells/spots colored by group.
#'
#' @details
#' This function visualizes spatial transcriptomics data with:
#' 1. Background tissue image (if available and show_image=TRUE)
#' 2. Spots/cells positioned according to their spatial coordinates
#' 3. Coloring based on the specified grouping variable
#'
#' The coordinates are transformed to match image orientation (flipped y-axis). The function
#' also supports cropping to specific regions using x_range and y_range parameters.
#'
#' The 'max' parameter is important for coordinate transformation and should match the image height.
#' Users can determine this by examining their image dimensions.
#'
#' @examples
#' # Load required packages
#' library(Seurat)
#' library(ggplot2)
#' library(dplyr)
#'
#' # Basic usage with default parameters
#' p1 <- stPlot(spatial_seurat_obj,
#'             group = "seurat_clusters",
#'             color_list = c("0" = "red", "1" = "blue", "2" = "green"))
#'
#' # Customize point size and transparency
#' p2 <- stPlot(spatial_seurat_obj,
#'             group = "seurat_clusters",
#'             color_list = c("0" = "red", "1" = "blue", "2" = "green"),
#'             point_size = 1.2,
#'             alpha = 0.8)
#'
#' # Add titles and labels
#' p3 <- stPlot(spatial_seurat_obj,
#'             group = "seurat_clusters",
#'             color_list = c("0" = "red", "1" = "blue", "2" = "green"),
#'             title = "Spatial Distribution of Cell Clusters",
#'             x_label = "X Coordinate (μm)",
#'             y_label = "Y Coordinate (μm)")
#'
#' # Zoom to a specific region
#' p4 <- stPlot(spatial_seurat_obj,
#'             group = "seurat_clusters",
#'             color_list = c("0" = "red", "1" = "blue", "2" = "green"),
#'             x_range = c(2000, 4000),
#'             y_range = c(1000, 3000))
#'
#' # Specify image name and max coordinate value
#' p5 <- stPlot(spatial_seurat_obj,
#'             group = "seurat_clusters",
#'             color_list = c("0" = "red", "1" = "blue", "2" = "green"),
#'             image_scale = "tissue_hires",
#'             max = 10000)
#'
#' # Hide background image to focus on spots
#' p6 <- stPlot(spatial_seurat_obj,
#'             group = "seurat_clusters",
#'             color_list = c("0" = "red", "1" = "blue", "2" = "green"),
#'             show_image = FALSE)
#'
#' # Display the plot
#' print(p1)
#'
#' # Save the plot with high resolution
#' ggsave("spatial_clusters.png", p1, width = 10, height = 8, dpi = 300)
#'
#' @export
stPlot <- function(object, group, color_list, point_size = 0.5, title = "",
                   x_label = "", y_label = "",
                   image_scale = "slice1.008um", alpha = 0.5, x_range = NULL,
                   y_range = NULL, max = NULL, show_image = TRUE) {
  # Load required packages
  require(ggplot2)
  require(dplyr)

  group_data <- object@meta.data[[group]]
  coordinates <- GetTissueCoordinates(object, scale = "hires")

  if (is.null(coordinates) || nrow(coordinates) == 0) {
    stop("坐标数据为空，无法生成图形")
  }

  Xmax <- if (is.null(max)) max(coordinates$y, na.rm = TRUE) else max

  plot_data <- data.frame(
    x = coordinates$y,
    y = Xmax - coordinates$x,
    cell = rownames(coordinates),
    group = group_data
  )

  rm(coordinates, group_data)
  gc()

  xlim_valid <- !is.null(x_range) && length(x_range) == 2
  ylim_valid <- !is.null(y_range) && length(y_range) == 2

  if (xlim_valid) {
    plot_data <- plot_data[plot_data$x >= x_range[1] & plot_data$x <= x_range[2], ]
  }

  if (ylim_valid) {
    plot_data <- plot_data[plot_data$y >= y_range[1] & plot_data$y <= y_range[2], ]
  }

  unique_groups <- unique(plot_data$group)
  color_mapping <- setNames(color_list[unique_groups], unique_groups)

  p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
    scale_color_manual(values = color_mapping) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "right",
      legend.justification = c(-2, 0.5),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22),

      panel.border = element_rect(colour = "black", fill = NA, size = 1),

      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(0.2, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 10, alpha = 1))) +
    labs(title = title, x = x_label, y = y_label)


  if (show_image) {
    image_data <- object@images[[image_scale]]
    if (!is.null(image_data)) {
      img <- image_data@image
      img_width <- ncol(img)
      img_height <- nrow(img)

      if (xlim_valid || ylim_valid) {


        x_crop <- c(
          max(1, floor(x_range[1])),
          min(img_width, ceiling(x_range[2]))
        )


        y_original_x_start <- Xmax - y_range[2]
        y_original_x_end <- Xmax - y_range[1]
        y_crop <- c(
          max(1, floor(y_original_x_start)),
          min(img_height, ceiling(y_original_x_end))
        )

        img_cropped <- img[y_crop[1]:y_crop[2], x_crop[1]:x_crop[2], , drop = FALSE]


        p <- p + annotation_raster(
          img_cropped,
          xmin = x_crop[1],
          xmax = x_crop[2],
          ymin = y_range[1],
          ymax = y_range[2]
        )
      } else {

        p <- p + annotation_raster(
          img,
          xmin = 0, xmax = img_width,
          ymin = 0, ymax = img_height
        )
      }

      rm(image_data, img)
      gc()
    }
  }


  p <- p +
    geom_point(size = point_size, alpha = alpha) +
    coord_cartesian(
      xlim = if (xlim_valid) x_range else c(NA_real_, NA_real_),
      ylim = if (ylim_valid) y_range else c(NA_real_, NA_real_),
      expand = FALSE
    )

  return(p)
}
