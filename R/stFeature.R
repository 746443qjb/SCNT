#' Visualize Gene/Feature Expression in Spatial Transcriptomics Data
#'
#' This function creates a spatial plot visualizing the expression of a single gene or feature
#' across the tissue section. Only cells expressing the feature above a minimum threshold
#' will be displayed, making it easier to visualize expression patterns without background noise.
#' The expression level is represented by a color gradient, and the plot can be overlaid
#' on the tissue histology image.
#'
#' @param object A Seurat object containing spatial transcriptomics data with tissue coordinates
#'   and image information.
#' @param feature Character string specifying the gene or feature name to visualize.
#'   Default is "gene".
#' @param color Vector of two colors defining the low-to-high expression gradient.
#'   Default is c("blue", "red").
#' @param point_size Numeric value for the size of points in the plot. Default is 0.5.
#' @param title Character string for the plot title. Default is empty.
#' @param x_label Character string for x-axis label. Default is empty.
#' @param y_label Character string for y-axis label. Default is empty.
#' @param max Numeric value specifying the maximum coordinate value (image height).
#'   If NULL, calculated automatically from data. Default is NULL.
#' @param scale Character string specifying the scale for tissue coordinates ("hires" or "lowres").
#'   Default is "hires".
#' @param image_scale Character string specifying which image in the Seurat object to use.
#'   This should match a name in object@images. Default is "slice1.008um".
#' @param alpha Numeric value between 0 and 1 for point transparency. Default is 0.5.
#' @param x_range Numeric vector of length 2 specifying the x-axis range to display.
#'   If NULL, shows the entire range. Default is NULL.
#' @param y_range Numeric vector of length 2 specifying the y-axis range to display.
#'   If NULL, shows the entire range. Default is NULL.
#' @param show_image Logical indicating whether to display the background tissue image.
#'   Default is TRUE.
#' @param min_expression Numeric threshold for minimum feature expression. Cells with expression
#'   at or below this value will not be displayed. Default is 0.
#'
#' @return A ggplot object representing the spatial feature plot with only expressing cells.
#'
#' @details
#' This function visualizes the expression of a single gene/feature across a tissue section with:
#' 1. Background tissue image (if available and show_image=TRUE)
#' 2. Spots/cells positioned according to their spatial coordinates
#' 3. Only cells expressing the feature above min_expression are displayed
#' 4. Coloring based on the expression level of the specified feature
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
#' # Basic usage with default parameters for a specific gene
#' p1 <- stFeature(spatial_seurat_obj,
#'                feature = "TPSB2")
#'
#' # Only show cells with expression > 1
#' p2 <- stFeature(spatial_seurat_obj,
#'                feature = "TPSB2",
#'                min_expression = 1)
#'
#' # Customize color gradient
#' p3 <- stFeature(spatial_seurat_obj,
#'                feature = "TPSB2",
#'                color = c("grey90", "darkred"))
#'
#' # Customize point size and transparency
#' p4 <- stFeature(spatial_seurat_obj,
#'                feature = "TPSB2",
#'                point_size = 1.2,
#'                alpha = 0.8)
#'
#' # Add title and specify coordinate scale
#' p5 <- stFeature(spatial_seurat_obj,
#'                feature = "TPSB2",
#'                title = "TPSB2 Expression in Tissue",
#'                scale = "hires")
#'
#' # Zoom to a specific region
#' p6 <- stFeature(spatial_seurat_obj,
#'                feature = "TPSB2",
#'                x_range = c(2000, 4000),
#'                y_range = c(1000, 3000))
#'
#' # Hide background image to focus on expression pattern
#' p7 <- stFeature(spatial_seurat_obj,
#'                feature = "TPSB2",
#'                show_image = FALSE)
#'
#' # Display the plot
#' print(p1)
#'
#' # Save the plot with high resolution
#' ggsave("spatial_feature.png", p1, width = 10, height = 8, dpi = 300)
#'
#' @export
stFeature <- function(object, feature = "gene", color = c("blue", "red"), point_size = 0.5, title = "",
                      x_label = "", y_label = "", max = NULL, scale = "hires",
                      image_scale = "slice1.008um", alpha = 0.5, x_range = NULL, y_range = NULL,
                      show_image = TRUE, min_expression = 0) {
  require(ggplot2)
  require(dplyr)


  expression_matrix <- GetAssayData(object, slot = "data")


  if (!feature %in% rownames(expression_matrix)) {
    stop(paste("Feature", feature, "not found in expression matrix"))
  }


  feature_data <- expression_matrix[feature, ]


  coordinates <- GetTissueCoordinates(object, scale = scale)
  if (is.null(coordinates) || nrow(coordinates) == 0) {
    stop("坐标数据为空，无法生成图形")
  }

  Xmax <- if (is.null(max)) max(coordinates$y, na.rm = TRUE) else max


  plot_data <- data.frame(
    x = coordinates$y,
    y = Xmax - coordinates$x,
    cell = rownames(coordinates)
  )


  plot_data$expression <- feature_data[plot_data$cell]


  plot_data$has_expression <- !is.na(plot_data$expression) & plot_data$expression > min_expression


  xlim_valid <- !is.null(x_range) && length(x_range) == 2
  ylim_valid <- !is.null(y_range) && length(y_range) == 2


  if (xlim_valid) {
    plot_data <- plot_data[plot_data$x >= x_range[1] & plot_data$x <= x_range[2], ]
  }

  if (ylim_valid) {
    plot_data <- plot_data[plot_data$y >= y_range[1] & plot_data$y <= y_range[2], ]
  }


  if (is.null(x_range)) x_range <- range(plot_data$x, na.rm = TRUE)
  if (is.null(y_range)) y_range <- range(plot_data$y, na.rm = TRUE)


  plot_data_expr <- plot_data[plot_data$has_expression, ]


  p <- ggplot() +
    scale_color_gradient(low = color[1], high = color[2]) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      legend.position = "right",
      legend.direction = "vertical",
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 25),
      legend.margin = margin(t = -22, b = -22),
      plot.title = element_text(size = 22, hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(0.2, "cm")
    ) +
    guides(color = guide_colorbar(title = "Expression", barwidth = 1)) +
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


  if (nrow(plot_data_expr) > 0) {
    p <- p + geom_point(
      data = plot_data_expr,
      aes(x = x, y = y, color = expression),
      size = point_size,
      alpha = alpha
    )
  } else {
    warning("No cells express ", feature, " above the threshold of ", min_expression)
  }


  p <- p + coord_cartesian(
    xlim = if (xlim_valid) x_range else c(NA_real_, NA_real_),
    ylim = if (ylim_valid) y_range else c(NA_real_, NA_real_),
    expand = FALSE
  )

  return(p)
}
