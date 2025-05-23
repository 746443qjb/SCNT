% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stFeature.R
\name{stFeature}
\alias{stFeature}
\title{Visualize Gene/Feature Expression in Spatial Transcriptomics Data}
\usage{
stFeature(
  object,
  feature = "gene",
  color = c("blue", "red"),
  point_size = 0.5,
  title = "",
  x_label = "",
  y_label = "",
  max = NULL,
  scale = "hires",
  image_scale = "slice1.008um",
  alpha = 0.5,
  x_range = NULL,
  y_range = NULL,
  show_image = TRUE,
  min_expression = 0
)
}
\arguments{
\item{object}{A Seurat object containing spatial transcriptomics data with tissue coordinates
and image information.}

\item{feature}{Character string specifying the gene or feature name to visualize.
Default is "gene".}

\item{color}{Vector of two colors defining the low-to-high expression gradient.
Default is c("blue", "red").}

\item{point_size}{Numeric value for the size of points in the plot. Default is 0.5.}

\item{title}{Character string for the plot title. Default is empty.}

\item{x_label}{Character string for x-axis label. Default is empty.}

\item{y_label}{Character string for y-axis label. Default is empty.}

\item{max}{Numeric value specifying the maximum coordinate value (image height).
If NULL, calculated automatically from data. Default is NULL.}

\item{scale}{Character string specifying the scale for tissue coordinates ("hires" or "lowres").
Default is "hires".}

\item{image_scale}{Character string specifying which image in the Seurat object to use.
This should match a name in object@images. Default is "slice1.008um".}

\item{alpha}{Numeric value between 0 and 1 for point transparency. Default is 0.5.}

\item{x_range}{Numeric vector of length 2 specifying the x-axis range to display.
If NULL, shows the entire range. Default is NULL.}

\item{y_range}{Numeric vector of length 2 specifying the y-axis range to display.
If NULL, shows the entire range. Default is NULL.}

\item{show_image}{Logical indicating whether to display the background tissue image.
Default is TRUE.}

\item{min_expression}{Numeric threshold for minimum feature expression. Cells with expression
at or below this value will not be displayed. Default is 0.}
}
\value{
A ggplot object representing the spatial feature plot with only expressing cells.
}
\description{
This function creates a spatial plot visualizing the expression of a single gene or feature
across the tissue section. Only cells expressing the feature above a minimum threshold
will be displayed, making it easier to visualize expression patterns without background noise.
The expression level is represented by a color gradient, and the plot can be overlaid
on the tissue histology image.
}
\details{
This function visualizes the expression of a single gene/feature across a tissue section with:
\enumerate{
\item Background tissue image (if available and show_image=TRUE)
\item Spots/cells positioned according to their spatial coordinates
\item Only cells expressing the feature above min_expression are displayed
\item Coloring based on the expression level of the specified feature
}

The coordinates are transformed to match image orientation (flipped y-axis). The function
also supports cropping to specific regions using x_range and y_range parameters.

The 'max' parameter is important for coordinate transformation and should match the image height.
Users can determine this by examining their image dimensions.
}
\examples{
# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)

# Basic usage with default parameters for a specific gene
p1 <- stFeature(spatial_seurat_obj,
               feature = "TPSB2")

# Only show cells with expression > 1
p2 <- stFeature(spatial_seurat_obj,
               feature = "TPSB2",
               min_expression = 1)

# Customize color gradient
p3 <- stFeature(spatial_seurat_obj,
               feature = "TPSB2",
               color = c("grey90", "darkred"))

# Customize point size and transparency
p4 <- stFeature(spatial_seurat_obj,
               feature = "TPSB2",
               point_size = 1.2,
               alpha = 0.8)

# Add title and specify coordinate scale
p5 <- stFeature(spatial_seurat_obj,
               feature = "TPSB2",
               title = "TPSB2 Expression in Tissue",
               scale = "hires")

# Zoom to a specific region
p6 <- stFeature(spatial_seurat_obj,
               feature = "TPSB2",
               x_range = c(2000, 4000),
               y_range = c(1000, 3000))

# Hide background image to focus on expression pattern
p7 <- stFeature(spatial_seurat_obj,
               feature = "TPSB2",
               show_image = FALSE)

# Display the plot
print(p1)

# Save the plot with high resolution
ggsave("spatial_feature.png", p1, width = 10, height = 8, dpi = 300)

}
