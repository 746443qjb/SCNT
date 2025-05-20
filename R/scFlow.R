#' Create Flow-like Plots for Gene Expression in Single-Cell/Spatial Data
#'
#' This function generates a series of flow-like plots to visualize the relationship between two genes
#' across different cell groups. Each plot shows a 2D density distribution with quadrants defined by
#' cutoff values, and displays the percentage of cells in each quadrant.
#'
#' @param obj A Seurat object containing gene expression data.
#' @param x_gene Character string specifying the gene name for the x-axis.
#' @param y_gene Character string specifying the gene name for the y-axis.
#' @param ref Optional character string specifying a reference gene. If provided, the function plots
#'   x_gene/ref and y_gene/ref. If NULL, raw expression values are used. Default is NULL.
#' @param group Character string specifying which metadata column to use for grouping cells.
#' @param x_cut Numeric value for the vertical cutoff line position.
#' @param y_cut Numeric value for the horizontal cutoff line position.
#' @param assays Character string specifying which assay to use for gene expression data. Default is "RNA".
#' @param color Vector of colors defining the density gradient. Default is c("white", "yellow", "orange", "red").
#' @param nrow Integer specifying the number of rows in the multi-plot layout. Default is 2.
#' @param ncol Integer specifying the number of columns in the multi-plot layout. Default is 3.
#' @param x_label_size Numeric value for the size of x-axis label. Default is 12.
#' @param y_label_size Numeric value for the size of y-axis label. Default is 12.
#' @param title_size Numeric value for the size of plot titles. Default is 14.
#' @param legend_size Numeric value for the size of legend text. Default is 10.
#' @param axis_text_size Numeric value for the size of axis text. Default is 10.
#' @param common_legend Logical indicating whether to use a common legend for all plots. Default is FALSE.
#' @param cutline_color Character string specifying the color for cutoff lines. Default is "red".
#' @param cutline_size Numeric value for the thickness of cutoff lines. Default is 0.8.
#' @param point_size Numeric value for the size of points. Default is 0.5.
#' @param point_color Character string specifying the color of points. Default is "black".
#' @param point_alpha Numeric value between 0 and 1 for point transparency. Default is 0.6.
#' @param theme Character string specifying the theme style ("default" or "dark"). Default is "default".
#' @param quad_text_size Numeric value for the size of quadrant percentage text. Default is 4.
#' @param quad_text_color Character string specifying the color of quadrant percentage text. Default is "black".
#' @param x_range Numeric vector of length 2 specifying the x-axis range. Default is NULL (auto-determined).
#' @param y_range Numeric vector of length 2 specifying the y-axis range. Default is NULL (auto-determined).
#'
#' @return A ggpubr arranged plot containing multiple flow plots, one for each cell group.
#'
#' @details
#' This function creates a series of flow-like plots with:
#' 1. 2D density visualization of gene expression or gene expression ratios
#' 2. Quadrants defined by user-specified cutoffs
#' 3. Percentage of cells in each quadrant
#'
#' When a reference gene is provided (`ref` parameter), the function normalizes expression values
#' of x_gene and y_gene by dividing by the reference gene expression. This is useful for visualizing
#' relative expression levels or gene ratios.
#'
#' The plots are arranged in a grid according to the `nrow` and `ncol` parameters, with one plot
#' for each unique value in the specified `group` column of the Seurat object's metadata.
#'
#' @examples
#' # Load required packages
#' library(Seurat)
#' library(dplyr)
#' library(ggplot2)
#' library(ggpubr)
#'
#' # Basic usage - comparing expression of two genes across cell clusters
#' p1 <- scflow(seurat_obj,
#'             x_gene = "CD4",
#'             y_gene = "CD8A",
#'             group = "seurat_clusters",
#'             x_cut = 1.0,
#'             y_cut = 1.0)
#'
#' # Using a reference gene to normalize expression (gene ratios)
#' p2 <- scflow(seurat_obj,
#'             x_gene = "MKI67",
#'             y_gene = "TOP2A",
#'             ref = "PCNA",     # Both genes will be normalized by PCNA expression
#'             group = "celltype",
#'             x_cut = 0.5,
#'             y_cut = 0.5)
#'
#' # Customizing appearance with dark theme and specific axis ranges
#' p3 <- scflow(seurat_obj,
#'             x_gene = "IL6",
#'             y_gene = "TNF",
#'             group = "treatment",
#'             x_cut = 1.5,
#'             y_cut = 1.0,
#'             theme = "dark",
#'             color = c("darkblue", "blue", "cyan", "yellow"),
#'             x_range = c(0, 5),
#'             y_range = c(0, 4))
#'
#' # Using a non-standard assay (e.g., for multi-modal data)
#' p4 <- scflow(seurat_obj,
#'             x_gene = "CD4",
#'             y_gene = "CD8A",
#'             group = "seurat_clusters",
#'             assays = "ADT",    # Using protein expression data instead of RNA
#'             x_cut = 2.0,
#'             y_cut = 2.0)
#'
#' # Display the plot
#' print(p1)
#'
#' # Save the plot with high resolution
#' ggsave("flow_plot.png", p1, width = 15, height = 10, dpi = 300)
#'
#' @export
scflow <- function(
    obj, x_gene, y_gene, ref = NULL,
    group, x_cut, y_cut, assays = "RNA",
    color = c("white", "yellow", "orange", "red"), nrow = 2, ncol = 3,
    x_label_size = 12, y_label_size = 12, title_size = 14, legend_size = 10, axis_text_size = 10,
    common_legend = FALSE,
    cutline_color = "red", cutline_size = 0.8,
    point_size = 0.5, point_color = "black", point_alpha = 0.6,
    theme = "default",
    quad_text_size = 4, quad_text_color = "black",
    x_range = NULL, y_range = NULL,
    skip_density = FALSE
) {

  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(tidyr)


  genes_to_check <- c(x_gene, y_gene)
  if (!is.null(ref)) {
    genes_to_check <- c(genes_to_check, ref)
  }


  if (!assays %in% names(obj@assays)) {
    stop(paste("Assay", assays, "not found in the Seurat object."))
  }


  if (!all(genes_to_check %in% rownames(obj@assays[[assays]]$data))) {
    stop("One or more specified genes not found in the data matrix.")
  }


  if (!group %in% colnames(obj@meta.data)) {
    stop(paste("Group column", group, "not found in meta.data"))
  }

  expr_data <- t(as.matrix(obj@assays[[assays]]$data[genes_to_check, ]))
  expr_data <- as.data.frame(expr_data) %>% rownames_to_column(var = "cellName")


  if (!is.null(ref)) {

    min_non_zero <- min(expr_data[[ref]][expr_data[[ref]] > 0], na.rm = TRUE)
    ref_values <- ifelse(expr_data[[ref]] <= 0, min_non_zero * 0.1, expr_data[[ref]])

    expr_data[[x_gene]] <- expr_data[[x_gene]] / ref_values
    expr_data[[y_gene]] <- expr_data[[y_gene]] / ref_values
    expr_data <- expr_data[, c("cellName", x_gene, y_gene)]
  }


  geneData <- expr_data
  geneData$group_label <- obj@meta.data[[group]][match(geneData$cellName, rownames(obj@meta.data))]
  geneData$group_label <- factor(geneData$group_label)


  group_order <- levels(geneData$group_label)


  theme_style <- if (theme == "dark") {
    theme(
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "white", fill = NA, linewidth = 1),
      axis.text = element_text(color = "white", size = axis_text_size),
      axis.title = element_text(color = "white", size = x_label_size),
      plot.title = element_text(color = "white", size = title_size, face = "bold", hjust = 0.5),
      legend.text = element_text(color = "white", size = legend_size),
      legend.title = element_text(color = "white", size = legend_size),
      legend.background = element_rect(fill = "black", color = "white"),
      legend.key = element_rect(fill = "black", color = NA)
    )
  } else {
    theme_minimal() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = x_label_size),
        plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
        legend.text = element_text(size = legend_size),
        legend.title = element_text(size = legend_size)
      )
  }

  custom_plot <- function(plot.data, x_gene, y_gene, x.cut, y.cut, title, x_range, y_range) {

    plot.data <- plot.data %>%
      filter(!is.na(.data[[x_gene]]) & !is.na(.data[[y_gene]]) &
               is.finite(.data[[x_gene]]) & is.finite(.data[[y_gene]]))

    if (nrow(plot.data) < 3) {
      warning("Insufficient valid data points for plotting density.")
      p <- ggplot(plot.data, aes(x = .data[[x_gene]], y = .data[[y_gene]])) +
        geom_point(size = point_size, color = point_color, alpha = point_alpha) +
        geom_vline(xintercept = x.cut, linetype = "dashed", color = cutline_color, linewidth = cutline_size) +
        geom_hline(yintercept = y.cut, linetype = "dashed", color = cutline_color, linewidth = cutline_size) +
        labs(title = title, x = x_gene, y = y_gene) +
        theme_style
      return(p)
    }

    if (is.null(x_range)) {
      x_range <- range(plot.data[[x_gene]], na.rm = TRUE)

      if (diff(x_range) < 1e-6) {
        x_center <- mean(x_range)
        x_range <- c(x_center - 0.5, x_center + 0.5)
      } else {
        x_buffer <- diff(x_range) * 0.05
        x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
      }
    }

    if (is.null(y_range)) {
      y_range <- range(plot.data[[y_gene]], na.rm = TRUE)

      if (diff(y_range) < 1e-6) {
        y_center <- mean(y_range)
        y_range <- c(y_center - 0.5, y_center + 0.5)
      } else {
        y_buffer <- diff(y_range) * 0.05
        y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
      }
    }


    plot.data <- plot.data %>%
      mutate(
        quadrant = case_when(
          .data[[x_gene]] < x.cut & .data[[y_gene]] < y.cut ~ "Q1",
          .data[[x_gene]] >= x.cut & .data[[y_gene]] < y.cut ~ "Q2",
          .data[[x_gene]] < x.cut & .data[[y_gene]] >= y.cut ~ "Q3",
          TRUE ~ "Q4"
        )
      )

    quadrant_counts <- plot.data %>%
      group_by(quadrant) %>%
      tally(name = "count") %>%
      mutate(percent = count / sum(count) * 100)


    p <- ggplot(plot.data, aes(x = .data[[x_gene]], y = .data[[y_gene]]))


    if (!skip_density) {
      tryCatch({
        p <- p +
          stat_density_2d(
            aes(fill = after_stat(level)),
            geom = "polygon",
            contour = TRUE,
            h = c(diff(x_range)/15, diff(y_range)/15),
            bins = 10
          ) +
          scale_fill_gradientn(colors = color, name = "Density")
      }, error = function(e) {
        warning(paste("Could not draw density contours:", e$message))
      })
    }

    p <- p +
      geom_point(size = point_size, color = point_color, alpha = point_alpha) +
      geom_vline(xintercept = x.cut, linetype = "dashed", color = cutline_color, linewidth = cutline_size) +
      geom_hline(yintercept = y.cut, linetype = "dashed", color = cutline_color, linewidth = cutline_size) +
      labs(
        title = title,
        x = x_gene,
        y = y_gene
      ) +
      xlim(x_range) +
      ylim(y_range) +
      theme_style

    for (quad in quadrant_counts$quadrant) {
      quad_data <- quadrant_counts %>% filter(quadrant == quad)
      x_pos <- ifelse(grepl("Q1|Q3", quad),
                      mean(x_range) - diff(x_range) * 0.3,
                      mean(x_range) + diff(x_range) * 0.3)
      y_pos <- ifelse(grepl("Q1|Q2", quad),
                      mean(y_range) - diff(y_range) * 0.3,
                      mean(y_range) + diff(y_range) * 0.3)

      p <- p +
        annotate("text", x = x_pos, y = y_pos,
                 label = paste0(round(quad_data$percent, 2), "%"),
                 color = quad_text_color, size = quad_text_size, fontface = "bold")
    }

    return(p)
  }

  gg.ls <- list()
  for (group in group_order) {
    current_data <- geneData %>% dplyr::filter(group_label == group)

    if (nrow(current_data) < 2) {
      gg.ls[[group]] <- ggplot() +
        labs(title = paste("No data or insufficient data for", group)) +
        theme_void()
      next
    }

    gg.ls[[group]] <- tryCatch({
      custom_plot(
        plot.data = current_data,
        x_gene = x_gene,
        y_gene = y_gene,
        x.cut = x_cut,
        y.cut = y_cut,
        title = group,
        x_range = x_range,
        y_range = y_range
      )
    }, error = function(e) {
      warning(paste("Error plotting group", group, ":", e$message))
      ggplot() +
        labs(title = paste("Error plotting", group, ":", e$message)) +
        theme_void()
    })
  }

  multi.page <- ggpubr::ggarrange(
    plotlist = gg.ls,
    nrow = nrow,
    ncol = ncol,
    common.legend = common_legend,
    legend = "right"
  )

  return(multi.page)
}
