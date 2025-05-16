#' Create Multiple Dot Plots for Gene Expression Across Cell Groups and Conditions
#'
#' This function generates a faceted dotplot visualization showing gene expression patterns
#' across different cell groups and experimental conditions. Each dot represents a
#' gene-group-condition combination, where the size indicates the percentage of cells
#' expressing the gene and the color intensity represents the average expression level.
#'
#' @param seurat_obj A Seurat object containing gene expression data and metadata.
#' @param markers Character vector of gene names to visualize.
#' @param group_var Character string specifying which metadata column contains cell group/type
#'   information (e.g., "celltype", "cluster"). Default is "celltype".
#' @param condition_var Character string specifying which metadata column contains condition
#'   information (e.g., "treatment", "timepoint"). Default is "condition".
#' @param assay_layer Character string specifying which expression data layer to use
#'   ("data", "counts", "scale.data"). Default is "data".
#' @param color_gradient Vector of colors defining the gradient for expression values.
#'   Default is c("blue", "yellow").
#' @param size_range Numeric vector of length 2 specifying the min and max dot sizes.
#'   Default is c(4, 15).
#' @param standardize Logical indicating whether to standardize expression values across
#'   groups for each gene. Default is TRUE.
#' @param scale_range Numeric vector of length 2 specifying the range to scale standardized
#'   expression values to. Only used when standardize=TRUE. Default is c(-2.5, 2.5).
#'
#' @return A ggplot object representing the faceted dotplot visualization.
#'
#' @details
#' This function creates a faceted dotplot with:
#' 1. Genes on the y-axis
#' 2. Conditions on the x-axis
#' 3. Cell groups as facets (columns)
#'
#' For each gene-group-condition combination, the function calculates:
#' - Percentage of cells expressing the gene (represented by dot size)
#' - Average expression level (represented by dot color)
#'
#' When standardize=TRUE, expression values are log-transformed, standardized across
#' groups for each gene, and then rescaled to the specified range. When standardize=FALSE,
#' the raw expression values are used directly.
#'
#' @examples
#' # Load required packages
#' library(Seurat)
#' library(dplyr)
#' library(ggplot2)
#' library(scales)
#'
#' # Basic usage with default parameters (standardized expression)
#' genes_of_interest <- c("IL6", "TNF", "IFNG", "IL1B", "IL10")
#' p1 <- scMultipleDot(seurat_obj,
#'                    markers = genes_of_interest,
#'                    group_var = "celltype",
#'                    condition_var = "treatment")
#'
#' # Using raw expression values (no standardization)
#' p2 <- scMultipleDot(seurat_obj,
#'                    markers = genes_of_interest,
#'                    group_var = "celltype",
#'                    condition_var = "treatment",
#'                    standardize = FALSE)
#'
#' # Custom color gradient and dot sizes
#' p3 <- scMultipleDot(seurat_obj,
#'                    markers = genes_of_interest,
#'                    color_gradient = c("navy", "white", "firebrick"),
#'                    size_range = c(2, 12))
#'
#' # Different grouping variables
#' p4 <- scMultipleDot(seurat_obj,
#'                    markers = genes_of_interest,
#'                    group_var = "seurat_clusters",
#'                    condition_var = "timepoint")
#'
#' # Display the plot
#' print(p1)
#'
#' # Save the plot with high resolution
#' ggsave("multiple_dotplot.png", p1, width = 15, height = 8, dpi = 300)
#'
#' @export
scMultipleDot <- function(
    seurat_obj,
    markers,
    group_var = "celltype",
    condition_var = "condition",
    assay_layer = "data",
    color_gradient = c("blue", "yellow"),
    size_range = c(4, 15),
    standardize = TRUE,
    scale_range = c(-2.5, 2.5)
) {
  # Load required packages
  require(dplyr)
  require(ggplot2)
  require(scales)

  if (!all(c(group_var, condition_var) %in% colnames(seurat_obj@meta.data))) {
    stop("Seurat 对象的 meta.data 中缺少 group_var 或 condition_var")
  }

  expression_matrix <- GetAssayData(seurat_obj, slot = assay_layer)
  group_info <- as.character(seurat_obj@meta.data[[group_var]])
  condition_info <- as.character(seurat_obj@meta.data[[condition_var]])

  group_levels <- unique(group_info)
  condition_levels <- unique(condition_info)

  dot_data <- data.frame(
    gene = character(),
    group = character(),
    condition = character(),
    avg_exp_raw = numeric(),
    pct_exp = numeric()
  )

  for (gene in markers) {
    if (!(gene %in% rownames(expression_matrix))) {
      print(paste("Gene not found:", gene))
      next
    }

    gene_exp <- expression_matrix[gene, ]

    for (g in group_levels) {
      for (c in condition_levels) {

        cells <- which(group_info == g & condition_info == c)

        if (length(cells) == 0) next


        avg_exp <- mean(gene_exp[cells], na.rm = TRUE)
        pct_exp <- mean(gene_exp[cells] > 0, na.rm = TRUE) * 100


        dot_data <- rbind(dot_data, data.frame(
          gene = gene,
          group = g,
          condition = c,
          avg_exp_raw = avg_exp,
          pct_exp = pct_exp
        ))
      }
    }
  }


  if (nrow(dot_data) == 0) {
    stop("dot_data is empty, check if genes exist in the expression matrix.")
  }


  if (standardize) {

    dot_data <- dot_data %>%
      group_by(gene) %>%
      mutate(

        avg_exp_scaled = scale(log1p(avg_exp_raw)),

        avg_exp_scaled = scales::rescale(
          avg_exp_scaled,
          to = scale_range,
          from = c(min(avg_exp_scaled), max(avg_exp_scaled))
        )
      ) %>%
      ungroup()
  } else {

    dot_data$avg_exp_scaled <- dot_data$avg_exp_raw
  }


  dot_data$gene <- factor(dot_data$gene, levels = markers)
  dot_data$group <- factor(dot_data$group, levels = group_levels)


  p <- ggplot(dot_data, aes(
    x = condition,
    y = gene,
    size = pct_exp,
    fill = avg_exp_scaled
  )) +
    geom_point(shape = 21, color = "black") +
    scale_size_continuous(range = size_range) +
    scale_fill_gradientn(colors = color_gradient) +
    scale_y_discrete(limits = markers) +
    facet_grid(cols = vars(group), scales = "fixed", space = "free") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 20, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 20, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 20, color = "black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18)
    ) +
    labs(
      x = "",
      y = "",
      fill = ifelse(standardize, "Scaled Expression", "Expression"),
      size = "Fraction of cells (%)"
    )

  return(p)
}
