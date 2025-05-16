#' Create Dotplot Visualizations for Single-Cell Gene Expression
#'
#' This function generates a dotplot visualization that displays gene expression patterns across different
#' cell types. Each dot represents a gene-celltype combination, where the size indicates the percentage of
#' cells expressing the gene and the color intensity represents the average expression level.
#'
#' @param seurat_obj A Seurat object containing gene expression data.
#' @param markers Character vector of gene names to visualize.
#' @param celltype_var Character string specifying which metadata column contains cell type information.
#'   Default is "celltype".
#' @param assay_layer Character string specifying which expression data layer to use ("data", "counts", "scale.data").
#'   Default is "data".
#' @param color Color palette for expression values. Default is c("#FFFFCC","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0").
#' @param standardize Logical indicating whether to standardize expression values across cell types for each gene.
#'   Default is TRUE.
#' @param scale_range Numeric vector of length 2 specifying the range to scale standardized expression values to.
#'   Only used when standardize=TRUE. Default is c(-1, 1).
#'
#' @return A ggplot object representing the dotplot visualization.
#'
#' @details
#' This function calculates two key metrics for each gene-celltype combination:
#' 1. Percentage of cells expressing the gene (represented by dot size)
#' 2. Average expression level (represented by dot color)
#'
#' When standardize=TRUE, the expression values are log-transformed, standardized across cell types
#' for each gene, and then rescaled to the specified range. When standardize=FALSE, the raw expression
#' values from the assay are used directly.
#'
#' @examples
#' # Load required packages
#' library(Seurat)
#' library(dplyr)
#' library(ggplot2)
#' library(scales)
#'
#' # Basic usage with default parameters
#' genes_of_interest <- c("CD8A", "CD4", "MS4A1", "CD14", "FCGR3A", "NKG7")
#' p1 <- scDot(seurat_obj, markers = genes_of_interest)
#'
#' # Using a different metadata column for cell types
#' p2 <- scDot(seurat_obj,
#'             markers = genes_of_interest,
#'             celltype_var = "seurat_clusters")
#'
#' # Without standardization, using raw expression values
#' p3 <- scDot(seurat_obj,
#'             markers = genes_of_interest,
#'             standardize = FALSE)
#'
#' # Using custom color palette and scale range
#' custom_colors <- c("navy", "white", "firebrick")
#' p4 <- scDot(seurat_obj,
#'             markers = genes_of_interest,
#'             color = custom_colors,
#'             scale_range = c(-2, 2))
#'
#' # Display the plot
#' print(p1)
#'
#' # Save the plot with high resolution
#' ggsave("dotplot_markers.png", p1, width = 12, height = 8, dpi = 300)
#'
#' @export
scDot <- function(seurat_obj, markers, celltype_var = "celltype", assay_layer = "data",
                  color = c("#FFFFCC","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0"),
                  standardize = TRUE, scale_range = c(-1, 1)) {
  require(dplyr)
  require(ggplot2)
  require(scales)


  celltypes <- seurat_obj@meta.data[[celltype_var]]


  expression_matrix <- GetAssayData(seurat_obj, slot = assay_layer)


  valid_markers <- intersect(markers, rownames(expression_matrix))
  if (length(valid_markers) == 0) {
    stop("None of the provided markers are present in the expression matrix.")
  }

  dot_data <- data.frame(
    gene = character(),
    celltype = character(),
    avg_exp_raw = numeric(),
    pct_exp = numeric()
  )

  for (gene in valid_markers) {
    gene_exp <- as.numeric(expression_matrix[gene, ])
    for (celltype in unique(celltypes)) {
      cell_indices <- celltypes == celltype


      if (standardize) {

        avg_exp_raw <- mean(expm1(gene_exp[cell_indices]))
      } else {

        avg_exp_raw <- mean(gene_exp[cell_indices])
      }

      pct_exp <- mean(gene_exp[cell_indices] > 0) * 100

      dot_data <- rbind(dot_data, data.frame(
        gene = gene,
        celltype = celltype,
        avg_exp_raw = avg_exp_raw,
        pct_exp = pct_exp
      ))
    }
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

  dot_data$gene <- factor(dot_data$gene, levels = rev(valid_markers))

  ggplot(dot_data, aes(x = celltype, y = gene, size = pct_exp, fill = avg_exp_scaled)) +
    geom_point(shape = 21,
               color = "black",
               stroke = 0.5) +
    scale_fill_gradientn(colors = color) +
    scale_size_continuous(range = c(4, 15)) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 25, color = "black", angle = 45, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 25, color = "black"),
      axis.ticks.y = element_blank(),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18)
    ) +
    labs(
      x = "Cell Type",
      y = "Gene",
      size = "Fraction of cells (%)",
      fill = "Expression"
    )
}
