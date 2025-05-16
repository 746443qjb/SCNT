#' GetPC
#'
#' This function determines the optimal number of Principal Components (PCs) to use in downstream analysis of a Seurat object.
#'
#' @param seurat A Seurat object that has been normalized, scaled, and processed with PCA or other dimensionality reduction methods.
#' @param reduction The dimensionality reduction method to use (e.g., "pca", "harmony"). Must exist in the Seurat object.
#' @param cum The cumulative percentage threshold to use for selecting PCs.
#' @param var The minimum variation percentage for selecting PCs.
#' @return The number of PCs that meet the given criteria.
#' @export
GetPC <- function(seurat, reduction = "pca", cum = 90, var = 5) {

  if (missing(seurat)) {
    stop("You must provide a Seurat object.")
  }
  if (!reduction %in% names(seurat@reductions)) {
    stop("The specified reduction method does not exist in the Seurat object.")
  }

  if (reduction == "pca") {
    stdevs <- seurat[[reduction]]@stdev
  } else {
    embeddings <- seurat[[reduction]]@cell.embeddings
    stdevs <- apply(embeddings, 2, sd)
  }

  if (is.null(stdevs) || length(stdevs) == 0) {
    stop("The specified reduction method does not contain valid standard deviations.")
  }
  pct <- stdevs / sum(stdevs) * 100
  cumu <- cumsum(pct)

  co1 <- which(cumu > cum & pct < var)[1]

  if (length(pct) > 1) {
    co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  } else {
    co2 <- NA
  }

  pcs <- min(co1, co2, na.rm = TRUE)

  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

  p <- ggplot(plot_df, aes(x = cumu, y = pct, label = rank, color = rank > pcs)) +
    geom_text(check_overlap = TRUE, size = 3.5, fontface = "bold") +
    geom_vline(xintercept = cum, color = "grey", linetype = "dashed") +
    geom_hline(yintercept = min(pct[pct > var], na.rm = TRUE), color = "grey", linetype = "dashed") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "skyblue")) +
    labs(
      title = "Optimal Principal Components Selection",
      x = "Cumulative Percentage",
      y = "Percentage of Variation"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  print(p)

  return(pcs)
}

# 示例用法：
# pcs <- GetPC(seurat = your_seurat_object, reduction = "pca", cum = 90, var = 5)

