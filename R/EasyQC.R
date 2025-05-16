#' EasyQC
#'
#' This function performs quality control on a Seurat object by calculating the percentages of mitochondrial, ribosomal, and hemoglobin genes and optionally plotting the QC metrics.
#'
#' @param seurat A Seurat object to be processed.
#' @param species Either 'human' or 'mouse', specifying which species the data comes from.
#' @param MT The mitochondrial gene filter threshold, default is 5.
#' @param RP The ribosomal gene filter threshold, default is 5.
#' @param HB The hemoglobin gene filter threshold, default is 1.
#' @param Feature_high The upper threshold for RNA feature count, default is 5000.
#' @param Feature_low The lower threshold for RNA feature count, default is 200.
#' @param Plot Logical value, whether to plot QC metrics (default TRUE).
#' @param color A vector of colors to be used in the plot, effective if Plot = TRUE. Default color set is provided.
#' @return A filtered Seurat object after quality control.
#' @export
EasyQC <- function(seurat, species = "human", MT = 5, RP = 5, HB = 1,
                     Feature_high = 5000, Feature_low = 200, Plot = TRUE,
                     color = c("#FF6F61", "#6B5B95", "#88B04B", "#98FB98", "#92A8D1",
                               "#955251", "#B565A7", "#009B77", "#F7CAC9", "#D65076",
                               "#45B8AC", "#EFC050", "#5B5EA6", "#9B2335", "#DFCFBE")) {

  if (missing(seurat)) {
    stop("You must provide a Seurat object.")
  }
  if (!(species %in% c("human", "mouse"))) {
    stop("Species must be either 'human' or 'mouse'.")
  }

  if (species == "human") {
    seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^MT-")
    seurat[["percent.rp"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^RP[sl]")
    seurat[["percent.hb"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^HB[^(p)]")
  } else if (species == "mouse") {
    seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^Mt-")
    seurat[["percent.rp"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^Rp[sl]")
    seurat[["percent.hb"]] <- Seurat::PercentageFeatureSet(seurat, pattern = "^Hb[^(p)]")
  }

  seurat <- subset(x = seurat,
                   subset = nFeature_RNA > Feature_low &
                     nFeature_RNA < Feature_high &
                     percent.mt < MT &
                     percent.rp < RP &
                     percent.hb < HB)

  if (Plot) {
    df <- Seurat::FetchData(seurat, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "percent.hb", "orig.ident"))
    df_long <- tidyr::pivot_longer(df, cols = c(nFeature_RNA, nCount_RNA, percent.mt, percent.rp, percent.hb),
                                   names_to = "indicator",
                                   values_to = "value")

    p <- ggpubr::ggviolin(df_long, x = "orig.ident", y = "value", fill = "indicator", palette = color,
                          add = "boxplot", add.params = list(fill = "white", error.plot = "linerange")) +
      facet_wrap(~indicator, scales = "fixed")

    print(p)
  }

  return(seurat)
}

# 示例用法：
# sce <- EasyQC(seurat = your_seurat_object, species = "human", Plot = TRUE)
