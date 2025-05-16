#' EasyDoublet
#'
#' This function identifies and removes doublets from a Seurat object using the DoubletFinder package.
#'
#' @param seurat A Seurat object to be processed.
#' @param PC The number of Principal Components (PCs) to use, which can be obtained from KidneyPC.
#' @param rate The rate of doublets per 1000 cells, default is 8.
#' @param select Either "high" or "low" to determine the method of doublet removal.
#' @return A filtered Seurat object without doublets.
#' @export
EasyDoublet <- function(seurat, PC, rate = 8, select = "high") {

  if (missing(seurat)) {
    stop("You must provide a Seurat object.")
  }
  if (missing(PC)) {
    stop("You must provide the number of PCs to use.")
  }
  if (!(select %in% c("high", "low"))) {
    stop("Select must be either 'high' or 'low'.")
  }

  sweep.res.list <- paramSweep(seurat, PCs = 1:PC, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), ]$pK))

  DoubletRate <- ncol(seurat) * rate * 1e-6
  nExp_poi <- round(DoubletRate * ncol(seurat))
  homotypic.prop <- modelHomotypic(seurat@meta.data$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  seurat <- doubletFinder(seurat, PCs = 1:PC, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(seurat@meta.data)[ncol(seurat@meta.data)] <- "doublet_low"
  seurat <- doubletFinder(seurat, PCs = 1:PC, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(seurat@meta.data)[ncol(seurat@meta.data)] <- "doublet_high"

  if (select == "high") {
    seurat <- subset(seurat, subset = doublet_high == "Singlet")
  } else if (select == "low") {
    seurat <- subset(seurat, subset = doublet_low == "Singlet")
  }

  return(seurat)
}

# 示例用法：
# seurat_filtered <- EasyDoublet(seurat = your_seurat_object, PC = pcs, rate = 8, select = "high")


