#' Read and Create Seurat Object from Visium Spatial Transcriptomics Data
#'
#' This function reads 10X Genomics Visium spatial transcriptomics data,
#' optionally supports binned data, loads the image, and returns a Seurat object.
#'
#' @param data.dir A string specifying the path to the Visium data directory. Can also point to binned output folder.
#' @param image.type Image resolution to load. Choose from `"lowres"` or `"hires"`. Default is `"lowres"`.
#' @param filename Name of the HDF5 count matrix file. Default is `"filtered_feature_bc_matrix.h5"`.
#' @param assay Name of the assay to use in the Seurat object. Default is `"Spatial"`.
#' @param slice Name of the tissue slice. Used as the image key in the Seurat object. Default is `"slice1"`.
#' @param bin.size Optional. If given (e.g., `c(16, 8)`), will read from binned_outputs folder using specified bin size.
#' @param filter.matrix Logical. Whether to filter the matrix to spots within tissue. Passed to `Read10X_Image()`. Default is `TRUE`.
#' @param to.upper Logical. Whether to convert gene names to uppercase. Default is `FALSE`.
#' @param image Optional. A custom image object (e.g., from `Read10X_Image`). If provided, it overrides image loading from disk.
#' @param ... Additional arguments passed to `Read10X_h5()`.
#'
#' @return A Seurat object containing the spatial transcriptomics data, including associated image.
#'
#' @importFrom Seurat Read10X_h5 Read10X_Image CreateSeuratObject
#' @export
ReadST <- function(data.dir, image.type = c("lowres", "hires"), filename = "filtered_feature_bc_matrix.h5",
                   assay = "Spatial", slice = "slice1", bin.size = NULL, filter.matrix = TRUE,
                   to.upper = FALSE, image = NULL, ...) {

  image.type <- match.arg(image.type)
  image.name <- if (image.type == "lowres") "tissue_lowres_image.png" else "tissue_hires_image.png"

  if (length(x = data.dir) > 1) {
    data.dir <- data.dir[1]
    warning(paste0("`data.dir` 期望单一值，但收到多个值 - ",
                   "继续使用第一个：'", data.dir, "'."), immediate. = TRUE)
  }
  if (!file.exists(data.dir)) {
    stop(paste0("没有找到文件或目录：", "'", data.dir, "'"))
  }

  if (is.null(bin.size) & file.exists(paste0(data.dir, "/binned_outputs"))) {
    bin.size <- c(16, 8)
  }

  if (!is.null(bin.size)) {
    bin.size.pretty <- paste0(sprintf("%03d", bin.size), "um")
    data.dirs <- paste0(data.dir, "/binned_outputs/", "square_", bin.size.pretty)
    assay.names <- paste0(assay, ".", bin.size.pretty)
    slice.names <- paste0(slice, ".", bin.size.pretty)
  }
  else {
    data.dirs <- data.dir
    assay.names <- assay
    slice.names <- slice
  }

  counts.paths <- lapply(data.dirs, file.path, filename)
  counts.list <- lapply(counts.paths, Read10X_h5, ...)

  if (to.upper) {
    rownames(counts) <- lapply(rownames(counts), toupper)
  }

  if (is.null(image)) {
    image.list <- mapply(Read10X_Image, file.path(data.dirs, "spatial"), image.name = image.name,
                         assay = assay.names, slice = slice.names, MoreArgs = list(filter.matrix = filter.matrix))
  }
  else {
    image.list <- c(image)
  }

  if (length(image.list) != length(counts.list)) {
    stop(paste0("图像的数量与计数矩阵的数量不匹配。请确保每个空间数据集都有对应的图像。"))
  }

  object.list <- mapply(CreateSeuratObject, counts.list, assay = assay.names)

  object.list <- mapply(function(.object, .image, .assay, .slice) {
    .image <- .image[Cells(.object)]
    .object[[.slice]] <- .image
    return(.object)
  }, object.list, image.list, assay.names, slice.names)

  object <- merge(object.list[[1]], y = object.list[-1])
  return(object)
}

# 示例用法：
# Visum data： sce=ReadST(data.dir=data_dir, image.type="hires")
# Visum data： sceHD=ReadST(data.dir=data_dir, image.type="hires"， bin.size = c("2","8","16")
