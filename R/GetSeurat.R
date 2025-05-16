#' GetSeurat
#'
#' This function converts an AnnData object stored in an H5AD file into a Seurat object.
#' It supports loading raw counts (if available) or normalized counts and includes metadata
#' and dimensionality reduction embeddings.
#'
#' @param h5ad_path A string specifying the file path to the H5AD file.
#' @param debug Logical. If TRUE, prints debugging information. Default is FALSE.
#' @return A Seurat object
#'
#' @export
GetSeurat <- function(h5ad_path, debug = FALSE) {
  library(reticulate)
  library(Seurat)
  library(Matrix)

  anndata <- import("anndata", convert = FALSE)

  adata <- anndata$read_h5ad(h5ad_path)

  use_raw <- FALSE
  var_names <- NULL
  if (!py_is_null_xptr(adata$raw)) {
    tryCatch({
      if (!is.null(adata$raw$X)) {
        use_raw <- TRUE
        var_names <- as.character(reticulate::py_to_r(adata$raw$var$index$to_list()))
      }
    }, error = function(e) {
      use_raw <- FALSE
    })
  }

  if (use_raw) {
    message("Using raw counts matrix from adata$raw$X")
    counts <- adata$raw$X
  } else {
    message("Using normalized counts matrix from adata$X")
    counts <- adata$X
    var_names <- as.character(reticulate::py_to_r(adata$var$index$to_list()))
  }

  counts <- reticulate::py_to_r(counts)

  if (!(inherits(counts, "dgCMatrix") || inherits(counts, "dgRMatrix"))) {
    counts <- as(as.matrix(counts), "dgCMatrix")
  }

  counts <- Matrix::t(counts)

  obs_names <- as.character(reticulate::py_to_r(adata$obs$index$to_list()))

  if (length(var_names) != nrow(counts)) {
    stop("Number of gene names does not match the number of rows in counts")
  }
  if (length(obs_names) != ncol(counts)) {
    stop("Number of cell names does not match the number of columns in counts")
  }

  rownames(counts) <- var_names
  colnames(counts) <- obs_names

  meta_data <- reticulate::py_to_r(adata$obs)
  rownames(meta_data) <- obs_names

  seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta_data, assay = "RNA")


  obsm_dict <- py_call(adata$obsm$as_dict)

  for (key in names(obsm_dict)) {
    embedding <- obsm_dict[[key]]
    embedding<- py_to_r(embedding)
    rownames(embedding) <- colnames(seurat_obj)
    colnames(embedding) <- paste0(key, "_", 1:ncol(embedding))
    fixed_key <- gsub("^X_", "", key)
    fixed_key <- paste0(fixed_key, "")
    seurat_obj[[fixed_key]] <- CreateDimReducObject(embeddings = embedding,
                                                    key = fixed_key, assay = "RNA")
  }
  return(seurat_obj)
}
# 示例用法：
# scnew=GetSeurat("sce.h5ad")
