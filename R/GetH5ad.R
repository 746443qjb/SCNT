#' GetH5ad
#'
#' This function converts a Seurat object to an H5AD file format.
#' Convert a Seurat object into a .h5ad file readable by Scanpy. For spatial transcriptomics data (mode = "st"),
#' this function includes spatial coordinates and metadata required by sc.pl.spatial() in Scanpy.
#' **Note: The image file (e.g. tissue_lowres_image.png) is not exported.** You must save it manually into a
#' folder named spatial/ in the same directory as the .h5ad file and match the path accordingly in Scanpy.
#'
#' @param seurat_obj A Seurat object to be converted.
#' @param output_path A character string specifying the path where the H5AD file will be saved.
#' @param mode Export mode: "sc" (single-cell) or "st" (spatial transcriptomics).
#' @param assay Name of the assay to use for expression data. Default is "RNA".
#' @param coord_cols Character vector specifying coordinate column names in coords dataframe. Default c("x", "y").
#' @param cell_col Character specifying the column name containing cell identifiers in coords. Default "cell".
#' @param scale Scale to use for tissue coordinates, e.g., "hires" or "lowres". Default is "hires".
#' @param external_coords Optional external coordinates dataframe. If provided, will be used instead of extracting coordinates from Seurat object.
#' @param debug Logical. If TRUE, prints debugging information. Default is FALSE.
#'
#' @details This function requires a working Python environment with anndata, numpy, and pandas libraries installed.
#' It is recommended to run joinlayers() on Seurat v5 objects before using this function to ensure proper compatibility.
#' This function supports dynamic addition of dimensionality reduction embeddings stored in the Seurat object.
#'
#' @section Python Configuration:
#' Before using this function, ensure Python is properly configured.
#'
#' @return No return value. The H5AD file will be saved at the specified location.
#' @export
GetH5ad <- function(seurat_obj, output_path, mode = c("sc", "st"), assay = "RNA",
                    coord_cols = c("x", "y"), cell_col = "cell", scale = "hires",
                    external_coords = NULL, debug = FALSE) {
  mode <- match.arg(mode)

  # Check if Python is properly configured
  if (!reticulate::py_available(initialize=TRUE)) {
    stop("Python environment not detected. Please configure Python before using this function.")
  }

  # Check if required Python packages are available
  tryCatch({
    anndata <- reticulate::import("anndata", delay_load = FALSE)
    np <- reticulate::import("numpy", delay_load = FALSE)
    pd <- reticulate::import("pandas", delay_load = FALSE)
  }, error = function(e) {
    stop("Required Python packages not found. Please install the following packages:\n",
         "  reticulate::py_install(c('anndata', 'numpy', 'pandas'))\n",
         "  OR manually: pip install anndata numpy pandas\n",
         "Error details: ", e$message)
  })

  # Expression matrix
  counts_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = 'counts')
  if (!inherits(counts_matrix, "matrix") && !inherits(counts_matrix, "dgCMatrix")) {
    counts_matrix <- as.matrix(counts_matrix)
  }

  original_cell_ids <- colnames(counts_matrix)
  gene_names <- rownames(counts_matrix)

  if (debug) {
    message("Original matrix dimensions: ", nrow(counts_matrix), " genes x ", ncol(counts_matrix), " cells")
    message("First 5 cell IDs: ", paste(head(original_cell_ids, 5), collapse=", "))
    message("First 5 gene names: ", paste(head(gene_names, 5), collapse=", "))
  }


  counts_matrix <- Matrix::t(counts_matrix)
  counts_matrix <- as(counts_matrix, "dgCMatrix")

  if (debug) {
    message("Transposed matrix dimensions: ", nrow(counts_matrix), " cells x ", ncol(counts_matrix), " genes")
    message("Row names (cells): ", paste(head(rownames(counts_matrix), 5), collapse=", "))
    message("Column names (genes): ", paste(head(colnames(counts_matrix), 5), collapse=", "))
  }


  meta_data <- seurat_obj@meta.data
  meta_data$barcode <- rownames(meta_data)

  adata <- anndata$AnnData(X = counts_matrix, obs = pd$DataFrame(meta_data))
  adata$var_names <- np$array(gene_names)


  if (length(seurat_obj@reductions) > 0) {
    for (reduction_name in names(seurat_obj@reductions)) {
      embeddings <- Seurat::Embeddings(seurat_obj, reduction = reduction_name)
      obsm_key <- paste0("X_", reduction_name)
      adata$obsm[obsm_key] <- np$array(embeddings)
    }
  }


  if (mode == "st") {

    if (is.null(external_coords)) {
      tryCatch({
        coords <- GetTissueCoordinates(seurat_obj, scale = scale)
        if (debug) {
          message("Successfully retrieved coordinates using GetTissueCoordinates()")
          message("Coordinate dimensions: ", nrow(coords), " x ", ncol(coords))
          message("Coordinate column names: ", paste(colnames(coords), collapse=", "))
        }
      }, error = function(e) {
        stop("Failed to get coordinates from Seurat object: ", e$message,
             "\nConsider providing external_coords parameter.")
      })
    } else {
      coords <- external_coords
      if (debug) {
        message("Using provided external coordinates")
        message("Coordinate dimensions: ", nrow(coords), " x ", ncol(coords))
        message("Coordinate column names: ", paste(colnames(coords), collapse=", "))
      }
    }

    if (!is.null(cell_col) && cell_col %in% colnames(coords)) {
      spatial_coords <- coords[, coord_cols, drop=FALSE]
      rownames(spatial_coords) <- coords[[cell_col]]
      if (debug) {
        message("Using cell IDs from column: ", cell_col)
        message("First 5 cell IDs from coordinates: ", paste(head(coords[[cell_col]], 5), collapse=", "))
      }
    } else {
      if (debug) {
        message("No cell_col found in coords, using row names as cell IDs")
        message("First 5 row names from coordinates: ", paste(head(rownames(coords), 5), collapse=", "))
      }
      spatial_coords <- coords[, coord_cols, drop=FALSE]
    }

    adata$obsm$`__setitem__`("spatial", np$array(as.matrix(spatial_coords)))

      image_name <- names(seurat_obj@images)[1]
      image_obj <- seurat_obj@images[[image_name]]
      scale_factors <- image_obj@scale.factors

      adata$uns[["spatial"]] <- reticulate::r_to_py(list(
        sample = list(
          images = list(),
          scalefactors = scale_factors,
          metadata = list()
        )
      ))

  }


  adata$write_h5ad(output_path)
  message("H5AD file successfully created at: ", output_path)
}

# 示例用法：
# GetH5ad(sce, "sce_sc.h5ad", assay = "Spatial",mode="sc")
# GetH5ad(sce, "sce_st.h5ad", assay = "Spatial",mode="st",scale="hires")
