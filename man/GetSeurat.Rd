% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetSeurat.R
\name{GetSeurat}
\alias{GetSeurat}
\title{GetSeurat}
\usage{
GetSeurat(h5ad_path, debug = FALSE)
}
\arguments{
\item{h5ad_path}{A string specifying the file path to the H5AD file.}

\item{debug}{Logical. If TRUE, prints debugging information. Default is FALSE.}
}
\value{
A Seurat object
}
\description{
This function converts an AnnData object stored in an H5AD file into a Seurat object.
It supports loading raw counts (if available) or normalized counts and includes metadata
and dimensionality reduction embeddings.
}
