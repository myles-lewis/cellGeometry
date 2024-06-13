
#' Reduce noise in single-cell data
#' 
#' Simple filter for removing noise in single-cell data.
#' 
#' @param cellmat Matrix of log2 mean gene expression in rows with cell types in
#'   columns.
#' @param noisefilter Numeric value below which mean log2 gene expression is reduced to 0.
#' @param noisefraction Numeric value. Maximum mean log2 gene expression across
#'   cell types is calculated and values in celltypes below this fraction are
#'   set to 0.
#' @returns Filtered mean gene expression matrix with genes in rows and cell
#'   types in columns.
#' @export
#'
reduceNoise <- function(cellmat,
                        noisefilter = 0.02,
                        noisefraction = 0.1) {
  genemax <- rowMaxs(cellmat)
  genelim <- pmin(genemax * noisefraction, noisefilter)
  cellmat[cellmat < genelim] <- 0
  cellmat
}
