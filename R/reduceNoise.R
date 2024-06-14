
#' Reduce noise in single-cell data
#' 
#' Simple filter for removing noise in single-cell data.
#' 
#' @param cellmat Matrix of log2 mean gene expression in rows with cell types in
#'   columns.
#' @param noisefilter Sets an upper bound for `noisefraction` cut-off below
#'   which gene expression is set to 0. Essentially gene expression above this
#'   level must be retained in the signature. Setting this higher can allow more
#'   suppression via noisefraction and can favour more highly expressed genes.
#' @param noisefraction Numeric value. Maximum mean log2 gene expression across
#'   cell types is calculated and values in celltypes below this fraction are
#'   set to 0. Set in conjunction with `noisefilter.` Note: if this is set too
#'   high (too close to 1), it can have a deleterious effect on deconvolution.
#' @returns Filtered mean gene expression matrix with genes in rows and cell
#'   types in columns.
#' @export
#'
reduceNoise <- function(cellmat,
                        noisefilter = 2,
                        noisefraction = 0.25) {
  genemax <- rowMaxs(cellmat)
  genelim <- pmin(genemax * noisefraction, noisefilter)
  cellmat[cellmat < genelim] <- 0
  cellmat
}
