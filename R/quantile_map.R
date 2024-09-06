
#' Quantile mapping function between two scRNA-Seq datasets
#' 
#' Quantile mapping to combine two scRNA-Seq datasets based on mapping either
#' the distribution of mean log2+1 gene expression in cell clusters to the
#' distribution of the 2nd dataset, or mapping the quantiles of one matrix of
#' gene expression (with genes in rows) to another.
#' 
#' @param x scRNA-Seq data whose distribution is to be mapped onto `y`: either a
#'   matrix of gene expression on log2+1 scale, or a 'cellMarkers' class object,
#'   in which case the `$genemeans` list element is extracted.
#' @param y Reference scRNA-Seq data: either a matrix of gene expression on
#'   log2+1 scale, or a 'cellMarkers' class object, in which case the
#'   `$genemeans` list element is extracted.
#' @param n Number of quantiles to split x and y
#' @param remove_noncoding Logical, whether to remove noncoding. This is a basic
#'   filter which looks at the gene names (rownames) in both matrices and
#'   removes genes containing "-" which are usually antisense or mitochondrial
#'   genes, or "." which are either pseudogenes or ribosomal genes.
#' @returns A function of form `FUN(x)` where `x` can be supplied as a numeric
#'   vector or matrix and the same type is returned.
#' @importFrom stats quantile
#' @export

quantile_map <- function(x, y, n = 1e4, remove_noncoding = TRUE) {
  if (inherits(x, "cellMarkers")) x <- x$genemeans
  if (inherits(y, "cellMarkers")) y <- y$genemeans
  common <- intersect(rownames(x), rownames(y))
  if (remove_noncoding) common <- common[!grepl("-|\\.", common)]
  message(length(common), " common genes")
  x <- as.vector(x[common, ])
  y <- as.vector(y[common, ])
  x <- x[x != 0]
  y <- y[y != 0]
  qx <- quantile(x, seq(0, 1, 1/n))
  qy <- quantile(y, seq(0, 1, 1/n))
  FUN <- approxfun(qx, qy, yleft = 0, rule = 2) |> suppressWarnings()
  approxfun.matrix <- function(x) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) {
      out <- FUN(as.vector(x))
      out <- matrix(out, nrow = nrow(x), dimnames = dimnames(x))
      return(out)
    }
    FUN(x)
  }
  approxfun.matrix
}
