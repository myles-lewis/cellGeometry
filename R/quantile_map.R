
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
#' @param n Number of quantiles to split `x` and `y`.
#' @param remove_noncoding Logical, whether to remove noncoding. This is a basic
#'   filter which looks at the gene names (rownames) in both matrices and
#'   removes genes containing "-" which are usually antisense or mitochondrial
#'   genes, or "." which are either pseudogenes or ribosomal genes.
#' @param knots Vector of quantile points for knots for fitting natural splines.
#' @param use_splines Logical whether to use "splines" to fit natural splines.
#'   If set to `FALSE`, interpolation is limited to the original range of `x`,
#'   i.e. it will clip for values > `max(x)`.
#' @returns A list object of class 'qqmap' containing:
#' \item{quantiles}{Dataframe containing matching quantiles of `x` and `y`}
#' \item{map}{A function of form `FUN(x)` where `x` can be supplied as a numeric
#'   vector or matrix and the same type is returned. The function converts given
#'   data points to the distribution of `y`.}
#' @details
#' The conversion uses the function [approxfun()] which uses interpolation. It
#' is not designed to perform stepwise (exact) quantile transformation. Zeroes
#' are removed before mapping each dataset to quantiles.
#' 
#' @seealso [approxfun()]
#' @importFrom stats quantile predict
#' @importFrom splines ns
#' @export

quantile_map <- function(x, y, n = 1e4, remove_noncoding = TRUE,
                         knots = c(0.25, 0.75, 0.85, 0.95, 0.97, 0.99, 0.999),
                         use_splines = TRUE) {
  xlab <- deparse(substitute(x))
  ylab <- deparse(substitute(y))
  if (inherits(x, "cellMarkers")) x <- x$genemeans
  if (inherits(y, "cellMarkers")) y <- y$genemeans
  common <- intersect(rownames(x), rownames(y))
  if (remove_noncoding) common <- common[!grepl("-|\\.", common)]
  message(length(common), " common genes")
  x <- as.vector(x[common, ])
  y <- as.vector(y[common, ])
  qi <- log10(seq(1, 10, length.out = n))  # better spacing
  qx <- quantile(x, qi)
  qy <- quantile(y, qi)
  df <- data.frame(qx, qy)
  if (use_splines) {
    kn <- quantile(qx, knots)
    fit <- lm(qy ~ ns(qx, knots = kn), df)
    qx <- seq(0, max(qx) *2, length.out = 1000)
    qy <- predict(fit, data.frame(qx))
    qy[qy < 0] <- 0
  }
  FUN <- approxfun(qx, qy, yleft = 0, rule = 2) |> suppressWarnings()
  FUN0 <- function(x) {
    nzero <- x != 0
    x[nzero] <- FUN(x[nzero])
    x
  }
  approxfun.matrix <- function(x) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) {
      out <- FUN0(as.vector(x))
      out <- matrix(out, nrow = nrow(x), dimnames = dimnames(x))
      return(out)
    }
    FUN0(x)
  }
  structure(list(quantiles = df, map = approxfun.matrix,
                 xlab = xlab, ylab = ylab), class = "qqmap")
}


#' Quantile-quantile plot
#' 
#' Produces a QQ plot showing the conversion function from the first dataset to
#' the second.
#' @param x A 'qqmap' class object created by [quantile_map()].
#' @param ... Optional plotting parameters passed to [plot()].
#' @returns No return value. Produces a QQ plot using base graphics with a red
#'   line showing the conversion function.
#' @importFrom graphics lines
#' @export
plot.qqmap <- function(x, ...) {
  new.args <- list(...)
  args <- list(x = x$quantiles$qx, y = x$quantiles$qy, cex = 0.5,
               xlab = x$xlab, ylab = x$ylab)
  if (length(new.args)) args[names(new.args)] <- new.args
  do.call("plot", args)
  
  xr <- par("usr")[1:2]
  xr[1] <- max(c(0, xr[1]))
  px <- seq(xr[1], xr[2], length.out = 1000)
  lines(px, x$map(px), col = "red")
}
