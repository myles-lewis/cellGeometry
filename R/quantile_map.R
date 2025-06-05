
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
#' @param remove_zeros Logical, whether to remove zeros from both datasets.
#'   This shifts the quantile relationships.
#' @param smooth Either "loess" or "lowess" which apply [loess()] or [lowess()]
#'   to smooth the QQ fitted line, or "ns" which uses natural splines via
#'   [ns()]. With any other value no smoothing is applied. With no smoothing or
#'   "loess/lowess", interpolation is limited to the original range of `x`, i.e.
#'   it will clip for values > `max(x)`.
#' @param span controls the degree of smoothing in [loess()] and [lowess()].
#' @param knots Vector of quantile points for knots for fitting natural splines.
#' @param respace Logical whether to respace quantile points so their x axis
#'   density is more even. Can help spline fitting.
#' @param silent Logical whether to suppress messages.
#' @returns A list object of class 'qqmap' containing:
#' \item{quantiles}{Dataframe containing matching quantiles of `x` and `y`}
#' \item{map}{A function of form `FUN(x)` where `x` can be supplied as a numeric
#'   vector or matrix and the same type is returned. The function converts given
#'   data points to the distribution of `y`.}
#' @details
#' The conversion uses the function [approxfun()] which uses interpolation. It
#' is not designed to perform stepwise (exact) quantile transformation of every
#' individual datapoint.
#' 
#' @seealso [approxfun()]
#' @importFrom stats quantile predict loess lowess
#' @importFrom splines ns
#' @export

quantile_map <- function(x, y, n = 1e4, remove_noncoding = TRUE,
                         remove_zeros = FALSE,
                         smooth = "loess",
                         span = 0.15,
                         knots = c(0.25, 0.75, 0.85, 0.95, 0.97, 0.99, 0.999),
                         respace = FALSE, silent = FALSE) {
  xlab <- deparse(substitute(x))
  ylab <- deparse(substitute(y))
  if (inherits(x, "cellMarkers")) x <- x$genemeans
  if (inherits(y, "cellMarkers")) y <- y$genemeans
  if (inherits(x, "data.frame")) x <- as.matrix(x)
  if (inherits(y, "data.frame")) y <- as.matrix(y)
  common <- intersect(rownames(x), rownames(y))
  if (remove_noncoding) common <- common[!grepl("-|\\.", common)]
  if (!silent) cat(length(common), "common genes\n")
  x <- x[common, ]
  y <- y[common, ]
  if (remove_zeros) {
    x <- as.vector(x)
    y <- as.vector(y)
    x <- x[x != 0]
    y <- y[y != 0]
  }
  qi <- log10(seq(1, 10, length.out = n))  # better spacing
  qx <- quantile(x, qi)
  qy <- quantile(y, qi)
  df <- data.frame(qx, qy)
  if (respace) {
    xseq <- seq(qx[1], qx[length(qx)], length.out = 1000)
    ind <- vapply(xseq, function(xi) {
      which(qx >= xi)[1]
    }, integer(1L))
    ind <- unique(ind)
    qx <- qx[ind]
    qy <- qy[ind]
    df <- data.frame(qx, qy)
  }
  if (smooth == "loess") {
    fit <- loess(qy ~ qx, span = span) |> suppressWarnings()
    qx <- seq(0, qx[length(qx)], length.out = 1000)
    qy <- try(predict(fit, data.frame(qx)), silent = TRUE)
    if (!inherits(qy, "try-error")) {
      qy[qy < 0] <- 0
    } else {
      # loess error
      cat("Error in loess. Switching to lowess.\n")
      smooth <- "lowess"
    }
  }
  if (smooth == "lowess") {
    fit <- lowess(df$qx, df$qy, f = span)
    qx <- fit$x
    qy <- fit$y
  } else if (smooth == "ns") {
    kn <- quantile(qx, knots)
    fit <- lm(qy ~ ns(qx, knots = kn), df)
    qx <- seq(0, qx[length(qx)] *2, length.out = 1000)
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
#' 
#' @param x A 'qqmap' class object created by [quantile_map()].
#' @param points Logical whether to show quantile points.
#' @param ... Optional plotting parameters passed to [plot()].
#' @returns No return value. Produces a QQ plot using base graphics with a red
#'   line showing the conversion function.
#' @importFrom graphics lines
#' @export
plot.qqmap <- function(x, points = TRUE, ...) {
  new.args <- list(...)
  if (points) {
    args <- list(x = x$quantiles$qx, y = x$quantiles$qy, cex = 0.5, las = 1,
                 xlab = x$xlab, ylab = x$ylab)
    if (length(new.args)) args[names(new.args)] <- new.args
    do.call("plot", args)
    
    xr <- par("usr")[1:2]
    xr[1] <- max(c(0, xr[1]))
    px <- seq(xr[1], xr[2], length.out = 1000)
    lines(px, x$map(px), col = "red", lwd = 1.5)
  } else {
    xr <- c(x$quantiles$qx[1], x$quantiles$qx[nrow(x$quantiles)])
    px <- seq(xr[1], xr[2], length.out = 1000)
    args = list(x = px, y = x$map(px), col = "red", type = "l", las = 1,
                lwd = 1.5, xlab = x$xlab, ylab = x$ylab)
    if (length(new.args)) args[names(new.args)] <- new.args
    do.call("plot", args)
  }
}


linear_qq <- function(x, y) {
  qm <- quantile_map(x, y)
  # xlim <- max(qm$quantiles$qx)
  xlim <- qm$quantiles$qx[nrow(qm$quantiles) -1]
  s <- seq(xlim / 2, xlim, length.out = 10)
  ratios <- qm$map(s) / s
  mean(ratios)
}
