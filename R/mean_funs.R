
#' Mean Objects
#'
#' Functions designed for use with [scmean()] to calculate mean gene expression
#' in each cell cluster across matrix rows.
#'
#' @param x A count matrix
#' @returns Numeric vector of mean values.
#'
#' `logmean` applies `log2(x+1)` then calculates `rowMeans`.
#'
#' `trimmean` applies a trimmed mean to each row of gene counts, excluding the
#' top and bottom 5% of values which helps to exclude outliers. Note, this needs
#' the `Rfast2` package to be installed. When `trimmean` is used with
#' [scmean()], `postFUN` is typically set to `log2s`. This simply applies
#' log2(x+1) after the trimmed mean of counts has been calculated.
#' @importFrom DelayedArray rowMeans
#' @export

logmean <- function(x) rowMeans(log2(x +1))

#' @rdname logmean
trimmean <- function(x) {
  tm <- Rfast2::rowTrimMean(x)
  names(tm) <- rownames(x)
  tm
}

#' @rdname logmean
log2s <- function(x) log2(x+1)
