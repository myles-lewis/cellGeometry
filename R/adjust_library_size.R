
#' Adjust count matrix by library size
#' 
#' Simple tool for adjusting raw count matrix by total library size. Library
#' size is calculated as column sums and columns are scaled to the median total
#' library size.
#' 
#' @param x Read count matrix with genes in rows and samples in columns.
#' @returns Matrix of adjusted read counts
#' @importFrom stats median
#' @export

adjust_library_size <- function(x) {
  cs <- colSums(x)
  a <- cs / median(cs)
  t(t(x) / a)
}
