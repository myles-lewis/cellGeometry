
#' Gene signature cosine similarity matrix
#'
#' Computes the cosine similarity matrix from the gene signature matrix of a
#' `cellMarkers` object or any matrix. Note that this function computes cosine
#' similarity between matrix columns, unlike [dist()] which computes the
#' distance metric between matrix rows.
#' 
#' @param x Either a `cellMarkers` class object or a matrix.
#' @param use_filter Logical whether to use filtered gene signature.
#' @returns A symmetric similarity matrix.
#' @export
#' 
cos_similarity <- function(x, use_filter = TRUE) {
  if (inherits(x, "cellMarkers")) {
    gs <- mk$geneset
    x <- if (use_filter) mk$genemeans_filtered[gs, ] else mk$genemeans[gs, ]
  }
  cos_sim(x)
}


# cosine similarity of columns of a matrix
cos_sim <- function(x) {
  sx <- scaleSphere(t(x))
  out <- tcrossprod(sx)
  diag(out) <- 1
  out
}
