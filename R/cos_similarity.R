
#' Gene signature cosine similarity matrix
#'
#' Computes the cosine similarity matrix from the gene signature matrix of a
#' `cellMarkers` object or any matrix. Note that this function computes cosine
#' similarity between matrix columns, unlike [dist()] which computes the
#' distance metric between matrix rows.
#' 
#' @param x Either a matrix or a 'cellMarkers' class or 'deconv' class object.
#' @param use_filter Logical whether to use filtered gene signature.
#' @returns A symmetric similarity matrix.
#' @export
#' 
cos_similarity <- function(x, use_filter = NULL) {
  if (inherits(x, "deconv")) {
    x <- x$mk
    if (is.null(use_filter)) use_filter <- x$call$use_filter
  }
  if (inherits(x, "cellMarkers")) {
    if (is.null(use_filter)) use_filter <- TRUE
    gs <- x$geneset
    x <- if (use_filter) x$genemeans_filtered[gs, ] else x$genemeans[gs, ]
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
