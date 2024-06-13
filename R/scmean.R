
#' Single-cell mean log gene expression across cell types
#' 
#' Workhorse function which takes as input a scRNA-Seq gene expression matrix
#' such as embedded in a Seurat object, calculates log2(counts +1) and averages
#' gene expression over a vector specifying cell subclasses or cell types. Very
#' large matrices are handled by slicing rows into blocks to avoid excess memory
#' requirements.
#' 
#' @param x matrix or sparse matrix of raw counts with genes in rows and cells
#'   in columns.
#' @param celltype a vector of cell subclasses or types whose length matches the
#'   number of columns in `x`. It is coerced to a factor. `NA` are tolerated and
#'   the matching columns in `x` are skipped.
#' @param big Logical, whether to invoke slicing of `x` into rows. This is
#'   invoked automatically if `x` is a large matrix with >2^31 elements.
#' @param verbose Logical, whether to print messages.
#' @returns a matrix of mean log2 gene expression across cell types with genes
#'   in rows and cell types in columns.
#' @export

scmean <- function(x, celltype, big = NULL, verbose = TRUE) {
  start0 <- Sys.time()
  if (!is.factor(celltype)) celltype <- factor(celltype)
  ok <- !is.na(celltype)
  dimx <- dim(x)
  if (dimx[2] != length(celltype)) stop("Incompatible dimensions")
  if (as.numeric(dimx[1]) * as.numeric(dimx[2]) > 2^31) big <- TRUE
  if (is.null(big) || !big) {
    # small matrix
    genemeans <- vapply(levels(celltype), function(i) {
      logmean(x[, which(celltype==i & ok)])
    }, numeric(dimx[1]))
    return(genemeans)
  }
  # large matrix
  s <- sliceIndex(dimx[1])
  genemeans <- vapply(levels(celltype), function(i) {
    if (verbose) cat(i, " ")
    start <- Sys.time()
    c_index <- which(celltype == i & ok)
    out <- lapply(s, function(j) {
      logmean(as.matrix(x[j, c_index])) |> suppressWarnings()
    })
    if (verbose) timer(start)
    unlist(out)
  }, numeric(dimx[1]))
  
  if (verbose) timer(start0, "Duration")
  genemeans
}

logmean <- function(x) rowMeans(log2(x +1))

sliceIndex <- function(nx, sliceSize = 2000) {
  if (is.null(sliceSize)) sliceSize <- nx
  s <- ceiling(nx / sliceSize)
  excess <- nx %% sliceSize
  lapply(seq_len(s), function(i) {
    if (i==s && excess != 0) return(seq_len(excess) + sliceSize * (i-1L))
    seq_len(sliceSize) + sliceSize * (i-1L)
  })
}

timer <- function(start, msg = NULL) {
  end <- Sys.time()
  if (is.null(msg)) {
    cat(paste0("(", format(end - start, digits = 3), ")\n"))
  } else {
    cat(msg, format(end - start, digits = 3), "\n")
  }
}
