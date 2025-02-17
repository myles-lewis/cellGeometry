
#' Single-cell apply a function to a matrix split by a factor
#' 
#' Workhorse function designed to handle large scRNA-Seq gene expression
#' matrices such as embedded Seurat matrices, and apply a function to columns of
#' the matrix split as a ragged array by an index factor, similar to [tapply()],
#' [by()] or [aggregate()]. Note that here the index is applied to columns as
#' these represent cells in the single-cell format, rather than rows as in
#' [aggregate()]. Very large matrices are handled by slicing rows into blocks to
#' avoid excess memory requirements.
#' 
#' @param x matrix, sparse matrix or DelayedMatrix of raw counts with genes in
#'   rows and cells in columns.
#' @param INDEX a factor whose length matches the number of columns in `x`. It
#'   is coerced to a factor. `NA` are tolerated and the matching columns in `x`
#'   are skipped.
#' @param FUN Function to be applied to each subblock of the matrix.
#' @param combine A function or a name of a function to apply to the list
#'   output to bind the final results together, e.g. 'cbind' or 'rbind' to
#'   return a matrix, or 'unlist' to return a vector.
#' @param combine2 A function or a name of a function to combine results after
#'   slicing. As the function is usually applied to blocks of 30000 genes or so,
#'   the result is usually a vector with an element per gene. Hence 'c' is the
#'   default function for combining vectors into a single longer vector. However
#'   if each gene returns a number of results (e.g. a vector or dataframe), then
#'   `combine2` could be set to 'rbind'.
#' @param progress Logical, whether to show progress.
#' @param sliceMem Max amount of memory in GB to allow for each subsetted count
#'   matrix object. When `x` is subsetted by each cell subclass, if the amount
#'   of memory would be above `sliceMem` then slicing is activated and the
#'   subsetted count matrix is divided into chunks and processed separately.
#'   The limit is just under 17.2 GB (2^34 / 1e9). At this level the subsetted
#'   matrix breaches the long vector limit (>2^31 elements).
#' @param cores Integer, number of cores to use for parallelisation using 
#'   `mclapply()`. Parallelisation is not available on windows. Warning:
#'   parallelisation increases the memory requirement by multiples of
#'   `sliceMem`.
#' @param ... Optional arguments passed to `FUN`.
#' @details
#' The limit on `sliceMem` is that the number of elements manipulated in each
#' block must be
#' kept below the long vector limit of 2^31 (around 2e9). Increasing `cores`
#' requires substantial amounts of spare RAM. `combine` works
#' in a similar way to `.combine` in `foreach()`; it works across the levels in
#' `INDEX`. `combine2` is nested and works across slices of genes (an inner
#' loop), so it is only invoked if slicing occurs which is when a matrix has a
#' larger memory footprint than `sliceMem`.
#' 
#' @returns By default returns a list, unless `combine` is invoked in which case
#'   the returned data type will depend on the functions specified by `FUN` and
#'   `combine`.
#' @seealso [scmean()] which applies a fixed function `logmean()` in a similar
#'   manner, and [slapply()] which applies a function to a big matrix with
#'   slicing but without splitting by an index factor.
#' @author Myles Lewis
#' @examples
#' # equivalent
#' m <- matrix(sample(0:100, 1000, replace = TRUE), nrow = 10)
#' cell_index <- sample(letters[1:5], 100, replace = TRUE)
#' o <- scmean(m, cell_index)
#' o2 <- scapply(m, cell_index, function(x) rowMeans(log2(x +1)),
#'               combine = "cbind")
#' identical(o, o2)
#' 
#' @importFrom mcprogress pmclapply
#' @export

scapply <- function(x, INDEX, FUN, combine = NULL, combine2 = "c",
                    progress = TRUE,
                    sliceMem = 16, cores = 1L, ...) {
  if (!is.factor(INDEX)) INDEX <- factor(INDEX)
  tab <- table(INDEX)
  ok <- !is.na(INDEX)
  dimx <- as.numeric(dim(x))
  if (dimx[2] != length(INDEX)) stop("Incompatible dimensions")
  if (sliceMem > 2^34 / 1e9) message("`sliceMem` is above the long vector limit")
  lev <- levels(INDEX)
  
  out <- pmclapply(seq_along(lev), function(i) {
    ind <- lev[i]
    c_index <- which(INDEX == ind & ok)
    n <- length(c_index) * dimx[1]
    bloc <- ceiling(n *8 / (sliceMem * 1e9))
    if (bloc == 1) {
      # unsliced
      xsub <- as.matrix(x[, c_index]) |> suppressWarnings()
      ret <- FUN(xsub, ...)
      xsub <- NULL
      return(ret)
    }
    # slice
    sliceSize <- ceiling(dimx[1] / bloc)
    s <- sliceIndex(dimx[1], sliceSize)
    out2 <- lapply(s, function(j) {
      xsub <- as.matrix(x[j, c_index]) |> suppressWarnings()
      ret <- FUN(xsub, ...)
      xsub <- NULL
      ret
    })
    if (!is.null(combine2)) out2 <- do.call(combine2, out2)
    out2
  }, mc.cores = cores, progress = progress)
  
  names(out) <- levels(INDEX)
  if (!is.null(combine)) out <- do.call(combine, out)
  out
}


#' Apply a function to a big matrix by slicing
#' 
#' Workhorse function ('slice apply') designed to handle large scRNA-Seq gene
#' expression matrices such as embedded Seurat matrices, and apply a function to
#' the whole matrix. Very large matrices are handled by slicing rows into blocks
#' to avoid excess memory requirements.
#' 
#' @param x matrix, sparse matrix or DelayedMatrix of raw counts with genes in
#'   rows and cells in columns.
#' @param FUN Function to be applied to each subblock of the matrix.
#' @param combine A function or a name of a function to combine results after
#'   slicing. As the function is usually applied to blocks of 30000 genes or so,
#'   the result is usually a vector with an element per gene. Hence 'c' is the
#'   default function for combining vectors into a single longer vector. However
#'   if each gene row returns a number of results (e.g. a vector or dataframe),
#'   then `combine` could be set to 'rbind'.
#' @param progress Logical, whether to show progress.
#' @param sliceMem Max amount of memory in GB to allow for each subsetted count
#'   matrix object. When `x` is subsetted by each cell subclass, if the amount
#'   of memory would be above `sliceMem` then slicing is activated and the
#'   subsetted count matrix is divided into chunks and processed separately.
#'   The limit is just under 17.2 GB (2^34 / 1e9). At this level the subsetted
#'   matrix breaches the long vector limit (>2^31 elements).
#' @param cores Integer, number of cores to use for parallelisation using 
#'   `mclapply()`. Parallelisation is not available on windows. Warning:
#'   parallelisation has increased memory requirements.
#' @param ... Optional arguments passed to `FUN`.
#' @details
#' The limit on `sliceMem` is that the number of elements manipulated in each
#' block must be kept below the long vector limit of 2^31 (around 2e9).
#' Increasing `cores` requires substantial amounts of spare RAM. `combine` works
#' in a similar way to `.combine` in `foreach()` across slices of genes; it is
#' only invoked if slicing occurs.
#' 
#' @returns The returned data type will depend on the functions specified by
#'   `FUN` and `combine`.
#' @seealso [scapply()]
#' @author Myles Lewis
#' @export

slapply <- function(x, FUN, combine = "c",
                    progress = TRUE,
                    sliceMem = 16, cores = 1L, ...) {
  dimx <- as.numeric(dim(x))
  if (sliceMem > 2^34 / 1e9) message("`sliceMem` is above the long vector limit")
  
  n <- dimx[1] * dimx[2]
  bloc <- ceiling(n *8 / (sliceMem * 1e9))
  if (bloc == 1) {
    # unsliced
    ret <- FUN(as.matrix(x), ...) |> suppressWarnings()
    return(ret)
  }
  # slice
  sliceSize <- ceiling(dimx[1] / bloc)
  s <- sliceIndex(dimx[1], sliceSize)
  out <- pmclapply(s, function(j) {
    xsub <- as.matrix(x[j, ]) |> suppressWarnings()
    ret <- FUN(xsub, ...)
    xsub <- NULL
    ret
  }, mc.cores = cores, progress = progress)
  if (!is.null(combine)) out <- do.call(combine, out)
  out
}
