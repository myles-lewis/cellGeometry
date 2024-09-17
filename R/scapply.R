
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
#' @param x matrix or sparse matrix of raw counts with genes in rows and cells
#'   in columns.
#' @param INDEX a factor whose length matches the number of columns in `x`. It
#'   is coerced to a factor. `NA` are tolerated and the matching columns in `x`
#'   are skipped.
#' @param FUN Function to be applied to each subblock of the matrix.
#' @param combine A function or a name of a function to apply to the list
#'   output to bind the final results together, e.g. 'cbind' or 'rbind' to
#'   return a matrix, or 'unlist' to return a vector.
#' @param combine2 A function or a name of a function to combine results after
#'   slicing, i.e. only invoked if `big` is `TRUE`. As the function is usually
#'   applied to blocks of 5000 genes or so, the result is usually a vector wih
#'   an element per gene. Hence 'c' is the default function for combining
#'   vectors into a single longer vector. However if each gene returns a number
#'   of results (e.g. a vector or dataframe), then `combine2` could be set to
#'   'rbind'.
#' @param big Logical, whether to invoke slicing of `x` into rows. This is
#'   invoked automatically if `x` is a large matrix with >2^31 elements.
#' @param verbose Logical, whether to show progress.
#' @param sliceSize Integer, number of rows of `x` to use in each slice if 
#'   `big = TRUE`.
#' @param cores Integer, number of cores to use for parallelisation using 
#'   `mclapply()`. Parallelisation is not available on windows. Warning:
#'   parallelisation has increased memory requirements.
#' @param ... Optional arguments passed to `FUN`.
#' @details
#' The limit on `sliceSize` is that the number of elements manipulated in each
#' block (i.e. `sliceSize` x number of cells in a given subclass/group) must be
#' kept below the long vector limit of 2^31 (around 2e9). Increasing `cores`
#' and/or `sliceSize` requires substantial amounts of spare RAM. `combine` works
#' in a similar way to `.combine` in `foreach()`; it works across the levels in
#' `INDEX`. `combine2` is nested and works across slices of genes (an inner
#' loop), so it is only invoked if `big` is `TRUE`.
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
#' @export

scapply <- function(x, INDEX, FUN, combine = NULL, combine2 = "c",
                    big = NULL, verbose = TRUE,
                    sliceSize = 5000L, cores = 1L, ...) {
  if (!is.factor(INDEX)) INDEX <- factor(INDEX)
  tab <- table(INDEX)
  if (any(tab * as.numeric(sliceSize) > 2^31))
    message("Warning: >2^31 matrix elements anticipated. `sliceSize` is too large")
  ok <- !is.na(INDEX)
  dimx <- dim(x)
  if (dimx[2] != length(INDEX)) stop("Incompatible dimensions")
  if (as.numeric(dimx[1]) * as.numeric(dimx[2]) > 2^31) big <- TRUE
  dimmax <- as.numeric(max(tab)) * sliceSize
  mem <- structure(dimmax * 8 * cores, class = "object_size")
  if (verbose & mem > 1e9)
    message("Absolute minimum memory headroom ", format(mem, units = "GB"))
  
  if (verbose) pb <- txtProgressBar2()
  lev <- levels(INDEX)
  if (is.null(big) || !big) {
    # small matrix
    out <- lapply(seq_along(lev), function(i) {
      ind <- lev[i]
      if (verbose) setTxtProgressBar(pb, i / length(lev))
      FUN(as.matrix(x[, which(INDEX==ind & ok)]), ...)
    })
  } else {
    # large matrix
    s <- sliceIndex(dimx[1], sliceSize)
    out <- lapply(seq_along(lev), function(i) {
      if (verbose) setTxtProgressBar(pb, i / length(lev))
      ind <- lev[i]
      c_index <- which(INDEX == ind & ok)
      # slicing inner loop
      out2 <- parallel::mclapply(s, function(j) {
        FUN(as.matrix(x[j, c_index]), ...) |> suppressWarnings()
      }, mc.cores = cores)
      if (!is.null(combine2)) out2 <- do.call(combine2, out2)
      out2
    })
  }
  names(out) <- levels(INDEX)
  if (!is.null(combine)) out <- do.call(combine, out)
  if (verbose) close(pb)
  
  out
}


#' Apply a function to a big matrix by slicing
#' 
#' Workhorse function ('slice apply') designed to handle large scRNA-Seq gene
#' expression matrices such as embedded Seurat matrices, and apply a function to
#' the whole matrix. Very large matrices are handled by slicing rows into blocks
#' to avoid excess memory requirements.
#' 
#' @param x matrix or sparse matrix of raw counts with genes in rows and cells
#'   in columns.
#' @param FUN Function to be applied to each subblock of the matrix.
#' @param combine A function or a name of a function to combine results after
#'   slicing, i.e. only invoked if `big` is `TRUE`. As the function is usually
#'   applied to blocks of 1000 genes or so, the result is usually a vector wih
#'   an element per gene. Hence 'c' is the default function for combining
#'   vectors into a single longer vector. However if each gene row returns a
#'   number of results (e.g. a vector or dataframe), then `combine` could be set
#'   to 'rbind'.
#' @param big Logical, whether to invoke slicing of `x` into rows. This is
#'   invoked automatically if `x` is a large matrix with >2^31 elements.
#' @param verbose Logical, whether to show progress.
#' @param sliceSize Integer, number of rows of `x` to use in each slice if 
#'   `big = TRUE`.
#' @param cores Integer, number of cores to use for parallelisation using 
#'   `mclapply()`. Parallelisation is not available on windows. Warning:
#'   parallelisation has increased memory requirements.
#' @param ... Optional arguments passed to `FUN`.
#' @details
#' The limit on `sliceSize` is that the number of elements manipulated in each
#' block (i.e. `sliceSize` x number of cells in a given subclass/group) must be
#' kept below the long vector limit of 2^31 (around 2e9). Increasing `cores`
#' and/or `sliceSize` requires substantial amounts of spare RAM. `combine` works
#' in a similar way to `.combine` in `foreach()` across slices of genes; it
#' is only invoked if `big` is `TRUE`.
#' 
#' @returns The returned data type will depend on the functions specified by
#'   `FUN` and `combine`.
#' @seealso [scapply()]
#' @author Myles Lewis
#' @export

slapply <- function(x, FUN, combine = "c",
                    big = NULL, verbose = TRUE,
                    sliceSize = 1000L, cores = 1L, ...) {
  start <- Sys.time()
  dimx <- dim(x)
  if (as.numeric(dimx[1]) * as.numeric(dimx[2]) > 2^31) big <- TRUE
  if (as.numeric(dimx[2]) * as.numeric(sliceSize) > 2^31)
    message("Warning: >2^31 matrix elements anticipated. `sliceSize` is too large")
  dimmax <- dimx[2] * sliceSize
  mem <- structure(dimmax * 8 * cores, class = "object_size")
  if (verbose & mem > 1e9)
    message("Absolute minimum memory headroom ", format(mem, units = "GB"))
  
  if (is.null(big) || !big) {
    # small matrix
    out <- FUN(as.matrix(x), ...)
  } else {
    # slice large matrix
    if (verbose) {
      if (cores == 1) pb <- txtProgressBar2() else mcProgressBar(0)
    }
    s <- sliceIndex(dimx[1], sliceSize)
    out <- parallel::mclapply(seq_along(s), function(j) {
      if (verbose) {
        if (cores == 1) {
          setTxtProgressBar(pb, j / length(s))
        } else mcProgressBar(j, length(s), cores)
      }
      ind <- s[[j]]
      FUN(as.matrix(x[ind, ]), ...) |> suppressWarnings()
    }, mc.cores = cores)
    if (!is.null(combine)) out <- do.call(combine, out)
    if (verbose) {
      if (cores == 1) close(pb) else closeProgress(start)
    }
  }
  
  out
}
