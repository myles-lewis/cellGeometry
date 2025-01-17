
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
#' @param FUN Function for applying mean. When applied to a matrix of count
#'   values, this must return a vector. Recommended options are `logmean` (the
#'   default) or `trimmean`.
#' @param postFUN Optional function to be applied to whole matrix after mean has
#'   been calculated, e.g. `log2s`.
#' @param big Logical, whether to invoke slicing of `x` into rows. This is
#'   invoked automatically if `x` is a large matrix with >2^31 elements.
#' @param verbose Logical, whether to print messages.
#' @param sliceSize Integer, number of rows of `x` to use in each slice if 
#'   `big = TRUE`.
#' @param cores Integer, number of cores to use for parallelisation using 
#'   `mclapply()`. Parallelisation is not available on windows. Warning:
#'   parallelisation has increased memory requirements.
#' @details
#' We find a significant speed up with `cores = 2`, which is almost twice as
#' fast as single core, but not much to be gained beyond this possibly due to
#' limits on memory traffic. The main speed up is in assigning the decompression
#' of a block from the sparse matrix to more than 1 core. Increasing `sliceSize`
#' also gives a speed up, but the limit on `sliceSize` is that the number of
#' elements manipulated in each block (i.e. `sliceSize` x number of cells in a
#' given subclass/group) must be kept below the long vector limit of 2^31
#' (around 2e9). Increasing `cores` and/or `sliceSize` requires substantial
#' amounts of spare RAM.
#' 
#' Mean functions which can be applied by setting `FUN` include `logmean` (the
#' default) which applies row means to log2(counts+1), or `trimmean` which
#' calculates the trimmed mean of the counts after top/bottom 5% of values have
#' been excluded. Alternatively `FUN = rowMeans` calculates the arithmetic mean
#' of counts.
#' 
#' If `FUN = trimmean` or `rowMeans`, `postFUN` needs to be set to `log2s` which
#' is a simple function which applies log2(x+1).
#' 
#' @returns a matrix of mean log2 gene expression across cell types with genes
#'   in rows and cell types in columns.
#' @seealso [scapply()] which is a more general version which can apply any
#'   function to the matrix. \code{\link[=logmean]{logmean}},  
#'   \code{\link[=trimmean]{trimmean}} are options for controlling the type of
#'   mean applied.
#' @author Myles Lewis
#' @importFrom parallel mclapply
#' @export

scmean <- function(x, celltype,
                   FUN = logmean, postFUN = NULL,
                   # FUN = trimmean,
                   # postFUN = function(x) log2(x +1),
                   big = NULL, verbose = TRUE,
                   sliceSize = 5000L, cores = 1L) {
  start0 <- Sys.time()
  if (!is.factor(celltype)) celltype <- factor(celltype)
  if (any(table(celltype) * as.numeric(sliceSize) > 2^31))
    message("Warning: >2^31 matrix elements anticipated. `sliceSize` is too large")
  ok <- !is.na(celltype)
  dimx <- dim(x)
  if (dimx[2] != length(celltype)) stop("Incompatible dimensions")
  if (as.numeric(dimx[1]) * as.numeric(dimx[2]) > 2^31) big <- TRUE
  
  if (is.null(big) || !big) {
    # small matrix
    genemeans <- vapply(levels(celltype), function(i) {
      FUN(as.matrix(x[, which(celltype==i & ok)])) |> suppressWarnings()
    }, numeric(dimx[1]))
    if (!is.null(postFUN)) genemeans <- postFUN(genemeans)
    return(genemeans)
  }
  
  # large matrix
  s <- sliceIndex(dimx[1], sliceSize)
  genemeans <- vapply(levels(celltype), function(i) {
    start <- Sys.time()
    c_index <- which(celltype == i & ok)
    if (verbose) cat(length(c_index), i, " ")
    out <- parallel::mclapply(s, function(j) {
      FUN(as.matrix(x[j, c_index])) |> suppressWarnings()
    }, mc.cores = cores)
    if (verbose) timer(start)
    unlist(out)
  }, numeric(dimx[1]))
  
  if (!is.null(postFUN)) genemeans <- postFUN(genemeans)
  if (verbose) timer(start0, "Duration")
  genemeans
}

#' Mean Objects
#'
#' Functions designed for use with [scmean()] to calculate mean gene expression
#' in each cell cluster.
#'
#' @param x A count matrix
#' @returns Numeric vector of mean values.
#'
#'   `logmean` applies `log2(x+1)` then calculates `rowMeans`.
#'
#'   `trimmean` applies a trimmed mean to each row of gene counts, excluding the
#'   top and bottom 5% of values which helps to exclude outliers. When `trimmean` is used with
#'   [scmean()], it is important to set `postFUN = log2s`. This simply applies
#'   log2(x+1) after the trimmed mean of counts has been calculated.
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

sliceIndex <- function(nx, sliceSize = 2000) {
  if (is.null(sliceSize)) sliceSize <- nx
  sliceSize <- as.integer(sliceSize)
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
