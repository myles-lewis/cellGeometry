
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
#' @param logfirst Logical whether log2 +1 is applied to counts first before
#'   mean is applied, or applied after the mean is calculated.
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
#' @returns a matrix of mean log2 gene expression across cell types with genes
#'   in rows and cell types in columns.
#' @seealso [scapply()] which is a more general version which can apply any
#'   function to the matrix.
#' @author Myles Lewis
#' @importFrom parallel mclapply
#' @export

scmean <- function(x, celltype, logfirst = TRUE, big = NULL, verbose = TRUE,
                   sliceSize = 5000L, cores = 1L) {
  start0 <- Sys.time()
  if (!is.factor(celltype)) celltype <- factor(celltype)
  if (any(table(celltype) * as.numeric(sliceSize) > 2^31))
    message("Warning: >2^31 matrix elements anticipated. `sliceSize` is too large")
  ok <- !is.na(celltype)
  dimx <- dim(x)
  if (dimx[2] != length(celltype)) stop("Incompatible dimensions")
  FUN <- if (logfirst) logmean else rowMeans
  if (as.numeric(dimx[1]) * as.numeric(dimx[2]) > 2^31) big <- TRUE
  if (is.null(big) || !big) {
    # small matrix
    genemeans <- vapply(levels(celltype), function(i) {
      FUN(as.matrix(x[, which(celltype==i & ok)])) |> suppressWarnings()
    }, numeric(dimx[1]))
    if (!logfirst) genemeans <- log2(genemeans +1)
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
  if (!logfirst) genemeans <- log2(genemeans +1)
  
  if (verbose) timer(start0, "Duration")
  genemeans
}

logmean <- function(x) rowMeans(log2(x +1))

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
