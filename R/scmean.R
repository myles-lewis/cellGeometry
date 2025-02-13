
#' Single-cell mean log gene expression across cell types
#' 
#' Workhorse function which takes as input a scRNA-Seq gene expression matrix
#' such as embedded in a Seurat object, calculates log2(counts +1) and averages
#' gene expression over a vector specifying cell subclasses or cell types. Very
#' large matrices are handled by slicing rows into blocks to avoid excess memory
#' requirements.
#' 
#' @param x matrix, sparse matrix or DelayedMatrix of raw counts with genes in
#'   rows and cells in columns.
#' @param celltype a vector of cell subclasses or types whose length matches the
#'   number of columns in `x`. It is coerced to a factor. `NA` are tolerated and
#'   the matching columns in `x` are skipped.
#' @param FUN Function for applying mean. When applied to a matrix of count
#'   values, this must return a vector. Recommended options are `logmean` (the
#'   default) or `trimmean`.
#' @param postFUN Optional function to be applied to whole matrix after mean has
#'   been calculated, e.g. `log2s`.
#' @param verbose Logical, whether to print messages.
#' @param sliceMem Max amount of memory in GB to allow for each subsetted count
#'   matrix object. When `x` is subsetted by each cell subclass, if the amount
#'   of memory would be above `sliceMem` then slicing is activated and the
#'   subsetted count matrix is divided into chunks and processed separately.
#'   This is indicated by addition of '...' in the timings. The limit is just
#'   under 17.2 GB (2^34 / 1e9). At this level the subsetted matrix breaches the
#'   long vector limit (>2^31 elements).
#' @param cores Integer, number of cores to use for parallelisation using 
#'   `mclapply()`. Parallelisation is not available on windows. Warning:
#'   parallelisation increases the memory requirement by multiples of
#'   `sliceMem`.
#' @details 
#' Mean functions which can be applied by setting `FUN` include `logmean` (the
#' default) which applies row means to log2(counts+1), or `trimmean` which
#' calculates the trimmed mean of the counts after top/bottom 5% of values have
#' been excluded. Alternatively `FUN = rowMeans` calculates the arithmetic mean
#' of counts.
#' 
#' If `FUN = trimmean` or `rowMeans`, `postFUN` needs to be set to `log2s` which
#' is a simple function which applies log2(x+1).
#' 
#' `sliceMem` can be set lower on machines with less RAM, but this will slow the
#' analysis down. `cores` increases the theoretical amount of memory required to
#' around `cores * sliceMem` in GB. For example on a 64 GB machine, we find a
#' significant speed increase with `cores = 3L`. Above this level, there is a
#' risk that memory swap will slow down processing.
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
                   verbose = TRUE,
                   sliceMem = 16, cores = 1L) {
  start0 <- Sys.time()
  if (!is.factor(celltype)) celltype <- factor(celltype)
  ok <- !is.na(celltype)
  dimx <- as.numeric(dim(x))
  if (dimx[2] != length(celltype)) stop("Incompatible dimensions")
  if (sliceMem > 2^34 / 1e9) message("`sliceMem` is above the long vector limit")
  
  # dynamic slicing
  genemeans <- mclapply(levels(celltype), function(i) {
    start <- Sys.time()
    c_index <- which(celltype == i & ok)
    n <- length(c_index) * dimx[1]
    bloc <- ceiling(n *8 / (sliceMem * 1e9))
    
    if (bloc == 1) {
      # unsliced
      xsub <- as.matrix(x[, c_index]) |> suppressWarnings()
      ret <- FUN(xsub)
      xsub <- NULL
      if (verbose) timer(start, paste0(length(c_index), " ", i, "  ("))
      return(ret)
    }
    
    # slice
    sliceSize <- ceiling(dimx[1] / bloc)
    s <- sliceIndex(dimx[1], sliceSize)
    out <- lapply(s, function(j) {
      xsub <- as.matrix(x[j, c_index]) |> suppressWarnings()
      ret <- FUN(xsub)
      xsub <- NULL
      ret
    })
    if (verbose) timer(start, paste0(length(c_index), " ", i, " ... ("))
    unlist(out)
  }, mc.cores = cores, mc.preschedule = FALSE)
  genemeans <- do.call("cbind", genemeans)
  colnames(genemeans) <- levels(celltype)
  
  if (!is.null(postFUN)) genemeans <- postFUN(genemeans)
  if (verbose) timer(start0, "Duration")
  genemeans
}

# xsub <- NULL is faster than rm(list="xsub")
# ought to reduce memory usage by mclapply


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

sliceIndex <- function(nx, sliceSize) {
  if (is.null(sliceSize)) sliceSize <- nx
  sliceSize <- as.integer(sliceSize)
  s <- ceiling(nx / sliceSize)
  excess <- nx %% sliceSize
  lapply(seq_len(s), function(i) {
    if (i==s && excess != 0) return(seq_len(excess) + sliceSize * (i-1L))
    seq_len(sliceSize) + sliceSize * (i-1L)
  })
}

#' @importFrom mcprogress cat_parallel
timer <- function(start, msg = NULL) {
  end <- Sys.time()
  tim <- format(end - start, digits = 3)
  if (is.null(msg)) {
    cat_parallel("(", tim, ")\n")
  } else {
    if (substr(msg, nchar(msg), nchar(msg)) == "(") {
      cat_parallel(msg, tim, ")\n")
    } else {
      cat_parallel(msg, " ", tim, "\n")
    }
  }
}
