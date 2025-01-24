# calculate both mean (log2(counts +1))
# and log2(arithmetic mean +1)
# returns list of 2 matrices

scmean2 <- function(x, celltype,
                    FUN = logmean, postFUN = NULL,
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
    gm <- lapply(levels(celltype), function(i) {
      xsub <- as.matrix(x[, which(celltype==i & ok)])
      m1 <- FUN(xsub) |> suppressWarnings()
      m2 <- rowMeans(xsub) |> suppressWarnings()
      cbind(m1, m2)
    })
    names(gm) <- levels(celltype)
    return(cols2mats(gm, postFUN))
  }
  
  # large matrix
  s <- sliceIndex(dimx[1], sliceSize)
  gm <- lapply(levels(celltype), function(i) {
    start <- Sys.time()
    c_index <- which(celltype == i & ok)
    if (verbose) cat(length(c_index), i, " ")
    out <- parallel::mclapply(s, function(j) {
      xsub <- as.matrix(x[j, c_index])
      m1 <- FUN(xsub) |> suppressWarnings()
      m2 <- rowMeans(xsub) |> suppressWarnings()
      cbind(m1, m2)
    }, mc.cores = cores)
    if (verbose) timer(start)
    do.call(rbind, out)
  })
  names(gm) <- levels(celltype)
  if (verbose) timer(start0, "Duration")
  cols2mats(gm, postFUN)
}

# separate log means and arith means
cols2mats <- function(gm, postFUN) {
  gm1 <- lapply(gm, function(x) x[, 1])
  gm1 <- do.call(cbind, gm1)
  colnames(gm1) <- names(gm)
  if (!is.null(postFUN)) gm1 <- postFUN(gm1)
  
  gm2 <- lapply(gm, function(x) x[, 2])
  gm2 <- do.call(cbind, gm2)
  colnames(gm2) <- names(gm)
  gm2 <- log2(gm2 +1)
  list(gm1, gm2)
}
