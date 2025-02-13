# calculate both mean (log2(counts +1))
# and log2(arithmetic mean +1)
# returns list of 2 matrices

scmean2 <- function(x, celltype,
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
  gm <- mclapply(levels(celltype), function(i) {
    start <- Sys.time()
    c_index <- which(celltype == i & ok)
    n <- length(c_index) * dimx[1]
    bloc <- ceiling(n *8 / (sliceMem * 1e9))
    
    if (bloc == 1) {
      # unsliced
      xsub <- as.matrix(x[, c_index]) |> suppressWarnings()
      m1 <- FUN(xsub) |> suppressWarnings()
      m2 <- rowMeans(xsub)
      xsub <- NULL
      if (verbose) timer(start, paste0(length(c_index), " ", i, "  ("))
      return(cbind(m1, m2))
    }
    
    # slice
    sliceSize <- ceiling(dimx[1] / bloc)
    s <- sliceIndex(dimx[1], sliceSize)
    out <- lapply(s, function(j) {
      xsub <- as.matrix(x[j, c_index]) |> suppressWarnings()
      m1 <- FUN(xsub) |> suppressWarnings()
      m2 <- rowMeans(xsub)
      xsub <- NULL
      cbind(m1, m2)
    })
    if (verbose) timer(start, paste0(length(c_index), " ", i, " ... ("))
    do.call(rbind, out)
  }, mc.cores = cores)
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
