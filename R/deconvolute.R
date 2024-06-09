
#' Deconvolute bulk RNA-Seq using single-cell RNA-Seq signature
#'
#' Deconvolution of bulk RNA-Seq using vector projection method with optional
#' compensation for spillover.
#'
#' @param mk object of class 'cellMarkers'
#' @param test matrix of bulk RNA-Seq counts to be deconvoluted
#' @param comp_amount either a single value from 0-1 for the amount of
#'   compensation or a numeric vector with the same length as the number of cell
#'   subclasses to deconvolute.
#' @param group_comp_amount either a single value from 0-1 for the amount of
#'   compensation for cell group analysis or a numeric vector with the same
#'   length as the number of cell groups to deconvolute.
#' @param adjust_comp logical, whether to optimise `comp_amount` to prevent
#'   negative cell proportion projections.
#' @param use_filter logical, whether to use denoised signature matrix
#' @param convert_bulk logical, whether to convert bulk RNA-Seq to scRNA-Seq
#'   scaling.
#' @returns a list object containing the cell proportions of each cell subclass
#'   in each sample in element `subclass`, and the cell proportions of each cell
#'   group in element `group`.
#' @importFrom matrixStats colMins
#' @importFrom stats optimise
#' @export
#'
deconvolute <- function(mk, test, comp_amount = 1,
                        group_comp_amount = 0,
                        adjust_comp = TRUE,
                        use_filter = TRUE,
                        convert_bulk = TRUE) {
  if (!inherits(mk, "cellMarkers")) stop ("Not a 'cellMarkers' class object")
  .call <- match.call()
  
  test <- as.matrix(test)
  
  # group first
  if (!is.null(mk$group_geneset)) {
    cellmat <- if (use_filter) {mk$groupmeans_filtered[mk$group_geneset, ]
    } else mk$groupmeans[mk$group_geneset, ]
    # cellmat <- sc2bulk(cellmat)
    logtest <- log2(test[mk$group_geneset, , drop = FALSE] +1)
    if (convert_bulk) logtest <- bulk2sc(logtest)
    gtest <- deconv_adjust(logtest, cellmat, group_comp_amount, adjust_comp)
  } else {
    gtest <- NULL
  }
  
  # cell subclasses
  cellmat <- if (use_filter) {mk$genemeans_filtered[mk$geneset, ]
  } else mk$genemeans[mk$geneset, ]
  # cellmat <- sc2bulk(cellmat)
  logtest2 <- log2(test[mk$geneset, , drop = FALSE] +1)
  if (convert_bulk) logtest2 <- bulk2sc(logtest2)
  atest <- deconv_adjust(logtest2, cellmat, comp_amount, adjust_comp)
  
  # subclass percent as nested percent of groups
  if (!is.null(gtest)) {
    pc2 <- lapply(levels(mk$cell_table), function(i) {
      pc1 <- gtest$percent[, i]
      ind <- mk$cell_table == i
      subclass_i <- atest$output[, ind, drop = FALSE]
      rs <- rowSums(subclass_i)
      subclass_i / rs * pc1
    })
    nest_percent <- do.call(cbind, pc2)
  } else nest_percent <- NULL
  
  out <- list(call = .call, mk = mk, subclass = atest, group = gtest,
              nest_percent = nest_percent, comp_amount = comp_amount)
  class(out) <- "deconv"
  out
}

deconv_adjust <- function(test, cellmat, comp_amount = 0,
                          adjust_comp = TRUE) {
  comp_amount <- rep_len(comp_amount, ncol(cellmat))
  names(comp_amount) <- colnames(cellmat)
  
  atest <- deconv(test, cellmat, comp_amount)
  if (any(atest$output < 0)) {
    if (adjust_comp) {
      message("optimising compensation")
      minout <- colMins(atest$output)
      w <- which(minout < 0)
      newcomps <- vapply(w, function(i) {
        f <- function(x) {
          newcomp <- comp_amount
          newcomp[i] <- x
          ntest <- deconv(test, cellmat, comp_amount = newcomp)
          min(ntest$output[, i], na.rm = TRUE)^2
        }
        xmin <- optimise(f, c(0, comp_amount[i]))
        xmin$minimum
      }, numeric(1))
      comp_amount[w] <- newcomps
      atest <- deconv(test, cellmat, comp_amount = comp_amount)
      atest$output[atest$output < 0] <- 0
    } else message("negative cell proportion projection detected")
  }
  atest
}

deconv <- function(test, cellmat, comp_amount = 0, equalWeight = FALSE) {
  if (!(length(comp_amount) %in% c(1, ncol(cellmat))))
    stop('comp_amount must be either single number or vector of length matching cellmat cols')
  if (!identical(rownames(test), rownames(cellmat)))
    stop('test and cell matrices must have same genes (rownames)')
  m_itself <- dotprod(cellmat, cellmat, equalWeight)
  rawcomp <- solve(m_itself)
  mixcomp <- solve(m_itself, t(comp_amount * diag(nrow(m_itself)) + (1-comp_amount) * t(m_itself)))
  output <- dotprod(test, cellmat, equalWeight) %*% mixcomp
  percent <- output / rowSums(output) * 100
  list(output = output, percent = percent, spillover = m_itself,
       compensation = mixcomp, rawcomp = rawcomp, comp_amount = comp_amount)
}


approxfun.matrix <- function(x, FUN) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    out <- vapply(seq_len(ncol(x)), function(i) {
      FUN(x[, i])
    }, numeric(nrow(x)))
    dimnames(out) <- dimnames(x)
    return(out)
  }
  FUN(x)
}

#' @importFrom stats approxfun
sc2bulk <- function(x) {
  sc2bulkfun <- approxfun(x = celseqfit$celseq, y = celseqfit$pred.bulk,
                          yleft = 0, rule = 2)
  approxfun.matrix(x, sc2bulkfun)
}


bulk2sc <- function(x) {
  bulk2scfun <- approxfun(x = celseqfit$pred.bulk, y = celseqfit$celseq,
                          yleft = 0, rule = 2)
  approxfun.matrix(x, bulk2scfun)
}
