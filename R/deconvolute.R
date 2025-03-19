
#' Deconvolute bulk RNA-Seq using single-cell RNA-Seq signature
#'
#' Deconvolution of bulk RNA-Seq using vector projection method with adjustable
#' compensation for spillover.
#'
#' @param mk object of class 'cellMarkers'. See [cellMarkers()].
#' @param test matrix of bulk RNA-Seq to be deconvoluted. We recommend raw
#'   counts as input, but normalised data can be provided, in which case set
#'   `log = FALSE`.
#' @param log Logical, whether to apply log2 +1 to count data in `test`. Set to
#'   `FALSE` if prenormalised bulk RNA-Seq data is provided.
#' @param count_space Logical, whether deconvolution is performed in count
#'   space (as opposed to log2 space). Signature and test revert to count scale
#'   by 2^ exponentiation during deconvolution.
#' @param comp_amount either a single value from 0-1 for the amount of
#'   compensation or a numeric vector with the same length as the number of cell
#'   subclasses to deconvolute.
#' @param group_comp_amount either a single value from 0-1 for the amount of
#'   compensation for cell group analysis or a numeric vector with the same
#'   length as the number of cell groups to deconvolute.
#' @param weights Optional vector of weights which affects how much each gene in
#'   the gene signature matrix affects the deconvolution.
#' @param adjust_comp logical, whether to optimise `comp_amount` to prevent
#'   negative cell proportion projections.
#' @param use_filter logical, whether to use denoised signature matrix.
#' @param arith_mean logical, whether to use arithmetic means (if available) for
#'   signature matrix. Mainly useful with pseudo-bulk simulation.
#' @param convert_bulk either "ref" to convert bulk RNA-Seq to scRNA-Seq scaling
#'   using reference data or "qqmap" using quantile mapping of the bulk to
#'   scRNA-Seq datasets, or "none" (or `FALSE`) for no conversion.
#' @param plot_comp logical, whether to analyse compensation values across
#'   subclasses.
#' @param IRW Logical, enables iterative reweighting of each gene in the gene
#'   signature matrix based on the absolute deviation of the residuals.
#' @param n_iter Number of iterations.
#' @param delta Regularisation term for the weighting function (to avoid
#'   division by zero).
#' @param bysample Logical, whether `comp_amount` is optimised per sample. This
#'   is a little slower.
#' @param verbose logical, whether to show additional information.
#' @returns A list object of S3 class 'deconv' containing:
#'   \item{call}{the matched call}
#'   \item{mk}{the original 'cellMarkers' class object}
#'   \item{subclass}{list object containing:
#'   \itemize{
#'       \item `output`, the amount of each subclass based purely on project 
#'       gene expression
#'       \item `percent`, the proportion of each subclass scaled as a percentage
#'       so that the total amount across all subclasses adds to 100%
#'       \item `spillover`, the spillover matrix
#'       \item `compensation`, the mixed final compensation matrix which
#'       incorporates `comp_amount`
#'       \item `rawcomp`, the original unadjusted compensation matrix
#'       \item `comp_amount`, the final values for the amount of compensation
#'       across each cell subclass after adjustment to prevent negative values
#'   }}
#'   \item{group}{similar list object to `subclass`, but with results for the 
#'   cell group analysis.}
#'   \item{nest_output}{alternative matrix of cell output results for each 
#'   subclass adjusted so that the cell outputs across subclasses are nested
#'   as a proportion of cell group outputs.}
#'   \item{nest_percent}{alternative matrix of cell proportion results for each 
#'   subclass adjusted so that the percentages across subclasses are nested
#'   within cell group percentages. The total percentage still adds to 100%.}
#'   \item{comp_amount}{original argument `comp_amount`}
#'   \item{comp_check}{optional list element returned when `plot_comp = TRUE`}
#' @seealso [cellMarkers()] [updateMarkers()]
#' @author Myles Lewis
#' @importFrom matrixStats colMins rowMins
#' @importFrom stats optimise
#' @export
#'
deconvolute <- function(mk, test, log = TRUE,
                        count_space = FALSE,
                        comp_amount = 1,
                        group_comp_amount = 0,
                        weights = NULL,
                        adjust_comp = TRUE,
                        use_filter = TRUE,
                        arith_mean = FALSE,
                        convert_bulk = "ref",
                        plot_comp = FALSE,
                        IRW = FALSE,
                        n_iter = 5,
                        delta = 0.01,
                        bysample = FALSE,
                        verbose = FALSE) {
  if (!inherits(mk, "cellMarkers")) stop("Not a 'cellMarkers' class object")
  .call <- match.call()
  
  test <- as.matrix(test)
  
  if (isTRUE(convert_bulk)) convert_bulk <- "ref"
  if (isFALSE(convert_bulk)) convert_bulk <- "none"
  if (convert_bulk == "qqmap") {
    message("Quantile map bulk to sc, ", appendLF = FALSE)
    qqmap <- quantile_map(log2(test +1), mk$genemeans, remove_zeros = TRUE)
  }
  bulk2scfun <- switch(convert_bulk, "ref" = bulk2sc, "qqmap" = qqmap$map)
  
  # group first
  if (!is.null(mk$group_geneset)) {
    cellmat <- if (use_filter) {mk$groupmeans_filtered[mk$group_geneset, ]
    } else mk$groupmeans[mk$group_geneset, ]
    if (!all(mk$group_geneset %in% rownames(test)))
      stop("test is missing some group signature genes")
    logtest <- test[mk$group_geneset, , drop = FALSE]
    if (log) logtest <- log2(logtest +1)
    if (convert_bulk != "none") logtest <- bulk2scfun(logtest)
    gtest <- deconv_adjust(logtest, cellmat, group_comp_amount, weights = NULL,
                           adjust_comp, count_space, bysample)
  } else {
    gtest <- NULL
  }
  
  # cell subclasses
  if (arith_mean) {
    cellmat <- if (use_filter) {mk$genemeans_filtered_ar[mk$geneset, ]
    } else mk$genemeans_ar[mk$geneset, ]
    if (is.null(cellmat)) stop("arithmetic mean not found")
  } else {
    cellmat <- if (use_filter) {mk$genemeans_filtered[mk$geneset, ]
    } else mk$genemeans[mk$geneset, ]
  }
  # cellmat <- sc2bulk(cellmat)
  if (!all(mk$geneset %in% rownames(test)))
    stop("test is missing some signature genes")
  logtest2 <- test[mk$geneset, , drop = FALSE]
  if (log) logtest2 <- log2(logtest2 +1)
  if (convert_bulk != "none") logtest2 <- bulk2scfun(logtest2)
  atest <- deconv_adjust_irw(logtest2, cellmat, comp_amount, weights,
                             adjust_comp, count_space, bysample,
                             IRW, n_iter, delta)
  
  if (verbose) {
    maxsp <- max_spill(atest$spillover)
    message("Max spillover ", format(maxsp, digits = 3))
    # message("Max/min compensation ", format(max_abs(atest$compensation), digits = 3))
  }
  
  # subclass nested within group output/percent
  if (!is.null(gtest)) {
    # subclass output as nested output of groups
    output2 <- lapply(levels(mk$cell_table), function(i) {
      output1 <- gtest$output[, i]
      ind <- mk$cell_table == i
      subclass_i <- atest$output[, ind, drop = FALSE]
      rs <- rowSums(subclass_i)
      subclass_i / rs * output1
    })
    nest_output <- do.call(cbind, output2)
    # subclass percent as nested percent of groups
    pc2 <- lapply(levels(mk$cell_table), function(i) {
      pc1 <- gtest$percent[, i]
      ind <- mk$cell_table == i
      subclass_i <- atest$output[, ind, drop = FALSE]
      rs <- rowSums(subclass_i)
      subclass_i / rs * pc1
    })
    nest_percent <- do.call(cbind, pc2)
  } else nest_output <- nest_percent <- NULL
  
  out <- list(call = .call, mk = mk, subclass = atest, group = gtest,
              nest_output = nest_output, nest_percent = nest_percent,
              comp_amount = comp_amount)
  if (convert_bulk == "qqmap") out$qqmap <- qqmap
  if (plot_comp) {
    message("analysing compensation")
    out$comp_check <- comp_check(logtest2, cellmat, comp_amount,
                                 weights, count_space)
  }
  class(out) <- "deconv"
  out
}

deconv_adjust_irw <- function(test, cellmat, comp_amount, weights,
                              adjust_comp, count_space, bysample, IRW, n_iter,
                              delta, ...) {
  if (!IRW) {
    return(deconv_adjust(test, cellmat, comp_amount, weights,
                         adjust_comp, count_space, bysample, ...))
  }
  fit1 <- fit <- deconv_adjust(test, cellmat, comp_amount, weights,
                               adjust_comp, count_space, bysample = FALSE,
                               resid = TRUE, verbose = FALSE, ...)
  
  for (i in seq_len(n_iter)) {
    abs_dev <- rowMeans(abs(fit$residuals))
    w <- 1 / pmax(abs_dev, delta)
    w <- w / mean(w)
    fit <- try(deconv_adjust(test, cellmat, comp_amount, weights = w,
                             adjust_comp, count_space, bysample = FALSE,
                             resid = TRUE, verbose = (i == n_iter), ...),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      warning(fit)
      fit1$weights <- w
      return(fit1)
    }
    fit1 <- fit
  }
  fit$weights <- w
  fit
}


deconv_adjust <- function(test, cellmat, comp_amount, weights,
                          adjust_comp, count_space, bysample, verbose = TRUE,
                          ...) {
  comp_amount <- rep_len(comp_amount, ncol(cellmat))
  names(comp_amount) <- colnames(cellmat)
  if (bysample) {
    return(deconv_adjust_bysample(test, cellmat, comp_amount, weights,
                                  adjust_comp, count_space))
  }
  
  atest <- deconv(test, cellmat, comp_amount, weights, count_space, ...)
  if (any(atest$output < 0)) {
    if (adjust_comp) {
      minout <- colMins(atest$output)
      w <- which(minout < 0)
      if (verbose) message("optimising compensation (", length(w), ")")
      newcomps <- vapply(w, function(i) {
        f <- function(x) {
          newcomp <- comp_amount
          newcomp[i] <- x
          ntest <- deconv(test, cellmat, comp_amount = newcomp, weights,
                          count_space)
          min(ntest$output[, i], na.rm = TRUE)^2
        }
        if (comp_amount[i] == 0) return(0)
        xmin <- optimise(f, c(0, comp_amount[i]))
        xmin$minimum
      }, numeric(1))
      comp_amount[w] <- newcomps
      atest <- deconv(test, cellmat, comp_amount = comp_amount, weights,
                      count_space, ...)
      # fix floating point errors
      atest$output[atest$output < 0] <- 0
      atest$percent[atest$percent < 0] <- 0
    } else if (verbose) message("negative cell proportion projection detected")
  }
  atest
}


deconv_adjust_bysample <- function(test, cellmat, comp_amount, weights,
                                   adjust_comp, count_space) {
  
  atest <- deconv(test, cellmat, comp_amount, weights, count_space)
  if (any(atest$output < 0)) {
    if (adjust_comp) {
      minout <- rowMins(atest$output)
      wr <- which(minout < 0)
      minout <- colMins(atest$output)
      w <- which(minout < 0)
      message("optimising compensation (", length(w), ")")
      btest <- lapply(wr, function(j) {
        testj <- test[, j, drop = FALSE]
        w <- which(atest$output[j, ] < 0)
        newcomps <- vapply(w, function(i) {
          f <- function(x) {
            newcomp <- comp_amount
            newcomp[i] <- x
            ntest <- quick_deconv(testj, cellmat, comp_amount = newcomp,
                                  weights, count_space)
            ntest[, i]^2
          }
          if (comp_amount[i] == 0) return(0)
          xmin <- optimise(f, c(0, comp_amount[i]))
          xmin$minimum
        }, numeric(1))
        comp_amount2 <- comp_amount
        comp_amount2[w] <- newcomps
        output <- quick_deconv(testj, cellmat, comp_amount = comp_amount2,
                               weights, count_space)
        list(output = output, comp_amount = comp_amount2)
      })
      newrows <- lapply(btest, function(x) x$output)
      newrows <- do.call(rbind, newrows)
      atest$output[wr, ] <- newrows
      comps <- lapply(btest, function(x) x$comp_amount)
      comps <- do.call(rbind, comps)
      atest$comp_amount <- rbind(atest$comp_amount, comps)
      atest$percent <- atest$output / rowSums(atest$output) * 100
      # fix floating point errors
      atest$output[atest$output < 0] <- 0
      atest$percent[atest$percent < 0] <- 0
    } else message("negative cell proportion projection detected")
  }
  atest
}


deconv <- function(test, cellmat, comp_amount = 0, weights = NULL,
                   count_space = FALSE, resid = FALSE) {
  if (!(length(comp_amount) %in% c(1, ncol(cellmat))))
    stop('comp_amount must be either single number or vector of length matching cellmat cols')
  if (!identical(rownames(test), rownames(cellmat)))
    stop('test and cell matrices must have same genes (rownames)')
  if (count_space) {
    test <- 2^test -1
    cellmat <- 2^cellmat -1
  }
  m_itself <- dotprod(cellmat, cellmat, weights)
  rawcomp <- solve(m_itself)
  mixcomp <- solve(m_itself, t(comp_amount * diag(nrow(m_itself)) + (1-comp_amount) * t(m_itself)))
  output <- dotprod(test, cellmat, weights) %*% mixcomp
  percent <- output / rowSums(output) * 100
  
  out <- list(output = output, percent = percent, spillover = m_itself,
              compensation = mixcomp, rawcomp = rawcomp,
              comp_amount = comp_amount)
  if (resid) out$residuals <- residuals_deconv(test, cellmat, output)
  
  out
}


quick_deconv <- function(test, cellmat, comp_amount, weights, count_space) {
  if (count_space) {
    test <- 2^test -1
    cellmat <- 2^cellmat -1
  }
  m_itself <- dotprod(cellmat, cellmat, weights)
  mixcomp <- solve(m_itself, t(comp_amount * diag(nrow(m_itself)) + (1-comp_amount) * t(m_itself)))
  dotprod(test, cellmat, weights) %*% mixcomp
}

approxfun.matrix <- function(x, FUN) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    out <- FUN(as.vector(x))
    out <- matrix(out, nrow = nrow(x), dimnames = dimnames(x))
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


comp_check <- function(test, cellmat, comp_amount, weights,
                       count_space) {
  comp_amount <- rep_len(comp_amount, ncol(cellmat))
  names(comp_amount) <- colnames(cellmat)
  px <- seq(0, 1, 0.05)
  
  out <- lapply(seq_len(ncol(cellmat)), function(i) {
    newcomp <- comp_amount
    vapply(px, function(ci) {
      newcomp[i] <- ci
      ntest <- deconv(test, cellmat, newcomp, weights, count_space)
      min(ntest$output[, i], na.rm = TRUE)
    }, numeric(1))
  })
  names(out) <- colnames(cellmat)
  out$px <- px
  out
}
