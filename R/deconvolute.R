
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
#' @param weight_method Optional. Choices include "equal" in which case weights
#'   are calculated so that each gene has equal weighting in the vector
#'   projection; "none"; or "irw" which enables iterative reweighting of genes
#'   based on residuals (see details). Setting this overrules any vector
#'   supplied by `weights`.
#' @param adjust_comp logical, whether to optimise `comp_amount` to prevent
#'   negative cell proportion projections.
#' @param use_filter logical, whether to use denoised signature matrix.
#' @param arith_mean logical, whether to use arithmetic means (if available) for
#'   signature matrix. Mainly useful with pseudo-bulk simulation.
#' @param convert_bulk either "ref" to convert bulk RNA-Seq to scRNA-Seq scaling
#'   using reference data or "qqmap" using quantile mapping of the bulk to
#'   scRNA-Seq datasets, or "none" (or `FALSE`) for no conversion.
#' @param check_comp logical, whether to analyse compensation values across
#'   subclasses.
#' @param npass Number of passes. If `npass` set to 2 or more this activates
#'   removal of genes with excess variance of the residuals.
#' @param var_cutoff Cutoff as Z-score for removing genes with high variance of
#'   residuals. Variances of residuals for each gene are first log2 transformed
#'   and then Z-score standardised.
#' @param verbose logical, whether to show messages.
#' @param cores Number of cores for parallelisation via `parallel::mclapply()`.
#' @details
#' Equal weighting of genes by setting `weight_method = "equal"` can help
#' devolution of subclusters whose signature genes have low expression. It is
#' enabled by default.
#' 
#' Multipass deconvolution can be activated by setting `npass` to 2 or higher.
#' This calculates the variance of the residuals across samples for each gene.
#' Genes whose variance of residuals are outliers based on Z-score
#' standardisation are removed during successive passes.
#' 
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
#'       \item `residuals`, residuals, that is gene expression minus fitted 
#'       values
#'       \item `resvar`, \eqn{s^2} the estimate of the gene expression variance 
#'       for each sample
#'       \item `se`, standard errors of cell counts
#'       \item `SE_method`, standard error method
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
#'   \item{comp_check}{optional list element returned when `check_comp = TRUE`}
#' @seealso [cellMarkers()] [updateMarkers()]
#' @author Myles Lewis
#' @importFrom matrixStats colMins rowMins
#' @importFrom stats optimise
#' @export
#'
deconvolute <- function(mk, test, log = TRUE,
                        count_space = TRUE,
                        comp_amount = 1,
                        group_comp_amount = 0,
                        weights = NULL,
                        weight_method = "equal",
                        adjust_comp = TRUE,
                        use_filter = TRUE,
                        arith_mean = FALSE,
                        convert_bulk = FALSE,
                        check_comp = FALSE,
                        npass = 1,
                        var_cutoff = 4,
                        verbose = TRUE, cores = 1L) {
  if (!inherits(mk, "cellMarkers")) stop("Not a 'cellMarkers' class object")
  .call <- match.call()
  weight_method <- match.arg(weight_method, c("none", "equal", "2pass"))
  test <- as.matrix(test)
  
  if (isTRUE(convert_bulk)) convert_bulk <- "ref"
  if (isFALSE(convert_bulk)) convert_bulk <- "none"
  if (convert_bulk == "qqmap") {
    if (verbose) cat("Quantile map bulk to sc, ")
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
                           adjust_comp, count_space, weight_method,
                           verbose = verbose, resid = FALSE)
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
  atest <- deconv_multipass(logtest2, cellmat, comp_amount, weights,
                            weight_method, adjust_comp, count_space,
                            var_cutoff, npass,
                            verbose, cores)
  
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
    nest_output[is.na(nest_output)] <- 0
    # subclass percent as nested percent of groups
    pc2 <- lapply(levels(mk$cell_table), function(i) {
      pc1 <- gtest$percent[, i]
      ind <- mk$cell_table == i
      subclass_i <- atest$output[, ind, drop = FALSE]
      rs <- rowSums(subclass_i)
      subclass_i / rs * pc1
    })
    nest_percent <- do.call(cbind, pc2)
    nest_percent[is.na(nest_percent)] <- 0
  } else nest_output <- nest_percent <- NULL
  
  out <- list(call = .call, mk = mk, subclass = atest, group = gtest,
              nest_output = nest_output, nest_percent = nest_percent,
              comp_amount = comp_amount)
  if (convert_bulk == "qqmap") out$qqmap <- qqmap
  if (check_comp) {
    if (verbose) cat("analysing compensation\n")
    out$comp_check <- comp_check(logtest2, cellmat, comp_amount,
                                 weights, count_space)
  }
  class(out) <- "deconv"
  out
}


deconv_multipass <- function(test, cellmat, comp_amount, weights, weight_method,
                             adjust_comp, count_space,
                             var_cutoff, npass,
                             verbose, cores) {
  fit <- deconv_adjust(test, cellmat, comp_amount, weights,
                       adjust_comp, count_space, weight_method,
                       cores, verbose)
  var.e <- if (count_space) log2(fit$var.e +1) else fit$var.e
  scale.var.e <- scale(var.e)[, 1]
  bigvar <- scale.var.e > var_cutoff
  if (verbose && npass == 1 && any(bigvar)) {
    cat("Detected genes with extreme residuals:",
        paste(names(var.e)[bigvar], collapse = ", "), "\n")
  }
  i <- 1
  while (any(bigvar) & i < npass) {
    i <- i +1
    if (verbose) cat("Pass", i, "- removed", paste(names(var.e)[bigvar],
                                                   collapse = ", "), "\n")
    # print(scale.var.e[bigvar], digits = 3)
    test <- test[!bigvar, , drop = FALSE]
    cellmat <- cellmat[!bigvar, , drop = FALSE]
    weights <- weights[!bigvar]
    fit <- deconv_adjust(test, cellmat, comp_amount, weights,
                         adjust_comp, count_space, weight_method,
                         cores, verbose)
    var.e <- if (count_space) log2(fit$var.e +1) else fit$var.e
    scale.var.e <- scale(var.e)[, 1]
    bigvar <- scale.var.e > var_cutoff
  }
  
  fit
}


deconv_adjust <- function(test, cellmat, comp_amount, weights,
                          adjust_comp, count_space,
                          weight_method = "", cores = 1L, verbose = TRUE,
                          resid = TRUE) {
  comp_amount <- rep_len(comp_amount, ncol(cellmat))
  names(comp_amount) <- colnames(cellmat)
  if (!identical(rownames(test), rownames(cellmat)))
    stop('test and cell matrices must have same genes (rownames)')
  if (count_space) {
    test <- 2^test -1
    cellmat <- 2^cellmat -1
  }
  if (!is.null(weights) && length(weights) != nrow(cellmat))
    stop("incorrect weights length")
  if (weight_method == "equal") weights <- equalweight(cellmat)
  if (any(nok <- weights == 0)) {
    test <- test[!nok, , drop = FALSE]
    cellmat <- cellmat[!nok, , drop = FALSE]
    weights <- weights[!nok] 
  }
  
  atest <- deconv(test, cellmat, comp_amount, weights)
  if (any(atest$output < 0)) {
    if (adjust_comp) {
      minout <- colMins(atest$output)
      w <- which(minout < 0)
      if (verbose) cat(paste0("optimising compensation (", length(w), ")\n"))
      
      newcomps <- pmclapply(seq_along(w), function(i) {
        wi <- w[i]
        f <- function(x) {
          newcomp <- comp_amount
          newcomp[wi] <- x
          ntest <- quick_deconv(test, cellmat, comp_amount = newcomp, weights)
          min(ntest[, wi], na.rm = TRUE)^2
        }
        if (comp_amount[wi] == 0) return(0)
        xmin <- optimise(f, c(0, comp_amount[wi]))
        xmin$minimum
      }, progress = verbose, mc.cores = cores)
      comp_amount[w] <- unlist(newcomps)
      atest <- deconv(test, cellmat, comp_amount = comp_amount, weights)
      # fix floating point errors
      atest$output[atest$output < 0] <- 0
      atest$percent[atest$percent < 0] <- 0
    } else if (verbose) cat("negative cell proportion projection detected\n")
  }
  if (resid) {
    atest$residuals <- r <- residuals_deconv(test, cellmat, atest$output)
    if (!is.null(weights)) {
      # adjust residuals & X by gene weights
      r <- r * weights
      cellmat <- cellmat * weights
    }
    rss <- colSums(r^2)
    rdf <- nrow(r) - ncol(cellmat)
    atest$resvar <- resvar <- rss/rdf
    Lv <- colSums(cellmat^2)
    diag_XTX <- diag(atest$compensation) / Lv
    atest$se1 <- sqrt(resvar %*% t(diag_XTX))
    # based on van Wieringen
    iXTX <- atest$compensation / Lv
    XTX <- crossprod(cellmat)
    diag2 <- diag(iXTX %*% XTX %*% t(iXTX))
    atest$se2 <- sqrt(resvar %*% t(diag2))
    # heteroscedasticity-consistent SE = HC0
    atest$se3 <- t(apply(r, 2, function(i) {
      XTXse <- crossprod(cellmat, i^2 * cellmat)
      sqrt(diag(iXTX %*% XTXse %*% t(iXTX)))
    }))
    # empirical Bayes estimate employing residuals row variance
    var.e <- matrixStats::rowVars(r)
    atest$se4 <- t(apply(r, 2, function(i) {
      vbsq <- (i^2 + var.e) /2
      XTXse <- crossprod(cellmat, vbsq * cellmat)
      sqrt(diag(iXTX %*% XTXse %*% t(iXTX)))
    }))
    # deploy residuals row variance
    XTXse <- crossprod(cellmat, var.e * cellmat)
    atest$se5 <-  sqrt(diag(iXTX %*% XTXse %*% t(iXTX)))
    atest$var.e <- var.e
    atest$mean.e <- rowMeans(r)
  }
  atest
}


deconv <- function(test, cellmat, comp_amount, weights) {
  m_itself <- dotprod(cellmat, cellmat, weights)
  rawcomp <- solve(m_itself)
  mixcomp <- solve(m_itself, t(comp_amount * diag(nrow(m_itself)) + (1-comp_amount) * t(m_itself)))
  output <- dotprod(test, cellmat, weights) %*% mixcomp
  percent <- output / rowSums(output) * 100
  
  list(output = output, percent = percent, spillover = m_itself,
       compensation = mixcomp, rawcomp = rawcomp,
       comp_amount = comp_amount)
}


quick_deconv <- function(test, cellmat, comp_amount, weights) {
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
  if (count_space) {
    test <- 2^test -1
    cellmat <- 2^cellmat -1
  }
  px <- seq(0, 1, 0.05)
  
  out <- lapply(seq_len(ncol(cellmat)), function(i) {
    newcomp <- comp_amount
    vapply(px, function(ci) {
      newcomp[i] <- ci
      ntest <- quick_deconv(test, cellmat, newcomp, weights)
      min(ntest[, i], na.rm = TRUE)
    }, numeric(1))
  })
  names(out) <- colnames(cellmat)
  out$px <- px
  out
}
