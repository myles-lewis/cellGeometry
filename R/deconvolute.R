
#' Deconvolute bulk RNA-Seq using single-cell RNA-Seq signature
#'
#' Deconvolution of bulk RNA-Seq using vector projection method with adjustable
#' compensation for spillover.
#'
#' @param mk object of class 'cellMarkers'. See [cellMarkers()].
#' @param test matrix of bulk RNA-Seq to be deconvoluted with genes in rows and
#'   samples in columns. We recommend raw counts as input, but normalised data
#'   can be provided, in which case set `logged_bulk = TRUE`.
#' @param logged_bulk Logical, whether log2 transformed bulk RNA-Seq data is
#'   used as input in `test`.
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
#' @param weight_method Optional. Choices include "none" or "equal" in which
#'   gene weights are calculated so that each gene has equal weighting in the
#'   vector projection; "equal" overrules any vector supplied by `weights`.
#' @param adjust_comp logical, whether to optimise `comp_amount` to prevent
#'   negative cell proportion projections.
#' @param use_filter logical, whether to use denoised signature matrix.
#' @param arith_mean logical, whether to use arithmetic means (if available) for
#'   signature matrix. Mainly useful with pseudo-bulk simulation.
#' @param convert_bulk either "ref" to convert bulk RNA-Seq to scRNA-Seq scaling
#'   using reference data or "qqmap" using quantile mapping of the bulk to
#'   scRNA-Seq datasets, or "none" (or `FALSE`) for no conversion.
#' @param check_comp logical, whether to analyse compensation values across
#'   subclasses. See [plot_comp()].
#' @param npass Number of passes. If `npass` set to 2 or more this activates
#'   removal of genes with excess variance of the residuals.
#' @param outlier_method Method for identifying outlying genes. Options are to
#'   use the variance of the residuals for each genes, Cook's distance or
#'   absolute Studentized residuals (see details).
#' @param outlier_cutoff Cutoff for removing genes which are outliers based on
#'   method selected by `outlier_method`.
#' @param outlier_quantile Controls quantile for the cutoff for identifying
#'   outliers for `outlier_method = "cook"` or `"rstudent"`.
#' @param verbose logical, whether to show messages.
#' @param cores Number of cores for parallelisation via `parallel::mclapply()`.
#' @details
#' Equal weighting of genes by setting `weight_method = "equal"` can help
#' devolution of subclusters whose signature genes have low expression. It is
#' enabled by default.
#' 
#' If a normalised (i.e. logged) bulk matrix is provided instead of raw counts,
#' then it is important that zero expression is true zero. For this reason we do
#' not recommend use of VST (variance stabilised transformed counts) which has a
#' variable offset.
#' 
#' Multipass deconvolution can be activated by setting `npass` to 2 or higher.
#' This is designed to remove genes which behave inconsistently due to noise in
#' either the sc or bulk datasets, which is increasingly likely if you have
#' larger signature geneset, i.e. if `nsubclass` is large. Or you may receive a
#' warning message "Detected genes with extreme residuals". Three methods are
#' available for identifying outlier genes (i.e. whose residuals are too noisy)
#' controlled by `outlier_method`:
#' - `var.e`, this calculates the variance of the residuals across samples for 
#' each gene. Genes whose variance of residuals are outliers based on Z-score
#' standardisation are removed during successive passes.
#' - `cooks`, this considers the deconvolution as if it were a regression and 
#' applies Cook's distance to the residuals and the hat matrix. This seems to be
#' the most stringent method (removes fewest genes).
#' - `rstudent`, externally Studentized residuals are used.
#' 
#' The cutoff specified by `outlier_cutoff` which is used to determine which
#' genes are outliers is very sensitive to the outlier method. With `var.e` the
#' variances are Z-score scaled. With Cook's distance it is typical to consider
#' a value of >1 as fairly strong indication of an outlier, while 0.5 is
#' considered a possible outlier. With Studentized residuals, these are expected
#' to be on a t distribution scale. However, since gene expression itself does
#' not derive from a normal distribution, the errors and residuals are not
#' normally distributed either, which probably explains the need for a very high
#' cut-off. In practice the choice of settings seems to be dataset dependent.
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
#'       \item `var.e`, variance of weighted residuals for each gene
#'       \item `weights`, vector of weights
#'       \item `resvar`, \eqn{s^2} the estimate of the gene expression variance 
#'       for each sample
#'       \item `se`, standard errors of cell counts
#'       \item `hat`, diagonal elements of the hat matrix
#'       \item `removed`, vector of outlying genes removed during successive 
#'       passes
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
#' @seealso [cellMarkers()] [updateMarkers()] [rstudent.deconv()]
#'   [cooks.distance.deconv()]
#' @author Myles Lewis
#' @importFrom matrixStats colMins rowQuantiles rowVars
#' @importFrom stats optimise
#' @export
#'
deconvolute <- function(mk, test,
                        logged_bulk = FALSE,
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
                        outlier_method = c("var.e", "cooks", "rstudent"),
                        outlier_cutoff = switch(outlier_method, var.e = 4,
                                                cooks = 1, rstudent = 10),
                        outlier_quantile = 0.9,
                        verbose = TRUE, cores = 1L) {
  if (!inherits(mk, "cellMarkers")) stop("Not a 'cellMarkers' class object")
  .call <- match.call()
  weight_method <- match.arg(weight_method, c("none", "equal"))
  outlier_method <- match.arg(outlier_method)
  test <- as.matrix(test)
  if (any(test < 0)) stop("`test` contains negative values")
  
  if (isTRUE(convert_bulk)) convert_bulk <- "ref"
  if (isFALSE(convert_bulk)) convert_bulk <- "none"
  if (convert_bulk == "qqmap") {
    if (verbose) message("Quantile map bulk to sc, ", appendLF = FALSE)
    logtest0 <- if (logged_bulk) test else log2(test +1)
    qqmap <- quantile_map(logtest0, mk$genemeans, remove_zeros = TRUE)
  }
  bulk2scfun <- switch(convert_bulk, "ref" = bulk2sc, "qqmap" = qqmap$map)
  
  # group first
  if (!is.null(mk$group_geneset)) {
    cellmat <- if (use_filter) {mk$groupmeans_filtered[mk$group_geneset, ]
    } else mk$groupmeans[mk$group_geneset, ]
    if (!all(mk$group_geneset %in% rownames(test)))
      stop("`test` is missing some group signature genes")
    logtest <- test[mk$group_geneset, , drop = FALSE]
    if (!logged_bulk) logtest <- log2(logtest +1)
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
    stop("`test` is missing some signature genes")
  logtest2 <- test[mk$geneset, , drop = FALSE]
  if (!logged_bulk) logtest2 <- log2(logtest2 +1)
  if (convert_bulk != "none") logtest2 <- bulk2scfun(logtest2)
  atest <- deconv_multipass(logtest2, cellmat, comp_amount, weights,
                            weight_method, adjust_comp, count_space, npass,
                            outlier_method, outlier_cutoff, outlier_quantile,
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
    if (verbose) message("analysing compensation")
    out$comp_check <- comp_check(logtest2, cellmat, comp_amount,
                                 weights, weight_method, count_space, cores)
  }
  class(out) <- "deconv"
  out
}


deconv_multipass <- function(test, cellmat, comp_amount, weights, weight_method,
                             adjust_comp, count_space, npass,
                             outlier_method, outlier_cutoff, outlier_quantile,
                             verbose, cores) {
  fit <- deconv_adjust(test, cellmat, comp_amount, weights,
                       adjust_comp, count_space, weight_method,
                       cores, verbose)
  metric <- outlier_metric(fit, outlier_method, outlier_quantile, count_space)
  outlier <- metric > outlier_cutoff
  if (verbose && npass == 1 && any(outlier)) {
    nm <- names(sort(metric[outlier], decreasing = TRUE))
    if (length(nm) > 20) nm <- c(nm[1:20], "...")
    message("Detected genes with extreme residuals: ",
            paste(nm, collapse = ", "))
  }
  i <- 1
  removed <- NULL
  while (any(outlier) & i < npass) {
    i <- i +1
    remove_genes <- sort(metric[outlier], decreasing = TRUE)
    removed <- c(removed, remove_genes)
    rem_names <- names(remove_genes)
    if (length(remove_genes) > 20) {
      rem_names <- c(rem_names[1:20],
                    paste0("... [", length(remove_genes), " genes]"))
    }
    if (verbose) message("Pass ", i, " - removed ",
                         paste(rem_names, collapse = ", "))
    test <- test[!outlier, , drop = FALSE]
    cellmat <- cellmat[!outlier, , drop = FALSE]
    weights <- weights[!outlier]
    if (nrow(cellmat) < ncol(cellmat)) stop("insufficient genes")
    fit <- deconv_adjust(test, cellmat, comp_amount, weights,
                         adjust_comp, count_space, weight_method,
                         cores, verbose)
    metric <- outlier_metric(fit, outlier_method, outlier_quantile, count_space)
    outlier <- metric > outlier_cutoff
  }
  fit$removed <- removed
  fit
}


outlier_metric <- function(fit, outlier_method, outlier_quantile, count_space) {
  metric <- switch(outlier_method,
                   cooks = cooks_distance_fit(fit),
                   rstudent = abs(rstudent_fit(fit)),
                   fit$var.e)
  if (outlier_method == "var.e") {
    if (count_space) metric <- log2(metric +1)
    return(scale(metric)[, 1])
  }
  metric[is.nan(metric)] <- 0
  rowQuantiles(metric, probs = outlier_quantile)
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
  oldcellmat <- cellmat
  oldtest <- test
  if (!is.null(weights)) {
    cellmat <- cellmat * weights
    test <- test * weights
  }
  
  m_itself <- dotprod(cellmat, cellmat)
  rawcomp <- solve(m_itself)
  atest <- deconv(test, cellmat, comp_amount, m_itself, rawcomp)
  if (any(atest$output < 0)) {
    if (adjust_comp) {
      minout <- colMins(atest$output)
      w <- which(minout < 0)
      if (verbose) message("optimising compensation (", length(w), ")")
      
      newcomps <- pmclapply(w, function(wi) {
        f <- function(x) {
          newcomp <- comp_amount
          newcomp[wi] <- x
          ntest <- quick_deconv(test, cellmat, newcomp, m_itself, rawcomp, wi)
          min(ntest, na.rm = TRUE)^2
        }
        if (comp_amount[wi] == 0) return(0)
        xmin <- optimise(f, c(0, comp_amount[wi]))
        xmin$minimum
      }, progress = verbose, mc.cores = cores)
      comp_amount[w] <- unlist(newcomps)
      atest <- deconv(test, cellmat, comp_amount, m_itself, rawcomp)
      # fix floating point errors
      if (any(z <- atest$output < 0)) {
        attr(atest$output, "min") <- min(atest$output)
        atest$output[z] <- 0
        atest$percent <- atest$output / rowSums(atest$output) * 100
      }
    } else if (verbose) message("negative cell proportion projection detected")
  }
  if (resid) {
    X <- cellmat
    atest$residuals <- r <- residuals_deconv(oldtest, oldcellmat, atest$output)
    # adjust residuals & X by gene weights
    if (!is.null(weights)) r <- r * weights
    # residuals row variance
    atest$var.e <- var.e <- rowVars(r)
    atest$weights <- weights
    rss <- colSums(r^2)
    rdf <- nrow(r) - ncol(X)
    atest$resvar <- rss/rdf
    # residuals row variance
    Lv <- colSums(X^2)
    iXTX <- atest$compensation / Lv
    XTXse <- crossprod(X, var.e * X)
    # var(beta) = (X' X)^-1 (X diag(e^2) X') (X' X)^-1
    atest$se <- sqrt(diag(iXTX %*% XTXse %*% t(iXTX)))
    # calculate hat matrix
    # X (X' X)^-1 X'
    atest$hat <- diag(X %*% iXTX %*% t(X))
  }
  atest
}


deconv <- function(test, cellmat, comp_amount, m_itself, rawcomp) {
  endcomp <- comp_amount * diag(nrow(m_itself)) + (1-comp_amount) * t(m_itself)
  mixcomp <- tcrossprod(rawcomp, endcomp)
  output <- dotprod(test, cellmat) %*% mixcomp
  percent <- output / rowSums(output) * 100
  
  list(output = output, percent = percent, spillover = m_itself,
       compensation = mixcomp, rawcomp = rawcomp,
       comp_amount = comp_amount)
}

# single column
quick_deconv <- function(test, cellmat, comp_amount, m_itself, rawcomp, wi) {
  endcomp <- (1-comp_amount[wi]) * m_itself[, wi]
  endcomp[wi] <- endcomp[wi] + comp_amount[wi]
  mixcomp <- rawcomp %*% endcomp
  dotprod(test, cellmat) %*% mixcomp
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


comp_check <- function(test, cellmat, comp_amount, weights, weight_method,
                       count_space, cores) {
  comp_amount <- rep_len(comp_amount, ncol(cellmat))
  names(comp_amount) <- colnames(cellmat)
  if (count_space) {
    test <- 2^test -1
    cellmat <- 2^cellmat -1
  }
  if (weight_method == "equal") weights <- equalweight(cellmat)
  if (any(nok <- weights == 0)) {
    test <- test[!nok, , drop = FALSE]
    cellmat <- cellmat[!nok, , drop = FALSE]
    weights <- weights[!nok]
  }
  if (!is.null(weights)) {
    cellmat <- cellmat * weights
    test <- test * weights
  }
  px <- seq(0, 1, 0.05)
  
  m_itself <- dotprod(cellmat, cellmat)
  rawcomp <- solve(m_itself)
  out <- mclapply(seq_len(ncol(cellmat)), function(i) {
    newcomp <- comp_amount
    vapply(px, function(ci) {
      newcomp[i] <- ci
      ntest <- quick_deconv(test, cellmat, newcomp, m_itself, rawcomp, i)
      min(ntest, na.rm = TRUE)
    }, numeric(1))
  }, mc.cores = cores)
  names(out) <- colnames(cellmat)
  out$px <- px
  out
}
