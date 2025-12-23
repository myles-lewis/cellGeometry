
cv_deconv <- function(test, cellmat, comp_amount, weights,
                      adjust_comp, count_space,
                      weight_method = "", lambda, cores = 1L, verbose = TRUE,
                      resid = TRUE) {
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
  
  # CV
  if (verbose) message("CV lambda")
  nfolds <- 10L
  lam_set <- lambda %||% c(0, 10^seq(-4, 0, length = 25))
  folds <- sample(rep_len(seq_len(nfolds), nrow(cellmat)))
  # for each fold do lambda set
  mse <- pmclapply(seq_len(nfolds), function(i) {
    train <- which(folds != i)
    deconv_lambda_set(test, cellmat, weights, comp_amount, adjust_comp, lam_set,
                      train)
  }, progress = verbose, mc.cores = cores)
  mse <- do.call(cbind, mse)
  cvm <- rowMeans(mse)
  cvsd <- matrixStats::rowSds(mse) # / sqrt(nfolds)
  mmse <- cbind(lambda = lam_set, cvm, cvsd)
  lambda.min <- lam_set[which.min(cvm)]
  cv1se <- min(cvm) + cvsd[which.min(cvm)]
  lambda.1se <- lam_set[min(which(cvm < cv1se))]
  list(mmse = mmse, lambda.min = lambda.min, lambda.1se = lambda.1se)
}


deconv_adjust_core <- function(test.cellmat, comp_amount, adjust_comp,
                               m_itself, lambda,
                               oldtest_test, oldcellmat_test, weights_test) {
  m_itself <- m_itself + diag(nrow(m_itself)) * lambda
  rawcomp <- solve(m_itself)
  atest <- deconv(test.cellmat, comp_amount, m_itself, rawcomp)
  if (any(atest$output < 0)) {
    if (adjust_comp) {
      minout <- colMins(atest$output)
      w <- which(minout < 0)
      newcomps <- lapply(w, function(wi) {
        if (comp_amount[wi] == 0) return(0)
        f <- function(x) {
          newcomp <- comp_amount
          newcomp[wi] <- x
          ntest <- quick_deconv(test.cellmat, newcomp, m_itself, rawcomp, wi)
          min(ntest, na.rm = TRUE)^2
        }
        xmin <- optimise(f, c(0, comp_amount[wi]))
        xmin$minimum
      })
      comp_amount[w] <- unlist(newcomps)
      atest <- deconv(test.cellmat, comp_amount, m_itself, rawcomp)
    }
  }
  
  # resvar on holdout
  r <- residuals_deconv(oldtest_test, oldcellmat_test, atest$output)
  # adjust residuals by gene weights
  if (!is.null(weights_test)) r <- r * weights_test
  # mse across all holdout bulk samples
  mean(r^2)
}


deconv_lambda_set <- function(test, cellmat, weights, comp_amount, adjust_comp,
                              lam_set, train) {
  oldcellmat <- cellmat
  oldtest <- test
  if (!is.null(weights)) {
    cellmat <- cellmat * weights
    test <- test * weights
  }
  oldtest_test <- oldtest[-train, , drop = FALSE]
  oldcellmat_test <- oldcellmat[-train, , drop = FALSE]
  
  # train-test
  m_itself <- dotprod(cellmat[train, ], cellmat[train, ])
  vapply(lam_set, function(i) {
    test.cellmat <- dotprod(test[train, ], cellmat[train, ])
    deconv_adjust_core(test.cellmat, comp_amount,
                       adjust_comp, m_itself, lambda = i,
                       oldtest_test, oldcellmat_test, weights[-train])
  }, numeric(1))
}
# returns vector of mse for each lambda


#' Plot deconvolution lambda cross-validation curve
#' 
#' Plots the cross-validation curve, and upper and lower standard error curves
#' as a function of the ridge parameter lambda.
#' 
#' @param fit Object of class 'deconv'
#' @param ... Optional arguments passed to `plot()`
#' @returns No return value. Plots the lambda CV curve using base graphics.
#' @importFrom graphics arrows
#' @export
plot_cv <- function(fit, ...) {
  if (!inherits(fit, "deconv")) stop("`fit` is not a 'deconv' class object")
  if (is.null(fit$subclass$cv)) stop("no cross-validation")
  xm <- fit$subclass$cv$mmse[, "lambda"]
  ok <- xm > 0
  xm <- xm[ok]
  ym <- fit$subclass$cv$mmse[ok, "cvm"]
  ysd <- fit$subclass$cv$mmse[ok, "cvsd"]
  ylo <- ym - ysd
  yhi <- ym + ysd
  
  args <- list(x = xm, y = ym,
               pch = 21, bg = "white",
               ylim = c(min(ylo), max(yhi)),
               log = "xy", xlab = "lambda", ylab = "MSE",
               panel.first = quote({
                 arrows(xm, ylo, xm, yhi, angle = 90, code = 3, length = 0.02)
               }))
  new.args <- list(...)
  if (length(new.args)) args[names(new.args)] <- new.args
  do.call(plot, args)
}
