
#' Regression Deletion Diagnostics
#' 
#' Functions for computing regression diagnostics including standardised or
#' Studentized residuals as well as Cook's distance.
#' 
#' Residuals are first adjusted for gene weights (if used). `rstandard` and
#' `rstudent` give standardized and Studentized residuals respectively.
#' Standardised residuals are calculated based on the hat matrix:
#' \deqn{H = X (X^T X)^{-1} X^T}
#' Leverage \eqn{h_{ii} = diag(H)} is used to standardise the residuals:
#' \deqn{t_i = \cfrac{\hat{\varepsilon_i}}{\hat{\sigma} \sqrt{1 - h_{ii}}}}
#' Studentized residuals are calculated based on excluding the \eqn{i} th case.
#' Note this corresponds to refitting the regression, but without recomputing
#' the non-negative compensation matrix. Cook's distance is calculated as:
#' \deqn{D_i = \cfrac{e_i^2}{ps^2} \left[\cfrac{h_{ii}}{(1 - h_{ii})^2} \right]}
#' where \eqn{p} is the number of predictors (cell subclasses) and \eqn{s^2} is
#' the mean squared error. In this model the intercept is not included.
#' 
#' @param model 'deconv' class object
#' @param ... retained for class compatibility
#' @returns Matrix of adjusted residuals or Cook's distance.
#' @seealso [stats::influence.measures()]
#' @importFrom stats rstudent rstandard cooks.distance
#' @export

rstudent.deconv <- function(model, ...) {
  r <- model$subclass$residuals
  weights <- model$subclass$weights
  if (!is.null(weights)) r <- r * weights
  hat <- hat(model)
  rdf <- nrow(r) - ncol(model$subclass$compensation)
  rss <- colSums(r^2)
  mse <- rss / nrow(r)
  std_res <- t(t(r / sqrt(1 - hat)) / sqrt(mse))
  # externally studentized
  std_res * sqrt((rdf - 1) / (rdf - std_res^2))
}

#' @rdname rstudent.deconv
#' @export
rstandard.deconv <- function(model, ...) {
  r <- model$subclass$residuals
  weights <- model$subclass$weights
  if (!is.null(weights)) r <- r * weights
  hat <- hat(model)
  rdf <- nrow(r) - ncol(model$subclass$compensation)
  rss <- colSums(r^2)
  mse <- rss / nrow(r)
  t(t(r / sqrt(1 - hat)) / sqrt(mse))
}


#' @rdname rstudent.deconv
#' @export
cooks.distance.deconv <- function(model, ...) {
  r <- model$subclass$residuals
  weights <- model$subclass$weights
  if (!is.null(weights)) r <- r * weights
  hat <- hat(model)
  p <- ncol(model$subclass$compensation)
  rdf <- nrow(r) - p
  rss <- colSums(r^2)
  mse <- rss / rdf
  (r / (1 - hat))^2 * outer(hat, mse, "/") / p
}


## functions for fit object
rstudent_fit <- function(fit) {
  r <- fit$residuals
  weights <- fit$weights
  if (!is.null(weights)) r <- r * weights
  hat <- hat_fit(fit)
  rdf <- nrow(r) - ncol(fit$compensation)
  rss <- colSums(r^2)
  mse <- rss / nrow(r)
  std_res <- t(t(r / sqrt(1 - hat)) / sqrt(mse))
  # externally studentized
  std_res * sqrt((rdf - 1) / (rdf - std_res^2))
}

cooks_distance_fit <- function(fit) {
  r <- fit$residuals
  weights <- fit$weights
  if (!is.null(weights)) r <- r * weights
  hat <- hat_fit(fit)
  p <- ncol(fit$compensation)
  rdf <- nrow(r) - p
  rss <- colSums(r^2)
  mse <- rss / rdf
  (r / (1 - hat))^2 * outer(hat, mse, "/") / p
}


hat <- function(model) {
  X <- model$subclass$X
  Lv <- colSums(X^2)
  iXTX <- model$subclass$compensation / Lv
  diag(X %*% iXTX %*% t(X))
}

hat_fit <- function(fit) {
  X <- fit$X
  Lv <- colSums(X^2)
  iXTX <- fit$compensation / Lv
  diag(X %*% iXTX %*% t(X))
}


#' Extract Deconvolution Residuals
#' 
#' Extracts residuals from a deconvolution model. As the model uses a reduced
#' signature gene set for deconvolution, in order to extract residuals for all
#' genes, these need to recalculated by supplying the bulk count matrix `test`.
#' 
#' @param object a 'deconv' class object
#' @param ... retained for class compatibility
#' @param test bulk gene expression matrix assumed to be in raw counts
#' @param arith_mean logical, whether to use arithmetic mean as gene signature
#' @param use_filter logical, whether to use denoised signature matrix
#' @returns Matrix of residuals.
#' @export
residuals.deconv <- function(object,
                             test = NULL,
                             arith_mean = FALSE, use_filter = FALSE, ...) {
  if (is.null(test)) return(object$subclass$residuals)
  # recalculate residuals
  if (is.null(arith_mean)) arith_mean <- object$call$arith_mean
  if (is.null(use_filter)) use_filter <- object$call$use_filter
  if (arith_mean) {
    cellmat <- if (use_filter) object$mk$genemeans_ar_filter else object$mk$genemeans_ar
    if (is.null(cellmat)) stop("arithmetic mean not available")
  } else {
    cellmat <- if (use_filter) object$mk$genemeans_filter else object$mk$genemeans
  }
  nm <- object$call$test
  .call <- match.call()
  if (nm != .call$test) {
    message("`", .call$test, "` does not match the bulk dataset `", nm,
            "` used in `", .call$object, "`")
  }
  geneset <- intersect(rownames(cellmat), rownames(test))
  cellmat <- cellmat[geneset, ]
  test <- test[geneset, ]
  cellmat <- 2^cellmat -1
  residuals_deconv(test, cellmat, object$subclass$output)
}


#' Standard errors of deconvoluted cell counts
#' 
#' Extracts standard errors of deconvoluted cell counts based on the row
#' variance of weighted residuals.
#' 
#' @param model a 'deconv' class object
#' @returns a vector of the standard errors of cell counts for each cell
#'   subclass
#' @seealso [deconvolute()]
#' @export
se <- function(model) {
  X <- model$subclass$X
  var.e <- model$subclass$var.e
  Lv <- colSums(X^2)
  iXTX <- model$subclass$compensation / Lv
  XTXse <- crossprod(X, var.e * X)
  # var(beta) = (X' X)^-1 (X diag(e^2) X') (X' X)^-1
  sqrt(diag(iXTX %*% XTXse %*% t(iXTX)))
}

