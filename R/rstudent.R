
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
#' @param type Specifies type of residuals, either raw or weighted. Ignored if
#'   `test` is specified.
#' @param test bulk gene expression matrix assumed to be in raw counts
#' @returns Matrix of residuals.
#' @export
residuals.deconv <- function(object, ...,
                             type = c("raw", "weight"),
                             test = NULL) {
  type <- match.arg(type)
  if (is.null(test)) {
    r <- object$subclass$residuals
    w <- object$subclass$weights
    if (type == "weight" && !is.null(w)) r <- r * w
    return(r)
  }
  # recalculate residuals
  arith_mean <- object$opt$arith_mean %||% FALSE
  use_filter <- object$opt$use_filter %||% TRUE
  count_space <- object$opt$count_space %||% TRUE
  cellmat <- get_cellmat(object$mk, arith_mean, use_filter, sub = TRUE)
  
  nm <- object$call$test
  .call <- match.call()
  if (nm != .call$test) {
    message("`", .call$test, "` does not match the bulk dataset `", nm,
            "` used in `", .call$object, "`")
  }
  geneset <- intersect(rownames(cellmat), rownames(test))
  cellmat <- cellmat[geneset, ]
  test <- test[geneset, ]
  if (count_space) cellmat <- 2^cellmat -1
  residuals_deconv(test, cellmat, object$subclass$output)
}


#' Standard errors of deconvoluted cell counts
#' 
#' Extracts standard errors of deconvoluted cell counts based on the row
#' variance of weighted residuals.
#' 
#' @param model a 'deconv' class object
#' @param type specifies standard error method. Default is "var.e", which we
#'   recommend as being the most accurate. The alternative options include OLS
#'   (ordinary least squares) and several heteroscedasticity-consistent methods,
#'   HC0, HC2 and HC3.
#' @returns with `type = "var.e"`, a vector of standard errors of cell counts
#'   for each cell subclass. For all other options, a matrix of SE for every
#'   bulk sample and every cell subclass.
#' @seealso [deconvolute()]
#' @export
se <- function(model, type = c("var.e", "OLS", "OLS2", "HC0", "HC2", "HC3")) {
  if (!inherits(model, "deconv")) stop("not a 'deconv' class object")
  type <- match.arg(type)
  X <- model$subclass$X
  var.e <- model$subclass$var.e
  Lv <- colSums(X^2)
  iXTX <- model$subclass$compensation / Lv
  resvar <- model$subclass$resvar
  XTXse <- crossprod(X, var.e * X)
  
  if (type == "var.e") {
    # var(beta) = (X' X)^-1 (X diag(e^2) X') (X' X)^-1
    se <- sqrt(diag(iXTX %*% XTXse %*% t(iXTX)))
    return(se)
  } else if (type == "OLS") {
    diag_XTX <- diag(model$subclass$compensation) / Lv
    se <- sqrt(resvar %*% t(diag_XTX))
  } else if (type == "OLS2") {
    # based on van Wieringen
    XTX <- crossprod(X)
    diag2 <- diag(iXTX %*% XTX %*% t(iXTX))
    se <- sqrt(resvar %*% t(diag2))
  } else {
    # heteroscedasticity-consistent SE, HC0-3
    r <- model$subclass$residuals
    weights <- model$subclass$weights
    if (!is.null(weights)) r <- r * weights
    hat <- hat(model)
    omega <- switch(type, "HC0" = function(i) i^2,
                    "HC2" = function(i) i^2 / (1 - hat),
                    "HC3" = function(i) i^2 / (1 - hat)^2)
    
    se <- t(apply(r, 2, function(i) {
      XTXse <- crossprod(X, omega(i) * X)
      sqrt(diag(iXTX %*% XTXse %*% t(iXTX)))
    }))
  }
  rownames(se) <- rownames(model$subclass$output)
  se
}


#' Compute condition number of deconvolution model
#' 
#' Computes the condition number of the spillover matrix from a deconvolution
#' model.
#' 
#' @param z a 'deconv' class object
#' @param ... arguments passed to [kappa()]
#' @returns The condition number.
#' @export
kappa.deconv <- function(z, ...) {
  kappa(z$subclass$spillover, ...)
}


#' Confidence Intervals for Deconvolution Models
#' 
#' Computes confidence intervals for fitted deconvolution models. Note that this
#' is anticipated to be most reliable when compensation values are close to 1
#' and lambda is close to 0.
#' 
#' @param object a fitted 'deconv' class model object.
#' @param parm for compatibility with S3 method. Not used.
#' @param level the confidence level required.
#' @param ... additional arguments for S3 compatibility.
#' @param type specifies standard error method, see [se()].
#' @returns List containing 2 matrices with lower and upper confidence
#'   intervals.
#' @seealso [se()]
#' @importFrom stats qnorm setNames
#' @export
confint.deconv <- function(object, parm, level = 0.95, ..., type = "var.e") {
  type <- match.arg(type, c("var.e", "OLS", "OLS2", "HC0", "HC2", "HC3"))
  se0 <- se(object, type)
  output <- object$subclass$output
  attr(output, "min") <- NULL
  if (type == "var.e") output <- t(output)
  a <- 1 - (1 - level)/2
  format_perc <- paste(format(c(1 - a, a), trim = TRUE, digits = 3,
                              scientific = FALSE), "%")
  fac <- qt(a, nrow(object$subclass$X))
  lci <- output - se0 * fac
  uci <- output + se0 * fac
  if (type == "var.e") {
    lci <- t(lci)
    uci <- t(uci)
  }
  setNames(list(lci, uci), format_perc)
}
