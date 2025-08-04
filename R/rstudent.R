
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
  hat <- model$subclass$hat
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
  hat <- model$subclass$hat
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
  hat <- model$subclass$hat
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
  hat <- fit$hat
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
  hat <- fit$hat
  p <- ncol(fit$compensation)
  rdf <- nrow(r) - p
  rss <- colSums(r^2)
  mse <- rss / rdf
  (r / (1 - hat))^2 * outer(hat, mse, "/") / p
}
