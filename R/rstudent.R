
#' Regression Deletion Diagnostics
#' 
#' Functions for computing regression diagnostics.
#' 
#' `rstandard` and `rstudent` give standardized and Studentized residuals
#' respectively.
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
