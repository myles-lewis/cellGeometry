
dotprod <- function(test, cellmat, weights = NULL) {
  msc <- cellmat
  if (!is.null(weights)) {
    if (length(weights) != nrow(cellmat)) stop("incorrect weights length")
    msc <- cellmat * weights
    test <- test * weights
  }
  md <- colSums(msc^2)
  t( t(crossprod(test, msc)) / md )
}

# t( t(t(test) %*% msc) / md )

comp_metric <- function(m) {
  m2 <- m - diag(nrow(m))
  mean(abs(m2))
}

max_spill <- function(m) {
  m2 <- m - diag(nrow(m))
  max(m2)
}

max_abs <- function(m) {
  if (abs(min(m)) > max(m)) return(min(m))
  max(m)
}

residuals_deconv <- function(test, cellmat, output, count_space) {
  pred <- tcrossprod(cellmat, output)
  # if (count_space) {
  #   pred[pred < 0] <- 0
  #   return(log2(test +1) - log2(pred +1))
  # }
  test - pred
}
