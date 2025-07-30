
dotprod <- function(test, cellmat, weights = NULL) {
  md <- colSums(cellmat^2)
  t( crossprod(cellmat, test) / md )
}

equalweight <- function(cellmat) {
  vecLength <- sqrt(rowSums(cellmat^2))
  1/vecLength
} 

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

residuals_deconv <- function(test, cellmat, output) {
  pred <- tcrossprod(cellmat, output)
  test - pred
}
