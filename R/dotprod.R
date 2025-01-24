
dotprod <- function(test, cellmat, equal_weight = FALSE) {
  msc <- cellmat
  if (equal_weight) {
    vecLength <- sqrt(rowSums(cellmat^2))
    msc <- cellmat / vecLength
    test <- test / vecLength
  }
  md <- colSums(msc^2)
  t( t(t(test) %*% msc) / md )
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
