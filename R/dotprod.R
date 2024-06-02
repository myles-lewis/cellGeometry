
dotprod <- function(test, cellmat, equalWeight = FALSE) {
  if (equalWeight) {
    vecLength <- sqrt(rowSums(cellmat^2))
    msc <- cellmat / vecLength
    test <- test / vecLength
  } else {
    msc <- cellmat
    geneScale <- 1
  }
  md <- colSums(msc^2)
  t( t(t(test) %*% msc) / md )
}
