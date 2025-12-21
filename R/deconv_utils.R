
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
