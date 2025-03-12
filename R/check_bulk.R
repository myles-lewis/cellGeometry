

bulk_stats <- function(bulkdata) {
  bulkmean <- rowMeans(log2(bulkdata +1))
  bulkvar <- matrixStats::rowVars(log2(bulkdata +1))
  posit <- rowSums(bulkdata > 1) / ncol(bulkdata)
  quant <- rank(bulkmean) / nrow(bulkdata)
  bulk_stat <- matrix(c(bulkmean, bulkvar, quant, posit),
                      nrow = nrow(bulkdata),
                      dimnames = list(rownames(bulkdata),
                                      c("mean", "var", "quantile", "prop_expr")))
}


check_bulk <- function(mk, bulk) {
  cm <- lapply(mk$best_angle, function(i) {
    genes <- rownames(i)[seq_len(mk$opt$nsubclass)]
    cormat <- cor(t(bulk[genes, ]))
    diag(cormat) <- NA
    rowMeans(cormat, na.rm = TRUE)
  })
}
