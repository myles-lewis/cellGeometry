
#' @export

diagnose <- function(mk) {
  nsubclass <- mk$call$nsubclass
  if (is.null(nsubclass)) nsubclass <- 5
  ngroup <- mk$call$ngroup
  if (is.null(ngroup)) ngroup <- 5
  
  no1 <- vapply(mk$best_angle, function(i) {
    ranks <- i[seq_len(nsubclass), "rank"]
    sum(ranks == 1)
  }, numeric(1))
  w <- which(no1 < 5)
  if (length(w) > 0) {
    cat("Weak subclass signatures:\n")
    cat(paste(paste0(colnames(mk$genemeans)[w], " ", no1[w], "/", nsubclass),
              collapse = "\n"))
    cat("\n")
  }
  
  no1 <- vapply(mk$group_angle, function(i) {
    ranks <- i[seq_len(ngroup), "rank"]
    sum(ranks == 1)
  }, numeric(1))
  w <- which(no1 < 5)
  if (length(w) > 0) {
    cat("Weak cell group signatures:\n")
    cat(paste(paste0(colnames(mk$groupmeans)[w], " ", no1[w], "/", ngroup),
              collapse = "\n"))
    cat("\n")
  }
  
  s1 <- mk$spillover
  smetric <- comp_metric(s1)
  cat("Standard mean spillover", format(smetric, digits = 3), "\n")
  
  m <- mk$genemeans_filtered[mk$geneset, ]
  s2 <- dotprod(m, m, equalWeight = TRUE)
  smetric <- comp_metric(s2)
  cat("Equal weight mean spillover", format(smetric, digits = 3), "\n")
  
}
