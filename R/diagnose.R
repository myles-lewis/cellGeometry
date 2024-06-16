
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
    cat(paste(paste0(format(colnames(mk$genemeans)[w], justify = "left"),
                     "  ", no1[w], "/", nsubclass),
              collapse = "\n"))
    cat("\n")
  }
  
  cat("\n")
  s1 <- mk$spillover
  spmax <- vapply(w, function(i) {
    x <- colnames(s1)[i]
    col_i <- s1[, i]
    col_i[i] <- 0
    wmax <- which.max(col_i)
    spmax <- colnames(s1)[wmax]
    c(spmax, format(max(col_i), digits = 3))
  }, character(2))
  s1w <- format(colnames(s1)[w], justify = "left")
  spmaxname <- format(spmax[1, ], justify = "left")
  s1set <- paste(s1w, " spills most into ", spmaxname, "", spmax[2, ])
  cat(paste(s1set, collapse = "\n"))
  cat("\n\n")
  
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
  
  
  smetric <- comp_metric(s1)
  cat("Standard mean spillover", format(smetric, digits = 3), "\n")
  
  m <- mk$genemeans_filtered[mk$geneset, ]
  s2 <- dotprod(m, m, equal_weight = TRUE)
  smetric <- comp_metric(s2)
  cat("Equal weight mean spillover", format(smetric, digits = 3), "\n")
  
}
