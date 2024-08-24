
#' Diagnostics for cellMarker signatures
#' 
#' Identifies cell subclasses or groups with weak signatures.
#' 
#' @param mk A 'cellMarkers' class object.
#' @returns No return value. Displays information about the cellMarkers
#'   signature showing cells subclasses with weak signatures and diagnostic
#'   information including which cell subclasses each problematic signature
#'   spills into.
#' @export

diagnose <- function(mk) {
  if (!inherits(mk, "cellMarkers")) stop ("Not a 'cellMarkers' class object")
  nsubclass <- mk$opt$nsubclass
  if (is.null(nsubclass)) nsubclass <- 5
  nsubclass <- max(nsubclass)
  ngroup <- mk$opt$ngroup
  if (is.null(ngroup)) ngroup <- 5
  ngroup <- max(ngroup)
  
  no1 <- vapply(mk$best_angle, function(i) {
    ranks <- i[seq_len(nsubclass), "rank"]
    sum(ranks == 1)
  }, numeric(1))
  w <- which(no1 < nsubclass)
  
  s1 <- mk$spillover
  if (length(w) > 0) {
    spmax <- vapply(w, function(i) {
      col_i <- s1[, i]
      col_i[i] <- 0
      wmax <- which.max(col_i)
      c(max(col_i), wmax)
    }, numeric(2))
    w <- w[order(spmax[1, ], decreasing = TRUE)]
    spmax <- spmax[, order(spmax[1, ], decreasing = TRUE), drop = FALSE]
    
    cat("Subclass signatures with fewer than", nsubclass,
        "first rank markers:\n")
    cat(paste(paste0(format(colnames(mk$genemeans)[w], justify = "left"),
                     "  ", no1[w], "/", nsubclass),
              collapse = "\n"))
    cat("\n\n")
    
    s1w <- format(colnames(s1)[w], justify = "left")
    spmaxname <- format(colnames(s1)[spmax[2, ]], justify = "left")
    s1set <- paste(s1w, " spills most into ", spmaxname, "",
                   format(spmax[1, ], digits = 3))
    cat(paste(s1set, collapse = "\n"))
    cat("\n\n")
  }
  
  no1 <- vapply(mk$group_angle, function(i) {
    ranks <- i[seq_len(ngroup), "rank"]
    sum(ranks == 1)
  }, numeric(1))
  w <- which(no1 < ngroup)
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
