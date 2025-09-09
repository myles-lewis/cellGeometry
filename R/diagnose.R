
#' Diagnostics for cellMarker signatures
#' 
#' Diagnostic tool which prints information for identifying cell subclasses or
#' groups with weak signatures.
#' 
#' @param object A 'cellMarkers' or 'deconv' class object.
#' @param group Character vector to focus on cell subclasses within a particular
#'   group or groups.
#' @param angle_cutoff Angle in degrees below which cell cluster vectors are
#'   considered to overlap too much. Range 0-90. See [cos_similarity()].
#' @param weak Number of 1st ranked genes for each cell cluster at which/below
#'   its gene set is considered weak.
#' @returns No return value. Prints information about the cellMarkers
#'   signature showing cells subclasses with weak signatures and diagnostic
#'   information including which cell subclasses each problematic signature
#'   spills into.
#' @export

diagnose <- function(object, group = NULL, angle_cutoff = 30, weak = 2) {
  if (inherits(object, "cellMarkers")) {
    mk <- object
    s1 <- mk$spillover
  } else if (inherits(object, "deconv")) {
    mk <- object$mk
    s1 <- object$subclass$spillover
  } else stop("incompatible object")
  if (!inherits(mk, "cellMarkers")) stop("Not a 'cellMarkers' class object")
  nsubclass <- mk$opt$nsubclass
  if (is.null(nsubclass)) nsubclass <- 5
  nsubclass <- max(nsubclass)
  ngroup <- mk$opt$ngroup
  if (is.null(ngroup)) ngroup <- 5
  ngroup <- max(ngroup)
  if (weak >= nsubclass) stop("`weak` is >= nsubclass")
  
  no1 <- vapply(mk$best_angle, function(i) {
    ranks <- i[seq_len(nsubclass), "rank"]
    sum(ranks == 1)
  }, numeric(1))
  ra <- rank_angle(cos_similarity(mk), angle_cutoff)
  
  ind <- TRUE
  if (!is.null(group)) {
    ind <- mk$cell_table %in% group
    if (sum(ind) == 0) stop("group not found")
    subcl <- colnames(mk$genemeans)[ind]
    ra_ind <- ra[, 1] %in% subcl | ra[, 2] %in% subcl
    ra <- ra[ra_ind, ]
  }
  w <- which(no1 <= weak & ind)
  w_ra <- NULL
  if (nrow(ra) > 0) {
    cat("Subclasses with vector overlap:\n")
    cat(paste(
      paste0(format(c("", as.character(ra[, 1]))), "   ",
             format(c("", as.character(ra[, 2]))), "   ",
             c("angle", format(ra$angle, digits = 3))),
      collapse = "\n"))
    cat("\n\n")
    w_ra <- unlist(lapply(ra[, 1:2], as.integer))
    w2 <- which(no1 == 0)
    s <- intersect(w2, w_ra)
    if (any(s)) {
      cat("Consider removing:\n")
      info(mk, s, no1)
    }
    s2 <- setdiff(intersect(w, w_ra), s)
    if (any(s2)) {
      cat("Consider fixing:\n")
      info(mk, s2, no1)
    }
  }
  
  if (length(w) > 0) {
    spmax <- vapply(w, function(i) {
      col_i <- s1[, i]
      col_i[i] <- 0
      wmax <- which.max(col_i)
      c(max(col_i), wmax)
    }, numeric(2))
    w <- w[order(spmax[1, ], decreasing = TRUE)]
    spmax <- spmax[, order(spmax[1, ], decreasing = TRUE), drop = FALSE]
    
    s3 <- setdiff(w, w_ra)
    if (any(s3)) {
      cat(paste0("Other weak subclass signatures (\u2264", weak,
                 " top rank markers):\n"))
      info(mk, s3, no1)
    }
    
    cat("Spillover:\n")
    s1w <- format(colnames(s1)[w], justify = "left")
    spmaxname <- format(colnames(s1)[spmax[2, ]], justify = "left")
    s1set <- paste(s1w, " spills most into ", spmaxname, "",
                   format(spmax[1, ], digits = 3))
    cat(paste(s1set, collapse = "\n"))
    cat("\n")
  }
  
  if (nrow(ra) == 0 & length(w) == 0) cat("No subclass signature issues\n")
  
  no1 <- vapply(mk$group_angle, function(i) {
    ranks <- i[seq_len(ngroup), "rank"]
    sum(ranks == 1)
  }, numeric(1))
  w <- which(no1 <= weak)
  if (length(w) > 0) {
    cat("\nWeak cell group signatures:\n")
    cat(paste(paste0(colnames(mk$groupmeans)[w], " ", no1[w], "/", ngroup),
              collapse = "\n"))
    cat("\n")
  }
  
  # smetric <- comp_metric(s1)
  # cat("Standard mean spillover", format(smetric, digits = 3), "\n")
  
  # m <- mk$genemeans_filtered[mk$geneset, ]
  # s2 <- dotprod(m, m)
  # smetric <- comp_metric(s2)
  # cat("Equal weight mean spillover", format(smetric, digits = 3), "\n")
  
}


info <- function(mk, w, no1) {
  nsubclass <- max(mk$opt$nsubclass)
  cat(paste(paste0(format(colnames(mk$genemeans)[w]), "   ",
                   format(as.vector(mk$subclass_table[w])),
                   "   ", no1[w], "/", nsubclass),
            collapse = "\n"))
  cat("\n\n")
}
