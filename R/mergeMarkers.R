
#' Merge cellMarker signatures
#' 
#' Takes 2 cellMarkers signatures, merges them and recalculates optimal gene
#' signatures.
#' 
#' @param mk1 The reference 'cellMarkers' class object.
#' @param mk2 A 'cellMarkers' class object containing cell signatures to merge
#'   into `mk1`.
#' @param remove_subclass Optional character vector of subclasses to remove when
#'   merging.
#' @param remove_group Optional character vector of cell groups to remove when
#'   merging.
#' @param transform Either "qq" which applies [quantile_map()] to `mk2` to
#'   quantile transform it onto the same distribution as `mk1`,
#'   "linear.qq", which determines the quantile transformation and then
#'   applies a linear approximation of this, "scale" which simply scales the
#'   gene expression by the value `scale`, or "none" for no transformation.
#' @param scale Numeric value determining the scaling factor for `mk2` if
#'   `transform` is set to "scale".
#' @param ... Optional arguments and settings passed to [updateMarkers()].
#' @returns A list object of S3 class 'cellMarkers'. See [cellMarkers()] for
#'   details. If `transform = "qq"` then an additional element `qqmerge` is
#'   returned containing the quantile mapping function between the 2 datasets.
#' @seealso [cellMarkers()] [updateMarkers()] [quantile_map()]
#' @export
mergeMarkers <- function(mk1, mk2,
                         remove_subclass = NULL,
                         remove_group = NULL,
                         transform = c("qq", "linear.qq", "scale", "none"),
                         scale = 1, ...) {
  .call <- match.call()
  if (!inherits(mk1, "cellMarkers")) stop("'mk1' is not a 'cellMarkers' object")
  if (!inherits(mk2, "cellMarkers")) stop("'mk2' is not a 'cellMarkers' object")
  xlab <- deparse(substitute(mk2))
  ylab <- deparse(substitute(mk1))
  
  transform <- match.arg(transform)
  qfun <- NULL
  if (transform == "qq") {
    message("Quantile transforming '", xlab, "'")
    qfun <- quantile_map(mk2, mk1, respace = TRUE) |> suppressMessages()
    mk2$genemeans <- qfun$map(mk2$genemeans)
    mk2$groupmeans <- qfun$map(mk2$groupmeans)
    if (!is.null(mk2$genemeans_ar)) mk2$genemeans_ar <- qfun$map(mk2$genemeans_ar)
    qfun$xlab <- xlab
    qfun$ylab <- ylab
  } else if (transform %in% c("scale", "linear.qq")) {
    if (transform == "linear.qq") {
      message("Quantile transformation, linear approximation")
      scale <- linear_qq(mk2, mk1)
    }
    mk2$genemeans <- mk2$genemeans * scale
    mk2$groupmeans <- mk2$groupmeans * scale
    mk2$genemeans_ar <- mk2$genemeans_ar * scale
  }
  
  common <- intersect(rownames(mk1$genemeans), rownames(mk2$genemeans))
  message(length(common), " common genes")
  gm1 <- mk1$genemeans[common, ]
  gm2 <- mk2$genemeans[common, ]
  genemeans <- cbind(gm1, gm2)
  gp1 <- mk1$groupmeans[common, ]
  gp2 <- mk2$groupmeans[common, ]
  groupmeans <- cbind(gp1, gp2)
  gm1_ar <- mk1$genemeans_ar[common, ]
  gm2_ar <- mk2$genemeans_ar[common, ]
  genemeans_ar <- cbind(gm1_ar, gm2_ar)
  
  # check for duplicate group
  dup <- duplicated(colnames(groupmeans))
  if (any(dup)) {
    dup_nm <- colnames(groupmeans)[dup]
    message("Duplicated groups: ", paste(dup_nm, collapse = ", "))
    nm <- paste0(dup_nm, ".1")
    colnames(groupmeans)[dup] <- nm
    w <- which(levels(mk2$cell_table) %in% dup_nm)
    levels(mk2$cell_table)[w] <- paste0(levels(mk2$cell_table)[w], ".1")
  }
  
  # remove subclass or group
  cell_table <- c(mk1$cell_table, mk2$cell_table)
  rem_subcl <- colnames(genemeans) %in% remove_subclass |
    cell_table %in% remove_group
  if (any(rem_subcl)) genemeans <- genemeans[, !rem_subcl]
  cell_table <- cell_table[!rem_subcl]
  old_levels <- levels(cell_table)
  cell_table <- droplevels(cell_table)
  gone <- setdiff(old_levels, levels(cell_table))
  remove_group <- c(remove_group, gone)
  if (any(remove_group %in% colnames(groupmeans))) {
    groupmeans <- groupmeans[, !colnames(groupmeans) %in% remove_group]
  }
  subclass_table <- c(mk1$subclass_table, mk2$subclass_table)
  subclass_table <- subclass_table[!rem_subcl]
  
  # check for duplicate subclass
  dup <- duplicated(colnames(genemeans))
  if (any(dup)) {
    dup_nm <- colnames(genemeans)[dup]
    message("Duplicated subclasses: ", paste(dup_nm, collapse = ", "))
    nm <- paste0(dup_nm, ".1")
    colnames(genemeans)[dup] <- nm
    colnames(genemeans_ar)[dup] <- nm
    names(cell_table)[dup] <- nm
    names(subclass_table)[dup] <- nm
  }
  
  mk1$genemeans <- genemeans
  mk1$groupmeans <- groupmeans
  mk1$cell_table <- cell_table
  mk1$subclass_table <- subclass_table
  mk1$genemeans_ar <- genemeans_ar
  if (!is.null(qfun)) mk1$qqmerge <- qfun
  
  updateMarkers(mk1, ...)
}
