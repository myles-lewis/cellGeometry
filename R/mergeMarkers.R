
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
#'   applies a linear approximation of this, or "scale" which simply scales the
#'   gene expression by the value `scale`.
#' @param scale Numeric value determining the scaling factor for `mk2` if
#'   `transform` is set to "scale".
#' @param ... Optional arguments and settings passed to [updateMarkers()].
#' @returns A list object of S3 class 'cellMarkers'. See [cellMarkers()] for
#'   details.
#' @seealso [cellMarkers()] [updateMarkers()] [quantile_map()]
#' @export
mergeMarkers <- function(mk1, mk2,
                         remove_subclass = NULL,
                         remove_group = NULL,
                         transform = c("qq", "linear.qq", "scale"),
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
    qfun <- quantile_map(mk2, mk1) |> suppressMessages()
    mk2$genemeans <- qfun$map(mk2$genemeans)
    mk2$groupmeans <- qfun$map(mk2$groupmeans)
    qfun$xlab <- xlab
    qfun$ylab <- ylab
  } else if (transform %in% c("scale", "linear.qq")) {
    if (transform == "linear.qq") {
      message("Quantile transformation, linear approximation")
      scale <- linear_qq(mk2, mk1)
    }
    mk2$genemeans <- mk2$genemeans * scale
    mk2$groupmeans <- mk2$groupmeans * scale
  }
  
  common <- intersect(rownames(mk1$genemeans), rownames(mk2$genemeans))
  message(length(common), " common genes")
  gm1 <- mk1$genemeans[common, ]
  gm2 <- mk2$genemeans[common, ]
  genemeans <- cbind(gm1, gm2)
  if (any(remove_subclass %in% colnames(genemeans))) {
    genemeans <- genemeans[, !colnames(genemeans) %in% remove_subclass]
  }
  
  cell_table <- c(mk1$cell_table, mk2$cell_table)
  cell_table <- cell_table[!names(cell_table) %in% remove_subclass]
  old_levels <- levels(cell_table)
  cell_table <- droplevels(cell_table)
  gone <- setdiff(old_levels, levels(cell_table))
  if (is.null(remove_group)) remove_group <- gone
  
  gp1 <- mk1$groupmeans[common, ]
  gp2 <- mk2$groupmeans[common, ]
  groupmeans <- cbind(gp1, gp2)
  if (any(remove_group %in% colnames(groupmeans))) {
    groupmeans <- groupmeans[, !colnames(groupmeans) %in% remove_group]
  }
  
  mk1$genemeans <- genemeans
  mk1$groupmeans <- groupmeans
  mk1$cell_table <- cell_table
  subclass_table <- c(mk1$subclass_table, mk2$subclass_table)
  subclass_table <- subclass_table[!names(subclass_table) %in% remove_subclass]
  mk1$subclass_table <- subclass_table
  if (!is.null(qfun)) mk1$qmap <- qfun
  
  updateMarkers(mk1, ...)
}
