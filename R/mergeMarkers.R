
#' @export
mergeMarkers <- function(mk1, mk2,
                         remove_subclass = NULL,
                         remove_group = NULL,
                         quantile_map = TRUE, ...) {
  .call <- match.call()
  if (!inherits(mk1, "cellMarkers")) stop("'mk1' is not a 'cellMarkers' object")
  if (!inherits(mk2, "cellMarkers")) stop("'mk2' is not a 'cellMarkers' object")
  
  if (quantile_map) {
    message("Quantile mapping")
    qfun <- quantile_map(mk2, mk1) |> suppressMessages()
    mk2$genemeans <- qfun$map(mk2$genemeans)
    mk2$groupmeans <- qfun$map(mk2$groupmeans)
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
  mk1$qmap <- qfun
  
  updateMarkers(mk1, ...)
}
