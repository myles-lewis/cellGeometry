
#' Compensation heatmap
#' 
#' Plots a heatmap of the compensation matrix for cell subclasses using 
#' ComplexHeatmap.
#' 
#' @param x object of class 'deconv' or a matrix of compensation values.
#' @param cell_table optional grouping vector to separate the heatmap rows and
#'   columns into groups.
#' @param text Logical whether to show values whose absolute value > `cutoff`.
#'   By default only shown for smaller matrices.
#' @param cutoff Absolute threshold for showing values.
#' @param fontsize Numeric value for font size for cell values when
#'   `text = TRUE`.
#' @param subset Character vector of groups to be subsetted.
#' @param ... optional arguments passed to [ComplexHeatmap::Heatmap()]
#' @returns No return value. Draws a ComplexHeatmap.
#' @export
#'
comp_heatmap <- function(x,
                         cell_table = NULL,
                         text = NULL,
                         cutoff = 0.2,
                         fontsize = 8,
                         subset = NULL, ...) {
  if (inherits(x, "deconv")) {
    if (is.null(cell_table)) cell_table <- x$mk$cell_table
    comp <- x$subclass$compensation
  } else {
    comp <- x
  }
  
  if (!is.null(subset) && !is.null(cell_table)) {
    s <- which(cell_table %in% subset)
    if (length(s) == 0) stop("no such subgroup")
    if (length(s) == 1) stop("subset too small")
    return(comp_heatmap(x = comp[s, s], text = text, cutoff = cutoff,
                        fontsize = fontsize, ...))
  }
  
  if (is.null(text)) text <- dim(comp)[1] <= 50
  layer_fun <- NULL
  if (text) {
    layer_fun <- function(j, i, x, y, width, height, fill) {
      v <- pindex(comp, i, j)
      l <- abs(v) > cutoff
      if (any(l)) {
        grid.text(sprintf("%.1f", v[l]), x[l], y[l],
                  gp = gpar(fontsize = fontsize))
      }
    }
  }
  
  dots <- list(...)
  args <- list(comp,
               cluster_rows = FALSE, row_split = cell_table,
               cluster_row_slices = FALSE, row_title = NULL,
               cluster_columns = FALSE, column_split = cell_table,
               cluster_column_slices = FALSE, column_title = NULL,
               column_names_rot = 90, column_names_gp = gpar(fontsize = 8),
               row_names_gp = gpar(fontsize = 8),
               layer_fun = layer_fun,
               heatmap_legend_param = list(title = 'compensation',
                                           legend_width = unit(6, "cm"),
                                           direction = "horizontal"))
  if (length(dots)) args[names(dots)] <- dots
  hm <- do.call(Heatmap, args) |> suppressMessages()
  draw(hm, heatmap_legend_side = "top")
}
