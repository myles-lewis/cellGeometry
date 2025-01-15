
#' Compensation heatmap
#' 
#' Plots a heatmap of the compensation matrix for cell subclasses using 
#' ComplexHeatmap.
#' 
#' @param x object of class 'deconv' or a matrix of compensation values.
#' @param cell_table optional grouping vector to separate the heatmap rows and
#'   columns into groups.
#' @param ... optional arguments passed to [ComplexHeatmap::Heatmap()]
#' @returns No return value. Draws a ComplexHeatmap.
#' @export
#'
comp_heatmap <- function(x, cell_table = NULL, ...) {
  if (inherits(x, "deconv")) {
    if (is.null(cell_table)) cell_table <- x$mk$cell_table
    comp <- x$subclass$compensation
  } else {
    comp <- x
  }
  
  hm2 <- Heatmap(comp,
                 cluster_rows = FALSE, row_split = cell_table,
                 cluster_row_slices = FALSE, row_title = NULL,
                 cluster_columns = FALSE, column_split = cell_table,
                 cluster_column_slices = FALSE, column_title = NULL,
                 column_names_rot = 75, column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(ifelse(abs(comp[i, j]) > 0.2, sprintf("%.1f",  comp[i, j]), ""),
                             x, y, gp = gpar(fontsize = 8))
                 },
                 heatmap_legend_param = list(title = 'compensation',
                                             legend_width = unit(6, "cm"),
                                             direction = "horizontal"),
                 ...) |> suppressMessages()
  draw(hm2, heatmap_legend_side = "top")
}
