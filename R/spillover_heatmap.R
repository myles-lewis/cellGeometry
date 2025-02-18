
#' Spillover heatmap
#' 
#' Produces a heatmap from a 'cellMarkers' or 'deconv' class object showing
#' estimated amount of spillover between cell subclasses. The amount that each
#' cell subclass's overall vector spillovers (projects) into other cell
#' subclasses' vectors is shown in each row. Thus the column gives an estimate
#' of how much the most influential (specific) genes for a cell subclass are
#' expressed in other cells.
#' 
#' @param x Either a 'cellMarkers' or 'deconv' class object or a spillover
#'   matrix.
#' @param col Vector of colours or colour mapping function passed to
#'   `ComplexHeatmap::Heatmap()`.
#' @param ... Optional arguments passed to [ComplexHeatmap::Heatmap()].
#' @returns No return value. Draws a heatmap using ComplexHeatmap.
#' @importFrom grid grid.text unit
#' @importFrom ComplexHeatmap draw
#' @importFrom circlize colorRamp2
#' @export

spillover_heatmap <- function(x,
                              col = colorRamp2(c(0, 0.5, 0.8, 1),
                                               c("#F4FAFF", "steelblue2", "purple", "red")),
                              ...) {
  cell_table <- NULL
  if (inherits(x, "deconv")) {
    cell_table <- x$mk$cell_table
    x <- x$subclass$spillover
  }
  if (inherits(x, 'cellMarkers')) {
    cell_table <- x$cell_table
    m_itself <- x$spillover
  } else m_itself <- x
  
  hm1 <- Heatmap(m_itself, col = col,
                 cluster_rows = FALSE, row_split = cell_table,
                 cluster_row_slices = FALSE, row_title = NULL,
                 cluster_columns = FALSE, column_split = cell_table,
                 cluster_column_slices = FALSE, column_title = NULL,
                 # column_names_rot = 75,
                 column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(ifelse(m_itself[i, j] > 0.5 & i!=j,
                                    sprintf("%.1f",  m_itself[i, j]), ""),
                             x, y, gp = gpar(fontsize = 8))
                 },
                 heatmap_legend_param = list(title = "spillover",
                                             legend_width = unit(6, "cm"),
                                             direction = "horizontal"),
                 ...) |> suppressMessages()
  draw(hm1, heatmap_legend_side = "top")
}
