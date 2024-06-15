
#' Spillover heatmap
#' 
#' Produces a heatmap from a 'cellMarkers' class object showing estimated amount
#' of spillover between cell subclasses. The amount that each cell subclass's
#' overall vector spillovers (projects) into other cell subclasses' vectors is
#' shown in each column.
#' 
#' @param x Either a cellMarkers class object or a spillover matrix.
#' @param ... Optional arguments passed to [ComplexHeatmap::Heatmap()].
#' @returns No return value. Draws a heatmap using ComplexHeatmap.
#' @importFrom grid grid.text unit
#' @importFrom ComplexHeatmap draw
#' @export

spillover_heatmap <- function(x, ...) {
  if (inherits(x, "deconv")) x <- x$mk
  if (inherits(x, 'cellMarkers')) {
    cell_table <- x$cell_table
    m_itself <- x$spillover
  } else {
    m_itself <- x
    cell_table <- NULL
  }
  hm1 <- Heatmap(m_itself,
                 cluster_rows = FALSE, row_split = cell_table,
                 cluster_row_slices = FALSE, row_title = NULL,
                 cluster_columns = FALSE, column_split = cell_table,
                 cluster_column_slices = FALSE, column_title = NULL,
                 column_names_rot = 75, column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(ifelse(m_itself[i, j] > 0.5 & i!=j,
                                    sprintf("%.1f",  m_itself[i, j]), ""),
                             x, y, gp = gpar(fontsize = 8))
                 },
                 heatmap_legend_param = list(title = "spillover",
                                             legend_width = unit(6, "cm"),
                                             direction = "horizontal"),
                 ...)
  draw(hm1, heatmap_legend_side = "top")
}
